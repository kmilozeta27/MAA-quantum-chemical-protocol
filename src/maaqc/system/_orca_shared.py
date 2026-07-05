"""stages/_orca_shared.py — Helpers ORCA compartidos por OptStage y SpStage.

Módulo PRIVADO (prefijo `_`): agrupa lo que en pipeline.py legacy eran
funciones sueltas a nivel de módulo (`run_orca`, `check_success`,
`resolve_opt_xyz`, `assign_cpu_lists`, etc.) para no duplicarlas entre
`opt.py` y `sp.py`. Portado 1:1 desde pipeline.py — mismo comportamiento,
mismos nombres, mismos templates de ORCA.

Ubicación pendiente de revisar: cuando se migren `system/affinity.py`
(CPU topology-aware, `ResourceProfile`) y `system/platform_detect.py`
según la skill maa-repackaging, `assign_cpu_lists`/`_check_cpu_budget`
deberían mudarse ahí (son específicos del ThinkPad T14, no de ORCA en
sí). Por ahora quedan acá para no inventar esa arquitectura sin que
kmilozeta la pida explícitamente.
"""

from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import threading
import time
from pathlib import Path
from typing import Optional

from ase.io import read as ase_read

# ─────────────────────────────────────────────
#  RESOLUCIÓN DEL BINARIO ORCA
# ─────────────────────────────────────────────

# Fallback si `orca` no está en PATH (pipeline.py legacy lo hardcodeaba
# así; ver python-kmilozeta skill — puede vivir en PATH vía .bashrc, en
# cuyo caso shutil.which() lo encuentra primero).
_ORCA_FALLBACK = Path("/home/kmilozeta/orca_6_1_1_avx2/orca")


def resolve_orca_bin() -> Path:
    found = shutil.which("orca")
    return Path(found) if found else _ORCA_FALLBACK


# ─────────────────────────────────────────────
#  NOHUP / GESTIÓN DE PROCESOS ACTIVOS (Ctrl-C)
# ─────────────────────────────────────────────

USE_NOHUP: bool = True   # ver docstring de pipeline.py: resiliencia SSH

_active_procs: list[subprocess.Popen] = []
_procs_lock = threading.Lock()


def _apply_nohup(cmd: list[str]) -> list[str]:
    return ["nohup"] + cmd if USE_NOHUP else cmd


def kill_active_procs() -> None:
    """Mata todos los subprocesos ORCA/taskset activos registrados.

    Expuesto (sin `_`) para que la futura CLI lo cablee a su propio
    manejador de SIGINT — este módulo no registra `signal.signal()` por
    su cuenta, para no pisar el manejador que instale `cli.py`/`__main__.py`.
    """
    with _procs_lock:
        procs = list(_active_procs)
    for p in procs:
        try:
            p.terminate()
        except OSError:
            pass
    for p in procs:
        try:
            p.wait(timeout=3)
        except subprocess.TimeoutExpired:
            try:
                p.kill()
            except OSError:
                pass


# ─────────────────────────────────────────────
#  ENTORNO OMP/MKL POR WORKER
# ─────────────────────────────────────────────

def build_omp_env(nprocs: int) -> dict:
    return {
        **os.environ,
        "OMP_NUM_THREADS": str(nprocs),
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        # Fix CPU híbrida P-core/E-core + taskset (ver pipeline.py:203-218
        # — hwloc "intersection without inclusion" en LEANSCF sin esto).
        "OMPI_MCA_hwloc_base_binding_policy": "none",
        "OMPI_MCA_rmaps_base_mapping_policy": "core",
        "OMPI_MCA_hwloc_base_use_hwthreads_as_cpus": "1",
    }


# ─────────────────────────────────────────────
#  TEMPLATES ORCA (opt / sp)
# ─────────────────────────────────────────────

TEMPLATE_OPT = """\
! wB97X-D3 def2-SVP CPCM(water) RIJCOSX TightSCF TightOpt
%method
  Grid 4
end
%scf
  MaxIter 300
end
%maxcore {MAXCORE}
%pal
  nprocs {NPROCS}
end

* xyz {CHARGE} {MULT}
{COORDS}
*
"""

TEMPLATE_SP = """\
! wB97X-D3 def2-TZVPD def2/J CPCM(water) RIJCOSX TightSCF defgrid3
%scf
  MaxIter 300
end
%maxcore {MAXCORE}
%pal
  nprocs {NPROCS}
end

* xyzfile {CHARGE} {MULT} opt_opt.xyz
"""


def get_template(step: str, templates_dir: Optional[Path] = None) -> str:
    """Devuelve el template inline. Si existe `<templates_dir>/<step>.inp.tmpl`,
    lo prioriza (mismo mecanismo de override que pipeline.py legacy)."""
    if templates_dir is not None:
        tmpl_file = templates_dir / f"{step}.inp.tmpl"
        if tmpl_file.exists():
            return tmpl_file.read_text()
    mapping = {"opt": TEMPLATE_OPT, "sp": TEMPLATE_SP}
    return mapping[step]


# ─────────────────────────────────────────────
#  VERIFICACIÓN DE ÉXITO / GEOMETRÍA
# ─────────────────────────────────────────────

def check_success(out_file: Path) -> bool:
    if not out_file.exists():
        return False
    try:
        text = out_file.read_text(errors="replace")
    except OSError:
        return False
    return "ORCA TERMINATED NORMALLY" in text


def resolve_opt_xyz(opt_dir: Path) -> Optional[Path]:
    """ORCA 6 escribe 'opt.xyz'; versiones previas usaban 'opt_opt.xyz'.
    Normaliza siempre a 'opt_opt.xyz' (copiando si hace falta) y la
    devuelve, o None si no se encuentra ninguna de las dos."""
    canonical = opt_dir / "opt_opt.xyz"
    if canonical.exists() and canonical.stat().st_size > 0:
        return canonical
    orca6_name = opt_dir / "opt.xyz"
    if orca6_name.exists() and orca6_name.stat().st_size > 0:
        shutil.copy(orca6_name, canonical)
        return canonical
    return None


def extract_sp_energy(out_file: Path) -> Optional[float]:
    if not out_file.exists():
        return None
    try:
        for line in out_file.read_text(errors="replace").splitlines():
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    return float(line.split()[-1])
                except ValueError:
                    pass
    except OSError:
        pass
    return None


def safe_copy(src: Path, dst: Path) -> bool:
    """Copia src -> dst solo si son archivos distintos (evita SameFileError
    cuando el origen ya tiene el nombre canónico)."""
    try:
        if src.resolve() == dst.resolve():
            return True
        shutil.copy(src, dst)
        return True
    except shutil.SameFileError:
        return True
    except OSError:
        return False


def atoms_to_xyz_block(atoms) -> str:
    lines = []
    for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        lines.append(f"  {symbol:<4s}  {pos[0]:14.8f}  {pos[1]:14.8f}  {pos[2]:14.8f}")
    return "\n".join(lines)


def read_all_conformers(xyz_file: Path) -> list:
    """Lee todos los frames de un multi-XYZ (CREST/CENSO). A diferencia del
    `read_conformers(xyz_file, nconf)` legacy, no trunca — la selección de
    CUÁLES confórmeros optimizar es responsabilidad de la etapa anterior
    (crest/censo), no de OPT (ver nota de diseño en opt.py)."""
    try:
        all_conf = ase_read(str(xyz_file), index=":")
    except Exception as e:
        raise RuntimeError(f"No se pudo leer {xyz_file}: {e}") from e
    return all_conf if isinstance(all_conf, list) else [all_conf]


def extract_conf_ids_from_xyz(source_xyz: Path) -> dict[int, str]:
    """{índice_0based → 'CONF#'} leyendo los comentarios extXYZ/plain del
    multi-XYZ. Preserva la trazabilidad CONF39 → conf_39 en vez de
    re-etiquetar conf_001, conf_002, ... (ver docstring original en
    pipeline.py, función `_extract_conf_ids_from_xyz`)."""
    result: dict[int, str] = {}
    try:
        raw_lines = source_xyz.read_text(errors="replace").splitlines()
    except OSError:
        return result

    frame_idx = 0
    line_idx = 0
    n_atoms: Optional[int] = None

    while line_idx < len(raw_lines):
        stripped = raw_lines[line_idx].strip()
        if stripped.isdigit():
            nat = int(stripped)
            if n_atoms is None:
                n_atoms = nat
            if nat == n_atoms and line_idx + 1 < len(raw_lines):
                comment = raw_lines[line_idx + 1].strip()
                m = re.search(r"(CONF\d+)=", comment)
                if m:
                    result[frame_idx] = m.group(1)
                else:
                    for tok in comment.split():
                        if tok.startswith("CONF"):
                            result[frame_idx] = tok
                            break
                frame_idx += 1
                line_idx += nat + 2
                continue
        line_idx += 1

    return result


# ─────────────────────────────────────────────
#  AFINIDAD DE CPU (ThinkPad T14, Core Ultra 7 — 14 CPUs lógicas)
# ─────────────────────────────────────────────

def assign_cpu_lists(batch_size: int, nprocs_per_job: int) -> list[list[int]]:
    """Ported 1:1 desde pipeline.py — topología real confirmada con
    `lscpu -e` (P-core 0-3 con SMT en pares, E-core 4-11, LP E-core 12-13).
    Ver docstring completo en pipeline.py:3259 para el razonamiento
    (evita el crash de hwloc "intersection without inclusion" en LEANSCF
    que produce asumir rangos contiguos == cores físicos independientes)."""
    P_CORE_BLOCKS = [[0, 1], [2, 3]]
    E_CORE_SINGLES = [4, 5, 6, 7, 8, 9, 10, 11]
    LP_CORE_SINGLES = [12, 13]

    lists: list[list[int]] = []

    if nprocs_per_job == 2:
        for block in P_CORE_BLOCKS:
            if len(lists) >= batch_size:
                break
            lists.append(list(block))
        e_pool = list(E_CORE_SINGLES)
        while len(lists) < batch_size and len(e_pool) >= 2:
            lists.append([e_pool.pop(0), e_pool.pop(0)])
        lp_pool = list(LP_CORE_SINGLES)
        while len(lists) < batch_size and len(lp_pool) >= 2:
            lists.append([lp_pool.pop(0), lp_pool.pop(0)])
    elif nprocs_per_job == 1:
        pool = list(E_CORE_SINGLES) + list(LP_CORE_SINGLES) + \
            [c for block in P_CORE_BLOCKS for c in block]
        for cpu in pool:
            if len(lists) >= batch_size:
                break
            lists.append([cpu])
    else:
        full_pool = [c for block in P_CORE_BLOCKS for c in block] + \
            E_CORE_SINGLES + LP_CORE_SINGLES
        idx = 0
        for _ in range(batch_size):
            chunk = full_pool[idx: idx + nprocs_per_job]
            if len(chunk) < nprocs_per_job:
                break
            lists.append(chunk)
            idx += nprocs_per_job

    return lists


def check_cpu_budget(batch_size: int, nprocs_per_job: int, log: logging.Logger) -> None:
    requested = batch_size * nprocs_per_job
    available = os.cpu_count() or 0
    if requested > available:
        log.warning(
            f"  ⚠  el lote pide {requested} cores ({batch_size}×{nprocs_per_job}) "
            f"pero la máquina solo tiene {available} — los rangos de taskset "
            f"se solaparán o excederán los cores reales. Reduce batch_size/nprocs."
        )


# ─────────────────────────────────────────────
#  RUNNER ORCA
# ─────────────────────────────────────────────

def run_orca(
    inp_file: Path, out_file: Path, log: logging.Logger,
    timeout_hours: float = 24.0,
    env: Optional[dict] = None,
    cpu_list: Optional[list[int]] = None,
) -> bool:
    """Lanza ORCA de forma síncrona. True si terminó bien. Si el .out ya
    existe y dice TERMINATED NORMALLY, omite la ejecución (resume gratis
    a nivel de archivo, independiente del resume vía state.json).

    Ported 1:1 desde pipeline.py:1039 — mismo manejo de timeout, taskset
    vía cpu_list, nohup + start_new_session para resiliencia SSH."""
    if check_success(out_file):
        log.info(f"  ↩  ya completado: {out_file.name} — saltando")
        return True

    if not inp_file.exists():
        log.error(f"  ✗  no existe el input: {inp_file}")
        return False

    orca_bin = resolve_orca_bin()
    if not orca_bin.exists():
        log.error(f"  ✗  no se encontró el ejecutable de ORCA en: {orca_bin}")
        return False

    cmd = [str(orca_bin), str(inp_file.resolve())]
    if cpu_list:
        cpu_str = ",".join(str(c) for c in cpu_list)
        cmd = ["taskset", "-c", cpu_str] + cmd
        log.info(f"  ▶  lanzando ORCA: {inp_file.name}  (cores {cpu_str})"
                 + ("  [nohup]" if USE_NOHUP else ""))
    else:
        log.info(f"  ▶  lanzando ORCA: {inp_file.name}"
                 + ("  [nohup]" if USE_NOHUP else ""))

    t0 = time.time()
    timeout = timeout_hours * 3600
    proc_env = env if env is not None else build_omp_env(4)

    try:
        with open(out_file, "w") as fout:
            proc = subprocess.Popen(
                _apply_nohup(cmd),
                stdout=fout,
                stderr=subprocess.STDOUT,
                cwd=inp_file.parent,
                env=proc_env,
                start_new_session=USE_NOHUP,
            )

        with _procs_lock:
            _active_procs.append(proc)

        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
            log.error(f"  ✗  {inp_file.name} TIMEOUT ({timeout_hours} h)")
            return False
        finally:
            with _procs_lock:
                try:
                    _active_procs.remove(proc)
                except ValueError:
                    pass

        elapsed = (time.time() - t0) / 60
        if check_success(out_file):
            log.info(f"  ✓  {inp_file.name} terminó en {elapsed:.1f} min")
            return True
        try:
            tail = out_file.read_text(errors="replace").splitlines()[-10:]
            log.debug(f"  últimas líneas de {out_file.name}:\n" +
                      "\n".join(f"    {l}" for l in tail))
        except OSError:
            pass
        log.error(f"  ✗  {inp_file.name} FALLÓ (returncode={proc.returncode}, "
                  f"{elapsed:.1f} min) — revisa {out_file}")
        return False
    except FileNotFoundError as e:
        if cpu_list and "taskset" in str(e):
            log.error("  ✗  no se encontró 'taskset' en PATH — instálalo "
                      "(paquete util-linux) o corre sin cpu_list")
        else:
            log.error(f"  ✗  no se pudo ejecutar ORCA — verifica: {orca_bin}")
        return False
    except Exception as e:
        log.error(f"  ✗  error inesperado al lanzar ORCA: {e}")
        return False
