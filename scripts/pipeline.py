#!/usr/bin/env python3
"""
pipeline.py — Pipeline ORCA automatizado para cálculos conformacionales
Flujo: CREST → Opt → NumFreq → SP → TDDFT

Uso (modo normal):
    python pipeline.py --mol molecules/mol_A.xyz --charge -1
    python pipeline.py --mol molecules/mol_A.xyz --charge 0 --start sp
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --nconf 5 --start tddft
    python pipeline.py --status                  # ver estado de todos los runs

Uso (modo expansión conformacional, --expand-crest):
    Diseñado para casos como el estado B donde la inversión energética GFN2→DFT
    exige explorar más confórmeros más allá de los primeros N ya calculados.
    El flujo es: extrae confs CREST[skip:skip+window] → OPT en paralelo (lotes
    de 3 × 4 cores) → filtro ΔE(DFT) < threshold kcal/mol → NUMFREQ → SP → TDDFT.

    python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest
    python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest \\
        --crest-skip 3 --crest-window 17 --dft-threshold 3.0 --opt-parallel 3
    python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest \\
        --start numfreq   # salta expand-OPT, retoma desde NUMFREQ
"""

import argparse
import datetime
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from ase.io import read, write

# ─────────────────────────────────────────────
#  CONFIGURACIÓN GLOBAL
# ─────────────────────────────────────────────

ORCA_BIN   = Path("/home/kmilozeta/orca_6_1_1_avx2/orca")
CREST_BIN  = Path("crest")          # asumido en PATH; ajusta si hace falta
WORK_ROOT  = Path("runs")           # directorio de salida
TMPL_DIR   = Path("templates")
MULT       = 1                      # multiplicidad siempre 1
NPROCS          = 4
NPROCS_NUMFREQ  = 6      # NumFreq usa más cores: cada desplazamiento es un SCF
                          # independiente → escala bien. Ajustable con --nprocs-numfreq.
MAXCORE         = 3500
MAXCORE_NUMFREQ = 2300   # memoria por core para NumFreq (6 cores × 2300 MB ≈ 13.8 GB)
CREST_TIMEOUT_HOURS = 16.0          # timeout CREST (horas)

# OMP_ENV se construye bajo demanda en build_omp_env() para que siempre
# refleje el valor actualizado de NPROCS (nunca queda obsoleto si el
# usuario pasa --nprocs desde CLI).
def build_omp_env(nprocs: int) -> dict:
    return {
        **os.environ,
        "OMP_NUM_THREADS": str(nprocs),
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }

# Alias global mantenido para compatibilidad con run_orca(env=None)
OMP_ENV: dict = build_omp_env(NPROCS)

STEPS = ["crest", "opt", "numfreq", "sp", "tddft"]

# ─────────────────────────────────────────────
#  CONSTANTES FÍSICAS / MODO EXPAND
# ─────────────────────────────────────────────

HARTREE_TO_KCAL     = 627.5094740631   # factor de conversión
EXPAND_BATCH_SIZE   = 3                # OPTs en paralelo (default; override con --total-cores)
EXPAND_DFT_THRESHOLD = 3.0             # kcal/mol (default filtro post-OPT)


# ─────────────────────────────────────────────
#  LOGGING
# ─────────────────────────────────────────────

def setup_logging(mol_dir: Path) -> logging.Logger:
    """
    Crea un logger con nombre único por molécula para evitar que dos
    instancias paralelas del script compartan handlers y dupliquen líneas.
    """
    logger_name = f"pipeline.{mol_dir.name}"
    log = logging.getLogger(logger_name)

    # Evitar añadir handlers duplicados si se llama más de una vez
    if log.handlers:
        return log

    log.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s",
                            datefmt="%H:%M:%S")
    # consola
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    # archivo
    mol_dir.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(mol_dir / "pipeline.log")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    log.addHandler(ch)
    log.addHandler(fh)
    return log


# ─────────────────────────────────────────────
#  HELPERS XYZ / ASE
# ─────────────────────────────────────────────

def atoms_to_xyz_block(atoms) -> str:
    """Devuelve bloque de coordenadas estilo ORCA (* xyz ... *)."""
    lines = []
    for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        lines.append(f"  {symbol:<4s}  {pos[0]:14.8f}  {pos[1]:14.8f}  {pos[2]:14.8f}")
    return "\n".join(lines)


def read_conformers(xyz_file: Path, nconf: int) -> list:
    """
    Lee hasta `nconf` confórmeros desde un multi-XYZ de CREST.
    Si el archivo tiene menos de nconf estructuras, devuelve los que haya
    sin lanzar excepción.
    """
    try:
        all_conf = read(str(xyz_file), index=":")
        # ASE puede devolver un solo Atoms en vez de lista si solo hay uno
        if not isinstance(all_conf, list):
            all_conf = [all_conf]
    except Exception as e:
        raise RuntimeError(f"No se pudo leer {xyz_file}: {e}") from e

    if len(all_conf) < nconf:
        # Advertir pero continuar con los disponibles
        import warnings
        warnings.warn(
            f"Se pidieron {nconf} confórmeros pero CREST solo generó "
            f"{len(all_conf)}. Se usarán todos los disponibles.",
            stacklevel=2,
        )
    return all_conf[:nconf]


# ─────────────────────────────────────────────
#  TEMPLATES
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

TEMPLATE_NUMFREQ = """\
! wB97X-D3 def2-SVP CPCM(water) RIJCOSX TightSCF NumFreq
%method
  Grid 4
end
%scf
  MaxIter 300
end
%maxcore {MAXCORE_NUMFREQ}
%pal
  nprocs {NPROCS_NUMFREQ}
end

* xyzfile {CHARGE} {MULT} opt_opt.xyz
"""

TEMPLATE_NUMFREQ_RESTART = """\
! wB97X-D3 def2-SVP CPCM(water) RIJCOSX TightSCF NumFreq
%method
  Grid 4
end
%scf
  MaxIter 300
end
%maxcore {MAXCORE_NUMFREQ}
%pal
  nprocs {NPROCS_NUMFREQ}
end
%freq
  Restart true
end

* xyzfile {CHARGE} {MULT} opt_opt.xyz
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

TEMPLATE_TDDFT = """\
! wB97X-D3 def2-TZVPD def2-TZVPD/C CPCM(water) RIJCOSX TightSCF
%pal
  nprocs {NPROCS}
end
%maxcore {MAXCORE}
%scf
  MaxIter 300
  MOInp "{GBW_PATH}"
end
%tddft
  NRoots    20
  MaxDim    200
  ETol      1e-6
  RTol      1e-6
  TDA       false
  Triplets  false
  PrintLevel 3
  DoCD      true
  DoNTO     true
end

* xyzfile {CHARGE} {MULT} sp_opt.xyz
"""


def get_template(step: str) -> str:
    """Devuelve el template inline. Si existe un .tmpl en templates/, lo prioriza."""
    tmpl_file = TMPL_DIR / f"{step}.inp.tmpl"
    if tmpl_file.exists():
        return tmpl_file.read_text()
    mapping = {
        "opt":             TEMPLATE_OPT,
        "numfreq":         TEMPLATE_NUMFREQ,
        "numfreq_restart": TEMPLATE_NUMFREQ_RESTART,
        "sp":              TEMPLATE_SP,
        "tddft":           TEMPLATE_TDDFT,
    }
    return mapping[step]


# ─────────────────────────────────────────────
#  VERIFICACIÓN DE ÉXITO
# ─────────────────────────────────────────────

def check_success(out_file: Path) -> bool:
    """Devuelve True si ORCA terminó normalmente."""
    if not out_file.exists():
        return False
    try:
        text = out_file.read_text(errors="replace")
    except OSError:
        return False
    return "ORCA TERMINATED NORMALLY" in text


def check_no_imaginary(out_file: Path) -> bool:
    """Devuelve True si no hay frecuencias imaginarias (NUMFREQ)."""
    if not out_file.exists():
        return True
    try:
        text = out_file.read_text(errors="replace")
    except OSError:
        return True
    for line in text.splitlines():
        if re.match(r"\s+\d+:\s+-\d", line):
            return False
    return True


def has_numfreq_restart(freq_dir: Path) -> bool:
    """
    Devuelve True si hay archivos de gradientes parciales de una NumFreq
    interrumpida. ORCA 6 guarda cada displacement como freq_D#####.res.Gradients.
    Con eso el Restart true puede retomar sin recalcular nada.
    """
    gradients = list(freq_dir.glob("freq_D*.res.Gradients"))
    return len(gradients) > 0


def extract_sp_energy(out_file: Path) -> Optional[float]:
    """Extrae la energía del SP en Hartree."""
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
    """
    Copia src → dst solo si son archivos distintos.
    Devuelve True si la copia se realizó o no era necesaria.
    Evita SameFileError cuando CREST ya escribe el nombre canónico.
    """
    try:
        if src.resolve() == dst.resolve():
            return True
        shutil.copy(src, dst)
        return True
    except shutil.SameFileError:
        return True   # mismo archivo, no hay nada que hacer
    except OSError as e:
        return False


def resolve_opt_xyz(opt_dir: Path) -> Optional[Path]:
    """
    ORCA 6 escribe la geometría optimizada como 'opt.xyz'.
    Versiones anteriores usaban 'opt_opt.xyz'.
    Esta función busca ambos nombres y devuelve el que exista,
    normalizando siempre a 'opt_opt.xyz' mediante una copia si hace falta.
    Devuelve la ruta a opt_opt.xyz, o None si no se encuentra ninguno.
    """
    canonical = opt_dir / "opt_opt.xyz"
    if canonical.exists() and canonical.stat().st_size > 0:
        return canonical
    # ORCA 6: buscar opt.xyz y normalizar
    orca6_name = opt_dir / "opt.xyz"
    if orca6_name.exists() and orca6_name.stat().st_size > 0:
        shutil.copy(orca6_name, canonical)
        return canonical
    return None


# ─────────────────────────────────────────────
#  RUNNER ORCA
# ─────────────────────────────────────────────

def run_orca(inp_file: Path, out_file: Path, log: logging.Logger,
             timeout_hours: float = 16.0,
             env: Optional[dict] = None) -> bool:
    """
    Lanza ORCA de forma síncrona. Devuelve True si terminó bien.
    Si el .out ya existe y tiene TERMINATED NORMALLY, omite la ejecución.
    `env` permite pasar un entorno de proceso personalizado; si es None
    se construye dinámicamente desde los globales NPROCS actuales (evita
    usar el OMP_ENV construido al inicio del módulo, que puede estar obsoleto).
    """
    if check_success(out_file):
        log.info(f"  ↩  ya completado: {out_file.name} — saltando")
        return True

    if not inp_file.exists():
        log.error(f"  ✗  no existe el input: {inp_file}")
        return False

    if not ORCA_BIN.exists():
        log.error(f"  ✗  no se encontró el ejecutable de ORCA en: {ORCA_BIN}")
        return False

    log.info(f"  ▶  lanzando ORCA: {inp_file.name}")
    t0 = time.time()
    try:
        with open(out_file, "w") as fout:
            proc = subprocess.run(
                [str(ORCA_BIN), str(inp_file.resolve())],
                stdout=fout,
                stderr=subprocess.STDOUT,
                cwd=inp_file.parent,
                env=env if env is not None else build_omp_env(NPROCS),
                timeout=timeout_hours * 3600,
            )
        elapsed = (time.time() - t0) / 60
        if check_success(out_file):
            log.info(f"  ✓  {inp_file.name} terminó en {elapsed:.1f} min")
            return True
        else:
            # Extraer las últimas líneas del .out para diagnóstico rápido
            try:
                tail = out_file.read_text(errors="replace").splitlines()[-10:]
                log.debug(f"  últimas líneas de {out_file.name}:\n" +
                          "\n".join(f"    {l}" for l in tail))
            except OSError:
                pass
            log.error(f"  ✗  {inp_file.name} FALLÓ (returncode={proc.returncode}, "
                      f"{elapsed:.1f} min) — revisa {out_file}")
            return False
    except subprocess.TimeoutExpired:
        log.error(f"  ✗  {inp_file.name} TIMEOUT ({timeout_hours} h)")
        return False
    except FileNotFoundError:
        log.error(f"  ✗  no se pudo ejecutar ORCA — verifica ORCA_BIN: {ORCA_BIN}")
        return False
    except Exception as e:
        log.error(f"  ✗  error inesperado al lanzar ORCA: {e}")
        return False


# ─────────────────────────────────────────────
#  ETAPAS DEL PIPELINE
# ─────────────────────────────────────────────

def stage_crest(mol_xyz: Path, mol_dir: Path, charge: int,
                log: logging.Logger) -> Optional[Path]:
    """
    Etapa 01_crest: corre CREST sobre la molécula de entrada.
    Devuelve la ruta al crest_conformers.xyz, o None si falla.
    """
    crest_dir = mol_dir / "01_crest"
    crest_dir.mkdir(parents=True, exist_ok=True)

    input_xyz      = crest_dir / "crest_input.xyz"
    conformers_xyz = crest_dir / "crest_conformers.xyz"

    # Skip si ya existe y no está vacío
    if conformers_xyz.exists() and conformers_xyz.stat().st_size > 0:
        log.info(f"[CREST] ya existe {conformers_xyz} — saltando")
        return conformers_xyz

    safe_copy(mol_xyz, input_xyz)

    log.info(f"[CREST] lanzando sobre {mol_xyz.name}  (carga={charge})")
    cmd = [
        str(CREST_BIN),
        str(input_xyz.name),
        "--T", str(NPROCS),
        "--chrg", str(charge),
        "--uhf", "0",
        "--alpb", "water",
    ]

    returncode = None
    t0 = time.time()
    try:
        with open(crest_dir / "crest.log", "w") as flog:
            proc = subprocess.run(
                cmd,
                stdout=flog,
                stderr=subprocess.STDOUT,
                cwd=crest_dir,
                env=build_omp_env(NPROCS),
                timeout=CREST_TIMEOUT_HOURS * 3600,
            )
        returncode = proc.returncode
    except subprocess.TimeoutExpired:
        log.error(f"[CREST] TIMEOUT ({CREST_TIMEOUT_HOURS} h) — proceso terminado")
        return None
    except FileNotFoundError:
        log.error(f"[CREST] no se encontró el ejecutable de CREST: {CREST_BIN}")
        return None
    except Exception as e:
        log.error(f"[CREST] error inesperado: {e}")
        return None

    elapsed = (time.time() - t0) / 60

    # CREST puede escribir el archivo con distintos nombres
    for candidate in ["crest_conformers.xyz", "crest_best.xyz"]:
        src = crest_dir / candidate
        if src.exists() and src.stat().st_size > 0:
            if not safe_copy(src, conformers_xyz):
                log.error(f"[CREST] error al copiar {src} → {conformers_xyz}")
                return None
            break

    if conformers_xyz.exists() and conformers_xyz.stat().st_size > 0:
        try:
            n = len(read(str(conformers_xyz), index=":"))
        except Exception:
            n = "?"
        log.info(f"[CREST] {n} confórmeros generados en {elapsed:.1f} min")
        return conformers_xyz
    else:
        log.error(f"[CREST] no se encontró conformers.xyz "
                  f"(returncode={returncode}) — revisa {crest_dir / 'crest.log'}")
        return None


def stage_opt(conformers_xyz: Path, mol_dir: Path,
              charge: int, nconf: int, log: logging.Logger) -> list[Path]:
    """
    Etapa 02_opt: optimización DFT de cada confórmero.
    Devuelve lista de directorios exitosos.
    """
    log.info(f"[OPT] optimizando los {nconf} mejores confórmeros")
    try:
        conformers = read_conformers(conformers_xyz, nconf)
    except RuntimeError as e:
        log.error(f"[OPT] no se pudo leer {conformers_xyz}: {e}")
        return []

    if not conformers:
        log.error("[OPT] el archivo de confórmeros está vacío")
        return []

    log.info(f"[OPT] {len(conformers)} confórmeros disponibles")
    template = get_template("opt")
    successful = []

    for i, atoms in enumerate(conformers, start=1):
        conf_id  = f"conf_{i:03d}"
        conf_dir = mol_dir / "02_opt" / conf_id
        conf_dir.mkdir(parents=True, exist_ok=True)

        inp_file = conf_dir / "opt.inp"
        out_file = conf_dir / "opt.out"

        if not check_success(out_file):
            coords   = atoms_to_xyz_block(atoms)
            inp_text = template.format(
                CHARGE=charge, MULT=MULT,
                NPROCS=NPROCS, MAXCORE=MAXCORE,
                COORDS=coords,
            )
            inp_file.write_text(inp_text)

        ok = run_orca(inp_file, out_file, log)
        if ok:
            # ORCA 6 escribe opt.xyz; versiones anteriores opt_opt.xyz.
            # Normalizar siempre a opt_opt.xyz para las etapas siguientes.
            xyz_out = resolve_opt_xyz(conf_dir)
            if xyz_out is None:
                log.warning(f"[OPT] {conf_id}: ORCA terminó pero no se generó "
                            f"geometría optimizada (opt.xyz / opt_opt.xyz) — excluido")
            else:
                successful.append(conf_dir)
        else:
            log.warning(f"[OPT] {conf_id} falló — excluido de etapas siguientes")

    log.info(f"[OPT] {len(successful)}/{len(conformers)} confórmeros optimizados")
    return successful


def stage_numfreq(opt_dirs: list[Path], mol_dir: Path,
                  charge: int, log: logging.Logger) -> list[Path]:
    """
    Etapa 03_numfreq: frecuencias numéricas.
    Usa NPROCS_NUMFREQ cores y MAXCORE_NUMFREQ MB/core (defaults: 6 / 2300 MB).
    Configurable vía --nprocs-numfreq / --maxcore-numfreq o --cores.

    Restart automático: si en la carpeta de un conformero hay archivos
    freq_D*.res.Gradients de una corrida interrumpida, escribe el inp con
    '%freq Restart true end' para que ORCA retome sin recalcular gradientes.
    Si el .out ya tiene TERMINATED NORMALLY, salta directamente.
    """
    log.info(f"[NUMFREQ] calculando frecuencias para {len(opt_dirs)} confórmeros")
    template_fresh   = get_template("numfreq")
    template_restart = TEMPLATE_NUMFREQ_RESTART
    successful = []

    for opt_dir in opt_dirs:
        conf_id = opt_dir.name
        opt_xyz = resolve_opt_xyz(opt_dir)

        if opt_xyz is None:
            log.warning(f"[NUMFREQ] {conf_id}: no se encontró opt_opt.xyz — saltando")
            continue

        freq_dir = mol_dir / "03_numfreq" / conf_id
        freq_dir.mkdir(parents=True, exist_ok=True)

        out_file = freq_dir / "freq.out"
        inp_file = freq_dir / "freq.inp"

        if not check_success(out_file):
            # Copiar geometría si aún no está en la carpeta de freq
            if not safe_copy(opt_xyz, freq_dir / "opt_opt.xyz"):
                log.warning(f"[NUMFREQ] {conf_id}: error copiando opt_opt.xyz — saltando")
                continue

            # Detectar restart automáticamente
            if has_numfreq_restart(freq_dir):
                n_done = len(list(freq_dir.glob("freq_D*.res.Gradients")))
                log.info(f"[NUMFREQ] {conf_id}: restart detectado "
                         f"({n_done} gradientes guardados) — activando Restart true")
                inp_text = template_restart.format(
                    CHARGE=charge, MULT=MULT,
                    NPROCS_NUMFREQ=NPROCS_NUMFREQ,
                    MAXCORE_NUMFREQ=MAXCORE_NUMFREQ,
                )
            else:
                inp_text = template_fresh.format(
                    CHARGE=charge, MULT=MULT,
                    NPROCS_NUMFREQ=NPROCS_NUMFREQ,
                    MAXCORE_NUMFREQ=MAXCORE_NUMFREQ,
                )
            inp_file.write_text(inp_text)

        ok = run_orca(inp_file, out_file, log)
        if not ok:
            log.warning(f"[NUMFREQ] {conf_id} falló — excluido")
            continue

        if check_no_imaginary(out_file):
            log.info(f"[NUMFREQ] {conf_id}: sin frecuencias imaginarias ✓")
        else:
            log.warning(f"[NUMFREQ] {conf_id}: frecuencias imaginarias detectadas "
                        f"— incluido de todos modos (revisar manualmente)")
        successful.append(freq_dir)

    log.info(f"[NUMFREQ] {len(successful)}/{len(opt_dirs)} completados")
    return successful


def stage_sp(freq_dirs: list[Path], mol_dir: Path,
             charge: int, log: logging.Logger) -> list[Path]:
    """
    Etapa 04_sp: single point de alta precisión.
    Devuelve directorios sp exitosos.
    """
    log.info(f"[SP] single point para {len(freq_dirs)} confórmeros")
    template   = get_template("sp")
    successful = []

    for freq_dir in freq_dirs:
        conf_id  = freq_dir.name
        opt_xyz  = resolve_opt_xyz(freq_dir)

        if opt_xyz is None:
            log.warning(f"[SP] {conf_id}: no se encontró opt_opt.xyz ni opt.xyz — saltando")
            continue

        sp_dir = mol_dir / "04_sp" / conf_id
        sp_dir.mkdir(parents=True, exist_ok=True)

        out_file = sp_dir / "sp.out"
        inp_file = sp_dir / "sp.inp"

        # Solo copiar geometría si el cálculo no está ya completo
        if not check_success(out_file):
            if not safe_copy(opt_xyz, sp_dir / "opt_opt.xyz"):
                log.warning(f"[SP] {conf_id}: error copiando opt_opt.xyz — saltando")
                continue
            inp_text = template.format(
                CHARGE=charge, MULT=MULT,
                NPROCS=NPROCS, MAXCORE=MAXCORE,
            )
            inp_file.write_text(inp_text)

        ok = run_orca(inp_file, out_file, log)
        if ok:
            successful.append(sp_dir)
        else:
            log.warning(f"[SP] {conf_id} falló — excluido")

    log.info(f"[SP] {len(successful)}/{len(freq_dirs)} completados")
    return successful


def stage_tddft(sp_dirs: list[Path], mol_dir: Path,
                charge: int, log: logging.Logger) -> list[Path]:
    """
    Etapa 05_tddft: TD-DFT con guess orbital del SP.
    """
    log.info(f"[TDDFT] calculando espectros para {len(sp_dirs)} confórmeros")
    template   = get_template("tddft")
    successful = []

    for sp_dir in sp_dirs:
        conf_id  = sp_dir.name
        gbw_file = sp_dir / "sp.gbw"
        opt_xyz  = resolve_opt_xyz(sp_dir)

        if not gbw_file.exists():
            log.warning(f"[TDDFT] {conf_id}: no se encontró sp.gbw — saltando")
            continue
        if opt_xyz is None:
            log.warning(f"[TDDFT] {conf_id}: no se encontró opt_opt.xyz ni opt.xyz — saltando")
            continue

        tddft_dir = mol_dir / "05_tddft" / conf_id
        tddft_dir.mkdir(parents=True, exist_ok=True)

        out_file = tddft_dir / "tddft.out"
        inp_file = tddft_dir / "tddft.inp"

        # Solo copiar archivos si el cálculo no está ya completo
        if not check_success(out_file):
            if not safe_copy(opt_xyz, tddft_dir / "sp_opt.xyz"):
                log.warning(f"[TDDFT] {conf_id}: error copiando sp_opt.xyz — saltando")
                continue
            if not safe_copy(gbw_file, tddft_dir / "sp.gbw"):
                log.warning(f"[TDDFT] {conf_id}: error copiando sp.gbw — saltando")
                continue
            inp_text = template.format(
                CHARGE=charge, MULT=MULT,
                NPROCS=NPROCS, MAXCORE=MAXCORE,
                GBW_PATH="sp.gbw",
            )
            inp_file.write_text(inp_text)

        ok = run_orca(inp_file, out_file, log)
        if ok:
            successful.append(tddft_dir)
        else:
            log.warning(f"[TDDFT] {conf_id} falló")

    log.info(f"[TDDFT] {len(successful)}/{len(sp_dirs)} completados")
    return successful


# ─────────────────────────────────────────────
#  EXPANSIÓN CONFORMACIONAL (--expand-crest)
# ─────────────────────────────────────────────

def find_opt_minimum(mol_dir: Path, log: logging.Logger) -> Optional[float]:
    """
    Escanea runs/<mol>/02_opt/conf_*/opt.out para encontrar la energía DFT
    mínima entre los confórmeros ya optimizados.

    Se usa como referencia E_ref en el filtro de expansión; si no hay opts
    previos exitosos devuelve None (el mínimo se determinará sobre el lote nuevo).
    """
    opt_root = mol_dir / "02_opt"
    if not opt_root.exists():
        return None

    min_e  : Optional[float] = None
    min_id : Optional[str]   = None

    for conf_dir in sorted(opt_root.glob("conf_*")):
        out_file = conf_dir / "opt.out"
        if not check_success(out_file):
            continue
        e = extract_sp_energy(out_file)
        if e is not None and (min_e is None or e < min_e):
            min_e  = e
            min_id = conf_dir.name

    if min_e is not None:
        log.info(f"[EXPAND] Mínimo DFT de opts previos: {min_id}  "
                 f"E = {min_e:.10f} Ha")
    else:
        log.info("[EXPAND] Sin opts previos con energía válida — "
                 "se usará el mínimo dinámico del lote actual")
    return min_e


def _opt_worker(task: dict) -> dict:
    """
    Worker para un único OPT, diseñado para ejecutarse en un hilo de
    ThreadPoolExecutor.  ORCA corre como subproceso externo (libera el GIL),
    por lo que el paralelismo real es efectivo.

    task debe contener:
        conf_id, conf_dir, atoms, charge, template, nprocs, maxcore, log

    Devuelve:
        dict con conf_id, conf_dir, energy (Ha o None), success, elapsed (min)
    """
    conf_id  : str           = task["conf_id"]
    conf_dir : Path          = task["conf_dir"]
    atoms                    = task["atoms"]
    charge   : int           = task["charge"]
    template : str           = task["template"]
    nprocs   : int           = task["nprocs"]
    maxcore  : int           = task["maxcore"]
    log      : logging.Logger = task["log"]

    inp_file = conf_dir / "opt.inp"
    out_file = conf_dir / "opt.out"

    conf_dir.mkdir(parents=True, exist_ok=True)

    if not check_success(out_file):
        coords   = atoms_to_xyz_block(atoms)
        inp_text = template.format(
            CHARGE=charge, MULT=MULT,
            NPROCS=nprocs, MAXCORE=maxcore,
            COORDS=coords,
        )
        inp_file.write_text(inp_text)

    t0      = time.time()
    # Env por worker: construido con el nprocs del task, no con el global.
    # Evita race conditions al modificar os.environ desde hilos concurrentes.
    worker_env = build_omp_env(nprocs)
    success = run_orca(inp_file, out_file, log, env=worker_env)
    elapsed = (time.time() - t0) / 60

    energy: Optional[float] = None
    if success:
        if resolve_opt_xyz(conf_dir) is None:
            log.warning(f"[EXPAND] {conf_id}: ORCA OK pero sin geometría optimizada")
            success = False
        else:
            energy = extract_sp_energy(out_file)

    return {
        "conf_id":  conf_id,
        "conf_dir": conf_dir,
        "energy":   energy,
        "success":  success,
        "elapsed":  elapsed,
    }


def stage_expand_opt(
    conformers_xyz : Path,
    mol_dir        : Path,
    charge         : int,
    skip           : int,
    window         : int,
    dft_threshold  : float,
    batch_size     : int,
    nprocs_per_job : int,
    maxcore_per_job: int,
    log            : logging.Logger,
) -> list[Path]:
    """
    Expansión conformacional desde CREST con OPT en paralelo y filtro ΔE DFT.

    Flujo detallado:
      1. Lee CREST conformers en el rango [skip : skip+window] (índices 0-based).
      2. Optimiza en lotes de `batch_size` con ThreadPoolExecutor.
         Cada OPT usa `nprocs_per_job` cores  →  batch_size × nprocs_per_job cores
         en paralelo. Los lotes son secuenciales; dentro de cada lote los OPTs
         corren simultáneamente.
      3. Determina E_ref = mínimo DFT global (opts previos ∪ lote actual).
      4. Descarta confórmeros con ΔE = E - E_ref > dft_threshold kcal/mol.
      5. Escribe trazabilidad en mol_dir/expand_discarded.log.

    Args:
        conformers_xyz  : ruta a 01_crest/crest_conformers.xyz
        mol_dir         : directorio raíz runs/<mol>/
        charge          : carga formal
        skip            : confórmeros CREST ya procesados (se omiten)
        window          : cuántos confórmeros adicionales intentar
        dft_threshold   : umbral ΔE en kcal/mol (recomendado 3.0)
        batch_size      : OPTs simultáneos por lote
        nprocs_per_job  : cores por job OPT (viene de --nprocs o calculado por --cores)
        maxcore_per_job : MB/core por job OPT
        log             : logger del pipeline principal

    Returns:
        Lista de directorios 02_opt/conf_XXX/ que pasan el filtro,
        listos para continuar en stage_numfreq.
    """
    log.info("=" * 60)
    log.info(f"[EXPAND] skip={skip}  window={window}  "
             f"threshold={dft_threshold} kcal/mol  "
             f"batch={batch_size} × {nprocs_per_job} cores = "
             f"{batch_size * nprocs_per_job} cores/lote")

    # ── Leer el multi-XYZ de CREST ──────────────────────────────
    try:
        all_confs = read(str(conformers_xyz), index=":")
        if not isinstance(all_confs, list):
            all_confs = [all_confs]
    except Exception as exc:
        log.error(f"[EXPAND] No se pudo leer {conformers_xyz}: {exc}")
        return []

    total = len(all_confs)
    log.info(f"[EXPAND] Total confórmeros en CREST: {total}")

    end_idx   = min(skip + window, total)
    new_confs = all_confs[skip:end_idx]

    if not new_confs:
        log.warning(f"[EXPAND] Rango [{skip}:{end_idx}] vacío — nada que procesar")
        return []

    log.info(f"[EXPAND] Procesando CREST[{skip}:{end_idx}] → "
             f"conf_{skip+1:03d} … conf_{end_idx:03d}  "
             f"({len(new_confs)} estructuras)")

    # ── Referencia energética de opts previos ────────────────────
    e_ref = find_opt_minimum(mol_dir, log)

    # ── Armar lista de tareas ────────────────────────────────────
    template  = get_template("opt")
    tasks_all = [
        {
            "conf_id":  f"conf_{skip + rel + 1:03d}",
            "conf_dir": mol_dir / "02_opt" / f"conf_{skip + rel + 1:03d}",
            "atoms":    atoms,
            "charge":   charge,
            "template": template,
            "nprocs":   nprocs_per_job,
            "maxcore":  maxcore_per_job,
            "log":      log,
        }
        for rel, atoms in enumerate(new_confs)
    ]

    n_lotes = -(-len(tasks_all) // batch_size)   # ceiling division
    log.info(f"[EXPAND] {len(tasks_all)} OPTs en {n_lotes} lote(s) de ≤{batch_size} "
             f"(CPUs disponibles: {os.cpu_count()})")
    log.info(f"[EXPAND] Cores por OPT: {nprocs_per_job}  "
             f"Cores en paralelo por lote: {batch_size * nprocs_per_job}  "
             f"RAM por lote ≈ {batch_size * nprocs_per_job * maxcore_per_job / 1024:.1f} GB")

    # ── Lotes paralelos ──────────────────────────────────────────
    opt_results: list[dict] = []

    for b_start in range(0, len(tasks_all), batch_size):
        batch    = tasks_all[b_start : b_start + batch_size]
        lote_num = b_start // batch_size + 1
        ids      = [t["conf_id"] for t in batch]
        log.info(f"[EXPAND] ── Lote {lote_num}/{n_lotes}: {ids}")

        with ThreadPoolExecutor(max_workers=len(batch)) as pool:
            futures = {pool.submit(_opt_worker, t): t["conf_id"] for t in batch}
            for fut in as_completed(futures):
                res = fut.result()
                opt_results.append(res)
                tag = "✓" if res["success"] else "✗"
                e_s = f"{res['energy']:.10f} Ha" if res["energy"] is not None else "—"
                log.info(f"[EXPAND]   {tag} {res['conf_id']}  "
                         f"E = {e_s}  ({res['elapsed']:.1f} min)")

    # ── Actualizar E_ref con el mínimo del lote actual ───────────
    ok_new = [r for r in opt_results if r["success"] and r["energy"] is not None]

    if ok_new:
        e_min_new = min(r["energy"] for r in ok_new)
        if e_ref is None or e_min_new < e_ref:
            e_ref = e_min_new
            log.info(f"[EXPAND] Nuevo mínimo DFT global (lote actual): "
                     f"{e_ref:.10f} Ha")

    if e_ref is None:
        log.error("[EXPAND] No se pudo determinar E_ref — imposible aplicar filtro")
        return []

    # ── Filtro ΔE < dft_threshold ────────────────────────────────
    log.info(f"[EXPAND] Filtro ΔE < {dft_threshold} kcal/mol  "
             f"(E_ref = {e_ref:.10f} Ha)")

    passing  : list[Path] = []
    discarded: list[str]  = []

    for r in sorted(ok_new, key=lambda x: x["energy"]):  # type: ignore[return-value]
        dE = (r["energy"] - e_ref) * HARTREE_TO_KCAL     # type: ignore[operator]
        if dE <= dft_threshold:
            log.info(f"[EXPAND]   ✓  {r['conf_id']}  ΔE = {dE:+.4f} kcal/mol → PASA")
            passing.append(r["conf_dir"])
        else:
            log.info(f"[EXPAND]   ✗  {r['conf_id']}  ΔE = {dE:+.4f} kcal/mol → DESCARTADO")
            discarded.append(r["conf_id"])

    failed = [r["conf_id"] for r in opt_results if not r["success"]]

    log.info(f"[EXPAND] ─────────────────────────────────────────────")
    log.info(f"[EXPAND]   OPTs intentados : {len(opt_results)}")
    log.info(f"[EXPAND]   Fallidos OPT     : {len(failed)}  "
             f"{'  '.join(failed) if failed else '—'}")
    log.info(f"[EXPAND]   Pasan filtro     : {len(passing)}")
    log.info(f"[EXPAND]   Descartados      : {len(discarded)}  "
             f"{'  '.join(discarded) if discarded else '—'}")
    log.info(f"[EXPAND] ─────────────────────────────────────────────")

    # ── Trazabilidad en disco (texto legible) ─────────────────────
    discard_log = mol_dir / "expand_discarded.log"
    e_ref_str = f"{e_ref:.10f} Ha" if e_ref is not None else "N/A"
    with open(discard_log, "a") as flog:
        flog.write(
            f"\n# {datetime.datetime.now().isoformat()}"
            f"  threshold={dft_threshold} kcal/mol"
            f"  E_ref={e_ref_str}\n"
        )
        for r in sorted(opt_results, key=lambda x: x["conf_id"]):
            if r["energy"] is not None and e_ref is not None:
                dE      = (r["energy"] - e_ref) * HARTREE_TO_KCAL
                verdict = "PASA" if dE <= dft_threshold else "DESCARTADO"
            else:
                dE      = float("nan")
                verdict = "FALLO_OPT"
            flog.write(
                f"  {r['conf_id']}  E={r['energy']}  "
                f"dE={dE:.4f} kcal/mol  {verdict}\n"
            )

    # ── Trazabilidad en disco (JSON estructurado, fuente de verdad) ─
    # Esto es lo que lee la recuperación de --start posteriores a opt,
    # para no tener que adivinar ni re-derivar qué confs pasaron el
    # filtro ΔE a partir de qué carpetas quedaron físicamente en 02_opt/.
    filter_result_path = mol_dir / "expand_filter_result.json"
    filter_result = {
        "timestamp":     datetime.datetime.now().isoformat(),
        "skip":          skip,
        "window":        window,
        "dft_threshold": dft_threshold,
        "e_ref_ha":      e_ref,
        "passing":       sorted(p.name for p in passing),
        "discarded":     sorted(discarded),
        "failed_opt":    sorted(failed),
    }
    with open(filter_result_path, "w") as fjson:
        json.dump(filter_result, fjson, indent=2)
    log.info(f"[EXPAND] Resultado del filtro guardado en {filter_result_path}")

    return passing


# ─────────────────────────────────────────────
#  STATUS GLOBAL
# ─────────────────────────────────────────────

def print_status():
    """Imprime estado de todos los cálculos en runs/."""
    runs = WORK_ROOT
    if not runs.exists():
        print("No se encontró el directorio 'runs/'. Aún no hay cálculos.")
        return

    print(f"\n{'─'*60}")
    print(f"  Estado de runs en: {runs.resolve()}")
    print(f"{'─'*60}")

    for mol_dir in sorted(runs.iterdir()):
        if not mol_dir.is_dir():
            continue
        print(f"\n  📁 {mol_dir.name}")
        for step_num in ["01_crest", "02_opt", "03_numfreq", "04_sp", "05_tddft"]:
            step_dir = mol_dir / step_num
            if not step_dir.exists():
                print(f"     {step_num:12s}  —  no iniciado")
                continue

            conf_dirs = [d for d in step_dir.iterdir()
                         if d.is_dir() and d.name.startswith("conf_")]
            if conf_dirs:
                ok = sum(1 for d in conf_dirs
                         if any(check_success(f) for f in d.glob("*.out")))
                print(f"     {step_num:12s}  ✓ {ok}/{len(conf_dirs)} confs")
            else:
                xyz = step_dir / "crest_conformers.xyz"
                if xyz.exists() and xyz.stat().st_size > 0:
                    try:
                        n = len(read(str(xyz), index=":"))
                        print(f"     {step_num:12s}  ✓ {n} confórmeros")
                    except Exception:
                        print(f"     {step_num:12s}  ⚠  XYZ presente pero ilegible")
                else:
                    log_f = step_dir / "crest.log"
                    status = "en progreso?" if log_f.exists() else "no iniciado"
                    print(f"     {step_num:12s}  — {status}")

    print(f"\n{'─'*60}\n")


# ─────────────────────────────────────────────
#  PUNTO DE ENTRADA
# ─────────────────────────────────────────────

def main():
    global NPROCS, MAXCORE, NPROCS_NUMFREQ, MAXCORE_NUMFREQ, OMP_ENV
    parser = argparse.ArgumentParser(
        description="Pipeline ORCA: CREST → Opt → NumFreq → SP → TDDFT",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos (modo normal):
  # Pipeline completo, 8 cores totales (usa 8 para OPT/SP/TDDFT y NumFreq)
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --cores 8

  # Retomar desde SP (opt y freq ya hechos)
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --start sp

  # Retomar NumFreq interrumpida (detección automática de restart)
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --start numfreq
  # → si hay freq_D*.res.Gradients en la carpeta, activa Restart true solo

  # Solo TDDFT, 3 confórmeros
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --nconf 3 --start tddft

  # Ver estado de todos los runs
  python pipeline.py --status

Ejemplos (modo expansión conformacional):
  # Expande B: 12 cores totales → 3 OPTs simultáneos × 4 cores cada uno
  python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest --cores 12

  # Control explícito: 2 OPTs simultáneos × 6 cores c/u (= 12 cores)
  python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest \\
      --cores 12 --opt-parallel 2

  # Ventana y threshold personalizados
  python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest \\
      --cores 12 --crest-skip 3 --crest-window 17 --dft-threshold 3.0

  # Retomar desde NUMFREQ (expand-OPT ya hecho)
  python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest \\
      --cores 8 --start numfreq
        """,
    )
    parser.add_argument("--mol",    type=Path, help="ruta al .xyz de entrada")
    parser.add_argument("--charge", type=int, default=0,
                        help="carga formal de la molécula (default: 0)")
    parser.add_argument("--nconf",  type=int, default=3,
                        help="número de confórmeros a procesar (default: 3)")
    parser.add_argument("--start",  choices=STEPS, default="crest",
                        help="etapa desde la que iniciar (default: crest)")
    parser.add_argument("--status", action="store_true",
                        help="mostrar estado de todos los runs y salir")

    # ── Recursos de cómputo ─────────────────────────────────────────────────
    res_group = parser.add_argument_group(
        "Recursos de cómputo",
        "Usa --cores para control sencillo (recomendado). Las opciones avanzadas\n"
        "permiten ajustar cada etapa por separado y sobreescriben --cores.",
    )
    res_group.add_argument(
        "--cores", type=int, default=None, metavar="N",
        help=(
            "Cores TOTALES disponibles para el pipeline. El script los distribuye\n"
            "automáticamente:\n"
            "  • Modo normal  → todos los cores van a cada job secuencial\n"
            "                   (OPT / NumFreq / SP / TDDFT)\n"
            "  • Modo expand  → se reparten entre --opt-parallel jobs en paralelo.\n"
            "                   Con --opt-parallel P: nprocs_por_job = cores // P\n"
            "                   (el resto de cores se ignora para no saturar RAM).\n"
            "  maxcore se recalcula como (RAM_total_MB × 0.85) // cores,\n"
            "  con RAM_total_MB = 28000 por defecto (ajustar en MAX_RAM_MB).\n"
            "  Ejemplo: --cores 12 --opt-parallel 3  →  4 cores × 3 OPTs paralelos."
        ),
    )
    res_group.add_argument(
        "--nprocs", type=int, default=None, metavar="N",
        help=(
            f"(avanzado) Cores por job de OPT / SP / TDDFT. "
            f"Sobreescribe la parte OPT de --cores. Default: {NPROCS}."
        ),
    )
    res_group.add_argument(
        "--maxcore", type=int, default=None, metavar="MB",
        help=(
            f"(avanzado) MB/core para OPT / SP / TDDFT. "
            f"Default: {MAXCORE}. RAM por job = --nprocs × MB."
        ),
    )
    res_group.add_argument(
        "--nprocs-numfreq", type=int, default=None, metavar="N",
        help=(
            f"(avanzado) Cores para NumFreq. "
            f"Sobreescribe la parte NumFreq de --cores. Default: {NPROCS_NUMFREQ}."
        ),
    )
    res_group.add_argument(
        "--maxcore-numfreq", type=int, default=None, metavar="MB",
        help=f"(avanzado) MB/core para NumFreq. Default: {MAXCORE_NUMFREQ}.",
    )

    # ── Modo expansión conformacional ──────────────────────────────────────
    parser.add_argument(
        "--expand-crest", action="store_true",
        help=(
            "Modo expansión: extrae confórmeros adicionales de CREST "
            "(más allá de los ya procesados), los optimiza en paralelo "
            "y continúa los que pasan el filtro ΔE por NUMFREQ → SP → TDDFT. "
            "Requiere que CREST ya haya completado (01_crest/ debe existir)."
        ),
    )
    parser.add_argument(
        "--crest-skip", type=int, default=3, metavar="N",
        help="confórmeros CREST ya procesados a omitir (default: 3)",
    )
    parser.add_argument(
        "--crest-window", type=int, default=17, metavar="N",
        help="confórmeros adicionales a intentar (default: 17, i.e. conf_4…20)",
    )
    parser.add_argument(
        "--dft-threshold", type=float, default=EXPAND_DFT_THRESHOLD, metavar="X",
        help=f"umbral ΔE(DFT) en kcal/mol para el filtro post-OPT "
             f"(default: {EXPAND_DFT_THRESHOLD})",
    )
    parser.add_argument(
        "--opt-parallel", type=int, default=EXPAND_BATCH_SIZE, metavar="N",
        help=(
            f"OPTs simultáneos por lote en modo expand (default: {EXPAND_BATCH_SIZE}). "
            f"Con --cores C, cada OPT usará C // N cores."
        ),
    )
    args = parser.parse_args()

    # ── Resolver recursos: --cores (simple) vs flags avanzados ───────────────
    # RAM total en MB asumida (ajusta si tienes diferente; 28 GB dejan margen al OS)
    MAX_RAM_MB = 28_000

    # Partir de los defaults globales
    nprocs         = NPROCS
    maxcore        = MAXCORE
    nprocs_numfreq = NPROCS_NUMFREQ
    maxcore_nf     = MAXCORE_NUMFREQ

    if args.cores is not None:
        # Modo simple: distribuir cores automáticamente.
        # NumFreq siempre usa nprocs_numfreq (default 6), no el mismo valor
        # que OPT/SP/TDDFT — escala mejor porque cada desplazamiento es un
        # SCF independiente. Si --cores < 6, se limita a lo disponible.
        c = args.cores

        if args.expand_crest:
            # En expand: nprocs_por_job = cores // opt_parallel
            # maxcore se calcula sobre los cores totales en paralelo para
            # no saturar RAM.
            p = args.opt_parallel
            nprocs = max(1, c // p)
            total_parallel_cores = p * nprocs   # puede ser < c si c no es divisible
            maxcore    = max(500, (MAX_RAM_MB * 85 // 100) // max(total_parallel_cores, 1))
            # NumFreq en expand usa todos los cores disponibles (no se paralelizan
            # multiples NumFreq simultaneamente en este modo)
            nprocs_numfreq = min(c, nprocs_numfreq)
            maxcore_nf = max(500, (MAX_RAM_MB * 85 // 100) // max(nprocs_numfreq, 1))
        else:
            nprocs     = c
            maxcore    = max(500, (MAX_RAM_MB * 85 // 100) // max(nprocs, 1))
            # NumFreq puede usar hasta nprocs_numfreq (default 6) o los cores
            # disponibles si --cores es menor que ese valor.
            nprocs_numfreq = min(c, nprocs_numfreq)
            maxcore_nf = max(500, (MAX_RAM_MB * 85 // 100) // max(nprocs_numfreq, 1))

    # Las flags avanzadas sobreescriben --cores si se pasan explícitamente
    if args.nprocs         is not None: nprocs         = args.nprocs
    if args.maxcore        is not None: maxcore        = args.maxcore
    if args.nprocs_numfreq is not None: nprocs_numfreq = args.nprocs_numfreq
    if args.maxcore_numfreq is not None: maxcore_nf    = args.maxcore_numfreq

    # Validación básica de RAM
    ram_opt = nprocs * maxcore / 1024
    ram_nf  = nprocs_numfreq * maxcore_nf / 1024
    if ram_opt > MAX_RAM_MB / 1024:
        print(f"ADVERTENCIA: OPT/SP/TDDFT requieren ~{ram_opt:.1f} GB "
              f"(> {MAX_RAM_MB/1024:.0f} GB asumidos). Considera reducir --nprocs o --maxcore.",
              file=sys.stderr)
    if ram_nf > MAX_RAM_MB / 1024:
        print(f"ADVERTENCIA: NumFreq requiere ~{ram_nf:.1f} GB "
              f"(> {MAX_RAM_MB/1024:.0f} GB asumidos). Considera reducir --nprocs-numfreq.",
              file=sys.stderr)

    # Actualizar globales (usados en las funciones stage_*)
    NPROCS          = nprocs
    MAXCORE         = maxcore
    NPROCS_NUMFREQ  = nprocs_numfreq
    MAXCORE_NUMFREQ = maxcore_nf
    OMP_ENV         = build_omp_env(NPROCS)

    if args.status:
        print_status()
        return

    if args.mol is None:
        parser.error("--mol es obligatorio salvo que uses --status")

    mol_xyz = args.mol.resolve()
    if not mol_xyz.exists():
        sys.exit(f"ERROR: no se encontró {mol_xyz}")

    if not mol_xyz.suffix.lower() == ".xyz":
        sys.exit(f"ERROR: el archivo de entrada debe ser un .xyz: {mol_xyz}")

    # Validar binarios antes de empezar
    if not ORCA_BIN.exists():
        sys.exit(f"ERROR: no se encontró el ejecutable de ORCA en: {ORCA_BIN}\n"
                 f"       Ajusta ORCA_BIN al inicio del script.")

    mol_name = mol_xyz.stem
    mol_dir  = WORK_ROOT / mol_name
    mol_dir.mkdir(parents=True, exist_ok=True)

    log = setup_logging(mol_dir)
    log.info("=" * 60)
    log.info(f"PIPELINE ORCA  —  molécula: {mol_name}")
    log.info(f"Carga: {args.charge}   Mult: {MULT}   Confs: {args.nconf}")
    log.info(f"Inicio desde etapa: {args.start.upper()}")
    log.info(f"ORCA: {ORCA_BIN}")
    log.info(f"Recursos OPT/SP/TDDFT : {NPROCS} cores / {MAXCORE} MB per core "
             f"→ {NPROCS * MAXCORE / 1024:.1f} GB por job")
    log.info(f"Recursos NumFreq       : {NPROCS_NUMFREQ} cores / {MAXCORE_NUMFREQ} MB per core "
             f"→ {NPROCS_NUMFREQ * MAXCORE_NUMFREQ / 1024:.1f} GB por job")
    if args.expand_crest:
        nprocs_per_expand = NPROCS
        total_expand_cores = args.opt_parallel * nprocs_per_expand
        total_expand_ram   = total_expand_cores * MAXCORE / 1024
        log.info(f"Expand batch           : {args.opt_parallel} × {nprocs_per_expand} = "
                 f"{total_expand_cores} cores en paralelo / RAM estimada ≈ {total_expand_ram:.1f} GB")
    log.info("=" * 60)

    start_idx = STEPS.index(args.start)

    # ════════════════════════════════════════════════════════════
    #  MODO EXPANSIÓN CONFORMACIONAL  (--expand-crest)
    # ════════════════════════════════════════════════════════════
    if args.expand_crest:
        log.info("MODO: Expansión conformacional (--expand-crest)")
        log.info(f"skip={args.crest_skip}  window={args.crest_window}  "
                 f"threshold={args.dft_threshold} kcal/mol  "
                 f"paralelo={args.opt_parallel}")

        conformers_xyz = mol_dir / "01_crest" / "crest_conformers.xyz"
        if not conformers_xyz.exists():
            sys.exit(
                f"ERROR: no existe {conformers_xyz}\n"
                f"       --expand-crest requiere que CREST haya completado "
                f"previamente (corre primero sin --expand-crest)."
            )

        # ── expand-OPT: lotes paralelos + filtro ────────────────
        if start_idx <= STEPS.index("opt"):
            expand_opt_dirs = stage_expand_opt(
                conformers_xyz=conformers_xyz,
                mol_dir=mol_dir,
                charge=args.charge,
                skip=args.crest_skip,
                window=args.crest_window,
                dft_threshold=args.dft_threshold,
                batch_size=args.opt_parallel,
                nprocs_per_job=NPROCS,
                maxcore_per_job=MAXCORE,
                log=log,
            )
            if not expand_opt_dirs:
                log.error("[EXPAND] Ningún confórmero nuevo pasó el filtro — abortando.")
                sys.exit(1)
        else:
            # Recuperar SOLO los confs que pasaron el filtro ΔE.
            # Fuente de verdad: expand_filter_result.json, escrito por
            # stage_expand_opt() en la corrida que hizo el OPT+filtro.
            # NO se reconstruye por índice de carpeta — eso ignoraría el
            # filtro ΔE y arrastraría confs ya DESCARTADOs (bug histórico).
            opt_root = mol_dir / "02_opt"
            if not opt_root.exists():
                sys.exit(
                    f"ERROR: no existe {opt_root}\n"
                    f"       Con --start {args.start} se esperan OPTs previos "
                    f"en 02_opt/."
                )

            filter_result_path = mol_dir / "expand_filter_result.json"
            if not filter_result_path.exists():
                sys.exit(
                    f"ERROR: no existe {filter_result_path}\n"
                    f"       --expand-crest --start {args.start} necesita el "
                    f"resultado del filtro ΔE de una corrida previa de "
                    f"stage_expand_opt (etapa OPT). Sin ese archivo no hay "
                    f"forma segura de saber qué confórmeros pasaron el "
                    f"filtro vs. cuáles fueron DESCARTADOs.\n"
                    f"       Corre primero con --start opt (o sin --start) "
                    f"para regenerar {filter_result_path.name}, o, si lo "
                    f"perdiste pero conservas el log, reconstrúyelo a mano "
                    f"a partir de expand_discarded.log."
                )

            with open(filter_result_path) as fjson:
                filter_result = json.load(fjson)

            passing_names = filter_result.get("passing", [])
            if not passing_names:
                sys.exit(
                    f"[EXPAND] {filter_result_path} no registra ningún "
                    f"confórmero PASA — nada que hacer."
                )

            expand_opt_dirs = sorted(
                opt_root / name for name in passing_names
                if (opt_root / name).is_dir()
            )
            missing = [n for n in passing_names if not (opt_root / n).is_dir()]
            if missing:
                log.warning(
                    f"[EXPAND] {len(missing)} conf(s) marcados PASA en "
                    f"{filter_result_path.name} pero sin carpeta en 02_opt/: "
                    f"{'  '.join(missing)}"
                )
            if not expand_opt_dirs:
                sys.exit(
                    f"[EXPAND] Ningún directorio de 02_opt/ coincide con los "
                    f"confs PASA registrados en {filter_result_path.name}."
                )

            log.info(
                f"[EXPAND] OPT saltado — {len(expand_opt_dirs)} dir(s) "
                f"recuperados de 02_opt/ según filtro ΔE registrado "
                f"({filter_result_path.name}): "
                f"{'  '.join(d.name for d in expand_opt_dirs)}"
            )

        # ── NUMFREQ (solo nuevos confs) ──────────────────────────
        if start_idx <= STEPS.index("numfreq"):
            freq_dirs = stage_numfreq(expand_opt_dirs, mol_dir, args.charge, log)
            if not freq_dirs:
                log.error("[EXPAND] Ningún NUMFREQ completado. Abortando.")
                sys.exit(1)
        else:
            opt_names = {d.name for d in expand_opt_dirs}
            freq_root = mol_dir / "03_numfreq"
            freq_dirs = (
                sorted(d for d in freq_root.glob("conf_*") if d.name in opt_names)
                if freq_root.exists() else []
            )
            log.info(f"[EXPAND] NUMFREQ saltado — "
                     f"{len(freq_dirs)} dir(s) recuperados de 03_numfreq/")

        # ── SP ───────────────────────────────────────────────────
        if start_idx <= STEPS.index("sp"):
            sp_dirs = stage_sp(freq_dirs, mol_dir, args.charge, log)
            if not sp_dirs:
                log.error("[EXPAND] Ningún SP completado. Abortando.")
                sys.exit(1)
        else:
            freq_names = {d.name for d in freq_dirs}
            sp_root = mol_dir / "04_sp"
            sp_dirs = (
                sorted(d for d in sp_root.glob("conf_*") if d.name in freq_names)
                if sp_root.exists() else []
            )
            log.info(f"[EXPAND] SP saltado — "
                     f"{len(sp_dirs)} dir(s) recuperados de 04_sp/")

        # ── TDDFT ────────────────────────────────────────────────
        if start_idx <= STEPS.index("tddft"):
            tddft_dirs = stage_tddft(sp_dirs, mol_dir, args.charge, log)
            if not tddft_dirs:
                log.warning("[EXPAND] Ningún TDDFT completado.")
            else:
                log.info(f"[EXPAND] {len(tddft_dirs)} espectros TD-DFT calculados ✓")
        else:
            log.info("[EXPAND] TDDFT no incluido en esta ejecución")

        log.info("=" * 60)
        log.info("PIPELINE (EXPAND) FINALIZADO")
        log.info(f"Resultados en: {mol_dir.resolve()}")
        log.info("=" * 60)
        return   # no continuar al flujo normal

    # ════════════════════════════════════════════════════════════
    #  MODO NORMAL
    # ════════════════════════════════════════════════════════════

    # ── Etapa 1: CREST ──────────────────────────────────────────
    conformers_xyz = mol_dir / "01_crest" / "crest_conformers.xyz"

    if start_idx == 0:
        conformers_xyz = stage_crest(mol_xyz, mol_dir, args.charge, log)
        if conformers_xyz is None:
            log.error("CREST falló. Abortando pipeline.")
            sys.exit(1)
    else:
        if not conformers_xyz.exists():
            sys.exit(f"ERROR: se pidió start={args.start} pero no existe "
                     f"{conformers_xyz}\nCorre primero con --start crest")
        log.info("[CREST] saltado (--start indica etapa posterior)")

    # ── Etapa 2: OPT ────────────────────────────────────────────
    if start_idx <= STEPS.index("opt"):
        opt_dirs = stage_opt(conformers_xyz, mol_dir, args.charge, args.nconf, log)
        if not opt_dirs:
            log.error("Ningún confórmero optimizado correctamente. Abortando.")
            sys.exit(1)
    else:
        opt_dirs = sorted((mol_dir / "02_opt").glob("conf_*")) \
                   if (mol_dir / "02_opt").exists() else []
        log.info(f"[OPT] saltado — encontrados {len(opt_dirs)} directorios previos")

    # ── Etapa 3: NUMFREQ ────────────────────────────────────────
    if start_idx <= STEPS.index("numfreq"):
        freq_dirs = stage_numfreq(opt_dirs, mol_dir, args.charge, log)
        if not freq_dirs:
            log.error("Ningún cálculo de frecuencias completado. Abortando.")
            sys.exit(1)
    else:
        freq_dirs = sorted((mol_dir / "03_numfreq").glob("conf_*")) \
                    if (mol_dir / "03_numfreq").exists() else []
        log.info(f"[NUMFREQ] saltado — encontrados {len(freq_dirs)} directorios previos")

    # ── Etapa 4: SP ─────────────────────────────────────────────
    if start_idx <= STEPS.index("sp"):
        sp_dirs = stage_sp(freq_dirs, mol_dir, args.charge, log)
        if not sp_dirs:
            log.error("Ningún SP completado. Abortando.")
            sys.exit(1)
    else:
        sp_dirs = sorted((mol_dir / "04_sp").glob("conf_*")) \
                  if (mol_dir / "04_sp").exists() else []
        log.info(f"[SP] saltado — encontrados {len(sp_dirs)} directorios previos")

    # ── Etapa 5: TDDFT ──────────────────────────────────────────
    if start_idx <= STEPS.index("tddft"):
        tddft_dirs = stage_tddft(sp_dirs, mol_dir, args.charge, log)
        if not tddft_dirs:
            log.warning("Ningún TDDFT completado.")
        else:
            log.info(f"[TDDFT] {len(tddft_dirs)} espectros calculados ✓")
    else:
        log.info("[TDDFT] no incluido en esta ejecución")

    log.info("=" * 60)
    log.info("PIPELINE FINALIZADO")
    log.info(f"Resultados en: {mol_dir.resolve()}")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
