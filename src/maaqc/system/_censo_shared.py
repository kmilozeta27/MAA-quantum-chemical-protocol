"""stages/_censo_shared.py — Helpers CENSO compartidos, usados por CensoStage.

Portado desde pipeline.py:
    _parse_prescreening_orca_sp / _parse_prescreening_xtb_energy   (~1234)
    _censo_prescreen_diagnose_dir / _censo_prescreen_recover        (~1276)
    stage_censo_prescreen                                           (~1438)
    parse_censo_screening_out / extract_conformers_by_id            (~1882)
    stage_censo_multitemp                                           (~2019)
    parse_boltzmann_weights / select_nconf_boltzmann                (~2108)
    select_nconf_auto                                               (~2277)

Cambios de firma respecto al legacy (todo lo demás — mensajes de log,
orden de búsqueda de outputs, criterios de recovery — es 1:1):
  - `charge`/`temperature`/`nprocs_crest`/`timeout_hours` ya no son
    constantes de módulo ni parámetros sueltos: vienen de `RecipeConfig`
    (`config.charge`, `config.censo.*`, `config.resources.nprocs_crest`).
  - `NCONF_MIN`/`NCONF_MAX` de `select_nconf_auto` ahora son parámetros
    (`config.censo.ncon_min`/`ncon_max`) en vez de constantes globales
    — para que un recipe YAML pueda ajustarlos sin tocar código.
  - Reutiliza `safe_copy`/`build_omp_env`/`USE_NOHUP`/`_apply_nohup`
    de `_orca_shared` en vez de duplicarlos (son genéricos de
    subprocess, no específicos de ORCA a pesar del nombre del módulo
    — si esto se vuelve confuso al migrar más stages,
    conviene extraerlos a `system/process.py`).
"""

from __future__ import annotations

import json
import logging
import math
import os
import re
import shutil
import subprocess
import time
from pathlib import Path
from typing import Optional

from ase.io import read as ase_read
from ase.io import write as ase_write

from maaqc.system._orca_shared import (
    USE_NOHUP,
    _apply_nohup,
    build_omp_env,
    safe_copy,
)

HARTREE_TO_KCAL = 627.5094740631

# Límites absolutos de select_nconf_auto si el recipe no los especifica
# (mismo default que NCONF_MIN/NCONF_MAX en pipeline.py legacy).
DEFAULT_NCONF_MIN = 3
DEFAULT_NCONF_MAX = 20

_WINDOW_SMALL = 2.0
_WINDOW_MEDIUM = 3.0
_WINDOW_LARGE = 5.0


# ─────────────────────────────────────────────
#  RESOLUCIÓN DE BINARIO (genérico — CENSO no vive en orca_6_1_1_avx2)
# ─────────────────────────────────────────────

def resolve_bin(name: str, conda_env: str = "sciencebitch") -> Path:
    """Ported desde pipeline.py:_resolve_bin. Busca en PATH y, si no,
    en <prefix del conda_env>/bin/<name> vía `conda env list --json`,
    sin necesitar que el env esté activado."""
    found = shutil.which(name)
    if found:
        return Path(found)
    try:
        result = subprocess.run(
            ["conda", "env", "list", "--json"],
            capture_output=True, text=True, timeout=15,
        )
        if result.returncode == 0:
            envs = json.loads(result.stdout).get("envs", [])
            for env_path in envs:
                if Path(env_path).name == conda_env:
                    candidate = Path(env_path) / "bin" / name
                    if candidate.is_file() and os.access(candidate, os.X_OK):
                        return candidate
    except Exception:
        pass
    return Path(name)


# ─────────────────────────────────────────────
#  PARSEO DE PRESCREENING (OrcaProc._sp / XtbProc._sp de CENSO)
# ─────────────────────────────────────────────

def _parse_prescreening_orca_sp(outfile: Path) -> Optional[float]:
    if not outfile.exists():
        return None
    try:
        text = outfile.read_text(errors="ignore")
    except OSError:
        return None
    if "ORCA TERMINATED NORMALLY" not in text:
        return None
    for line in text.splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                return float(line.split()[-1])
            except (ValueError, IndexError):
                return None
    return None


def _parse_prescreening_xtb_energy(outfile: Path) -> Optional[float]:
    if not outfile.exists():
        return None
    try:
        text = outfile.read_text(errors="ignore")
    except OSError:
        return None
    for line in text.splitlines():
        if "| TOTAL ENERGY" in line:
            try:
                return float(line.split()[3])
            except (ValueError, IndexError):
                return None
    return None


# ─────────────────────────────────────────────
#  DIAGNÓSTICO / RECOVERY DE CARPETA CENSO PARCIAL
# ─────────────────────────────────────────────

def _censo_prescreen_diagnose_dir(censo_dir: Path, log: logging.Logger) -> dict:
    search_roots = [censo_dir]
    prescreen_subdir = censo_dir / "0_PRESCREENING"
    if prescreen_subdir.is_dir():
        search_roots.append(prescreen_subdir)

    conf_dirs = sorted(
        {
            d for root in search_roots if root.is_dir()
            for d in root.iterdir() if d.is_dir() and d.name.startswith("CONF")
        },
        key=lambda d: int(re.search(r"\d+", d.name).group())
        if re.search(r"\d+", d.name) else 0,
    )
    n_total = len(conf_dirs)
    n_sp_ok = 0
    n_gsolv_ok = 0
    n_full = 0

    for cd in conf_dirs:
        sp_ok = _parse_prescreening_orca_sp(cd / "sp" / "sp.out") is not None
        gas_ok = _parse_prescreening_xtb_energy(cd / "xtb_gsolv" / "gas.out") is not None
        sol_ok = _parse_prescreening_xtb_energy(cd / "xtb_gsolv" / "solv.out") is not None
        if sp_ok:
            n_sp_ok += 1
        if gas_ok and sol_ok:
            n_gsolv_ok += 1
        if sp_ok and gas_ok and sol_ok:
            n_full += 1

    can_recover = (n_full == n_total > 0)
    log.debug(
        f"[CENSO-RECOVER]  diagnóstico {censo_dir.name}: "
        f"{n_total} CONF dirs | sp_ok={n_sp_ok} | gsolv_ok={n_gsolv_ok} | "
        f"full={n_full} | can_recover={can_recover}"
    )
    return {
        "conf_dirs": conf_dirs,
        "n_total": n_total,
        "n_sp_ok": n_sp_ok,
        "n_gsolv_ok": n_gsolv_ok,
        "n_full": n_full,
        "can_recover": can_recover,
    }


def _censo_prescreen_recover(censo_dir: Path, log: logging.Logger) -> Optional[Path]:
    diag = _censo_prescreen_diagnose_dir(censo_dir, log)
    if diag["n_total"] == 0:
        log.debug("[CENSO-RECOVER]  no hay carpetas CONF* — nada que recuperar")
        return None

    data: dict[str, dict] = {}
    missing: list[str] = []

    for cd in diag["conf_dirs"]:
        name = cd.name
        energy = _parse_prescreening_orca_sp(cd / "sp" / "sp.out")
        if energy is None:
            missing.append(name)
            continue
        energy_gas = _parse_prescreening_xtb_energy(cd / "xtb_gsolv" / "gas.out")
        energy_solv = _parse_prescreening_xtb_energy(cd / "xtb_gsolv" / "solv.out")
        if energy_gas is None or energy_solv is None:
            missing.append(name)
            continue
        gsolv = energy_solv - energy_gas
        grrho = 0.0
        gtot = energy + gsolv + grrho
        data[name] = {"energy": energy, "gsolv": gsolv, "grrho": grrho, "gtot": gtot}

    if missing:
        log.warning(
            f"[CENSO-RECOVER]  {len(missing)}/{diag['n_total']} confórmeros "
            f"sin datos completos: {', '.join(missing[:10])}"
            f"{'…' if len(missing) > 10 else ''}"
        )
        return None
    if not data:
        log.warning("[CENSO-RECOVER]  ningún confórmero recuperable")
        return None

    recovered_json = censo_dir / "0_PRESCREENING_recovered.json"
    dump = {"part_name": "prescreening", "data": data}
    try:
        recovered_json.write_text(json.dumps(dump, indent=4))
    except OSError as e:
        log.error(f"[CENSO-RECOVER]  no se pudo escribir {recovered_json}: {e}")
        return None

    log.info(f"[CENSO-RECOVER]  ✓  JSON reconstruido → {recovered_json.name} ({len(data)} confórmeros)")
    return recovered_json


# ─────────────────────────────────────────────
#  PARSEO DE 1_SCREENING.out
# ─────────────────────────────────────────────

def parse_censo_screening_out(screening_out: Path) -> Optional[list[str]]:
    """Lista de CONF# ordenados por Gtot (más estable primero), o None."""
    if not screening_out.exists():
        return None
    conf_ids: list[tuple[float, str]] = []
    try:
        lines = screening_out.read_text(errors="replace").splitlines()
    except Exception:
        return None
    for line in lines:
        line = line.strip()
        m = re.match(
            r'^(CONF\d+)\s+'
            r'[\-\d.]+\s+'
            r'[\-\d.]+\s+'
            r'[\-\d.]+\s+'
            r'([\-\d.]+)\s+'
            r'([\-\d.]+)\s+'
            r'[\d.]+',
            line,
        )
        if m:
            conf_ids.append((float(m.group(2)), m.group(1)))
    if not conf_ids:
        return None
    conf_ids.sort(key=lambda x: x[0])
    return [cid for _, cid in conf_ids]


def parse_boltzmann_weights(screening_out: Path) -> dict[str, float]:
    if not screening_out.exists():
        return {}
    result: dict[str, float] = {}
    try:
        for line in screening_out.read_text(errors="replace").splitlines():
            line = line.strip()
            m = re.match(
                r'^(CONF\d+)\s+'
                r'[\-\d.]+\s+'
                r'[\-\d.]+\s+'
                r'[\-\d.]+\s+'
                r'[\-\d.]+\s+'
                r'[\-\d.]+\s+'
                r'([\d.]+)',
                line,
            )
            if m:
                result[m.group(1)] = float(m.group(2))
    except Exception:
        return {}
    return result


def extract_conformers_by_id(
    source_xyz: Path, conf_ids: list[str], output_xyz: Path, log: logging.Logger,
) -> bool:
    try:
        all_frames = ase_read(str(source_xyz), index=":")
        if not isinstance(all_frames, list):
            all_frames = [all_frames]
    except Exception as e:
        log.warning(f"[CENSO] no se pudo leer {source_xyz}: {e}")
        return False

    frame_index: dict[str, object] = {}
    try:
        raw_lines = source_xyz.read_text(errors="replace").splitlines()
        n_atoms = len(all_frames[0]) if all_frames else 0
        fi, li = 0, 0
        while li < len(raw_lines) and fi < len(all_frames):
            line = raw_lines[li].strip()
            if line.isdigit() and int(line) == n_atoms:
                comment = raw_lines[li + 1].strip() if li + 1 < len(raw_lines) else ""
                token = comment.split()[0] if comment.split() else ""
                if token.startswith("CONF"):
                    frame_index[token] = all_frames[fi]
                fi += 1
                li += n_atoms + 2
            else:
                li += 1
    except Exception:
        for atoms in all_frames:
            comment = atoms.info.get("comment", "") or atoms.info.get("description", "")
            token = str(comment).split()[0] if comment else ""
            if token.startswith("CONF"):
                frame_index[token] = atoms

    selected = []
    missing = []
    for cid in conf_ids:
        if cid in frame_index:
            selected.append(frame_index[cid])
        else:
            missing.append(cid)

    if missing:
        log.warning(
            f"[CENSO]  ⚠  {len(missing)} conf(s) del screening no encontrados "
            f"en {source_xyz.name}: {missing[:5]}{'…' if len(missing) > 5 else ''}"
        )
    if not selected:
        return False
    try:
        ase_write(str(output_xyz), selected)
        return True
    except Exception as e:
        log.warning(f"[CENSO] error escribiendo ensemble filtrado: {e}")
        return False


# ─────────────────────────────────────────────
#  SELECCIÓN BOLTZMANN / AUTO (post-CENSO, pre-OPT)
# ─────────────────────────────────────────────

def select_nconf_boltzmann(
    conformers_xyz: Path, screening_out: Path, log: logging.Logger,
    threshold: float = 5.0, nconf_min: int = DEFAULT_NCONF_MIN,
) -> tuple[int, list[str]]:
    weights = parse_boltzmann_weights(screening_out)
    if not weights:
        log.warning(f"[NCONF-BOLTZMANN] no se pudo leer {screening_out.name} — cayendo en select_nconf_auto")
        return 0, []

    xyz_conf_ids: list[str] = []
    try:
        lines = conformers_xyz.read_text(errors="replace").splitlines()
        n_atoms = None
        li = 0
        while li < len(lines):
            l = lines[li].strip()
            if l.isdigit():
                nat = int(l)
                if n_atoms is None:
                    n_atoms = nat
                if nat == n_atoms and li + 1 < len(lines):
                    comment = lines[li + 1].strip()
                    m = re.search(r'(CONF\d+)=', comment)
                    if m:
                        xyz_conf_ids.append(m.group(1))
                    else:
                        for tok in comment.split():
                            if tok.startswith("CONF"):
                                xyz_conf_ids.append(tok)
                                break
                    li += nat + 2
                    continue
            li += 1
    except Exception as e:
        log.warning(f"[NCONF-BOLTZMANN] no se pudo leer {conformers_xyz.name}: {e} — cayendo en select_nconf_auto")
        return 0, []

    if not xyz_conf_ids:
        log.warning(f"[NCONF-BOLTZMANN] no se encontraron tokens CONF# en {conformers_xyz.name} — cayendo en select_nconf_auto")
        return 0, []

    above = [(cid, weights[cid]) for cid in xyz_conf_ids if weights.get(cid, 0.0) >= threshold]
    below = [(cid, weights[cid]) for cid in xyz_conf_ids if 0 < weights.get(cid, 0.0) < threshold]
    no_wt = [cid for cid in xyz_conf_ids if cid not in weights]
    above.sort(key=lambda x: -x[1])

    selected = [cid for cid, _ in above]
    if len(selected) < nconf_min:
        below.sort(key=lambda x: -x[1])
        extras = [cid for cid, _ in below[: nconf_min - len(selected)]]
        selected += extras
        log.warning(
            f"[NCONF-BOLTZMANN] solo {len(above)} conf(s) ≥ {threshold}% — "
            f"añadiendo {len(extras)} conf(s) extra para garantizar nconf_min={nconf_min}: {extras}"
        )

    nconf = len(selected)
    above_str = ", ".join(f"{cid}({weights[cid]:.1f}%)" for cid, _ in above)
    log.info(f"[NCONF-BOLTZMANN] umbral={threshold}% → {len(above)} conf(s) seleccionados: {above_str}")
    if no_wt:
        log.warning(f"[NCONF-BOLTZMANN] {len(no_wt)} conf(s) sin peso en {screening_out.name}: {no_wt} — descartados")
    log.info(f"[NCONF-BOLTZMANN] → nconf={nconf} (nconf_min={nconf_min})")
    return nconf, selected


def select_nconf_auto(
    conformers_xyz: Path, log: logging.Logger,
    nconf_min: int = DEFAULT_NCONF_MIN, nconf_max: int = DEFAULT_NCONF_MAX,
) -> int:
    try:
        all_conf = ase_read(str(conformers_xyz), index=":")
        if not isinstance(all_conf, list):
            all_conf = [all_conf]
        n_avail = len(all_conf)
    except Exception as e:
        log.warning(f"[NCONF-AUTO] no se pudo leer {conformers_xyz}: {e} — usando nconf={nconf_min}")
        return nconf_min

    if n_avail == 0:
        log.warning(f"[NCONF-AUTO] ensemble vacío — usando nconf={nconf_min}")
        return nconf_min

    if n_avail <= nconf_max:
        nconf = max(nconf_min, n_avail)
        log.info(f"[NCONF-AUTO] {n_avail} confórmeros disponibles ≤ {nconf_max} — todos pasan a OPT → nconf={nconf}")
        return nconf

    energies_hartree: list[float] = []
    try:
        raw_lines = conformers_xyz.read_text(errors="replace").splitlines()
        n_at = len(all_conf[0]) if all_conf else 0
        fi, li = 0, 0
        while li < len(raw_lines) and fi < len(all_conf):
            l = raw_lines[li].strip()
            if l.isdigit() and int(l) == n_at:
                comment = raw_lines[li + 1].strip() if li + 1 < len(raw_lines) else ""
                for tok in comment.split():
                    try:
                        energies_hartree.append(float(tok))
                        break
                    except ValueError:
                        continue
                fi += 1
                li += n_at + 2
            else:
                li += 1
    except Exception:
        for atoms in all_conf:
            e = atoms.info.get("energy") or atoms.info.get("Energy")
            if e is None:
                comment = atoms.info.get("comment", "") or atoms.info.get("description", "")
                for token in str(comment).split():
                    try:
                        e = float(token)
                        break
                    except ValueError:
                        continue
            if e is not None:
                try:
                    energies_hartree.append(float(e))
                except (TypeError, ValueError):
                    pass

    if len(energies_hartree) >= 2:
        e_min, e_max = min(energies_hartree), max(energies_hartree)
        spread_kcal = (e_max - e_min) * HARTREE_TO_KCAL
        if spread_kcal <= _WINDOW_SMALL:
            nconf, tier = 5, f"spread={spread_kcal:.2f} kcal/mol ≤ {_WINDOW_SMALL} (paisaje plano)"
        elif spread_kcal <= _WINDOW_MEDIUM:
            nconf, tier = 10, f"spread={spread_kcal:.2f} kcal/mol ≤ {_WINDOW_MEDIUM} (moderado)"
        elif spread_kcal <= _WINDOW_LARGE:
            nconf, tier = 15, f"spread={spread_kcal:.2f} kcal/mol ≤ {_WINDOW_LARGE} (amplio)"
        else:
            nconf, tier = 20, f"spread={spread_kcal:.2f} kcal/mol > {_WINDOW_LARGE} (rugoso)"
        log.info(f"[NCONF-AUTO] {n_avail} confórmeros disponibles — {tier} → nconf={nconf}")
    else:
        if n_avail <= nconf_max * 2:
            nconf = min(n_avail, nconf_max)
            log.warning(
                f"[NCONF-AUTO] {n_avail} confórmeros disponibles — energías no legibles "
                f"(ensemble post-CENSO); tomando nconf={nconf} (min(n_avail, nconf_max))"
            )
        else:
            ratio = math.log(n_avail) / math.log(100)
            nconf = int(nconf_min + (nconf_max - nconf_min) * min(ratio, 1.0))
            log.warning(f"[NCONF-AUTO] {n_avail} confórmeros disponibles — energías no legibles; heurística de tamaño → nconf={nconf}")

    return max(nconf_min, min(nconf_max, nconf))


# ─────────────────────────────────────────────
#  EJECUCIÓN DE CENSO (prescreening + screening, fail-open)
# ─────────────────────────────────────────────

def run_censo_prescreen(
    conformers_xyz: Path, censo_dir: Path, charge: int, unpaired: int,
    temperature: float, nprocs_crest: int, timeout_hours: float,
    log: logging.Logger,
) -> Path:
    """Ported desde pipeline.py::stage_censo_prescreen (~1438-1875).

    Fail-open: cualquier fallo (CENSO ausente, timeout, sin ensemble
    reconocible) devuelve `conformers_xyz` sin tocar — el pipeline sigue
    con el ensemble CREST crudo, igual que el legacy.
    """
    censo_dir.mkdir(parents=True, exist_ok=True)
    output_xyz = censo_dir / "censo_ensemble_sorted.xyz"

    screening_candidates = [
        censo_dir / "1_SCREENING.out",
        censo_dir / "CENSO" / "1_SCREENING.out",
        censo_dir / "prescreening" / "1_SCREENING.out",
    ]
    has_screening = any(parse_censo_screening_out(p) for p in screening_candidates)

    if output_xyz.exists() and output_xyz.stat().st_size > 0 and has_screening:
        log.info(f"[CENSO]  ↩  screening ya completado — saltando ({output_xyz.name})")
        return output_xyz

    if output_xyz.exists() and output_xyz.stat().st_size > 0 and not has_screening:
        log.warning(
            f"[CENSO]  ⚠  {output_xyz.name} existe pero proviene solo de prescreening "
            f"— se considera incompleto. Relanzando CENSO con --screening."
        )

    # ── Recovery desde carpeta parcial ──────────────────────────────────
    diag = _censo_prescreen_diagnose_dir(censo_dir, log)
    if diag["n_total"] > 0:
        log.info(
            f"[CENSO]  ⚙  carpeta previa con {diag['n_total']} CONF dirs "
            f"(sp_ok={diag['n_sp_ok']}, gsolv_ok={diag['n_gsolv_ok']}, completos={diag['n_full']})"
        )
        prev_json = censo_dir / "0_PRESCREENING_recovered.json"
        if prev_json.exists() and prev_json.stat().st_size > 0:
            recovered_json: Optional[Path] = prev_json
        elif diag["can_recover"]:
            recovered_json = _censo_prescreen_recover(censo_dir, log)
        elif diag["n_full"] > 0:
            n_incomplete = diag["n_total"] - diag["n_full"]
            log.warning(f"[CENSO]  ⚙  {diag['n_full']}/{diag['n_total']} completas — eliminando {n_incomplete} incompleta(s)")
            for cd in diag["conf_dirs"]:
                sp_ok = _parse_prescreening_orca_sp(cd / "sp" / "sp.out") is not None
                gas_ok = _parse_prescreening_xtb_energy(cd / "xtb_gsolv" / "gas.out") is not None
                sol_ok = _parse_prescreening_xtb_energy(cd / "xtb_gsolv" / "solv.out") is not None
                if not (sp_ok and gas_ok and sol_ok):
                    try:
                        shutil.rmtree(cd)
                    except OSError as e:
                        log.warning(f"[CENSO]  no se pudo eliminar {cd.name}: {e}")
            diag = _censo_prescreen_diagnose_dir(censo_dir, log)
            recovered_json = _censo_prescreen_recover(censo_dir, log) if diag["can_recover"] else None
        else:
            recovered_json = None

        if recovered_json is not None:
            local_input = censo_dir / conformers_xyz.name
            if safe_copy(conformers_xyz, local_input):
                censo_bin = resolve_bin("censo")
                cmd_reload = [
                    str(censo_bin), "--prescreening", "--screening",
                    "-i", local_input.name, "-c", str(charge), "-u", str(unpaired),
                    "--solvent", "water", "--maxcores", str(nprocs_crest),
                    "-T", f"{temperature:.2f}", "--reload", str(recovered_json),
                ]
                log.info(f"[CENSO-RECOVER]  ▶  relanzando CENSO con --reload ({recovered_json.name})")
                try:
                    with open(censo_dir / "censo_recover.log", "w") as flog:
                        proc = subprocess.run(
                            _apply_nohup(cmd_reload), stdout=flog, stderr=subprocess.STDOUT,
                            cwd=censo_dir, env=build_omp_env(nprocs_crest),
                            timeout=timeout_hours * 3600, start_new_session=USE_NOHUP,
                        )
                    for name in ("censo_ensemble_sorted.xyz", "censo_ensemble.xyz",
                                 "crest_conformers_sorted.xyz", "crest_conformers.xyz.sorted"):
                        p = censo_dir / name
                        if p.exists() and p.stat().st_size > 0 and p != local_input:
                            if safe_copy(p, output_xyz):
                                log.info(f"[CENSO-RECOVER]  ✓  ensemble recuperado con --reload → {output_xyz.name}")
                                return output_xyz
                            break
                    log.warning(f"[CENSO-RECOVER]  ⚠  --reload no produjo ensemble (returncode={proc.returncode})")
                except subprocess.TimeoutExpired:
                    log.error("[CENSO-RECOVER]  ✗  TIMEOUT en --reload")
                except Exception as e:
                    log.error(f"[CENSO-RECOVER]  ✗  error en --reload: {e}")

        log.warning(f"[CENSO]  🗑  eliminando {len(diag['conf_dirs'])} subcarpetas CONF* para relanzar limpio")
        for cd in diag["conf_dirs"]:
            try:
                shutil.rmtree(cd)
            except OSError as e:
                log.warning(f"[CENSO]  no se pudo eliminar {cd.name}: {e}")

    # ── Recovery desde 1_SCREENING.out existente (sin relanzar CENSO) ──
    for sc in screening_candidates:
        conf_ids = parse_censo_screening_out(sc)
        if conf_ids:
            log.info(f"[CENSO]  📄  {sc.relative_to(censo_dir)} con {len(conf_ids)} conf(s) — reconstruyendo sin relanzar CENSO")
            geom_sources = [conformers_xyz]
            for cand in (sc.parent / "1_SCREENING.xyz", sc.parent.parent / "1_SCREENING.xyz"):
                if cand.exists() and cand.stat().st_size > 0:
                    geom_sources.append(cand)
            for geom_src in geom_sources:
                if extract_conformers_by_id(geom_src, conf_ids, output_xyz, log):
                    log.info(f"[CENSO]  ✓  ensemble reconstruido desde {sc.name} → {output_xyz.name}")
                    return output_xyz
            log.warning(f"[CENSO]  ⚠  no se pudo reconstruir desde {sc.name} — relanzando CENSO")
            break

    # ── Lanzamiento limpio ───────────────────────────────────────────────
    local_input = censo_dir / conformers_xyz.name
    if not safe_copy(conformers_xyz, local_input):
        log.warning(f"[CENSO] no se pudo copiar {conformers_xyz.name} — usando ensemble CREST sin filtrar")
        return conformers_xyz

    censo_bin = resolve_bin("censo")
    cmd = [
        str(censo_bin), "--prescreening", "--screening",
        "-i", local_input.name, "-c", str(charge), "-u", str(unpaired),
        "--solvent", "water", "--maxcores", str(nprocs_crest),
        "-T", f"{temperature:.2f}",
    ]
    log.info(f"[CENSO]  ▶  lanzando pre-screening sobre {conformers_xyz.name} (carga={charge}, T={temperature:.2f} K)")
    t0 = time.time()
    try:
        with open(censo_dir / "censo.log", "w") as flog:
            proc = subprocess.run(
                _apply_nohup(cmd), stdout=flog, stderr=subprocess.STDOUT,
                cwd=censo_dir, env=build_omp_env(nprocs_crest),
                timeout=timeout_hours * 3600, start_new_session=USE_NOHUP,
            )
        returncode = proc.returncode
    except subprocess.TimeoutExpired:
        log.error(f"[CENSO]  ✗  TIMEOUT ({timeout_hours} h) — usando ensemble CREST sin filtrar")
        return conformers_xyz
    except FileNotFoundError:
        log.error(f"[CENSO]  ✗  no se encontró el ejecutable ({censo_bin}) — usando ensemble CREST sin filtrar")
        return conformers_xyz
    except Exception as e:
        log.error(f"[CENSO]  ✗  error inesperado: {e} — usando ensemble CREST sin filtrar")
        return conformers_xyz

    elapsed = (time.time() - t0) / 60

    # ── Diagnóstico temprano: parsear censo.log para errores de config ──────
    censo_log_path = censo_dir / "censo.log"
    censo_log_snippet = ""
    if censo_log_path.exists():
        try:
            tail = censo_log_path.read_text(errors="replace").splitlines()
            error_lines = [l for l in tail if any(kw in l for kw in (
                "validation error", "Path for", "is required", "is not set",
                "UserWarning", "No context found",
            ))]
            if error_lines:
                censo_log_snippet = " | ".join(error_lines[-4:])
                log.warning(
                    f"[CENSO]  ✗  CENSO terminó con errores de configuración "
                    f"(returncode={returncode}) — posible causa: ~/.censo2rc incompleto "
                    f"o falta 'tm'/'orca' en [paths]. Detalle: {censo_log_snippet}"
                )
        except Exception:
            pass

    # Reconstrucción preferente desde .out con energías (Gtot)
    post_out_candidates = [
        censo_dir / "1_SCREENING.out", censo_dir / "CENSO" / "1_SCREENING.out",
        censo_dir / "prescreening" / "1_SCREENING.out",
        censo_dir / "0_PRESCREENING.out", censo_dir / "CENSO" / "0_PRESCREENING.out",
        censo_dir / "prescreening" / "0_PRESCREENING.out",
    ]
    for po in post_out_candidates:
        conf_ids = parse_censo_screening_out(po)
        if not conf_ids:
            continue
        geom_sources = [conformers_xyz]
        for xyz_name in ("1_SCREENING.xyz", "0_PRESCREENING.xyz"):
            for cand in (po.parent / xyz_name, po.parent.parent / xyz_name):
                if cand.exists() and cand.stat().st_size > 0:
                    geom_sources.append(cand)
        for geom_src in geom_sources:
            if extract_conformers_by_id(geom_src, conf_ids, output_xyz, log):
                log.info(
                    f"[CENSO]  ✓  pre-screening terminó en {elapsed:.1f} min — "
                    f"{len(conf_ids)} confórmero(s) sobreviven (desde {po.relative_to(censo_dir)})"
                )
                return output_xyz

    # Búsqueda robusta de output (fallback sin energías parseables)
    fixed_candidates = [
        "censo_ensemble_sorted.xyz", "censo_ensemble.xyz",
        "crest_conformers_sorted.xyz", "crest_conformers.xyz.sorted",
    ]
    candidate_paths = [censo_dir / c for c in fixed_candidates]
    candidate_paths += sorted(
        (p for p in censo_dir.rglob("*.xyz") if p != local_input and p.name != output_xyz.name),
        key=lambda p: p.stat().st_mtime if p.exists() else 0, reverse=True,
    )
    for p in candidate_paths:
        if not p.exists() or p.stat().st_size == 0 or p == local_input:
            continue
        if not safe_copy(p, output_xyz):
            log.warning(f"[CENSO] no se pudo normalizar {p.name} — usando ensemble CREST sin filtrar")
            return conformers_xyz
        log.info(f"[CENSO]  ✓  pre-screening terminó en {elapsed:.1f} min ({p.relative_to(censo_dir)})")
        return output_xyz

    if not censo_log_snippet:
        log.warning(
            f"[CENSO]  ✗  no se encontró un ensemble filtrado reconocible "
            f"(returncode={returncode}) en {censo_dir} — revisa {censo_log_path} "
            f"a mano. Usando ensemble CREST sin filtrar."
        )
    return conformers_xyz


def trange_to_temperatures(t_ini: float, t_fin: float, n_points: int) -> list[float]:
    """Genera `n_points` temperaturas equiespaciadas entre t_ini y t_fin
    (sin numpy — linspace manual). n_points==1 -> [t_ini]."""
    if n_points <= 1:
        return [t_ini]
    step = (t_fin - t_ini) / (n_points - 1)
    return [t_ini + i * step for i in range(n_points)]


def run_censo_multitemp(
    conformers_xyz: Path, base_censo_dir: Path, charge: int, unpaired: int,
    temperatures: list[float], nprocs_crest: int, timeout_hours: float,
    log: logging.Logger,
) -> dict[float, Path]:
    """Ported desde pipeline.py::stage_censo_multitemp (~2019)."""
    base_censo_dir.mkdir(parents=True, exist_ok=True)
    results: dict[float, Path] = {}

    log.info(f"[CENSO-MULTITEMP]  ▶  {len(temperatures)} temperaturas: {temperatures[0]:.2f} K … {temperatures[-1]:.2f} K")

    for temp in temperatures:
        subdir_name = f"T{temp:.2f}K"
        log.info(f"[CENSO-MULTITEMP]  ── T = {temp:.2f} K ({subdir_name})")
        result_xyz = run_censo_prescreen(
            conformers_xyz=conformers_xyz, censo_dir=base_censo_dir / subdir_name,
            charge=charge, unpaired=unpaired, temperature=temp,
            nprocs_crest=nprocs_crest, timeout_hours=timeout_hours, log=log,
        )
        results[temp] = result_xyz

    summary_path = base_censo_dir / "multitemp_summary.json"
    summary_data = {
        "temperatures_K": temperatures,
        "results": {f"{temp:.2f}": str(path) for temp, path in results.items()},
        "n_temperatures": len(temperatures),
        "ensemble_input": str(conformers_xyz),
    }
    try:
        with open(summary_path, "w") as f:
            json.dump(summary_data, f, indent=2)
        log.info(f"[CENSO-MULTITEMP]  ✓  resumen guardado en {summary_path.name}")
    except OSError as e:
        log.warning(f"[CENSO-MULTITEMP]  ⚠  no se pudo guardar resumen JSON: {e}")

    n_ok = sum(1 for p in results.values() if p != conformers_xyz)
    log.info(f"[CENSO-MULTITEMP]  ✓  completado — {n_ok}/{len(temperatures)} con ensemble filtrado")
    return results
