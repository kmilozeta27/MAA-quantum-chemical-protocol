"""Adc2Stage — ADC2 sobre el confórmero de menor energía SP. OPCIONAL.

Migrado desde pipeline.py:
  - stage_adc2()          líneas 2972-3135
  - extract_adc2_lambda() líneas 929-998
  - load_exp_lambda()     líneas 270-351

No corre por defecto: requiere config.adc2.mandatory = True (equivalente
a --adc2 explícito en el legacy) o config.entry == "sp-only" (equivalente
a --start adc2 standalone, hoy manejado por
orchestrator/dag.py::STANDALONE_ENTRIES).

`requires = ["tddft"]` se mantiene ESTÁTICO — ver dag.py para cómo se
resuelve el caso --entry sp-only sin romper el contrato de Stage. Este
archivo no necesita saber nada sobre --entry: is_complete()/inputs()
siempre validan lo mismo (04_sp/ del confórmero dominante); es dag.py
quien decide si además exige que 'tddft' esté completo.
"""

from __future__ import annotations

import datetime
import json
import logging
import re
import shutil
from pathlib import Path
from typing import Optional

from maaqc.config.schema import RecipeConfig
from maaqc.stages.base import Stage, StageResult

# TODO(migración): estas utilidades ya existen en pipeline.py y presumiblemente
# ya se migraron junto con opt.py/sp.py/tddft.py (que según el chat ya están
# migrados). Ajustar el import real cuando se confirme dónde quedaron
# (candidato: maaqc.io.orca_utils). Se listan explícitas acá para dejar clara
# la dependencia y no duplicar lógica de parsing ORCA.
from maaqc.io.orca_utils import (  # type: ignore[import-not-found]
    check_success,
    extract_sp_energy,
    get_template,
    resolve_opt_xyz,
    run_orca,
    safe_copy,
)

log = logging.getLogger(__name__)

HARTREE_TO_KCAL = 627.5094740631
MULT = 1  # multiplicidad siempre 1 (idéntico a pipeline.py:173)

# ── λmax experimental de fallback (pipeline.py:252-265) ──────────────────
# Orden de búsqueda idéntico al legacy: 1) config.adc2.exp_lambda_nm (nuevo,
# reemplaza molecules.yml) 2) este diccionario 3) None.
EXP_LAMBDA_FALLBACK: dict[str, float] = {
    "A": 335.0, "B": 335.0, "C": 335.0, "D": 335.0, "R": 335.0,
    "SA": 334.0, "SB": 334.0, "SC": 334.0, "SD": 334.0, "SR": 334.0,
}


def load_exp_lambda(mol_name: str, config: RecipeConfig) -> Optional[float]:
    """λmax experimental (nm) para mol_name.

    Orden de búsqueda (adaptado de pipeline.py:load_exp_lambda, que leía
    molecules.yml en texto plano — reemplazado acá por el YAML de recipe,
    que ya pasa por schema/pydantic):
      1. config.adc2.exp_lambda_nm   (declarado explícito en el recipe)
      2. EXP_LAMBDA_FALLBACK         (compatibilidad con runs legacy)
      3. None (stage_adc2 omite Δλ pero no falla)
    """
    exp_from_config = getattr(getattr(config, "adc2", None), "exp_lambda_nm", None)
    if exp_from_config is not None:
        return float(exp_from_config)

    if mol_name in EXP_LAMBDA_FALLBACK:
        log.debug(
            "[EXP_LAMBDA] '%s' no está en config.adc2.exp_lambda_nm — "
            "usando fallback (%.1f nm)", mol_name, EXP_LAMBDA_FALLBACK[mol_name],
        )
        return EXP_LAMBDA_FALLBACK[mol_name]

    log.warning(
        "[EXP_LAMBDA] λmax experimental no disponible para '%s'. "
        "Añádelo en config.adc2.exp_lambda_nm para obtener Δλ.", mol_name,
    )
    return None


def extract_adc2_lambda(out_file: Path) -> Optional[tuple[float, float, float]]:
    """Extrae λmax (nm), E (eV) y fosc de un output ADC(2)/ORCA (nroots 1,
    S0→S1). Migrado sin cambios de lógica desde pipeline.py:929-998.

    Busca la tabla "SPECTRUM FOR RIGHT TRANSITION MOMENTS" /
    "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" y toma la
    primera fila de datos (válido para nroots=1).
    """
    if not out_file.exists():
        return None

    _data_re = re.compile(
        r"^\s*\S+\s+->\s+\S+"
        r"\s+([\d.]+)"
        r"\s+[\d.]+\s+([\d.]+)"
        r"\s+([\d.]+)"
    )
    _right_header = "SPECTRUM FOR RIGHT TRANSITION MOMENTS"
    _abs_header = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"

    try:
        text = out_file.read_text(errors="replace")
    except OSError:
        return None

    right_pos = text.find(_right_header)
    if right_pos == -1:
        return None

    abs_pos = text.find(_abs_header, right_pos)
    if abs_pos == -1:
        return None

    segment = text[abs_pos:]
    for line in segment.splitlines():
        m = _data_re.match(line)
        if m:
            try:
                energy_ev = float(m.group(1))
                lambda_nm = float(m.group(2))
                fosc = float(m.group(3))
                return (lambda_nm, energy_ev, fosc)
            except ValueError:
                return None

    return None


def _select_dominant_conformer(sp_root: Path) -> Optional[tuple[float, Path]]:
    """Confórmero dominante: menor energía SP entre los 04_sp/conf_* con
    sp.out exitoso. Factoriza el bloque de selección de
    pipeline.py:stage_adc2 (líneas 3015-3037) para reusarlo en
    is_complete()/inputs()/run() sin duplicar lógica.
    """
    energies: list[tuple[float, Path]] = []
    for conf_dir in sorted(sp_root.glob("conf_*")):
        out_file = conf_dir / "sp.out"
        if not check_success(out_file):
            continue
        e = extract_sp_energy(out_file)
        if e is not None:
            energies.append((e, conf_dir))
    if not energies:
        return None
    energies.sort(key=lambda x: x[0])
    return energies[0]


class Adc2Stage(Stage):
    name = "adc2"
    requires = ["tddft"]  # resuelto por dag.py según --entry, no acá

    def is_complete(self, workdir: Path) -> bool:
        sp_root = workdir / "04_sp"
        if not sp_root.exists():
            return False
        dominant = _select_dominant_conformer(sp_root)
        if dominant is None:
            return False
        _, sp_dir = dominant
        adc2_out = workdir / "07_adc2" / sp_dir.name / "adc2.out"
        return check_success(adc2_out)

    def inputs(self, workdir: Path) -> list[Path]:
        # El input real es "04_sp/ con al menos un confórmero exitoso", pero
        # el contrato pide list[Path]: se valida la existencia del directorio
        # raíz acá; la selección del confórmero dominante (y su geometría)
        # se resuelve dentro de run()/is_complete(), donde si falta se
        # reporta un error específico en vez de un Path genérico.
        return [workdir / "04_sp"]

    def outputs(self, workdir: Path) -> list[Path]:
        return [workdir / "07_adc2" / "adc2_results.json"]

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        """Migrado desde pipeline.py:stage_adc2 (líneas 2972-3135), adaptado
        a la firma Stage.run(config, workdir) — mol_dir/charge/mol_name/log
        del legacy vienen ahora de config y workdir.
        """
        mol_name = getattr(config, "molecule_name", None) or workdir.name
        charge = config.charge
        nprocs_adc2 = getattr(config.adc2, "nprocs", 4)
        maxcore_adc2 = getattr(config.adc2, "maxcore_mb", 5500)

        sp_root = workdir / "04_sp"
        if not sp_root.exists():
            return StageResult(ok=False, message="No existe 04_sp/ — corre primero SP.")

        dominant = _select_dominant_conformer(sp_root)
        if dominant is None:
            return StageResult(
                ok=False,
                message="Ningún confórmero SP exitoso con energía válida.",
            )
        e_min, sp_dir = dominant
        conf_id = sp_dir.name
        log.info(
            "[ADC2] confórmero dominante: %s (E_SP = %.10f Ha)", conf_id, e_min,
        )

        geom_src = resolve_opt_xyz(sp_dir)
        if geom_src is None:
            return StageResult(
                ok=False,
                message=(
                    f"No se encontró opt_opt.xyz ni opt.xyz en {sp_dir} — "
                    f"el confórmero dominante ({conf_id}) tiene sp.out pero "
                    f"le falta la geometría en 04_sp/."
                ),
            )

        adc2_root = workdir / "07_adc2"
        adc2_dir = adc2_root / conf_id
        adc2_dir.mkdir(parents=True, exist_ok=True)

        out_file = adc2_dir / "adc2.out"
        inp_file = adc2_dir / "adc2.inp"

        if not check_success(out_file):
            if not safe_copy(geom_src, adc2_dir / "adc2_geom.xyz"):
                return StageResult(
                    ok=False,
                    message=f"{conf_id}: error copiando geometría.",
                )
            template = get_template("adc2")
            inp_text = template.format(
                CHARGE=charge, MULT=MULT,
                NPROCS_ADC2=nprocs_adc2, MAXCORE_ADC2=maxcore_adc2,
            )
            inp_file.write_text(inp_text)

        ok = run_orca(inp_file, out_file, log)
        if not ok:
            return StageResult(ok=False, message=f"{conf_id}: ADC2 falló.")

        log.info("[ADC2] %s completado ✓ → %s", conf_id, adc2_dir)

        result = extract_adc2_lambda(out_file)
        record: dict = {}
        if result is not None:
            lambda_nm, energy_ev, fosc = result
            lambda_exp = load_exp_lambda(mol_name, config)
            delta_lambda = (lambda_nm - lambda_exp) if lambda_exp is not None else None

            record = {
                "mol_name": mol_name,
                "conf_id": conf_id,
                "lambda_nm": round(lambda_nm, 2),
                "energy_eV": round(energy_ev, 6),
                "fosc": round(fosc, 8),
                "lambda_exp_nm": lambda_exp,
                "delta_lambda_nm": round(delta_lambda, 2) if delta_lambda is not None else None,
                "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
            }
            results_json = adc2_root / "adc2_results.json"
            try:
                results_json.write_text(
                    json.dumps(record, indent=2, ensure_ascii=False) + "\n"
                )
                log.info("[ADC2] resultados guardados → %s", results_json)
            except OSError as e:
                log.warning("[ADC2] no se pudo escribir %s: %s", results_json, e)
        else:
            log.warning(
                "[ADC2] %s: cálculo OK pero no se pudo extraer λmax "
                "del espectro RIGHT — revisar %s", conf_id, out_file.name,
            )

        return StageResult(ok=True, message=f"ADC2 completado en {adc2_dir}", data=record)
