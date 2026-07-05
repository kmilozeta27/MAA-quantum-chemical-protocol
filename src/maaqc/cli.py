"""cli.py — CLI de subcomandos `maaqc` (run / resume / collect / report / doctor).

Cablea, sin reinventar interfaces:

    config/schema.py        -> RecipeConfig, load_recipe()
    stages/base.py           -> Stage, StageResult
    orchestrator/dag.py      -> STAGE_ORDER, resolve_stages_for_entry(),
                                 is_standalone_entry(), validate_stage_ready()
    orchestrator/state.py    -> resolve_resume_index(), needs_rerun(),
                                 record_stage_result()

Reemplaza los flags sueltos de pipeline.py legacy (--start adc2,
--skipnumfreq, --expand-crest) por subcomandos + --entry + recipe YAML,
según el contrato acordado en la skill maa-repackaging.

Registro de stages (`STAGE_CLASS_NAMES`): la migración de pipeline.py es
incremental, stage por stage (ver skill, sección 7 — "orden de
migración"). Este módulo NO asume que los 7 stages ya existen: importa
cada uno de forma perezosa y, si `maaqc.stages.<name>` todavía no fue
migrado, lo reporta como "pendiente" en vez de reventar con un
ImportError críptico al armar el parser.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

from maaqc.config.schema import RecipeConfig, load_recipe
from maaqc.orchestrator import dag, state
from maaqc.stages.base import Stage

# Nombre de clase esperado por módulo de stage (confirmado por dag.py,
# que ya nombra "SpStage" y "Adc2Stage" en su propio docstring).
STAGE_CLASS_NAMES: dict[str, str] = {
    "crest": "CrestStage",
    "censo": "CensoStage",
    "opt": "OptStage",
    "numfreq": "NumFreqStage",
    "sp": "SpStage",
    "tddft": "TddftStage",
    "adc2": "Adc2Stage",
}


# ─────────────────────────────────────────────
#  REGISTRO DE STAGES (import perezoso, tolera migración parcial)
# ─────────────────────────────────────────────

def _stage_has_stub_methods(stage: Stage, probe_dir: Path) -> Optional[str]:
    """None si `stage` tiene lógica real; si no, el motivo por el que se
    lo considera un stub sin migrar todavía.

    Detecta el patrón usado en TODA la migración (crest.py, censo.py,
    numfreq.py, tddft.py): clases que YA implementan los 4 métodos
    abstractos de `Stage` (por eso Python las deja instanciar sin
    quejarse — ABC solo exige que el método EXISTA, no que tenga
    cuerpo real) pero cuyo cuerpo es únicamente `raise
    NotImplementedError(...)`. Antes de este chequeo, `maaqc doctor`
    reportaba estos 4 stages como "[OK] migrado" y `maaqc run` recién
    reventaba a mitad de la cadena, en el primer `is_complete()` real.

    Se prueba `is_complete()` con un `probe_dir` vacío (no `workdir` del
    run real) porque es el único método del contrato que no debería
    tener efectos secundarios para un directorio inexistente en ninguno
    de los stages ya migrados (opt/sp/adc2 simplemente devuelven False).
    Solo se captura `NotImplementedError` — cualquier otra excepción es
    un bug real y debe propagarse, no esconderse acá.
    """
    try:
        stage.is_complete(probe_dir)
    except NotImplementedError as exc:
        detail = str(exc) or "is_complete() sin implementar"
        return f"stub sin migrar ({detail})"
    return None


def _load_stage_registry() -> tuple[dict[str, Stage], dict[str, str]]:
    """Instancia los stages ya migrados a `maaqc.stages.<name>`.

    Devuelve (registro, pendientes):
        registro:   {nombre: instancia de Stage} para los ya migrados
                     CON lógica real (ver `_stage_has_stub_methods`).
        pendientes: {nombre: motivo} para los que faltan — import
                     falló, la clase no existe todavía, o existe pero
                     sigue siendo un stub (`raise NotImplementedError`).

    No asume que los 7 stages existen — eso rompería `maaqc doctor` y
    el `--help` mismo mientras la migración de pipeline.py está en
    curso (ver skill, orden de migración sugerido: opt/sp primero).
    """
    registry: dict[str, Stage] = {}
    pending: dict[str, str] = {}
    probe_dir = Path("/nonexistent-maaqc-stage-probe")
    for name, class_name in STAGE_CLASS_NAMES.items():
        try:
            module = __import__(f"maaqc.stages.{name}", fromlist=[class_name])
        except ImportError as exc:
            pending[name] = f"módulo maaqc.stages.{name} no encontrado ({exc})"
            continue
        stage_cls = getattr(module, class_name, None)
        if stage_cls is None:
            pending[name] = f"{class_name} no existe en maaqc.stages.{name}"
            continue
        instance = stage_cls()
        stub_reason = _stage_has_stub_methods(instance, probe_dir)
        if stub_reason is not None:
            pending[name] = stub_reason
            continue
        registry[name] = instance
    return registry, pending


def _require_stages(
    stage_names: list[str], registry: dict[str, Stage], pending: dict[str, str],
) -> Optional[str]:
    """None si todos los `stage_names` están migrados; si no, un mensaje
    de error listando qué falta (para abortar el comando con claridad)."""
    missing = [n for n in stage_names if n not in registry]
    if not missing:
        return None
    detalles = "; ".join(f"{n}: {pending.get(n, 'no migrado')}" for n in missing)
    return f"Stages sin migrar todavía, no se puede continuar: {detalles}"


# ─────────────────────────────────────────────
#  CONFIG: carga + override de --entry por CLI
# ─────────────────────────────────────────────

def _load_config(config_path: Path, entry_override: Optional[str]) -> RecipeConfig:
    config = load_recipe(config_path)
    if entry_override is not None and entry_override != config.entry:
        config = config.model_copy(update={"entry": entry_override})
    return config


# ─────────────────────────────────────────────
#  EJECUCIÓN DE LA CADENA (compartida por run/resume)
# ─────────────────────────────────────────────

def _execute_chain(
    workdir: Path,
    config: RecipeConfig,
    start_index: int,
    registry: dict[str, Stage],
    pending: dict[str, str],
) -> int:
    stage_names = dag.resolve_stages_for_entry(config.entry)

    error = _require_stages(stage_names, registry, pending)
    if error:
        print(f"[maaqc] {error}", file=sys.stderr)
        return 1

    if start_index >= len(stage_names):
        print("[maaqc] Nada que correr: toda la cadena ya está al día.")
        return 0

    for name in stage_names[start_index:]:
        stage = registry[name]

        ready, reason = dag.validate_stage_ready(stage, workdir, config.entry, registry)
        if not ready:
            print(f"[maaqc] Abortando en '{name}': {reason}", file=sys.stderr)
            return 1

        print(f"[maaqc] Corriendo stage '{name}'...")
        result = stage.run(config, workdir)
        if not result.ok:
            print(f"[maaqc] Stage '{name}' falló: {result.message}", file=sys.stderr)
            return 1

        state.record_stage_result(
            workdir, name, config, stage.outputs(workdir), message=result.message,
        )
        print(f"[maaqc] Stage '{name}' completo. {result.message}".rstrip())

    print("[maaqc] Cadena completa.")
    return 0


# ─────────────────────────────────────────────
#  SUBCOMANDOS
# ─────────────────────────────────────────────

def cmd_run(args: argparse.Namespace) -> int:
    workdir = Path(args.workdir) if args.workdir else Path(args.input).resolve().parent
    config = _load_config(Path(args.config), args.entry)
    if args.reference is not None:
        config = config.model_copy(update={"reference": args.reference})

    registry, pending = _load_stage_registry()
    stage_names = dag.resolve_stages_for_entry(config.entry)

    error = _require_stages(stage_names, registry, pending)
    if error:
        print(f"[maaqc] {error}", file=sys.stderr)
        return 1

    start_index = state.resolve_resume_index(workdir, stage_names, registry, config)
    return _execute_chain(workdir, config, start_index, registry, pending)


def cmd_resume(args: argparse.Namespace) -> int:
    workdir = Path(args.workdir)
    config = _load_config(Path(args.config), args.entry)

    registry, pending = _load_stage_registry()
    stage_names = dag.resolve_stages_for_entry(config.entry)

    if args.from_stage is not None:
        if args.from_stage not in stage_names:
            print(
                f"[maaqc] '--from {args.from_stage}' no pertenece a la cadena "
                f"resuelta para entry='{config.entry}': {stage_names}",
                file=sys.stderr,
            )
            return 1
        start_index = stage_names.index(args.from_stage)
    else:
        error = _require_stages(stage_names, registry, pending)
        if error:
            print(f"[maaqc] {error}", file=sys.stderr)
            return 1
        start_index = state.resolve_resume_index(workdir, stage_names, registry, config)

    return _execute_chain(workdir, config, start_index, registry, pending)


def cmd_collect(args: argparse.Namespace) -> int:
    # TODO: cablear analysis/results_collector.py una vez migrado desde
    # results_collector.py legacy (mismo patrón que load_recipe: no
    # inventar la interfaz acá, esperar el módulo real). Estados
    # esperados vía --states A,B,C,D,R (ver protonation_cycle.py).
    print(
        "[maaqc] 'collect' pendiente: falta migrar "
        "analysis/results_collector.py.",
        file=sys.stderr,
    )
    return 1


def cmd_report(args: argparse.Namespace) -> int:
    # TODO: mismo caso que collect — pendiente de un módulo de reporte
    # real (analysis/ o un report.py nuevo) antes de fijar la interfaz.
    print(
        "[maaqc] 'report' pendiente: falta definir el módulo de reporte.",
        file=sys.stderr,
    )
    return 1


DOCTOR_BINARIES: dict[str, list[str]] = {
    "orca": ["orca", "--version"],
    "crest": ["crest", "--version"],
    "xtb": ["xtb", "--version"],
}


def _check_binary(name: str, version_cmd: list[str]) -> tuple[bool, str]:
    path = shutil.which(version_cmd[0])
    if path is None:
        return False, "no encontrado en PATH"
    try:
        proc = subprocess.run(
            version_cmd, capture_output=True, text=True, timeout=10,
        )
        first_line = (proc.stdout or proc.stderr).strip().splitlines()
        version = first_line[0] if first_line else "(sin salida de versión)"
    except (OSError, subprocess.TimeoutExpired) as exc:
        return True, f"encontrado en {path}, pero no se pudo consultar versión ({exc})"
    return True, f"{path} — {version}"


def cmd_doctor(args: argparse.Namespace) -> int:
    print("[maaqc] Chequeando binarios...")
    all_ok = True
    for name, version_cmd in DOCTOR_BINARIES.items():
        ok, detail = _check_binary(name, version_cmd)
        status = "OK" if ok else "FALTA"
        print(f"  [{status}] {name}: {detail}")
        all_ok = all_ok and ok

    print("[maaqc] Chequeando stages migrados...")
    registry, pending = _load_stage_registry()
    for name in dag.STAGE_ORDER:
        if name in registry:
            print(f"  [OK] stage '{name}' migrado")
        else:
            print(f"  [PENDIENTE] stage '{name}': {pending[name]}")

    return 0 if all_ok else 1


# ─────────────────────────────────────────────
#  PARSER
# ─────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="maaqc",
        description="MAA Quantum Characterization — CREST/CENSO/OPT/NumFreq/SP/TDDFT/ADC2",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    entry_choices = ["smiles", "xyz", "crest-ensemble", "sp-only", "calibrated"]

    p_run = subparsers.add_parser("run", help="Corre (o retoma) la cadena completa")
    p_run.add_argument("input", help="Molécula de entrada (SMILES, .xyz, etc.)")
    p_run.add_argument("--config", required=True, help="Recipe YAML (ej. publication.yml)")
    p_run.add_argument("--workdir", default=None, help="Directorio de trabajo (default: junto al input)")
    p_run.add_argument("--entry", choices=entry_choices, default=None, help="Override de config.entry")
    p_run.add_argument("--reference", default=None, help="Para --entry calibrated (ej. shinorine)")
    p_run.set_defaults(func=cmd_run)

    p_resume = subparsers.add_parser("resume", help="Retoma una corrida existente en --workdir")
    p_resume.add_argument("--workdir", required=True)
    p_resume.add_argument("--config", required=True)
    p_resume.add_argument("--entry", choices=entry_choices, default=None)
    p_resume.add_argument(
        "--from", dest="from_stage", default=None, choices=list(dag.STAGE_ORDER),
        help="Forzar retoma desde este stage (ignora resume automático vía state.json)",
    )
    p_resume.set_defaults(func=cmd_resume)

    p_collect = subparsers.add_parser("collect", help="Consolida resultados entre estados de protonación")
    p_collect.add_argument("--workdir", required=True)
    p_collect.add_argument("--states", default="A,B,C,D,R", help="Estados a consolidar, separados por coma")
    p_collect.add_argument("--output", default=None)
    p_collect.set_defaults(func=cmd_collect)

    p_report = subparsers.add_parser("report", help="Genera un reporte de la corrida")
    p_report.add_argument("--workdir", required=True)
    p_report.add_argument("--format", choices=["html", "md", "json"], default="html")
    p_report.set_defaults(func=cmd_report)

    p_doctor = subparsers.add_parser("doctor", help="Chequea binarios (ORCA/CREST/xtb) y stages migrados")
    p_doctor.set_defaults(func=cmd_doctor)

    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
