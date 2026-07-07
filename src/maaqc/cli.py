"""cli.py — CLI de subcomandos `maaqc` (run / resume / collect / report / doctor).

Cablea, sin reinventar interfaces:

    config/schema.py        -> RecipeConfig, load_recipe()
    stages/base.py           -> Stage, StageResult
    orchestrator/dag.py      -> STAGE_ORDER, resolve_stages_for_entry(),
                                 is_standalone_entry(), validate_stage_ready()
    orchestrator/state.py    -> resolve_resume_index(), needs_rerun(),
                                 record_stage_result()
    analysis/results_collector.py -> collect_results(), print_summary(),
                                 compare_protonation_states(),
                                 print_comparison(), export_results(),
                                 export_comparison() (ya migrado 1:1 del
                                 legacy, ver docstring del módulo)
    analysis/adc2_compare.py -> collect_adc2_results(), print_adc2_summary(),
                                 compare_adc2_tddft(),
                                 print_adc2_tddft_comparison(),
                                 export_adc2_comparison() (ya migrado 1:1
                                 del legacy; results_collector.py ya tenía
                                 export_adc2_mol() migrado, pero cmd_collect
                                 no lo invocaba porque no había de dónde
                                 sacar df_adc2 — este módulo lo resuelve)
    analysis/report.py       -> build_report_data(), write_report()
                                 (interfaz NUEVA, no existía en el legacy;
                                 definida acá porque `report` no tenía
                                 contrato previo)

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
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Optional

from maaqc.analysis import adc2_compare as adc2
from maaqc.analysis import report as report_mod
from maaqc.analysis import results_collector as rc
from maaqc.config.schema import RecipeConfig, load_recipe
from maaqc.orchestrator import dag, state
from maaqc.stages.base import Stage
from maaqc.system.platform_detect import resolve_bin

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


def _resolve_states(run_root: Path, states_arg: str, purpose: str) -> Optional[list[str]]:
    """Parsea `--states a,b,c`, avisa (sin abortar) de los que no tienen
    subdir bajo `run_root`, y devuelve la lista final. `None` si no queda
    ningún estado usable (el caller debe abortar el comando en ese caso).

    Compartido por `collect` y `report`: ambos reciben el mismo
    `run_root` (dir que contiene un subdir por estado, ej. runs/A,
    runs/B) y el mismo criterio de qué estados existen de verdad.
    """
    states = [s.strip() for s in states_arg.split(",") if s.strip()]
    if not states:
        print(f"[maaqc] '--states' vacío, no hay nada que {purpose}.", file=sys.stderr)
        return None

    missing = [s for s in states if not (run_root / s).exists()]
    if missing:
        print(
            f"[maaqc] Aviso: no se encontró subdir para los estados {missing} "
            f"bajo {run_root}; se omiten.",
            file=sys.stderr,
        )
        states = [s for s in states if s not in missing]

    if not states:
        print(f"[maaqc] Ningún estado disponible para {purpose}.", file=sys.stderr)
        return None
    return states


def _resolve_reference(ref_arg: Optional[str], states: list[str]) -> Optional[str]:
    """None si no hay --reference o si no pertenece a --states (con
    aviso); si no, el nombre de estado a usar como referencia."""
    if ref_arg is None:
        return None
    if ref_arg not in states:
        print(
            f"[maaqc] Aviso: --reference '{ref_arg}' no está en --states "
            f"{states}; se continúa sin referencia.",
            file=sys.stderr,
        )
        return None
    return ref_arg


def cmd_collect(args: argparse.Namespace) -> int:
    run_root = Path(args.workdir)
    if not run_root.exists():
        print(f"[maaqc] '--workdir {run_root}' no existe.", file=sys.stderr)
        return 1

    states = _resolve_states(run_root, args.states, "consolidar")
    if states is None:
        return 1
    ref_mol = _resolve_reference(args.reference, states)

    # Recolección + impresión por estado individual (mol_dir = run_root/
    # <estado>, mismo patrón que Stage.outputs(workdir)). Un estado con
    # error no aborta a los demás; se reporta y se sigue.
    any_error = False
    for state_name in states:
        mol_dir = run_root / state_name
        try:
            df_conf, df_tddft = rc.collect_results(mol_dir)
        except FileNotFoundError as exc:
            print(f"[maaqc] {state_name}: {exc}", file=sys.stderr)
            any_error = True
            continue
        rc.print_summary(df_conf, df_tddft, mol_name=state_name)

        # ADC2 (opcional, por estado): arma df_adc2 vía
        # adc2_compare.collect_adc2_results() — vacío si el estado todavía
        # no corrió 07_adc2/, lo cual es normal (ADC2 es la etapa final).
        # Se imprime igual que print_summary, sin depender de --export,
        # para que `maaqc collect` sin --export ya muestre ADC2 si existe
        # (mismo criterio que confórmeros/TDDFT arriba). export_adc2_mol()
        # ya estaba migrado en results_collector.py pero nada lo invocaba
        # porque no había de dónde sacar df_adc2 hasta ahora.
        df_adc2 = adc2.collect_adc2_results(mol_dir, mol_name=state_name)
        if not df_adc2.empty:
            adc2.print_adc2_summary(df_adc2, mol_name=state_name, df_conf=df_conf)

        if args.export:
            rc.export_results(mol_dir, df_conf, df_tddft, mol_name=state_name, fmt=args.export)
            if not df_adc2.empty:
                rc.export_adc2_mol(mol_dir, df_adc2, mol_name=state_name, fmt=args.export)

    # Comparación entre estados. compare_protonation_states redescubre
    # TODOS los subdirs con 04_sp/ bajo run_root; se filtra después a los
    # --states pedidos para no comparar estados que el usuario no pidió
    # aunque existan en disco (ej. una corrida vieja abandonada al lado).
    df_cmp, _dfs_conf, _dfs_tddft = rc.compare_protonation_states(run_root, ref_mol=ref_mol)
    if not df_cmp.empty:
        df_cmp = df_cmp[df_cmp["molecule"].isin(states)].reset_index(drop=True)
    rc.print_comparison(df_cmp, ref_mol=ref_mol)

    # Comparación ADC2 vs TDDFT (el benchmark central del paper). A
    # diferencia de compare_protonation_states, compare_adc2_tddft() ya
    # acepta `mols` explícito, así que se le pasan directo los --states
    # resueltos en vez de filtrar después.
    df_adc2_cmp = adc2.compare_adc2_tddft(run_root, mols=states)
    if not df_adc2_cmp.empty:
        adc2.print_adc2_tddft_comparison(df_adc2_cmp)

    if args.export:
        rc.export_comparison(run_root, df_cmp, fmt=args.export)
        if not df_adc2_cmp.empty:
            adc2.export_adc2_comparison(run_root, df_adc2_cmp, fmt=args.export)

    return 1 if any_error else 0


def cmd_report(args: argparse.Namespace) -> int:
    run_root = Path(args.workdir)
    if not run_root.exists():
        print(f"[maaqc] '--workdir {run_root}' no existe.", file=sys.stderr)
        return 1

    states = _resolve_states(run_root, args.states, "reportar")
    if states is None:
        return 1
    ref_mol = _resolve_reference(args.reference, states)

    # Interfaz de `report` (nueva, no migrada del legacy: pipeline.py no
    # tenía generación de reportes, solo print_summary/print_comparison a
    # stdout). Definida en analysis/report.py, reutilizando los DataFrames
    # de results_collector.py sin reimplementar collect/compare.
    data = report_mod.build_report_data(run_root, states, ref_mol=ref_mol)
    try:
        out_path = report_mod.write_report(data, run_root, args.format)
    except ValueError as exc:
        print(f"[maaqc] {exc}", file=sys.stderr)
        return 1

    print(f"[maaqc] Reporte generado: {out_path}")
    return 0


DOCTOR_BINARIES: dict[str, list[str]] = {
    "orca": ["orca", "--version"],
    "crest": ["crest", "--version"],
    "censo": ["censo", "--version"],
    "xtb": ["xtb", "--version"],
}


# Instrucciones de instalación por binario faltante, mostradas por
# `maaqc doctor` cuando `_check_binary` no lo encuentra. crest/censo/xtb
# tienen paquete conda-forge (mismo comando en ThinkPad Ubuntu y en el
# MacBook Air institucional, porque `conda`/`mamba` no dependen del SO);
# ORCA es la excepción real — requiere licencia académica y se
# descarga a mano, no hay gestor de paquetes que lo resuelva. Por eso
# para ORCA el hint no es un comando, sino los pasos + la env var que
# `resolve_bin` ya sabe leer (MAAQC_ORCA_BIN, ver platform_detect.py).
def _install_hint(name: str) -> str:
    is_macos = platform.system() == "Darwin"

    if name == "orca":
        return (
            "ORCA no tiene gestor de paquetes (requiere licencia académica):\n"
            "      1. Registrate y descargá el tarball en https://orcaforum.kofo.mpg.de/\n"
            "         (sección 'ORCA downloads', elegí la versión para "
            + ("macOS" if is_macos else "Linux")
            + ").\n"
            "      2. Descomprimilo donde te resulte cómodo (no requiere instalación).\n"
            "      3. Si el ejecutable `orca` no queda en tu PATH, fijalo explícito:\n"
            "         export MAAQC_ORCA_BIN=/ruta/a/orca_6_1_1/orca\n"
            "         (resolve_bin() respeta esta variable antes que cualquier búsqueda)."
        )

    if name in ("crest", "censo", "xtb"):
        # Paquete conda-forge disponible para los tres — mismo comando
        # en ambas máquinas. Se sugiere el env activo si hay uno, o
        # 'sciencebitch' (env de referencia de kmilo) como default.
        return f"conda install -c conda-forge {name}   # (o: mamba install -c conda-forge {name})"

    return "sin instrucciones de instalación registradas para este binario."


def _check_binary(name: str, version_cmd: list[str]) -> tuple[bool, str]:
    # resolve_bin busca, en orden: override MAAQC_<NAME>_BIN -> PATH
    # actual -> TODOS los envs de conda (no uno hardcodeado) -> rutas
    # típicas de instalación manual (~, /opt, /usr/local, /Applications).
    # Esto reemplaza el shutil.which() sobre el PATH crudo, que fallaba
    # para crest (vive en el env base de conda, no en el env `maaqc`
    # activo) y para orca (instalación manual fuera del PATH).
    resolved = resolve_bin(version_cmd[0])
    if not (resolved.is_absolute() and resolved.is_file()):
        return False, "no encontrado en PATH, en ningún env de conda, ni en rutas típicas"
    version_cmd = [str(resolved)] + version_cmd[1:]
    try:
        proc = subprocess.run(
            version_cmd, capture_output=True, text=True, timeout=10,
        )
        first_line = (proc.stdout or proc.stderr).strip().splitlines()
        version = first_line[0] if first_line else "(sin salida de versión)"
    except (OSError, subprocess.TimeoutExpired) as exc:
        return True, f"encontrado en {resolved}, pero no se pudo consultar versión ({exc})"
    return True, f"{resolved} — {version}"


def cmd_doctor(args: argparse.Namespace) -> int:
    print("[maaqc] Chequeando binarios...")
    all_ok = True
    for name, version_cmd in DOCTOR_BINARIES.items():
        ok, detail = _check_binary(name, version_cmd)
        status = "OK" if ok else "FALTA"
        print(f"  [{status}] {name}: {detail}")
        if not ok:
            print(f"      -> instalar: {_install_hint(name)}")
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
    p_collect.add_argument("--workdir", required=True, help="Dir que contiene un subdir por estado (ej. runs/A, runs/B)")
    p_collect.add_argument("--states", default="A,B,C,D,R", help="Estados a consolidar, separados por coma")
    p_collect.add_argument(
        "--reference", default=None,
        help="Estado a usar como referencia para ΔE (ej. R). Debe estar en --states",
    )
    p_collect.add_argument(
        "--export", choices=["csv", "xlsx", "both"], default=None,
        help="Además de imprimir, exportar resultados por estado y la comparación (export_results/export_comparison)",
    )
    p_collect.set_defaults(func=cmd_collect)

    p_report = subparsers.add_parser("report", help="Genera un reporte de la corrida")
    p_report.add_argument("--workdir", required=True, help="Dir que contiene un subdir por estado (ej. runs/A, runs/B)")
    p_report.add_argument("--states", default="A,B,C,D,R", help="Estados a incluir en el reporte, separados por coma")
    p_report.add_argument(
        "--reference", default=None,
        help="Estado a usar como referencia para ΔE (ej. R). Debe estar en --states",
    )
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
