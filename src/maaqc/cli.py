"""CLI de maaqc — subcomandos estilo git/CENSO.

Reemplaza los flags sueltos del pipeline legacy:

    pipeline.py --mol A.xyz --charge -1                  -> maaqc run
    pipeline.py --mol A.xyz --charge -1 --start opt       -> maaqc resume --from opt
    pipeline.py --mol A.xyz --charge -1 --start sp        -> maaqc resume --from sp
    pipeline.py --mol A.xyz --charge -1 --start adc2      -> maaqc run --entry sp-only
    pipeline.py --status                                  -> maaqc doctor / maaqc resume --list
    results_collector.py (su propio main())               -> maaqc collect / maaqc report

TODO (en orden, ver docs/troubleshooting.md "Orden de migración"):
    1. `doctor`   — trivial, migra check_required_bins()/_resolve_bin()
                    de pipeline.py:93-158 casi sin cambios.
    2. `run`      — depende de que orchestrator/dag.py y al menos un
                    Stage real (probablemente StageOpt) estén migrados.
    3. `resume`   — depende de orchestrator/state.py (generaliza
                    recover_prescreening.py + _censo_prescreen_recover
                    de pipeline.py:1342-1437).
    4. `collect` / `report` — envuelven analysis/results_collector.py
                    una vez migrado desde el results_collector.py legacy.

No se implementa lógica real todavía: cada subcomando levanta
NotImplementedError con un mensaje que indica qué migrar antes de
poder usarlo, para que `maaqc --help` ya sea navegable durante la
migración incremental.
"""

from __future__ import annotations

import argparse
import sys


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="maaqc",
        description="MAA Quantum Characterization pipeline (CREST -> CENSO -> OPT -> NumFreq -> SP -> TDDFT -> ADC2)",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # ── maaqc run ────────────────────────────────────────────────────
    p_run = subparsers.add_parser("run", help="Corre el pipeline desde cero (o desde un --entry específico)")
    p_run.add_argument("molecule", type=str, help="Archivo de entrada (.xyz, .smi, según --entry)")
    p_run.add_argument("--config", type=str, required=True, metavar="RECIPE.yml",
                        help="Recipe YAML (ver src/maaqc/config/presets/)")
    p_run.add_argument(
        "--entry",
        choices=["smiles", "xyz", "crest-ensemble", "sp-only", "calibrated"],
        default="xyz",
        help=(
            "Punto de entrada al pipeline. 'crest-ensemble' reemplaza "
            "--expand-crest; 'sp-only' reemplaza --start adc2; "
            "'calibrated' usa shifts empíricos (--reference)."
        ),
    )
    p_run.add_argument("--reference", type=str, default=None,
                        help="Molécula de referencia para --entry calibrated (ej. shinorine)")
    p_run.set_defaults(func=_cmd_run)

    # ── maaqc resume ─────────────────────────────────────────────────
    p_resume = subparsers.add_parser("resume", help="Retoma un run existente desde el último stage completo (o uno explícito)")
    p_resume.add_argument("--from", dest="from_stage", type=str, default=None,
                           choices=["crest", "censo", "opt", "numfreq", "sp", "tddft", "adc2"],
                           help="Forzar reanudación desde este stage (default: autodetectar vía orchestrator/state.py)")
    p_resume.add_argument("--list", action="store_true",
                           help="Solo lista el estado de todos los runs (reemplaza pipeline.py --status)")
    p_resume.set_defaults(func=_cmd_resume)

    # ── maaqc collect ────────────────────────────────────────────────
    p_collect = subparsers.add_parser("collect", help="Recolecta y consolida resultados de uno o más estados de protonación")
    p_collect.add_argument("--states", type=str, required=True, metavar="A,B,C,D,R",
                            help="Lista separada por comas de estados/moléculas a recolectar")
    p_collect.set_defaults(func=_cmd_collect)

    # ── maaqc report ─────────────────────────────────────────────────
    p_report = subparsers.add_parser("report", help="Genera reporte (tablas, espectros, comparación ADC2 vs TDDFT)")
    p_report.add_argument("--format", choices=["html", "md", "csv"], default="html")
    p_report.set_defaults(func=_cmd_report)

    # ── maaqc doctor ─────────────────────────────────────────────────
    p_doctor = subparsers.add_parser("doctor", help="Verifica binarios ORCA/CREST/CENSO, versiones y PATH")
    p_doctor.set_defaults(func=_cmd_doctor)

    return parser


def _cmd_run(args: argparse.Namespace) -> int:
    raise NotImplementedError(
        "maaqc run: pendiente de orchestrator/dag.py + al menos un Stage migrado. "
        "Ver docs/troubleshooting.md, sección 'Orden de migración'."
    )


def _cmd_resume(args: argparse.Namespace) -> int:
    raise NotImplementedError(
        "maaqc resume: pendiente de orchestrator/state.py "
        "(generaliza recover_prescreening.py y pipeline.py:_censo_prescreen_recover)."
    )


def _cmd_collect(args: argparse.Namespace) -> int:
    raise NotImplementedError(
        "maaqc collect: pendiente de migrar analysis/results_collector.py:collect_results()."
    )


def _cmd_report(args: argparse.Namespace) -> int:
    raise NotImplementedError(
        "maaqc report: pendiente de migrar analysis/results_collector.py (print_summary, "
        "compare_protonation_states, plot_*)."
    )


def _cmd_doctor(args: argparse.Namespace) -> int:
    raise NotImplementedError(
        "maaqc doctor: migrar pipeline.py:_resolve_bin/check_required_bins (líneas 93-158) "
        "a system/platform_detect.py — es el subcomando más simple, buen primer PR de CLI real."
    )


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv if argv is not None else sys.argv[1:])
    try:
        sys.exit(args.func(args) or 0)
    except NotImplementedError as e:
        parser.exit(status=1, message=f"[no implementado todavía] {e}\n")


if __name__ == "__main__":
    main()
