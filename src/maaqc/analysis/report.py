"""Generación de reportes (HTML/MD/JSON) de una corrida `maaqc collect`/`report`.

Interfaz NUEVA — no migrada del legacy. `pipeline.py`/`results_collector.py`
(legacy) solo tenían `print_summary`/`print_comparison` a stdout; nunca
existió un formato de reporte serializable. Este módulo lo define desde
cero, siguiendo el mismo criterio que el resto de `analysis/`:

  - Recibe `run_root: Path` explícito (mismo patrón que
    `compare_protonation_states`), nunca una constante global.
  - NO reimplementa recolección/comparación: reutiliza
    `results_collector.collect_results()` / `.compare_protonation_states()`
    tal cual están, y solo los serializa en distintos formatos.

Contrato:
    build_report_data(run_root, states, ref_mol=None) -> ReportData
    render_html(data) -> str
    render_markdown(data) -> str
    render_json(data) -> str
    write_report(data, workdir, fmt) -> Path   # ruta del archivo escrito

Si más adelante se necesitan reportes por-molécula individual (no solo la
comparación entre estados), extender `ReportData` con lo que haga falta en
vez de introducir una segunda clase paralela.

ADC2 (ver skill maa-repackaging, tarea B1): SÍ entra en `ReportData`, como
sección nueva (`df_adc2_cmp`), en vez de quedar fuera de alcance de
`report`. Se decidió así porque `compare_adc2_tddft()` ya devuelve un
DataFrame por-estado con la misma forma que `df_cmp` (una fila por
molécula/estado) — es una comparación más entre estados, no un artefacto
por-confórmero como `dfs_conf`/`dfs_tddft`, así que encaja en el mismo
nivel que la sección "Comparación entre estados" ya existente en vez de
justificar una salida aparte de `collect`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import pandas as pd

from . import adc2_compare as adc2
from . import results_collector as rc


@dataclass
class ReportData:
    run_root: Path
    states: list[str]
    ref_mol: Optional[str]
    df_cmp: pd.DataFrame
    dfs_conf: dict[str, pd.DataFrame] = field(default_factory=dict)
    dfs_tddft: dict[str, pd.DataFrame] = field(default_factory=dict)
    # Comparación ADC2 vs TDDFT (compare_adc2_tddft), una fila por estado.
    # Vacío (DataFrame sin filas) si ningún estado corrió 07_adc2/ todavía
    # — ADC2 es la etapa final de la cadena, así que esto es normal en
    # corridas parciales, no un error.
    df_adc2_cmp: pd.DataFrame = field(default_factory=pd.DataFrame)


# Columnas de detalle por confórmero que se muestran en MD/HTML/JSON.
# Mismo subconjunto que usa print_summary() para la tabla base, para que
# el reporte no diverja visualmente de lo que ya se ve con `maaqc collect`.
_CONFORMER_DISPLAY_COLS = (
    "conformer", "E_hartree", "dE_kcal", "boltz_298K",
    "has_imaginary", "n_imaginary",
)


def build_report_data(
    run_root: Path, states: list[str], ref_mol: Optional[str] = None,
) -> ReportData:
    """Arma los datos del reporte a partir de results_collector.py.

    `states` ya viene resuelto/validado por el caller (cli.py:
    `_resolve_states`) — acá no se re-valida contra el filesystem.
    """
    df_cmp, dfs_conf, dfs_tddft = rc.compare_protonation_states(run_root, ref_mol=ref_mol)
    if not df_cmp.empty:
        df_cmp = df_cmp[df_cmp["molecule"].isin(states)].reset_index(drop=True)
    dfs_conf = {k: v for k, v in dfs_conf.items() if k in states}
    dfs_tddft = {k: v for k, v in dfs_tddft.items() if k in states}

    # compare_adc2_tddft() ya acepta `mols` explícito (mismo patrón que
    # collect_results(workdir) con workdir explícito), así que se le pasan
    # directo los `states` resueltos en vez de filtrar después como con
    # df_cmp arriba.
    df_adc2_cmp = adc2.compare_adc2_tddft(run_root, mols=states)

    return ReportData(
        run_root=run_root, states=states, ref_mol=ref_mol,
        df_cmp=df_cmp, dfs_conf=dfs_conf, dfs_tddft=dfs_tddft,
        df_adc2_cmp=df_adc2_cmp,
    )


def _conformer_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in _CONFORMER_DISPLAY_COLS if c in df.columns]


def render_markdown(data: ReportData) -> str:
    lines = [f"# Reporte maaqc — {data.run_root.name}", ""]
    header = f"Estados: {', '.join(data.states)}"
    if data.ref_mol:
        header += f"  (ref: {data.ref_mol})"
    lines += [header, "", "## Comparación entre estados", ""]

    if data.df_cmp.empty:
        lines.append("_Sin datos de comparación._")
    else:
        lines.append("```")
        lines.append(data.df_cmp.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        lines.append("```")
    lines.append("")

    lines += ["## Comparación ADC2 vs TDDFT", ""]
    if data.df_adc2_cmp.empty:
        lines.append("_Sin datos ADC2 (ningún estado corrió 07_adc2/ todavía)._")
    else:
        lines.append("```")
        lines.append(data.df_adc2_cmp.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
        lines.append("```")
    lines.append("")

    for state_name, df in data.dfs_conf.items():
        lines.append(f"## Estado {state_name}")
        lines.append("")
        if df.empty:
            lines.append("_Sin confórmeros SP._")
        else:
            cols = _conformer_cols(df)
            lines.append("```")
            lines.append(df[cols].to_string(index=False, float_format=lambda x: f"{x:.4f}"))
            lines.append("```")
        lines.append("")

    return "\n".join(lines)


def render_html(data: ReportData) -> str:
    header = f"Estados: {', '.join(data.states)}"
    if data.ref_mol:
        header += f" (ref: {data.ref_mol})"

    parts = [
        "<!doctype html><html><head><meta charset='utf-8'>",
        f"<title>maaqc report — {data.run_root.name}</title>",
        "<style>body{font-family:sans-serif;margin:2rem} "
        "table{border-collapse:collapse} "
        "td,th{border:1px solid #ccc;padding:4px 8px;text-align:right}</style>",
        "</head><body>",
        f"<h1>Reporte maaqc — {data.run_root.name}</h1>",
        f"<p>{header}</p>",
        "<h2>Comparación entre estados</h2>",
    ]

    if data.df_cmp.empty:
        parts.append("<p><em>Sin datos de comparación.</em></p>")
    else:
        parts.append(data.df_cmp.to_html(index=False, float_format=lambda x: f"{x:.3f}"))

    parts.append("<h2>Comparación ADC2 vs TDDFT</h2>")
    if data.df_adc2_cmp.empty:
        parts.append("<p><em>Sin datos ADC2 (ningún estado corrió 07_adc2/ todavía).</em></p>")
    else:
        parts.append(data.df_adc2_cmp.to_html(index=False, float_format=lambda x: f"{x:.3f}"))

    for state_name, df in data.dfs_conf.items():
        parts.append(f"<h2>Estado {state_name}</h2>")
        if df.empty:
            parts.append("<p><em>Sin confórmeros SP.</em></p>")
            continue
        cols = _conformer_cols(df)
        parts.append(df[cols].to_html(index=False, float_format=lambda x: f"{x:.4f}"))

    parts.append("</body></html>")
    return "\n".join(parts)


def render_json(data: ReportData) -> str:
    payload = {
        "run_root": str(data.run_root),
        "states": data.states,
        "reference": data.ref_mol,
        "comparison": data.df_cmp.to_dict(orient="records"),
        "comparison_adc2_tddft": data.df_adc2_cmp.to_dict(orient="records"),
        "conformers": {
            state_name: (
                df.drop(columns=["_thermochem"], errors="ignore")
                .to_dict(orient="records")
            )
            for state_name, df in data.dfs_conf.items()
        },
    }
    # default=str cubre np.float64/NaN/Path residuales sin reventar el dump.
    return json.dumps(payload, indent=2, default=str)


_RENDERERS = {
    "html": (render_html, ".html"),
    "md": (render_markdown, ".md"),
    "json": (render_json, ".json"),
}


def write_report(data: ReportData, workdir: Path, fmt: str) -> Path:
    """Renderiza `data` en `fmt` y lo escribe como `workdir/report.<ext>`.

    Levanta ValueError si `fmt` no es uno de los soportados — el caller
    (cli.py) ya restringe `--format` con `choices=`, así que esto es
    defensivo (por si `write_report` se llama desde otro lado sin pasar
    por argparse).
    """
    if fmt not in _RENDERERS:
        raise ValueError(f"formato de reporte no soportado: {fmt!r}")
    render_fn, ext = _RENDERERS[fmt]
    content = render_fn(data)
    out_path = workdir / f"report{ext}"
    out_path.write_text(content, encoding="utf-8")
    return out_path
