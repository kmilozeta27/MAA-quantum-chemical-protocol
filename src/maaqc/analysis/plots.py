"""Visualización: espectros UV-Vis, comparaciones, distribución de Boltzmann.

Migrar desde results_collector.py:
  - _build_spectrum()                 líneas 1539-1559
  - plot_spectrum()                   líneas 1560-1602
  - plot_comparison()                 líneas 1603-1718
  - plot_adc2_tddft_comparison()      líneas 3130-3259

y desde protonation_cycle.py:
  - plot_distribution()               líneas 1919-1984
  - plot_distribution_5()             líneas 1985-2049
  - plot_conformer_boltzmann()        líneas 2050-2108

Decisión (Issue C5): matplotlib queda como dependencia OPCIONAL, no dura.
En pyproject.toml de maaqc:

    [project.optional-dependencies]
    plots = ["matplotlib>=3.7"]

y se preserva el try/except perezoso del legacy dentro de cada función
(no se mueve el import al tope del módulo), para que `maaqc collect` /
`maaqc report --format json` sigan funcionando en corridas headless sin
X11/matplotlib instalado. `maaqc doctor` puede reportar si el extra
`plots` está disponible, igual que reporta binarios ORCA/CREST/xtb.

Estado de migración:
  - plot_spectrum(), plot_comparison(), plot_adc2_tddft_comparison():
    pendientes, dependen de C3 (results_collector.py, no incluido acá).
  - plot_distribution(), plot_distribution_5(), plot_conformer_boltzmann():
    migradas desde protonation_cycle.py:1919-2108, lógica sin modificar.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


def plot_spectrum(df_tddft: pd.DataFrame, df_conf: pd.DataFrame, mol_name: str, **kwargs):
    raise NotImplementedError("Migrar desde results_collector.py:plot_spectrum (1560-1602)")


def plot_comparison(*args, **kwargs):
    raise NotImplementedError("Migrar desde results_collector.py:plot_comparison (1603-1718)")


def plot_adc2_tddft_comparison(*args, **kwargs):
    raise NotImplementedError(
        "Migrar desde results_collector.py:plot_adc2_tddft_comparison (3130-3259)"
    )


def plot_distribution(
    df_dist: pd.DataFrame,
    pKa1: Optional[float],
    pKa2: Optional[float],
    output_path: Optional[Path] = None,
    ph_highlight: Optional[float] = None,
):
    """
    Diagrama de distribución de especies (fracciones molares vs pH).
    Guarda como PNG si se provee output_path; si no, muestra en pantalla.

    Migrado desde protonation_cycle.py:1919-1984, lógica sin modificar.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: pip install maaqc[plots]")
        return

    has_dianion = "f_dianion" in df_dist.columns
    pH_vals = df_dist["pH"].values

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = ["#2196F3", "#4CAF50", "#F44336"]  # azul, verde, rojo
    labels = ["AH₂ (q=0)", "AH⁻ (q=−1)", "A²⁻ (q=−2)"]
    cols = ["f_neutral", "f_anion"] + (["f_dianion"] if has_dianion else [])

    for col, label, color in zip(cols, labels, colors):
        if col in df_dist.columns:
            ax.plot(pH_vals, df_dist[col].values * 100, label=label,
                    color=color, linewidth=2.0)

    # Líneas verticales en pKa
    if pKa1:
        ax.axvline(pKa1, color="#2196F3", linestyle="--", alpha=0.6,
                   label=f"pKa₁ = {pKa1:.2f}")
    if pKa2:
        ax.axvline(pKa2, color="#4CAF50", linestyle="--", alpha=0.6,
                   label=f"pKa₂ = {pKa2:.2f}")

    # pH oceánico antártico
    if ph_highlight:
        ax.axvline(ph_highlight, color="gray", linestyle=":", alpha=0.8,
                   label=f"pH = {ph_highlight} (referencia)")
        ax.axvspan(ph_highlight - 0.1, ph_highlight + 0.1,
                   alpha=0.08, color="gray")

    ax.set_xlabel("pH", fontsize=12)
    ax.set_ylabel("Fracción molar (%)", fontsize=12)
    ax.set_title("Distribución de estados de protonación\n"
                 "(Henderson-Hasselbalch, ciclo de Born)", fontsize=11)
    ax.set_ylim(-2, 102)
    ax.set_xlim(pH_vals[0], pH_vals[-1])
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g%%"))
    ax.legend(fontsize=9, loc="center right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Gráfico guardado en: {output_path}")
    else:
        plt.show()
    plt.close(fig)


def plot_distribution_5(
    df_dist: pd.DataFrame,
    pKa0: Optional[float],
    pKa1: Optional[float],
    pKa2: Optional[float],
    output_path: Optional[Path] = None,
    ph_highlight: Optional[float] = None,
):
    """Diagrama de distribución de 5 especies (catión, R, D, anión, dianión) vs pH.

    Migrado desde protonation_cycle.py:1985-2049, lógica sin modificar.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: pip install maaqc[plots]")
        return

    pH_vals = df_dist["pH"].values
    fig, ax = plt.subplots(figsize=(9, 5.5))

    series = []
    if "f_cation" in df_dist.columns:
        series.append(("f_cation", "A⁺ catión (q=+1)", "#9C27B0", "-"))
    series.append(("f_R", "R neutro (q=0)", "#2196F3", "-"))
    series.append(("f_D", "D zwitterión (q=0)", "#00BCD4", "--"))
    if "f_anion" in df_dist.columns:
        series.append(("f_anion", "B⁻ anión (q=−1)", "#4CAF50", "-"))
    if "f_dianion" in df_dist.columns:
        series.append(("f_dianion", "C²⁻ dianión (q=−2)", "#F44336", "-"))

    for col, label, color, ls in series:
        if col in df_dist.columns:
            ax.plot(pH_vals, df_dist[col].values * 100, label=label,
                    color=color, linewidth=2.0, linestyle=ls)

    if pKa0 is not None:
        ax.axvline(pKa0, color="#9C27B0", linestyle=":", alpha=0.6, label=f"pKa₀ = {pKa0:.2f}")
    if pKa1 is not None:
        ax.axvline(pKa1, color="#4CAF50", linestyle=":", alpha=0.6, label=f"pKa₁ = {pKa1:.2f}")
    if pKa2 is not None:
        ax.axvline(pKa2, color="#F44336", linestyle=":", alpha=0.6, label=f"pKa₂ = {pKa2:.2f}")

    if ph_highlight is not None:
        ax.axvline(ph_highlight, color="gray", linestyle=":", alpha=0.8,
                   label=f"pH = {ph_highlight} (referencia)")
        ax.axvspan(ph_highlight - 0.1, ph_highlight + 0.1, alpha=0.08, color="gray")

    ax.set_xlabel("pH", fontsize=12)
    ax.set_ylabel("Fracción molar (%)", fontsize=12)
    ax.set_title("Distribución de estados de protonación\n"
                 "catión ⇌ [neutro ⇌ zwitterión] ⇌ anión ⇌ dianión", fontsize=11)
    ax.set_ylim(-2, 102)
    ax.set_xlim(pH_vals[0], pH_vals[-1])
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g%%"))
    ax.legend(fontsize=8, loc="center right", ncol=1)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Gráfico guardado en: {output_path}")
    else:
        plt.show()
    plt.close(fig)


def plot_conformer_boltzmann(data_list: list[dict], output_path: Optional[Path] = None):
    """
    Gráfico de barras apiladas: distribución conformacional (298 K) por estado.
    Muestra los confórmeros con mayor peso en cada estado de protonación.

    Migrado desde protonation_cycle.py:2050-2108, lógica sin modificar.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
    except ImportError:
        return

    fig, axes = plt.subplots(1, len(data_list), figsize=(5 * len(data_list), 4),
                             sharey=True)
    if len(data_list) == 1:
        axes = [axes]

    for ax, data in zip(axes, data_list):
        df = data["df"]
        if df.empty or "boltz_298K" not in df.columns:
            ax.set_title(data["mol_name"])
            ax.text(0.5, 0.5, "Sin datos", ha="center", va="center",
                    transform=ax.transAxes)
            continue

        # Filtrar confs con boltz > 1%
        col_b = "boltz_298K"
        sub = df[df[col_b] > 0.01].sort_values(col_b, ascending=False)

        colors = cm.viridis(np.linspace(0.2, 0.85, len(sub)))
        bars = ax.bar(range(len(sub)), sub[col_b].values * 100,
                      color=colors, width=0.7)
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(sub["conformer"].tolist(), rotation=45, ha="right",
                           fontsize=7)
        ax.set_title(data["mol_name"], fontsize=9)
        ax.set_ylabel("Boltzmann @ 298 K (%)" if ax == axes[0] else "")
        ax.set_ylim(0, 105)

        # Anotar porcentaje sobre cada barra
        for bar, pct in zip(bars, sub[col_b].values * 100):
            if pct > 5:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                        f"{pct:.1f}%", ha="center", fontsize=7)

    fig.suptitle("Distribución conformacional (Boltzmann, 298 K)", fontsize=11)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Gráfico conformacional guardado en: {output_path}")
    else:
        plt.show()
    plt.close(fig)
