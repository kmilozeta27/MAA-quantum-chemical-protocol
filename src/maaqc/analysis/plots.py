"""Visualización: espectros UV-Vis, comparaciones, distribución de Boltzmann.

Migrado desde results_collector.py:
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

CAMBIOS DE FIRMA respecto al legacy (A2, alineado con el patrón ya
adoptado en results_collector.py/C2 y adc2_compare.py/A1, donde el
WORK_ROOT global se eliminó):
  - plot_spectrum(df_tddft, df_conf, mol_name, sigma, x_range,
    output_path=None)
        Antes escribía directo a WORK_ROOT / mol_name / "uvvis_spectrum.png"
        y siempre hacía plt.show(). Ahora sigue el mismo contrato que
        plot_distribution()/plot_distribution_5(): si se pasa output_path
        guarda el PNG ahí (el caller decide la ruta, típicamente
        workdir / "uvvis_spectrum.png"); si no, muestra en pantalla.
  - plot_comparison(dfs_conf, dfs_tddft, sigma, x_range, output_path=None)
        Antes escribía a WORK_ROOT / "comparison_uvvis.png". Mismo cambio:
        output_path explícito (típicamente run_root / "comparison_uvvis.png"),
        opcional.
  - plot_adc2_tddft_comparison(df_cmp, dfs_conf, dfs_tddft, sigma, x_range,
    output_path=None)
        Antes escribía a WORK_ROOT / "comparison_adc2_tddft.png". Mismo
        cambio. El shape de df_cmp es el que retorna
        analysis/adc2_compare.py:compare_adc2_tddft() (A1): columnas
        molecule, lambda_tddft_nm, lambda_adc2_nm, delta_nm, delta_eV,
        n_confs_tddft, n_confs_adc2 — no se tocó la lógica de graficado,
        que ya consumía exactamente esas columnas en el legacy.
  - Ninguna de las tres vuelve a hacer plt.show() incondicional cuando se
    guarda a archivo (antes el legacy hacía savefig() Y show() siempre;
    show() en un run headless sin backend interactivo revienta o cuelga).
    Ahora: si output_path se especifica, solo savefig(); si no, solo show().

Estado de migración:
  - plot_spectrum(), plot_comparison(), plot_adc2_tddft_comparison():
    migradas desde results_collector.py:1539-1718, 3130-3259. Lógica de
    graficado sin modificar; solo se adaptó la escritura de salida (ver
    "CAMBIOS DE FIRMA" arriba) y se eliminó la dependencia del TEMP_LABELS
    global de results_collector.py (ahora se importa desde ..constants,
    centralizado en D0).
  - plot_distribution(), plot_distribution_5(), plot_conformer_boltzmann():
    migradas desde protonation_cycle.py:1919-2108, lógica sin modificar.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..constants import TEMP_LABELS


def _build_spectrum(
    df_tddft: pd.DataFrame,
    df_conf: pd.DataFrame,
    sigma: float,
    x: np.ndarray,
    temp_label: str = "298K",
) -> np.ndarray:
    """Construye el espectro Gaussiano-suavizado Boltzmann-promediado.

    Migrado desde results_collector.py:_build_spectrum (1539-1559),
    lógica sin modificar.
    """
    col = f"boltz_{temp_label}"
    if col not in df_conf.columns:
        col = "boltz_298K"

    y = np.zeros_like(x)
    if df_tddft.empty:
        return y

    df = df_tddft.merge(df_conf[["conformer", col]], on="conformer", how="left")
    df[col] = df[col].fillna(1.0 / max(len(df_tddft), 1))

    for _, row in df.iterrows():
        gauss = np.exp(-0.5 * ((x - row["wavelength_nm"]) / sigma) ** 2)
        y += row["fosc"] * row[col] * gauss
    return y


def plot_spectrum(
    df_tddft: pd.DataFrame,
    df_conf: pd.DataFrame,
    mol_name: str,
    sigma: float = 20.0,
    x_range: tuple = (200, 500),
    output_path: Optional[Path] = None,
):
    """Espectro UV-Vis individual: Boltzmann-promediado a múltiples temperaturas.

    Migrado desde results_collector.py:plot_spectrum (1560-1602), lógica
    de graficado sin modificar. Cambio de firma: output_path explícito
    en vez de WORK_ROOT / mol_name / "uvvis_spectrum.png" (ver docstring
    del módulo).
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: pip install maaqc[plots]")
        return

    if df_tddft.empty:
        print("Sin datos TDDFT para graficar.")
        return

    x = np.linspace(x_range[0], x_range[1], 1500)
    colors = ["#0077bb", "#33bbee", "#ee7733", "#cc3311"]

    fig, ax = plt.subplots(figsize=(9, 4.5))

    for lbl, color in zip(TEMP_LABELS, colors):
        y = _build_spectrum(df_tddft, df_conf, sigma, x, temp_label=lbl)
        ax.plot(x, y, lw=1.8, label=f"T = {lbl}", color=color)

    # Palillos de excitaciones individuales (ponderados, 298 K)
    col298 = "boltz_298K"
    df_sticks = df_tddft.merge(df_conf[["conformer", col298]], on="conformer", how="left")
    df_sticks[col298] = df_sticks[col298].fillna(0)
    ax.vlines(df_sticks["wavelength_nm"], 0,
              df_sticks["fosc"] * df_sticks[col298],
              colors="gray", alpha=0.35, lw=0.8)

    ax.set_xlabel("Longitud de onda (nm)", fontsize=12)
    ax.set_ylabel("Intensidad (u.a.)", fontsize=12)
    ax.set_title(f"UV-Vis  —  {mol_name}  (σ = {sigma} nm)", fontsize=13)
    ax.legend(fontsize=10)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Espectro guardado en: {output_path}")
    else:
        plt.show()
    plt.close(fig)


def plot_comparison(
    dfs_conf: dict,
    dfs_tddft: dict,
    sigma: float = 20.0,
    x_range: tuple = (200, 500),
    output_path: Optional[Path] = None,
):
    """
    Superpone espectros UV-Vis (298 K) de todos los estados de protonación,
    normalizados al máximo individual para facilitar comparación de forma.

    Migrado desde results_collector.py:plot_comparison (1603-1718), lógica
    de graficado sin modificar. Cambio de firma: output_path explícito en
    vez de WORK_ROOT / "comparison_uvvis.png" (ver docstring del módulo).
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: pip install maaqc[plots]")
        return

    x = np.linspace(x_range[0], x_range[1], 1500)

    # Paleta colorblind-friendly (Paul Tol)
    palette = [
        "#4477AA", "#EE6677", "#228833", "#CCBB44",
        "#66CCEE", "#AA3377", "#BBBBBB"
    ]

    fig, axes = plt.subplots(2, 1, figsize=(10, 8),
                             gridspec_kw={"height_ratios": [3, 1]})
    ax_spec, ax_bar = axes

    mol_names = sorted(dfs_conf.keys())
    lambda_maxs = []

    for i, mol in enumerate(mol_names):
        df_t = dfs_tddft.get(mol, pd.DataFrame())
        df_c = dfs_conf.get(mol, pd.DataFrame())
        if df_t.empty or df_c.empty:
            continue

        y = _build_spectrum(df_t, df_c, sigma, x, temp_label="298K")
        y_norm = y / y.max() if y.max() > 0 else y

        color = palette[i % len(palette)]
        ax_spec.plot(x, y_norm, lw=2, label=mol.replace("_", " "), color=color)

        # Marcar λ_max
        lmax = x[np.argmax(y_norm)]
        lambda_maxs.append((mol, lmax, color))
        ax_spec.axvline(lmax, color=color, lw=0.7, ls="--", alpha=0.6)

    ax_spec.set_xlabel("Longitud de onda (nm)", fontsize=12)
    ax_spec.set_ylabel("Intensidad normalizada", fontsize=12)
    ax_spec.set_title(f"Comparación UV-Vis — Estados de protonación  "
                      f"(σ = {sigma} nm, T = 298 K)", fontsize=13)
    ax_spec.legend(fontsize=9, loc="upper right")

    # Gráfico de barras de λ_max
    mols_bar  = [m for m, _, _ in lambda_maxs]
    lmax_bar  = [l for _, l, _ in lambda_maxs]
    cols_bar  = [c for _, _, c in lambda_maxs]
    y_pos = range(len(mols_bar))
    ax_bar.barh(y_pos, lmax_bar, color=cols_bar, alpha=0.85)
    ax_bar.set_yticks(list(y_pos))
    ax_bar.set_yticklabels([m.replace("_", " ") for m in mols_bar], fontsize=9)
    ax_bar.set_xlabel("λ_max (nm)", fontsize=11)
    ax_bar.set_title("λ_max Boltzmann-promediado (S1, 298 K)", fontsize=11)

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Espectro comparativo guardado en: {output_path}")
    else:
        plt.show()
    plt.close(fig)


def plot_adc2_tddft_comparison(
    df_cmp: pd.DataFrame,
    dfs_conf: dict,
    dfs_tddft: dict,
    sigma: float = 20.0,
    x_range: tuple = (200, 500),
    output_path: Optional[Path] = None,
):
    """
    Grafica la comparación ADC2 vs TDDFT mostrando para cada estado:
      - Espectro TDDFT (gaussiano, normalizado)
      - Línea vertical λ_max TDDFT
      - Línea vertical λ_max ADC2
      - Flechas mostrando el desplazamiento Δ

    Panel inferior: barra de Δ(nm) por estado de protonación.

    Migrado desde results_collector.py:plot_adc2_tddft_comparison
    (3130-3259), lógica de graficado sin modificar. df_cmp tiene el shape
    que retorna analysis/adc2_compare.py:compare_adc2_tddft() (A1):
    columnas molecule, lambda_tddft_nm, lambda_adc2_nm, delta_nm,
    delta_eV, n_confs_tddft, n_confs_adc2 — coincide 1:1 con lo que este
    plot ya esperaba en el legacy. Cambio de firma: output_path explícito
    en vez de WORK_ROOT / "comparison_adc2_tddft.png" (ver docstring del
    módulo).
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: pip install maaqc[plots]")
        return

    if df_cmp.empty:
        return

    # Filtrar solo filas con datos completos
    valid = df_cmp.dropna(subset=["lambda_tddft_nm", "lambda_adc2_nm"])
    if valid.empty:
        print("  Sin pares TDDFT+ADC2 para graficar.")
        return

    palette = [
        "#4477AA", "#EE6677", "#228833", "#CCBB44",
        "#66CCEE", "#AA3377", "#BBBBBB"
    ]

    fig, axes = plt.subplots(
        2, 1, figsize=(11, 9),
        gridspec_kw={"height_ratios": [3, 1.2]}
    )
    ax_spec, ax_delta = axes

    x = np.linspace(x_range[0], x_range[1], 1500)

    for i, (_, row) in enumerate(valid.iterrows()):
        mol   = row["molecule"]
        color = palette[i % len(palette)]

        df_t = dfs_tddft.get(mol, pd.DataFrame())
        df_c = dfs_conf.get(mol, pd.DataFrame())

        # Espectro TDDFT
        if not df_t.empty and not df_c.empty:
            y = _build_spectrum(df_t, df_c, sigma, x, temp_label="298K")
            y_norm = y / y.max() if y.max() > 0 else y
            ax_spec.plot(x, y_norm, lw=1.8, alpha=0.7,
                         label=mol.replace("_", " "), color=color)

        # λ_max TDDFT
        lmax_td = row["lambda_tddft_nm"]
        ax_spec.axvline(lmax_td, color=color, lw=1.5, ls="--", alpha=0.85,
                        label=f"λ_TDDFT {mol.replace('_', ' ')} = {lmax_td:.0f} nm")

        # λ_max ADC2
        lmax_a2 = row["lambda_adc2_nm"]
        ax_spec.axvline(lmax_a2, color=color, lw=2.0, ls=":",
                        label=f"λ_ADC2  {mol.replace('_', ' ')} = {lmax_a2:.0f} nm")

        # Flecha TDDFT → ADC2
        y_arrow = 0.85 - i * 0.12
        if 0 < y_arrow < 1:
            ax_spec.annotate(
                "", xy=(lmax_a2, y_arrow), xytext=(lmax_td, y_arrow),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.4)
            )
            mid = (lmax_td + lmax_a2) / 2
            ax_spec.text(mid, y_arrow + 0.03,
                         f"Δ={row['delta_nm']:+.1f}nm",
                         ha="center", va="bottom", fontsize=8, color=color)

    ax_spec.set_xlabel("Longitud de onda (nm)", fontsize=12)
    ax_spec.set_ylabel("Intensidad TD-DFT (norm.)", fontsize=12)
    ax_spec.set_title(
        f"TD-DFT vs ADC(2) — λ_max por estado de protonación  "
        f"(σ = {sigma} nm, T = 298 K)", fontsize=12
    )
    # Leyenda compacta: solo las líneas de espectro
    handles = [
        plt.Line2D([0], [0], color=palette[i % len(palette)], lw=1.8, alpha=0.7,
                   label=row["molecule"].replace("_", " "))
        for i, (_, row) in enumerate(valid.iterrows())
    ]
    handles += [
        plt.Line2D([0], [0], color="gray", lw=1.5, ls="--", label="λ_max TD-DFT"),
        plt.Line2D([0], [0], color="gray", lw=2.0, ls=":", label="λ_max ADC(2)"),
    ]
    ax_spec.legend(handles=handles, fontsize=8, loc="upper right")

    # Panel inferior: Δ(nm) por estado
    mols_bar  = valid["molecule"].tolist()
    delta_bar = valid["delta_nm"].tolist()
    colors_bar = [palette[i % len(palette)] for i in range(len(mols_bar))]
    y_pos = range(len(mols_bar))

    bars = ax_delta.barh(list(y_pos), delta_bar, color=colors_bar, alpha=0.85)
    ax_delta.axvline(0, color="black", lw=0.8, ls="-")

    # Añadir valor numérico en cada barra
    for bar, val in zip(bars, delta_bar):
        xpos = val + 0.3 if val >= 0 else val - 0.3
        ha   = "left" if val >= 0 else "right"
        ax_delta.text(xpos, bar.get_y() + bar.get_height() / 2,
                      f"{val:+.1f} nm", va="center", ha=ha, fontsize=9)

    ax_delta.set_yticks(list(y_pos))
    ax_delta.set_yticklabels([m.replace("_", " ") for m in mols_bar], fontsize=9)
    ax_delta.set_xlabel("Δλ = λ_TDDFT − λ_ADC2  (nm)", fontsize=11)
    ax_delta.set_title("Desplazamiento TD-DFT vs ADC(2)  "
                       "[positivo = TDDFT redshifted]", fontsize=10)

    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Plot ADC2 vs TDDFT guardado en: {output_path}")
    else:
        plt.show()
    plt.close(fig)


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
