"""Comparación ADC2 vs TDDFT (el benchmark central del paper).

Migrado 1:1 desde results_collector.py (legacy):
  - collect_adc2_results()          líneas 2822-2849  → collect_adc2_results()
  - compute_adc2_lambda_max()       líneas 2850-2885  → compute_adc2_lambda_max()
  - compare_adc2_tddft()            líneas 2886-2968  → compare_adc2_tddft()
  - _boltz_weighted_lmax()          líneas 2969-2985  → _boltz_weighted_lmax()
  - print_adc2_summary()            líneas 2986-3019  → print_adc2_summary()
  - print_adc2_tddft_comparison()   líneas 3020-3116  → print_adc2_tddft_comparison()
  - export_adc2_comparison()        líneas 3117-3129  → export_adc2_comparison()

Depende de io/parsing.py:parse_adc2_excitations y de analysis/boltzmann.py
(pesos Boltzmann para el λmax ponderado) — no se reimplementó ninguno de
los dos, solo se importan.

CAMBIOS DE FIRMA respecto al legacy (para alinear con el patrón ya
adoptado en results_collector.py / C2, donde WORK_ROOT global se eliminó):
  - collect_adc2_results(mol_name)                → collect_adc2_results(workdir, mol_name=None)
        `workdir` es el dir de la molécula/estado (equivalente al viejo
        WORK_ROOT / mol_name), igual que export_adc2_mol(workdir, ...) en
        results_collector.py. `mol_name` queda opcional, solo para
        etiquetar filas si se desea (por defecto usa workdir.name).
  - compare_adc2_tddft(mols, temp_label)          → compare_adc2_tddft(run_root, mols=None, temp_label="298K")
        `run_root` es el dir que contiene los subdirs de cada estado de
        protonación (equivalente al viejo WORK_ROOT), igual que
        compare_protonation_states(run_root, ...) en results_collector.py.
        Internamente se llama collect_results(run_root / mol) (nueva firma
        con workdir explícito) en vez de collect_results(mol) (legacy).
  - export_adc2_comparison(df_cmp, fmt)           → export_adc2_comparison(run_root, df_cmp, fmt)
        Mismo patrón que export_comparison(run_root, df_cmp, fmt) ya
        migrado: el archivo de salida se escribe en run_root/, no en un
        WORK_ROOT global.
  - compute_adc2_lambda_max() y _boltz_weighted_lmax() no tocan filesystem
    en el legacy tampoco, así que quedan con la misma firma.
  - print_adc2_summary() y print_adc2_tddft_comparison() son solo
    impresión sobre DataFrames ya cargados, misma firma; se ajustó el
    mensaje de "ruta esperada" en print_adc2_summary() para no asumir un
    prefijo "runs/" fijo (ya no existe ese global).

⚠️ TODO(kmilozeta): ADC2_SUBDIR ("07_adc2") era una constante de módulo
local en el results_collector.py legacy (línea 2707), y en pipeline.py
aparece siempre como literal "07_adc2" (nunca se centralizó ahí tampoco).
La dejo definida localmente acá por ahora, igual que en el legacy —a
diferencia de las constantes de boltzmann.py (TEMPS_K, etc.) esta no se
comparte con otro módulo migrado hasta ahora, así que no hay riesgo de
duplicación/desincronización. Si en algún momento aparece un
`stages/censo.py` o similar que también necesite conocer el nombre de la
subcarpeta ADC2, mover esto a `..constants` en ese momento.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from ..io.parsing import parse_adc2_excitations
from .results_collector import collect_results

# Subcarpeta de la etapa ADC2 dentro del workdir de cada molécula/estado
# (ver nota TODO arriba sobre por qué queda local a este módulo).
ADC2_SUBDIR = "07_adc2"


# ─────────────────────────────────────────────
#  RECOLECCIÓN Y λ_max ADC2
# ─────────────────────────────────────────────

def collect_adc2_results(workdir: Path, mol_name: Optional[str] = None) -> pd.DataFrame:
    """
    Recopila todas las excitaciones ADC(2) de los confórmeros de una molécula.

    Estructura esperada:
        <workdir>/07_adc2/conf_<id>/adc2.out

    Parámetros
    ----------
    workdir   : directorio de la molécula/estado (antes WORK_ROOT / mol_name).
    mol_name  : etiqueta opcional para logging; por defecto workdir.name.

    Retorna DataFrame con columnas:
        conformer, state, E_eV, E_cm1, wl_nm, fosc, fosc_right, fosc_left, source
    O DataFrame vacío si no hay archivos.
    """
    adc2_dir = workdir / ADC2_SUBDIR
    if not adc2_dir.exists():
        return pd.DataFrame()

    rows = []
    for conf_dir in sorted(adc2_dir.glob("conf_*")):
        if not conf_dir.is_dir():
            continue
        conf_id = conf_dir.name
        adc2_out = conf_dir / "adc2.out"
        for exc in parse_adc2_excitations(adc2_out):
            rows.append({"conformer": conf_id, **exc})

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def compute_adc2_lambda_max(
    df_adc2: pd.DataFrame,
    df_conf: pd.DataFrame,
    temp_label: str = "298K",
) -> Optional[float]:
    """
    Calcula λ_max ADC(2) Boltzmann-promediado (igual que en TDDFT).

    Usa la columna boltz_{temp_label} de df_conf.
    Sólo considera el estado S1 (state == 1).

    Retorna λ_max en nm o None si no hay datos suficientes.
    """
    if df_adc2.empty or df_conf.empty:
        return None

    col_b = f"boltz_{temp_label}"
    if col_b not in df_conf.columns:
        return None

    s1 = df_adc2[df_adc2["state"] == 1].copy()
    if s1.empty:
        return None

    s1 = s1.merge(df_conf[["conformer", col_b]], on="conformer", how="left")
    w = s1[col_b].fillna(0)
    if w.sum() <= 0:
        return None

    return float((s1["wl_nm"] * w).sum() / w.sum())


# ─────────────────────────────────────────────
#  COMPARACIÓN ADC2 vs TDDFT (ΔTDDFT-ADC2)
# ─────────────────────────────────────────────

def compare_adc2_tddft(
    run_root: Path,
    mols: Optional[list[str]] = None,
    temp_label: str = "298K",
) -> pd.DataFrame:
    """
    Compara λ_max(TD-DFT) vs λ_max(ADC2) entre estados de protonación.

    Para cada molécula disponible bajo run_root/ que tenga datos TD-DFT Y
    ADC2, calcula:
        lambda_tddft_nm   : λ_max Boltzmann-promediado TD-DFT (nm)
        lambda_adc2_nm    : λ_max Boltzmann-promediado ADC2 (nm)
        delta_nm          : λ_TDDFT - λ_ADC2 (nm)
                            > 0 → TDDFT sobreestima λ (redshift)
                            < 0 → TDDFT subestima λ (blueshift)
        delta_eV          : conversión a eV  (Δ positivo = TDDFT demasiado bajo en energía)
        n_confs_adc2      : número de confórmeros con ADC2 completado
        n_confs_tddft     : número de confórmeros con TDDFT disponible

    También imprime un análisis de sistematicidad (vía
    print_adc2_tddft_comparison):
        - ¿Es Δ constante entre estados? (std_delta → corrección sistemática viable)
        - Corrección empírica recomendada: shift = mean(Δ) ± std(Δ)

    Parámetros
    ----------
    run_root    : directorio que contiene los subdirs de cada estado de
                  protonación (antes WORK_ROOT).
    mols        : lista de nombres de mol (subdirs de run_root). None → todos
                  los que tengan 04_sp/.
    temp_label  : temperatura para el Boltzmann (default "298K").

    Retorna DataFrame con una fila por molécula.
    """
    if mols is None:
        mols = sorted(p.name for p in run_root.iterdir()
                      if p.is_dir() and (p / "04_sp").exists())

    rows = []
    for mol in mols:
        mol_dir = run_root / mol

        # ── Cargar TDDFT ───────────────────────────────────────────────────
        try:
            df_conf, df_tddft = collect_results(mol_dir)
        except Exception as e:
            warnings.warn(f"[ADC2 vs TDDFT] {mol}: error al cargar TDDFT — {e}")
            df_conf = pd.DataFrame()
            df_tddft = pd.DataFrame()

        lmax_tddft: Optional[float] = None
        n_confs_tddft = 0
        if not df_tddft.empty and not df_conf.empty:
            lmax_tddft = df_conf.attrs.get("lambda_max_boltz_298K") \
                        if temp_label == "298K" else \
                        _boltz_weighted_lmax(df_tddft, df_conf, temp_label)
            n_confs_tddft = int(df_tddft["conformer"].nunique())

        # ── Cargar ADC2 ────────────────────────────────────────────────────
        df_adc2 = collect_adc2_results(mol_dir, mol_name=mol)
        lmax_adc2: Optional[float] = None
        n_confs_adc2 = 0
        if not df_adc2.empty and not df_conf.empty:
            lmax_adc2 = compute_adc2_lambda_max(df_adc2, df_conf, temp_label)
            n_confs_adc2 = int(df_adc2["conformer"].nunique())

        # ── Δ = TDDFT − ADC2 ──────────────────────────────────────────────
        delta_nm: Optional[float] = None
        delta_eV: Optional[float] = None
        if lmax_tddft is not None and lmax_adc2 is not None and lmax_adc2 > 0:
            delta_nm = lmax_tddft - lmax_adc2
            # Conversión eV: E = hc/λ → ΔE ≈ -hc·Δλ/λ² (diferencia de energías)
            # Calculamos directamente como diferencia de energías en eV:
            ev_tddft = 1239.84193 / lmax_tddft   # eV
            ev_adc2 = 1239.84193 / lmax_adc2     # eV
            delta_eV = ev_tddft - ev_adc2        # >0 → TDDFT subestima E (redshift)

        rows.append({
            "molecule": mol,
            "lambda_tddft_nm": lmax_tddft,
            "lambda_adc2_nm": lmax_adc2,
            "delta_nm": delta_nm,        # TDDFT − ADC2 en nm
            "delta_eV": delta_eV,        # TDDFT − ADC2 en eV
            "n_confs_tddft": n_confs_tddft,
            "n_confs_adc2": n_confs_adc2,
        })

    return pd.DataFrame(rows)


def _boltz_weighted_lmax(
    df_tddft: pd.DataFrame,
    df_conf: pd.DataFrame,
    temp_label: str,
) -> Optional[float]:
    """Helper: λ_max TDDFT Boltzmann-promediado a temperatura arbitraria."""
    col_b = f"boltz_{temp_label}"
    if col_b not in df_conf.columns:
        return None
    s1 = df_tddft[df_tddft["state"] == 1].copy()
    if s1.empty:
        return None
    s1 = s1.merge(df_conf[["conformer", col_b]], on="conformer", how="left")
    w = s1[col_b].fillna(0)
    return float((s1["wavelength_nm"] * w).sum() / w.sum()) if w.sum() > 0 else None


# ─────────────────────────────────────────────
#  IMPRESIÓN
# ─────────────────────────────────────────────

def print_adc2_summary(df_adc2: pd.DataFrame, mol_name: str,
                        df_conf: pd.DataFrame):
    """Imprime tabla de excitaciones ADC2 de una molécula."""
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  ADC(2) EXCITACIONES — {mol_name}")
    print(sep)

    if df_adc2.empty:
        print(f"  Sin datos ADC2 para '{mol_name}'.")
        print(f"  Ruta esperada: <workdir>/{ADC2_SUBDIR}/conf_*/adc2.out")
        print(f"{sep}\n")
        return

    # S1 por confórmero con pesos Boltzmann
    col_b = "boltz_298K"
    s1 = df_adc2[df_adc2["state"] == 1].copy()
    if not s1.empty and not df_conf.empty and col_b in df_conf.columns:
        s1 = s1.merge(df_conf[["conformer", col_b]], on="conformer", how="left")
        disp_cols = ["conformer", "wl_nm", "E_eV", "fosc", col_b]
        disp_cols = [c for c in disp_cols if c in s1.columns]
        print("\n  ── S1 por confórmero (ADC2, ponderado 298 K) ──")
        print(s1[disp_cols].to_string(index=False,
              float_format=lambda x: f"{x:.4f}"))

    lmax = compute_adc2_lambda_max(df_adc2, df_conf, "298K")
    if lmax:
        print(f"\n  λ_max ADC2 Boltzmann-promediado (298 K, S1): {lmax:.1f} nm")

    n = df_adc2["conformer"].nunique()
    print(f"  Confórmeros ADC2 completados: {n}")
    print(f"{sep}\n")


def print_adc2_tddft_comparison(df_cmp: pd.DataFrame):
    """
    Imprime la tabla comparativa ADC2 vs TDDFT con análisis de sistematicidad.
    """
    sep = "═" * 80
    print(f"\n{sep}")
    print(f"  COMPARACIÓN ADC(2) vs TD-DFT — λ_max Boltzmann-promediado (298 K)")
    print(sep)

    if df_cmp.empty:
        print("  Sin datos para comparar.")
        return

    # ── Tabla por estado de protonación ─────────────────────────────────────
    print(f"\n  {'Molécula':<25} {'λ_TDDFT(nm)':>12} {'λ_ADC2(nm)':>12} "
          f"{'Δ(nm)':>10} {'Δ(eV)':>10} {'N_ADC2':>7}")
    print(f"  {'─'*25} {'─'*12} {'─'*12} {'─'*10} {'─'*10} {'─'*7}")

    valid_rows = []
    for _, row in df_cmp.iterrows():
        tddft = row.get("lambda_tddft_nm")
        adc2 = row.get("lambda_adc2_nm")
        delta = row.get("delta_nm")
        devv = row.get("delta_eV")
        n = row.get("n_confs_adc2", 0)

        flag = ""
        if tddft is None or adc2 is None:
            flag = "  ⚠ datos incompletos"

        print(f"  {row['molecule']:<25} "
              f"{'—':>12}" if tddft is None else
              f"  {row['molecule']:<25} {tddft:>12.1f}", end="")
        print(f" {adc2:>12.1f}" if adc2 is not None else f" {'—':>12}", end="")
        print(f" {delta:>+10.1f}" if delta is not None else f" {'—':>10}", end="")
        print(f" {devv:>+10.4f}" if devv is not None else f" {'—':>10}", end="")
        print(f" {int(n):>7}{flag}")

        if delta is not None:
            valid_rows.append(row)

    if not valid_rows:
        print(f"\n  Sin pares TDDFT+ADC2 completos para análisis de sistematicidad.")
        print(f"{sep}\n")
        return

    # ── Análisis de sistematicidad ───────────────────────────────────────────
    deltas_nm = [r["delta_nm"] for r in valid_rows]
    deltas_eV = [r["delta_eV"] for r in valid_rows]

    mean_nm = float(np.mean(deltas_nm))
    std_nm = float(np.std(deltas_nm, ddof=1)) if len(deltas_nm) > 1 else 0.0
    mean_eV = float(np.mean(deltas_eV))
    std_eV = float(np.std(deltas_eV, ddof=1)) if len(deltas_eV) > 1 else 0.0

    print(f"\n  ── Sistematicidad del desplazamiento TDDFT − ADC2 ──")
    print(f"  N pares válidos       : {len(deltas_nm)}")
    print(f"  Δ medio (nm)          : {mean_nm:+.2f} ± {std_nm:.2f}")
    print(f"  Δ medio (eV)          : {mean_eV:+.4f} ± {std_eV:.4f}")

    # Evaluación de sistematicidad: si std < |mean|/3 el desplazamiento es sistemático
    systematic = (std_nm < abs(mean_nm) / 3.0) if mean_nm != 0 else (std_nm < 2.0)

    if systematic and len(deltas_nm) > 1:
        direction = "redshift" if mean_nm > 0 else "blueshift"
        print(f"\n  ✓  Desplazamiento SISTEMÁTICO ({direction}):")
        print(f"     TD-DFT predice λ_max {abs(mean_nm):.1f} nm más {direction.replace('shift','ed')}"
              f" que ADC2 (CV = {std_nm/abs(mean_nm)*100:.1f}%)")
        print(f"\n  ► Corrección empírica recomendada:")
        print(f"     λ_corr(nm) = λ_TDDFT(nm) − ({mean_nm:+.1f} nm)")
        print(f"     E_corr(eV) = E_TDDFT(eV) − ({mean_eV:+.4f} eV)")
        print(f"     Incertidumbre: ±{std_nm:.1f} nm  (±{std_eV:.4f} eV)")
    elif len(deltas_nm) == 1:
        print(f"\n  ⚠  Solo 1 par: imposible evaluar sistematicidad.")
        print(f"     Desplazamiento observado: {mean_nm:+.1f} nm ({mean_eV:+.4f} eV)")
    else:
        print(f"\n  ⚠  Desplazamiento NO sistemático (std/|mean| = "
              f"{std_nm/abs(mean_nm)*100:.1f}% > 33%).")
        print(f"     La corrección empírica tiene alta incertidumbre (±{std_nm:.1f} nm).")
        print(f"     Considerar revisar geometrías o ampliar ensemble conformacional.")

    # ADC2 inter-estados: comparación relativa
    if len(valid_rows) > 1:
        lmax_adc2_vals = [r["lambda_adc2_nm"] for r in valid_rows]
        mols_valid = [r["molecule"] for r in valid_rows]
        print(f"\n  ── Comparación λ_max ADC2 entre estados de protonación ──")
        adc2_min = min(lmax_adc2_vals)
        adc2_max = max(lmax_adc2_vals)
        mol_min = mols_valid[lmax_adc2_vals.index(adc2_min)]
        mol_max = mols_valid[lmax_adc2_vals.index(adc2_max)]
        print(f"  Rango ADC2: {adc2_min:.1f} – {adc2_max:.1f} nm")
        print(f"  Estado con λ_max ADC2 menor: {mol_min} ({adc2_min:.1f} nm)")
        print(f"  Estado con λ_max ADC2 mayor: {mol_max} ({adc2_max:.1f} nm)")

    print(f"\n{sep}\n")


# ─────────────────────────────────────────────
#  EXPORTACIÓN
# ─────────────────────────────────────────────

def export_adc2_comparison(run_root: Path, df_cmp: pd.DataFrame, fmt: str):
    """Exporta la tabla comparativa ADC2 vs TDDFT a CSV y/o XLSX."""
    if df_cmp.empty:
        return
    out = run_root / "comparison_adc2_tddft"
    if fmt in ("csv", "both"):
        df_cmp.to_csv(out.with_suffix(".csv"), index=False)
        print(f"  CSV ADC2 vs TDDFT → {out.with_suffix('.csv')}")
    if fmt in ("xlsx", "both"):
        df_cmp.to_excel(out.with_suffix(".xlsx"), index=False)
        print(f"  XLSX ADC2 vs TDDFT → {out.with_suffix('.xlsx')}")
