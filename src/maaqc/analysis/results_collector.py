"""Recolección, comparación e impresión/exportación de resultados por estado
de protonación.

Migrado 1:1 desde results_collector.py (legacy):
  - collect_results()              líneas 1229-1313
  - print_summary()                líneas 1320-1379
  - compare_protonation_states()   líneas 1386-1484
  - print_comparison()             líneas 1487-1532
  - export_results()               líneas 3260-3301
  - export_adc2_mol()              líneas 3304-3322
  - export_comparison()            líneas 3326-3335

Depende de C2 (analysis/boltzmann.py) terminado: los nombres de función
cambiaron de CamelCase a snake_case (`compute_G_total` → `compute_g_total`,
`compute_boltzmann_multiT` → `compute_boltzmann_multi_t`); las llamadas
abajo ya usan los nombres nuevos.

⚠️ TODO(kmilozeta) — dos dependencias sin resolver, no incluidas en el
alcance de este issue, que hay que verificar antes de que este módulo
corra:

1. `analyze_crest_to_dft()` / `print_crest_dft_analysis()` (legacy
   results_collector.py:524-642): `collect_results()` y `print_summary()`
   las llaman pero no están en la lista de funciones migradas a
   `io/parsing.py` que se mencionó en sesiones previas, ni tienen issue
   propio en las notas que tengo. Asumí un import desde `.crest_dft`
   (módulo nuevo, no confirmado) — hay que decidir dónde viven
   realmente (¿`analysis/crest_dft.py`? ¿parte de `censo.py` /
   `CensoStage`, dado que ahí es donde se compara ranking CREST vs DFT?)
   antes de que el import de abajo funcione.

2. `WORK_ROOT`: en el legacy es un `Path("runs")` global fijo. Acá lo
   importo desde `..constants` igual que en boltzmann.py, pero en la
   arquitectura nueva (orchestrator/state.py con workdir por
   molécula/estado) probablemente `collect_results` y
   `compare_protonation_states` deberían recibir el workdir como
   parámetro en vez de derivarlo de una constante global — no lo cambié
   acá para no mezclar refactor de arquitectura con migración 1:1, pero
   quedó pendiente de decidir.
"""

from __future__ import annotations

import warnings
from typing import Optional

import numpy as np
import pandas as pd

from ..constants import HARTREE_TO_KCAL, T_STD, WORK_ROOT
from ..io.parsing import (
    parse_imaginary_freqs,
    parse_opt_energy,
    parse_sp_energy,
    parse_tddft_excitations,
    parse_thermochem,
)
from .boltzmann import (
    compute_boltzmann_multi_t,
    compute_g_total,
    print_boltzmann_analysis,
)
from .crest_dft import analyze_crest_to_dft, print_crest_dft_analysis  # TODO: confirmar ubicación real


# ─────────────────────────────────────────────
#  RECOLECCIÓN DE RESULTADOS
# ─────────────────────────────────────────────

def collect_results(mol_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Devuelve (df_conf, df_tddft) con todos los campos enriquecidos:
      df_conf  — una fila por confórmero
      df_tddft — una fila por excitación TDDFT de cada confórmero
    """
    mol_dir = WORK_ROOT / mol_name
    if not mol_dir.exists():
        raise FileNotFoundError(f"No se encontró {mol_dir}")

    conf_rows: list[dict] = []
    tddft_rows: list[dict] = []

    sp_dirs = sorted((mol_dir / "04_sp").glob("conf_*")) \
              if (mol_dir / "04_sp").exists() else []

    for sp_dir in sp_dirs:
        conf_id  = sp_dir.name
        sp_out   = sp_dir / "sp.out"
        e_sp     = parse_sp_energy(sp_out)

        # Energía DFT post-opt (de la optimización, no del SP)
        opt_out  = mol_dir / "02_opt" / conf_id / "opt.out"
        e_opt    = parse_opt_energy(opt_out)

        # Frecuencias imaginarias
        freq_out   = mol_dir / "03_numfreq" / conf_id / "freq.out"
        imag_freqs = parse_imaginary_freqs(freq_out)

        conf_rows.append({
            "conformer":      conf_id,
            "E_hartree":      e_sp,       # SP de alta precisión (def2-TZVPD)
            "E_opt_hartree":  e_opt,      # energía al final de la Opt (def2-SVP)
            "has_imaginary":  len(imag_freqs) > 0,
            "n_imaginary":    len(imag_freqs),
            "imag_freqs_cm1": str(imag_freqs) if imag_freqs else "—",
        })

        # TDDFT
        tddft_out = mol_dir / "05_tddft" / conf_id / "tddft.out"
        for exc in parse_tddft_excitations(tddft_out):
            tddft_rows.append({"conformer": conf_id, **exc})

    df_conf = pd.DataFrame(conf_rows)

    if df_conf.empty:
        print(f"⚠  No se encontraron resultados SP para '{mol_name}'")
        return df_conf, pd.DataFrame(tddft_rows)

    # ── Correcciones térmicas ORCA (NumFreq) por confórmero ────────────────
    # compute_boltzmann_multi_t() necesita la columna "thermo_ok" (creada por
    # compute_g_total) para usar el ciclo de Gibbs compuesto; sin este paso
    # caía siempre al modo "electronic_only" (ΔE pura) aunque el NumFreq
    # se hubiera parseado bien.
    _thermochem_cache: list[Optional[dict]] = []
    for _, row in df_conf.iterrows():
        conf_id  = row["conformer"]
        freq_out = mol_dir / "03_numfreq" / conf_id / "freq.out"
        tc = parse_thermochem(freq_out, T=T_STD)
        _thermochem_cache.append(tc if tc.get("ok") else None)
    df_conf["_thermochem"] = _thermochem_cache
    df_conf = compute_g_total(df_conf)

    # Boltzmann multi-temperatura
    df_conf = compute_boltzmann_multi_t(df_conf)

    # Análisis CREST→DFT
    df_conf = analyze_crest_to_dft(mol_dir, df_conf)

    df_tddft = pd.DataFrame(tddft_rows)

    # λ_max Boltzmann-promediado (S1, 298 K) — adjunto como atributo del DF
    if not df_tddft.empty and "boltz_298K" in df_conf.columns:
        s1 = df_tddft[df_tddft["state"] == 1].copy()
        s1 = s1.merge(df_conf[["conformer", "boltz_298K"]], on="conformer", how="left")
        w = s1["boltz_298K"].fillna(0)
        if w.sum() > 0:
            df_conf.attrs["lambda_max_boltz_298K"] = float(
                (s1["wavelength_nm"] * w).sum() / w.sum()
            )
            df_conf.attrs["f_max_boltz_298K"] = float(
                (s1["fosc"] * w).sum() / w.sum()
            )

    return df_conf, df_tddft


# ─────────────────────────────────────────────
#  IMPRESIÓN DE RESULTADOS
# ─────────────────────────────────────────────

def print_summary(df_conf: pd.DataFrame, df_tddft: pd.DataFrame, mol_name: str):
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  RESULTADOS: {mol_name}")
    print(sep)

    if df_conf.empty:
        print("  Sin datos de confórmeros.")
        return

    # ── Tabla de confórmeros base ──────────────────────────────────────────
    print("\n  ── Confórmeros (SP def2-TZVPD / CPCM) ──")
    base_cols = ["conformer", "E_hartree", "dE_kcal", "boltz_298K",
                 "has_imaginary", "n_imaginary"]
    base_cols = [c for c in base_cols if c in df_conf.columns]
    print(df_conf[base_cols].to_string(
        index=False, float_format=lambda x: f"{x:.6f}"
    ))

    lmax = df_conf.attrs.get("lambda_max_boltz_298K")
    if lmax:
        fmax = df_conf.attrs.get("f_max_boltz_298K", float("nan"))
        print(f"\n  λ_max Boltzmann-promediado (298 K, S1): {lmax:.1f} nm  "
              f"(f_avg = {fmax:.5f})")

    # ── Análisis CREST→DFT ────────────────────────────────────────────────
    print_crest_dft_analysis(df_conf)

    # ── Boltzmann multi-temperatura ────────────────────────────────────────
    print_boltzmann_analysis(df_conf)

    # ── Excitaciones TDDFT ────────────────────────────────────────────────
    if not df_tddft.empty:
        print("\n  ── Excitaciones S1 por confórmero ──")
        s1 = df_tddft[df_tddft["state"] == 1].copy()
        if not s1.empty:
            s1 = s1.merge(
                df_conf[["conformer", "boltz_298K", "dE_kcal"]],
                on="conformer", how="left"
            )
            disp_cols = ["conformer", "wavelength_nm", "E_eV", "fosc",
                         "boltz_298K", "dE_kcal"]
            disp_cols = [c for c in disp_cols if c in s1.columns]
            print(s1[disp_cols].to_string(
                index=False, float_format=lambda x: f"{x:.4f}"
            ))

        # Top 5 excitaciones por fosc ponderada
        if "boltz_298K" in df_conf.columns:
            top = df_tddft.merge(
                df_conf[["conformer", "boltz_298K"]], on="conformer", how="left"
            ).copy()
            top["f_weighted"] = top["fosc"] * top["boltz_298K"].fillna(0)
            top5 = top.nlargest(5, "f_weighted")[
                ["conformer", "state", "wavelength_nm", "E_eV", "fosc", "f_weighted"]
            ]
            print("\n  ── Top 5 excitaciones por intensidad ponderada (298 K) ──")
            print(top5.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print(f"\n{sep}\n")


# ─────────────────────────────────────────────
#  COMPARACIÓN ENTRE ESTADOS DE PROTONACIÓN
# ─────────────────────────────────────────────

def compare_protonation_states(
    ref_mol: Optional[str] = None,
) -> pd.DataFrame:
    """
    Recopila resultados de todos los estados de protonación disponibles en runs/
    y construye una tabla comparativa.

    La tabla incluye por mol:
      - Número de confórmeros con SP completado
      - Confórmero de mínima energía
      - Distribución de peso Boltzmann (confs > 1%, > 5%, > 20%)
      - N_eff a 298 K y 250 K
      - λ_max Boltzmann-promediado (S1, 298 K)
      - f_avg (298 K)
      - E_min (Hartree) — para comparar entre estados si mismo nivel de teoría

    Si ref_mol se especifica, calcula también ΔE relativa entre estados
    (solo significativa si tienen la misma carga y nivel de teoría).

    Retorna DataFrame de comparación.
    """
    mols = sorted(p.name for p in WORK_ROOT.iterdir()
                  if p.is_dir() and (p / "04_sp").exists())

    if not mols:
        print("No se encontraron resultados SP en ninguna molécula bajo runs/")
        return pd.DataFrame()

    rows = []
    dfs_conf: dict[str, pd.DataFrame] = {}
    dfs_tddft: dict[str, pd.DataFrame] = {}

    for mol in mols:
        try:
            df_c, df_t = collect_results(mol)
        except Exception as e:
            print(f"  ⚠  {mol}: error al leer resultados — {e}")
            continue

        dfs_conf[mol]  = df_c
        dfs_tddft[mol] = df_t

        if df_c.empty:
            continue

        valid = df_c["E_hartree"].notna()
        n_confs = valid.sum()
        e_min   = df_c.loc[valid, "E_hartree"].min() if n_confs > 0 else None
        best    = df_c.loc[df_c["E_hartree"].idxmin(), "conformer"] \
                  if n_confs > 0 else "—"

        # Distribución conformacional (298 K)
        col298 = "boltz_298K"
        n_1pct  = int((df_c[col298] > 0.01).sum()) if col298 in df_c.columns else None
        n_5pct  = int((df_c[col298] > 0.05).sum()) if col298 in df_c.columns else None
        n_20pct = int((df_c[col298] > 0.20).sum()) if col298 in df_c.columns else None
        n_eff   = df_c.attrs.get("N_eff_298K")

        # N_eff a 250 K
        col250 = "boltz_250K"
        if col250 in df_c.columns:
            p250 = df_c.loc[valid, col250].clip(lower=1e-15)
            s250 = -(p250 * np.log(p250)).sum()
            n_eff_250 = np.exp(s250)
        else:
            n_eff_250 = None

        lmax = df_c.attrs.get("lambda_max_boltz_298K")
        favg = df_c.attrs.get("f_max_boltz_298K")

        # Inversiones CREST→DFT
        n_inv = None
        if "rank_change" in df_c.columns:
            n_inv = int((df_c["rank_change"].abs() > 0).sum())

        rows.append({
            "molecule":       mol,
            "n_confs_sp":     n_confs,
            "best_conformer": best,
            "E_min_hartree":  e_min,
            "N_eff_298K":     n_eff,
            "N_eff_250K":     n_eff_250,
            "n_confs_>1%":    n_1pct,
            "n_confs_>5%":    n_5pct,
            "n_confs_>20%":   n_20pct,
            "lambda_max_nm":  lmax,
            "f_avg_298K":     favg,
            "n_rank_inversions": n_inv,
        })

    df_cmp = pd.DataFrame(rows)

    if ref_mol and ref_mol in df_cmp["molecule"].values:
        e_ref = df_cmp.loc[df_cmp["molecule"] == ref_mol, "E_min_hartree"].values[0]
        if e_ref is not None:
            df_cmp["dE_vs_ref_kcal"] = (df_cmp["E_min_hartree"] - e_ref) * HARTREE_TO_KCAL
            df_cmp["dE_vs_ref_kJ"]   = df_cmp["dE_vs_ref_kcal"] * 4.184

    return df_cmp, dfs_conf, dfs_tddft


def print_comparison(df_cmp: pd.DataFrame, ref_mol: Optional[str] = None):
    sep = "─" * 80
    print(f"\n{sep}")
    print(f"  COMPARACIÓN DE ESTADOS DE PROTONACIÓN  "
          f"({'ref: ' + ref_mol if ref_mol else 'sin referencia'})")
    print(sep)

    if df_cmp.empty:
        print("  Sin datos para comparar.")
        return

    # Tabla principal
    core_cols = ["molecule", "n_confs_sp", "N_eff_298K", "N_eff_250K",
                 "n_confs_>1%", "lambda_max_nm", "f_avg_298K",
                 "n_rank_inversions"]
    if "dE_vs_ref_kcal" in df_cmp.columns:
        core_cols.insert(4, "dE_vs_ref_kcal")
    core_cols = [c for c in core_cols if c in df_cmp.columns]

    print("\n  ── Resumen por estado ──")
    print(df_cmp[core_cols].to_string(
        index=False, float_format=lambda x: f"{x:.3f}"
    ))

    # Destacar quién tiene mayor peso UV (λ_max más alto)
    if "lambda_max_nm" in df_cmp.columns and df_cmp["lambda_max_nm"].notna().any():
        idx_max = df_cmp["lambda_max_nm"].idxmax()
        mol_max = df_cmp.loc[idx_max, "molecule"]
        lmax    = df_cmp.loc[idx_max, "lambda_max_nm"]
        print(f"\n  Estado con λ_max más alto: {mol_max}  ({lmax:.1f} nm)")

    # N_eff: cuál tiene ensemble más desordenado (más conformacionalmente flexible)
    if "N_eff_298K" in df_cmp.columns and df_cmp["N_eff_298K"].notna().any():
        idx_neff = df_cmp["N_eff_298K"].idxmax()
        mol_neff = df_cmp.loc[idx_neff, "molecule"]
        neff_v   = df_cmp.loc[idx_neff, "N_eff_298K"]
        print(f"  Estado más conformacionalmente flexible (298 K): {mol_neff}  "
              f"(N_eff = {neff_v:.2f})")

    if "dE_vs_ref_kcal" in df_cmp.columns:
        print(f"\n  ΔE relativas (mínimo de cada estado vs {ref_mol}):")
        print("  NOTA: ΔE entre estados de protonación con distinta carga no son")
        print("        termodinámicamente directamente comparables sin correcciones")
        print("        de solvatación absoluta y energía del protón.")

    print(f"\n{sep}\n")


# ─────────────────────────────────────────────
#  EXPORTACIÓN
# ─────────────────────────────────────────────

def export_results(df_conf: pd.DataFrame, df_tddft: pd.DataFrame,
                   mol_name: str, fmt: str):
    mol_dir = WORK_ROOT / mol_name

    # Añadir columnas de Boltzmann para todas las temperaturas a df_tddft
    if not df_tddft.empty and not df_conf.empty:
        boltz_cols = [c for c in df_conf.columns if c.startswith("boltz_")]
        df_tddft = df_tddft.merge(
            df_conf[["conformer"] + boltz_cols + ["dE_kcal"]],
            on="conformer", how="left"
        )

    if fmt in ("csv", "both"):
        df_conf.to_csv(mol_dir / "conformers.csv", index=False)
        df_tddft.to_csv(mol_dir / "tddft_excitations.csv", index=False)
        print(f"  CSVs guardados en {mol_dir}/")

    if fmt in ("xlsx", "both"):
        xlsx_path = mol_dir / f"{mol_name}_results.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            df_conf.to_excel(writer, sheet_name="Conformers", index=False)
            df_tddft.to_excel(writer, sheet_name="TDDFT", index=False)

            # Hoja resumen CREST→DFT
            crest_cols = [c for c in df_conf.columns
                          if c in ("conformer", "rank_crest", "rank_dft",
                                   "rank_change", "dE_xtb_kcal",
                                   "dE_kcal", "dE_reranking_kcal")]
            if len(crest_cols) > 2:
                df_conf[crest_cols].to_excel(
                    writer, sheet_name="CREST_vs_DFT", index=False
                )

            # Hoja multi-T Boltzmann
            boltz_cols = ["conformer"] + [c for c in df_conf.columns
                                           if c.startswith("boltz_")]
            if len(boltz_cols) > 1:
                df_conf[boltz_cols].to_excel(
                    writer, sheet_name="Boltzmann_multiT", index=False
                )

        print(f"  Excel guardado en: {xlsx_path}")


def export_adc2_mol(df_adc2: "pd.DataFrame", mol_name: str, fmt: str):
    """Exporta excitaciones ADC2 de una molecula al CSV/XLSX del estado."""
    if df_adc2.empty:
        return
    mol_dir = WORK_ROOT / mol_name
    if fmt in ("csv", "both"):
        out = mol_dir / "adc2_excitations.csv"
        df_adc2.to_csv(out, index=False)
        print(f"  ADC2 CSV -> {out}")
    if fmt in ("xlsx", "both"):
        xlsx_path = mol_dir / f"{mol_name}_results.xlsx"
        mode = "a" if xlsx_path.exists() else "w"
        kw   = {"if_sheet_exists": "replace"} if mode == "a" else {}
        try:
            with pd.ExcelWriter(xlsx_path, engine="openpyxl", mode=mode, **kw) as w:
                df_adc2.to_excel(w, sheet_name="ADC2", index=False)
            print(f"  ADC2 -> {xlsx_path} (hoja: ADC2)")
        except Exception as exc:
            warnings.warn(f"[ADC2 export] {exc}")


def export_comparison(df_cmp: pd.DataFrame, fmt: str):
    if df_cmp.empty:
        return
    out = WORK_ROOT / "comparison_protonation_states"
    if fmt in ("csv", "both"):
        df_cmp.to_csv(out.with_suffix(".csv"), index=False)
        print(f"  CSV comparativo guardado en: {out.with_suffix('.csv')}")
    if fmt in ("xlsx", "both"):
        df_cmp.to_excel(out.with_suffix(".xlsx"), index=False)
        print(f"  Excel comparativo guardado en: {out.with_suffix('.xlsx')}")
