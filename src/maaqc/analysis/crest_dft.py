"""Comparación de ranking conformacional CREST (GFN2-xTB) vs DFT.

Migrado 1:1 desde results_collector.py (legacy):
  - analyze_crest_to_dft()       líneas 524-587
  - print_crest_dft_analysis()   líneas 588-647

Vive en analysis/ (no en stages/censo.py) porque es análisis post-pipeline
sobre resultados ya calculados: lee energías xTB de `01_crest/crest_conformers.xyz`
(vía io/parsing.parse_crest_energies) y las contrasta con las energías DFT
en el DataFrame de confórmeros que produce results_collector.collect_results().
CensoStage opera antes en el DAG (pre-screening del ensemble); no comparte
inputs ni momento del flujo con este módulo.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from ..constants import HARTREE_TO_KCAL
from ..io.parsing import parse_crest_energies


def analyze_crest_to_dft(mol_dir: Path, df_conf: pd.DataFrame) -> pd.DataFrame:
    """
    Compara el orden conformacional de CREST (GFN2-xTB) con el orden DFT post-opt.

    Añade al DataFrame de confórmeros:
      - E_xtb_hartree    : energía xTB del confórmero en CREST (Hartree)
      - dE_xtb_kcal      : ΔE(xTB) respecto al mínimo xTB (kcal/mol)
      - rank_crest        : rango CREST (1 = más estable por xTB)
      - rank_dft          : rango DFT (1 = más estable por SP)
      - rank_change       : rank_dft - rank_crest (0 = no hubo inversión)
      - dE_reranking_kcal : ΔE_DFT - ΔE_xTB (cuánto cambió la brecha energética)

    Retorna df_conf enriquecido.
    """
    conformers_xyz = mol_dir / "01_crest" / "crest_conformers.xyz"
    xtb_energies = parse_crest_energies(conformers_xyz)

    if not xtb_energies:
        warnings.warn("[CREST→DFT] No se pudieron leer energías xTB de CREST.")
        return df_conf

    # Mapear conf_id → índice CREST (conf_001 = índice 0, etc.)
    xtb_map: dict[str, float | None] = {}
    for conf_id in df_conf["conformer"]:
        try:
            idx = int(conf_id.split("_")[1]) - 1   # "conf_003" → 2
        except (IndexError, ValueError):
            continue
        xtb_map[conf_id] = xtb_energies[idx] if idx < len(xtb_energies) else None

    df = df_conf.copy()
    df["E_xtb_hartree"] = df["conformer"].map(xtb_map)

    valid_xtb = df["E_xtb_hartree"].notna()
    valid_dft = df["E_hartree"].notna()

    if valid_xtb.any():
        e_xtb_min = df.loc[valid_xtb, "E_xtb_hartree"].min()
        df["dE_xtb_kcal"] = (df["E_xtb_hartree"] - e_xtb_min) * HARTREE_TO_KCAL

    # Rankings (1 = más estable)
    if valid_xtb.any():
        df.loc[valid_xtb, "rank_crest"] = (
            df.loc[valid_xtb, "E_xtb_hartree"]
            .rank(method="min")
            .astype("Int64")
        )
    if valid_dft.any():
        df.loc[valid_dft, "rank_dft"] = (
            df.loc[valid_dft, "E_hartree"]
            .rank(method="min")
            .astype("Int64")
        )

    if "rank_crest" in df.columns and "rank_dft" in df.columns:
        df["rank_change"] = df["rank_dft"] - df["rank_crest"]

    # Cuánto cambió la brecha energética relativa al pasar de xTB a DFT
    if "dE_xtb_kcal" in df.columns and "dE_kcal" in df.columns:
        df["dE_reranking_kcal"] = df["dE_kcal"] - df["dE_xtb_kcal"]

    return df


def print_crest_dft_analysis(df: pd.DataFrame) -> None:
    """Imprime tabla de análisis de reordenamiento CREST→DFT."""
    needed = ["conformer", "rank_crest", "rank_dft", "rank_change",
              "dE_xtb_kcal", "dE_kcal", "dE_reranking_kcal"]
    cols = [c for c in needed if c in df.columns]
    if len(cols) < 3:
        print("  (energías xTB no disponibles — omitiendo análisis CREST→DFT)")
        return

    sub = df[cols].copy()
    inversions = df["rank_change"].abs() if "rank_change" in df.columns else pd.Series(dtype=float)
    n_inv = (inversions > 0).sum() if not inversions.empty else 0

    print("\n  ── Reordenamiento conformacional CREST (GFN2-xTB) → DFT ──")
    print(f"  Inversiones de rango detectadas: {n_inv}/{len(df)}")

    # Cabecera
    header_map = {
        "conformer":         "conf",
        "rank_crest":        "rank_xTB",
        "rank_dft":          "rank_DFT",
        "rank_change":       "Δrank",
        "dE_xtb_kcal":       "ΔE_xTB(kcal)",
        "dE_kcal":           "ΔE_DFT(kcal)",
        "dE_reranking_kcal": "ΔΔE(kcal)",
    }
    print("\n  " + "  ".join(f"{header_map.get(c, c):>13}" for c in cols))
    print("  " + "─" * (15 * len(cols)))
    for _, row in sub.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            if c == "conformer":
                vals.append(f"{v:>13}")
            elif pd.isna(v):
                vals.append(f"{'—':>13}")
            elif isinstance(v, (int, np.integer)):
                flag = " ◀" if c == "rank_change" and abs(int(v)) > 0 else ""
                vals.append(f"{int(v):>12}{flag}")
            else:
                vals.append(f"{float(v):>13.4f}")
        print("  " + "  ".join(vals))

    if "dE_reranking_kcal" in df.columns:
        max_shift = df["dE_reranking_kcal"].abs().max()
        print(f"\n  Máximo desplazamiento energético xTB→DFT: {max_shift:.3f} kcal/mol")
        if max_shift > 2.0:
            print("  ⚠  La función de energía xTB difiere significativamente de DFT")
            print("     en este sistema. Recomendable ampliar el muestreo CREST.")
