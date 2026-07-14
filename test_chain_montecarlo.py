#!/usr/bin/env python3
"""
test_chain_montecarlo.py — Sanity tests para run_chain_montecarlo() y
_boltzmann_5state_vectorized() (Fase 2)

No requiere archivos ORCA reales: mockea compute_pka() con un side_effect
que responde distinto según qué par (data_acid, data_base) se le pasó,
para simular los tres tramos de la cadena (A→R/D, R/D→B, B→C).

Uso:
    conda activate sciencebitch
    python test_chain_montecarlo.py
"""

from __future__ import annotations
from unittest.mock import patch

import numpy as np
import pandas as pd

import protonation_cycle as pc


# ── Estados sintéticos de la cadena (solo necesitan un 'id' para que el
#    side_effect de compute_pka() los reconozca) ──────────────────────────
DATA_A       = {"id": "A"}
DATA_NEUTRAL = {"id": "neutral"}
DATA_B       = {"id": "B"}
DATA_C       = {"id": "C"}

# pKa's deterministicos sintéticos para cada tramo (valores arbitrarios,
# pero con pKa1 < pKa2 como correspondería físicamente a un diácido)
PKA0_CENTRAL, PKA1_CENTRAL, PKA2_CENTRAL = 3.0, 7.5, 12.0


def _dG_aq_from_pka(pka: float, T: float = pc.T_STD) -> float:
    return pka * (pc.R_KCAL * T * pc.LN10)


def make_chain_side_effect(with_terms: bool = True):
    """Fabrica un side_effect para mockear compute_pka() según el par de
    estados recibido, devolviendo un base_result determinístico distinto
    por tramo (pKa0: A→neutral, pKa1: neutral→B, pKa2: B→C)."""

    def _side_effect(data_acid, data_base, T=pc.T_STD, method="auto",
                      gcorr_source="harmonic", verbose=False):
        pair = (data_acid.get("id"), data_base.get("id"))
        pka_map = {
            ("A", "neutral"): PKA0_CENTRAL,
            ("neutral", "B"): PKA1_CENTRAL,
            ("B", "C"):       PKA2_CENTRAL,
        }
        if pair not in pka_map:
            return {"ok": False, "method": None, "warnings": [f"par no reconocido: {pair}"]}

        pka = pka_map[pair]
        result = {
            "ok": True,
            "method": "direct",
            "T": T,
            "dG_aq_kcal": _dG_aq_from_pka(pka, T),
            "pKa": pka,
            "warnings": [],
        }
        if with_terms:
            result["E_sp_aq_acid_Eh"] = -500.0
            result["E_sp_aq_base_Eh"] = -499.5
        else:
            result["E_sp_aq_acid_Eh"] = None
            result["E_sp_aq_base_Eh"] = None
        return result

    return _side_effect


def test_end_to_end_chain_runs_ok():
    pH_values = list(np.arange(0, 14.1, 1.0))
    with patch.object(pc, "compute_pka", side_effect=make_chain_side_effect()):
        chain = pc.run_chain_montecarlo(
            DATA_A, DATA_NEUTRAL, DATA_B, DATA_C,
            f_R=0.9, f_D=0.1, pH_values=pH_values,
            n_samples=5_000, seed=42,
        )

    assert chain["ok"], chain["warnings"]
    for k in ("mc_pka0", "mc_pka1", "mc_pka2"):
        assert chain[k] is not None and chain[k]["ok"], k

    df = chain["df_dist_bands"]
    assert not df.empty
    expected_cols = {
        "pH",
        "f_neutral_macro_lo", "f_neutral_macro_median", "f_neutral_macro_hi",
        "f_R_lo", "f_R_median", "f_R_hi",
        "f_D_lo", "f_D_median", "f_D_hi",
        "f_cation_lo", "f_cation_median", "f_cation_hi",
        "f_anion_lo", "f_anion_median", "f_anion_hi",
        "f_dianion_lo", "f_dianion_median", "f_dianion_hi",
    }
    assert expected_cols.issubset(set(df.columns)), set(df.columns)
    print(f"OK  test_end_to_end_chain_runs_ok: df_dist_bands shape={df.shape}, "
          f"columnas OK")


def test_vectorized_matches_scalar_at_zero_noise():
    """Con sigma=0 en todo, cada muestra MC es idéntica al valor
    determinístico -> la mediana de _boltzmann_5state_vectorized() debe
    coincidir con boltzmann_5state() evaluado en los pKa deterministicos."""
    pH_values = list(np.arange(0, 14.1, 0.5))
    f_R, f_D = 0.8, 0.2

    with patch.object(pc, "compute_pka", side_effect=make_chain_side_effect()):
        chain = pc.run_chain_montecarlo(
            DATA_A, DATA_NEUTRAL, DATA_B, DATA_C,
            f_R=f_R, f_D=f_D, pH_values=pH_values,
            n_samples=2_000, seed=1,
            sigma_elec_kcal=0.0, sigma_gcorr_kcal=0.0, sigma_hplus_kcal=0.0,
        )
    assert chain["ok"]
    df_mc = chain["df_dist_bands"]

    df_ref = pc.boltzmann_5state(
        PKA0_CENTRAL, PKA1_CENTRAL, PKA2_CENTRAL, f_R, f_D, pH_values,
    )

    for frac_ref, frac_mc in [
        ("f_neutral_macro", "f_neutral_macro_median"),
        ("f_R", "f_R_median"),
        ("f_D", "f_D_median"),
        ("f_cation", "f_cation_median"),
        ("f_anion", "f_anion_median"),
        ("f_dianion", "f_dianion_median"),
    ]:
        np.testing.assert_allclose(
            df_mc[frac_mc].to_numpy(), df_ref[frac_ref].to_numpy(),
            atol=1e-9, err_msg=f"mismatch en {frac_ref}",
        )
        # con sigma=0 las bandas deben colapsar (lo == hi == median)
        lo_col = frac_mc.replace("_median", "_lo")
        hi_col = frac_mc.replace("_median", "_hi")
        np.testing.assert_allclose(df_mc[lo_col].to_numpy(), df_mc[frac_mc].to_numpy(), atol=1e-9)
        np.testing.assert_allclose(df_mc[hi_col].to_numpy(), df_mc[frac_mc].to_numpy(), atol=1e-9)

    print("OK  test_vectorized_matches_scalar_at_zero_noise: "
          "mediana MC == boltzmann_5state() determinístico, bandas colapsadas")


def test_shared_hplus_correlates_pka1_pka2():
    """El shared_hplus_draw debe producir una correlación positiva fuerte
    entre pKa1_samples y pKa2_samples (mismo desplazamiento de la
    constante del protón afecta a ambos con el mismo signo), en
    contraste con generarlos de forma independiente."""
    with patch.object(pc, "compute_pka", side_effect=make_chain_side_effect()):
        chain = pc.run_chain_montecarlo(
            DATA_A, DATA_NEUTRAL, DATA_B, DATA_C,
            f_R=1.0, f_D=0.0, pH_values=[7.0],
            n_samples=20_000, seed=7,
            sigma_elec_kcal=0.3, sigma_gcorr_kcal=0.2, sigma_hplus_kcal=2.5,
        )
    assert chain["ok"]
    pka1 = chain["mc_pka1"]["pka_samples"]
    pka2 = chain["mc_pka2"]["pka_samples"]
    corr_shared = np.corrcoef(pka1, pka2)[0, 1]

    # Referencia: mismos sigmas pero SIN compartir el draw de H+ (cada
    # tramo con su propio rng, sin shared_hplus_draw)
    with patch.object(pc, "compute_pka", side_effect=make_chain_side_effect()):
        mc1_indep = pc.compute_pka_montecarlo(
            DATA_NEUTRAL, DATA_B, n_samples=20_000,
            sigma_elec_kcal=0.3, sigma_gcorr_kcal=0.2, sigma_hplus_kcal=2.5,
            rng=np.random.default_rng(101),
        )
        mc2_indep = pc.compute_pka_montecarlo(
            DATA_B, DATA_C, n_samples=20_000,
            sigma_elec_kcal=0.3, sigma_gcorr_kcal=0.2, sigma_hplus_kcal=2.5,
            rng=np.random.default_rng(202),
        )
    corr_indep = np.corrcoef(mc1_indep["pka_samples"], mc2_indep["pka_samples"])[0, 1]

    assert corr_shared > 0.5, corr_shared
    assert abs(corr_indep) < 0.05, corr_indep
    print(f"OK  test_shared_hplus_correlates_pka1_pka2: "
          f"corr(compartido)={corr_shared:.3f}  corr(independiente)={corr_indep:.3f}")


def test_missing_cation_state_omits_branch():
    """Si data_A es None (no hay catión, como en el caso real cuando
    --protonated se omite), pKa0 debe omitirse y f_cation no debe aparecer
    en las bandas — igual que hace boltzmann_5state()."""
    with patch.object(pc, "compute_pka", side_effect=make_chain_side_effect()):
        chain = pc.run_chain_montecarlo(
            None, DATA_NEUTRAL, DATA_B, DATA_C,
            f_R=1.0, f_D=0.0, pH_values=[1.0, 7.0, 13.0],
            n_samples=2_000, seed=5,
        )
    assert chain["ok"]
    assert chain["mc_pka0"] is None
    assert chain["mc_pka1"] is not None and chain["mc_pka1"]["ok"]
    df = chain["df_dist_bands"]
    assert "f_cation_median" not in df.columns
    assert "f_anion_median" in df.columns
    assert "f_dianion_median" in df.columns
    print("OK  test_missing_cation_state_omits_branch: f_cation ausente, resto presente")


def test_bands_are_ordered_lo_le_median_le_hi():
    pH_values = list(np.arange(0, 14.1, 1.0))
    with patch.object(pc, "compute_pka", side_effect=make_chain_side_effect()):
        chain = pc.run_chain_montecarlo(
            DATA_A, DATA_NEUTRAL, DATA_B, DATA_C,
            f_R=0.9, f_D=0.1, pH_values=pH_values,
            n_samples=5_000, seed=13,
        )
    df = chain["df_dist_bands"]
    for frac in ["f_neutral_macro", "f_R", "f_D", "f_cation", "f_anion", "f_dianion"]:
        lo, med, hi = df[f"{frac}_lo"], df[f"{frac}_median"], df[f"{frac}_hi"]
        assert (lo <= med + 1e-12).all(), frac
        assert (med <= hi + 1e-12).all(), frac
        assert (lo >= -1e-9).all() and (hi <= 1 + 1e-9).all(), frac
    print("OK  test_bands_are_ordered_lo_le_median_le_hi")


def test_macrostate_states_without_separate_terms():
    """El neutro/zwitterión suele venir como macroestado combinado (sin
    E_sp_aq por separado, ver data_neutral_macro en _run_extended); la
    cadena completa no debe romperse por eso, solo debe caer al fallback de
    compute_pka_montecarlo (sigma combinado) para los tramos que lo tocan."""
    with patch.object(pc, "compute_pka", side_effect=make_chain_side_effect(with_terms=False)):
        chain = pc.run_chain_montecarlo(
            DATA_A, DATA_NEUTRAL, DATA_B, DATA_C,
            f_R=1.0, f_D=0.0, pH_values=[7.0],
            n_samples=2_000, seed=9,
        )
    assert chain["ok"]
    assert any("Sin desglose E_sp_aq/G_corr" in w for w in chain["warnings"])
    print("OK  test_macrostate_states_without_separate_terms: fallback propagado sin romper la cadena")


def test_reproducibility_with_seed():
    with patch.object(pc, "compute_pka", side_effect=make_chain_side_effect()):
        chain1 = pc.run_chain_montecarlo(
            DATA_A, DATA_NEUTRAL, DATA_B, DATA_C,
            f_R=1.0, f_D=0.0, pH_values=[7.0],
            n_samples=1_000, seed=321,
        )
        chain2 = pc.run_chain_montecarlo(
            DATA_A, DATA_NEUTRAL, DATA_B, DATA_C,
            f_R=1.0, f_D=0.0, pH_values=[7.0],
            n_samples=1_000, seed=321,
        )
    np.testing.assert_allclose(
        chain1["mc_pka0"]["pka_samples"], chain2["mc_pka0"]["pka_samples"],
    )
    np.testing.assert_allclose(
        chain1["df_dist_bands"]["f_anion_median"].to_numpy(),
        chain2["df_dist_bands"]["f_anion_median"].to_numpy(),
    )
    print("OK  test_reproducibility_with_seed: misma semilla -> mismos resultados")


if __name__ == "__main__":
    tests = [
        test_end_to_end_chain_runs_ok,
        test_vectorized_matches_scalar_at_zero_noise,
        test_shared_hplus_correlates_pka1_pka2,
        test_missing_cation_state_omits_branch,
        test_bands_are_ordered_lo_le_median_le_hi,
        test_macrostate_states_without_separate_terms,
        test_reproducibility_with_seed,
    ]
    for t in tests:
        t()
    print(f"\nTodos los tests de Fase 2 pasaron ({len(tests)}/{len(tests)}).")
