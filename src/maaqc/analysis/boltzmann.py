"""Ponderación Boltzmann multi-temperatura y G_total compuesto.

Migrado 1:1 desde results_collector.py:
  - compute_G_total()               líneas 648-764  → compute_g_total()
  - compute_boltzmann_multiT()      líneas 765-848  → compute_boltzmann_multi_t()
  - print_boltzmann_analysis()      líneas 849-932  → print_boltzmann_analysis()

Revisado el posible solapamiento con protonation_cycle.py:1586-1703
(`boltzmann_protonation_states`, `boltzmann_5state`): NO hay duplicación.
Esas dos funciones calculan la distribución de MACROESTADOS de
protonación (A⁺/R/B⁻/C²⁻) en función del pH vía Henderson-Hasselbalch
encadenado, usando pKa's ya calculados — pertenecen a `pka_cycle.py`
(ya migrado según el contexto de la sesión). Las funciones de este
módulo, en cambio, ponderan CONFÓRMEROS dentro de un único estado de
protonación, usando G_total(T) por confórmero. No comparten inputs
(pKa vs. G_total por confórmero) ni granularidad (macroestado vs.
confórmero), así que no hay nada que fusionar ni deduplicar.

NO incluido acá (fuera de alcance de este módulo, viven en pipeline.py /
stages/censo.py):
  - parse_boltzmann_weights()       líneas 2108-2137
  - select_nconf_boltzmann()        líneas 2138-2276
Estas seleccionan/truncan confórmeros para el CENSO, no calculan pesos
— ver el issue abierto de CensoStage debiendo escribir el ensemble
pre-filtrado en el path que espera OptStage._source_xyz().

⚠️ TODO(kmilozeta): las constantes TEMPS_K / TEMP_LABELS / T_STD /
HARTREE_TO_KCAL / KB_KCAL se importan acá desde un módulo compartido
asumido (`..constants`). Si `pka_cycle.py` ya las definió/importó desde
otro lugar, ajustar este import antes de usar — a propósito no dupliqué
las constantes localmente, para que este archivo y pka_cycle.py no
terminen con dos copias que se puedan desincronizar.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from ..constants import (
    HARTREE_TO_KCAL,
    KB_KCAL,
    TEMP_LABELS,
    TEMPS_K,
    T_STD,
)


def compute_g_total(df_conf: pd.DataFrame) -> pd.DataFrame:
    """
    Aplica el ciclo termodinámico compuesto para calcular G_total de cada
    confórmero, usando las temperaturas definidas en TEMPS_K.

    Ciclo implementado (estándar en literatura de MAAs, Ho & Coote 2010):

        G_total_i(T) = E_SP(TZVPD/CPCM) + G_corr_i(T, SVP/CPCM)

    donde G_corr_i(T) = H_corr_i - T·S_total_i
    y H_corr, S_total provienen del NumFreq a def2-SVP/CPCM.

    IMPORTANTE — temperatura del NumFreq:
    El NumFreq de ORCA se corre a 298.15 K por defecto. Para las otras
    temperaturas (250 K, 274.15 K, 310 K) se recalcula G_corr(T) usando
    la aproximación rígido-rotor/oscilador-armónico (RRHO):

        G_corr(T) = H_corr(T_ref) - T · S_total(T_ref)

    Esto es una aproximación válida cuando |T - T_ref| < 60 K (Grimme 2012),
    adecuada para nuestro rango 250–310 K. Para rigor máximo se debería
    relanzar el NumFreq a cada T, pero el costo computacional no lo justifica
    en este sistema.

    Columnas añadidas al DataFrame:
        G_corr_298K       : corrección de Gibbs a 298.15 K (Hartree)
        G_total_{T}K      : energía libre de Gibbs compuesta a cada T (Hartree)
        dG_{T}K_kcal      : ΔG relativo al mínimo a cada T (kcal/mol)
        H_corr            : corrección entálpica (Hartree, independiente de T)
        S_total_EhK       : entropía total (Eh/K)
        ZPE               : energía de punto cero (Hartree)
        thermo_ok         : True si el NumFreq fue parseado correctamente
        thermo_warning    : descripción de problemas encontrados
    """
    df = df_conf.copy()

    # Columnas nuevas con valores por defecto
    df["G_corr_298K"] = np.nan
    df["H_corr"] = np.nan
    df["S_total_EhK"] = np.nan
    df["ZPE"] = np.nan
    df["thermo_ok"] = False
    df["thermo_warning"] = ""

    for label in TEMP_LABELS:
        df[f"G_total_{label}"] = np.nan
        df[f"dG_{label}_kcal"] = np.nan

    for idx, row in df.iterrows():
        tc = row.get("_thermochem", None)
        if tc is None or not isinstance(tc, dict):
            df.at[idx, "thermo_warning"] = "sin datos termoquímicos"
            continue

        if not tc.get("ok", False):
            issues = []
            for field in ("ZPE", "H_corr", "S_total", "G_corr"):
                if tc.get(field) is None:
                    issues.append(field)
            df.at[idx, "thermo_warning"] = f"campos faltantes: {', '.join(issues)}"
            # Intentar uso parcial si al menos H_corr y S_total están
            if tc.get("H_corr") is None or tc.get("S_total") is None:
                continue

        e_sp = row.get("E_hartree")
        if e_sp is None or pd.isna(e_sp):
            df.at[idx, "thermo_warning"] = (
                df.at[idx, "thermo_warning"] + " | E_SP faltante"
            ).lstrip(" | ")
            continue

        H_corr = tc["H_corr"]
        S_total = tc["S_total"]  # Eh/K
        ZPE = tc.get("ZPE")

        df.at[idx, "H_corr"] = H_corr
        df.at[idx, "S_total_EhK"] = S_total
        df.at[idx, "ZPE"] = ZPE
        df.at[idx, "thermo_ok"] = True

        # G_corr a cada T  (aproximación RRHO, válida ±60 K respecto a T_ref)
        for T, label in zip(TEMPS_K, TEMP_LABELS):
            G_corr_T = H_corr - T * S_total
            G_total = e_sp + G_corr_T
            df.at[idx, f"G_total_{label}"] = G_total

        # G_corr de referencia a 298.15 K
        df.at[idx, "G_corr_298K"] = H_corr - T_STD * S_total

    # ΔG relativo al mínimo para cada temperatura
    for T, label in zip(TEMPS_K, TEMP_LABELS):
        col_g = f"G_total_{label}"
        col_dg = f"dG_{label}_kcal"
        valid_g = df["thermo_ok"] & df[col_g].notna()
        if valid_g.any():
            g_min = df.loc[valid_g, col_g].min()
            df.loc[valid_g, col_dg] = (df.loc[valid_g, col_g] - g_min) * HARTREE_TO_KCAL

    # ΔE electrónica pura (de SP) — mantenida para comparación
    valid_e = df["E_hartree"].notna()
    if valid_e.any():
        e_min = df.loc[valid_e, "E_hartree"].min()
        df["dE_kcal"] = (df["E_hartree"] - e_min) * HARTREE_TO_KCAL

    # Advertencia global si hay confórmeros sin termoquímica
    n_missing = (~df["thermo_ok"]).sum()
    if n_missing > 0:
        warnings.warn(
            f"[THERMOCHEM] {n_missing}/{len(df)} confórmeros sin correcciones "
            f"térmicas válidas. El Boltzmann de esos confs se omitirá. "
            f"Detalles en columna 'thermo_warning'."
        )

    return df


def compute_boltzmann_multi_t(df_conf: pd.DataFrame) -> pd.DataFrame:
    """
    Calcula pesos de Boltzmann a todas las temperaturas en TEMPS_K.

    Usa ΔG = G_total(T) - G_min(T) cuando las correcciones térmicas están
    disponibles (ciclo compuesto completo). Si thermo_ok es False para algún
    confórmero, ese confórmero queda excluido del ensemble y su peso es NaN.

    Fallback: si NINGÚN confórmero tiene termoquímica, usa ΔE electrónica
    pura con advertencia explícita en df.attrs["boltzmann_mode"].

    Métricas adicionales calculadas para cada T:
        N_eff(T)      = exp(-Σ p_i · ln(p_i))  (número efectivo de confs)
        S_mix(T)/R    = -Σ p_i · ln(p_i)       (entropía de mezcla en unidades R)

    Columna de conveniencia:
        boltz_298K    = alias de boltz_298K (label "298K"), para compatibilidad
                        con el resto del código.
    """
    df = df_conf.copy()

    # Determinar si el ciclo completo está disponible
    has_thermo = "thermo_ok" in df.columns and df["thermo_ok"].any()

    if not has_thermo:
        warnings.warn(
            "[BOLTZMANN] No hay correcciones térmicas disponibles. "
            "Usando solo energía electrónica SP (ΔE). "
            "Los pesos Boltzmann NO incluyen ZPE ni entropía conformacional. "
            "Reportar como 'Boltzmann weighting based on electronic energies'."
        )
        df.attrs["boltzmann_mode"] = "electronic_only"
        valid = df["E_hartree"].notna()
        e_min = df.loc[valid, "E_hartree"].min()
        df.loc[valid, "dE_kcal"] = (df.loc[valid, "E_hartree"] - e_min) * HARTREE_TO_KCAL

        for T, label in zip(TEMPS_K, TEMP_LABELS):
            col = f"boltz_{label}"
            raw = np.exp(-df.loc[valid, "dE_kcal"] / (KB_KCAL * T))
            df.loc[valid, col] = raw / raw.sum()
    else:
        df.attrs["boltzmann_mode"] = "gibbs_composite"
        for T, label in zip(TEMPS_K, TEMP_LABELS):
            col_g = f"G_total_{label}"
            col_dg = f"dG_{label}_kcal"
            col_b = f"boltz_{label}"

            valid = df["thermo_ok"] & df[col_g].notna()
            if not valid.any():
                warnings.warn(f"[BOLTZMANN] Sin G_total válido a {label}.")
                continue

            dG = df.loc[valid, col_dg]
            raw = np.exp(-dG / (KB_KCAL * T))
            df.loc[valid, col_b] = raw / raw.sum()

    # ── Métricas del ensemble ─────────────────────────────────────────────
    for T, label in zip(TEMPS_K, TEMP_LABELS):
        col_b = f"boltz_{label}"
        if col_b not in df.columns:
            continue
        valid_b = df[col_b].notna()
        if not valid_b.any():
            continue
        p = df.loc[valid_b, col_b].clip(lower=1e-15)
        s_mix = -(p * np.log(p)).sum()  # S_mix / R (adimensional)
        n_eff = np.exp(s_mix)
        df.attrs[f"S_mix_{label}"] = float(s_mix)
        df.attrs[f"N_eff_{label}"] = float(n_eff)
        df.attrs[f"n_above_1pct_{label}"] = int((df[col_b] > 0.01).sum())

    # ── Alias de conveniencia para 298K ──────────────────────────────────
    # La etiqueta "298K" corresponde al índice 2 en TEMP_LABELS
    lbl_298 = "298K"
    col_298 = f"boltz_{lbl_298}"
    if col_298 in df.columns:
        df["boltz_298K"] = df[col_298]
        df.attrs["S_mix_298K"] = df.attrs.get(f"S_mix_{lbl_298}")
        df.attrs["N_eff_298K"] = df.attrs.get(f"N_eff_{lbl_298}")
        df.attrs["n_above_1pct_298K"] = df.attrs.get(f"n_above_1pct_{lbl_298}")

    return df


def print_boltzmann_analysis(df: pd.DataFrame):
    """Imprime tabla de pesos Boltzmann multi-temperatura con modo termodinámico."""
    mode = df.attrs.get("boltzmann_mode", "unknown")
    mode_str = {
        "gibbs_composite": "ΔG compuesto [E_SP(TZVPD) + G_corr(SVP/CPCM)]",
        "electronic_only": "ΔE electrónica pura — SIN ZPE/entropía ⚠",
        "unknown": "modo desconocido",
    }.get(mode, mode)

    print(f"\n  ── Pesos de Boltzmann por temperatura ──")
    print(f"  Modo: {mode_str}")

    avail_labels = [lbl for lbl in TEMP_LABELS if f"boltz_{lbl}" in df.columns]
    if not avail_labels:
        print("  (sin columnas boltz_* disponibles)")
        return

    # Encabezado con temperatura en K
    header = f"\n  {'conf':>10}  "
    header += "  ".join(f"{lbl:>8}" for lbl in avail_labels)
    if mode == "gibbs_composite":
        header += "   ΔG_298(kcal)"
    else:
        header += "   ΔE(kcal)"
    print(header)
    print("  " + "─" * (14 + 10 * len(avail_labels) + 16))

    for _, row in df.iterrows():
        thermo_flag = ""
        if "thermo_ok" in df.columns and not row.get("thermo_ok", True):
            thermo_flag = " [no G_corr]"

        vals = [f"  {row['conformer']:>10}"]
        for lbl in avail_labels:
            col = f"boltz_{lbl}"
            v = row.get(col, np.nan)
            if pd.isna(v):
                vals.append(f"  {'—':>8}")
            else:
                marker = " ★" if v > 0.50 else (" ·" if v > 0.10 else "  ")
                vals.append(f"  {v:6.4f}{marker}")

        # ΔG o ΔE según modo
        dg_col = "dG_298K_kcal" if mode == "gibbs_composite" else "dE_kcal"
        dg_val = row.get(dg_col, np.nan)
        if pd.isna(dg_val):
            vals.append(f"   {'—':>12}")
        else:
            vals.append(f"   {dg_val:+12.4f}")

        print("".join(vals) + thermo_flag)

    # Métricas del ensemble
    print()
    for lbl in avail_labels:
        s = df.attrs.get(f"S_mix_{lbl}")
        n = df.attrs.get(f"N_eff_{lbl}")
        n1 = df.attrs.get(f"n_above_1pct_{lbl}")
        if s is not None:
            print(f"  {lbl:>5}:  S_mix/R = {s:.4f}  │  N_eff = {n:.2f}  │  "
                  f"confs >1% = {n1}/{len(df)}")

    # Sensibilidad a temperatura: ¿cambia el conf dominante entre extremos?
    b_low = f"boltz_{TEMP_LABELS[0]}"  # 250K
    b_high = f"boltz_{TEMP_LABELS[-1]}"  # 310K
    if b_low in df.columns and b_high in df.columns:
        valid_low = df[b_low].notna()
        valid_high = df[b_high].notna()
        if valid_low.any() and valid_high.any():
            dom_low = df.loc[df[b_low].idxmax(), "conformer"]
            dom_high = df.loc[df[b_high].idxmax(), "conformer"]
            if dom_low != dom_high:
                print(f"\n  ⚠  Confórmero dominante cambia con T:")
                print(f"     {TEMP_LABELS[0]} → {dom_low}   │   "
                      f"{TEMP_LABELS[-1]} → {dom_high}")
            else:
                print(f"\n  Confórmero dominante estable entre "
                      f"{TEMP_LABELS[0]}–{TEMP_LABELS[-1]}: {dom_low}")
