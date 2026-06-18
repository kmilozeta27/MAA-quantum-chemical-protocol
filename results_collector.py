#!/usr/bin/env python3
"""
results_collector.py — Recolector y analizador de resultados del pipeline ORCA
                        con ciclo termodinámico completo y pKa teórico

═══════════════════════════════════════════════════════════════════════════════
FUNDAMENTO TERMODINÁMICO
═══════════════════════════════════════════════════════════════════════════════

Ciclo de Gibbs compuesto (Isobaric-isoentropic approach):

  G_i = E_SP(def2-TZVPD) + G_corr(def2-SVP, CPCM)

donde G_corr = H_corr + ZPE - T·S_total (extraída del NumFreq).
Este esquema —SP de alta precisión con corrección térmica de geometría
optimizada— es estándar en la literatura de MAAs y química computacional
de solutos en agua.

Poblaciones de Boltzmann:

  p_i(T) = exp(-ΔG_i / RT) / Σ_j exp(-ΔG_j / RT)

calculadas a T ∈ {250.0, 274.15, 298.15, 310.0} K.
  250.0 K  — mínimo oceánico Antártico (WAP surface, winter)
  274.15 K — temperatura media oceánico Antártico superficial (objetivo FAE)
  298.15 K — condición estándar IUPAC
  310.0 K  — temperatura fisiológica humana

Ciclo termodinámico de pKa (método de Born / ciclo termoquímico):

  AH(aq) → A⁻(aq) + H⁺(aq)       ΔG_aq

  ΔG_aq = ΔG_gas + ΔΔG_solv

donde:
  ΔG_gas   = G(A⁻, gas) + G(H⁺, gas) - G(AH, gas)
             con G(H⁺, gas) = -6.28 kcal/mol  (Tissandier et al. 1998)
  ΔΔG_solv = ΔG_solv(A⁻) + ΔG_solv(H⁺) - ΔG_solv(AH)
             con ΔG_solv(H⁺) = -265.9 kcal/mol (Kelly et al. 2006)
             y ΔG_solv extraída de CPCM como diferencia E(CPCM) - E(gas)

  pKa = ΔG_aq / (RT · ln 10)     a T = 298.15 K

NOTA: Para el pKa se requiere un par de cálculos SP adicionales en fase gas
(sin CPCM) para AH y A⁻, cuyas rutas se configuran en PKA_GAS_SP_DIR.
Si esos archivos no existen, el módulo de pKa se omite con advertencia.

═══════════════════════════════════════════════════════════════════════════════
MÓDULOS DE ANÁLISIS
═══════════════════════════════════════════════════════════════════════════════

  1. Correcciones térmicas (ZPE, H, S, G_corr) desde NumFreq ORCA 6
  2. G_total = E_SP(TZVPD) + G_corr(SVP) — ciclo compuesto rigoroso
  3. Boltzmann multi-temperatura basado en ΔG (no ΔE)
  4. Reordenamiento conformacional CREST (GFN2-xTB) → DFT post-opt
  5. Excitaciones TD-DFT: λ (nm), f, E (eV) + λ_max Boltzmann-promediado
  6. pKa teórico via ciclo termodinámico Born (requiere SP en fase gas)
  7. Comparación entre estados de protonación (todas las mols en runs/)

═══════════════════════════════════════════════════════════════════════════════
USO
═══════════════════════════════════════════════════════════════════════════════

  python results_collector.py --mol mglyval               # un estado
  python results_collector.py --mol mglyval --export xlsx
  python results_collector.py --mol mglyval --plot
  python results_collector.py --compare                   # todos los estados
  python results_collector.py --compare --plot
  python results_collector.py --compare --ref mglyval_neutral
  python results_collector.py --pka --acid mglyval_neutral --base mglyval_anion

Para el pKa, cada mol debe tener un SP en fase gas bajo:
  runs/<mol>/06_sp_gas/conf_<best>/sp_gas.out
(calcular con el mismo funcional/base que el SP acuoso pero sin CPCM)
"""

import argparse
import re
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

try:
    from ase.io import read as ase_read
    HAS_ASE = True
except ImportError:
    HAS_ASE = False

# ─────────────────────────────────────────────
#  CONSTANTES
# ─────────────────────────────────────────────

WORK_ROOT       = Path("runs")
HARTREE_TO_KCAL = 627.5094740631
HARTREE_TO_KJ   = HARTREE_TO_KCAL * 4.184
HARTREE_TO_EV   = 27.211396132
KB_KCAL         = 0.0019872041     # kcal/(mol·K)
R_KCAL          = 1.9872041e-3     # kcal/(mol·K)  — constante de gas ideal
LN10            = 2.302585093

# Temperaturas de análisis:
#   250.0 K  — mínimo oceánico Antártico WAP (invierno)
#   274.15 K — temperatura media superficie oceáno Antártico (objetivo FAE)
#   298.15 K — condición estándar IUPAC
#   310.0 K  — temperatura fisiológica humana
TEMPS_K     = [250.0, 274.15, 298.15, 310.0]
TEMP_LABELS = ["250K", "274K", "298K", "310K"]
T_STD       = 298.15               # K — temperatura de referencia para pKa

# ── Constantes del ciclo de pKa (todas a 298.15 K) ──────────────────────────
# G(H⁺, gas) = -6.28 kcal/mol  [Tissandier et al., J. Phys. Chem. A 1998]
# ΔG_solv(H⁺) = -265.9 kcal/mol [Kelly, Cramer & Truhlar, J. Phys. Chem. B 2006]
# Corrección de estado estándar de fase gas → solución (1 atm → 1 M):
#   δG_conc = RT·ln(RT/P°V°) = RT·ln(24.5) ≈ +1.894 kcal/mol a 298.15 K
G_HPLUS_GAS      = -6.28           # kcal/mol  Tissandier 1998
DG_SOLV_HPLUS    = -265.9          # kcal/mol  Kelly 2006 (absoluto, escala absoluta)
DG_CONC_CORR     = 1.894           # kcal/mol  corrección 1 atm → 1 M (RT·ln(24.5))

# Directorio opcional para SPs en fase gas (para el ciclo de pKa)
# Estructura esperada: runs/<mol>/06_sp_gas/<conf_id>/sp_gas.out
PKA_GAS_SP_SUBDIR = "06_sp_gas"


# ─────────────────────────────────────────────
#  PARSERS — individuales y robustos
# ─────────────────────────────────────────────

def parse_sp_energy(out_file: Path) -> Optional[float]:
    """Extrae la última energía FINAL SINGLE POINT ENERGY del .out de SP."""
    if not out_file.exists():
        return None
    energy = None
    for line in out_file.read_text(errors="replace").splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                energy = float(line.split()[-1])
            except ValueError:
                pass
    return energy


def parse_opt_energy(out_file: Path) -> Optional[float]:
    """
    Extrae la energía del último ciclo de optimización (paso final convergido).
    ORCA reporta 'FINAL SINGLE POINT ENERGY' también durante Opt, así que
    tomamos la última ocurrencia, que corresponde al mínimo convergido.
    """
    return parse_sp_energy(out_file)


def parse_imaginary_freqs(freq_out: Path) -> list[float]:
    """Devuelve lista de frecuencias imaginarias (cm⁻¹, valores negativos)."""
    imag = []
    if not freq_out.exists():
        return imag
    for line in freq_out.read_text(errors="replace").splitlines():
        m = re.match(r"\s+\d+:\s+(-\d+\.\d+)\s+cm", line)
        if m:
            imag.append(float(m.group(1)))
    return imag


def parse_tddft_excitations(tddft_out: Path) -> list[dict]:
    """
    Parsea excitaciones del bloque 'ABSORPTION SPECTRUM VIA TRANSITION
    ELECTRIC DIPOLE MOMENTS' del output TDDFT de ORCA 6.
    Formato: state_label -> final_label   E_eV  E_cm-1  wl_nm  fosc ...
    """
    excitations = []
    if not tddft_out.exists():
        return excitations

    pattern = re.compile(
        r"^\s*\S+\s*->\s*(\d+)-\S+\s+"
        r"([\d.]+)\s+"        # E_eV
        r"[\d.]+\s+"          # E_cm-1 (ignorado)
        r"([\d.]+)\s+"        # wavelength_nm
        r"([\d.]+)"           # fosc
    )
    in_block = False
    for line in tddft_out.read_text(errors="replace").splitlines():
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
            in_block = True
            continue
        if in_block:
            if "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" in line:
                break
            m = pattern.match(line)
            if m:
                excitations.append({
                    "state":         int(m.group(1)),
                    "E_eV":          float(m.group(2)),
                    "wavelength_nm": float(m.group(3)),
                    "fosc":          float(m.group(4)),
                })
    return excitations


def parse_crest_energies(conformers_xyz: Path) -> list[Optional[float]]:
    """
    Lee las energías GFN2-xTB embebidas por CREST en el multi-XYZ.
    CREST escribe la energía en Hartree en la línea de comentario (línea 2
    de cada frame), típicamente como un número flotante negativo solo o
    precedido de 'Etot='.
    Devuelve lista con la energía de cada frame (None si no parseable).
    """
    if not conformers_xyz.exists():
        return []

    energies: list[Optional[float]] = []

    # Intentar primero con ASE (más robusto para multi-XYZ)
    if HAS_ASE:
        try:
            frames = ase_read(str(conformers_xyz), index=":")
            if not isinstance(frames, list):
                frames = [frames]
            for atoms in frames:
                comment = atoms.info.get("comment", "")
                e = _parse_energy_from_comment(comment)
                energies.append(e)
            return energies
        except Exception:
            pass  # fallback al parser manual

    # Fallback: parser manual línea a línea
    text = conformers_xyz.read_text(errors="replace")
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            natoms = int(line)
        except ValueError:
            i += 1
            continue
        if i + 1 < len(lines):
            comment = lines[i + 1]
            e = _parse_energy_from_comment(comment)
        else:
            e = None
        energies.append(e)
        i += natoms + 2  # saltar el bloque completo
    return energies


def _parse_energy_from_comment(comment: str) -> Optional[float]:
    """Extrae energía de la línea de comentario de un frame XYZ de CREST."""
    comment = comment.strip()
    # Patrón "Etot= -123.456789"
    m = re.search(r"Etot[=\s]+([-\d.]+)", comment)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    # Número flotante negativo solo (formato CREST clásico)
    m = re.match(r"^\s*([-\d.]+)\s*$", comment)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    # xTB escribe a veces: "-123.456789 Eh"
    m = re.search(r"([-\d.]+)\s*Eh", comment)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    return None


# ─────────────────────────────────────────────
#  PARSER TERMOQUÍMICA (NUMFREQ ORCA 6)
# ─────────────────────────────────────────────

def parse_thermochem(freq_out: Path, T: float = T_STD) -> dict:
    """
    Extrae las correcciones termoquímicas del bloque THERMOCHEMISTRY de ORCA 6.

    ORCA 6 reporta un bloque completo a la temperatura especificada en el input
    (por defecto 298.15 K). Este parser extrae todas las contribuciones
    relevantes para el ciclo termodinámico compuesto:

        G_total_i = E_SP(TZVPD/CPCM) + G_corr(SVP/CPCM)

    donde G_corr = H_corr + ZPE - T·S_total se obtiene como:
        G_corr = G_total_freq - E_elec_freq
    siendo G_total_freq la 'Total Gibbs free energy' del NumFreq
    y E_elec_freq la energía electrónica del mismo NumFreq (def2-SVP).

    Retorna dict con las siguientes claves (todas en Hartree salvo indicación):
        E_elec_freq    : energía electrónica del NumFreq (referencia del ciclo)
        ZPE            : corrección de punto cero
        H_corr         : corrección entálpica total (ZPE + contribuciones térmicas)
        S_total        : entropía total (Hartree/K = Eh/K)
        G_corr         : H_corr - T·S_total  (corrección de Gibbs)
        G_total_freq   : E_elec_freq + G_corr (Gibbs total a nivel SVP)
        T_parsed       : temperatura a la que ORCA calculó la termoquímica (K)
        n_imag         : número de frecuencias imaginarias detectadas
        imag_freqs_cm1 : lista de valores imaginarios (cm⁻¹, negativos)
        ok             : True si todos los campos críticos fueron encontrados

    Si el archivo no existe o el bloque no se puede parsear, retorna un dict
    con ok=False y None en todos los campos numéricos.

    Notas sobre el formato ORCA 6:
    ─────────────────────────────
    El bloque de termoquímica en ORCA 6 (NumFreq + CPCM) tiene esta estructura:

        THERMOCHEMISTRY AT 298.15K
        Temperature         ...    298.15 K
        ...
        Zero point energy                   ...      0.123456 Eh
        Total thermal correction            ...      0.000789 Eh     (ZPE ya incluido)
        Total Enthalpy                      ...    -500.123456 Eh
        Final entropy term                  ...      0.045678 Eh/K  (= T·S / T, i.e. S)
        Total Gibbs free energy             ...    -500.234567 Eh

    ORCA reporta 'Final entropy term' como S (no T·S). Para obtener G:
        G = H - T·S  →  G_corr = H_corr - T * S_total
    donde H_corr = Total Enthalpy - E_elec, S_total = Final entropy term.
    """
    null = {
        "E_elec_freq":    None,
        "ZPE":            None,
        "H_corr":         None,
        "S_total":        None,
        "G_corr":         None,
        "G_total_freq":   None,
        "T_parsed":       None,
        "n_imag":         0,
        "imag_freqs_cm1": [],
        "ok":             False,
    }

    if not freq_out.exists():
        return null

    try:
        text = freq_out.read_text(errors="replace")
    except OSError:
        return null

    lines = text.splitlines()

    # ── 1. Frecuencias imaginarias ────────────────────────────────────────
    imag_freqs: list[float] = []
    for line in lines:
        m = re.match(r"\s+\d+:\s+(-\d+\.\d+)\s+cm", line)
        if m:
            imag_freqs.append(float(m.group(1)))

    # ── 2. Energía electrónica del NumFreq (FINAL SINGLE POINT ENERGY) ───
    # Tomamos la PRIMERA ocurrencia: corresponde al SP inicial de la freq,
    # antes de los desplazamientos numéricos.
    E_elec_freq: Optional[float] = None
    for line in lines:
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                E_elec_freq = float(line.split()[-1])
                break
            except ValueError:
                pass

    # ── 3. Bloque THERMOCHEMISTRY ─────────────────────────────────────────
    # Buscar el bloque más reciente (última ocurrencia, por si hay más de uno)
    thermo_start = -1
    for idx, line in enumerate(lines):
        if re.search(r"THERMOCHEMISTRY AT\s+[\d.]+\s*K", line, re.IGNORECASE):
            thermo_start = idx

    if thermo_start < 0:
        warnings.warn(
            f"[THERMOCHEM] No se encontró bloque THERMOCHEMISTRY en {freq_out}. "
            f"¿El cálculo terminó normalmente? ¿Es un NumFreq?"
        )
        result = dict(null)
        result.update({
            "E_elec_freq":    E_elec_freq,
            "n_imag":         len(imag_freqs),
            "imag_freqs_cm1": imag_freqs,
        })
        return result

    # Extraer temperatura del encabezado del bloque
    T_parsed: Optional[float] = None
    m = re.search(r"THERMOCHEMISTRY AT\s+([\d.]+)\s*K",
                  lines[thermo_start], re.IGNORECASE)
    if m:
        try:
            T_parsed = float(m.group(1))
        except ValueError:
            pass

    # Parsear líneas del bloque hasta el siguiente separador o fin de archivo
    ZPE:          Optional[float] = None
    H_total:      Optional[float] = None   # Total Enthalpy (Eh)  = E_elec + H_corr
    S_total:      Optional[float] = None   # Final entropy term   (Eh/K)
    G_total_freq: Optional[float] = None   # Total Gibbs free energy (Eh)

    # Patrones robustos para ORCA 6 (espaciado variable, unidades variables)
    _pat_zpe    = re.compile(r"Zero point energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_h      = re.compile(r"Total [Ee]nthalpy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_s      = re.compile(r"Final entropy term\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_g      = re.compile(r"Total Gibbs free energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)

    for line in lines[thermo_start:]:
        if ZPE is None:
            m = _pat_zpe.search(line)
            if m:
                try: ZPE = float(m.group(1))
                except ValueError: pass
        if H_total is None:
            m = _pat_h.search(line)
            if m:
                try: H_total = float(m.group(1))
                except ValueError: pass
        if S_total is None:
            m = _pat_s.search(line)
            if m:
                try: S_total = float(m.group(1))
                except ValueError: pass
        if G_total_freq is None:
            m = _pat_g.search(line)
            if m:
                try: G_total_freq = float(m.group(1))
                except ValueError: pass

    # ── 4. Derivar correcciones ───────────────────────────────────────────
    H_corr:  Optional[float] = None
    G_corr:  Optional[float] = None

    if E_elec_freq is not None and H_total is not None:
        H_corr = H_total - E_elec_freq

    if H_corr is not None and S_total is not None:
        # S_total viene en Eh/K  →  T·S en Eh
        G_corr = H_corr - T * S_total

    # Verificación de consistencia: G_corr debe reproducir G_total_freq ± 1e-6 Eh
    if (G_corr is not None and G_total_freq is not None
            and E_elec_freq is not None):
        G_reconstructed = E_elec_freq + G_corr
        discrepancy = abs(G_reconstructed - G_total_freq)
        if discrepancy > 1e-5:  # > 6 kcal/mol: algo muy mal
            warnings.warn(
                f"[THERMOCHEM] Discrepancia G reconstruida vs ORCA: "
                f"{discrepancy*HARTREE_TO_KCAL:.4f} kcal/mol en {freq_out}. "
                f"Revisar el output manualmente."
            )

    ok = all(x is not None for x in [ZPE, H_corr, S_total, G_corr, G_total_freq])

    return {
        "E_elec_freq":    E_elec_freq,
        "ZPE":            ZPE,
        "H_corr":         H_corr,
        "S_total":        S_total,          # Eh/K
        "G_corr":         G_corr,           # Eh, a temperatura T
        "G_total_freq":   G_total_freq,     # Eh
        "T_parsed":       T_parsed,
        "n_imag":         len(imag_freqs),
        "imag_freqs_cm1": imag_freqs,
        "ok":             ok,
    }


def parse_gas_phase_sp(sp_gas_out: Path) -> Optional[float]:
    """
    Extrae la energía de un SP en fase gas (sin CPCM).
    Usado en el ciclo de pKa para calcular ΔG_solv de cada especie.
    Retorna la última 'FINAL SINGLE POINT ENERGY' en Hartree, o None.
    """
    return parse_sp_energy(sp_gas_out)


# ─────────────────────────────────────────────
#  ANÁLISIS CREST → DFT
# ─────────────────────────────────────────────

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
    xtb_map: dict[str, Optional[float]] = {}
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


def print_crest_dft_analysis(df: pd.DataFrame):
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

    print(f"\n  ── Reordenamiento conformacional CREST (GFN2-xTB) → DFT ──")
    print(f"  Inversiones de rango detectadas: {n_inv}/{len(df)}")

    def fmt(x):
        if pd.isna(x):
            return "    —   "
        return f"{x:+8.3f}" if isinstance(x, float) else f"{x:>8}"

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
    print("\n  " + "  ".join(f"{header_map.get(c,c):>13}" for c in cols))
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


# ─────────────────────────────────────────────
#  CICLO TERMODINÁMICO COMPUESTO + BOLTZMANN
# ─────────────────────────────────────────────

def compute_G_total(df_conf: pd.DataFrame) -> pd.DataFrame:
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
    df["G_corr_298K"]     = np.nan
    df["H_corr"]          = np.nan
    df["S_total_EhK"]     = np.nan
    df["ZPE"]             = np.nan
    df["thermo_ok"]       = False
    df["thermo_warning"]  = ""

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

        H_corr  = tc["H_corr"]
        S_total = tc["S_total"]   # Eh/K
        ZPE     = tc.get("ZPE")

        df.at[idx, "H_corr"]       = H_corr
        df.at[idx, "S_total_EhK"]  = S_total
        df.at[idx, "ZPE"]          = ZPE
        df.at[idx, "thermo_ok"]    = True

        # G_corr a cada T  (aproximación RRHO, válida ±60 K respecto a T_ref)
        for T, label in zip(TEMPS_K, TEMP_LABELS):
            G_corr_T = H_corr - T * S_total
            G_total  = e_sp + G_corr_T
            df.at[idx, f"G_total_{label}"] = G_total

        # G_corr de referencia a 298.15 K
        df.at[idx, "G_corr_298K"] = H_corr - T_STD * S_total

    # ΔG relativo al mínimo para cada temperatura
    n_ok = df["thermo_ok"].sum()
    for T, label in zip(TEMPS_K, TEMP_LABELS):
        col_g   = f"G_total_{label}"
        col_dg  = f"dG_{label}_kcal"
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


def compute_boltzmann_multiT(df_conf: pd.DataFrame) -> pd.DataFrame:
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
            col_g  = f"G_total_{label}"
            col_dg = f"dG_{label}_kcal"
            col_b  = f"boltz_{label}"

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
        s_mix  = -(p * np.log(p)).sum()     # S_mix / R (adimensional)
        n_eff  = np.exp(s_mix)
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
        "gibbs_composite":  "ΔG compuesto [E_SP(TZVPD) + G_corr(SVP/CPCM)]",
        "electronic_only":  "ΔE electrónica pura — SIN ZPE/entropía ⚠",
        "unknown":          "modo desconocido",
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
            v   = row.get(col, np.nan)
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
        s   = df.attrs.get(f"S_mix_{lbl}")
        n   = df.attrs.get(f"N_eff_{lbl}")
        n1  = df.attrs.get(f"n_above_1pct_{lbl}")
        if s is not None:
            print(f"  {lbl:>5}:  S_mix/R = {s:.4f}  │  N_eff = {n:.2f}  │  "
                  f"confs >1% = {n1}/{len(df)}")

    # Sensibilidad a temperatura: ¿cambia el conf dominante entre extremos?
    b_low  = f"boltz_{TEMP_LABELS[0]}"    # 250K
    b_high = f"boltz_{TEMP_LABELS[-1]}"   # 310K
    if b_low in df.columns and b_high in df.columns:
        valid_low  = df[b_low].notna()
        valid_high = df[b_high].notna()
        if valid_low.any() and valid_high.any():
            dom_low  = df.loc[df[b_low].idxmax(),  "conformer"]
            dom_high = df.loc[df[b_high].idxmax(), "conformer"]
            if dom_low != dom_high:
                print(f"\n  ⚠  Confórmero dominante cambia con T:")
                print(f"     {TEMP_LABELS[0]} → {dom_low}   │   "
                      f"{TEMP_LABELS[-1]} → {dom_high}")
            else:
                print(f"\n  Confórmero dominante estable entre "
                      f"{TEMP_LABELS[0]}–{TEMP_LABELS[-1]}: {dom_low}")


# ─────────────────────────────────────────────
#  CICLO TERMODINÁMICO DE pKa
# ─────────────────────────────────────────────

def compute_pka(
    mol_acid: str,
    mol_base: str,
    T: float = T_STD,
) -> dict:
    """
    Calcula el pKa teórico de la reacción AH(aq) → A⁻(aq) + H⁺(aq) usando
    el ciclo termodinámico de Born (ciclo termoquímico de dos-ramas).

    ════════════════════════════════════════════════════════════════════════
    CICLO DE BORN (esquema de dos niveles)
    ════════════════════════════════════════════════════════════════════════

      AH(gas) ──────────────────────────────────────► A⁻(gas) + H⁺(gas)
         │                  ΔG_gas                          │
     ΔG_solv(AH)                              ΔG_solv(A⁻) + ΔG_solv(H⁺)
         │                                                  │
         ▼                                                  ▼
      AH(aq)  ──────────────────────────────────────► A⁻(aq) + H⁺(aq)
                              ΔG_aq

    ΔG_gas   = G_total(A⁻, gas) + G(H⁺, gas) - G_total(AH, gas)
    ΔG_solv  = E_SP(CPCM) - E_SP(gas)   para cada especie
    ΔG_aq    = ΔG_gas + ΔΔG_solv
    pKa      = ΔG_aq / (RT · ln10)

    ════════════════════════════════════════════════════════════════════════
    NIVEL DE TEORÍA
    ════════════════════════════════════════════════════════════════════════

    G_total en fase gas se aproxima como:
        G_gas(X) = E_SP_gas(X) + G_corr(X)
    donde G_corr viene del NumFreq (que usa CPCM). Esta es la aproximación
    de Ho & Coote (Theor. Chem. Acc. 2010): la corrección térmica es casi
    independiente del solvente para moléculas rígidas.

    ΔG_solv(X) = E_SP_aq(X, CPCM) - E_SP_gas(X)
    usando el mismo funcional/base en ambos (ωB97X-D3/def2-TZVPD).

    ════════════════════════════════════════════════════════════════════════
    CONSTANTES DE REFERENCIA (todas a 298.15 K)
    ════════════════════════════════════════════════════════════════════════

    G(H⁺, gas)    = -6.28 kcal/mol    [Tissandier et al. 1998]
    ΔG_solv(H⁺)   = -265.9 kcal/mol   [Kelly, Cramer & Truhlar 2006]
    δG_conc       = +1.894 kcal/mol    [corrección 1 atm → 1 M, RT·ln(24.5)]

    ════════════════════════════════════════════════════════════════════════
    PARÁMETROS
    ════════════════════════════════════════════════════════════════════════

    mol_acid : nombre del directorio en runs/ para AH (forma protonada)
    mol_base : nombre del directorio en runs/ para A⁻ (forma desprotonada)
    T        : temperatura de cálculo (default 298.15 K)

    ════════════════════════════════════════════════════════════════════════
    ARCHIVOS REQUERIDOS
    ════════════════════════════════════════════════════════════════════════

    Para cada mol (acid y base):
      · runs/<mol>/04_sp/<best_conf>/sp.out          — SP acuoso (CPCM)
      · runs/<mol>/06_sp_gas/<best_conf>/sp_gas.out  — SP fase gas (sin CPCM)
      · runs/<mol>/03_numfreq/<best_conf>/freq.out   — correcciones térmicas

    El confórmero representativo es el de menor G_total a 298 K.
    Si no existe 06_sp_gas, el ciclo de pKa no puede calcularse.

    ════════════════════════════════════════════════════════════════════════
    RETORNA
    ════════════════════════════════════════════════════════════════════════

    dict con:
        pKa              : valor teórico (float o None si falta algún dato)
        dG_aq_kcal       : ΔG_aq en kcal/mol
        dG_gas_kcal      : ΔG_gas en kcal/mol
        dG_solv_acid     : ΔG_solv(AH) en kcal/mol
        dG_solv_base     : ΔG_solv(A⁻) en kcal/mol
        dG_solv_hplus    : ΔG_solv(H⁺) constante de literatura
        G_total_acid_gas : G_total(AH) en fase gas (Hartree)
        G_total_base_gas : G_total(A⁻) en fase gas (Hartree)
        conf_acid        : confórmero representativo de AH
        conf_base        : confórmero representativo de A⁻
        T                : temperatura usada (K)
        ok               : True si el pKa se calculó correctamente
        warnings         : lista de advertencias
    """
    result: dict = {
        "pKa":              None,
        "dG_aq_kcal":       None,
        "dG_gas_kcal":      None,
        "dG_solv_acid":     None,
        "dG_solv_base":     None,
        "dG_solv_hplus":    DG_SOLV_HPLUS,
        "G_total_acid_gas": None,
        "G_total_base_gas": None,
        "conf_acid":        None,
        "conf_base":        None,
        "T":                T,
        "ok":               False,
        "warnings":         [],
    }

    def _warn(msg: str):
        result["warnings"].append(msg)
        warnings.warn(f"[pKa] {msg}")

    # ── Cargar resultados de ambas formas ─────────────────────────────────
    try:
        df_acid, _ = collect_results(mol_acid)
        df_base, _ = collect_results(mol_base)
    except Exception as e:
        _warn(f"Error al cargar resultados: {e}")
        return result

    if df_acid.empty or df_base.empty:
        _warn("DataFrames vacíos para ácido o base.")
        return result

    # ── Confórmero representativo: mínimo G_total a 298.15 K ─────────────
    lbl_298 = "298K"
    col_g298 = f"G_total_{lbl_298}"

    def _best_conf(df: pd.DataFrame, mol: str) -> Optional[str]:
        if col_g298 in df.columns and df["thermo_ok"].any():
            valid = df["thermo_ok"] & df[col_g298].notna()
            if valid.any():
                return df.loc[df.loc[valid, col_g298].idxmin(), "conformer"]
        # Fallback: mínima E_SP
        if "E_hartree" in df.columns and df["E_hartree"].notna().any():
            _warn(f"[{mol}] Usando E_SP (sin G_corr) para seleccionar confórmero.")
            return df.loc[df["E_hartree"].idxmin(), "conformer"]
        return None

    conf_acid = _best_conf(df_acid, mol_acid)
    conf_base = _best_conf(df_base, mol_base)
    result["conf_acid"] = conf_acid
    result["conf_base"] = conf_base

    if conf_acid is None or conf_base is None:
        _warn("No se pudo determinar el confórmero representativo.")
        return result

    # ── Energías SP acuosas (CPCM) ───────────────────────────────────────
    acid_sp_aq_file = WORK_ROOT / mol_acid / "04_sp" / conf_acid / "sp.out"
    base_sp_aq_file = WORK_ROOT / mol_base / "04_sp" / conf_base / "sp.out"

    E_acid_aq = parse_sp_energy(acid_sp_aq_file)
    E_base_aq = parse_sp_energy(base_sp_aq_file)

    if E_acid_aq is None:
        _warn(f"No se encontró SP acuoso de AH en {acid_sp_aq_file}")
    if E_base_aq is None:
        _warn(f"No se encontró SP acuoso de A⁻ en {base_sp_aq_file}")

    # ── Energías SP en fase gas ───────────────────────────────────────────
    acid_sp_gas_file = (WORK_ROOT / mol_acid / PKA_GAS_SP_SUBDIR
                        / conf_acid / "sp_gas.out")
    base_sp_gas_file = (WORK_ROOT / mol_base / PKA_GAS_SP_SUBDIR
                        / conf_base / "sp_gas.out")

    E_acid_gas = parse_gas_phase_sp(acid_sp_gas_file)
    E_base_gas = parse_gas_phase_sp(base_sp_gas_file)

    if E_acid_gas is None:
        _warn(
            f"SP en fase gas de AH no encontrado en {acid_sp_gas_file}. "
            f"Corre: orca sp_gas.inp > sp_gas.out (sin CPCM, mismo funcional/base)."
        )
    if E_base_gas is None:
        _warn(
            f"SP en fase gas de A⁻ no encontrado en {base_sp_gas_file}. "
            f"Corre: orca sp_gas.inp > sp_gas.out (sin CPCM, mismo funcional/base)."
        )

    if E_acid_gas is None or E_base_gas is None:
        return result  # no hay suficiente para continuar

    # ── Correcciones térmicas (G_corr) ────────────────────────────────────
    acid_freq_file = WORK_ROOT / mol_acid / "03_numfreq" / conf_acid / "freq.out"
    base_freq_file = WORK_ROOT / mol_base / "03_numfreq" / conf_base / "freq.out"

    tc_acid = parse_thermochem(acid_freq_file, T=T)
    tc_base = parse_thermochem(base_freq_file, T=T)

    G_corr_acid: Optional[float] = None
    G_corr_base: Optional[float] = None

    if tc_acid.get("ok"):
        G_corr_acid = tc_acid["G_corr"]
    else:
        _warn(
            f"Correcciones térmicas de AH no disponibles o incompletas "
            f"({acid_freq_file}). Se usará G_gas ≈ E_gas (sin ZPE ni entropía). "
            f"El pKa resultante tendrá mayor incertidumbre."
        )

    if tc_base.get("ok"):
        G_corr_base = tc_base["G_corr"]
    else:
        _warn(
            f"Correcciones térmicas de A⁻ no disponibles o incompletas "
            f"({base_freq_file}). Se usará G_gas ≈ E_gas."
        )

    # ── G_total en fase gas para cada especie ─────────────────────────────
    #   G_gas = E_SP_gas + G_corr(RRHO, de NumFreq con CPCM)
    #   La aproximación Ho & Coote justifica usar G_corr del cálculo CPCM
    #   incluso para la fase gas (la corrección térmica es casi la misma).
    G_acid_gas = E_acid_gas + (G_corr_acid if G_corr_acid is not None else 0.0)
    G_base_gas = E_base_gas + (G_corr_base if G_corr_base is not None else 0.0)

    result["G_total_acid_gas"] = G_acid_gas
    result["G_total_base_gas"] = G_base_gas

    # ── ΔG_gas (en kcal/mol) ─────────────────────────────────────────────
    #   ΔG_gas = G(A⁻, gas) + G(H⁺, gas) - G(AH, gas)
    dG_gas_hartree = G_base_gas - G_acid_gas   # en Hartree (H⁺ se suma aparte)
    dG_gas_kcal    = dG_gas_hartree * HARTREE_TO_KCAL + G_HPLUS_GAS
    result["dG_gas_kcal"] = dG_gas_kcal

    # ── ΔG_solv para cada especie ─────────────────────────────────────────
    #   ΔG_solv(X) = E_SP_aq(X, CPCM) - E_SP_gas(X)
    #   En kcal/mol. Incluye corrección de concentración estándar δG_conc.
    if E_acid_aq is not None:
        dG_solv_acid = (E_acid_aq - E_acid_gas) * HARTREE_TO_KCAL + DG_CONC_CORR
        result["dG_solv_acid"] = dG_solv_acid
    else:
        _warn("No se puede calcular ΔG_solv(AH) sin E_SP acuoso.")
        return result

    if E_base_aq is not None:
        dG_solv_base = (E_base_aq - E_base_gas) * HARTREE_TO_KCAL + DG_CONC_CORR
        result["dG_solv_base"] = dG_solv_base
    else:
        _warn("No se puede calcular ΔG_solv(A⁻) sin E_SP acuoso.")
        return result

    # ── ΔG_aq (ciclo de Born completo) ───────────────────────────────────
    #   ΔG_aq = ΔG_gas + ΔG_solv(A⁻) + ΔG_solv(H⁺) - ΔG_solv(AH)
    dG_aq_kcal = (
        dG_gas_kcal
        + dG_solv_base
        + DG_SOLV_HPLUS
        - dG_solv_acid
    )
    result["dG_aq_kcal"] = dG_aq_kcal

    # ── pKa ──────────────────────────────────────────────────────────────
    pKa = dG_aq_kcal / (R_KCAL * T * LN10)
    result["pKa"]  = pKa
    result["ok"]   = True

    return result


def print_pka_analysis(pka_result: dict, mol_acid: str, mol_base: str):
    """Imprime el análisis completo del ciclo termodinámico de pKa."""
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  ANÁLISIS DE pKa TEÓRICO")
    print(f"  AH = {mol_acid}   │   A⁻ = {mol_base}")
    print(f"  T  = {pka_result['T']:.2f} K")
    print(sep)

    if not pka_result["ok"]:
        print("\n  ✗  El ciclo de pKa no pudo completarse.")
        for w in pka_result["warnings"]:
            print(f"     ⚠  {w}")
        return

    print(f"\n  Confórmero representativo AH : {pka_result['conf_acid']}")
    print(f"  Confórmero representativo A⁻ : {pka_result['conf_base']}")

    print(f"\n  ── Contribuciones al ciclo de Born ──")
    print(f"  {'ΔG_gas':>30} = {pka_result['dG_gas_kcal']:+10.3f} kcal/mol")
    print(f"  {'ΔG_solv(AH)':>30} = {pka_result['dG_solv_acid']:+10.3f} kcal/mol")
    print(f"  {'ΔG_solv(A⁻)':>30} = {pka_result['dG_solv_base']:+10.3f} kcal/mol")
    print(f"  {'ΔG_solv(H⁺) [Kelly 2006]':>30} = {pka_result['dG_solv_hplus']:+10.3f} kcal/mol")
    print(f"  {'─'*48}")
    print(f"  {'ΔG_aq':>30} = {pka_result['dG_aq_kcal']:+10.3f} kcal/mol")
    print(f"\n  pKa (teórico, {pka_result['T']:.2f} K) = {pka_result['pKa']:.2f}")
    print(f"\n  Nota: incertidumbre típica ±1–2 unidades de pKa para ciclos")
    print(f"  Born/CPCM a nivel ωB97X-D3/def2-TZVPD (Ho & Coote 2010).")

    if pka_result["warnings"]:
        print(f"\n  Advertencias:")
        for w in pka_result["warnings"]:
            print(f"    ⚠  {w}")

    print(f"\n{sep}\n")


# ─────────────────────────────────────────────
#  COLECTOR PRINCIPAL
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

    # Boltzmann multi-temperatura
    df_conf = compute_boltzmann_multiT(df_conf)

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
#  ESPECTROS UV-Vis
# ─────────────────────────────────────────────

def _build_spectrum(df_tddft: pd.DataFrame, df_conf: pd.DataFrame,
                    sigma: float, x: np.ndarray,
                    temp_label: str = "298K") -> np.ndarray:
    """Construye el espectro Gaussiano-suavizado Boltzmann-promediado."""
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


def plot_spectrum(df_tddft: pd.DataFrame, df_conf: pd.DataFrame, mol_name: str,
                  sigma: float = 20.0, x_range: tuple = (200, 500)):
    """Espectro UV-Vis individual: Boltzmann-promediado a múltiples temperaturas."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
    except ImportError:
        print("matplotlib no disponible — saltando plot")
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

    out_png = WORK_ROOT / mol_name / "uvvis_spectrum.png"
    plt.savefig(out_png, dpi=150)
    print(f"  Espectro guardado en: {out_png}")
    plt.show()


def plot_comparison(dfs_conf: dict, dfs_tddft: dict,
                    sigma: float = 20.0, x_range: tuple = (200, 500)):
    """
    Superpone espectros UV-Vis (298 K) de todos los estados de protonación,
    normalizados al máximo individual para facilitar comparación de forma.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib no disponible — saltando plot")
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
    out_png = WORK_ROOT / "comparison_uvvis.png"
    plt.savefig(out_png, dpi=150)
    print(f"  Espectro comparativo guardado en: {out_png}")
    plt.show()


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


# ─────────────────────────────────────────────
#  PUNTO DE ENTRADA
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Recolecta y analiza resultados ORCA: energías, "
                    "Boltzmann multi-T, reordenamiento CREST→DFT, "
                    "y comparación entre estados de protonación.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  python results_collector.py --mol mglyval
  python results_collector.py --mol mglyval --export xlsx --plot
  python results_collector.py --compare
  python results_collector.py --compare --plot --ref mglyval_neutral
        """,
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--mol",     metavar="NAME",
                      help="analizar un estado (nombre del directorio en runs/)")
    mode.add_argument("--compare", action="store_true",
                      help="comparar todos los estados de protonación en runs/")

    parser.add_argument("--ref",    metavar="NAME", default=None,
                        help="estado de referencia para ΔE relativa (solo con --compare)")
    parser.add_argument("--export", choices=["csv", "xlsx", "both"], default=None,
                        help="exportar tablas")
    parser.add_argument("--plot",   action="store_true",
                        help="generar espectro(s) UV-Vis")
    parser.add_argument("--sigma",  type=float, default=20.0,
                        help="anchura Gaussiana en nm para los plots (default: 20)")
    parser.add_argument("--xmin",   type=float, default=200.0,
                        help="límite inferior del eje λ en nm (default: 200)")
    parser.add_argument("--xmax",   type=float, default=500.0,
                        help="límite superior del eje λ en nm (default: 500)")
    args = parser.parse_args()

    x_range = (args.xmin, args.xmax)

    if args.mol:
        df_conf, df_tddft = collect_results(args.mol)
        print_summary(df_conf, df_tddft, args.mol)

        if args.export:
            export_results(df_conf, df_tddft, args.mol, args.export)
        if args.plot:
            plot_spectrum(df_tddft, df_conf, args.mol,
                          sigma=args.sigma, x_range=x_range)

    elif args.compare:
        result = compare_protonation_states(ref_mol=args.ref)
        if isinstance(result, tuple):
            df_cmp, dfs_conf, dfs_tddft = result
        else:
            df_cmp = result
            dfs_conf = dfs_tddft = {}

        print_comparison(df_cmp, ref_mol=args.ref)

        # También imprimir resumen individual de cada estado
        for mol in sorted(dfs_conf.keys()):
            print_summary(dfs_conf[mol], dfs_tddft.get(mol, pd.DataFrame()), mol)

        if args.export:
            export_comparison(df_cmp, args.export)
            for mol in sorted(dfs_conf.keys()):
                if not dfs_conf[mol].empty:
                    export_results(dfs_conf[mol],
                                   dfs_tddft.get(mol, pd.DataFrame()),
                                   mol, args.export)
        if args.plot:
            plot_comparison(dfs_conf, dfs_tddft,
                            sigma=args.sigma, x_range=x_range)


if __name__ == "__main__":
    main()
