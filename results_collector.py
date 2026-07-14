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

⚠️ Nota de parsing: ORCA 6 imprime 'Final entropy term' como T·S YA
EVALUADO a la T del bloque (no como S puro en Eh/K). parse_thermochem()
deriva S_total = TS_total / T_parsed antes de usarlo en cualquier fórmula
de G_corr; ver docstring de esa función para el detalle.

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
        TS_total       : término T·S a T_parsed, tal cual lo imprime ORCA
                         ('Final entropy term'), en Eh
        S_total        : entropía pura derivada, S = TS_total / T_parsed (Eh/K)
        G_corr         : H_corr - TS_total  (corrección de Gibbs a T_parsed)
        G_total_freq   : E_elec_freq + G_corr (Gibbs total a nivel SVP)
        T_parsed       : temperatura a la que ORCA calculó la termoquímica (K)
        n_imag         : número de frecuencias imaginarias detectadas
        imag_freqs_cm1 : lista de valores imaginarios (cm⁻¹, negativos)
        ok             : True si todos los campos críticos fueron encontrados

    Si el archivo no existe o el bloque no se puede parsear, retorna un dict
    con ok=False y None en todos los campos numéricos.

    Notas sobre el formato ORCA 6:
    ─────────────────────────────
    El bloque de termoquímica en ORCA 6 (NumFreq + CPCM, con Quasi-RRHO) tiene
    esta estructura:

        THERMOCHEMISTRY AT 298.15K
        Zero point energy                   ...      0.123456 Eh
        Total thermal correction                     0.000789 Eh     (ZPE ya incluido)
        Total Enthalpy                      ...    -500.123456 Eh
        Final entropy term                  ...      0.045678 Eh   48.42 kcal/mol
        Final Gibbs free energy             ...    -500.234567 Eh

    ⚠️ IMPORTANTE: 'Final entropy term' de ORCA 6 NO es la entropía S en
    Eh/K. Es el término T·S ya evaluado a la temperatura del bloque
    (T_ref = 298.15 K por defecto), expresado en Eh. Se puede confirmar
    dividiendo el valor en kcal/mol que ORCA imprime al lado: por ejemplo
    0.07716168 Eh × 627.509 kcal/mol/Eh ≈ 48.42 kcal/mol, que es exactamente
    el valor que ORCA muestra en esa misma línea — si fuera S puro (Eh/K),
    al multiplicar por T=298.15 quedaría ~23 Eh, muy por encima del rango
    físico razonable para esta corrección.

    Por lo tanto:
        G = H - T·S  →  G_corr = H_corr - TS_ref
    donde H_corr = Total Enthalpy - E_elec, TS_ref = 'Final entropy term'
    (ya es T·S a T_ref, NO multiplicar por T de nuevo).

    Si se necesita S puro en Eh/K (p.ej. para recalcular G_corr a otra T),
    se debe derivar como S = TS_ref / T_ref.
    """
    null = {
        "E_elec_freq":    None,
        "ZPE":            None,
        "H_corr":         None,
        "TS_total":       None,
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
    TS_total:     Optional[float] = None   # Final entropy term = T·S a T_parsed (Eh)
    G_total_freq: Optional[float] = None   # Total/Final Gibbs free energy (Eh)

    # Patrones robustos para ORCA 6 (espaciado variable, unidades variables)
    _pat_zpe    = re.compile(r"Zero point energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_h      = re.compile(r"Total [Ee]nthalpy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_s      = re.compile(r"Final entropy term\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    # ORCA 6.1.1 con Quasi-RRHO (Quasi RRHO ... True) imprime esta línea como
    # "Final Gibbs free energy", no "Total Gibbs free energy". Sin Quasi-RRHO
    # (o en otras versiones/anotaciones) puede aparecer la otra forma.
    # Aceptamos ambas para que el parser sea robusto a los dos modos.
    _pat_g      = re.compile(r"(?:Total|Final) Gibbs free energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)

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
        if TS_total is None:
            m = _pat_s.search(line)
            if m:
                try: TS_total = float(m.group(1))
                except ValueError: pass
        if G_total_freq is None:
            m = _pat_g.search(line)
            if m:
                try: G_total_freq = float(m.group(1))
                except ValueError: pass

    # ── 4. Derivar correcciones ───────────────────────────────────────────
    H_corr:  Optional[float] = None
    G_corr:  Optional[float] = None
    S_total: Optional[float] = None   # entropía pura derivada (Eh/K)

    if E_elec_freq is not None and H_total is not None:
        H_corr = H_total - E_elec_freq

    if TS_total is not None and T_parsed is not None and T_parsed != 0:
        # TS_total ('Final entropy term') es T·S evaluado a T_parsed (Eh).
        # Despejamos S puro para poder reescalar a la T solicitada por el
        # caller (parámetro T), que puede diferir de T_parsed.
        S_total = TS_total / T_parsed

    if H_corr is not None and S_total is not None:
        # G_corr a la temperatura T solicitada (no necesariamente T_parsed).
        # Si T == T_parsed, esto reproduce exactamente G_corr = H_corr - TS_total.
        G_corr = H_corr - T * S_total

    # Verificación de consistencia (solo válida si T == T_parsed, ya que
    # G_total_freq de ORCA está fijo a T_parsed):
    if (G_corr is not None and G_total_freq is not None
            and E_elec_freq is not None and T_parsed is not None
            and abs(T - T_parsed) < 1e-6):
        G_reconstructed = E_elec_freq + G_corr
        discrepancy = abs(G_reconstructed - G_total_freq)
        if discrepancy > 1e-5:  # > 6 kcal/mol: algo muy mal
            warnings.warn(
                f"[THERMOCHEM] Discrepancia G reconstruida vs ORCA: "
                f"{discrepancy*HARTREE_TO_KCAL:.4f} kcal/mol en {freq_out}. "
                f"Revisar el output manualmente."
            )

    ok = all(x is not None for x in [ZPE, H_corr, TS_total, G_corr, G_total_freq])

    return {
        "E_elec_freq":    E_elec_freq,
        "ZPE":            ZPE,
        "H_corr":         H_corr,
        "TS_total":       TS_total,          # Eh, T·S a T_parsed (tal cual ORCA)
        "S_total":        S_total,           # Eh/K, derivado = TS_total / T_parsed
        "G_corr":         G_corr,            # Eh, a la temperatura T solicitada
        "G_total_freq":   G_total_freq,      # Eh
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

    # ── Correcciones térmicas ORCA (NumFreq) por confórmero ────────────────
    # compute_boltzmann_multiT() necesita la columna "thermo_ok" (creada por
    # compute_G_total) para usar el ciclo de Gibbs compuesto; sin este paso
    # caía siempre al modo "electronic_only" (ΔE pura) aunque el NumFreq
    # se hubiera parseado bien.
    _thermochem_cache: list[Optional[dict]] = []
    for _, row in df_conf.iterrows():
        conf_id  = row["conformer"]
        freq_out = mol_dir / "03_numfreq" / conf_id / "freq.out"
        tc = parse_thermochem(freq_out, T=T_STD)
        _thermochem_cache.append(tc if tc.get("ok") else None)
    df_conf["_thermochem"] = _thermochem_cache
    df_conf = compute_G_total(df_conf)

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

        # Marcar λ_max de S1 (Boltzmann-promediado, 298 K) — NO el argmax del
        # espectro combinado, que mezcla todos los estados (S1, S2, S3...)
        # ponderados por fosc y puede quedar dominado por un estado superior
        # de mayor fuerza de oscilador. Usamos el mismo valor que va al xlsx
        # (df_conf.attrs["lambda_max_boltz_298K"], calculado solo con S1) para
        # que la figura y la tabla exportada sean consistentes.
        lmax = df_c.attrs.get("lambda_max_boltz_298K")
        if lmax is None:
            # Fallback: recalcular filtrando a S1 si el atributo no está disponible
            s1 = df_t[df_t["state"] == 1].copy()
            if not s1.empty:
                s1 = s1.merge(df_c[["conformer", "boltz_298K"]], on="conformer", how="left")
                w = s1["boltz_298K"].fillna(0)
                lmax = float((s1["wavelength_nm"] * w).sum() / w.sum()) if w.sum() > 0 else None

        if lmax is not None:
            lambda_maxs.append((mol, lmax, color))
            ax_spec.axvline(lmax, color=color, lw=0.7, ls="--", alpha=0.6)

    ax_spec.set_xlabel("Longitud de onda (nm)", fontsize=12)
    ax_spec.set_ylabel("Intensidad normalizada", fontsize=12)
    ax_spec.set_title(f"Comparación UV-Vis — Estados de protonación  "
                      f"(σ = {sigma} nm, T = 298 K; curva = todos los estados, "
                      f"líneas punteadas = λ_max de S1)", fontsize=13)
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


# ═══════════════════════════════════════════════════════════════════════════════
#  APÉNDICE GOODVIBES — corrección cuasi-armónica Grimme mRRHO
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Fundamento:
#  ORCA calcula G usando el tratamiento rígido-rotor/oscilador-armónico (RRHO)
#  estándar, que sobreestima la contribución entrópica de modos de baja
#  frecuencia (< ~100 cm⁻¹), produciendo G demasiado negativo.
#  GoodVibes aplica la corrección cuasi-RRHO de Grimme (J. Chem. Theory Comput.
#  2012, 8, 4165): por debajo del cutoff f₀ (default 100 cm⁻¹), la entropía
#  vibracional se sustituye por una interpolación entre el oscilador armónico
#  y un rotor libre, suavizando el error sistemático.
#
#  Lo que se compara en este apéndice:
#
#    G_ORCA(i) = E_SP(TZVPD) + G_corr_ORCA(RRHO, SVP/CPCM)
#    G_GV(i)   = E_SP(TZVPD) + G_corr_GV(mRRHO Grimme, SVP/CPCM)
#
#    ΔG_corr(i)   = G_corr_GV - G_corr_ORCA           [Eh]  — diferencia pura
#    ΔΔG_rel(i,T) = ΔG_GV_rel(i,T) - ΔG_ORCA_rel(i,T) [kcal/mol] — impacto en ranking
#    Δboltz(i,T)  = boltz_GV(i,T) - boltz_ORCA(i,T)    — impacto en Boltzmann
#
#  Flujo:
#    run_goodvibes_on_mol()  → invoca `python -m goodvibes` como subprocess,
#                              guarda JSON en runs/<mol>/<conf>/goodvibes/,
#                              retorna dict conf_id → GV_thermochem
#    compute_gv_appendix()  → construye DataFrame de comparación
#    print_gv_appendix()    → tabla formateada al estilo del script
#    export_gv_appendix()   → guarda CSV/XLSX del apéndice en runs/<mol>/
#
# ═══════════════════════════════════════════════════════════════════════════════

import json
import subprocess
import sys


# ── Parámetros GoodVibes (modificables) ─────────────────────────────────────
GV_FREQ_CUTOFF = 100.0    # cm⁻¹ — cutoff para corrección mRRHO Grimme
GV_QS          = "grimme" # esquema quasi-armónico: "grimme" o "truhlar"
GV_SUBDIR      = "goodvibes"  # subcarpeta creada dentro de cada conf_xxx

# Entorno conda donde está instalado GoodVibes.
# Si se deja en None, el script lo detecta automáticamente buscando en los
# entornos conocidos de miniforge3/miniconda3/anaconda3.
GV_CONDA_ENV   = "sciencebitch"


def _resolve_gv_python() -> str:
    """
    Devuelve la ruta absoluta al intérprete Python que tiene goodvibes instalado.

    Estrategia (en orden de prioridad):
      1. Si GV_CONDA_ENV está definido, buscar ese entorno en los prefijos
         estándar de miniforge3, miniconda3 y anaconda3 en /home/kmilozeta.
      2. Intentar `python -m goodvibes --version` con el ejecutable actual
         (sys.executable) por si el entorno ya lo incluye.
      3. Buscar `goodvibes` en el PATH como script de entrada.
      4. Fallback: devolver sys.executable (fallará con el error habitual).
    """
    import shutil

    # ── 1. Entorno conda explícito ────────────────────────────────────────
    if GV_CONDA_ENV:
        home = Path.home()
        for base in ["miniforge3", "miniconda3", "anaconda3", "conda"]:
            candidate = home / base / "envs" / GV_CONDA_ENV / "bin" / "python"
            if candidate.exists():
                # Verificar que goodvibes esté realmente disponible allí
                try:
                    r = subprocess.run(
                        [str(candidate), "-c", "import goodvibes"],
                        capture_output=True, timeout=10
                    )
                    if r.returncode == 0:
                        return str(candidate)
                except Exception:
                    pass

    # ── 2. Python actual ──────────────────────────────────────────────────
    try:
        r = subprocess.run(
            [sys.executable, "-c", "import goodvibes"],
            capture_output=True, timeout=10
        )
        if r.returncode == 0:
            return sys.executable
    except Exception:
        pass

    # ── 3. goodvibes en PATH (instalado como script de consola) ──────────
    gv_script = shutil.which("goodvibes")
    if gv_script:
        # Derivar el Python del mismo directorio bin/
        py_candidate = Path(gv_script).parent / "python"
        if py_candidate.exists():
            return str(py_candidate)

    # ── 4. Fallback ───────────────────────────────────────────────────────
    return sys.executable


def _probe_gv_capabilities(python_exe: str) -> dict:
    """
    Ejecuta `python -m goodvibes --help` para detectar qué flags están
    disponibles en la versión instalada.

    Retorna un dict con:
        supports_json_format  : bool — existe --output_format json (v4+)
        supports_csv          : bool — existe --csv (v3.x)
        version               : str  — número de versión o "?"
    """
    caps = {"supports_json_format": False, "supports_csv": False, "version": "?"}
    try:
        r = subprocess.run(
            [python_exe, "-m", "goodvibes", "--help"],
            capture_output=True, text=True, timeout=15
        )
        help_text = r.stdout + r.stderr

        # Detectar flags de salida disponibles
        if "--output_format" in help_text or "output_format" in help_text:
            caps["supports_json_format"] = True
        if "--csv" in help_text:
            caps["supports_csv"] = True

        # Intentar versión
        r2 = subprocess.run(
            [python_exe, "-m", "goodvibes", "--version"],
            capture_output=True, text=True, timeout=10
        )
        for line in (r2.stdout + r2.stderr).splitlines():
            m = re.search(r"([\d]+\.[\d]+\.[\d]+|v[\d.]+)", line)
            if m:
                caps["version"] = m.group(1).lstrip("v")
                break
    except Exception:
        pass
    return caps


def _parse_goodvibes_dat(dat_path: Path) -> dict:
    """
    Parsea el archivo de texto tabulado que GoodVibes genera siempre
    (Goodvibes_output.dat o similar) cuando no se solicita salida JSON.

    GoodVibes imprime una tabla con el formato (columnas separadas por espacios):
        Structure  E_SPC  ZPE  H  T.S  T.qh-S  G(T)  qh-G(T)
    o bien (versión 3.x, con nombres ligeramente distintos):
        Structure  E  ZPE  H  T.S  T.qh-S  G(T)  qh-G(T)

    Las unidades son Hartree para energías y Hartree/K para S si se reporta.
    Las líneas de datos comienzan con el nombre del archivo (sin extensión).

    Retorna el mismo formato que _parse_goodvibes_json():
        H_GV, S_GV (Eh/K), TS_GV (Eh), G_GV (Eh), ZPE_GV (Eh),
        E_elec_GV (Eh), T_GV (K), G_corr_GV (Eh), ok, warning
    """
    null = {
        "H_GV": None, "S_GV": None, "TS_GV": None, "G_GV": None,
        "ZPE_GV": None, "E_elec_GV": None, "T_GV": None,
        "G_corr_GV": None, "ok": False, "warning": "dat no parseado",
    }

    if not dat_path.exists():
        null["warning"] = f"Archivo .dat no encontrado: {dat_path}"
        return null

    try:
        text = dat_path.read_text(errors="replace")
    except OSError as exc:
        null["warning"] = f"No se pudo leer {dat_path}: {exc}"
        return null

    lines = text.splitlines()

    # ── Detectar temperatura del encabezado ──────────────────────────────
    T_GV: Optional[float] = None
    for line in lines:
        m = re.search(r"Temp\s*=\s*([\d.]+)\s*K", line, re.I)
        if not m:
            m = re.search(r"at\s+([\d.]+)\s*K", line, re.I)
        if m:
            try:
                T_GV = float(m.group(1))
                break
            except ValueError:
                pass
    if T_GV is None:
        T_GV = 298.15  # asumir condición estándar si no se detecta

    # ── Detectar línea de encabezado de columnas ─────────────────────────
    # GoodVibes imprime una línea con "Structure" y los nombres de columnas.
    # Buscamos esa línea y la siguiente que sea numérica (datos).
    header_idx = -1
    for i, line in enumerate(lines):
        if re.search(r"\bStructure\b", line, re.I):
            header_idx = i
            break

    if header_idx < 0:
        null["warning"] = f"No se encontró encabezado 'Structure' en {dat_path}"
        return null

    header = lines[header_idx].split()

    # Normalizar nombres de columna a minúsculas para búsqueda robusta
    def _col(names: list[str]) -> Optional[int]:
        """Devuelve el índice (0-based) de la primera columna que coincida."""
        for name in names:
            nl = name.lower()
            for j, h in enumerate(header):
                if h.lower() == nl:
                    return j
        return None

    # Columnas posibles según versión de GoodVibes
    idx_e    = _col(["E_SPC", "E_SP", "E", "E_DFT"])
    idx_zpe  = _col(["ZPE"])
    idx_h    = _col(["H"])
    idx_ts   = _col(["T.S", "TS", "T*S"])       # T·S a T_GV (Eh)
    idx_tqs  = _col(["T.qh-S", "T.QH-S"])       # T·qhS (corrección GoodVibes)
    idx_g    = _col(["G(T)", "G"])               # G RRHO (sin corrección quasi)
    idx_qhg  = _col(["qh-G(T)", "qh-G", "QH-G(T)"])  # G quasi-RRHO (el que queremos)

    # Buscar primera línea de datos después del encabezado
    # NOTA: GoodVibes v4 antepone un símbolo de estado (p.ej. "o", "x", "?")
    # ANTES del nombre de la estructura, produciendo una columna extra al inicio
    # que NO existe en el encabezado. Ejemplo real:
    #   encabezado : Structure   E   ZPE   H   qh-H   T.S   T.qh-S   G(T)   qh-G(T)
    #   datos      : o   freq   -1220.967288   0.403359   ...
    # Los índices derivados del encabezado (base 0) deben desplazarse +1 para
    # las líneas de datos cuando se detecta este símbolo inicial.
    E_elec: Optional[float] = None
    ZPE:    Optional[float] = None
    H_GV:   Optional[float] = None
    TS_GV:  Optional[float] = None   # T·qhS (Eh) — con corrección cuasi-RRHO
    G_GV:   Optional[float] = None   # qh-G(T) absoluto (Eh)

    for line in lines[header_idx + 1:]:
        parts = line.split()
        if len(parts) < 3:
            continue

        # Detectar offset: si parts[0] es un símbolo de 1 carácter (o/x/?)
        # y parts[1] no es numérico (es el nombre del archivo), hay offset=1.
        # Si parts[0] no es numérico directamente (es el nombre), offset=0.
        offset = 0
        if len(parts[0]) == 1 and not parts[0].lstrip('-').replace('.','',1).isdigit():
            # parts[0] es símbolo de estado → parts[1] es nombre, datos desde parts[2]
            try:
                float(parts[2])
                offset = 1
            except (ValueError, IndexError):
                continue
        else:
            # parts[0] es nombre → datos desde parts[1]
            try:
                float(parts[1])
                offset = 0
            except (ValueError, IndexError):
                continue

        def _get_col(idx: Optional[int]) -> Optional[float]:
            if idx is None:
                return None
            real_idx = idx + offset
            if real_idx >= len(parts):
                return None
            try:
                return float(parts[real_idx])
            except ValueError:
                return None

        E_elec = _get_col(idx_e)
        ZPE    = _get_col(idx_zpe)
        H_GV   = _get_col(idx_h)
        TS_GV  = _get_col(idx_tqs) or _get_col(idx_ts)  # preferir quasi-RRHO
        G_GV   = _get_col(idx_qhg) or _get_col(idx_g)   # preferir quasi-RRHO
        break  # solo necesitamos la primera fila de datos

    if G_GV is None:
        null["warning"] = (
            f"No se pudo extraer G(T) del .dat ({dat_path}). "
            f"Columnas detectadas: {header}"
        )
        return null

    S_GV: Optional[float] = None
    if TS_GV is not None and T_GV > 0:
        S_GV = TS_GV / T_GV  # Eh/K

    G_corr_GV: Optional[float] = None
    if G_GV is not None and E_elec is not None:
        G_corr_GV = G_GV - E_elec

    return {
        "H_GV":      H_GV,
        "S_GV":      S_GV,        # Eh/K
        "TS_GV":     TS_GV,       # Eh (T·qhS)
        "G_GV":      G_GV,        # Eh absoluto (qh-G(T))
        "ZPE_GV":    ZPE,         # Eh
        "E_elec_GV": E_elec,      # Eh
        "T_GV":      T_GV,        # K
        "G_corr_GV": G_corr_GV,   # Eh (puede ser None si no hay E_elec)
        "ok":        True,
        "warning":   "",
    }


def run_goodvibes_on_mol(mol_name: str) -> dict[str, dict]:
    """
    Invoca GoodVibes sobre todos los freq.out de 03_numfreq/<conf_*>/freq.out
    para la molécula mol_name.

    Flujo para cada confórmero:
      1. Resuelve el Python correcto (entorno sciencebitch o auto-detección).
      2. Detecta las capacidades de la versión instalada (--output_format json
         vs parseo del .dat de salida clásico).
      3. Crea runs/<mol>/<conf>/goodvibes/ y ejecuta GoodVibes desde allí.
      4. Parsea JSON (v4+) o el .dat tabulado (v3.x) según disponibilidad.

    Retorna dict:
        { "conf_001": {"H_GV": float, "S_GV": float, "G_GV": float,
                       "ZPE_GV": float, "T_GV": float, "ok": bool,
                       "json_path": Path, "warning": str}, ... }

    Si GoodVibes no está disponible o falla para un confórmero, "ok" = False
    y "warning" describe el error — el resto del script no se interrumpe.
    """
    mol_dir  = WORK_ROOT / mol_name
    freq_dir = mol_dir / "03_numfreq"

    if not freq_dir.exists():
        warnings.warn(f"[GoodVibes] {freq_dir} no existe — omitiendo {mol_name}.")
        return {}

    conf_dirs = sorted(freq_dir.glob("conf_*"))
    if not conf_dirs:
        warnings.warn(f"[GoodVibes] Sin subdirectorios conf_* en {freq_dir}.")
        return {}

    # ── Resolver Python y capacidades UNA SOLA VEZ para todos los conf ────
    python_exe = _resolve_gv_python()
    caps       = _probe_gv_capabilities(python_exe)

    # Verificar que goodvibes sea invocable con el ejecutable resuelto
    try:
        probe = subprocess.run(
            [python_exe, "-c", "import goodvibes"],
            capture_output=True, timeout=10
        )
        if probe.returncode != 0:
            warnings.warn(
                f"[GoodVibes] No se puede importar goodvibes con {python_exe}.\n"
                f"  Instalar con: conda activate sciencebitch && pip install goodvibes\n"
                f"  stderr: {probe.stderr.strip()[:200]}"
            )
            return {c.name: {
                "H_GV": None, "S_GV": None, "G_GV": None,
                "ZPE_GV": None, "T_GV": None, "ok": False,
                "json_path": None,
                "warning": f"goodvibes no disponible en {python_exe}",
            } for c in conf_dirs}
    except Exception as exc:
        warnings.warn(f"[GoodVibes] Error al verificar goodvibes: {exc}")
        return {}

    print(f"  [GoodVibes v{caps['version']}] Python: {python_exe}")
    use_json = caps["supports_json_format"]
    print(f"  [GoodVibes] Modo salida: {'JSON (--output_format json)' if use_json else '.dat tabulado (clásico)'}")

    results: dict[str, dict] = {}

    for conf_dir in conf_dirs:
        conf_id  = conf_dir.name
        freq_out = conf_dir / "freq.out"

        null_entry = {
            "H_GV":    None, "S_GV":   None, "G_GV":  None,
            "ZPE_GV":  None, "T_GV":   None,
            "ok":      False, "json_path": None, "warning": "",
        }

        if not freq_out.exists():
            null_entry["warning"] = f"freq.out no encontrado en {conf_dir}"
            results[conf_id] = null_entry
            continue

        # Crear subcarpeta goodvibes/ dentro del directorio del confórmero
        gv_dir = conf_dir / GV_SUBDIR
        gv_dir.mkdir(exist_ok=True)

        # ── Construir comando según capacidades de la versión ─────────────
        # GoodVibes siempre necesita la ruta al .out; corremos desde gv_dir
        # para que los archivos de salida queden ahí con nombres predecibles.
        cmd_base = [
            python_exe, "-m", "goodvibes",
            str(freq_out.resolve()),
            "-q",
            "--qs", GV_QS,
            "-f", str(GV_FREQ_CUTOFF),
            "-t", "298.15",
        ]

        if use_json:
            # v4+: --output_format json → escribe <stem>_goodvibes.json en cwd
            json_path = gv_dir / f"{freq_out.stem}_goodvibes.json"
            cmd = cmd_base + ["--output_format", "json"]
        else:
            # v3.x: sin flag de formato → escribe Goodvibes_output.dat en cwd
            json_path = None
            cmd = cmd_base  # salida .dat clásica

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120,
                cwd=str(gv_dir),
            )
        except subprocess.TimeoutExpired:
            null_entry["warning"] = "GoodVibes timeout (>120 s)"
            results[conf_id] = null_entry
            continue
        except FileNotFoundError:
            null_entry["warning"] = (
                f"Ejecutable no encontrado: {python_exe}\n"
                "  Instalar con: conda activate sciencebitch && pip install goodvibes"
            )
            results[conf_id] = null_entry
            continue

        if proc.returncode != 0:
            null_entry["warning"] = (
                f"GoodVibes retornó código {proc.returncode}. "
                f"stderr: {proc.stderr[:300]}"
            )
            results[conf_id] = null_entry
            continue

        # ── Parsear la salida según el formato producido ───────────────────
        if use_json and json_path is not None:
            if not json_path.exists():
                # v4+ con --output_format json pero JSON no generado →
                # intentar fallback con el .dat que puede haber escrito de todas formas
                dat_candidates = list(gv_dir.glob("*.dat")) + list(gv_dir.glob("Goodvibes_output*"))
                if dat_candidates:
                    null_entry["warning"] = (
                        f"JSON no generado por GoodVibes v4 en {json_path}; "
                        f"usando .dat de fallback."
                    )
                    gv_data = _parse_goodvibes_dat(dat_candidates[0])
                    if gv_data["ok"]:
                        gv_data["json_path"] = dat_candidates[0]
                        results[conf_id] = gv_data
                    else:
                        null_entry["warning"] += f" {gv_data['warning']}"
                        results[conf_id] = null_entry
                else:
                    null_entry["warning"] = (
                        f"JSON no generado por GoodVibes en {json_path}. "
                        f"stdout: {proc.stdout[:300]}"
                    )
                    results[conf_id] = null_entry
                continue

            try:
                gv_data = _parse_goodvibes_json(json_path, freq_out.name)
            except Exception as exc:
                null_entry["warning"] = f"Error parseando JSON de GoodVibes: {exc}"
                results[conf_id] = null_entry
                continue

            gv_data["json_path"] = json_path
            results[conf_id] = gv_data

        else:
            # Parsear .dat tabulado (v3.x o cuando JSON no está disponible)
            dat_candidates = (
                list(gv_dir.glob("Goodvibes_output*"))
                + list(gv_dir.glob("*.dat"))
            )
            if not dat_candidates:
                null_entry["warning"] = (
                    f"GoodVibes no generó ningún archivo de salida en {gv_dir}. "
                    f"stdout: {proc.stdout[:400]}"
                )
                results[conf_id] = null_entry
                continue

            dat_path = dat_candidates[0]
            try:
                gv_data = _parse_goodvibes_dat(dat_path)
            except Exception as exc:
                null_entry["warning"] = f"Error parseando .dat de GoodVibes: {exc}"
                results[conf_id] = null_entry
                continue

            if not gv_data["ok"]:
                null_entry["warning"] = gv_data.get("warning", "Parseo .dat falló")
                results[conf_id] = null_entry
                continue

            gv_data["json_path"] = dat_path  # reaprovechamos el campo para la ruta
            results[conf_id] = gv_data

    return results


def _parse_goodvibes_json(json_path: Path, filename: str) -> dict:
    """
    Parsea el JSON producido por GoodVibes v4.x y extrae H, S, G (mRRHO), ZPE.

    El JSON de GoodVibes v4 tiene la forma:
        {
          "timestamp": "...",
          "version": "4.x",
          "molecules": [
            {
              "filename": "freq.out",
              "t":  298.15,
              "zpe": <float Eh>,
              "h":  <float Eh>,    ← entalpía total (E_elec + H_corr_GV)
              "ts": <float Eh>,    ← T·S (en Eh)
              "g":  <float Eh>,    ← G total (E_elec + G_corr_GV)
              ...
            }
          ]
        }

    Nota: GoodVibes reporta H y G como energías absolutas que incluyen E_elec.
    Para obtener G_corr_GV (compatible con el ciclo compuesto del script):
        G_corr_GV = G_GV - E_elec_GV
    donde E_elec_GV es la energía electrónica que GoodVibes lee del .out
    (campo "sp" o "e_sp" en el JSON, según versión).
    """
    with open(json_path) as f:
        data = json.load(f)

    # GoodVibes v4 agrupa bajo "molecules" (lista) o directamente como dict
    molecules = data.get("molecules") or data.get("data") or []

    # Buscar la entrada que corresponde a nuestro freq.out
    entry = None
    stem  = Path(filename).stem  # "freq"

    for mol in molecules:
        fn = mol.get("filename", mol.get("file", ""))
        if Path(fn).stem == stem or fn == filename:
            entry = mol
            break

    # Fallback: si solo hay una entrada, usarla directamente
    if entry is None and len(molecules) == 1:
        entry = molecules[0]

    if entry is None:
        raise ValueError(
            f"No se encontró entrada para '{filename}' en {json_path}. "
            f"Claves disponibles: {[m.get('filename','?') for m in molecules]}"
        )

    # Extraer campos — GoodVibes v4 usa minúsculas
    def _get(*keys):
        for k in keys:
            v = entry.get(k)
            if v is not None:
                try:
                    return float(v)
                except (TypeError, ValueError):
                    pass
        return None

    T_GV   = _get("t", "temp", "temperature")
    ZPE_GV = _get("zpe", "ZPE")
    H_GV   = _get("h", "H", "enthalpy")           # absoluto (Eh)
    TS_GV  = _get("ts", "tds", "t*s")             # T·S en Eh
    G_GV   = _get("g", "G", "gibbs", "qh_g")      # absoluto (Eh), con mRRHO
    E_elec = _get("sp", "e_sp", "e_elec",
                  "electronic_energy", "scf_energy")

    # Algunas versiones reportan S directamente (en Eh/K), otras T·S
    S_GV = None
    if TS_GV is not None and T_GV is not None and T_GV > 0:
        S_GV = TS_GV / T_GV   # Eh/K
    else:
        S_GV = _get("s", "S", "entropy")

    # G_corr_GV = G_total_GV - E_elec
    # Si GoodVibes no reporta E_elec en el JSON, lo dejamos None;
    # el apéndice usará G_GV directamente para ΔG_rel (diferencias).
    G_corr_GV = None
    if G_GV is not None and E_elec is not None:
        G_corr_GV = G_GV - E_elec

    ok = G_GV is not None

    return {
        "H_GV":      H_GV,
        "S_GV":      S_GV,       # Eh/K
        "TS_GV":     TS_GV,      # Eh
        "G_GV":      G_GV,       # Eh absoluto (E_elec + G_corr_GV)
        "G_corr_GV": G_corr_GV,  # Eh — corrección pura (None si no hay E_elec)
        "ZPE_GV":    ZPE_GV,     # Eh
        "E_elec_GV": E_elec,     # Eh — energía electrónica leída por GoodVibes
        "T_GV":      T_GV,       # K
        "ok":        ok,
        "warning":   "" if ok else "G_GV no encontrado en JSON",
    }


def compute_gv_appendix(
    mol_name: str,
    df_conf: pd.DataFrame,
    gv_results: dict[str, dict],
) -> pd.DataFrame:
    """
    Construye el DataFrame del apéndice GoodVibes, comparando G_ORCA vs G_GV
    para cada confórmero a todas las temperaturas.

    Requiere que df_conf tenga las columnas producidas por compute_G_total()
    (G_total_{label}, dG_{label}_kcal, boltz_{label}, H_corr, S_total_EhK).
    Si compute_G_total() no fue llamado previamente (el flujo normal de
    collect_results no lo llama), este módulo lo invoca internamente.

    Columnas del DataFrame resultante:
        conformer          — id del confórmero
        E_SP_hartree       — energía SP ORCA (def2-TZVPD, Eh)
        — ORCA (RRHO) —
        G_corr_ORCA_Eh     — corrección G de ORCA a 298 K (Eh)
        H_corr_ORCA_Eh     — corrección H de ORCA (Eh)
        S_ORCA_EhK         — entropía ORCA (Eh/K)
        G_total_ORCA_Eh    — G compuesto ORCA a 298 K (Eh)
        dG_ORCA_kcal       — ΔG relativo ORCA a 298 K (kcal/mol)
        boltz_ORCA_298K    — peso Boltzmann ORCA a 298 K
        — GoodVibes (mRRHO Grimme) —
        G_GV_abs_Eh        — G absoluto GoodVibes (Eh, incluye E_elec)
        G_corr_GV_Eh       — corrección G de GV (Eh, None si sin E_elec en JSON)
        S_GV_EhK           — entropía quasi-RRHO GV (Eh/K)
        ZPE_GV_Eh          — ZPE de GoodVibes (Eh)
        G_total_GV_Eh      — G compuesto GV a 298 K = E_SP + G_corr_GV (Eh)
        dG_GV_kcal         — ΔG relativo GV a 298 K (kcal/mol)
        boltz_GV_298K      — peso Boltzmann GV a 298 K
        — Diferencias —
        dG_corr_Eh         — G_corr_GV - G_corr_ORCA (Eh)
        dG_corr_kcal       — ídem en kcal/mol
        dS_EhK             — S_GV - S_ORCA (Eh/K)  ← origen del Δ
        ddG_rel_kcal       — ΔG_GV_rel - ΔG_ORCA_rel a 298 K (kcal/mol)
        d_boltz_298K       — boltz_GV - boltz_ORCA a 298 K
        rank_ORCA          — ranking conformacional ORCA (1=más estable)
        rank_GV            — ranking conformacional GV
        rank_change_GV     — rank_GV - rank_ORCA (0 = sin inversión)
        gv_ok              — True si GoodVibes procesó este confórmero
        gv_warning         — mensaje de advertencia si gv_ok=False
    """
    mol_dir = WORK_ROOT / mol_name

    # ── Asegurar correcciones térmicas ORCA en df_conf ────────────────────
    # Si compute_G_total ya fue llamado, G_total_298K existirá.
    # Si no (flujo normal), lo recalculamos internamente aquí.
    if "G_total_298K" not in df_conf.columns or "H_corr" not in df_conf.columns:
        # Cargar termoquímica ORCA para cada confórmero
        df_work = df_conf.copy()
        _thermochem_cache: list[Optional[dict]] = []
        for _, row in df_work.iterrows():
            conf_id  = row["conformer"]
            freq_out = mol_dir / "03_numfreq" / conf_id / "freq.out"
            tc = parse_thermochem(freq_out, T=T_STD)
            _thermochem_cache.append(tc if tc.get("ok") else None)
        df_work["_thermochem"] = _thermochem_cache
        df_work = compute_G_total(df_work)
        df_work = compute_boltzmann_multiT(df_work)
    else:
        df_work = df_conf.copy()

    rows: list[dict] = []

    for _, row in df_work.iterrows():
        conf_id  = row["conformer"]
        gv       = gv_results.get(conf_id, {})
        gv_ok    = gv.get("ok", False)

        e_sp = row.get("E_hartree")

        # ── ORCA RRHO ──────────────────────────────────────────────────────
        g_corr_orca = row.get("G_corr_298K")     # Eh
        h_corr_orca = row.get("H_corr")           # Eh
        s_orca      = row.get("S_total_EhK")      # Eh/K
        g_total_orca = row.get("G_total_298K")    # Eh
        dg_orca      = row.get("dG_298K_kcal")    # kcal/mol
        boltz_orca   = row.get("boltz_298K")

        # ── GoodVibes mRRHO ────────────────────────────────────────────────
        g_gv_abs   = gv.get("G_GV")              # Eh absoluto
        g_corr_gv  = gv.get("G_corr_GV")         # Eh (puede ser None)
        s_gv       = gv.get("S_GV")              # Eh/K
        zpe_gv     = gv.get("ZPE_GV")            # Eh
        e_elec_gv  = gv.get("E_elec_GV")         # Eh

        # G compuesto GV = E_SP(TZVPD) + G_corr_GV
        # Si GoodVibes no reporta E_elec separado, usamos G_GV directamente
        # como proxy (válido para diferencias relativas entre confórmeros).
        g_total_gv: Optional[float] = None
        if g_corr_gv is not None and e_sp is not None:
            g_total_gv = e_sp + g_corr_gv          # ciclo compuesto correcto
        elif g_gv_abs is not None and e_elec_gv is not None and e_sp is not None:
            # Reconstruir G_corr a partir del absoluto
            g_corr_gv  = g_gv_abs - e_elec_gv
            g_total_gv = e_sp + g_corr_gv
        elif g_gv_abs is not None:
            # Último recurso: usar G absoluto de GV para los ΔG (pérdida del SP TZVPD)
            g_total_gv = g_gv_abs

        # ── Diferencias ────────────────────────────────────────────────────
        dg_corr_eh   = None
        dg_corr_kcal = None
        ds_ehk       = None

        if g_corr_gv is not None and g_corr_orca is not None:
            dg_corr_eh   = g_corr_gv - g_corr_orca
            dg_corr_kcal = dg_corr_eh * HARTREE_TO_KCAL
        if s_gv is not None and s_orca is not None:
            ds_ehk = s_gv - s_orca

        rows.append({
            "conformer":        conf_id,
            "E_SP_hartree":     e_sp,
            # ORCA
            "G_corr_ORCA_Eh":   g_corr_orca,
            "H_corr_ORCA_Eh":   h_corr_orca,
            "S_ORCA_EhK":       s_orca,
            "G_total_ORCA_Eh":  g_total_orca,
            "dG_ORCA_kcal":     dg_orca,
            "boltz_ORCA_298K":  boltz_orca,
            # GoodVibes
            "G_GV_abs_Eh":      g_gv_abs,
            "G_corr_GV_Eh":     g_corr_gv,
            "S_GV_EhK":         s_gv,
            "ZPE_GV_Eh":        zpe_gv,
            "G_total_GV_Eh":    g_total_gv,
            # (dG_GV_kcal y boltz_GV_298K se calculan abajo tras tener todos los g_total)
            "gv_ok":            gv_ok,
            "gv_warning":       gv.get("warning", ""),
            # Diferencias parciales (las que no dependen del mínimo del ensemble)
            "dG_corr_Eh":       dg_corr_eh,
            "dG_corr_kcal":     dg_corr_kcal,
            "dS_EhK":           ds_ehk,
        })

    df_gv = pd.DataFrame(rows)

    if df_gv.empty:
        return df_gv

    # ── ΔG relativo y Boltzmann para GV (necesita el mínimo del ensemble) ──
    valid_gv = df_gv["G_total_GV_Eh"].notna() & df_gv["gv_ok"]

    if valid_gv.any():
        g_min_gv = df_gv.loc[valid_gv, "G_total_GV_Eh"].min()
        df_gv["dG_GV_kcal"] = np.nan
        df_gv.loc[valid_gv, "dG_GV_kcal"] = (
            (df_gv.loc[valid_gv, "G_total_GV_Eh"] - g_min_gv) * HARTREE_TO_KCAL
        )

        # Boltzmann GV a 298 K
        dg_gv = df_gv.loc[valid_gv, "dG_GV_kcal"]
        raw   = np.exp(-dg_gv / (KB_KCAL * T_STD))
        df_gv["boltz_GV_298K"] = np.nan
        df_gv.loc[valid_gv, "boltz_GV_298K"] = raw / raw.sum()
    else:
        df_gv["dG_GV_kcal"]    = np.nan
        df_gv["boltz_GV_298K"] = np.nan

    # ── ΔΔG relativo y Δboltz ─────────────────────────────────────────────
    df_gv["ddG_rel_kcal"] = df_gv["dG_GV_kcal"] - df_gv["dG_ORCA_kcal"]
    df_gv["d_boltz_298K"] = df_gv["boltz_GV_298K"] - df_gv["boltz_ORCA_298K"]

    # ── Rankings ──────────────────────────────────────────────────────────
    valid_orca = df_gv["G_total_ORCA_Eh"].notna()
    if valid_orca.any():
        df_gv.loc[valid_orca, "rank_ORCA"] = (
            df_gv.loc[valid_orca, "G_total_ORCA_Eh"]
            .rank(method="min").astype("Int64")
        )
    if valid_gv.any():
        df_gv.loc[valid_gv, "rank_GV"] = (
            df_gv.loc[valid_gv, "G_total_GV_Eh"]
            .rank(method="min").astype("Int64")
        )
    if "rank_ORCA" in df_gv.columns and "rank_GV" in df_gv.columns:
        df_gv["rank_change_GV"] = df_gv["rank_GV"] - df_gv["rank_ORCA"]

    return df_gv


def print_gv_appendix(df_gv: pd.DataFrame, mol_name: str):
    """
    Imprime el apéndice GoodVibes con el mismo estilo que el resto del script.
    Muestra tres secciones:
      1. Correcciones absolutas ORCA vs GV (G_corr, S, ZPE)
      2. Impacto en ΔG relativo y pesos Boltzmann a 298 K
      3. Cambios de ranking conformacional (si los hay)
    """
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  APÉNDICE GOODVIBES — Corrección cuasi-armónica mRRHO (Grimme)")
    print(f"  Molécula : {mol_name}")
    print(f"  Método   : -q --qs grimme -f {GV_FREQ_CUTOFF:.0f} cm⁻¹ (v{_gv_version()})")
    print(f"  Ciclo    : G_GV = E_SP(TZVPD/CPCM) + G_corr_GV(mRRHO, SVP/CPCM)")
    print(sep)

    if df_gv.empty:
        print("  Sin datos de GoodVibes.")
        return

    # Confórmeros sin GV
    failed = df_gv[~df_gv["gv_ok"]]
    if not failed.empty:
        print(f"\n  ⚠  {len(failed)} confórmero(s) sin datos GoodVibes:")
        for _, r in failed.iterrows():
            print(f"     {r['conformer']:>10}  →  {r['gv_warning']}")

    ok = df_gv[df_gv["gv_ok"]].copy()
    if ok.empty:
        print("  Sin confórmeros válidos de GoodVibes para comparar.")
        return

    # ── Sección 1: Correcciones absolutas ────────────────────────────────
    print(f"\n  ── 1. Correcciones absolutas (Eh) ── "
          f"[ORCA=RRHO  |  GV=mRRHO Grimme, f₀={GV_FREQ_CUTOFF:.0f} cm⁻¹]")
    print(f"\n  {'conf':>10}  {'G_corr_ORCA':>14}  {'G_corr_GV':>12}  "
          f"{'ΔG_corr(Eh)':>13}  {'ΔG_corr(kcal)':>14}  "
          f"{'S_ORCA(Eh/K)':>13}  {'S_GV(Eh/K)':>11}  {'ΔS(Eh/K)':>10}")
    print("  " + "─" * 100)

    for _, r in ok.iterrows():
        def _f(v, fmt=".6f"):
            return f"{v:{fmt}}" if pd.notna(v) else "      —     "

        dg_sign = "+" if (r.get("dG_corr_Eh") or 0) >= 0 else ""
        ds_sign = "+" if (r.get("dS_EhK") or 0) >= 0 else ""
        print(
            f"  {r['conformer']:>10}  "
            f"{_f(r.get('G_corr_ORCA_Eh')):>14}  "
            f"{_f(r.get('G_corr_GV_Eh')):>12}  "
            f"{dg_sign}{_f(r.get('dG_corr_Eh'), '.6f'):>12}  "
            f"{dg_sign}{_f(r.get('dG_corr_kcal'), '.4f'):>13}  "
            f"{_f(r.get('S_ORCA_EhK')):>13}  "
            f"{_f(r.get('S_GV_EhK')):>11}  "
            f"{ds_sign}{_f(r.get('dS_EhK'), '.6f'):>9}"
        )

    # Estadística del Δ(G_corr)
    dg_vals = ok["dG_corr_kcal"].dropna()
    if not dg_vals.empty:
        print(f"\n  Δ(G_corr) medio  : {dg_vals.mean():+.4f} kcal/mol")
        print(f"  Δ(G_corr) rango  : [{dg_vals.min():+.4f}, {dg_vals.max():+.4f}] kcal/mol")
        print(f"  Interpretación   : GV{'aumenta' if dg_vals.mean()>0 else 'reduce'} G_corr "
              f"{'(mRRHO menos negativo → G más alto)' if dg_vals.mean()>0 else '(mRRHO más negativo → G más bajo)'}")

    # ── Sección 2: Impacto en ΔG relativo y Boltzmann ────────────────────
    print(f"\n  ── 2. Impacto en ΔG relativo y pesos Boltzmann (T = 298.15 K) ──")
    print(f"\n  {'conf':>10}  {'ΔG_ORCA':>9}  {'ΔG_GV':>8}  {'ΔΔG_rel':>9}  "
          f"{'boltz_ORCA':>11}  {'boltz_GV':>10}  {'Δboltz':>8}  {'dominante':>9}")
    print("  " + "─" * 85)

    for _, r in ok.sort_values("dG_GV_kcal").iterrows():
        dg_o   = r.get("dG_ORCA_kcal")
        dg_gv  = r.get("dG_GV_kcal")
        ddg    = r.get("ddG_rel_kcal")
        bo     = r.get("boltz_ORCA_298K")
        bg     = r.get("boltz_GV_298K")
        db     = r.get("d_boltz_298K")

        def _fv(v, fmt):
            return f"{v:{fmt}}" if pd.notna(v) else "   —   "

        dom_gv = " ★" if pd.notna(bg) and bg > 0.50 else ("  ·" if pd.notna(bg) and bg > 0.10 else "   ")
        print(
            f"  {r['conformer']:>10}  "
            f"{_fv(dg_o,  '+9.3f')}  "
            f"{_fv(dg_gv, '+8.3f')}  "
            f"{_fv(ddg,   '+9.4f')}  "
            f"{_fv(bo,    '11.4f')}  "
            f"{_fv(bg,    '10.4f')}  "
            f"{_fv(db,    '+8.4f')}"
            f"{dom_gv}"
        )

    # Resumen del impacto sobre Boltzmann
    ddg_vals = ok["ddG_rel_kcal"].dropna()
    db_vals  = ok["d_boltz_298K"].dropna()
    if not ddg_vals.empty:
        print(f"\n  |ΔΔG_rel| máximo : {ddg_vals.abs().max():.4f} kcal/mol")
        if ddg_vals.abs().max() < 0.1:
            print("  → Impacto negligible en los ΔG relativos (< 0.1 kcal/mol)")
        elif ddg_vals.abs().max() < 0.5:
            print("  → Impacto pequeño pero no nulo en los ΔG relativos (0.1–0.5 kcal/mol)")
        else:
            print("  → Impacto significativo en los ΔG relativos (> 0.5 kcal/mol) ⚠")
    if not db_vals.empty:
        print(f"  |Δboltz| máximo  : {db_vals.abs().max():.4f} "
              f"({'cambio relevante en ensemble' if db_vals.abs().max() > 0.05 else 'cambio menor en ensemble'})")

    # ── Sección 3: Rankings ───────────────────────────────────────────────
    if "rank_change_GV" in df_gv.columns:
        inversions = df_gv["rank_change_GV"].dropna()
        n_inv = (inversions.abs() > 0).sum()
        print(f"\n  ── 3. Ranking conformacional: ORCA vs GoodVibes ──")
        print(f"  Inversiones detectadas: {n_inv}/{len(ok)}")

        if n_inv > 0:
            print(f"\n  {'conf':>10}  {'rank_ORCA':>10}  {'rank_GV':>8}  {'Δrank':>7}")
            print("  " + "─" * 42)
            for _, r in df_gv[df_gv["rank_change_GV"].abs() > 0].iterrows():
                rch = r.get("rank_change_GV", 0)
                flag = " ◀ inversión"
                print(
                    f"  {r['conformer']:>10}  "
                    f"{str(r.get('rank_ORCA','—')):>10}  "
                    f"{str(r.get('rank_GV','—')):>8}  "
                    f"{int(rch):>+7}{flag}"
                )
        else:
            print("  → El ranking conformacional es idéntico con ORCA-RRHO y GV-mRRHO.")

    print(f"\n  Archivos GV guardados en:")
    for conf_id in df_gv["conformer"]:
        p = WORK_ROOT / mol_name / "03_numfreq" / conf_id / GV_SUBDIR
        if p.exists():
            # El nombre del JSON depende del stem del .out (normalmente freq_goodvibes.json)
            json_files = list(p.glob("*_goodvibes.json"))
            fname = json_files[0].name if json_files else "*.json"
            print(f"    {p}/{fname}")

    print(f"\n{sep}\n")


def _gv_version() -> str:
    """Devuelve la versión de GoodVibes instalada, o '?' si no se puede determinar."""
    try:
        python_exe = _resolve_gv_python()
        result = subprocess.run(
            [python_exe, "-m", "goodvibes", "--version"],
            capture_output=True, text=True, timeout=10
        )
        for line in (result.stdout + result.stderr).splitlines():
            m = re.search(r"([\d]+\.[\d]+\.[\d]+|v[\d.]+)", line)
            if m:
                return m.group(1).lstrip("v")
    except Exception:
        pass
    return "?"


def export_gv_appendix(df_gv: pd.DataFrame, mol_name: str, fmt: str):
    """
    Guarda el DataFrame del apéndice GoodVibes en runs/<mol>/.
    Formatos: "csv", "xlsx" o "both".
    Si fmt incluye xlsx, lo agrega como hoja extra del Excel existente si ya
    existe, o crea un archivo dedicado goodvibes_appendix.xlsx.
    """
    if df_gv.empty:
        return

    mol_dir   = WORK_ROOT / mol_name
    base_name = f"{mol_name}_goodvibes_appendix"

    # Columnas a exportar (en orden lógico)
    export_cols = [
        "conformer", "E_SP_hartree",
        "G_corr_ORCA_Eh", "G_corr_GV_Eh", "dG_corr_Eh", "dG_corr_kcal",
        "S_ORCA_EhK", "S_GV_EhK", "dS_EhK", "ZPE_GV_Eh",
        "G_total_ORCA_Eh", "G_total_GV_Eh",
        "dG_ORCA_kcal", "dG_GV_kcal", "ddG_rel_kcal",
        "boltz_ORCA_298K", "boltz_GV_298K", "d_boltz_298K",
        "rank_ORCA", "rank_GV", "rank_change_GV",
        "gv_ok", "gv_warning",
    ]
    export_cols = [c for c in export_cols if c in df_gv.columns]
    df_out = df_gv[export_cols].copy()

    if fmt in ("csv", "both"):
        csv_path = mol_dir / f"{base_name}.csv"
        df_out.to_csv(csv_path, index=False)
        print(f"  GoodVibes appendix CSV → {csv_path}")

    if fmt in ("xlsx", "both"):
        # Intentar agregar hoja al Excel existente
        xlsx_main = mol_dir / f"{mol_name}_results.xlsx"
        xlsx_gv   = mol_dir / f"{base_name}.xlsx"

        target = xlsx_main if xlsx_main.exists() else xlsx_gv
        mode   = "a" if xlsx_main.exists() else "w"
        kw     = {"if_sheet_exists": "replace"} if mode == "a" else {}

        try:
            with pd.ExcelWriter(target, engine="openpyxl", mode=mode, **kw) as writer:
                df_out.to_excel(writer, sheet_name="GoodVibes_mRRHO", index=False)
            print(f"  GoodVibes appendix XLSX → {target} (hoja: GoodVibes_mRRHO)")
        except Exception as exc:
            # Si falla al agregar al Excel existente, escribir archivo separado
            warnings.warn(f"[GoodVibes export] No se pudo agregar a {target}: {exc}. "
                          f"Escribiendo {xlsx_gv}.")
            with pd.ExcelWriter(xlsx_gv, engine="openpyxl") as writer:
                df_out.to_excel(writer, sheet_name="GoodVibes_mRRHO", index=False)
            print(f"  GoodVibes appendix XLSX → {xlsx_gv}")


# ═══════════════════════════════════════════════════════════════════════════════
#  MÓDULO ADC(2) — extracción de λ_max y comparación con TD-DFT
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Fundamento:
#  ADC(2) (Algebraic Diagrammatic Construction, segundo orden) es un método
#  ab initio para estados excitados más preciso que TD-DFT pero también más
#  costoso. La comparación sistemática entre λ_max(TD-DFT) y λ_max(ADC2)
#  permite evaluar si el error del funcional es sistemático (desplazamiento
#  constante) y si puede corregirse empíricamente aplicando un Δλ = λ_TDDFT - λ_ADC2
#  uniforme al espectro UV-Vis calculado con TD-DFT.
#
#  Estructura de archivos esperada:
#    runs/<mol>/07_adc2/conf_<id>/adc2.out
#
#  Por cada conformero se parsea UN root (nroots 1), extrayendo:
#    - E_eV     : energía de excitación vertical (eV)
#    - wl_nm    : longitud de onda (nm)
#    - fosc     : fuerza de oscilador (del bloque LEFT-RIGHT, el más equilibrado)
#
#  Luego se calcula λ_max ADC2 Boltzmann-promediada al igual que en TDDFT y
#  se compara estado a estado con λ_max TDDFT.
#
# ═══════════════════════════════════════════════════════════════════════════════

ADC2_SUBDIR = "07_adc2"   # subcarpeta dentro de runs/<mol>/


def parse_adc2_excitations(adc2_out: Path) -> list[dict]:
    """
    Extrae la(s) excitación(es) del bloque ADC(2) de un output ORCA 6.

    Busca el bloque "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    del espectro LEFT-RIGHT (el más equilibrado entre los tres que imprime ORCA
    para ADC2: RIGHT, LEFT y LEFT-RIGHT). Si no encuentra el bloque LEFT-RIGHT,
    usa el bloque RIGHT como fallback.

    Formato de la línea de datos (ORCA 6):
        0-1A  ->  1-1A    3.464792   27945.4   357.8   0.416153266   ...

    Retorna lista de dicts con:
        state      : int  — número de root (1, 2, …)
        E_eV       : float — energía de excitación (eV)
        E_cm1      : float — energía en cm⁻¹
        wl_nm      : float — longitud de onda (nm)
        fosc       : float — fuerza de oscilador (LEFT-RIGHT si disponible)
        fosc_right : float — fuerza de oscilador RIGHT (siempre disponible)
        fosc_left  : float — fuerza de oscilador LEFT (si disponible)
        source     : str  — qué bloque se usó como fosc principal
    """
    if not adc2_out.exists():
        return []

    try:
        text = adc2_out.read_text(errors="replace")
    except OSError:
        return []

    lines = text.splitlines()

    # Patrón para líneas de excitación del bloque de absorción:
    # "  0-1A  ->  1-1A    3.464792   27945.4   357.8   0.416153266  ..."
    pat_line = re.compile(
        r"^\s*(\d+)-\S+\s*->\s*(\d+)-\S+\s+"   # initial_root -> final_root
        r"([\d.]+)\s+"                           # E_eV
        r"([\d.]+)\s+"                           # E_cm-1
        r"([\d.]+)\s+"                           # wl_nm
        r"([\d.]+)"                              # fosc
    )

    # Buscar los tres bloques de absorción que imprime ADC2.
    # ORCA 6 los imprime en este orden: RIGHT, LEFT, LEFT-RIGHT.
    # Cada bloque está precedido por su título y la cabecera de columnas.
    block_tags = {
        "right":      "SPECTRUM FOR RIGHT TRANSITION MOMENTS",
        "left":       "SPECTRUM FOR LEFT TRANSITION MOMENTS",
        "left_right": "SPECTRUM FOR LEFT-RIGHT TRANSITION MOMENTS",
    }
    # La línea de datos de absorción viene tras el separador "---" y la cabecera;
    # marcamos el inicio del bloque con "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC"
    abs_header = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"

    # Recolectar excitaciones por tipo de bloque
    collected: dict[str, list[dict]] = {"right": [], "left": [], "left_right": []}

    current_block: Optional[str] = None
    in_abs_block   = False

    for line in lines:
        # Detectar a qué bloque de momentos de transición entramos
        for key, tag in block_tags.items():
            if tag in line:
                current_block = key
                in_abs_block  = False
                break

        # Dentro del bloque correcto, esperar el encabezado de absorción
        if current_block and not in_abs_block and abs_header in line:
            in_abs_block = True
            continue

        if in_abs_block and current_block:
            m = pat_line.match(line)
            if m:
                root_initial = int(m.group(1))
                root_final   = int(m.group(2))
                # El estado excitado es root_final (convenio ORCA: "0->1" = S1)
                state_n = root_final if root_final > 0 else 1
                collected[current_block].append({
                    "state":  state_n,
                    "E_eV":   float(m.group(3)),
                    "E_cm1":  float(m.group(4)),
                    "wl_nm":  float(m.group(5)),
                    "fosc":   float(m.group(6)),
                })
            # Salir del bloque si encontramos una línea vacía seguida de "---"
            # o el comienzo de otro bloque (se detectará en el bucle externo)

    # Construir lista final cruzando los tres bloques
    # Usamos LEFT-RIGHT como fosc principal; si no existe, RIGHT.
    primary   = collected["left_right"] or collected["right"]
    right_map = {e["state"]: e["fosc"] for e in collected["right"]}
    left_map  = {e["state"]: e["fosc"] for e in collected["left"]}

    result = []
    for exc in primary:
        s = exc["state"]
        result.append({
            "state":      s,
            "E_eV":       exc["E_eV"],
            "E_cm1":      exc["E_cm1"],
            "wl_nm":      exc["wl_nm"],
            "fosc":       exc["fosc"],
            "fosc_right": right_map.get(s, float("nan")),
            "fosc_left":  left_map.get(s, float("nan")),
            "source":     "left_right" if collected["left_right"] else "right",
        })
    return result


def collect_adc2_results(mol_name: str) -> pd.DataFrame:
    """
    Recopila todas las excitaciones ADC(2) de los confórmeros de una molécula.

    Estructura esperada:
        runs/<mol>/07_adc2/conf_<id>/adc2.out

    Retorna DataFrame con columnas:
        conformer, state, E_eV, E_cm1, wl_nm, fosc, fosc_right, fosc_left, source
    O DataFrame vacío si no hay archivos.
    """
    mol_dir   = WORK_ROOT / mol_name
    adc2_dir  = mol_dir / ADC2_SUBDIR
    if not adc2_dir.exists():
        return pd.DataFrame()

    rows = []
    for conf_dir in sorted(adc2_dir.glob("conf_*")):
        if not conf_dir.is_dir():
            continue
        conf_id  = conf_dir.name
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
    w  = s1[col_b].fillna(0)
    if w.sum() <= 0:
        return None

    return float((s1["wl_nm"] * w).sum() / w.sum())


# ─────────────────────────────────────────────
#  COMPARACIÓN ADC2 vs TDDFT (ΔTDDFT-ADC2)
# ─────────────────────────────────────────────

def compare_adc2_tddft(
    mols: Optional[list[str]] = None,
    temp_label: str = "298K",
) -> pd.DataFrame:
    """
    Compara λ_max(TD-DFT) vs λ_max(ADC2) entre estados de protonación.

    Para cada molécula disponible en runs/ que tenga datos TD-DFT Y ADC2,
    calcula:
        lambda_tddft_nm   : λ_max Boltzmann-promediado TD-DFT (nm)
        lambda_adc2_nm    : λ_max Boltzmann-promediado ADC2 (nm)
        delta_nm          : λ_TDDFT - λ_ADC2 (nm)
                            > 0 → TDDFT sobreestima λ (redshift)
                            < 0 → TDDFT subestima λ (blueshift)
        delta_eV          : conversión a eV  (Δ positivo = TDDFT demasiado bajo en energía)
        n_confs_adc2      : número de confórmeros con ADC2 completado
        n_confs_tddft     : número de confórmeros con TDDFT disponible

    También imprime un análisis de sistematicidad:
        - ¿Es Δ constante entre estados? (std_delta → corrección sistemática viable)
        - Corrección empírica recomendada: shift = mean(Δ) ± std(Δ)

    Parámetros
    ----------
    mols        : lista de nombres de mol (directorios en runs/). None → todos.
    temp_label  : temperatura para el Boltzmann (default "298K").

    Retorna DataFrame con una fila por molécula.
    """
    if mols is None:
        mols = sorted(p.name for p in WORK_ROOT.iterdir()
                      if p.is_dir() and (p / "04_sp").exists())

    rows = []
    for mol in mols:
        # ── Cargar TDDFT ───────────────────────────────────────────────────
        try:
            df_conf, df_tddft = collect_results(mol)
        except Exception as e:
            warnings.warn(f"[ADC2 vs TDDFT] {mol}: error al cargar TDDFT — {e}")
            df_conf  = pd.DataFrame()
            df_tddft = pd.DataFrame()

        lmax_tddft: Optional[float] = None
        n_confs_tddft = 0
        if not df_tddft.empty and not df_conf.empty:
            lmax_tddft    = df_conf.attrs.get("lambda_max_boltz_298K") \
                            if temp_label == "298K" else \
                            _boltz_weighted_lmax(df_tddft, df_conf, temp_label)
            n_confs_tddft = int(df_tddft["conformer"].nunique())

        # ── Cargar ADC2 ────────────────────────────────────────────────────
        df_adc2      = collect_adc2_results(mol)
        lmax_adc2: Optional[float] = None
        n_confs_adc2 = 0
        if not df_adc2.empty and not df_conf.empty:
            lmax_adc2    = compute_adc2_lambda_max(df_adc2, df_conf, temp_label)
            n_confs_adc2 = int(df_adc2["conformer"].nunique())

        # ── Δ = TDDFT − ADC2 ──────────────────────────────────────────────
        delta_nm: Optional[float] = None
        delta_eV: Optional[float] = None
        if lmax_tddft is not None and lmax_adc2 is not None and lmax_adc2 > 0:
            delta_nm = lmax_tddft - lmax_adc2
            # Conversión eV: E = hc/λ → ΔE ≈ -hc·Δλ/λ² (diferencia de energías)
            # Calculamos directamente como diferencia de energías en eV:
            ev_tddft = 1239.84193 / lmax_tddft   # eV
            ev_adc2  = 1239.84193 / lmax_adc2    # eV
            delta_eV = ev_tddft - ev_adc2         # >0 → TDDFT subestima E (redshift)

        rows.append({
            "molecule":        mol,
            "lambda_tddft_nm": lmax_tddft,
            "lambda_adc2_nm":  lmax_adc2,
            "delta_nm":        delta_nm,        # TDDFT − ADC2 en nm
            "delta_eV":        delta_eV,        # TDDFT − ADC2 en eV
            "n_confs_tddft":   n_confs_tddft,
            "n_confs_adc2":    n_confs_adc2,
        })

    return pd.DataFrame(rows)


def _boltz_weighted_lmax(
    df_tddft: pd.DataFrame,
    df_conf:  pd.DataFrame,
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
    w  = s1[col_b].fillna(0)
    return float((s1["wavelength_nm"] * w).sum() / w.sum()) if w.sum() > 0 else None


def print_adc2_summary(df_adc2: pd.DataFrame, mol_name: str,
                       df_conf: pd.DataFrame):
    """Imprime tabla de excitaciones ADC2 de una molécula."""
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  ADC(2) EXCITACIONES — {mol_name}")
    print(sep)

    if df_adc2.empty:
        print(f"  Sin datos ADC2 para '{mol_name}'.")
        print(f"  Ruta esperada: runs/{mol_name}/{ADC2_SUBDIR}/conf_*/adc2.out")
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
        adc2  = row.get("lambda_adc2_nm")
        delta = row.get("delta_nm")
        devv  = row.get("delta_eV")
        n     = row.get("n_confs_adc2", 0)

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
    std_nm  = float(np.std(deltas_nm, ddof=1)) if len(deltas_nm) > 1 else 0.0
    mean_eV = float(np.mean(deltas_eV))
    std_eV  = float(np.std(deltas_eV, ddof=1)) if len(deltas_eV) > 1 else 0.0

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
        mols_valid     = [r["molecule"]        for r in valid_rows]
        print(f"\n  ── Comparación λ_max ADC2 entre estados de protonación ──")
        adc2_min  = min(lmax_adc2_vals)
        adc2_max  = max(lmax_adc2_vals)
        mol_min   = mols_valid[lmax_adc2_vals.index(adc2_min)]
        mol_max   = mols_valid[lmax_adc2_vals.index(adc2_max)]
        print(f"  Rango ADC2: {adc2_min:.1f} – {adc2_max:.1f} nm")
        print(f"  Estado con λ_max ADC2 menor: {mol_min} ({adc2_min:.1f} nm)")
        print(f"  Estado con λ_max ADC2 mayor: {mol_max} ({adc2_max:.1f} nm)")

    print(f"\n{sep}\n")


def export_adc2_comparison(df_cmp: pd.DataFrame, fmt: str):
    """Exporta la tabla comparativa ADC2 vs TDDFT a CSV y/o XLSX."""
    if df_cmp.empty:
        return
    out = WORK_ROOT / "comparison_adc2_tddft"
    if fmt in ("csv", "both"):
        df_cmp.to_csv(out.with_suffix(".csv"), index=False)
        print(f"  CSV ADC2 vs TDDFT → {out.with_suffix('.csv')}")
    if fmt in ("xlsx", "both"):
        df_cmp.to_excel(out.with_suffix(".xlsx"), index=False)
        print(f"  XLSX ADC2 vs TDDFT → {out.with_suffix('.xlsx')}")


def plot_adc2_tddft_comparison(df_cmp: pd.DataFrame,
                                dfs_conf: dict, dfs_tddft: dict,
                                sigma: float = 20.0,
                                x_range: tuple = (200, 500)):
    """
    Grafica la comparación ADC2 vs TDDFT mostrando para cada estado:
      - Espectro TDDFT (gaussiano, normalizado)
      - Línea vertical λ_max TDDFT
      - Línea vertical λ_max ADC2
      - Flechas mostrando el desplazamiento Δ

    Panel inferior: barra de Δ(nm) por estado de protonación.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("matplotlib no disponible — saltando plot ADC2")
        return

    if df_cmp.empty:
        return

    # Filtrar solo filas con datos completos
    valid = df_cmp.dropna(subset=["lambda_tddft_nm", "lambda_adc2_nm"])
    if valid.empty:
        print("  Sin pares TDDFT+ADC2 para graficar.")
        return

    n_mols = len(valid)
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
    out_png = WORK_ROOT / "comparison_adc2_tddft.png"
    plt.savefig(out_png, dpi=150)
    print(f"  Plot ADC2 vs TDDFT guardado en: {out_png}")
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
  python results_collector.py --mol mglyval --goodvibes
  python results_collector.py --mol mglyval --goodvibes --export csv
  python results_collector.py --compare
  python results_collector.py --compare --plot --ref mglyval_neutral
        """,
    )

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--mol",     metavar="NAME",
                      help="analizar un estado (nombre del directorio en runs/)")
    mode.add_argument("--compare", action="store_true",
                      help="comparar todos los estados de protonacion en runs/")

    parser.add_argument("--ref",    metavar="NAME", default=None,
                        help="estado de referencia para dE relativa (solo con --compare)")
    parser.add_argument("--export", choices=["csv", "xlsx", "both"], default=None,
                        help="exportar tablas")
    parser.add_argument("--plot",   action="store_true",
                        help="generar espectro(s) UV-Vis")
    parser.add_argument("--goodvibes", action="store_true",
                        help="correr GoodVibes (mRRHO Grimme) sobre cada "
                             "freq.out y comparar contra el G_corr de ORCA "
                             "(solo con --mol). Imprime y exporta el apendice "
                             "de comparacion a runs/<mol>/.")
    parser.add_argument("--adc2",   action="store_true",
                        help="extraer lambda_max ADC(2) y comparar con TD-DFT. "
                             "Con --mol: analiza ese estado. "
                             "Con --compare: compara todos los estados y evalua "
                             "la sistematicidad del desplazamiento TDDFT-ADC2.")
    parser.add_argument("--sigma",  type=float, default=20.0,
                        help="anchura Gaussiana en nm para los plots (default: 20)")
    parser.add_argument("--xmin",   type=float, default=200.0,
                        help="limite inferior del eje lambda en nm (default: 200)")
    parser.add_argument("--xmax",   type=float, default=500.0,
                        help="limite superior del eje lambda en nm (default: 500)")
    args = parser.parse_args()

    x_range = (args.xmin, args.xmax)

    if args.mol:
        df_conf, df_tddft = collect_results(args.mol)
        print_summary(df_conf, df_tddft, args.mol)

        if args.adc2:
            df_adc2 = collect_adc2_results(args.mol)
            print_adc2_summary(df_adc2, args.mol, df_conf)
            if args.export:
                export_adc2_mol(df_adc2, args.mol, args.export)
            if not df_adc2.empty:
                df_cmp_adc2 = compare_adc2_tddft(mols=[args.mol])
                print_adc2_tddft_comparison(df_cmp_adc2)
                if args.export:
                    export_adc2_comparison(df_cmp_adc2, args.export)

        if args.export:
            export_results(df_conf, df_tddft, args.mol, args.export)
        if args.plot:
            plot_spectrum(df_tddft, df_conf, args.mol,
                          sigma=args.sigma, x_range=x_range)
        if args.goodvibes:
            print(f"\n{chr(9472)*74}")
            print(f"  GOODVIBES (mRRHO Grimme) vs ORCA (RRHO) -- {args.mol}")
            print(f"{chr(9472)*74}")
            gv_results = run_goodvibes_on_mol(args.mol)
            df_gv = compute_gv_appendix(args.mol, df_conf, gv_results)
            print(df_gv.to_string(index=False))
            if args.export:
                export_gv_appendix(df_gv, args.mol, args.export)

    elif args.compare:
        result = compare_protonation_states(ref_mol=args.ref)
        if isinstance(result, tuple):
            df_cmp, dfs_conf, dfs_tddft = result
        else:
            df_cmp = result
            dfs_conf = dfs_tddft = {}

        print_comparison(df_cmp, ref_mol=args.ref)

        for mol in sorted(dfs_conf.keys()):
            print_summary(dfs_conf[mol], dfs_tddft.get(mol, pd.DataFrame()), mol)

        if args.adc2:
            mols_available = sorted(dfs_conf.keys()) or None
            df_cmp_adc2 = compare_adc2_tddft(mols=mols_available)
            print_adc2_tddft_comparison(df_cmp_adc2)
            for mol in sorted(dfs_conf.keys()):
                df_c   = dfs_conf.get(mol, pd.DataFrame())
                df_adc2 = collect_adc2_results(mol)
                print_adc2_summary(df_adc2, mol, df_c)
                if args.export and not df_adc2.empty:
                    export_adc2_mol(df_adc2, mol, args.export)
            if args.export:
                export_adc2_comparison(df_cmp_adc2, args.export)
            if args.plot:
                plot_adc2_tddft_comparison(
                    df_cmp_adc2, dfs_conf, dfs_tddft,
                    sigma=args.sigma, x_range=x_range
                )

        if args.export:
            export_comparison(df_cmp, args.export)
            for mol in sorted(dfs_conf.keys()):
                if not dfs_conf[mol].empty:
                    export_results(dfs_conf[mol],
                                   dfs_tddft.get(mol, pd.DataFrame()),
                                   mol, args.export)
        if args.plot and not args.adc2:
            plot_comparison(dfs_conf, dfs_tddft,
                            sigma=args.sigma, x_range=x_range)


if __name__ == "__main__":
    main()
