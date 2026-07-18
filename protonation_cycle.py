#!/usr/bin/env python3
"""
protonation_cycle.py — Ciclo termodinámico de pKa y distribución de Boltzmann
                        entre estados de protonación con distinta carga y N_átomos

═══════════════════════════════════════════════════════════════════════════════
FUNDAMENTO TERMODINÁMICO
═══════════════════════════════════════════════════════════════════════════════

Compara N estados de protonación de la misma molécula esencial:
  AH₂  (q=0,  N   átomos)   →  estado neutro
  AH⁻  (q=−1, N−1 átomos)   →  primer anión (primera desprotonación)
  A²⁻  (q=−2, N−2 átomos)   →  segundo anión (segunda desprotonación)

El problema central: NO se pueden comparar energías DFT brutas entre estados
con distinta carga/composición. La solución es el ciclo termodinámico de Born
(Ho & Coote, Theor. Chem. Acc. 2010), que incluye explícitamente la energía
del protón transferido.

CICLO DE BORN para AH₂ → AH⁻ + H⁺  (pKa₁):

  AH₂(gas) ─────────────────────────────────► AH⁻(gas) + H⁺(gas)
      │               ΔG_gas                        │
  ΔG_solv(AH₂)                         ΔG_solv(AH⁻) + ΔG_solv(H⁺)
      │                                              │
      ▼                                              ▼
  AH₂(aq) ─────────────────────────────────► AH⁻(aq) + H⁺(aq)
                          ΔG_aq

  ΔG_aq = ΔG_gas + ΔG_solv(AH⁻) + ΔG_solv(H⁺) − ΔG_solv(AH₂)
  pKa   = ΔG_aq / (RT · ln10)

Lo mismo para el segundo pKa usando AH⁻ y A²⁻.

ENERGÍAS COMPUESTAS (ciclo Ho & Coote):

  G_total(X) = E_SP(TZVPD, CPCM) + G_corr(SVP, CPCM)
  G_gas(X)   = E_SP(TZVPD, gas)  + G_corr(SVP, CPCM)    [aproximación RRHO]
  ΔG_solv(X) = E_SP(TZVPD, CPCM) − E_SP(TZVPD, gas) + δG_conc

  δG_conc = +1.894 kcal/mol  [corrección 1 atm → 1 M a 298.15 K]

CONSTANTES DEL PROTÓN (298.15 K):

  G(H⁺, gas)   = −6.28 kcal/mol     [Tissandier et al. J. Phys. Chem. A 1998]
  ΔG_solv(H⁺)  = −265.9 kcal/mol    [Kelly, Cramer & Truhlar J. Phys. Chem. B 2006]

BOLTZMANN DE ESTADOS DE PROTONACIÓN:

  Se calcula a partir de los pKa usando Henderson-Hasselbalch encadenado.
  La fracción de cada especie depende del pH, no de comparar G's directamente.

  Para diacido con pKa₁ < pKa₂:
    f(AH₂) = 1 / D
    f(AH⁻) = 10^(pH − pKa₁) / D
    f(A²⁻) = 10^(2·pH − pKa₁ − pKa₂) / D
    D       = 1 + 10^(pH − pKa₁) + 10^(2·pH − pKa₁ − pKa₂)

BOLTZMANN DE CONFÓRMEROS (dentro de cada estado):

  Previo al ciclo de pKa, cada estado necesita su G_Boltzmann (energía libre
  efectiva del ensemble conformacional), calculada desde el confórmero más
  estable o como promedio de Boltzmann explícito.

═══════════════════════════════════════════════════════════════════════════════
ESTRUCTURA DE DIRECTORIOS (igual que pipeline.py / results_collector.py)
═══════════════════════════════════════════════════════════════════════════════

  /home/kmilozeta/PAPER1/
  ├── MGLYVAL/
  │   └── runs/
  │       ├── R/                         ← estado de referencia/neutro
  │       │   ├── 03_numfreq/conf_001/freq.out
  │       │   ├── 04_sp/conf_001/sp.out          ← SP acuoso (CPCM)
  │       │   └── 04b_sp_gas/conf_001/sp_gas.out ← SP fase gas (sin CPCM) [NUEVO]
  │       ├── A/                         ← estado de protonación A
  │       │   └── (mismo árbol)
  │       ├── B/                         ← estado de protonación B
  │       │   └── (mismo árbol)
  │       └── C/                         ← estado de protonación C
  │           └── (mismo árbol)
  ├── SHINORINE/
  │   └── runs/
  │       └── (mismo árbol)
  └── (otras moléculas...)

WORK_ROOT se construye, por default, como <cwd>/<MOL>/runs, usando el nombre
de molécula pasado en --mol (p. ej. --mol MGLYVAL corriendo desde ~/PAPER1
da ~/PAPER1/MGLYVAL/runs). Se puede sobreescribir del todo con
--work-root /ruta/a/runs si necesitas apuntar a otro lado (--mol se ignora
en ese caso). Los nombres de estado que se pasan en
--neutral/--anion/--dianion son los nombres de carpeta dentro de runs/
(p. ej. "R", "A", "B", "C"), no rutas completas.

Para el ciclo de pKa se requieren los SP en fase gas (04b_sp_gas).
Si no existen, el módulo de pKa se omite con advertencia clara.

═══════════════════════════════════════════════════════════════════════════════
USO
═══════════════════════════════════════════════════════════════════════════════

  # Caso más común: tres estados, pKa1 + pKa2 + distribución
  python protonation_cycle.py \\
      --mol MGLYVAL \\
      --neutral  R \\
      --anion    B \\
      --dianion  C

  # Solo pKa1 (R → B), para SHINORINE
  python protonation_cycle.py \\
      --mol SHINORINE --neutral R --anion B

  # Evaluar distribución a pH específico además del rango
  python protonation_cycle.py \\
      --mol MGLYVAL --neutral R --anion B --dianion C \\
      --ph 8.1

  # Exportar resultados
  python protonation_cycle.py \\
      --mol MGLYVAL --neutral R --anion B --dianion C \\
      --export xlsx --plot

  # Modo detallado (muestra contribuciones del ciclo de Born paso a paso)
  python protonation_cycle.py \\
      --mol MGLYVAL --neutral R --anion B --dianion C \\
      --verbose
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# ─────────────────────────────────────────────
#  CONSTANTES  (idénticas a results_collector.py)
# ─────────────────────────────────────────────

WORK_ROOT        = Path.cwd() / "runs"  # default; sobreescribible con --work-root
HARTREE_TO_KCAL  = 627.5094740631
HARTREE_TO_KJ    = HARTREE_TO_KCAL * 4.184
R_KCAL           = 1.9872041e-3        # kcal/(mol·K)
LN10             = 2.302585093
T_STD            = 298.15              # K

# Temperaturas de análisis (contexto FAE antártico)
TEMPS_K     = [250.0, 274.15, 298.15, 310.0]
TEMP_LABELS = ["250K", "274K", "298K", "310K"]

# G(H⁺, aq) — DOS CONVENCIONES DISTINTAS, NO INTERCAMBIABLES:
#
# 1) Para el CICLO DE BORN (gas + ΔG_solv separados):
#      G(H⁺,gas) = -6.28 kcal/mol      [Tissandier et al. 1998]
#      ΔG_solv(H⁺) = -265.9 kcal/mol   [Kelly, Cramer & Truhlar 2006]
#      G(H⁺,aq) implícito = -272.18 kcal/mol
#    Esta convención está calibrada para combinarse con E_SP en FASE GAS
#    pura (sin CPCM). Úsala solo si G_gas_acid/G_gas_base vienen de SP sin
#    solvente implícito.
#
# 2) Para el MÉTODO DIRECTO (solo energías ya en CPCM, sin SP fase gas):
#    Las energías E_SP_CPCM(X) ya incluyen el efecto de solvatación de la
#    molécula vía el modelo de continuo. Usar G(H⁺,aq) de Tissandier/Kelly
#    aquí sobreestima la estabilización del protón en ~8-10 kcal/mol, porque
#    esa constante fue derivada para cerrar un ciclo gas+solv explícito, no
#    para emparejarse con energías que ya son "solución implícita".
#    El valor correcto para esquemas de fase condensada / esquema directo
#    es el de Pliego & Riveros, J. Phys. Chem. A 2002, 106, 7434:
#      G(H⁺,aq) = -264.61 kcal/mol  (escala molar, 298.15 K, calibrado con
#                                     ciclos de fase condensada, consistente
#                                     con Zhan & Dixon 2001 ≈ -262.4 kcal/mol)
G_HPLUS_GAS   = -6.28    # kcal/mol   G(H⁺, gas)       Tissandier et al. 1998
DG_SOLV_HPLUS = -265.9   # kcal/mol   ΔG_solv(H⁺)      Kelly et al. 2006
DG_CONC_CORR  = 1.894    # kcal/mol   δG_conc 1 atm→1M  RT·ln(24.5) a 298.15 K

# G(H⁺, aq) PARA EL CICLO DE BORN (gas + solv explícitos, NO usar con el método directo)
G_HPLUS_AQ_BORN   = G_HPLUS_GAS + DG_SOLV_HPLUS   # kcal/mol  ≈ -272.18

# G(H⁺, aq) PARA EL MÉTODO DIRECTO (energías ya en CPCM) — Pliego & Riveros 2002
G_HPLUS_AQ_DIRECT = -264.61   # kcal/mol   [Pliego & Riveros, JPCA 2002, 106, 7434]

# Alias retro-compatible: NO usar en código nuevo, mantenido solo por si algo
# externo importa G_HPLUS_AQ directamente.
G_HPLUS_AQ = G_HPLUS_AQ_DIRECT

# Subdirectorios dentro de cada mol en runs/
SUBDIR_SP      = "04_sp"
SUBDIR_FREQ    = "03_numfreq"
SUBDIR_SP_GAS  = "04b_sp_gas"  # SPs sin CPCM para el ciclo de pKa

# ─────────────────────────────────────────────
#  GOODVIBES  (corrección cuasi-armónica opcional)
# ─────────────────────────────────────────────
# G_corr estándar (parse_thermochem) usa el oscilador armónico puro para
# TODOS los modos vibracionales, incluidos los rotores internos de baja
# frecuencia (<100 cm⁻¹) — típicos en M-Gly:Val (cadenas laterales Gly/Val,
# grupos OH/NH libres). Ahí la aproximación armónica sobreestima la entropía
# y por ende puede desviar bastante el G_corr. GoodVibes (Grimme, Head-Gordon)
# corrige esto con un tratamiento cuasi-armónico (mRRHO) para esos modos.
# Estos valores son sobreescritos por --gv-cutoff/--gv-qh/--gv-qs en main().
GV_CUTOFF = 100.0      # cm⁻¹, frequency cutoff (default de GoodVibes)
GV_QH     = False      # aplicar también la corrección cuasi-armónica de
                        # Head-Gordon a la entalpía (si no, H se deja armónica)
GV_QS     = "grimme"   # método de entropía cuasi-armónica: 'grimme' o 'truhlar'

# Entorno conda donde está instalado GoodVibes. Igual que en results_collector.py:
# si se deja en None, se salta el paso 1 de _resolve_gv_python() y se prueba
# directamente con sys.executable / PATH.
GV_CONDA_ENV = "sciencebitch"

# Caché del intérprete resuelto (evita repetir la búsqueda en disco + los
# subprocess de verificación en cada llamada a run_goodvibes(), que se invoca
# una vez por confórmero — potencialmente decenas de veces por corrida).
_GV_PYTHON_CACHE: Optional[str] = None

# pH de referencia oceánico antártico (WAP superficial, verano)
PH_ANTARCTIC   = 8.1
PH_STANDARD    = 7.0


# ─────────────────────────────────────────────
#  PARSERS  (extraídos de results_collector.py,
#            sin dependencias cruzadas)
# ─────────────────────────────────────────────

def parse_sp_energy(out_file: Path) -> Optional[float]:
    """
    Extrae la última 'FINAL SINGLE POINT ENERGY' del .out de ORCA.
    Retorna energía en Hartree, o None si no se encuentra.
    """
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


def parse_theory_level(out_file: Path) -> Optional[str]:
    """
    Extrae la línea de keywords ORCA ('! funcional base CPCM ...') del .out,
    tal como aparece ecoizada al inicio del cálculo (línea '|  1> ! ...').

    Se usa solo para diagnóstico de consistencia entre estados — NO se usa
    para ningún cálculo numérico. Retorna None si no se encuentra.
    """
    if not out_file.exists():
        return None
    pat = re.compile(r"^\s*\|\s*1>\s*!\s*(.+)$")
    for line in out_file.read_text(errors="replace").splitlines():
        m = pat.match(line)
        if m:
            return m.group(1).strip()
    return None


def parse_charge_mult_natoms(out_file: Path) -> dict:
    """
    Extrae carga total, multiplicidad y número de átomos de un .out de ORCA.

    Se usa para validar que los estados de protonación pasados por
    --protonated/--neutral/--anion/--dianion/--zwitterion realmente difieren
    entre sí en la forma esperada (Δq=-1/ΔN_átomos=-1 para una desprotonación,
    Δq=0/ΔN_átomos=0 para un par de tautómeros) — ver
    _validate_state_charges(). Sin este chequeo, un mapeo de carpetas
    desalineado (p. ej. --anion apuntando al estado equivocado, o un
    best_conf que no corresponde a la geometría que se cree) se procesa sin
    quejarse y produce un pKa numéricamente válido pero físicamente absurdo.

    Prioriza la línea de eco del input ('|  N> * xyz Q M'), que es
    exactamente lo que se envió a ORCA. Si no aparece (p. ej. geometría
    leída de un archivo externo en vez de un bloque xyz inline), cae a los
    patrones de la sección de configuración SCF ('Total Charge ... Charge
    ... N', 'Multiplicity ... Mult ... N').

    n_atoms se cuenta a partir del primer bloque
    'CARTESIAN COORDINATES (ANGSTROEM)' del archivo.

    Retorna {"charge": int|None, "mult": int|None, "n_atoms": int|None}.
    """
    null = {"charge": None, "mult": None, "n_atoms": None}
    if not out_file.exists():
        return null
    try:
        lines = out_file.read_text(errors="replace").splitlines()
    except OSError:
        return null

    charge: Optional[int] = None
    mult: Optional[int] = None
    n_atoms: Optional[int] = None

    # 1) Eco del input: "|  N> * xyz Q M"
    pat_echo = re.compile(r"^\s*\|\s*\d+>\s*\*\s*xyz\s+(-?\d+)\s+(\d+)", re.IGNORECASE)
    for line in lines:
        m = pat_echo.match(line)
        if m:
            charge, mult = int(m.group(1)), int(m.group(2))
            break

    # 2) Fallback: bloque de configuración SCF ("Total Charge ... N")
    if charge is None:
        pat_q = re.compile(r"Total Charge\s+Charge\s+\.+\s*(-?\d+)", re.IGNORECASE)
        for line in lines:
            m = pat_q.search(line)
            if m:
                charge = int(m.group(1))
                break
    if mult is None:
        pat_m = re.compile(r"Multiplicity\s+Mult\s+\.+\s*(\d+)", re.IGNORECASE)
        for line in lines:
            m = pat_m.search(line)
            if m:
                mult = int(m.group(1))
                break

    # 3) Número de átomos: primer bloque CARTESIAN COORDINATES (ANGSTROEM)
    for idx, line in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line.upper():
            j = idx + 2  # saltar encabezado y línea separadora '---...'
            count = 0
            while j < len(lines) and lines[j].strip():
                count += 1
                j += 1
            n_atoms = count if count > 0 else None
            break

    return {"charge": charge, "mult": mult, "n_atoms": n_atoms}


def _validate_state_charges(
    data_a: dict,
    data_b: dict,
    dq_expected: int,
    dn_expected: int,
    label_a: str,
    label_b: str,
) -> list[str]:
    """
    Verifica que data_b difiera de data_a exactamente en dq_expected unidades
    de carga y dn_expected átomos. Ejemplos:
      - compute_pka(acid, base):              dq_expected=-1, dn_expected=-1
      - compute_tautomer_equilibrium(R, D):    dq_expected=0,  dn_expected=0

    No aborta nada por sí sola: retorna una lista de mensajes. Los mensajes
    que empiezan con '❌' son fatales — el llamador (compute_pka /
    compute_tautomer_equilibrium) los trata como error y NO calcula el
    pKa/equilibrio, en vez de seguir de largo con un resultado numérico sin
    sentido físico. Los mensajes sin '❌' son informativos (p. ej. cuando no
    se pudo parsear la carga/N° de átomos de alguno de los dos estados) y
    solo se agregan a las advertencias.
    """
    msgs: list[str] = []
    qa, qb = data_a.get("charge"), data_b.get("charge")
    na, nb = data_a.get("n_atoms"), data_b.get("n_atoms")

    if qa is None or qb is None:
        msgs.append(
            f"Could not verify charge between '{label_a}' and '{label_b}' "
            f"(charge not determined from the ORCA .out for at least one of "
            f"the two states). Manually verify that the expected charge "
            f"relation (Δq={dq_expected:+d}) is correct."
        )
    elif (qb - qa) != dq_expected:
        msgs.append(
            f"❌ INCONSISTENT CHARGE: '{label_a}' (q={qa:+d}) → '{label_b}' "
            f"(q={qb:+d}) implies Δq={qb - qa:+d}, but Δq={dq_expected:+d} "
            f"was expected. Check whether the folder mapping "
            f"(--protonated/--neutral/--anion/--dianion/--zwitterion) "
            f"corresponds to the state you think it is."
        )

    if na is None or nb is None:
        msgs.append(
            f"Could not verify atom count between '{label_a}' and '{label_b}' "
            f"(atom count not determined from the ORCA .out for at least one "
            f"of the two states). Manually verify that the expected relation "
            f"(ΔN={dn_expected:+d}) is correct."
        )
    elif (nb - na) != dn_expected:
        msgs.append(
            f"❌ INCONSISTENT ATOM COUNT: '{label_a}' (N={na}) → "
            f"'{label_b}' (N={nb}) implies ΔN={nb - na:+d}, but ΔN={dn_expected:+d} "
            f"was expected. Check whether the folder mapping "
            f"corresponds to the state you think it is, or whether the "
            f"representative conformer (best_conf) corresponds to the "
            f"expected geometry."
        )

    return msgs


def parse_thermochem(freq_out: Path, T: float = T_STD) -> dict:
    """
    Extrae correcciones termoquímicas del bloque THERMOCHEMISTRY de ORCA 6.

    Lógica idéntica a results_collector.py:
      - 'Final entropy term' es T·S a T_parsed (en Eh), NO es S puro.
      - S_total = TS_total / T_parsed   (Eh/K)
      - G_corr(T) = H_corr − T · S_total

    Retorna dict con claves:
        E_elec_freq, ZPE, H_corr, TS_total, S_total,
        G_corr, G_total_freq, T_parsed, n_imag, ok
    Todos los valores en Hartree (salvo T_parsed en K y n_imag entero).
    """
    null = {
        "E_elec_freq": None, "ZPE": None,       "H_corr": None,
        "TS_total":    None, "S_total": None,    "G_corr": None,
        "G_total_freq": None, "T_parsed": None,
        "n_imag": 0,         "ok": False,
    }

    if not freq_out.exists():
        return null

    try:
        lines = freq_out.read_text(errors="replace").splitlines()
    except OSError:
        return null

    # Frecuencias imaginarias
    n_imag = sum(
        1 for l in lines
        if re.match(r"\s+\d+:\s+(-\d+\.\d+)\s+cm", l)
    )

    # Energía electrónica (primera ocurrencia = SP del NumFreq)
    E_elec_freq: Optional[float] = None
    for line in lines:
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                E_elec_freq = float(line.split()[-1])
                break
            except ValueError:
                pass

    # Bloque THERMOCHEMISTRY (última ocurrencia)
    thermo_start = -1
    for idx, line in enumerate(lines):
        if "THERMOCHEMISTRY AT" in line.upper():
            thermo_start = idx

    if thermo_start < 0:
        result = dict(null)
        result["E_elec_freq"] = E_elec_freq
        result["n_imag"] = n_imag
        return result

    # Temperatura del bloque
    T_parsed: Optional[float] = None
    m = re.search(r"THERMOCHEMISTRY AT\s+([\d.]+)\s*K",
                  lines[thermo_start], re.IGNORECASE)
    if m:
        try:
            T_parsed = float(m.group(1))
        except ValueError:
            pass

    # Patrones ORCA 6
    _pat_zpe = re.compile(r"Zero point energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_h   = re.compile(r"Total [Ee]nthalpy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_s   = re.compile(r"Final entropy term\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_g   = re.compile(
        r"(?:Total|Final) Gibbs free energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I
    )

    ZPE = H_total = TS_total = G_total_freq = None
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

    # Derivar correcciones
    H_corr  = (H_total - E_elec_freq) if (H_total is not None and E_elec_freq is not None) else None
    S_total = (TS_total / T_parsed)   if (TS_total is not None and T_parsed and T_parsed != 0) else None
    G_corr  = (H_corr - T * S_total)  if (H_corr is not None and S_total is not None) else None

    ok = all(x is not None for x in [ZPE, H_corr, TS_total, G_corr, G_total_freq])

    return {
        "E_elec_freq":  E_elec_freq,
        "ZPE":          ZPE,
        "H_corr":       H_corr,
        "TS_total":     TS_total,
        "S_total":      S_total,       # Eh/K
        "G_corr":       G_corr,        # Eh, a temperatura T
        "G_total_freq": G_total_freq,  # Eh
        "T_parsed":     T_parsed,
        "n_imag":       n_imag,
        "ok":           ok,
    }


def _resolve_gv_python(force: bool = False) -> str:
    """
    Devuelve la ruta absoluta al intérprete Python que tiene goodvibes instalado.

    Idéntico en espíritu a _resolve_gv_python() de results_collector.py (misma
    estrategia en el mismo orden de prioridad), para que ambos scripts se
    comporten igual sea cual sea el shell/entorno desde el que se invoquen:

      1. Si GV_CONDA_ENV está definido, buscar ese entorno en los prefijos
         estándar de miniforge3, miniconda3, anaconda3 y conda en $HOME.
      2. Intentar `python -c "import goodvibes"` con el ejecutable actual
         (sys.executable), por si el entorno activo ya lo incluye.
      3. Buscar `goodvibes` en el PATH como script de consola y derivar el
         Python de su mismo directorio bin/.
      4. Fallback: sys.executable (fallará más adelante con un warning claro
         en vez de un `shutil.which` silencioso).

    El resultado se cachea en _GV_PYTHON_CACHE porque run_goodvibes() se llama
    una vez por confórmero — resolver esto en disco + subprocess en cada
    llamada sería lento e innecesario.
    """
    global _GV_PYTHON_CACHE
    if _GV_PYTHON_CACHE is not None and not force:
        return _GV_PYTHON_CACHE

    if GV_CONDA_ENV:
        home = Path.home()
        for base in ["miniforge3", "miniconda3", "anaconda3", "conda"]:
            candidate = home / base / "envs" / GV_CONDA_ENV / "bin" / "python"
            if candidate.exists():
                try:
                    r = subprocess.run(
                        [str(candidate), "-c", "import goodvibes"],
                        capture_output=True, timeout=10,
                    )
                    if r.returncode == 0:
                        _GV_PYTHON_CACHE = str(candidate)
                        return _GV_PYTHON_CACHE
                except (subprocess.TimeoutExpired, OSError):
                    pass

    try:
        r = subprocess.run(
            [sys.executable, "-c", "import goodvibes"],
            capture_output=True, timeout=10,
        )
        if r.returncode == 0:
            _GV_PYTHON_CACHE = sys.executable
            return _GV_PYTHON_CACHE
    except (subprocess.TimeoutExpired, OSError):
        pass

    gv_script = shutil.which("goodvibes")
    if gv_script:
        py_candidate = Path(gv_script).parent / "python"
        if py_candidate.exists():
            _GV_PYTHON_CACHE = str(py_candidate)
            return _GV_PYTHON_CACHE

    _GV_PYTHON_CACHE = sys.executable
    return _GV_PYTHON_CACHE


def run_goodvibes(
    freq_out: Path,
    T: float = T_STD,
    cutoff: float = GV_CUTOFF,
    qs: str = GV_QS,
    qh: bool = GV_QH,
) -> dict:
    """
    Corre GoodVibes UNA vez sobre un freq.out de ORCA y extrae de la MISMA
    tabla dos juegos de valores paralelos:

      • "armónico puro"     → columnas T.S / G(T)       (nunca amortiguadas
                               por el cutoff — así las define GoodVibes; T.S
                               es el RRHO estándar, sin mRRHO)
      • "cuasi-armónico"    → columnas T.qh-S / qh-G(T)  (mRRHO de Grimme,
                               amortiguado con `cutoff`, default 100 cm⁻¹)

    Esto es clave para que la comparación "armónico vs. GoodVibes" del
    script sea limpia: antes, la pasada "armónica" salía de parse_thermochem()
    (bloque THERMOCHEMISTRY de ORCA, que ya viene con `Quasi RRHO = True` y
    su propio cutoff interno ~1 cm⁻¹) mientras que la pasada "GoodVibes"
    salía de un cálculo aparte — dos fuentes de datos distintas, no
    comparables limpiamente. Ahora ambas pasadas vienen del mismo cálculo de
    frecuencias, la misma detección de simetría, etc.; la única diferencia
    real entre ellas es si se aplica o no el amortiguamiento de Grimme.

    Invocación (idéntica a la de terminal):

        python -m goodvibes.GoodVibes <freq.out> \\
            --qs grimme --freq_cutoff 100 -t 298.15

    Se parsea la tabla de stdout (no `--json`/`--output_format`, cuyo nombre
    de archivo de salida varía entre versiones de GoodVibes) — más simple,
    algo más lento pero robusto a esos cambios de flags.

    Retorna dict con:
        scf_energy, ZPE, n_imag, ok, warning
        G_corr_harmonic, H_corr_harmonic, S_total_harmonic, G_total_harmonic
        G_corr_qh,       H_corr_qh,       S_total_qh,       G_total_qh
    más alias retro-compatibles G_corr/H_corr/S_total/qh_gibbs_free_energy
    (= los *_qh, que es lo que antes usaba la pasada "goodvibes").
    """
    null = {
        "G_corr": None, "H_corr": None, "S_total": None,
        "scf_energy": None, "ZPE": None, "qh_gibbs_free_energy": None,
        "G_corr_harmonic": None, "H_corr_harmonic": None,
        "S_total_harmonic": None, "G_total_harmonic": None,
        "G_corr_qh": None, "H_corr_qh": None,
        "S_total_qh": None, "G_total_qh": None,
        "n_imag": 0, "ok": False, "warning": None,
    }

    if not freq_out.exists():
        return {**null, "warning": f"No existe {freq_out}"}

    python_exe = _resolve_gv_python()

    try:
        r = subprocess.run(
            [python_exe, "-c", "import goodvibes"],
            capture_output=True, timeout=10,
        )
        if r.returncode != 0:
            return {**null, "warning": (
                f"No se puede importar goodvibes con {python_exe}. "
                f"Instalar con: conda activate {GV_CONDA_ENV or '<env>'} && "
                f"pip install goodvibes --break-system-packages"
            )}
    except (subprocess.TimeoutExpired, OSError) as e:
        return {**null, "warning": f"No se pudo verificar goodvibes con {python_exe}: {e}"}

    cmd = [
        python_exe, "-m", "goodvibes.GoodVibes", str(freq_out),
        "--qs", qs, "--freq_cutoff", str(cutoff), "-t", str(T),
    ]
    if qh:
        cmd.append("--qh")

    try:
        with tempfile.TemporaryDirectory() as tmp:
            proc = subprocess.run(
                cmd, cwd=tmp, capture_output=True, text=True, timeout=180,
            )
    except (subprocess.TimeoutExpired, OSError) as e:
        return {**null, "warning": (
            f"Fallo al ejecutar goodvibes ({python_exe} -m goodvibes.GoodVibes): {e}"
        )}

    if proc.returncode != 0:
        stderr_tail = (proc.stderr or "").strip()[-400:]
        return {**null, "warning": (
            f"goodvibes ({python_exe}) terminó con returncode={proc.returncode}. "
            f"{stderr_tail}"
        )}

    stdout = proc.stdout or ""
    parsed = _parse_goodvibes_table(stdout, freq_out.stem)
    if parsed is None:
        stdout_tail = stdout.strip()[-600:]
        return {**null, "warning": (
            f"No se pudo parsear la tabla de salida de goodvibes para "
            f"{freq_out.name}. Salida (cola):\n{stdout_tail}"
        )}

    scf_energy = parsed.get("E")
    ZPE        = parsed.get("ZPE")
    G_qh       = parsed.get("qh-G(T)")
    G_harm     = parsed.get("G(T)")
    if scf_energy is None or G_qh is None or G_harm is None:
        return {**null, "warning": (
            "goodvibes no reportó columnas E / G(T) / qh-G(T) en la tabla "
            "(freq.out posiblemente incompleto o mal formado)."
        )}

    # Entalpía armónica pura: siempre la columna H estándar.
    H_harm = parsed.get("H")
    H_corr_harmonic = (H_harm - scf_energy) if H_harm is not None else None

    # Entalpía cuasi-armónica: columna qh-H (Head-Gordon) si --qh está
    # activo; si no, GoodVibes no amortigua la entalpía y
    # H_corr_qh == H_corr_harmonic.
    H_qh = parsed.get("qh-H") if "qh-H" in parsed else H_harm
    H_corr_qh = (H_qh - scf_energy) if H_qh is not None else None

    # T.S / T.qh-S ya vienen multiplicadas por T; se recupera S dividiendo por T.
    T_S    = parsed.get("T.S")
    T_qh_S = parsed.get("T.qh-S")
    S_total_harmonic = (T_S / T) if T_S is not None else None
    S_total_qh       = (T_qh_S / T) if T_qh_S is not None else None

    n_imag = 0
    m = re.search(r"(\d+)\s+imaginary frequenc", stdout, re.IGNORECASE)
    if m:
        n_imag = int(m.group(1))

    G_corr_harmonic = G_harm - scf_energy
    G_corr_qh       = G_qh - scf_energy

    return {
        "scf_energy":         scf_energy,
        "ZPE":                ZPE,
        "n_imag":             n_imag,
        "ok":                 True,
        "warning":            None,
        "G_corr_harmonic":    G_corr_harmonic,
        "H_corr_harmonic":    H_corr_harmonic,
        "S_total_harmonic":   S_total_harmonic,
        "G_total_harmonic":   G_harm,
        "G_corr_qh":          G_corr_qh,
        "H_corr_qh":          H_corr_qh,
        "S_total_qh":         S_total_qh,
        "G_total_qh":         G_qh,
        # alias retro-compatibles (= *_qh, comportamiento previo de "goodvibes")
        "G_corr":               G_corr_qh,
        "H_corr":               H_corr_qh,
        "S_total":              S_total_qh,
        "qh_gibbs_free_energy": G_qh,
    }


def _parse_goodvibes_table(stdout: str, structure_stem: str) -> Optional[dict]:
    """
    Parsea la tabla de resultados que GoodVibes imprime por stdout, del
    estilo:

          Structure              E        ZPE              H        T.S     T.qh-S           G(T)        qh-G(T)
        ────────────────────────────────────────────────────────────────────────────────────────────────────────
             freq        -1220.971714   0.402267   -1220.542612   0.085076   0.079189   -1220.627688   -1220.621801
        ────────────────────────────────────────────────────────────────────────────────────────────────────────

    Devuelve un dict {nombre_columna: valor_float} para la fila cuya
    columna "Structure" coincide con `structure_stem` (nombre del archivo
    sin extensión; en la convención del pipeline siempre es "freq", ya que
    el archivo es conf_NNN/freq.out), o None si la tabla no aparece o no se
    pudo parsear.

    NOTA: GoodVibes antepone un marcador de estado de una sola columna
    (aparece como un bullet/checkmark, p. ej. una "o" suelta) ANTES del
    nombre de estructura en cada fila de datos:

        o     freq        -1220.967288   0.403345   ...

    Si se asumiera que la estructura es siempre el primer token, esa "o"
    correría el conteo de columnas en uno y la fila se descartaría por
    "no calzar" con el header — silenciosamente, sin fallar con un error
    claro. Por eso acá se busca `structure_stem` en la posición que sea
    dentro de los tokens, y se toman los valores inmediatamente después.
    """
    lines = stdout.splitlines()
    header_tokens = None
    header_idx = None
    for i, line in enumerate(lines):
        if "Structure" in line and "G(T)" in line:
            header_tokens = re.split(r"\s{2,}", line.strip())
            header_idx = i
            break
    if header_tokens is None:
        return None

    n_values = len(header_tokens) - 1  # columnas numéricas tras "Structure"

    for line in lines[header_idx + 1:]:
        stripped = line.strip()
        if not stripped or set(stripped) <= {"─", "-", "_"}:
            continue
        tokens = stripped.split()
        if structure_stem not in tokens:
            continue
        idx = tokens.index(structure_stem)
        value_tokens = tokens[idx + 1:]
        if len(value_tokens) != n_values:
            continue
        values: dict = {"Structure": structure_stem}
        for key, val in zip(header_tokens[1:], value_tokens):
            try:
                values[key] = float(val)
            except ValueError:
                values[key] = val
        return values

    return None


def get_thermochem(
    freq_out: Path,
    T: float = T_STD,
    source: str = "harmonic",
    gv_cutoff: float = GV_CUTOFF,
    gv_qs: str = GV_QS,
    gv_qh: bool = GV_QH,
) -> dict:
    """
    Wrapper único para obtener G_corr, con dos fuentes intercambiables que
    ahora comparten la MISMA corrida de GoodVibes sobre el mismo freq.out
    (ver run_goodvibes()) — solo difieren en qué columnas de esa tabla usan:

      source="harmonic"  → columnas T.S / G(T) de GoodVibes
                            (RRHO puro, NUNCA amortiguado por el cutoff)
      source="goodvibes" → columnas T.qh-S / qh-G(T) de GoodVibes
                            (mRRHO de Grimme, amortiguado con gv_cutoff)

    Antes, "harmonic" salía de parse_thermochem() (bloque THERMOCHEMISTRY
    de ORCA, que ya trae `Quasi RRHO = True` con su propio cutoff interno
    ~1 cm⁻¹) mientras que "goodvibes" salía de un cálculo aparte — dos
    fuentes de datos distintas, no comparables. parse_thermochem() se deja
    definida arriba solo como utilidad de referencia/inspección manual, ya
    no la usa esta función.

    Retorna siempre el shape de parse_thermochem() (E_elec_freq, ZPE,
    H_corr, TS_total, S_total, G_corr, G_total_freq, T_parsed, n_imag, ok),
    más "source" y "warning" para trazabilidad. Esto permite que
    collect_mol_data(), _state_G_total_direct() y compute_pka() usen la
    misma interfaz sin importar la fuente.
    """
    if source not in ("harmonic", "goodvibes"):
        raise ValueError(f"source debe ser 'harmonic' o 'goodvibes', no {source!r}")

    gv = run_goodvibes(freq_out, T=T, cutoff=gv_cutoff, qs=gv_qs, qh=gv_qh)

    if source == "harmonic":
        S_total  = gv.get("S_total_harmonic")
        TS_total = (T * S_total) if S_total is not None else None
        return {
            "E_elec_freq":  gv.get("scf_energy"),
            "ZPE":          gv.get("ZPE"),
            "H_corr":       gv.get("H_corr_harmonic"),
            "TS_total":     TS_total,
            "S_total":      S_total,
            "G_corr":       gv.get("G_corr_harmonic"),
            "G_total_freq": gv.get("G_total_harmonic"),
            "T_parsed":     T,
            "n_imag":       gv.get("n_imag", 0),
            "ok":           gv.get("ok", False),
            "source":       "harmonic",
            "warning":      gv.get("warning"),
        }

    # source == "goodvibes"
    S_total  = gv.get("S_total_qh")
    TS_total = (T * S_total) if S_total is not None else None
    return {
        "E_elec_freq":  gv.get("scf_energy"),
        "ZPE":          gv.get("ZPE"),
        "H_corr":       gv.get("H_corr_qh"),
        "TS_total":     TS_total,
        "S_total":      S_total,
        "G_corr":       gv.get("G_corr_qh"),
        "G_total_freq": gv.get("G_total_qh"),
        "T_parsed":     T,
        "n_imag":       gv.get("n_imag", 0),
        "ok":           gv.get("ok", False),
        "source":       "goodvibes",
        "warning":      gv.get("warning"),
    }


# ─────────────────────────────────────────────
#  COLECTOR POR MOL  (versión ligera, sin TDDFT
#  ni CREST reranking — solo lo necesario para
#  el ciclo de pKa)
# ─────────────────────────────────────────────

def collect_mol_data(mol_name: str, verbose: bool = False,
                      gcorr_source: str = "harmonic") -> dict:
    """
    Recolecta los datos termodinámicos de todos los confórmeros de una molécula.

    Para cada confórmero con SP acuoso:
      - E_sp_aq   : energía SP con CPCM (Hartree)
      - E_sp_gas  : energía SP sin CPCM (Hartree, de 04b_sp_gas — opcional)
      - G_corr(T) : corrección de Gibbs de NumFreq, re-evaluada a cada T (Hartree)
      - G_total(T): E_sp_aq + G_corr(T)    ciclo compuesto Ho & Coote

    Retorna dict con:
        mol_name        : str
        conformers      : list[dict]  — datos por confórmero
        df              : pd.DataFrame
        best_conf_298K  : str — id del confórmero con menor G_total a 298.15 K
        G_boltz(T)      : dict {label: G_boltzmann en Hartree} — ensemble efectivo
        ok              : bool
        warnings        : list[str]
    """
    mol_dir = WORK_ROOT / mol_name
    warns: list[str] = []

    if not mol_dir.exists():
        return {"mol_name": mol_name, "ok": False,
                "warnings": [f"Directory not found: {mol_dir}"],
                "conformers": [], "df": pd.DataFrame()}

    sp_dirs = sorted((mol_dir / SUBDIR_SP).glob("conf_*")) \
              if (mol_dir / SUBDIR_SP).exists() else []

    if not sp_dirs:
        return {"mol_name": mol_name, "ok": False,
                "warnings": [f"No SP conformers in {mol_dir / SUBDIR_SP}"],
                "conformers": [], "df": pd.DataFrame()}

    rows: list[dict] = []

    for sp_dir in sp_dirs:
        conf_id = sp_dir.name

        # SP acuoso (CPCM)
        sp_aq_file  = sp_dir / "sp.out"
        e_sp_aq     = parse_sp_energy(sp_aq_file)

        # SP fase gas (sin CPCM) — para el ciclo de pKa
        sp_gas_file = mol_dir / SUBDIR_SP_GAS / conf_id / "sp_gas.out"
        e_sp_gas    = parse_sp_energy(sp_gas_file)

        # Correcciones térmicas — leídas a T_STD, S_total re-escalado.
        # gcorr_source="harmonic" → columnas T.S/G(T) de GoodVibes (RRHO
        # puro, sin amortiguar); "goodvibes" → T.qh-S/qh-G(T) (mRRHO
        # cuasi-armónico). Ambas vienen de la misma corrida de GoodVibes
        # (ver get_thermochem()/run_goodvibes()), así que un fallo de
        # GoodVibes afecta a las dos fuentes por igual.
        freq_file = mol_dir / SUBDIR_FREQ / conf_id / "freq.out"
        tc = get_thermochem(freq_file, T=T_STD, source=gcorr_source,
                            gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)
        if not tc["ok"] and tc.get("warning"):
            warns.append(f"[{mol_name}/{conf_id}] GoodVibes ({gcorr_source}): {tc['warning']}")

        # Carga/multiplicidad/N° átomos — para validación cruzada entre
        # estados (ver _validate_state_charges). Se intenta primero desde el
        # SP acuoso y, si falla, desde el freq.out del mismo confórmero.
        cmn = parse_charge_mult_natoms(sp_aq_file)
        if cmn["charge"] is None and cmn["n_atoms"] is None:
            cmn = parse_charge_mult_natoms(freq_file)

        row: dict = {
            "conformer":   conf_id,
            "E_sp_aq":     e_sp_aq,
            "E_sp_gas":    e_sp_gas,
            "thermo_ok":   tc["ok"],
            "n_imag":      tc["n_imag"],
            "H_corr":      tc["H_corr"],
            "S_total_EhK": tc["S_total"],
            "ZPE":         tc["ZPE"],
            "charge":      cmn["charge"],
            "mult":        cmn["mult"],
            "n_atoms":     cmn["n_atoms"],
        }

        # G_total compuesto a cada temperatura
        for T, lbl in zip(TEMPS_K, TEMP_LABELS):
            if e_sp_aq is not None and tc["ok"]:
                G_corr_T  = tc["H_corr"] - T * tc["S_total"]
                row[f"G_total_{lbl}"] = e_sp_aq + G_corr_T
            else:
                row[f"G_total_{lbl}"] = None

        rows.append(row)

        if verbose:
            gas_str = f"{e_sp_gas:.8f}" if e_sp_gas is not None else "—"
            aq_str  = f"{e_sp_aq:.8f}"  if e_sp_aq is not None  else "—"
            print(f"    {conf_id}: E_aq={aq_str}  E_gas={gas_str}  "
                  f"thermo={'✓' if tc['ok'] else '✗'}  n_imag={tc['n_imag']}")

    df = pd.DataFrame(rows)

    # Pesos de Boltzmann conformacional a cada T
    G_boltz: dict[str, Optional[float]] = {}
    for T, lbl in zip(TEMPS_K, TEMP_LABELS):
        col = f"G_total_{lbl}"
        valid = df[col].notna() & df["thermo_ok"]
        if valid.any():
            G_vals = df.loc[valid, col].values.astype(float)
            G_min  = G_vals.min()
            dG     = G_vals - G_min                           # en Hartree
            dG_kcal= dG * HARTREE_TO_KCAL
            weights = np.exp(-dG_kcal / (R_KCAL * T))
            weights /= weights.sum()
            df.loc[valid, f"boltz_{lbl}"] = weights
            df.loc[valid, f"dG_{lbl}_kcal"] = dG_kcal

            # G_Boltzmann del ensemble (free energy efectivo):
            # -RT·ln(Σ exp(-G_i/RT)) = G_min - RT·ln(Σ exp(-ΔG_i/RT))
            ln_Z   = np.log(np.sum(np.exp(-dG_kcal / (R_KCAL * T))))
            G_boltz[lbl] = G_min - (R_KCAL * T * ln_Z) / HARTREE_TO_KCAL
        else:
            G_boltz[lbl] = None
            if not df["thermo_ok"].any():
                warns.append(
                    f"[{mol_name}] No valid thermal corrections. "
                    f"Conformational Boltzmann distribution not available."
                )

    # Confórmero representativo (mínimo G_total a 298.15 K)
    col_298 = "G_total_298K"
    best_conf = None
    if col_298 in df.columns:
        valid = df[col_298].notna() & df["thermo_ok"]
        if valid.any():
            best_conf = df.loc[df.loc[valid, col_298].idxmin(), "conformer"]

    # Advertencia sobre frecuencias imaginarias
    n_imag_any = (df["n_imag"] > 0).sum()
    if n_imag_any > 0:
        warns.append(
            f"[{mol_name}] {n_imag_any} conformer(s) with imaginary frequencies. "
            f"Check geometries before using them in the pKa cycle."
        )

    # Warning if gas-phase SP is missing
    n_missing_gas = df["E_sp_gas"].isna().sum()
    if n_missing_gas > 0:
        warns.append(
            f"[{mol_name}] {n_missing_gas}/{len(df)} conformers without gas-phase SP "
            f"(expected at runs/{mol_name}/{SUBDIR_SP_GAS}/conf_NNN/sp_gas.out). "
            f"The pKa cycle requires the gas-phase SP of the representative conformer."
        )

    # ── Carga / multiplicidad / N° átomos a nivel de estado ────────────────
    # Todos los confórmeros de un mismo estado deben compartir carga, mult.
    # y N° de átomos (son el mismo sistema químico, solo geometrías
    # distintas). Si no es así, algo está mal en runs/<mol_name>/ (p. ej.
    # confórmeros de dos estados de protonación mezclados en la misma
    # carpeta) y NO se puede confiar en ningún valor único para la
    # validación cruzada con estados vecinos (_validate_state_charges).
    charges_seen  = sorted(set(df["charge"].dropna().astype(int)))
    mults_seen    = sorted(set(df["mult"].dropna().astype(int)))
    natoms_seen   = sorted(set(df["n_atoms"].dropna().astype(int)))

    state_charge  = charges_seen[0] if len(charges_seen) == 1 else None
    state_mult    = mults_seen[0]   if len(mults_seen)   == 1 else None
    state_n_atoms = natoms_seen[0]  if len(natoms_seen)  == 1 else None

    if len(charges_seen) > 1:
        warns.append(
            f"[{mol_name}] Inconsistent charge across conformers: "
            f"{charges_seen}. Cannot validate against neighboring states; "
            f"check runs/{mol_name}/{SUBDIR_SP}/conf_*/sp.out."
        )
    if len(natoms_seen) > 1:
        warns.append(
            f"[{mol_name}] Inconsistent atom count across conformers: "
            f"{natoms_seen}. Cannot validate against neighboring states; "
            f"check runs/{mol_name}/{SUBDIR_SP}/conf_*/sp.out."
        )
    if not charges_seen or not natoms_seen:
        warns.append(
            f"[{mol_name}] Could not determine charge and/or atom count "
            f"from any .out file. Cross-validation of charge/atom count "
            f"against neighboring states will be disabled for this "
            f"state (only an informational warning will be issued, not "
            f"an error)."
        )

    return {
        "mol_name":     mol_name,
        "conformers":   rows,
        "df":           df,
        "best_conf":    best_conf,
        "G_boltz":      G_boltz,      # Hartree, ensemble efectivo por T
        "charge":       state_charge,
        "mult":         state_mult,
        "n_atoms":      state_n_atoms,
        "ok":           len(rows) > 0,
        "warnings":     warns,
        "gcorr_source": gcorr_source,
    }


# ─────────────────────────────────────────────
#  ENERGÍA LIBRE REPRESENTATIVA DE UN ESTADO (método directo)
# ─────────────────────────────────────────────

def _state_G_total_direct(data: dict, T: float = T_STD,
                           gcorr_source: str = "harmonic",
                           use_ensemble: bool = False) -> dict:
    """
    G_total_CPCM del estado, para el método directo. Funciona con dos tipos
    de entrada:

      1. Un dict "normal" salido de collect_mol_data(): por default
         (use_ensemble=False) se relee freq.out del confórmero representativo
         (best_conf) y se recalcula G_corr a T (E_SP_aq(best_conf) + G_corr).

         Si use_ensemble=True, en vez de eso se usa data["G_boltz"][T_label]
         — el G_Boltzmann del ENSEMBLE conformacional completo, ya calculado
         por collect_mol_data() (línea ~1005-1008) — evitando descartar la
         población de confórmeros que no son el mínimo. Si no hay una entrada
         válida de G_boltz para la T pedida, cae automáticamente al
         comportamiento de best_conf y lo deja registrado en 'warn' (no falla
         silenciosamente).

      2. Un dict "macroestado" (data["is_macrostate"] = True), construido por
         compute_tautomer_equilibrium(): ya trae G_total_direct_Eh precalculado
         (combinación de Boltzmann entre neutro y zwitterión) y se usa tal cual,
         sin volver a tocar archivos. use_ensemble no aplica aquí (el macro-
         estado ya es, por construcción, una combinación de Boltzmann).

    Retorna dict: G_total_Eh, conf, theory_sp, theory_freq, ok, warn (str|None).
    """
    if data.get("is_macrostate"):
        return {
            "G_total_Eh":  data["G_total_direct_Eh"],
            "conf":        data.get("conf_label", "macrostate"),
            "theory_sp":   data.get("theory_label"),
            "theory_freq": data.get("theory_label"),
            "ok":          True,
            "warn":        None,
        }

    mol  = data["mol_name"]
    conf = data.get("best_conf")
    if conf is None:
        return {"ok": False, "warn": f"No representative conformer for {mol}."}

    freq_file   = WORK_ROOT / mol / SUBDIR_FREQ / conf / "freq.out"
    theory_sp   = parse_theory_level(WORK_ROOT / mol / SUBDIR_SP / conf / "sp.out")
    theory_freq = parse_theory_level(freq_file)

    warn_fallback = None
    if use_ensemble:
        T_label = next(
            (lbl for T_i, lbl in zip(TEMPS_K, TEMP_LABELS) if np.isclose(T_i, T, atol=1e-3)),
            None,
        )
        G_ens = data.get("G_boltz", {}).get(T_label) if T_label is not None else None
        if G_ens is not None:
            return {
                "G_total_Eh":  G_ens,
                "conf":        f"ensemble Boltzmann ({conf}=best_conf)",
                "theory_sp":   theory_sp,
                "theory_freq": theory_freq,
                "ok":          True,
                "warn":        None,
            }
        warn_fallback = (
            f"use_ensemble=True but no valid G_boltz at T={T:.2f} K "
            f"for {mol} (T_label={T_label!r}); using the representative "
            f"conformer (best_conf={conf}) instead."
        )

    df = data["df"]
    row = df[df["conformer"] == conf].iloc[0]
    E_aq = row["E_sp_aq"]
    if E_aq is None or pd.isna(E_aq):
        return {"ok": False, "warn": f"No aqueous SP for {mol}/{conf}."}

    tc = get_thermochem(freq_file, T=T, source=gcorr_source,
                        gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)
    G_corr = tc["G_corr"] if tc["ok"] else None
    warn = None
    if not tc["ok"]:
        extra = f" ({tc['warning']})" if tc.get("warning") else ""
        warn = f"G_corr ({gcorr_source}) not available for {mol}/{conf}. Using G≈E_SP instead (no ZPE/entropy).{extra}"

    if warn_fallback:
        warn = f"{warn_fallback} {warn}" if warn else warn_fallback

    G_total = E_aq + (G_corr if G_corr is not None else 0.0)

    return {
        "G_total_Eh":  G_total,
        "conf":        conf,
        "theory_sp":   theory_sp,
        "theory_freq": theory_freq,
        "ok":          True,
        "warn":        warn,
    }


# ─────────────────────────────────────────────
#  EQUILIBRIO DE TAUTÓMEROS  (neutro ⇌ zwitterión, sin H⁺)
# ─────────────────────────────────────────────

def compute_tautomer_equilibrium(
    data_R: dict,
    data_D: dict,
    T: float = T_STD,
    verbose: bool = False,
    gcorr_source: str = "harmonic",
) -> dict:
    """
    Compara el neutro (R) y el zwitterión (D) como TAUTÓMEROS: mismo q=0,
    mismo número de átomos (solo se desplaza un H interno, p. ej. de COOH a
    N). A diferencia del ciclo de pKa, esto NO involucra H⁺ de bulk — es una
    comparación directa de G_total_CPCM entre confórmeros representativos.

      ΔG_taut = G_total(D) − G_total(R)     [kcal/mol]
      K_taut  = [D]/[R] = exp(−ΔG_taut / RT)
      f_R     = 1 / (1 + K_taut)
      f_D     = K_taut / (1 + K_taut)

    También construye el "macroestado neutro" combinado (R+D) vía Boltzmann,
    que se usa como base/ácido en compute_pka() para los pasos con el catión
    y con el anión:

      G_macro = G_min − RT·ln( Σ exp(−(G_i − G_min)/RT) )   i ∈ {R, D}

    Nota: f_R y f_D son independientes del pH (es un equilibrio intramolecular),
    así que se calculan una sola vez aquí y se reusan en toda la distribución.
    """
    mol_R = data_R["mol_name"]
    mol_D = data_D["mol_name"]
    warns: list[str] = []
    result: dict = {
        "mol_R": mol_R, "mol_D": mol_D, "T": T, "ok": False, "warnings": warns,
        "conf_R": None, "conf_D": None, "G_R_Eh": None, "G_D_Eh": None,
        "dG_taut_kcal": None, "K_taut": None, "f_R": None, "f_D": None,
        "G_macro_Eh": None, "theory_R": None, "theory_D": None,
    }

    def _warn(msg: str):
        warns.append(msg)
        warnings.warn(f"[tautómero {mol_R}⇌{mol_D}] {msg}")

    # ── Validación de carga / N° de átomos entre R y D ──────────────────────
    # R y D son tautómeros: mismo q, mismo N° de átomos (solo se desplaza un
    # H interno). Si difieren, no son tautómeros del mismo par ácido-base y
    # el "equilibrio" que se calcularía no tendría sentido físico.
    charge_msgs = _validate_state_charges(
        data_R, data_D, dq_expected=0, dn_expected=0,
        label_a=mol_R, label_b=mol_D,
    )
    for m in charge_msgs:
        _warn(m)
    if any(m.startswith("❌") for m in charge_msgs):
        _warn(
            f"Aborting tautomer equilibrium ({mol_R}⇌{mol_D}) due to "
            f"charge/atom-count inconsistency."
        )
        return result

    gR = _state_G_total_direct(data_R, T=T, gcorr_source=gcorr_source)
    gD = _state_G_total_direct(data_D, T=T, gcorr_source=gcorr_source)

    if not gR.get("ok"):
        _warn(gR.get("warn") or f"Could not obtain G_total for {mol_R}.")
        return result
    if not gD.get("ok"):
        _warn(gD.get("warn") or f"Could not obtain G_total for {mol_D}.")
        return result
    if gR.get("warn"):
        warns.append(gR["warn"])
    if gD.get("warn"):
        warns.append(gD["warn"])

    if gR["theory_sp"] and gD["theory_sp"] and gR["theory_sp"] != gD["theory_sp"]:
        _warn(
            f"INCONSISTENT THEORY LEVEL in the aqueous SP between "
            f"{mol_R} ('{gR['theory_sp']}') and {mol_D} ('{gD['theory_sp']}'). "
            f"The tautomer comparison is not valid if basis set/functional/CPCM differ."
        )

    G_R = gR["G_total_Eh"]
    G_D = gD["G_total_Eh"]
    result["conf_R"], result["conf_D"] = gR["conf"], gD["conf"]
    result["G_R_Eh"], result["G_D_Eh"] = G_R, G_D
    result["theory_R"], result["theory_D"] = gR["theory_sp"], gD["theory_sp"]

    dG_taut_kcal = (G_D - G_R) * HARTREE_TO_KCAL
    K_taut = np.exp(-dG_taut_kcal / (R_KCAL * T))
    f_R = 1.0 / (1.0 + K_taut)
    f_D = K_taut / (1.0 + K_taut)

    # Macroestado combinado (Boltzmann de dos "confórmeros" R y D)
    G_min = min(G_R, G_D)
    dG_R_kcal = (G_R - G_min) * HARTREE_TO_KCAL
    dG_D_kcal = (G_D - G_min) * HARTREE_TO_KCAL
    lnZ = np.log(np.exp(-dG_R_kcal / (R_KCAL * T)) + np.exp(-dG_D_kcal / (R_KCAL * T)))
    G_macro = G_min - (R_KCAL * T * lnZ) / HARTREE_TO_KCAL

    result.update({
        "dG_taut_kcal": dG_taut_kcal, "K_taut": K_taut,
        "f_R": f_R, "f_D": f_D, "G_macro_Eh": G_macro, "ok": True,
    })

    if verbose:
        _print_tautomer_cycle(result)

    return result


def _print_tautomer_cycle(r: dict):
    """Imprime el equilibrio de tautómeros neutro⇌zwitterión."""
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  TAUTOMER EQUILIBRIUM (no H⁺) — {r['mol_R']} ⇌ {r['mol_D']}")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. R (neutral)     : {r['conf_R']}")
    print(f"  Conf. D (zwitterion)  : {r['conf_D']}")
    print(f"  SP level (R) : {r['theory_R'] or '—'}")
    print(f"  SP level (D) : {r['theory_D'] or '—'}")
    print()
    print(f"  {'G_total_CPCM(R)  [Eh]':<40} {r['G_R_Eh']:>14.8f}")
    print(f"  {'G_total_CPCM(D)  [Eh]':<40} {r['G_D_Eh']:>14.8f}")
    print(f"  {'─'*40}  {'─'*14}")
    print(f"  {'ΔG_taut = G(D) − G(R)  [kcal/mol]':<40} {r['dG_taut_kcal']:>+14.3f}")
    print(f"  {'K_taut = [D]/[R]':<40} {r['K_taut']:>14.4g}")
    print()
    print(f"  f(R, neutral)     = {r['f_R']*100:6.2f} %")
    print(f"  f(D, zwitterion)  = {r['f_D']*100:6.2f} %")
    dominant = "R (neutral)" if r["f_R"] >= r["f_D"] else "D (zwitterion)"
    print(f"  Dominant tautomer: {dominant}")
    print(f"\n  Neutral macrostate G (Boltzmann R+D) = {r['G_macro_Eh']:.8f} Eh")
    print(f"  (this combined value is what enters the pKa cycle as the")
    print(f"   'neutral/zwitterion' state versus the cation and the anion)")
    if r["warnings"]:
        print(f"\n  Warnings:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


# ─────────────────────────────────────────────
#  CICLO TERMODINÁMICO DE pKa
# ─────────────────────────────────────────────

def compute_pka(
    data_acid: dict,
    data_base: dict,
    T: float = T_STD,
    verbose: bool = False,
    method: str = "auto",
    gcorr_source: str = "harmonic",
) -> dict:
    """
    Calcula el pKa de AH → A⁻ + H⁺.

    Dos métodos disponibles:

    1. "born"   — Ciclo de Born completo (Ho & Coote 2010). Requiere SP en
                  fase gas (04b_sp_gas) además del SP acuoso (CPCM):
                    ΔG_solv(X) = E_SP(CPCM) − E_SP(gas) + δG_conc
                    ΔG_aq = ΔG_gas + ΔG_solv(A⁻) + ΔG_solv(H⁺) − ΔG_solv(AH)
                  Usa G(H⁺,aq) implícito de Tissandier 1998 + Kelly 2006
                  (−272.18 kcal/mol), calibrado para cerrar el ciclo gas+solv
                  explícito.

    2. "direct" — Método directo solo-CPCM. No requiere SP fase gas; usa
                  directamente G_total ya calculado con CPCM (E_SP_CPCM + G_corr)
                  y la energía libre del protón de Pliego & Riveros 2002
                  (G_HPLUS_AQ_DIRECT = −264.61 kcal/mol):
                    ΔG_aq = G_total_CPCM(A⁻) − G_total_CPCM(AH) + G(H⁺, aq)
                  IMPORTANTE: esta constante NO es la misma que la de "born".
                  Las energías E_SP_CPCM ya incluyen el efecto de solvatación
                  de la molécula vía el modelo de continuo; usar la constante
                  de Tissandier/Kelly aquí sobreestima la estabilización del
                  protón en ~8-10 kcal/mol porque esa constante está calibrada
                  para emparejarse con energías en fase gas pura, no con
                  energías ya en solución implícita.

    method="auto" (default): usa "born" si hay SP fase gas disponibles para
    ambos confórmeros representativos; si no, cae automáticamente a "direct"
    y lo deja registrado en result["method"] / result["warnings"].

    Usa el confórmero representativo (mínimo G_total a 298 K) de cada estado.
    Retorna dict con todas las contribuciones al ciclo y el pKa final.
    """
    mol_acid = data_acid["mol_name"]
    mol_base = data_base["mol_name"]
    warns: list[str] = []
    result: dict = {
        "mol_acid": mol_acid, "mol_base": mol_base,
        "T": T, "pKa": None, "ok": False, "warnings": warns,
        "method": None,
        # contribuciones (ciclo de Born)
        "dG_gas_kcal": None, "dG_solv_acid_kcal": None,
        "dG_solv_base_kcal": None, "dG_solv_hplus_kcal": DG_SOLV_HPLUS,
        "dG_aq_kcal": None,
        # contribuciones (método directo)
        "G_total_acid_Eh": None, "G_total_base_Eh": None,
        "G_hplus_aq_kcal": None,
        # datos intermedios
        "conf_acid": None, "conf_base": None,
        "G_gas_acid_Eh": None, "G_gas_base_Eh": None,
        "E_sp_aq_acid_Eh": None, "E_sp_aq_base_Eh": None,
        "E_sp_gas_acid_Eh": None, "E_sp_gas_base_Eh": None,
        "G_corr_acid_Eh": None, "G_corr_base_Eh": None,
        "theory_sp_acid": None, "theory_sp_base": None,
        "theory_freq_acid": None, "theory_freq_base": None,
    }

    def _warn(msg: str):
        warns.append(msg)
        warnings.warn(f"[pKa {mol_acid}→{mol_base}] {msg}")

    if method not in ("auto", "born", "direct"):
        raise ValueError(f"method debe ser 'auto', 'born' o 'direct', no {method!r}")

    # ── Validación de carga / N° de átomos entre ácido y base ──────────────
    # Un pKa modela AH → A⁻ + H⁺: la base SIEMPRE debe tener Δq=-1 y
    # ΔN_átomos=-1 respecto al ácido. Este chequeo corre ANTES de cualquier
    # cálculo, para atrapar de inmediato un mapeo de carpetas desalineado
    # (p. ej. --anion apuntando al estado equivocado, o un best_conf que no
    # corresponde a la geometría esperada) en vez de dejar que el script
    # calcule un pKa numéricamente válido pero físicamente sin sentido.
    charge_msgs = _validate_state_charges(
        data_acid, data_base, dq_expected=-1, dn_expected=-1,
        label_a=mol_acid, label_b=mol_base,
    )
    for m in charge_msgs:
        _warn(m)
    if any(m.startswith("❌") for m in charge_msgs):
        _warn(
            f"Aborting pKa calculation ({mol_acid}→{mol_base}) due to "
            f"charge/atom-count inconsistency. Check the state mapping "
            f"before continuing."
        )
        return result

    # ── Rama simplificada: uno de los dos lados es un macroestado combinado ─
    # (neutro+zwitterión, construido por compute_tautomer_equilibrium).
    #
    # method='born' SÍ se puede resolver aquí, y de forma EXACTA, sin ningún
    # SP de fase gas para la mezcla. La razón es la misma identidad algebraica
    # ya documentada arriba de compute_pka_montecarlo() para el caso de una
    # sola especie: con las definiciones que usa este script,
    #
    #   G_gas(X)    = E_gas(X) + G_corr(X)
    #   ΔG_solv(X)  = E_aq(X) − E_gas(X) + δG_conc
    #
    # TODOS los términos en E_gas se cancelan al construir ΔG_aq del ciclo de
    # Born, y lo que queda es exactamente la resta de energías CPCM (idéntica
    # a la del método directo) más la constante G(H⁺,aq) correspondiente:
    #
    #   ΔG_aq(born)   = G_total_CPCM(base) − G_total_CPCM(acid) + G_HPLUS_AQ_BORN
    #   ΔG_aq(direct) = G_total_CPCM(base) − G_total_CPCM(acid) + G_HPLUS_AQ_DIRECT
    #
    # Esa cancelación se extiende sin cambios al macroestado combinado: si se
    # define el "G_gas" de la mezcla con la MISMA entropía de mezcla (mismos
    # pesos f_R/f_D) que ya usa G_macro_Eh en compute_tautomer_equilibrium, el
    # término de entropía de mezcla aparece igual en G_macro_aq y en
    # G_macro_gas y se cancela exactamente igual que G_corr en el caso de una
    # sola especie. No hace falta calcular ni siquiera construir ese
    # G_macro_gas explícitamente — el resultado final es el mismo de arriba,
    # con G_total_CPCM(macro) = G_macro_Eh ya calculado por
    # compute_tautomer_equilibrium(). Por eso tampoco hacía falta el SP de
    # fase gas para el lado que NO es macroestado: esos términos también se
    # cancelan ahí.
    if data_acid.get("is_macrostate") or data_base.get("is_macrostate"):
        # use_ensemble=True: para el lado que NO es macroestado (A o B), usa
        # el G_Boltzmann del ensemble conformacional completo en vez de solo
        # el best_conf — ver nota #3 del review. Para el lado que SÍ es
        # macroestado, _state_G_total_direct ignora use_ensemble (ya entra
        # por la rama is_macrostate de arriba).
        g_acid = _state_G_total_direct(data_acid, T=T, gcorr_source=gcorr_source,
                                       use_ensemble=True)
        g_base = _state_G_total_direct(data_base, T=T, gcorr_source=gcorr_source,
                                       use_ensemble=True)

        result["conf_acid"] = g_acid.get("conf")
        result["conf_base"] = g_base.get("conf")
        result["theory_sp_acid"] = g_acid.get("theory_sp")
        result["theory_sp_base"] = g_base.get("theory_sp")

        if not g_acid.get("ok"):
            _warn(g_acid.get("warn") or f"Could not obtain G_total for {mol_acid}.")
            return result
        if not g_base.get("ok"):
            _warn(g_base.get("warn") or f"Could not obtain G_total for {mol_base}.")
            return result
        if g_acid.get("warn"):
            warns.append(g_acid["warn"])
        if g_base.get("warn"):
            warns.append(g_base["warn"])

        if (g_acid["theory_sp"] and g_base["theory_sp"]
                and g_acid["theory_sp"] != g_base["theory_sp"]):
            _warn(
                f"INCONSISTENT THEORY LEVEL between {mol_acid} "
                f"('{g_acid['theory_sp']}') and {mol_base} ('{g_base['theory_sp']}')."
            )

        G_total_acid = g_acid["G_total_Eh"]
        G_total_base = g_base["G_total_Eh"]
        result["G_total_acid_Eh"] = G_total_acid
        result["G_total_base_Eh"] = G_total_base

        if method == "born":
            use_born_macro = True
        elif method == "direct":
            use_born_macro = False
        else:  # "auto"
            # Exacto sin SP de fase gas, por la derivación algebraica de
            # arriba: el ciclo de Born se reduce a la misma resta CPCM que
            # 'direct', solo con la constante G_HPLUS_AQ_BORN en vez de
            # G_HPLUS_AQ_DIRECT. No hay razón para que 'auto' degrade a
            # 'direct' aquí — a diferencia de la rama simple, donde 'auto'
            # SÍ depende de si existen los SP de fase gas.
            use_born_macro = True

        if use_born_macro:
            result["G_hplus_aq_kcal"] = G_HPLUS_AQ_BORN
            result["method"] = "born"
            if method == "auto":
                _warn(
                    f"method='auto' resolved to 'born' for the combined "
                    f"neutral+zwitterion macrostate "
                    f"({mol_acid if data_acid.get('is_macrostate') else mol_base}): "
                    f"no gas-phase SP is needed here (see derivation note "
                    f"above), so 'born' is used exactly, same as it would be "
                    f"for a single species with gas-phase SPs available."
                )
            else:
                _warn(
                    f"method='born' with a combined neutral+zwitterion macrostate "
                    f"({mol_acid if data_acid.get('is_macrostate') else mol_base}): "
                    f"no gas-phase SP is used (or needed) for the mixture. ΔG_aq "
                    f"is obtained via the algebraic reduction of the Born cycle "
                    f"to G_total_CPCM(base) − G_total_CPCM(acid) + G_HPLUS_AQ_BORN "
                    f"— the same reduction that already applies to a single "
                    f"species in this script, since the E_gas terms cancel "
                    f"identically (see the derivation note above "
                    f"compute_pka_montecarlo())."
                )
        else:
            result["G_hplus_aq_kcal"] = G_HPLUS_AQ_DIRECT
            result["method"] = "direct"
            if method == "direct":
                _warn(
                    f"method='direct' explicitly requested for the combined "
                    f"neutral+zwitterion macrostate "
                    f"({mol_acid if data_acid.get('is_macrostate') else mol_base})."
                )

        dG_aq_kcal = (
            (G_total_base - G_total_acid) * HARTREE_TO_KCAL
            + result["G_hplus_aq_kcal"]
        )
        result["dG_aq_kcal"] = dG_aq_kcal
        result["pKa"] = dG_aq_kcal / (R_KCAL * T * LN10)
        result["ok"] = True

        if verbose:
            if result["method"] == "born":
                _print_born_cycle_macrostate(result)
            else:
                _print_direct_cycle(result)

        return result

    # ── Obtener confórmero representativo de cada estado ──────────────────
    conf_acid = data_acid.get("best_conf")
    conf_base = data_base.get("best_conf")
    result["conf_acid"] = conf_acid
    result["conf_base"] = conf_base

    if conf_acid is None:
        _warn(f"No representative conformer for {mol_acid}.")
        return result
    if conf_base is None:
        _warn(f"No representative conformer for {mol_base}.")
        return result

    df_acid = data_acid["df"]
    df_base = data_base["df"]

    row_acid = df_acid[df_acid["conformer"] == conf_acid].iloc[0]
    row_base = df_base[df_base["conformer"] == conf_base].iloc[0]

    # ── Verificación de consistencia de nivel de teoría ────────────────────
    # Compara las líneas de keywords ORCA del SP acuoso y del freq.out entre
    # acid y base. Un desajuste aquí (distinta base, funcional, o grid) invalida
    # la resta de energías de ~1000+ Eh para sacar una diferencia de ~0.4 Eh,
    # independientemente de qué convención de G(H⁺) se use después.
    sp_acid_file  = WORK_ROOT / mol_acid / SUBDIR_SP / conf_acid / "sp.out"
    sp_base_file  = WORK_ROOT / mol_base / SUBDIR_SP / conf_base / "sp.out"
    freq_acid_lvl = parse_theory_level(WORK_ROOT / mol_acid / SUBDIR_FREQ / conf_acid / "freq.out")
    freq_base_lvl = parse_theory_level(WORK_ROOT / mol_base / SUBDIR_FREQ / conf_base / "freq.out")
    sp_acid_lvl   = parse_theory_level(sp_acid_file)
    sp_base_lvl   = parse_theory_level(sp_base_file)

    result["theory_sp_acid"]   = sp_acid_lvl
    result["theory_sp_base"]   = sp_base_lvl
    result["theory_freq_acid"] = freq_acid_lvl
    result["theory_freq_base"] = freq_base_lvl

    if sp_acid_lvl and sp_base_lvl and sp_acid_lvl != sp_base_lvl:
        _warn(
            f"INCONSISTENT THEORY LEVEL in the aqueous SP between "
            f"{mol_acid} ('{sp_acid_lvl}') and {mol_base} ('{sp_base_lvl}'). "
            f"The pKa cycle is NOT valid if the basis set/functional/CPCM "
            f"differ between states — check before reporting."
        )
    if freq_acid_lvl and freq_base_lvl and freq_acid_lvl != freq_base_lvl:
        _warn(
            f"INCONSISTENT THEORY LEVEL in NumFreq between "
            f"{mol_acid} ('{freq_acid_lvl}') and {mol_base} ('{freq_base_lvl}')."
        )

    # ── Energías SP acuosas (CPCM) ────────────────────────────────────────
    E_aq_acid = row_acid["E_sp_aq"]
    E_aq_base = row_base["E_sp_aq"]
    result["E_sp_aq_acid_Eh"] = E_aq_acid
    result["E_sp_aq_base_Eh"] = E_aq_base

    if E_aq_acid is None or pd.isna(E_aq_acid):
        _warn(f"No aqueous SP for {mol_acid}/{conf_acid}.")
        return result
    if E_aq_base is None or pd.isna(E_aq_base):
        _warn(f"No aqueous SP for {mol_base}/{conf_base}.")
        return result

    # ── Energías SP fase gas (pueden faltar) ───────────────────────────────
    E_gas_acid = row_acid["E_sp_gas"]
    E_gas_base = row_base["E_sp_gas"]
    result["E_sp_gas_acid_Eh"] = E_gas_acid
    result["E_sp_gas_base_Eh"] = E_gas_base

    have_gas = not (
        E_gas_acid is None or pd.isna(E_gas_acid) or
        E_gas_base is None or pd.isna(E_gas_base)
    )

    # ── Decidir método ──────────────────────────────────────────────────────
    if method == "born" and not have_gas:
        _warn(
            f"method='born' requested but gas-phase SPs are missing for "
            f"{mol_acid}/{conf_acid} and/or {mol_base}/{conf_base}. "
            f"Create runs/{mol_acid}/{SUBDIR_SP_GAS}/{conf_acid}/sp_gas.inp "
            f"and the equivalent for {mol_base} (same functional/basis, no CPCM)."
        )
        return result

    use_born = have_gas and method in ("auto", "born")

    if method == "auto" and not have_gas:
        _warn(
            f"No gas-phase SP for {mol_acid}/{conf_acid} and/or {mol_base}/{conf_base}. "
            f"Using the direct method (CPCM only) instead: "
            f"ΔG_aq = G_total_CPCM(base) − G_total_CPCM(acid) + G(H⁺,aq). "
            f"For the full Born cycle, create the gas-phase SPs in "
            f"{SUBDIR_SP_GAS}/conf_NNN/sp_gas.out."
        )

    # ── Correcciones térmicas G_corr a temperatura T (necesarias en ambos métodos) ─
    freq_acid = WORK_ROOT / mol_acid / SUBDIR_FREQ / conf_acid / "freq.out"
    freq_base = WORK_ROOT / mol_base / SUBDIR_FREQ / conf_base / "freq.out"
    tc_acid = get_thermochem(freq_acid, T=T, source=gcorr_source,
                             gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)
    tc_base = get_thermochem(freq_base, T=T, source=gcorr_source,
                             gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)

    G_corr_acid: Optional[float] = None
    G_corr_base: Optional[float] = None

    if tc_acid["ok"]:
        G_corr_acid = tc_acid["G_corr"]
    else:
        extra = f" ({tc_acid['warning']})" if tc_acid.get("warning") else ""
        _warn(
            f"G_corr ({gcorr_source}) not available for {mol_acid}/{conf_acid}. "
            f"Using G ≈ E_SP instead (no ZPE/entropy; pKa has higher uncertainty).{extra}"
        )
    if tc_base["ok"]:
        G_corr_base = tc_base["G_corr"]
    else:
        extra = f" ({tc_base['warning']})" if tc_base.get("warning") else ""
        _warn(
            f"G_corr ({gcorr_source}) not available for {mol_base}/{conf_base}. "
            f"Using G ≈ E_SP.{extra}"
        )

    result["G_corr_acid_Eh"] = G_corr_acid
    result["G_corr_base_Eh"] = G_corr_base

    if use_born:
        # ── G_total en fase gas ───────────────────────────────────────────
        # G_gas(X) = E_SP_gas(X) + G_corr(X)
        # Aproximación Ho & Coote: G_corr del cálculo CPCM es válido en gas.
        G_gas_acid = E_gas_acid + (G_corr_acid if G_corr_acid is not None else 0.0)
        G_gas_base = E_gas_base + (G_corr_base if G_corr_base is not None else 0.0)
        result["G_gas_acid_Eh"] = G_gas_acid
        result["G_gas_base_Eh"] = G_gas_base

        # ── ΔG_gas ─────────────────────────────────────────────────────────
        # ΔG_gas = G(A⁻, gas) + G(H⁺, gas) − G(AH, gas)
        dG_gas_kcal = (G_gas_base - G_gas_acid) * HARTREE_TO_KCAL + G_HPLUS_GAS
        result["dG_gas_kcal"] = dG_gas_kcal

        # ── ΔG_solv por especie ───────────────────────────────────────────
        # ΔG_solv(X) = E_SP(CPCM) − E_SP(gas) + δG_conc
        dG_solv_acid = (E_aq_acid - E_gas_acid) * HARTREE_TO_KCAL + DG_CONC_CORR
        dG_solv_base = (E_aq_base - E_gas_base) * HARTREE_TO_KCAL + DG_CONC_CORR
        result["dG_solv_acid_kcal"] = dG_solv_acid
        result["dG_solv_base_kcal"] = dG_solv_base

        # ── ΔG_aq (ciclo de Born completo) ────────────────────────────────
        # ΔG_aq = ΔG_gas + ΔG_solv(A⁻) + ΔG_solv(H⁺) − ΔG_solv(AH)
        dG_aq_kcal = (
            dG_gas_kcal
            + dG_solv_base
            + DG_SOLV_HPLUS
            - dG_solv_acid
        )
        result["dG_aq_kcal"] = dG_aq_kcal
        result["method"] = "born"

    else:
        # ── Método directo: solo energías CPCM ───────────────────────────
        # G_total_CPCM(X) = E_SP(CPCM) + G_corr(X)
        G_total_acid = E_aq_acid + (G_corr_acid if G_corr_acid is not None else 0.0)
        G_total_base = E_aq_base + (G_corr_base if G_corr_base is not None else 0.0)
        result["G_total_acid_Eh"] = G_total_acid
        result["G_total_base_Eh"] = G_total_base
        result["G_hplus_aq_kcal"] = G_HPLUS_AQ_DIRECT

        # ΔG_aq = G_total_CPCM(A⁻) − G_total_CPCM(AH) + G(H⁺, aq)
        dG_aq_kcal = (G_total_base - G_total_acid) * HARTREE_TO_KCAL + G_HPLUS_AQ_DIRECT
        result["dG_aq_kcal"] = dG_aq_kcal
        result["method"] = "direct"

    # ── pKa ──────────────────────────────────────────────────────────────
    pKa = result["dG_aq_kcal"] / (R_KCAL * T * LN10)
    result["pKa"] = pKa
    result["ok"]  = True

    if verbose:
        if result["method"] == "born":
            _print_born_cycle(result)
        else:
            _print_direct_cycle(result)

    return result


# ─────────────────────────────────────────────
#  MONTE CARLO — PROPAGACIÓN DE INCERTIDUMBRE EN pKa (FASE 1)
# ─────────────────────────────────────────────
#
# NOTA DE DERIVACIÓN (por qué el ruido se aplica solo a E_aq y G_corr, no a
# E_gas por separado, en el método "born"):
#
#   dG_gas       = (G_gas_base − G_gas_acid)·H2K + G_HPLUS_GAS
#   dG_solv_X    = (E_aq_X − E_gas_X)·H2K + DG_CONC_CORR
#   dG_aq(born)  = dG_gas + dG_solv_base + DG_SOLV_HPLUS − dG_solv_acid
#
#   Sustituyendo G_gas_X = E_gas_X + G_corr_X y expandiendo, TODOS los
#   términos en E_gas_acid/E_gas_base se cancelan algebraicamente (y también
#   DG_CONC_CORR, que aparece una vez con cada signo). Lo que queda es
#   idéntico en forma al método directo:
#
#   dG_aq(born) = (E_aq_base − E_aq_acid + G_corr_base − G_corr_acid)·H2K
#                 + (G_HPLUS_GAS + DG_SOLV_HPLUS)      # = G_HPLUS_AQ_BORN
#
#   dG_aq(direct) = (E_aq_base − E_aq_acid + G_corr_base − G_corr_acid)·H2K
#                   + G_HPLUS_AQ_DIRECT
#
#   O sea: dado cómo compute_pka() arma el ciclo de Born hoy, "born" y
#   "direct" tienen EXACTAMENTE la misma sensibilidad al ruido en E_aq/G_corr;
#   solo difieren en qué constante de G(H⁺,aq) se suma. Por eso el modelo de
#   ruido de abajo es el mismo para ambos métodos, y solo cambia el punto de
#   partida (base_result["dG_aq_kcal"]). Si algún día cambia la fórmula de
#   compute_pka() (p. ej. DG_CONC_CORR deja de cancelarse), esta nota y el
#   test test_born_and_direct_same_noise_sensitivity() dejan de ser válidos
#   y hay que revisar esta función.
#
# FUENTES DE INCERTIDUMBRE MODELADAS (ver plan, sección 1):
#   - sigma_elec_kcal  : energía electrónica SP (funcional/base), por estado,
#                        independiente entre ácido y base.
#   - sigma_gcorr_kcal : corrección térmica G_corr (ZPE+TS), por estado,
#                        independiente entre ácido y base.
#   - sigma_hplus_kcal : G(H⁺,aq) — constante compartida por TODA la cadena
#                        de pKa's. Por eso admite shared_hplus_draw: al
#                        correr varios pKa's de la misma cadena (pKa0, pKa1,
#                        pKa2) hay que pasarles el MISMO array de muestras
#                        para esta constante, o la incertidumbre de
#                        cantidades derivadas (Δ entre pKa's, fracciones de
#                        especie) queda sobreestimada.

def compute_pka_montecarlo(
    data_acid: dict,
    data_base: dict,
    T: float = T_STD,
    method: str = "auto",
    gcorr_source: str = "harmonic",
    n_samples: int = 10_000,
    sigma_elec_kcal: float = 0.5,
    sigma_gcorr_kcal: float = 0.4,
    sigma_hplus_kcal: float = 2.5,
    shared_hplus_draw: Optional[np.ndarray] = None,
    rng: Optional[np.random.Generator] = None,
) -> dict:
    """
    Propaga incertidumbre en el pKa de AH → A⁻ + H⁺ vía Monte Carlo,
    reusando compute_pka() como fuente única de verdad para el valor
    determinístico (mismo mapeo de estados, misma validación de carga/N°
    de átomos, mismo método born/direct/auto).

    No vuelve a leer archivos ni reimplementa el ciclo termodinámico: toma
    dG_aq_kcal ya calculado por compute_pka() como centro de la distribución
    y le suma perturbaciones gaussianas sobre los términos que compute_pka()
    ya expone en su dict de resultado (E_sp_aq_acid/base_Eh — si están
    disponibles — o, si no, sobre un sigma combinado, ver el caso de
    macroestados más abajo).

    Parameters
    ----------
    data_acid, data_base : dict
        Mismos dicts que recibiría compute_pka() (salida de
        collect_mol_data() o de compute_tautomer_equilibrium()).
    n_samples : int
        Número de muestras Monte Carlo.
    sigma_elec_kcal, sigma_gcorr_kcal : float
        Desviación estándar asumida para la energía electrónica SP y la
        corrección térmica, respectivamente, por estado (kcal/mol).
    sigma_hplus_kcal : float
        Desviación estándar asumida para G(H⁺,aq) — dispersión típica entre
        valores de literatura (Pliego-Riveros, Kelly-Cramer-Truhlar,
        Zhan-Dixon); ver plan, sección 1.
    shared_hplus_draw : np.ndarray, opcional
        Si se provee (shape (n_samples,)), se usa como la perturbación de
        G(H⁺,aq) en vez de samplear una nueva — para correlacionar el mismo
        error sistemático entre varios pKa's de una misma cadena (ver
        run_chain_montecarlo(), fase 2).
    rng : np.random.Generator, opcional
        Generador de números aleatorios; si no se provee, usa uno sin
        semilla fija (no reproducible entre corridas).

    Returns
    -------
    dict con:
        ok               : bool
        base_result      : dict — salida de compute_pka() (determinístico)
        method_used      : 'born' | 'direct' (resuelto, nunca 'auto')
        pka_samples      : np.ndarray, shape (n_samples,)
        pka_mean         : float
        pka_std          : float (ddof=1)
        pka_median       : float
        pka_ci95         : tuple (lo, hi) — percentiles 2.5/97.5
        hplus_draw       : np.ndarray — la perturbación de G(H⁺,aq) usada
                            (para reusar como shared_hplus_draw en otro pKa)
        warnings         : list[str]
    """
    base_result = compute_pka(
        data_acid, data_base, T=T, method=method,
        gcorr_source=gcorr_source, verbose=False,
    )

    out: dict = {
        "ok": False,
        "base_result": base_result,
        "method_used": base_result.get("method"),
        "n_samples": n_samples,
        "pka_samples": None,
        "pka_mean": None, "pka_std": None, "pka_median": None,
        "pka_ci95": None,
        "hplus_draw": None,
        "sigma_elec_kcal": sigma_elec_kcal,
        "sigma_gcorr_kcal": sigma_gcorr_kcal,
        "sigma_hplus_kcal": sigma_hplus_kcal,
        "warnings": list(base_result.get("warnings", [])),
    }

    if not base_result.get("ok"):
        out["warnings"].append(
            "compute_pka() determinístico falló (ok=False); no se puede "
            "construir una distribución Monte Carlo encima de un ciclo "
            "inválido. Revisa los warnings de base_result."
        )
        return out

    if n_samples < 2:
        out["warnings"].append(f"n_samples={n_samples} demasiado bajo; se requieren ≥2.")
        return out

    if shared_hplus_draw is not None and len(shared_hplus_draw) != n_samples:
        raise ValueError(
            f"shared_hplus_draw tiene {len(shared_hplus_draw)} muestras, "
            f"pero n_samples={n_samples}. Deben coincidir para poder sumarlo "
            f"término a término con el resto de las perturbaciones."
        )

    if rng is None:
        rng = np.random.default_rng()

    d_hplus = (
        shared_hplus_draw if shared_hplus_draw is not None
        else rng.normal(0.0, sigma_hplus_kcal, size=n_samples)
    )

    base_dG_aq = base_result["dG_aq_kcal"]

    # ── ¿Tenemos E_sp_aq por separado (caso normal) o solo G_total combinado
    #    (caso macroestado neutro+zwitterión, ver compute_tautomer_equilibrium)? ──
    have_separate_terms = (
        base_result.get("E_sp_aq_acid_Eh") is not None and
        base_result.get("E_sp_aq_base_Eh") is not None
    )

    if have_separate_terms:
        # Ruido independiente en energía electrónica y corrección térmica,
        # por estado (ácido y base son cálculos distintos → errores no
        # correlacionados entre sí). Ver nota de derivación arriba: esto
        # es válido tanto para method_used='born' como 'direct'.
        d_elec_acid  = rng.normal(0.0, sigma_elec_kcal,  size=n_samples)
        d_elec_base  = rng.normal(0.0, sigma_elec_kcal,  size=n_samples)
        d_gcorr_acid = rng.normal(0.0, sigma_gcorr_kcal, size=n_samples)
        d_gcorr_base = rng.normal(0.0, sigma_gcorr_kcal, size=n_samples)

        dG_aq_samples = (
            base_dG_aq
            + (d_elec_base - d_elec_acid)
            + (d_gcorr_base - d_gcorr_acid)
            + d_hplus
        )
    else:
        # Macroestado combinado (o, por alguna razón, sin desglose): no hay
        # forma de separar contribución electrónica de térmica dentro de
        # G_total_Eh ya mezclado con pesos de Boltzmann internos. Se usa un
        # sigma combinado en cuadratura como aproximación conservadora y se
        # deja registrado en warnings — esto es una simplificación de Fase 1,
        # no un tratamiento riguroso del error del ensemble.
        out["warnings"].append(
            "Sin desglose E_sp_aq/G_corr por separado (probablemente un "
            "macroestado neutro+zwitterión); se usa un sigma combinado "
            "sqrt(sigma_elec² + sigma_gcorr²) por estado en vez de separar "
            "las dos fuentes de error."
        )
        sigma_total = float(np.hypot(sigma_elec_kcal, sigma_gcorr_kcal))
        d_total_acid = rng.normal(0.0, sigma_total, size=n_samples)
        d_total_base = rng.normal(0.0, sigma_total, size=n_samples)

        dG_aq_samples = (
            base_dG_aq
            + (d_total_base - d_total_acid)
            + d_hplus
        )

    pka_samples = dG_aq_samples / (R_KCAL * T * LN10)
    ci_lo, ci_hi = np.percentile(pka_samples, [2.5, 97.5])

    out.update({
        "ok": True,
        "pka_samples": pka_samples,
        "pka_mean":   float(np.mean(pka_samples)),
        "pka_std":    float(np.std(pka_samples, ddof=1)),
        "pka_median": float(np.median(pka_samples)),
        "pka_ci95":   (float(ci_lo), float(ci_hi)),
        "hplus_draw": d_hplus,
    })
    return out


def _print_born_cycle(r: dict):
    """Imprime las contribuciones detalladas del ciclo de Born."""
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  BORN CYCLE — {r['mol_acid']} → {r['mol_base']} + H⁺")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. AH  : {r['conf_acid']}")
    print(f"  Conf. A⁻  : {r['conf_base']}")
    print(f"  SP level (AH)  : {r['theory_sp_acid'] or '—'}")
    print(f"  SP level (A⁻)  : {r['theory_sp_base'] or '—'}")
    print()
    print(f"  {'Term':<40} {'Value':>12}")
    print(f"  {'─'*40}  {'─'*12}")

    def fmt(v): return f"{v:+12.3f}" if v is not None else "       N/A"

    print(f"  {'G_gas(AH)  [Eh]':<40} {r['G_gas_acid_Eh']:>12.8f}")
    print(f"  {'G_gas(A⁻)  [Eh]':<40} {r['G_gas_base_Eh']:>12.8f}")
    print(f"  {'G(H⁺, gas)  [kcal/mol]':<40} {G_HPLUS_GAS:>12.3f}")
    print()
    print(f"  {'ΔG_gas  [kcal/mol]':<40} {fmt(r['dG_gas_kcal'])}")
    print(f"  {'ΔG_solv(AH)  [kcal/mol]':<40} {fmt(r['dG_solv_acid_kcal'])}")
    print(f"  {'ΔG_solv(A⁻)  [kcal/mol]':<40} {fmt(r['dG_solv_base_kcal'])}")
    print(f"  {'ΔG_solv(H⁺) [Kelly 2006]  [kcal/mol]':<40} {DG_SOLV_HPLUS:>12.3f}")
    print(f"  {'─'*40}  {'─'*12}")
    print(f"  {'ΔG_aq  [kcal/mol]':<40} {fmt(r['dG_aq_kcal'])}")
    print(f"  {'ΔG_aq  [kJ/mol]':<40} {r['dG_aq_kcal']*4.184:>12.3f}" if r['dG_aq_kcal'] else "")
    print()
    pka_str = f"{r['pKa']:.2f}" if r['pKa'] is not None else "N/A"
    print(f"  pKa = ΔG_aq / (RT·ln10)  =  {pka_str}")
    print(f"  (typical uncertainty ±1–2 units; Born/CPCM/TZVPD)")
    if r["warnings"]:
        print(f"\n  Warnings:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


def _print_direct_cycle(r: dict):
    """Imprime las contribuciones del método directo (solo CPCM, sin fase gas)."""
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  DIRECT METHOD (CPCM only) — {r['mol_acid']} → {r['mol_base']} + H⁺")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. AH  : {r['conf_acid']}")
    print(f"  Conf. A⁻  : {r['conf_base']}")
    print(f"  SP level (AH)  : {r['theory_sp_acid'] or '—'}")
    print(f"  SP level (A⁻)  : {r['theory_sp_base'] or '—'}")
    print()
    print(f"  {'G_total_CPCM(AH)  [Eh]':<40} {r['G_total_acid_Eh']:>12.8f}")
    print(f"  {'G_total_CPCM(A⁻)  [Eh]':<40} {r['G_total_base_Eh']:>12.8f}")
    print(f"  {'G(H⁺, aq)  [kcal/mol]':<40} {r['G_hplus_aq_kcal']:>12.3f}"
          f"   [Pliego & Riveros 2002, JPCA 106, 7434]")
    print(f"  {'─'*40}  {'─'*12}")
    dG = r["dG_aq_kcal"]
    print(f"  {'ΔG_aq  [kcal/mol]':<40} {dG:+12.3f}" if dG is not None else "")
    print(f"  {'ΔG_aq  [kJ/mol]':<40} {dG*4.184:>12.3f}" if dG is not None else "")
    print()
    pka_str = f"{r['pKa']:.2f}" if r['pKa'] is not None else "N/A"
    print(f"  pKa = ΔG_aq / (RT·ln10)  =  {pka_str}")
    print(f"  (direct method: valid if the solvation error of the theory")
    print(f"  level is consistent across species; typical uncertainty ±1–2")
    print(f"  units, similar to the full Born cycle)")
    if r["warnings"]:
        print(f"\n  Warnings:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


def _print_born_cycle_macrostate(r: dict):
    """
    Imprime el ciclo de Born cuando uno de los dos lados es un macroestado
    combinado neutro+zwitterión. No hay G_gas/ΔG_solv por especie que mostrar
    (no se calculan ni se necesitan — ver la nota en compute_pka()): en vez
    de eso, se muestra explícitamente la reducción algebraica a la resta de
    energías CPCM con la constante G_HPLUS_AQ_BORN.
    """
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  BORN CYCLE (macrostate, algebraic reduction) — "
          f"{r['mol_acid']} → {r['mol_base']} + H⁺")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. AH  : {r['conf_acid']}")
    print(f"  Conf. A⁻  : {r['conf_base']}")
    print(f"  SP level (AH)  : {r['theory_sp_acid'] or '—'}")
    print(f"  SP level (A⁻)  : {r['theory_sp_base'] or '—'}")
    print()
    print(f"  No explicit gas-phase SP is used for the mixture — the Born")
    print(f"  cycle reduces algebraically to the CPCM-only energy difference")
    print(f"  (same E_gas-cancellation already present for a single species")
    print(f"  in this script; see the derivation note above")
    print(f"  compute_pka_montecarlo()).")
    print()
    print(f"  {'G_total_CPCM(AH)  [Eh]':<40} {r['G_total_acid_Eh']:>12.8f}")
    print(f"  {'G_total_CPCM(A⁻)  [Eh]':<40} {r['G_total_base_Eh']:>12.8f}")
    print(f"  {'G(H⁺, aq)  [kcal/mol]':<40} {r['G_hplus_aq_kcal']:>12.3f}"
          f"   [Tissandier 1998 + Kelly 2006]")
    print(f"  {'─'*40}  {'─'*12}")
    dG = r["dG_aq_kcal"]
    print(f"  {'ΔG_aq  [kcal/mol]':<40} {dG:+12.3f}" if dG is not None else "")
    print(f"  {'ΔG_aq  [kJ/mol]':<40} {dG*4.184:>12.3f}" if dG is not None else "")
    print()
    pka_str = f"{r['pKa']:.2f}" if r['pKa'] is not None else "N/A"
    print(f"  pKa = ΔG_aq / (RT·ln10)  =  {pka_str}")
    print(f"  (typical uncertainty ±1–2 units; Born/CPCM/TZVPD)")
    if r["warnings"]:
        print(f"\n  Warnings:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


# ─────────────────────────────────────────────
#  DISTRIBUCIÓN DE BOLTZMANN ENTRE ESTADOS
# ─────────────────────────────────────────────

def boltzmann_protonation_states(
    pKa1: Optional[float],
    pKa2: Optional[float],
    pH_values: list[float],
    T: float = T_STD,
) -> pd.DataFrame:
    """
    Calcula la distribución de especies entre estados de protonación
    a múltiples valores de pH usando Henderson-Hasselbalch encadenado.

    Para un diácido con pKa1 (AH₂ → AH⁻) y pKa2 (AH⁻ → A²⁻):

      f(AH₂) = 1 / D
      f(AH⁻) = 10^(pH − pKa1) / D
      f(A²⁻) = 10^(2·pH − pKa1 − pKa2) / D
      D = 1 + 10^(pH − pKa1) + 10^(2·pH − pKa1 − pKa2)

    Si solo hay pKa1 (sin pKa2 / sin dianión):
      f(AH)  = 1 / (1 + 10^(pH − pKa1))
      f(A⁻)  = 10^(pH − pKa1) / (1 + 10^(pH − pKa1))

    Retorna DataFrame con columnas:
      pH, f_neutral, f_anion, f_dianion (si pKa2 disponible)
    """
    rows = []
    for pH in pH_values:
        if pKa1 is not None and pKa2 is not None:
            a0 = 1.0
            a1 = 10 ** (pH - pKa1)
            a2 = 10 ** (2 * pH - pKa1 - pKa2)
            D  = a0 + a1 + a2
            rows.append({
                "pH":        pH,
                "f_neutral": a0 / D,
                "f_anion":   a1 / D,
                "f_dianion": a2 / D,
            })
        elif pKa1 is not None:
            a0 = 1.0
            a1 = 10 ** (pH - pKa1)
            D  = a0 + a1
            rows.append({
                "pH":        pH,
                "f_neutral": a0 / D,
                "f_anion":   a1 / D,
            })
        else:
            rows.append({"pH": pH})

    return pd.DataFrame(rows)


def boltzmann_5state(
    pKa0: Optional[float],
    pKa1: Optional[float],
    pKa2: Optional[float],
    f_R: Optional[float],
    f_D: Optional[float],
    pH_values: list[float],
) -> pd.DataFrame:
    """
    Distribución de Henderson-Hasselbalch encadenada para la cadena completa:

        A⁺ (catión)  ⇌  [R (neutro) ⇌ D (zwitterión)]  ⇌  B⁻ (anión)  ⇌  C²⁻ (dianión)
              pKa0                  (f_R / f_D, sin pH)        pKa1          pKa2

    pKa0 = A → neutro/zwitterión + H⁺
    pKa1 = neutro/zwitterión → B + H⁺
    pKa2 = B → C + H⁺

    f_R, f_D: fracción de R y D DENTRO del macroestado neutro (constante con
    el pH, viene de compute_tautomer_equilibrium). Si no se proveen, se asume
    f_R = 1, f_D = 0 (todo el macroestado se reporta como "neutro").

    Cualquier pKa que sea None se trata como "ese estado no existe" (se omite
    su rama de la cadena, igual que boltzmann_protonation_states).

    Retorna DataFrame con columnas: pH, f_cation, f_R, f_D, f_neutral_macro,
    f_anion, f_dianion (las que apliquen).
    """
    if f_R is None or f_D is None:
        f_R, f_D = 1.0, 0.0

    def _pow10(exp):
        # clip para evitar overflow cuando un pKa está muy lejos del rango de pH
        # (con pKa's físicamente razonables esto nunca se activa)
        return 10.0 ** np.clip(exp, -300.0, 300.0)

    rows = []
    for pH in pH_values:
        r0 = _pow10(pKa0 - pH) if pKa0 is not None else 0.0
        r1 = _pow10(pH - pKa1) if pKa1 is not None else 0.0
        r2 = (_pow10(2 * pH - pKa1 - pKa2)
              if (pKa1 is not None and pKa2 is not None) else 0.0)

        D = 1.0 + r0 + r1 + r2
        f_neutral_macro = 1.0 / D
        row = {
            "pH":               pH,
            "f_neutral_macro":  f_neutral_macro,
            "f_R":              f_neutral_macro * f_R,
            "f_D":              f_neutral_macro * f_D,
        }
        if pKa0 is not None:
            row["f_cation"] = r0 / D
        if pKa1 is not None:
            row["f_anion"] = r1 / D
        if pKa1 is not None and pKa2 is not None:
            row["f_dianion"] = r2 / D
        rows.append(row)

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────
#  MONTE CARLO — CADENA COMPLETA DE 5 ESTADOS (FASE 2)
# ─────────────────────────────────────────────

def _boltzmann_5state_vectorized(
    pKa0_samples: Optional[np.ndarray],
    pKa1_samples: Optional[np.ndarray],
    pKa2_samples: Optional[np.ndarray],
    f_R: Optional[float],
    f_D: Optional[float],
    pH_values: list[float],
) -> pd.DataFrame:
    """
    Versión vectorizada de boltzmann_5state() sobre un eje adicional de
    muestras Monte Carlo. En vez de llamar boltzmann_5state() una vez por
    cada una de las N muestras (con N=10 000 y ~140 valores de pH serían
    ~1.4M evaluaciones en loop de Python puro), arma arrays 2D
    (n_samples, n_pH) con broadcasting de numpy y reduce a percentiles por
    columna de pH al final.

    Reproduce exactamente la misma fórmula que boltzmann_5state() (mismo
    r0/r1/r2/D) — si esa función cambia, esta hay que actualizarla en
    paralelo; no hay una única fuente de verdad aquí porque una versión
    escalar con loop sobre N=10⁴ muestras sería demasiado lenta para uso
    interactivo.

    Parameters
    ----------
    pKa0_samples, pKa1_samples, pKa2_samples : np.ndarray o None
        Arrays de muestras MC (mismo n_samples cada uno, generados con el
        mismo shared_hplus_draw si vienen de la misma cadena — ver
        run_chain_montecarlo()). None si ese pKa no existe en la cadena
        (mismo criterio que boltzmann_5state: se omite esa rama).
    f_R, f_D : float, opcional
        Igual que en boltzmann_5state(); default 1.0/0.0 si no se proveen.
    pH_values : list[float]

    Returns
    -------
    pd.DataFrame con columna 'pH' y, para cada fracción que aplique
    (f_neutral_macro, f_R, f_D, f_cation, f_anion, f_dianion), tres
    columnas: '{frac}_lo' (percentil 2.5), '{frac}_median', '{frac}_hi'
    (percentil 97.5). DataFrame vacío si no hay ninguna muestra disponible.
    """
    if f_R is None or f_D is None:
        f_R, f_D = 1.0, 0.0

    have_pka0 = pKa0_samples is not None
    have_pka1 = pKa1_samples is not None
    have_pka2 = have_pka1 and (pKa2_samples is not None)

    if not have_pka0 and not have_pka1:
        return pd.DataFrame()

    n_samples = (
        len(pKa0_samples) if have_pka0 else
        len(pKa1_samples)
    )
    for name, arr in [("pKa0", pKa0_samples), ("pKa1", pKa1_samples), ("pKa2", pKa2_samples)]:
        if arr is not None and len(arr) != n_samples:
            raise ValueError(
                f"{name}_samples tiene {len(arr)} muestras, se esperaban "
                f"{n_samples} (todas las ramas de la cadena deben compartir "
                f"el mismo n_samples)."
            )

    pH_arr = np.asarray(pH_values, dtype=float)             # shape (M,)

    def _pow10_2d(exp_2d):
        return 10.0 ** np.clip(exp_2d, -300.0, 300.0)

    zeros_2d = np.zeros((n_samples, len(pH_arr)))

    r0 = (_pow10_2d(pKa0_samples[:, None] - pH_arr[None, :])
          if have_pka0 else zeros_2d)
    r1 = (_pow10_2d(pH_arr[None, :] - pKa1_samples[:, None])
          if have_pka1 else zeros_2d)
    r2 = (_pow10_2d(2 * pH_arr[None, :] - pKa1_samples[:, None] - pKa2_samples[:, None])
          if have_pka2 else zeros_2d)

    D = 1.0 + r0 + r1 + r2
    f_neutral_macro = 1.0 / D
    f_R_arr = f_neutral_macro * f_R
    f_D_arr = f_neutral_macro * f_D

    def _summ(arr_2d):
        lo, med, hi = np.percentile(arr_2d, [2.5, 50.0, 97.5], axis=0)
        return lo, med, hi

    cols: dict = {"pH": pH_arr}

    for label, arr_2d in [
        ("f_neutral_macro", f_neutral_macro),
        ("f_R", f_R_arr),
        ("f_D", f_D_arr),
    ]:
        lo, med, hi = _summ(arr_2d)
        cols[f"{label}_lo"], cols[f"{label}_median"], cols[f"{label}_hi"] = lo, med, hi

    if have_pka0:
        lo, med, hi = _summ(r0 / D)
        cols["f_cation_lo"], cols["f_cation_median"], cols["f_cation_hi"] = lo, med, hi
    if have_pka1:
        lo, med, hi = _summ(r1 / D)
        cols["f_anion_lo"], cols["f_anion_median"], cols["f_anion_hi"] = lo, med, hi
    if have_pka2:
        lo, med, hi = _summ(r2 / D)
        cols["f_dianion_lo"], cols["f_dianion_median"], cols["f_dianion_hi"] = lo, med, hi

    return pd.DataFrame(cols)


def run_chain_montecarlo(
    data_A: Optional[dict],
    data_neutral_macro: Optional[dict],
    data_B: Optional[dict],
    data_C: Optional[dict],
    f_R: Optional[float],
    f_D: Optional[float],
    pH_values: list[float],
    T: float = T_STD,
    method: str = "auto",
    gcorr_source: str = "harmonic",
    n_samples: int = 10_000,
    sigma_elec_kcal: float = 0.5,
    sigma_gcorr_kcal: float = 0.4,
    sigma_hplus_kcal: float = 2.5,
    seed: Optional[int] = None,
) -> dict:
    """
    Propaga incertidumbre por Monte Carlo a través de TODA la cadena de 5
    estados (A⁺ ⇌ [R⇌D] ⇌ B⁻ ⇌ C²⁻), preservando la correlación física de
    G(H⁺,aq) entre pKa0/pKa1/pKa2, y devuelve bandas de confianza para las
    fracciones de especie vs. pH.

    Mismos pares acid/base que usa _run_extended() para la cadena
    determinística:
        pKa0 = compute_pka(data_A, data_neutral_macro, ...)
        pKa1 = compute_pka(data_neutral_macro, data_B,  ...)
        pKa2 = compute_pka(data_B, data_C, ...)

    Cualquier par con un miembro en None se omite (mismo criterio que
    _run_extended: pKa0 requiere data_A, pKa2 requiere data_B y data_C).

    Punto clave de diseño: la perturbación de G(H⁺,aq) (shared_hplus_draw)
    se genera UNA sola vez acá y se pasa idéntica a las tres llamadas de
    compute_pka_montecarlo. Esto hace que, en las N muestras generadas,
    pKa0/pKa1/pKa2 se muevan de forma correlacionada cuando la constante del
    protón se desplaza — que es lo físicamente correcto, porque es la MISMA
    constante en los tres pasos de desprotonación. El ruido de energía
    electrónica/térmica (sigma_elec_kcal, sigma_gcorr_kcal) en cambio usa
    generadores independientes por tramo (rng.spawn), porque A, R/D, B y C
    son cálculos ORCA distintos con errores no correlacionados entre sí.

    Returns
    -------
    dict con:
        ok             : bool
        mc_pka0/1/2    : dict o None — salida completa de
                         compute_pka_montecarlo() para cada tramo presente
        df_dist_bands  : pd.DataFrame — salida de
                         _boltzmann_5state_vectorized() (columnas
                         {frac}_lo/_median/_hi por cada fracción aplicable)
        warnings       : list[str]
    """
    rng = np.random.default_rng(seed)

    # Constante compartida por toda la cadena — un solo draw, reusado en
    # las tres llamadas de compute_pka_montecarlo de abajo.
    shared_hplus_draw = rng.normal(0.0, sigma_hplus_kcal, size=n_samples)

    # Sub-generadores independientes para el ruido de E_elec/G_corr de cada
    # tramo (no deben compartirse entre tramos: cada uno involucra cálculos
    # ORCA distintos).
    rng_pka0, rng_pka1, rng_pka2 = rng.spawn(3)

    out: dict = {
        "ok": False,
        "n_samples": n_samples,
        "mc_pka0": None, "mc_pka1": None, "mc_pka2": None,
        "df_dist_bands": None,
        "warnings": [],
    }

    mc_pka0 = mc_pka1 = mc_pka2 = None

    if data_A is not None and data_neutral_macro is not None:
        mc_pka0 = compute_pka_montecarlo(
            data_A, data_neutral_macro, T=T, method=method, gcorr_source=gcorr_source,
            n_samples=n_samples,
            sigma_elec_kcal=sigma_elec_kcal, sigma_gcorr_kcal=sigma_gcorr_kcal,
            sigma_hplus_kcal=sigma_hplus_kcal,
            shared_hplus_draw=shared_hplus_draw, rng=rng_pka0,
        )
        out["warnings"].extend(mc_pka0["warnings"])

    if data_neutral_macro is not None and data_B is not None:
        mc_pka1 = compute_pka_montecarlo(
            data_neutral_macro, data_B, T=T, method=method, gcorr_source=gcorr_source,
            n_samples=n_samples,
            sigma_elec_kcal=sigma_elec_kcal, sigma_gcorr_kcal=sigma_gcorr_kcal,
            sigma_hplus_kcal=sigma_hplus_kcal,
            shared_hplus_draw=shared_hplus_draw, rng=rng_pka1,
        )
        out["warnings"].extend(mc_pka1["warnings"])

    if data_B is not None and data_C is not None:
        mc_pka2 = compute_pka_montecarlo(
            data_B, data_C, T=T, method=method, gcorr_source=gcorr_source,
            n_samples=n_samples,
            sigma_elec_kcal=sigma_elec_kcal, sigma_gcorr_kcal=sigma_gcorr_kcal,
            sigma_hplus_kcal=sigma_hplus_kcal,
            shared_hplus_draw=shared_hplus_draw, rng=rng_pka2,
        )
        out["warnings"].extend(mc_pka2["warnings"])

    out["mc_pka0"], out["mc_pka1"], out["mc_pka2"] = mc_pka0, mc_pka1, mc_pka2

    if (mc_pka0 is None or not mc_pka0["ok"]) and (mc_pka1 is None or not mc_pka1["ok"]):
        out["warnings"].append(
            "No hay suficientes tramos válidos para construir la cadena MC "
            "(se necesita al menos A+neutro/zwitterión, o neutro/zwitterión+B, "
            "con compute_pka() determinístico exitoso)."
        )
        return out

    pKa0_samples = mc_pka0["pka_samples"] if (mc_pka0 and mc_pka0["ok"]) else None
    pKa1_samples = mc_pka1["pka_samples"] if (mc_pka1 and mc_pka1["ok"]) else None
    pKa2_samples = mc_pka2["pka_samples"] if (mc_pka2 and mc_pka2["ok"]) else None

    df_bands = _boltzmann_5state_vectorized(
        pKa0_samples, pKa1_samples, pKa2_samples, f_R, f_D, pH_values,
    )

    out["ok"] = True
    out["df_dist_bands"] = df_bands
    return out


def print_mc_summary(chain_mc: Optional[dict], pka_labels: list) -> None:
    """
    Imprime tabla comparativa pKa determinístico vs. Monte Carlo (media±std,
    CI95) — análoga a print_pka_summary(), pero para la salida de
    run_chain_montecarlo().

    pka_labels: lista de tuplas (key, resultado_dict_o_None, etiqueta_str),
    p. ej. [("pKa0", pka0_result, "pKa₀ (A→R/D)"), ("pKa1", pka1_result,
    "pKa₁ (R/D→B)"), ("pKa2", pka2_result, "pKa₂ (B→C)")] — key debe
    coincidir con mc_pka0/mc_pka1/mc_pka2 dentro de chain_mc.
    """
    sep = "═" * 74
    print(f"\n{sep}")
    print("  PROPAGACIÓN DE INCERTIDUMBRE — MONTE CARLO")
    print(sep)

    if chain_mc is None or not chain_mc.get("ok"):
        print("\n  ⚠  El Monte Carlo no pudo ejecutarse.")
        for w in (chain_mc or {}).get("warnings", []):
            print(f"    ⚠  {w}")
        print(f"{sep}\n")
        return

    mc_map = {
        "pKa0": chain_mc.get("mc_pka0"),
        "pKa1": chain_mc.get("mc_pka1"),
        "pKa2": chain_mc.get("mc_pka2"),
    }

    print(f"\n  N = {chain_mc['n_samples']:,} muestras  "
          f"(G(H⁺,aq) compartido entre pKa0/pKa1/pKa2 — ver nota abajo)")
    print(f"\n  {'':<22} {'determinístico':>14} {'media MC':>10} "
          f"{'± std':>8} {'CI 95%':>18}")
    print(f"  {'-'*22} {'-'*14} {'-'*10} {'-'*8} {'-'*18}")

    any_row = False
    for key, det, label in pka_labels:
        mc = mc_map.get(key)
        if det is None and mc is None:
            continue
        any_row = True
        det_v = f"{det['pKa']:.2f}" if (det and det.get("ok")) else "—"
        if mc is not None and mc.get("ok"):
            mean_v = f"{mc['pka_mean']:.2f}"
            std_v  = f"±{mc['pka_std']:.2f}"
            lo, hi = mc["pka_ci95"]
            ci_v   = f"[{lo:.2f}, {hi:.2f}]"
        else:
            mean_v = std_v = ci_v = "—"
        print(f"  {label:<22} {det_v:>14} {mean_v:>10} {std_v:>8} {ci_v:>18}")

    if not any_row:
        print("  (sin pKa's disponibles para Monte Carlo)")

    sample_mc = next((mc for mc in mc_map.values() if mc is not None and mc.get("ok")), None)
    if sample_mc is not None:
        print(f"\n  σ asumidos: E_elec=±{sample_mc['sigma_elec_kcal']:.2f}  "
              f"G_corr=±{sample_mc['sigma_gcorr_kcal']:.2f}  "
              f"G(H⁺,aq)=±{sample_mc['sigma_hplus_kcal']:.2f} kcal/mol")

    if chain_mc["warnings"]:
        print(f"\n  ⚠  Advertencias:")
        for w in chain_mc["warnings"]:
            print(f"     ⚠  {w}")

    print(f"\n  Nota: G(H⁺,aq) es la misma constante física en los tres pasos")
    print(f"  de desprotonación, así que su perturbación se comparte entre")
    print(f"  pKa0/pKa1/pKa2 en cada muestra MC (no se samplea independiente")
    print(f"  por tramo) — esto correlaciona el sesgo sistemático entre pKa's")
    print(f"  y evita inflar artificialmente la incertidumbre de cantidades")
    print(f"  derivadas como Δ(pKa2−pKa1) o la fracción de especie a un pH")
    print(f"  dado. Las bandas de confianza de la distribución de especies")
    print(f"  quedan en df_dist_bands (percentiles 2.5/50/97.5).")
    print(f"{sep}\n")


# ─────────────────────────────────────────────
#  IMPRESIÓN DE RESULTADOS
# ─────────────────────────────────────────────

def print_mol_summary(data: dict, verbose: bool = False):
    """Imprime resumen del ensemble conformacional de una molécula."""
    mol = data["mol_name"]
    df  = data["df"]
    sep = "─" * 74

    print(f"\n{sep}")
    print(f"  ESTADO: {mol}")
    print(sep)

    if df.empty:
        print("  ⚠  Sin datos.")
        return

    q, mult, n_at = data.get("charge"), data.get("mult"), data.get("n_atoms")
    q_str  = f"{q:+d}"  if q  is not None else "?"
    m_str  = f"{mult}"  if mult is not None else "?"
    n_str  = f"{n_at}"  if n_at is not None else "?"
    print(f"  Carga: {q_str}   Multiplicidad: {m_str}   N° átomos: {n_str}")

    # Tabla de confórmeros
    cols_base = ["conformer", "E_sp_aq", "thermo_ok", "n_imag"]
    cols_g    = [c for c in df.columns if c.startswith("G_total_")]
    cols_b    = [c for c in df.columns if c.startswith("boltz_")]
    cols_dg   = [c for c in df.columns if c.startswith("dG_")]

    print(f"\n  Confórmeros SP (acuoso, CPCM):")
    disp = df[cols_base + cols_g[:2] + cols_b[:2] + cols_dg[:2]].copy()
    print(disp.to_string(index=False, float_format=lambda x: f"{x:.5f}"))

    print(f"\n  Confórmero representativo (min G_total @ 298K): {data['best_conf']}")
    print(f"\n  G_Boltzmann del ensemble (energía libre efectiva, Hartree):")
    for lbl, G in data["G_boltz"].items():
        if G is not None:
            print(f"    {lbl}: {G:.8f} Eh")
        else:
            print(f"    {lbl}: N/A")

    for w in data["warnings"]:
        print(f"\n  ⚠  {w}")


def _print_single_pka(pka: dict, label: str):
    """Imprime un único resultado de pKa con su etiqueta (pKa₀, pKa₁, pKa₂, ...)."""
    if pka is None:
        return
    if pka["ok"]:
        print(f"\n  {label}  ({pka['mol_acid']} → {pka['mol_base']} + H⁺)"
              f"   [método: {pka['method']}]")
        if pka["method"] == "born" and pka.get("dG_gas_kcal") is not None:
            print(f"  ΔG_gas        = {pka['dG_gas_kcal']:+.3f} kcal/mol")
            print(f"  ΔG_solv(ácido) = {pka['dG_solv_acid_kcal']:+.3f} kcal/mol")
            print(f"  ΔG_solv(base)  = {pka['dG_solv_base_kcal']:+.3f} kcal/mol")
            print(f"  ΔG_solv(H⁺)   = {DG_SOLV_HPLUS:+.3f} kcal/mol  [Kelly 2006]")
        elif pka["method"] == "born":
            # Macroestado: ciclo de Born resuelto vía la reducción algebraica
            # (sin SP de fase gas para la mezcla) — ver compute_pka().
            print(f"  G(H⁺, aq)     = {pka['G_hplus_aq_kcal']:+.3f} kcal/mol "
                  f"[Tissandier 1998 + Kelly 2006; macroestado, reducción algebraica]")
        else:
            print(f"  G(H⁺, aq)     = {pka['G_hplus_aq_kcal']:+.3f} kcal/mol "
                  f"[Pliego & Riveros 2002]")
        print(f"  ΔG_aq         = {pka['dG_aq_kcal']:+.3f} kcal/mol")
        print(f"\n  ━━►  {label}  =  {pka['pKa']:.2f}  (T = {pka['T']:.2f} K)")
        if pka["warnings"]:
            print(f"\n  ⚠  Advertencias (revisar antes de reportar):")
            for w in pka["warnings"]:
                print(f"     ⚠  {w}")
    else:
        print(f"\n  {label}: NO calculado.")
        for w in pka["warnings"]:
            print(f"    ⚠  {w}")


def print_pka_summary(pka_list: list):
    """
    Imprime resumen de todos los pKa's calculados.

    pka_list: lista de tuplas (resultado_dict_o_None, etiqueta_str), en orden,
    p. ej. [(pka0, "pKa₀ (A→R/D)"), (pka1, "pKa₁ (R/D→B)"), (pka2, "pKa₂ (B→C)")]
    """
    sep = "═" * 74
    print(f"\n{sep}")
    print(f"  RESULTADOS FINALES — CICLOS DE pKa")
    print(sep)

    for pka, label in pka_list:
        if pka is not None:
            _print_single_pka(pka, label)

    print(f"\n  Nota: incertidumbre típica ±1–2 unidades de pKa para métodos")
    print(f"  Born/CPCM a nivel ωB97X-D3/def2-TZVPD (Ho & Coote 2010). El método")
    print(f"  'direct' usa solo energías CPCM (sin SP fase gas) con G(H⁺,aq) de")
    print(f"  Pliego & Riveros 2002 — esta constante NO es intercambiable con la")
    print(f"  de Tissandier/Kelly usada en 'born', que está calibrada para SP en")
    print(f"  fase gas pura. Cada método usa su propia referencia del protón.")
    print(f"  Los pasos que involucran el macroestado neutro/zwitterión combinado")
    print(f"  resuelven method='auto' a 'born' (exacto, sin SP fase gas — ver")
    print(f"  advertencias arriba). Con --method direct explícito usan 'direct'.")
    print(f"{sep}\n")


def _cross_validate_ph_scale(
    data_acid: dict,
    data_base: dict,
    result: dict,
    T: float = T_STD,
    gcorr_source: str = "harmonic",
) -> dict:
    """
    Diagnóstico de consistencia de escala de energía libre del protón.

    Problema que resuelve: con method='auto', pKa₁ (R/D→B, macroestado) y
    pKa₂ (B→C, simple) pueden resolver a métodos distintos. pKa₁ resuelve a
    'born' siempre que involucre el macroestado (exacto ahí, sin necesitar
    SP fase gas — ver rama is_macrostate en compute_pka), mientras que pKa₂
    depende de si existen SP fase gas para B y C. 'direct' y 'born' usan
    constantes G(H⁺,aq) distintas y NO intercambiables (Pliego & Riveros
    −264.61 kcal/mol vs. Tissandier/Kelly −272.18 kcal/mol) — sumar un pKa
    'direct' con un pKa 'born' en la misma cadena mezcla dos escalas de
    referencia del protón. Bajo method='auto' esto ahora se trata como error
    duro en el caller (ver chequeo tras _cross_validate_ph_scale más abajo),
    no solo como diagnóstico.

    Esta función NO reemplaza el resultado ya calculado (`result`); solo
    lo audita. Si `result['method'] == 'born'`, recalcula el mismo par
    ácido/base forzando method='direct' (reusando compute_pka, que ya
    implementa ambas rutas) y reporta el delta entre ambos pKa como
    diagnóstico. Si `result['method']` ya es 'direct', no hay nada que
    cross-validar (ambos pasos ya comparten la misma escala) y se retorna
    ok=False sin advertencia.
    """
    diag: dict = {
        "ok": False, "applicable": False,
        "pKa_born": None, "pKa_direct": None,
        "delta_pKa": None, "delta_dG_kcal": None,
        "warnings": [],
    }

    if result is None or not result.get("ok"):
        return diag

    if result.get("method") != "born":
        # Ya está en 'direct': mismo pie que pKa₁/pKa₀, nada que cruzar.
        return diag

    diag["applicable"] = True
    diag["pKa_born"] = result["pKa"]

    direct_result = compute_pka(
        data_acid, data_base, T=T, verbose=False,
        method="direct", gcorr_source=gcorr_source,
    )

    if not direct_result.get("ok"):
        diag["warnings"].append(
            "No se pudo recalcular con method='direct' para cross-validar "
            "contra el resultado 'born' (revisar warnings del sub-cálculo)."
        )
        return diag

    diag["ok"] = True
    diag["pKa_direct"] = direct_result["pKa"]
    diag["delta_pKa"] = result["pKa"] - direct_result["pKa"]
    diag["delta_dG_kcal"] = result["dG_aq_kcal"] - direct_result["dG_aq_kcal"]
    return diag


def _print_ph_scale_diagnostic(diag: dict, label: str = "pKa₂ (B→C)"):
    """Imprime el resultado de _cross_validate_ph_scale(), si aplica."""
    if not diag.get("applicable"):
        return

    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  DIAGNÓSTICO — consistencia de escala del protón ({label})")
    print(f"  {sep}")

    if not diag["ok"]:
        print(f"  ⚠  No se pudo completar el cross-check:")
        for w in diag["warnings"]:
            print(f"     {w}")
        print(f"  {sep}\n")
        return

    print(f"  {'method=\"born\"  -> pKa':<40} {diag['pKa_born']:>10.2f}")
    print(f"  {'method=\"direct\" (forzado) -> pKa':<40} {diag['pKa_direct']:>10.2f}")
    print(f"  {'Delta pKa (born - direct)':<40} {diag['delta_pKa']:>+10.2f}")
    print(f"  {'Delta dG_aq (born - direct)  [kcal/mol]':<40} {diag['delta_dG_kcal']:>+10.2f}")
    print()
    if abs(diag["delta_pKa"]) > 2.0:
        print(f"  \u26a0  Delta pKa > 2 unidades: born y direct discrepan mas de lo")
        print(f"     tipicamente esperable entre metodos. Si este resultado se")
        print(f"     va a combinar con pKa1 (que siempre usa 'direct'), evaluar")
        print(f"     reportar tambien la version 'direct' de este paso para que")
        print(f"     toda la cadena comparta la misma referencia de G(H+,aq).")
    else:
        print(f"  Delta dentro de lo esperable entre convenciones de G(H+,aq);")
        print(f"  no invalida el resultado, pero recordar que 'direct' es la")
        print(f"  escala que comparte el resto de la cadena (pKa0/pKa1).")
    print(f"  {sep}\n")


def _pka_temperature_scan(
    steps: list,
    gcorr_source: str = "harmonic",
    method: str = "auto",
) -> dict:
    """
    Escanea pKa en las 4 temperaturas de TEMPS_K/TEMP_LABELS (mismas T a las
    que collect_mol_data() ya tabula G_boltz por ensemble conformacional).

    Nota #8 del review: T estaba hardcodeada a T_STD=298.15 K en toda la
    orquestación (_run_classic/_run_extended), pese a que el pipeline ya
    calcula G a 4 temperaturas para el Boltzmann conformacional. Esta
    función reusa esa infraestructura vía compute_pka(T=...) para reportar
    cómo cambia cada pKa con T, sin fijar una sola T "de fábrica" — relevante
    para MAAs en agua de mar antártica (~274 K), no solo a 298.15 K estándar.

    steps: lista de tuplas (label, data_acid, data_base). Si data_acid o
    data_base es None, ese step se omite (igual que en el flujo principal).

    Retorna {label: {T_label: pka_result_dict}}.
    """
    scan: dict = {}
    for label, data_acid, data_base in steps:
        if data_acid is None or data_base is None:
            continue
        scan[label] = {}
        for T, lbl in zip(TEMPS_K, TEMP_LABELS):
            scan[label][lbl] = compute_pka(
                data_acid, data_base, T=T, verbose=False,
                method=method, gcorr_source=gcorr_source,
            )
    return scan


def _print_pka_temperature_scan(scan: dict, reported_T: float):
    """Imprime tabla pKa vs. T para cada step de _pka_temperature_scan()."""
    if not scan:
        return

    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  ESCANEO DE TEMPERATURA — pKa(T)")
    print(f"  {sep}")
    print(f"  {'Paso':<18}" + "".join(f"{lbl:>12}" for lbl in TEMP_LABELS))
    print(f"  {'─'*18}" + "".join(f"{'─'*12}" for _ in TEMP_LABELS))

    for label, by_T in scan.items():
        row = f"  {label:<18}"
        for T, lbl in zip(TEMPS_K, TEMP_LABELS):
            r = by_T.get(lbl)
            marker = "*" if np.isclose(T, reported_T, atol=1e-3) else " "
            if r and r.get("ok"):
                row += f"{r['pKa']:>11.2f}{marker}"
            else:
                row += f"{'N/A':>11}{marker}"
        print(row)

    print(f"\n  * = T reportada en el resto del output (--temp, default "
          f"{T_STD:.2f} K).")
    print(f"  Nota: el pKa incluye T tanto en G_corr (vía GoodVibes/NumFreq)")
    print(f"  como en RT·ln10 al final del ciclo; el efecto no es solo el")
    print(f"  denominador. No confundir con --goodvibes (harmonic vs.")
    print(f"  quasi-harmonic), que es un eje ortogonal a esta T.")
    print(f"  {sep}\n")


def _extended_temperature_scan(
    data_A: Optional[dict], data_R: dict, data_D: dict,
    data_B: Optional[dict], data_C: Optional[dict],
    gcorr_source: str = "harmonic",
    method: str = "auto",
) -> dict:
    """
    Análogo a _pka_temperature_scan() pero para la cadena extendida de 5
    estados, donde pKa₀/pKa₁ dependen del macroestado neutro/zwitterión
    (data_neutral_macro), que a su vez está "horneado" a una única T dentro
    de taut_result["G_macro_Eh"]. No alcanza con pasarle T distinta a
    compute_pka() para esos dos pasos — hay que reconstruir el equilibrio
    de tautómeros y el macroestado en cada T del escaneo.

    Retorna {label: {T_label: pka_result_dict}}, mismo formato que
    _pka_temperature_scan(), más "f_D" (fracción de zwitterión por T) para
    contexto, ya que también depende de T.
    """
    scan: dict = {"f_D": {}}

    for T, lbl in zip(TEMPS_K, TEMP_LABELS):
        taut_T = compute_tautomer_equilibrium(data_R, data_D, T=T, verbose=False,
                                              gcorr_source=gcorr_source)
        if not taut_T["ok"]:
            gR_T = _state_G_total_direct(data_R, T=T, gcorr_source=gcorr_source)
            if not gR_T.get("ok"):
                scan["f_D"][lbl] = None
                continue
            taut_T = {
                "f_R": 1.0, "f_D": 0.0,
                "G_macro_Eh": gR_T["G_total_Eh"],
                "theory_R": gR_T["theory_sp"], "theory_D": None,
            }

        scan["f_D"][lbl] = taut_T["f_D"]

        macro_T = {
            "mol_name":          "neutro/zwitterión",
            "is_macrostate":     True,
            "G_total_direct_Eh": taut_T["G_macro_Eh"],
            "conf_label":        f"macroestado @ {lbl}",
            "theory_label":      taut_T["theory_R"] or taut_T["theory_D"],
            "charge":            data_R.get("charge"),
            "n_atoms":           data_R.get("n_atoms"),
        }

        if data_A is not None:
            scan.setdefault("pKa0", {})[lbl] = compute_pka(
                data_A, macro_T, T=T, verbose=False,
                method=method, gcorr_source=gcorr_source,
            )
        if data_B is not None:
            scan.setdefault("pKa1", {})[lbl] = compute_pka(
                macro_T, data_B, T=T, verbose=False,
                method=method, gcorr_source=gcorr_source,
            )
        if data_B is not None and data_C is not None:
            scan.setdefault("pKa2", {})[lbl] = compute_pka(
                data_B, data_C, T=T, verbose=False,
                method=method, gcorr_source=gcorr_source,
            )

    return scan


def _print_extended_temperature_scan(scan: dict, reported_T: float):
    """Imprime tabla pKa(T) + f_D(T) para la cadena extendida de 5 estados."""
    pka_labels = [k for k in ("pKa0", "pKa1", "pKa2") if k in scan]
    if not pka_labels:
        return

    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  ESCANEO DE TEMPERATURA — pKa(T) y f_D(T)  (cadena extendida)")
    print(f"  {sep}")
    print(f"  {'Paso':<18}" + "".join(f"{lbl:>12}" for lbl in TEMP_LABELS))
    print(f"  {'─'*18}" + "".join(f"{'─'*12}" for _ in TEMP_LABELS))

    row = f"  {'f_D (zwitt.)':<18}"
    for lbl in TEMP_LABELS:
        fD = scan["f_D"].get(lbl)
        row += f"{fD*100:>10.1f}% " if fD is not None else f"{'N/A':>11} "
    print(row)

    for label in pka_labels:
        row = f"  {label:<18}"
        for T, lbl in zip(TEMPS_K, TEMP_LABELS):
            r = scan[label].get(lbl)
            marker = "*" if np.isclose(T, reported_T, atol=1e-3) else " "
            if r and r.get("ok"):
                row += f"{r['pKa']:>11.2f}{marker}"
            else:
                row += f"{'N/A':>11}{marker}"
        print(row)

    print(f"\n  * = T reportada en el resto del output (--temp, default "
          f"{T_STD:.2f} K).")
    print(f"  f_D depende de T (equilibrio R⇌D) y se re-cierra el macroestado")
    print(f"  neutro/zwitterión en cada T — no es solo re-evaluar RT·ln10 con")
    print(f"  un G_macro fijo. pKa₂ (B→C) no usa el macroestado y puede diferir")
    print(f"  en método ('born' vs 'direct') de pKa₀/pKa₁; ver diagnóstico de")
    print(f"  cross-validación de escala del protón más arriba si aplica.")
    print(f"  {sep}\n")


def print_distribution_table(
    df_dist: pd.DataFrame,
    pKa1: Optional[float],
    pKa2: Optional[float],
    ph_highlight: Optional[float] = None,
):
    """Imprime tabla de distribución de especies por pH."""
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  DISTRIBUCIÓN DE ESTADOS DE PROTONACIÓN  "
          f"(pKa₁ = {pKa1:.2f}" + (f", pKa₂ = {pKa2:.2f}" if pKa2 else "") + ")")
    print(sep)

    has_dianion = "f_dianion" in df_dist.columns

    # Cabecera
    if has_dianion:
        print(f"\n  {'pH':>6}  {'f(AH₂)':>10}  {'f(AH⁻)':>10}  {'f(A²⁻)':>10}  {'Dominante'}")
        print(f"  {'─'*6}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*12}")
    else:
        print(f"\n  {'pH':>6}  {'f(AH)':>10}  {'f(A⁻)':>10}  {'Dominante'}")
        print(f"  {'─'*6}  {'─'*10}  {'─'*10}  {'─'*12}")

    for _, row in df_dist.iterrows():
        pH   = row["pH"]
        fn   = row.get("f_neutral", 0.0)
        fa   = row.get("f_anion",   0.0)
        fd   = row.get("f_dianion", 0.0) if has_dianion else None

        fracs = {"neutral": fn, "anion": fa}
        if has_dianion:
            fracs["dianion"] = fd
        dominant = max(fracs, key=fracs.get)

        marker = " ◀" if (ph_highlight is not None and abs(pH - ph_highlight) < 0.05) else ""

        if has_dianion:
            print(f"  {pH:6.2f}  {fn:10.4f}  {fa:10.4f}  {fd:10.4f}  "
                  f"{dominant}{marker}")
        else:
            print(f"  {pH:6.2f}  {fn:10.4f}  {fa:10.4f}  {dominant}{marker}")

    # Puntos de cambio de especie dominante
    print(f"\n  Punto isoelectrónico (dominancia igual):")
    if pKa1:
        print(f"    AH₂ = AH⁻  →  pH = pKa₁ = {pKa1:.2f}")
    if pKa2:
        print(f"    AH⁻ = A²⁻  →  pH = pKa₂ = {pKa2:.2f}")
    if pKa1 and pKa2:
        pH_PI = (pKa1 + pKa2) / 2
        print(f"    Punto isoeléctrico (pI) ≈ (pKa₁ + pKa₂)/2 = {pH_PI:.2f}")

    print(sep)


def print_distribution_table_5(
    df_dist: pd.DataFrame,
    pKa0: Optional[float],
    pKa1: Optional[float],
    pKa2: Optional[float],
    ph_highlight: Optional[float] = None,
):
    """Imprime tabla de distribución de 5 especies (catión, R, D, anión, dianión)."""
    sep = "─" * 92
    print(f"\n{sep}")
    hdr = "  DISTRIBUCIÓN DE ESTADOS DE PROTONACIÓN  ("
    parts = []
    if pKa0 is not None: parts.append(f"pKa₀ = {pKa0:.2f}")
    if pKa1 is not None: parts.append(f"pKa₁ = {pKa1:.2f}")
    if pKa2 is not None: parts.append(f"pKa₂ = {pKa2:.2f}")
    print(hdr + ", ".join(parts) + ")")
    print(sep)

    has_cation  = "f_cation"  in df_dist.columns
    has_anion   = "f_anion"   in df_dist.columns
    has_dianion = "f_dianion" in df_dist.columns

    cols   = ["pH"]
    labels = [f"{'pH':>6}"]
    if has_cation:  labels.append(f"{'f(A⁺)':>9}")
    labels.append(f"{'f(R)':>9}")
    labels.append(f"{'f(D)':>9}")
    if has_anion:   labels.append(f"{'f(B⁻)':>9}")
    if has_dianion: labels.append(f"{'f(C²⁻)':>9}")
    labels.append("Dominante")
    print(f"\n  {'  '.join(labels)}")
    print(f"  {'─'*(len(''.join(labels)) + 2*(len(labels)-1))}")

    for _, row in df_dist.iterrows():
        pH = row["pH"]
        fracs = {"R": row.get("f_R", 0.0), "D": row.get("f_D", 0.0)}
        vals = [f"{pH:6.2f}"]
        if has_cation:
            fracs["A⁺"] = row.get("f_cation", 0.0)
            vals.append(f"{row.get('f_cation', 0.0):9.4f}")
        vals.append(f"{row.get('f_R', 0.0):9.4f}")
        vals.append(f"{row.get('f_D', 0.0):9.4f}")
        if has_anion:
            fracs["B⁻"] = row.get("f_anion", 0.0)
            vals.append(f"{row.get('f_anion', 0.0):9.4f}")
        if has_dianion:
            fracs["C²⁻"] = row.get("f_dianion", 0.0)
            vals.append(f"{row.get('f_dianion', 0.0):9.4f}")
        dominant = max(fracs, key=fracs.get)
        marker = " ◀" if (ph_highlight is not None and abs(pH - ph_highlight) < 0.05) else ""
        print(f"  {'  '.join(vals)}  {dominant}{marker}")

    print(f"\n  Nota: f(R) y f(D) son las fracciones DENTRO del macroestado neutro")
    print(f"  (equilibrio de tautómeros, independiente del pH); f(R)+f(D) = f(neutro/zwitterión total).")
    if pKa0 is not None and pKa1 is not None:
        pH_PI = (pKa0 + pKa1) / 2
        print(f"\n  pI aproximado ≈ (pKa₀ + pKa₁)/2 = {pH_PI:.2f}")
    print(sep)


# ─────────────────────────────────────────────
#  VISUALIZACIÓN (matplotlib)
# ─────────────────────────────────────────────

def plot_distribution(
    df_dist: pd.DataFrame,
    pKa1: Optional[float],
    pKa2: Optional[float],
    output_path: Optional[Path] = None,
    ph_highlight: Optional[float] = None,
    chain_mc: Optional[dict] = None,
):
    """
    Diagrama de distribución de especies (fracciones molares vs pH).
    Guarda como PNG si se provee output_path; si no, muestra en pantalla.

    chain_mc : dict de run_chain_montecarlo(), opcional (--montecarlo). En
        _run_classic() se construye con data_A=None, f_R=1.0, f_D=0.0, así
        que la columna 'f_neutral_macro' de df_dist_bands es equivalente a
        'f_neutral' acá (f_R=1 → f_R_macro == f_neutral_macro).
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: conda install matplotlib")
        return

    has_dianion = "f_dianion" in df_dist.columns
    pH_vals = df_dist["pH"].values

    df_bands = None
    if chain_mc is not None and chain_mc.get("ok"):
        cand = chain_mc.get("df_dist_bands")
        if cand is not None and not cand.empty and len(cand) == len(pH_vals):
            df_bands = cand
        elif cand is not None and not cand.empty:
            print("  ⚠  df_dist_bands tiene distinto largo que df_dist "
                  "(pH_all no coincide) — se omiten las bandas MC en el gráfico.")

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = ["#2196F3", "#4CAF50", "#F44336"]  # azul, verde, rojo
    labels = ["AH₂ (q=0)", "AH⁻ (q=−1)", "A²⁻ (q=−2)"]
    cols   = ["f_neutral", "f_anion"] + (["f_dianion"] if has_dianion else [])
    # mapeo a columnas de df_dist_bands (run_chain_montecarlo siempre usa
    # nomenclatura de 5 estados internamente)
    band_cols = {"f_neutral": "f_neutral_macro", "f_anion": "f_anion", "f_dianion": "f_dianion"}

    for col, label, color in zip(cols, labels, colors):
        if col in df_dist.columns:
            bcol = band_cols.get(col)
            if df_bands is not None and bcol and f"{bcol}_lo" in df_bands.columns and f"{bcol}_hi" in df_bands.columns:
                ax.fill_between(pH_vals,
                                 df_bands[f"{bcol}_lo"].values * 100,
                                 df_bands[f"{bcol}_hi"].values * 100,
                                 color=color, alpha=0.15, linewidth=0, zorder=1)
            ax.plot(pH_vals, df_dist[col].values * 100, label=label,
                    color=color, linewidth=2.0, zorder=2)

    if df_bands is not None:
        ax.plot([], [], color="gray", alpha=0.3, linewidth=8,
                label="banda MC 95% (percentil 2.5–97.5)")

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
    chain_mc: Optional[dict] = None,
):
    """
    Diagrama de distribución de 5 especies (catión, R, D, anión, dianión) vs pH.

    chain_mc : dict de run_chain_montecarlo(), opcional (--montecarlo). Si se
        provee y trae df_dist_bands no vacío, se dibuja una banda sombreada
        de percentil 2.5–97.5 alrededor de cada curva (misma columna de
        pH que df_dist, ya que ambas se generan con el mismo pH_all — ver
        _run_extended()/_run_classic()).
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: conda install matplotlib")
        return

    pH_vals = df_dist["pH"].values

    df_bands = None
    if chain_mc is not None and chain_mc.get("ok"):
        cand = chain_mc.get("df_dist_bands")
        if cand is not None and not cand.empty and len(cand) == len(pH_vals):
            df_bands = cand
        elif cand is not None and not cand.empty:
            print("  ⚠  df_dist_bands tiene distinto largo que df_dist "
                  "(pH_all no coincide) — se omiten las bandas MC en el gráfico.")

    fig, ax = plt.subplots(figsize=(9, 5.5))

    series = []
    if "f_cation" in df_dist.columns:
        series.append(("f_cation", "A⁺ catión (q=+1)", "#9C27B0", "-"))
    series.append(("f_R", "R neutro (q=0)",       "#2196F3", "-"))
    series.append(("f_D", "D zwitterión (q=0)",   "#00BCD4", "--"))
    if "f_anion" in df_dist.columns:
        series.append(("f_anion", "B⁻ anión (q=−1)",   "#4CAF50", "-"))
    if "f_dianion" in df_dist.columns:
        series.append(("f_dianion", "C²⁻ dianión (q=−2)", "#F44336", "-"))

    for col, label, color, ls in series:
        if col in df_dist.columns:
            if df_bands is not None and f"{col}_lo" in df_bands.columns and f"{col}_hi" in df_bands.columns:
                ax.fill_between(pH_vals,
                                 df_bands[f"{col}_lo"].values * 100,
                                 df_bands[f"{col}_hi"].values * 100,
                                 color=color, alpha=0.15, linewidth=0, zorder=1)
            ax.plot(pH_vals, df_dist[col].values * 100, label=label,
                    color=color, linewidth=2.0, linestyle=ls, zorder=2)

    if df_bands is not None:
        ax.plot([], [], color="gray", alpha=0.3, linewidth=8,
                label="banda MC 95% (percentil 2.5–97.5)")

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


def plot_pka_histograms(chain_mc: Optional[dict], pka_labels: list,
                        output_path: Optional[Path] = None):
    """
    Histograma de las muestras Monte Carlo de cada pKa presente en la
    cadena (chain_mc['mc_pka0']/['mc_pka1']/['mc_pka2']), con el pKa
    determinístico (compute_pka(), línea vertical) y la media/CI95 de la
    distribución MC superpuestos. Complementa las bandas de confianza de
    plot_distribution_5()/plot_distribution() con la vista "por pKa" que
    pide la Fase 4 del plan (histograma de cada pKa individual).

    pka_labels : misma estructura que print_mc_summary() — lista de
        tuplas (key, base_result, label) con key en {'pKa0','pKa1','pKa2'}
        y base_result el resultado determinístico de compute_pka() (para
        la línea vertical de referencia).
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: conda install matplotlib")
        return

    if chain_mc is None or not chain_mc.get("ok"):
        return

    key_to_mc = {"pKa0": "mc_pka0", "pKa1": "mc_pka1", "pKa2": "mc_pka2"}
    panels = []
    for key, base_result, label in pka_labels:
        mc = chain_mc.get(key_to_mc.get(key))
        if mc is not None and mc.get("ok") and mc.get("pka_samples") is not None:
            panels.append((key, base_result, label, mc))

    if not panels:
        return

    fig, axes = plt.subplots(1, len(panels), figsize=(4.5 * len(panels), 4),
                             sharey=False)
    if len(panels) == 1:
        axes = [axes]

    for ax, (key, base_result, label, mc) in zip(axes, panels):
        samples = mc["pka_samples"]
        ax.hist(samples, bins=60, color="#2196F3", alpha=0.75, edgecolor="none")

        pka_det = base_result["pKa"] if (base_result and base_result.get("ok")) else None
        if pka_det is not None:
            ax.axvline(pka_det, color="black", linestyle="-", linewidth=1.5,
                       label=f"determinístico = {pka_det:.2f}")

        ax.axvline(mc["pka_mean"], color="#F44336", linestyle="--", linewidth=1.5,
                   label=f"media MC = {mc['pka_mean']:.2f}")
        lo95, hi95 = mc["pka_ci95"]
        ax.axvspan(lo95, hi95, color="#F44336", alpha=0.08,
                   label=f"CI95 = [{lo95:.2f}, {hi95:.2f}]")

        ax.set_title(label, fontsize=10)
        ax.set_xlabel("pKa")
        ax.set_ylabel("Frecuencia" if ax == axes[0] else "")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, alpha=0.25)

    fig.suptitle(f"Distribución Monte Carlo de pKa (N={chain_mc['n_samples']:,} muestras)",
                fontsize=11)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Histogramas de pKa (MC) guardados en: {output_path}")
    else:
        plt.show()
    plt.close(fig)


# ─────────────────────────────────────────────
#  EXPORTACIÓN
# ─────────────────────────────────────────────

def export_results(
    pka_list: list,
    df_dist: pd.DataFrame,
    data_list: list[dict],
    fmt: str,
    taut_result: Optional[dict] = None,
    chain_mc: Optional[dict] = None,
):
    """
    Exporta resultados del ciclo de pKa y distribución de especies.

    pka_list: lista de tuplas (resultado_dict_o_None, etiqueta_str), p. ej.
      [(pka0, "pKa0"), (pka1, "pKa1"), (pka2, "pKa2")]
    taut_result: dict de compute_tautomer_equilibrium(), opcional (se agrega
      como fila/hoja adicional si se provee).
    chain_mc: dict de run_chain_montecarlo(), opcional (--montecarlo). Si se
      provee y ok=True, agrega las muestras crudas de pKa0/pKa1/pKa2
      (mc_samples.csv / hoja MC_samples) y las bandas de confianza de la
      distribución de especies (mc_distribution_bands.csv / hoja
      MC_distribution_bands) — útil para reproducir el análisis de
      incertidumbre sin recorrer todo el pipeline.
    """
    out_dir = WORK_ROOT
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── DataFrame del ciclo de Born / directo ────────────────────────────
    cycle_rows = []
    for pka_r, label in pka_list:
        if pka_r is None:
            continue
        cycle_rows.append({
            "reaction":          label,
            "mol_acid":          pka_r.get("mol_acid"),
            "mol_base":          pka_r.get("mol_base"),
            "method":            pka_r.get("method"),
            "conf_acid":         pka_r.get("conf_acid"),
            "conf_base":         pka_r.get("conf_base"),
            "T_K":               pka_r.get("T"),
            "dG_gas_kcal":       pka_r.get("dG_gas_kcal"),
            "dG_solv_acid_kcal": pka_r.get("dG_solv_acid_kcal"),
            "dG_solv_base_kcal": pka_r.get("dG_solv_base_kcal"),
            "dG_solv_Hplus_kcal": DG_SOLV_HPLUS if pka_r.get("method") == "born" else None,
            "G_hplus_aq_kcal":   pka_r.get("G_hplus_aq_kcal"),
            "dG_aq_kcal":        pka_r.get("dG_aq_kcal"),
            "dG_aq_kJ":          (pka_r.get("dG_aq_kcal") or 0) * 4.184,
            "pKa":               pka_r.get("pKa"),
            "ok":                pka_r.get("ok"),
        })

    df_cycle = pd.DataFrame(cycle_rows)

    # ── DataFrame del equilibrio de tautómeros (opcional) ────────────────
    df_taut = None
    if taut_result is not None and taut_result.get("ok"):
        df_taut = pd.DataFrame([{
            "mol_R":         taut_result.get("mol_R"),
            "mol_D":         taut_result.get("mol_D"),
            "conf_R":        taut_result.get("conf_R"),
            "conf_D":        taut_result.get("conf_D"),
            "T_K":           taut_result.get("T"),
            "G_R_Eh":        taut_result.get("G_R_Eh"),
            "G_D_Eh":        taut_result.get("G_D_Eh"),
            "dG_taut_kcal":  taut_result.get("dG_taut_kcal"),
            "K_taut":        taut_result.get("K_taut"),
            "f_R":           taut_result.get("f_R"),
            "f_D":           taut_result.get("f_D"),
            "G_macro_Eh":    taut_result.get("G_macro_Eh"),
        }])

    # ── DataFrame de distribución conformacional ─────────────────────────
    conf_rows = []
    for d in data_list:
        df_c = d["df"]
        for _, row in df_c.iterrows():
            r = {"mol": d["mol_name"]}
            r.update(row.to_dict())
            conf_rows.append(r)
    df_confs = pd.DataFrame(conf_rows)

    # ── DataFrame de muestras crudas Monte Carlo (opcional) ───────────────
    df_mc_samples = None
    df_mc_bands = None
    if chain_mc is not None and chain_mc.get("ok"):
        mc_cols = {}
        for key, colname in [("mc_pka0", "pKa0"), ("mc_pka1", "pKa1"), ("mc_pka2", "pKa2")]:
            mc = chain_mc.get(key)
            if mc is not None and mc.get("ok"):
                mc_cols[colname] = mc["pka_samples"]
        if mc_cols:
            df_mc_samples = pd.DataFrame(mc_cols)
        if chain_mc.get("df_dist_bands") is not None and not chain_mc["df_dist_bands"].empty:
            df_mc_bands = chain_mc["df_dist_bands"]

    if fmt in ("csv", "both"):
        df_cycle.to_csv(out_dir / "pka_cycle.csv", index=False)
        df_dist.to_csv(out_dir / "protonation_distribution.csv", index=False)
        df_confs.to_csv(out_dir / "conformer_boltzmann.csv", index=False)
        if df_taut is not None:
            df_taut.to_csv(out_dir / "tautomer_equilibrium.csv", index=False)
        if df_mc_samples is not None:
            df_mc_samples.to_csv(out_dir / "mc_samples.csv", index=False)
        if df_mc_bands is not None:
            df_mc_bands.to_csv(out_dir / "mc_distribution_bands.csv", index=False)
        print(f"  CSVs guardados en: {out_dir}/")

    if fmt in ("xlsx", "both"):
        xlsx_path = out_dir / "protonation_cycle_results.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            if not df_cycle.empty:
                df_cycle.to_excel(writer, sheet_name="pKa_cycle", index=False)
            if df_taut is not None:
                df_taut.to_excel(writer, sheet_name="Tautomer_equilibrium", index=False)
            df_dist.to_excel(writer, sheet_name="Species_distribution", index=False)
            if not df_confs.empty:
                # columnas relevantes
                keep = ["mol", "conformer", "E_sp_aq", "thermo_ok", "n_imag",
                        "G_total_298K", "boltz_298K", "dG_298K_kcal"]
                keep = [c for c in keep if c in df_confs.columns]
                df_confs[keep].to_excel(writer, sheet_name="Conformer_Boltzmann", index=False)
            if df_mc_samples is not None:
                df_mc_samples.to_excel(writer, sheet_name="MC_samples", index=False)
            if df_mc_bands is not None:
                df_mc_bands.to_excel(writer, sheet_name="MC_distribution_bands", index=False)
        print(f"  Excel guardado en: {xlsx_path}")


# ─────────────────────────────────────────────
#  PUNTO DE ENTRADA
# ─────────────────────────────────────────────

def main():
    # Debe ir ANTES de cualquier lectura de estas variables en esta función
    # (p. ej. como default= en add_argument más abajo). Python exige que la
    # declaración `global` preceda todo uso del nombre dentro del mismo
    # scope, o produce SyntaxError: "name used prior to global declaration".
    global WORK_ROOT, GV_CUTOFF, GV_QS, GV_QH

    parser = argparse.ArgumentParser(
        description="Ciclo termodinámico de pKa y distribución de Boltzmann "
                    "entre estados de protonación (Born cycle, Ho & Coote 2010).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  # Los tres estados clásicos (pKa1 + pKa2 + distribución):
  python protonation_cycle.py \\
      --mol MGLYVAL --neutral R --anion B --dianion C

  # Solo pKa1:
  python protonation_cycle.py --mol MGLYVAL --neutral R --anion B

  # Cadena completa de 5 estados con catión y zwitterión (tautómero de R):
  python protonation_cycle.py \\
      --mol MGLYVAL \\
      --neutral R --anion B --dianion C --protonated A --zwitterion D \\
      --ph 7 --method direct --export xlsx --plot --verbose

  # Con exportación, gráficos y pH de referencia antártico:
  python protonation_cycle.py \\
      --mol SHINORINE --neutral R --anion B --dianion C \\
      --ph 8.1 --export xlsx --plot --verbose
        """,
    )

    parser.add_argument("--neutral",    metavar="DIR",
                        help="Carpeta en runs/ para el neutro R (q=0), p. ej. 'R'")
    parser.add_argument("--anion",      metavar="DIR",
                        help="Carpeta en runs/ para el anión B (q=−1), p. ej. 'B'")
    parser.add_argument("--dianion",    metavar="DIR",
                        help="Carpeta en runs/ para el dianión C (q=−2), p. ej. 'C'")
    parser.add_argument("--protonated", metavar="DIR",
                        help="Carpeta en runs/ para el catión protonado A (q=+1), "
                             "p. ej. 'A'. OPCIONAL: solo tiene efecto junto con "
                             "--zwitterion. Si se omite, no se calcula pKa₀ y la "
                             "cadena arranca en el macroestado neutro/zwitterión.")
    parser.add_argument("--zwitterion", metavar="DIR",
                        help="Carpeta en runs/ para el zwitterión D — TAUTÓMERO de "
                             "--neutral (mismo q=0, mismos átomos), p. ej. 'D'. "
                             "Activa la cadena extendida; requiere también "
                             "--neutral. --protonated es opcional dentro de esta "
                             "cadena.")
    parser.add_argument("--mol", metavar="MOL", default=None,
                        help="Nombre de la carpeta de la molécula dentro del "
                             "directorio actual, p. ej. 'MGLYVAL' o 'SHINORINE'. "
                             "Se usa para construir la ruta por default: "
                             "<cwd>/<MOL>/runs. Ignorado si se pasa --work-root.")
    parser.add_argument("--work-root", metavar="DIR", default=None,
                        help="Ruta explícita a la carpeta runs/ que contiene los "
                             "estados (default: <cwd>/<MOL>/runs usando --mol, "
                             "p. ej. ~/PAPER1/MGLYVAL/runs). Si se pasa, tiene "
                             "prioridad sobre --mol.")
    parser.add_argument("--ph",       type=float, default=PH_ANTARCTIC,
                        help=f"pH de referencia para la distribución "
                             f"(default: {PH_ANTARCTIC}, agua de mar antártica)")
    parser.add_argument("--ph-range", nargs=2, type=float, default=[0.0, 14.0],
                        metavar=("pH_MIN", "pH_MAX"),
                        help="Rango de pH para el diagrama (default: 0 14)")
    parser.add_argument("--ph-step",  type=float, default=0.1,
                        help="Paso del eje de pH (default: 0.1)")
    parser.add_argument("--method",   choices=["auto", "born", "direct"], default="auto",
                        help="Método para ΔG_aq: 'born' (ciclo completo, requiere "
                             "SP fase gas), 'direct' (solo CPCM, sin SP fase gas), "
                             "'auto' usa born si hay SP gas y si no cae a direct "
                             "(default: auto). Los pasos que involucran el "
                             "macroestado neutro/zwitterión (con --protonated/"
                             "--zwitterion) siempre usan 'direct'.")
    parser.add_argument("--temp",     type=float, default=T_STD, choices=TEMPS_K,
                        help=f"Temperatura T (K) para el pKa reportado (default: "
                             f"{T_STD} K, estándar). Restringido a {TEMPS_K} porque "
                             f"G_boltz (usado en los pasos con macroestado "
                             f"neutro/zwitterión, --protonated/--zwitterion) solo "
                             f"está tabulado a esas 4 T; pedir otra T degrada "
                             f"silenciosamente esos pasos a un solo confórmero en "
                             f"vez del ensemble Boltzmann. Útil para reportar el "
                             f"pKa a T de agua de mar antártica (274.15 K) en vez "
                             f"del estándar de 298.15 K. Ver también el escaneo de "
                             f"las 4 T que se imprime siempre como diagnóstico.")
    parser.add_argument("--export",   choices=["csv", "xlsx", "both"], default=None,
                        help="Exportar resultados")
    parser.add_argument("--plot",     action="store_true",
                        help="Generar gráficos de distribución y Boltzmann conformacional")
    parser.add_argument("--verbose",  action="store_true",
                        help="Mostrar contribuciones detalladas del ciclo de Born")
    parser.add_argument("--goodvibes", action="store_true",
                        help="Corre el ciclo COMPLETO dos veces: una con G_corr "
                             "armónico estándar de NumFreq (RRHO, tal cual lo "
                             "reporta ORCA) y otra con la corrección cuasi-armónica "
                             "de GoodVibes (mRRHO de Grimme para los modos de baja "
                             "frecuencia), y al final compara cuánto cambian los "
                             "pKa's y las fracciones de especie. Requiere el "
                             "binario 'goodvibes' instalado "
                             "(pip install goodvibes --break-system-packages).")
    parser.add_argument("--gv-cutoff", type=float, default=GV_CUTOFF, metavar="CM-1",
                        help=f"Frequency cutoff de GoodVibes en cm⁻¹ "
                             f"(default: {GV_CUTOFF:.0f})")
    parser.add_argument("--gv-qs", choices=["grimme", "truhlar"], default=GV_QS,
                        help=f"Método de interpolación de entropía cuasi-armónica "
                             f"de GoodVibes (default: {GV_QS})")
    parser.add_argument("--gv-qh", action="store_true",
                        help="Además de la entropía, aplica también la corrección "
                             "cuasi-armónica de Head-Gordon a la entalpía "
                             "(default: solo la entropía, la más relevante)")
    parser.add_argument("--montecarlo", action="store_true",
                        help="Propaga incertidumbre por Monte Carlo a través de "
                             "toda la cadena de pKa's (ver --mc-*), preservando "
                             "la correlación física de G(H⁺,aq) entre pKa0/pKa1/"
                             "pKa2. Imprime media±std/CI95 junto al valor "
                             "determinístico y, si además se pasa --export, "
                             "agrega las muestras crudas y las bandas de "
                             "confianza de la distribución de especies.")
    parser.add_argument("--mc-samples", type=int, default=10_000, metavar="N",
                        help="Número de muestras Monte Carlo (default: 10000)")
    parser.add_argument("--mc-sigma-elec", type=float, default=0.5, metavar="KCAL",
                        help="σ asumida para la energía electrónica SP por "
                             "estado, en kcal/mol (default: 0.5, típico para "
                             "ωB97X-D3/def2-TZVPD)")
    parser.add_argument("--mc-sigma-gcorr", type=float, default=0.4, metavar="KCAL",
                        help="σ asumida para la corrección térmica G_corr por "
                             "estado, en kcal/mol (default: 0.4; mayor si hay "
                             "n_imag>0 o gcorr_source='harmonic' en vez de mRRHO)")
    parser.add_argument("--mc-sigma-hplus", type=float, default=2.5, metavar="KCAL",
                        help="σ asumida para G(H⁺,aq) — dispersión típica entre "
                             "valores de literatura (Pliego-Riveros −264.61, "
                             "Kelly-Cramer-Truhlar −265.9, Zhan-Dixon ≈−262.4; "
                             "default: 2.5). Esta es sistemática y se comparte "
                             "entre los tres pKa's de la cadena, no se muestrea "
                             "independiente por tramo.")
    parser.add_argument("--mc-seed", type=int, default=None, metavar="N",
                        help="Semilla para reproducibilidad del Monte Carlo "
                             "(default: sin fijar)")
    args = parser.parse_args()

    GV_CUTOFF = args.gv_cutoff
    GV_QS     = args.gv_qs
    GV_QH     = args.gv_qh
    if args.work_root:
        WORK_ROOT = Path(args.work_root).expanduser().resolve()
    elif args.mol:
        WORK_ROOT = (Path.cwd() / args.mol / "runs").resolve()
    else:
        parser.error(
            "Debes indicar --mol MOL (p. ej. --mol MGLYVAL) o, alternativamente, "
            "--work-root /ruta/a/runs. Ya no se asume ./runs en el directorio actual, "
            "porque ahora conviven varias moléculas (MGLYVAL, SHINORINE, etc.) bajo "
            "el mismo PAPER1/."
        )
    if not WORK_ROOT.exists():
        parser.error(
            f"No existe la carpeta runs/: {WORK_ROOT}\n"
            f"Verifica que la molécula '{args.mol}' exista bajo {Path.cwd()} "
            f"(se esperaba {Path.cwd() / (args.mol or '<MOL>') / 'runs'}), "
            f"o usa --work-root /ruta/a/runs para indicar la correcta."
        )

    # Verificar que se pasó al menos un par de estados
    if args.neutral is None and args.anion is None:
        parser.error("Debes proporcionar al menos --neutral y --anion.")

    # El catión protonado A⁺ (--protonated) es OPCIONAL dentro de la cadena
    # extendida: si no se pasa, simplemente se omite pKa₀ y la distribución
    # arranca en el macroestado neutro/zwitterión. Lo único que activa el
    # modo extendido (con tautómero R⇌D) es --zwitterion.
    use_extended = args.zwitterion is not None
    if use_extended:
        missing = [name for name, val in
                   [("--zwitterion", args.zwitterion), ("--neutral", args.neutral)]
                   if val is None]
        if missing:
            parser.error(
                "La cadena extendida (zwitterión) requiere --zwitterion y "
                "--neutral simultáneamente. Falta(n): " + ", ".join(missing) + ". "
                "--protonated (catión A⁺, q=+1) es opcional: si se omite, "
                "simplemente no se calcula pKa₀ y la cadena arranca en el "
                "macroestado neutro/zwitterión."
            )
    elif args.protonated is not None:
        parser.error(
            "--protonated (catión A⁺) solo tiene sentido junto con "
            "--zwitterion (cadena extendida); para la cadena clásica de 3 "
            "estados no se usa. Agrega --zwitterion o quita --protonated."
        )

    print("\n" + "═" * 74)
    print("  PROTONATION CYCLE  —  FAE-MAA / kmilozeta")
    print("═" * 74)
    if args.mol:
        print(f"  molécula : {args.mol}")
    print(f"  runs/ : {WORK_ROOT}")

    if args.goodvibes:
        _run_comparison(args, use_extended)
    elif use_extended:
        _run_extended(args)
    else:
        _run_classic(args)


def _run_classic(args, gcorr_source: str = "harmonic"):
    """
    Cadena original de 3 estados: neutro → anión → dianión.

    gcorr_source controla de dónde sale G_corr para TODOS los estados:
      "harmonic"  → RRHO puro de NumFreq (default, comportamiento original)
      "goodvibes" → mRRHO cuasi-armónico de GoodVibes (ver --goodvibes)

    Retorna un dict con los resultados clave (pKa1, pKa2, df_dist, ...) para
    que _run_comparison() pueda diferenciar entre las dos pasadas.
    """
    data_neutral = data_anion = data_dianion = None

    if args.neutral:
        print(f"\n▸ Cargando estado neutro: {args.neutral}")
        data_neutral = collect_mol_data(args.neutral, verbose=args.verbose,
                                        gcorr_source=gcorr_source)
        print_mol_summary(data_neutral, verbose=args.verbose)

    if args.anion:
        print(f"\n▸ Cargando anión: {args.anion}")
        data_anion = collect_mol_data(args.anion, verbose=args.verbose,
                                      gcorr_source=gcorr_source)
        print_mol_summary(data_anion, verbose=args.verbose)

    if args.dianion:
        print(f"\n▸ Cargando dianión: {args.dianion}")
        data_dianion = collect_mol_data(args.dianion, verbose=args.verbose,
                                        gcorr_source=gcorr_source)
        print_mol_summary(data_dianion, verbose=args.verbose)

    pka1_result = pka2_result = None

    if data_neutral and data_anion:
        print(f"\n▸ Calculando pKa₁  ({args.neutral} → {args.anion} + H⁺)  @ T={args.temp:.2f} K")
        pka1_result = compute_pka(data_neutral, data_anion,
                                  T=args.temp, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if data_anion and data_dianion:
        print(f"\n▸ Calculando pKa₂  ({args.anion} → {args.dianion} + H⁺)  @ T={args.temp:.2f} K")
        pka2_result = compute_pka(data_anion, data_dianion,
                                  T=args.temp, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if pka1_result is None and pka2_result is None:
        print("\n  ⚠  No se pudo calcular ningún pKa. "
              "Revisa los directorios y los SP fase gas.")
        sys.exit(1)

    temp_scan = _pka_temperature_scan(
        [("pKa1", data_neutral, data_anion), ("pKa2", data_anion, data_dianion)],
        gcorr_source=gcorr_source, method=args.method,
    )
    _print_pka_temperature_scan(temp_scan, reported_T=args.temp)

    print_pka_summary([(pka1_result, "pKa₁"), (pka2_result, "pKa₂")])

    pKa1 = pka1_result["pKa"] if (pka1_result and pka1_result["ok"]) else None
    pKa2 = pka2_result["pKa"] if (pka2_result and pka2_result["ok"]) else None

    if pKa1 is not None:
        pH_range = np.arange(args.ph_range[0], args.ph_range[1] + args.ph_step / 2,
                             args.ph_step)
        pH_highlights = [args.ph, PH_STANDARD]
        if pKa1: pH_highlights.append(round(pKa1, 1))
        if pKa2: pH_highlights.append(round(pKa2, 1))
        pH_all = sorted(set(list(pH_range) + pH_highlights))
        pH_all = [p for p in pH_all if args.ph_range[0] <= p <= args.ph_range[1]]

        df_dist = boltzmann_protonation_states(pKa1, pKa2, pH_all)

        pH_display = np.arange(args.ph_range[0], args.ph_range[1] + 0.5, 0.5)
        pH_display_list = sorted(set(
            [round(p, 1) for p in pH_display] + pH_highlights
        ))
        pH_display_list = [p for p in pH_display_list
                           if args.ph_range[0] <= p <= args.ph_range[1]]
        df_display = boltzmann_protonation_states(pKa1, pKa2, pH_display_list)
        print_distribution_table(df_display, pKa1, pKa2, ph_highlight=args.ph)

        print(f"\n  ── Fracciones a pH de referencia ──")
        for pH_ref, label in [
            (args.ph,     f"pH {args.ph} (agua de mar antártica)"),
            (PH_STANDARD, f"pH {PH_STANDARD} (estándar)"),
        ]:
            row = boltzmann_protonation_states(pKa1, pKa2, [pH_ref])
            if not row.empty:
                fn = row["f_neutral"].iloc[0] * 100 if "f_neutral" in row else None
                fa = row["f_anion"].iloc[0]   * 100 if "f_anion"   in row else None
                fd = row["f_dianion"].iloc[0] * 100 if "f_dianion" in row else None
                parts = []
                if fn is not None: parts.append(f"f(AH₂)={fn:.1f}%")
                if fa is not None: parts.append(f"f(AH⁻)={fa:.1f}%")
                if fd is not None: parts.append(f"f(A²⁻)={fd:.1f}%")
                print(f"    {label}: {', '.join(parts)}")

        # ── Monte Carlo (cadena de 3 estados: no hay catión A ni macroestado,
        #    se reusa run_chain_montecarlo() con data_A=None, f_R=1/f_D=0) ──
        chain_mc = None
        if args.montecarlo:
            print(f"\n▸ Propagando incertidumbre Monte Carlo "
                  f"(N={args.mc_samples:,} muestras)...")
            chain_mc = run_chain_montecarlo(
                None, data_neutral, data_anion, data_dianion,
                f_R=1.0, f_D=0.0, pH_values=pH_all,
                T=args.temp, method=args.method, gcorr_source=gcorr_source,
                n_samples=args.mc_samples,
                sigma_elec_kcal=args.mc_sigma_elec,
                sigma_gcorr_kcal=args.mc_sigma_gcorr,
                sigma_hplus_kcal=args.mc_sigma_hplus,
                seed=args.mc_seed,
            )
            print_mc_summary(chain_mc, [
                ("pKa1", pka1_result, "pKa₁"),
                ("pKa2", pka2_result, "pKa₂"),
            ])
    else:
        df_dist = pd.DataFrame()
        chain_mc = None

    data_list = [d for d in [data_neutral, data_anion, data_dianion] if d is not None]

    if args.export and not df_dist.empty:
        export_results([(pka1_result, "pKa1"), (pka2_result, "pKa2")],
                       df_dist, data_list, args.export, chain_mc=chain_mc)

    if args.plot and pKa1 is not None:
        out_dir = WORK_ROOT
        plot_distribution(df_dist, pKa1, pKa2,
                          output_path=out_dir / "species_distribution.png",
                          ph_highlight=args.ph, chain_mc=chain_mc)
        plot_conformer_boltzmann(data_list,
                                 output_path=out_dir / "conformer_boltzmann.png")
        if args.montecarlo:
            plot_pka_histograms(chain_mc, [
                ("pKa1", pka1_result, "pKa₁"),
                ("pKa2", pka2_result, "pKa₂"),
            ], output_path=out_dir / "pka_mc_histograms.png")

    print("\nTerminado.\n")

    return {
        "gcorr_source": gcorr_source,
        "pKa1": pKa1, "pKa2": pKa2,
        "df_dist": df_dist,
        "data_list": data_list,
        "chain_mc": chain_mc,
    }


def _run_extended(args, gcorr_source: str = "harmonic"):
    """
    Cadena de 5 estados:
      A⁺ (catión) ⇌ [R (neutro) ⇌ D (zwitterión)] ⇌ B⁻ (anión) ⇌ C²⁻ (dianión)

    gcorr_source: ver docstring de _run_classic(). Retorna dict con los
    resultados clave (pKa0/1/2, f_R, f_D, df_dist, ...) para _run_comparison().
    """
    data_A = None
    if args.protonated:
        print(f"\n▸ Cargando catión protonado: {args.protonated}")
        data_A = collect_mol_data(args.protonated, verbose=args.verbose,
                                  gcorr_source=gcorr_source)
        print_mol_summary(data_A, verbose=args.verbose)
    else:
        print("\n▸ Catión protonado A⁺ no proporcionado (--protonated omitido): "
              "se omite pKa₀, la cadena arranca en el macroestado neutro/zwitterión.")

    print(f"\n▸ Cargando estado neutro: {args.neutral}")
    data_R = collect_mol_data(args.neutral, verbose=args.verbose,
                              gcorr_source=gcorr_source)
    print_mol_summary(data_R, verbose=args.verbose)

    print(f"\n▸ Cargando zwitterión (tautómero de {args.neutral}): {args.zwitterion}")
    data_D = collect_mol_data(args.zwitterion, verbose=args.verbose,
                              gcorr_source=gcorr_source)
    print_mol_summary(data_D, verbose=args.verbose)

    data_B = data_C = None
    if args.anion:
        print(f"\n▸ Cargando anión: {args.anion}")
        data_B = collect_mol_data(args.anion, verbose=args.verbose,
                                  gcorr_source=gcorr_source)
        print_mol_summary(data_B, verbose=args.verbose)

    if args.dianion:
        print(f"\n▸ Cargando dianión: {args.dianion}")
        data_C = collect_mol_data(args.dianion, verbose=args.verbose,
                                  gcorr_source=gcorr_source)
        print_mol_summary(data_C, verbose=args.verbose)

    # ── Equilibrio de tautómeros R ⇌ D (sin H⁺) ────────────────────────────
    print(f"\n▸ Calculando equilibrio de tautómeros  ({args.neutral} ⇌ {args.zwitterion})  @ T={args.temp:.2f} K")
    taut_result = compute_tautomer_equilibrium(data_R, data_D, T=args.temp, verbose=args.verbose,
                                               gcorr_source=gcorr_source)

    if not taut_result["ok"]:
        # No se pudo cerrar el equilibrio R⇌D (típicamente porque D no tiene
        # confórmero representativo). Esto NO debe abortar todo el script:
        # se cae a usar R solo como estado "neutro" (f_R=100%, f_D=0%) y se
        # sigue con la cadena de pKa para los estados que sí tienen datos.
        print(f"\n  ⚠  No se pudo calcular el equilibrio de tautómeros "
              f"({args.neutral}⇌{args.zwitterion}). Se continúa usando solo "
              f"'{args.neutral}' como estado neutro "
              f"(f_{args.neutral}=100%, f_{args.zwitterion}=0%).")
        for w in taut_result["warnings"]:
            print(f"    ⚠  {w}")

        gR_fallback = _state_G_total_direct(data_R, T=args.temp, gcorr_source=gcorr_source)
        if not gR_fallback.get("ok"):
            print(f"\n  ⚠  Tampoco se pudo obtener G_total para "
                  f"'{args.neutral}'. Abortando.")
            sys.exit(1)

        taut_result = {
            "mol_R": args.neutral, "mol_D": args.zwitterion, "T": args.temp,
            "ok": True,
            "warnings": taut_result["warnings"] + [
                f"Macroestado neutro construido solo con '{args.neutral}' "
                f"(sin datos válidos para '{args.zwitterion}')."
            ],
            "conf_R": gR_fallback["conf"], "conf_D": None,
            "G_R_Eh": gR_fallback["G_total_Eh"], "G_D_Eh": None,
            "dG_taut_kcal": None, "K_taut": None,
            "f_R": 1.0, "f_D": 0.0,
            "G_macro_Eh": gR_fallback["G_total_Eh"],
            "theory_R": gR_fallback["theory_sp"], "theory_D": None,
        }
        print(f"    → G_total({args.neutral}, confórmero {gR_fallback['conf']}) "
              f"= {gR_fallback['G_total_Eh']:.8f} Eh  usado como G_macro")
    elif not args.verbose:
        _print_tautomer_cycle(taut_result)

    if taut_result["conf_D"] is not None:
        conf_label = (f"{args.neutral}:{taut_result['conf_R']} + "
                      f"{args.zwitterion}:{taut_result['conf_D']} (Boltzmann)")
        mol_label  = f"{args.neutral}/{args.zwitterion}"
    else:
        conf_label = f"{args.neutral}:{taut_result['conf_R']} (sin {args.zwitterion})"
        mol_label  = args.neutral

    data_neutral_macro = {
        "mol_name":          mol_label,
        "is_macrostate":     True,
        "G_total_direct_Eh": taut_result["G_macro_Eh"],
        "conf_label":        conf_label,
        "theory_label":      taut_result["theory_R"] or taut_result["theory_D"],
        # R y D son tautómeros (mismo q, mismo N° átomos, ya validado en
        # compute_tautomer_equilibrium), así que cualquiera de los dos sirve
        # como representante del macroestado para la validación cruzada de
        # compute_pka() frente al catión (A) y al anión (B).
        "charge":            data_R.get("charge"),
        "n_atoms":           data_R.get("n_atoms"),
    }

    # ── Cadena de pKa: A→[R/D], [R/D]→B, B→C ───────────────────────────────
    pka0_result = pka1_result = pka2_result = None

    if data_A:
        print(f"\n▸ Calculando pKa₀  ({args.protonated} → neutro/zwitterión + H⁺)  @ T={args.temp:.2f} K")
        pka0_result = compute_pka(data_A, data_neutral_macro,
                                  T=args.temp, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if data_B:
        print(f"\n▸ Calculando pKa₁  (neutro/zwitterión → {args.anion} + H⁺)  @ T={args.temp:.2f} K")
        pka1_result = compute_pka(data_neutral_macro, data_B,
                                  T=args.temp, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if data_B and data_C:
        print(f"\n▸ Calculando pKa₂  ({args.anion} → {args.dianion} + H⁺)  @ T={args.temp:.2f} K")
        pka2_result = compute_pka(data_B, data_C,
                                  T=args.temp, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

        # Nota #2 del review: pKa1 siempre usa 'direct' (macroestado), pero
        # pKa2 puede resolver a 'born' con method='auto'. Ambos métodos usan
        # constantes G(H+,aq) distintas y no intercambiables — este chequeo
        # cuantifica la discrepancia cuando aplica, sin alterar pka2_result.
        ph_scale_diag = _cross_validate_ph_scale(
            data_B, data_C, pka2_result, T=args.temp, gcorr_source=gcorr_source,
        )
        if ph_scale_diag.get("applicable"):
            _print_ph_scale_diagnostic(ph_scale_diag, label="pKa₂ (B→C)")

        # Con method='auto', pKa1 (macroestado) y pKa2 (simple) pueden ahora
        # resolver a métodos distintos ('born' vs 'direct') dependiendo de si
        # existen SP de fase gas para B y C. Antes esto quedaba escondido
        # porque pKa1 SIEMPRE caía a 'direct' (bug ya corregido arriba); con
        # el fix, mezclar 'born' (G_HPLUS_AQ_BORN) y 'direct'
        # (G_HPLUS_AQ_DIRECT) en la misma cadena de pKa's suma dos escalas de
        # G(H+,aq) no intercambiables sin que nada lo señale. Se aborta en
        # vez de solo advertir.
        if (args.method == "auto" and pka1_result is not None
                and pka1_result.get("ok") and pka2_result.get("ok")
                and pka1_result.get("method") != pka2_result.get("method")):
            raise SystemExit(
                f"\n❌ INCONSISTENT G(H+,aq) REFERENCE within the same pKa "
                f"chain under method='auto':\n"
                f"   pKa₁ ([R/D]→{args.anion}) resolved to "
                f"'{pka1_result.get('method')}'\n"
                f"   pKa₂ ({args.anion}→{args.dianion}) resolved to "
                f"'{pka2_result.get('method')}'\n"
                f"   These use different, non-interchangeable G(H+,aq) "
                f"constants (Pliego-Riveros {G_HPLUS_AQ_DIRECT:.2f} kcal/mol "
                f"vs. Tissandier/Kelly {G_HPLUS_AQ_BORN:.2f} kcal/mol).\n"
                f"   Fix: either provide gas-phase SPs for all species in the "
                f"chain (so 'auto' resolves to 'born' everywhere), or force "
                f"a single method explicitly with --method direct or "
                f"--method born."
            )

    # Nota #8 del review: reconstruye tautómero+macroestado y recalcula
    # pKa0/1/2 en las 4 T de TEMPS_K (no solo re-evalúa RT·ln10 a T fija).
    ext_temp_scan = _extended_temperature_scan(
        data_A, data_R, data_D, data_B, data_C,
        gcorr_source=gcorr_source, method=args.method,
    )
    _print_extended_temperature_scan(ext_temp_scan, reported_T=args.temp)

    print_pka_summary([
        (pka0_result, "pKa₀ (A→R/D)"),
        (pka1_result, "pKa₁ (R/D→B)"),
        (pka2_result, "pKa₂ (B→C)"),
    ])

    pKa0 = pka0_result["pKa"] if (pka0_result and pka0_result["ok"]) else None
    pKa1 = pka1_result["pKa"] if (pka1_result and pka1_result["ok"]) else None
    pKa2 = pka2_result["pKa"] if (pka2_result and pka2_result["ok"]) else None
    f_R, f_D = taut_result["f_R"], taut_result["f_D"]

    if pKa0 is None and pKa1 is None:
        print("\n  ⚠  No se pudo calcular ningún pKa de la cadena principal. "
              "Revisa los directorios.")
        sys.exit(1)

    # ── Distribución de especies ──────────────────────────────────────────
    pH_range = np.arange(args.ph_range[0], args.ph_range[1] + args.ph_step / 2, args.ph_step)
    pH_highlights = [args.ph, PH_STANDARD]
    for pk in (pKa0, pKa1, pKa2):
        if pk is not None:
            pH_highlights.append(round(pk, 1))
    pH_all = sorted(set(list(pH_range) + pH_highlights))
    pH_all = [p for p in pH_all if args.ph_range[0] <= p <= args.ph_range[1]]

    df_dist = boltzmann_5state(pKa0, pKa1, pKa2, f_R, f_D, pH_all)

    pH_display = np.arange(args.ph_range[0], args.ph_range[1] + 0.5, 0.5)
    pH_display_list = sorted(set([round(p, 1) for p in pH_display] + pH_highlights))
    pH_display_list = [p for p in pH_display_list if args.ph_range[0] <= p <= args.ph_range[1]]
    df_display = boltzmann_5state(pKa0, pKa1, pKa2, f_R, f_D, pH_display_list)
    print_distribution_table_5(df_display, pKa0, pKa1, pKa2, ph_highlight=args.ph)

    print(f"\n  ── Fracciones a pH de referencia ──")
    for pH_ref, label in [
        (args.ph,     f"pH {args.ph} (agua de mar antártica)"),
        (PH_STANDARD, f"pH {PH_STANDARD} (estándar)"),
    ]:
        row = boltzmann_5state(pKa0, pKa1, pKa2, f_R, f_D, [pH_ref])
        if not row.empty:
            r0 = row.iloc[0]
            parts = []
            if "f_cation" in row: parts.append(f"f(A⁺)={r0['f_cation']*100:.1f}%")
            parts.append(f"f(R)={r0['f_R']*100:.1f}%")
            parts.append(f"f(D)={r0['f_D']*100:.1f}%")
            if "f_anion" in row: parts.append(f"f(B⁻)={r0['f_anion']*100:.1f}%")
            if "f_dianion" in row: parts.append(f"f(C²⁻)={r0['f_dianion']*100:.1f}%")
            print(f"    {label}: {', '.join(parts)}")

    # ── Monte Carlo (cadena completa de 5 estados, G(H⁺,aq) correlacionado
    #    entre pKa0/pKa1/pKa2 — ver run_chain_montecarlo()) ─────────────────
    chain_mc = None
    if args.montecarlo:
        print(f"\n▸ Propagando incertidumbre Monte Carlo por toda la cadena "
              f"(N={args.mc_samples:,} muestras)...")
        chain_mc = run_chain_montecarlo(
            data_A, data_neutral_macro, data_B, data_C,
            f_R=f_R, f_D=f_D, pH_values=pH_all,
            T=args.temp, method=args.method, gcorr_source=gcorr_source,
            n_samples=args.mc_samples,
            sigma_elec_kcal=args.mc_sigma_elec,
            sigma_gcorr_kcal=args.mc_sigma_gcorr,
            sigma_hplus_kcal=args.mc_sigma_hplus,
            seed=args.mc_seed,
        )
        print_mc_summary(chain_mc, [
            ("pKa0", pka0_result, "pKa₀ (A→R/D)"),
            ("pKa1", pka1_result, "pKa₁ (R/D→B)"),
            ("pKa2", pka2_result, "pKa₂ (B→C)"),
        ])

    # ── Exportar ─────────────────────────────────────────────────────────
    data_list = [d for d in [data_A, data_R, data_D, data_B, data_C] if d is not None]

    if args.export and not df_dist.empty:
        export_results(
            [(pka0_result, "pKa0"), (pka1_result, "pKa1"), (pka2_result, "pKa2")],
            df_dist, data_list, args.export, taut_result=taut_result,
            chain_mc=chain_mc,
        )

    # ── Gráficos ─────────────────────────────────────────────────────────
    if args.plot:
        out_dir = WORK_ROOT
        plot_distribution_5(
            df_dist, pKa0, pKa1, pKa2,
            output_path=out_dir / "species_distribution.png",
            ph_highlight=args.ph,
            chain_mc=chain_mc,
        )
        plot_conformer_boltzmann(
            data_list, output_path=out_dir / "conformer_boltzmann.png",
        )
        if args.montecarlo:
            plot_pka_histograms(chain_mc, [
                ("pKa0", pka0_result, "pKa₀ (A→R/D)"),
                ("pKa1", pka1_result, "pKa₁ (R/D→B)"),
                ("pKa2", pka2_result, "pKa₂ (B→C)"),
            ], output_path=out_dir / "pka_mc_histograms.png")

    print("\nTerminado.\n")

    return {
        "gcorr_source": gcorr_source,
        "pKa0": pKa0, "pKa1": pKa1, "pKa2": pKa2,
        "f_R": f_R, "f_D": f_D,
        "df_dist": df_dist,
        "data_list": data_list,
        "chain_mc": chain_mc,
    }


def _run_comparison(args, use_extended: bool):
    """
    Corre el ciclo COMPLETO dos veces con distinta fuente de G_corr —
    armónico estándar (NumFreq/RRHO) vs. cuasi-armónico (GoodVibes/mRRHO) —
    y compara al final cuánto cambian los pKa's y la distribución de especies.

    Activado con --goodvibes. Usa _run_classic() o _run_extended() según
    corresponda, sin duplicar lógica: solo cambia gcorr_source entre pasadas.
    """
    gv_python = _resolve_gv_python()
    try:
        check = subprocess.run(
            [gv_python, "-c", "import goodvibes"],
            capture_output=True, timeout=10,
        )
    except (subprocess.TimeoutExpired, OSError) as e:
        print(f"\n  ⚠  No se pudo verificar goodvibes con {gv_python}: {e}\n")
        sys.exit(1)
    if check.returncode != 0:
        print(f"\n  ⚠  No se puede importar goodvibes con {gv_python}.")
        print(f"     Instalar con: conda activate {GV_CONDA_ENV or '<env>'} && "
              f"pip install goodvibes --break-system-packages")
        print("     (o, en un entorno conda: conda install -c conda-forge goodvibes)\n")
        sys.exit(1)

    runner = _run_extended if use_extended else _run_classic

    print("\n" + "═" * 74)
    print("  PASADA 1/2  —  G_corr ARMÓNICO ESTÁNDAR (NumFreq, RRHO)")
    print("═" * 74)
    res_harm = runner(args, gcorr_source="harmonic")

    qh_tag = ", +qh (Head-Gordon)" if args.gv_qh else ""
    print("\n" + "═" * 74)
    print(f"  PASADA 2/2  —  G_corr CUASI-ARMÓNICO (GoodVibes, "
          f"qs={args.gv_qs}, cutoff={args.gv_cutoff:.0f} cm⁻¹{qh_tag})")
    print("═" * 74)
    res_gv = runner(args, gcorr_source="goodvibes")

    _print_goodvibes_comparison(res_harm, res_gv, extended=use_extended, args=args)


def _print_goodvibes_comparison(res_harm: dict, res_gv: dict,
                                 extended: bool, args) -> None:
    """
    Imprime una tabla comparativa entre la pasada armónica (NumFreq) y la
    cuasi-armónica (GoodVibes): Δ en cada pKa y Δ en las fracciones de
    especie a los pH de referencia (--ph y pH 7.0 estándar).
    """
    sep = "═" * 74
    print(f"\n{sep}")
    print("  COMPARACIÓN  —  G_corr ARMÓNICO (NumFreq) vs. CUASI-ARMÓNICO (GoodVibes)")
    print(sep)

    if res_harm is None or res_gv is None:
        print("\n  ⚠  Una de las dos pasadas no produjo resultados; no se puede comparar.\n")
        return

    pka_keys = ([("pKa0", "pKa₀ (A→R/D)")] if extended else []) + [
        ("pKa1", "pKa₁"), ("pKa2", "pKa₂"),
    ]

    print(f"\n  {'':<20} {'armónico':>12} {'GoodVibes':>12} {'Δ (GV−arm)':>12}")
    print(f"  {'-'*20} {'-'*12} {'-'*12} {'-'*12}")
    any_pka = False
    for key, label in pka_keys:
        v1, v2 = res_harm.get(key), res_gv.get(key)
        if v1 is None and v2 is None:
            continue
        any_pka = True
        v1s = f"{v1:.2f}" if v1 is not None else "—"
        v2s = f"{v2:.2f}" if v2 is not None else "—"
        ds  = f"{(v2 - v1):+.2f}" if (v1 is not None and v2 is not None) else "—"
        print(f"  {label:<20} {v1s:>12} {v2s:>12} {ds:>12}")
    if not any_pka:
        print("  (sin pKa's calculados en ninguna de las dos pasadas)")

    # ── Fracciones de especie a los pH de referencia ────────────────────
    def _row_at_ph(df: Optional[pd.DataFrame], ph: float):
        if df is None or df.empty or "pH" not in df.columns:
            return None
        sub = df[np.isclose(df["pH"].astype(float), ph, atol=1e-6)]
        return sub.iloc[0] if not sub.empty else None

    for ph_ref, ph_label in [
        (args.ph,     f"pH {args.ph} (referencia)"),
        (PH_STANDARD, f"pH {PH_STANDARD} (estándar)"),
    ]:
        row_h = _row_at_ph(res_harm.get("df_dist"), ph_ref)
        row_g = _row_at_ph(res_gv.get("df_dist"), ph_ref)
        if row_h is None or row_g is None:
            continue
        frac_cols = [c for c in row_h.index if c.startswith("f_")]
        if not frac_cols:
            continue
        print(f"\n  Fracciones de especie a {ph_label}:")
        print(f"  {'especie':<14} {'armónico':>10} {'GoodVibes':>10} {'Δ (pp)':>10}")
        for c in frac_cols:
            fh = float(row_h[c]) * 100 if pd.notna(row_h[c]) else None
            fg = float(row_g[c]) * 100 if (c in row_g and pd.notna(row_g[c])) else None
            if fh is None and fg is None:
                continue
            fhs = f"{fh:.2f}%" if fh is not None else "—"
            fgs = f"{fg:.2f}%" if fg is not None else "—"
            dps = f"{(fg - fh):+.2f}" if (fh is not None and fg is not None) else "—"
            print(f"  {c:<14} {fhs:>10} {fgs:>10} {dps:>10}")

    print(f"\n  Nota: la aproximación armónica pura sobreestima la entropía de los")
    print(f"  modos de baja frecuencia (rotores internos, <{args.gv_cutoff:.0f} cm⁻¹) —")
    print(f"  típicos en las cadenas laterales Gly/Val y grupos OH/NH libres de")
    print(f"  M-Gly:Val. Diferencias de pKa >0.3-0.5 unidades, o de fracción de")
    print(f"  especie >5 puntos porcentuales entre ambas columnas, son razón para")
    print(f"  preferir el valor de GoodVibes (mRRHO) en el paper en vez del RRHO puro.")
    print(f"{sep}\n")


if __name__ == "__main__":
    main()
