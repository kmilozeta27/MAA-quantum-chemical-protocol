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
  │       │   └── 06_sp_gas/conf_001/sp_gas.out  ← SP fase gas (sin CPCM) [NUEVO]
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

Para el ciclo de pKa se requieren los SP en fase gas (06_sp_gas).
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
SUBDIR_SP_GAS  = "06_sp_gas"   # SPs sin CPCM para el ciclo de pKa

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
            f"No se pudo verificar la carga entre '{label_a}' y '{label_b}' "
            f"(carga no determinada desde el .out de ORCA en al menos uno de "
            f"los dos estados). Verificar manualmente que la relación de "
            f"carga esperada (Δq={dq_expected:+d}) sea correcta."
        )
    elif (qb - qa) != dq_expected:
        msgs.append(
            f"❌ CARGA INCONSISTENTE: '{label_a}' (q={qa:+d}) → '{label_b}' "
            f"(q={qb:+d}) implica Δq={qb - qa:+d}, pero se esperaba "
            f"Δq={dq_expected:+d}. Revisar si el mapeo de carpetas "
            f"(--protonated/--neutral/--anion/--dianion/--zwitterion) "
            f"corresponde al estado que crees que es."
        )

    if na is None or nb is None:
        msgs.append(
            f"No se pudo verificar N° de átomos entre '{label_a}' y '{label_b}' "
            f"(N° de átomos no determinado desde el .out de ORCA en al menos "
            f"uno de los dos estados). Verificar manualmente que la relación "
            f"esperada (ΔN={dn_expected:+d}) sea correcta."
        )
    elif (nb - na) != dn_expected:
        msgs.append(
            f"❌ N° DE ÁTOMOS INCONSISTENTE: '{label_a}' (N={na}) → "
            f"'{label_b}' (N={nb}) implica ΔN={nb - na:+d}, pero se esperaba "
            f"ΔN={dn_expected:+d}. Revisar si el mapeo de carpetas "
            f"corresponde al estado que crees que es, o si el confórmero "
            f"representativo (best_conf) corresponde a la geometría "
            f"esperada."
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
    Corre GoodVibes (`python -m goodvibes`, con el intérprete resuelto por
    _resolve_gv_python()) sobre un freq.out de ORCA y extrae la corrección de
    Gibbs CUASI-ARMÓNICA (mRRHO de Grimme para la entropía, y opcionalmente
    Head-Gordon para la entalpía si qh=True).

    A diferencia de parse_thermochem() (RRHO puro, como lo reporta ORCA
    directamente), GoodVibes reinterpola los modos de baja frecuencia
    (< cutoff, default 100 cm⁻¹) como rotores libres en vez de osciladores
    armónicos, evitando la sobreestimación de S (y por ende de −TS) que la
    aproximación armónica pura produce ahí.

    NOTA sobre la invocación: se usa `[python_exe, "-m", "goodvibes", ...]`
    en vez de depender de que el binario `goodvibes` esté en el PATH del
    proceso actual (`shutil.which("goodvibes")` se usaba antes, pero eso
    falla silenciosamente en shells no interactivas — SSH sin `conda activate`,
    cron, nohup en otra sesión — donde el PATH no incluye el bin/ del entorno
    sciencebitch aunque GoodVibes esté perfectamente instalado ahí). Misma
    estrategia que ya usa results_collector.py.

    Retorna dict con claves compatibles con el resto del pipeline:
        G_corr, H_corr, S_total (Eh, Eh/K), n_imag, ok, warning
    más los campos crudos de GoodVibes (scf_energy, qh_gibbs_free_energy).
    """
    null = {
        "G_corr": None, "H_corr": None, "S_total": None,
        "scf_energy": None, "qh_gibbs_free_energy": None,
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

    with tempfile.TemporaryDirectory() as tmp:
        json_path = Path(tmp) / "gv_result.json"
        cmd = [
            python_exe, "-m", "goodvibes",
            "--temp", str(T), "-f", str(cutoff), "--qs", qs,
        ]
        if qh:
            cmd.append("--qh")
        cmd += ["--json", str(json_path), str(freq_out)]

        try:
            proc = subprocess.run(
                cmd, cwd=tmp, capture_output=True, text=True, timeout=180,
            )
        except (subprocess.TimeoutExpired, OSError) as e:
            return {**null, "warning": f"Fallo al ejecutar goodvibes ({python_exe} -m goodvibes): {e}"}

        if not json_path.exists():
            stderr_tail = (proc.stderr or "").strip()[-400:]
            return {**null, "warning": (
                f"goodvibes ({python_exe}) no generó JSON (returncode={proc.returncode}). "
                f"{stderr_tail}"
            )}

        try:
            payload = json.loads(json_path.read_text())
        except (OSError, json.JSONDecodeError) as e:
            return {**null, "warning": f"No se pudo leer el JSON de goodvibes: {e}"}

    results = payload.get("results") or []
    if not results:
        return {**null, "warning": "goodvibes no devolvió resultados para este archivo."}

    thermo = results[0].get("thermo") or {}
    scf_energy = thermo.get("scf_energy")
    qh_gibbs   = thermo.get("qh_gibbs_free_energy")
    # Entalpía: si --qh está activo, usa la entalpía cuasi-armónica (Head-Gordon);
    # si no, GoodVibes deja la entalpía en su forma armónica estándar.
    enthalpy   = thermo.get("qh_enthalpy") if qh else thermo.get("enthalpy")
    # Entropía: siempre cuasi-armónica (mRRHO), que es el objeto del análisis.
    entropy    = thermo.get("qh_entropy")
    im_freqs   = thermo.get("im_frequency_wn") or []
    n_imag     = len(im_freqs) if isinstance(im_freqs, list) else 0

    if scf_energy is None or qh_gibbs is None:
        return {**null, "warning": (
            "goodvibes no reportó scf_energy/qh_gibbs_free_energy "
            "(freq.out posiblemente incompleto o mal formado)."
        )}

    G_corr = qh_gibbs - scf_energy
    H_corr = (enthalpy - scf_energy) if enthalpy is not None else None

    return {
        "G_corr":               G_corr,
        "H_corr":               H_corr,
        "S_total":              entropy,      # Eh/K
        "scf_energy":           scf_energy,
        "qh_gibbs_free_energy": qh_gibbs,
        "n_imag":               n_imag,
        "ok":                   True,
        "warning":              None,
    }


def get_thermochem(
    freq_out: Path,
    T: float = T_STD,
    source: str = "harmonic",
    gv_cutoff: float = GV_CUTOFF,
    gv_qs: str = GV_QS,
    gv_qh: bool = GV_QH,
) -> dict:
    """
    Wrapper único para obtener G_corr, con dos fuentes intercambiables:

      source="harmonic"  → parse_thermochem()  (RRHO puro, tal cual ORCA)
      source="goodvibes"  → run_goodvibes()      (mRRHO cuasi-armónico)

    Retorna siempre el shape de parse_thermochem() (E_elec_freq, ZPE,
    H_corr, TS_total, S_total, G_corr, G_total_freq, T_parsed, n_imag, ok),
    más "source" y "warning" para trazabilidad. Esto permite que
    collect_mol_data(), _state_G_total_direct() y compute_pka() usen la
    misma interfaz sin importar la fuente.
    """
    if source not in ("harmonic", "goodvibes"):
        raise ValueError(f"source debe ser 'harmonic' o 'goodvibes', no {source!r}")

    if source == "harmonic":
        tc = parse_thermochem(freq_out, T=T)
        tc["source"] = "harmonic"
        tc["warning"] = None
        return tc

    gv = run_goodvibes(freq_out, T=T, cutoff=gv_cutoff, qs=gv_qs, qh=gv_qh)
    TS_total = (T * gv["S_total"]) if gv.get("S_total") is not None else None
    return {
        "E_elec_freq":  gv.get("scf_energy"),
        "ZPE":          None,   # GoodVibes no expone ZPE por separado en G_corr
        "H_corr":       gv.get("H_corr"),
        "TS_total":     TS_total,
        "S_total":      gv.get("S_total"),
        "G_corr":       gv.get("G_corr"),
        "G_total_freq": gv.get("qh_gibbs_free_energy"),
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
      - E_sp_gas  : energía SP sin CPCM (Hartree, de 06_sp_gas — opcional)
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
                "warnings": [f"Directorio no encontrado: {mol_dir}"],
                "conformers": [], "df": pd.DataFrame()}

    sp_dirs = sorted((mol_dir / SUBDIR_SP).glob("conf_*")) \
              if (mol_dir / SUBDIR_SP).exists() else []

    if not sp_dirs:
        return {"mol_name": mol_name, "ok": False,
                "warnings": [f"Sin confórmeros SP en {mol_dir / SUBDIR_SP}"],
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
        # gcorr_source="harmonic" → RRHO puro de NumFreq (ORCA); "goodvibes"
        # → mRRHO cuasi-armónico (ver get_thermochem()).
        freq_file = mol_dir / SUBDIR_FREQ / conf_id / "freq.out"
        tc = get_thermochem(freq_file, T=T_STD, source=gcorr_source,
                            gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)
        if gcorr_source == "goodvibes" and not tc["ok"] and tc.get("warning"):
            warns.append(f"[{mol_name}/{conf_id}] GoodVibes: {tc['warning']}")

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
                    f"[{mol_name}] Sin correcciones térmicas válidas. "
                    f"Boltzmann conformacional no disponible."
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
            f"[{mol_name}] {n_imag_any} confórmero(s) con frecuencias imaginarias. "
            f"Revisar geometrías antes de usar en el ciclo de pKa."
        )

    # Advertencia si falta SP gas
    n_missing_gas = df["E_sp_gas"].isna().sum()
    if n_missing_gas > 0:
        warns.append(
            f"[{mol_name}] {n_missing_gas}/{len(df)} confórmeros sin SP fase gas "
            f"(esperado en runs/{mol_name}/{SUBDIR_SP_GAS}/conf_NNN/sp_gas.out). "
            f"El ciclo de pKa requiere SP gas del confórmero representativo."
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
            f"[{mol_name}] Carga inconsistente entre confórmeros: "
            f"{charges_seen}. No se puede validar contra estados vecinos; "
            f"revisar runs/{mol_name}/{SUBDIR_SP}/conf_*/sp.out."
        )
    if len(natoms_seen) > 1:
        warns.append(
            f"[{mol_name}] N° de átomos inconsistente entre confórmeros: "
            f"{natoms_seen}. No se puede validar contra estados vecinos; "
            f"revisar runs/{mol_name}/{SUBDIR_SP}/conf_*/sp.out."
        )
    if not charges_seen or not natoms_seen:
        warns.append(
            f"[{mol_name}] No se pudo determinar carga y/o N° de átomos "
            f"desde ningún .out. La validación cruzada de carga/N° átomos "
            f"frente a estados vecinos quedará deshabilitada para este "
            f"estado (solo se emitirá una advertencia informativa, no un "
            f"error)."
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
                           gcorr_source: str = "harmonic") -> dict:
    """
    G_total_CPCM (E_SP_aq + G_corr) del confórmero representativo de un estado,
    para el método directo. Funciona con dos tipos de entrada:

      1. Un dict "normal" salido de collect_mol_data(): se relee freq.out del
         confórmero representativo (best_conf) y se recalcula G_corr a T.

      2. Un dict "macroestado" (data["is_macrostate"] = True), construido por
         compute_tautomer_equilibrium(): ya trae G_total_direct_Eh precalculado
         (combinación de Boltzmann entre neutro y zwitterión) y se usa tal cual,
         sin volver a tocar archivos.

    Retorna dict: G_total_Eh, conf, theory_sp, theory_freq, ok, warn (str|None).
    """
    if data.get("is_macrostate"):
        return {
            "G_total_Eh":  data["G_total_direct_Eh"],
            "conf":        data.get("conf_label", "macroestado"),
            "theory_sp":   data.get("theory_label"),
            "theory_freq": data.get("theory_label"),
            "ok":          True,
            "warn":        None,
        }

    mol  = data["mol_name"]
    conf = data.get("best_conf")
    if conf is None:
        return {"ok": False, "warn": f"No hay confórmero representativo para {mol}."}

    df = data["df"]
    row = df[df["conformer"] == conf].iloc[0]
    E_aq = row["E_sp_aq"]
    if E_aq is None or pd.isna(E_aq):
        return {"ok": False, "warn": f"Sin SP acuoso para {mol}/{conf}."}

    freq_file = WORK_ROOT / mol / SUBDIR_FREQ / conf / "freq.out"
    tc = get_thermochem(freq_file, T=T, source=gcorr_source,
                        gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)
    G_corr = tc["G_corr"] if tc["ok"] else None
    warn = None
    if not tc["ok"]:
        extra = f" ({tc['warning']})" if tc.get("warning") else ""
        warn = f"G_corr ({gcorr_source}) no disponible para {mol}/{conf}. Se usa G≈E_SP (sin ZPE/entropía).{extra}"

    G_total = E_aq + (G_corr if G_corr is not None else 0.0)
    theory_sp   = parse_theory_level(WORK_ROOT / mol / SUBDIR_SP / conf / "sp.out")
    theory_freq = parse_theory_level(freq_file)

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
            f"Abortando equilibrio de tautómeros ({mol_R}⇌{mol_D}) por "
            f"inconsistencia de carga/N° de átomos."
        )
        return result

    gR = _state_G_total_direct(data_R, T=T, gcorr_source=gcorr_source)
    gD = _state_G_total_direct(data_D, T=T, gcorr_source=gcorr_source)

    if not gR.get("ok"):
        _warn(gR.get("warn") or f"No se pudo obtener G_total para {mol_R}.")
        return result
    if not gD.get("ok"):
        _warn(gD.get("warn") or f"No se pudo obtener G_total para {mol_D}.")
        return result
    if gR.get("warn"):
        warns.append(gR["warn"])
    if gD.get("warn"):
        warns.append(gD["warn"])

    if gR["theory_sp"] and gD["theory_sp"] and gR["theory_sp"] != gD["theory_sp"]:
        _warn(
            f"NIVEL DE TEORÍA INCONSISTENTE en el SP acuoso entre "
            f"{mol_R} ('{gR['theory_sp']}') y {mol_D} ('{gD['theory_sp']}'). "
            f"La comparación de tautómeros no es válida si difieren base/funcional/CPCM."
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
    print(f"  EQUILIBRIO DE TAUTÓMEROS (sin H⁺) — {r['mol_R']} ⇌ {r['mol_D']}")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. R (neutro)      : {r['conf_R']}")
    print(f"  Conf. D (zwitterión)  : {r['conf_D']}")
    print(f"  Nivel SP(R)  : {r['theory_R'] or '—'}")
    print(f"  Nivel SP(D)  : {r['theory_D'] or '—'}")
    print()
    print(f"  {'G_total_CPCM(R)  [Eh]':<40} {r['G_R_Eh']:>14.8f}")
    print(f"  {'G_total_CPCM(D)  [Eh]':<40} {r['G_D_Eh']:>14.8f}")
    print(f"  {'─'*40}  {'─'*14}")
    print(f"  {'ΔG_taut = G(D) − G(R)  [kcal/mol]':<40} {r['dG_taut_kcal']:>+14.3f}")
    print(f"  {'K_taut = [D]/[R]':<40} {r['K_taut']:>14.4g}")
    print()
    print(f"  f(R, neutro)      = {r['f_R']*100:6.2f} %")
    print(f"  f(D, zwitterión)  = {r['f_D']*100:6.2f} %")
    dominant = "R (neutro)" if r["f_R"] >= r["f_D"] else "D (zwitterión)"
    print(f"  Tautómero dominante: {dominant}")
    print(f"\n  G_macroestado neutro (Boltzmann R+D) = {r['G_macro_Eh']:.8f} Eh")
    print(f"  (este valor combinado es el que entra al ciclo de pKa como")
    print(f"   estado 'neutro/zwitterión' frente al catión y al anión)")
    if r["warnings"]:
        print(f"\n  Advertencias:")
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
                  fase gas (06_sp_gas) además del SP acuoso (CPCM):
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
            f"Abortando cálculo de pKa ({mol_acid}→{mol_base}) por "
            f"inconsistencia de carga/N° de átomos. Revisar el mapeo de "
            f"estados antes de continuar."
        )
        return result

    # ── Rama simplificada: uno de los dos lados es un macroestado combinado ─
    # (neutro+zwitterión, construido por compute_tautomer_equilibrium). No hay
    # forma físicamente limpia de combinar dos SP de fase gas de especies
    # distintas en un ciclo de Born, así que este caso SIEMPRE usa el método
    # directo, sin importar lo que pida --method.
    if data_acid.get("is_macrostate") or data_base.get("is_macrostate"):
        if method == "born":
            _warn(
                f"method='born' solicitado, pero uno de los estados "
                f"({mol_acid if data_acid.get('is_macrostate') else mol_base}) "
                f"es un macroestado combinado neutro+zwitterión. El ciclo de Born "
                f"completo requiere SP fase gas de una única especie; no está "
                f"definido para una mezcla. Se usa el método directo."
            )

        g_acid = _state_G_total_direct(data_acid, T=T, gcorr_source=gcorr_source)
        g_base = _state_G_total_direct(data_base, T=T, gcorr_source=gcorr_source)

        result["conf_acid"] = g_acid.get("conf")
        result["conf_base"] = g_base.get("conf")
        result["theory_sp_acid"] = g_acid.get("theory_sp")
        result["theory_sp_base"] = g_base.get("theory_sp")

        if not g_acid.get("ok"):
            _warn(g_acid.get("warn") or f"No se pudo obtener G_total para {mol_acid}.")
            return result
        if not g_base.get("ok"):
            _warn(g_base.get("warn") or f"No se pudo obtener G_total para {mol_base}.")
            return result
        if g_acid.get("warn"):
            warns.append(g_acid["warn"])
        if g_base.get("warn"):
            warns.append(g_base["warn"])

        if (g_acid["theory_sp"] and g_base["theory_sp"]
                and g_acid["theory_sp"] != g_base["theory_sp"]):
            _warn(
                f"NIVEL DE TEORÍA INCONSISTENTE entre {mol_acid} "
                f"('{g_acid['theory_sp']}') y {mol_base} ('{g_base['theory_sp']}')."
            )

        G_total_acid = g_acid["G_total_Eh"]
        G_total_base = g_base["G_total_Eh"]
        result["G_total_acid_Eh"] = G_total_acid
        result["G_total_base_Eh"] = G_total_base
        result["G_hplus_aq_kcal"] = G_HPLUS_AQ_DIRECT

        dG_aq_kcal = (G_total_base - G_total_acid) * HARTREE_TO_KCAL + G_HPLUS_AQ_DIRECT
        result["dG_aq_kcal"] = dG_aq_kcal
        result["method"] = "direct"
        result["pKa"] = dG_aq_kcal / (R_KCAL * T * LN10)
        result["ok"] = True

        if verbose:
            _print_direct_cycle(result)

        return result

    # ── Obtener confórmero representativo de cada estado ──────────────────
    conf_acid = data_acid.get("best_conf")
    conf_base = data_base.get("best_conf")
    result["conf_acid"] = conf_acid
    result["conf_base"] = conf_base

    if conf_acid is None:
        _warn(f"No hay confórmero representativo para {mol_acid}.")
        return result
    if conf_base is None:
        _warn(f"No hay confórmero representativo para {mol_base}.")
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
            f"NIVEL DE TEORÍA INCONSISTENTE en el SP acuoso entre "
            f"{mol_acid} ('{sp_acid_lvl}') y {mol_base} ('{sp_base_lvl}'). "
            f"El ciclo de pKa NO es válido si la base/funcional/CPCM difieren "
            f"entre estados — revisar antes de reportar."
        )
    if freq_acid_lvl and freq_base_lvl and freq_acid_lvl != freq_base_lvl:
        _warn(
            f"NIVEL DE TEORÍA INCONSISTENTE en NumFreq entre "
            f"{mol_acid} ('{freq_acid_lvl}') y {mol_base} ('{freq_base_lvl}')."
        )

    # ── Energías SP acuosas (CPCM) ────────────────────────────────────────
    E_aq_acid = row_acid["E_sp_aq"]
    E_aq_base = row_base["E_sp_aq"]
    result["E_sp_aq_acid_Eh"] = E_aq_acid
    result["E_sp_aq_base_Eh"] = E_aq_base

    if E_aq_acid is None or pd.isna(E_aq_acid):
        _warn(f"Sin SP acuoso para {mol_acid}/{conf_acid}.")
        return result
    if E_aq_base is None or pd.isna(E_aq_base):
        _warn(f"Sin SP acuoso para {mol_base}/{conf_base}.")
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
            f"method='born' solicitado pero faltan SP fase gas para "
            f"{mol_acid}/{conf_acid} y/o {mol_base}/{conf_base}. "
            f"Crear runs/{mol_acid}/{SUBDIR_SP_GAS}/{conf_acid}/sp_gas.inp "
            f"y el equivalente para {mol_base} (mismo funcional/base, sin CPCM)."
        )
        return result

    use_born = have_gas and method in ("auto", "born")

    if method == "auto" and not have_gas:
        _warn(
            f"Sin SP fase gas para {mol_acid}/{conf_acid} y/o {mol_base}/{conf_base}. "
            f"Se usa el método directo (solo CPCM) en su lugar: "
            f"ΔG_aq = G_total_CPCM(base) − G_total_CPCM(acid) + G(H⁺,aq). "
            f"Para el ciclo de Born completo, crear los SP fase gas en "
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
            f"G_corr ({gcorr_source}) no disponible para {mol_acid}/{conf_acid}. "
            f"Se usa G ≈ E_SP (sin ZPE/entropía; pKa con mayor incertidumbre).{extra}"
        )
    if tc_base["ok"]:
        G_corr_base = tc_base["G_corr"]
    else:
        extra = f" ({tc_base['warning']})" if tc_base.get("warning") else ""
        _warn(
            f"G_corr ({gcorr_source}) no disponible para {mol_base}/{conf_base}. "
            f"Se usa G ≈ E_SP.{extra}"
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


def _print_born_cycle(r: dict):
    """Imprime las contribuciones detalladas del ciclo de Born."""
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  CICLO DE BORN — {r['mol_acid']} → {r['mol_base']} + H⁺")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. AH  : {r['conf_acid']}")
    print(f"  Conf. A⁻  : {r['conf_base']}")
    print(f"  Nivel SP(AH)  : {r['theory_sp_acid'] or '—'}")
    print(f"  Nivel SP(A⁻)  : {r['theory_sp_base'] or '—'}")
    print()
    print(f"  {'Término':<40} {'Valor':>12}")
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
    print(f"  (incertidumbre típica ±1–2 unidades; Born/CPCM/TZVPD)")
    if r["warnings"]:
        print(f"\n  Advertencias:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


def _print_direct_cycle(r: dict):
    """Imprime las contribuciones del método directo (solo CPCM, sin fase gas)."""
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  MÉTODO DIRECTO (solo CPCM) — {r['mol_acid']} → {r['mol_base']} + H⁺")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. AH  : {r['conf_acid']}")
    print(f"  Conf. A⁻  : {r['conf_base']}")
    print(f"  Nivel SP(AH)  : {r['theory_sp_acid'] or '—'}")
    print(f"  Nivel SP(A⁻)  : {r['theory_sp_base'] or '—'}")
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
    print(f"  (método directo: válido si el error de solvatación del nivel de")
    print(f"  teoría es consistente entre especies; incertidumbre típica ±1–2")
    print(f"  unidades, similar al ciclo de Born completo)")
    if r["warnings"]:
        print(f"\n  Advertencias:")
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
        if pka["method"] == "born":
            print(f"  ΔG_gas        = {pka['dG_gas_kcal']:+.3f} kcal/mol")
            print(f"  ΔG_solv(ácido) = {pka['dG_solv_acid_kcal']:+.3f} kcal/mol")
            print(f"  ΔG_solv(base)  = {pka['dG_solv_base_kcal']:+.3f} kcal/mol")
            print(f"  ΔG_solv(H⁺)   = {DG_SOLV_HPLUS:+.3f} kcal/mol  [Kelly 2006]")
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
    print(f"  siempre usan 'direct' (ver advertencias arriba si pediste 'born').")
    print(f"{sep}\n")


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
):
    """
    Diagrama de distribución de especies (fracciones molares vs pH).
    Guarda como PNG si se provee output_path; si no, muestra en pantalla.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: conda install matplotlib")
        return

    has_dianion = "f_dianion" in df_dist.columns
    pH_vals = df_dist["pH"].values

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = ["#2196F3", "#4CAF50", "#F44336"]  # azul, verde, rojo
    labels = ["AH₂ (q=0)", "AH⁻ (q=−1)", "A²⁻ (q=−2)"]
    cols   = ["f_neutral", "f_anion"] + (["f_dianion"] if has_dianion else [])

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
    """Diagrama de distribución de 5 especies (catión, R, D, anión, dianión) vs pH."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: conda install matplotlib")
        return

    pH_vals = df_dist["pH"].values
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


# ─────────────────────────────────────────────
#  EXPORTACIÓN
# ─────────────────────────────────────────────

def export_results(
    pka_list: list,
    df_dist: pd.DataFrame,
    data_list: list[dict],
    fmt: str,
    taut_result: Optional[dict] = None,
):
    """
    Exporta resultados del ciclo de pKa y distribución de especies.

    pka_list: lista de tuplas (resultado_dict_o_None, etiqueta_str), p. ej.
      [(pka0, "pKa0"), (pka1, "pKa1"), (pka2, "pKa2")]
    taut_result: dict de compute_tautomer_equilibrium(), opcional (se agrega
      como fila/hoja adicional si se provee).
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

    if fmt in ("csv", "both"):
        df_cycle.to_csv(out_dir / "pka_cycle.csv", index=False)
        df_dist.to_csv(out_dir / "protonation_distribution.csv", index=False)
        df_confs.to_csv(out_dir / "conformer_boltzmann.csv", index=False)
        if df_taut is not None:
            df_taut.to_csv(out_dir / "tautomer_equilibrium.csv", index=False)
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
                             "p. ej. 'A'. Requiere también --zwitterion.")
    parser.add_argument("--zwitterion", metavar="DIR",
                        help="Carpeta en runs/ para el zwitterión D — TAUTÓMERO de "
                             "--neutral (mismo q=0, mismos átomos), p. ej. 'D'. "
                             "Requiere también --protonated y --neutral.")
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

    use_extended = args.protonated is not None or args.zwitterion is not None
    if use_extended:
        missing = [name for name, val in
                   [("--protonated", args.protonated), ("--zwitterion", args.zwitterion),
                    ("--neutral", args.neutral)] if val is None]
        if missing:
            parser.error(
                "La cadena extendida (catión + zwitterión) requiere "
                "--protonated, --zwitterion y --neutral simultáneamente. "
                f"Falta(n): {', '.join(missing)}."
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
        print(f"\n▸ Calculando pKa₁  ({args.neutral} → {args.anion} + H⁺)")
        pka1_result = compute_pka(data_neutral, data_anion,
                                  T=T_STD, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if data_anion and data_dianion:
        print(f"\n▸ Calculando pKa₂  ({args.anion} → {args.dianion} + H⁺)")
        pka2_result = compute_pka(data_anion, data_dianion,
                                  T=T_STD, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if pka1_result is None and pka2_result is None:
        print("\n  ⚠  No se pudo calcular ningún pKa. "
              "Revisa los directorios y los SP fase gas.")
        sys.exit(1)

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
    else:
        df_dist = pd.DataFrame()

    data_list = [d for d in [data_neutral, data_anion, data_dianion] if d is not None]

    if args.export and not df_dist.empty:
        export_results([(pka1_result, "pKa1"), (pka2_result, "pKa2")],
                       df_dist, data_list, args.export)

    if args.plot and pKa1 is not None:
        out_dir = WORK_ROOT
        plot_distribution(df_dist, pKa1, pKa2,
                          output_path=out_dir / "species_distribution.png",
                          ph_highlight=args.ph)
        plot_conformer_boltzmann(data_list,
                                 output_path=out_dir / "conformer_boltzmann.png")

    print("\nTerminado.\n")

    return {
        "gcorr_source": gcorr_source,
        "pKa1": pKa1, "pKa2": pKa2,
        "df_dist": df_dist,
        "data_list": data_list,
    }


def _run_extended(args, gcorr_source: str = "harmonic"):
    """
    Cadena de 5 estados:
      A⁺ (catión) ⇌ [R (neutro) ⇌ D (zwitterión)] ⇌ B⁻ (anión) ⇌ C²⁻ (dianión)

    gcorr_source: ver docstring de _run_classic(). Retorna dict con los
    resultados clave (pKa0/1/2, f_R, f_D, df_dist, ...) para _run_comparison().
    """
    print(f"\n▸ Cargando catión protonado: {args.protonated}")
    data_A = collect_mol_data(args.protonated, verbose=args.verbose,
                              gcorr_source=gcorr_source)
    print_mol_summary(data_A, verbose=args.verbose)

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
    print(f"\n▸ Calculando equilibrio de tautómeros  ({args.neutral} ⇌ {args.zwitterion})")
    taut_result = compute_tautomer_equilibrium(data_R, data_D, T=T_STD, verbose=args.verbose,
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

        gR_fallback = _state_G_total_direct(data_R, T=T_STD, gcorr_source=gcorr_source)
        if not gR_fallback.get("ok"):
            print(f"\n  ⚠  Tampoco se pudo obtener G_total para "
                  f"'{args.neutral}'. Abortando.")
            sys.exit(1)

        taut_result = {
            "mol_R": args.neutral, "mol_D": args.zwitterion, "T": T_STD,
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

    print(f"\n▸ Calculando pKa₀  ({args.protonated} → neutro/zwitterión + H⁺)")
    pka0_result = compute_pka(data_A, data_neutral_macro,
                              T=T_STD, verbose=args.verbose, method=args.method,
                              gcorr_source=gcorr_source)

    if data_B:
        print(f"\n▸ Calculando pKa₁  (neutro/zwitterión → {args.anion} + H⁺)")
        pka1_result = compute_pka(data_neutral_macro, data_B,
                                  T=T_STD, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if data_B and data_C:
        print(f"\n▸ Calculando pKa₂  ({args.anion} → {args.dianion} + H⁺)")
        pka2_result = compute_pka(data_B, data_C,
                                  T=T_STD, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

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

    # ── Exportar ─────────────────────────────────────────────────────────
    data_list = [d for d in [data_A, data_R, data_D, data_B, data_C] if d is not None]

    if args.export and not df_dist.empty:
        export_results(
            [(pka0_result, "pKa0"), (pka1_result, "pKa1"), (pka2_result, "pKa2")],
            df_dist, data_list, args.export, taut_result=taut_result,
        )

    # ── Gráficos ─────────────────────────────────────────────────────────
    if args.plot:
        out_dir = WORK_ROOT
        plot_distribution_5(
            df_dist, pKa0, pKa1, pKa2,
            output_path=out_dir / "species_distribution.png",
            ph_highlight=args.ph,
        )
        plot_conformer_boltzmann(
            data_list, output_path=out_dir / "conformer_boltzmann.png",
        )

    print("\nTerminado.\n")

    return {
        "gcorr_source": gcorr_source,
        "pKa0": pKa0, "pKa1": pKa1, "pKa2": pKa2,
        "f_R": f_R, "f_D": f_D,
        "df_dist": df_dist,
        "data_list": data_list,
    }


def _run_comparison(args, use_extended: bool):
    """
    Corre el ciclo COMPLETO dos veces con distinta fuente de G_corr —
    armónico estándar (NumFreq/RRHO) vs. cuasi-armónico (GoodVibes/mRRHO) —
    y compara al final cuánto cambian los pKa's y la distribución de especies.

    Activado con --goodvibes. Usa _run_classic() o _run_extended() según
    corresponda, sin duplicar lógica: solo cambia gcorr_source entre pasadas.
    """
    if shutil.which("goodvibes") is None:
        print("\n  ⚠  No se encontró el binario 'goodvibes' en PATH.")
        print("     Instálalo con: pip install goodvibes --break-system-packages")
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
