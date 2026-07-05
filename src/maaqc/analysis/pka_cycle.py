"""analysis/pka_cycle.py — Ciclo termodinámico de pKa (fuente de verdad).

Migrado desde protonation_cycle.py. Confirmado con kmilozeta: esta es la
ÚNICA implementación en uso; el compute_pka de results_collector.py
(ciclo de Born puro, sin fallback directo ni _validate_state_charges)
queda descartado — no se migra ni se referencia desde acá.

`compute_pka()` NO es una función standalone: consume los dicts que arma
`collect_mol_data()` (lee 04_sp/, 03_numfreq/, opcionalmente 06_sp_gas/ y
GoodVibes), y usa `_validate_state_charges()` / `_state_G_total_direct()`
como piezas internas. Por eso se migra el núcleo físico completo del
módulo original — no solo las líneas 1174-1509 de compute_pka — excepto:
  - argparse / main() / exportación a xlsx / plots (quedan en
    protonation_cycle.py legacy, fuera de alcance de esta migración).
  - boltzmann_protonation_states() y lo que sigue después de la línea
    ~1600 (distribución de especies por pH): no forma parte de
    compute_pka, se migra aparte si se necesita.

Se incluye compute_tautomer_equilibrium() porque es el único productor de
dicts con is_macrostate=True, que compute_pka() sabe consumir
explícitamente (rama de macroestado combinado neutro+zwitterión).
"""

from __future__ import annotations

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
#  CONSTANTES  (idénticas a protonation_cycle.py)
# ─────────────────────────────────────────────

WORK_ROOT = Path.cwd() / "runs"  # TODO: reemplazar por workdir explícito del Stage al integrar
HARTREE_TO_KCAL = 627.5094740631
R_KCAL = 1.9872041e-3      # kcal/(mol·K)
LN10 = 2.302585093
T_STD = 298.15             # K

# G(H⁺, aq) — dos convenciones NO intercambiables (ver docstring original,
# protonation_cycle.py:160-193):
G_HPLUS_GAS = -6.28        # kcal/mol   Tissandier et al. 1998
DG_SOLV_HPLUS = -265.9     # kcal/mol   Kelly, Cramer & Truhlar 2006
DG_CONC_CORR = 1.894       # kcal/mol   δG_conc 1 atm→1M a 298.15 K
G_HPLUS_AQ_BORN = G_HPLUS_GAS + DG_SOLV_HPLUS          # ≈ -272.18, solo ciclo de Born
G_HPLUS_AQ_DIRECT = -264.61  # kcal/mol  Pliego & Riveros, JPCA 2002, 106, 7434 — método directo (§3.5 del paper)

SUBDIR_SP = "04_sp"
SUBDIR_FREQ = "03_numfreq"
SUBDIR_SP_GAS = "06_sp_gas"

# ── GoodVibes (corrección cuasi-armónica opcional) ────────────────────────
GV_CUTOFF = 100.0
GV_QH = False
GV_QS = "grimme"
GV_CONDA_ENV = "sciencebitch"
_GV_PYTHON_CACHE: Optional[str] = None


# ─────────────────────────────────────────────
#  PARSING DE .out DE ORCA
# ─────────────────────────────────────────────

def parse_sp_energy(out_file: Path) -> Optional[float]:
    """Energía SP (Hartree) desde 'FINAL SINGLE POINT ENERGY'."""
    if not out_file.exists():
        return None
    try:
        for line in out_file.read_text(errors="replace").splitlines():
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    return float(line.split()[-1])
                except ValueError:
                    pass
    except OSError:
        pass
    return None


def parse_theory_level(out_file: Path) -> Optional[str]:
    """Línea de keywords ORCA ('! funcional base CPCM ...'), solo para
    diagnóstico de consistencia entre estados."""
    if not out_file.exists():
        return None
    pat = re.compile(r"^\s*\|\s*1>\s*!\s*(.+)$")
    for line in out_file.read_text(errors="replace").splitlines():
        m = pat.match(line)
        if m:
            return m.group(1).strip()
    return None


def parse_charge_mult_natoms(out_file: Path) -> dict:
    """Carga total, multiplicidad y N° de átomos desde un .out de ORCA.
    Usado por _validate_state_charges() para atrapar un mapeo de carpetas
    desalineado antes de calcular ningún pKa."""
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

    pat_echo = re.compile(r"^\s*\|\s*\d+>\s*\*\s*xyz\s+(-?\d+)\s+(\d+)", re.IGNORECASE)
    for line in lines:
        m = pat_echo.match(line)
        if m:
            charge, mult = int(m.group(1)), int(m.group(2))
            break

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

    for idx, line in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line.upper():
            j = idx + 2
            count = 0
            while j < len(lines) and lines[j].strip():
                count += 1
                j += 1
            n_atoms = count if count > 0 else None
            break

    return {"charge": charge, "mult": mult, "n_atoms": n_atoms}


def _validate_state_charges(
    data_a: dict, data_b: dict, dq_expected: int, dn_expected: int,
    label_a: str, label_b: str,
) -> list[str]:
    """Verifica que data_b difiera de data_a en dq_expected de carga y
    dn_expected átomos. Mensajes con '❌' son fatales (el llamador aborta);
    el resto son informativos."""
    msgs: list[str] = []
    qa, qb = data_a.get("charge"), data_b.get("charge")
    na, nb = data_a.get("n_atoms"), data_b.get("n_atoms")

    if qa is None or qb is None:
        msgs.append(
            f"No se pudo verificar la carga entre '{label_a}' y '{label_b}'. "
            f"Verificar manualmente Δq={dq_expected:+d}."
        )
    elif (qb - qa) != dq_expected:
        msgs.append(
            f"❌ CARGA INCONSISTENTE: '{label_a}' (q={qa:+d}) → '{label_b}' "
            f"(q={qb:+d}) implica Δq={qb - qa:+d}, se esperaba Δq={dq_expected:+d}."
        )

    if na is None or nb is None:
        msgs.append(
            f"No se pudo verificar N° de átomos entre '{label_a}' y '{label_b}'. "
            f"Verificar manualmente ΔN={dn_expected:+d}."
        )
    elif (nb - na) != dn_expected:
        msgs.append(
            f"❌ N° DE ÁTOMOS INCONSISTENTE: '{label_a}' (N={na}) → "
            f"'{label_b}' (N={nb}) implica ΔN={nb - na:+d}, se esperaba "
            f"ΔN={dn_expected:+d}."
        )

    return msgs


# ─────────────────────────────────────────────
#  TERMOQUÍMICA (RRHO armónico + GoodVibes cuasi-armónico)
# ─────────────────────────────────────────────

def parse_thermochem(freq_out: Path, T: float = T_STD) -> dict:
    """Correcciones termoquímicas del bloque THERMOCHEMISTRY de ORCA 6
    (RRHO puro). 'Final entropy term' es T·S, no S puro."""
    null = {
        "E_elec_freq": None, "ZPE": None, "H_corr": None,
        "TS_total": None, "S_total": None, "G_corr": None,
        "G_total_freq": None, "T_parsed": None,
        "n_imag": 0, "ok": False,
    }
    if not freq_out.exists():
        return null
    try:
        lines = freq_out.read_text(errors="replace").splitlines()
    except OSError:
        return null

    n_imag = sum(1 for l in lines if re.match(r"\s+\d+:\s+(-\d+\.\d+)\s+cm", l))

    E_elec_freq: Optional[float] = None
    for line in lines:
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                E_elec_freq = float(line.split()[-1])
                break
            except ValueError:
                pass

    thermo_start = -1
    for idx, line in enumerate(lines):
        if "THERMOCHEMISTRY AT" in line.upper():
            thermo_start = idx

    if thermo_start < 0:
        result = dict(null)
        result["E_elec_freq"] = E_elec_freq
        result["n_imag"] = n_imag
        return result

    T_parsed: Optional[float] = None
    m = re.search(r"THERMOCHEMISTRY AT\s+([\d.]+)\s*K", lines[thermo_start], re.IGNORECASE)
    if m:
        try:
            T_parsed = float(m.group(1))
        except ValueError:
            pass

    _pat_zpe = re.compile(r"Zero point energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_h = re.compile(r"Total [Ee]nthalpy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_s = re.compile(r"Final entropy term\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_g = re.compile(r"(?:Total|Final) Gibbs free energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)

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

    H_corr = (H_total - E_elec_freq) if (H_total is not None and E_elec_freq is not None) else None
    S_total = (TS_total / T_parsed) if (TS_total is not None and T_parsed) else None
    G_corr = (H_corr - T * S_total) if (H_corr is not None and S_total is not None) else None

    ok = all(x is not None for x in [ZPE, H_corr, TS_total, G_corr, G_total_freq])

    return {
        "E_elec_freq": E_elec_freq, "ZPE": ZPE, "H_corr": H_corr,
        "TS_total": TS_total, "S_total": S_total, "G_corr": G_corr,
        "G_total_freq": G_total_freq, "T_parsed": T_parsed,
        "n_imag": n_imag, "ok": ok,
    }


def _resolve_gv_python(force: bool = False) -> str:
    """Ruta al intérprete Python con goodvibes instalado. Cacheada."""
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
    freq_out: Path, T: float = T_STD, cutoff: float = GV_CUTOFF,
    qs: str = GV_QS, qh: bool = GV_QH,
) -> dict:
    """Corre GoodVibes sobre un freq.out y extrae G_corr cuasi-armónico
    (mRRHO de Grimme, opcionalmente Head-Gordon en la entalpía si qh=True)."""
    null = {
        "G_corr": None, "H_corr": None, "S_total": None,
        "scf_energy": None, "qh_gibbs_free_energy": None,
        "n_imag": 0, "ok": False, "warning": None,
    }
    if not freq_out.exists():
        return {**null, "warning": f"No existe {freq_out}"}

    python_exe = _resolve_gv_python()

    try:
        r = subprocess.run([python_exe, "-c", "import goodvibes"], capture_output=True, timeout=10)
        if r.returncode != 0:
            return {**null, "warning": (
                f"No se puede importar goodvibes con {python_exe}. "
                f"conda activate {GV_CONDA_ENV or '<env>'} && pip install goodvibes --break-system-packages"
            )}
    except (subprocess.TimeoutExpired, OSError) as e:
        return {**null, "warning": f"No se pudo verificar goodvibes con {python_exe}: {e}"}

    with tempfile.TemporaryDirectory() as tmp:
        json_path = Path(tmp) / "gv_result.json"
        cmd = [python_exe, "-m", "goodvibes", "--temp", str(T), "-f", str(cutoff), "--qs", qs]
        if qh:
            cmd.append("--qh")
        cmd += ["--json", str(json_path), str(freq_out)]

        try:
            proc = subprocess.run(cmd, cwd=tmp, capture_output=True, text=True, timeout=180)
        except (subprocess.TimeoutExpired, OSError) as e:
            return {**null, "warning": f"Fallo al ejecutar goodvibes: {e}"}

        if not json_path.exists():
            stderr_tail = (proc.stderr or "").strip()[-400:]
            return {**null, "warning": f"goodvibes no generó JSON (returncode={proc.returncode}). {stderr_tail}"}

        try:
            payload = json.loads(json_path.read_text())
        except (OSError, json.JSONDecodeError) as e:
            return {**null, "warning": f"No se pudo leer el JSON de goodvibes: {e}"}

    results = payload.get("results") or []
    if not results:
        return {**null, "warning": "goodvibes no devolvió resultados para este archivo."}

    thermo = results[0].get("thermo") or {}
    scf_energy = thermo.get("scf_energy")
    qh_gibbs = thermo.get("qh_gibbs_free_energy")
    enthalpy = thermo.get("qh_enthalpy") if qh else thermo.get("enthalpy")
    entropy = thermo.get("qh_entropy")
    im_freqs = thermo.get("im_frequency_wn") or []
    n_imag = len(im_freqs) if isinstance(im_freqs, list) else 0

    if scf_energy is None or qh_gibbs is None:
        return {**null, "warning": "goodvibes no reportó scf_energy/qh_gibbs_free_energy."}

    G_corr = qh_gibbs - scf_energy
    H_corr = (enthalpy - scf_energy) if enthalpy is not None else None

    return {
        "G_corr": G_corr, "H_corr": H_corr, "S_total": entropy,
        "scf_energy": scf_energy, "qh_gibbs_free_energy": qh_gibbs,
        "n_imag": n_imag, "ok": True, "warning": None,
    }


def get_thermochem(
    freq_out: Path, T: float = T_STD, source: str = "harmonic",
    gv_cutoff: float = GV_CUTOFF, gv_qs: str = GV_QS, gv_qh: bool = GV_QH,
) -> dict:
    """Wrapper único: source="harmonic" (RRHO/ORCA) o "goodvibes" (mRRHO).
    Mismo shape de salida en ambos casos."""
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
        "E_elec_freq": gv.get("scf_energy"), "ZPE": None,
        "H_corr": gv.get("H_corr"), "TS_total": TS_total,
        "S_total": gv.get("S_total"), "G_corr": gv.get("G_corr"),
        "G_total_freq": gv.get("qh_gibbs_free_energy"), "T_parsed": T,
        "n_imag": gv.get("n_imag", 0), "ok": gv.get("ok", False),
        "source": "goodvibes", "warning": gv.get("warning"),
    }


# ─────────────────────────────────────────────
#  COLECTOR POR MOL (versión ligera para el ciclo de pKa)
# ─────────────────────────────────────────────

TEMPS_K = [250.0, 274.15, 298.15, 310.0]
TEMP_LABELS = ["250K", "274K", "298K", "310K"]


def collect_mol_data(mol_name: str, verbose: bool = False, gcorr_source: str = "harmonic") -> dict:
    """Recolecta datos termodinámicos de todos los confórmeros de una
    molécula: E_sp_aq, E_sp_gas (opcional), G_corr(T), G_total(T). Produce
    el dict que compute_pka() espera como data_acid/data_base."""
    mol_dir = WORK_ROOT / mol_name
    warns: list[str] = []

    if not mol_dir.exists():
        return {"mol_name": mol_name, "ok": False,
                "warnings": [f"Directorio no encontrado: {mol_dir}"],
                "conformers": [], "df": pd.DataFrame()}

    sp_dirs = sorted((mol_dir / SUBDIR_SP).glob("conf_*")) if (mol_dir / SUBDIR_SP).exists() else []
    if not sp_dirs:
        return {"mol_name": mol_name, "ok": False,
                "warnings": [f"Sin confórmeros SP en {mol_dir / SUBDIR_SP}"],
                "conformers": [], "df": pd.DataFrame()}

    rows: list[dict] = []
    for sp_dir in sp_dirs:
        conf_id = sp_dir.name
        sp_aq_file = sp_dir / "sp.out"
        e_sp_aq = parse_sp_energy(sp_aq_file)

        sp_gas_file = mol_dir / SUBDIR_SP_GAS / conf_id / "sp_gas.out"
        e_sp_gas = parse_sp_energy(sp_gas_file)

        freq_file = mol_dir / SUBDIR_FREQ / conf_id / "freq.out"
        tc = get_thermochem(freq_file, T=T_STD, source=gcorr_source,
                             gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)
        if gcorr_source == "goodvibes" and not tc["ok"] and tc.get("warning"):
            warns.append(f"[{mol_name}/{conf_id}] GoodVibes: {tc['warning']}")

        cmn = parse_charge_mult_natoms(sp_aq_file)
        if cmn["charge"] is None and cmn["n_atoms"] is None:
            cmn = parse_charge_mult_natoms(freq_file)

        row: dict = {
            "conformer": conf_id, "E_sp_aq": e_sp_aq, "E_sp_gas": e_sp_gas,
            "thermo_ok": tc["ok"], "n_imag": tc["n_imag"],
            "H_corr": tc["H_corr"], "S_total_EhK": tc["S_total"], "ZPE": tc["ZPE"],
            "charge": cmn["charge"], "mult": cmn["mult"], "n_atoms": cmn["n_atoms"],
        }

        for T, lbl in zip(TEMPS_K, TEMP_LABELS):
            if e_sp_aq is not None and tc["ok"]:
                G_corr_T = tc["H_corr"] - T * tc["S_total"]
                row[f"G_total_{lbl}"] = e_sp_aq + G_corr_T
            else:
                row[f"G_total_{lbl}"] = None

        rows.append(row)

    df = pd.DataFrame(rows)

    G_boltz: dict[str, Optional[float]] = {}
    for T, lbl in zip(TEMPS_K, TEMP_LABELS):
        col = f"G_total_{lbl}"
        valid = df[col].notna() & df["thermo_ok"]
        if valid.any():
            G_vals = df.loc[valid, col].values.astype(float)
            G_min = G_vals.min()
            dG_kcal = (G_vals - G_min) * HARTREE_TO_KCAL
            weights = np.exp(-dG_kcal / (R_KCAL * T))
            weights /= weights.sum()
            df.loc[valid, f"boltz_{lbl}"] = weights
            df.loc[valid, f"dG_{lbl}_kcal"] = dG_kcal
            ln_Z = np.log(np.sum(np.exp(-dG_kcal / (R_KCAL * T))))
            G_boltz[lbl] = G_min - (R_KCAL * T * ln_Z) / HARTREE_TO_KCAL
        else:
            G_boltz[lbl] = None
            if not df["thermo_ok"].any():
                warns.append(f"[{mol_name}] Sin correcciones térmicas válidas.")

    col_298 = "G_total_298K"
    best_conf = None
    if col_298 in df.columns:
        valid = df[col_298].notna() & df["thermo_ok"]
        if valid.any():
            best_conf = df.loc[df.loc[valid, col_298].idxmin(), "conformer"]

    n_imag_any = (df["n_imag"] > 0).sum()
    if n_imag_any > 0:
        warns.append(f"[{mol_name}] {n_imag_any} confórmero(s) con frecuencias imaginarias.")

    n_missing_gas = df["E_sp_gas"].isna().sum()
    if n_missing_gas > 0:
        warns.append(
            f"[{mol_name}] {n_missing_gas}/{len(df)} confórmeros sin SP fase gas "
            f"(esperado en runs/{mol_name}/{SUBDIR_SP_GAS}/conf_NNN/sp_gas.out)."
        )

    charges_seen = sorted(set(df["charge"].dropna().astype(int)))
    mults_seen = sorted(set(df["mult"].dropna().astype(int)))
    natoms_seen = sorted(set(df["n_atoms"].dropna().astype(int)))

    state_charge = charges_seen[0] if len(charges_seen) == 1 else None
    state_mult = mults_seen[0] if len(mults_seen) == 1 else None
    state_n_atoms = natoms_seen[0] if len(natoms_seen) == 1 else None

    if len(charges_seen) > 1:
        warns.append(f"[{mol_name}] Carga inconsistente entre confórmeros: {charges_seen}.")
    if len(natoms_seen) > 1:
        warns.append(f"[{mol_name}] N° de átomos inconsistente entre confórmeros: {natoms_seen}.")
    if not charges_seen or not natoms_seen:
        warns.append(f"[{mol_name}] No se pudo determinar carga y/o N° de átomos.")

    return {
        "mol_name": mol_name, "conformers": rows, "df": df,
        "best_conf": best_conf, "G_boltz": G_boltz,
        "charge": state_charge, "mult": state_mult, "n_atoms": state_n_atoms,
        "ok": len(rows) > 0, "warnings": warns, "gcorr_source": gcorr_source,
    }


# ─────────────────────────────────────────────
#  ENERGÍA LIBRE REPRESENTATIVA DE UN ESTADO (método directo)
# ─────────────────────────────────────────────

def _state_G_total_direct(data: dict, T: float = T_STD, gcorr_source: str = "harmonic") -> dict:
    """G_total_CPCM del confórmero representativo, para el método directo.
    Soporta dicts normales (de collect_mol_data) y macroestados combinados
    (data["is_macrostate"]=True, de compute_tautomer_equilibrium)."""
    if data.get("is_macrostate"):
        return {
            "G_total_Eh": data["G_total_direct_Eh"],
            "conf": data.get("conf_label", "macroestado"),
            "theory_sp": data.get("theory_label"),
            "theory_freq": data.get("theory_label"),
            "ok": True, "warn": None,
        }

    mol = data["mol_name"]
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
        warn = f"G_corr ({gcorr_source}) no disponible para {mol}/{conf}.{extra}"

    G_total = E_aq + (G_corr if G_corr is not None else 0.0)
    theory_sp = parse_theory_level(WORK_ROOT / mol / SUBDIR_SP / conf / "sp.out")
    theory_freq = parse_theory_level(freq_file)

    return {
        "G_total_Eh": G_total, "conf": conf, "theory_sp": theory_sp,
        "theory_freq": theory_freq, "ok": True, "warn": warn,
    }


# ─────────────────────────────────────────────
#  EQUILIBRIO DE TAUTÓMEROS (único productor de macroestados combinados)
# ─────────────────────────────────────────────

def compute_tautomer_equilibrium(
    data_R: dict, data_D: dict, T: float = T_STD,
    verbose: bool = False, gcorr_source: str = "harmonic",
) -> dict:
    """Compara neutro (R) y zwitterión (D) como tautómeros (mismo q, mismo
    N° de átomos) y construye el macroestado combinado que compute_pka()
    usa como 'neutro' frente a catión/anión."""
    mol_R, mol_D = data_R["mol_name"], data_D["mol_name"]
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

    charge_msgs = _validate_state_charges(data_R, data_D, dq_expected=0, dn_expected=0,
                                           label_a=mol_R, label_b=mol_D)
    for m in charge_msgs:
        _warn(m)
    if any(m.startswith("❌") for m in charge_msgs):
        _warn(f"Abortando equilibrio de tautómeros ({mol_R}⇌{mol_D}) por inconsistencia de carga/N° de átomos.")
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
        _warn(f"NIVEL DE TEORÍA INCONSISTENTE entre {mol_R} y {mol_D}.")

    G_R, G_D = gR["G_total_Eh"], gD["G_total_Eh"]
    result["conf_R"], result["conf_D"] = gR["conf"], gD["conf"]
    result["G_R_Eh"], result["G_D_Eh"] = G_R, G_D
    result["theory_R"], result["theory_D"] = gR["theory_sp"], gD["theory_sp"]

    dG_taut_kcal = (G_D - G_R) * HARTREE_TO_KCAL
    K_taut = np.exp(-dG_taut_kcal / (R_KCAL * T))
    f_R = 1.0 / (1.0 + K_taut)
    f_D = K_taut / (1.0 + K_taut)

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
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  EQUILIBRIO DE TAUTÓMEROS (sin H⁺) — {r['mol_R']} ⇌ {r['mol_D']}")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. R (neutro)      : {r['conf_R']}")
    print(f"  Conf. D (zwitterión)  : {r['conf_D']}")
    print(f"  {'G_total_CPCM(R)  [Eh]':<40} {r['G_R_Eh']:>14.8f}")
    print(f"  {'G_total_CPCM(D)  [Eh]':<40} {r['G_D_Eh']:>14.8f}")
    print(f"  {'ΔG_taut = G(D) − G(R)  [kcal/mol]':<40} {r['dG_taut_kcal']:>+14.3f}")
    print(f"  f(R)={r['f_R']*100:6.2f} %   f(D)={r['f_D']*100:6.2f} %")
    print(f"  G_macroestado neutro (Boltzmann R+D) = {r['G_macro_Eh']:.8f} Eh")
    if r["warnings"]:
        print("  Advertencias:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


# ─────────────────────────────────────────────
#  CICLO TERMODINÁMICO DE pKa  (fuente de verdad)
# ─────────────────────────────────────────────

def compute_pka(
    data_acid: dict, data_base: dict, T: float = T_STD, verbose: bool = False,
    method: str = "auto", gcorr_source: str = "harmonic",
) -> dict:
    """Calcula el pKa de AH → A⁻ + H⁺.

    method="born"   — ciclo de Born completo (Ho & Coote 2010), requiere
                       SP fase gas (06_sp_gas) además del SP acuoso (CPCM).
    method="direct" — solo CPCM, usa G(H⁺,aq) de Pliego & Riveros 2002
                       (G_HPLUS_AQ_DIRECT), consistente con §3.5 del paper.
    method="auto"   — "born" si hay SP gas para ambos estados; si no, cae
                       a "direct" automáticamente (registrado en
                       result["method"] / result["warnings"]).
    """
    mol_acid = data_acid["mol_name"]
    mol_base = data_base["mol_name"]
    warns: list[str] = []
    result: dict = {
        "mol_acid": mol_acid, "mol_base": mol_base,
        "T": T, "pKa": None, "ok": False, "warnings": warns, "method": None,
        "dG_gas_kcal": None, "dG_solv_acid_kcal": None,
        "dG_solv_base_kcal": None, "dG_solv_hplus_kcal": DG_SOLV_HPLUS,
        "dG_aq_kcal": None,
        "G_total_acid_Eh": None, "G_total_base_Eh": None, "G_hplus_aq_kcal": None,
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

    charge_msgs = _validate_state_charges(
        data_acid, data_base, dq_expected=-1, dn_expected=-1,
        label_a=mol_acid, label_b=mol_base,
    )
    for m in charge_msgs:
        _warn(m)
    if any(m.startswith("❌") for m in charge_msgs):
        _warn(f"Abortando cálculo de pKa ({mol_acid}→{mol_base}) por inconsistencia de carga/N° de átomos.")
        return result

    # ── Rama macroestado combinado (neutro+zwitterión) ──────────────────
    if data_acid.get("is_macrostate") or data_base.get("is_macrostate"):
        if method == "born":
            _warn(
                f"method='born' solicitado, pero uno de los estados es un "
                f"macroestado combinado. El ciclo de Born no está definido "
                f"para una mezcla; se usa el método directo."
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

        if g_acid["theory_sp"] and g_base["theory_sp"] and g_acid["theory_sp"] != g_base["theory_sp"]:
            _warn(f"NIVEL DE TEORÍA INCONSISTENTE entre {mol_acid} y {mol_base}.")

        G_total_acid, G_total_base = g_acid["G_total_Eh"], g_base["G_total_Eh"]
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

    # ── Confórmero representativo de cada estado ────────────────────────
    conf_acid = data_acid.get("best_conf")
    conf_base = data_base.get("best_conf")
    result["conf_acid"], result["conf_base"] = conf_acid, conf_base

    if conf_acid is None:
        _warn(f"No hay confórmero representativo para {mol_acid}.")
        return result
    if conf_base is None:
        _warn(f"No hay confórmero representativo para {mol_base}.")
        return result

    df_acid, df_base = data_acid["df"], data_base["df"]
    row_acid = df_acid[df_acid["conformer"] == conf_acid].iloc[0]
    row_base = df_base[df_base["conformer"] == conf_base].iloc[0]

    sp_acid_file = WORK_ROOT / mol_acid / SUBDIR_SP / conf_acid / "sp.out"
    sp_base_file = WORK_ROOT / mol_base / SUBDIR_SP / conf_base / "sp.out"
    freq_acid_lvl = parse_theory_level(WORK_ROOT / mol_acid / SUBDIR_FREQ / conf_acid / "freq.out")
    freq_base_lvl = parse_theory_level(WORK_ROOT / mol_base / SUBDIR_FREQ / conf_base / "freq.out")
    sp_acid_lvl = parse_theory_level(sp_acid_file)
    sp_base_lvl = parse_theory_level(sp_base_file)

    result["theory_sp_acid"], result["theory_sp_base"] = sp_acid_lvl, sp_base_lvl
    result["theory_freq_acid"], result["theory_freq_base"] = freq_acid_lvl, freq_base_lvl

    if sp_acid_lvl and sp_base_lvl and sp_acid_lvl != sp_base_lvl:
        _warn(f"NIVEL DE TEORÍA INCONSISTENTE en el SP acuoso entre {mol_acid} y {mol_base}.")
    if freq_acid_lvl and freq_base_lvl and freq_acid_lvl != freq_base_lvl:
        _warn(f"NIVEL DE TEORÍA INCONSISTENTE en NumFreq entre {mol_acid} y {mol_base}.")

    E_aq_acid, E_aq_base = row_acid["E_sp_aq"], row_base["E_sp_aq"]
    result["E_sp_aq_acid_Eh"], result["E_sp_aq_base_Eh"] = E_aq_acid, E_aq_base

    if E_aq_acid is None or pd.isna(E_aq_acid):
        _warn(f"Sin SP acuoso para {mol_acid}/{conf_acid}.")
        return result
    if E_aq_base is None or pd.isna(E_aq_base):
        _warn(f"Sin SP acuoso para {mol_base}/{conf_base}.")
        return result

    E_gas_acid, E_gas_base = row_acid["E_sp_gas"], row_base["E_sp_gas"]
    result["E_sp_gas_acid_Eh"], result["E_sp_gas_base_Eh"] = E_gas_acid, E_gas_base

    have_gas = not (
        E_gas_acid is None or pd.isna(E_gas_acid) or
        E_gas_base is None or pd.isna(E_gas_base)
    )

    if method == "born" and not have_gas:
        _warn(
            f"method='born' solicitado pero faltan SP fase gas para "
            f"{mol_acid}/{conf_acid} y/o {mol_base}/{conf_base}. Crear "
            f"runs/{mol_acid}/{SUBDIR_SP_GAS}/{conf_acid}/sp_gas.inp y el "
            f"equivalente para {mol_base}."
        )
        return result

    use_born = have_gas and method in ("auto", "born")

    if method == "auto" and not have_gas:
        _warn(
            f"Sin SP fase gas para {mol_acid}/{conf_acid} y/o "
            f"{mol_base}/{conf_base}. Se usa el método directo (solo CPCM)."
        )

    freq_acid = WORK_ROOT / mol_acid / SUBDIR_FREQ / conf_acid / "freq.out"
    freq_base = WORK_ROOT / mol_base / SUBDIR_FREQ / conf_base / "freq.out"
    tc_acid = get_thermochem(freq_acid, T=T, source=gcorr_source,
                              gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)
    tc_base = get_thermochem(freq_base, T=T, source=gcorr_source,
                              gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)

    G_corr_acid = tc_acid["G_corr"] if tc_acid["ok"] else None
    if not tc_acid["ok"]:
        extra = f" ({tc_acid['warning']})" if tc_acid.get("warning") else ""
        _warn(f"G_corr ({gcorr_source}) no disponible para {mol_acid}/{conf_acid}.{extra}")

    G_corr_base = tc_base["G_corr"] if tc_base["ok"] else None
    if not tc_base["ok"]:
        extra = f" ({tc_base['warning']})" if tc_base.get("warning") else ""
        _warn(f"G_corr ({gcorr_source}) no disponible para {mol_base}/{conf_base}.{extra}")

    result["G_corr_acid_Eh"], result["G_corr_base_Eh"] = G_corr_acid, G_corr_base

    if use_born:
        G_gas_acid = E_gas_acid + (G_corr_acid if G_corr_acid is not None else 0.0)
        G_gas_base = E_gas_base + (G_corr_base if G_corr_base is not None else 0.0)
        result["G_gas_acid_Eh"], result["G_gas_base_Eh"] = G_gas_acid, G_gas_base

        dG_gas_kcal = (G_gas_base - G_gas_acid) * HARTREE_TO_KCAL + G_HPLUS_GAS
        result["dG_gas_kcal"] = dG_gas_kcal

        dG_solv_acid = (E_aq_acid - E_gas_acid) * HARTREE_TO_KCAL + DG_CONC_CORR
        dG_solv_base = (E_aq_base - E_gas_base) * HARTREE_TO_KCAL + DG_CONC_CORR
        result["dG_solv_acid_kcal"], result["dG_solv_base_kcal"] = dG_solv_acid, dG_solv_base

        dG_aq_kcal = dG_gas_kcal + dG_solv_base + DG_SOLV_HPLUS - dG_solv_acid
        result["dG_aq_kcal"] = dG_aq_kcal
        result["method"] = "born"
    else:
        G_total_acid = E_aq_acid + (G_corr_acid if G_corr_acid is not None else 0.0)
        G_total_base = E_aq_base + (G_corr_base if G_corr_base is not None else 0.0)
        result["G_total_acid_Eh"], result["G_total_base_Eh"] = G_total_acid, G_total_base
        result["G_hplus_aq_kcal"] = G_HPLUS_AQ_DIRECT

        dG_aq_kcal = (G_total_base - G_total_acid) * HARTREE_TO_KCAL + G_HPLUS_AQ_DIRECT
        result["dG_aq_kcal"] = dG_aq_kcal
        result["method"] = "direct"

    pKa = result["dG_aq_kcal"] / (R_KCAL * T * LN10)
    result["pKa"] = pKa
    result["ok"] = True

    if verbose:
        if result["method"] == "born":
            _print_born_cycle(result)
        else:
            _print_direct_cycle(result)

    return result


def _print_born_cycle(r: dict):
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  CICLO DE BORN — {r['mol_acid']} → {r['mol_base']} + H⁺")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. AH  : {r['conf_acid']}   Conf. A⁻  : {r['conf_base']}")

    def fmt(v): return f"{v:+12.3f}" if v is not None else "       N/A"

    print(f"  {'ΔG_gas  [kcal/mol]':<40} {fmt(r['dG_gas_kcal'])}")
    print(f"  {'ΔG_solv(AH)  [kcal/mol]':<40} {fmt(r['dG_solv_acid_kcal'])}")
    print(f"  {'ΔG_solv(A⁻)  [kcal/mol]':<40} {fmt(r['dG_solv_base_kcal'])}")
    print(f"  {'ΔG_solv(H⁺) [Kelly 2006]  [kcal/mol]':<40} {DG_SOLV_HPLUS:>12.3f}")
    print(f"  {'ΔG_aq  [kcal/mol]':<40} {fmt(r['dG_aq_kcal'])}")
    pka_str = f"{r['pKa']:.2f}" if r['pKa'] is not None else "N/A"
    print(f"  pKa = ΔG_aq / (RT·ln10)  =  {pka_str}  (incertidumbre típica ±1–2)")
    if r["warnings"]:
        print("  Advertencias:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


def _print_direct_cycle(r: dict):
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  MÉTODO DIRECTO (solo CPCM) — {r['mol_acid']} → {r['mol_base']} + H⁺")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. AH  : {r['conf_acid']}   Conf. A⁻  : {r['conf_base']}")
    print(f"  {'G(H⁺, aq)  [kcal/mol]':<40} {r['G_hplus_aq_kcal']:>12.3f}"
          f"   [Pliego & Riveros 2002, JPCA 106, 7434]")
    dG = r["dG_aq_kcal"]
    print(f"  {'ΔG_aq  [kcal/mol]':<40} {dG:+12.3f}" if dG is not None else "")
    pka_str = f"{r['pKa']:.2f}" if r['pKa'] is not None else "N/A"
    print(f"  pKa = ΔG_aq / (RT·ln10)  =  {pka_str}  (incertidumbre típica ±1–2)")
    if r["warnings"]:
        print("  Advertencias:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")
