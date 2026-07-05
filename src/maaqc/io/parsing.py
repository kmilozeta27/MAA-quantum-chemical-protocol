"""Parsers puros de outputs ORCA/xTB/CREST. Sin side effects, sin logging.

Módulo que elimina la duplicación detectada entre results_collector.py y
protonation_cycle.py. Cada función documenta su origen exacto.

Dos discrepancias reales (no solo duplicaciones) que existían entre
variantes legacy ya se resolvieron explícitamente, con la decisión
documentada en cada función en vez de un merge silencioso:
  - parse_sp_energy ahora exige 'ORCA TERMINATED NORMALLY' (antes solo
    parse_orca_sp_energy lo exigía); parse_orca_sp_energy quedó como
    alias.
  - extract_adc2_lambda ahora es un wrapper de parse_adc2_excitations,
    así que ambas usan la misma prioridad de bloque (LEFT-RIGHT > RIGHT).
"""

from __future__ import annotations

import re
import warnings
from pathlib import Path
from typing import Optional

try:
    from ase.io import read as ase_read
    HAS_ASE = True
except ImportError:
    HAS_ASE = False

HARTREE_TO_KCAL = 627.5094740631
T_STD = 298.15  # K — temperatura de referencia para pKa


def parse_sp_energy(out_file: Path) -> Optional[float]:
    """Extrae la última FINAL SINGLE POINT ENERGY del .out de SP.

    Portado de results_collector.py:136, idéntica carácter por carácter a
    protonation_cycle.py:235. También idéntica en lógica a
    pipeline.py:extract_sp_energy (913-926); se pliega acá como única
    fuente de verdad.

    DECISIÓN (resuelve discrepancia (c) con parse_orca_sp_energy, antes
    sin resolver): ahora EXIGE que 'ORCA TERMINATED NORMALLY' esté en el
    texto antes de aceptar cualquier energía -- si el cálculo no terminó
    normalmente, retorna None en vez de la última energía parcial
    encontrada. Se mantiene el barrido de "última ocurrencia gana" (no el
    "primera ocurrencia + return" de la variante vieja de recovery) porque
    parse_opt_energy depende de este mismo bucle para tomar el paso final
    convergido de un .out de OPT con múltiples 'FINAL SINGLE POINT ENERGY'.

    ⚠️ Implicación para el pipeline en vivo: cualquier llamador que hoy
    lea un .out todavía en progreso (ORCA corriendo, sin terminar) para
    monitoreo/progreso empezará a recibir None en vez de la última energía
    parcial. Si `pipeline.py` tiene ese uso en algún punto (a confirmar al
    migrar `stage_sp`/`stage_opt`), va a necesitar su propio chequeo de
    progreso en vivo separado de este parser.
    """
    if not out_file.exists():
        return None
    try:
        text = out_file.read_text(errors="replace")
    except OSError:
        return None

    if "ORCA TERMINATED NORMALLY" not in text:
        return None

    energy = None
    for line in text.splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                energy = float(line.split()[-1])
            except ValueError:
                pass
    return energy


def parse_opt_energy(out_file: Path) -> Optional[float]:
    """Extrae la energía del último ciclo de optimización (paso final
    convergido). ORCA reporta 'FINAL SINGLE POINT ENERGY' también durante
    Opt, así que tomamos la última ocurrencia.

    Portado de results_collector.py:150 (alias directo de parse_sp_energy).
    """
    return parse_sp_energy(out_file)


def parse_imaginary_freqs(freq_out: Path) -> list[float]:
    """Devuelve lista de frecuencias imaginarias (cm⁻¹, valores negativos).

    Portado de results_collector.py:159.
    """
    imag: list[float] = []
    if not freq_out.exists():
        return imag
    for line in freq_out.read_text(errors="replace").splitlines():
        m = re.match(r"\s+\d+:\s+(-\d+\.\d+)\s+cm", line)
        if m:
            imag.append(float(m.group(1)))
    return imag


def parse_tddft_excitations(tddft_out: Path) -> list[dict]:
    """Parsea excitaciones del bloque 'ABSORPTION SPECTRUM VIA TRANSITION
    ELECTRIC DIPOLE MOMENTS' del output TDDFT de ORCA 6.

    Portado de results_collector.py:171.
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
    """Lee las energías GFN2-xTB embebidas por CREST en el multi-XYZ.

    Portado de results_collector.py:207.
    """
    if not conformers_xyz.exists():
        return []

    energies: list[Optional[float]] = []

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
    """Extrae energía de la línea de comentario de un frame XYZ de CREST.

    Portado de results_collector.py:258.
    """
    comment = comment.strip()
    m = re.search(r"Etot[=\s]+([-\d.]+)", comment)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    m = re.match(r"^\s*([-\d.]+)\s*$", comment)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    m = re.search(r"([-\d.]+)\s*Eh", comment)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    return None


def parse_thermochem(freq_out: Path, T: float = T_STD) -> dict:
    """Extrae las correcciones termoquímicas del bloque THERMOCHEMISTRY
    de ORCA 6.

        G_total_i = E_SP(TZVPD/CPCM) + G_corr(SVP/CPCM)
        G_corr = H_corr - T·S_total

    ⚠️ 'Final entropy term' de ORCA 6 NO es la entropía S en Eh/K: es el
    término T·S ya evaluado a T_parsed, expresado en Eh.

    Retorna dict con: E_elec_freq, ZPE, H_corr, TS_total, S_total, G_corr,
    G_total_freq, T_parsed, n_imag, imag_freqs_cm1, ok.

    Portado de results_collector.py:289. protonation_cycle.py:406 tiene la
    MISMA lógica de cálculo (confirmado línea por línea) pero con una
    diferencia real de forma: su dict de retorno NO incluye la clave
    `imag_freqs_cm1` (solo `n_imag`), y su detección de bloque
    THERMOCHEMISTRY usa `"THERMOCHEMISTRY AT" in line.upper()` en vez del
    regex de acá -- equivalentes en la práctica para outputs reales de
    ORCA 6, pero no literalmente el mismo código. Se adopta esta versión
    (results_collector.py) como fuente de verdad por ser superset estricto;
    protonation_cycle.py debe importar desde acá al migrar.
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

    imag_freqs: list[float] = []
    for line in lines:
        m = re.match(r"\s+\d+:\s+(-\d+\.\d+)\s+cm", line)
        if m:
            imag_freqs.append(float(m.group(1)))

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

    T_parsed: Optional[float] = None
    m = re.search(r"THERMOCHEMISTRY AT\s+([\d.]+)\s*K",
                  lines[thermo_start], re.IGNORECASE)
    if m:
        try:
            T_parsed = float(m.group(1))
        except ValueError:
            pass

    ZPE:          Optional[float] = None
    H_total:      Optional[float] = None
    TS_total:     Optional[float] = None
    G_total_freq: Optional[float] = None

    _pat_zpe = re.compile(r"Zero point energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_h   = re.compile(r"Total [Ee]nthalpy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_s   = re.compile(r"Final entropy term\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_g   = re.compile(r"(?:Total|Final) Gibbs free energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)

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

    H_corr:  Optional[float] = None
    G_corr:  Optional[float] = None
    S_total: Optional[float] = None

    if E_elec_freq is not None and H_total is not None:
        H_corr = H_total - E_elec_freq

    if TS_total is not None and T_parsed is not None and T_parsed != 0:
        S_total = TS_total / T_parsed

    if H_corr is not None and S_total is not None:
        G_corr = H_corr - T * S_total

    if (G_corr is not None and G_total_freq is not None
            and E_elec_freq is not None and T_parsed is not None
            and abs(T - T_parsed) < 1e-6):
        G_reconstructed = E_elec_freq + G_corr
        discrepancy = abs(G_reconstructed - G_total_freq)
        if discrepancy > 1e-5:
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
        "TS_total":       TS_total,
        "S_total":        S_total,
        "G_corr":         G_corr,
        "G_total_freq":   G_total_freq,
        "T_parsed":       T_parsed,
        "n_imag":         len(imag_freqs),
        "imag_freqs_cm1": imag_freqs,
        "ok":             ok,
    }


def parse_gas_phase_sp(sp_gas_out: Path) -> Optional[float]:
    """Extrae la energía de un SP en fase gas (sin CPCM), usado en el
    ciclo de pKa para ΔG_solv de cada especie.

    Portado de results_collector.py:511 (alias directo de parse_sp_energy).
    """
    return parse_sp_energy(sp_gas_out)


def parse_adc2_excitations(adc2_out: Path) -> list[dict]:
    """Extrae la(s) excitación(es) del bloque ADC(2) de un output ORCA 6.

    Prioriza el bloque LEFT-RIGHT (el más equilibrado de los tres que
    imprime ORCA: RIGHT, LEFT, LEFT-RIGHT); si no existe, usa RIGHT.

    Retorna lista de dicts: state, E_eV, E_cm1, wl_nm, fosc, fosc_right,
    fosc_left, source.

    Portado de results_collector.py:2710. Es la única fuente de verdad
    para leer excitaciones ADC(2): extract_adc2_lambda (más abajo) es
    ahora un wrapper delgado sobre esta función, ya no tiene su propio
    parseo del bloque RIGHT.
    """
    if not adc2_out.exists():
        return []

    try:
        text = adc2_out.read_text(errors="replace")
    except OSError:
        return []

    lines = text.splitlines()

    pat_line = re.compile(
        r"^\s*(\d+)-\S+\s*->\s*(\d+)-\S+\s+"
        r"([\d.]+)\s+"
        r"([\d.]+)\s+"
        r"([\d.]+)\s+"
        r"([\d.]+)"
    )

    block_tags = {
        "right":      "SPECTRUM FOR RIGHT TRANSITION MOMENTS",
        "left":       "SPECTRUM FOR LEFT TRANSITION MOMENTS",
        "left_right": "SPECTRUM FOR LEFT-RIGHT TRANSITION MOMENTS",
    }
    abs_header = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"

    collected: dict[str, list[dict]] = {"right": [], "left": [], "left_right": []}

    current_block: Optional[str] = None
    in_abs_block   = False

    for line in lines:
        for key, tag in block_tags.items():
            if tag in line:
                current_block = key
                in_abs_block  = False
                break

        if current_block and not in_abs_block and abs_header in line:
            in_abs_block = True
            continue

        if in_abs_block and current_block:
            m = pat_line.match(line)
            if m:
                root_final = int(m.group(2))
                state_n = root_final if root_final > 0 else 1
                collected[current_block].append({
                    "state":  state_n,
                    "E_eV":   float(m.group(3)),
                    "E_cm1":  float(m.group(4)),
                    "wl_nm":  float(m.group(5)),
                    "fosc":   float(m.group(6)),
                })

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


def parse_theory_level(out_file: Path) -> Optional[str]:
    """Extrae la línea de keywords ORCA ('! funcional base CPCM ...') tal
    como aparece ecoizada al inicio del cálculo. Solo diagnóstico, no se
    usa en ningún cálculo numérico.

    Portado de protonation_cycle.py:252.
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
    """Extrae carga total, multiplicidad y número de átomos de un .out
    de ORCA -- usado para validar consistencia entre estados de
    protonación antes de calcular un pKa.

    Retorna {"charge": int|None, "mult": int|None, "n_atoms": int|None}.

    Portado de protonation_cycle.py:270.
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


def extract_adc2_lambda(out_file: Path) -> Optional[tuple[float, float, float]]:
    """Extrae λmax (nm), E (eV) y fosc de un output ADC(2)/ORCA
    (nroots=1, S0→S1).

    Devuelve (lambda_nm, energy_eV, fosc) del primer estado (S1), o None
    si no se pudo parsear.

    RESUELTO: antes esta función leía únicamente el bloque RIGHT
    (pipeline.py:929 original), mientras que parse_adc2_excitations
    priorizaba LEFT-RIGHT. Ahora es un wrapper delgado sobre
    parse_adc2_excitations, así que ambas usan la misma prioridad de
    bloque (LEFT-RIGHT si existe, si no RIGHT) y el mismo fosc principal
    -- para nroots=1 solo cambia de nombre/forma de retorno (tupla en vez
    de dict), no de fuente de datos.
    """
    excitations = parse_adc2_excitations(out_file)
    if not excitations:
        return None
    exc = excitations[0]  # nroots=1 -> único estado, S0->S1
    return (exc["wl_nm"], exc["E_eV"], exc["fosc"])


def parse_orca_sp_energy(outfile: Path) -> Optional[float]:
    """Alias de parse_sp_energy, mantenido por compatibilidad con el
    nombre usado en recover_prescreening.py.

    RESUELTO: antes esta función era la única de las dos que exigía
    'ORCA TERMINATED NORMALLY'; ahora parse_sp_energy exige exactamente
    lo mismo, así que dejaron de ser dos comportamientos distintos.
    Se conserva como nombre separado solo para no romper el sitio de
    llamada en recover_prescreening.py al migrarlo -- si preferís
    eliminarla del todo y actualizar ese call site para usar
    parse_sp_energy directamente, es un cambio de una línea allá.
    """
    return parse_sp_energy(outfile)


def parse_xtb_total_energy(outfile: Path) -> Optional[float]:
    """Replica el parseo de XtbProc._sp: busca '| TOTAL ENERGY', token[3].

    Portado de recover_prescreening.py.
    """
    if not outfile.exists():
        return None
    try:
        text = outfile.read_text(errors="ignore")
    except OSError:
        return None

    for line in text.splitlines():
        if "| TOTAL ENERGY" in line:
            try:
                return float(line.split()[3])
            except (ValueError, IndexError):
                return None
    return None
