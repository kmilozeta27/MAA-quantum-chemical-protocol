"""Wrapper de GoodVibes (correcciones termoquímicas alternativas a mano).

No estaba contemplado como módulo propio en el esqueleto original de la
skill maa-repackaging (que solo listaba boltzmann/pka_cycle/monte_carlo/
parse_uvvis bajo analysis/) — se agrega acá porque el legacy tiene
suficiente lógica propia de GoodVibes como para no mezclarla con
results_collector.py ni con pka_cycle.py.

FUENTE CANÓNICA (Issue C1, primera mitad, ya migrada): `_resolve_gv_python()`
y `run_goodvibes()` vienen de `pka_cycle.py:282-397` — es la versión que
ya estaba en uso activo para el ciclo de pKa (llamada indirectamente vía
`get_thermochem(..., source="goodvibes")`). `pka_cycle.py` importa estas
dos funciones desde acá en vez de definirlas localmente.

Issue C1, segunda mitad (esta migración) — RESUELTO:
Había 3 candidatas para la parte "a nivel de molécula/appendix" de la
integración con GoodVibes:

  1. `pka_cycle.py:282-397` (fuente canónica, ver arriba) — opera sobre
     UN freq.out, invoca `python -m goodvibes --json <path> <freq_out>`
     (flag `--json` de GoodVibes v4+ apuntando a un archivo explícito en
     un tempdir), y devuelve un dict simplificado
     (G_corr/H_corr/S_total/scf_energy/qh_gibbs_free_energy/n_imag/ok).
     Pensada para consumirse una vez por confórmero desde el ciclo de pKa.

  2. `protonation_cycle.py:_resolve_gv_python()` (525-587) y
     `run_goodvibes()` (588-708) — comparación línea por línea contra (1):
     MISMA lógica, mismo orden de prioridad de resolución de intérprete,
     mismo flag `--json` + tempdir, mismo dict de retorno. Las únicas
     diferencias son cosméticas (docstrings más largos explicando el
     porqué de `python -m goodvibes` en vez de `shutil.which`, y algunos
     mensajes de warning con más contexto). Es, en la práctica, una
     copia-pega de (1) con comentarios extra, no una función distinta.
     → DESCARTADA explícitamente: no se migra nada de acá. La versión
     canónica (1), ya migrada en este archivo, cubre el mismo caso de uso.
     `protonation_cycle.py` puede apuntar a importar desde acá si en algún
     momento se quiere deduplicar esa copia en el legacy.

  3. `results_collector.py:_resolve_gv_python()` (1719-1772) y
     `run_goodvibes_on_mol()` (1981-2181) — a diferencia de (1) y (2),
     ESTA NO es una copia: opera sobre TODOS los confórmeros de una
     molécula a la vez (`03_numfreq/conf_*/freq.out`), prueba capacidades
     de la versión de GoodVibes instalada (`_probe_gv_capabilities`,
     JSON vs .dat tabulado según versión), y parsea ambos formatos
     posibles (`_parse_goodvibes_json` para v4+, `_parse_goodvibes_dat`
     como fallback para v3.x o cuando el JSON no se generó). Es la base
     del "apéndice GoodVibes" del paper: comparación G_corr/S ORCA-RRHO
     vs GV-mRRHO por confórmero, con Δ, rankings y pesos Boltzmann
     (`compute_gv_appendix`, `print_gv_appendix`, `export_gv_appendix`,
     más `_gv_version` como utilidad compartida). Su `_resolve_gv_python`
     propia (1719-1772) sí es funcionalmente idéntica a la canónica (1) en
     lógica (mismo orden de prioridad) pero le falta el parámetro `force`
     y el cacheo en `_GV_PYTHON_CACHE` — se descarta también, sin migrar,
     y todo lo de abajo reutiliza la `_resolve_gv_python()` canónica de
     este mismo archivo.
     → MIGRADA tal cual (con `WORK_ROOT` reemplazado por `runs_root`
     explícito, ver Issue C0, y `GV_FREQ_CUTOFF`/`GV_QS` unificados con
     las constantes ya existentes acá — ver notas en la sección de
     constantes más abajo).

TODO pendiente (fuera de alcance de esta migración): `protonation_cycle.py`
también tiene `_print_goodvibes_comparison()` (2721-2795), que compara
pKa's armónicos vs GoodVibes a nivel de ciclo de pKa completo (no de
apéndice por confórmero). Es conceptualmente más cercana a `pka_cycle.py`
que a este módulo — si se migra, debería vivir ahí o en un módulo de
reporting de pKa, no acá.
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

from ..constants import HARTREE_TO_KCAL, KB_KCAL, T_STD

# ─────────────────────────────────────────────
#  CONSTANTES PROPIAS DE GOODVIBES
#  T_STD/HARTREE_TO_KCAL/KB_KCAL ya no se definen acá (vivían duplicadas,
#  idénticas a las de pka_cycle.py/parsing.py) — se importan desde
#  ..constants, fuente única.
#
#  GV_CUTOFF: results_collector.py tenía esta misma constante bajo el
#  nombre `GV_FREQ_CUTOFF` (mismo valor, 100.0 cm⁻¹) — se unifica al
#  nombre que ya existía acá (`GV_CUTOFF`) en vez de mantener dos
#  constantes idénticas con nombres distintos; `GV_FREQ_CUTOFF` no
#  sobrevive a esta migración.
#  GV_QS: mismo nombre y mismo valor ("grimme") en las 3 fuentes — sin
#  cambios.
#  GV_SUBDIR: nueva acá (no existía antes de esta migración) — viene tal
#  cual de results_collector.py, subcarpeta creada dentro de cada
#  conf_xxx/ para los archivos de salida de GoodVibes.
#  SUBDIR_FREQ: duplica intencionalmente el valor de SUBDIR_FREQ definido
#  en pka_cycle.py ("03_numfreq") — no se importa desde ahí para evitar un
#  ciclo de imports (pka_cycle.py ya importa DE este módulo).
# ─────────────────────────────────────────────

GV_CUTOFF = 100.0
GV_QH = False
GV_QS = "grimme"
GV_CONDA_ENV = "sciencebitch"
GV_SUBDIR = "goodvibes"
SUBDIR_FREQ = "03_numfreq"
_GV_PYTHON_CACHE: Optional[str] = None


def _resolve_gv_python(force: bool = False) -> str:
    """Ruta al intérprete Python con goodvibes instalado. Cacheada.

    Migrado desde pka_cycle.py:282-323 (fuente canónica elegida para esta
    migración porque era la versión efectivamente en uso por el ciclo de
    pKa). Sin cambios de lógica respecto al original. Las variantes de
    `protonation_cycle.py` (525-587) y `results_collector.py` (1719-1772)
    son funcionalmente equivalentes a esta y se descartaron explícitamente
    (ver docstring del módulo, Issue C1) — todo el resto de este archivo,
    incluidas las funciones a nivel de molécula migradas más abajo, reusa
    esta única implementación.
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
    freq_out: Path, T: float = T_STD, cutoff: float = GV_CUTOFF,
    qs: str = GV_QS, qh: bool = GV_QH,
) -> dict:
    """Corre GoodVibes sobre un freq.out y extrae G_corr cuasi-armónico
    (mRRHO de Grimme, opcionalmente Head-Gordon en la entalpía si qh=True).

    Migrado desde pka_cycle.py:326-395, sin cambios de lógica. Es la
    función que `pka_cycle.get_thermochem(..., source="goodvibes")` llama
    -- pka_cycle.py importa esto desde acá en vez de tenerlo local.

    Opera sobre UN freq.out. Para procesar todos los confórmeros de una
    molécula de una sola vez (y comparar contra las correcciones RRHO de
    ORCA), ver `run_goodvibes_on_mol()` más abajo — es una función
    distinta con su propio formato de invocación a GoodVibes, no un
    wrapper de esta.
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


# ═══════════════════════════════════════════════════════════════════════
#  APÉNDICE GOODVIBES POR MOLÉCULA (Issue C1, segunda mitad)
#  Migrado desde results_collector.py:1773-2709. A diferencia de
#  run_goodvibes() (arriba), estas funciones operan sobre TODOS los
#  confórmeros de una molécula y construyen la comparación ORCA-RRHO vs
#  GV-mRRHO que se reporta en el paper como apéndice.
# ═══════════════════════════════════════════════════════════════════════


def _probe_gv_capabilities(python_exe: str) -> dict:
    """
    Ejecuta `python -m goodvibes --help` para detectar qué flags están
    disponibles en la versión instalada.

    Retorna un dict con:
        supports_json_format  : bool — existe --output_format json (v4+)
        supports_csv          : bool — existe --csv (v3.x)
        version               : str  — número de versión o "?"

    Migrado desde results_collector.py:1773-1811, sin cambios de lógica.
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

    Migrado desde results_collector.py:1812-1978, sin cambios de lógica.
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
        if len(parts[0]) == 1 and not parts[0].lstrip('-').replace('.', '', 1).isdigit():
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

    Migrado desde results_collector.py:2182-2281, sin cambios de lógica.
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
            f"Claves disponibles: {[m.get('filename', '?') for m in molecules]}"
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


def run_goodvibes_on_mol(
    mol_name: str, runs_root: Path, T: float = T_STD,
) -> dict[str, dict]:
    """
    Invoca GoodVibes sobre todos los freq.out de <runs_root>/<mol_name>/
    03_numfreq/<conf_*>/freq.out.

    Flujo para cada confórmero:
      1. Resuelve el Python correcto (entorno sciencebitch o auto-detección),
         reusando la `_resolve_gv_python()` canónica de este módulo.
      2. Detecta las capacidades de la versión instalada (--output_format json
         vs parseo del .dat de salida clásico).
      3. Crea <runs_root>/<mol_name>/03_numfreq/<conf>/goodvibes/ y ejecuta
         GoodVibes desde allí.
      4. Parsea JSON (v4+) o el .dat tabulado (v3.x) según disponibilidad.

    Retorna dict:
        { "conf_001": {"H_GV": float, "S_GV": float, "G_GV": float,
                       "ZPE_GV": float, "T_GV": float, "ok": bool,
                       "json_path": Path, "warning": str}, ... }

    Si GoodVibes no está disponible o falla para un confórmero, "ok" = False
    y "warning" describe el error — el resto del script no se interrumpe.

    Migrado desde results_collector.py:1981-2179 (Issue C1, segunda mitad).
    Cambios respecto al original, ninguno silencioso:
      - `WORK_ROOT` (global del legacy) reemplazado por `runs_root: Path`
        explícito (Issue C0, mismo patrón que pka_cycle.py).
      - `GV_FREQ_CUTOFF`/`GV_QS` locales del legacy unificados con las
        constantes ya existentes en este módulo (`GV_CUTOFF`, `GV_QS`).
      - Temperatura hardcodeada `"-t", "298.15"` en el comando del legacy
        ahora es el parámetro `T: float = T_STD`, siguiendo el mismo
        patrón que `run_goodvibes()` (arriba) en vez de hardcodear.
      - `_resolve_gv_python()` propia del legacy descartada; se usa la
        canónica del módulo (ver docstring del módulo, Issue C1).
    """
    mol_dir  = runs_root / mol_name
    freq_dir = mol_dir / SUBDIR_FREQ

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
                f"  Instalar con: conda activate {GV_CONDA_ENV} && pip install goodvibes\n"
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
            "-f", str(GV_CUTOFF),
            "-t", str(T),
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
                f"  Instalar con: conda activate {GV_CONDA_ENV} && pip install goodvibes"
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


def compute_gv_appendix(
    mol_name: str,
    runs_root: Path,
    df_conf: pd.DataFrame,
    gv_results: dict[str, dict],
) -> pd.DataFrame:
    """
    Construye el DataFrame del apéndice GoodVibes, comparando G_ORCA vs G_GV
    para cada confórmero a todas las temperaturas.

    Requiere que df_conf tenga las columnas producidas por compute_g_total()
    (G_total_{label}, dG_{label}_kcal, boltz_{label}, H_corr, S_total_EhK).
    Si compute_g_total() no fue llamado previamente (el flujo normal de
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

    Migrado desde results_collector.py:2284-2465 (Issue C1, segunda mitad).
    `WORK_ROOT` reemplazado por `runs_root: Path` explícito (Issue C0).

    Migrado desde results_collector.py:2284-2465 (Issue C1, segunda mitad).
    `WORK_ROOT` reemplazado por `runs_root: Path` explícito (Issue C0).

    `compute_g_total()` y `compute_boltzmann_multi_t()` viven en
    `analysis/boltzmann.py` (migración C2, que las renombró de
    `compute_G_total`/`compute_boltzmann_multiT` a snake_case). Fix
    (2026-07-08): este import usaba todavía los nombres camelCase viejos
    y rompía con ImportError en cuanto se llegaba a esta rama —
    detectado por
    `tests/test_analysis/test_goodvibes.py::test_compute_gv_appendix_boltzmann_import_names_bug`.
    Corregido a los nombres reales; no se reimplementa esa lógica acá
    para no duplicarla.
    """
    try:
        from .boltzmann import compute_boltzmann_multi_t, compute_g_total
    except ImportError as exc:
        raise ImportError(
            "compute_gv_appendix() necesita compute_g_total() y "
            "compute_boltzmann_multi_t() desde analysis/boltzmann.py, que "
            "todavía no parece estar migrado (o vive en otra ruta). "
            "Migrar results_collector.py -> analysis/boltzmann.py primero, "
            "o ajustar este import si el módulo quedó en otro lugar."
        ) from exc

    from ..io.parsing import parse_thermochem

    mol_dir = runs_root / mol_name

    # ── Asegurar correcciones térmicas ORCA en df_conf ────────────────────
    # Si compute_g_total ya fue llamado, G_total_298K existirá.
    # Si no (flujo normal), lo recalculamos internamente aquí.
    if "G_total_298K" not in df_conf.columns or "H_corr" not in df_conf.columns:
        # Cargar termoquímica ORCA para cada confórmero
        df_work = df_conf.copy()
        _thermochem_cache: list[Optional[dict]] = []
        for _, row in df_work.iterrows():
            conf_id  = row["conformer"]
            freq_out = mol_dir / SUBDIR_FREQ / conf_id / "freq.out"
            tc = parse_thermochem(freq_out, T=T_STD)
            _thermochem_cache.append(tc if tc.get("ok") else None)
        df_work["_thermochem"] = _thermochem_cache
        df_work = compute_g_total(df_work)
        df_work = compute_boltzmann_multi_t(df_work)
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


def _gv_version() -> str:
    """Devuelve la versión de GoodVibes instalada, o '?' si no se puede determinar.

    Migrado desde results_collector.py:2611-2625, sin cambios de lógica
    (ya usaba la `_resolve_gv_python()` canónica en el legacy también).
    """
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


def print_gv_appendix(df_gv: pd.DataFrame, mol_name: str, runs_root: Path) -> None:
    """
    Imprime el apéndice GoodVibes con el mismo estilo que el resto del script.
    Muestra tres secciones:
      1. Correcciones absolutas ORCA vs GV (G_corr, S, ZPE)
      2. Impacto en ΔG relativo y pesos Boltzmann a 298 K
      3. Cambios de ranking conformacional (si los hay)

    Migrado desde results_collector.py:2468-2609. `WORK_ROOT` reemplazado
    por `runs_root: Path` explícito (Issue C0) -- por eso esta función
    ahora recibe `runs_root`, que el original no necesitaba al depender
    del global.
    """
    sep = "─" * 74
    print(f"\n{sep}")
    print("  APÉNDICE GOODVIBES — Corrección cuasi-armónica mRRHO (Grimme)")
    print(f"  Molécula : {mol_name}")
    print(f"  Método   : -q --qs grimme -f {GV_CUTOFF:.0f} cm⁻¹ (v{_gv_version()})")
    print("  Ciclo    : G_GV = E_SP(TZVPD/CPCM) + G_corr_GV(mRRHO, SVP/CPCM)")
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
          f"[ORCA=RRHO  |  GV=mRRHO Grimme, f₀={GV_CUTOFF:.0f} cm⁻¹]")
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
        print(f"  Interpretación   : GV{'aumenta' if dg_vals.mean() > 0 else 'reduce'} G_corr "
              f"{'(mRRHO menos negativo → G más alto)' if dg_vals.mean() > 0 else '(mRRHO más negativo → G más bajo)'}")

    # ── Sección 2: Impacto en ΔG relativo y Boltzmann ────────────────────
    print("\n  ── 2. Impacto en ΔG relativo y pesos Boltzmann (T = 298.15 K) ──")
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
        print("\n  ── 3. Ranking conformacional: ORCA vs GoodVibes ──")
        print(f"  Inversiones detectadas: {n_inv}/{len(ok)}")

        if n_inv > 0:
            print(f"\n  {'conf':>10}  {'rank_ORCA':>10}  {'rank_GV':>8}  {'Δrank':>7}")
            print("  " + "─" * 42)
            for _, r in df_gv[df_gv["rank_change_GV"].abs() > 0].iterrows():
                rch = r.get("rank_change_GV", 0)
                flag = " ◀ inversión"
                print(
                    f"  {r['conformer']:>10}  "
                    f"{str(r.get('rank_ORCA', '—')):>10}  "
                    f"{str(r.get('rank_GV', '—')):>8}  "
                    f"{int(rch):>+7}{flag}"
                )
        else:
            print("  → El ranking conformacional es idéntico con ORCA-RRHO y GV-mRRHO.")

    print("\n  Archivos GV guardados en:")
    for conf_id in df_gv["conformer"]:
        p = runs_root / mol_name / SUBDIR_FREQ / conf_id / GV_SUBDIR
        if p.exists():
            # El nombre del JSON depende del stem del .out (normalmente freq_goodvibes.json)
            json_files = list(p.glob("*_goodvibes.json"))
            fname = json_files[0].name if json_files else "*.json"
            print(f"    {p}/{fname}")

    print(f"\n{sep}\n")


def export_gv_appendix(df_gv: pd.DataFrame, mol_name: str, runs_root: Path, fmt: str) -> None:
    """
    Guarda el DataFrame del apéndice GoodVibes en <runs_root>/<mol_name>/.
    Formatos: "csv", "xlsx" o "both".
    Si fmt incluye xlsx, lo agrega como hoja extra del Excel existente si ya
    existe, o crea un archivo dedicado goodvibes_appendix.xlsx.

    Migrado desde results_collector.py:2628-2679. `WORK_ROOT` reemplazado
    por `runs_root: Path` explícito (Issue C0).
    """
    if df_gv.empty:
        return

    mol_dir   = runs_root / mol_name
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
