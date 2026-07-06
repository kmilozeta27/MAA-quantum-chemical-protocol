"""Wrapper de GoodVibes (correcciones termoquímicas alternativas a mano).

No estaba contemplado como módulo propio en el esqueleto original de la
skill maa-repackaging (que solo listaba boltzmann/pka_cycle/monte_carlo/
parse_uvvis bajo analysis/) — se agrega acá porque el legacy tiene
suficiente lógica propia de GoodVibes como para no mezclarla con
results_collector.py ni con pka_cycle.py.

FUENTE CANÓNICA (Issue C1, ya migrado): `_resolve_gv_python()` y
`run_goodvibes()` vienen de `pka_cycle.py:282-397` — es la versión que
ya estaba en uso activo para el ciclo de pKa (llamada indirectamente vía
`get_thermochem(..., source="goodvibes")`). `pka_cycle.py` ahora importa
estas dos funciones desde acá en vez de definirlas localmente.

Pendiente (segunda mitad de Issue C1, todavía sin resolver): revisar si
lo siguiente, migrado desde legacy, aporta algo que la versión de
pka_cycle.py no tenga, o si es redundante y se descarta. NO se migra
nada de esto todavía -- son notas de dónde está en el legacy, no código.

Desde results_collector.py:
  - _resolve_gv_python()        líneas 1719-1772  (variante propia, corre
                                 aparte de la de pka_cycle.py -- comparar
                                 antes de decidir si se unifican)
  - _probe_gv_capabilities()    líneas 1773-1811
  - _parse_goodvibes_dat()      líneas 1812-1877
  - run_goodvibes_on_mol()      líneas 1981-2181
  - _parse_goodvibes_json()     líneas 2182-2235
  - compute_gv_appendix()       líneas 2284-2467
  - print_gv_appendix()         líneas 2468-2610
  - _gv_version()               líneas 2611-2627
  - export_gv_appendix()        líneas 2628-2709

y desde protonation_cycle.py (revisar solapamiento con lo de arriba
antes de migrar -- ambos scripts invocan GoodVibes por separado):
  - _resolve_gv_python()        líneas 525-587  (otra copia/variante más,
                                 ya van 3 candidatas: pka_cycle.py [ahora
                                 canónica acá], results_collector.py y
                                 esta -- ninguna de las dos últimas se
                                 migró todavía)
  - run_goodvibes()             líneas 588-708
  - _print_goodvibes_comparison() líneas 2721-2795
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Optional

# ─────────────────────────────────────────────
#  CONSTANTES  (idénticas a las que tenía pka_cycle.py antes de migrar)
# ─────────────────────────────────────────────

T_STD = 298.15             # K
GV_CUTOFF = 100.0
GV_QH = False
GV_QS = "grimme"
GV_CONDA_ENV = "sciencebitch"
_GV_PYTHON_CACHE: Optional[str] = None


def _resolve_gv_python(force: bool = False) -> str:
    """Ruta al intérprete Python con goodvibes instalado. Cacheada.

    Migrado desde pka_cycle.py:282-323 (fuente canónica elegida para esta
    migración porque era la versión efectivamente en uso por el ciclo de
    pKa). Sin cambios de lógica respecto al original.
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
    -- pka_cycle.py ahora importa esto desde acá en vez de tenerlo local.
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


def run_goodvibes_on_mol(mol_name: str) -> dict[str, dict]:
    raise NotImplementedError(
        "Migrar desde results_collector.py:run_goodvibes_on_mol (1981-2181); "
        "revisar solapamiento con protonation_cycle.py:run_goodvibes (588-708) "
        "antes de decidir si son la misma función o dos usos distintos. "
        "OJO: esta es una función DISTINTA de run_goodvibes() (arriba), que "
        "ya está migrada -- esta variante '_on_mol' opera a nivel de "
        "molécula/conjunto de confórmeros, no de un único freq.out."
    )
