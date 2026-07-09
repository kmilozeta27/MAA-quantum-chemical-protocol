"""Tests de contrato para maaqc.analysis.goodvibes.

Alcance: _parse_goodvibes_dat() y _parse_goodvibes_json() son funciones
puras de parseo de texto/JSON (no llaman a subprocess ni a GoodVibes de
verdad), así que son testeables sin mockear nada externo -- se
construyen archivos .dat/.json sintéticos pero con el formato real que
GoodVibes v3.x/v4.x produce (documentado en las docstrings de cada
función).

run_goodvibes(), run_goodvibes_on_mol() y _resolve_gv_python() SÍ llaman
a subprocess/goodvibes real y no se testean acá (necesitarían el binario
goodvibes instalado o mockear subprocess.run, fuera de alcance de esta
tarea de docstring+contrato mínimo).

Incluye además un test de regresión (test_compute_gv_appendix_*) para el
bug ya identificado en la migración: compute_gv_appendix() (línea ~856)
importa `compute_G_total, compute_boltzmann_multiT` (camelCase) desde
.boltzmann, pero boltzmann.py migró a snake_case
(compute_g_total, compute_boltzmann_multi_t). Este test debe FALLAR
hasta que se corrija el import en goodvibes.py -- es intencional, mismo
patrón que test_goodvibes_symbols_used_by_pka_cycle_exist() en
test_imports.py.
"""

from __future__ import annotations

import importlib
import json
import tempfile
from pathlib import Path

import pytest

goodvibes = importlib.import_module("maaqc.analysis.goodvibes")


# ─────────────────────────────────────────────
#  _parse_goodvibes_dat
# ─────────────────────────────────────────────

# Valores compartidos entre los tests de offset=0 y offset=1, para que
# ambos caminos de parseo (con/sin símbolo de estado al inicio) deban
# producir exactamente el mismo resultado numérico.
_E_ELEC = -1220.967288
_ZPE = 0.403359
_H = -1220.550000
_TS = 0.030000
_TQS = 0.028000
_G = -1220.580000
_QHG = -1220.578000
_T = 298.15


def test_parse_goodvibes_dat_v4_with_state_symbol_offset(tmp_path):
    """GoodVibes v4 antepone un símbolo de estado ('o', 'x', '?') antes
    del nombre de archivo en la línea de datos, lo que desplaza todos
    los índices +1 respecto al encabezado -- este es el caso
    explícitamente documentado en el comentario de la función."""
    dat_path = tmp_path / "Goodvibes_output.dat"
    dat_path.write_text(
        "GoodVibes v4.1 -- Temp = 298.15 K\n"
        "   Structure       E          ZPE          H         qh-H"
        "        T.S       T.qh-S        G(T)       qh-G(T)\n"
        f"   o   freq   {_E_ELEC}   {_ZPE}   {_H}   {_H - 0.001}"
        f"   {_TS}   {_TQS}   {_G}   {_QHG}\n"
    )

    result = goodvibes._parse_goodvibes_dat(dat_path)

    assert result["ok"] is True
    assert result["T_GV"] == pytest.approx(_T)
    assert result["E_elec_GV"] == pytest.approx(_E_ELEC)
    assert result["ZPE_GV"] == pytest.approx(_ZPE)
    assert result["G_GV"] == pytest.approx(_QHG)  # prefiere qh-G(T) sobre G(T)
    assert result["TS_GV"] == pytest.approx(_TQS)  # prefiere T.qh-S sobre T.S
    assert result["S_GV"] == pytest.approx(_TQS / _T)
    assert result["G_corr_GV"] == pytest.approx(_QHG - _E_ELEC)


def test_parse_goodvibes_dat_v3_without_state_symbol(tmp_path):
    """GoodVibes v3.x clásico: sin símbolo de estado, datos desde
    parts[1] (offset=0). Con los mismos valores numéricos que el test
    v4 de arriba, para verificar que ambos caminos de parseo coinciden."""
    dat_path = tmp_path / "Goodvibes_output.dat"
    dat_path.write_text(
        "GoodVibes v3.2 -- at 298.15 K\n"
        "   Structure       E          ZPE          H"
        "        T.S       T.qh-S        G(T)       qh-G(T)\n"
        f"   freq   {_E_ELEC}   {_ZPE}   {_H}"
        f"   {_TS}   {_TQS}   {_G}   {_QHG}\n"
    )

    result = goodvibes._parse_goodvibes_dat(dat_path)

    assert result["ok"] is True
    assert result["G_GV"] == pytest.approx(_QHG)
    assert result["G_corr_GV"] == pytest.approx(_QHG - _E_ELEC)


def test_parse_goodvibes_dat_defaults_temperature_when_not_found(tmp_path):
    """Si no se detecta 'Temp = X K' ni 'at X K' en el archivo, debe
    asumir 298.15 K en vez de fallar."""
    dat_path = tmp_path / "Goodvibes_output.dat"
    dat_path.write_text(
        "GoodVibes v3.2 -- sin línea de temperatura reconocible\n"
        "   Structure       E          ZPE          H"
        "        T.S       T.qh-S        G(T)       qh-G(T)\n"
        f"   freq   {_E_ELEC}   {_ZPE}   {_H}"
        f"   {_TS}   {_TQS}   {_G}   {_QHG}\n"
    )

    result = goodvibes._parse_goodvibes_dat(dat_path)

    assert result["T_GV"] == pytest.approx(298.15)


def test_parse_goodvibes_dat_missing_file_returns_warning(tmp_path):
    result = goodvibes._parse_goodvibes_dat(tmp_path / "no_existe.dat")

    assert result["ok"] is False
    assert "no encontrado" in result["warning"]


def test_parse_goodvibes_dat_missing_header_returns_warning(tmp_path):
    dat_path = tmp_path / "sin_header.dat"
    dat_path.write_text("esto no tiene ninguna tabla de GoodVibes adentro\n")

    result = goodvibes._parse_goodvibes_dat(dat_path)

    assert result["ok"] is False
    assert "Structure" in result["warning"]


# ─────────────────────────────────────────────
#  _parse_goodvibes_json
# ─────────────────────────────────────────────

def test_parse_goodvibes_json_v4_single_molecule(tmp_path):
    json_path = tmp_path / "freq_goodvibes.json"
    json_path.write_text(json.dumps({
        "version": "4.1",
        "molecules": [{
            "filename": "freq.out", "t": _T, "zpe": _ZPE,
            "h": _H, "ts": _TQS, "g": _QHG, "sp": _E_ELEC,
        }],
    }))

    result = goodvibes._parse_goodvibes_json(json_path, "freq.out")

    assert result["ok"] is True
    assert result["T_GV"] == pytest.approx(_T)
    assert result["G_GV"] == pytest.approx(_QHG)
    assert result["E_elec_GV"] == pytest.approx(_E_ELEC)
    assert result["S_GV"] == pytest.approx(_TQS / _T)
    assert result["G_corr_GV"] == pytest.approx(_QHG - _E_ELEC)


def test_parse_goodvibes_json_matches_correct_entry_by_filename(tmp_path):
    """Con varias moléculas en el mismo JSON (GoodVibes corrido sobre un
    batch), debe elegir la entrada cuyo filename coincide con el pedido,
    no simplemente la primera."""
    json_path = tmp_path / "batch_goodvibes.json"
    json_path.write_text(json.dumps({
        "molecules": [
            {"filename": "otro_conformero.out", "t": _T, "g": -999.0, "sp": -998.0},
            {"filename": "freq.out", "t": _T, "g": _QHG, "sp": _E_ELEC},
        ],
    }))

    result = goodvibes._parse_goodvibes_json(json_path, "freq.out")

    assert result["G_GV"] == pytest.approx(_QHG)
    assert result["E_elec_GV"] == pytest.approx(_E_ELEC)


def test_parse_goodvibes_json_falls_back_to_single_entry(tmp_path):
    """Si hay una sola entrada y el nombre no calza exactamente (p.ej.
    GoodVibes normalizó el filename de otra forma), se usa esa única
    entrada como fallback en vez de fallar."""
    json_path = tmp_path / "freq_goodvibes.json"
    json_path.write_text(json.dumps({
        "molecules": [{"filename": "algo_distinto.out", "t": _T, "g": _QHG, "sp": _E_ELEC}],
    }))

    result = goodvibes._parse_goodvibes_json(json_path, "freq.out")

    assert result["G_GV"] == pytest.approx(_QHG)


def test_parse_goodvibes_json_no_match_multiple_entries_raises_valueerror(tmp_path):
    """Con más de una entrada y ninguna coincidencia por nombre, no hay
    fallback razonable -- debe levantar ValueError explícito en vez de
    adivinar cuál usar."""
    json_path = tmp_path / "batch_goodvibes.json"
    json_path.write_text(json.dumps({
        "molecules": [
            {"filename": "conf_a.out", "t": _T, "g": -1.0, "sp": -2.0},
            {"filename": "conf_b.out", "t": _T, "g": -3.0, "sp": -4.0},
        ],
    }))

    with pytest.raises(ValueError, match="No se encontró entrada"):
        goodvibes._parse_goodvibes_json(json_path, "freq.out")


def test_parse_goodvibes_json_prefers_ts_over_direct_s():
    """Si TS y T son ambos parseables, S_GV se deriva de TS/T aunque el
    JSON también traiga una clave 's' directa -- documenta la prioridad
    real del código (rama `if TS_GV is not None and T_GV is not None`)."""
    entry = {"t": 100.0, "ts": 10.0, "s": 999.0, "g": -5.0, "sp": -6.0}
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "x_goodvibes.json"
        p.write_text(json.dumps({"molecules": [{**entry, "filename": "x.out"}]}))
        result = goodvibes._parse_goodvibes_json(p, "x.out")

    assert result["S_GV"] == pytest.approx(10.0 / 100.0)  # TS/T, no el 999.0 de "s"


# ─────────────────────────────────────────────
#  Regresión: bug conocido de nombres camelCase/snake_case
#  (ver "Top of mind" -- goodvibes.py::compute_gv_appendix())
# ─────────────────────────────────────────────

def test_compute_gv_appendix_boltzmann_import_names_bug():
    """compute_gv_appendix() (goodvibes.py, línea ~856) hace:

        from .boltzmann import compute_G_total, compute_boltzmann_multiT

    pero boltzmann.py fue migrado a snake_case
    (compute_g_total, compute_boltzmann_multi_t, ver Completed migrations
    / C2). Este test debe FALLAR hasta que se corrija esa línea de
    import en goodvibes.py -- es la misma clase de bug que D0 (nombres
    viejos sobreviviendo a una migración de convención), solo que en
    goodvibes.py en vez de boltzmann.py, y en sentido inverso (goodvibes
    todavía usa el nombre viejo).

    Fix sugerido: cambiar el import a
        from .boltzmann import compute_g_total, compute_boltzmann_multi_t
    y renombrar los dos usos locales (líneas ~883-884) a snake_case.
    """
    boltzmann = importlib.import_module("maaqc.analysis.boltzmann")

    assert hasattr(boltzmann, "compute_G_total"), (
        "maaqc.analysis.boltzmann.compute_G_total no existe -- fue "
        "renombrado a compute_g_total en la migración C2. "
        "goodvibes.py:~856 (compute_gv_appendix) todavía importa con el "
        "nombre camelCase viejo; actualizar ese import a snake_case."
    )
    assert hasattr(boltzmann, "compute_boltzmann_multiT"), (
        "maaqc.analysis.boltzmann.compute_boltzmann_multiT no existe -- "
        "fue renombrado a compute_boltzmann_multi_t en la migración C2. "
        "goodvibes.py:~856 (compute_gv_appendix) todavía importa con el "
        "nombre camelCase viejo; actualizar ese import a snake_case."
    )
