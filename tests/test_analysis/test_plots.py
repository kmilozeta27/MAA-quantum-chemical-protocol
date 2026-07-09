"""Tests de contrato para maaqc.analysis.plots.

matplotlib es dependencia OPCIONAL (Issue C5) -- todos los tests que
generan figuras usan pytest.importorskip("matplotlib") para no romper
en runs sin el extra `maaqc[plots]` instalado (headless CI, `maaqc doctor`
sin X11, etc.), consistente con el try/except perezoso que el propio
módulo implementa función por función.

_build_spectrum() es la excepción: es pura numpy/pandas (no importa
matplotlib), así que se testea sin el importorskip.

TEMP_LABELS se lee de maaqc.constants en vez de hardcodearse acá, porque
plot_spectrum() itera exactamente esa lista -- si algún día cambia
(se agrega/quita una temperatura), este test debe seguir siendo válido
sin tocarlo.
"""

from __future__ import annotations

import importlib

import numpy as np
import pandas as pd
import pytest

plots = importlib.import_module("maaqc.analysis.plots")
constants = importlib.import_module("maaqc.constants")

REAL_ADC2_WL_NM = 357.8
REAL_TDDFT_WL_NM = 299.2


# ─────────────────────────────────────────────
#  Fixtures compartidas
# ─────────────────────────────────────────────

def _make_df_conf(conformers_boltz298: dict) -> pd.DataFrame:
    """df_conf mínimo con boltz_298K y boltz_{label} = boltz_298K para
    cada label de TEMP_LABELS (suficiente para que _build_spectrum /
    plot_spectrum no exploten buscando columnas inexistentes)."""
    rows = []
    for conf, w in conformers_boltz298.items():
        row = {"conformer": conf, "boltz_298K": w}
        for lbl in constants.TEMP_LABELS:
            row.setdefault(f"boltz_{lbl}", w)
        rows.append(row)
    return pd.DataFrame(rows)


def _make_df_tddft(wavelength_nm: float, fosc: float, conformer: str = "conf_001") -> pd.DataFrame:
    return pd.DataFrame([
        {"conformer": conformer, "state": 1, "wavelength_nm": wavelength_nm, "fosc": fosc},
    ])


# ─────────────────────────────────────────────
#  _build_spectrum (pura, sin matplotlib)
# ─────────────────────────────────────────────

def test_build_spectrum_peaks_near_real_tddft_wavelength():
    """El espectro gaussiano-suavizado debe tener su máximo cerca del
    wavelength_nm real de tddft.out (299.2 nm), con sigma chico para que
    el pico no se difumine fuera de tolerancia."""
    df_tddft = _make_df_tddft(REAL_TDDFT_WL_NM, fosc=0.77)
    df_conf = _make_df_conf({"conf_001": 1.0})
    x = np.linspace(200, 500, 1500)

    y = plots._build_spectrum(df_tddft, df_conf, sigma=5.0, x=x, temp_label="298K")

    assert y.shape == x.shape
    peak_x = x[np.argmax(y)]
    assert peak_x == pytest.approx(REAL_TDDFT_WL_NM, abs=1.0)


def test_build_spectrum_empty_tddft_returns_zeros():
    x = np.linspace(200, 500, 100)
    y = plots._build_spectrum(pd.DataFrame(), pd.DataFrame(), sigma=20.0, x=x, temp_label="298K")
    assert y.shape == x.shape
    assert np.all(y == 0)


def test_build_spectrum_missing_boltz_column_falls_back_to_298K():
    """Si se pide un temp_label sin columna boltz_{label} en df_conf, cae
    a boltz_298K en vez de reventar (ver `col = 'boltz_298K'` de fallback
    en la función)."""
    df_tddft = _make_df_tddft(REAL_TDDFT_WL_NM, fosc=0.77)
    df_conf = pd.DataFrame([{"conformer": "conf_001", "boltz_298K": 1.0}])
    x = np.linspace(200, 500, 500)

    y = plots._build_spectrum(df_tddft, df_conf, sigma=10.0, x=x, temp_label="999K_no_existe")

    assert not np.all(y == 0)


# ─────────────────────────────────────────────
#  plot_spectrum / plot_comparison / plot_adc2_tddft_comparison
# ─────────────────────────────────────────────

def test_plot_spectrum_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    df_tddft = _make_df_tddft(REAL_TDDFT_WL_NM, fosc=0.77)
    df_conf = _make_df_conf({"conf_001": 1.0})
    out = tmp_path / "uvvis_spectrum.png"

    plots.plot_spectrum(df_tddft, df_conf, mol_name="R_test", output_path=out)

    assert out.exists() and out.stat().st_size > 0


def test_plot_spectrum_empty_tddft_does_not_write_file(tmp_path, capsys):
    pytest.importorskip("matplotlib")
    out = tmp_path / "should_not_exist.png"

    plots.plot_spectrum(pd.DataFrame(), pd.DataFrame(), mol_name="R_test", output_path=out)

    assert not out.exists()
    assert "Sin datos TDDFT" in capsys.readouterr().out


def test_plot_comparison_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    dfs_conf = {
        "R": _make_df_conf({"conf_001": 1.0}),
        "B": _make_df_conf({"conf_001": 1.0}),
    }
    dfs_tddft = {
        "R": _make_df_tddft(REAL_TDDFT_WL_NM, 0.77),
        "B": _make_df_tddft(304.0, 0.6),
    }
    out = tmp_path / "comparison_uvvis.png"

    plots.plot_comparison(dfs_conf, dfs_tddft, output_path=out)

    assert out.exists() and out.stat().st_size > 0


def test_plot_adc2_tddft_comparison_writes_png(tmp_path):
    """df_cmp con la forma exacta que retorna
    adc2_compare.compare_adc2_tddft() (columnas molecule, lambda_tddft_nm,
    lambda_adc2_nm, delta_nm, delta_eV, n_confs_tddft, n_confs_adc2),
    usando los dos números reales de esta sesión."""
    pytest.importorskip("matplotlib")
    df_cmp = pd.DataFrame([{
        "molecule": "R",
        "lambda_tddft_nm": REAL_TDDFT_WL_NM,
        "lambda_adc2_nm": REAL_ADC2_WL_NM,
        "delta_nm": REAL_TDDFT_WL_NM - REAL_ADC2_WL_NM,
        "delta_eV": 0.0,
        "n_confs_tddft": 1,
        "n_confs_adc2": 1,
    }])
    dfs_conf = {"R": _make_df_conf({"conf_001": 1.0})}
    dfs_tddft = {"R": _make_df_tddft(REAL_TDDFT_WL_NM, 0.77)}
    out = tmp_path / "comparison_adc2_tddft.png"

    plots.plot_adc2_tddft_comparison(df_cmp, dfs_conf, dfs_tddft, output_path=out)

    assert out.exists() and out.stat().st_size > 0


def test_plot_adc2_tddft_comparison_empty_df_cmp_does_not_write(tmp_path):
    pytest.importorskip("matplotlib")
    out = tmp_path / "should_not_exist.png"

    plots.plot_adc2_tddft_comparison(pd.DataFrame(), {}, {}, output_path=out)

    assert not out.exists()


# ─────────────────────────────────────────────
#  plot_distribution / plot_distribution_5
# ─────────────────────────────────────────────

def test_plot_distribution_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    ph = np.linspace(0, 14, 50)
    df_dist = pd.DataFrame({
        "pH": ph,
        "f_neutral": 1 / (1 + 10 ** (ph - 1.34)),
        "f_anion": 1 - 1 / (1 + 10 ** (ph - 1.34)),
    })
    out = tmp_path / "distribution.png"

    plots.plot_distribution(df_dist, pKa1=1.34, pKa2=24.75, output_path=out, ph_highlight=8.1)

    assert out.exists() and out.stat().st_size > 0


def test_plot_distribution_5_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    ph = np.linspace(0, 14, 50)
    df_dist = pd.DataFrame({
        "pH": ph,
        "f_R": np.full_like(ph, 0.4),
        "f_D": np.full_like(ph, 0.1),
        "f_anion": np.full_like(ph, 0.5),
    })
    out = tmp_path / "distribution_5.png"

    plots.plot_distribution_5(df_dist, pKa0=None, pKa1=1.34, pKa2=24.75,
                               output_path=out, ph_highlight=8.1)

    assert out.exists() and out.stat().st_size > 0


# ─────────────────────────────────────────────
#  plot_conformer_boltzmann
# ─────────────────────────────────────────────

def test_plot_conformer_boltzmann_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    data_list = [
        {"mol_name": "R", "df": pd.DataFrame([
            {"conformer": "conf_001", "boltz_298K": 0.85},
            {"conformer": "conf_002", "boltz_298K": 0.15},
        ])},
        {"mol_name": "B", "df": pd.DataFrame([
            {"conformer": "conf_001", "boltz_298K": 1.0},
        ])},
    ]
    out = tmp_path / "conformer_boltzmann.png"

    plots.plot_conformer_boltzmann(data_list, output_path=out)

    assert out.exists() and out.stat().st_size > 0


def test_plot_conformer_boltzmann_handles_empty_df_gracefully(tmp_path):
    """Un estado sin datos (df vacío o sin boltz_298K) no debe reventar
    el gráfico completo -- debe mostrar 'Sin datos' en ese panel y seguir."""
    pytest.importorskip("matplotlib")
    data_list = [
        {"mol_name": "R", "df": pd.DataFrame()},
        {"mol_name": "B", "df": pd.DataFrame([{"conformer": "conf_001", "boltz_298K": 1.0}])},
    ]
    out = tmp_path / "conformer_boltzmann_partial.png"

    plots.plot_conformer_boltzmann(data_list, output_path=out)

    assert out.exists() and out.stat().st_size > 0
