"""Tests de contrato para maaqc.analysis.adc2_compare.

A diferencia de test_imports.py (que solo prueba "importa sin explotar"),
estos tests ejercitan las funciones puras de este módulo con datos —
en particular con los dos números reales extraídos de los outputs ORCA
subidos en la sesión de la Tarea A:

  - adc2.out  → bloque "SPECTRUM FOR LEFT-RIGHT TRANSITION MOMENTS"
                (la prioridad correcta post-fix del bug LEFT/RIGHT):
                transición 0-1A -> 1-1A, wl_nm = 357.8, fosc = 0.416153266
  - tddft.out → estado más brillante (0-1A -> 1-1A) de los 15 calculados:
                wavelength_nm = 299.2, fosc = 0.772302746

Alcance: solo se testean compute_adc2_lambda_max() y
_boltz_weighted_lmax(), que son funciones puras sobre DataFrames (no
tocan filesystem ni llaman a collect_results()/parse_adc2_excitations()).
compare_adc2_tddft() en sí depende de results_collector.collect_results()
y de un layout de directorio real (07_adc2/conf_*/adc2.out, 04_sp/, etc.)
que no está migrado/disponible en esta sesión -- se deja fuera de
alcance a propósito en vez de mockear un contrato que no se ha visto.
"""

from __future__ import annotations

import importlib

import pandas as pd
import pytest

adc2_compare = importlib.import_module("maaqc.analysis.adc2_compare")

# Números reales, ver docstring del módulo.
REAL_ADC2_WL_NM = 357.8
REAL_ADC2_FOSC_LEFT_RIGHT = 0.416153266
REAL_TDDFT_WL_NM = 299.2
REAL_TDDFT_FOSC = 0.772302746


def test_compute_adc2_lambda_max_single_conformer_matches_real_output():
    """Con un solo confórmero (peso Boltzmann = 1.0), el λ_max ADC2
    Boltzmann-promediado debe ser exactamente el wl_nm reportado por
    ORCA en el bloque LEFT-RIGHT de adc2.out -- no hay promedio que
    hacer con un solo dato."""
    df_adc2 = pd.DataFrame([
        {"conformer": "conf_001", "state": 1, "wl_nm": REAL_ADC2_WL_NM,
         "fosc": REAL_ADC2_FOSC_LEFT_RIGHT},
    ])
    df_conf = pd.DataFrame([
        {"conformer": "conf_001", "boltz_298K": 1.0},
    ])

    lmax = adc2_compare.compute_adc2_lambda_max(df_adc2, df_conf, "298K")

    assert lmax == pytest.approx(REAL_ADC2_WL_NM)


def test_compute_adc2_lambda_max_weighted_average_two_conformers():
    """Con dos confórmeros, el resultado debe ser el promedio ponderado
    por boltz_298K, no un promedio simple -- esto es lo que distingue
    'Boltzmann-promediado' de solo tomar el S1 más bajo en energía."""
    df_adc2 = pd.DataFrame([
        {"conformer": "conf_001", "state": 1, "wl_nm": 300.0, "fosc": 0.5},
        {"conformer": "conf_002", "state": 1, "wl_nm": 320.0, "fosc": 0.5},
    ])
    df_conf = pd.DataFrame([
        {"conformer": "conf_001", "boltz_298K": 0.7},
        {"conformer": "conf_002", "boltz_298K": 0.3},
    ])

    lmax = adc2_compare.compute_adc2_lambda_max(df_adc2, df_conf, "298K")

    expected = (300.0 * 0.7 + 320.0 * 0.3) / (0.7 + 0.3)
    assert lmax == pytest.approx(expected)


def test_compute_adc2_lambda_max_ignores_non_s1_states():
    """Solo state == 1 (S1) debe contar -- un S2 con peso alto no debe
    contaminar el λ_max reportado."""
    df_adc2 = pd.DataFrame([
        {"conformer": "conf_001", "state": 1, "wl_nm": REAL_ADC2_WL_NM, "fosc": 0.4},
        {"conformer": "conf_001", "state": 2, "wl_nm": 180.0, "fosc": 0.9},
    ])
    df_conf = pd.DataFrame([{"conformer": "conf_001", "boltz_298K": 1.0}])

    lmax = adc2_compare.compute_adc2_lambda_max(df_adc2, df_conf, "298K")

    assert lmax == pytest.approx(REAL_ADC2_WL_NM)


@pytest.mark.parametrize("df_adc2,df_conf", [
    (pd.DataFrame(), pd.DataFrame([{"conformer": "conf_001", "boltz_298K": 1.0}])),
    (pd.DataFrame([{"conformer": "conf_001", "state": 1, "wl_nm": 300.0}]), pd.DataFrame()),
])
def test_compute_adc2_lambda_max_returns_none_on_empty_input(df_adc2, df_conf):
    assert adc2_compare.compute_adc2_lambda_max(df_adc2, df_conf, "298K") is None


def test_compute_adc2_lambda_max_returns_none_when_boltz_column_missing():
    """Si df_conf no tiene la columna boltz_{temp_label} (p.ej. se pidió
    una temperatura para la que no se corrió el ciclo Boltzmann), la
    función no debe reventar -- debe devolver None."""
    df_adc2 = pd.DataFrame([{"conformer": "conf_001", "state": 1, "wl_nm": 300.0}])
    df_conf = pd.DataFrame([{"conformer": "conf_001", "boltz_274K": 1.0}])

    assert adc2_compare.compute_adc2_lambda_max(df_adc2, df_conf, "298K") is None


def test_boltz_weighted_lmax_single_conformer_matches_real_tddft_output():
    """Mismo razonamiento que el test ADC2 de arriba, pero para
    _boltz_weighted_lmax() (TDDFT). Nótese que la columna acá es
    'wavelength_nm', no 'wl_nm' como en el DataFrame de ADC2 -- son
    nombres de columna distintos entre los dos módulos, a propósito."""
    df_tddft = pd.DataFrame([
        {"conformer": "conf_001", "state": 1,
         "wavelength_nm": REAL_TDDFT_WL_NM, "fosc": REAL_TDDFT_FOSC},
    ])
    df_conf = pd.DataFrame([{"conformer": "conf_001", "boltz_298K": 1.0}])

    lmax = adc2_compare._boltz_weighted_lmax(df_tddft, df_conf, "298K")

    assert lmax == pytest.approx(REAL_TDDFT_WL_NM)


def test_boltz_weighted_lmax_returns_none_when_all_weights_zero():
    df_tddft = pd.DataFrame([
        {"conformer": "conf_001", "state": 1, "wavelength_nm": 300.0, "fosc": 0.5},
    ])
    df_conf = pd.DataFrame([{"conformer": "conf_001", "boltz_298K": 0.0}])

    assert adc2_compare._boltz_weighted_lmax(df_tddft, df_conf, "298K") is None


def test_real_adc2_vs_tddft_delta_matches_manual_calculation():
    """Sanity check end-to-end (sin pasar por compare_adc2_tddft(), que
    necesita filesystem real): con los dos números reales de esta sesión,
    el delta_nm que compare_adc2_tddft() calcularía (λ_TDDFT − λ_ADC2)
    es un blueshift de TDDFT de ~58.6 nm respecto a ADC2 -- vale la pena
    dejarlo anclado en un test para no perder de vista esta cifra si en
    algún momento compare_adc2_tddft() se refactoriza."""
    delta_nm = REAL_TDDFT_WL_NM - REAL_ADC2_WL_NM
    assert delta_nm == pytest.approx(-58.6, abs=0.05)
