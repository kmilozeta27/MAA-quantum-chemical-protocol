"""Comparación ADC2 vs TDDFT (el benchmark central del paper).

Migrar desde results_collector.py:
  - collect_adc2_results()          líneas 2822-2849
  - compute_adc2_lambda_max()       líneas 2850-2885
  - compare_adc2_tddft()            líneas 2886-2968
  - _boltz_weighted_lmax()          líneas 2969-2985
  - print_adc2_summary()            líneas 2986-3019
  - print_adc2_tddft_comparison()   líneas 3020-3116
  - export_adc2_comparison()        líneas 3117-3129

Depende de io/parsing.py:parse_adc2_excitations y de boltzmann.py
(pesos Boltzmann para el λmax ponderado) — no reimplementar ninguno.
"""

from __future__ import annotations

import pandas as pd


def compare_adc2_tddft(*args, **kwargs) -> pd.DataFrame:
    raise NotImplementedError("Migrar desde results_collector.py:compare_adc2_tddft (2886-2968)")
