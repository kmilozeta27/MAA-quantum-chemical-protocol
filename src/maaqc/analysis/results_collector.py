"""Consolidación de resultados por estado de protonación (`maaqc collect`).

Migrar desde results_collector.py (legacy), la parte que NO es parsing
puro (eso va a io/parsing.py) ni Boltzmann/pKa/GoodVibes/ADC2/plots (eso
va a los otros módulos de analysis/):

  - analyze_crest_to_dft()              líneas 524-587
  - print_crest_dft_analysis()          líneas 588-647
  - collect_results()                   líneas 1229-1319
  - print_summary()                     líneas 1320-1385
  - compare_protonation_states()        líneas 1386-1486
  - print_comparison()                  líneas 1487-1538
  - export_results() / export_comparison() / export_adc2_mol()
                                         líneas 3260-3341

Este módulo depende de io/parsing.py (parsers) y de boltzmann.py
(compute_G_total, compute_boltzmann_multiT) — NO debe reimplementar
ninguno de los dos.
"""

from __future__ import annotations

import pandas as pd


def collect_results(mol_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    raise NotImplementedError("Migrar desde results_collector.py:collect_results (1229-1319)")


def compare_protonation_states(*args, **kwargs):
    raise NotImplementedError(
        "Migrar desde results_collector.py:compare_protonation_states (1386-1486)"
    )
