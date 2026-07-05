"""Visualización: espectros UV-Vis, comparaciones, distribución de Boltzmann.

Migrar desde results_collector.py:
  - _build_spectrum()                 líneas 1539-1559
  - plot_spectrum()                   líneas 1560-1602
  - plot_comparison()                 líneas 1603-1718
  - plot_adc2_tddft_comparison()      líneas 3130-3259

y desde protonation_cycle.py:
  - plot_distribution()               líneas 1919-1984
  - plot_distribution_5()             líneas 1985-2049
  - plot_conformer_boltzmann()        líneas 2050-2108

Nota: el legacy importa matplotlib de forma perezosa (try/except dentro
de cada función) para no romper corridas headless sin backend gráfico.
matplotlib es dependencia dura en pyproject.toml de maaqc -- si se
quiere preservar ese fallback opcional, mover matplotlib a un extra
`[project.optional-dependencies] plots` y mantener el try/except acá.
"""

from __future__ import annotations

import pandas as pd


def plot_spectrum(df_tddft: pd.DataFrame, df_conf: pd.DataFrame, mol_name: str, **kwargs):
    raise NotImplementedError("Migrar desde results_collector.py:plot_spectrum (1560-1602)")
