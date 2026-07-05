"""Ponderación Boltzmann multi-temperatura y G_total compuesto.

Migrar desde results_collector.py:
  - compute_G_total()               líneas 648-764
  - compute_boltzmann_multiT()      líneas 765-848
  - print_boltzmann_analysis()      líneas 849-932

y desde protonation_cycle.py (revisar solapamiento antes de migrar,
puede haber la misma duplicación ya detectada entre parsers):
  - boltzmann_protonation_states()  líneas 1586-1637
  - boltzmann_5state()              líneas 1638-1703

y desde pipeline.py (selección de confórmeros basada en pesos Boltzmann,
no cálculo de pesos en sí — posible dependencia cruzada a resolver):
  - parse_boltzmann_weights()       líneas 2108-2137
  - select_nconf_boltzmann()        líneas 2138-2276
"""

from __future__ import annotations

import pandas as pd


def compute_g_total(df_conf: pd.DataFrame) -> pd.DataFrame:
    raise NotImplementedError("Migrar desde results_collector.py:compute_G_total (648-764)")


def compute_boltzmann_multi_t(df_conf: pd.DataFrame) -> pd.DataFrame:
    raise NotImplementedError(
        "Migrar desde results_collector.py:compute_boltzmann_multiT (765-848)"
    )
