"""Wrapper de GoodVibes (correcciones termoquímicas alternativas a mano).

No estaba contemplado como módulo propio en el esqueleto original de la
skill maa-repackaging (que solo listaba boltzmann/pka_cycle/monte_carlo/
parse_uvvis bajo analysis/) — se agrega acá porque el legacy tiene
suficiente lógica propia de GoodVibes como para no mezclarla con
results_collector.py ni con pka_cycle.py.

Migrar desde results_collector.py:
  - _resolve_gv_python()        líneas 1719-1772
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
  - _resolve_gv_python()        líneas 525-587  (otra copia/variante)
  - run_goodvibes()             líneas 588-708
  - _print_goodvibes_comparison() líneas 2721-2795
"""

from __future__ import annotations

from pathlib import Path


def run_goodvibes_on_mol(mol_name: str) -> dict[str, dict]:
    raise NotImplementedError(
        "Migrar desde results_collector.py:run_goodvibes_on_mol (1981-2181); "
        "revisar solapamiento con protonation_cycle.py:run_goodvibes (588-708) "
        "antes de decidir si son la misma función o dos usos distintos"
    )
