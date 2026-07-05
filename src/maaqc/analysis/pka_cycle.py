"""Ciclo termodinámico de pKa (esquema directo CPCM, Pliego & Riveros 2002)
y equilibrio de tautómeros.

ATENCIÓN: hay DOS implementaciones de compute_pka en el legacy, y no son
la misma duplicación exacta que los parsers de io/parsing.py — antes de
migrar, confirmar con kmilo si son variantes con propósito distinto o
si una es la evolución de la otra:

  - results_collector.py:compute_pka()     líneas 933-1187
  - protonation_cycle.py:compute_pka()     líneas 1174-1509  (bastante más
                                             grande: incluye
                                             _validate_state_charges,
                                             constante de Pliego & Riveros,
                                             y soporte para GoodVibes vía
                                             _run_classic/_run_extended/
                                             _run_comparison en su propio
                                             main())

Resto a migrar desde protonation_cycle.py:
  - compute_tautomer_equilibrium()   líneas 1030-1066
  - _validate_state_charges()        líneas 344-405  (validación de
                                       consistencia de carga/N_átomos
                                       entre estados -- IMPORTANTE:
                                       preservar los mensajes '❌' que
                                       marcan error fatal vs. advertencia
                                       informativa, no colapsar a un
                                       solo nivel de severidad)
  - _state_G_total_direct()          líneas 966-1029
  - get_thermochem() / collect_mol_data()   líneas 709-965

Depende de io/parsing.py (parse_sp_energy, parse_thermochem,
parse_charge_mult_natoms) y de boltzmann.py — no reimplementar.
"""

from __future__ import annotations

from typing import Optional


def compute_pka(*args, **kwargs) -> dict:
    raise NotImplementedError(
        "Decidir primero cuál de las dos implementaciones (results_collector.py:933 "
        "vs protonation_cycle.py:1174) es la fuente de verdad antes de migrar"
    )


def compute_tautomer_equilibrium(*args, **kwargs) -> dict:
    raise NotImplementedError(
        "Migrar desde protonation_cycle.py:compute_tautomer_equilibrium (1030-1066)"
    )


def validate_state_charges(
    data_a: dict, data_b: dict, dq_expected: int, dn_expected: int,
    label_a: str, label_b: str,
) -> list[str]:
    raise NotImplementedError(
        "Migrar desde protonation_cycle.py:_validate_state_charges (344-405), "
        "preservando distinción '❌' (fatal) vs. informativo"
    )
