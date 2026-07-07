"""Test de humo para maaqc.analysis.

Motivación (D0): boltzmann.py quedó con un import roto a ..constants
durante varias sesiones de migración porque ningún test importaba
siquiera el paquete analysis/ — el error solo aparecía al ejecutar el
pipeline completo. Este archivo existe únicamente para que "el módulo
carga sin ImportError/AttributeError" quede cubierto por CI, sin esperar
a un test funcional de cada stage.

Alcance deliberadamente mínimo: NO prueba que las funciones hagan lo que
dicen (para eso hacen falta fixtures de datos ORCA/CENSO reales, fuera
de alcance acá). Solo prueba que:
  1. el módulo importa sin reventar (sintaxis, imports circulares,
     nombres que no existen en ..constants, etc.)
  2. los símbolos públicos que otros módulos migrados ya importan
     explícitamente (p.ej. pka_cycle.py importa GV_CUTOFF/GV_QH/GV_QS/
     run_goodvibes desde goodvibes.py) siguen existiendo con ese nombre.

Lista de módulos: se mantiene explícita (no se descubre con pkgutil)
a propósito, para que agregar un módulo nuevo a analysis/ sin agregarlo
acá se note en review — un test que "descubre solo" los módulos nunca
puede fallar por un módulo faltante en la lista.

Estado de cada módulo al momento de escribir este test (ver conversación
de migración maa-repackaging para detalle línea por línea):
  - pka_cycle.py     : migrado completo
  - boltzmann.py     : migrado completo
  - goodvibes.py     : migrado completo (parte 1 del Issue C1; la
                       segunda mitad -- run_goodvibes_on_mol -- sigue
                       siendo NotImplementedError a propósito)
  - adc2_compare.py  : stub, compare_adc2_tddft() -> NotImplementedError
  - monte_carlo.py   : stub vacío -- ni siquiera tiene funciones
                       (ver docstring del propio archivo: puede que ni
                       corresponda a código real, a confirmar con
                       kmilozeta). Igual debe importar limpio.
  - parse_uvvis.py   : stub, extract_lambda_max() -> NotImplementedError
  - plots.py         : parcialmente migrado (plot_distribution,
                       plot_distribution_5, plot_conformer_boltzmann
                       reales; plot_spectrum, plot_comparison,
                       plot_adc2_tddft_comparison -> NotImplementedError)
"""

from __future__ import annotations

import importlib

import pytest

# Nombres relativos a maaqc.analysis. Si se agrega un módulo nuevo al
# paquete y no aparece acá, este test NO lo cubre -- agregarlo a mano.
ANALYSIS_MODULES = [
    "pka_cycle",
    "boltzmann",
    "goodvibes",
    "adc2_compare",
    "monte_carlo",
    "parse_uvvis",
    "plots",
]


@pytest.mark.parametrize("module_name", ANALYSIS_MODULES)
def test_analysis_module_imports(module_name: str) -> None:
    """Cada módulo de analysis/ debe importar sin lanzar excepciones.

    Parametrizado (en vez de un loop con un solo test) para que si un
    módulo falla, pytest reporte cuál puntualmente en vez de abortar
    todo el archivo en el primer import roto -- eso es exactamente lo
    que hubiera atrapado el bug de D0 sin ambigüedad sobre cuál módulo
    era el culpable.
    """
    importlib.import_module(f"maaqc.analysis.{module_name}")


def test_analysis_package_itself_imports() -> None:
    """El paquete analysis/ (su __init__.py) también debe cargar solo."""
    importlib.import_module("maaqc.analysis")


def test_goodvibes_symbols_used_by_pka_cycle_exist() -> None:
    """pka_cycle.py importa estos 4 nombres explícitamente desde
    goodvibes.py (from .goodvibes import GV_CUTOFF, GV_QH, GV_QS,
    run_goodvibes). Si alguno se renombra o se borra en goodvibes.py sin
    actualizar pka_cycle.py, ESTE test debe fallar acá -- más claro que
    esperar a que falle el import de pka_cycle.py y haya que adivinar
    cuál de los 4 símbolos fue.
    """
    goodvibes = importlib.import_module("maaqc.analysis.goodvibes")
    for symbol in ("GV_CUTOFF", "GV_QH", "GV_QS", "run_goodvibes"):
        assert hasattr(goodvibes, symbol), (
            f"maaqc.analysis.goodvibes.{symbol} no existe; "
            f"pka_cycle.py depende de este nombre exacto."
        )


def test_constants_symbols_used_by_boltzmann_exist() -> None:
    """boltzmann.py importa estos 5 nombres desde ..constants. Mismo
    razonamiento que el test anterior, para el otro lado del bug de D0
    (el import que efectivamente estaba roto)."""
    constants = importlib.import_module("maaqc.constants")
    for symbol in ("HARTREE_TO_KCAL", "KB_KCAL", "TEMP_LABELS", "TEMPS_K", "T_STD"):
        assert hasattr(constants, symbol), (
            f"maaqc.constants.{symbol} no existe; boltzmann.py depende "
            f"de este nombre exacto."
        )
