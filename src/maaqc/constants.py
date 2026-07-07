"""Constantes físicas y numéricas compartidas por todo el paquete maaqc.

Fuente única de verdad. Antes de esta migración estas constantes estaban
duplicadas de forma idéntica en pka_cycle.py, parsing.py y goodvibes.py
(HARTREE_TO_KCAL, T_STD), solo definidas en pka_cycle.py sin exportar
(TEMPS_K, TEMP_LABELS), o repetidas bajo dos nombres distintos para la
misma cantidad física (R_KCAL en pka_cycle.py / KB_KCAL en boltzmann.py).

Historial de cada constante (para trazabilidad, no borrar):

  HARTREE_TO_KCAL : 627.5094740631
      Vivía idéntica en io/parsing.py:29 y analysis/pka_cycle.py:98.

  T_STD : 298.15  (K)
      Vivía idéntica en io/parsing.py:30, analysis/pka_cycle.py:101 y
      analysis/goodvibes.py (T_STD local).

  TEMPS_K / TEMP_LABELS
      Vivían solo en analysis/pka_cycle.py:192-193, sin exportar.
      boltzmann.py ya las importaba desde ..constants (import adelantado,
      a la espera de que este archivo existiera).

  KB_KCAL : 1.9872041e-3  (kcal·mol⁻¹·K⁻¹)
      MISMA constante que pka_cycle.py llamaba R_KCAL (línea 99, usada en
      líneas 266, 270, 271, 427, 434, 435, 564, 681). kB expresado por mol
      (kB·N_A) es numéricamente la constante de los gases R — no son dos
      constantes, es una sola cantidad física con dos nombres históricos
      distintos según qué módulo la escribió primero. Se resuelve el
      duplicado a favor de KB_KCAL (nombre que boltzmann.py ya usaba) y
      pka_cycle.py se actualiza para llamarla igual — ver notas de edición
      en pka_cycle.py. No se deja alias R_KCAL acá a propósito: mantenerlo
      solo reintroduciría el problema original (dos nombres válidos para
      lo mismo, con riesgo de que alguien defina un tercer valor bajo uno
      de los dos nombres sin tocar el otro).
"""

from __future__ import annotations

from pathlib import Path

# Directorio raíz de corridas por DEFAULT (legacy: Path("runs") global fijo).
#
# NO es equivalente al WORK_ROOT que analysis/results_collector.py importaba
# antes de esta migración. Ese uso quedó reemplazado: collect_results() y
# compare_protonation_states() ahora reciben workdir/run_root como parámetro
# explícito (mismo patrón que orchestrator/state.py y el contrato de Stage:
# is_complete(workdir), outputs(workdir), etc.), en vez de derivarlo de una
# constante de módulo — necesario para correr varios estados de protonación
# (A/B/C/D/R) en paralelo con workdirs distintos.
#
# DEFAULT_RUN_ROOT existe solo como valor por default para el CLI (cli.py):
# si el usuario no pasa --workdir/--run-root explícito, el CLI resuelve acá
# y pasa el Path resultante a las funciones de análisis. Ningún módulo bajo
# analysis/ debería importar esta constante directamente.
DEFAULT_RUN_ROOT: Path = Path("runs")

# Conversión de energía
HARTREE_TO_KCAL: float = 627.5094740631

# Temperatura de referencia estándar (K)
T_STD: float = 298.15

# Rango de temperaturas para el análisis multi-T (Boltzmann, pKa)
TEMPS_K: list[float] = [250.0, 274.15, 298.15, 310.0]
TEMP_LABELS: list[str] = ["250K", "274K", "298K", "310K"]

# Constante de Boltzmann expresada por mol (kcal·mol⁻¹·K⁻¹).
# Numéricamente igual a la constante de los gases R, porque kB·N_A = R.
# Nombre unificado: pka_cycle.py la llamaba R_KCAL antes de esta migración;
# ver docstring arriba.
KB_KCAL: float = 1.9872041e-3
