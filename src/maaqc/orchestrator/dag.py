"""DAG de dependencias entre stages + lógica de resume/entry points.

Reemplaza:
  - pipeline.py:_run_one()          líneas 3425-3849 (decide qué stages
                                     correr según --start/--skipnumfreq/etc.)
  - pipeline.py:_split_argv_by_mol() líneas 3850-3888 (modo multi-molécula:
                                     esto se vuelve una cola de runs del
                                     orchestrator, no parsing de argv)

Responsable de resolver los dos TODOs de diseño dejados pendientes en
stages/sp.py y stages/adc2.py: `requires` no es siempre una lista
estática, depende de config (skip_numfreq) o del --entry elegido
(sp-only ignora todo lo previo a 04_sp/). La resolución de qué stages
correr para un --entry dado debe vivir acá, como una función que mapea
(entry, config) -> lista de Stage a ejecutar en orden — no como
lógica ad hoc dentro de cada Stage.is_complete().

STAGE_ORDER y STAGE_REGISTRY reemplazan la lista plana STEPS =
["crest", "opt", "numfreq", "sp", "tddft", "adc2"] de pipeline.py
(línea 218) -- nótese que esa lista legacy ya ni siquiera incluye
"censo", otro síntoma de que el orden real vivía implícito en el código
de _run_one() en vez de en una estructura de datos explícita.
"""

from __future__ import annotations

from maaqc.config.schema import RecipeConfig
from maaqc.stages.adc2 import Adc2Stage
from maaqc.stages.base import Stage
from maaqc.stages.censo import CensoStage
from maaqc.stages.crest import CrestStage
from maaqc.stages.numfreq import NumfreqStage
from maaqc.stages.opt import OptStage
from maaqc.stages.sp import SpStage
from maaqc.stages.tddft import TddftStage

STAGE_REGISTRY: dict[str, type[Stage]] = {
    "crest": CrestStage,
    "censo": CensoStage,
    "opt": OptStage,
    "numfreq": NumfreqStage,
    "sp": SpStage,
    "tddft": TddftStage,
    "adc2": Adc2Stage,
}

STAGE_ORDER: list[str] = ["crest", "censo", "opt", "numfreq", "sp", "tddft", "adc2"]


def resolve_stages_for_entry(entry: str, config: RecipeConfig) -> list[str]:
    """Mapea --entry + config a la lista de stages a ejecutar, en orden.

    TODO: implementar los 5 casos de pipeline.py (líneas 150-161 del
    esqueleto de la skill):
      - smiles          -> cadena completa desde crest
      - xyz             -> cadena completa desde crest
      - crest-ensemble  -> salta crest (ya hay confórmeros), empieza en censo
      - sp-only         -> ignora crest/censo/opt/numfreq/tddft, solo
                            requiere 04_sp/ (equivalente a --start adc2)
      - calibrated      -> como xyz/crest-ensemble pero usa config.reference
                            para shifts empíricos (ej. shinorine) en adc2

    Además debe aplicar config.skip_numfreq (salta "numfreq" del resultado).
    """
    raise NotImplementedError


def run_pipeline(config: RecipeConfig, workdir) -> None:
    """Corre (o retoma) el pipeline completo según resolve_stages_for_entry()."""
    raise NotImplementedError
