"""CensoStage — CENSO prescreening (obligatorio entre CREST y OPT).

Migrar desde pipeline.py:
  - stage_censo_prescreen()      líneas 1438-1881  (modo normal, single T)
  - stage_censo_multitemp()      líneas 2019-2107  (modo --multitemp, corre
                                                     CENSO N veces, un
                                                     subdirectorio T<T>K/
                                                     por temperatura)
  - parse_censo_screening_out()  líneas 1882-1932
  - extract_conformers_by_id()   líneas 1933-2018
  - select_nconf_auto()          líneas 2277-2414
  - select_nconf_boltzmann()     líneas 2138-2276  (usa --boltzmann-threshold)

multitemp es un MODO de este stage (config.censo.multitemp), no un stage
aparte: ambos producen el mismo tipo de output (ensemble filtrado listo
para OPT), solo difieren en cuántas veces corre CENSO y dónde queda cada
resultado.

is_complete() ES el caso de uso central de orchestrator/state.py: no
debe confiar en que CENSO haya escrito 0_PRESCREENING.json (solo lo hace
si el proceso terminó sin interrupciones — ver recover_prescreening.py
original y pipeline.py:_censo_prescreen_recover, líneas 1342-1437).
En su lugar debe reconstruir el estado por confórmero parseando
directamente sp/sp.out + xtb_gsolv/{gas,solv}.out, delegando en
orchestrator.state.rebuild_censo_state().
"""

from __future__ import annotations

from pathlib import Path

from maaqc.config.schema import RecipeConfig
from maaqc.stages.base import Stage, StageResult


class CensoStage(Stage):
    name = "censo"
    requires = ["crest"]

    def is_complete(self, workdir: Path) -> bool:
        # TODO: delegar en orchestrator.state.rebuild_censo_state(workdir)
        # en vez de buscar 0_PRESCREENING.json directamente.
        raise NotImplementedError(
            "Migrar desde pipeline.py:stage_censo_prescreen + "
            "_censo_prescreen_recover (líneas 1342-1881)"
        )

    def inputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def outputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        raise NotImplementedError
