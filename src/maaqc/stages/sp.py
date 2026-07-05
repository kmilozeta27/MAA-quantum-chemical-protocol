"""SpStage — single point fase gas (etapa OBLIGATORIA del pipeline).

Migrar desde pipeline.py:
  - stage_sp()            líneas 2672-2747
  - _sp_worker()          líneas 3201-3258
  - stage_sp_nroot()      líneas 2881-2971  (modo --nroot: SP+NRoots1 sobre
                                              top-N confórmeros por E_OPT;
                                              es un MODO de este stage, no
                                              un stage aparte -- análogo a
                                              multitemp en censo.py)
  - extract_sp_energy()   líneas 913-928

TODO de diseño (no resuelto todavía, discutir antes de migrar): con
--skipnumfreq, SP depende directo de OptStage; sin --skipnumfreq depende
de NumfreqStage. Eso significa que `requires` no puede ser una lista
estática de clase si se quiere mantener el contrato Stage tal cual está
en base.py. Opciones: (a) requires siempre incluye ambos y
is_complete()/run() deciden cuál input usar según config.skip_numfreq,
o (b) requires se vuelve un método en vez de atributo de clase. Definir
esto en orchestrator/dag.py antes de migrar sp.py en serio.
"""

from __future__ import annotations

from pathlib import Path

from maaqc.config.schema import RecipeConfig
from maaqc.stages.base import Stage, StageResult


class SpStage(Stage):
    name = "sp"
    requires = ["numfreq"]  # TODO: ver nota de diseño arriba (skip_numfreq)

    def is_complete(self, workdir: Path) -> bool:
        raise NotImplementedError("Migrar desde pipeline.py:stage_sp (líneas 2672-2747)")

    def inputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def outputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        raise NotImplementedError
