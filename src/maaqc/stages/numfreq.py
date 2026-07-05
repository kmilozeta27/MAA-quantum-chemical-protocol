"""NumfreqStage — frecuencias vibracionales numéricas (termoquímica).

Migrar desde pipeline.py:
  - stage_numfreq()          líneas 2601-2671
  - has_numfreq_restart()    líneas 903-912

Nota: config.skip_numfreq (= --skipnumfreq legacy) hace que el
orchestrator salte este stage completo (OPT -> SP -> TDDFT directo),
no que NumfreqStage se ejecute como no-op. Esa decisión vive en
orchestrator/dag.py, no acá.
"""

from __future__ import annotations

from pathlib import Path

from maaqc.config.schema import RecipeConfig
from maaqc.stages.base import Stage, StageResult


class NumfreqStage(Stage):
    name = "numfreq"
    requires = ["opt"]

    def is_complete(self, workdir: Path) -> bool:
        raise NotImplementedError("Migrar desde pipeline.py:stage_numfreq (líneas 2601-2671)")

    def inputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def outputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        raise NotImplementedError
