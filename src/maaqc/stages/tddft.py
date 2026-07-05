"""TddftStage — TD-DFT sobre los confórmeros que completaron SP.

Migrar desde pipeline.py:
  - stage_tddft()       líneas 2814-2880
  - _tddft_worker()     líneas 2748-2813

Corre después de SP salvo que config.entry sea posterior en la cadena
(ej. --entry sp-only detiene antes de llegar acá).
"""

from __future__ import annotations

from pathlib import Path

from maaqc.config.schema import RecipeConfig
from maaqc.stages.base import Stage, StageResult


class TddftStage(Stage):
    name = "tddft"
    requires = ["sp"]

    def is_complete(self, workdir: Path) -> bool:
        raise NotImplementedError("Migrar desde pipeline.py:stage_tddft (líneas 2814-2880)")

    def inputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def outputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        raise NotImplementedError
