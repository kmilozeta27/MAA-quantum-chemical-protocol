"""OptStage — optimización de geometría DFT sobre confórmeros seleccionados.

Migrar desde pipeline.py:
  - stage_opt()          líneas 2462-2600
  - _opt_worker()        líneas 3136-3200  (worker de ThreadPoolExecutor)
  - resolve_opt_xyz()    líneas 1016-1038
  - check_success() / check_no_imaginary()  líneas 878-902

Candidato sugerido como PRIMER stage a migrar completo (junto con sp.py):
es el más autocontenido, sin las complicaciones de recovery de censo.py
ni el árbol de modos de adc2.py. Buen lugar para validar el contrato
Stage corriendo legacy y maaqc en paralelo sobre el mismo input y
comparando outputs (energías OPT, geometrías) antes de confiar en la
migración.
"""

from __future__ import annotations

from pathlib import Path

from maaqc.config.schema import RecipeConfig
from maaqc.stages.base import Stage, StageResult


class OptStage(Stage):
    name = "opt"
    requires = ["censo"]

    def is_complete(self, workdir: Path) -> bool:
        raise NotImplementedError("Migrar desde pipeline.py:stage_opt (líneas 2462-2600)")

    def inputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def outputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        raise NotImplementedError
