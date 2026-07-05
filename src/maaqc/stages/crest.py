"""CrestStage — búsqueda conformacional GFN2-xTB.

Migrar desde pipeline.py:stage_crest() (líneas 1152-1233 en el legacy).
No tiene stages previos (requires = []).

Puntos a preservar al migrar:
  - Uso de CREST_BIN resuelto vía _resolve_bin() -> system/platform_detect.py.
  - Timeout CREST_TIMEOUT_HOURS (24h default).
  - Energías embebidas en el multi-XYZ de salida se parsean después con
    io/parsing.py::parse_crest_energies (NO reimplementar el parseo acá).
"""

from __future__ import annotations

from pathlib import Path

from maaqc.config.schema import RecipeConfig
from maaqc.stages.base import Stage, StageResult


class CrestStage(Stage):
    name = "crest"
    requires: list[str] = []

    def is_complete(self, workdir: Path) -> bool:
        raise NotImplementedError("Migrar desde pipeline.py:stage_crest (líneas 1152-1233)")

    def inputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def outputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        raise NotImplementedError
