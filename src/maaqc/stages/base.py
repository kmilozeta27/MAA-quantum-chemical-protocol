"""Contrato común que todo stage debe implementar.

Esto es lo que reemplaza flags hardcodeados como --start adc2,
--skipnumfreq, --expand-crest por detección automática de qué ya está
satisfecho en disco. Cada stage_* de pipeline.py (legacy) se convierte
en una subclase de Stage:

    stage_crest            -> CrestStage
    stage_censo_prescreen  -> CensoStage   (+ stage_censo_multitemp como
                                             modo del mismo stage, no un
                                             stage aparte -- ver censo.py)
    stage_opt              -> OptStage
    stage_numfreq          -> NumfreqStage
    stage_sp               -> SpStage      (stage_sp_nroot es un modo,
                                             no un stage nuevo -- ver sp.py)
    stage_tddft            -> TddftStage
    stage_adc2             -> Adc2Stage

`is_complete()` es la pieza que generaliza recover_prescreening.py y
pipeline.py:_censo_prescreen_recover — cada stage debe poder reconstruir
su propio estado leyendo los .out en disco, en vez de depender de un
JSON consolidado que el programa externo (CENSO) solo escribe si el
proceso terminó sin interrupciones.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from maaqc.config.schema import RecipeConfig


@dataclass
class StageResult:
    """Resultado de correr (o verificar) un stage."""

    stage_name: str
    success: bool
    outputs: list[Path] = field(default_factory=list)
    # Datos estructurados específicos del stage (ej. energías por confórmero,
    # confórmero de menor energía para ADC2, etc.) — cada stage define su
    # propia forma; el orchestrator solo los pasa entre stages sin interpretarlos.
    data: dict[str, Any] = field(default_factory=dict)
    message: str = ""


class Stage(ABC):
    """Interfaz común para todos los stages del pipeline.

    Contrato acordado en la skill maa-repackaging — no introducir un
    patrón ad hoc por stage; si un stage nuevo necesita algo que esta
    interfaz no cubre, extenderla acá (afecta a todos los stages) en
    vez de resolverlo solo en la subclase.
    """

    name: str
    requires: list[str]  # nombres de stages previos requeridos, ej. ["censo"]

    @abstractmethod
    def is_complete(self, workdir: Path) -> bool:
        """True si este stage ya produjo outputs válidos y reusables.

        Debe reconstruir el estado leyendo archivos reales en disco
        (.out, .xyz, etc.), no confiar únicamente en un marcador o JSON
        consolidado que un proceso externo pudo no haber escrito si fue
        interrumpido (ver recover_prescreening.py para el caso CENSO).
        """
        raise NotImplementedError

    @abstractmethod
    def inputs(self, workdir: Path) -> list[Path]:
        """Archivos que este stage necesita que ya existan antes de correr."""
        raise NotImplementedError

    @abstractmethod
    def outputs(self, workdir: Path) -> list[Path]:
        """Archivos que este stage produce cuando corre exitosamente."""
        raise NotImplementedError

    @abstractmethod
    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        """Ejecuta el stage. Asume que requires ya están is_complete()."""
        raise NotImplementedError
