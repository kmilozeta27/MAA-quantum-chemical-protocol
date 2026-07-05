"""stages/base.py — Contrato común para todos los stages del pipeline.

Decisión de diseño (ver discusión adc2.py / dag.py): `requires` se queda
como atributo ESTÁTICO de clase en los 7 stages, sin excepción — incluido
Adc2Stage, que en apariencia necesitaría requires dinámico según
config.entry. Quien resuelve "qué validar según el --entry" es
`orchestrator/dag.py`, no `Stage`. Esto evita reabrir este contrato cada
vez que aparezca un --entry nuevo con el patrón "saltar cadena completa".
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from maaqc.config.schema import RecipeConfig


@dataclass
class StageResult:
    ok: bool
    message: str = ""
    data: dict[str, Any] = field(default_factory=dict)


class Stage(ABC):
    """Interfaz que todo stage debe implementar.

    Atributos de clase:
        name:     nombre corto del stage (ej. "adc2"), usado en logs,
                  state.json y en `requires` de otros stages.
        requires: nombres de stages previos que deben estar is_complete()
                  en la cadena NORMAL. Lista estática — no evaluar config
                  acá. Casos "a veces existe, a veces no" (--entry
                  standalone) se resuelven en dag.py, no acá.
    """

    name: str
    requires: list[str] = []

    @abstractmethod
    def is_complete(self, workdir: Path) -> bool:
        """True si este stage ya terminó exitosamente para `workdir`."""
        ...

    @abstractmethod
    def inputs(self, workdir: Path) -> list[Path]:
        """Archivos/directorios que este stage necesita leer para correr.

        Se valida SIEMPRE, sin importar el --entry (ver dag.py). Esta es
        la fuente de verdad para "¿tengo lo que necesito en disco?",
        independiente de si los stages previos de la cadena normal
        corrieron o no en este run.
        """
        ...

    @abstractmethod
    def outputs(self, workdir: Path) -> list[Path]:
        """Artefactos que este stage produce cuando termina bien."""
        ...

    @abstractmethod
    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        ...
