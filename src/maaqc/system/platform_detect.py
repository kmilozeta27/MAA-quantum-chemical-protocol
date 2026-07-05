"""Detección de plataforma (ThinkPad Ubuntu vs MacBook Air institucional)
y resolución robusta de binarios (ORCA/CREST/CENSO).

Migrar desde pipeline.py:
  - _resolve_bin()          líneas 93-129
  - check_required_bins()   líneas 132-158

check_required_bins() es el cuerpo de `maaqc doctor` (ver cli.py) — es
el subcomando más simple del CLI nuevo, buen primer PR de CLI real una
vez migrado esto.

TODO: agregar detección explícita de macOS (setup.py ya tenía
auto-detección para el MacBook Air institucional, usuario `ieii`) —
el legacy resuelve esto en setup.py, no en pipeline.py; confirmar con
kmilo si ese setup.py también debe migrarse acá o si queda fuera del
alcance de este paquete (es infra de instalación, no del pipeline).
"""

from __future__ import annotations

import shutil
from pathlib import Path


def resolve_bin(name: str, conda_env: str = "sciencebitch") -> Path:
    raise NotImplementedError("Migrar desde pipeline.py:_resolve_bin (líneas 93-129)")


def check_required_bins(required: dict[str, Path]) -> list[tuple[str, Path]]:
    """Retorna lista de (nombre, path) de binarios faltantes (vacía si todo ok).

    Migrar desde pipeline.py:check_required_bins (líneas 132-158), pero
    retornando la lista en vez de hacer log.error + SystemExit directo
    — eso pasa a ser responsabilidad de `maaqc doctor` en cli.py.
    """
    raise NotImplementedError
