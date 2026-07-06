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

import json
import os
import shutil
import subprocess
from pathlib import Path


def resolve_bin(name: str, conda_env: str = "sciencebitch") -> Path:
    """
    Resuelve un ejecutable a una ruta absoluta de forma robusta.

    Antes, CENSO_BIN/CREST_BIN eran solo el nombre del comando (p.ej. Path("censo")),
    así que su existencia dependía 100% de que el PATH del proceso que lanzó
    pipeline.py ya tuviera el env conda correcto activado. Si el script se
    invocaba desde un shell/cron/nohup sin 'conda activate sciencebitch', el
    fallo solo se descubría al llamar subprocess.run() — a veces horas después
    (p.ej. tras 70 min de CREST), con un FileNotFoundError tardío.

    Busca en orden:
      1. PATH actual del proceso (shutil.which) — caso normal, env ya activado.
      2. <prefix de conda_env>/bin/<name> — vía `conda env list --json`,
         sin necesitar que el env esté activado.
    Si no se encuentra en ningún lado, devuelve Path(name) tal cual (para no
    romper compatibilidad); el caller debe loggear/abortar temprano con
    check_required_bins().
    """
    found = shutil.which(name)
    if found:
        return Path(found)
    try:
        result = subprocess.run(
            ["conda", "env", "list", "--json"],
            capture_output=True, text=True, timeout=15,
        )
        if result.returncode == 0:
            envs = json.loads(result.stdout).get("envs", [])
            for env_path in envs:
                if Path(env_path).name == conda_env:
                    candidate = Path(env_path) / "bin" / name
                    if candidate.is_file() and os.access(candidate, os.X_OK):
                        return candidate
    except Exception:
        pass
    return Path(name)


def check_required_bins(required: dict[str, Path]) -> list[tuple[str, Path]]:
    """Retorna lista de (nombre, path) de binarios faltantes (vacía si todo ok).

    Migrado desde pipeline.py:check_required_bins (líneas 132-158), pero
    retornando la lista en vez de hacer log.error + SystemExit directo
    — eso pasa a ser responsabilidad de `maaqc doctor` en cli.py (decidir
    si abortar, qué mensaje mostrar al usuario, etc.).

    `required` es un dict {nombre_binario: path_resuelto}, típicamente
    construido con resolve_bin(), p.ej.:
        check_required_bins({"CREST": resolve_bin("crest"),
                              "CENSO": resolve_bin("censo")})
    """
    missing: list[tuple[str, Path]] = []
    for bin_name, bin_path in required.items():
        resolved = shutil.which(str(bin_path)) or (
            str(bin_path) if bin_path.is_absolute() and bin_path.is_file() else None
        )
        if resolved is None:
            missing.append((bin_name, bin_path))
    return missing
