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
from typing import Optional


# Ubicaciones extra donde buscar si el binario no está en PATH ni en
# ningún env de conda. ORCA en particular casi nunca se instala vía
# conda (requiere licencia, se descarga como tarball y se descomprime
# donde a cada uno le resulte cómodo) — por eso necesita una lista de
# rutas típicas, a diferencia de crest/censo que sí suelen vivir en
# algún env.
_EXTRA_SEARCH_DIRS: list[Path] = [
    Path.home(),
    Path("/opt"),
    Path("/usr/local"),
    Path("/usr/local/bin"),
    Path("/Applications"),
]


def _env_var_override(name: str) -> Optional[Path]:
    """Permite fijar la ruta de un binario a mano vía variable de entorno,
    p.ej. MAAQC_ORCA_BIN=/Users/ieii/orca_6_1_1/orca. Tiene prioridad
    sobre cualquier otra búsqueda: si kmilo la define, se respeta tal
    cual (sin validar existencia acá — eso lo hace check_required_bins)."""
    env_key = f"MAAQC_{name.upper()}_BIN"
    value = os.environ.get(env_key)
    return Path(value).expanduser() if value else None


def _search_all_conda_envs(name: str) -> Optional[Path]:
    """Busca `<env>/bin/<name>` en TODOS los envs de conda, no en uno
    hardcodeado. Antes esto solo miraba un env fijo ("sciencebitch"),
    que no es el env real de kmilo en ninguna de sus dos máquinas — por
    eso `crest` (instalado en el env base, /opt/anaconda3/bin/crest) no
    se encontraba al correr `maaqc doctor` desde el env `maaqc`: conda
    saca el bin de base del PATH al activar otro env con nombre."""
    try:
        result = subprocess.run(
            ["conda", "env", "list", "--json"],
            capture_output=True, text=True, timeout=15,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0:
        return None
    try:
        envs = json.loads(result.stdout).get("envs", [])
    except json.JSONDecodeError:
        return None
    for env_path in envs:
        candidate = Path(env_path) / "bin" / name
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return candidate
    return None


def _search_extra_dirs(name: str, max_depth: int = 3) -> Optional[Path]:
    """Búsqueda superficial (BFS acotado a max_depth) en ubicaciones
    típicas de instalación manual — relevante sobre todo para ORCA, que
    no vive en conda. No sigue symlinks para evitar ciclos."""
    for base in _EXTRA_SEARCH_DIRS:
        if not base.is_dir():
            continue
        frontier = [base]
        depth = 0
        while frontier and depth < max_depth:
            next_frontier: list[Path] = []
            for d in frontier:
                try:
                    for child in d.iterdir():
                        if child.is_symlink():
                            continue
                        if child.is_dir():
                            next_frontier.append(child)
                        elif child.name == name and os.access(child, os.X_OK):
                            return child
                except (PermissionError, OSError):
                    continue
            frontier = next_frontier
            depth += 1
    return None


def resolve_bin(name: str) -> Path:
    """
    Resuelve un ejecutable a una ruta absoluta de forma robusta.

    Antes, CENSO_BIN/CREST_BIN eran solo el nombre del comando (p.ej. Path("censo")),
    así que su existencia dependía 100% de que el PATH del proceso que lanzó
    pipeline.py ya tuviera el env conda correcto activado. Si el script se
    invocaba desde un shell/cron/nohup sin el env correcto activado, el
    fallo solo se descubría al llamar subprocess.run() — a veces horas después
    (p.ej. tras 70 min de CREST), con un FileNotFoundError tardío.

    Busca en orden:
      1. Override explícito vía MAAQC_<NAME>_BIN (variable de entorno).
      2. PATH actual del proceso (shutil.which) — caso normal, env ya activado.
      3. TODOS los envs de conda vía `conda env list --json` (no un
         nombre de env hardcodeado), sin necesitar que estén activados.
      4. Ubicaciones típicas de instalación manual (~, /opt, /usr/local,
         /Applications), profundidad acotada — relevante para ORCA.
    Si no se encuentra en ningún lado, devuelve Path(name) tal cual (para no
    romper compatibilidad); el caller debe loggear/abortar temprano con
    check_required_bins().
    """
    override = _env_var_override(name)
    if override is not None:
        return override

    found = shutil.which(name)
    if found:
        return Path(found)

    from_conda = _search_all_conda_envs(name)
    if from_conda is not None:
        return from_conda

    from_extra = _search_extra_dirs(name)
    if from_extra is not None:
        return from_extra

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

    Un path se considera "encontrado" si es absoluto y apunta a un
    archivo ejecutable (caso normal: ya vino resuelto por resolve_bin),
    o si sigue siendo un nombre pelado que todavía resuelve por PATH
    (shutil.which) — cubre el caso de que alguien pase un Path("nombre")
    sin pasar por resolve_bin primero.
    """
    missing: list[tuple[str, Path]] = []
    for bin_name, bin_path in required.items():
        found = (
            bin_path.is_absolute()
            and bin_path.is_file()
            and os.access(bin_path, os.X_OK)
        ) or shutil.which(str(bin_path)) is not None
        if not found:
            missing.append((bin_name, bin_path))
    return missing
