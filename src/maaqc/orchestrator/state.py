"""Estado explícito por molécula/estado de protonación (state.json).

Reemplaza el reload manual de CENSO (--reload archivo.json, que exige
que TODOS los confórmeros de la ensemble actual estén presentes) por un
mecanismo de resume robusto ante configs cambiantes: hash de config +
checksums de outputs por stage.

Esto GENERALIZA dos mecanismos de recovery hoy independientes y ad hoc:

  1. recover_prescreening.py (script standalone) — reconstruye
     0_PRESCREENING.json parseando sp/sp.out + xtb_gsolv/{gas,solv}.out
     a mano, confórmero por confórmero, porque CENSO solo escribe el
     JSON consolidado si el proceso terminó sin interrupciones.

  2. pipeline.py:_censo_prescreen_recover (líneas 1342-1437) — matching
     por energía DFT sobre 1_SCREENING.out cuando los tokens de string
     esperados no coinciden (otro síntoma del mismo problema: CENSO no
     tiene "skip si el .out ya existe").

En maaqc, rebuild_censo_state() debe ser el comportamiento NORMAL de
CensoStage.is_complete() (ver stages/censo.py), no una herramienta que
hay que acordarse de correr manualmente cuando algo se corta.

Funciones a portar directo desde recover_prescreening.py (lógica de
parseo ya validada contra outputs reales, no reescribir):
  - parse_orca_sp_energy()      -> aquí o en io/parsing.py (ver TODO abajo)
  - parse_xtb_total_energy()    -> aquí o en io/parsing.py (ver TODO abajo)
  - natural_conf_sort_key()

TODO de diseño: parse_orca_sp_energy/parse_xtb_total_energy son parsers
puros de archivos ORCA/xTB — probablemente deberían vivir en
io/parsing.py junto a parse_sp_energy (que hace casi lo mismo, ver
duplicación ya detectada entre results_collector.py y
protonation_cycle.py) y no duplicarse de nuevo acá. rebuild_censo_state()
en este módulo debería importar esos parsers, no reimplementarlos.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class StageState:
    stage_name: str
    complete: bool
    config_hash: str
    output_checksums: dict[str, str]
    data: dict[str, Any]


def load_state(workdir: Path) -> dict[str, StageState]:
    """Lee state.json del workdir. Vacío si no existe (run nuevo)."""
    raise NotImplementedError("Definir formato de state.json antes de implementar")


def save_state(workdir: Path, state: dict[str, StageState]) -> None:
    raise NotImplementedError


def rebuild_censo_state(censo_dir: Path) -> dict[str, dict[str, float]]:
    """Reconstruye el estado de CENSO prescreening parseando .out en disco.

    Migrar la lógica de recover_prescreening.py::recover() acá,
    reutilizando parsers de io/parsing.py en vez de reimplementarlos.
    Debe reportar también qué confórmeros no se pudieron recuperar
    (comportamiento actual: lista `missing` + `errors` con razón).
    """
    raise NotImplementedError(
        "Migrar recover_prescreening.py::recover() + parse_orca_sp_energy/"
        "parse_xtb_total_energy, delegando el parseo a io/parsing.py"
    )
