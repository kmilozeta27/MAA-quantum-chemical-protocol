"""Lectura/escritura de archivos .xyz multi-frame (confórmeros CREST/CENSO).

Migrado desde pipeline.py (vía system/_orca_shared.py, donde vivían
duplicadas mientras se portaban opt.py/sp.py — ver CHANGELOG):
  - atoms_to_xyz_block()          líneas 687-694 (legacy)
  - read_all_conformers()         líneas 695-856 (legacy: read_conformers)
  - extract_conf_ids_from_xyz()   líneas 2415-2461 (legacy)

Depende de ase.io (read/write) tal como el legacy — mantener esa
dependencia, no reemplazarla por un parser manual salvo que ASE falle
en algún caso concreto que se documente acá.

Fuente única de verdad
-----------------------
Este módulo es la fuente única de verdad para estas 3 funciones.
`system/_orca_shared.py` (usado por OptStage/SpStage) las importa desde
acá en vez de definir su propia copia — antes había dos copias
idénticas + una tercera (`extract_conf_ids_from_xyz`) que solo estaba
documentada, no declarada, en este archivo. Si algo necesita I/O de
XYZ, importa de `maaqc.io.xyz`, no de `maaqc.system._orca_shared`.

Decisión de diseño: por qué NO se migró `read_conformers(xyz_file, nconf)`
---------------------------------------------------------------------------
El legacy truncaba a `nconf` confórmeros dentro de la misma función de
lectura. Esa firma se descarta deliberadamente (no es un stub pendiente,
es un "no vamos a portar esto tal cual"): en la arquitectura actual la
responsabilidad de decidir CUÁNTOS confórmeros pasan a OPT ya no es de
la etapa de lectura, es de CENSO (`select_nconf_auto` /
`select_nconf_boltzmann` en `system/_censo_shared.py`, invocadas desde
`stages/censo.py`). Mezclar de nuevo el truncamiento acá reintroduciría
la misma responsabilidad en dos lugares. Por eso `read_all_conformers()`
lee siempre el ensemble completo; quien la llama (típicamente
`OptStage`) ya recibe un XYZ pre-filtrado por CENSO y no necesita
truncar de nuevo.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional

from ase.io import read as ase_read


def atoms_to_xyz_block(atoms) -> str:
    """Formatea un objeto `ase.Atoms` como bloque de coordenadas para
    un input de ORCA (sin línea de conteo ni comentario, solo
    `símbolo  x  y  z` por átomo, alineado a 4 decimales... 8 en
    realidad, ver formato abajo). Portado 1:1 desde pipeline.py:687-694."""
    lines = []
    for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        lines.append(f"  {symbol:<4s}  {pos[0]:14.8f}  {pos[1]:14.8f}  {pos[2]:14.8f}")
    return "\n".join(lines)


def read_all_conformers(xyz_file: Path) -> list:
    """Lee todos los frames de un multi-XYZ (CREST/CENSO) como lista de
    `ase.Atoms`. A diferencia del `read_conformers(xyz_file, nconf)`
    legacy, NO trunca — ver nota de diseño en el docstring del módulo."""
    try:
        all_conf = ase_read(str(xyz_file), index=":")
    except Exception as e:
        raise RuntimeError(f"No se pudo leer {xyz_file}: {e}") from e
    return all_conf if isinstance(all_conf, list) else [all_conf]


def extract_conf_ids_from_xyz(source_xyz: Path) -> dict[int, str]:
    """{índice_0based → 'CONF#'} leyendo los comentarios extXYZ/plain del
    multi-XYZ. Preserva la trazabilidad CONF39 → conf_39 en vez de
    re-etiquetar conf_001, conf_002, ... (ver docstring original en
    pipeline.py, función `_extract_conf_ids_from_xyz`, líneas 2415-2461)."""
    result: dict[int, str] = {}
    try:
        raw_lines = source_xyz.read_text(errors="replace").splitlines()
    except OSError:
        return result

    frame_idx = 0
    line_idx = 0
    n_atoms: Optional[int] = None

    while line_idx < len(raw_lines):
        stripped = raw_lines[line_idx].strip()
        if stripped.isdigit():
            nat = int(stripped)
            if n_atoms is None:
                n_atoms = nat
            if nat == n_atoms and line_idx + 1 < len(raw_lines):
                comment = raw_lines[line_idx + 1].strip()
                m = re.search(r"(CONF\d+)=", comment)
                if m:
                    result[frame_idx] = m.group(1)
                else:
                    for tok in comment.split():
                        if tok.startswith("CONF"):
                            result[frame_idx] = tok
                            break
                frame_idx += 1
                line_idx += nat + 2
                continue
        line_idx += 1

    return result
