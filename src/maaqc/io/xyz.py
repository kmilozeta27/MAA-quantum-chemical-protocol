"""Lectura/escritura de archivos .xyz multi-frame (confórmeros CREST/CENSO).

Migrar desde pipeline.py:
  - atoms_to_xyz_block()          líneas 687-694
  - read_conformers()             líneas 695-856
  - _extract_conf_ids_from_xyz()  líneas 2415-2461

Depende de ase.io (read/write) tal como el legacy — mantener esa
dependencia, no reemplazarla por un parser manual salvo que ASE falle
en algún caso concreto que se documente acá.
"""

from __future__ import annotations

from pathlib import Path


def read_conformers(xyz_file: Path, nconf: int) -> list:
    raise NotImplementedError("Migrar desde pipeline.py:read_conformers (líneas 695-856)")


def atoms_to_xyz_block(atoms) -> str:
    raise NotImplementedError("Migrar desde pipeline.py:atoms_to_xyz_block (líneas 687-694)")
