"""Tests para maaqc.io.xyz — usa tests/fixtures/glycine.xyz (1 frame)."""

from pathlib import Path

import pytest

from maaqc.io.xyz import (
    atoms_to_xyz_block,
    extract_conf_ids_from_xyz,
    read_all_conformers,
)

FIXTURE = Path(__file__).parent.parent / "fixtures" / "glycine.xyz"


def test_read_all_conformers_returns_list():
    confs = read_all_conformers(FIXTURE)
    assert isinstance(confs, list)
    assert len(confs) == 1
    assert len(confs[0]) == 10  # 10 átomos, ver header del .xyz


def test_read_all_conformers_missing_file_raises():
    with pytest.raises(RuntimeError):
        read_all_conformers(Path("/no/existe.xyz"))


def test_atoms_to_xyz_block_format():
    confs = read_all_conformers(FIXTURE)
    block = atoms_to_xyz_block(confs[0])
    lines = block.splitlines()
    assert len(lines) == 10
    # símbolo + 3 floats con 8 decimales, sin línea de conteo/comentario
    first = lines[0].split()
    assert first[0] == "N"
    assert len(first) == 4
    float(first[1])  # no debe lanzar


def test_extract_conf_ids_from_xyz_no_conf_tag():
    # glycine.xyz no trae 'CONF#=' en el comentario -> dict vacío, no crash
    result = extract_conf_ids_from_xyz(FIXTURE)
    assert result == {}


def test_extract_conf_ids_from_xyz_missing_file_returns_empty():
    result = extract_conf_ids_from_xyz(Path("/no/existe.xyz"))
    assert result == {}
