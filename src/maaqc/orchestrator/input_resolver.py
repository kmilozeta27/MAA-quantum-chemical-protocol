"""orchestrator/input_resolver.py — Materialización del input inicial de
la cadena, según `config.entry`.

Contexto (issue "materialización del input inicial de la cadena", Problema 1)
-------------------------------------------------------------------------
`cmd_run`/`cmd_resume` nunca escribían el archivo que `CrestStage` espera
(`workdir/input.xyz`, ver `stages/crest.py`, pregunta abierta #1 de su
docstring). `args.input` se leía en el argparse pero solo se usaba para
calcular el `workdir` default — el valor en sí se perdía.

Este módulo resuelve exactamente eso: dado `args.input` (+ `args.ensemble`
para `--entry crest-ensemble`) y el `config.entry` ya resuelto, deja en
`workdir` lo que el primer stage de la cadena resuelta por
`orchestrator/dag.py::resolve_stages_for_entry()` necesita encontrar.

Deliberadamente NO vive acá:
  - Orden de stages / qué stage es "el primero" → orchestrator/dag.py.
  - Parsing de argparse → cli.py (que sigue delgado; esto es lógica de
    negocio de materialización de inputs, no de parsing de flags).

Contrato de error
-----------------
Sigue el mismo estilo que `cli.py::_require_stages`: devuelve `None` si
todo salió bien, o un `str` con el mensaje de error a imprimir en stderr
(el caller — `cmd_run` — decide el `return 1`, esta función no aborta el
proceso ni lanza excepciones para casos de validación esperables).

Idempotencia
------------
No pisa un `input.xyz` (o `crest_conformers.xyz`) ya materializado si el
contenido coincide con el que se copiaría — evitar invalidar
`resume`/`state.json` en corridas repetidas de `cmd_run` sobre un workdir
que ya tiene progreso.

Pendiente (Problema 3, fuera de alcance de este módulo por ahora)
------------------------------------------------------------------
La rama `smiles` requiere una conversión SMILES→3D. `rdkit` todavía no es
dependencia del proyecto (ver pyproject.toml) — se importa acá de forma
perezosa y, si falta, se devuelve un mensaje de error claro en vez de un
`ImportError` crudo. Confirmar con kmilo si `rdkit` entra como dependencia
dura o como extra opcional (`[project.optional-dependencies] smiles`)
antes de cerrar esa rama.
"""

from __future__ import annotations

import filecmp
from pathlib import Path
from typing import Optional

from maaqc.config.schema import RecipeConfig
from maaqc.stages.crest import CrestStage
from maaqc.system._orca_shared import safe_copy

# Mismos nombres que ya usan CrestStage/OptStage — no reinventar convención.
_INPUT_XYZ_NAME = CrestStage.INPUT_FILENAME          # "input.xyz" (workdir/)
_CREST_STAGE_DIR = CrestStage.STAGE_DIR              # "01_crest"
_CREST_OUTPUT_NAME = CrestStage.OUTPUT_FILENAME       # "crest_conformers.xyz"


def _same_content(src: Path, dst: Path) -> bool:
    """True si `dst` ya existe y su contenido es idéntico al de `src`."""
    if not dst.exists():
        return False
    try:
        return filecmp.cmp(src, dst, shallow=False)
    except OSError:
        return False


def _materialize_from_path(src_arg: str, dest: Path) -> Optional[str]:
    """Copia un .xyz ya existente (`src_arg`) a `dest`, idempotente."""
    src = Path(src_arg)
    if not src.is_file():
        return f"no se encontró el archivo de entrada: {src}"
    if _same_content(src, dest):
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not safe_copy(src, dest):
        return f"no se pudo copiar {src} → {dest}"
    return None


def _materialize_from_smiles(smiles: str, dest: Path) -> Optional[str]:
    """Genera una geometría 3D desde un SMILES y la escribe en `dest`.

    Pendiente de confirmación (Problema 3): `rdkit` no es dependencia
    todavía. Import perezoso a propósito, para que las ramas xyz /
    crest-ensemble / sp-only funcionen sin rdkit instalado.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return (
            "entry='smiles' requiere el paquete 'rdkit' (todavía no es "
            "dependencia de maaqc — ver Problema 3 del issue de "
            "input_resolver). Instalalo manualmente por ahora, o usá "
            "--entry xyz con una geometría ya generada."
        )

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f"SMILES inválido, rdkit no pudo parsearlo: {smiles!r}"

    mol = Chem.AddHs(mol)
    embed_status = AllChem.EmbedMolecule(mol, randomSeed=0xC0FFEE)
    if embed_status != 0:
        return f"rdkit no pudo generar una conformación 3D para: {smiles!r}"
    AllChem.MMFFOptimizeMolecule(mol)

    xyz_block = Chem.MolToXYZBlock(mol)
    if _same_content_str(xyz_block, dest):
        return None
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(xyz_block)
    return None


def _same_content_str(text: str, dest: Path) -> bool:
    if not dest.exists():
        return False
    try:
        return dest.read_text() == text
    except OSError:
        return False


def materialize_entry_input(
    input_arg: str,
    workdir: Path,
    config: RecipeConfig,
    ensemble_arg: Optional[str] = None,
) -> Optional[str]:
    """Prepara los archivos que el primer stage de la cadena resuelta para
    `config.entry` espera encontrar en `workdir`.

    Devuelve `None` si todo quedó materializado correctamente, o un `str`
    con el mensaje de error a imprimir en stderr (mismo estilo que
    `cli.py::_require_stages` — no lanza excepciones para errores de
    validación esperables como archivos faltantes o flags requeridos).
    """
    entry = config.entry

    if entry == "xyz":
        return _materialize_from_path(input_arg, workdir / _INPUT_XYZ_NAME)

    if entry == "smiles":
        return _materialize_from_smiles(input_arg, workdir / _INPUT_XYZ_NAME)

    if entry == "crest-ensemble":
        if not ensemble_arg:
            return (
                "--entry crest-ensemble requiere --ensemble "
                "(ruta a un multi-XYZ de confórmeros ya generado)"
            )
        dest = workdir / _CREST_STAGE_DIR / _CREST_OUTPUT_NAME
        return _materialize_from_path(ensemble_arg, dest)

    if entry == "calibrated":
        if not config.reference:
            return "--entry calibrated requiere config.reference (ej. --reference shinorine)"
        # Mismo tratamiento que xyz/smiles según la forma del input.
        if input_arg.lower().endswith(".xyz"):
            return _materialize_from_path(input_arg, workdir / _INPUT_XYZ_NAME)
        return _materialize_from_smiles(input_arg, workdir / _INPUT_XYZ_NAME)

    if entry == "sp-only":
        # Standalone: Adc2Stage.inputs() valida directo desde 04_sp/, no
        # hay nada que materializar acá.
        return None

    return f"entry desconocido, no se sabe qué input materializar: {entry!r}"
