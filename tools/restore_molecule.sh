#!/usr/bin/env bash
# =============================================================================
#  restore_molecule.sh — vuelve a traer una carpeta de molécula al computador.
#
#  Uso:
#     tools/restore_molecule.sh <MOL>             # solo texto, desde GitHub
#     tools/restore_molecule.sh <MOL> <DESTINO>   # + binarios desde el disco externo
#
#  Ejemplo:
#     tools/restore_molecule.sh DG /media/kmilozeta/BACKUP1/PAPER1
# =============================================================================
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

[[ $# -ge 1 ]] || die "uso: $0 <MOL> [DESTINO]"
MOL="${1%/}"
SRC="${2:-}"

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"

git sparse-checkout add -- "$MOL"
echo "== $MOL restaurado desde git ($(git ls-files -- "$MOL" | wc -l) archivos)"

if [[ -n "$SRC" ]]; then
    [[ -d "$SRC/$MOL" ]] || die "no existe '$SRC/$MOL' (¿disco montado?)"
    # --ignore-existing: nunca pisa lo que viene de git
    rsync -a --ignore-existing --info=progress2 "$SRC/$MOL/" "$MOL/"
    echo "== binarios copiados desde $SRC/$MOL ($(du -sh "$MOL" | cut -f1) en total)"
fi
