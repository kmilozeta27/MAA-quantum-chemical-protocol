#!/usr/bin/env bash
# =============================================================================
#  archive_molecule.sh — respalda una carpeta de molécula en disco(s) externo(s),
#  sube a GitHub los resultados de texto y la saca del computador.
#
#  Uso:
#     tools/archive_molecule.sh <MOL> <DESTINO1> [DESTINO2 ...]
#
#  Ejemplo (doble backup):
#     tools/archive_molecule.sh DG /media/kmilozeta/BACKUP1/PAPER1 /media/kmilozeta/BACKUP2/PAPER1
#
#  Pasos:
#    1. rsync COMPLETO de <MOL>/ (binarios incluidos) a cada destino
#    2. verifica que la copia esté completa (tamaño+fecha; --checksum con VERIFY=checksum)
#    3. git add <MOL> (respeta .gitignore → solo texto), commit y push
#    4. comprueba que GitHub tenga exactamente lo mismo que el commit local
#    5. saca <MOL> del sparse-checkout (sigue en GitHub, ya no en el disco local)
#    6. pide confirmación y borra lo que queda (binarios ignorados) del disco local
#
#  Para traerla de vuelta: tools/restore_molecule.sh <MOL> [DESTINO]
# =============================================================================
set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

[[ $# -ge 2 ]] || die "uso: $0 <MOL> <DESTINO1> [DESTINO2 ...]"
MOL="${1%/}"; shift
DESTS=("$@")

ROOT="$(git rev-parse --show-toplevel)"
cd "$ROOT"
[[ -d "$MOL" ]] || die "no existe la carpeta '$MOL' en $ROOT"
[[ "$MOL" != */* && "$MOL" != .* ]] || die "MOL debe ser una carpeta de primer nivel (ej. DG)"

BRANCH="$(git rev-parse --abbrev-ref HEAD)"
SIZE_KB=$(du -sk "$MOL" | cut -f1)
echo "== $MOL: $(du -sh "$MOL" | cut -f1) en disco local, rama $BRANCH"

# --- 1-2. copiar y verificar en cada destino ---------------------------------
for D in "${DESTS[@]}"; do
    [[ -d "$D" ]] || die "destino '$D' no existe (¿disco montado?)"
    [[ -w "$D" ]] || die "sin permiso de escritura en '$D'"
    FREE_KB=$(df -Pk "$D" | awk 'NR==2 {print $4}')
    EXIST_KB=$( [[ -d "$D/$MOL" ]] && du -sk "$D/$MOL" | cut -f1 || echo 0 )
    (( FREE_KB + EXIST_KB > SIZE_KB )) || die "espacio insuficiente en '$D'"

    echo "-- copiando a $D/$MOL"
    rsync -a --info=progress2 "$MOL/" "$D/$MOL/"

    echo "-- verificando $D/$MOL"
    VFLAGS="-a -n -i"
    [[ "${VERIFY:-}" == "checksum" ]] && VFLAGS="$VFLAGS --checksum"
    DIFF=$(rsync $VFLAGS "$MOL/" "$D/$MOL/")
    [[ -z "$DIFF" ]] || { echo "$DIFF" | head; die "la copia en '$D' no coincide"; }
    echo "   OK ($(du -sh "$D/$MOL" | cut -f1))"
done

# --- 3-4. subir el registro de texto a GitHub ---------------------------------
git add -- "$MOL"
if ! git diff --cached --quiet -- "$MOL"; then
    git commit -q -m "results: archive $MOL (backup en disco externo)" -- "$MOL"
fi
git push -q origin "$BRANCH"
git fetch -q origin "$BRANCH"
[[ "$(git rev-parse HEAD)" == "$(git rev-parse "origin/$BRANCH")" ]] \
    || die "origin/$BRANCH no coincide con HEAD; no se borra nada"
[[ -z "$(git status --porcelain -- "$MOL")" ]] \
    || die "quedan cambios sin commitear en $MOL; no se borra nada"
echo "== GitHub al día: $(git ls-files -- "$MOL" | wc -l) archivos de $MOL en origin/$BRANCH"

# --- 5. sacar del sparse-checkout --------------------------------------------
if [[ "$(git config core.sparseCheckout || true)" == "true" ]]; then
    mapfile -t KEEP < <(git sparse-checkout list | grep -vx "$MOL" || true)
    git sparse-checkout set -- "${KEEP[@]}"
    echo "== $MOL fuera del sparse-checkout (sigue en GitHub)"
fi

# --- 6. borrar lo que quede (binarios ignorados) -----------------------------
if [[ -d "$MOL" ]]; then
    echo "Quedan $(du -sh "$MOL" | cut -f1) en ./$MOL (binarios ya respaldados en: ${DESTS[*]})."
    read -r -p "Escribe '$MOL' para borrar ./$MOL del computador: " ANS
    if [[ "$ANS" == "$MOL" ]]; then
        rm -rf -- "$MOL"
        echo "== ./$MOL borrado del disco local."
    else
        echo "No se borró nada."
    fi
fi
