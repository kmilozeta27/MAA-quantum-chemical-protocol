#!/usr/bin/env bash
# watch_pipeline.sh — monitoreo automático del pipeline FAE-MAA
#
# Uso:
#   ./watch_pipeline.sh <mol_name>
#   ./watch_pipeline.sh mol_B          # ej: runs/mol_B/
#   ./watch_pipeline.sh mol_B --step opt   # forzar etapa específica
#
# Sin --step: detecta automáticamente la etapa activa leyendo qué
# archivos .out están siendo escritos en este momento (via lsof o mtime).

set -euo pipefail

# ── Colores ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GRN='\033[0;32m'; YEL='\033[1;33m'
CYN='\033[0;36m'; BLD='\033[1m';    RST='\033[0m'

# ── Args ──────────────────────────────────────────────────────────────────────
MOL_NAME="${1:-}"
FORCE_STEP="${3:-}"   # --step <nombre> → solo monitorea esa etapa

if [[ -z "$MOL_NAME" ]]; then
    echo -e "${RED}Uso: $0 <mol_name> [--step crest|opt|numfreq|sp|tddft|adc2]${RST}"
    echo -e "     Ejemplo: $0 mol_B"
    echo -e "              $0 mol_B --step numfreq"
    exit 1
fi

# Directorio raíz del pipeline (relativo al repo)
RUNS_DIR="${RUNS_DIR:-runs}"
MOL_DIR="$RUNS_DIR/$MOL_NAME"
PIPELINE_LOG="$MOL_DIR/pipeline.log"

if [[ ! -d "$MOL_DIR" ]]; then
    echo -e "${YEL}[watch] Esperando a que se cree $MOL_DIR ...${RST}"
    while [[ ! -d "$MOL_DIR" ]]; do sleep 2; done
fi

# ── Helpers ───────────────────────────────────────────────────────────────────

# Retorna el primer archivo .out activo (modificado en últimos 30 s) en $1/
find_active_out() {
    local dir="$1"
    # Archivos .out o .log modificados en los últimos 60 segundos
    find "$dir" -maxdepth 3 \( -name "*.out" -o -name "*.log" \) \
        -newer <(date -d '60 seconds ago' +%Y%m%d%H%M%S 2>/dev/null \
                 || date -v-60S +%Y%m%d%H%M%S 2>/dev/null \
                 || echo "19700101000000") \
        2>/dev/null | head -1
}

# Detecta qué etapa está corriendo ahora basándose en los directorios existentes
# y cuál tiene el .out más recientemente modificado.
detect_current_step() {
    local best_file="" best_time=0

    for stage_dir in \
        "$MOL_DIR/01_crest" \
        "$MOL_DIR/01b_censo_prescreen" \
        "$MOL_DIR/02_opt" \
        "$MOL_DIR/03_numfreq" \
        "$MOL_DIR/04_sp" \
        "$MOL_DIR/05_tddft" \
        "$MOL_DIR/06_adc2" \
        "$MOL_DIR/06_sp_gas"
    do
        [[ -d "$stage_dir" ]] || continue
        while IFS= read -r f; do
            mtime=$(stat -c %Y "$f" 2>/dev/null || stat -f %m "$f" 2>/dev/null || echo 0)
            if (( mtime > best_time )); then
                best_time=$mtime; best_file="$f"
            fi
        done < <(find "$stage_dir" -maxdepth 3 \( -name "*.out" -o -name "crest.log" \) 2>/dev/null)
    done

    echo "$best_file"
}

# Devuelve todos los .out de una etapa dada (para multi-conformer)
get_stage_outs() {
    local step="$1"
    case "$step" in
        crest)   find "$MOL_DIR/01_crest"         -name "crest.log"  2>/dev/null ;;
        censo)   find "$MOL_DIR/01b_censo_prescreen" \( -name "*.log" -o -name "*.out" \) 2>/dev/null ;;
        opt)     find "$MOL_DIR/02_opt"     -maxdepth 2 -name "opt.out"    | sort 2>/dev/null ;;
        numfreq) find "$MOL_DIR/03_numfreq" -maxdepth 2 -name "freq.out"   | sort 2>/dev/null ;;
        sp)      find "$MOL_DIR/04_sp"      -maxdepth 2 -name "sp.out"     | sort 2>/dev/null ;;
        tddft)   find "$MOL_DIR/05_tddft"   -maxdepth 2 -name "tddft.out"  | sort 2>/dev/null ;;
        adc2)    find "$MOL_DIR/06_adc2"    -maxdepth 2 -name "adc2.out"   | sort 2>/dev/null ;;
        sp_gas)  find "$MOL_DIR/06_sp_gas"  -maxdepth 2 -name "sp_gas.out" | sort 2>/dev/null ;;
        *)       echo "" ;;
    esac
}

# tail -f con encabezado coloreado; termina si el proceso ORCA muere
tail_watch() {
    local logfile="$1" label="$2"
    echo -e "\n${CYN}${BLD}━━━  $label  ━━━${RST}"
    echo -e "${YEL}  archivo: $logfile${RST}\n"

    if [[ ! -f "$logfile" ]]; then
        echo -e "${YEL}[watch] esperando a que aparezca $logfile ...${RST}"
        local waited=0
        while [[ ! -f "$logfile" ]]; do
            sleep 2; (( waited += 2 ))
            (( waited > 600 )) && { echo -e "${RED}[watch] timeout esperando $logfile${RST}"; return 1; }
        done
    fi

    # tail -f; se interrumpe con Ctrl-C o cuando watch_pipeline detecta
    # que la etapa terminó (el proceso principal cambia de archivo)
    tail -f "$logfile" &
    TAIL_PID=$!
    echo "$TAIL_PID"   # devuelve PID para que el caller pueda matarlo
}

kill_tail() {
    [[ -n "${TAIL_PID:-}" ]] && kill "$TAIL_PID" 2>/dev/null && wait "$TAIL_PID" 2>/dev/null || true
    TAIL_PID=""
}

# ── Banner inicial ─────────────────────────────────────────────────────────────
echo -e "${BLD}${GRN}"
echo "╔══════════════════════════════════════════════╗"
echo "║       watch_pipeline — FAE-MAA Monitor       ║"
echo "╚══════════════════════════════════════════════╝"
echo -e "${RST}"
echo -e "  Molécula : ${BLD}$MOL_NAME${RST}"
echo -e "  Dir      : ${BLD}$MOL_DIR${RST}"
echo -e "  Log      : ${BLD}$PIPELINE_LOG${RST}"
echo ""

# ── Modo forzado (--step) ─────────────────────────────────────────────────────
if [[ -n "$FORCE_STEP" ]]; then
    echo -e "${YEL}[watch] modo manual — etapa: $FORCE_STEP${RST}"
    outs=( $(get_stage_outs "$FORCE_STEP") )
    if [[ ${#outs[@]} -eq 0 ]]; then
        echo -e "${YEL}[watch] no hay archivos .out para '$FORCE_STEP' aún, esperando...${RST}"
        while [[ ${#outs[@]} -eq 0 ]]; do
            sleep 3
            outs=( $(get_stage_outs "$FORCE_STEP") )
        done
    fi
    # Si hay múltiples conformers, hacer tail de todos en paralelo
    if [[ ${#outs[@]} -eq 1 ]]; then
        tail -f "${outs[0]}"
    else
        echo -e "${CYN}[watch] ${#outs[@]} conformers — mostrando el activo más reciente${RST}"
        tail -f "${outs[@]}"
    fi
    exit 0
fi

# ── Modo automático ────────────────────────────────────────────────────────────
# Primero: mostrar el pipeline.log en background para ver estado global
echo -e "${CYN}[watch] Mostrando pipeline.log (log maestro) + archivo activo de cada etapa${RST}"
echo -e "${YEL}  Pulsa Ctrl-C para salir${RST}\n"

TAIL_PID=""
LAST_FILE=""
LAST_STEP=""

# Esperar a que el pipeline.log exista
if [[ ! -f "$PIPELINE_LOG" ]]; then
    echo -e "${YEL}[watch] esperando pipeline.log ...${RST}"
    while [[ ! -f "$PIPELINE_LOG" ]]; do sleep 2; done
fi

# Arrancar tail del pipeline.log
echo -e "\n${BLD}${GRN}━━━━━━━━━━━━━━━  PIPELINE LOG (maestro)  ━━━━━━━━━━━━━━━${RST}\n"
tail -f "$PIPELINE_LOG" &
PIPELINE_TAIL=$!

# Bucle principal: detectar etapa activa y hacer tail del .out correcto
while true; do
    active_file=$(detect_current_step)

    if [[ -n "$active_file" && "$active_file" != "$LAST_FILE" ]]; then
        # Etapa nueva detectada
        kill_tail

        # Extraer label legible del path
        step_label=$(echo "$active_file" | grep -oP '(0[1-6][^/]*)' | head -1 || echo "etapa")
        conf_label=$(echo "$active_file" | grep -oP 'conf_\d+' || echo "")
        full_label="$step_label${conf_label:+ / $conf_label}"

        echo -e "\n${BLD}${CYN}━━━━━━━━━━━━━━━  $full_label  ━━━━━━━━━━━━━━━${RST}"
        echo -e "${YEL}  → $active_file${RST}\n"

        tail -f "$active_file" &
        TAIL_PID=$!
        LAST_FILE="$active_file"
    fi

    sleep 8
done
