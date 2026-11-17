#!/bin/bash
# Preoptimizacion GFN2-xTB/ALPB(agua) de cada xyz de start/ y fixed_start/ (1 hilo c/u, 12 en paralelo).
# Uso: bash run_xtb_opt.sh   (~ unos minutos). Resultados en xtb_opt/<nombre>/xtbopt.xyz
cd "$(dirname "$0")"; mkdir -p xtb_opt
declare -A Q=( [A_fixH]=1 [D_regen]=0 [B_fixH]=-1 )
while IFS=$'\t' read -r f q n; do Q[${f%.xyz}]=$q; done < start/charges.tsv
job() { n=$1; q=$2; src=$(ls start/$n.xyz fixed_start/$n.xyz 2>/dev/null | head -1)
  mkdir -p xtb_opt/$n; cp $src xtb_opt/$n/in.xyz; cd xtb_opt/$n
  OMP_NUM_THREADS=1 xtb in.xyz --gfn 2 --alpb water --chrg $q --opt tight > xtb.out 2>&1; }
export -f job
for k in "${!Q[@]}"; do echo "$k ${Q[$k]}"; done | xargs -P 12 -L1 bash -c 'job $0 $1'
