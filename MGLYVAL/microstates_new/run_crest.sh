#!/bin/bash
# CREST (iMTD-GC, GFN2/ALPB agua) sobre las semillas preoptimizadas. NO LANZADO todavia.
# Mismas opciones que maaqc.py (--T N --chrg q --uhf 0 --alpb water; topologia de referencia por defecto,
# asi CREST no cambia tautomero y las semillas s00/s01/s10/s11 conservan su isomeria C=N / orientacion N-R).
# Seleccion por defecto: las 2 semillas de menor energia xTB de cada microestado (2 x 7 = 14 corridas).
# Con ALL=1 corre las 28. Estimacion: 48 atomos, 12 hilos, ~1.5-4 h por corrida (14 corridas ~1-2 dias en serie).
# Uso: nohup bash run_crest.sh > crest_all.log 2>&1 &
cd "$(dirname "$0")"; NT=${NT:-12}
python3 - <<'PY' > crest_jobs.txt
import re, collections
Q = {l.split("\t")[0][:-4]: l.split("\t")[1] for l in open("start/charges.tsv")}
rows = [l.split("\t") for l in open("xtb_opt/summary.tsv")]
g = collections.defaultdict(list)
for r in rows:
    n = r[0]
    if n in Q: g[n.rsplit("_s", 1)[0]].append((float(r[1]), n))
import os
k = 99 if os.environ.get("ALL") else 2
for m, v in sorted(g.items()):
    for e, n in sorted(v)[:k]: print(n, Q[n])
PY
while read -r n q; do
  # reanudable: salta corridas terminadas; una corrida incompleta se aparta (no se borra) y se repite
  if grep -qs "CREST terminated normally" crest/$n/crest.log; then echo "skip $n (terminada)"; continue; fi
  if [ -d crest/$n ]; then mv crest/$n crest/${n}.interrumpida_$(date +%Y%m%d_%H%M%S); fi
  mkdir -p crest/$n; cp xtb_opt/$n/xtbopt.xyz crest/$n/in.xyz
  ( cd crest/$n && crest in.xyz --T $NT --chrg $q --uhf 0 --alpb water > crest.log 2>&1 )
done < crest_jobs.txt
