#!/usr/bin/env python3
"""Analiza crest/*/crest_conformers.xyz: nº conformeros, ventana de energia, H sobre N8/N9/O3-O6 (transferencias de proton),
clase E/Z (dihedros C=N) y Boltzmann combinado de las dos semillas de cada microestado. Escribe crest_summary.tsv"""
import os, re, glob, collections, numpy as np
exec(open("build_microstates.py").read().split("# ---- bases")[0])
KT = 0.593  # kcal/mol a 298.15 K
H2K = 627.5095
def frames(path):
    L = open(path).read().split("\n"); i = 0; out = []
    while i < len(L) and L[i].strip():
        n = int(L[i]); E = float(L[i+1].split()[0]); out.append((E, [l.split() for l in L[i+2:i+2+n]])); i += n + 2
    return out
def feats(atoms):
    s = [a[0] for a in atoms]; x = np.array([[float(v) for v in a[1:4]] for a in atoms])
    nb = bonds(s, x)
    hs = lambda i: len([j for j in nb[i] if s[j] == "H"])
    H = "".join(str(hs(i)) for i in (7, 8, 3, 4, 5, 6))   # N8 N9 O3 O4 O5 O6
    rn = lambda c, n_: [j for j in nb[c] if j != n_ and s[j] == "C"][0]
    d7 = dih(x, 16, 7, 12, rn(12, 7)); d9 = dih(x, 21, 8, 13, rn(13, 8))
    cls = lambda d: "a" if abs(d) > 90 else "b"
    return H, cls(d7) + cls(d9)
rows = []; per = collections.defaultdict(list)
for d in sorted(glob.glob("crest/*/")):
    n = os.path.basename(d.rstrip("/"))
    if "interrumpida" in n or not os.path.exists(d + "crest_conformers.xyz"): continue
    fr = frames(d + "crest_conformers.xyz"); E0 = fr[0][0]
    ft = [feats(a) for _, a in fr]
    Hc = collections.Counter(h for h, _ in ft); Ec = collections.Counter(e for _, e in ft)
    w3 = sum(1 for E, _ in fr if (E - E0) * H2K <= 3)
    rows.append((n, len(fr), E0, w3, dict(Hc), dict(Ec)))
    micro = n.rsplit("_s", 1)[0]
    for (E, _), (h, e) in zip(fr, ft): per[micro].append((E, n, h, e))
with open("crest_summary.tsv", "w") as fh:
    for r in rows:
        print("%-12s N=%4d Emin=%.6f Eh  <=3kcal:%3d  Hsites(N8 N9 O3 O4 O5 O6)=%s  EZ=%s" % r); fh.write("\t".join(map(str, r)) + "\n")
print("\n== Microestados (dos semillas combinadas; xTB/ALPB, 298.15 K) ==")
for m, v in sorted(per.items()):
    Emin = min(e for e, *_ in v)
    w = np.array([np.exp(-(e - Emin) * H2K / KT) for e, *_ in v]); w /= w.sum()
    by = collections.defaultdict(float); byH = collections.defaultdict(float); byE = collections.defaultdict(float)
    for wi, (e, n, h, ez) in zip(w, v): by[n] += wi; byH[h] += wi; byE[ez] += wi
    best = min(v)[1]
    print(m, "Emin=%.6f (semilla %s)" % (Emin, best), "| peso por semilla:", {k: round(x, 3) for k, x in by.items()},
          "| H:", {k: round(x, 3) for k, x in byH.items()}, "| E/Z:", {k: round(x, 3) for k, x in byE.items()})
