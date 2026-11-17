#!/usr/bin/env python3
"""Resume xtb_opt/*: energia, distancias clave, H sobre O/N y dihedros C=N (E/Z). Escribe xtb_opt/summary.tsv"""
import os, glob, re, numpy as np
exec(open("build_microstates.py").read().split("# ---- bases")[0])
rows = []
for d in sorted(glob.glob("xtb_opt/*/")):
    n = os.path.basename(d.rstrip("/"))
    if not os.path.exists(d + "xtbopt.xyz"): rows.append((n, "FALLO")); continue
    s, x = rd(d + "xtbopt.xyz")
    E = float(re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+)", open(d + "xtb.out").read())[-1])
    nb = bonds(s, x)
    hs = lambda i: sorted(j for j in nb[i] if s[j] == "H")
    Hon = {k: len(hs(i)) for k, i in (("N8", 7), ("N9", 8), ("O3", 3), ("O4", 4), ("O5", 5), ("O6", 6))}
    dd = lambda i, j: np.linalg.norm(x[i] - x[j])
    # dihedro C16-N7-C12-C(ring vecino) y C21-N8-C13-C(ring)
    rn = lambda c, n_: [j for j in nb[c] if j != n_ and s[j] == "C"][0]
    d7 = dih(x, 16, 7, 12, rn(12, 7)); d9 = dih(x, 21, 8, 13, rn(13, 8))
    rows.append((n, E, "Hsobre=" + " ".join(f"{k}:{v}" for k, v in Hon.items()),
                 f"C12-N7={dd(12,7):.2f} C13-N9={dd(13,8):.2f}",
                 f"C20-O3/O4={dd(20,3):.2f}/{dd(20,4):.2f} C23-O5/O6={dd(23,5):.2f}/{dd(23,6):.2f}",
                 f"dih7={d7:.0f} dih9={d9:.0f}"))
with open("xtb_opt/summary.tsv", "w") as fh:
    for r in rows: fh.write("\t".join(map(str, r)) + "\n"); print("\t".join(map(str, r)))
