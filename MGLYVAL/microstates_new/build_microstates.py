#!/usr/bin/env python3
"""Genera xyz de partida de D', R', C' y 010/100 (cribado) para MGLYVAL.

Solo usa numpy (no hay RDKit en este entorno). Parte de ../molecules/{A,B}.xyz
sin modificarlos. Indices 0-based; la numeracion 1-based del plan es +1:
  N8 -> idx 7 (Val), N9 -> idx 8 (Gly); H31 en N8, H46 en N9;
  Val-COOH: C20 (O3,O4), H47 en O3; Gly-COOH: C23 (O5,O6), H48 en O5.
Semillas E/Z: sXY, X=1 si se invierte la orientacion del sustituyente del N8 (idx 7),
Y=1 igual para N9 (idx 8). La inversion gira 180 grados el subarbol (Val o Gly)
alrededor del eje C(anillo)-N y, si el N lleva H, intercambia H y C.
"""
import numpy as np, os
HERE = os.path.dirname(os.path.abspath(__file__))
MOL = os.path.join(HERE, "..", "molecules")

def rd(f):
    L = open(f).read().split("\n"); n = int(L[0])
    a = [l.split() for l in L[2:2 + n]]
    return [x[0] for x in a], np.array([[float(v) for v in x[1:4]] for x in a])

def wr(f, s, x, comment):
    with open(f, "w") as fh:
        fh.write(f"{len(s)}\n{comment}\n")
        for e, p in zip(s, x):
            fh.write(f"{e:2s} {p[0]:16.10f} {p[1]:16.10f} {p[2]:16.10f}\n")

def bonds(s, x):
    r = {"H": .31, "C": .76, "N": .71, "O": .66}
    n = len(s); nb = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(x[i] - x[j]) < 1.25 * (r[s[i]] + r[s[j]]):
                nb[i].append(j); nb[j].append(i)
    return nb

def subtree(nb, start, block):
    seen = {block, start}; st = [start]
    while st:
        i = st.pop()
        for j in nb[i]:
            if j not in seen: seen.add(j); st.append(j)
    seen.discard(block); return sorted(seen)

def rot(x, idx, p0, ax, ang):
    ax = ax / np.linalg.norm(ax); c, s_ = np.cos(ang), np.sin(ang)
    v = x[idx] - p0
    x[idx] = p0 + v * c + np.cross(ax, v) * s_ + np.outer(v @ ax, ax) * (1 - c)

def fix_OH(x, o, h, d=0.97):
    u = x[h] - x[o]; x[h] = x[o] + d * u / np.linalg.norm(u)

def fix_NH(x, n, h, d=1.01):
    u = x[h] - x[n]; x[h] = x[n] + d * u / np.linalg.norm(u)

def flip(s, x, N, ring_c, side_c, hN=None):
    """Invierte orientacion del sustituyente (side_c) en el N indicado."""
    nb = bonds(s, x)
    sub = subtree(nb, side_c, N)
    old = x[side_c].copy()
    rot(x, sub, x[N], x[N] - x[ring_c], np.pi)
    if hN is not None:                      # amidinio: H va a la posicion antigua del C
        u = old - x[N]; x[hN] = x[N] + 1.01 * u / np.linalg.norm(u)
    return x

def drop(s, x, idxs):
    k = [i for i in range(len(s)) if i not in idxs]
    return [s[i] for i in k], x[k]

def dih(x, a, b, c, d):
    b0, b1, b2 = x[a] - x[b], x[c] - x[b], x[d] - x[c]
    b1 /= np.linalg.norm(b1)
    v = b0 - (b0 @ b1) * b1; w = b2 - (b2 @ b1) * b1
    return np.degrees(np.arctan2(np.cross(b1, v) @ w, v @ w))

# ---- bases saneadas (solo se corrige la distancia X-H; el resto lo cura xtb) ----
sA, A = rd(os.path.join(MOL, "A.xyz")); fix_OH(A, 3, 47); fix_OH(A, 5, 48)
sB, B = rd(os.path.join(MOL, "B.xyz")); fix_NH(B, 8, 46)
out = {}
# A_fixH (q=+1), D regenerado (q=0, Gly-COO-, Val-COOH), B_fixH (q=-1)
wr(f"{HERE}/fixed_start/A_fixH.xyz", sA, A, "A.xyz con O4-H=0.97 (q=+1)")
s, x = drop(sA, A, [48]); wr(f"{HERE}/fixed_start/D_regen.xyz", s, x, "D regenerado desde A: Gly-COO-, Val-COOH (q=0)")
wr(f"{HERE}/fixed_start/B_fixH.xyz", sB, B, "B.xyz con N9-H=1.01 (q=-1)")

# definiciones: (nombre, base, indices H a quitar, carga)
H_VALN, H_GLYN, H_VALO, H_GLYO = 31, 46, 47, 48   # N8-H, N9-H, Val-COOH, Gly-COOH (en A)
defs = [
    ("Dp",       "A", [H_VALO],                0, "D' (101): Val-COO-, Gly-COOH, amidinio"),
    ("Rp",       "A", [H_VALN],                0, "R': imina, H en N9 (N8 imina), 2xCOOH"),
    ("010_N8H",  "A", [H_GLYO, H_GLYN],       -1, "010: Gly-COO-, Val-COOH; N9 imina, H en N8"),
    ("010_N9H",  "A", [H_GLYO, H_VALN],       -1, "010: Gly-COO-, Val-COOH; N8 imina, H en N9"),
    ("100_N8H",  "A", [H_VALO, H_GLYN],       -1, "100: Gly-COOH, Val-COO-; N9 imina, H en N8"),
    ("100_N9H",  "A", [H_VALO, H_VALN],       -1, "100: Gly-COOH, Val-COO-; N8 imina, H en N9"),
]
for name, base, rem, q, desc in defs:
    for fx in (0, 1):
        for fy in (0, 1):
            s, x = list(sA), A.copy()
            hN7 = H_VALN if H_VALN not in rem else None
            hN8 = H_GLYN if H_GLYN not in rem else None
            if fx: x = flip(s, x, 7, 12, 16, hN7)
            if fy: x = flip(s, x, 8, 13, 21, hN8)
            s2, x2 = drop(s, x, rem)
            fn = f"{name}_s{fx}{fy}.xyz"
            wr(f"{HERE}/start/{fn}", s2, x2, f"{desc}; q={q}; semilla s{fx}{fy}")
            out[fn] = (q, len(s2))
# C' (q=-2): desde B, N9-H presente y N8 desprotonado
for fx in (0, 1):
    for fy in (0, 1):
        s, x = list(sB), B.copy()
        if fx: x = flip(s, x, 7, 12, 16, None)   # N8 sin H -> solo gira
        if fy: x = flip(s, x, 8, 13, 21, 46)
        # quitar H31 DESPUES de girar para que la lista siga siendo valida
        s2, x2 = drop(s, x, [H_VALN])
        fn = f"Cp_s{fx}{fy}.xyz"
        wr(f"{HERE}/start/{fn}", s2, x2, f"C': dianion, H en N9, N8 imina; q=-2; semilla s{fx}{fy}")
        out[fn] = (-2, len(s2))
with open(f"{HERE}/start/charges.tsv", "w") as fh:
    for k, (q, n) in sorted(out.items()): fh.write(f"{k}\t{q}\t{n}\n")
print(len(out), "xyz escritos")
