# MGLYVAL/microstates_new (P1-13)

Microestados nuevos del plan (`docs/CALIBRATION_PLAN.md`), generados con numpy a partir de `../molecules/{A,B}.xyz`
(sin modificarlos; RDKit no esta instalado, si lo esta xtb 6.7.1, crest 3.0.2, obabel, ase).
Numeracion 1-based del plan: N8 = Val (idx 7), N9 = Gly (idx 8). Val-COOH = C21 (O4,O5), Gly-COOH = C24 (O6,O7).

| Carpeta/archivo | Contenido |
|---|---|
| `build_microstates.py` | genera todo; reproducible |
| `fixed_start/` | `A_fixH` (O4-H 0.87 -> 0.97), `B_fixH` (N9-H 1.17 -> 1.01), `D_regen` (desde A, Gly-COO-/Val-COOH). Solo reparan X-H; las longitudes C-O y C-N de A/B/D siguen mal hasta optimizar |
| `start/` | 28 xyz: `Dp` (D', q=0), `Rp` (R', q=0), `Cp` (C', q=-2), `010_N8H/N9H`, `100_N8H/N9H` (q=-1; el sufijo indica el N que conserva H). Notacion 010 = Gly,Val,N: 010 = Gly-COO-, Val-COOH; 100 = Gly-COOH, Val-COO- |
| `start/charges.tsv` | carga y n.o de atomos de cada xyz |
| `xtb_opt/` | `run_xtb_opt.sh`: GFN2/ALPB(agua) --opt tight de los 31 xyz (32 s); `summary.tsv` via `analyze_xtb.py` |
| `run_crest.sh` | CREST listo, NO lanzado (ver estimacion dentro) |

Semillas `sXY`: X=1 invierte la orientacion del sustituyente (Val) en N8, Y=1 el de Gly en N9 (giro 180 grados del subarbol
alrededor del eje C-N; si el N lleva H, se intercambian H y C). Los dihedros C(sust)-N-C-C(anillo) en `summary.tsv`
confirman que cada semilla se mantuvo tras la optimizacion (~+160/170 grados frente a ~-30/-5 grados para N8; ~-10/-20 frente a ~+160/-180 para N9).
Las 28 optimizaciones terminaron con el tautomero/protonacion pretendidos (comprobado por conectividad).

Pendiente: CREST, DFT (P1-14), asignar nombre E/Z CIP formal a las semillas, y restaurar/optimizar A/R/D/C.
