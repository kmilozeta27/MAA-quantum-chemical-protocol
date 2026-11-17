# TODO — PAPER1 (maaqc / MAAs, Marine Drugs)

> Archivo vivo. Se actualiza al terminar o cambiar cada tarea: marcar `[x]`, anotar fecha y resultado en el **Registro**, y mover lo bloqueado a **Decisiones abiertas**.
> Origen: inventario local (agente `reader`/Haiku) + plan verificado (agente `planner`/Opus), 2026-11-17.
> Leyenda: `[ ]` pendiente · `[~]` en curso · `[x]` hecho · `[!]` bloqueado. Carril A = pKa (8 hilos, ≤10 GB). Carril B = gold standard CC (4 hilos, ~18 GB).

## Prioridad (recomendación)

**1º calibración de pKa; el gold standard corre en paralelo como carril secundario.**

- La calibración es la tesis del paper. La clase que decide la especiación a pH 8.1 (B→C de MGLYVAL, `amino_dianion`) tiene n=1 y error crudo +6.0 pKa (v3, SMD/direct). Además está mal emparejada: B→C es la desprotonación de un amidinio vinílogo, no de un NH₃⁺ sp³.
- El gold standard está más avanzado de lo que dice el README: ya hay STEOM-DLPNO-CCSD sobre la molécula completa (SVP/SVPD/ma-SVP/CPCM) = 3.41–3.62 eV frente a 4.20 eV de ωB97X-D3 (exp. 3.70 eV / 335 nm). Falta solo el incremento TZVP (caro).
- No hay dependencia dura entre ambos: B (−1) domina a pH 7 con cualquier calibrador, así que el benchmark sobre B no espera a la calibración.
- Nunca correr a la vez un NumFreq de 12 núcleos / 25.7 GB y el carril B. Un solo job MDCI a la vez (los `.tmp` llenan el disco).

## Estado real (verificado, corrige el README)

| Sistema | Estado real |
|---|---|
| AP (APA/APR) | 6 etapas OK, 1 conformero cada uno. Único calibrador `pyridinium`. |
| DG 4DG_A/T | OPT, FREQ, SP, TDDFT. Sin SP_GAS/SP_SMD. |
| DG 4DG_B | 7 OPT, 3/4 NumFreq, sin SP ni TDDFT. |
| MGLYVAL B | CREST→CENSO (33 sobrevivientes)→6 conf: OPT/NumFreq/SP 6/6. conf_004 (−8.7 cm⁻¹) y conf_005 (−17.6 cm⁻¹) con imaginarias. TDDFT rehecho a mano en `conf_5`/`conf_6`. |
| MGLYVAL A/R/D/C | SP CPCM según `pipeline.log`; `pka` nunca corrió con éxito en este layout. |
| 07_adc2 | STEOM SVP/SVPD/ma-SVP/CPCM terminados; ADC(2) cortado por timeout (f=0.037, no es el estado brillante); maTZVP y bt-PNO-EOM (gold_01) fallaron. |

Falsedades del README a corregir: "every conformer converged", "TD-DFT (10 compounds) ✅ Done", "≈8 unidades" (es 6 con SMD/direct; 8 era born/CPCM), "complete run record" de DG.

---

## Fase 0 — Higiene y bugs (1–2 días, sin ORCA)

- [~] **P0-1** (2026-11-17: datos localizados en el disco TOSHIBA `PAPER1/` y `PAPER1.2/`, ver docs/ARCHIVO_DISCO.md; falta restaurar a una ruta aparte con `--base`, sin tocar el sparse-checkout) Restaurar outputs de calibradores y MGLYVAL A/R/D/C/B desde discos externos (`tools/restore_molecule.sh`, ver con `git show HEAD:tools/restore_molecule.sh`). **No tocar el sparse-checkout** (ver memoria). Hecho: `python maaqc.py calibrate --base <ruta> --spsmd --method direct --gcorr goodvibes` reproduce `calib_v3.csv` con |Δ|<0.01; si no, documentar diferencia v2/v3.
- [x] **P0-2** (2026-11-17, ver Registro) Bug de datos: `collect` omite `conf_5`/`conf_6` (globs en `maaqc.py:~10146` y `~10395`); `B_results.xlsx` pierde el 39 % del peso de Boltzmann (incl. el conformero dominante, 36 %). Resolver id de conformero numéricamente. Hecho: hoja TDDFT con 6 conf y pesos suman 1.0. *Verificado solo con fixture sintético (`test_collect_tddft_includes_all_conformer_dir_names`); falta regenerar `B_results.xlsx` real cuando 04_sp/02_opt/03_numfreq de B estén en disco (P0-1).*
- [x] **P0-3** (2026-11-17) CENSO no encontrado → hoy sigue con CREST sin filtrar (`maaqc.py:~3052/3074/3130`). Convertir en error fatal salvo flag explícito. *Hecho: `CensoUnavailableError`; flag `--allow-no-censo` (comando `run`).*
- [ ] **P0-4** Imaginarias conf_004/conf_005: reoptimizar con TightOpt + desplazamiento por el modo (~2×5 h NumFreq) o documentar umbral (|ν|<20 cm⁻¹ con qRRHO). Recomendado: reoptimizar.
- [x] **P0-5** (2026-11-17) Sacar `EMPIRICAL_CALIBRATION` obsoleto (`maaqc.py:~11232`) del camino por defecto; exigir `--calibration-file` o cargar `calibration_v3.json`. Marcar `calibration/calibration.csv` como legacy. *Hecho: `--empirical-correction` exige `--calibration-file` o `calibration/calibration_v3.json` (si no, error); escape `--allow-builtin-calibration`; legacy anotado en docs/METHODS.md (el CSV no se tocó). La tabla embebida sigue en el código, solo fuera del camino por defecto.*
- [x] **P0-6** (2026-11-17) Corregir README y METHODS (lista de falsedades arriba; gold standard ya tiene STEOM; tabla real por clase). *Hecho solo README.md y docs/METHODS.md; el manuscrito .docx no se tocó. La tabla por clase usa los errores crudos de `calibration_v3.json` (leído con `git show`).*
- [x] **P0-7** (2026-11-17) Tests nuevos (fixtures sintéticos en `tmp_path`; archivos reales con `needs(...)` para skip limpio en sparse):
  - [x] `test_collect_tddft_includes_all_conformer_dir_names`
  - [x] `test_extract_steom_vs_cis_spectrum` (CCSD_CPCM.out → 3.549 eV; CCSD_maTZVP.out → None)
  - [x] `test_adc2_requires_normal_termination`
  - [x] `test_fit_empirical_offset_n1_rmse_nan`, `test_fit_auto_switches_to_linear_n3_span1`, `test_loo_rmse_two_points`
  - [x] `test_refuses_cross_class_without_force`, `test_calibration_file_overrides_builtin`
  - [x] `test_speciation_fractions_sum_to_one_and_monotonic`
  - [x] `test_mc_shared_proton_draw_reduces_dpKa_variance`, `test_mc_reproducible_with_seed`
  - [x] `test_boltzmann_weights_multiT_sum_one`, `test_censo_missing_aborts`
  - [x] `test_nto_selection_dual_criterion`, `test_parse_multiwfn_esp_log`
  - Hecho: `python -m pytest test_maaqc.py -q` → **54 passed, 0 failed, 0 skipped** (2026-11-17). Pendiente: confirmar CI `.github/workflows/tests.yml` (no verificado en esta corrida). Cambio de comportamiento: `extract_adc2_lambda` / `parse_adc2_excitations` ahora devuelven None / [] si el `.out` no terminó normalmente (el `adc2.out` real de conf_006 está cortado por timeout).
- [~] **P0-8** (2026-11-17, solo diagnóstico, sin cambios de código) Confirmar que `_setup_run_resources()` respeta límites de hilos/RAM (NumFreq por defecto 12 núcleos / 25.7 GB, `maaqc.py:528-531`); ajustar para moléculas <20 átomos. *Diagnóstico: sí hay flags (`--cores N`, `--nprocs`, `--maxcore`, `--nprocs-numfreq`, `--maxcore-numfreq`); `ResourceProfile` (maaqc.py ~1004-1250) reparte: CREST/CENSO min(c,12); OPT 3 jobs×4; NumFreq 1 job×min(c,12) con maxcore=RAM·0.85/nprocs; SP 3×4; TDDFT 2×6 con maxcore/1.4; ADC2 ≤4 cores, ≥5000 MB/core. RAM = total del sistema (psutil), `--cores` solo limita cores, no RAM; no hay flag de RAM total ni escalado por número de átomos → ajuste para <20 átomos pendiente (decisión: usar --nprocs-numfreq/--maxcore-numfreq manual).*

## Fase 1 — Prioridad A: calibración de pKa (carril A)

Diseño completo, tablas de calibradores y criterios de aceptación: [docs/CALIBRATION_PLAN.md](docs/CALIBRATION_PLAN.md). Sitios a calibrar: COOH de Gly (i), COOH de Val (ii) y N–H del amidinio vinílogo (iii). Macro-pKa: M1 (+1→0, 1.er COOH), M2 (0→−1, 2.º COOH), M3 (−1→−2, N–H con 2×COO⁻).

**Sin ORCA (hacer primero):**
- [~] **P1-9** (2026-11-17: texto en docs/METHODS.md §11, en español; falta crear el tag al commitear, traducir o aceptar la mezcla de idioma y reconciliar §6 con la regla de pendiente congelada del plan: n≥5, rango ≥3, ΔAICc>2) Congelar especificación (SMD primario, qh-GoodVibes, offset por clase, criterios de aceptación) en METHODS + tag git `pka-spec-v1`. Hecho: tag existe y METHODS lo cita.
- [x] **P1-10** (2026-11-17, docs/CALIBRATION_LOO.md: carboxyl RMSE_LOO = 0.198; repetir tras el cambio de Asp de P1-11) Recalcular offsets y LOO por clase desde `git show HEAD:calibration/calib_v3.csv`; tabla error crudo vs carga del producto. Hecho: carboxyl RMSE_LOO ≈ 0.20.
- [~] **P1-11** (2026-11-17: constante `PKA_EXP` en maaqc.py, valores verificados contra Williams/Evans [W]; Asp ahora 2.09/3.86/9.82; falta contraste con IUPAC/NIST SRD 46 y verificar los pKa de los calibradores nuevos) Unificar pKa experimentales en una fuente (IUPAC/NIST SRD 46); corregir Asp en `CALIBRATION_SETS` (~l.16908) y `EMPIRICAL_CALIBRATION` (~l.11232) (código usa 1.88/3.65/9.60, [W] da 2.09/3.86/9.82).
- [~] **P1-12** (2026-11-17: hecho en maaqc.py con 60 tests passed; sin cablear a la CLI (`pka`/`calibrate`); faltan bootstrap conformacional, incertidumbre de pendiente en MC lineal y fracciones f_R/f_D con N microestados) Código: N microestados por carga con degeneración (`_build_neutral_macrostate` / `compute_tautomer_equilibrium` aceptan solo 2 tautómeros), LOO en modo lineal, criterio AICc, MC con error de calibración (hoy no lo hace), tipos nuevos en `EMPIRICAL_TYPE_META`. Tests: `test_loo_linear`, `test_macrostate_n_micro_degeneracy`, `test_mc_calibrated_no_hplus`, `test_double_difference_class`.
- [~] **P1-13** (2026-11-17: 28 xyz con semillas + xtb preoptimizado en MGLYVAL/microstates_new/, CREST listo en run_crest.sh pero sin lanzar; A.xyz, D.xyz y B.xyz arrancaban distorsionados, ver README de esa carpeta) Generar xyz y CREST (con semillas E/Z) de D′ (Val-COO⁻/Gly-COOH), R′ (imina N9–H), C′ (dianión N9–H) y 010/100 (cribado). Revisar además que A.xyz, D.xyz y B.xyz no tengan geometrías de partida mal armadas (O–H 0.87 Å, N9–H 1.17 Å).

**Carril A, con ORCA (8 hilos, ≤10 GB):**
- [ ] **P1-14** DFT de D′ y C′ (prioridad), luego R′; verificar A/R/D/C restaurados (P0-1). ~2.5 d por microestado.
- [ ] **P1-15** Clase dianión: Asp-AC + IDA + MeIDA (reemplaza P1-1; n=4). 2.5–3.5 d.
- [ ] **P1-16** Amidinio +1→0: 2-AP, benzamidina, acetamidina (+4-DMAP opcional) con AP (refina P1-2). 1–1.5 d.
- [ ] **P1-17** Amonio +1→0: metilamina, etanolamina, Tris, éster de glicina (para doble diferencia). 1–1.5 d.
- [ ] **P1-18** carboxyl: creatina, betaína, β-Ala, GABA. 1.5–2 d.
- [ ] **P1-19** Test de aditividad: nicotínico, isonicotínico, piridina (reemplaza P1-3; arginina fuera del ajuste). 1–1.5 d.
- [ ] **P1-20** Validación externa: sarcosina, N,N-dimetilglicina, glicilglicina, fenilguanidina (no entran al ajuste). 1.5 d.
- [ ] **P1-4′** `pka` de MGLYVAL completo con macroestados de N microestados y MC. Dep: P1-9, P1-12, P1-14. Hecho: M1/M2/M3 crudos y calibrados con IC; fracciones a pH 7.0 y 8.1.
- [ ] **P1-5′** Sensibilidad de M3: dianión solo / doble diferencia / amidinio +1→0 / crudo; fracción de C a pH 8.1.
- [ ] **P1-6** Validación global (MAE, RMSE, LOO por clase y global; CPCM vs SMD; armónico vs qh). Tabla SI vía `calibrate --csv`.
- [ ] **P1-21** (opcional, alto valor) Porphyra-334 como validación externa MAA contra pKa 2.92 ± 0.11 de RMN (Chuang 2014). 5–8 d.
- [ ] **P1-7** DG: completar NumFreq de 4DG_B y SP; SP_GAS/SP_SMD para A/B/T. Prioridad menor.
- [ ] **P1-8** (opcional, caro) Cluster-continuo con aguas explícitas en M3. Solo si sobra tiempo.

Orden: P1-9 → P1-12 → P1-13 → P1-15 → P1-16/17 → P1-14 → P1-18 → P1-19 → P1-4′ → P1-5′ → P1-20 → P1-6. Carril A sin porphyra-334 ≈ 3 semanas de cómputo. Los NumFreq de MGLYVAL (47–49 át., 12 núcleos / 25.7 GB por defecto) no pueden coincidir con un MDCI del carril B.

## Fase 2 — Prioridad B: gold standard de excitaciones (carril B)

- [ ] **P2-1** Fijar modelo de cromóforo truncado: catión (+1), amidinio con COO⁻→CH₃/H; generar `cromoforo_31at.xyz` desde conf_6; resolver la inconsistencia gold_01/03 (q=0, `cromoforo_S1.xyz`) vs gold_02 (q=+1, archivo inexistente). Hecho: STEOM-DLPNO SVP gas truncado vs completo (3.624 eV) |Δ|≤0.1 eV y NTO equivalente. ~1–3 h.
- [ ] **P2-2** Corregir gold_01: `%moinp` lee un `.gbw` con el mismo nombre que el job (falla en GUESS); usar `guess_SVP.gbw` o quitar MOREAD. Correr bt-PNO-EOM-CCSD y STEOM-DLPNO def2-TZVP en el truncado. Hecho: terminación normal, NRoots≥3, carácter activo ≥96 %. 1–2 días, `%maxcore 4500`×4, scratch en SSD con 50–100 GB libres.
- [ ] **P2-3** Corrección de base aditiva (focal point): E(completo, SVPD, CPCM) + [E(trunc, TZVP(D)) − E(trunc, SVPD)]; falta STEOM SVPD+CPCM en molécula completa. Hecho: S1 de referencia con incertidumbre.
- [ ] **P2-4** Benchmark de funcionales sobre conf_6, TZVPD CPCM, 15 raíces: B3LYP, CAM-B3LYP, M06-2X, ωB97M-V, ωB2PLYP o SCS-ωPBEPP86 (PBE0 y ωB97X-D3 ya hechos). ~1–2 h c/u, 6 hilos; intercalar en carril A. Hecho: tabla ΔE(funcional − STEOM) y Δ vs exp.
- [ ] **P2-5** Validar en segunda geometría: STEOM SVP en conf_004 (23.6 %) y, si P1-5 lo pide, en C. Hecho: dispersión conformacional del error TD-DFT <0.05 eV.
- [ ] **P2-6** ADC(2): decidir (recomendado: retirar del paper; STEOM ya lo cubre). Si se conserva: NRoots≥3 + CPCM sobre el truncado.
- [ ] **P2-7** Comparación con experimento: ensemble de Boltzmann de 6 conf (tras P0-2) + corrección vertical→máximo de banda (justificar 0.1–0.2 eV). Hecho: λmax crudo y corregido por offset STEOM vs 335 nm.
- [ ] **P2-8** Aclarar fallo de `CCSD_maTZVP` (¿disco/scratch lleno en MDCI `.tmp`?) y limpiar `sel_nr*` sin `.out` de `07_adc2/conf_6/`.

## Fase 3 — Manuscrito y cierre

- [ ] **P3-1** Tablas SI finales: calibradores (versión final), sensibilidad B→C, benchmark de funcionales, focal point STEOM.
- [ ] **P3-2** Extender a otras MAAs (shinorine, porphyra-334, …) solo con especificación congelada; si no hay tiempo, acotar el paper a MGLYVAL + 4DG.
- [ ] **P3-3** Comparación con estructuras de docking (10 single points, ~1 día).
- [~] **P3-5** (2026-11-17) Manuscrito actualizado según carpeta local → `~/Descargas/MAAqc manuscript draft v2.docx` (original intacto). Editados: 2.2, 2.3 (microestados, Tabla 4, extensión planificada, sesgo B→C), 2.5 (datos STEOM preliminares), 3.1, 4.5–4.7, sección 5, tabla "Current status". Pendiente: resultados numéricos finales, Tabla 3/4 con calibradores nuevos, referencias (Chuang 2014, Williams–Evans), abstract y Tabla 1 ("nine MAAs", λmax) no verificables con la carpeta local.
- [~] **P3-6** (2026-11-17) Respaldo final en disco EC IEII: `/run/media/kmilozeta/EC IEII/KMILO/MAAqc_BACKUP_FINAL/` (01_calculos, 02_repositorio_y_scripts con bundle git, 03_papers_MAAs, 04_manuscrito, 90_duplicados_para_revisar). Hecho: fusión TOSHIBA+portátil sin duplicados, bundle git verificado, 0 pérdidas detectadas por tamaño/fecha. Pendiente: verificación checksum completa (en curso), segunda copia en otro disco, borrar duplicados de 90_, reparar TOSHIBA con chkdsk, y recién entonces limpiar el portátil (no se borró nada local).
- [x] **P3-7** (2026-11-17) Papers de Zotero copiados y ordenados en `docs/papers/` (164 PDFs, índice, BibTeX); falta revisar clasificación y descargar los PDFs que Zotero registra sin archivo (Dreuw 2015, Nishida 2021, etc.).
- [ ] **P3-4** Reproducibilidad: archivar `.inp`/outputs de texto de calibradores en git, `--mc-seed` fijo, etiqueta git de la especificación congelada.

---

## Orden de ejecución sugerido

1. Fase 0 completa (P0-1 a P0-3 primero: desbloquean todo).
2. Carril A: ver orden en la Fase 1 (P1-9 → P1-12 → P1-13 → …).
3. Carril B en paralelo: P2-1 → P2-2 → P2-3 → P2-5; P2-4 en los huecos del carril A.
4. P2-7, luego Fase 3.

## Decisiones abiertas

- [ ] **D1′** Valor primario de M3 (B→C): doble diferencia si pasa el test de aditividad; si no, clase dianión sola; las demás como sensibilidad. "B domina a pH 8.1" solo si P(pKa M3 > 9.1) ≥ 0.95 en todas las variantes.
- [ ] **D5** SMD como primario (recomendado) o CPCM. **D6** Incluir porphyra-334 (1–2 semanas). **D7** Excluir guanidinios del ajuste (recomendado). **D8** Fuente única de pKa: IUPAC/NIST SRD 46.
- [ ] **D2** Carga del cromóforo truncado: recomendado catión (+1), validado con P2-1.
- [ ] **D3** Imaginarias conf_004/005: reoptimizar (recomendado) o aceptar con qRRHO.
- [ ] **D4** ADC(2): retirar del paper (recomendado).

## Riesgos

- El offset de calibración absorbe errores sistemáticos de solvatación; con n<3 la pendiente no es estimable → reportar LOO, no solo RMSE.
- Solo existe un pKa experimental de una MAA (COOH de porphyra-334, 2.92 ± 0.11, RMN en D₂O; Chuang 2014); ninguno para el N–H amidinio: validación indirecta para M3.
- λmax experimental sin pH ni especie rastreados (aviso del README).
- STEOM-DLPNO: error típico ±0.1–0.2 eV; el truncado puede perder la polarización de los COO⁻.
- Costo/disco del TZVP en MDCI.
- Entre v2 y v3 de la calibración (mismo día) cambian Gly y Val sin explicación → reproducibilidad.

## No verificado (pendiente de comprobar)

- Outputs de calibradores y de MGLYVAL B 01–04c y A/R/D/C (fuera de disco y de git); si los SP_GAS/SP_SMD de A/R/D/C terminaron normal.
- ~~Si `pytest` pasa hoy~~ → 54 passed (2026-11-17).
- "Percentage Active Character" de las corridas STEOM existentes.
- Causa exacta del fallo de CCSD_maTZVP.
- ~~Flags de hilos/RAM en `run`~~ → existen, ver P0-8.
- pKa de literatura de los calibradores nuevos propuestos.
- pH y solvente de los λmax del README.
- Contenido de `tools/restore_molecule.sh`.

## Registro de cambios

- 2026-11-17 — Archivo creado a partir del inventario y plan de los agentes `reader` y `planner`. Ninguna tarea ejecutada aún.
- 2026-11-17 — Fase 0 parcial (P0-2, P0-3, P0-5, P0-6 docs, P0-7, diagnóstico P0-8). `maaqc.py`: ids numéricos de conformero (`_conf_num`, `_ensemble_conf_dirs`, `_find_conf_dir`, `_canon_conf_id`) en `collect_results`, ESP y ADC2; `CensoUnavailableError` + `--allow-no-censo`; `resolve_calibration_source` + `--allow-builtin-calibration`; extractores ADC2 exigen terminación normal; el dispatcher `main()` ahora imprime los mensajes de `SystemExit(str)` (antes se perdían, rc=1 mudo). `test_maaqc.py`: +19 tests, `pytest` = 54 passed. README/METHODS corregidos. Sin ORCA, sin commit. P0-1, P0-4 sin tocar.
- 2026-11-17 — Fase 1 reescrita con el diseño de calibración de tres sitios (docs/CALIBRATION_PLAN.md): P1-1/2/3 reemplazados por P1-15/16/19; añadidos P1-9 a P1-21. Manuscrito pendiente de actualizar (P3-5).
- 2026-11-17 — Manuscrito v2 generado (P3-5 en curso); docs/CALIBRATION_PLAN.md creado.
- 2026-11-17 — Inventario del disco TOSHIBA y plan de archivo (docs/ARCHIVO_DISCO.md); snapshot del portátil en MAAqc_ARCHIVO verificado; papers de Zotero en docs/papers.
- 2026-11-17 — Respaldo final consolidado en EC IEII/KMILO/MAAqc_BACKUP_FINAL; papers de filosofía retirados del repo (copiados a 05_otros_no_proyecto).
- 2026-11-17 — P1-9 a P1-13 con 4 agentes en paralelo. P1-10 hecha; P1-9, 11, 12 y 13 parciales (ver notas en cada tarea). `maaqc.py`: `PKA_EXP`, `compute_micro_equilibrium`, AICc, LOO lineal, MC calibrado, `apply_double_difference_correction`; pytest = 60 passed. Sin commit ni tag.
