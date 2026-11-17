# Plan de calibración de pKa para MGLYVAL (tres sitios)

Diseño del 2026-11-17 (agente `planner`, revisado). La carpeta local es la fuente de verdad. Las tareas con ID `P1-x` viven en [../TODO.md](../TODO.md); este documento guarda el razonamiento. Los pKa experimentales marcados **NV** no están verificados en fuente primaria.

## 1. Sitios y microestados

Sitios a calibrar: **(i)** COOH de Gly, **(ii)** COOH de Val, **(iii)** N–H del amidinio vinílogo que une el ciclo con los sustituyentes. Asignación por conectividad de `MGLYVAL/molecules/*.xyz`: el carboxilo C21 es el de Val (unido a N8) y C24 es el de Gly (unido a N9).

| xyz | q | Gly-COOH | Val-COOH | N8–H | N9–H | Especie |
|---|---|---|---|---|---|---|
| A | +1 | H | H | sí | sí | amidinio, 2×COOH |
| R | 0 | H | H | sí | no (C=N9 1.29 Å) | imina neutra |
| D | 0 | COO⁻ | H | sí | sí | zwitterión |
| B | −1 | COO⁻ | COO⁻ | sí | sí | amidinio dicarboxilato |
| C | −2 | COO⁻ | COO⁻ | sí | no | imina dianión (tautómero N8–H) |

La imagen del usuario coincide con A, R, D y B.

**Microestados que faltan** (notación Gly, Val, N con 1 = protonado):

- D′ (101): Val-COO⁻ / Gly-COOH, amidinio, q=0. **Obligatorio.**
- R′: imina con H en N9, q=0. Necesario para cerrar el macroestado neutro.
- C′: dianión con H en N9, q=−2. **Obligatorio.**
- 010 y 100: imina + un COOH, q=−1, dos tautómeros cada uno. Solo cribado CREST/CENSO; pasan a DFT si quedan a ≤3 kcal/mol del mínimo de su carga.
- Isómeros E/Z del C=N y orientaciones de los N-sustituyentes del amidinio: CREST/GFN2 no cruza esas barreras, hay que sembrarlos explícitamente y comprobar con qué isómero quedaron los 6 conformeros de B.

**Macro-pKa que se reportan:**

- **pKaM1** (+1→0): A → {D, D′, R, R′}. Primer COOH con el amidinio presente.
- **pKaM2** (0→−1): {D, D′, …} → {B, 010, 100}. Segundo COOH.
- **pKaM3** (−1→−2): {B, …} → {C, C′}. N–H del amidinio con 2×COO⁻.

**Micro-pKa (diagnóstico):** A→D (Gly), A→D′ (Val), D→B, D′→B, A→R. Cierre de ciclo: pK(A→D)+pK(D→B) = pK(A→D′)+pK(D′→B); la discrepancia tras calibrar es el error interno.

**Cuidado con la clase de pKaM1.** Hoy el código acepta exactamente dos tautómeros en el macroestado neutro (R y D). A→R, como microscópico, es la desprotonación del amidinio y es minoritaria: Chuang et al. 2014 informan que el N imínico sigue protonado en medio ácido. La clase `carboxyl` solo es correcta si el camino dominante es A→D/D′; hay que comprobarlo tras el cálculo.

**Riesgo en las geometrías de partida.** `A.xyz` y `D.xyz` dicen "CREST best conformer" pero parecen armadas a mano (D con coordenadas redondeadas y C–O ≈ 1.43 Å; A con O4–H = 0.87 Å; B.xyz con N9–H = 1.17 Å). Las salidas DFT de A/R/D/C no están en disco, así que no se puede confirmar que se optimizaron desde una estructura sana.

El pKa de MGLYVAL **nunca se ha calculado** en el layout actual (`git show HEAD:calibration/pka_mglyval_raw.log` dice "Sin datos" en los 5 estados). El "18.02" citado en un comentario de `maaqc.py` viene de la ronda antigua born/CPCM.

## 2. Calibradores

Fuentes: **[W]** Williams/Evans, compilación de pKa ACS Org. Div. (4/2022); **[F]** Fitch et al., Protein Sci. 2015, 24, 752; **[Ch]** Chuang et al., Mar. Drugs 2014. Recomendación: declarar una sola compilación primaria (IUPAC / NIST SRD 46) para todos los valores. Inconsistencias a resolver: Asp en el código (McMurry) 1.88/3.65/9.60 frente a [W] 2.09/3.86/9.82; 4-AP 9.17 frente a 8.96 "termodinámico" en la misma compilación.

**Sitios (i)/(ii).** Gly/Val zwitteriónicos sirven para la electrostática (+1→0), pero en la MAA la carga está deslocalizada y el N es sp², así que se espera un pKa algo mayor que el de Gly (2.34). Concuerda con porphyra-334: 2.92 ± 0.11 (RMN, D₂O, 25 °C) [Ch]. Mejores análogos: creatina 2.67, betaína 1.84, sarcosina 2.23, N,N-dimetilglicina 2.08 (todos [W]). N-acetilglicina e hipurato **no** sirven para M1 (vecino neutro: transición 0→−1, otra clase). Glicociamina: solo valores predichos, excluir. Para M2 (segundo COOH): IDA 2.98 y MeIDA 2.81, más Asp β 3.86 y Glu γ 4.25. El rango físico de la clase +1→0 es ~2.4 unidades, por eso el modelo es solo offset.

**Sitio (iii).** Los pKa experimentales por encima de 12 no son fiables (arginina 12.48 [W] frente a 13.8 [F]); no usar guanidinios en el ajuste. Diseño de **doble diferencia**, que respeta la protección cross-class:

`pKa3,corr = pKa3,calc − δ(amonio, −1→−2) − [δ(amidinio, +1→0) − δ(amonio, +1→0)]`

El término de carga sale de calibradores con la misma carga; la corrección química, de un par con la misma carga (+1→0). Supone que ambos errores se suman: se comprueba con el test de aditividad (piridinio con COO⁻).

| # | Molécula | Estados | Át. | pKa exp. (fuente) | Sitio / clase | Estado |
|---|---|---|---|---|---|---|
| 1 | Glicina | +1/0/−1 | 11 | 2.34; 9.60 [W] | carboxyl | hecho (v3) |
| 2 | Valina | +1/0/−1 | 20 | 2.32; 9.62 [W] | carboxyl | hecho |
| 3 | Asp | +1/0/−1 (+AC −2) | 17 | 2.09/3.86/9.82 [W] | carboxyl, carboxyl_2, dianión | falta AC |
| 4 | Glu | hasta −2 | 20 | 2.19/4.25/9.67 [W] | idem | hecho |
| 5 | Creatina | +1/0 | 19 | 2.67 [W] | carboxyl (análogo de Gly) | nuevo |
| 6 | Betaína | +1/0 | 20 | 1.84 [W] | carboxyl | nuevo |
| 7 | β-Alanina | +1/0 | 14 | 3.60 [W] | carboxyl (rango) | nuevo |
| 8 | GABA | +1/0 | 17 | 4.23 [W] | carboxyl (rango) | nuevo |
| 9 | IDA | +1/0/−1/−2 | 17 | 2.98; 9.89 [W]; pKa1 ≈1.87 NV | carboxyl_2 y dianión | nuevo |
| 10 | MeIDA | 0/−1/−2 | 19 | 2.81; 10.18 [W] | carboxyl_2 y dianión | nuevo |
| 11 | 4-AP | +1/0 | 14 | 9.17 [W ref. 41] | amidinio | hecho (único) |
| 12 | 2-AP | +1/0 | 14 | 6.86 [W] | amidinio | nuevo |
| 13 | Benzamidina | +1/0 | 18 | 11.6 [W] | amidinio | nuevo |
| 14 | Acetamidina | +1/0 | 11 | 12.52 [W] | amidinio | nuevo |
| 15 | 4-DMAP | +1/0 | 20 | 9.6 (NV) | amidinio | opcional |
| 16 | Metilamina | +1/0 | 7 | 10.62 [W] | amonio | nuevo |
| 17 | Etanolamina | +1/0 | 12 | 9.50 [W] | amonio | nuevo |
| 18 | Tris | +1/0 | 20 | 8.10 [W] | amonio | nuevo |
| 19 | Éster de glicina | +1/0 | 16 | 7.75 (éster NV) | amonio | nuevo |
| 20–22 | Nicotínico, isonicotínico, piridina | +1/0/−1 | 11–14 | 4.77; 4.90; 5.17 [W] | test de aditividad | nuevo |
| H1–H4 | Sarcosina, N,N-DMGly, glicilglicina, fenilguanidina | +1/0 | 14–18 | 2.23; 2.08; 3.06; 10.88 [W] | **validación externa** | — |
| H5 | Arginina | 0/−1 | 26 | 12.48 [W] / 13.8 [F] | solo sensibilidad | — |
| X | **Porphyra-334** | +1/0/−1 | 46 | 2.92 ± 0.11 [Ch] | **validación externa MAA** (M1/M2) | opcional, alto valor |

IDA y las amidinas tienen microestados equivalentes: la degeneración da ±RT ln 2 = 0.30 unidades y se cuenta explícitamente. En MGLYVAL g=1.

## 3. Diseño estadístico (declarado antes de mirar MAAs)

1. **Modelo por clase:** offset por defecto. Pendiente libre solo con n≥5, rango de pKa calculado ≥3 y AICc(lineal) < AICc(offset) − 2. En la práctica solo `amidinio` (6.9–12.5) podría usar pendiente.
2. **LOO en todos los modos** (hoy el modo lineal devuelve NaN).
3. **Validación externa:** H1–H4 y porphyra-334 no entran en el ajuste.
4. **Criterios de aceptación:** RMSE_LOO ≤ 1.0 por clase (≤ 0.5 para carboxyl); |error medio externo| ≤ 1.0; cierre de ciclo de micro-pKa ≤ 0.5; aditividad |[δ(piridinio, 0→−1) − δ(piridinio, +1→0)] − [δ(amonio, 0→−1) − δ(amonio, +1→0)]| ≤ 1.0 (si falla, la doble diferencia queda solo como sensibilidad).
5. **Ya calculable con `git show HEAD:calibration/calib_v3.csv`:** carboxyl n=4 offset −1.36, RMSE_LOO ≈ 0.20; carboxyl_2 n=2, RMSE_LOO ≈ 0.91; amino n=2, ≈ 1.15; amino_dianion y pyridinium n=1 (LOO indefinido).
6. **Monte Carlo con error de calibración:** por paso calibrado, `b ~ b̂ + t_{n−1}·s_LOO·√(1+1/n)`, con el mismo sorteo por clase compartido entre todos los usos de esa clase. En pasos calibrados no se muestrea σ_H⁺ ni el σ_elec/σ_gcorr genérico (el residuo ya los contiene); se suma la incertidumbre conformacional (bootstrap sobre conformeros). **El código actual no hace esto** (`run_chain_montecarlo` ~l.16128 recibe estados crudos).
7. **Protección cross-class:** se mantiene; nuevos tipos en `EMPIRICAL_TYPE_META` (`amidinium`, `ammonium`, `ammonium_dianion`, `pyridinium_anion`). La doble diferencia es un tipo compuesto explícito, no `--force-charge-class`.

## 4. Especificación congelada (calibradores y MAA idénticos)

| Elemento | Valor |
|---|---|
| Muestreo | CREST 3.0.2 (GFN2/ALPB-agua) → CENSO 3.0.9 (CENSO ausente = error fatal). Peso ≥5 %, ≥3 conformeros; semillas E/Z explícitas |
| Opt / freq | ωB97X-D3/def2-SVP CPCM(agua) RIJCOSX TightSCF TightOpt; NumFreq. Sin imaginarias \|ν\|>20 cm⁻¹ (si las hay, reoptimizar) |
| SP primario | ωB97X-D3/def2-TZVPD **SMD(agua)** sobre geometría CPCM; CPCM como sensibilidad |
| G térmica | GoodVibes qRRHO (Grimme, 100 cm⁻¹), H armónica, 298.15 K |
| Método | direct: ΔG = G_aq(base) − G_aq(ácido) + G(H⁺,aq) |
| G(H⁺,aq) | −270.286 kcal/mol (−6.28 + 1.894 − 265.9) |
| Microestados | G_macro = −RT ln ΣᵢΣ_conf exp(−G/RT) sobre todos los microestados de la misma carga, con degeneración gᵢ |
| Temperatura | Calibración solo a 298.15 K |

SMD como primario porque v3 ya está en SMD y está parametrizado para iones; sigue siendo una elección ciega porque el pKa de MGLYVAL aún no se calculó. Congelar con tag `pka-spec-v1` antes de P1-4′.

## 5. Decisiones y riesgos

- **D1′** Valor primario de pKaM3: doble diferencia si pasa aditividad; si no, clase dianión sola. "B domina a pH 8.1" solo se acepta si P(pKaM3 > 9.1) ≥ 0.95 en *todas* las variantes.
- **D5** SMD primario. **D6** Incluir porphyra-334 (1–2 semanas). **D7** Excluir guanidinios del ajuste. **D8** Fuente única de pKa (IUPAC / NIST SRD 46).
- Riesgos: clase dianión de rango estrecho (9.7–10.2); aditividad carga+química puede fallar; E/Z y tautómeros pueden mover los macro-pKa 0.3–1 unidad; fuerza iónica de los datos (IDA I=0.03) frente a I=0 del cálculo; pKa de porphyra-334 medido en D₂O con conversión empírica; ~10 días de máquina para los cuatro microestados nuevos de MGLYVAL.
- **No verificado:** 4-DMAP 9.6, IDA pKa1 1.87, N-acetilglicina 3.67, hipurato, N-acetilvalina, N-fenilglicina, glicociamina, qué éster de glicina corresponde a 7.75, nota de fuerza iónica de IDA, asignación exacta del 2.92 en [Ch], tiempos de cómputo (estimados), Hessianas analíticas con RIJCOSX+CPCM en ORCA 6.1.1.

Fuentes: [Williams/Evans](https://organicchemistrydata.org/hansreich/resources/pka/pka_data/pka-compilation-williams.pdf), [Fitch 2015](https://pmc.ncbi.nlm.nih.gov/articles/PMC4420524), [Chuang 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4178487).
