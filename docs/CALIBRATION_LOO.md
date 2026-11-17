# Offsets y LOO por clase (tarea P1-10)

Recalculado el 2026-11-17 desde `git show HEAD:calibration/calib_v3.csv` (10 puntos, goodvibes, SMD, directo). Script: `calib_loo.py` en el scratchpad de la sesión (pandas/numpy, no versionado).

## Método

- Error crudo = pKa_calc − pKa_exp (columna `error` del CSV). Offset = media de (pKa_exp − pKa_calc), con el mismo signo que `intercept` en `maaqc.py` (~l.11517); es decir, offset = −error medio.
- RMSE = desviación estándar poblacional (divisor n) de los residuos tras aplicar el offset, igual que `maaqc.py` (~l.11519).
- RMSE_LOO: se quita el punto k, se recalcula el offset con los restantes y se evalúa el residuo del punto k. n = 1 → NaN. n = 2 → el offset de cada fold es el sesgo del otro punto, de modo que RMSE_LOO = |b1 − b2| (el doble del RMSE).
- RMSE_crudo = raíz del error cuadrático medio sin corregir (sin offset).
- Carga del producto: carga de la especie desprotonada del paso (la clase define el cambio de carga).

## Tabla por clase

| Clase | Paso (carga) | n | Offset medio | Error crudo medio | RMSE (con offset) | RMSE_LOO | RMSE crudo |
|---|---|---|---|---|---|---|---|
| carboxyl | +1→0 | 4 | −1.364 | +1.364 | 0.149 | **0.198** | 1.372 |
| carboxyl_2 | 0→−1 | 2 | −1.858 | +1.858 | 0.454 | 0.908 | 1.912 |
| amino | 0→−1 | 2 | −2.566 | +2.566 | 0.576 | 1.153 | 2.630 |
| amino_dianion | −1→−2 | 1 | −6.009 | +6.009 | NaN | NaN | 6.009 |
| pyridinium | +1→0 | 1 | −1.575 | +1.575 | NaN | NaN | 1.575 |

Puntos individuales (error crudo): carboxyl: Asp 1.135, Glu 1.507, Gly 1.484, Val 1.331; carboxyl_2: Asp 2.312, Glu 1.404; amino: Gly 1.990, Val 3.143; amino_dianion: Glu 6.009; pyridinium: AP 1.575.

## Error crudo frente a carga del producto

| Carga del producto | Clases | n | Error crudo medio | Rango |
|---|---|---|---|---|
| 0 | carboxyl, pyridinium | 5 | +1.41 | 1.14 a 1.58 |
| −1 | carboxyl_2, amino | 4 | +2.21 | 1.40 a 3.14 |
| −2 | amino_dianion | 1 | +6.01 | 6.01 |

| Clase | Carga producto | Error crudo medio |
|---|---|---|
| carboxyl (+1→0) | 0 | +1.36 |
| pyridinium (+1→0) | 0 | +1.58 |
| carboxyl_2 (0→−1) | −1 | +1.86 |
| amino (0→−1) | −1 | +2.57 |
| amino_dianion (−1→−2) | −2 | +6.01 |

El error sobreestima el pKa siempre (signo positivo) y crece con la carga negativa del producto: ~1.4, ~2.2 y ~6.0 unidades para carga 0, −1 y −2. Es consistente con una estabilización insuficiente de aniones en SMD. Con solo un punto en −2 (Glu, además de otro tipo químico) el valor de −2 no permite separar el efecto de carga del efecto químico. Dentro de la misma carga el error depende además de la clase: amino (−1) frente a carboxyl_2 (−1) difieren en 0.7.

## Contraste con las fuentes

- `calibration_v3.json` no guarda offset, RMSE ni LOO por clase: solo `points` (con `error`) y `EMPIRICAL_CALIBRATION` (tabla de calc/exp). Los valores de este documento se derivan de los mismos puntos; los errores individuales coinciden con los de `docs/METHODS.md` (carboxyl +1.14 a +1.51; carboxyl_2 +2.31 y +1.40; amino +1.99 y +3.14; amino_dianion +6.01; pyridinium +1.58).
- `README.md` y `docs/METHODS.md` no citan cifras de RMSE ni LOO; el valor 0.20 aparece solo como criterio en `TODO.md` (P1-10). METHODS.md solo afirma que n = 1 da RMSE indefinido sin validación cruzada, lo que se confirma (amino_dianion, pyridinium).
- Criterio de TODO: carboxyl RMSE_LOO ≈ 0.20. Resultado 0.198: **se cumple**. El RMSE sin LOO es 0.149.
- Reproduce la fórmula de `maaqc.py` (bucle LOO en ~l.11523 a 11527) en modo offset.

## Notas y advertencias

- Con n = 2, RMSE_LOO (0.91 y 1.15) es poco informativo: es |b1 − b2| y se debe leer como dispersión entre los dos calibradores, no como incertidumbre de predicción.
- Los pKa experimentales de Asp en el CSV (1.88 / 3.65) son los de McMurry; la compilación [W] da 2.09 / 3.86 (P1-11). Con 2.09 el error de Asp en carboxyl bajaría a 0.93 y el RMSE_LOO de carboxyl cambiaría; habrá que rehacer este cálculo tras P1-11.
- Los pKa calculados provienen del esquema A→R/D; la validez de la clase carboxyl para MGLYVAL depende de que A→D/D′ sea el camino dominante (ver `CALIBRATION_PLAN.md`).
