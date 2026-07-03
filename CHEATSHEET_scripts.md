# Cheatsheet — pipeline.py / protonation_cycle.py / results_collector.py

Guía rápida de uso para no enredarse con las opciones. `conda activate sciencebitch` antes de correr cualquiera de los tres.

---

## 1. `pipeline.py` — CREST → CENSO → OPT → NumFreq → SP → TDDFT → ADC2

CREST y CENSO (pre-screening) son **siempre obligatorios** al inicio de una molécula nueva. El nº de confórmeros que pasa a OPT se decide solo (mín 3, máx 20), salvo que uses `--boltzmann-threshold`.

### 1.1 Flags básicos

| Flag | Default | Qué hace |
|---|---|---|
| `--mol RUTA.xyz` | — | Input. Obligatorio salvo con `--status`. |
| `--charge N` | `0` | Carga formal. |
| `--start {crest,opt,numfreq,sp,tddft,adc2}` | `crest` | Etapa desde la que arranca. `crest` y `opt` son equivalentes al inicio (ambos corren CREST+CENSO). Usa `opt/numfreq/sp/tddft` solo para **retomar** un run con etapas previas ya en disco. |
| `--status` | — | Imprime estado de todos los runs y sale (no requiere `--mol`). |
| `--skipnumfreq` | — | Salta NumFreq: SP/TDDFT usan la geometría de `02_opt/` directo. **Sin correcciones térmicas** → no sirve para el ciclo pKa (`06_sp_gas`). |

⚠️ `--skipnumfreq` es **incompatible** con `--start numfreq/sp/tddft` (ya que esos asumen que numfreq ya corrió con correcciones térmicas). Solo combínalo con `--start crest` o `--start opt`.

### 1.2 Recursos de cómputo

Con `--cores` el script autodistribuye por etapa. Los flags avanzados (`--nprocs`, `--maxcore`, etc.) **sobreescriben** lo que calcule `--cores`. Sin ningún flag de recursos, autodetecta la topología real del ThinkPad.

| Flag | Nivel | Afecta |
|---|---|---|
| `--cores N` | simple (recomendado) | Reparte cores entre OPT/NumFreq/SP/TDDFT automáticamente. `maxcore ≈ (RAM_total×0.85)/cores`. |
| `--nprocs N` | avanzado | Cores por job de OPT/SP/TDDFT (default interno `4`). |
| `--maxcore MB` | avanzado | MB/core para OPT/SP/TDDFT (default `3500`). RAM/job = `nprocs × maxcore`. |
| `--nprocs-numfreq N` | avanzado | Cores para NumFreq (default `8`). |
| `--maxcore-numfreq MB` | avanzado | MB/core para NumFreq (default `2800`). |

TDDFT y ADC2 (`nprocs`/`maxcore`) **no son configurables por CLI**, siempre vienen del `ResourceProfile` autodetectado.

**Precedencia:** flags avanzados explícitos > `--cores N` > autodetección pura.

### 1.3 Modo multi-temperatura CENSO (`--multitemp`)

Corre CENSO varias veces (una por temperatura), cada ensemble en `01b_censo_prescreen/T<temp>K/`. El pipeline sigue a OPT usando el ensemble de `--temperature` (default `298.15`).

| Flag | Default | Nota |
|---|---|---|
| `--multitemp` | off | Activa el modo. |
| `--trange T_INI T_FIN N` | `274 298 24` (solo si `--multitemp`) | Rango en K, N puntos equiespaciados incl. extremos. Si `N=1`, usa solo `T_INI`. **Solo tiene efecto junto con `--multitemp`** (error si lo pasas sin él). |
| `--temperature T` | `298.15` | Sin `--multitemp`: única T usada por CENSO. Con `--multitemp`: T de referencia para elegir el ensemble que pasa a OPT (no tiene que coincidir exacto con `--trange`; si no coincide, se añade al rango con warning). |

⚠️ `--multitemp` es **incompatible** con `--nroot` (corre `--nroot` aparte, sobre un run ya completo).

### 1.4 Selección de confórmeros para OPT

| Flag | Default | Nota |
|---|---|---|
| `--boltzmann-threshold PCT` | `5.0` | Toma confórmeros con peso Boltzmann ≥ PCT% desde `1_SCREENING.out` de CENSO. Si deja menos que el mínimo, se completa automáticamente con `select_nconf_auto`. Si `1_SCREENING.out` no existe, cae a `select_nconf_auto` igual. |

*(No hay flag para desactivarlo del todo — si no te sirve el behavior, es el fallback automático el que manda cuando falta el `.out`.)*

### 1.5 Modo SP+NRoots1 (`--nroot`)

Corre un SP con TD-DFT mínimo (NRoots 1) sobre los mejores confórmeros por energía DFT de OPT. **Requiere que OPT ya haya terminado.** Escribe en `06_sp_nroot/<conf_id>/`.

| Flag | Default | Nota |
|---|---|---|
| `--nroot` | off | Activa el modo (standalone: corre y termina, no continúa el pipeline normal). |
| `--nroot-top N` | `5` | Top-N confórmeros por E_OPT. |

### 1.6 `--start adc2` (caso especial)

**No** retoma la cadena normal. Solo necesita `04_sp/` (con `.out` + `opt_opt.xyz` por confórmero) — ignora por completo `01_crest/`, `02_opt/`, `03_numfreq/`, `05_tddft/` aunque no existan. Detecta automáticamente el confórmero de menor energía SP y corre ADC2 solo sobre ese.

Pensado para tu flujo real: como ADC2 necesita ~400 GB, respaldas todo en disco externo y borras lo pesado, dejando solo `04_sp/`.

```bash
python pipeline.py --mol molecules/mol_A.xyz --charge -1 --start adc2
```

### 1.7 Combinaciones típicas

```bash
# Pipeline completo desde cero, 8 cores totales
python pipeline.py --mol molecules/mglyval.xyz --charge -1 --cores 8

# Retomar desde OPT (CREST+CENSO ya en disco)
python pipeline.py --mol molecules/mglyval.xyz --charge -1 --start opt

# Retomar desde SP (opt y freq ya hechos)
python pipeline.py --mol molecules/mglyval.xyz --charge -1 --start sp

# Retomar NumFreq interrumpida (detecta restart si hay freq_D*.res.Gradients)
python pipeline.py --mol molecules/mglyval.xyz --charge -1 --start numfreq

# OPT → SP → TDDFT directo, sin correcciones térmicas
python pipeline.py --mol molecules/mol_A.xyz --charge -1 --skipnumfreq

# Multitemp con defaults (274–298 K, 24 puntos), OPT usa 298.15 K
python pipeline.py --mol molecules/mol_A.xyz --charge -1 --multitemp

# Multitemp con rango explícito
python pipeline.py --mol molecules/mol_A.xyz --charge -1 \
    --multitemp --trange 274 298 24 --temperature 298.15

# Selección de confórmeros por peso Boltzmann ≥5% en vez de auto
python pipeline.py --mol molecules/mol_A.xyz --charge -1 --boltzmann-threshold 5.0

# SP+NRoots1 sobre top-3 confórmeros por E_OPT (requiere 02_opt/ completo)
python pipeline.py --mol molecules/mglyval.xyz --charge -1 --nroot --nroot-top 3

# ADC2 standalone (solo 04_sp/ en disco)
python pipeline.py --mol molecules/mglyval.xyz --charge -1 --start adc2

# Ver estado de todos los runs
python pipeline.py --status
```

### 1.8 Combinaciones inválidas (el script las rechaza)

- `--trange` sin `--multitemp` → error.
- `--multitemp` junto con `--nroot` → error.
- `--skipnumfreq` junto con `--start numfreq/sp/tddft` → error (solo vale con `--start crest` u `opt`).
- Correr sin `--mol` y sin `--status` → error.
- `--mol` cuyo archivo no exista o no sea `.xyz` → error.

---

## 2. `protonation_cycle.py` — pKa (ciclo de Born) + Boltzmann entre estados

Compara AH₂ (q=0) → AH⁻ (q=−1) → A²⁻ (q=−2) usando el ciclo termodinámico de Ho & Coote. Requiere que existan en `runs/<estado>/`: `03_numfreq/`, `04_sp/` (SP acuoso CPCM) y, para el método `born`, `06_sp_gas/` (SP fase gas).

Los nombres que pasas en `--neutral/--anion/--dianion` son **nombres de carpeta dentro de `runs/`** (p. ej. `R`, `A`, `B`, `C`), no rutas.

### 2.1 Flags

| Flag | Default | Nota |
|---|---|---|
| `--neutral DIR` | — | Carpeta de AH₂ (q=0). |
| `--anion DIR` | — | Carpeta de AH⁻ (q=−1). |
| `--dianion DIR` | — | Carpeta de A²⁻ (q=−2). Opcional: sin este solo se calcula pKa₁. |
| `--ph VALOR` | pH antártico (`PH_ANTARCTIC`) | pH puntual de referencia para la distribución. |
| `--ph-range MIN MAX` | `0 14` | Rango de pH para el diagrama. |
| `--ph-step PASO` | `0.1` | Paso del eje pH. |
| `--method {auto,born,direct}` | `auto` | `born`: ciclo completo (requiere SP fase gas). `direct`: solo CPCM, sin fase gas. `auto`: usa `born` si hay SP gas, si no cae a `direct`. |
| `--export {csv,xlsx,both}` | ninguno | Exportar resultados. |
| `--plot` | off | Gráficos de distribución + Boltzmann conformacional. |
| `--verbose` | off | Muestra contribuciones paso a paso del ciclo de Born. |

⚠️ Se requiere **al menos** `--neutral` y `--anion` (error si faltan).

### 2.2 Combinaciones típicas

```bash
# Caso más común: 3 estados, pKa1 + pKa2 + distribución
python protonation_cycle.py --neutral R --anion B --dianion C

# Solo pKa1 (R → B)
python protonation_cycle.py --neutral R --anion B

# Distribución a un pH específico además del rango completo
python protonation_cycle.py --neutral R --anion B --dianion C --ph 8.1

# Exportar resultados con gráficos
python protonation_cycle.py --neutral R --anion B --dianion C --export xlsx --plot

# Modo detallado (contribuciones del ciclo de Born paso a paso)
python protonation_cycle.py --neutral R --anion B --dianion C --verbose

# Todo junto, pH antártico explícito
python protonation_cycle.py \
    --neutral R --anion B --dianion C \
    --ph 8.1 --export xlsx --plot --verbose
```

---

## 3. `results_collector.py` — recolección de resultados + comparación entre estados

Modo obligatorio y **mutuamente excluyente**: `--mol NAME` (un estado) **o** `--compare` (todos los estados en `runs/`). Uno de los dos es requerido.

### 3.1 Flags

| Flag | Requiere | Default | Nota |
|---|---|---|---|
| `--mol NAME` | — | — | Analiza un solo estado (carpeta en `runs/`). |
| `--compare` | — | — | Compara todos los estados de protonación en `runs/`. |
| `--ref NAME` | solo con `--compare` | ninguno | Estado de referencia para ΔE relativa. |
| `--export {csv,xlsx,both}` | — | ninguno | Exportar tablas. |
| `--plot` | — | off | Genera espectro(s) UV-Vis. |
| `--goodvibes` | solo con `--mol` | off | Corre GoodVibes (mRRHO Grimme) sobre cada `freq.out` y compara contra el G_corr de ORCA (RRHO). Exporta el apéndice a `runs/<mol>/`. |
| `--adc2` | — | off | Extrae λmax ADC(2) y compara con TD-DFT. Con `--mol`: analiza ese estado. Con `--compare`: compara todos los estados y evalúa sistematicidad del desplazamiento TDDFT−ADC2. |
| `--sigma NM` | — | `20.0` | Anchura Gaussiana (nm) para los plots. |
| `--xmin NM` | — | `200.0` | Límite inferior eje λ. |
| `--xmax NM` | — | `500.0` | Límite superior eje λ. |

Nota de comportamiento: con `--compare --adc2 --plot`, el plot de ADC2-vs-TDDFT reemplaza al plot de comparación normal (no se generan ambos).

### 3.2 Combinaciones típicas

```bash
# Un solo estado
python results_collector.py --mol mglyval
python results_collector.py --mol mglyval --export xlsx
python results_collector.py --mol mglyval --plot

# GoodVibes vs ORCA para ese estado
python results_collector.py --mol mglyval --goodvibes
python results_collector.py --mol mglyval --goodvibes --export csv

# ADC2 vs TDDFT para un solo estado
python results_collector.py --mol mglyval --adc2
python results_collector.py --mol mglyval --adc2 --export xlsx

# Comparar todos los estados
python results_collector.py --compare
python results_collector.py --compare --plot
python results_collector.py --compare --ref mglyval_neutral

# Comparar ADC2 vs TDDFT en todos los estados
python results_collector.py --compare --adc2 --plot
```

---

## 4. Flujo end-to-end recomendado (un estado de protonación nuevo)

```bash
conda activate sciencebitch

# 1. Pipeline completo hasta ADC2 (usa todos los cores del ThinkPad, 12 razonable)
nohup python pipeline.py --mol molecules/mol_R.xyz --charge -1 --cores 12 &

# 2. Si necesitas el ciclo pKa, corre también el SP en fase gas (fuera de pipeline.py,
#    ver docstring de protonation_cycle.py: va en runs/<estado>/06_sp_gas/)

# 3. Revisar resultados de ese estado
python results_collector.py --mol mol_R --export xlsx --plot

# 4. Repetir 1–3 para cada estado (A, B, C, D, R...)

# 5. Comparar todos los estados entre sí
python results_collector.py --compare --ref mol_R --plot --adc2

# 6. Ciclo de pKa + distribución de Boltzmann entre estados
python protonation_cycle.py --neutral R --anion B --dianion C --export xlsx --plot --verbose
```

## 5. Reglas rápidas para no enredarte

- **¿Molécula nueva?** → sin `--start` (equivale a `crest`), obligatorio `--mol` y `--charge`.
- **¿Retomando un run interrumpido?** → mira qué carpeta llegó más lejos (`01_crest`, `02_opt`, `03_numfreq`, `04_sp`, `05_tddft`) y usa `--start` con la **siguiente** etapa.
- **¿Ya borraste todo salvo `04_sp/` para liberar espacio antes de ADC2?** → `--start adc2`, nada más.
- **¿Quieres Boltzmann multi-T de CENSO?** → `--multitemp` (nunca `--nroot` a la vez).
- **¿Quieres saltarte NumFreq porque no te importa el ciclo pKa esta vez?** → `--skipnumfreq`, pero solo arrancando desde `crest` u `opt`.
- **¿pKa?** → primero necesitas `06_sp_gas/` en cada estado (fase gas, sin CPCM), luego `protonation_cycle.py`.
- **¿Comparar todo al final?** → `results_collector.py --compare`.
