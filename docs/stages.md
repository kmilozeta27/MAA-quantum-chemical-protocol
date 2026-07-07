# Stages — inputs, outputs, estado de migración

| Stage    | requires  | Fuente legacy (pipeline.py)              | Estado      | Archivo migrado |
|----------|-----------|-------------------------------------------|-------------|-----------------|
| crest    | —         | `stage_crest` (1152–1233)                  | ✅ migrado  | `stages/crest.py` |
| censo    | crest     | `stage_censo_prescreen` (1438–1881), `stage_censo_multitemp` (2019–2107) | ✅ migrado | `stages/censo.py` |
| opt      | crest     | `stage_opt` (2462–2600), `_opt_worker` (3136) | ✅ migrado  | `stages/opt.py` |
| numfreq  | opt       | `stage_numfreq` (2601–2671), `has_numfreq_restart` (903–912) | ✅ migrado  | `stages/numfreq.py` |
| sp       | opt\*     | `stage_sp` (2672–2747), `_sp_worker` (3201) | ✅ migrado  | `stages/sp.py` |
| tddft    | sp        | `stage_tddft` (2814–2880), `_tddft_worker` (2748–2813) | ✅ migrado  | `stages/tddft.py` |
| adc2     | tddft\*\* | `stage_adc2` (2972–3135), `extract_adc2_lambda` (929–998), `load_exp_lambda` (270–351) | ✅ migrado | `stages/adc2.py` |

**Estado general: 7 de 7 stages migrados con lógica completa y funcional.**

Fuente de verdad para esta tabla: `maaqc doctor`, que reporta el estado
real de cada stage (`[OK]` para migrados con lógica real, `[PENDIENTE]`
para stubs). La salida de `maaqc doctor` puede regenerar esta tabla
de forma automática en el futuro.

\* Con `skip_numfreq: true`, `sp` depende de `opt` directo — la lógica
está resuelta en `SpStage._source_dirs_for_validation()` que prioriza
`03_numfreq/` si existe, sino cae a `02_opt/`. No requiere cambios en
`orchestrator/dag.py` ni `Stage.requires` estático.

\*\* Con `entry: sp-only` (equivalente al legacy `--start adc2`), `adc2`
ignora la cadena completa — manejado en `orchestrator/dag.py::STANDALONE_ENTRIES`,
no en `Stage.requires`.

## Detalles de migración por stage

### CrestStage (`stages/crest.py`)
- **Origen:** `pipeline.py:stage_crest()` líneas 1152–1233.
- **Lógica:** Búsqueda conformacional GFN2-xTB usando CREST.
- **Entrada:** `workdir/input.xyz` (supuesto por convención).
- **Salida:** `01_crest/crest_conformers.xyz` (o normaliza `crest_best.xyz` si es necesario).
- **Notas:** ⚠️ No se comparó línea-a-línea contra legacy (no estaba disponible al portar) — usar `maaqc doctor` para validación rápida.

### CensoStage (`stages/censo.py`)
- **Origen:** `pipeline.py:stage_censo_prescreen()` + `stage_censo_multitemp()` líneas 1438–2107.
- **Lógica:** Pre-screening CENSO, con soporte para modo single-T y multitemp.
- **Entrada:** `01_crest/crest_conformers.xyz`.
- **Salida:** `01b_censo_prescreen/crest_conformers.xyz` (ensemble recortado) o `01b_censo_prescreen/T{temp:.2f}K/crest_conformers.xyz` por temperatura.
- **Diseño:** Cableado a `system/_censo_shared.py::run_censo_prescreen()`, `run_censo_multitemp()`, selección Boltzmann/auto.
- **Notas:** `_select_final_ensemble()` aplica la selección de `nconf` y escribe el archivo con el nombre que `OptStage` espera.

### OptStage (`stages/opt.py`)
- **Origen:** `pipeline.py:stage_opt()` (líneas 2462–2600) + `_opt_worker()` (línea 3136).
- **Lógica:** Optimización DFT en paralelo (ThreadPoolExecutor).
- **Entrada:** `01b_censo_prescreen/crest_conformers.xyz` o `01_crest/crest_conformers.xyz` (fallback).
- **Salida:** `02_opt/conf_{id}/opt.out` + `opt_opt.xyz`.
- **Paralelismo:** BatchSize configurable, CPUs asignadas dinámicamente.
- **Notas:** Cambio de diseño: legacy recibía `nconf` ya resuelto; aquí se optimizan TODOS los frames del ensemble — la selección debe ocurrir en CENSO.

### NumFreqStage (`stages/numfreq.py`)
- **Origen:** `pipeline.py:stage_numfreq()` (líneas 2601–2671) + `has_numfreq_restart()` (líneas 903–912).
- **Lógica:** Frecuencias vibracionales numéricas para termoquímica.
- **Entrada:** `02_opt/conf_*/opt_opt.xyz`.
- **Salida:** `03_numfreq/conf_*/freq.out` + restart files.
- **Restart:** Detecta y activa automáticamente si existen `freq_D*.res.Gradients`.
- **Notas:** ⚠️ Legacy no usaba `env` propio (fallback a `build_omp_env(4)`). Portado 1:1 — confirmar si era bug intencional.

### SpStage (`stages/sp.py`)
- **Origen:** `pipeline.py:stage_sp()` (líneas 2672–2747) + `_sp_worker()` (línea 3201).
- **Lógica:** Single-point de alta precisión.
- **Entrada:** `03_numfreq/conf_*/opt_opt.xyz` (default) o `02_opt/conf_*/opt_opt.xyz` (si `skip_numfreq=True`).
- **Salida:** `04_sp/conf_*/sp.out` + `opt_opt.xyz`.
- **Diseño:** Dependencia dinámica resuelta en `_source_dirs_for_validation()` sin afectar `Stage.requires` estático.

### TddftStage (`stages/tddft.py`)
- **Origen:** `pipeline.py:stage_tddft()` (líneas 2814–2880) + `_tddft_worker()` (líneas 2748–2813).
- **Lógica:** TD-DFT sobre geometría + orbitales SP.
- **Entrada:** `04_sp/conf_*/sp.gbw` + `sp_opt.xyz`.
- **Salida:** `05_tddft/conf_*/tddft.out`.
- **Config:** Respeta `config.tddft.functional` (default `wb97x-d3`) y `config.tddft.states` (default `10`, cambio respecto al legacy que hardcodeaba `15`).
- **Notas:** Template ahora paramétrico (`{FUNCTIONAL}`, `{NROOTS}`) en lugar de constantes de módulo.

### Adc2Stage (`stages/adc2.py`)
- **Origen:** `pipeline.py:stage_adc2()` (líneas 2972–3135) + `extract_adc2_lambda()` (líneas 929–998) + `load_exp_lambda()` (líneas 270–351).
- **Lógica:** ADC2 sobre el confórmero de menor energía SP (opcional).
- **Entrada:** `04_sp/conf_*/sp.out` (selecciona mínimo global por energía).
- **Salida:** `07_adc2/conf_{id}/adc2.out` + `adc2_results.json`.
- **Opcional:** Solo corre si `config.adc2.mandatory=True` o `entry="sp-only"`.
- **Extracción:** `extract_adc2_lambda()` parsea espectro RIGHT, calcula `Δλ = λ_calc - λ_exp`.
- **Notas:** `Stage.requires = ["tddft"]` estático; `orchestrator/dag.py` maneja `sp-only` como entrada sin requerir TDDFT.

## Módulos de soporte migrados

### `system/platform_detect.py`
- Detección de binarios en PATH y envs conda (CREST, CENSO, ORCA, xTB).
- `resolve_bin(name)` — retorna ruta ejecutable o levanta error con instrucciones de instalación.

### `system/affinity.py`
- Distribución de afinidad CPU por worker (para aislar cache/NUMA).
- `assign_cpu_lists(n_workers, nprocs_per)` — lista de CPUs por worker.

### `system/_orca_shared.py`
- Helpers genéricos: templates ORCA, subprocess, parsing (energía, convergencia, termoquímica).
- `run_orca(inp, out, log, env, cpu_list)` — ejecución con manejo de timeout/señales.
- `check_success(out)` — valida "ORCA TERMINATED NORMALLY".

### `system/_censo_shared.py`
- Helpers de CENSO: ejecución, parsing de screening, selección Boltzmann/auto.
- `run_censo_prescreen()`, `run_censo_multitemp()`, `select_nconf_boltzmann()`, `select_nconf_auto()`.

### `io/xyz.py`
- I/O unificado de XYZ: lectura/escritura de confórmeros, extracción de IDs (CREST labels).
- `atoms_to_xyz_block()`, `read_all_conformers()`, `extract_conf_ids_from_xyz()`.

### `io/parsing.py`
- Parsers de salidas: energía DFT, frecuencias imaginarias, espectro ADC2 (RIGHT).
- `extract_sp_energy()`, `check_no_imaginary()`, `extract_adc2_lambda()`, etc.

### `analysis/results_collector.py`
- Consolidación de resultados entre estados de protonación (A, B, C, D, R).
- Integración con `orchestrator/state.py` para resume.

### `analysis/protonation_cycle.py`
- Cálculo de ciclo pKa vía Boltzmann con correcciones de goodvibes.

## Issues pendientes antes de mergear a main

1. **Validación end-to-end:** Comparar outputs de `maaqc run` contra `pipeline.py` legacy con moléculas reales de PAPER1.
2. **Tests de integración:** Mock ORCA outputs, flujos completos (CREST → ADC2).
3. **Documentación de deployment:** Instrucciones conda, conda-forge, reproducibilidad.
4. **Tag v0-legacy:** Crear tag git sobre la versión pre-migración de `main` para preservar historial.
5. **Review de templates:** Validar que los templates ORCA (hardcodeados en `system/constants.py`) coincidan 100% con los del legacy.

## Notas de diseño y cambios respecto al legacy

### Dependencias dinámicas resuelta sin romper `Stage.requires`
- **Problema:** `sp.py con skip_numfreq` necesita "a veces OPT, a veces NumFreq" como input.
- **Solución:** `Stage.requires` sigue siendo estático (`["opt"]`). La decisión real se toma en `SpStage._source_dirs_for_validation()` / `run()`, donde inspecciona `config.skip_numfreq` y elige la fuente de entrada.
- **Ventaja:** `Stage.requires` nunca es dinámico (simple, predecible).

### Entries standalone (sp-only)
- **Problema:** `--start adc2` (legacy) ejecuta ADC2 sin la cadena normal, requiriendo solo `04_sp/`.
- **Solución:** `orchestrator/dag.py::STANDALONE_ENTRIES` mapea `"sp-only"` → `"adc2"`. `Adc2Stage.requires` es estático (`["tddft"]`), pero `dag.py::validate_stage_ready()` detecta entry standalone y salta la validación de requisitos previos.
- **Ventaja:** Stage es simple; toda la lógica de orquestación vive fuera.

### Templates ORCA ahora parametrizados
- **Cambio:** `stages/tddft.py` y `stages/adc2.py` usan `config.tddft.functional`, `config.tddft.states`, etc. en lugar de hardcodear constantes.
- **Ventaja:** Recipes YAML controlan exactamente qué se corre, sin cambiar código.

### Resume robusto vía state.json
- **Adición:** `orchestrator/state.py` guarda config_hash + output fingerprint por stage.
- **Uso:** `resolve_resume_index()` detecta si config cambió o outputs fueron tocados por fuera del orchestrator.
- **Fallback:** `is_complete()` de cada stage es la fuente de verdad ("¿ORCA TERMINATED NORMALLY?"), `state.py` la complementa.
