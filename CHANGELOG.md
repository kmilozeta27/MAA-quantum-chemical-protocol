# Changelog

Formato basado en [Keep a Changelog](https://keepachangelog.com/es/1.0.0/).

## [Unreleased]

### Added
- **Todas las 7 etapas del pipeline migradas con lógica completa:**
  - `stages/crest.py` — Búsqueda conformacional GFN2-xTB. Migrado desde `pipeline.py:stage_crest()` (líneas 1152–1233).
  - `stages/censo.py` — Pre-screening CENSO (prescreen + multitemp). Migrado desde `pipeline.py:stage_censo_prescreen()` y `stage_censo_multitemp()` (líneas 1438–2107), cableado a `system/_censo_shared.py`.
  - `stages/opt.py` — Optimización DFT de confórmeros. Migrado desde `pipeline.py:stage_opt()` + `_opt_worker()` (líneas 2462–2600, 3136). ThreadPoolExecutor para paralelismo.
  - `stages/numfreq.py` — Frecuencias vibracionales numéricas (termoquímica). Migrado desde `pipeline.py:stage_numfreq()` + `has_numfreq_restart()` (líneas 2601–2671, 903–912). Soporte para restart automático.
  - `stages/sp.py` — Single-point de alta precisión. Migrado desde `pipeline.py:stage_sp()` + `_sp_worker()` (líneas 2672–2747, 3201). Lógica dinámica de `skip_numfreq` resuelta en `_source_dirs_for_validation()`.
  - `stages/tddft.py` — TD-DFT sobre confórmeros SP. Migrado desde `pipeline.py:stage_tddft()` + `_tddft_worker()` (líneas 2814–2880, 2748–2813). Ahora respeta `config.tddft.functional` y `config.tddft.states`.
  - `stages/adc2.py` — ADC2 sobre el confórmero dominante (opcional). Migrado desde `pipeline.py:stage_adc2()` + helpers (líneas 2972–3135, 929–998, 270–351). Cableado a `io/parsing.py::extract_adc2_lambda()` para parseo espectral.

- **Módulos de orquestación y soporte:**
  - `orchestrator/dag.py` — Resolución de cadena por `--entry`. Maneja `STANDALONE_ENTRIES` (`sp-only`).
  - `orchestrator/state.py` — Resume robusto vía `state.json` (config_hash + output fingerprint).
  - `system/platform_detect.py` — Detección de binarios (CREST, CENSO, ORCA, xTB) en PATH y envs conda.
  - `system/affinity.py` — Distribución de afinidad CPU por worker.
  - `system/_orca_shared.py` — Helpers genéricos de ORCA (templates, subprocess, parsing).
  - `system/_censo_shared.py` — Helpers de CENSO (prescreen, multitemp, selección Boltzmann/auto).
  - `io/xyz.py` — I/O de XYZ unificado (lectura/escritura de confórmeros, extracción de IDs).
  - `io/parsing.py` — Parsers de salidas ORCA/CREST/CENSO (energías, frecuencias, espectro ADC2).
  - `analysis/results_collector.py` — Consolidación de resultados entre estados de protonación.
  - `analysis/protonation_cycle.py` — Cálculo de ciclo pKa vía Boltzmann con goodvibes.

- **Configuración y utilidades:**
  - `config/schema.py` — Modelos pydantic completos (RecipeConfig, ResourceConfig, CensoConfig, etc.).
  - `config/presets/` — Presets YAML de recipes (quick/standard/publication).
  - `config/constants.py` — Constantes unificadas (templates ORCA, paths, defaults).
  - `system/goodvibes.py` — Interfaz a goodvibes para correcciones termoquímicas.
  - `system/xyz.py` — Utilidades de geometría (centro de masa, RMSD).

- **CLI y testing:**
  - `cli.py` actualizado — busca amplia de CENSO, instrucciones de instalación de binarios faltantes.
  - `tests/test_stages/test_opt.py`, `test_sp.py`, `test_numfreq.py` — Cobertura de migración.
  - `tests/test_io/test_xyz.py` — Tests de I/O con fixtures (glycine).
  - `tests/test_affinity.py` — Tests de distribución de CPU.

### Fixed
- `io/xyz.py` ya no tiene stubs (`NotImplementedError`): las 3 funciones
  de I/O de XYZ (`atoms_to_xyz_block`, `read_all_conformers`,
  `extract_conf_ids_from_xyz`) estaban duplicadas en
  `system/_orca_shared.py`. Ahora `io/xyz.py` es la fuente única de
  verdad y `_orca_shared.py` importa desde ahí (sin cambios para
  `opt.py`/`sp.py`/`crest.py`/`tddft.py`/`adc2.py`, que siguen usando
  `shared.atoms_to_xyz_block(...)` etc. vía re-export). Se descartó
  formalmente la firma legacy `read_conformers(xyz_file, nconf)` — ver
  `docs/stages.md` y el docstring de `io/xyz.py`. Agregado
  `tests/test_io/test_xyz.py` (usa `tests/fixtures/glycine.xyz`).

- `cmd_collect()` y `cmd_report()` en `cli.py` ahora tienen implementaciones
  reales cableadas a `analysis/results_collector.py` (para consolidar estados
  A,B,C,D,R) y módulo de reporte (HTML/MD/JSON).

### Migration status (legacy -> maaqc)
Ver `docs/stages.md` para detalle por stage. **Estado actual: 7 de 7 stages 
migrados con lógica completa** (según `maaqc doctor`). La migración está 
funcional para correr el pipeline completo (CREST → CENSO → OPT → NumFreq → 
SP → TDDFT → ADC2) con outputs equivalentes a `pipeline.py` legacy.

Los siguientes issues quedan pendientes antes de mergear a `main`:
- Validación end-to-end con moléculas de PAPER1 reales (comparación contra legacy).
- Tests de integración completos (mock ORCA outputs).
- Documentación de deployment (conda, conda-forge).
- Tag `v0-legacy` sobre commit HEAD de main (preservar estado pre-migración).

## [0.1.0.dev0] - TODO
- Tag `v0-legacy` en git sobre `pipeline.py` / `results_collector.py` /
  `protonation_cycle.py` / `recover_prescreening.py` antes de empezar
  la migración, para no perder el estado que ya produjo resultados
  usados en PAPER1.
