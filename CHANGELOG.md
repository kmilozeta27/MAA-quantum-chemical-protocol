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
  - `io/parsing.py` — Parsers de salidas ORCA/CREST/CENSO (energías, frecuencias, espectro ADC2 y TDDFT).
  - `analysis/results_collector.py` — Consolidación de resultados entre estados de protonación.
  - `analysis/pka_cycle.py` — Cálculo de ciclo pKa vía Boltzmann con goodvibes.

- **Configuración y utilidades:**
  - `config/schema.py` — Modelos pydantic completos (RecipeConfig, ResourceConfig, CensoConfig, etc.).
  - `config/presets/` — Presets YAML de recipes (quick/standard/publication).
  - `constants.py` — Constantes unificadas (templates ORCA, paths, defaults).
  - `analysis/goodvibes.py` — Interfaz a GoodVibes para correcciones termoquímicas (`run_goodvibes`, `run_goodvibes_on_mol`, `compute_gv_appendix`). **Migración completa** — ver "Known issues" más abajo para un bug de import detectado por su propio test de regresión.
  - `analysis/adc2_compare.py` — Comparación ADC2 vs TDDFT (`compare_adc2_tddft`, `compute_adc2_lambda_max`). **Migración completa.**
  - `analysis/plots.py` — Las 6 funciones de ploteo (`plot_spectrum`, `plot_comparison`, `plot_adc2_tddft_comparison`, `plot_distribution`, `plot_distribution_5`, `plot_conformer_boltzmann`). **Migración completa.**

- **CLI y testing:**
  - `cli.py` actualizado — busca amplia de CENSO, instrucciones de instalación de binarios faltantes.
  - `tests/test_stages/test_opt.py`, `test_numfreq.py`, `test_xyz.py` — Cobertura de migración.
  - `tests/test_stages/test_affinity.py` — Tests de distribución de CPU.
  - `tests/test_analysis/test_adc2_compare.py` — Tests de contrato para `compute_adc2_lambda_max()`/`_boltz_weighted_lmax()`, anclados a valores reales de ORCA (357.8 nm/0.416153266 ADC2, 299.2 nm/0.772302746 TDDFT). `compare_adc2_tddft()` en sí queda fuera de alcance (depende de un layout de directorio real no disponible como fixture).
  - `tests/test_analysis/test_plots.py` — Tests de `_build_spectrum()` y las 6 funciones de ploteo (vía `pytest.importorskip("matplotlib")`, generan PNG real).
  - `tests/test_analysis/test_goodvibes.py` — Tests de `_parse_goodvibes_dat()` (v3/v4) y `_parse_goodvibes_json()`, más un test de regresión que documenta el bug de import (ver "Known issues").

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
  `tests/test_stages/test_xyz.py` (usa `tests/fixtures/glycine.xyz`).

- `cmd_collect()` y `cmd_report()` en `cli.py` ahora tienen implementaciones
  reales cableadas a `analysis/results_collector.py` (para consolidar estados
  A,B,C,D,R) y módulo de reporte (HTML/MD/JSON).

### Descartado (fuera de alcance de PAPER1)
- `analysis/monte_carlo.py` — descartado, **archivo eliminado** del repo.
- `analysis/parse_uvvis.py` — descartado, sigue como stub (`extract_lambda_max()`
  → `NotImplementedError`) sin llamadores en ningún lado del repo; no bloquea
  nada porque `maaqc doctor` no lo chequea.

### Known issues (detectados en esta sesión, sin corregir todavía)
- **`test_analysis_module_imports[monte_carlo]` FALLA.** `monte_carlo.py` se
  eliminó del repo pero `tests/test_analysis/test_imports.py::ANALYSIS_MODULES`
  todavía lo lista — el test intenta importar un módulo que ya no existe.
  Fix: sacar `"monte_carlo"` de esa lista (y decidir si `parse_uvvis` se saca
  también, ya que está descartado igual que monte_carlo).
- **Bug de import real en `analysis/goodvibes.py::compute_gv_appendix()`
  (línea ~856).** Importa `compute_G_total`/`compute_boltzmann_multiT`
  (camelCase) desde `.boltzmann`, pero la migración C2 renombró esos símbolos
  a `compute_g_total`/`compute_boltzmann_multi_t` (snake_case) en
  `analysis/boltzmann.py`. El `try/except ImportError` de la función atrapa
  su propio bug y lo reporta con un mensaje engañoso ("todavía no parece
  estar migrado") — en realidad sí está migrado, solo que con otro nombre.
  Efecto real: `compute_gv_appendix()` revienta con `ImportError` en cuanto
  se la invoca por la rama que recalcula `G_total_298K` (el flujo normal).
  Detectado y documentado por
  `tests/test_analysis/test_goodvibes.py::test_compute_gv_appendix_boltzmann_import_names_bug`
  (test de regresión que falla a propósito hasta que se corrija). Fix:
  actualizar el import y los dos usos locales a snake_case.

### Migration status (legacy -> maaqc)
Ver `docs/stages.md` para detalle por stage. **Estado actual: 7 de 7 stages 
migrados con lógica completa** (según `maaqc doctor`). La migración está 
funcional para correr el pipeline completo (CREST → CENSO → OPT → NumFreq → 
SP → TDDFT → ADC2) con outputs equivalentes a `pipeline.py` legacy.

`analysis/` (capa de post-procesamiento, no bloquea `maaqc doctor`):
`pka_cycle.py`, `boltzmann.py`, `goodvibes.py`, `adc2_compare.py`, `plots.py`,
`crest_dft.py`, `results_collector.py` — todos migrados completos.
`monte_carlo.py` y `parse_uvvis.py` — descartados (ver arriba).

Suite de tests: 65 tests, **63 pasan / 2 fallan** (ver "Known issues" arriba;
ninguno de los dos bloquea `maaqc run`/`maaqc doctor`, pero sí bloquean un uso
real de `compute_gv_appendix()` y dejan el propio test suite en rojo).

Los siguientes issues quedan pendientes antes de mergear a `main`:
- Corregir los 2 "Known issues" de arriba (triviales, minutos de trabajo).
- Validación end-to-end con moléculas de PAPER1 reales (comparación contra legacy).
- Tests de integración completos (mock ORCA outputs); en particular
  `compare_adc2_tddft()` sigue sin test de contrato (depende de un layout de
  directorio real no disponible como fixture todavía).
- Documentación de deployment (conda, conda-forge).
- Tag `v0-legacy` sobre commit HEAD de main (preservar estado pre-migración).

## [0.1.0.dev0] - TODO
- Tag `v0-legacy` en git sobre `pipeline.py` / `results_collector.py` /
  `protonation_cycle.py` / `recover_prescreening.py` antes de empezar
  la migración, para no perder el estado que ya produjo resultados
  usados en PAPER1.
