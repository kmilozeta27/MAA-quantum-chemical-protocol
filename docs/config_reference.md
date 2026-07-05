# Config reference — recipes YAML

Cada campo, su default y de dónde viene en el legacy. Ver
`src/maaqc/config/schema.py` para los modelos pydantic reales (fuente
de verdad); esta tabla es documentación derivada, no autoritativa —
si diverge del schema, el schema gana.

## Nivel raíz (`RecipeConfig`)

| Campo         | Default | Legacy equivalente                              |
|---------------|---------|--------------------------------------------------|
| `extends`     | `null`  | (nuevo — herencia de presets, no existía)         |
| `entry`       | `xyz`   | `--start` / `--expand-crest` / `--start adc2`     |
| `reference`   | `null`  | (nuevo — para `entry: calibrated`, ej. shinorine) |
| `skip_numfreq`| `false` | `--skipnumfreq`                                   |
| `charge`      | `null`  | `--charge` (obligatorio en la práctica)           |
| `multiplicity`| `1`     | `MULT` (constante, siempre 1)                     |

## `resources`

| Campo             | Default | Legacy (pipeline.py línea) |
|-------------------|---------|------------------------------|
| `nprocs_crest`    | 12      | `NPROCS_CREST` (180)         |
| `nprocs`          | 4       | `NPROCS` (181)               |
| `nprocs_numfreq`  | 8       | `NPROCS_NUMFREQ` (182)       |
| `nprocs_tddft`    | 6       | `NPROCS_TDDFT` (183)         |
| `nprocs_adc2`     | 4       | `NPROCS_ADC2` (186)          |
| `maxcore`         | 3500    | `MAXCORE` (184)              |
| `maxcore_numfreq` | 2800    | `MAXCORE_NUMFREQ` (185)      |
| `maxcore_tddft`   | 2000    | `MAXCORE_TDDFT` (187, sic — ver nota) |
| `maxcore_adc2`    | 5500    | `MAXCORE_ADC2` (187)         |
| `opt_batch_size`  | 3       | `OPT_BATCH_SIZE` (189)       |
| `sp_batch_size`   | 3       | `SP_BATCH_SIZE` (190)        |
| `tddft_batch_size`| 2       | `TDDFT_BATCH_SIZE` (191)     |
| `cores_total`     | `null`  | `--cores`                    |
| `max_ram_mb`      | 28000   | comentario de `--cores` en argparse (línea ~3992) |

Nota: `maxcore_tddft`/`maxcore_adc2` — verificar contra pipeline.py real
al migrar `resources.py`, la numeración de línea exacta de estas dos
constantes puede haberse desplazado entre ediciones del archivo.

## `censo`

| Campo                 | Default            | Legacy |
|-----------------------|--------------------|--------|
| `temperature`         | 298.15             | `--temperature` |
| `multitemp`           | `false`            | `--multitemp` |
| `trange`              | `[274.0, 298.0, 24]` | `--trange` (default cuando `--multitemp` está activo) |
| `boltzmann_threshold` | `null`             | `--boltzmann-threshold` (null = `select_nconf_auto`) |
| `ncon_min` / `ncon_max` | 3 / 20           | límites hardcodeados en `select_nconf_auto` (pipeline.py:2277-2414) |
| `timeout_hours`       | 24.0               | `CENSO_TIMEOUT_HOURS` (175) |

## `tddft`

| Campo        | Default      | Legacy |
|--------------|--------------|--------|
| `functional` | `wb97x-d3`   | hardcodeado en templates/ (confirmar al migrar `stages/tddft.py`) |
| `states`     | 10           | hardcodeado en templates/ (confirmar) |

## `adc2`

| Campo                 | Default        | Legacy |
|-----------------------|----------------|--------|
| `mandatory`           | `false`        | `--adc2` (opt-in explícito) |
| `basis`               | `def2-TZVP`    | hardcodeado en `stage_adc2` (confirmar línea exacta al migrar) |
| `exp_lambda_fallback` | `{}`           | `EXP_LAMBDA_FALLBACK` (222-241) |
| `molecules_yml`       | `null`         | `_MOLECULES_YML` (259) |

## `nroot`

| Campo     | Default | Legacy |
|-----------|---------|--------|
| `enabled` | `false` | `--nroot` |
| `top_n`   | 5       | `--nroot-top` |
