# Stages — inputs, outputs, estado de migración

| Stage    | requires  | Fuente legacy (pipeline.py)              | Estado      |
|----------|-----------|-------------------------------------------|-------------|
| crest    | —         | `stage_crest` (1152-1233)                  | no migrado |
| censo    | crest     | `stage_censo_prescreen` (1438-1881), `stage_censo_multitemp` (2019-2107) | no migrado |
| opt      | censo     | `stage_opt` (2462-2600)                    | no migrado |
| numfreq  | opt       | `stage_numfreq` (2601-2671)                | no migrado |
| sp       | numfreq\* | `stage_sp` (2672-2747), `stage_sp_nroot` (2881-2971) | no migrado |
| tddft    | sp        | `stage_tddft` (2814-2880)                  | no migrado |
| adc2     | tddft\*\* | `stage_adc2` (2972-3135)                   | no migrado |

\* Con `skip_numfreq: true`, `sp` depende de `opt` directo — ver TODO de
diseño en `src/maaqc/stages/sp.py` (requires dinámico según config,
pendiente de resolver en `orchestrator/dag.py`).

\*\* Con `entry: sp-only` (equivalente al legacy `--start adc2`), `adc2`
ignora la cadena completa y solo exige `04_sp/` — ver TODO de diseño en
`src/maaqc/stages/adc2.py`.

## Orden de migración sugerido

1. `system/affinity.py` + `system/platform_detect.py` — autocontenidos,
   sin dependencia de parsing ORCA.
2. `io/parsing.py` — portar todos los parsers duplicados (ver el header
   de ese archivo para el mapeo exacto legacy -> maaqc).
3. `stages/opt.py` y `stages/sp.py` — los stages más simples, correr en
   paralelo con `pipeline.py` legacy sobre el mismo input y comparar
   outputs antes de confiar en la migración.
4. `orchestrator/state.py` — generaliza `recover_prescreening.py`.
5. `stages/censo.py` — el más complejo (recovery, multitemp).
6. `analysis/` completo, apoyado en `io/parsing.py` ya unificado.
7. `cli.py` real sobre el orchestrator ya migrado.

## Notas por stage

Ver los docstrings de cada archivo en `src/maaqc/stages/` — cada uno
documenta líneas exactas de origen en el legacy y TODOs de diseño
pendientes (no se documentan dos veces acá para evitar que diverjan).
