# Stages — inputs, outputs, estado de migración

| Stage    | requires  | Fuente legacy (pipeline.py)              | Estado      |
|----------|-----------|-------------------------------------------|-------------|
| crest    | —         | `stage_crest` (1152-1233)                  | migrado |
| censo    | crest     | `stage_censo_prescreen` (1438-1881), `stage_censo_multitemp` (2019-2107) | pendiente (stub, `NotImplementedError`) |
| opt      | crest     | `stage_opt` (2462-2600)                    | migrado |
| numfreq  | opt       | `stage_numfreq` (2601-2671)                | pendiente (stub, `NotImplementedError`) |
| sp       | numfreq\* | `stage_sp` (2672-2747), `stage_sp_nroot` (2881-2971) | migrado |
| tddft    | sp        | `stage_tddft` (2814-2880)                  | migrado |
| adc2     | tddft\*\* | `stage_adc2` (2972-3135)                   | migrado |

Fuente de verdad para esta tabla: `maaqc doctor`, que reporta el estado
real de cada stage (`[OK]` para migrados con lógica real, `[PENDIENTE]`
para stubs). Regenerar esta tabla desde ahí en vez de mantenerla a mano
para evitar que vuelva a divergir.

\* Con `skip_numfreq: true`, `sp` depende de `opt` directo — ver TODO de
diseño en `src/maaqc/stages/sp.py` (requires dinámico según config,
pendiente de resolver en `orchestrator/dag.py`).

\*\* Con `entry: sp-only` (equivalente al legacy `--start adc2`), `adc2`
ignora la cadena completa y solo exige `04_sp/` — ver TODO de diseño en
`src/maaqc/stages/adc2.py`.

## Orden de migración sugerido

1. `system/affinity.py` + `system/platform_detect.py` — autocontenidos,
   sin dependencia de parsing ORCA. ✅ hecho.
2. `io/parsing.py` — portar todos los parsers duplicados (ver el header
   de ese archivo para el mapeo exacto legacy -> maaqc). ✅ hecho.
3. `stages/opt.py` y `stages/sp.py` — los stages más simples, correr en
   paralelo con `pipeline.py` legacy sobre el mismo input y comparar
   outputs antes de confiar en la migración. ✅ hecho (junto con
   `crest.py`, `tddft.py` y `adc2.py`, ver tabla arriba).
4. `orchestrator/state.py` — generaliza `recover_prescreening.py`.
5. `stages/censo.py` — el más complejo (recovery, multitemp). ⬅ en curso.
6. `analysis/` completo, apoyado en `io/parsing.py` ya unificado.
7. `cli.py` real sobre el orchestrator ya migrado.

## Duplicación resuelta: `io/xyz.py`

`io/xyz.py` tenía 2 stubs (`NotImplementedError`) aunque las 3 funciones
correspondientes ya estaban migradas y en uso activo, viviendo
duplicadas en `system/_orca_shared.py` (`atoms_to_xyz_block`,
`read_all_conformers`, `extract_conf_ids_from_xyz` — esta última ni
siquiera declarada como stub). Resuelto: `io/xyz.py` es ahora la fuente
única de verdad; `_orca_shared.py` importa desde ahí. La firma legacy
`read_conformers(xyz_file, nconf)` se descartó formalmente (no se migra
tal cual) — el truncamiento a `nconf` ahora es responsabilidad de CENSO
(`select_nconf_auto`/`select_nconf_boltzmann`), no de la lectura de XYZ.
Ver docstring de `io/xyz.py` para el detalle de la decisión.

## Notas por stage

Ver los docstrings de cada archivo en `src/maaqc/stages/` — cada uno
documenta líneas exactas de origen en el legacy y TODOs de diseño
pendientes (no se documentan dos veces acá para evitar que diverjan).
