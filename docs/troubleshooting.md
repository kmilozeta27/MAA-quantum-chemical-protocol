# Troubleshooting

Conocimiento acumulado del pipeline legacy, sistematizado. Completar a
medida que se migra cada stage — no reescribir desde cero, portar el
conocimiento real ya validado.

## CENSO 3.x

- **CENSO nunca reusa cálculos previos.** No tiene mecanismo de "skip si
  el .out ya existe": si el proceso se interrumpe, CENSO no escribe
  `0_PRESCREENING.json` (solo lo hace al final, cuando TODOS los
  confórmeros terminan en la misma ejecución). El único punto de reuso
  es `--reload archivo.json`, que exige que TODOS los confórmeros de la
  ensemble actual estén presentes en ese JSON — no acepta reload
  parcial.
  → En maaqc esto se resuelve en `orchestrator/state.py::rebuild_censo_state()`,
    generalizando lo que hoy hace `recover_prescreening.py`.

- Ver `pipeline.py:_censo_prescreen_diagnose_dir` (1276-1341) para el
  diagnóstico de por qué un directorio CENSO quedó incompleto antes de
  intentar recovery.

## CPU híbrida (ThinkPad T14, Intel Core Ultra 7, P-core/E-core)

- `hwloc` (usado internamente por el `mpirun` de ORCA en la fase
  LEANSCF) lee la topología directo de `/sys/.../cluster` sin filtrar
  por el cpuset heredado de `taskset`. En una CPU híbrida esto produce
  cpusets de cluster L1d que no son subconjunto coherente del cpuset
  restringido ("intersection without inclusion"), y aunque hwloc dice
  que ignora el error y continúa, el reparto de trabajo entre ranks MPI
  en LEANSCF queda corrupto y aborta justo al iniciar la formación de
  Fock.
  → Fix: forzar binding por núcleo en vez de por topología vía
    `OMPI_MCA_hwloc_base_binding_policy=none` +
    `OMPI_MCA_rmaps_base_mapping_policy=core` +
    `OMPI_MCA_hwloc_base_use_hwthreads_as_cpus=1`
    (ver `system/affinity.py::build_omp_env`, migrado desde
    `pipeline.py` líneas 197-233). **No revertir este fix al migrar** —
    es producto de debugging real contra un fallo silencioso y tardío.

## Binarios no encontrados a mitad de un run largo

- Antes: `CENSO_BIN`/`CREST_BIN` eran solo el nombre del comando, así
  que su existencia dependía 100% de que el PATH del proceso que lanzó
  `pipeline.py` tuviera el env conda correcto activado. Si se invocaba
  desde un shell/cron/nohup sin `conda activate sciencebitch`, el fallo
  solo se descubría al llamar `subprocess.run()` — a veces horas
  después (p. ej. tras 70 min de CREST).
  → Fix: `_resolve_bin()` busca primero en el PATH actual
    (`shutil.which`) y si falla, en `<prefix del env>/bin/<nombre>` vía
    `conda env list --json`, sin necesitar que el env esté activado.
    → En maaqc: `system/platform_detect.py::resolve_bin` +
      `maaqc doctor` para chequear esto ANTES de lanzar CREST.

## CENSO `--reload` y matching por energía

- Cuando los tokens de string esperados en `1_SCREENING.out` no
  coinciden (cambios de formato entre versiones 3.x de CENSO),
  `pipeline.py:_censo_prescreen_recover` (1342-1437) cae a matching por
  energía DFT en vez de parseo de string. Esta es la generalización que
  `orchestrator/state.py` debe absorber como comportamiento normal, no
  como fallback de emergencia.

## GoodVibes

- Dos resoluciones independientes del intérprete Python de GoodVibes
  (`results_collector.py:_resolve_gv_python` y
  `protonation_cycle.py:_resolve_gv_python`) — revisar si son idénticas
  antes de migrar a `analysis/goodvibes.py`, puede ser otra duplicación
  del mismo tipo que los parsers de ORCA.

## Disco

- ADC2 necesita ~400GB de disco. Flujo real: respaldar todo en disco
  externo y borrar los archivos pesados del computador antes de correr
  ADC2, dejando solo `04_sp/` disponible localmente. Esto es lo que
  justifica que `entry: sp-only` (`--start adc2` legacy) ignore por
  completo `crest/censo/opt/numfreq/tddft` aunque no existan en disco.
