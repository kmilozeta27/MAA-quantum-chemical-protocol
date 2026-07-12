# Parche: coordinación de slots de CPU en OrcaProc._sp()/opt() (CENSO)

## Problema que resuelve

CENSO (`src/censo/processing/orca_processor.py`) llama a ORCA vía
`self._make_call(call, outputpath, jobdir)` **sin pasar `env=`** en `_sp()`
(línea ~508) y en `opt()` (línea ~655). Eso hace que `_make_call()` (en
`processor.py`) caiga siempre a `env=env or ENVIRON`, donde `ENVIRON` es
`os.environ.copy()` congelado UNA sola vez al importar `params.py` — el
mismo diccionario para TODAS las llamadas paralelas, sin importar cuántos
"workers"/hilos de Dask las disparen.

Confirmado con `ps --forest` real: cuando CENSO corre 2+ SP en paralelo
(cada uno invocando su propio `mpirun -np N`), los `mpirun` heredan el
mismo `MAAQC_PE_LIST` — si `pipeline.py` lo fijara a nivel del proceso
`censo` padre — resultando en pe-lists idénticos para invocaciones
concurrentes que no coordinan entre sí. Esto ya fue mitigado del lado de
`pipeline.py` con `set_pe_list=False` (ver historial en
`build_omp_env()`), pero eso deja el binding librado al kernel — sin
garantía de reparto limpio entre 2 slots de 4 cores dedicados.

Este parche resuelve el caso concreto: **2 slots fijos de 4 cores cada
uno, sin que los SP paralelos de CENSO se disputen cores entre sí**,
asignando cada conformero de forma determinística a un slot.

## Mecanismo

1. `pipeline.py` sigue confinando el árbol completo de CENSO con
   `taskset -c <pool completo>` (sin cambios ahí).
2. `pipeline.py` fija una variable de entorno nueva, `MAAQC_CENSO_SLOTS`,
   con el pool particionado en N grupos de M PUs cada uno, separados por
   `|`. Ejemplo para 2 slots de 4 PUs: `"0,1,4,5|2,3,6,7"`.
3. Este parche, en `orca_processor.py`, intercepta `_sp()` y `opt()`:
   antes de llamar a `_make_call`, arma un `env` explícito (copia de
   `ENVIRON` + el `MAAQC_PE_LIST` correspondiente al slot que le toca a
   `job.conf.name`), y se lo pasa a `_make_call(..., env=env)`.
4. El slot se deriva con un hash estable de `job.conf.name` módulo N —
   determinístico entre corridas (mismo conformero → mismo slot siempre),
   sin necesitar coordinación en tiempo real entre hilos de Dask (que
   sería mucho más invasivo de implementar de forma segura).

## Por qué esto es suficiente sin coordinación en tiempo real

Con N=2 slots y una distribución razonable de conformeros (docenas típico
en un ensemble CREST/CENSO), un hash módulo 2 reparte aproximadamente
mitad/mitad. No es un scheduler dinámico — dos conformeros del mismo slot
que casualmente corran al mismo tiempo sí competirían por esas 4 PUs —
pero esto es enormemente mejor que el estado actual (pe-list idéntico
para TODOS, o binding libre del kernel sin ninguna partición). Si en el
futuro se necesita perfectamente 0 colisión también en ese caso raro,
haría falta un lock file real (ver sección "Extensión futura" abajo).

## Aplicación del parche

Backup primero, SIEMPRE:

```bash
cp /home/kmilozeta/CENSO/src/censo/processing/orca_processor.py \
   /home/kmilozeta/CENSO/src/censo/processing/orca_processor.py.bak-pre-slots-patch
```

Aplicar los cambios de `orca_processor.py.patch` (ver ese archivo, diff
unificado) con:

```bash
cd /home/kmilozeta/CENSO
patch -p1 --dry-run < /ruta/a/orca_processor.py.patch   # simular primero
patch -p1 < /ruta/a/orca_processor.py.patch             # aplicar de verdad
```

O aplicar manualmente los 2 bloques marcados `# MAAQC-PATCH` en el archivo
`orca_processor.py.patched` (versión completa ya parcheada, para copiar/
comparar con diff visual si preferís no usar `patch`).

## IMPORTANTE — este parche se pierde en cualquier reinstalación de CENSO

Como está aplicado directamente sobre el paquete instalado en
`sciencebitch` (no es un fork ni un paquete separado), CUALQUIER
`pip install -U censo` / reinstalación / `conda update` en ese entorno
va a sobreescribir `orca_processor.py` y perder el parche silenciosamente
— sin error, sin aviso. `pipeline.py doctor` debería (pendiente,
Paso siguiente) verificar al arrancar que el parche sigue presente antes
de lanzar CENSO, y advertir si no lo está.

## Extensión futura (si 2 slots no alcanza / se necesita 0 colisión real)

Reemplazar el hash estático por un lock file real en
`/tmp/maaqc_censo_slots.lock` con `fcntl.flock`, donde cada llamada a
`_sp()`/`opt()` adquiere el primer slot libre y lo libera al terminar
(try/finally alrededor de `_make_call`). Más robusto pero más invasivo —
no implementado en este parche porque no fue pedido y añade superficie de
fallo (locks huérfanos si el proceso muere abruptamente).
