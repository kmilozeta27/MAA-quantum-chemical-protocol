#!/usr/bin/env python3
"""
pipeline.py — Pipeline ORCA automatizado para cálculos conformacionales
Flujo: CREST → CENSO (prescreen) → Opt → NumFreq → SP → TDDFT → [ADC2]

El pipeline siempre empieza en CREST. CENSO pre-screening es obligatorio
entre CREST y OPT. El número de confórmeros que pasan a OPT se determina
automáticamente a partir del panorama conformacional post-CENSO (mín 3, máx 20).

SP (fase gas) es la etapa OBLIGATORIA: el pipeline solo se da por
FINALIZADO cuando al menos un confórmero completa SP. ADC2 es OPCIONAL
y NO corre por defecto — debe pedirse explícitamente con --adc2 (para
incluirla al final del pipeline normal) o de forma standalone con
--start adc2 (ver nota más abajo).

Carga automática por nombre de archivo:
    --charge ahora es OPCIONAL. Si se omite, se infiere del sufijo del
    nombre del .xyz según la convención de estados de protonación:
        ...A → +1   ...B → -1   ...C → -2   ...D → 0   ...R → 0
    (p.ej. "SB.xyz" -> carga -1, "P334R.xyz" -> carga 0). Si el nombre
    no sigue esa convención, se avisa y se usa 0 por defecto — pasa
    --charge explícito para forzar otro valor o silenciar el aviso.

Uso (modo normal):
    python pipeline.py --mol molecules/mol_A.xyz              # carga +1 inferida por sufijo "A"
    python pipeline.py --mol molecules/mol_A.xyz --charge 1   # equivalente, forma explícita
    python pipeline.py --mol molecules/mol_A.xyz --charge 0 --start opt
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --start sp
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --skipnumfreq   # OPT → SP → TDDFT directo
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --nroot         # SP+NRoots1 top-5 por E_OPT
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --nroot --nroot-top 3  # top-3
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --adc2         # pipeline normal + ADC2 al final (opcional)
    python pipeline.py --status                  # ver estado de todos los runs

Uso (temperatura personalizada):
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --temperature 298.15
    # → CENSO corre a 298.15 K (default cuando no se usa --multitemp)

Uso (modo multi-temperatura, --multitemp):
    # Análisis Boltzmann a múltiples temperaturas. CENSO corre N veces —
    # una por temperatura — guardando cada ensemble en un subdirectorio
    # separado (01b_censo_prescreen/T<temp>K/). El pipeline selecciona los
    # confórmeros para OPT usando el ensemble de --temperature (default 298.15 K).
    # --trange es OPCIONAL: default = 274 298 24 (24 puntos, 274 K … 298 K).
    #
    # Equivale al comando que usaría CENSO directamente:
    #   censo -i A.xyz --charge -1 --multitemp --trange 274 298 24
    #         --temperature 298.15 --solvent h2o
    #         (prog/func se leen de ~/.censo2rc; --prog eliminado en CENSO 2.x)
    #
    # Con todos los defaults (equivalente al comando de arriba):
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --multitemp
    #
    # Con trange explícito:
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 \\
        --multitemp --trange 274 298 24 --temperature 298.15

Nota: --start crest y --start opt son equivalentes al inicio (ambos requieren
      que CREST+CENSO corran). Solo usar --start opt/numfreq/sp/tddft/adc2
      para retomar un pipeline con CREST+CENSO ya completados.

Nota (SP obligatorio / ADC2 opcional): en el pipeline normal, SP (fase gas)
      es la etapa mínima requerida — el pipeline se aborta si SP no completa
      ningún confórmero, y solo se declara "PIPELINE FINALIZADO" una vez que
      SP terminó. TDDFT corre después de SP salvo que --start sea posterior.
      ADC2 en cambio ya NO corre automáticamente: debe pedirse explícitamente
      con --adc2 para incluirla al final del pipeline normal.

Nota (--start adc2): caso especial, NO retoma la cadena normal. Solo requiere
      que exista 04_sp/ (con sus .out y opt_opt.xyz por confórmero) — ignora
      por completo 01_crest/, 02_opt/, 03_numfreq/ y 05_tddft/, incluso si no
      existen en disco. Identifica automáticamente el confórmero de menor
      energía SP y corre ADC2 solo sobre ese. Pensado para el flujo real de
      trabajo: como ADC2 necesita ~400GB de disco, Kmilo respalda todo en un
      disco externo y borra los archivos pesados del computador antes de
      correrlo, dejando únicamente 04_sp/ disponible.
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --start adc2

────────────────────────────────────────────────────────────────────────
CHANGELOG — historial de bugs y parches relevantes
────────────────────────────────────────────────────────────────────────
Bitácora de problemas ya resueltos que motivaron partes no obvias del
código. Se mantiene aquí, consolidada, para no interrumpir la lectura
de la lógica activa; cada función afectada deja solo una referencia
corta a la entrada correspondiente.

[MPI-1] pe-list compartido entre sub-jobs paralelos de CENSO (jul 2026)
    Con bind-to hwthread ya corregido (ver MPI-2), CENSO seguía saturando
    solo 2 cores físicos (4/12 hilos) al correr sus SP internos en
    paralelo, con throttling térmico medido (102°C). Causa: CENSO congela
    ENVIRON una sola vez al lanzar el proceso `censo` padre, y cada
    `mpirun` interno que CENSO invoca hereda ESE MISMO entorno — los N
    sub-jobs paralelos que CENSO lanza colapsaban sobre el mismo
    MAAQC_PE_LIST fijo, todos pidiendo los mismos primeros cores.
    Mitigación aplicada: build_omp_env(set_pe_list=False) para los
    llamados desde el camino de CENSO — taskset sigue confinando el
    árbol completo al pool asignado, pero MAAQC_PE_LIST no se fija, así
    que el kernel reparte libremente dentro del pool (no es pinning
    óptimo, pero es estrictamente mejor que la colisión forzada previa).
    Solución completa instalada además: MAAQC_CENSO_SLOTS + un parche
    local en CENSO (processing/orca_processor.py del paquete instalado,
    fuera de este repo — se pierde en cualquier reinstalación/update de
    CENSO) que particiona el pool en N slots y deriva el slot de cada
    confórmero con un hash estable de su nombre. Si los SP de CENSO
    vuelven a mostrar colisión parcial, lo primero a revisar es si ese
    parche sigue presente.

[MPI-2] namespace viejo vs PRTE, y pe-list como índice de CORE no de PU
    (jul 2026)
    Las variables OMPI_MCA_hwloc_base_binding_policy /
    rmaps_base_mapping_policy son el namespace viejo (pre-5.0, runtime
    ORTE); Open MPI ≥5.0 migró a PRTE y las renombró bajo el prefijo
    PRTE_MCA_ — pero fijar la política de mapping como "default" del DVM
    tampoco sirve (PRTE la rechaza con "unsupported-default-policy"). La
    única vía que PRTE acepta para un pe-list explícito es --map-by en
    la línea de comandos de mpirun/prterun, que ORCA invoca internamente
    sin exponer flags. Solución: un wrapper de `mpirun` (symlink en PATH
    antes que el real, ver scripts/mpirun-wrapper.sh) que lee
    MAAQC_PE_LIST y le inyecta `--map-by pe-list=<lista>:ORDERED:hwtcpus
    --bind-to hwthread`.
    Un segundo bug independiente: ORCA seguía usando solo 6/12 cores en
    batches que deberían repartirse en 12. Diagnóstico con `mpirun
    --display-map --do-not-launch` + `lstopo`: PRTE interpreta CUALQUIER
    lista numérica pasada a pe-list/rankfile como ÍNDICE DE CORE, nunca
    como PU/hilo de hardware (confirmado: `Cpu Type: CORE` en el mapa) —
    independiente del flag --bind-to. taskset -c, en cambio, sí usa PU
    directa (confirmado vía /proc/<pid>/status). Como MAAQC_PE_LIST se
    armaba con la misma lista de PUs que taskset, cada PU terminaba
    resolviendo a un core físico distinto de lo esperado (offset +2
    constante en la zona E-core; agrupamiento inesperado en la zona
    P-core con SMT). Fix: MAAQC_PE_LIST se arma traduciendo cada PU a su
    índice de Core vía la tabla PU_TO_CORE (derivada de `lstopo`);
    taskset sigue recibiendo la lista de PUs cruda — esa ruta nunca tuvo
    el bug.

[CENSO-BUDGET] OMP_NUM_THREADS al tamaño del pool completo, no del
    sub-job (ago 2026)
    Todos los call-sites de CENSO armaban su env con
    OMP_NUM_THREADS = NPROCS_CREST (el pool completo, ej. 12), pero
    CENSO decide internamente cuántos sub-jobs correr en paralelo según
    maxcores // omp_min (con --omp-min default 4 y --maxcores 12, decide
    3 sub-jobs a la vez). Con OMP_NUM_THREADS=12 en vez de ~4 por
    sub-job, los 3 sub-jobs concurrentes pedían 12 hilos CADA UNO → 36
    hilos compitiendo por 12 cores confinados por taskset (3x
    sobre-suscripción real). Fix: fijar --omp-min explícitamente y usar
    ESE mismo valor — nunca maxcores — como OMP_NUM_THREADS
    (censo_threads_per_job()).
    Bug relacionado, hallado al revisar por qué la sobre-suscripción
    persistía con 2+ instancias del pipeline corriendo a la vez:
    _acquire_cpu_pool() puede devolver un pool MÁS CORTO que lo pedido
    (incluso vacío) si otra instancia ya reservó PUs, pero los
    call-sites de CENSO seguían usando el valor NOMINAL pedido (no
    len(cpu_pool) real) tanto para --maxcores como para
    censo_threads_per_job(). Fix: _reconcile_censo_budget() ajusta
    (maxcores, omp_min) al tamaño real de cpu_pool después de
    reservarlo; si el pool queda vacío, aborta el lanzamiento de CENSO
    en vez de correrlo sin confinar.

[CENSO-LOCK] relanzamientos duplicados sobre el mismo censo_dir
    (10 jul 2026, run A/MGLYVAL)
    stage_censo_prescreen() decidía si "recuperar" un run mirando solo
    el contenido de censo_dir (carpetas CONF* con .out parseables), pero
    eso no indica si CENSO sigue corriendo — las carpetas CONF*
    sobreviven después de que CENSO termina un part y avanza al
    siguiente. Con el proceso original todavía vivo y ya en screening,
    una segunda invocación del pipeline vio "192 CONF dirs completas" y
    disparó un --reload que relanzó un SEGUNDO proceso CENSO desde
    CONF1, descartando el progreso de screening del primero — y lo mismo
    volvió a pasar con un tercero: 3 CENSOs sucesivos rehaciendo
    prescreening desde cero sin que ninguno llegara a terminar
    screening. Fix: un lockfile propio del pipeline
    (.censo_pipeline.lock, no confundir con censo.log) que contiene el
    PID real del proceso CENSO mientras corre; stage_censo_prescreen()
    se niega a tocar censo_dir si ya hay un CENSO vivo ahí, sin importar
    el contenido de las carpetas CONF*.

[BIN-RESOLVE] CENSO_BIN/CREST_BIN dependían del PATH del proceso lanzador
    Antes, CENSO_BIN/CREST_BIN eran solo el nombre del comando, así que
    su existencia dependía 100% de que el shell/cron/nohup que lanzó
    pipeline.py ya tuviera 'conda activate sciencebitch' hecho — si no,
    el fallo solo se descubría al llamar subprocess.run(), a veces
    horas después (p.ej. tras 70 min de CREST). Fix: _resolve_bin()
    busca primero en el PATH actual (shutil.which) y si no lo encuentra
    cae a <prefix del env conda>/bin/<nombre> vía `conda env list
    --json`, sin necesitar que el env esté activado. check_required_bins()
    además verifica ambos binarios ANTES de lanzar CREST, para no
    descubrir un "censo: command not found" después de gastar cómputo.
────────────────────────────────────────────────────────────────────────
"""

import argparse
import datetime
import fcntl
import json
import logging
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from ase.io import read, write

# ─────────────────────────────────────────────
#  CONFIGURACIÓN GLOBAL
# ─────────────────────────────────────────────

def _resolve_bin(name: str, conda_env: str = "sciencebitch") -> Path:
    """
    Resuelve un ejecutable a una ruta absoluta de forma robusta.

    Antes, CENSO_BIN/CREST_BIN eran solo el nombre del comando (p.ej. Path("censo")),
    así que su existencia dependía 100% de que el PATH del proceso que lanzó
    pipeline.py ya tuviera el env conda correcto activado. Si el script se
    invocaba desde un shell/cron/nohup sin 'conda activate sciencebitch', el
    fallo solo se descubría al llamar subprocess.run() — a veces horas después
    (p.ej. tras 70 min de CREST), con un FileNotFoundError tardío.

    Busca en orden:
      1. PATH actual del proceso (shutil.which) — caso normal, env ya activado.
      2. <prefix de conda_env>/bin/<name> — vía `conda env list --json`,
         sin necesitar que el env esté activado.
    Si no se encuentra en ningún lado, devuelve Path(name) tal cual (para no
    romper compatibilidad); el caller debe loggear/abortar temprano con
    check_required_bins().
    """
    found = shutil.which(name)
    if found:
        return Path(found)
    try:
        result = subprocess.run(
            ["conda", "env", "list", "--json"],
            capture_output=True, text=True, timeout=15,
        )
        if result.returncode == 0:
            envs = json.loads(result.stdout).get("envs", [])
            for env_path in envs:
                if Path(env_path).name == conda_env:
                    candidate = Path(env_path) / "bin" / name
                    if candidate.is_file() and os.access(candidate, os.X_OK):
                        return candidate
    except Exception:
        pass
    return Path(name)


def check_required_bins(log: logging.Logger) -> None:
    """
    Verifica al inicio del pipeline (antes de lanzar CREST) que los binarios
    requeridos existan y sean ejecutables, con mensaje de diagnóstico claro.
    Evita descubrir un 'censo: command not found' después de 70+ min de CREST.
    """
    missing = []
    for bin_name, bin_path in (("CREST", CREST_BIN), ("CENSO", CENSO_BIN)):
        resolved = shutil.which(str(bin_path)) or (
            str(bin_path) if bin_path.is_absolute() and bin_path.is_file() else None
        )
        if resolved is None:
            missing.append((bin_name, bin_path))
        else:
            log.debug(f"[CHECK-BIN] {bin_name} → {resolved}")
    if missing:
        for bin_name, bin_path in missing:
            log.error(
                f"[CHECK-BIN]  ✗  {bin_name} no encontrado ({bin_path}). "
                f"¿'conda activate sciencebitch' antes de lanzar pipeline.py? "
                f"¿está instalado en ese env (pip/conda list | grep -i {bin_name.lower()})?"
            )
        raise SystemExit(
            "ERROR: faltan ejecutables requeridos: "
            + ", ".join(n for n, _ in missing)
            + ". Pipeline abortado ANTES de CREST para no desperdiciar horas de cómputo."
        )


ORCA_BIN   = Path(os.environ.get("PIPELINE_ORCA_BIN",
                                  "/home/kmilozeta/orca_6_1_1_avx2/orca"))
CREST_BIN  = _resolve_bin("crest")  # resuelto vía PATH o conda env sciencebitch
CENSO_BIN  = _resolve_bin("censo")  # resuelto vía PATH o conda env sciencebitch


def _find_paper_root(start: Optional[Path] = None) -> Path:
    """
    Resuelve la carpeta raíz del proyecto (equivalente a lo que antes era
    PAPER_ROOT hardcodeado a /home/kmilozeta/PAPER1) sin asumir ninguna
    ruta fija. Orden de resolución:

      1. Variable de entorno PIPELINE_ROOT, si está seteada.
      2. Subir desde `start` (o desde el cwd si no se pasa nada) buscando
         un marcador de raíz de proyecto: un directorio que contenga
         alguno de README.md, LICENSE, tasks.json o una carpeta scripts/.
      3. Si no se encuentra ningún marcador, se usa el cwd tal cual
         (fallback razonable: correr el pipeline desde la raíz del
         proyecto sigue funcionando sin configurar nada).

    Esto reemplaza la ruta fija PAPER_ROOT = Path("/home/kmilozeta/PAPER1"):
    ahora el proyecto puede vivir en cualquier máquina/ruta sin editar el
    script, y --status recorre la raíz real detectada.
    """
    env_root = os.environ.get("PIPELINE_ROOT")
    if env_root:
        return Path(env_root).expanduser().resolve()

    here = (start or Path.cwd()).resolve()
    markers = {"README.md", "LICENSE", "tasks.json", "scripts"}
    for candidate in [here, *here.parents]:
        if any((candidate / m).exists() for m in markers):
            return candidate
    return here


PAPER_ROOT = _find_paper_root()   # carpeta madre del proyecto (auto-detectada)
WORK_ROOT  = PAPER_ROOT                        # alias, usado solo por print_status()
TMPL_DIR   = Path("templates")
MULT       = 1                      # multiplicidad siempre 1
CREST_TIMEOUT_HOURS = 24.0          # timeout CREST (horas)
CENSO_TIMEOUT_HOURS = 24.0          # timeout CENSO pre-screening (horas)

# ── Convención de carga por sufijo de estado de protonación ────────────
# En todos los proyectos de la serie MAA, el nombre del archivo .xyz
# (su stem, sin importar el prefijo de la molécula: A, SA, PA, mol_A...)
# termina en una letra que codifica el estado de protonación:
#
#   ...A  → catión      (+1)
#   ...B  → anión simple (-1)
#   ...C  → dianión      (-2)
#   ...D  → neutro       ( 0)
#   ...R  → neutro       ( 0)   (forma de referencia)
#
# Esto permite inferir --charge automáticamente a partir del nombre del
# archivo, sin tener que pasarlo a mano en cada corrida. --charge sigue
# aceptándose como override explícito si se quiere forzar otra cosa.
STATE_CHARGE_MAP: dict[str, int] = {
    "A": 1,
    "B": -1,
    "C": -2,
    "D": 0,
    "R": 0,
}


def infer_charge_from_name(mol_name: str) -> Optional[int]:
    """
    Infiere la carga formal a partir del nombre de la molécula (stem del
    .xyz), usando STATE_CHARGE_MAP sobre la última letra del nombre.

    Reglas:
      - Se mira solo la ÚLTIMA letra alfabética del stem (p.ej. "SB" → 'B',
        "mol_A" → 'A', "P334R" → 'R', "MR" → 'R'). No importa el prefijo.
      - Si el nombre no termina en una letra reconocida en
        STATE_CHARGE_MAP (p.ej. termina en dígito, o es algo como
        "MGLYVAL" sin sufijo de estado), devuelve None: no se puede
        inferir con confianza y el caller debe pedir --charge explícito.
      - No es case-sensitive: "b" y "B" son equivalentes.

    Devuelve la carga inferida (int) o None si no se pudo inferir.
    """
    if not mol_name:
        return None
    last_char = mol_name[-1].upper()
    if not last_char.isalpha():
        return None
    return STATE_CHARGE_MAP.get(last_char)


def resolve_charge(mol_name: str, explicit_charge: Optional[int],
                    log: Optional[logging.Logger] = None) -> int:
    """
    Resuelve la carga a usar para `mol_name`, combinando lo que el usuario
    pasó explícitamente por --charge (si lo pasó) con la inferencia
    automática por sufijo (infer_charge_from_name).

    Prioridad:
      1. Si el usuario pasó --charge explícitamente (explicit_charge no es
         None), se usa ese valor tal cual — es un override manual y
         siempre gana. Si además hay una carga inferida y NO coincide,
         se loggea un warning (probable error de tipeo) pero se respeta
         igual lo que pidió el usuario.
      2. Si no se pasó --charge, se usa la carga inferida por sufijo.
      3. Si no se pasó --charge y no se pudo inferir nada, se cae al
         default histórico (0) y se avisa por log — igual que el
         comportamiento previo del script (default=0 en argparse).
    """
    inferred = infer_charge_from_name(mol_name)

    if explicit_charge is not None:
        if inferred is not None and inferred != explicit_charge:
            msg = (f"--charge {explicit_charge} fue especificado manualmente, pero "
                   f"el nombre '{mol_name}' sugiere carga {inferred:+d} por su "
                   f"sufijo de estado de protonación. Se respeta --charge "
                   f"{explicit_charge} (override manual) — verifica que sea "
                   f"intencional.")
            log.warning(msg) if log else print(f"AVISO: {msg}", file=sys.stderr)
        return explicit_charge

    if inferred is not None:
        msg = f"Carga inferida automáticamente para '{mol_name}': {inferred:+d} (sufijo '{mol_name[-1].upper()}')"
        log.info(msg) if log else print(msg)
        return inferred

    msg = (f"No se pudo inferir la carga de '{mol_name}' por su nombre "
           f"(no termina en A/B/C/D/R). Usando default 0 — pasa --charge "
           f"explícitamente si no es correcto.")
    log.warning(msg) if log else print(f"AVISO: {msg}", file=sys.stderr)
    return 0

# Ver [CENSO-BUDGET] en el CHANGELOG: OMP_NUM_THREADS debe reflejar los
# hilos que le tocan a CADA sub-job paralelo de CENSO, no el pool
# completo — CENSO congela ENVIRON una sola vez y lo hereda tal cual en
# cada sub-job concurrente que lanza.
CENSO_OMP_MIN = 4    # mismo valor que se pasa a --omp-min en el CLI de CENSO


def censo_threads_per_job(maxcores: int, omp_min: int = CENSO_OMP_MIN) -> int:
    """
    Replica (aproximada) de cómo CENSO reparte maxcores entre sub-jobs
    paralelos dado --omp-min: n_parallel = maxcores // omp_min sub-jobs
    corriendo a la vez, cada uno con maxcores // n_parallel hilos OpenMP.
    Este es el valor que hay que usar como OMP_NUM_THREADS del env que
    hereda TODO el árbol de CENSO -- no maxcores -- para que cada sub-job
    concurrente pida la cantidad de hilos que efectivamente le corresponde
    y no el pool completo (ver comentario junto a CENSO_OMP_MIN arriba).
    """
    n_parallel = max(1, maxcores // max(omp_min, 1))
    return max(1, maxcores // n_parallel)


def censo_n_parallel_jobs(maxcores: int, omp_min: int = CENSO_OMP_MIN) -> int:
    """
    Cuántos sub-jobs corre CENSO en paralelo dado maxcores/omp_min (ver
    get_cluster() en censo/parallel.py: threads_per_worker = ncores //
    ompmin, y Dask solo admite tantas tareas concurrentes como el
    presupuesto de recursos "CPU" alcance). Usado para dimensionar
    MAAQC_CENSO_SLOTS -- necesitamos exactamente un slot de CPUs por
    sub-job concurrente posible, ni más ni menos.
    """
    return max(1, maxcores // max(omp_min, 1))


# Ver [CENSO-BUDGET] en el CHANGELOG (segunda mitad): reconcilia
# (--maxcores, --omp-min) contra el tamaño REAL del cpu_pool reservado —
# que puede ser más chico que lo pedido, incluso vacío, si otra instancia
# de pipeline.py ya tomó PUs — en vez de asumir que siempre coincide.
def _reconcile_censo_budget(requested_cores: int, omp_min: int,
                             cpu_pool: list[int], log: logging.Logger,
                             tag: str) -> tuple[int, int]:
    """
    Ajusta (effective_maxcores, effective_omp_min) al tamaño REAL de
    cpu_pool (devuelto por _acquire_cpu_pool()), que puede ser menor que
    requested_cores si otra instancia de pipeline.py ya reservó PUs.

    Devuelve (0, omp_min) si cpu_pool quedó vacío — el caller debe
    interpretar maxcores == 0 como señal de abortar el lanzamiento de
    CENSO (liberando el lock) en vez de correrlo sin taskset.
    """
    n_reserved = len(cpu_pool)
    if n_reserved == 0:
        log.error(
            f"[{tag}] pool de CPUs vacío (todo reservado por otra "
            f"instancia de pipeline.py) — abortando lanzamiento de CENSO "
            f"en vez de correrlo sin confinar"
        )
        return 0, omp_min
    if n_reserved < requested_cores:
        log.warning(
            f"[{tag}] pool reducido: pedidos {requested_cores}, "
            f"reservados {n_reserved} (otra instancia tiene el resto) — "
            f"ajustando --maxcores/--omp-min/OMP_NUM_THREADS al valor real "
            f"en vez de sobre-suscribir el taskset confinado"
        )
    return n_reserved, min(omp_min, n_reserved)


# DEPRECADO (2026-07): reparto manual de slots de cores para los SP/OPT
# paralelos internos de CENSO (ver MAAQC_CENSO_SLOTS junto a
# build_omp_env() y partition_cpu_pool_into_slots()). Ya no se usa en
# ningún llamado a stage_censo_prescreen() — seguía mostrando colisión
# parcial entre sub-jobs pese al parche local en orca_processor.py, así
# que se abandonó a favor de dejar que el kernel reparta libremente
# dentro del pool confinado por taskset. build_omp_env() y
# partition_cpu_pool_into_slots() se dejan intactas por si se retoma el
# enfoque más adelante, pero ningún caller les pasa censo_slot_groups hoy.
MAAQC_CENSO_N_SLOTS = 2

# Los defaults de recursos (NPROCS, MAXCORE, etc.) se inicializan más abajo
# por ResourceProfile.detect() tras leer la topología real del sistema.
# Estos valores son placeholders que main() sobreescribe siempre.
NPROCS_CREST    = 12   # CREST/CENSO: xTB escala bien, priorizar más cores
NPROCS          = 4
NPROCS_NUMFREQ  = 8
NPROCS_TDDFT    = 6
MAXCORE         = 3500
MAXCORE_NUMFREQ = 2800
MAXCORE_TDDFT   = 2000
NPROCS_ADC2     = 4
MAXCORE_ADC2    = 5500
# Batch sizes para modo normal (sobreescritos por ResourceProfile en main())
OPT_BATCH_SIZE  = 3
SP_BATCH_SIZE   = 3
TDDFT_BATCH_SIZE = 2

# ── Topología REAL de este ThinkPad T14 (Core Ultra 7), confirmada con
#    `lscpu -e`. Elevada a nivel de módulo porque tanto assign_cpu_lists()
#    (pinning por-job, usado por ORCA) como assign_cpu_pool() (pinning de
#    un único proceso con paralelismo interno propio, usado por CENSO)
#    necesitan partir de la misma tabla CPU lógica → core físico. Es
#    también la tabla que build_omp_env() usa para traducir PU→Core al
#    armar MAAQC_PE_LIST (ver [MPI-2] en el CHANGELOG).
#
#    CPU lógica  →  CORE físico   →  tipo
#    0, 1        →  CORE 0        →  P-core (2 hilos SMT, mismo L2)
#    2, 3        →  CORE 1        →  P-core (2 hilos SMT, mismo L2)
#    4           →  CORE 2        →  E-core
#    5           →  CORE 3        →  E-core
#    6           →  CORE 4        →  E-core
#    7           →  CORE 5        →  E-core
#    8           →  CORE 6        →  E-core
#    9           →  CORE 7        →  E-core
#    10          →  CORE 8        →  E-core
#    11          →  CORE 9        →  E-core
#    12          →  CORE 10       →  LP E-core (max 400 MHz — evitar)
#    13          →  CORE 11       →  LP E-core (max ~1.3 GHz — evitar)
P_CORE_BLOCKS   = [[0, 1], [2, 3]]                 # CORE 0, CORE 1 (SMT)
E_CORE_SINGLES  = [4, 5, 6, 7, 8, 9, 10, 11]        # CORE 2-9
LP_CORE_SINGLES = [12, 13]                          # CORE 10-11 (lentos)

# ── Traducción PU (CPU lógica OS) → índice de Core, para PRTE ─────────────
# Ver [MPI-2] en el CHANGELOG: PRTE interpreta los números de
# `pe-list=`/`rankfile` como ÍNDICE DE CORE, nunca como PU/hilo de
# hardware — taskset -c, en cambio, sí usa PU directa, por eso solo lo
# que se manda a MAAQC_PE_LIST necesita esta traducción.
#
# Como los 2 P-cores tienen 2 PUs c/u (SMT) y el resto de los cores tiene
# 1 PU c/u, la tabla no es un offset constante. IMPORTANTE: esta
# traducción NO deduplica — si cpu_list trae los 2 hilos SMT de un mismo
# P-core (ej. [0, 1]), ambos deben traducir al MISMO índice de Core, para
# que `--bind-to hwthread` + `:ORDERED` reparta un hilo SMT distinto a
# cada rank dentro del mismo core. Deduplicar acá rompería esa asignación.
CORE_TO_PUS: dict[int, list[int]] = {
    0: [0, 1], 1: [2, 3],                        # P-cores, SMT
    2: [4], 3: [5], 4: [6], 5: [7],
    6: [8], 7: [9], 8: [10], 9: [11],            # E-cores
    10: [12], 11: [13],                          # LP E-cores
}
PU_TO_CORE: dict[int, int] = {
    pu: core for core, pus in CORE_TO_PUS.items() for pu in pus
}

# OMP_ENV se construye bajo demanda en build_omp_env() para que siempre
# refleje el valor actualizado de NPROCS (nunca queda obsoleto si el
# usuario pasa --nprocs desde CLI). Ver [MPI-1] y [MPI-2] en el CHANGELOG
# del módulo para el porqué de MAAQC_PE_LIST, set_pe_list y
# MAAQC_CENSO_SLOTS — el resumen operativo es:
#   - taskset -c (aplicado por el caller) siempre recibe PUs crudas.
#   - MAAQC_PE_LIST (leído por el wrapper de mpirun) necesita índices de
#     CORE, no de PU, así que se traduce vía PU_TO_CORE.
#   - set_pe_list=False (usado por el camino de CENSO) omite
#     MAAQC_PE_LIST para no forzar la misma lista fija sobre N sub-jobs
#     paralelos que CENSO lanza con el mismo ENVIRON heredado.
#   - MAAQC_CENSO_SLOTS es leído por un parche instalado aparte en CENSO
#     (fuera de este repo) para repartir esos sub-jobs por hash de
#     confórmero.
def build_omp_env(nprocs: int, cpu_list: Optional[list[int]] = None,
                   set_pe_list: bool = True,
                   censo_slot_groups: Optional[list[list[int]]] = None) -> dict:
    env = {
        **os.environ,
        "OMP_NUM_THREADS": str(nprocs),
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
    }
    # Limpiar las variables del namespace viejo si venían heredadas del
    # entorno del proceso padre — PRTE las ignora silenciosamente, pero
    # dejarlas sueltas es ruido y puede confundir a quien lea el entorno
    # de un proceso ORCA con `/proc/<pid>/environ` en el futuro.
    for stale_var in (
        "OMPI_MCA_hwloc_base_binding_policy",
        "OMPI_MCA_rmaps_base_mapping_policy",
        "OMPI_MCA_hwloc_base_use_hwthreads_as_cpus",
    ):
        env.pop(stale_var, None)

    if cpu_list and set_pe_list:
        # `cpu_list` son PUs/CPUs lógicas OS — la MISMA lista que se usa
        # para `taskset -c` en run_orca()/_run_censo_locked(). Pero PRTE
        # espera índices de CORE en pe-list, no PUs (ver historial arriba
        # y PU_TO_CORE junto a P_CORE_BLOCKS), así que hay que traducir
        # antes de setear la env var. Sin deduplicar: si dos PUs de
        # cpu_list caen en el mismo core físico (hilos SMT de un P-core),
        # el índice de ese core debe repetirse una vez por PU/rank para
        # que PRTE reparta un hilo SMT distinto a cada rank dentro del
        # mismo core (--bind-to hwthread + :ORDERED en el wrapper).
        pe_list = [PU_TO_CORE[pu] for pu in cpu_list]
        env["MAAQC_PE_LIST"] = ",".join(str(c) for c in pe_list)
    elif cpu_list and not set_pe_list:
        # cpu_list se usa solo para taskset (aplicado por el caller, ver
        # _run_censo_locked), NO para MAAQC_PE_LIST — ver Bug #1 arriba.
        # Si el proceso padre heredó un MAAQC_PE_LIST de una invocación
        # anterior (no debería, pero por si acaso), lo limpiamos para que
        # el wrapper de mpirun no lo use por accidente.
        env.pop("MAAQC_PE_LIST", None)

    if censo_slot_groups:
        # Cada grupo es una lista de PUs crudas (sin traducir a Core —
        # esa traducción la hace el parche de CENSO en el momento de
        # construir el MAAQC_PE_LIST real de cada llamada, reutilizando
        # la misma tabla PU_TO_CORE duplicada localmente allá).
        env["MAAQC_CENSO_SLOTS"] = "|".join(
            ",".join(str(pu) for pu in group) for group in censo_slot_groups
        )
    else:
        env.pop("MAAQC_CENSO_SLOTS", None)

    return env

# Alias global mantenido para compatibilidad con run_orca(env=None)
OMP_ENV: dict = build_omp_env(NPROCS)

STEPS = ["crest", "opt", "numfreq", "sp", "tddft", "adc2"]

# ─────────────────────────────────────────────
#  REGISTRO DE λmax EXPERIMENTAL
# ─────────────────────────────────────────────
#
# Asocia el nombre de carpeta de corrida (= mol_xyz.stem, p.ej. "A", "SA")
# con el λmax experimental en nm.  Se consulta en stage_adc2 para calcular
# Δλ = λ_ADC2 − λ_exp y registrarlo en el log y en adc2_results.json.
#
# Fuente canónica: molecules.yml en el directorio de trabajo (junto al script).
# Si el archivo no existe o no contiene la clave, se usa el diccionario
# de fallback EXP_LAMBDA_FALLBACK definido aquí.
#
# Formato molecules.yml:
#   exp_lambda:
#     A:   335.0   # M-Gly:Val estado A (carboxilato libre)
#     B:   335.0   # M-Gly:Val estado B
#     C:   335.0   # M-Gly:Val estado C
#     D:   335.0   # M-Gly:Val estado D (zwitterion)
#     R:   335.0   # M-Gly:Val estado R (neutro tautómero)
#     SA:  334.0   # Shinorine estado A
#     SB:  334.0   # Shinorine estado B
#     SR:  334.0   # Shinorine estado R
#
# Para agregar una molécula nueva sin tocar el código:
#   echo "  MiMol: 320.0" >> molecules.yml

EXP_LAMBDA_FALLBACK: dict[str, float] = {
    # M-Gly:Val — λmax experimental ~ 335 nm (UV-Vis en agua, Antarctic MAAs)
    "A":  335.0,
    "B":  335.0,
    "C":  335.0,
    "D":  335.0,
    "R":  335.0,
    # Shinorine — λmax experimental ~ 334 nm
    "SA": 334.0,
    "SB": 334.0,
    "SC": 334.0,
    "SD": 334.0,
    "SR": 334.0,
}

_MOLECULES_YML = Path("molecules.yml")   # ruta relativa al CWD donde corre el script


def load_exp_lambda(mol_name: str, log: logging.Logger) -> Optional[float]:
    """
    Devuelve el λmax experimental (nm) para mol_name.

    Orden de búsqueda:
      1. molecules.yml  →  sección exp_lambda:<mol_name>
      2. EXP_LAMBDA_FALLBACK  (diccionario en este script)
      3. None  (no disponible; stage_adc2 omite Δλ pero no falla)

    El YAML se lee sin PyYAML para no añadir dependencia: solo se parsean
    líneas "  clave: valor" dentro del bloque "exp_lambda:", lo que es
    suficiente para el formato plano que usamos.
    """
    # ── Intentar leer molecules.yml ──────────────────────────────
    if _MOLECULES_YML.exists():
        try:
            text = _MOLECULES_YML.read_text(errors="replace")
            in_block = False
            for line in text.splitlines():
                stripped = line.strip()
                if stripped.startswith("exp_lambda:"):
                    in_block = True
                    continue
                if in_block:
                    # Salir del bloque si aparece una clave de nivel superior (sin indentación)
                    if line and not line[0].isspace():
                        in_block = False
                        continue
                    # Parsear "  clave: valor"
                    if ":" in stripped and not stripped.startswith("#"):
                        key, _, val = stripped.partition(":")
                        key = key.strip()
                        val = val.split("#")[0].strip()   # ignorar comentarios inline
                        if key == mol_name:
                            try:
                                return float(val)
                            except ValueError:
                                log.warning(
                                    f"[EXP_LAMBDA] molecules.yml: valor inválido "
                                    f"para '{mol_name}': '{val}'"
                                )
        except OSError as e:
            log.warning(f"[EXP_LAMBDA] No se pudo leer molecules.yml: {e}")

    # ── Fallback: diccionario en el script ───────────────────────
    if mol_name in EXP_LAMBDA_FALLBACK:
        log.debug(
            f"[EXP_LAMBDA] '{mol_name}' no encontrado en molecules.yml "
            f"— usando fallback ({EXP_LAMBDA_FALLBACK[mol_name]:.1f} nm)"
        )
        return EXP_LAMBDA_FALLBACK[mol_name]

    log.warning(
        f"[EXP_LAMBDA] λmax experimental no disponible para '{mol_name}'. "
        f"Añádelo en molecules.yml bajo 'exp_lambda:' para obtener Δλ."
    )
    return None


# ─────────────────────────────────────────────
#  NOHUP — lanzamiento inmune a HUP
# ─────────────────────────────────────────────

# Controla si todos los subprocesos de cómputo (CREST, CENSO, ORCA)
# se lanzan con `nohup`. Cuando está activo, el proceso hijo ignora
# SIGHUP y sobrevive al cierre del terminal o a una desconexión SSH.
#
# Valor por defecto: True (comportamiento recomendado para cálculos largos).
#
# Cómo cambiarlo sin tocar el código:
#   — En el futuro se leerá desde environment.yml o pipeline_config.local.yml.
#     Por ahora, la única forma de desactivarlo es poner USE_NOHUP = False
#     aquí, o sobreescribir la variable antes de importar el módulo.
#
# Implementación:
#   _apply_nohup(cmd) antepone ["nohup"] al comando si USE_NOHUP es True.
#   Para subprocess.Popen (ORCA) se combina con start_new_session=True para
#   asegurar que el proceso hijo quede en su propio grupo de procesos y no
#   reciba SIGHUP del terminal padre.
USE_NOHUP: bool = True


def _apply_nohup(cmd: list[str]) -> list[str]:
    """
    Antepone 'nohup' al comando si USE_NOHUP está activo.

    Ejemplo:
        ["orca", "calc.inp"]  →  ["nohup", "orca", "calc.inp"]

    No modifica el comando si USE_NOHUP es False, lo que permite
    recuperar el comportamiento anterior sin tocar las llamadas a
    subprocess en el resto del script.
    """
    if USE_NOHUP:
        return ["nohup"] + cmd
    return cmd


# ─────────────────────────────────────────────
#  GESTIÓN DE SEÑALES / CTRL-C
# ─────────────────────────────────────────────

_active_procs: list[subprocess.Popen] = []   # procesos ORCA activos
_procs_lock   = threading.Lock()

def _signal_process_group(p: subprocess.Popen, sig: int) -> None:
    """Envía `sig` a TODO el grupo de procesos de p, no solo a su PID líder.

    Los procesos ORCA se lanzan con start_new_session=True, lo que los
    convierte en líderes de su propio grupo (pgid == pid). Cuando el
    comando real es mpirun (ORCA paralelo), mpirun lanza sus workers como
    hijos suyos dentro de ese mismo grupo. p.terminate()/p.kill() solo
    señalizan el PID de mpirun: si mpirun no alcanza a reenviar la señal
    a sus workers antes de morir (p.ej. porque le llega un SIGKILL a los
    3s), esos workers quedan huérfanos y siguen vivos consumiendo CPU —
    exactamente el síntoma de "Ctrl-C corta el pipeline pero los procesos
    siguen vivos". os.killpg() evita esto señalizando el grupo completo.
    """
    try:
        pgid = os.getpgid(p.pid)
        os.killpg(pgid, sig)
    except (ProcessLookupError, PermissionError, OSError):
        # El proceso ya no existe, o no se pudo resolver el pgid (p.ej. si
        # algún día USE_NOHUP/start_new_session pasa a False) — fallback al
        # comportamiento anterior de señalizar solo el PID directo.
        try:
            if sig == signal.SIGKILL:
                p.kill()
            else:
                p.terminate()
        except OSError:
            pass


def _kill_active_procs() -> None:
    """Mata todos los subprocesos ORCA/taskset registrados en _active_procs,
    incluyendo a los workers hijos que hayan lanzado (p.ej. mpirun), vía
    señales a todo el process group."""
    with _procs_lock:
        procs = list(_active_procs)
    for p in procs:
        _signal_process_group(p, signal.SIGTERM)
    # Espera breve y SIGKILL si siguen vivos
    for p in procs:
        try:
            p.wait(timeout=3)
        except subprocess.TimeoutExpired:
            _signal_process_group(p, signal.SIGKILL)

def _sigint_handler(signum, frame) -> None:  # noqa: ANN001
    print("\n[pipeline] Ctrl-C recibido — terminando procesos ORCA activos…",
          file=sys.stderr)
    _kill_active_procs()
    sys.exit(130)   # código estándar para Ctrl-C (128 + SIGINT)

signal.signal(signal.SIGINT,  _sigint_handler)
signal.signal(signal.SIGTERM, _sigint_handler)

# ─────────────────────────────────────────────
#  CONSTANTES FÍSICAS
# ─────────────────────────────────────────────

HARTREE_TO_KCAL     = 627.5094740631   # factor de conversión


# ─────────────────────────────────────────────
#  AUTODETECCIÓN DE RECURSOS
# ─────────────────────────────────────────────

class ResourceProfile:
    """
    Detecta la topología real del sistema y calcula los parámetros de
    recursos óptimos para cada etapa del pipeline ORCA.

    Reglas de distribución por tipo de cálculo
    ──────────────────────────────────────────
    CREST → xTB no tiene cuello de botella MPI como ORCA.
            El flag --T N reparte exploraciones iMTD-GC independientes entre
            núcleos de forma casi perfectamente lineal hasta saturar el BW de
            memoria (~12 cores en el T14 Gen5).
            Se usa MAX_CREST_NPROCS (default 12) para priorizar throughput.

    OPT  → 3 confórmeros en paralelo × 4 cores/job = 12 cores totales.
            El paralelismo entre confórmeros es perfecto (no hay comunicación
            entre jobs); el MPI interno de cada job sigue con su cuello de
            botella pero a 4 cores la eficiencia se mantiene alta.
            OPT_BATCH_SIZE = 3, MAX_OPT_NPROCS_PER_JOB = 4.

    NUMFREQ → paralelismo ideal: cada desplazamiento es un SCF independiente.
              Usa hasta MAX_NUMFREQ_NPROCS en un único job secuencial.
              No hay comunicación MPI entre desplazamientos → escala lineal.

    SP   → mismo esquema que OPT: 3 × 4 cores = 12 cores totales.
           SP_BATCH_SIZE = 3, MAX_SP_NPROCS = 4.

    TDDFT → 2 confórmeros en paralelo × 6 cores/job = 12 cores totales.
            Con NRoots 15 el consumo de RAM por job es mayor que SP,
            por eso se usan 2 jobs en vez de 3.
            Se aplica RAM_PENALTY_TDDFT sobre el maxcore calculado.
            TDDFT_BATCH_SIZE = 2, TDDFT_NPROCS_PER_JOB = 6.

    ADC2  → muy costoso en memoria (transformación integral 4-índice).
            Se limita explícitamente a NPROCS_ADC2_MAX cores y se asigna
            RAM generosa por core.

    maxcore (MB/core) se calcula como:
        floor(RAM_usable_MB / total_cores_activos / SAFETY)
    donde RAM_usable_MB = RAM_total × 0.85 (margen para SO y buffers).

    Parámetros ajustables al inicio de la clase:
        SAFETY_FACTOR          — fracción de RAM reservada para el SO
        MAX_CREST_NPROCS       — límite de cores para CREST/CENSO (xTB)
        MAX_NUMFREQ_NPROCS     — límite de cores para NUMFREQ
        MAX_OPT_NPROCS_PER_JOB — cores/job en OPT
        OPT_BATCH_SIZE         — jobs paralelos en OPT
        MAX_SP_NPROCS          — cores/job en SP
        SP_BATCH_SIZE          — jobs paralelos en SP modo normal
        TDDFT_NPROCS_PER_JOB   — cores/job en TDDFT
        TDDFT_BATCH_SIZE       — jobs paralelos en TDDFT
        RAM_PENALTY_TDDFT      — divisor extra de maxcore para TDDFT
        NPROCS_ADC2_MAX        — máximo cores para ADC2
        MAXCORE_ADC2_MIN       — mínimo MB/core para ADC2
    """

    # ── Parámetros de política (ajustar aquí si cambia el hardware) ──────
    SAFETY_FACTOR          = 0.85   # fracción de RAM que puede usar ORCA
    MAX_CREST_NPROCS       = 12    # CREST/xTB paraleliza bien (--T N): priorizar hasta 12 cores
    MAX_NUMFREQ_NPROCS     = 12    # NUMFREQ: escala lineal hasta saturar BW de memoria (~12 cores)
    MAX_OPT_NPROCS_PER_JOB = 4     # OPT: 3 jobs × 4 cores = 12 cores en paralelo
    OPT_BATCH_SIZE         = 3     # OPT modo normal: siempre 3 confórmeros en paralelo
    MAX_SP_NPROCS          = 4     # SP: 3 jobs × 4 cores = 12 cores en paralelo
    SP_BATCH_SIZE          = 3     # SP modo normal: siempre 3 confórmeros en paralelo
    TDDFT_NPROCS_PER_JOB   = 6     # TDDFT: 2 jobs × 6 cores = 12 cores en paralelo
    TDDFT_BATCH_SIZE       = 2     # TDDFT: 2 confórmeros en paralelo
    RAM_PENALTY_TDDFT      = 1.4   # TDDFT consume ~40% más RAM que SP por NRoots
    NPROCS_ADC2_MAX        = 4     # ADC2: más cores no ayudan, solo gastan RAM
    MAXCORE_ADC2_MIN       = 5000  # MB/core mínimo para ADC2 (transformación 4-índice)

    def __init__(
        self,
        total_cores : int,
        ram_mb      : int,
        batch_size  : int = 1,
    ) -> None:
        """
        Args:
            total_cores  : CPUs lógicas disponibles (de os.cpu_count() o --cores).
            ram_mb       : RAM total en MB (de psutil o estimada).
            batch_size   : jobs simultáneos por lote (default 1).
        """
        self.total_cores = total_cores
        self.ram_mb      = ram_mb
        self.batch_size  = batch_size
        self._compute()

    # ── Interfaz pública ────────────────────────────────────────────────

    @classmethod
    def detect(cls, batch_size: int = 1) -> "ResourceProfile":
        """
        Construye un ResourceProfile leyendo la topología real del sistema.

        Usa psutil para la RAM (más fiable que /proc/meminfo en algunos
        kernels). Si psutil no está disponible, estima 28 GB (valor
        conservador para el ThinkPad T14 Gen5 con 32 GB).
        """
        total_cores = os.cpu_count() or 4

        try:
            import psutil
            ram_mb = psutil.virtual_memory().total // (1024 * 1024)
        except ImportError:
            # Fallback: leer /proc/meminfo directamente
            try:
                with open("/proc/meminfo") as f:
                    for line in f:
                        if line.startswith("MemTotal"):
                            ram_mb = int(line.split()[1]) // 1024
                            break
                    else:
                        ram_mb = 28_000
            except OSError:
                ram_mb = 28_000

        return cls(total_cores=total_cores, ram_mb=ram_mb, batch_size=batch_size)

    @classmethod
    def detect_with_batch(cls, total_cores: int, batch_size: int, ram_mb: int = 0) -> "ResourceProfile":
        """
        Variante que recibe total_cores explícito (de --cores)
        y batch_size (de --parallel); recalcula el perfil de recursos.
        """
        if ram_mb <= 0:
            try:
                import psutil
                ram_mb = psutil.virtual_memory().total // (1024 * 1024)
            except (ImportError, Exception):
                try:
                    with open("/proc/meminfo") as f:
                        for line in f:
                            if line.startswith("MemTotal"):
                                ram_mb = int(line.split()[1]) // 1024
                                break
                        else:
                            ram_mb = 28_000
                except OSError:
                    ram_mb = 28_000

        return cls(total_cores=total_cores, ram_mb=ram_mb, batch_size=batch_size)

    def summary(self) -> str:
        """Devuelve un bloque de texto legible para el log de arranque."""
        ram_opt    = self.nprocs_opt    * self.maxcore_opt    / 1024
        ram_nf     = self.nprocs_numfreq * self.maxcore_numfreq / 1024
        ram_sp     = self.nprocs_sp     * self.maxcore_sp     / 1024
        ram_tddft  = self.nprocs_tddft  * self.maxcore_tddft  / 1024
        ram_adc2   = self.nprocs_adc2   * self.maxcore_adc2   / 1024
        ram_total  = self.ram_mb / 1024
        ram_usable = ram_total * self.SAFETY_FACTOR

        opt_batch  = self.batch_size if self.batch_size > 1 else self.OPT_BATCH_SIZE
        sp_batch   = self.batch_size if self.batch_size > 1 else self.SP_BATCH_SIZE

        return (
            f"  Sistema         : {self.total_cores} CPUs lógicas  "
            f"RAM {ram_total:.0f} GB total / {ram_usable:.0f} GB usable\n"
            f"  CREST/CENSO     : {self.nprocs_crest} cores  (xTB --T, escala lineal)\n"
            f"  OPT             : {opt_batch} jobs × {self.nprocs_opt} cores / {self.maxcore_opt} MB  "
            f"→ {ram_opt:.1f} GB/job\n"
            f"  SP              : {sp_batch} jobs × {self.nprocs_sp} cores / {self.maxcore_sp} MB  "
            f"→ {ram_sp:.1f} GB/job\n"
            f"  TDDFT           : {self.TDDFT_BATCH_SIZE} jobs × {self.nprocs_tddft} cores / {self.maxcore_tddft} MB  "
            f"→ {ram_tddft:.1f} GB/job  (×{self.RAM_PENALTY_TDDFT} RAM penalty, NRoots 15)\n"
            f"  NumFreq         : 1 job × {self.nprocs_numfreq} cores / {self.maxcore_numfreq} MB  "
            f"→ {ram_nf:.1f} GB  (escala lineal)\n"
            f"  ADC2            : 1 job × {self.nprocs_adc2} cores / {self.maxcore_adc2} MB  "
            f"→ {ram_adc2:.1f} GB"
        )

    # ── Cálculo interno ──────────────────────────────────────────────────

    def _compute(self) -> None:
        """Calcula todos los parámetros de recursos según las reglas de política."""
        c   = self.total_cores
        p   = self.batch_size
        ram = int(self.ram_mb * self.SAFETY_FACTOR)

        # ── CREST / CENSO ─────────────────────────────────────────────────
        # xTB paralleliza bien (--T N): cada núcleo corre una exploración
        # iMTD-GC independiente → escala casi lineal hasta saturar BW de memoria.
        # Se usa MAX_CREST_NPROCS para priorizar throughput sin desperdiciar cores
        # que beneficiarían a ORCA (que escala peor).
        self.nprocs_crest = min(c, self.MAX_CREST_NPROCS)

        # ── OPT ──────────────────────────────────────────────────────────
        # Siempre OPT_BATCH_SIZE jobs × MAX_OPT_NPROCS_PER_JOB cores.
        effective_opt_batch = p if p > 1 else self.OPT_BATCH_SIZE
        if p <= 1:
            self.nprocs_opt = min(c // max(self.OPT_BATCH_SIZE, 1),
                                  self.MAX_OPT_NPROCS_PER_JOB)
        else:
            self.nprocs_opt = min(max(1, c // p), self.MAX_OPT_NPROCS_PER_JOB)

        active_opt_cores = self.nprocs_opt * effective_opt_batch
        self.maxcore_opt = max(500, ram // max(active_opt_cores, 1))

        # ── NUMFREQ ──────────────────────────────────────────────────────
        # Siempre 1 job; escala lineal hasta MAX_NUMFREQ_NPROCS (saturaación BW mem).
        self.nprocs_numfreq  = min(c, self.MAX_NUMFREQ_NPROCS)
        self.maxcore_numfreq = max(500, ram // max(self.nprocs_numfreq, 1))

        # ── SP ───────────────────────────────────────────────────────────
        # Siempre SP_BATCH_SIZE jobs × MAX_SP_NPROCS cores.
        if p <= 1:
            self.nprocs_sp = min(c // max(self.SP_BATCH_SIZE, 1),
                                 self.MAX_SP_NPROCS)
        else:
            self.nprocs_sp = min(max(1, c // p), self.MAX_SP_NPROCS)

        effective_sp_batch = p if p > 1 else self.SP_BATCH_SIZE
        active_sp_cores = self.nprocs_sp * effective_sp_batch
        self.maxcore_sp = max(500, ram // max(active_sp_cores, 1))

        # ── TDDFT ────────────────────────────────────────────────────────
        # Siempre TDDFT_BATCH_SIZE jobs × TDDFT_NPROCS_PER_JOB cores.
        # maxcore reducido por RAM_PENALTY_TDDFT para acomodar vectores de
        # excitación (NRoots 15).
        self.nprocs_tddft  = min(c // max(self.TDDFT_BATCH_SIZE, 1),
                                 self.TDDFT_NPROCS_PER_JOB)
        active_tddft_cores = self.nprocs_tddft * self.TDDFT_BATCH_SIZE
        raw_maxcore_tddft  = max(500, ram // max(active_tddft_cores, 1))
        self.maxcore_tddft = max(500, int(raw_maxcore_tddft / self.RAM_PENALTY_TDDFT))

        # ── ADC2 ─────────────────────────────────────────────────────────
        # Costoso en memoria; se limita a NPROCS_ADC2_MAX y se prioriza RAM/core.
        self.nprocs_adc2  = min(c, self.NPROCS_ADC2_MAX)
        self.maxcore_adc2 = max(
            self.MAXCORE_ADC2_MIN,
            ram // max(self.nprocs_adc2, 1),
        )


# ─────────────────────────────────────────────
#  LOGGING
# ─────────────────────────────────────────────

def setup_logging(mol_dir: Path) -> logging.Logger:
    """
    Crea un logger con nombre único por molécula para evitar que dos
    instancias paralelas del script compartan handlers y dupliquen líneas.
    """
    # mol_dir es <grupo>/runs/<mol_name>/ — el nombre de la molécula es
    # directamente mol_dir.name (no mol_dir.parent.name, que sería "runs").
    logger_name = f"pipeline.{mol_dir.name}"
    log = logging.getLogger(logger_name)

    # Evitar añadir handlers duplicados si se llama más de una vez
    if log.handlers:
        return log

    log.setLevel(logging.DEBUG)
    fmt = logging.Formatter("%(asctime)s  %(levelname)-8s  %(message)s",
                            datefmt="%H:%M:%S")
    # consola
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(fmt)
    # archivo
    mol_dir.mkdir(parents=True, exist_ok=True)
    fh = logging.FileHandler(mol_dir / "pipeline.log")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)

    log.addHandler(ch)
    log.addHandler(fh)
    return log


# ─────────────────────────────────────────────
#  HELPERS XYZ / ASE
# ─────────────────────────────────────────────

def atoms_to_xyz_block(atoms) -> str:
    """Devuelve bloque de coordenadas estilo ORCA (* xyz ... *)."""
    lines = []
    for symbol, pos in zip(atoms.get_chemical_symbols(), atoms.get_positions()):
        lines.append(f"  {symbol:<4s}  {pos[0]:14.8f}  {pos[1]:14.8f}  {pos[2]:14.8f}")
    return "\n".join(lines)


def read_conformers(xyz_file: Path, nconf: int) -> list:
    """
    Lee hasta `nconf` confórmeros desde un multi-XYZ de CREST.
    Si el archivo tiene menos de nconf estructuras, devuelve los que haya
    sin lanzar excepción.
    """
    try:
        all_conf = read(str(xyz_file), index=":")
        # ASE puede devolver un solo Atoms en vez de lista si solo hay uno
        if not isinstance(all_conf, list):
            all_conf = [all_conf]
    except Exception as e:
        raise RuntimeError(f"No se pudo leer {xyz_file}: {e}") from e

    if len(all_conf) < nconf:
        # Advertir pero continuar con los disponibles
        import warnings
        warnings.warn(
            f"Se pidieron {nconf} confórmeros pero CREST solo generó "
            f"{len(all_conf)}. Se usarán todos los disponibles.",
            stacklevel=2,
        )
    return all_conf[:nconf]


# ─────────────────────────────────────────────
#  TEMPLATES
# ─────────────────────────────────────────────

TEMPLATE_OPT = """\
! wB97X-D3 def2-SVP CPCM(water) RIJCOSX TightSCF TightOpt
%method
  Grid 4
end
%scf
  MaxIter 300
end
%maxcore {MAXCORE}
%pal
  nprocs {NPROCS}
end

* xyz {CHARGE} {MULT}
{COORDS}
*
"""

TEMPLATE_NUMFREQ = """\
! wB97X-D3 def2-SVP CPCM(water) RIJCOSX TightSCF NumFreq
%method
  Grid 4
end
%scf
  MaxIter 300
end
%maxcore {MAXCORE_NUMFREQ}
%pal
  nprocs {NPROCS_NUMFREQ}
end

* xyzfile {CHARGE} {MULT} opt_opt.xyz
"""

TEMPLATE_NUMFREQ_RESTART = """\
! wB97X-D3 def2-SVP CPCM(water) RIJCOSX TightSCF NumFreq
%method
  Grid 4
end
%scf
  MaxIter 300
end
%maxcore {MAXCORE_NUMFREQ}
%pal
  nprocs {NPROCS_NUMFREQ}
end
%freq
  Restart true
end

* xyzfile {CHARGE} {MULT} opt_opt.xyz
"""

TEMPLATE_SP = """\
! wB97X-D3 def2-TZVPD def2/J CPCM(water) RIJCOSX TightSCF defgrid3
%scf
  MaxIter 300
end
%maxcore {MAXCORE}
%pal
  nprocs {NPROCS}
end

* xyzfile {CHARGE} {MULT} opt_opt.xyz
"""

TEMPLATE_SP_GAS = """\
! wB97X-D3 def2-TZVPD def2/J RIJCOSX TightSCF defgrid3
%scf
  MaxIter 300
end
%maxcore {MAXCORE}
%pal
  nprocs {NPROCS}
end

* xyzfile {CHARGE} {MULT} opt_opt.xyz
"""

TEMPLATE_SP_SMD = """\
! wB97X-D3 def2-TZVPD def2/J SMD(water) RIJCOSX TightSCF defgrid3
%scf
  MaxIter 300
end
%maxcore {MAXCORE}
%pal
  nprocs {NPROCS}
end

* xyzfile {CHARGE} {MULT} opt_opt.xyz
"""

TEMPLATE_SP_NROOT = """\
! wB97X-D3 def2-TZVPD def2/J CPCM(water) RIJCOSX TightSCF defgrid3
%scf
  MaxIter 300
end
%tddft
  NRoots    1
  MaxDim    20
  ETol      1e-6
  RTol      1e-6
  TDA       false
  Triplets  false
  PrintLevel 3
end
%maxcore {MAXCORE}
%pal
  nprocs {NPROCS}
end

* xyzfile {CHARGE} {MULT} opt_opt.xyz
"""

TEMPLATE_TDDFT = """\
! wB97X-D3 def2-TZVPD def2-TZVPD/C CPCM(water) RIJCOSX TightSCF
%pal
  nprocs {NPROCS}
end
%maxcore {MAXCORE}
%scf
  MaxIter 300
  MOInp "{GBW_PATH}"
end
%tddft
  NRoots    15
  MaxDim    150
  ETol      1e-6
  RTol      1e-6
  TDA       false
  Triplets  false
  PrintLevel 3
  DoCD      true
  DoNTO     true
end

* xyzfile {CHARGE} {MULT} sp_opt.xyz
"""

TEMPLATE_ADC2 = """\
! ADC2 def2-SVP def2-SVP/C RIJCOSX TightSCF
%scf
  MaxIter 500
end
%maxcore {MAXCORE_ADC2}
%pal
  nprocs {NPROCS_ADC2}
end
%mdci
  nroots 1
  trafotype trafo_ri
  maxiter 50
  MaxCore {MAXCORE_ADC2}
end

* xyzfile {CHARGE} {MULT} adc2_geom.xyz
"""


def get_template(step: str) -> str:
    """Devuelve el template inline. Si existe un .tmpl en templates/, lo prioriza."""
    tmpl_file = TMPL_DIR / f"{step}.inp.tmpl"
    if tmpl_file.exists():
        return tmpl_file.read_text()
    mapping = {
        "opt":             TEMPLATE_OPT,
        "numfreq":         TEMPLATE_NUMFREQ,
        "numfreq_restart": TEMPLATE_NUMFREQ_RESTART,
        "sp":              TEMPLATE_SP,
        "sp_gas":          TEMPLATE_SP_GAS,
        "sp_smd":          TEMPLATE_SP_SMD,
        "sp_nroot":        TEMPLATE_SP_NROOT,
        "tddft":           TEMPLATE_TDDFT,
        "adc2":            TEMPLATE_ADC2,
    }
    return mapping[step]


# ─────────────────────────────────────────────
#  VERIFICACIÓN DE ÉXITO
# ─────────────────────────────────────────────

def check_success(out_file: Path) -> bool:
    """Devuelve True si ORCA terminó normalmente."""
    if not out_file.exists():
        return False
    try:
        text = out_file.read_text(errors="replace")
    except OSError:
        return False
    return "ORCA TERMINATED NORMALLY" in text


def check_no_imaginary(out_file: Path) -> bool:
    """Devuelve True si no hay frecuencias imaginarias (NUMFREQ)."""
    if not out_file.exists():
        return True
    try:
        text = out_file.read_text(errors="replace")
    except OSError:
        return True
    for line in text.splitlines():
        if re.match(r"\s+\d+:\s+-\d", line):
            return False
    return True


def has_numfreq_restart(freq_dir: Path) -> bool:
    """
    Devuelve True si hay archivos de gradientes parciales de una NumFreq
    interrumpida. ORCA 6 guarda cada displacement como freq_D#####.res.Gradients.
    Con eso el Restart true puede retomar sin recalcular nada.
    """
    gradients = list(freq_dir.glob("freq_D*.res.Gradients"))
    return len(gradients) > 0


def extract_sp_energy(out_file: Path) -> Optional[float]:
    """Extrae la energía del SP en Hartree."""
    if not out_file.exists():
        return None
    try:
        for line in out_file.read_text(errors="replace").splitlines():
            if "FINAL SINGLE POINT ENERGY" in line:
                try:
                    return float(line.split()[-1])
                except ValueError:
                    pass
    except OSError:
        pass
    return None


def extract_adc2_lambda(out_file: Path) -> Optional[tuple[float, float, float]]:
    """
    Extrae λmax (nm), E (eV) y fosc de un output ADC(2)/ORCA (nroots 1, S0→S1).

    Busca el bloque "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    dentro de la sección "SPECTRUM FOR RIGHT TRANSITION MOMENTS" — que es el
    referente estándar para ADC(2) en ORCA — y lee la primera fila de datos
    (única raíz calculada, transición 0-1A → 1-1A).

    Formato esperado de la línea de datos:
        0-1A  ->  1-1A    <eV>   <cm-1>   <nm>   <fosc>   ...

    Devuelve:
        (lambda_nm, energy_eV, fosc)  si se encontraron los datos
        None                          si no se pudo parsear

    Notas:
    - El output ADC(2) de ORCA incluye tres tablas (RIGHT, LEFT, LEFT-RIGHT).
      Solo se extrae la tabla RIGHT (primera que aparece), ignorando las demás.
    - Se toma la primera línea de datos encontrada: válido cuando nroots=1.
    - El patrón de línea de cabecera es robusto a cambios menores de espaciado.
    """
    if not out_file.exists():
        return None

    # Patrón de la línea de datos: al menos 7 columnas separadas por espacios.
    # Columnas: label_from -> label_to  eV  cm-1  nm  fosc  D2  DX  DY  DZ
    # Ejemplo:
    #   0-1A  ->  1-1A    3.464792   27945.4   357.8   0.394893694   4.65206  ...
    _data_re = re.compile(
        r"^\s*\S+\s+->\s+\S+"           # etiquetas (e.g. "0-1A  ->  1-1A")
        r"\s+([\d.]+)"                   # grupo 1: energía en eV
        r"\s+[\d.]+\s+([\d.]+)"         # cm-1 (ignorado), grupo 2: nm
        r"\s+([\d.]+)"                   # grupo 3: fosc
    )
    _right_header = "SPECTRUM FOR RIGHT TRANSITION MOMENTS"
    _abs_header   = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"

    try:
        text = out_file.read_text(errors="replace")
    except OSError:
        return None

    # Localizar el bloque RIGHT (primera aparición) y buscar la tabla de absorción
    right_pos = text.find(_right_header)
    if right_pos == -1:
        return None

    # Dentro del bloque RIGHT buscar la línea de cabecera de la tabla
    abs_pos = text.find(_abs_header, right_pos)
    if abs_pos == -1:
        return None

    # Avanzar línea a línea desde la cabecera hasta la primera fila de datos
    # (hay una línea de sub-cabecera con "(eV) (cm-1) ..." y la línea de "---")
    segment = text[abs_pos:]
    for line in segment.splitlines():
        m = _data_re.match(line)
        if m:
            try:
                energy_ev  = float(m.group(1))
                lambda_nm  = float(m.group(2))
                fosc       = float(m.group(3))
                return (lambda_nm, energy_ev, fosc)
            except ValueError:
                return None

    return None


def safe_copy(src: Path, dst: Path) -> bool:
    """
    Copia src → dst solo si son archivos distintos.
    Devuelve True si la copia se realizó o no era necesaria.
    Evita SameFileError cuando CREST ya escribe el nombre canónico.
    """
    try:
        if src.resolve() == dst.resolve():
            return True
        shutil.copy(src, dst)
        return True
    except shutil.SameFileError:
        return True   # mismo archivo, no hay nada que hacer
    except OSError as e:
        return False


def _parse_charge_mult_from_sp_inp(sp_inp: Path) -> Optional[tuple[int, int]]:
    """
    Extrae (charge, mult) directo de la línea '* xyzfile {CHARGE} {MULT} ...'
    ya escrita en un sp.inp existente (04_sp/<conf_id>/sp.inp).

    Usado por SP_GAS para garantizar que usa exactamente la misma carga/mult
    que el SP CPCM que le dio origen, sin importar el modo de invocación
    (evita desalineación con args.charge en modo batch, donde cada corrida
    — B, C, D, R — tiene su propia carga).
    """
    m = re.search(r'^\*\s*xyzfile\s+(-?\d+)\s+(\d+)', sp_inp.read_text(), re.MULTILINE)
    return (int(m.group(1)), int(m.group(2))) if m else None


def resolve_opt_xyz(opt_dir: Path) -> Optional[Path]:
    """
    ORCA 6 escribe la geometría optimizada como 'opt.xyz'.
    Versiones anteriores usaban 'opt_opt.xyz'.
    Esta función busca ambos nombres y devuelve el que exista,
    normalizando siempre a 'opt_opt.xyz' mediante una copia si hace falta.
    Devuelve la ruta a opt_opt.xyz, o None si no se encuentra ninguno.
    """
    canonical = opt_dir / "opt_opt.xyz"
    if canonical.exists() and canonical.stat().st_size > 0:
        return canonical
    # ORCA 6: buscar opt.xyz y normalizar
    orca6_name = opt_dir / "opt.xyz"
    if orca6_name.exists() and orca6_name.stat().st_size > 0:
        shutil.copy(orca6_name, canonical)
        return canonical
    return None


# ─────────────────────────────────────────────
#  RUNNER ORCA
# ─────────────────────────────────────────────

def run_orca(inp_file: Path, out_file: Path, log: logging.Logger,
             timeout_hours: float = 24.0,
             env: Optional[dict] = None,
             cpu_list: Optional[list[int]] = None) -> bool:
    """
    Lanza ORCA de forma síncrona. Devuelve True si terminó bien.
    Si el .out ya existe y tiene TERMINATED NORMALLY, omite la ejecución.
    `env` permite pasar un entorno de proceso personalizado; si es None
    se construye dinámicamente desde los globales NPROCS actuales (evita
    usar el OMP_ENV construido al inicio del módulo, que puede estar obsoleto).

    `cpu_list`: lista explícita de CPUs lógicas (no necesariamente
    contiguas) a las que confinar el proceso vía `taskset -c cpu1,cpu2,...`.
    Se usa una lista en vez de un rango lo-hi porque la topología real de
    esta máquina (CPU híbrida P-core/E-core, ver `assign_cpu_lists`) no
    siempre permite asignar rangos contiguos sin partir un core físico
    SMT entre dos jobs distintos — confinar 2 ranks MPI a los 2 hilos SMT
    del MISMO P-core (ej. CPUs 0,1) no es un problema en sí, pero mezclar
    eso con asunciones de "cores contiguos = cores físicos distintos"
    llevó a un crash de hwloc en LEANSCF (ver historial: "intersection
    without inclusion" + abort en la formación de Fock).
    """
    if check_success(out_file):
        log.info(f"  ↩  ya completado: {out_file.name} — saltando")
        return True

    if not inp_file.exists():
        log.error(f"  ✗  no existe el input: {inp_file}")
        return False

    if not ORCA_BIN.exists():
        log.error(f"  ✗  no se encontró el ejecutable de ORCA en: {ORCA_BIN}")
        return False

    cmd = [str(ORCA_BIN), str(inp_file.resolve())]
    if cpu_list:
        cpu_str = ",".join(str(c) for c in cpu_list)
        cmd = ["taskset", "-c", cpu_str] + cmd
        log.info(f"  ▶  lanzando ORCA: {inp_file.name}  (cores {cpu_str})"
                 + ("  [nohup]" if USE_NOHUP else ""))
    else:
        log.info(f"  ▶  lanzando ORCA: {inp_file.name}"
                 + ("  [nohup]" if USE_NOHUP else ""))

    t0       = time.time()
    timeout  = timeout_hours * 3600
    proc_env = env if env is not None else build_omp_env(NPROCS)

    # start_new_session=True coloca el proceso en un nuevo grupo de procesos,
    # evitando que reciba SIGHUP al cerrarse el terminal padre — complementa
    # el nohup del comando cuando USE_NOHUP es True.  Con USE_NOHUP=False se
    # mantiene el comportamiento original: Ctrl-C propaga al grupo y
    # _sigint_handler puede terminar los jobs sin wait() pendiente.
    try:
        with open(out_file, "w") as fout:
            proc = subprocess.Popen(
                _apply_nohup(cmd),
                stdout=fout,
                stderr=subprocess.STDOUT,
                cwd=inp_file.parent,
                env=proc_env,
                start_new_session=USE_NOHUP,
            )

        # Registrar para que _sigint_handler pueda matar el proceso
        with _procs_lock:
            _active_procs.append(proc)

        try:
            proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
            log.error(f"  ✗  {inp_file.name} TIMEOUT ({timeout_hours} h)")
            return False
        finally:
            with _procs_lock:
                try:
                    _active_procs.remove(proc)
                except ValueError:
                    pass

        elapsed = (time.time() - t0) / 60
        if check_success(out_file):
            log.info(f"  ✓  {inp_file.name} terminó en {elapsed:.1f} min")
            return True
        else:
            # Extraer las últimas líneas del .out para diagnóstico rápido
            try:
                tail = out_file.read_text(errors="replace").splitlines()[-10:]
                log.debug(f"  últimas líneas de {out_file.name}:\n" +
                          "\n".join(f"    {l}" for l in tail))
            except OSError:
                pass
            log.error(f"  ✗  {inp_file.name} FALLÓ (returncode={proc.returncode}, "
                      f"{elapsed:.1f} min) — revisa {out_file}")
            return False
    except FileNotFoundError as e:
        if cpu_list and "taskset" in str(e):
            log.error(f"  ✗  no se encontró 'taskset' en PATH — instálalo "
                      f"(paquete util-linux) o llama a run_orca sin cpu_list")
        else:
            log.error(f"  ✗  no se pudo ejecutar ORCA — verifica ORCA_BIN: {ORCA_BIN}")
        return False
    except Exception as e:
        log.error(f"  ✗  error inesperado al lanzar ORCA: {e}")
        return False


# ─────────────────────────────────────────────
#  ETAPAS DEL PIPELINE
# ─────────────────────────────────────────────

def stage_crest(mol_xyz: Path, mol_dir: Path, charge: int,
                log: logging.Logger, autoprotonate: bool = False,
                autodeprotonate: bool = False) -> Optional[Path]:
    """
    Etapa 01_crest: corre CREST sobre la molécula de entrada.

    Modo normal: muestreo conformacional estándar.
    Modo autoprotonate/autodeprotonate: usa el modo --protonate/--deprotonate
    de CREST para generar y rankear sitios candidatos de (des)protonación,
    en vez del muestreo conformacional normal. Se mantienen --T/--chrg/
    --uhf/--alpb water igual que en modo normal — solo cambia el flag de
    modo. El resultado se normaliza a 01_crest/crest_conformers.xyz para
    que el resto del pipeline lo use sin cambios.

    ATENCIÓN: las estructuras resultantes de estos dos modos tienen carga
    real = charge ± 1 (CREST lo maneja internamente sumando/restando el
    protón sobre el --chrg de partida), pero esta función NO ajusta
    `charge` para el resto del pipeline — eso queda a criterio del caller
    (ver --autoprotonate/--autodeprotonate en el CLI).

    Devuelve la ruta al crest_conformers.xyz, o None si falla.
    """
    crest_dir = mol_dir / "01_crest"
    crest_dir.mkdir(parents=True, exist_ok=True)

    input_xyz      = crest_dir / "crest_input.xyz"
    conformers_xyz = crest_dir / "crest_conformers.xyz"

    # Skip si ya existe y no está vacío
    if conformers_xyz.exists() and conformers_xyz.stat().st_size > 0:
        log.info(f"[CREST] ya existe {conformers_xyz} — saltando")
        return conformers_xyz

    safe_copy(mol_xyz, input_xyz)

    if autoprotonate:
        modo_str = "  [--autoprotonate → crest --protonate]"
    elif autodeprotonate:
        modo_str = "  [--autodeprotonate → crest --deprotonate]"
    else:
        modo_str = ""

    log.info(f"[CREST] lanzando sobre {mol_xyz.name}  (carga={charge})"
             + modo_str
             + ("  [nohup]" if USE_NOHUP else ""))
    cmd = [str(CREST_BIN), str(input_xyz.name)]
    if autoprotonate:
        cmd.append("--protonate")
    elif autodeprotonate:
        cmd.append("--deprotonate")
    cmd += [
        "--T", str(NPROCS_CREST),
        "--chrg", str(charge),
        "--uhf", "0",
        "--alpb", "water",
    ]

    returncode = None
    t0 = time.time()
    try:
        with open(crest_dir / "crest.log", "w") as flog:
            proc = subprocess.run(
                _apply_nohup(cmd),
                stdout=flog,
                stderr=subprocess.STDOUT,
                cwd=crest_dir,
                env=build_omp_env(NPROCS_CREST),
                timeout=CREST_TIMEOUT_HOURS * 3600,
                start_new_session=USE_NOHUP,
            )
        returncode = proc.returncode
    except subprocess.TimeoutExpired:
        log.error(f"[CREST] TIMEOUT ({CREST_TIMEOUT_HOURS} h) — proceso terminado")
        return None
    except FileNotFoundError:
        log.error(f"[CREST] no se encontró el ejecutable de CREST: {CREST_BIN}")
        return None
    except Exception as e:
        log.error(f"[CREST] error inesperado: {e}")
        return None

    elapsed = (time.time() - t0) / 60

    # CREST puede escribir el archivo con distintos nombres. En modo
    # --protonate/--deprotonate el ensemble final se llama protonated.xyz /
    # deprotonated.xyz en vez de crest_conformers.xyz — se prioriza ese
    # nombre pero se mantienen los candidatos normales como fallback.
    if autoprotonate:
        output_candidates = ["protonated.xyz", "crest_conformers.xyz", "crest_best.xyz"]
    elif autodeprotonate:
        output_candidates = ["deprotonated.xyz", "crest_conformers.xyz", "crest_best.xyz"]
    else:
        output_candidates = ["crest_conformers.xyz", "crest_best.xyz"]

    for candidate in output_candidates:
        src = crest_dir / candidate
        if src.exists() and src.stat().st_size > 0:
            if not safe_copy(src, conformers_xyz):
                log.error(f"[CREST] error al copiar {src} → {conformers_xyz}")
                return None
            break

    if conformers_xyz.exists() and conformers_xyz.stat().st_size > 0:
        try:
            n = len(read(str(conformers_xyz), index=":"))
        except Exception:
            n = "?"
        etiqueta = ("sitios de protonación candidatos" if autoprotonate else
                    "sitios de desprotonación candidatos" if autodeprotonate else
                    "confórmeros")
        log.info(f"[CREST] {n} {etiqueta} generados en {elapsed:.1f} min")
        if autoprotonate or autodeprotonate:
            carga_real = charge + 1 if autoprotonate else charge - 1
            signo = "+ 1" if autoprotonate else "- 1"
            log.warning(
                f"[CREST] las estructuras de {conformers_xyz.name} tienen "
                f"carga REAL = {carga_real}  (--charge={charge} {signo}). "
                f"Las etapas siguientes de ESTE run seguirán usando "
                f"--charge={charge} sin ajustar — si vas a continuar el "
                f"pipeline sobre la estructura elegida, arranca un run "
                f"nuevo con --charge {carga_real}."
            )
        return conformers_xyz
    else:
        log.error(f"[CREST] no se encontró conformers.xyz "
                  f"(returncode={returncode}) — revisa {crest_dir / 'crest.log'}")
        return None


# ─────────────────────────────────────────────
#  RECUPERACIÓN DE PRESCREENING PARCIAL
# ─────────────────────────────────────────────

def _parse_prescreening_orca_sp(outfile: Path) -> Optional[float]:
    """
    Replica OrcaProc._sp de CENSO: busca 'FINAL SINGLE POINT ENERGY'
    en sp/sp.out y exige que el archivo contenga 'ORCA TERMINATED NORMALLY'.
    """
    if not outfile.exists():
        return None
    try:
        text = outfile.read_text(errors="ignore")
    except OSError:
        return None
    if "ORCA TERMINATED NORMALLY" not in text:
        return None
    for line in text.splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                return float(line.split()[-1])
            except (ValueError, IndexError):
                return None
    return None


def _parse_prescreening_xtb_energy(outfile: Path) -> Optional[float]:
    """
    Replica XtbProc._sp de CENSO: busca '| TOTAL ENERGY' y devuelve
    el token[3] (energía xTB en Hartree).
    """
    if not outfile.exists():
        return None
    try:
        text = outfile.read_text(errors="ignore")
    except OSError:
        return None
    for line in text.splitlines():
        if "| TOTAL ENERGY" in line:
            try:
                return float(line.split()[3])
            except (ValueError, IndexError):
                return None
    return None


def _censo_prescreen_diagnose_dir(censo_dir: Path, log: logging.Logger) -> dict:
    """
    Inspecciona censo_dir (01b_censo_prescreen/) y devuelve un dict con
    el estado de los cálculos de prescreening:

        {
          "conf_dirs"  : [Path, ...],   # carpetas CONF* encontradas
          "n_total"    : int,
          "n_sp_ok"    : int,           # con sp/sp.out + TERMINATED NORMALLY
          "n_gsolv_ok" : int,           # con xtb_gsolv/ completos
          "n_full"     : int,           # con los tres archivos correctos
          "can_recover": bool,          # True si n_full == n_total > 0
        }

    No lanza excepciones; si no hay carpetas CONF* devuelve todo en cero.
    """
    # CENSO 3.x escribe las carpetas CONF* dentro de censo_dir/0_PRESCREENING/
    # (ver censo.log: "...using XtbProc in 0_PRESCREENING/CONF141/xtb_gsolv...").
    # Antes solo se miraba censo_dir directamente y nunca se encontraba nada
    # que recuperar. Buscamos en censo_dir y, si existe, en censo_dir/0_PRESCREENING/.
    _search_roots = [censo_dir]
    _prescreen_subdir = censo_dir / "0_PRESCREENING"
    if _prescreen_subdir.is_dir():
        _search_roots.append(_prescreen_subdir)

    conf_dirs = sorted(
        {
            d for root in _search_roots if root.is_dir()
            for d in root.iterdir() if d.is_dir() and d.name.startswith("CONF")
        },
        key=lambda d: int(re.search(r"\d+", d.name).group())
                      if re.search(r"\d+", d.name) else 0,
    )
    n_total    = len(conf_dirs)
    n_sp_ok    = 0
    n_gsolv_ok = 0
    n_full     = 0

    for cd in conf_dirs:
        sp_ok  = _parse_prescreening_orca_sp(cd / "sp" / "sp.out") is not None
        gas_ok = _parse_prescreening_xtb_energy(cd / "xtb_gsolv" / "gas.out") is not None
        sol_ok = _parse_prescreening_xtb_energy(cd / "xtb_gsolv" / "solv.out") is not None
        if sp_ok:
            n_sp_ok += 1
        if gas_ok and sol_ok:
            n_gsolv_ok += 1
        if sp_ok and gas_ok and sol_ok:
            n_full += 1

    can_recover = (n_full == n_total > 0)

    log.debug(
        f"[CENSO-RECOVER]  diagnóstico {censo_dir.name}: "
        f"{n_total} CONF dirs | sp_ok={n_sp_ok} | gsolv_ok={n_gsolv_ok} | "
        f"full={n_full} | can_recover={can_recover}"
    )
    return {
        "conf_dirs"  : conf_dirs,
        "n_total"    : n_total,
        "n_sp_ok"    : n_sp_ok,
        "n_gsolv_ok" : n_gsolv_ok,
        "n_full"     : n_full,
        "can_recover": can_recover,
    }


def _censo_prescreen_recover(
    censo_dir: Path,
    log: logging.Logger,
) -> Optional[Path]:
    """
    Reconstruye 0_PRESCREENING_recovered.json desde los .out ya calculados
    en censo_dir (01b_censo_prescreen/), sin relanzar nada.

    Reproduce exactamente la estructura que espera EnsembleData.read_output
    de CENSO 3.0.x:
        {
          "part_name": "prescreening",
          "data": {
            "CONF1": {"energy": ..., "gsolv": ..., "grrho": 0.0, "gtot": ...},
            ...
          }
        }

    Devuelve la ruta al JSON si tuvo éxito con ≥1 confórmero, o None si
    hay confórmeros sin datos suficientes (en ese caso el llamador debe
    eliminar la carpeta y relanzar CENSO desde cero).

    Diferencia frente al script standalone recover_prescreening.py:
    • No escribe a stdout ni llama sys.exit(); usa el logger del pipeline.
    • Devuelve None en vez de sys.exit() cuando faltan datos.
    • Escribe 0_PRESCREENING_recovered.json junto al censo.log existente.
    """
    diag = _censo_prescreen_diagnose_dir(censo_dir, log)

    if diag["n_total"] == 0:
        log.debug("[CENSO-RECOVER]  no hay carpetas CONF* — nada que recuperar")
        return None

    data: dict[str, dict] = {}
    missing: list[str]    = []

    for cd in diag["conf_dirs"]:
        name = cd.name

        energy = _parse_prescreening_orca_sp(cd / "sp" / "sp.out")
        if energy is None:
            missing.append(name)
            log.debug(f"[CENSO-RECOVER]  {name}: sp.out ausente o incompleto")
            continue

        energy_gas  = _parse_prescreening_xtb_energy(cd / "xtb_gsolv" / "gas.out")
        energy_solv = _parse_prescreening_xtb_energy(cd / "xtb_gsolv" / "solv.out")

        if energy_gas is None or energy_solv is None:
            faltante = "gas.out" if energy_gas is None else "solv.out"
            missing.append(name)
            log.debug(f"[CENSO-RECOVER]  {name}: {faltante} ausente o sin parsear")
            continue

        gsolv = energy_solv - energy_gas
        grrho = 0.0
        gtot  = energy + gsolv + grrho

        data[name] = {
            "energy": energy,
            "gsolv" : gsolv,
            "grrho" : grrho,
            "gtot"  : gtot,
        }

    if missing:
        log.warning(
            f"[CENSO-RECOVER]  {len(missing)}/{diag['n_total']} confórmeros "
            f"sin datos completos: {', '.join(missing[:10])}"
            f"{'…' if len(missing) > 10 else ''}"
        )
        # CENSO exige el set COMPLETO en --reload; si faltan datos no podemos usar recovery.
        return None

    if not data:
        log.warning("[CENSO-RECOVER]  ningún confórmero recuperable")
        return None

    recovered_json = censo_dir / "0_PRESCREENING_recovered.json"
    dump = {
        "part_name": "prescreening",
        "data"     : data,
    }
    try:
        recovered_json.write_text(json.dumps(dump, indent=4))
    except OSError as e:
        log.error(f"[CENSO-RECOVER]  no se pudo escribir {recovered_json}: {e}")
        return None

    log.info(
        f"[CENSO-RECOVER]  ✓  JSON reconstruido → {recovered_json.name} "
        f"({len(data)} confórmeros)"
    )
    return recovered_json


# ─────────────────────────────────────────────
#  LOCK de CENSO — evita relanzamientos duplicados sobre el mismo censo_dir
# ─────────────────────────────────────────────
#
# Ver [CENSO-LOCK] en el CHANGELOG: el contenido de censo_dir (carpetas
# CONF*) no indica si CENSO sigue vivo, así que "recuperar" un run
# mirando solo eso podía relanzar un CENSO duplicado sobre uno que ya
# estaba en curso. El lockfile de abajo (con el PID real del proceso) es
# la única fuente de verdad sobre si CENSO sigue corriendo.
def _censo_lock_file(censo_dir: Path) -> Path:
    return censo_dir / ".censo_pipeline.lock"


def _censo_lock_owner(censo_dir: Path, log: logging.Logger) -> Optional[int]:
    """
    Devuelve el PID de un proceso CENSO activo sobre censo_dir, o None si
    no hay ninguno corriendo (lock ausente, obsoleto, o de un PID reciclado
    por el SO para otro programa).

    IMPORTANTE: esta es la unica fuente de verdad sobre si CENSO sigue
    vivo. El estado de las carpetas CONF*/*.out NO sirve para eso (ver
    comentario arriba) -- por eso el lock se escribe con el PID real del
    Popen, no se infiere del contenido del directorio.
    """
    lock = _censo_lock_file(censo_dir)
    if not lock.exists():
        return None
    try:
        pid = int(lock.read_text().strip())
    except (OSError, ValueError):
        log.warning(f"[CENSO-LOCK]  {lock.name} ilegible -- se ignora y se elimina")
        try:
            lock.unlink()
        except OSError:
            pass
        return None

    try:
        os.kill(pid, 0)   # no mata nada, solo comprueba si el PID existe
    except ProcessLookupError:
        log.debug(f"[CENSO-LOCK]  lock obsoleto (PID {pid} ya no existe) -- limpiando")
        try:
            lock.unlink()
        except OSError:
            pass
        return None
    except PermissionError:
        # el PID existe pero pertenece a otro usuario -- asumir vivo por seguridad
        return pid

    # PID vivo: si psutil esta disponible, confirmar que en efecto es un
    # proceso de censo y no un PID reciclado por el SO para otro programa
    # (evita falsos positivos que bloquearian el pipeline indefinidamente).
    try:
        import psutil
        cmdline = " ".join(psutil.Process(pid).cmdline()).lower()
        if "censo" not in cmdline:
            log.debug(
                f"[CENSO-LOCK]  PID {pid} vivo pero no parece ser censo "
                f"(cmdline no contiene 'censo') -- lock obsoleto, limpiando"
            )
            try:
                lock.unlink()
            except OSError:
                pass
            return None
    except Exception:
        pass  # sin psutil o proceso inaccesible: confiar en la sola existencia del PID

    return pid


def _run_censo_locked(cmd: list[str], censo_dir: Path, log_file: Path,
                       env: dict, timeout_hours: float,
                       cpu_list: Optional[list[int]] = None) -> subprocess.CompletedProcess:
    """
    Reemplazo drop-in de subprocess.run() para lanzar CENSO: misma firma de
    excepciones (subprocess.TimeoutExpired, FileNotFoundError) pero ademas
    escribe .censo_pipeline.lock con el PID real apenas arranca el proceso
    y lo borra siempre al terminar (exito, timeout, o error), via Popen en
    vez de run() para poder conocer el PID antes de bloquear en wait().

    Tambien registra el proceso en _active_procs para que Ctrl-C
    (_sigint_handler) lo mate igual que a los procesos ORCA.

    `cpu_list`: confina el proceso CENSO (y por herencia TODOS los xtb
    que lanza internamente vía su propio scheduler Dask) a estas CPUs
    lógicas mediante `taskset -c`. Sin esto, CENSO calcula correctamente
    cuántos sub-jobs correr en paralelo y con cuántos hilos OpenMP cada
    uno (respetando --maxcores), pero el kernel queda libre de repartir
    esos hilos entre CUALQUIERA de las CPUs lógicas de la máquina —
    incluyendo aterrizar hilos de sub-jobs *distintos* sobre el mismo
    par SMT de un P-core, o sobre los LP E-cores lentos. Resultado
    medido: varios procesos xtb paralelos compitiendo por los mismos
    cores físicos, con solo ~6/12 cores realmente aprovechados aunque
    CENSO cree estar usando 12. Ver assign_cpu_pool() para el porqué
    de por qué este pinning es a nivel de proceso completo y no
    por-sub-job (a diferencia de run_orca()/assign_cpu_lists()).
    """
    if cpu_list:
        cpu_str = ",".join(str(c) for c in cpu_list)
        cmd = ["taskset", "-c", cpu_str] + cmd

    lock = _censo_lock_file(censo_dir)
    with open(log_file, "w") as flog:
        proc = subprocess.Popen(
            _apply_nohup(cmd), stdout=flog, stderr=subprocess.STDOUT,
            cwd=censo_dir, env=env, start_new_session=USE_NOHUP,
        )
        try:
            lock.write_text(str(proc.pid))
        except OSError:
            pass
        with _procs_lock:
            _active_procs.append(proc)
        try:
            returncode = proc.wait(timeout=timeout_hours * 3600)
        except subprocess.TimeoutExpired:
            proc.kill()
            try:
                proc.wait(timeout=30)
            except subprocess.TimeoutExpired:
                pass
            raise
        finally:
            with _procs_lock:
                if proc in _active_procs:
                    _active_procs.remove(proc)
            try:
                lock.unlink(missing_ok=True)
            except OSError:
                pass
    return subprocess.CompletedProcess(cmd, returncode)


def stage_censo_prescreen(conformers_xyz: Path, mol_dir: Path,
                           charge: int, log: logging.Logger,
                           temperature: float = 274.15,
                           subdir: Optional[str] = None) -> Path:
    """
    Etapa 01b_censo: pre-screening DFT liviano con CENSO antes del OPT
    costoso (ωB97X-D3/def2-SVP). Descarta confórmeros de alta energía con
    un nivel barato (funcional/base/threshold definidos en ~/.censo2rc,
    no por CLI en CENSO 3.x).

    Fail-open: si CENSO no está, falla, o no produce un ensemble
    reconocible, se devuelve conformers_xyz sin tocar y el pipeline
    sigue exactamente como antes (no bloquea nada).

    Args:
        conformers_xyz : ruta al ensemble CREST (multi-XYZ).
        mol_dir        : directorio raíz del run (<grupo>/runs/<mol_name>/).
        charge         : carga formal de la molécula.
        log            : logger del pipeline.
        temperature    : temperatura para CENSO en K (default 274.15 = 1 °C Antártico).
        subdir         : si se especifica, CENSO corre en 01b_censo_prescreen/<subdir>/
                         (usado en modo --multitemp para aislar cada temperatura).
                         Si es None, usa 01b_censo_prescreen/ directamente.

    Devuelve la ruta al ensemble filtrado (multi-XYZ), o conformers_xyz
    original si el pre-screening no se pudo aplicar.
    """
    base_censo_dir = mol_dir / "01b_censo_prescreen"
    censo_dir = (base_censo_dir / subdir) if subdir else base_censo_dir
    censo_dir.mkdir(parents=True, exist_ok=True)

    # ── Guardia de concurrencia: ¿ya hay un CENSO vivo sobre censo_dir? ──────
    # Esto va ANTES de cualquier diagnóstico o lectura de CONF*/*.out: el
    # bug real (ver comentario sobre _censo_lock_owner más arriba) era que
    # esta función juzgaba si CENSO estaba "atascado" mirando el contenido
    # del directorio, que sigue existiendo aunque CENSO esté vivo y ya haya
    # avanzado a screening. Si hay un PID vivo, no tocamos nada — ni
    # diagnóstico, ni --reload, ni limpieza de carpetas — y devolvemos el
    # ensemble sin filtrar tal cual llegó, para que el pipeline no bloquee
    # ni piense que falló (fail-open, igual que el resto de esta función).
    _owner_pid = _censo_lock_owner(censo_dir, log)
    if _owner_pid is not None:
        log.warning(
            f"[CENSO]  ⏳  ya hay un proceso CENSO activo sobre {censo_dir.name}/ "
            f"(PID {_owner_pid}) — no se toca nada para no pisarle el trabajo "
            f"(este era el bug: relanzar --reload sobre un run que seguía vivo "
            f"y ya había avanzado a screening). Vuelve a invocar el pipeline "
            f"más tarde, cuando termine."
        )
        return conformers_xyz

    output_xyz = censo_dir / "censo_ensemble_sorted.xyz"

    # ── Verificación de que el ensemble proviene de SCREENING, no solo de
    # prescreening. output_xyz es un nombre genérico que históricamente se
    # llenaba con cualquier .xyz encontrado (incluido el de prescreening);
    # ahora exigimos evidencia de que 1_SCREENING.out existe y es parseable
    # antes de considerar el run completo. Si solo hay prescreening (p.ej.
    # de un run anterior con --prescreening sin --screening), se trata como
    # incompleto y se relanza CENSO con --screening incluido.
    _screening_out_candidates = [
        censo_dir / "1_SCREENING.out",
        censo_dir / "CENSO" / "1_SCREENING.out",
        censo_dir / "prescreening" / "1_SCREENING.out",
    ]
    _screening_conf_ids: Optional[list[str]] = None
    for _p in _screening_out_candidates:
        _parsed = parse_censo_screening_out(_p)
        if _parsed:
            _screening_conf_ids = _parsed
            break
    _has_screening = _screening_conf_ids is not None

    # ── Reconstrucción DIRECTA desde 1_SCREENING.out + 1_SCREENING.xyz ───────
    # FIX (2026-07-14): si el SCREENING (no solo el prescreening) ya terminó
    # — es decir, _has_screening es True —, NO hace falta invocar CENSO de
    # nuevo bajo ninguna forma (ni limpio, ni --reload): las energías
    # Gtot/ΔGtot/Boltzmann y las geometrías de los confórmeros sobrevivientes
    # YA están en disco (1_SCREENING.out + 1_SCREENING.xyz). Este bloque se
    # ejecuta ANTES de la verificación de output_xyz y ANTES del bloque de
    # diagnóstico/--reload sobre 0_PRESCREENING/ más abajo.
    #
    # Bug real que esto corrige: con prescreening Y screening ya 100%
    # completos en disco, el flujo anterior igual llegaba al bloque de
    # "recuperación de carpeta parcial de CENSO" (_censo_prescreen_diagnose_dir
    # sobre 0_PRESCREENING/), que diagnostica sobre las carpetas CONF* del
    # PRESCREENING — no sabe nada de que el SCREENING ya terminó — y podía
    # terminar relanzando CENSO (limpio o con --reload) para un trabajo que
    # ya estaba hecho y disponible en 1_SCREENING.out. Ahora, apenas se
    # confirma que 1_SCREENING.out es parseable, se intenta reconstruir el
    # ensemble directamente desde los archivos existentes y se retorna sin
    # tocar CENSO en absoluto.
    if _has_screening:
        _screening_out_path = next(
            (p for p in _screening_out_candidates if parse_censo_screening_out(p)),
            None,
        )
        # Fuentes de geometría en orden de prioridad:
        #   1. conformers_xyz original de CREST (comentarios "CONF7  -E")
        #   2. 1_SCREENING.xyz generado por CENSO (comentarios = solo CONF#)
        _geom_sources = [conformers_xyz]
        if _screening_out_path is not None:
            for _cand in (
                _screening_out_path.parent / "1_SCREENING.xyz",
                _screening_out_path.parent.parent / "1_SCREENING.xyz",
            ):
                if _cand.exists() and _cand.stat().st_size > 0 and _cand not in _geom_sources:
                    _geom_sources.append(_cand)
        _rebuilt_direct = False
        for _geom_src in _geom_sources:
            if extract_conformers_by_id(_geom_src, _screening_conf_ids, output_xyz, log):
                log.info(
                    f"[CENSO]  ✓  screening ya completo en disco "
                    f"({(_screening_out_path.relative_to(censo_dir) if _screening_out_path else '1_SCREENING.out')}) "
                    f"— {len(_screening_conf_ids)} confórmero(s), ensemble reconstruido "
                    f"SIN invocar CENSO (geometrías: {_geom_src.name}) → {output_xyz.name}"
                )
                _rebuilt_direct = True
                break
            log.debug(
                f"[CENSO]  fuente {_geom_src.name} no contiene los CONF# de "
                f"1_SCREENING.out — probando siguiente fuente de geometrías"
            )
        if _rebuilt_direct:
            return output_xyz
        log.warning(
            f"[CENSO]  ⚠  1_SCREENING.out parseable ({len(_screening_conf_ids)} conf(s)) "
            f"pero ninguna fuente de geometría disponible ({', '.join(g.name for g in _geom_sources)}) "
            f"contiene esos CONF# — se continúa con la verificación de "
            f"{output_xyz.name} / recovery de prescreening más abajo"
        )

    # ── Verificación cruzada output_xyz ↔ 1_SCREENING.out / conformers_xyz ──
    # El bug real: "existe output_xyz" y "existe 1_SCREENING.out parseable"
    # se verificaban por separado, sin comprobar que ambos describen el
    # MISMO screening. Un censo_ensemble_sorted.xyz viejo (de un run
    # anterior, --reload fallido, --multitemp sin subdir, etc.) que quedó
    # en este mismo censo_dir podía convivir con un 1_SCREENING.out fresco
    # y válido, y esta función los trataba como pareja sin más evidencia.
    # Ahora exigimos que output_xyz tenga el mismo n° de átomos por frame
    # que el ensemble CREST de este run (conformers_xyz, ya garantizado
    # correcto porque viene del mismo mol_dir), y que al menos el CONF#
    # de mayor peso Boltzmann del .out exista como frame real en output_xyz.
    # Si algo no cuadra, no confiamos en el archivo viejo: se re-corre CENSO.
    if output_xyz.exists() and output_xyz.stat().st_size > 0 and _has_screening:
        _stale_reason = None
        try:
            _crest_atoms = read(str(conformers_xyz), index="0")
            _expected_n = len(_crest_atoms)
            _out_atoms = read(str(output_xyz), index=":")
            if not isinstance(_out_atoms, list):
                _out_atoms = [_out_atoms]
            if not _out_atoms:
                _stale_reason = f"{output_xyz.name} no tiene frames"
            elif len(_out_atoms[0]) != _expected_n:
                _stale_reason = (
                    f"{output_xyz.name} tiene {len(_out_atoms[0])} átomos/frame, "
                    f"pero el ensemble CREST de este run ({conformers_xyz.name}) "
                    f"tiene {_expected_n} — no pueden ser del mismo screening"
                )
            else:
                _out_text = output_xyz.read_text(errors="replace")
                _top_conf = _screening_conf_ids[0] if _screening_conf_ids else None
                if _top_conf and _top_conf not in _out_text:
                    _stale_reason = (
                        f"el confórmero de mayor peso según 1_SCREENING.out "
                        f"({_top_conf}) no aparece en {output_xyz.name}"
                    )
        except Exception as e:
            _stale_reason = f"no se pudo verificar consistencia ({e})"

        if _stale_reason is None:
            log.info(f"[CENSO]  ↩  screening ya completado — saltando ({output_xyz.name})")
            return output_xyz
        else:
            log.warning(
                f"[CENSO]  ⚠  {output_xyz.name} existe y hay un 1_SCREENING.out "
                f"parseable, pero NO son consistentes entre sí: {_stale_reason}. "
                f"Esto indica un {output_xyz.name} obsoleto de otra corrida "
                f"(posiblemente otro estado de protonación) que quedó en "
                f"{censo_dir} y no fue sobrescrito. Se lo descarta y se "
                f"relanza CENSO desde cero para regenerarlo correctamente."
            )
            try:
                output_xyz.unlink()
            except OSError:
                pass
            _has_screening = False

    if output_xyz.exists() and output_xyz.stat().st_size > 0 and not _has_screening:
        log.warning(
            f"[CENSO]  ⚠  {output_xyz.name} existe pero proviene solo de "
            f"prescreening (no se encontró 1_SCREENING.out parseable) — "
            f"se considera incompleto. Relanzando CENSO con --screening."
        )

    # ── Recuperación de carpeta parcial de CENSO (prescreening interrumpido) ──
    # Si existe la carpeta 01b_censo_prescreen/ con carpetas CONF* internas
    # (generadas por un run de CENSO interrumpido), intentamos reconstruir
    # 0_PRESCREENING_recovered.json y relanzar CENSO con --reload en lugar
    # de correr todo desde cero.
    #
    # Flujo:
    #   1. Detectar si hay carpetas CONF* en censo_dir.
    #   2. Si hay un 0_PRESCREENING_recovered.json previo, reutilizarlo
    #      directamente (permite reintentos tras un --reload que falló).
    #   3. Diagnosticar si los .out de todas las CONF* son parseables.
    #   4a. Si SÍ → reconstruir JSON y relanzar con --reload.
    #   4b. Si NO (datos incompletos) → borrar carpetas corruptas y dejar
    #       que CENSO corra desde cero más abajo.
    #
    # Nota: la carpeta censo_dir en sí se conserva (contiene el censo.log
    # original y el crest_conformers.xyz copiado); solo se eliminan las
    # subcarpetas CONF* que impiden a CENSO correr limpio.
    _diag = _censo_prescreen_diagnose_dir(censo_dir, log)
    if _diag["n_total"] > 0:
        log.info(
            f"[CENSO]  ⚙  carpeta previa detectada con {_diag['n_total']} "
            f"CONF dirs (sp_ok={_diag['n_sp_ok']}, gsolv_ok={_diag['n_gsolv_ok']}, "
            f"completos={_diag['n_full']})"
        )

        # ── ¿Ya existe un JSON de recovery de una ejecución anterior? ──────
        _prev_json = censo_dir / "0_PRESCREENING_recovered.json"
        if _prev_json.exists() and _prev_json.stat().st_size > 0:
            log.info(
                f"[CENSO]  ♻  reutilizando {_prev_json.name} existente "
                f"— relanzando CENSO con --reload"
            )
            _recovered_json: Optional[Path] = _prev_json
        elif _diag["can_recover"]:
            # Todos los confórmeros tienen datos completos → reconstruir
            log.info(
                f"[CENSO]  🔧  todos los .out son parseables "
                f"— reconstruyendo {censo_dir.name}/0_PRESCREENING_recovered.json"
            )
            _recovered_json = _censo_prescreen_recover(censo_dir, log)
        elif _diag["n_full"] > 0:
            # Recovery parcial: algunos CONF* completos, otros no.
            # CENSO --reload exige el set COMPLETO (no parcial), así que
            # no podemos mezclar incompletos en el JSON. En vez de descartar
            # los n_full ya calculados, borramos SOLO las carpetas incompletas
            # y reconstruimos el JSON sobre el subconjunto consistente que
            # queda (n_full == nuevo n_total). Los pocos incompletos se
            # recalculan al relanzar CENSO limpio sobre el ensemble reducido.
            _n_incomplete = _diag["n_total"] - _diag["n_full"]
            log.warning(
                f"[CENSO]  ⚙  {_diag['n_full']}/{_diag['n_total']} CONF dirs "
                f"completas — eliminando {_n_incomplete} incompleta(s) y "
                f"reconstruyendo recovery sobre el resto en vez de recalcular todo"
            )
            for _cd in _diag["conf_dirs"]:
                sp_ok  = _parse_prescreening_orca_sp(_cd / "sp" / "sp.out") is not None
                gas_ok = _parse_prescreening_xtb_energy(_cd / "xtb_gsolv" / "gas.out") is not None
                sol_ok = _parse_prescreening_xtb_energy(_cd / "xtb_gsolv" / "solv.out") is not None
                if not (sp_ok and gas_ok and sol_ok):
                    try:
                        shutil.rmtree(_cd)
                        log.debug(f"[CENSO]  rm -rf {_cd.name} (incompleta)")
                    except OSError as _e:
                        log.warning(f"[CENSO]  no se pudo eliminar {_cd.name}: {_e}")
            _diag = _censo_prescreen_diagnose_dir(censo_dir, log)
            if _diag["can_recover"]:
                _recovered_json = _censo_prescreen_recover(censo_dir, log)
            else:
                log.warning(
                    "[CENSO]  ⚠  diagnóstico tras limpieza parcial sigue "
                    "inconsistente — recovery no es posible"
                )
                _recovered_json = None
        else:
            # Datos incompletos → el recovery no es posible; limpiar y relanzar
            _recovered_json = None

        if _recovered_json is not None:
            # ── Relanzar CENSO con --reload sobre el JSON reconstruido ──────
            local_input = censo_dir / conformers_xyz.name
            _reload_attempted = False  # default: sin esto, referenciar
                                        # _reload_attempted más abajo lanzaría
                                        # NameError si safe_copy falla aquí
            if not safe_copy(conformers_xyz, local_input):
                log.warning(
                    f"[CENSO-RECOVER]  no se pudo copiar {conformers_xyz.name} "
                    f"— abortando recovery, relanzando limpio"
                )
            else:
                # ago 2026: reservar el pool y reconciliarlo ANTES de armar
                # cmd_reload — mismo fix que en el lanzamiento principal, ver
                # _reconcile_censo_budget(). Si el pool queda vacío, no se
                # intenta el --reload: se deja _reload_attempted=False para
                # que el flujo caiga directo a la limpieza + relanzamiento
                # normal más abajo (que reconcilia el pool de nuevo, fresco,
                # sin dejar _diag en un estado inconsistente).
                _recover_cpu_pool = _acquire_cpu_pool(NPROCS_CREST, log, "CENSO-RECOVER")
                _eff_maxcores, _eff_omp_min = _reconcile_censo_budget(
                    NPROCS_CREST, CENSO_OMP_MIN, _recover_cpu_pool, log, "CENSO-RECOVER"
                )
                _reload_attempted = _eff_maxcores > 0
                if not _reload_attempted:
                    log.error(f"[CENSO-RECOVER]  ✗  sin CPUs disponibles — "
                              f"se omite --reload, se limpiará la carpeta y "
                              f"CENSO correrá desde cero")
            if _reload_attempted:
                cmd_reload = [
                    str(CENSO_BIN),
                    # FIX (bug del 10 jul): NO pasar "--prescreening" aquí.
                    # _recovered_json ya contiene el part "prescreening"
                    # COMPLETO para los n_full == n_total confórmeros (así
                    # lo garantiza _censo_prescreen_recover, que devuelve
                    # None si falta algún dato). Pero si además le decimos
                    # a CENSO "--prescreening" en el CLI, CENSO ejecuta ese
                    # part igual — el --reload no lo hace tratarlo como ya
                    # resuelto, solo le da un punto de partida — y recalcula
                    # xtb_gsolv + sp para las 192 CONF* desde cero. Eso era
                    # exactamente lo que producía "Reloaded results..."
                    # seguido inmediatamente de "Running xtb_gsolv... CONF1"
                    # en censo.log: tres corridas completas de prescreening
                    # (~2.5 h cada una) sin avanzar nunca a screening.
                    # Con solo "--screening" + "--reload", CENSO toma el
                    # prescreening ya resuelto como dado y arranca directo
                    # en screening.
                    "--screening",
                    "-i", local_input.name,
                    "-c", str(charge),
                    "-u", "0",
                    "--solvent", "water",
                    "--maxcores", str(_eff_maxcores),   # ago 2026: reconciliado
                                                         # contra len(_recover_cpu_pool)
                                                         # real, no NPROCS_CREST nominal
                    "--omp-min", str(_eff_omp_min),     # ídem
                    "-T", f"{temperature:.2f}",
                    "--reload", str(_recovered_json),
                ]
                log.info(
                    f"[CENSO-RECOVER]  ▶  relanzando CENSO con --reload "
                    f"({_recovered_json.name}), solo --screening "
                    f"(prescreening se toma del reload, no se recalcula)"
                    + ("  [nohup]" if USE_NOHUP else "")
                )
                try:
                    # 2026-07: mismo cambio que en el lanzamiento principal
                    # (ver comentario junto a censo_cpu_pool más abajo en
                    # stage_censo_prescreen) — se abandonó el reparto manual
                    # vía MAAQC_CENSO_SLOTS; _recover_cpu_pool solo fija el
                    # límite exterior vía taskset (necesario para coordinar
                    # con otras instancias de pipeline.py corriendo en
                    # paralelo), y el kernel reparte libremente los
                    # sub-jobs internos de CENSO dentro de ese pool.
                    #
                    # ago 2026: OMP_NUM_THREADS del env ya NO es NPROCS_CREST
                    # (el pool completo) sino censo_threads_per_job(...) —
                    # ver CENSO_OMP_MIN arriba para el porqué (sin esto, cada
                    # sub-job concurrente de CENSO pedía el pool completo de
                    # hilos para sí mismo, sobre-suscribiendo 3x el pool
                    # confinado por taskset).
                    _recover_cpu_pool = _acquire_cpu_pool(NPROCS_CREST, log, "CENSO-RECOVER")
                    _proc = _run_censo_locked(
                        cmd_reload, censo_dir,
                        censo_dir / "censo_recover.log",
                        build_omp_env(censo_threads_per_job(NPROCS_CREST),
                                       cpu_list=_recover_cpu_pool,
                                       set_pe_list=False),
                        CENSO_TIMEOUT_HOURS,
                        cpu_list=_recover_cpu_pool,
                    )
                    # Buscar el ensemble de salida (mismo esquema que abajo)
                    _fixed = [
                        "censo_ensemble_sorted.xyz",
                        "censo_ensemble.xyz",
                        "crest_conformers_sorted.xyz",
                        "crest_conformers.xyz.sorted",
                    ]
                    for _name in _fixed:
                        _p = censo_dir / _name
                        if _p.exists() and _p.stat().st_size > 0 and _p != local_input:
                            if not safe_copy(_p, output_xyz):
                                break
                            try:
                                _n = len(read(str(output_xyz), index=":"))
                            except Exception:
                                _n = "?"
                            log.info(
                                f"[CENSO-RECOVER]  ✓  ensemble recuperado con --reload "
                                f"→ {output_xyz.name} ({_n} conf(s))"
                            )
                            return output_xyz
                    # No encontró ensemble → recovery falló; continuar al lanzamiento limpio
                    log.warning(
                        f"[CENSO-RECOVER]  ⚠  --reload no produjo un ensemble "
                        f"reconocible (returncode={_proc.returncode}) — "
                        f"revisa censo_recover.log. Se limpiará la carpeta y CENSO "
                        f"correrá desde cero."
                    )
                except subprocess.TimeoutExpired:
                    log.error("[CENSO-RECOVER]  ✗  TIMEOUT en --reload")
                except Exception as _e:
                    log.error(f"[CENSO-RECOVER]  ✗  error en --reload: {_e}")
                finally:
                    _cpu_lock_release()

        # ── Limpieza de subcarpetas CONF* corruptas / no recuperables ───────
        # Solo se borran las carpetas CONF* individuales; el directorio
        # censo_dir (y archivos sueltos como censo.log, JSON, input XYZ)
        # se conservan para diagnóstico.
        log.warning(
            f"[CENSO]  🗑  eliminando {len(_diag['conf_dirs'])} subcarpetas "
            f"CONF* corruptas/incompletas en {censo_dir.name}/ para relanzar "
            f"CENSO desde cero"
        )
        for _cd in _diag["conf_dirs"]:
            try:
                shutil.rmtree(_cd)
                log.debug(f"[CENSO]  rm -rf {_cd.name}")
            except OSError as _e:
                log.warning(f"[CENSO]  no se pudo eliminar {_cd.name}: {_e}")
        log.info("[CENSO]  ✓  limpieza completada — relanzando CENSO limpio")

    # ── Intento de recuperación desde 1_SCREENING.out existente ─────────────
    # Si CENSO ya corrió manualmente (o en un run previo) y dejó el archivo
    # de resultados del pre-screening, parsear directamente en vez de
    # relanzar CENSO (que puede no estar en PATH o estar en otro entorno).
    # Busca en las ubicaciones más habituales que CENSO 2.x / 3.x produce.
    _screening_candidates = [
        censo_dir / "1_SCREENING.out",
        censo_dir / "CENSO" / "1_SCREENING.out",
        censo_dir / "prescreening" / "1_SCREENING.out",
    ]
    for _sc in _screening_candidates:
        _conf_ids = parse_censo_screening_out(_sc)
        if _conf_ids:
            log.info(
                f"[CENSO]  📄  encontrado {_sc.relative_to(censo_dir)} "
                f"con {len(_conf_ids)} conf(s) — reconstruyendo ensemble "
                f"filtrado sin relanzar CENSO"
            )
            # Fuentes de geometría en orden de prioridad:
            #   1. conformers_xyz original de CREST (comentarios = solo energía)
            #   2. 1_SCREENING.xyz generado por CENSO (comentarios = CONF#)
            #      CENSO renombra los frames con CONF# al leerlos; es el único
            #      archivo que contiene los tokens CONF# como comentario XYZ.
            _geom_sources = [conformers_xyz]
            for _screening_xyz_candidate in [
                _sc.parent / "1_SCREENING.xyz",
                _sc.parent.parent / "1_SCREENING.xyz",
            ]:
                if _screening_xyz_candidate.exists() and _screening_xyz_candidate.stat().st_size > 0:
                    _geom_sources.append(_screening_xyz_candidate)

            _rebuilt = False
            for _geom_src in _geom_sources:
                if extract_conformers_by_id(_geom_src, _conf_ids, output_xyz, log):
                    log.info(
                        f"[CENSO]  ✓  ensemble reconstruido desde {_sc.name} "
                        f"(geometrías: {_geom_src.name}) "
                        f"→ {output_xyz.name} ({len(_conf_ids)} conf(s))"
                    )
                    _rebuilt = True
                    break
                log.debug(
                    f"[CENSO]  fuente {_geom_src.name} no contiene los CONF# — "
                    f"probando siguiente fuente de geometrías"
                )

            if _rebuilt:
                return output_xyz
            log.warning(
                f"[CENSO]  ⚠  no se pudo reconstruir desde {_sc.name} "
                f"— continuando con intento de relanzar CENSO"
            )
            break
    # ─────────────────────────────────────────────────────────────────────────

    # CENSO escribe en cwd; copiar el ensemble de entrada ahí y referenciarlo
    # por nombre relativo evita ambigüedad sobre rutas absolutas en -i.
    local_input = censo_dir / conformers_xyz.name
    if not safe_copy(conformers_xyz, local_input):
        log.warning(f"[CENSO] no se pudo copiar {conformers_xyz.name} a {censo_dir} "
                    f"— usando ensemble CREST sin filtrar")
        return conformers_xyz

    log.info(f"[CENSO]  ▶  lanzando pre-screening sobre {conformers_xyz.name} "
             f"(carga={charge}, T={temperature:.2f} K)"
             + ("  [nohup]" if USE_NOHUP else ""))
    t0 = time.time()
    try:
        # ago 2026: reservar el pool ANTES de armar --maxcores/--omp-min —
        # ver _reconcile_censo_budget() arriba. Sin esto, --maxcores usaba
        # siempre NPROCS_CREST nominal aunque censo_cpu_pool viniera más
        # chico (otra instancia de pipeline.py con PUs ya reservados),
        # causando que CENSO programe más sub-jobs concurrentes / más
        # hilos por sub-job de los que el taskset confinado realmente
        # soporta.
        censo_cpu_pool = _acquire_cpu_pool(NPROCS_CREST, log, "CENSO")
        eff_maxcores, eff_omp_min = _reconcile_censo_budget(
            NPROCS_CREST, CENSO_OMP_MIN, censo_cpu_pool, log, "CENSO"
        )
        if eff_maxcores == 0:
            _cpu_lock_release()
            log.error(f"[CENSO]  ✗  sin CPUs disponibles — "
                      f"usando ensemble CREST sin filtrar")
            return conformers_xyz

        cmd = [
            str(CENSO_BIN),
            "--prescreening",
            "--screening",               # FIX: faltaba — sin esto CENSO 3.x corre
                                          # SOLO prescreening y nunca llega a
                                          # screening (1_SCREENING.out nunca existe).
            "-i", local_input.name,
            "-c", str(charge),
            "-u", "0",                  # MULT global siempre 1 → singlete, capa cerrada
            "--solvent", "water",
            # prog/func/basis se leen de ~/.censo2rc ([prescreening] prog=xtb, func=gfn2, smgsolv=alpb)
            # --prog fue eliminado en CENSO 2.x; ya no se pasa por CLI
            "--maxcores", str(eff_maxcores),    # ago 2026: reconciliado contra
                                                 # len(censo_cpu_pool) real, no
                                                 # NPROCS_CREST nominal — ver
                                                 # _reconcile_censo_budget()
            "--omp-min", str(eff_omp_min),      # ídem — antes CENSO_OMP_MIN fijo,
                                                 # ahora acotado al pool real
            "-T", f"{temperature:.2f}",             # temperatura configurable (default 274.15 K = 1°C Antártico)
        ]
        # 2026-07: se abandonó el intento de repartir manualmente los
        # sub-jobs internos de CENSO (MAAQC_CENSO_SLOTS + parche local en
        # orca_processor.py) — seguía mostrando colisión parcial entre
        # sub-jobs (Cpus_allowed_list solapados) pese a todas las capas de
        # fix anteriores. En vez de seguir forzando la asignación a mano,
        # censo_cpu_pool solo fija el LÍMITE EXTERIOR del árbol completo
        # de CENSO vía taskset (esto sí es necesario: coordina con otras
        # instancias de pipeline.py corriendo en paralelo, vía
        # CPU_LOCK_FILE). Dentro de ese pool, el reparto entre los propios
        # sub-jobs de CENSO (Dask + xtb --parallel) queda 100% a cargo del
        # scheduler del kernel — set_pe_list=False evita que se le imponga
        # un MAAQC_PE_LIST compartido entre sub-jobs (esa era la causa raíz
        # original del solapamiento, ver Bug #1 junto a build_omp_env()).
        #
        # ago 2026 (segunda vuelta): censo_threads_per_job() ya usaba
        # hilos-por-sub-job en vez del pool completo, pero seguía
        # calculándolo sobre el valor NOMINAL pedido (NPROCS_CREST),
        # ignorando que censo_cpu_pool puede venir más chico si otra
        # instancia de pipeline.py ya reservó PUs — ver
        # _reconcile_censo_budget() arriba, que ahora fija
        # (eff_maxcores, eff_omp_min) al tamaño REAL del pool antes de
        # armar --maxcores/--omp-min/OMP_NUM_THREADS.
        proc = _run_censo_locked(
            cmd, censo_dir, censo_dir / "censo.log",
            build_omp_env(censo_threads_per_job(eff_maxcores, eff_omp_min),
                           cpu_list=censo_cpu_pool,
                           set_pe_list=False),
            CENSO_TIMEOUT_HOURS,
            cpu_list=censo_cpu_pool,
        )
        returncode = proc.returncode
    except subprocess.TimeoutExpired:
        log.error(f"[CENSO]  ✗  TIMEOUT ({CENSO_TIMEOUT_HOURS} h) — "
                  f"usando ensemble CREST sin filtrar")
        return conformers_xyz
    except FileNotFoundError as e:
        if "taskset" in str(e):
            log.error(f"  ✗  no se encontró 'taskset' en PATH — instálalo "
                      f"(paquete util-linux) o llama a stage_censo_prescreen "
                      f"sin cpu pinning")
        else:
            log.error(f"[CENSO]  ✗  no se encontró el ejecutable ({CENSO_BIN}) — "
                      f"¿'conda activate sciencebitch'? Usando ensemble CREST sin filtrar")
        return conformers_xyz
    except Exception as e:
        log.error(f"[CENSO]  ✗  error inesperado: {e} — usando ensemble CREST sin filtrar")
        return conformers_xyz
    finally:
        _cpu_lock_release()

    elapsed = (time.time() - t0) / 60

    # ── Diagnóstico temprano: parsear censo.log para errores de config ──────
    censo_log_path = censo_dir / "censo.log"
    _censo_log_snippet = ""
    if censo_log_path.exists():
        try:
            _tail = censo_log_path.read_text(errors="replace").splitlines()
            # Buscar indicadores de error de configuración de CENSO 3.x
            _error_lines = [l for l in _tail if any(kw in l for kw in (
                "validation error", "Path for", "is required", "is not set",
                "UserWarning", "No context found",
            ))]
            if _error_lines:
                _censo_log_snippet = " | ".join(_error_lines[-4:])
                log.warning(
                    f"[CENSO]  ✗  CENSO terminó con errores de configuración "
                    f"(returncode={returncode}) — posible causa: ~/.censo2rc incompleto "
                    f"o falta 'tm'/'orca' en [paths]. Detalle: {_censo_log_snippet}"
                )
        except Exception:
            pass

    # ── Reconstrucción preferente desde el .out con energías (Gtot) ─────────
    # CENSO escribe las energías reales (Gtot/ΔGtot/Boltzmann weight) en
    # archivos .out (1_SCREENING.out si corrió screening, o 0_PRESCREENING.out
    # si solo corrió prescreening). Los .xyz que CENSO deja (censo_ensemble_*.xyz,
    # 1_SCREENING.xyz, etc.) tienen comentarios "CONF7" SIN energía numérica,
    # por lo que select_nconf_auto() no puede leerlas si se usan directamente.
    # Por eso, antes de copiar cualquier .xyz crudo, intentamos reconstruir el
    # ensemble vía parse_censo_screening_out() + extract_conformers_by_id(),
    # que sí preserva el orden por Gtot y permite recuperar energías después.
    _post_out_candidates = [
        censo_dir / "1_SCREENING.out",
        censo_dir / "CENSO" / "1_SCREENING.out",
        censo_dir / "prescreening" / "1_SCREENING.out",
        censo_dir / "0_PRESCREENING.out",
        censo_dir / "CENSO" / "0_PRESCREENING.out",
        censo_dir / "prescreening" / "0_PRESCREENING.out",
    ]
    for _po in _post_out_candidates:
        _conf_ids = parse_censo_screening_out(_po)
        if not _conf_ids:
            continue
        _geom_sources = [conformers_xyz]
        for _xyz_name in ("1_SCREENING.xyz", "0_PRESCREENING.xyz"):
            for _cand in (_po.parent / _xyz_name, _po.parent.parent / _xyz_name):
                if _cand.exists() and _cand.stat().st_size > 0:
                    _geom_sources.append(_cand)
        for _geom_src in _geom_sources:
            if extract_conformers_by_id(_geom_src, _conf_ids, output_xyz, log):
                log.info(
                    f"[CENSO]  ✓  pre-screening terminó en {elapsed:.1f} min — "
                    f"{len(_conf_ids)} confórmero(s) sobreviven, ensemble "
                    f"reconstruido desde {_po.relative_to(censo_dir)} "
                    f"(geometrías: {_geom_src.name}, ordenado por Gtot)"
                )
                return output_xyz
        log.debug(
            f"[CENSO]  {_po.name} parseado pero no se pudo reconstruir "
            f"geometría — probando candidatos .xyz crudos"
        )

    # ── Búsqueda robusta de output: nombres fijos + glob de respaldo ─────────
    # CENSO 2/3 puede escribir en censo_dir directamente o en subdirectorios
    # (CENSO/, prescreening/, etc.) con nombres que varían entre versiones.
    # NOTA: este fallback se usa solo si no había .out con energías parseables
    # arriba; estos .xyz traen comentarios "CONF#" sin energía numérica.
    _fixed_candidates = [
        "censo_ensemble_sorted.xyz",
        "censo_ensemble.xyz",
        "crest_conformers_sorted.xyz",
        "crest_conformers.xyz.sorted",
    ]
    # Candidatos fijos primero (orden de prioridad)
    _candidate_paths = [censo_dir / c for c in _fixed_candidates]
    # Luego cualquier .xyz en subdirectorios que no sea el input copiado
    _candidate_paths += sorted(
        (p for p in censo_dir.rglob("*.xyz")
         if p != local_input and p.name != output_xyz.name),
        key=lambda p: p.stat().st_mtime if p.exists() else 0,
        reverse=True,  # el más reciente primero
    )

    for p in _candidate_paths:
        if not p.exists() or p.stat().st_size == 0:
            continue
        if p == local_input:
            continue  # ignorar el input que copiamos nosotros
        if not safe_copy(p, output_xyz):
            log.warning(f"[CENSO] no se pudo normalizar {p.name} → {output_xyz.name} "
                        f"— usando ensemble CREST sin filtrar")
            return conformers_xyz
        try:
            n = len(read(str(output_xyz), index=":"))
        except Exception:
            n = "?"
        log.info(f"[CENSO]  ✓  pre-screening terminó en {elapsed:.1f} min — "
                 f"{n} confórmero(s) sobreviven ({p.relative_to(censo_dir)})")
        return output_xyz

    # Si llegamos aquí, CENSO no produjo ningún ensemble utilizable
    if not _censo_log_snippet:
        # No hubo error de config detectado → informar genéricamente
        log.warning(
            f"[CENSO]  ✗  no se encontró un ensemble filtrado reconocible "
            f"(returncode={returncode}) en {censo_dir} — revisa {censo_log_path} "
            f"a mano. Usando ensemble CREST sin filtrar."
        )
    return conformers_xyz


# ─────────────────────────────────────────────
#  PARSEO DE 1_SCREENING.out DE CENSO
# ─────────────────────────────────────────────

def parse_censo_screening_out(screening_out: Path) -> list[str] | None:
    """
    Parsea el archivo 1_SCREENING.out de CENSO y devuelve la lista de
    identificadores de confórmeros que sobrevivieron el pre-screening,
    ordenados por Gtot (más estable primero).

    Formato esperado (CENSO 2.x / 3.x):
        CONF#      E (DFT)   ΔGsolv   GmRRHO   Gtot   ΔGtot   Boltzmann weight
        CONF7     ...                           ...    0.00      41.71
        CONF290   ...                           ...    0.11      34.76
        ...

    Devuelve lista de strings tipo ["CONF7", "CONF290", ...] o None si
    no se pudo parsear ninguna línea de datos.
    """
    if not screening_out.exists():
        return None

    conf_ids: list[tuple[float, str]] = []  # (Gtot, conf_id)

    try:
        lines = screening_out.read_text(errors="replace").splitlines()
    except Exception:
        return None

    for line in lines:
        line = line.strip()
        # Línea de datos: empieza con CONF seguido de dígitos
        m = re.match(
            r'^(CONF\d+)\s+'           # CONF#
            r'[\-\d.]+\s+'             # E (DFT)
            r'[\-\d.]+\s+'             # ΔGsolv
            r'[\-\d.]+\s+'             # GmRRHO
            r'([\-\d.]+)\s+'           # Gtot  ← usamos este para ordenar
            r'([\-\d.]+)\s+'           # ΔGtot
            r'[\d.]+',                 # Boltzmann weight
            line,
        )
        if m:
            conf_id = m.group(1)       # "CONF7"
            gtot    = float(m.group(2))
            conf_ids.append((gtot, conf_id))

    if not conf_ids:
        return None

    # Ordenar por Gtot (ya vienen ordenados en el .out, pero por seguridad)
    conf_ids.sort(key=lambda x: x[0])
    return [cid for _, cid in conf_ids]


def extract_conformers_by_id(
    source_xyz: Path,
    conf_ids: list[str],
    output_xyz: Path,
    log: logging.Logger,
) -> bool:
    """
    Extrae del multi-XYZ de CREST los frames cuyo comentario contiene
    uno de los identificadores en conf_ids (e.g. "CONF7"), y los escribe
    en output_xyz en el mismo orden que conf_ids.

    CREST numera los frames secuencialmente en el comentario de cada bloque
    XYZ (ej: "CONF7  -1218.034197"). Busca el primer token del comentario
    como clave.

    Devuelve True si se escribió al menos un confórmero.
    """
    try:
        all_frames = read(str(source_xyz), index=":")
        if not isinstance(all_frames, list):
            all_frames = [all_frames]
    except Exception as e:
        log.warning(f"[CENSO] no se pudo leer {source_xyz}: {e}")
        return False

    # Construir índice: conf_id → atoms
    # ASE solo popula atoms.info para formato extXYZ (key=value en comentario).
    # Para XYZ simples (comentario = "CONF7" o "-76.301"), atoms.info queda
    # vacío y hay que leer el comentario directamente del archivo de texto.
    frame_index: dict[str, object] = {}
    try:
        _raw_lines = source_xyz.read_text(errors="replace").splitlines()
        _n_atoms   = len(all_frames[0]) if all_frames else 0
        _fi = 0   # índice en all_frames
        _li = 0   # índice en _raw_lines
        while _li < len(_raw_lines) and _fi < len(all_frames):
            _line = _raw_lines[_li].strip()
            if _line.isdigit() and int(_line) == _n_atoms:
                _comment = _raw_lines[_li + 1].strip() if _li + 1 < len(_raw_lines) else ""
                _token   = _comment.split()[0] if _comment.split() else ""
                if _token.startswith("CONF"):
                    frame_index[_token] = all_frames[_fi]
                _fi += 1
                _li += _n_atoms + 2
            else:
                _li += 1
    except Exception:
        # Fallback: intentar con atoms.info (funciona si es extXYZ)
        for atoms in all_frames:
            comment = (
                atoms.info.get("comment", "")
                or atoms.info.get("description", "")
            )
            token = str(comment).split()[0] if comment else ""
            if token.startswith("CONF"):
                frame_index[token] = atoms

    selected = []
    missing  = []
    for cid in conf_ids:
        if cid in frame_index:
            selected.append(frame_index[cid])
        else:
            missing.append(cid)

    if missing:
        log.warning(
            f"[CENSO]  ⚠  {len(missing)} conf(s) del screening no encontrados "
            f"en {source_xyz.name}: {missing[:5]}{'…' if len(missing) > 5 else ''}"
        )

    if not selected:
        return False

    try:
        write(str(output_xyz), selected)
        return True
    except Exception as e:
        log.warning(f"[CENSO] error escribiendo ensemble filtrado: {e}")
        return False


# ─────────────────────────────────────────────
#  CENSO STANDALONE (--censo-standalone)
# ─────────────────────────────────────────────
#
# Motivación (ago 2026): Kmilo seguía invocando `censo` pelado a mano en la
# terminal para correr un ensemble suelto (fuera del flujo CREST->pipeline),
# perdiendo con eso TODO el aislamiento de cores que stage_censo_prescreen()
# ya tiene resuelto (taskset -c sobre un pool que evita mezclar SMT/E-core/
# LP-core, OMP_NUM_THREADS=maxcores + MKL/OPENBLAS_NUM_THREADS=1, lockfile
# para no relanzar sobre un CENSO que ya sigue vivo). Resultado medido: los
# mismos síntomas que ya se habían diagnosticado y arreglado una vez —
# varios sub-jobs de CENSO compitiendo por los mismos cores físicos en vez
# de repartirse el pool completo, con throttling térmico de por medio.
#
# run_censo_standalone() NO reimplementa nada nuevo: arma el mismo `cmd` y
# usa exactamente el mismo trío _acquire_cpu_pool() + build_omp_env(...,
# set_pe_list=False) + _run_censo_locked() que ya usan los dos call-sites
# reales de CENSO (lanzamiento limpio y --reload, ver stage_censo_prescreen
# más arriba). A diferencia de stage_censo_prescreen(), a propósito NO
# hace diagnóstico/recovery de corridas parciales previas — es un modo para
# invocaciones manuales puntuales sobre un ensemble ya filtrado o crudo, no
# para el flujo automatizado completo del pipeline (para eso, seguir usando
# el pipeline normal, que sí hace recovery).
def run_censo_standalone(input_xyz: Path, charge: int, out_dir: Path,
                          temperature: float, maxcores: int,
                          log: logging.Logger,
                          omp_min: int = CENSO_OMP_MIN) -> bool:
    """
    Corre CENSO (--prescreening --screening) sobre `input_xyz` dentro de
    `out_dir`, con el mismo aislamiento de cores que usa el pipeline
    completo para CENSO. Devuelve True si CENSO terminó con returncode 0
    (independientemente de si se pudo localizar un ensemble de salida
    reconocible — eso solo afecta el nivel del log final).
    """
    owner_pid = _censo_lock_owner(out_dir, log)
    if owner_pid is not None:
        log.error(
            f"[CENSO-STANDALONE]  ✗  ya hay un proceso CENSO activo sobre "
            f"{out_dir} (PID {owner_pid}) — no se relanza para no pisarle "
            f"el trabajo (este era exactamente el bug de relanzamientos "
            f"duplicados ya diagnosticado antes). Espera a que termine o "
            f"revisa si el lock quedó huérfano."
        )
        return False

    local_input = out_dir / input_xyz.name
    if not safe_copy(input_xyz, local_input):
        log.error(f"[CENSO-STANDALONE]  ✗  no se pudo copiar {input_xyz.name} "
                  f"a {out_dir}")
        return False

    # ago 2026: reservar y reconciliar el pool ANTES de armar --maxcores/
    # --omp-min — ver _reconcile_censo_budget(). Antes, --maxcores usaba
    # siempre el `maxcores` pedido aunque cpu_pool viniera más chico (otra
    # instancia de pipeline.py con PUs ya reservados).
    cpu_pool = _acquire_cpu_pool(maxcores, log, "CENSO-STANDALONE")
    eff_maxcores, eff_omp_min = _reconcile_censo_budget(
        maxcores, omp_min, cpu_pool, log, "CENSO-STANDALONE"
    )
    if eff_maxcores == 0:
        _cpu_lock_release()
        log.error(f"[CENSO-STANDALONE]  ✗  sin CPUs disponibles — abortando")
        return False

    cmd = [
        str(CENSO_BIN),
        "--prescreening", "--screening",   # sin --screening, CENSO 3.x se
                                            # detiene tras prescreening
        "-i", local_input.name,
        "-c", str(charge),
        "-u", "0",                        # MULT global siempre 1 (singlete)
        "--solvent", "water",
        "--maxcores", str(eff_maxcores),  # ago 2026: reconciliado contra
                                           # len(cpu_pool) real, no el
                                           # `maxcores` nominal pedido
        "--omp-min", str(eff_omp_min),    # ídem
        "-T", f"{temperature:.2f}",
    ]

    log.info(f"[CENSO-STANDALONE]  ▶  lanzando sobre {input_xyz.name} "
             f"(carga={charge}, T={temperature:.2f} K, maxcores={eff_maxcores}, "
             f"omp-min={eff_omp_min}) en {out_dir}"
             + ("  [nohup]" if USE_NOHUP else ""))

    try:
        # OMP_NUM_THREADS = hilos-por-sub-job (censo_threads_per_job), NO
        # maxcores — ver CENSO_OMP_MIN arriba. Con maxcores el pool
        # completo, cada uno de los N sub-jobs paralelos que CENSO corre a
        # la vez pedía el pool completo para sí mismo (CENSO congela
        # ENVIRON una sola vez y no reparte OMP_NUM_THREADS entre
        # sub-jobs), sobre-suscribiendo el pool confinado por taskset.
        proc = _run_censo_locked(
            cmd, out_dir, out_dir / "censo.log",
            build_omp_env(censo_threads_per_job(eff_maxcores, eff_omp_min),
                           cpu_list=cpu_pool, set_pe_list=False),
            CENSO_TIMEOUT_HOURS,
            cpu_list=cpu_pool,
        )
    except subprocess.TimeoutExpired:
        log.error(f"[CENSO-STANDALONE]  ✗  TIMEOUT ({CENSO_TIMEOUT_HOURS} h)")
        return False
    except FileNotFoundError as e:
        if "taskset" in str(e):
            log.error("[CENSO-STANDALONE]  ✗  no se encontró 'taskset' en PATH "
                      "(paquete util-linux)")
        else:
            log.error(f"[CENSO-STANDALONE]  ✗  no se encontró el ejecutable "
                      f"({CENSO_BIN}) — ¿'conda activate sciencebitch'?")
        return False
    finally:
        _cpu_lock_release()

    if proc.returncode != 0:
        log.error(f"[CENSO-STANDALONE]  ✗  CENSO terminó con returncode="
                  f"{proc.returncode} — revisa {out_dir / 'censo.log'}")
        return False

    for _name in ("censo_ensemble_sorted.xyz", "censo_ensemble.xyz",
                  "crest_conformers_sorted.xyz", "crest_conformers.xyz.sorted"):
        _p = out_dir / _name
        if _p.exists() and _p.stat().st_size > 0:
            log.info(f"[CENSO-STANDALONE]  ✓  ensemble filtrado → {_p}")
            return True

    log.warning(
        f"[CENSO-STANDALONE]  ⚠  CENSO terminó (returncode=0) pero no se "
        f"encontró un ensemble de salida reconocible en {out_dir} — revisa "
        f"censo.log / 1_SCREENING.out a mano."
    )
    return True


def _run_censo_standalone_cli(args) -> None:
    """Punto de entrada de --censo-standalone: valida args, prepara out_dir
    y logger, y delega en run_censo_standalone(). Sale del proceso con
    código 0/1 según éxito."""
    input_xyz = args.censo_standalone.resolve()
    if not input_xyz.is_file():
        print(f"ERROR: no existe el ensemble {input_xyz}")
        raise SystemExit(1)

    out_dir = (args.out.resolve() if args.out
               else input_xyz.parent / f"{input_xyz.stem}_censo")
    out_dir.mkdir(parents=True, exist_ok=True)

    log = setup_logging(out_dir)
    check_required_bins(log)

    charge = resolve_charge(input_xyz.stem, args.charge, log)

    maxcores = args.cores if args.cores else NPROCS_CREST
    ok = run_censo_standalone(
        input_xyz=input_xyz, charge=charge, out_dir=out_dir,
        temperature=args.temperature, maxcores=maxcores, log=log,
        omp_min=args.omp_min,
    )
    raise SystemExit(0 if ok else 1)


# ─────────────────────────────────────────────
#  CENSO MULTI-TEMPERATURA
# ─────────────────────────────────────────────

def stage_censo_multitemp(
    conformers_xyz: Path,
    mol_dir: Path,
    charge: int,
    log: logging.Logger,
    temperatures: list[float],
) -> dict[float, Path]:
    """
    Corre stage_censo_prescreen() para cada temperatura en `temperatures`.

    Cada temperatura produce su propio subdirectorio dentro de
    01b_censo_prescreen/:

        01b_censo_prescreen/
        └── T274.15K/
        │   ├── censo.log
        │   └── censo_ensemble_sorted.xyz
        └── T283.00K/
            ├── censo.log
            └── censo_ensemble_sorted.xyz
        ...

    El nombre del subdirectorio tiene formato "T{temp:.2f}K" redondeado
    a 2 decimales (ej. "T274.15K", "T298.00K").

    Devuelve un dict {temperatura_K → ruta_ensemble}. Si CENSO falla para
    una temperatura se registra la ruta del ensemble CREST sin filtrar
    (fail-open, consistente con stage_censo_prescreen()).

    Los resultados se guardan además en:
        01b_censo_prescreen/multitemp_summary.json
    para facilitar análisis posterior con results_collector.py.
    """
    base_censo_dir = mol_dir / "01b_censo_prescreen"
    base_censo_dir.mkdir(parents=True, exist_ok=True)

    results: dict[float, Path] = {}

    log.info(
        f"[CENSO-MULTITEMP]  ▶  {len(temperatures)} temperaturas: "
        f"{temperatures[0]:.2f} K … {temperatures[-1]:.2f} K"
    )

    for temp in temperatures:
        subdir_name = f"T{temp:.2f}K"
        log.info(f"[CENSO-MULTITEMP]  ── T = {temp:.2f} K ({subdir_name})")

        result_xyz = stage_censo_prescreen(
            conformers_xyz=conformers_xyz,
            mol_dir=mol_dir,
            charge=charge,
            log=log,
            temperature=temp,
            subdir=subdir_name,
        )
        results[temp] = result_xyz

    # ── Resumen JSON ─────────────────────────────────────────────────────────
    summary_path = base_censo_dir / "multitemp_summary.json"
    summary_data = {
        "temperatures_K": temperatures,
        "results": {
            f"{temp:.2f}": str(path)
            for temp, path in results.items()
        },
        "n_temperatures": len(temperatures),
        "ensemble_input": str(conformers_xyz),
    }
    try:
        with open(summary_path, "w") as f:
            json.dump(summary_data, f, indent=2)
        log.info(f"[CENSO-MULTITEMP]  ✓  resumen guardado en {summary_path.name}")
    except OSError as e:
        log.warning(f"[CENSO-MULTITEMP]  ⚠  no se pudo guardar resumen JSON: {e}")

    n_ok = sum(1 for p in results.values() if p != conformers_xyz)
    log.info(
        f"[CENSO-MULTITEMP]  ✓  completado — {n_ok}/{len(temperatures)} "
        f"temperaturas con ensemble filtrado, "
        f"{len(temperatures) - n_ok} usando CREST sin filtrar (fail-open)"
    )

    return results


# ─────────────────────────────────────────────
#  FILTRO BOLTZMANN POST-CENSO
# ─────────────────────────────────────────────

def parse_boltzmann_weights(screening_out: Path) -> dict[str, float]:
    """
    Parsea los pesos de Boltzmann de 1_SCREENING.out de CENSO.

    Devuelve dict {conf_id → peso%} (ej. {"CONF7": 41.71, "CONF290": 34.76}).
    Devuelve dict vacío si el archivo no existe o no se pudo parsear.
    """
    if not screening_out.exists():
        return {}
    result: dict[str, float] = {}
    try:
        for line in screening_out.read_text(errors="replace").splitlines():
            line = line.strip()
            m = re.match(
                r'^(CONF\d+)\s+'
                r'[\-\d.]+\s+'   # E(DFT)
                r'[\-\d.]+\s+'   # ΔGsolv
                r'[\-\d.]+\s+'   # GmRRHO
                r'[\-\d.]+\s+'   # Gtot
                r'[\-\d.]+\s+'   # ΔGtot
                r'([\d.]+)',        # Boltzmann weight %
                line,
            )
            if m:
                result[m.group(1)] = float(m.group(2))
    except Exception:
        return {}
    return result


def select_nconf_boltzmann(
    conformers_xyz: Path,
    screening_out:  Path,
    log:            logging.Logger,
    threshold:      float = 5.0,
    nconf_min:      int   = 3,          # = NCONF_MIN (definido más abajo)
) -> tuple[int, list[str]]:
    """
    Selecciona confórmeros del ensemble post-CENSO con peso Boltzmann ≥ threshold%.

    Lee los pesos de Boltzmann de 1_SCREENING.out (calculados por CENSO sobre
    geometrías GFN2-xTB), filtra por umbral y devuelve (nconf, conf_ids).

    El ensemble censo_ensemble_sorted.xyz puede tener comentarios en formato
    extXYZ ("CONF7=T") o plain ("CONF7"); se extrae el token CONF# en ambos casos.

    Args:
        conformers_xyz : ruta al ensemble post-CENSO (censo_ensemble_sorted.xyz).
        screening_out  : ruta al 1_SCREENING.out con los pesos Boltzmann.
        log            : logger del pipeline.
        threshold      : porcentaje mínimo (default 5.0 → ≥5%).
        nconf_min      : mínimo de confórmeros garantizado (default NCONF_MIN=3).

    Returns:
        (nconf, conf_ids) donde conf_ids es la lista ordenada por peso Boltzmann
        descendente. Si no se puede aplicar el filtro, devuelve (0, []) para que
        el llamador caiga en el modo automático.
    """
    # ── Leer pesos Boltzmann ─────────────────────────────────────────────────
    weights = parse_boltzmann_weights(screening_out)
    if not weights:
        log.warning(
            f"[NCONF-BOLTZMANN] no se pudo leer {screening_out.name} "
            f"— cayendo en select_nconf_auto"
        )
        return 0, []

    # ── Extraer CONF# del ensemble XYZ ──────────────────────────────────────
    # Soporta comentarios extXYZ ("Properties=... CONF7=T pbc=...") y
    # plain ("CONF7" o "CONF7  -1218.03").
    xyz_conf_ids: list[str] = []
    try:
        _lines   = conformers_xyz.read_text(errors="replace").splitlines()
        _n_atoms = None
        _li      = 0
        while _li < len(_lines):
            _l = _lines[_li].strip()
            if _l.isdigit():
                _nat = int(_l)
                if _n_atoms is None:
                    _n_atoms = _nat
                if _nat == _n_atoms and _li + 1 < len(_lines):
                    _comment = _lines[_li + 1].strip()
                    # extXYZ: CONF7=T
                    _m = re.search(r'(CONF\d+)=', _comment)
                    if _m:
                        xyz_conf_ids.append(_m.group(1))
                    else:
                        # plain: primer token que empiece con CONF
                        for _tok in _comment.split():
                            if _tok.startswith("CONF"):
                                xyz_conf_ids.append(_tok)
                                break
                    _li += _nat + 2
                    continue
            _li += 1
    except Exception as e:
        log.warning(
            f"[NCONF-BOLTZMANN] no se pudo leer {conformers_xyz.name}: {e} "
            f"— cayendo en select_nconf_auto"
        )
        return 0, []

    if not xyz_conf_ids:
        log.warning(
            f"[NCONF-BOLTZMANN] no se encontraron tokens CONF# en "
            f"{conformers_xyz.name} — cayendo en select_nconf_auto"
        )
        return 0, []

    # ── Filtrar por umbral Boltzmann ─────────────────────────────────────────
    above  = [(cid, weights[cid]) for cid in xyz_conf_ids if weights.get(cid, 0.0) >= threshold]
    below  = [(cid, weights[cid]) for cid in xyz_conf_ids if 0 < weights.get(cid, 0.0) < threshold]
    no_wt  = [cid for cid in xyz_conf_ids if cid not in weights]

    # Ordenar por peso descendente
    above.sort(key=lambda x: -x[1])

    # Garantizar nconf_min: si el umbral deja < nconf_min, completar con los
    # mejores por debajo del umbral hasta llegar a nconf_min.
    selected = [cid for cid, _ in above]
    if len(selected) < nconf_min:
        below.sort(key=lambda x: -x[1])
        extras = [cid for cid, _ in below[: nconf_min - len(selected)]]
        selected += extras
        log.warning(
            f"[NCONF-BOLTZMANN] solo {len(above)} conf(s) ≥ {threshold}% — "
            f"añadiendo {len(extras)} conf(s) extra para garantizar "
            f"nconf_min={nconf_min}: {extras}"
        )

    nconf = len(selected)

    # ── Log resumen ──────────────────────────────────────────────────────────
    above_str = ", ".join(
        f"{cid}({weights[cid]:.1f}%)" for cid, _ in above
    )
    log.info(
        f"[NCONF-BOLTZMANN] umbral={threshold}% → "
        f"{len(above)} conf(s) seleccionados: {above_str}"
    )
    if below:
        log.info(
            f"[NCONF-BOLTZMANN] {len(below)} conf(s) por debajo del umbral "
            f"(mayor: {below[0][0]} {below[0][1]:.2f}%) — descartados"
        )
    if no_wt:
        log.warning(
            f"[NCONF-BOLTZMANN] {len(no_wt)} conf(s) sin peso en "
            f"{screening_out.name}: {no_wt} — descartados"
        )
    log.info(f"[NCONF-BOLTZMANN] → nconf={nconf} (nconf_min={nconf_min})")
    return nconf, selected


# ─────────────────────────────────────────────
#  AUTO-SELECCIÓN DE CONFÓRMEROS PARA OPT
# ─────────────────────────────────────────────

# Límites absolutos (nunca se violan)
NCONF_MIN = 3
NCONF_MAX = 20

# Energías de GFN2-xTB en kcal/mol para definir ventanas conformacionales
_WINDOW_SMALL  = 2.0   # paisaje plano  → cobertura generosa
_WINDOW_MEDIUM = 3.0   # rango de decisión intermedio
_WINDOW_LARGE  = 5.0   # paisaje amplio → cota conservadora


def select_nconf_auto(conformers_xyz: Path, log: logging.Logger) -> int:
    """
    Decide automáticamente cuántos confórmeros pasan a OPT a partir del
    panorama conformacional del ensemble post-CENSO (energías GFN2-xTB o
    relativas al mínimo si no se pueden extraer).

    Heurística por número de confórmeros disponibles y spread de energías:

        n_avail ≤ NCONF_MIN            → usa todos (= NCONF_MIN)
        n_avail ≤ NCONF_MAX            → usa todos (el ensemble ya es pequeño)
        n_avail > NCONF_MAX:
            spread ≤ _WINDOW_SMALL     → 5  (paisaje muy plano, muchos conf. equivalentes)
            spread ≤ _WINDOW_MEDIUM    → 10 (paisaje moderado)
            spread ≤ _WINDOW_LARGE     → 15 (paisaje amplio)
            spread > _WINDOW_LARGE     → 20 (máximo, paisaje muy rugoso)

    El spread se calcula con las energías del comentario de cada frame XYZ
    (campo típico de CREST/CENSO). Si no se pueden leer, se usa n_avail como
    proxy y se elige la cota media conservadora.

    Siempre devuelve un valor en [NCONF_MIN, NCONF_MAX].
    """
    # ── Leer ensemble ────────────────────────────────────────────────────────
    try:
        all_conf = read(str(conformers_xyz), index=":")
        if not isinstance(all_conf, list):
            all_conf = [all_conf]
        n_avail = len(all_conf)
    except Exception as e:
        log.warning(f"[NCONF-AUTO] no se pudo leer {conformers_xyz}: {e} "
                    f"— usando nconf={NCONF_MIN}")
        return NCONF_MIN

    if n_avail == 0:
        log.warning("[NCONF-AUTO] ensemble vacío — usando nconf={NCONF_MIN}")
        return NCONF_MIN

    # Caso trivial: todos los disponibles caben dentro del límite
    if n_avail <= NCONF_MAX:
        nconf = max(NCONF_MIN, n_avail)
        log.info(f"[NCONF-AUTO] {n_avail} confórmeros disponibles ≤ {NCONF_MAX} "
                 f"— todos pasan a OPT → nconf={nconf}")
        return nconf

    # ── Extraer energías del comentario de cada frame ────────────────────────
    # CREST escribe la energía en Hartree en la línea de comentario (línea 2
    # de cada bloque XYZ). ASE solo la expone en atoms.info para formato
    # extXYZ (key=value); para XYZ simples hay que leer el texto directamente.
    #
    # El ensemble post-CENSO (censo_ensemble_sorted.xyz) tiene comentarios
    # "CONF7", "CONF1", etc. (sin valor numérico): en ese caso no se pueden
    # extraer energías y n_avail se usa como proxy (ver rama else abajo).
    energies_hartree: list[float] = []
    try:
        _raw_lines = conformers_xyz.read_text(errors="replace").splitlines()
        _n_at = len(all_conf[0]) if all_conf else 0
        _fi = 0
        _li = 0
        while _li < len(_raw_lines) and _fi < len(all_conf):
            _l = _raw_lines[_li].strip()
            if _l.isdigit() and int(_l) == _n_at:
                _comment = _raw_lines[_li + 1].strip() if _li + 1 < len(_raw_lines) else ""
                # Buscar primer token numérico en el comentario
                for _tok in _comment.split():
                    try:
                        energies_hartree.append(float(_tok))
                        break
                    except ValueError:
                        continue
                _fi += 1
                _li += _n_at + 2
            else:
                _li += 1
    except Exception:
        # Fallback: atoms.info (funciona para extXYZ)
        for atoms in all_conf:
            e = atoms.info.get("energy") or atoms.info.get("Energy")
            if e is None:
                comment = atoms.info.get("comment", "") or atoms.info.get("description", "")
                for token in str(comment).split():
                    try:
                        e = float(token)
                        break
                    except ValueError:
                        continue
            if e is not None:
                try:
                    energies_hartree.append(float(e))
                except (TypeError, ValueError):
                    pass

    # ── Calcular spread y decidir nconf ─────────────────────────────────────
    if len(energies_hartree) >= 2:
        e_min = min(energies_hartree)
        e_max = max(energies_hartree)
        spread_kcal = (e_max - e_min) * HARTREE_TO_KCAL

        if spread_kcal <= _WINDOW_SMALL:
            nconf = 5
            tier  = f"spread={spread_kcal:.2f} kcal/mol ≤ {_WINDOW_SMALL} (paisaje plano)"
        elif spread_kcal <= _WINDOW_MEDIUM:
            nconf = 10
            tier  = f"spread={spread_kcal:.2f} kcal/mol ≤ {_WINDOW_MEDIUM} (moderado)"
        elif spread_kcal <= _WINDOW_LARGE:
            nconf = 15
            tier  = f"spread={spread_kcal:.2f} kcal/mol ≤ {_WINDOW_LARGE} (amplio)"
        else:
            nconf = 20
            tier  = f"spread={spread_kcal:.2f} kcal/mol > {_WINDOW_LARGE} (rugoso)"

        log.info(
            f"[NCONF-AUTO] {n_avail} confórmeros disponibles — {tier} "
            f"→ nconf={nconf} (límites: [{NCONF_MIN}, {NCONF_MAX}])"
        )
    else:
        # Sin energías legibles: si el ensemble ya es suficientemente pequeño
        # (viene de CENSO pre-screening), tomar el máximo permitido o todos.
        # Heurística logarítmica solo para ensembles grandes sin energías.
        import math
        if n_avail <= NCONF_MAX * 2:
            # Ensemble post-CENSO sin energías numéricas (comentarios = CONF#):
            # tomar min(n_avail, NCONF_MAX) — no tiene sentido restringir más.
            nconf = min(n_avail, NCONF_MAX)
            log.warning(
                f"[NCONF-AUTO] {n_avail} confórmeros disponibles — energías no legibles "
                f"(ensemble post-CENSO); tomando nconf={nconf} (min(n_avail, NCONF_MAX))"
            )
        else:
            ratio = math.log(n_avail) / math.log(100)  # 100 conf ≈ nconf_max
            nconf = int(NCONF_MIN + (NCONF_MAX - NCONF_MIN) * min(ratio, 1.0))
            log.warning(
                f"[NCONF-AUTO] {n_avail} confórmeros disponibles — energías no legibles; "
                f"usando heurística de tamaño → nconf={nconf}"
            )

    return max(NCONF_MIN, min(NCONF_MAX, nconf))


def _extract_conf_ids_from_xyz(source_xyz: Path) -> dict[int, str]:
    """
    Lee el multi-XYZ de CREST/CENSO y devuelve un dict {índice_0based → CONF#}.

    Soporta comentarios extXYZ ("CONF7=T ...") y plain ("CONF7  -1218.034197").
    Si un frame no tiene token CONF# en el comentario, se omite del dict.

    Usado por stage_opt para propagar la identificación
    CREST original (CONF39, CONF29, etc.) a los directorios de trabajo ORCA,
    evitando el re-etiquetado conf_001, conf_002, ... que borra la trazabilidad.
    """
    result: dict[int, str] = {}
    try:
        raw_lines = source_xyz.read_text(errors="replace").splitlines()
    except OSError:
        return result

    frame_idx = 0
    line_idx  = 0
    n_atoms: Optional[int] = None

    while line_idx < len(raw_lines):
        stripped = raw_lines[line_idx].strip()
        if stripped.isdigit():
            nat = int(stripped)
            if n_atoms is None:
                n_atoms = nat
            if nat == n_atoms and line_idx + 1 < len(raw_lines):
                comment = raw_lines[line_idx + 1].strip()
                # extXYZ: "CONF7=T ..."
                m = re.search(r"(CONF\d+)=", comment)
                if m:
                    result[frame_idx] = m.group(1)
                else:
                    # plain: primer token que empiece con CONF
                    for tok in comment.split():
                        if tok.startswith("CONF"):
                            result[frame_idx] = tok
                            break
                frame_idx += 1
                line_idx  += nat + 2
                continue
        line_idx += 1

    return result


def stage_opt(conformers_xyz: Path, mol_dir: Path,
              charge: int, nconf: int, log: logging.Logger,
              conf_ids: Optional[list[str]] = None,
              expected_natoms: Optional[int] = None) -> list[Path]:
    """
    Etapa 02_opt: optimización DFT de cada confórmero.
    `nconf` viene de select_nconf_auto() o select_nconf_boltzmann().

    Args:
        conf_ids : lista explícita de CONF# a optimizar (modo --boltzmann-threshold).
                   Si se especifica, se ignora nconf y se seleccionan los frames
                   correspondientes del ensemble. Si es None, se toman los primeros
                   nconf frames (comportamiento clásico).
        expected_natoms : número de átomos que DEBE tener cada confórmero,
                   derivado de la molécula base (mol_xyz) que dio origen a
                   este mol_dir/charge. Si se especifica y no coincide con
                   el ensemble leído, se aborta con RuntimeError en vez de
                   seguir adelante silenciosamente.

                   Esta validación existe porque conformers_xyz, mol_dir y
                   charge llegan aquí como tres valores independientes: nada
                   impide en el nivel de orquestación (p.ej. al iterar sobre
                   varios estados de protonación A/B/C/D/R) pasar el
                   ensemble de un estado con el charge de otro. Ni
                   read_conformers() ni el modo Boltzmann verifican eso por
                   sí solos — ambos toman lo que se les da y confían en que
                   ya está alineado. expected_natoms es el candado que
                   faltaba: si no coincide, es evidencia directa de un cruce
                   de estado/charge, y ORCA fallaría en el SCF (o, peor,
                   convergería a una geometría sin sentido químico) en vez
                   de avisar aquí.

    Corre OPT_BATCH_SIZE confórmeros en paralelo (ThreadPoolExecutor),
    cada uno con NPROCS cores y MAXCORE MB. El paralelismo entre confórmeros
    es perfecto (no hay comunicación entre jobs).

    Devuelve lista de directorios exitosos.
    """
    if conf_ids:
        # Modo Boltzmann: extraer frames por CONF# desde el ensemble
        # Soporta comentarios extXYZ ("CONF7=T") y plain ("CONF7")
        _selection_mode = f"Boltzmann ({len(conf_ids)} conf(s): {', '.join(conf_ids)})"
        log.info(f"[OPT] optimizando {len(conf_ids)} confórmeros ({_selection_mode})")
        try:
            _all = read(str(conformers_xyz), index=":")
            if not isinstance(_all, list):
                _all = [_all]
        except Exception as e:
            log.error(f"[OPT] no se pudo leer {conformers_xyz}: {e}")
            return []
        # Construir índice CONF# → atoms (igual que extract_conformers_by_id)
        _frame_idx: dict[str, object] = {}
        try:
            _raw = conformers_xyz.read_text(errors="replace").splitlines()
            _nat = len(_all[0]) if _all else 0
            _fi, _li = 0, 0
            while _li < len(_raw) and _fi < len(_all):
                _l = _raw[_li].strip()
                if _l.isdigit() and int(_l) == _nat:
                    _comment = _raw[_li + 1].strip() if _li + 1 < len(_raw) else ""
                    _m = re.search(r"(CONF\d+)=", _comment)
                    _tok = _m.group(1) if _m else (_comment.split()[0] if _comment.split() else "")
                    if _tok.startswith("CONF"):
                        _frame_idx[_tok] = _all[_fi]
                    _fi += 1
                    _li += _nat + 2
                else:
                    _li += 1
        except Exception:
            pass
        conformers = [_frame_idx[cid] for cid in conf_ids if cid in _frame_idx]
        _missing = [cid for cid in conf_ids if cid not in _frame_idx]
        if _missing:
            log.warning(f"[OPT] CONF# no encontrados en ensemble: {_missing}")
        # En modo Boltzmann, conf_ids ya tiene los CONF# en orden de peso descendente;
        # se usan directamente como nombres de directorio para trazabilidad completa.
        _crest_labels: list[str] = [cid for cid in conf_ids if cid in _frame_idx]
    else:
        # Modo clásico: primeros nconf frames
        log.info(f"[OPT] optimizando {nconf} confórmeros (seleccionados automáticamente)")
        try:
            conformers = read_conformers(conformers_xyz, nconf)
        except RuntimeError as e:
            log.error(f"[OPT] no se pudo leer {conformers_xyz}: {e}")
            return []
        # Extraer los CONF# del XYZ para los primeros nconf frames
        _xyz_id_map = _extract_conf_ids_from_xyz(conformers_xyz)
        _crest_labels = [_xyz_id_map.get(i, f"CONF{i+1}") for i in range(len(conformers))]
        if any(lbl.startswith("CONF") and not lbl[4:].isdigit() for lbl in _crest_labels):
            log.warning("[OPT] algunos frames no tienen CONF# en el comentario; "
                        "usando CONF<n> como fallback")

    if not conformers:
        log.error("[OPT] el archivo de confórmeros está vacío")
        return []

    # ── Validación estado/charge: el candado que faltaba ─────────────────
    # conformers_xyz, mol_dir y charge llegan aquí como argumentos
    # independientes. Si en el nivel de orquestación se cruzó el ensemble
    # de un estado de protonación con el charge de otro (p.ej. al iterar
    # A/B/C/D/R), nada dentro de esta función lo detectaría por sí solo:
    # read_conformers() hace slicing posicional puro y el modo Boltzmann
    # solo empareja por CONF#, ninguno de los dos sabe qué charge
    # "debería" tener el ensemble que se le pasó. Si el llamador conoce el
    # conteo de átomos esperado (derivado de mol_xyz para este estado),
    # lo verificamos aquí y abortamos con un mensaje claro en vez de dejar
    # que ORCA falle en el SCF — o, peor, converja sin sentido químico.
    _actual_natoms = len(conformers[0])
    if any(len(c) != _actual_natoms for c in conformers):
        log.error(
            f"[OPT] los confórmeros leídos de {conformers_xyz} no tienen "
            f"todos el mismo número de átomos (esperado {_actual_natoms} "
            f"según el primer frame). Ensemble corrupto o mezclado."
        )
        return []
    if expected_natoms is not None and _actual_natoms != expected_natoms:
        raise RuntimeError(
            f"[OPT] DESALINEACIÓN estado/charge detectada: "
            f"conformers_xyz={conformers_xyz} tiene {_actual_natoms} átomos "
            f"por frame, pero se esperaban {expected_natoms} átomos para "
            f"mol_dir={mol_dir} con charge={charge}. "
            f"Esto es exactamente el tipo de cruce que produce geometrías "
            f"como 'opt.inp de 46 átomos con charge de un estado de 49 "
            f"átomos': el conformers_xyz de un estado de protonación se "
            f"pasó junto con el charge de otro. Verifica en el nivel de "
            f"orquestación (el ciclo que itera sobre los estados A/B/C/D/R) "
            f"que conformers_xyz, mol_dir y charge correspondan todos al "
            f"mismo estado antes de llamar a stage_opt()."
        )

    log.info(f"[OPT] {len(conformers)} confórmeros disponibles — "
             f"batch {OPT_BATCH_SIZE}x{NPROCS} cores")
    log.info(f"[OPT] IDs CREST → dirs: "
             + "  ".join(f"{lbl}→conf_{lbl[4:]}" for lbl in _crest_labels))
    template   = get_template("opt")
    successful = []

    tasks_all = [
        {
            "conf_id":  f"conf_{crest_lbl[4:]}",      # CONF39 → conf_39
            "conf_dir": mol_dir / "02_opt" / f"conf_{crest_lbl[4:]}",
            "atoms":    atoms,
            "charge":   charge,
            "template": template,
            "nprocs":   NPROCS,
            "maxcore":  MAXCORE,
            "log":      log,
        }
        for crest_lbl, atoms in zip(_crest_labels, conformers)
    ]

    n_lotes_est = -(-len(tasks_all) // OPT_BATCH_SIZE)   # ceiling division (estimado)
    _check_cpu_budget(OPT_BATCH_SIZE, NPROCS, log)

    # BUGFIX: assign_cpu_lists() puede devolver MENOS listas que el batch_size
    # pedido cuando no hay suficientes cores físicos independientes (está
    # documentado en su docstring). El código anterior hacía zip(batch,
    # cpu_lists) y lanzaba max_workers=len(batch) de todas formas — los
    # tasks que quedaban fuera del zip nunca recibían "cpu_list" y
    # run_orca() los lanzaba SIN taskset (mpirun "pelado", compitiendo con
    # los jobs sí pineados). Ahora el tamaño real del lote se recorta al
    # número de cpu_lists efectivamente disponibles, y solo se lanzan esos.
    idx      = 0
    lote_num = 0
    while idx < len(tasks_all):
        want      = min(OPT_BATCH_SIZE, len(tasks_all) - idx)
        cpu_lists = _acquire_cpu_lists(want, NPROCS, log, "OPT")
        if not cpu_lists:
            log.error("[OPT] no hay suficientes cores físicos independientes "
                       "para asignar ni un solo job — abortando etapa OPT")
            break
        batch    = tasks_all[idx: idx + len(cpu_lists)]
        lote_num += 1
        ids      = [t["conf_id"] for t in batch]

        for t, cl in zip(batch, cpu_lists):
            t["cpu_list"] = cl

        log.info(f"[OPT] ── Lote {lote_num}/~{n_lotes_est}: {ids}  "
                 f"(cores: {[','.join(map(str, cl)) for cl in cpu_lists]})")

        with ThreadPoolExecutor(max_workers=len(batch)) as pool:
            futures = {pool.submit(_opt_worker, t): t["conf_id"] for t in batch}
            try:
                for fut in as_completed(futures):
                    res = fut.result()
                    tag = "✓" if res["success"] else "✗"
                    log.info(f"[OPT]   {tag} {res['conf_id']}  ({res['elapsed']:.1f} min)")
                    if res["success"]:
                        xyz_out = resolve_opt_xyz(res["conf_dir"])
                        if xyz_out is None:
                            log.warning(f"[OPT] {res['conf_id']}: ORCA OK pero sin "
                                        f"geometría optimizada — excluido")
                        else:
                            successful.append(res["conf_dir"])
                    else:
                        log.warning(f"[OPT] {res['conf_id']} falló — excluido de etapas siguientes")
            except KeyboardInterrupt:
                log.warning("[OPT] Ctrl-C — cancelando lote, esperando workers…")
                pool.shutdown(wait=False, cancel_futures=True)
                _kill_active_procs()
                _cpu_lock_release()
                sys.exit(130)

        _cpu_lock_release()
        idx += len(cpu_lists)

    log.info(f"[OPT] {len(successful)}/{len(conformers)} confórmeros optimizados")
    return successful


def stage_numfreq(opt_dirs: list[Path], mol_dir: Path,
                  charge: int, log: logging.Logger) -> list[Path]:
    """
    Etapa 03_numfreq: frecuencias numéricas.
    Usa NPROCS_NUMFREQ cores y MAXCORE_NUMFREQ MB/core (defaults: 6 / 2300 MB).
    Configurable vía --nprocs-numfreq / --maxcore-numfreq o --cores.

    Restart automático: si en la carpeta de un conformero hay archivos
    freq_D*.res.Gradients de una corrida interrumpida, escribe el inp con
    '%freq Restart true end' para que ORCA retome sin recalcular gradientes.
    Si el .out ya tiene TERMINATED NORMALLY, salta directamente.
    """
    log.info(f"[NUMFREQ] calculando frecuencias para {len(opt_dirs)} confórmeros")
    template_fresh   = get_template("numfreq")
    template_restart = TEMPLATE_NUMFREQ_RESTART
    successful = []

    for opt_dir in opt_dirs:
        conf_id = opt_dir.name
        opt_xyz = resolve_opt_xyz(opt_dir)

        if opt_xyz is None:
            log.warning(f"[NUMFREQ] {conf_id}: no se encontró opt_opt.xyz — saltando")
            continue

        freq_dir = mol_dir / "03_numfreq" / conf_id
        freq_dir.mkdir(parents=True, exist_ok=True)

        out_file = freq_dir / "freq.out"
        inp_file = freq_dir / "freq.inp"

        if not check_success(out_file):
            # Copiar geometría si aún no está en la carpeta de freq
            if not safe_copy(opt_xyz, freq_dir / "opt_opt.xyz"):
                log.warning(f"[NUMFREQ] {conf_id}: error copiando opt_opt.xyz — saltando")
                continue

            # Detectar restart automáticamente
            if has_numfreq_restart(freq_dir):
                n_done = len(list(freq_dir.glob("freq_D*.res.Gradients")))
                log.info(f"[NUMFREQ] {conf_id}: restart detectado "
                         f"({n_done} gradientes guardados) — activando Restart true")
                inp_text = template_restart.format(
                    CHARGE=charge, MULT=MULT,
                    NPROCS_NUMFREQ=NPROCS_NUMFREQ,
                    MAXCORE_NUMFREQ=MAXCORE_NUMFREQ,
                )
            else:
                inp_text = template_fresh.format(
                    CHARGE=charge, MULT=MULT,
                    NPROCS_NUMFREQ=NPROCS_NUMFREQ,
                    MAXCORE_NUMFREQ=MAXCORE_NUMFREQ,
                )
            inp_file.write_text(inp_text)

        ok = run_orca(inp_file, out_file, log)
        if not ok:
            log.warning(f"[NUMFREQ] {conf_id} falló — excluido")
            continue

        if check_no_imaginary(out_file):
            log.info(f"[NUMFREQ] {conf_id}: sin frecuencias imaginarias ✓")
        else:
            log.warning(f"[NUMFREQ] {conf_id}: frecuencias imaginarias detectadas "
                        f"— incluido de todos modos (revisar manualmente)")
        successful.append(freq_dir)

    log.info(f"[NUMFREQ] {len(successful)}/{len(opt_dirs)} completados")
    return successful


def stage_sp(freq_dirs: list[Path], mol_dir: Path,
             charge: int, log: logging.Logger,
             batch_size: int = 0,
             nprocs_per_job: Optional[int] = None,
             maxcore_per_job: Optional[int] = None) -> list[Path]:
    """
    Etapa 04_sp: single point de alta precisión.
    Devuelve directorios sp exitosos.

    batch_size=0  → usa SP_BATCH_SIZE global (modo normal: 3 jobs × 4 cores).
    batch_size>0  → valor explícito (usado con --parallel N).
    nprocs_per_job / maxcore_per_job sobreescriben los globals NPROCS / MAXCORE
    para cada worker paralelo; si son None se usan los globals vigentes.
    """
    _batch   = batch_size if batch_size > 0 else SP_BATCH_SIZE
    _nprocs  = nprocs_per_job  if nprocs_per_job  is not None else NPROCS
    _maxcore = maxcore_per_job if maxcore_per_job is not None else MAXCORE

    log.info(f"[SP] single point para {len(freq_dirs)} confórmeros "
             f"(batch={_batch}, {_nprocs} cores/job)")
    template   = get_template("sp")
    successful = []

    n_lotes = -(-len(freq_dirs) // _batch)   # ceiling division
    log.info(f"[SP] {len(freq_dirs)} SPs en {n_lotes} lote(s) de ≤{_batch} "
             f"(CPUs activos/lote: {_batch * _nprocs})")

    tasks_all = [
        {
            "conf_id":  d.name,
            "freq_dir": d,
            "charge":   charge,
            "template": template,
            "nprocs":   _nprocs,
            "maxcore":  _maxcore,
            "mol_dir":  mol_dir,
            "log":      log,
        }
        for d in freq_dirs
    ]

    _check_cpu_budget(_batch, _nprocs, log)

    # BUGFIX: ver comentario equivalente en stage_opt — assign_cpu_lists()
    # puede devolver menos listas que _batch; recortamos el lote real al
    # número de cpu_lists disponibles para que ningún job se lance sin
    # taskset asignado.
    idx      = 0
    lote_num = 0
    while idx < len(tasks_all):
        want      = min(_batch, len(tasks_all) - idx)
        cpu_lists = _acquire_cpu_lists(want, _nprocs, log, "SP")
        if not cpu_lists:
            log.error("[SP] no hay suficientes cores físicos independientes "
                       "para asignar ni un solo job — abortando etapa SP")
            break
        batch    = tasks_all[idx: idx + len(cpu_lists)]
        lote_num += 1
        ids      = [t["conf_id"] for t in batch]

        for t, cl in zip(batch, cpu_lists):
            t["cpu_list"] = cl

        log.info(f"[SP] ── Lote {lote_num}/~{n_lotes}: {ids}  "
                 f"(cores: {[','.join(map(str, cl)) for cl in cpu_lists]})")

        with ThreadPoolExecutor(max_workers=len(batch)) as pool:
            futures = {pool.submit(_sp_worker, t): t["conf_id"] for t in batch}
            try:
                for fut in as_completed(futures):
                    res = fut.result()
                    tag = "✓" if res["success"] else "✗"
                    log.info(f"[SP]   {tag} {res['conf_id']}  ({res['elapsed']:.1f} min)")
                    if res["success"] and res["sp_dir"] is not None:
                        successful.append(res["sp_dir"])
                    else:
                        log.warning(f"[SP] {res['conf_id']} falló — excluido")
            except KeyboardInterrupt:
                log.warning("[SP] Ctrl-C — cancelando lote, esperando workers…")
                pool.shutdown(wait=False, cancel_futures=True)
                _kill_active_procs()
                _cpu_lock_release()
                sys.exit(130)

        _cpu_lock_release()
        idx += len(cpu_lists)

    log.info(f"[SP] {len(successful)}/{len(freq_dirs)} completados")
    return successful


def stage_sp_gas(sp_dirs: list[Path], mol_dir: Path, log: logging.Logger,
                  batch_size: int = 0,
                  nprocs_per_job: Optional[int] = None,
                  maxcore_per_job: Optional[int] = None) -> list[Path]:
    """
    Etapa 04b_sp_gas (OPCIONAL — requiere --spgas): single point fase-gas
    (idéntico al SP con CPCM(water) pero sin el bloque CPCM), sobre los
    mismos confórmeros y con la misma geometría/charge/mult que un
    04_sp/<conf_id>/ ya exitoso.

    Devuelve directorios 04b_sp_gas/<conf_id>/ exitosos.

    batch_size=0  → usa SP_BATCH_SIZE global (comparte recursos con SP normal).
    nprocs_per_job / maxcore_per_job sobreescriben los globals NPROCS / MAXCORE
    para cada worker paralelo; si son None se usan los globals vigentes.
    """
    _batch   = batch_size if batch_size > 0 else SP_BATCH_SIZE
    _nprocs  = nprocs_per_job  if nprocs_per_job  is not None else NPROCS
    _maxcore = maxcore_per_job if maxcore_per_job is not None else MAXCORE

    log.info(f"[SP_GAS] single point fase-gas para {len(sp_dirs)} confórmeros "
             f"(batch={_batch}, {_nprocs} cores/job)")
    template   = get_template("sp_gas")
    successful = []

    n_lotes = -(-len(sp_dirs) // _batch)   # ceiling division
    log.info(f"[SP_GAS] {len(sp_dirs)} SPs fase-gas en {n_lotes} lote(s) de ≤{_batch} "
             f"(CPUs activos/lote: {_batch * _nprocs})")

    tasks_all = [
        {
            "conf_id":  d.name,
            "sp_dir":   d,
            "template": template,
            "nprocs":   _nprocs,
            "maxcore":  _maxcore,
            "mol_dir":  mol_dir,
            "log":      log,
        }
        for d in sp_dirs
    ]

    _check_cpu_budget(_batch, _nprocs, log)

    idx      = 0
    lote_num = 0
    while idx < len(tasks_all):
        want      = min(_batch, len(tasks_all) - idx)
        cpu_lists = _acquire_cpu_lists(want, _nprocs, log, "SP_GAS")
        if not cpu_lists:
            log.error("[SP_GAS] no hay suficientes cores físicos independientes "
                       "para asignar ni un solo job — abortando etapa SP_GAS")
            break
        batch    = tasks_all[idx: idx + len(cpu_lists)]
        lote_num += 1
        ids      = [t["conf_id"] for t in batch]

        for t, cl in zip(batch, cpu_lists):
            t["cpu_list"] = cl

        log.info(f"[SP_GAS] ── Lote {lote_num}/~{n_lotes}: {ids}  "
                 f"(cores: {[','.join(map(str, cl)) for cl in cpu_lists]})")

        with ThreadPoolExecutor(max_workers=len(batch)) as pool:
            futures = {pool.submit(_sp_gas_worker, t): t["conf_id"] for t in batch}
            try:
                for fut in as_completed(futures):
                    res = fut.result()
                    tag = "✓" if res["success"] else "✗"
                    log.info(f"[SP_GAS]   {tag} {res['conf_id']}  ({res['elapsed']:.1f} min)")
                    if res["success"] and res["sp_gas_dir"] is not None:
                        successful.append(res["sp_gas_dir"])
                    else:
                        log.warning(f"[SP_GAS] {res['conf_id']} falló — excluido")
            except KeyboardInterrupt:
                log.warning("[SP_GAS] Ctrl-C — cancelando lote, esperando workers…")
                pool.shutdown(wait=False, cancel_futures=True)
                _kill_active_procs()
                _cpu_lock_release()
                sys.exit(130)

        _cpu_lock_release()
        idx += len(cpu_lists)

    log.info(f"[SP_GAS] {len(successful)}/{len(sp_dirs)} completados")
    return successful


def _find_screening_out(mol_dir: Path) -> Optional[Path]:
    """
    Busca 1_SCREENING.out de CENSO para mol_dir, probando las mismas rutas
    candidatas que usa stage_censo_prescreen() (base, CENSO/, prescreening/).
    Devuelve la primera candidata parseable, o None si ninguna lo es.
    """
    censo_dir = mol_dir / "01b_censo_prescreen"
    candidates = [
        censo_dir / "1_SCREENING.out",
        censo_dir / "CENSO" / "1_SCREENING.out",
        censo_dir / "prescreening" / "1_SCREENING.out",
    ]
    for p in candidates:
        if parse_boltzmann_weights(p):
            return p
    return None


def _filter_dirs_by_boltzmann(
    conf_dirs: list[Path],
    mol_dir:   Path,
    log:       logging.Logger,
    threshold: float,
    stage_tag: str,
) -> list[Path]:
    """
    Filtra conf_dirs (p.ej. 04_sp/conf_*) dejando solo los CONF# cuyo peso
    Boltzmann (de 1_SCREENING.out de CENSO, ver _find_screening_out) sea
    ESTRICTAMENTE MAYOR que threshold%. Es un filtro adicional por-estado:
    cada mol_dir/run (cada estado del ensemble, p.ej. cada especie de
    protonación) se filtra de forma independiente, ya que conf_dirs y
    mol_dir siempre corresponden a un único run/estado.

    Fail-open: si no hay 1_SCREENING.out parseable (p.ej. run sin CENSO
    screening) o si el filtro no deja ningún confórmero, se cae de vuelta
    a conf_dirs sin filtrar, con un warning — para no bloquear runs
    antiguos o corridos sin este dato.
    """
    screening_out = _find_screening_out(mol_dir)
    if screening_out is None:
        log.warning(
            f"[{stage_tag}] no se encontró 1_SCREENING.out parseable en "
            f"{mol_dir.name}/01b_censo_prescreen/ — no se puede filtrar por "
            f"Boltzmann (>{threshold}%), se usan los {len(conf_dirs)} "
            f"confórmero(s) disponibles sin filtrar"
        )
        return conf_dirs

    weights = parse_boltzmann_weights(screening_out)
    kept, dropped, no_wt = [], [], []
    for d in conf_dirs:
        w = weights.get(d.name)
        if w is None:
            no_wt.append(d.name)
        elif w > threshold:
            kept.append(d)
        else:
            dropped.append((d.name, w))

    if not kept:
        log.warning(
            f"[{stage_tag}] ningún confórmero de {mol_dir.name} supera "
            f">{threshold}% Boltzmann — cayendo de vuelta a los "
            f"{len(conf_dirs)} confórmero(s) sin filtrar"
        )
        return conf_dirs

    kept_str = ", ".join(f"{d.name}({weights[d.name]:.1f}%)" for d in kept)
    log.info(
        f"[{stage_tag}] filtro Boltzmann >{threshold}% ({mol_dir.name}): "
        f"{len(kept)}/{len(conf_dirs)} seleccionados: {kept_str}"
    )
    if dropped:
        log.info(
            f"[{stage_tag}] {len(dropped)} conf(s) descartados por debajo "
            f"del umbral: {[f'{cid}({w:.1f}%)' for cid, w in dropped]}"
        )
    if no_wt:
        log.warning(
            f"[{stage_tag}] {len(no_wt)} conf(s) sin peso en "
            f"{screening_out.name}, descartados: {no_wt}"
        )
    return kept


def stage_sp_smd(sp_dirs: list[Path], mol_dir: Path, log: logging.Logger,
                  batch_size: int = 0,
                  nprocs_per_job: Optional[int] = None,
                  maxcore_per_job: Optional[int] = None,
                  boltzmann_min_pct: float = 20.0) -> list[Path]:
    """
    Etapa 04c_sp_smd (OPCIONAL — requiere --spsmd): single point con
    SMD(water) (idéntico al SP con CPCM(water) pero usando SMD en vez de
    CPCM), sobre los mismos confórmeros y con la misma geometría/charge/mult
    que un 04_sp/<conf_id>/ ya exitoso.

    Antes de correr, filtra sp_dirs quedándose SOLO con los confórmeros
    cuyo peso Boltzmann (de 1_SCREENING.out de CENSO, propio de este
    mol_dir/estado) sea > boltzmann_min_pct (default 20%) — es decir, los
    mejores confórmeros de este estado, no todo el ensemble que llegó a SP.
    Fail-open: si no hay dato Boltzmann disponible o el filtro no deja
    ningún confórmero, corre sobre sp_dirs sin filtrar (ver
    _filter_dirs_by_boltzmann).

    Devuelve directorios 04c_sp_smd/<conf_id>/ exitosos.

    batch_size=0  → usa SP_BATCH_SIZE global (comparte recursos con SP normal).
    nprocs_per_job / maxcore_per_job sobreescriben los globals NPROCS / MAXCORE
    para cada worker paralelo; si son None se usan los globals vigentes.
    """
    sp_dirs = _filter_dirs_by_boltzmann(
        sp_dirs, mol_dir, log, threshold=boltzmann_min_pct, stage_tag="SP_SMD"
    )

    _batch   = batch_size if batch_size > 0 else SP_BATCH_SIZE
    _nprocs  = nprocs_per_job  if nprocs_per_job  is not None else NPROCS
    _maxcore = maxcore_per_job if maxcore_per_job is not None else MAXCORE

    log.info(f"[SP_SMD] single point SMD(water) para {len(sp_dirs)} confórmeros "
             f"(batch={_batch}, {_nprocs} cores/job)")
    template   = get_template("sp_smd")
    successful = []

    n_lotes = -(-len(sp_dirs) // _batch)   # ceiling division
    log.info(f"[SP_SMD] {len(sp_dirs)} SPs SMD en {n_lotes} lote(s) de ≤{_batch} "
             f"(CPUs activos/lote: {_batch * _nprocs})")

    tasks_all = [
        {
            "conf_id":  d.name,
            "sp_dir":   d,
            "template": template,
            "nprocs":   _nprocs,
            "maxcore":  _maxcore,
            "mol_dir":  mol_dir,
            "log":      log,
        }
        for d in sp_dirs
    ]

    _check_cpu_budget(_batch, _nprocs, log)

    idx      = 0
    lote_num = 0
    while idx < len(tasks_all):
        want      = min(_batch, len(tasks_all) - idx)
        cpu_lists = _acquire_cpu_lists(want, _nprocs, log, "SP_SMD")
        if not cpu_lists:
            log.error("[SP_SMD] no hay suficientes cores físicos independientes "
                       "para asignar ni un solo job — abortando etapa SP_SMD")
            break
        batch    = tasks_all[idx: idx + len(cpu_lists)]
        lote_num += 1
        ids      = [t["conf_id"] for t in batch]

        for t, cl in zip(batch, cpu_lists):
            t["cpu_list"] = cl

        log.info(f"[SP_SMD] ── Lote {lote_num}/~{n_lotes}: {ids}  "
                 f"(cores: {[','.join(map(str, cl)) for cl in cpu_lists]})")

        with ThreadPoolExecutor(max_workers=len(batch)) as pool:
            futures = {pool.submit(_sp_smd_worker, t): t["conf_id"] for t in batch}
            try:
                for fut in as_completed(futures):
                    res = fut.result()
                    tag = "✓" if res["success"] else "✗"
                    log.info(f"[SP_SMD]   {tag} {res['conf_id']}  ({res['elapsed']:.1f} min)")
                    if res["success"] and res["sp_smd_dir"] is not None:
                        successful.append(res["sp_smd_dir"])
                    else:
                        log.warning(f"[SP_SMD] {res['conf_id']} falló — excluido")
            except KeyboardInterrupt:
                log.warning("[SP_SMD] Ctrl-C — cancelando lote, esperando workers…")
                pool.shutdown(wait=False, cancel_futures=True)
                _kill_active_procs()
                _cpu_lock_release()
                sys.exit(130)

        _cpu_lock_release()
        idx += len(cpu_lists)

    log.info(f"[SP_SMD] {len(successful)}/{len(sp_dirs)} completados")
    return successful


def stage_spsmd_batch(group_dir: Path, log: logging.Logger,
                       boltzmann_min_pct: float = 20.0) -> None:
    """
    Modo batch de --spsmd: recorre TODAS las corridas dentro de una carpeta
    de grupo (<grupo>/runs/*/04_sp/conf_*) y calcula el SP_SMD que falte
    para cada confórmero con SP CPCM ya exitoso, sin tocar CREST/OPT/
    NUMFREQ/TDDFT de ninguna. Como stage_sp_smd usa check_success() para
    saltar los ya hechos, correr el mismo comando dos veces es seguro
    (solo procesa los faltantes).

    Cada mol_dir es un run/estado independiente, así que el filtro
    Boltzmann (> boltzmann_min_pct, ver stage_sp_smd) se aplica por
    separado dentro de cada uno.

    Se invoca cuando --spsmd se pasa con --mol apuntando a una carpeta de
    grupo (ej. MGLYVAL) en vez de a un .xyz específico.
    """
    runs_dir = group_dir / "runs"
    if not runs_dir.is_dir():
        sys.exit(f"ERROR: no se encontró {runs_dir} — "
                  f"¿'{group_dir.name}' es una carpeta de grupo válida?")

    total_ok, total_pend = 0, 0
    for mol_dir in sorted(runs_dir.iterdir()):
        if not mol_dir.is_dir():
            continue
        sp_root = mol_dir / "04_sp"
        if not sp_root.exists():
            continue
        sp_dirs = [d for d in sorted(sp_root.glob("conf_*"))
                   if check_success(d / "sp.out")]
        if not sp_dirs:
            continue

        log.info(f"[SP_SMD] {mol_dir.name}: {len(sp_dirs)} confórmeros con SP CPCM completo")
        done = stage_sp_smd(sp_dirs, mol_dir, log, boltzmann_min_pct=boltzmann_min_pct)
        total_ok += len(done)
        total_pend += len(sp_dirs) - len(done)

    log.info(f"[SP_SMD] batch finalizado: {total_ok} completados, {total_pend} pendientes/fallidos")


def stage_spgas_batch(group_dir: Path, log: logging.Logger) -> None:
    """
    Modo batch de --spgas: recorre TODAS las corridas dentro de una carpeta
    de grupo (<grupo>/runs/*/04_sp/conf_*) y calcula el SP_GAS que falte
    para cada confórmero con SP CPCM ya exitoso, sin tocar CREST/OPT/
    NUMFREQ/TDDFT de ninguna. Como stage_sp_gas usa check_success() para
    saltar los ya hechos, correr el mismo comando dos veces es seguro
    (solo procesa los faltantes).

    Se invoca cuando --spgas se pasa con --mol apuntando a una carpeta de
    grupo (ej. MGLYVAL) en vez de a un .xyz específico.
    """
    runs_dir = group_dir / "runs"
    if not runs_dir.is_dir():
        sys.exit(f"ERROR: no se encontró {runs_dir} — "
                  f"¿'{group_dir.name}' es una carpeta de grupo válida?")

    total_ok, total_pend = 0, 0
    for mol_dir in sorted(runs_dir.iterdir()):
        if not mol_dir.is_dir():
            continue
        sp_root = mol_dir / "04_sp"
        if not sp_root.exists():
            continue
        sp_dirs = [d for d in sorted(sp_root.glob("conf_*"))
                   if check_success(d / "sp.out")]
        if not sp_dirs:
            continue

        log.info(f"[SP_GAS] {mol_dir.name}: {len(sp_dirs)} confórmeros con SP CPCM completo")
        done = stage_sp_gas(sp_dirs, mol_dir, log)
        total_ok += len(done)
        total_pend += len(sp_dirs) - len(done)

    log.info(f"[SP_GAS] batch finalizado: {total_ok} completados, {total_pend} pendientes/fallidos")


def _tddft_worker(task: dict) -> dict:
    """
    Worker para un único TDDFT, diseñado para ejecutarse en un hilo de
    ThreadPoolExecutor. ORCA corre como subproceso externo (libera el GIL).

    task debe contener:
        conf_id, sp_dir, charge, template, nprocs, maxcore, mol_dir, log
    Opcional:
        cpu_list: list[int] de CPUs lógicas exclusivas para este worker.

    Devuelve:
        dict con conf_id, tddft_dir (Path o None), success, elapsed (min)
    """
    conf_id  : str            = task["conf_id"]
    sp_dir   : Path           = task["sp_dir"]
    charge   : int            = task["charge"]
    template : str            = task["template"]
    nprocs   : int            = task["nprocs"]
    maxcore  : int            = task["maxcore"]
    mol_dir  : Path           = task["mol_dir"]
    log      : logging.Logger = task["log"]
    cpu_list : Optional[list[int]] = task.get("cpu_list")

    gbw_file = sp_dir / "sp.gbw"
    opt_xyz  = resolve_opt_xyz(sp_dir)

    if not gbw_file.exists():
        log.warning(f"[TDDFT] {conf_id}: no se encontró sp.gbw — saltando")
        return {"conf_id": conf_id, "tddft_dir": None, "success": False, "elapsed": 0.0}
    if opt_xyz is None:
        log.warning(f"[TDDFT] {conf_id}: no se encontró opt_opt.xyz ni opt.xyz — saltando")
        return {"conf_id": conf_id, "tddft_dir": None, "success": False, "elapsed": 0.0}

    tddft_dir = mol_dir / "05_tddft" / conf_id
    tddft_dir.mkdir(parents=True, exist_ok=True)

    out_file = tddft_dir / "tddft.out"
    inp_file = tddft_dir / "tddft.inp"

    if not check_success(out_file):
        if not safe_copy(opt_xyz, tddft_dir / "sp_opt.xyz"):
            log.warning(f"[TDDFT] {conf_id}: error copiando sp_opt.xyz — saltando")
            return {"conf_id": conf_id, "tddft_dir": None, "success": False, "elapsed": 0.0}
        if not safe_copy(gbw_file, tddft_dir / "sp.gbw"):
            log.warning(f"[TDDFT] {conf_id}: error copiando sp.gbw — saltando")
            return {"conf_id": conf_id, "tddft_dir": None, "success": False, "elapsed": 0.0}
        inp_text = template.format(
            CHARGE=charge, MULT=MULT,
            NPROCS=nprocs, MAXCORE=maxcore,
            GBW_PATH="sp.gbw",
        )
        inp_file.write_text(inp_text)

    t0         = time.time()
    worker_env = build_omp_env(nprocs, cpu_list=cpu_list)
    success    = run_orca(inp_file, out_file, log, env=worker_env, cpu_list=cpu_list)
    elapsed    = (time.time() - t0) / 60

    return {
        "conf_id":   conf_id,
        "tddft_dir": tddft_dir if success else None,
        "success":   success,
        "elapsed":   elapsed,
    }


def stage_tddft(sp_dirs: list[Path], mol_dir: Path,
                charge: int, log: logging.Logger) -> list[Path]:
    """
    Etapa 05_tddft: TD-DFT (NRoots 15) con guess orbital del SP.

    Corre TDDFT_BATCH_SIZE confórmeros en paralelo (2 jobs × 6 cores = 12 cores).
    MAXCORE_TDDFT ya incluye el factor RAM_PENALTY_TDDFT aplicado en ResourceProfile.
    """
    log.info(f"[TDDFT] calculando espectros para {len(sp_dirs)} confórmeros — "
             f"batch {TDDFT_BATCH_SIZE}×{NPROCS_TDDFT} cores (NRoots 15)")
    template   = get_template("tddft")
    successful = []

    tasks_all = [
        {
            "conf_id":  d.name,
            "sp_dir":   d,
            "charge":   charge,
            "template": template,
            "nprocs":   NPROCS_TDDFT,
            "maxcore":  MAXCORE_TDDFT,
            "mol_dir":  mol_dir,
            "log":      log,
        }
        for d in sp_dirs
    ]

    n_lotes_est = -(-len(tasks_all) // TDDFT_BATCH_SIZE)   # ceiling division (estimado)
    _check_cpu_budget(TDDFT_BATCH_SIZE, NPROCS_TDDFT, log)

    # BUGFIX: ver comentario equivalente en stage_opt/stage_sp — assign_cpu_lists()
    # puede devolver menos listas que TDDFT_BATCH_SIZE; recortamos el lote real al
    # número de cpu_lists disponibles para que ningún job se lance sin taskset
    # asignado (el mpirun "pelado" visto en ps -eLf).
    idx      = 0
    lote_num = 0
    while idx < len(tasks_all):
        want      = min(TDDFT_BATCH_SIZE, len(tasks_all) - idx)
        cpu_lists = _acquire_cpu_lists(want, NPROCS_TDDFT, log, "TDDFT")
        if not cpu_lists:
            log.error("[TDDFT] no hay suficientes cores físicos independientes "
                       "para asignar ni un solo job — abortando etapa TDDFT")
            break
        batch    = tasks_all[idx: idx + len(cpu_lists)]
        lote_num += 1
        ids      = [t["conf_id"] for t in batch]

        for t, cl in zip(batch, cpu_lists):
            t["cpu_list"] = cl

        log.info(f"[TDDFT] ── Lote {lote_num}/~{n_lotes_est}: {ids}  "
                 f"(cores: {[','.join(map(str, cl)) for cl in cpu_lists]})")

        with ThreadPoolExecutor(max_workers=len(batch)) as pool:
            futures = {pool.submit(_tddft_worker, t): t["conf_id"] for t in batch}
            try:
                for fut in as_completed(futures):
                    res = fut.result()
                    tag = "✓" if res["success"] else "✗"
                    log.info(f"[TDDFT]   {tag} {res['conf_id']}  ({res['elapsed']:.1f} min)")
                    if res["success"] and res["tddft_dir"] is not None:
                        successful.append(res["tddft_dir"])
                    else:
                        log.warning(f"[TDDFT] {res['conf_id']} falló")
            except KeyboardInterrupt:
                log.warning("[TDDFT] Ctrl-C — cancelando lote, esperando workers…")
                pool.shutdown(wait=False, cancel_futures=True)
                _kill_active_procs()
                _cpu_lock_release()
                sys.exit(130)

        _cpu_lock_release()
        idx += len(cpu_lists)

    log.info(f"[TDDFT] {len(successful)}/{len(sp_dirs)} completados")
    return successful


# ─────────────────────────────────────────────
#  ETAPA SP_NROOT  (--nroot)
# ─────────────────────────────────────────────

def stage_sp_nroot(mol_dir: Path, charge: int, log: logging.Logger,
                   n_top: int = 5) -> list[Path]:
    """
    Etapa 06_sp_nroot: SP con NRoots 1 (TD-DFT mínimo) sobre los N mejores
    confórmeros según energía DFT de OPT (02_opt/conf_*/opt.out).

    Selección: ordena todos los conf_* exitosos de 02_opt/ por energía SCF
    (FINAL SINGLE POINT ENERGY del OPT) y toma los n_top más bajos.

    Escribe los inputs como 'sp_nroot.inp' en 06_sp_nroot/<conf_id>/.
    Copia la geometría optimizada (opt_opt.xyz) de 02_opt/<conf_id>/.

    Devuelve lista de directorios 06_sp_nroot/<conf_id>/ completados con éxito.
    """
    opt_root = mol_dir / "02_opt"
    if not opt_root.exists():
        log.error("[SP_NROOT] No existe 02_opt/ — corre primero las etapas OPT.")
        return []

    # ── Recopilar energías OPT ────────────────────────────────────
    energies: list[tuple[float, Path]] = []
    for conf_dir in sorted(opt_root.glob("conf_*")):
        out_file = conf_dir / "opt.out"
        if not check_success(out_file):
            continue
        e = extract_sp_energy(out_file)
        if e is not None:
            energies.append((e, conf_dir))
        else:
            log.warning(f"[SP_NROOT] {conf_dir.name}: OPT completado pero sin energía parseable — excluido")

    if not energies:
        log.error("[SP_NROOT] Ningún confórmero OPT exitoso con energía válida.")
        return []

    # Ordenar por energía ascendente y tomar los n_top mejores
    energies.sort(key=lambda x: x[0])
    selected = energies[:n_top]

    e_ref = selected[0][0]
    log.info(f"[SP_NROOT] {len(energies)} conf(s) OPT disponibles → "
             f"seleccionando top-{len(selected)} por E(OPT):")
    for rank, (e, conf_dir) in enumerate(selected, start=1):
        dE = (e - e_ref) * HARTREE_TO_KCAL
        log.info(f"[SP_NROOT]   #{rank:2d}  {conf_dir.name}  "
                 f"E = {e:.10f} Ha  ΔE = {dE:+.4f} kcal/mol")

    # ── Calcular SP+NRoots 1 ─────────────────────────────────────
    template   = get_template("sp_nroot")
    nroot_root = mol_dir / "06_sp_nroot"
    nroot_root.mkdir(parents=True, exist_ok=True)
    successful = []

    for e, opt_dir in selected:
        conf_id  = opt_dir.name
        opt_xyz  = resolve_opt_xyz(opt_dir)

        if opt_xyz is None:
            log.warning(f"[SP_NROOT] {conf_id}: no se encontró opt_opt.xyz — saltando")
            continue

        nroot_dir = nroot_root / conf_id
        nroot_dir.mkdir(parents=True, exist_ok=True)

        out_file = nroot_dir / "sp_nroot.out"
        inp_file = nroot_dir / "sp_nroot.inp"

        if not check_success(out_file):
            if not safe_copy(opt_xyz, nroot_dir / "opt_opt.xyz"):
                log.warning(f"[SP_NROOT] {conf_id}: error copiando opt_opt.xyz — saltando")
                continue
            inp_text = template.format(
                CHARGE=charge, MULT=MULT,
                NPROCS=NPROCS, MAXCORE=MAXCORE,
            )
            inp_file.write_text(inp_text)

        ok = run_orca(inp_file, out_file, log)
        if ok:
            successful.append(nroot_dir)
        else:
            log.warning(f"[SP_NROOT] {conf_id} falló")

    log.info(f"[SP_NROOT] {len(successful)}/{len(selected)} completados")
    return successful


# ─────────────────────────────────────────────
#  ETAPA ADC2  (séptima etapa, OPCIONAL — requiere --adc2 o --start adc2)
# ─────────────────────────────────────────────

def stage_adc2(
    mol_dir: Path,
    charge: int,
    log: logging.Logger,
    mol_name: str = "",
) -> Optional[Path]:
    """
    Etapa 07_adc2: validación ADC2 (NRoots 1) sobre un único confórmero —
    el dominante, definido como el de menor energía SP (FINAL SINGLE POINT
    ENERGY en 04_sp/conf_*/sp.out).

    La geometría se toma directamente de 04_sp/<conf_id>/opt_opt.xyz (copiada
    ahí por la propia etapa SP desde 02_opt/03_numfreq). ADC2 NO depende de
    TDDFT ni de ninguna otra etapa posterior a SP — así --start adc2 solo
    necesita que 04_sp/ exista, sin requerir 05_tddft ni el resto del
    pipeline (relevante porque ADC2 necesita ~400GB de disco, por lo que
    Kmilo suele respaldar y borrar todo lo demás antes de correrlo).

    Requiere únicamente que 04_sp/<conf_id>/ exista (out + opt_opt.xyz) para
    el confórmero dominante. Si falta la geometría, aborta esa etapa con
    error explicado.

    mol_name: nombre de la corrida (mol_dir.name si no se pasa explícitamente).
              Se usa para buscar λmax experimental en molecules.yml /
              EXP_LAMBDA_FALLBACK y calcular Δλ = λ_ADC2 − λ_exp.

    Devuelve el directorio 07_adc2/<conf_id>/ si el cálculo terminó bien,
    o None si falló o no se pudo seleccionar confórmero.

    Efectos secundarios:
      - Escribe 07_adc2/adc2_results.json con los valores numéricos de la
        transición S0→S1 (λ_nm, E_eV, fosc, λ_exp_nm, delta_lambda_nm).
        El archivo se sobreescribe en cada llamada exitosa con parseo correcto.
    """
    # Resolver mol_name desde mol_dir si no se pasó explícitamente
    if not mol_name:
        mol_name = mol_dir.name

    sp_root = mol_dir / "04_sp"
    if not sp_root.exists():
        log.error("[ADC2] No existe 04_sp/ — corre primero la etapa SP.")
        return None

    # ── Seleccionar confórmero dominante: menor energía SP ───────
    energies: list[tuple[float, Path]] = []
    for conf_dir in sorted(sp_root.glob("conf_*")):
        out_file = conf_dir / "sp.out"
        if not check_success(out_file):
            continue
        e = extract_sp_energy(out_file)
        if e is not None:
            energies.append((e, conf_dir))
        else:
            log.warning(f"[ADC2] {conf_dir.name}: SP completado pero sin energía "
                        f"parseable — excluido")

    if not energies:
        log.error("[ADC2] Ningún confórmero SP exitoso con energía válida. Abortando etapa.")
        return None

    energies.sort(key=lambda x: x[0])
    e_min, sp_dir = energies[0]
    conf_id = sp_dir.name

    log.info(f"[ADC2] {len(energies)} conf(s) SP disponibles — "
             f"confórmero dominante: {conf_id} (E_SP = {e_min:.10f} Ha)")
    if len(energies) > 1:
        runner_up = energies[1]
        dE = (runner_up[0] - e_min) * HARTREE_TO_KCAL
        log.info(f"[ADC2] siguiente más estable: {runner_up[1].name} "
                 f"(ΔE = +{dE:.4f} kcal/mol respecto al dominante)")

    # ── Localizar geometría: directamente desde 04_sp/<conf_id>/ ─
    # (opt_opt.xyz ya fue copiada ahí por la etapa SP; no depende de TDDFT)
    geom_src = resolve_opt_xyz(sp_dir)
    if geom_src is None:
        log.error(
            f"[ADC2] No se encontró opt_opt.xyz ni opt.xyz en {sp_dir} — "
            f"el confórmero dominante ({conf_id}) tiene sp.out pero le falta "
            f"la geometría en 04_sp/. Vuelve a correr la etapa SP para este "
            f"confórmero."
        )
        return None

    # ── Preparar y correr ADC2 ────────────────────────────────────
    adc2_root = mol_dir / "07_adc2"
    adc2_dir  = adc2_root / conf_id
    adc2_dir.mkdir(parents=True, exist_ok=True)

    out_file = adc2_dir / "adc2.out"
    inp_file = adc2_dir / "adc2.inp"

    if not check_success(out_file):
        if not safe_copy(geom_src, adc2_dir / "adc2_geom.xyz"):
            log.error(f"[ADC2] {conf_id}: error copiando geometría — abortando etapa.")
            return None
        template = get_template("adc2")
        inp_text = template.format(
            CHARGE=charge, MULT=MULT,
            NPROCS_ADC2=NPROCS_ADC2, MAXCORE_ADC2=MAXCORE_ADC2,
        )
        inp_file.write_text(inp_text)

    ok = run_orca(inp_file, out_file, log)
    if ok:
        log.info(f"[ADC2] {conf_id} completado ✓ → {adc2_dir}")

        # ── Extraer λmax S0→S1 de la única raíz calculada ────────
        result = extract_adc2_lambda(out_file)
        if result is not None:
            lambda_nm, energy_ev, fosc = result

            # Buscar λmax experimental y calcular Δλ
            lambda_exp = load_exp_lambda(mol_name, log)
            if lambda_exp is not None:
                delta_lambda = lambda_nm - lambda_exp
                sign = "+" if delta_lambda >= 0 else ""
                log.info(
                    f"[ADC2] λmax(S0→S1) = {lambda_nm:.1f} nm  "
                    f"({energy_ev:.4f} eV)  fosc = {fosc:.6f}  │  "
                    f"λexp = {lambda_exp:.1f} nm  │  "
                    f"Δλ = {sign}{delta_lambda:.1f} nm"
                )
            else:
                delta_lambda = None
                log.info(
                    f"[ADC2] λmax(S0→S1) = {lambda_nm:.1f} nm  "
                    f"({energy_ev:.4f} eV)  fosc = {fosc:.6f}"
                )

            # ── Guardar resultados en JSON ────────────────────────
            results_json = adc2_root / "adc2_results.json"
            record: dict = {
                "mol_name":        mol_name,
                "conf_id":         conf_id,
                "lambda_nm":       round(lambda_nm, 2),
                "energy_eV":       round(energy_ev, 6),
                "fosc":            round(fosc, 8),
                "lambda_exp_nm":   lambda_exp,
                "delta_lambda_nm": (
                    round(delta_lambda, 2) if delta_lambda is not None else None
                ),
                "timestamp": datetime.datetime.now().isoformat(timespec="seconds"),
            }
            try:
                results_json.write_text(
                    json.dumps(record, indent=2, ensure_ascii=False) + "\n"
                )
                log.info(f"[ADC2] resultados guardados → {results_json}")
            except OSError as e:
                log.warning(f"[ADC2] no se pudo escribir {results_json}: {e}")

        else:
            log.warning(
                f"[ADC2] {conf_id}: cálculo OK pero no se pudo extraer λmax "
                f"del espectro RIGHT — revisar {out_file.name}"
            )

        return adc2_dir

    log.error(f"[ADC2] {conf_id} falló.")
    return None


def _opt_worker(task: dict) -> dict:
    """
    Worker para un único OPT, diseñado para ejecutarse en un hilo de
    ThreadPoolExecutor.  ORCA corre como subproceso externo (libera el GIL),
    por lo que el paralelismo real es efectivo.

    task debe contener:
        conf_id, conf_dir, atoms, charge, template, nprocs, maxcore, log
    Opcional:
        cpu_list: list[int] de CPUs lógicas exclusivas para este worker,
                  asignado por el orquestador del lote vía assign_cpu_lists()
                  para evitar que OpenMPI distribuya mal los procesos entre
                  workers concurrentes.

    Devuelve:
        dict con conf_id, conf_dir, energy (Ha o None), success, elapsed (min)
    """
    conf_id  : str           = task["conf_id"]
    conf_dir : Path          = task["conf_dir"]
    atoms                    = task["atoms"]
    charge   : int           = task["charge"]
    template : str           = task["template"]
    nprocs   : int           = task["nprocs"]
    maxcore  : int           = task["maxcore"]
    log      : logging.Logger = task["log"]
    cpu_list : Optional[list[int]] = task.get("cpu_list")

    inp_file = conf_dir / "opt.inp"
    out_file = conf_dir / "opt.out"

    conf_dir.mkdir(parents=True, exist_ok=True)

    if not check_success(out_file):
        coords   = atoms_to_xyz_block(atoms)
        inp_text = template.format(
            CHARGE=charge, MULT=MULT,
            NPROCS=nprocs, MAXCORE=maxcore,
            COORDS=coords,
        )
        inp_file.write_text(inp_text)

    t0      = time.time()
    # Env por worker: construido con el nprocs del task, no con el global.
    # Evita race conditions al modificar os.environ desde hilos concurrentes.
    worker_env = build_omp_env(nprocs, cpu_list=cpu_list)
    success = run_orca(inp_file, out_file, log, env=worker_env, cpu_list=cpu_list)
    elapsed = (time.time() - t0) / 60

    energy: Optional[float] = None
    if success:
        if resolve_opt_xyz(conf_dir) is None:
            log.warning(f"[OPT] {conf_id}: ORCA OK pero sin geometría optimizada")
            success = False
        else:
            energy = extract_sp_energy(out_file)

    return {
        "conf_id":  conf_id,
        "conf_dir": conf_dir,
        "energy":   energy,
        "success":  success,
        "elapsed":  elapsed,
    }


def _sp_worker(task: dict) -> dict:
    """
    Worker para un único SP, diseñado para ejecutarse en un hilo de
    ThreadPoolExecutor.  ORCA corre como subproceso externo (libera el GIL),
    por lo que el paralelismo real es efectivo.

    task debe contener:
        conf_id, freq_dir, charge, template, nprocs, maxcore, mol_dir, log
    Opcional:
        cpu_list: list[int] de CPUs lógicas exclusivas para este worker.

    Devuelve:
        dict con conf_id, sp_dir (Path o None), success, elapsed (min)
    """
    conf_id  : str            = task["conf_id"]
    freq_dir : Path           = task["freq_dir"]
    charge   : int            = task["charge"]
    template : str            = task["template"]
    nprocs   : int            = task["nprocs"]
    maxcore  : int            = task["maxcore"]
    mol_dir  : Path           = task["mol_dir"]
    log      : logging.Logger = task["log"]
    cpu_list : Optional[list[int]] = task.get("cpu_list")

    opt_xyz = resolve_opt_xyz(freq_dir)
    if opt_xyz is None:
        log.warning(f"[SP] {conf_id}: no se encontró opt_opt.xyz ni opt.xyz — saltando")
        return {"conf_id": conf_id, "sp_dir": None, "success": False, "elapsed": 0.0}

    sp_dir = mol_dir / "04_sp" / conf_id
    sp_dir.mkdir(parents=True, exist_ok=True)

    out_file = sp_dir / "sp.out"
    inp_file = sp_dir / "sp.inp"

    if not check_success(out_file):
        if not safe_copy(opt_xyz, sp_dir / "opt_opt.xyz"):
            log.warning(f"[SP] {conf_id}: error copiando opt_opt.xyz — saltando")
            return {"conf_id": conf_id, "sp_dir": None, "success": False, "elapsed": 0.0}
        inp_text = template.format(
            CHARGE=charge, MULT=MULT,
            NPROCS=nprocs, MAXCORE=maxcore,
        )
        inp_file.write_text(inp_text)

    t0         = time.time()
    worker_env = build_omp_env(nprocs, cpu_list=cpu_list)
    success    = run_orca(inp_file, out_file, log, env=worker_env, cpu_list=cpu_list)
    elapsed    = (time.time() - t0) / 60

    return {
        "conf_id": conf_id,
        "sp_dir":  sp_dir if success else None,
        "success": success,
        "elapsed": elapsed,
    }


def _sp_gas_worker(task: dict) -> dict:
    """
    Worker para un único SP fase-gas (sin CPCM), diseñado para ejecutarse
    en un hilo de ThreadPoolExecutor. Mismo patrón que _sp_worker, pero:
      - parte de un 04_sp/<conf_id>/ ya exitoso (sp_dir), no de un freq_dir.
      - charge/mult se leen de sp_dir/"sp.inp" (no de task["charge"]), para
        garantizar que coincide exactamente con el SP CPCM que le dio origen.
      - la geometría se copia desde sp_dir/"opt_opt.xyz" (la misma que usó
        el SP CPCM), no desde 02_opt/ ni 03_numfreq/ directamente — así el
        SP_GAS queda 100% trazable al mismo SP.

    task debe contener:
        conf_id, sp_dir, template, nprocs, maxcore, mol_dir, log
    Opcional:
        cpu_list: list[int] de CPUs lógicas exclusivas para este worker.

    Devuelve:
        dict con conf_id, sp_gas_dir (Path o None), success, elapsed (min)
    """
    conf_id  : str            = task["conf_id"]
    sp_dir   : Path           = task["sp_dir"]
    template : str            = task["template"]
    nprocs   : int            = task["nprocs"]
    maxcore  : int            = task["maxcore"]
    mol_dir  : Path           = task["mol_dir"]
    log      : logging.Logger = task["log"]
    cpu_list : Optional[list[int]] = task.get("cpu_list")

    sp_inp = sp_dir / "sp.inp"
    charge_mult = _parse_charge_mult_from_sp_inp(sp_inp)
    if charge_mult is None:
        log.warning(f"[SP_GAS] {conf_id}: no se pudo leer charge/mult de "
                    f"{sp_inp} — saltando")
        return {"conf_id": conf_id, "sp_gas_dir": None, "success": False, "elapsed": 0.0}
    charge, mult = charge_mult

    opt_xyz = sp_dir / "opt_opt.xyz"
    if not opt_xyz.exists():
        log.warning(f"[SP_GAS] {conf_id}: no se encontró {opt_xyz} — saltando")
        return {"conf_id": conf_id, "sp_gas_dir": None, "success": False, "elapsed": 0.0}

    sp_gas_dir = mol_dir / "04b_sp_gas" / conf_id
    sp_gas_dir.mkdir(parents=True, exist_ok=True)

    out_file = sp_gas_dir / "sp_gas.out"
    inp_file = sp_gas_dir / "sp_gas.inp"

    if not check_success(out_file):
        if not safe_copy(opt_xyz, sp_gas_dir / "opt_opt.xyz"):
            log.warning(f"[SP_GAS] {conf_id}: error copiando opt_opt.xyz — saltando")
            return {"conf_id": conf_id, "sp_gas_dir": None, "success": False, "elapsed": 0.0}
        inp_text = template.format(
            CHARGE=charge, MULT=mult,
            NPROCS=nprocs, MAXCORE=maxcore,
        )
        inp_file.write_text(inp_text)

    t0         = time.time()
    worker_env = build_omp_env(nprocs, cpu_list=cpu_list)
    success    = run_orca(inp_file, out_file, log, env=worker_env, cpu_list=cpu_list)
    elapsed    = (time.time() - t0) / 60

    return {
        "conf_id":    conf_id,
        "sp_gas_dir": sp_gas_dir if success else None,
        "success":    success,
        "elapsed":    elapsed,
    }


def _sp_smd_worker(task: dict) -> dict:
    """
    Worker para un único SP con SMD(water), diseñado para ejecutarse en un
    hilo de ThreadPoolExecutor. Idéntico patrón que _sp_gas_worker, pero:
      - en vez de fase-gas (sin bloque de solvente), usa SMD(water) en lugar
        del CPCM(water) del SP original.
      - parte de un 04_sp/<conf_id>/ ya exitoso (sp_dir), no de un freq_dir.
      - charge/mult se leen de sp_dir/"sp.inp" (no de task["charge"]), para
        garantizar que coincide exactamente con el SP CPCM que le dio origen.
      - la geometría se copia desde sp_dir/"opt_opt.xyz" (la misma que usó
        el SP CPCM), no desde 02_opt/ ni 03_numfreq/ directamente — así el
        SP_SMD queda 100% trazable al mismo SP.

    task debe contener:
        conf_id, sp_dir, template, nprocs, maxcore, mol_dir, log
    Opcional:
        cpu_list: list[int] de CPUs lógicas exclusivas para este worker.

    Devuelve:
        dict con conf_id, sp_smd_dir (Path o None), success, elapsed (min)
    """
    conf_id  : str            = task["conf_id"]
    sp_dir   : Path           = task["sp_dir"]
    template : str            = task["template"]
    nprocs   : int            = task["nprocs"]
    maxcore  : int            = task["maxcore"]
    mol_dir  : Path           = task["mol_dir"]
    log      : logging.Logger = task["log"]
    cpu_list : Optional[list[int]] = task.get("cpu_list")

    sp_inp = sp_dir / "sp.inp"
    charge_mult = _parse_charge_mult_from_sp_inp(sp_inp)
    if charge_mult is None:
        log.warning(f"[SP_SMD] {conf_id}: no se pudo leer charge/mult de "
                    f"{sp_inp} — saltando")
        return {"conf_id": conf_id, "sp_smd_dir": None, "success": False, "elapsed": 0.0}
    charge, mult = charge_mult

    opt_xyz = sp_dir / "opt_opt.xyz"
    if not opt_xyz.exists():
        log.warning(f"[SP_SMD] {conf_id}: no se encontró {opt_xyz} — saltando")
        return {"conf_id": conf_id, "sp_smd_dir": None, "success": False, "elapsed": 0.0}

    sp_smd_dir = mol_dir / "04c_sp_smd" / conf_id
    sp_smd_dir.mkdir(parents=True, exist_ok=True)

    out_file = sp_smd_dir / "sp_smd.out"
    inp_file = sp_smd_dir / "sp_smd.inp"

    if not check_success(out_file):
        if not safe_copy(opt_xyz, sp_smd_dir / "opt_opt.xyz"):
            log.warning(f"[SP_SMD] {conf_id}: error copiando opt_opt.xyz — saltando")
            return {"conf_id": conf_id, "sp_smd_dir": None, "success": False, "elapsed": 0.0}
        inp_text = template.format(
            CHARGE=charge, MULT=mult,
            NPROCS=nprocs, MAXCORE=maxcore,
        )
        inp_file.write_text(inp_text)

    t0         = time.time()
    worker_env = build_omp_env(nprocs, cpu_list=cpu_list)
    success    = run_orca(inp_file, out_file, log, env=worker_env, cpu_list=cpu_list)
    elapsed    = (time.time() - t0) / 60

    return {
        "conf_id":    conf_id,
        "sp_smd_dir": sp_smd_dir if success else None,
        "success":    success,
        "elapsed":    elapsed,
    }


# ── Coordinación de CPUs ENTRE PROCESOS (2026-07-12) ──────────────────────
#
# assign_cpu_lists()/assign_cpu_pool() reparten cores correctamente DENTRO
# de un mismo proceso (todos los workers de un ThreadPoolExecutor comparten
# el mismo cálculo), pero cada invocación de `python3 pipeline.py` es un
# proceso independiente que arranca su cálculo SIEMPRE desde el mismo pool
# estático (E_CORE_SINGLES = [4,5,6,7,8,9,10,11], empezando por el frente).
#
# Confirmado con `ps` + `/proc/<pid>/status`: corriendo SP_GAS (--cores 8)
# en paralelo con otro SP (--cores 4) ya en curso, ambos procesos reclamaron
# PUs 4-9 (oversuscritos 3-5x cada uno) mientras 10-11 nunca se usaron —
# "7 cores en total" en vez de los ~12 esperados entre ambos. No hay nada
# que le avise a un proceso que el otro ya reservó ciertos PUs.
#
# Fix: lock file compartido en /tmp (sobrevive solo mientras el sistema no
# reinicia, lo cual es deseable — descarta reservas viejas automáticamente
# en cada boot) con reservas por PID, protegido con flock para acceso
# atómico. Cada proceso, antes de calcular su asignación, lee qué PUs
# tienen reservados OTROS PIDs vivos (podando los de PIDs muertos) y los
# excluye del pool disponible; después de calcular, registra su propia
# reserva. Fail-open: si el lock file no se puede leer/escribir por
# cualquier motivo, el pipeline sigue funcionando igual que antes (con
# riesgo de colisión, pero sin abortar).
CPU_LOCK_FILE = Path("/tmp/maaqc_cpu_lock.json")


def _cpu_lock_busy_cpus() -> set[int]:
    """
    Devuelve el conjunto de CPUs lógicas reservadas AHORA MISMO por OTROS
    procesos vivos (pipeline.py u otros que usen este lock), leyendo y
    depurando CPU_LOCK_FILE. Nunca incluye la reserva del propio proceso.
    """
    if not CPU_LOCK_FILE.exists():
        return set()
    my_pid = str(os.getpid())
    try:
        with open(CPU_LOCK_FILE, "r+") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            try:
                raw = f.read().strip()
                data = json.loads(raw) if raw else {}
                alive: dict[str, list[int]] = {}
                for pid_str, cpus in data.items():
                    if pid_str == my_pid:
                        continue
                    try:
                        os.kill(int(pid_str), 0)
                    except (OSError, ValueError):
                        continue   # PID muerto o inválido — se poda
                    alive[pid_str] = cpus
                # Reescribir depurado (poda de PIDs muertos), preservando
                # la propia reserva si existía.
                if my_pid in data:
                    alive[my_pid] = data[my_pid]
                f.seek(0)
                json.dump(alive, f)
                f.truncate()
            finally:
                fcntl.flock(f, fcntl.LOCK_UN)
        alive.pop(my_pid, None)
        return {c for cpus in alive.values() for c in cpus}
    except (OSError, ValueError, json.JSONDecodeError):
        return set()   # fail-open


def _cpu_lock_reserve(cpus: list[int]) -> None:
    """Registra los CPUs que el proceso actual (PID) usará en este lote."""
    try:
        CPU_LOCK_FILE.touch(exist_ok=True)
        with open(CPU_LOCK_FILE, "r+") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            try:
                raw = f.read().strip()
                data = json.loads(raw) if raw else {}
                data[str(os.getpid())] = cpus
                f.seek(0)
                json.dump(data, f)
                f.truncate()
            finally:
                fcntl.flock(f, fcntl.LOCK_UN)
    except (OSError, ValueError, json.JSONDecodeError):
        pass   # fail-open: sin coordinación, pero el pipeline sigue


def _cpu_lock_release() -> None:
    """Libera la reserva de CPUs del proceso actual (llamar al terminar cada lote)."""
    if not CPU_LOCK_FILE.exists():
        return
    try:
        with open(CPU_LOCK_FILE, "r+") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            try:
                raw = f.read().strip()
                data = json.loads(raw) if raw else {}
                data.pop(str(os.getpid()), None)
                f.seek(0)
                json.dump(data, f)
                f.truncate()
            finally:
                fcntl.flock(f, fcntl.LOCK_UN)
    except (OSError, ValueError, json.JSONDecodeError):
        pass


def _acquire_cpu_lists(want: int, nprocs_per_job: int,
                        log: logging.Logger, tag: str) -> list[list[int]]:
    """
    Wrapper de assign_cpu_lists() que coordina con OTRAS instancias de
    pipeline.py corriendo en paralelo vía CPU_LOCK_FILE: excluye los PUs
    que otro proceso vivo ya tiene reservados, calcula la asignación sobre
    el resto, y registra la reserva propia bajo el PID actual.

    El caller es responsable de llamar _cpu_lock_release() cuando el lote
    termina (éxito, fallo o Ctrl-C) para liberar los PUs para otros procesos.
    """
    busy_others = _cpu_lock_busy_cpus()
    cpu_lists = assign_cpu_lists(want, nprocs_per_job, exclude=busy_others)
    if busy_others:
        log.info(f"[{tag}] excluyendo PUs reservados por otra(s) instancia(s) "
                 f"de pipeline.py: {sorted(busy_others)}")
    if cpu_lists:
        _cpu_lock_reserve([c for cl in cpu_lists for c in cl])
    return cpu_lists


def _acquire_cpu_pool(n_cores: int, log: logging.Logger, tag: str) -> list[int]:
    """
    Wrapper de assign_cpu_pool() — análogo a _acquire_cpu_lists() pero para
    el pool "plano" que usa CENSO (un único proceso de larga duración, no
    un lote de workers). Excluye los PUs que otra instancia de pipeline.py
    ya tiene reservados y registra la reserva propia bajo el PID actual.

    El caller es responsable de llamar _cpu_lock_release() cuando CENSO
    termina (éxito, timeout, error o Ctrl-C) para liberar los PUs.
    """
    busy_others = _cpu_lock_busy_cpus()
    cpu_pool = assign_cpu_pool(n_cores, exclude=busy_others)
    if busy_others:
        log.info(f"[{tag}] excluyendo PUs reservados por otra(s) instancia(s) "
                 f"de pipeline.py: {sorted(busy_others)}")
    if cpu_pool:
        _cpu_lock_reserve(cpu_pool)
    return cpu_pool


def assign_cpu_lists(batch_size: int, nprocs_per_job: int,
                      exclude: Optional[set[int]] = None) -> list[list[int]]:
    """
    Calcula listas explícitas de CPUs lógicas (no necesariamente
    contiguas) para cada worker de un lote, basadas en la topología REAL
    de este ThinkPad T14 (Core Ultra 7), confirmada con `lscpu -e`:

        CPU lógica  →  CORE físico   →  tipo
        0, 1        →  CORE 0        →  P-core (2 hilos SMT, mismo L2)
        2, 3        →  CORE 1        →  P-core (2 hilos SMT, mismo L2)
        4           →  CORE 2        →  E-core
        5           →  CORE 3        →  E-core
        6           →  CORE 4        →  E-core
        7           →  CORE 5        →  E-core
        8           →  CORE 6        →  E-core
        9           →  CORE 7        →  E-core
        10          →  CORE 8        →  E-core
        11          →  CORE 9        →  E-core
        12          →  CORE 10       →  LP E-core (max 400 MHz — evitar)
        13          →  CORE 11       →  LP E-core (max ~1.3 GHz — evitar)

    Por qué NO se puede usar aritmética ingenua (i*n, (i+1)*n-1) ni rangos
    lo-hi contiguos: los pares (0,1) y (2,3) son hilos SMT del MISMO core
    físico, no dos cores distintos. Confinar 2 ranks MPI ahí no es
    problemático en sí mismo, pero asumir que "rango contiguo == cores
    físicos independientes" en una CPU híbrida hizo que hwloc (usado
    internamente por el mpirun que ORCA lanza para LEANSCF) leyera una
    topología de cluster L1d/L2 inconsistente con el cpuset heredado de
    taskset, produciendo "intersection without inclusion" y abortando el
    cálculo justo al iniciar la formación de Fock (visto en conf_008).

    Estrategia: cada worker recibe un conjunto explícito de CPUs lógicas
    que corresponden a unidades físicas reales — 1 P-core completo (sus 2
    hilos SMT juntos, nunca repartidos entre dos workers distintos) o
    E-cores individuales. Los LP E-cores se usan solo como último recurso
    porque corren a relojes mucho más bajos y desbalancean el lote
    (todos los workers de un batch esperan al más lento).

    nprocs_per_job determina cuántas CPUs lógicas recibe cada worker:
      - == 2: preferentemente 2 E-cores individuales (físicamente
        independientes); los P-core blocks (2 hilos SMT de UN solo core
        físico) se usan solo como último recurso — ver BUGFIX abajo.
      - == 1: 1 E-core (preferido) o, si no alcanzan, 1 hilo SMT de P-core.
      - otro valor (ej. 6 para NumFreq): reparto secuencial sobre el pool
        completo en el mismo orden de prioridad.

    Devuelve una lista de listas de CPUs, una por worker, en el mismo
    orden que `batch`. Si no hay suficientes CPUs físicas para satisfacer
    batch_size workers, la lista devuelta es más corta que batch_size
    (el caller debe manejarlo o confiar en _check_cpu_budget para advertir
    antes de llegar aquí).

    `exclude`, si se pasa, es un conjunto de PUs a NO usar — normalmente
    los PUs que otra instancia de pipeline.py ya tiene reservados AHORA
    MISMO (ver _acquire_cpu_lists() / CPU_LOCK_FILE). Un P_CORE_BLOCK se
    descarta por completo si cualquiera de sus 2 PUs está excluido (no
    tiene sentido partir un bloque SMT entre dos procesos distintos).
    """
    excl        = exclude or set()
    p_blocks    = [b for b in P_CORE_BLOCKS if not (set(b) & excl)]
    e_singles   = [c for c in E_CORE_SINGLES if c not in excl]
    lp_singles  = [c for c in LP_CORE_SINGLES if c not in excl]

    lists: list[list[int]] = []

    if nprocs_per_job == 2:
        # BUGFIX (2026-07-12, descubierto corriendo SP_GAS con --cores 8,
        # pero afecta a CUALQUIER etapa con nprocs_per_job==2 — SP/OPT
        # incluidos, simplemente no se había dado la combinación de
        # --cores que fuerza nprocs=2 antes de ahora):
        #
        # La idea original era dar cada job de 2 ranks (mpirun -np 2) un
        # P-core completo (sus 2 hilos SMT), asumiendo que
        # --map-by pe-list=<core>,<core>:ORDERED:hwtcpus --bind-to hwthread
        # (inyectado por el wrapper de mpirun vía MAAQC_PE_LIST) reparte
        # cada rank a un hilo SMT distinto del MISMO índice de Core.
        #
        # En la práctica, PRTE round-robin-mapea por CORE por defecto: con
        # pe-list="0,0" (2 entradas repitiendo el mismo índice de Core),
        # PRTE ve un único objeto Core disponible para 2 ranks y aborta con
        # "prte-rmaps-rr:not-enough-cpus" en Startup — el wrapper actual no
        # fuerza --map-by a granularidad hwthread, solo --bind-to hwthread
        # (que solo controla el binding DESPUÉS de mapear, no la unidad de
        # mapeo). Confirmado con conf_001/conf_002 (cores 0,1 y 2,3 — ambos
        # P_CORE_BLOCKS) fallando, mientras conf_003 (cores 4,5 — 2 E-cores
        # físicamente independientes) corrió sin problema.
        #
        # Fix: para nprocs_per_job==2 priorizar 2 E-cores físicamente
        # independientes (mismo patrón ya usado para nprocs_per_job>=3, ver
        # comentario más abajo) y dejar los P_CORE_BLOCKS (SMT de 1 core)
        # como último recurso, con advertencia — hasta que el wrapper de
        # mpirun soporte --map-by hwthread explícito para este caso.
        e_pool = list(e_singles)
        while len(lists) < batch_size and len(e_pool) >= 2:
            lists.append([e_pool.pop(0), e_pool.pop(0)])

        for block in p_blocks:
            if len(lists) >= batch_size:
                break
            lists.append(list(block))

        lp_pool = list(lp_singles)
        while len(lists) < batch_size and len(lp_pool) >= 2:
            lists.append([lp_pool.pop(0), lp_pool.pop(0)])

    elif nprocs_per_job == 1:
        pool = list(e_singles) + list(lp_singles) + \
               [c for block in p_blocks for c in block]
        for cpu in pool:
            if len(lists) >= batch_size:
                break
            lists.append([cpu])

    else:
        # BUGFIX: para nprocs_per_job >= 3 NO se puede meter los P_CORE_BLOCKS
        # combinados en el mismo pool secuencial que los E-cores. Cada bloque
        # P-core son 2 CPUs lógicas (SMT) respaldadas por 1 solo core físico;
        # los dos bloques juntos ([0,1,2,3]) son 4 CPUs lógicas pero solo 2
        # cores físicos. Trocear ese pool en chunks de nprocs_per_job=4 le
        # asignaba a un job 4 "cores" que en realidad son 2 físicos con HT,
        # produciendo oversubscription 2x (el patrón de utilización 66/33
        # observado: ese job compite consigo mismo por 2 cores reales
        # mientras los otros jobs, en E-cores reales, sí rinden ~100%).
        # Para trabajos anchos (nprocs_per_job >= 3) solo se reparten cores
        # físicamente independientes: E-cores primero, LP E-cores como
        # último recurso. Los P-core blocks quedan excluidos aquí — no hay
        # forma de darle 3+ hilos a un job y que sigan siendo físicamente
        # reales dentro de un bloque de 2 hilos SMT.
        full_pool = list(e_singles) + list(lp_singles)
        idx = 0
        for _ in range(batch_size):
            chunk = full_pool[idx: idx + nprocs_per_job]
            if len(chunk) < nprocs_per_job:
                break
            lists.append(chunk)
            idx += nprocs_per_job

    return lists


def assign_cpu_pool(n_cores: int, exclude: Optional[set[int]] = None) -> list[int]:
    """
    Calcula una única lista plana de CPUs lógicas — no una lista por
    worker como assign_cpu_lists() — para confinar vía `taskset -c` un
    proceso PADRE que reparte su propio presupuesto de cores entre
    sub-procesos internamente (caso CENSO).

    Por qué esto es distinto de assign_cpu_lists():
    assign_cpu_lists() reparte cores ENTRE jobs hermanos que el pipeline
    lanza por separado (ej. 3 SP en paralelo), y por eso a partir de
    nprocs_per_job >= 3 excluye los P_CORE_BLOCKS: un solo job de ese
    tamaño no puede recibir 2 hilos SMT del mismo core físico y seguir
    siendo "N cores físicos independientes" para su propio MPI/OpenMP.

    CENSO es distinto: el proceso `censo` que lanzamos aquí NO es un
    único cálculo monolítico — es un scheduler (Dask) que decide él
    mismo cuántos conformeros corren en paralelo y con cuántos hilos
    OpenMP cada uno (ver censo/parallel.py: get_cluster()/set_omp()),
    y llama a xtb con `--parallel N` por cada sub-job. Lo único que
    hace falta es confinar el ÁRBOL COMPLETO de procesos (censo +
    todos los xtb que lanza) a un conjunto fijo de N cores físicamente
    reales — sin esto, el kernel es libre de repartir los hilos OpenMP
    de sub-jobs *distintos* sobre el mismo par SMT de un P-core (o
    sobre los LP E-cores lentos), que es exactamente el patrón medido
    de "varios procesos compitiendo por los mismos cores, 6/12 cores
    aprovechados" — mismo problema estructural que el histórico bug de
    ORCA (ver comentario en assign_cpu_lists), pero a nivel de proceso
    completo en vez de por-job.

    Por eso aquí SÍ se pueden incluir los P_CORE_BLOCKS completos como
    unidades atómicas (2 hilos SMT juntos, nunca partidos): no importa
    que sean "solo" 1 core físico con 2 hilos, porque es el propio
    scheduler interno de CENSO (Dask + xtb --parallel) el que decide
    cómo subdividir el presupuesto entre sub-jobs — nuestro taskset solo
    fija el límite exterior, igual que --maxcores fija el límite lógico.

    Orden de prioridad: P-core blocks completos primero (mayor
    frecuencia), luego E-cores, y LP E-cores solo si no alcanza con lo
    anterior (se avisa por log en el caller si esto ocurre).

    `exclude`, si se pasa, son PUs reservados por OTRA instancia de
    pipeline.py AHORA MISMO (ver _acquire_cpu_lists()/CPU_LOCK_FILE) —
    se descartan del pool antes de repartir. Igual que en
    assign_cpu_lists(), un P_CORE_BLOCK se descarta completo si
    cualquiera de sus 2 PUs está excluido.

    Devuelve una lista de hasta n_cores CPUs lógicas. Si n_cores excede
    las CPUs físicas reales de la máquina (12 en este ThinkPad: 2
    P-cores + 8 E-cores, sin contar LP), la lista incluye también los
    LP E-cores como último recurso y, si aún así no alcanza, se trunca
    a las CPUs lógicas disponibles en el sistema.
    """
    excl = exclude or set()
    p_blocks   = [b for b in P_CORE_BLOCKS if not (set(b) & excl)]
    e_singles  = [c for c in E_CORE_SINGLES if c not in excl]
    lp_singles = [c for c in LP_CORE_SINGLES if c not in excl]

    pool: list[int] = []
    for block in p_blocks:
        if len(pool) >= n_cores:
            break
        pool.extend(block)

    for cpu in e_singles:
        if len(pool) >= n_cores:
            break
        pool.append(cpu)

    for cpu in lp_singles:
        if len(pool) >= n_cores:
            break
        pool.append(cpu)

    total_logical = os.cpu_count() or len(pool)
    return [c for c in pool if c < total_logical][:n_cores]


def partition_cpu_pool_into_slots(pool: list[int], n_slots: int) -> list[list[int]]:
    """
    Particiona un pool ya calculado (ej. por assign_cpu_pool()) en
    `n_slots` grupos contiguos de PUs, cada uno destinado a un sub-job
    paralelo de CENSO (ver MAAQC_CENSO_SLOTS junto a build_omp_env()).

    Contexto (julio 2026): CENSO congela su entorno una sola vez
    (ENVIRON en params.py) y no expone ningún mecanismo nativo para que
    sub-jobs paralelos (2+ SP simultáneos, cada uno con su propio
    mpirun) reciban listas de cores distintas — confirmado leyendo
    censo/parallel.py (LocalCluster con n_workers=1, un solo pool de
    hilos) y processing/orca_processor.py (_sp()/opt() llaman a
    _make_call() sin env= explícito, cayendo siempre a ENVIRON global).
    Un parche local en orca_processor.py (fuera de este repo, ver
    censo_patch/README.md) intercepta esas llamadas y lee
    MAAQC_CENSO_SLOTS para asignar cada conformero a un slot de forma
    determinística (hash estable de job.conf.name módulo n_slots).

    Esta función solo arma la partición del lado de pipeline.py — no
    reparte P_CORE_BLOCKS/E_CORE_SINGLES con ninguna lógica especial de
    prioridad (a diferencia de assign_cpu_lists()/assign_cpu_pool()):
    simplemente corta `pool` (ya ordenado por assign_cpu_pool() con
    P-cores primero) en `n_slots` trozos de tamaño lo más parejo
    posible, preservando el orden. Si `pool` no es divisible
    exactamente, los primeros slots reciben una PU extra antes que los
    últimos (ej. pool de 10 PUs en 3 slots → tamaños [4, 3, 3], no
    [3, 3, 4]) — irrelevante en la práctica porque el caso de uso real
    es n_slots=2 con pool par (8 PUs → [4, 4]).

    Devuelve una lista de `n_slots` listas de PUs. Si `pool` tiene menos
    elementos que `n_slots`, algunos slots devueltos pueden quedar
    vacíos (lista vacía) — el caller debe decidir si eso es aceptable
    (build_omp_env() con censo_slot_groups conteniendo un grupo vacío
    simplemente no fija MAAQC_PE_LIST para ese slot en el parche de
    CENSO, ver fail-open en _maaqc_build_slot_env()).
    """
    if n_slots <= 0:
        return []
    if not pool:
        return [[] for _ in range(n_slots)]

    base_size, remainder = divmod(len(pool), n_slots)
    slots: list[list[int]] = []
    idx = 0
    for i in range(n_slots):
        size = base_size + (1 if i < remainder else 0)
        slots.append(pool[idx: idx + size])
        idx += size
    return slots


def _check_cpu_budget(batch_size: int, nprocs_per_job: int, log: logging.Logger) -> None:
    """Advierte (no bloquea) si el lote pide más cores de los físicamente disponibles."""
    requested = batch_size * nprocs_per_job
    available = os.cpu_count() or 0
    if requested > available:
        log.warning(
            f"  ⚠  el lote pide {requested} cores ({batch_size}×{nprocs_per_job}) "
            f"pero la máquina solo tiene {available} — los rangos de taskset "
            f"se solaparán o excederán los cores reales. Reduce --parallel o "
            f"--nprocs-per-job."
        )




# ─────────────────────────────────────────────
#  STATUS GLOBAL
# ─────────────────────────────────────────────

def print_status():
    """Imprime estado de todos los cálculos en PAPER_ROOT/<grupo>/runs/<mol_name>/.

    PAPER_ROOT se auto-detecta (ver _find_paper_root): usa
    $PIPELINE_ROOT si está seteada, o sube desde el cwd buscando
    README.md/LICENSE/tasks.json/scripts/. Si el resultado no es el
    esperado, corre desde la raíz del proyecto o exporta PIPELINE_ROOT.
    """
    if not WORK_ROOT.exists():
        print(f"No se encontró el directorio '{WORK_ROOT}/'. Aún no hay cálculos.")
        return

    print(f"\n{'─'*60}")
    print(f"  Estado de runs en: {WORK_ROOT.resolve()}")
    print(f"  (raíz auto-detectada — usa $PIPELINE_ROOT para forzar otra)")
    print(f"{'─'*60}")

    # Estructura: PAPER_ROOT/<grupo>/runs/<mol_name>/<etapas>
    # (p.ej. PAPER1/MGLYVAL/runs/D/01_crest/, .../D/02_opt/, etc.)
    for group_dir in sorted(WORK_ROOT.iterdir()):
        if not group_dir.is_dir():
            continue
        runs_dir = group_dir / "runs"
        if not runs_dir.is_dir():
            continue
        for mol_dir in sorted(runs_dir.iterdir()):
            if not mol_dir.is_dir():
                continue
            print(f"\n  📁 {group_dir.name}/{mol_dir.name}")
            for step_num in ["01_crest", "02_opt", "03_numfreq", "04_sp", "04b_sp_gas", "04c_sp_smd", "05_tddft", "06_sp_nroot", "07_adc2"]:
                step_dir = mol_dir / step_num
                if not step_dir.exists():
                    print(f"     {step_num:12s}  —  no iniciado")
                    continue

                conf_dirs = [d for d in step_dir.iterdir()
                             if d.is_dir() and d.name.startswith("conf_")]
                if conf_dirs:
                    ok = sum(1 for d in conf_dirs
                             if any(check_success(f) for f in d.glob("*.out")))
                    print(f"     {step_num:12s}  ✓ {ok}/{len(conf_dirs)} confs")
                else:
                    xyz = step_dir / "crest_conformers.xyz"
                    if xyz.exists() and xyz.stat().st_size > 0:
                        try:
                            n = len(read(str(xyz), index=":"))
                            print(f"     {step_num:12s}  ✓ {n} confórmeros")
                        except Exception:
                            print(f"     {step_num:12s}  ⚠  XYZ presente pero ilegible")
                    else:
                        log_f = step_dir / "crest.log"
                        status = "en progreso?" if log_f.exists() else "no iniciado"
                        print(f"     {step_num:12s}  — {status}")

    print(f"\n{'─'*60}\n")


# ─────────────────────────────────────────────
#  PUNTO DE ENTRADA
# ─────────────────────────────────────────────

def _setup_run_resources(args):
    """
    Valida --multitemp/--trange/--skipnumfreq, detecta/resuelve los
    recursos de cómputo (cores/RAM por etapa, vía ResourceProfile) y
    actualiza los globales usados por las funciones stage_*.

    Se separó de _run_one() para poder compartirlo entre una corrida
    normal (una molécula) y el modo ciclo de protonación automático
    (--autoprotonate + --autodeprotonate juntos, ver
    _run_protonation_cycle), que corre tres pipelines completos (R/A/B)
    con los MISMOS recursos de cómputo detectados una sola vez.

    Devuelve (censo_temperatures, resources) donde resources es un dict
    con los valores ya resueltos (nprocs_crest, nprocs, maxcore, ...)
    para poder loggearlos en el resumen de cada corrida.
    """
    global NPROCS_CREST, NPROCS, MAXCORE, NPROCS_NUMFREQ, MAXCORE_NUMFREQ, \
           NPROCS_TDDFT, MAXCORE_TDDFT, NPROCS_ADC2, MAXCORE_ADC2, \
           OPT_BATCH_SIZE, SP_BATCH_SIZE, TDDFT_BATCH_SIZE, OMP_ENV

    # ── Validar --multitemp / --trange ───────────────────────────────────────
    if args.trange is not None and not args.multitemp:
        sys.exit(
            "ERROR: --trange solo tiene efecto junto con --multitemp."
        )
    if args.multitemp and args.nroot:
        sys.exit(
            "ERROR: --multitemp es incompatible con --nroot. "
            "Corre --nroot como paso independiente sobre un run ya completado."
        )
    # En modo --multitemp: --trange es opcional; default = 274 298 24
    # En modo normal: temperatura única = --temperature (default 298.15)
    if args.multitemp:
        t_ini, t_fin, n_pts = args.trange if args.trange is not None else (274.0, 298.0, 24)
        n_pts = int(round(n_pts))
        if n_pts < 1:
            sys.exit("ERROR: --trange N debe ser ≥ 1.")
        if n_pts == 1:
            censo_temperatures = [t_ini]
        else:
            import numpy as _np
            censo_temperatures = list(_np.linspace(t_ini, t_fin, n_pts))
        # Validar que --temperature cae dentro del rango generado
        ref_temp = args.temperature
        temp_set = {round(t, 8) for t in censo_temperatures}
        # Buscar la temperatura más cercana en el rango generado
        closest = min(censo_temperatures, key=lambda t: abs(t - ref_temp))
        if abs(closest - ref_temp) > 1e-4:
            # La temperatura de referencia no está exactamente en el linspace;
            # añadirla al rango para que CENSO la calcule igualmente.
            censo_temperatures = sorted(set(censo_temperatures + [ref_temp]))
            import warnings as _w
            _w.warn(
                f"--temperature {ref_temp} K no coincide exactamente con ningún "
                f"punto de --trange; se añade al rango. "
                f"Punto más cercano era {closest:.4f} K (Δ={abs(closest - ref_temp):.4f} K).",
                stacklevel=1,
            )
    else:
        # Modo normal: temperatura única
        censo_temperatures = [args.temperature]

    if args.skipnumfreq and args.start in ("numfreq", "sp", "tddft"):
        sys.exit(
            f"ERROR: --skipnumfreq es incompatible con --start {args.start}. "
            f"--skipnumfreq ya implica saltar NUMFREQ desde OPT; usa solo "
            f"--start opt (o --start crest) junto con --skipnumfreq."
        )

    # ── Resolver recursos: autodetección + flags CLI ──────────────────────────
    #
    # Lógica de precedencia (mayor a menor):
    #   1. Flags avanzados explícitos (--nprocs, --maxcore, --nprocs-numfreq, ...)
    #   2. --cores N  → ResourceProfile.detect_with_batch(N, --parallel)
    #   3. Sin --cores → ResourceProfile.detect()  (lee topología real del sistema)
    #
    # ResourceProfile calcula por etapa:
    #   OPT     : OPT_BATCH_SIZE × MAX_OPT_NPROCS_PER_JOB  (3×4 = 12 cores)
    #   NUMFREQ : 1 × MAX_NUMFREQ_NPROCS                   (1×12 = 12 cores)
    #   SP      : SP_BATCH_SIZE × MAX_SP_NPROCS             (3×4 = 12 cores)
    #   TDDFT   : TDDFT_BATCH_SIZE × TDDFT_NPROCS_PER_JOB  (2×6 = 12 cores)
    #   ADC2    : fijo ≤ NPROCS_ADC2_MAX + RAM generosa por core

    if args.cores is not None:
        rp = ResourceProfile.detect_with_batch(
            total_cores = args.cores,
            batch_size  = 1,
        )
    else:
        rp = ResourceProfile.detect(batch_size=1)

    # Extraer valores del perfil autodetectado
    nprocs_crest   = rp.nprocs_crest
    nprocs         = rp.nprocs_opt
    maxcore        = rp.maxcore_opt
    nprocs_numfreq = rp.nprocs_numfreq
    maxcore_nf     = rp.maxcore_numfreq

    # Sobreescritura por flags manuales (mayor precedencia)
    if args.nprocs is not None:
        nprocs = args.nprocs
    if args.maxcore is not None:
        maxcore = args.maxcore
    if args.nprocs_numfreq is not None:
        nprocs_numfreq = args.nprocs_numfreq
    if args.maxcore_numfreq is not None:
        maxcore_nf = args.maxcore_numfreq

    # Valores derivados para TDDFT y ADC2 (no expuestos en CLI; vienen siempre del perfil)
    nprocs_tddft   = rp.nprocs_tddft
    maxcore_tddft  = rp.maxcore_tddft
    nprocs_adc2    = rp.nprocs_adc2
    maxcore_adc2   = rp.maxcore_adc2
    opt_batch_size = rp.OPT_BATCH_SIZE
    sp_batch_size  = rp.SP_BATCH_SIZE
    tddft_batch    = rp.TDDFT_BATCH_SIZE

    # Actualizar globales (usados en las funciones stage_*)
    NPROCS_CREST     = nprocs_crest
    NPROCS           = nprocs
    MAXCORE          = maxcore
    NPROCS_NUMFREQ   = nprocs_numfreq
    MAXCORE_NUMFREQ  = maxcore_nf
    NPROCS_TDDFT     = nprocs_tddft
    MAXCORE_TDDFT    = maxcore_tddft
    NPROCS_ADC2      = nprocs_adc2
    MAXCORE_ADC2     = maxcore_adc2
    OPT_BATCH_SIZE   = opt_batch_size
    SP_BATCH_SIZE    = sp_batch_size
    TDDFT_BATCH_SIZE = tddft_batch
    OMP_ENV          = build_omp_env(NPROCS)

    resources = dict(
        nprocs_crest=nprocs_crest, nprocs=nprocs, maxcore=maxcore,
        nprocs_numfreq=nprocs_numfreq, maxcore_nf=maxcore_nf,
        nprocs_tddft=nprocs_tddft, maxcore_tddft=maxcore_tddft,
        nprocs_adc2=nprocs_adc2, maxcore_adc2=maxcore_adc2,
        opt_batch_size=opt_batch_size, sp_batch_size=sp_batch_size,
        tddft_batch=tddft_batch,
    )
    return censo_temperatures, resources


def _run_one(args, parser):
    """
    Corre el pipeline completo (CREST -> CENSO -> OPT -> NumFreq -> SP -> TDDFT -> ADC2)
    para UNA molecula/carga/configuracion ya resuelta en `args`.

    Se separo de main() para poder encadenar varias moleculas en un solo
    comando (--mol repetido, ver _split_argv_by_mol), corriendo el pipeline
    consecutivamente para cada una en el orden en que se escribieron en
    la terminal.
    """
    censo_temperatures, _res = _setup_run_resources(args)
    nprocs_crest   = _res["nprocs_crest"]
    nprocs         = _res["nprocs"]
    maxcore        = _res["maxcore"]
    nprocs_numfreq = _res["nprocs_numfreq"]
    maxcore_nf     = _res["maxcore_nf"]
    nprocs_tddft   = _res["nprocs_tddft"]
    maxcore_tddft  = _res["maxcore_tddft"]
    opt_batch_size = _res["opt_batch_size"]
    sp_batch_size  = _res["sp_batch_size"]
    tddft_batch    = _res["tddft_batch"]

    if args.status:
        print_status()
        return

    if args.mol is None:
        parser.error("--mol es obligatorio salvo que uses --status")

    mol_xyz = args.mol.resolve()
    if not mol_xyz.exists():
        sys.exit(f"ERROR: no se encontró {mol_xyz}")

    # ── Modo batch SP_GAS sobre carpeta de grupo (--spgas + --mol carpeta) ──
    # No es un .xyz: --mol apunta directo a la carpeta de grupo (ej. MGLYVAL),
    # y se recorren TODOS los runs/ dentro de ella calculando el SP_GAS que
    # falte, sin correr ningún otro paso del pipeline.
    if args.spgas and mol_xyz.is_dir():
        group_dir = mol_xyz
        log = setup_logging(group_dir)
        log.info("=" * 60)
        log.info(f"MODO: SP_GAS batch (--spgas) sobre carpeta de grupo → {group_dir.name}")
        log.info("=" * 60)
        stage_spgas_batch(group_dir, log)
        log.info("=" * 60)
        log.info("SP_GAS BATCH FINALIZADO")
        log.info("=" * 60)
        return

    # ── Modo batch SP_SMD sobre carpeta de grupo (--spsmd + --mol carpeta) ──
    # Mismo patrón que el modo batch de --spgas: --mol apunta directo a la
    # carpeta de grupo (ej. MGLYVAL), y se recorren TODOS los runs/ dentro
    # de ella calculando el SP_SMD que falte, sin correr ningún otro paso.
    if args.spsmd and mol_xyz.is_dir():
        group_dir = mol_xyz
        log = setup_logging(group_dir)
        log.info("=" * 60)
        log.info(f"MODO: SP_SMD batch (--spsmd) sobre carpeta de grupo → {group_dir.name}")
        log.info("=" * 60)
        stage_spsmd_batch(group_dir, log, boltzmann_min_pct=args.spsmd_boltzmann_threshold)
        log.info("=" * 60)
        log.info("SP_SMD BATCH FINALIZADO")
        log.info("=" * 60)
        return

    if not mol_xyz.suffix.lower() == ".xyz":
        sys.exit(f"ERROR: el archivo de entrada debe ser un .xyz: {mol_xyz}")

    # Validar binarios antes de empezar
    if not ORCA_BIN.exists():
        sys.exit(f"ERROR: no se encontró el ejecutable de ORCA en: {ORCA_BIN}\n"
                 f"       Ajusta ORCA_BIN al inicio del script.")

    mol_name = mol_xyz.stem

    # Convención de carpetas: <grupo>/molecules/<mol_name>.xyz  →  <grupo>/runs/<mol_name>/
    # (p.ej. PAPER1/MGLYVAL/molecules/D.xyz  →  PAPER1/MGLYVAL/runs/D/).
    # Si el .xyz no vive dentro de una carpeta "molecules/", se usa su carpeta
    # contenedora directamente como "grupo" (fallback razonable, evita romper
    # invocaciones con otra estructura de carpetas).
    if mol_xyz.parent.name == "molecules":
        group_dir = mol_xyz.parent.parent
    else:
        group_dir = mol_xyz.parent

    mol_dir = group_dir / "runs" / mol_name
    mol_dir.mkdir(parents=True, exist_ok=True)

    log = setup_logging(mol_dir)

    # ── Resolución de carga (auto por sufijo, o --charge explícito) ──────
    # args.charge llega en None si el usuario no pasó --charge: se infiere
    # aquí y se sobrescribe args.charge con el valor resuelto, de modo que
    # TODO el resto de _run_one()/_run_downstream_from_crest() (que ya
    # asumían args.charge como int) siga funcionando sin más cambios.
    args.charge = resolve_charge(mol_name, args.charge, log)

    log.info("=" * 60)
    log.info(f"PIPELINE ORCA  —  molécula: {mol_name}")
    log.info(f"Carga: {args.charge}   Mult: {MULT}")
    log.info(f"Flujo: CREST → CENSO (prescreen) → OPT (nconf auto) → NumFreq → SP → TDDFT → ADC2")
    log.info(f"Inicio desde etapa: {args.start.upper()}")
    log.info(f"ORCA: {ORCA_BIN}")
    # ── Temperatura / Multitemp ──────────────────────────────────────────────
    if args.multitemp:
        t_ini_str = f"{censo_temperatures[0]:.2f} K"
        t_fin_str = f"{censo_temperatures[-1]:.2f} K"
        log.info(
            f"Temperatura CENSO      : MULTITEMP — {len(censo_temperatures)} puntos "
            f"[{t_ini_str} … {t_fin_str}]  "
            f"(OPT usa T={args.temperature:.2f} K como referencia)"
        )
    else:
        log.info(f"Temperatura CENSO      : {censo_temperatures[0]:.2f} K")
    log.info(f"Recursos CREST/CENSO   : {nprocs_crest} cores  (xTB --T, escala lineal)")
    log.info(f"Recursos OPT           : {opt_batch_size}×{nprocs} cores / {maxcore} MB  "
             f"→ {nprocs * maxcore / 1024:.1f} GB/job")
    log.info(f"Recursos SP            : {sp_batch_size}×{nprocs} cores / {maxcore} MB  "
             f"→ {nprocs * maxcore / 1024:.1f} GB/job")
    log.info(f"Recursos TDDFT         : {tddft_batch}×{nprocs_tddft} cores / {maxcore_tddft} MB  "
             f"→ {nprocs_tddft * maxcore_tddft / 1024:.1f} GB/job  (NRoots 15)")
    log.info(f"Recursos NumFreq       : 1×{nprocs_numfreq} cores / {maxcore_nf} MB  "
             f"→ {nprocs_numfreq * maxcore_nf / 1024:.1f} GB")
    log.info("=" * 60)

    # ── Verificación temprana de binarios requeridos ─────────────────────
    # Antes, un 'censo' no encontrado solo se descubría dentro de
    # stage_censo_prescreen() vía FileNotFoundError — después de que CREST
    # ya hubiera corrido (p.ej. 70+ min de cómputo). check_required_bins()
    # aborta aquí, antes de gastar tiempo de cómputo, con mensaje claro
    # sobre activar el conda env 'sciencebitch'.
    check_required_bins(log)

    start_idx = STEPS.index(args.start)

    # ════════════════════════════════════════════════════════════
    #  MODO SP + NROOTS 1  (--nroot)
    # ════════════════════════════════════════════════════════════
    if args.nroot:
        log.info("MODO: SP + NRoots 1 (--nroot)")
        log.info(f"Top-N confórmeros por E_OPT: {args.nroot_top}")
        nroot_dirs = stage_sp_nroot(mol_dir, args.charge, log, n_top=args.nroot_top)
        if not nroot_dirs:
            log.error("[SP_NROOT] Ningún cálculo completado.")
            sys.exit(1)
        log.info("=" * 60)
        log.info("PIPELINE (SP_NROOT) FINALIZADO")
        log.info(f"Resultados en: {(mol_dir / '06_sp_nroot').resolve()}")
        log.info("=" * 60)
        return

    # ════════════════════════════════════════════════════════════
    #  MODO ADC2 standalone (--start adc2)
    # ════════════════════════════════════════════════════════════
    # ADC2 solo depende de 04_sp/ (energías + opt_opt.xyz de cada confórmero).
    # Se salta CREST, CENSO, OPT, NUMFREQ, TDDFT por completo — ni siquiera
    # se verifica que existan — porque en el flujo real de Kmilo esos
    # archivos ya fueron borrados del disco (ADC2 requiere ~400GB libres,
    # así que respalda todo en un disco externo y limpia el compu antes de
    # correr esta etapa). Lo único indispensable es 04_sp/.
    if args.start == "adc2":
        log.info("MODO: ADC2 standalone (--start adc2)")
        log.info("[ADC2] saltando CREST/CENSO/OPT/NUMFREQ/TDDFT — "
                 "solo se usa 04_sp/ (confórmero de menor energía)")
        sp_root = mol_dir / "04_sp"
        if not sp_root.exists():
            sys.exit(
                f"ERROR: --start adc2 requiere que exista {sp_root}\n"
                f"       (los .out y opt_opt.xyz de cada confórmero SP). "
                f"Corre primero la etapa SP."
            )
        adc2_dir = stage_adc2(mol_dir, args.charge, log, mol_name=mol_name)
        if adc2_dir is None:
            log.error("[ADC2] No se pudo completar la validación ADC2.")
            sys.exit(1)
        log.info("=" * 60)
        log.info(f"[ADC2] validación completada ✓ → {adc2_dir}")
        log.info("PIPELINE (ADC2 standalone) FINALIZADO")
        log.info(f"Resultados en: {(mol_dir / '07_adc2').resolve()}")
        log.info("=" * 60)
        return

    # ════════════════════════════════════════════════════════════
    #  MODO NORMAL
    # ════════════════════════════════════════════════════════════

    # ── Etapa 1: CREST ──────────────────────────────────────────
    # CREST siempre corre cuando start_idx == 0 (default).
    # Si se usa --start opt/numfreq/sp/tddft/adc2, se asume que
    # 01_crest/crest_conformers.xyz ya existe.
    conformers_xyz = mol_dir / "01_crest" / "crest_conformers.xyz"

    if start_idx == 0:
        conformers_xyz = stage_crest(
            mol_xyz, mol_dir, args.charge, log,
            autoprotonate=args.autoprotonate,
            autodeprotonate=args.autodeprotonate,
        )
        if conformers_xyz is None:
            log.error("CREST falló. Abortando pipeline.")
            sys.exit(1)
    else:
        if not conformers_xyz.exists():
            sys.exit(
                f"ERROR: se pidió start={args.start} pero no existe "
                f"{conformers_xyz}\n"
                f"       Corre primero desde el inicio (sin --start) para "
                f"generar el ensemble CREST+CENSO."
            )
        if args.autoprotonate or args.autodeprotonate:
            log.warning(
                "[CREST] --autoprotonate/--autodeprotonate no tienen efecto "
                "aquí: CREST no corre en esta invocación "
                f"(--start {args.start} indica etapa posterior)."
            )
        log.info("[CREST] saltado (--start indica etapa posterior)")

    _run_downstream_from_crest(
        args, conformers_xyz, mol_dir, args.charge, mol_xyz, mol_name,
        log, censo_temperatures, start_idx,
    )


def _run_downstream_from_crest(args, conformers_xyz, mol_dir, charge,
                                natoms_ref_xyz, mol_name, log,
                                censo_temperatures, start_idx):
    """
    Corre CENSO (prescreen) → OPT → NumFreq → SP → [SP_GAS] → [SP_SMD] →
    TDDFT → [ADC2] sobre un ensemble CREST ya generado (conformers_xyz).

    Se separó del cuerpo de _run_one() para poder reutilizarlo tanto en
    una corrida normal (una molécula/carga) como en el modo ciclo de
    protonación automático (--autoprotonate + --autodeprotonate juntos,
    ver _run_protonation_cycle), que llama a esta función tres veces —
    una por estado de protonación (R/A/B) — con su propio mol_dir/charge/
    natoms_ref_xyz, pero compartiendo el resto de opciones de `args`
    (--skipnumfreq, --adc2, --spgas, --spsmd, --multitemp,
    --boltzmann-threshold, --temperature, etc.) tal cual.

    natoms_ref_xyz es la molécula base usada SOLO para derivar
    expected_natoms (sanity check de que conformers_xyz corresponde a
    ESTE estado/carga) — en una corrida normal es el mismo mol_xyz de
    entrada; en modo ciclo, para los estados A/B es la estructura
    representativa (des)protonada, no la molécula neutra original.
    """
    # ── Etapa 1b: CENSO pre-screening (siempre obligatorio) ──────
    # Corre entre CREST y OPT. Fail-open: si CENSO no está disponible
    # o falla, el pipeline continúa con el ensemble CREST sin filtrar.
    # En retoma (start > opt) se omite si ya existe el ensemble filtrado.
    #
    # Modo --multitemp: corre CENSO para cada temperatura en censo_temperatures
    # y guarda los ensembles en subdirectorios separados (T<temp>K/).
    # El pipeline continúa con el ensemble de --temperature (default 298.15 K)
    # para seleccionar los confórmeros que pasan a OPT; los demás ensembles
    # quedan disponibles para análisis Boltzmann multi-temperatura con
    # results_collector.py.
    if start_idx <= STEPS.index("opt"):
        if args.multitemp:
            multitemp_results = stage_censo_multitemp(
                conformers_xyz=conformers_xyz,
                mol_dir=mol_dir,
                charge=charge,
                log=log,
                temperatures=censo_temperatures,
            )
            # Seleccionar el ensemble cuya temperatura coincide con --temperature.
            # Como censo_temperatures puede contener floats de linspace con
            # imprecisión numérica, buscamos el punto más cercano.
            ref_temp = args.temperature
            opt_temp = min(multitemp_results.keys(), key=lambda t: abs(t - ref_temp))
            conformers_xyz = multitemp_results[opt_temp]
            log.info(
                f"[MULTITEMP] ensemble de referencia para OPT: "
                f"T={opt_temp:.2f} K (--temperature={ref_temp:.2f} K) "
                f"→ {conformers_xyz}"
            )
        else:
            conformers_xyz = stage_censo_prescreen(
                conformers_xyz, mol_dir, charge, log,
                temperature=censo_temperatures[0],
            )

    # ── Etapa 2: OPT ────────────────────────────────────────────
    # El número de confórmeros se determina según el modo activo:
    #   --boltzmann-threshold PCT : toma confs con peso Boltzmann ≥ PCT%
    #                               (requiere 1_SCREENING.out de CENSO)
    #   modo automático (default) : select_nconf_auto basado en spread de energías
    if start_idx <= STEPS.index("opt"):
        boltz_conf_ids: list[str] = []
        if args.boltzmann_threshold is not None:
            # Buscar 1_SCREENING.out en las ubicaciones estándar de CENSO
            _censo_prescreen_dir = mol_dir / "01b_censo_prescreen"
            _screening_candidates = [
                _censo_prescreen_dir / "1_SCREENING.out",
                _censo_prescreen_dir / "CENSO" / "1_SCREENING.out",
                _censo_prescreen_dir / "prescreening" / "1_SCREENING.out",
            ]
            _screening_out = next(
                (p for p in _screening_candidates if p.exists()), None
            )
            if _screening_out:
                nconf, boltz_conf_ids = select_nconf_boltzmann(
                    conformers_xyz, _screening_out, log,
                    threshold=args.boltzmann_threshold,
                )
            else:
                log.debug(
                    f"[NCONF-BOLTZMANN] 1_SCREENING.out no disponible en "
                    f"{_censo_prescreen_dir} — usando select_nconf_auto"
                )
                nconf = 0

            if nconf == 0:
                nconf = select_nconf_auto(conformers_xyz, log)
                boltz_conf_ids = []
        else:
            nconf = select_nconf_auto(conformers_xyz, log)

        # expected_natoms se deriva de natoms_ref_xyz: es la molécula base
        # de ESTE estado de protonación (A/B/C/D/R), la misma que arrancó
        # CREST para este mol_dir/charge (o, en modo ciclo automático, la
        # estructura representativa (des)protonada de ese estado). Si
        # conformers_xyz apunta al ensemble de otro estado, este número no
        # va a coincidir, y stage_opt() lo detecta antes de escribir
        # ningún opt.inp.
        try:
            _mol_atoms = read(str(natoms_ref_xyz))
            _expected_natoms = len(_mol_atoms)
        except Exception as e:
            log.warning(
                f"[OPT] no se pudo leer {natoms_ref_xyz} para derivar "
                f"expected_natoms ({e}); se omite la validación estado/carga."
            )
            _expected_natoms = None

        opt_dirs = stage_opt(
            conformers_xyz, mol_dir, charge, nconf, log,
            conf_ids=boltz_conf_ids if boltz_conf_ids else None,
            expected_natoms=_expected_natoms,
        )
        if not opt_dirs:
            log.error("Ningún confórmero optimizado correctamente. Abortando.")
            sys.exit(1)
    else:
        opt_dirs = sorted((mol_dir / "02_opt").glob("conf_*")) \
                   if (mol_dir / "02_opt").exists() else []
        log.info(f"[OPT] saltado — encontrados {len(opt_dirs)} directorios previos")

    # ── Etapa 3: NUMFREQ ────────────────────────────────────────
    if args.skipnumfreq:
        log.warning(
            "[NUMFREQ] SALTADO por --skipnumfreq — SP/TDDFT usarán "
            "geometría de 02_opt/ directamente. Sin correcciones térmicas "
            "de NumFreq en el ciclo termodinámico de pKa."
        )
        freq_dirs = opt_dirs
    elif start_idx <= STEPS.index("numfreq"):
        freq_dirs = stage_numfreq(opt_dirs, mol_dir, charge, log)
        if not freq_dirs:
            log.error("Ningún cálculo de frecuencias completado. Abortando.")
            sys.exit(1)
    else:
        freq_dirs = sorted((mol_dir / "03_numfreq").glob("conf_*")) \
                    if (mol_dir / "03_numfreq").exists() else []
        log.info(f"[NUMFREQ] saltado — encontrados {len(freq_dirs)} directorios previos")

    # ── Etapa 4: SP (fase gas) — OBLIGATORIA ─────────────────────
    # SP es la etapa mínima requerida para dar el pipeline por finalizado.
    # Si no se completa ningún confórmero, se aborta antes de llegar a
    # "PIPELINE FINALIZADO" (ver sys.exit(1) abajo).
    if start_idx <= STEPS.index("sp"):
        sp_dirs = stage_sp(freq_dirs, mol_dir, charge, log)
        if not sp_dirs:
            log.error("Ningún SP completado. Abortando "
                      "(SP fase gas es obligatorio para finalizar el pipeline).")
            sys.exit(1)
    else:
        sp_dirs = sorted((mol_dir / "04_sp").glob("conf_*")) \
                  if (mol_dir / "04_sp").exists() else []
        log.info(f"[SP] saltado — encontrados {len(sp_dirs)} directorios previos")
        if not sp_dirs:
            log.error("[SP] no hay directorios previos disponibles. Abortando "
                      "(SP fase gas es obligatorio para finalizar el pipeline).")
            sys.exit(1)

    # ── Etapa 4b: SP_GAS (OPCIONAL — requiere --spgas) ───────────
    # Idéntico al SP (Etapa 4) pero sin CPCM(water). No bloquea
    # "PIPELINE FINALIZADO" — es aditivo, igual que --adc2.
    if args.spgas:
        log.info(f"[SP_GAS] calculando SP fase-gas (sin CPCM) para "
                 f"{len(sp_dirs)} confórmeros de 04_sp/")
        spgas_dirs = stage_sp_gas(sp_dirs, mol_dir, log)
        if not spgas_dirs:
            log.warning("[SP_GAS] ningún cálculo completado.")
        else:
            log.info(f"[SP_GAS] {len(spgas_dirs)}/{len(sp_dirs)} completados ✓")
    else:
        log.info("[SP_GAS] no incluido en esta ejecución (opcional — usa --spgas)")

    # ── Etapa 4c: SP_SMD (OPCIONAL — requiere --spsmd) ────────────
    # Idéntico al SP (Etapa 4) pero con SMD(water) en vez de CPCM(water).
    # No bloquea "PIPELINE FINALIZADO" — es aditivo, igual que --spgas/--adc2.
    if args.spsmd:
        log.info(f"[SP_SMD] calculando SP con SMD(water) para "
                 f"{len(sp_dirs)} confórmeros de 04_sp/ "
                 f"(filtro Boltzmann > {args.spsmd_boltzmann_threshold}%)")
        spsmd_dirs = stage_sp_smd(sp_dirs, mol_dir, log,
                                   boltzmann_min_pct=args.spsmd_boltzmann_threshold)
        if not spsmd_dirs:
            log.warning("[SP_SMD] ningún cálculo completado.")
        else:
            log.info(f"[SP_SMD] {len(spsmd_dirs)}/{len(sp_dirs)} completados ✓")
    else:
        log.info("[SP_SMD] no incluido en esta ejecución (opcional — usa --spsmd)")

    # ── Etapa 5: TDDFT ──────────────────────────────────────────
    if start_idx <= STEPS.index("tddft"):
        tddft_dirs = stage_tddft(sp_dirs, mol_dir, charge, log)
        if not tddft_dirs:
            log.warning("Ningún TDDFT completado.")
        else:
            log.info(f"[TDDFT] {len(tddft_dirs)} espectros calculados ✓")
    else:
        log.info("[TDDFT] no incluido en esta ejecución")

    # ── Etapa 6: ADC2 (OPCIONAL — requiere --adc2) ───────────────
    # ADC2 ya NO corre por defecto. Solo se ejecuta si el usuario lo pide
    # explícitamente con --adc2 (o de forma standalone con --start adc2,
    # manejado en su propio bloque más arriba). El pipeline se considera
    # FINALIZADO con solo SP (obligatorio) completado; ADC2 es un extra.
    if args.adc2 and start_idx <= STEPS.index("adc2"):
        adc2_dir = stage_adc2(mol_dir, charge, log, mol_name=mol_name)
        if adc2_dir is None:
            log.warning("[ADC2] No se pudo completar la validación ADC2.")
        else:
            log.info(f"[ADC2] validación completada ✓ → {adc2_dir}")
    else:
        log.info("[ADC2] no incluido en esta ejecución (opcional — usa --adc2 para activarla)")

    log.info("=" * 60)
    log.info("PIPELINE FINALIZADO  (requisito obligatorio SP fase-gas: OK)")
    log.info(f"Resultados en: {mol_dir.resolve()}")
    log.info("=" * 60)




def _split_argv_by_mol(argv):
    """
    Divide argv en:
      - preamble: tokens ANTES del primer --mol (flags "globales" que
        aplican a todas las moleculas: --status, --cores, --nprocs,
        --maxcore, --nprocs-numfreq, --maxcore-numfreq, etc.)
      - mol_blocks: lista de bloques. Cada bloque empieza en un token
        "--mol" (o "--mol=...") y contiene todos los flags que le
        siguen hasta el proximo "--mol" (o el final de argv).

    Esto permite encadenar varias moleculas en un solo comando, donde
    cada una lleva sus propias especificaciones (--charge, --start,
    --skipnumfreq, --nroot, --multitemp, --boltzmann-threshold, etc.)
    respetando el orden en que se escriben en la terminal, por ejemplo:

        python pipeline.py --mol molecules/B.xyz --charge -1 \\
                            --mol molecules/C.xyz --charge -2

    Aqui --charge -1 pertenece solo a B.xyz y --charge -2 solo a C.xyz.
    Un flag repetido dentro de un mismo bloque (p.ej. --cores puesto de
    nuevo despues de un --mol) sobreescribe el valor del preamble para
    esa molecula en particular, ya que argparse siempre se queda con la
    ultima ocurrencia de un flag.
    """
    def _is_mol_flag(tok):
        return tok == "--mol" or tok.startswith("--mol=")

    mol_indices = [i for i, tok in enumerate(argv) if _is_mol_flag(tok)]
    if not mol_indices:
        return argv, []

    preamble = argv[:mol_indices[0]]
    blocks = []
    for j, start in enumerate(mol_indices):
        end = mol_indices[j + 1] if j + 1 < len(mol_indices) else len(argv)
        blocks.append(argv[start:end])
    return preamble, blocks


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Pipeline ORCA: CREST → CENSO (prescreen) → Opt → NumFreq → SP → TDDFT → ADC2\n"
            "CREST y CENSO pre-screening son obligatorios. El número de confórmeros para OPT\n"
            "se determina automáticamente (mín 3, máx 20) según el panorama conformacional."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos (modo normal):
  # Pipeline completo desde cero, 8 cores totales
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --cores 8

  # Retomar desde OPT (CREST+CENSO ya hechos en 01_crest/ y 01b_censo_prescreen/)
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --start opt

  # Retomar desde SP (opt y freq ya hechos)
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --start sp

  # Retomar NumFreq interrumpida (detección automática de restart)
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --start numfreq
  # → si hay freq_D*.res.Gradients en la carpeta, activa Restart true solo

  # SP+NRoots 1 sobre los 5 mejores confórmeros por E_OPT (requiere 02_opt/ completo)
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --nroot
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --nroot --nroot-top 3

  # ADC2 standalone: SOLO requiere 04_sp/ (ignora crest/opt/numfreq/tddft
  # aunque no existan en disco). Útil cuando ya se borró todo salvo 04_sp/
  # para liberar espacio (ADC2 necesita ~400GB libres).
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --start adc2

  # Pipeline normal + ADC2 al final (ADC2 es OPCIONAL: sin --adc2 el
  # pipeline se da por finalizado apenas SP, que es obligatorio, completa).
  python pipeline.py --mol molecules/mglyval.xyz --charge -1 --adc2

  # Ver estado de todos los runs
  python pipeline.py --status

  # Correr SOLO CENSO sobre un ensemble ya existente (reemplaza invocar
  # `censo` a mano en la terminal — arma OMP_NUM_THREADS/MKL/OPENBLAS y
  # taskset igual que el pipeline completo, evita saturación de cores):
  python pipeline.py --censo-standalone runs/A/01_crest/crest_conformers.xyz \\
                      --charge -1 --temperature 274.15 --cores 12

Ejemplos (modo multi-molécula, --mol repetido):
  # Corre B.xyz (carga -1) y luego C.xyz (carga -2), consecutivamente.
  # Cada molécula toma las especificaciones que la siguen en la terminal
  # hasta el próximo --mol.
  python pipeline.py --mol molecules/B.xyz --charge -1 \\
                      --mol molecules/C.xyz --charge -2

  # Igual, pero C.xyz además salta NumFreq y B.xyz retoma desde OPT:
  python pipeline.py --mol molecules/B.xyz --charge -1 --start opt \\
                      --mol molecules/C.xyz --charge -2 --skipnumfreq

  # Flags globales (--cores, etc.) puestos ANTES del primer --mol se
  # heredan en todas las moléculas de la cola:
  python pipeline.py --cores 8 \\
                      --mol molecules/B.xyz --charge -1 \\
                      --mol molecules/C.xyz --charge -2
        """,
    )
    parser.add_argument("--mol",    type=Path, help="ruta al .xyz de entrada")
    parser.add_argument(
        "--charge", type=int, default=None,
        help=(
            "carga formal de la molécula. OPCIONAL: si se omite, se "
            "infiere automáticamente del nombre del archivo según el "
            "sufijo de estado de protonación (…A=+1, …B=-1, …C=-2, "
            "…D=0, …R=0). Pasa --charge explícito para forzar otro "
            "valor o si el nombre no sigue esa convención."
        ),
    )
    _protonation_group = parser.add_mutually_exclusive_group()
    _protonation_group.add_argument(
        "--autoprotonate", action="store_true",
        help=(
            "Solo afecta la etapa CREST inicial (01_crest/, requiere "
            "--start crest, el default): en vez de muestreo conformacional "
            "normal, corre CREST en modo --protonate para generar y "
            "rankear sitios de protonación candidatos. Se mantienen "
            "--T/--chrg/--uhf 0/--alpb water igual que el CREST normal — "
            "solo cambia el flag de modo. El resultado (protonated.xyz) se "
            "copia a 01_crest/crest_conformers.xyz. ATENCIÓN: las "
            "estructuras resultantes tienen carga REAL = --charge + 1; "
            "este flag NO propaga ese cambio a las etapas posteriores "
            "(OPT/SP/...), que siguen usando --charge tal cual — arranca "
            "un run nuevo con el --charge correcto sobre la estructura "
            "elegida para continuar el pipeline."
        ),
    )
    _protonation_group.add_argument(
        "--autodeprotonate", action="store_true",
        help=(
            "Análogo a --autoprotonate pero con el modo --deprotonate de "
            "CREST (sitios de desprotonación candidatos), misma mecánica "
            "(--T/--chrg/--uhf 0/--alpb water sin cambios, salida "
            "normalizada a 01_crest/crest_conformers.xyz). ATENCIÓN: carga "
            "REAL resultante = --charge - 1; tampoco se propaga a las "
            "etapas posteriores de este run."
        ),
    )
    parser.add_argument("--start",  choices=STEPS, default="crest",
                        help="etapa desde la que iniciar (default: crest). "
                             "SP (fase gas) es obligatorio: el pipeline no se "
                             "da por finalizado sin completarlo. "
                             "ADC2 es opcional en el pipeline normal (usa --adc2 "
                             "para incluirla); --start adc2 es especial: solo "
                             "requiere 04_sp/, ignora todas las demás etapas y "
                             "sus directorios.")
    parser.add_argument("--status", action="store_true",
                        help="mostrar estado de todos los runs y salir")
    parser.add_argument(
        "--censo-standalone", type=Path, default=None, metavar="ENSEMBLE_XYZ",
        help=(
            "Corre SOLO CENSO (--prescreening --screening) sobre un ensemble "
            "multi-XYZ ya existente (p.ej. crest_conformers.xyz), con el "
            "MISMO aislamiento de cores (taskset + OMP_NUM_THREADS/MKL/"
            "OPENBLAS + lockfile anti-relanzamiento) que usa el pipeline "
            "normal para CENSO — evita el problema de invocar `censo` a "
            "mano sin esas variables, que hace que sature mal los cores y "
            "arriesgue pisar una corrida anterior sobre el mismo directorio. "
            "Ignora --mol/--start/--multitemp/etc. Usa --charge (obligatorio "
            "si != 0), --temperature y --cores (total de cores para CENSO; "
            "default: mismo NPROCS_CREST del pipeline). Salida en "
            "<ENSEMBLE_XYZ>_censo/ junto al xyz de entrada, salvo --out."
        ),
    )
    parser.add_argument(
        "--out", type=Path, default=None, metavar="DIR",
        help="(solo con --censo-standalone) directorio de salida.",
    )
    parser.add_argument(
        "--omp-min", type=int, default=CENSO_OMP_MIN, metavar="N",
        help=(
            "(solo con --censo-standalone) hilos OpenMP mínimos por "
            f"sub-job de CENSO — se pasa tal cual a `censo --omp-min` Y se "
            f"usa para calcular OMP_NUM_THREADS del entorno (nunca el pool "
            f"completo de --cores; ver CENSO_OMP_MIN en el código para el "
            f"porqué). Default: {CENSO_OMP_MIN}, igual que el default de CENSO."
        ),
    )
    parser.add_argument(
        "--skipnumfreq", action="store_true",
        help=(
            "Salta NUMFREQ y pasa directo de OPT a SP/TDDFT. "
            "Usa la geometría de 02_opt/ "
            "directamente como input de SP, sin correcciones térmicas de NumFreq."
        ),
    )
    parser.add_argument(
        "--adc2", action="store_true",
        help=(
            "Incluye la etapa ADC2 (07_adc2/) al final del pipeline normal. "
            "OPCIONAL: sin este flag, el pipeline se da por FINALIZADO al "
            "completar SP (04_sp/, fase gas), que es la etapa obligatoria "
            "y requisito mínimo de finalización. TDDFT corre igualmente "
            "salvo que --start sea posterior a esa etapa; ADC2 en cambio "
            "requiere pedirse explícitamente con --adc2 (o de forma "
            "standalone con --start adc2, que sigue funcionando igual)."
        ),
    )
    parser.add_argument(
        "--spgas", action="store_true",
        help=(
            "Incluye la etapa SP fase-gas (04b_sp_gas/) — OPCIONAL, solo se "
            "calcula si se pide explícitamente. Idéntico al SP con CPCM(water) "
            "pero sin CPCM (misma geometría, mismo charge/mult, leídos de "
            "04_sp/<conf_id>/opt_opt.xyz y sp.inp). "
            "Si --mol apunta a un .xyz: se agrega como paso extra al final de "
            "esa corrida. Si --mol apunta a una carpeta de grupo (ej. MGLYVAL, "
            "sin molecules/): calcula todos los sp_gas faltantes de TODOS los "
            "runs/ dentro de esa carpeta, sin correr ningún otro paso."
        ),
    )
    parser.add_argument(
        "--spsmd", action="store_true",
        help=(
            "Incluye la etapa SP con SMD(water) (04c_sp_smd/) — OPCIONAL, "
            "solo se calcula si se pide explícitamente. Idéntico al SP con "
            "CPCM(water) pero usando SMD(water) en vez de CPCM(water) (misma "
            "geometría, mismo charge/mult, leídos de 04_sp/<conf_id>/"
            "opt_opt.xyz y sp.inp). "
            "Si --mol apunta a un .xyz: se agrega como paso extra al final de "
            "esa corrida. Si --mol apunta a una carpeta de grupo (ej. MGLYVAL, "
            "sin molecules/): calcula todos los sp_smd faltantes de TODOS los "
            "runs/ dentro de esa carpeta, sin correr ningún otro paso."
        ),
    )

    # ── Recursos de cómputo ─────────────────────────────────────────────────
    res_group = parser.add_argument_group(
        "Recursos de cómputo",
        "Usa --cores para control sencillo (recomendado). Las opciones avanzadas\n"
        "permiten ajustar cada etapa por separado y sobreescriben --cores.",
    )
    res_group.add_argument(
        "--cores", type=int, default=None, metavar="N",
        help=(
            "Cores TOTALES disponibles para el pipeline. El script los distribuye\n"
            "automáticamente entre etapas (OPT / NumFreq / SP / TDDFT).\n"
            "  maxcore se recalcula como (RAM_total_MB × 0.85) // cores,\n"
            "  con RAM_total_MB = 28000 por defecto (ajustar en MAX_RAM_MB)."
        ),
    )
    res_group.add_argument(
        "--nprocs", type=int, default=None, metavar="N",
        help=(
            f"(avanzado) Cores por job de OPT / SP / TDDFT. "
            f"Sobreescribe la parte OPT de --cores. Default: {NPROCS}."
        ),
    )
    res_group.add_argument(
        "--maxcore", type=int, default=None, metavar="MB",
        help=(
            f"(avanzado) MB/core para OPT / SP / TDDFT. "
            f"Default: {MAXCORE}. RAM por job = --nprocs × MB."
        ),
    )
    res_group.add_argument(
        "--nprocs-numfreq", type=int, default=None, metavar="N",
        help=(
            f"(avanzado) Cores para NumFreq. "
            f"Sobreescribe la parte NumFreq de --cores. Default: {NPROCS_NUMFREQ}."
        ),
    )
    res_group.add_argument(
        "--maxcore-numfreq", type=int, default=None, metavar="MB",
        help=f"(avanzado) MB/core para NumFreq. Default: {MAXCORE_NUMFREQ}.",
    )

    # ── Modo SP + NRoots 1 ─────────────────────────────────────────────────
    parser.add_argument(
        "--nroot", action="store_true",
        help=(
            "Modo SP+NRoots 1: calcula un SP con TD-DFT mínimo (NRoots 1) "
            "sobre los mejores confórmeros según energía DFT de OPT (02_opt/). "
            "Los inputs se escriben en 06_sp_nroot/<conf_id>/sp_nroot.inp. "
            "Requiere que OPT haya completado."
        ),
    )
    parser.add_argument(
        "--nroot-top", type=int, default=5, metavar="N",
        help=(
            "Número de confórmeros (top-N por E_OPT) sobre los que correr "
            "el SP+NRoots 1 con --nroot (default: 5)."
        ),
    )

    # ── Modo multi-temperatura (CENSO) ─────────────────────────────────────
    parser.add_argument(
        "--multitemp", action="store_true",
        help=(
            "Modo multi-temperatura: corre CENSO pre-screening a múltiples "
            "temperaturas definidas por --trange T_ini T_fin N_puntos. "
            "Cada temperatura produce su propio ensemble en "
            "01b_censo_prescreen/T<temp>K/. "
            "--trange es opcional (default: 274 298 24)."
        ),
    )
    parser.add_argument(
        "--trange", type=float, nargs=3, metavar=("T_INI", "T_FIN", "N"),
        default=None,
        help=(
            "Rango de temperaturas para --multitemp. "
            "Formato: T_INI T_FIN N  (en K, N = número de puntos, incluye extremos). "
            "Default cuando --multitemp está activo: 274 298 24 "
            "(24 puntos equiespaciados de 274 K a 298 K). "
            "Si N=1, se usa solo T_INI."
        ),
    )
    parser.add_argument(
        "--temperature", type=float, default=298.15, metavar="T",
        help=(
            "Temperatura de referencia para CENSO: "
            "en modo normal (sin --multitemp) es la única temperatura usada; "
            "en modo --multitemp es la temperatura cuyo ensemble CENSO se usa "
            "para seleccionar los confórmeros que pasan a OPT. "
            "En K. Default: 298.15 (25 °C)."
        ),
    )

    parser.add_argument(
        "--boltzmann-threshold", type=float, default=5.0, metavar="PCT",
        help=(
            "Umbral de peso Boltzmann (%%) para seleccionar confórmeros que pasan "
            "a OPT. Requiere que CENSO haya generado 1_SCREENING.out con los pesos. "
            "Ejemplo: --boltzmann-threshold 5.0 selecciona confs con peso ≥ 5%% "
            "(típicamente CONF7=41.7%%, CONF290=34.8%%, CONF81=5.9%% para SB). "
            "Si no se especifica, se usa select_nconf_auto (comportamiento por defecto). "
            "Si el umbral deja < NCONF_MIN confórmeros, se completan automáticamente."
        ),
    )
    parser.add_argument(
        "--spsmd-boltzmann-threshold", type=float, default=20.0, metavar="PCT",
        help=(
            "Umbral de peso Boltzmann (%%) para filtrar qué confórmeros de "
            "04_sp/ se optimizan con --spsmd (SP_SMD): solo corre sobre los "
            "confórmeros con peso Boltzmann > PCT en 1_SCREENING.out de CENSO "
            "(por estado/run). Default: 20.0. Si no hay 1_SCREENING.out "
            "disponible o el filtro no deja ningún confórmero, se cae de "
            "vuelta a todos los confórmeros de 04_sp/ (fail-open)."
        ),
    )


    raw_argv = sys.argv[1:]

    # --censo-standalone es global y de mayor prioridad, igual que --status:
    # si aparece, se ignora todo lo demas (--mol, --start, --multitemp...) y
    # solo se corre CENSO sobre el ensemble dado, vía _run_censo_standalone_cli
    # (que sale del proceso por su cuenta con SystemExit).
    if "--censo-standalone" in raw_argv:
        cs_args = parser.parse_args(raw_argv)
        if cs_args.censo_standalone:
            _run_censo_standalone_cli(cs_args)
            return

    # --status es global y de mayor prioridad: si aparece en cualquier
    # parte del comando, se ignora todo lo demas (incluyendo --mol
    # repetidos) y solo se imprime el estado de los runs.
    if "--status" in raw_argv:
        status_args = parser.parse_args(raw_argv)
        if status_args.status:
            print_status()
            return

    preamble, mol_blocks = _split_argv_by_mol(raw_argv)

    if not mol_blocks:
        # Sin --mol: comportamiento clasico de un solo run
        # (parser.error si --mol era obligatorio y no se dio).
        args = parser.parse_args(raw_argv)
        _run_one(args, parser)
        return

    if len(mol_blocks) == 1:
        args = parser.parse_args(preamble + mol_blocks[0])
        _run_one(args, parser)
        return

    # ================================================================
    #  MODO MULTI-MOLECULA: --mol repetido -> correr consecutivamente
    # ================================================================
    # Cada bloque (--mol ... hasta el proximo --mol) se parsea como
    # preamble + bloque, de modo que las especificaciones (--charge,
    # --start, --skipnumfreq, --nroot, --multitemp, etc.) se apliquen
    # segun el orden en que se escribieron en la terminal, y los flags
    # globales (--cores, --status, etc.) puestos ANTES del primer --mol
    # se hereden en todas las moleculas salvo que se sobreescriban
    # dentro del bloque de una molecula en particular.
    n = len(mol_blocks)
    print(f"\n{'='*60}")
    print(f"MODO MULTI-MOLECULA: {n} moleculas en cola (ejecucion secuencial)")
    print(f"{'='*60}")

    resumen = []
    for i, block in enumerate(mol_blocks, start=1):
        run_argv = preamble + block
        try:
            args = parser.parse_args(run_argv)
        except SystemExit:
            print(f"[{i}/{n}] ERROR: argumentos invalidos en este bloque -- se omite.")
            resumen.append((f"run_{i}", "ERROR: argumentos invalidos"))
            continue

        etiqueta = args.mol.name if args.mol else f"run_{i}"
        carga_str = args.charge if args.charge is not None else "auto"
        print(f"\n{'#'*60}")
        print(f"[{i}/{n}] {etiqueta}  (carga={carga_str}, start={args.start})")
        print(f"{'#'*60}\n")

        try:
            _run_one(args, parser)
            resumen.append((etiqueta, "OK"))
        except SystemExit as e:
            codigo = e.code if isinstance(e.code, int) else 1
            print(f"\n[{i}/{n}] {etiqueta}: termino con error (sys.exit={codigo}). "
                  f"Continuando con la siguiente molecula...\n")
            resumen.append((etiqueta, f"ERROR (sys.exit={codigo})"))
        except Exception as e:
            print(f"\n[{i}/{n}] {etiqueta}: excepcion no controlada: {e!r}. "
                  f"Continuando con la siguiente molecula...\n")
            resumen.append((etiqueta, f"EXCEPCION: {e}"))

    print(f"\n{'='*60}")
    print("RESUMEN MULTI-MOLECULA")
    print(f"{'='*60}")
    for etiqueta, estado in resumen:
        print(f"  {etiqueta:30s} -> {estado}")
    print(f"{'='*60}\n")

if __name__ == "__main__":
    main()
