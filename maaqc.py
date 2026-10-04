#!/usr/bin/env python3
"""
maaqc.py — MAA Quantum Chemistry, herramienta única
═══════════════════════════════════════════════════════════════════════════════

Unifica los cuatro scripts del protocolo (pipeline / results_collector /
protonation_cycle / extract_tddft_cubes) en un solo ejecutable con
subcomandos. Las llamadas a ORCA, orca_plot, CREST, CENSO, GoodVibes y
Multiwfn son idénticas a las de los scripts originales, igual que los
flags de cada subcomando.

    maaqc.py run      …   etapas de cómputo (ex pipeline.py)
    maaqc.py status       estado de todos los runs (atajo de `run --status`)
    maaqc.py collect  …   análisis de resultados (ex results_collector.py)
    maaqc.py pka      …   ciclo de protonación y pKa (ex protonation_cycle.py)
    maaqc.py cubes    …   cubos HOMO/LUMO y NTO (ex extract_tddft_cubes.py)

Equivalencias directas:

    python pipeline.py --mol m.xyz --charge -1 --cores 8
        → python maaqc.py run --mol m.xyz --charge -1 --cores 8

    python results_collector.py --mol MGLYVAL --state D --export xlsx
        → python maaqc.py collect --mol MGLYVAL --state D --export xlsx

    python protonation_cycle.py --mol MGLYVAL --neutral R --anion B --dianion C
        → python maaqc.py pka --mol MGLYVAL --neutral R --anion B --dianion C

    python extract_tddft_cubes.py --root PAPER1 --grid 80
        → python maaqc.py cubes --root PAPER1 --grid 80

NOTAS DE LA UNIFICACIÓN
───────────────────────────────────────────────────────────────────────────────
Colisiones de nombres resueltas (19 detectadas por AST):

  · 12 constantes eran byte-idénticas en varios archivos
    (HARTREE_TO_KCAL, HARTREE_TO_KJ, R_KCAL, LN10, T_STD, TEMPS_K,
    TEMP_LABELS, G_HPLUS_GAS, DG_SOLV_HPLUS, DG_CONC_CORR, GV_CONDA_ENV,
    GV_QS) → definición única en el NÚCLEO COMÚN de abajo.

  · WORK_ROOT existía en tres archivos con TRES semánticas distintas
    (PAPER_ROOT / <MOL>/runs / cwd/runs) → separadas en PIPE_WORK_ROOT,
    RC_WORK_ROOT y PC_WORK_ROOT. Esta era la colisión que rompía en
    silencio al concatenar los archivos.

  · 5 funciones homónimas pero divergentes. Se conserva la versión de
    protonation_cycle.py (la más evolucionada, y la que el test de Monte
    Carlo parchea por nombre) y la de results_collector.py pasa a rc_*:
        compute_pka        → rc_compute_pka        (nombres vs dicts de estado)
        parse_thermochem   → rc_parse_thermochem   (220 vs 117 líneas)
        parse_sp_energy    → rc_parse_sp_energy
        export_results     → rc_export_results
        _resolve_gv_python → _rc_resolve_gv_python

  · main() ×4 → _main_run / _main_collect / _main_pka / _main_cubes,
    despachadas por el CLI del final.

numpy, pandas y ase.io se cargan de forma perezosa: `maaqc.py status`
ya no paga el arranque de pandas/ase, que era ~1 s de los scripts viejos.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import fcntl
import json
import logging
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time
import warnings

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field as _dc_field
from datetime import datetime as _dt_datetime
from pathlib import Path
from typing import Optional


# ══════════════════════════════════════════════════════════════════════════════
#  §  NÚCLEO COMÚN
#     Constantes compartidas e imports perezosos
# ══════════════════════════════════════════════════════════════════════════════

# ── Conversión de unidades ───────────────────────────────────────────────────
HARTREE_TO_KCAL = 627.5094740631
HARTREE_TO_KJ   = HARTREE_TO_KCAL * 4.184
R_KCAL          = 1.9872041e-3     # kcal/(mol·K) — constante de gas ideal
LN10            = 2.302585093

# ── Temperaturas de análisis ─────────────────────────────────────────────────
#   250.0 K  — mínimo oceánico Antártico WAP (invierno)
#   274.15 K — temperatura media superficie océano Antártico (objetivo FAE)
#   298.15 K — condición estándar IUPAC
#   310.0 K  — temperatura fisiológica humana
TEMPS_K     = [250.0, 274.15, 298.15, 310.0]
TEMP_LABELS = ["250K", "274K", "298K", "310K"]
T_STD       = 298.15               # K — temperatura de referencia para pKa

# ── Constantes del ciclo de pKa (todas a 298.15 K) ───────────────────────────
G_HPLUS_GAS   = -6.28      # kcal/mol  G(H⁺, gas)        Tissandier et al. 1998
DG_SOLV_HPLUS = -265.9     # kcal/mol  ΔG_solv(H⁺)       Kelly et al. 2006
DG_CONC_CORR  = 1.894      # kcal/mol  δG_conc 1 atm→1 M  RT·ln(24.5)

# ── GoodVibes ────────────────────────────────────────────────────────────────
GV_CONDA_ENV = "sciencebitch"
GV_QS        = "grimme"    # entropía cuasi-armónica: 'grimme' o 'truhlar'


class _LazyModule:
    """Proxy que importa el módulo real en el primer acceso a un atributo.

    Mantiene `np`, `pd` y matplotlib fuera del arranque: un `maaqc.py status`
    o un `--help` no tienen por qué cargar el stack científico completo.
    Es transparente para isinstance(), anotaciones (hay `from __future__
    import annotations`) y acceso normal a atributos.
    """

    __slots__ = ("_name", "_mod")

    def __init__(self, name: str):
        object.__setattr__(self, "_name", name)
        object.__setattr__(self, "_mod", None)

    def _load(self):
        mod = object.__getattribute__(self, "_mod")
        if mod is None:
            import importlib
            mod = importlib.import_module(object.__getattribute__(self, "_name"))
            object.__setattr__(self, "_mod", mod)
        return mod

    def __getattr__(self, item):
        return getattr(self._load(), item)

    def __dir__(self):
        return dir(self._load())

    def __repr__(self):
        name = object.__getattribute__(self, "_name")
        state = "cargado" if object.__getattribute__(self, "_mod") else "perezoso"
        return f"<_LazyModule {name!r} ({state})>"


np = _LazyModule("numpy")
pd = _LazyModule("pandas")


def read(*args, **kwargs):
    """ase.io.read diferido (lo usa la etapa CREST/OPT)."""
    from ase.io import read as _read
    return _read(*args, **kwargs)


def write(*args, **kwargs):
    """ase.io.write diferido."""
    from ase.io import write as _write
    return _write(*args, **kwargs)


# ══════════════════════════════════════════════════════════════════════════════
#  §  PIPELINE
#     Etapas de cómputo: CREST → CENSO → OPT → NumFreq → SP → TDDFT → ADC2
# ══════════════════════════════════════════════════════════════════════════════

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
PIPE_WORK_ROOT  = PAPER_ROOT                        # alias, usado solo por print_status()
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
      - Los dígitos finales se ignoran ("AB2" → 'B', "GB2" → 'B').
      - Si el nombre no termina en una letra reconocida en
        STATE_CHARGE_MAP (p.ej. es algo como
        "MGLYVAL" sin sufijo de estado), devuelve None: no se puede
        inferir con confianza y el caller debe pedir --charge explícito.
      - No es case-sensitive: "b" y "B" son equivalentes.

    Devuelve la carga inferida (int) o None si no se pudo inferir.
    """
    if not mol_name:
        return None
    # Tautómeros / variantes numeradas (AB2, GB2, ...): se ignoran los dígitos
    # finales y se mira la última letra. Antes 'AB2' devolvía None y caía al
    # default 0 (¡carga equivocada para un anión!).
    stem = mol_name.rstrip("0123456789")
    if not stem:
        return None
    last_char = stem[-1].upper()
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

R_GAS_KCAL          = 1.987204e-3      # cte. de los gases, kcal/(mol·K)


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


# ─────────────────────────────────────────────
#  BOLTZMANN DESDE ENERGÍAS SP PRE-EXISTENTES
#  (nivel DFT del propio SP, no xTB/DFT de prescreening de CENSO — usado
#  como criterio principal de selección para --spsmd, ver stage_sp_smd)
# ─────────────────────────────────────────────

def compute_boltzmann_from_energies(
    energies:    dict[str, float],
    temperature: float = 298.15,
) -> dict[str, float]:
    """
    Calcula pesos de Boltzmann (%) a partir de energías electrónicas
    (Hartree) ya calculadas — típicamente FINAL SINGLE POINT ENERGY de
    04_sp/conf_*/sp.out — en vez de reusar el peso Boltzmann de CENSO
    (que viene de xTB/DFT de prescreening, un nivel de teoría distinto).

    w_i = exp(-ΔE_i / R·T) / Σ_j exp(-ΔE_j / R·T),  ΔE_i = (E_i − E_min) en kcal/mol

    Args:
        energies    : {conf_id → energía electrónica en Hartree}.
        temperature : temperatura en K (default 298.15; usa la misma que
                      --temperature para ser consistente con el resto del
                      pipeline, p.ej. 274.15 K en runs Antárticos).

    Devuelve {conf_id → peso%}, ordenable por valor descendente. Dict
    vacío si `energies` está vacío.
    """
    if not energies:
        return {}

    e_min = min(energies.values())
    boltz_factors = {
        cid: math.exp(-((e - e_min) * HARTREE_TO_KCAL) / (R_GAS_KCAL * temperature))
        for cid, e in energies.items()
    }
    total = sum(boltz_factors.values())
    if total <= 0:
        return {cid: 0.0 for cid in energies}
    return {cid: 100.0 * f / total for cid, f in boltz_factors.items()}


def extract_sp_boltzmann_weights(
    sp_dirs:     list[Path],
    log:         logging.Logger,
    temperature: float = 298.15,
    stage_tag:   str   = "SP_SMD",
) -> dict[str, float]:
    """
    Extrae los pesos de Boltzmann de un conjunto de confórmeros (p.ej. todos
    los que llegaron a 04_sp/ con SP exitoso) usando sus energías de SP
    PRE-EXISTENTES (FINAL SINGLE POINT ENERGY de <conf_dir>/sp.out) en vez
    del peso Boltzmann de 1_SCREENING.out de CENSO.

    Esto da el peso Boltzmann al nivel DFT real del SP (el mismo con el
    que luego se pesan espectros/observables) y sirve como criterio de
    selección más fiable que el de CENSO, cuyo nivel de teoría (xTB/DFT de
    prescreening) puede subestimar la brecha energética entre confórmeros
    y dejar pasar de más.

    Fail-open por confórmero: si sp.out falta, no terminó bien, o no tiene
    energía parseable, ese confórmero se excluye del cálculo (con warning)
    en vez de abortar toda la etapa.

    Devuelve {conf_id → peso%} (vacío si ningún confórmero tiene energía
    SP válida).
    """
    energies: dict[str, float] = {}
    for d in sp_dirs:
        out_file = d / "sp.out"
        if not check_success(out_file):
            log.warning(f"[{stage_tag}-BOLTZMANN] {d.name}: sp.out ausente o "
                        f"incompleto — excluido del cálculo Boltzmann")
            continue
        e = extract_sp_energy(out_file)
        if e is None:
            log.warning(f"[{stage_tag}-BOLTZMANN] {d.name}: sp.out sin "
                        f"energía parseable — excluido del cálculo Boltzmann")
            continue
        energies[d.name] = e

    if not energies:
        log.warning(f"[{stage_tag}-BOLTZMANN] ningún confórmero con energía "
                    f"SP válida — no se pudo calcular Boltzmann desde SP")
        return {}

    weights = compute_boltzmann_from_energies(energies, temperature=temperature)

    ranked = sorted(weights.items(), key=lambda x: -x[1])
    ranked_str = ", ".join(f"{cid}({w:.1f}%)" for cid, w in ranked)
    log.info(f"[{stage_tag}-BOLTZMANN] pesos Boltzmann desde energías SP "
              f"preexistentes (T={temperature:.2f} K), "
              f"{len(weights)} conf(s): {ranked_str}")

    return weights


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


def _filter_dirs_by_sp_weights(
    conf_dirs: list[Path],
    weights:   dict[str, float],
    threshold: float,
    log:       logging.Logger,
    stage_tag: str,
) -> list[Path]:
    """
    Filtra conf_dirs dejando solo los que superan `threshold`% de peso
    Boltzmann calculado desde energías de SP PRE-EXISTENTES (ver
    extract_sp_boltzmann_weights) — a diferencia de _filter_dirs_by_boltzmann,
    que usa el peso de 1_SCREENING.out de CENSO (nivel xTB/DFT de
    prescreening, distinto y potencialmente menos fiable que el propio SP).

    Fail-open: si `weights` está vacío o ningún confórmero supera el
    umbral, devuelve conf_dirs sin filtrar con un warning.
    """
    if not weights:
        log.warning(f"[{stage_tag}] sin pesos Boltzmann de SP disponibles — "
                    f"no se puede filtrar, se usan los {len(conf_dirs)} sin filtrar")
        return conf_dirs

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
        log.warning(f"[{stage_tag}] ningún confórmero supera >{threshold}% "
                    f"Boltzmann (SP) — cayendo de vuelta a los {len(conf_dirs)} sin filtrar")
        return conf_dirs

    kept_str = ", ".join(f"{d.name}({weights[d.name]:.1f}%)" for d in kept)
    log.info(f"[{stage_tag}] filtro Boltzmann (SP, energías preexistentes) "
              f">{threshold}%: {len(kept)}/{len(conf_dirs)} seleccionados: {kept_str}")
    if dropped:
        log.info(f"[{stage_tag}] {len(dropped)} conf(s) descartados por debajo "
                  f"del umbral: {[f'{cid}({w:.1f}%)' for cid, w in dropped]}")
    if no_wt:
        log.warning(f"[{stage_tag}] {len(no_wt)} conf(s) sin peso SP calculado, "
                    f"descartados: {no_wt}")
    return kept


def _collect_spsmd_selection(
    sp_dirs:           list[Path],
    mol_dir:           Path,
    log:               logging.Logger,
    boltzmann_min_pct: float = 20.0,
    temperature:       float = 298.15,
) -> list[Path]:
    """
    Selecciona, para un mol_dir/estado dado, qué confórmeros de sp_dirs
    (04_sp/conf_*/ ya exitosos) deben pasar a SP_SMD, y persiste el
    resultado en 04c_sp_smd/sp_boltzmann_weights.json. NO ejecuta nada —
    solo decide (ver _run_spsmd_execution para la ejecución).

    El peso se calcula PRIMERO desde las energías de SP PRE-EXISTENTES
    (FINAL SINGLE POINT ENERGY de 04_sp/<conf_id>/sp.out, ver
    extract_sp_boltzmann_weights y _filter_dirs_by_sp_weights) — al nivel
    DFT real del propio SP, no al nivel xTB/DFT de prescreening de CENSO,
    que puede subestimar la brecha energética entre confórmeros y dejar
    pasar de más (fail-open) cuando ninguno supera el umbral a ese nivel.

    Solo si NINGÚN confórmero tiene energía de SP parseable (caso límite,
    no debería pasar con sp_dirs ya exitosos) se cae de vuelta al filtro
    Boltzmann de CENSO (1_SCREENING.out, ver _filter_dirs_by_boltzmann) como
    respaldo; y si tampoco hay dato ahí, se corre sobre sp_dirs sin filtrar.

    Devuelve la lista de directorios 04_sp/<conf_id>/ seleccionados.
    """
    sp_boltzmann_weights = extract_sp_boltzmann_weights(
        sp_dirs, log, temperature=temperature, stage_tag="SP_SMD"
    )

    if sp_boltzmann_weights:
        selected_dirs = _filter_dirs_by_sp_weights(
            sp_dirs, sp_boltzmann_weights, threshold=boltzmann_min_pct,
            log=log, stage_tag="SP_SMD",
        )
    else:
        log.warning("[SP_SMD] sin energías de SP válidas para Boltzmann — "
                    "cayendo al filtro Boltzmann de CENSO (1_SCREENING.out)")
        selected_dirs = _filter_dirs_by_boltzmann(
            sp_dirs, mol_dir, log, threshold=boltzmann_min_pct, stage_tag="SP_SMD"
        )

    if sp_boltzmann_weights:
        smd_root = mol_dir / "04c_sp_smd"
        smd_root.mkdir(parents=True, exist_ok=True)
        weights_path = smd_root / "sp_boltzmann_weights.json"
        try:
            weights_path.write_text(json.dumps(
                {
                    "temperature_K":     temperature,
                    "source":            "04_sp/<conf_id>/sp.out (FINAL SINGLE POINT ENERGY)",
                    "threshold_pct":     boltzmann_min_pct,
                    "weights_pct":       sp_boltzmann_weights,
                    "selected_conf_ids": [d.name for d in selected_dirs],
                },
                indent=2,
            ))
        except OSError as e:
            log.warning(f"[SP_SMD-BOLTZMANN] no se pudo guardar {weights_path.name}: {e}")

    return selected_dirs


def _run_spsmd_execution(
    tasks_all:       list[dict],
    log:             logging.Logger,
    batch_size:      int = 0,
    nprocs_per_job:  Optional[int] = None,
    maxcore_per_job: Optional[int] = None,
) -> list[Path]:
    """
    Ejecuta una lista de espera (queue) de tareas SP_SMD ya SELECCIONADAS
    (ver _collect_spsmd_selection), potencialmente de varios mol_dir /
    estados de protonación distintos a la vez, en lotes de ≤batch_size
    jobs en paralelo (batch_size=0 → SP_BATCH_SIZE global, típicamente
    3 jobs × 4 cores = 12 cores).

    Armar la queue completa ANTES de llamar a esta función (ver
    stage_spsmd_batch) es lo que evita que un estado con un único
    confórmero seleccionado dispare un lote de 1 job en solitario
    (desperdiciando 8 de 12 cores): se junta con lo seleccionado de los
    demás estados y se lanza todo en el mismo lote de 3.

    Cada dict de tasks_all debe tener: conf_id, sp_dir, template, nprocs,
    maxcore, mol_dir, log. Se le agrega cpu_list por lote.

    Devuelve la lista de directorios 04c_sp_smd/<conf_id>/ exitosos.
    """
    if not tasks_all:
        log.info("[SP_SMD] cola vacía — nada que ejecutar")
        return []

    _batch = batch_size if batch_size > 0 else SP_BATCH_SIZE

    n_lotes = -(-len(tasks_all) // _batch)   # ceiling division
    _nprocs_first = tasks_all[0]["nprocs"]
    log.info(f"[SP_SMD] {len(tasks_all)} SPs SMD en cola, {n_lotes} lote(s) de "
             f"≤{_batch} (CPUs activos/lote: {_batch * _nprocs_first})")

    successful: list[Path] = []
    _check_cpu_budget(_batch, _nprocs_first, log)

    idx      = 0
    lote_num = 0
    while idx < len(tasks_all):
        want      = min(_batch, len(tasks_all) - idx)
        nprocs_want = tasks_all[idx]["nprocs"]
        cpu_lists = _acquire_cpu_lists(want, nprocs_want, log, "SP_SMD")
        if not cpu_lists:
            log.error("[SP_SMD] no hay suficientes cores físicos independientes "
                       "para asignar ni un solo job — abortando etapa SP_SMD")
            break
        batch    = tasks_all[idx: idx + len(cpu_lists)]
        lote_num += 1
        labels   = [f"{t['mol_dir'].name}/{t['conf_id']}" for t in batch]

        for t, cl in zip(batch, cpu_lists):
            t["cpu_list"] = cl

        log.info(f"[SP_SMD] ── Lote {lote_num}/~{n_lotes}: {labels}  "
                 f"(cores: {[','.join(map(str, cl)) for cl in cpu_lists]})")

        with ThreadPoolExecutor(max_workers=len(batch)) as pool:
            futures = {
                pool.submit(_sp_smd_worker, t): f"{t['mol_dir'].name}/{t['conf_id']}"
                for t in batch
            }
            try:
                for fut in as_completed(futures):
                    res   = fut.result()
                    label = futures[fut]
                    tag   = "✓" if res["success"] else "✗"
                    log.info(f"[SP_SMD]   {tag} {label}  ({res['elapsed']:.1f} min)")
                    if res["success"] and res["sp_smd_dir"] is not None:
                        successful.append(res["sp_smd_dir"])
                    else:
                        log.warning(f"[SP_SMD] {label} falló — excluido")
            except KeyboardInterrupt:
                log.warning("[SP_SMD] Ctrl-C — cancelando lote, esperando workers…")
                pool.shutdown(wait=False, cancel_futures=True)
                _kill_active_procs()
                _cpu_lock_release()
                sys.exit(130)

        _cpu_lock_release()
        idx += len(cpu_lists)

    log.info(f"[SP_SMD] {len(successful)}/{len(tasks_all)} completados")
    return successful


def stage_sp_smd(sp_dirs: list[Path], mol_dir: Path, log: logging.Logger,
                  batch_size: int = 0,
                  nprocs_per_job: Optional[int] = None,
                  maxcore_per_job: Optional[int] = None,
                  boltzmann_min_pct: float = 20.0,
                  temperature: float = 298.15) -> list[Path]:
    """
    Etapa 04c_sp_smd (OPCIONAL — requiere --spsmd): single point con
    SMD(water) (idéntico al SP con CPCM(water) pero usando SMD en vez de
    CPCM), sobre los mismos confórmeros y con la misma geometría/charge/mult
    que un 04_sp/<conf_id>/ ya exitoso, para UN SOLO mol_dir/estado.

    Selecciona primero qué confórmeros corren (ver
    _collect_spsmd_selection: filtro Boltzmann de SP PRE-EXISTENTE,
    respaldo en CENSO) y luego los ejecuta en lotes (ver
    _run_spsmd_execution). Para el modo batch multi-estado que arma UNA
    sola cola con la selección de TODOS los estados antes de ejecutar
    (evitando lotes de 1 job), ver stage_spsmd_batch.

    Devuelve directorios 04c_sp_smd/<conf_id>/ exitosos.

    batch_size=0  → usa SP_BATCH_SIZE global (comparte recursos con SP normal).
    nprocs_per_job / maxcore_per_job sobreescriben los globals NPROCS / MAXCORE
    para cada worker paralelo; si son None se usan los globals vigentes.
    temperature   : T (K) para el cálculo Boltzmann de SP (default 298.15;
                    pasar --temperature para runs Antárticos a 274.15 K).
    """
    selected_dirs = _collect_spsmd_selection(
        sp_dirs, mol_dir, log,
        boltzmann_min_pct=boltzmann_min_pct, temperature=temperature,
    )

    _nprocs  = nprocs_per_job  if nprocs_per_job  is not None else NPROCS
    _maxcore = maxcore_per_job if maxcore_per_job is not None else MAXCORE
    template = get_template("sp_smd")

    log.info(f"[SP_SMD] single point SMD(water) para {len(selected_dirs)} "
             f"confórmero(s) de {mol_dir.name}")

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
        for d in selected_dirs
    ]

    return _run_spsmd_execution(tasks_all, log, batch_size=batch_size,
                                 nprocs_per_job=_nprocs, maxcore_per_job=_maxcore)


def stage_spsmd_batch(group_dir: Path, log: logging.Logger,
                       boltzmann_min_pct: float = 20.0,
                       temperature: float = 298.15,
                       batch_size: int = 0,
                       nprocs_per_job: Optional[int] = None,
                       maxcore_per_job: Optional[int] = None) -> None:
    """
    Modo batch de --spsmd: recorre TODAS las corridas (estados de
    protonación) dentro de una carpeta de grupo (<grupo>/runs/*/04_sp/conf_*).

    Funciona en DOS FASES para no desperdiciar capacidad de cómputo:

    1) SELECCIÓN: para cada mol_dir/estado, decide qué confórmeros pasan a
       SP_SMD (ver _collect_spsmd_selection — filtro Boltzmann de SP
       PRE-EXISTENTE, respaldo en CENSO) y los junta TODOS en una única
       lista de espera (queue), sin importar de qué estado vienen.
    2) EJECUCIÓN: recién con la queue completa arma lotes de ≤batch_size
       jobs en paralelo (default SP_BATCH_SIZE, típicamente 3 jobs × 4
       cores = 12 cores) y los lanza (ver _run_spsmd_execution).

    Así, un estado con un único confórmero seleccionado (p.ej. porque solo
    uno supera el umbral Boltzmann) NO dispara un lote de 1 job en
    solitario: se sigue con el siguiente estado, se suma su selección a la
    cola, y todos se lanzan juntos en el mismo lote de 3.

    Como stage_sp_smd/_run_spsmd_execution usan check_success() para
    saltar los ya hechos, correr el mismo comando dos veces es seguro
    (solo procesa los faltantes).

    Se invoca cuando --spsmd se pasa con --mol apuntando a una carpeta de
    grupo (ej. MGLYVAL) en vez de a un .xyz específico.
    """
    runs_dir = group_dir / "runs"
    if not runs_dir.is_dir():
        sys.exit(f"ERROR: no se encontró {runs_dir} — "
                  f"¿'{group_dir.name}' es una carpeta de grupo válida?")

    _nprocs  = nprocs_per_job  if nprocs_per_job  is not None else NPROCS
    _maxcore = maxcore_per_job if maxcore_per_job is not None else MAXCORE
    template = get_template("sp_smd")

    # ── Fase 1: selección por estado — arma la cola completa ─────────────
    tasks_all: list[dict] = []
    n_estados = 0
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

        n_estados += 1
        log.info(f"[SP_SMD] {mol_dir.name}: {len(sp_dirs)} confórmeros con SP CPCM completo")
        selected_dirs = _collect_spsmd_selection(
            sp_dirs, mol_dir, log,
            boltzmann_min_pct=boltzmann_min_pct, temperature=temperature,
        )
        tasks_all.extend({
            "conf_id":  d.name,
            "sp_dir":   d,
            "template": template,
            "nprocs":   _nprocs,
            "maxcore":  _maxcore,
            "mol_dir":  mol_dir,
            "log":      log,
        } for d in selected_dirs)

    log.info(f"[SP_SMD] selección completa: {len(tasks_all)} confórmero(s) "
             f"en cola de {n_estados} estado(s) — lanzando en batch")

    # ── Fase 2: ejecuta la cola completa en lotes ─────────────────────────
    successful = _run_spsmd_execution(
        tasks_all, log, batch_size=batch_size,
        nprocs_per_job=_nprocs, maxcore_per_job=_maxcore,
    )

    total_ok   = len(successful)
    total_pend = len(tasks_all) - total_ok
    log.info(f"[SP_SMD] batch finalizado: {total_ok} completados, {total_pend} pendientes/fallidos")

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
    if not PIPE_WORK_ROOT.exists():
        print(f"No se encontró el directorio '{PIPE_WORK_ROOT}/'. Aún no hay cálculos.")
        return

    print(f"\n{'─'*60}")
    print(f"  Estado de runs en: {PIPE_WORK_ROOT.resolve()}")
    print(f"  (raíz auto-detectada — usa $PIPELINE_ROOT para forzar otra)")
    print(f"{'─'*60}")

    # Estructura: PAPER_ROOT/<grupo>/runs/<mol_name>/<etapas>
    # (p.ej. PAPER1/MGLYVAL/runs/D/01_crest/, .../D/02_opt/, etc.)
    for group_dir in sorted(PIPE_WORK_ROOT.iterdir()):
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
        stage_spsmd_batch(group_dir, log,
                           boltzmann_min_pct=args.spsmd_boltzmann_threshold,
                           temperature=args.temperature)
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
                                   boltzmann_min_pct=args.spsmd_boltzmann_threshold,
                                   temperature=args.temperature)
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


def _main_run():
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


# ══════════════════════════════════════════════════════════════════════════════
#  §  COLLECTOR
#     Recolección y análisis: Boltzmann multi-T, GoodVibes, ESP/Multiwfn, ADC2 vs TD-DFT
# ══════════════════════════════════════════════════════════════════════════════

try:
    from ase.io import read as ase_read
    HAS_ASE = True
except ImportError:
    HAS_ASE = False

# ─────────────────────────────────────────────
#  CONSTANTES
# ─────────────────────────────────────────────

# Raíz del proyecto donde viven las carpetas de cada molécula, p.ej.:
#   PROJECT_ROOT/MGLYVAL/runs/<state>/04_sp/conf_xxx/...
# Se puede sobrescribir con --project-root; por defecto usa la ruta real
# del sistema de kmilozeta.
PROJECT_ROOT    = Path("/home/kmilozeta/PAPER1")

# WORK_ROOT se recalcula en main() a partir de --mol como
#   PROJECT_ROOT / <MOL> / "runs"
# y cada "state" (antes llamado mol_name en todo el código) es un
# subdirectorio de WORK_ROOT: WORK_ROOT / <state> / 04_sp / conf_xxx / ...
RC_WORK_ROOT       = Path("runs")
HARTREE_TO_EV   = 27.211396132
KB_KCAL         = 0.0019872041     # kcal/(mol·K)

# Temperaturas de análisis:
#   250.0 K  — mínimo oceánico Antártico WAP (invierno)
#   274.15 K — temperatura media superficie oceáno Antártico (objetivo FAE)
#   298.15 K — condición estándar IUPAC
#   310.0 K  — temperatura fisiológica humana

# ── Constantes del ciclo de pKa (todas a 298.15 K) ──────────────────────────
# G(H⁺, gas) = -6.28 kcal/mol  [Tissandier et al., J. Phys. Chem. A 1998]
# ΔG_solv(H⁺) = -265.9 kcal/mol [Kelly, Cramer & Truhlar, J. Phys. Chem. B 2006]
# Corrección de estado estándar de fase gas → solución (1 atm → 1 M):
#   δG_conc = RT·ln(RT/P°V°) = RT·ln(24.5) ≈ +1.894 kcal/mol a 298.15 K

# Directorio opcional para SPs en fase gas (para el ciclo de pKa)
# Estructura esperada: runs/<mol>/06_sp_gas/<conf_id>/sp_gas.out
PKA_GAS_SP_SUBDIR = "06_sp_gas"

# ── Constantes Multiwfn / ESP de superficie ─────────────────────────────────
# Ejecutables (rutas completas — mismo criterio que ORCA: nunca confiar en PATH)
ORCA_2MKL_BIN   = "/home/kmilozeta/orca_6_1_1_avx2/orca_2mkl"
MULTIWFN_BIN    = "/home/kmilozeta/multiwfn/Multiwfn_2026.6.2_bin_Linux/Multiwfn"

# Subcarpeta del pipeline donde vive la función de onda de cada confórmero
# (el .gbw usado para el análisis ESP es el de la etapa TDDFT: runs/<mol>/05_tddft/conf_XXX/)
ESP_SOURCE_SUBDIR = "05_tddft"

# Subcarpeta creada dentro de cada conf_xxx para no ensuciar 05_tddft/
ESP_SUBDIR        = "esp_multiwfn"

# Timeouts generosos: orca_2mkl es rápido, pero el análisis de superficie ESP
# (mapeo + búsqueda de extremos sobre miles de puntos) puede tardar varios
# minutos por confórmero según el tamaño de la molécula y la densidad de malla.
ORCA_2MKL_TIMEOUT = 120     # s
MULTIWFN_TIMEOUT  = 21600   # s — 6 h por confórmero (ESP puede demorar bastante)

# Multiwfn respeta OMP_NUM_THREADS para paralelizar el análisis de superficie
# ESP (mapeo de malla + búsqueda de extremos). Por defecto hereda el entorno
# del proceso padre, que puede estar limitado a 2. Ajustar según núcleos
# disponibles en la máquina (ver también compquim-kmilozeta: MKL/OpenBLAS a 1
# para evitar sobre-suscripción cuando Multiwfn ya reparte con OMP).
MULTIWFN_NTHREADS = 8


# ─────────────────────────────────────────────
#  PARSERS — individuales y robustos
# ─────────────────────────────────────────────

def rc_parse_sp_energy(out_file: Path) -> Optional[float]:
    """Extrae la última energía FINAL SINGLE POINT ENERGY del .out de SP."""
    if not out_file.exists():
        return None
    energy = None
    for line in out_file.read_text(errors="replace").splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                energy = float(line.split()[-1])
            except ValueError:
                pass
    return energy


def parse_opt_energy(out_file: Path) -> Optional[float]:
    """
    Extrae la energía del último ciclo de optimización (paso final convergido).
    ORCA reporta 'FINAL SINGLE POINT ENERGY' también durante Opt, así que
    tomamos la última ocurrencia, que corresponde al mínimo convergido.
    """
    return rc_parse_sp_energy(out_file)


def parse_imaginary_freqs(freq_out: Path) -> list[float]:
    """Devuelve lista de frecuencias imaginarias (cm⁻¹, valores negativos)."""
    imag = []
    if not freq_out.exists():
        return imag
    for line in freq_out.read_text(errors="replace").splitlines():
        m = re.match(r"\s+\d+:\s+(-\d+\.\d+)\s+cm", line)
        if m:
            imag.append(float(m.group(1)))
    return imag


def parse_tddft_excitations(tddft_out: Path) -> list[dict]:
    """
    Parsea excitaciones del bloque 'ABSORPTION SPECTRUM VIA TRANSITION
    ELECTRIC DIPOLE MOMENTS' del output TDDFT de ORCA 6.
    Formato: state_label -> final_label   E_eV  E_cm-1  wl_nm  fosc ...
    """
    excitations = []
    if not tddft_out.exists():
        return excitations

    pattern = re.compile(
        r"^\s*\S+\s*->\s*(\d+)-\S+\s+"
        r"([\d.]+)\s+"        # E_eV
        r"[\d.]+\s+"          # E_cm-1 (ignorado)
        r"([\d.]+)\s+"        # wavelength_nm
        r"([\d.]+)"           # fosc
    )
    in_block = False
    for line in tddft_out.read_text(errors="replace").splitlines():
        if "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS" in line:
            in_block = True
            continue
        if in_block:
            if "ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS" in line:
                break
            m = pattern.match(line)
            if m:
                excitations.append({
                    "state":         int(m.group(1)),
                    "E_eV":          float(m.group(2)),
                    "wavelength_nm": float(m.group(3)),
                    "fosc":          float(m.group(4)),
                })
    return excitations


def parse_crest_energies(conformers_xyz: Path) -> list[Optional[float]]:
    """
    Lee las energías GFN2-xTB embebidas por CREST en el multi-XYZ.
    CREST escribe la energía en Hartree en la línea de comentario (línea 2
    de cada frame), típicamente como un número flotante negativo solo o
    precedido de 'Etot='.
    Devuelve lista con la energía de cada frame (None si no parseable).
    """
    if not conformers_xyz.exists():
        return []

    energies: list[Optional[float]] = []

    # Intentar primero con ASE (más robusto para multi-XYZ)
    if HAS_ASE:
        try:
            frames = ase_read(str(conformers_xyz), index=":")
            if not isinstance(frames, list):
                frames = [frames]
            for atoms in frames:
                comment = atoms.info.get("comment", "")
                e = _parse_energy_from_comment(comment)
                energies.append(e)
            return energies
        except Exception:
            pass  # fallback al parser manual

    # Fallback: parser manual línea a línea
    text = conformers_xyz.read_text(errors="replace")
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        try:
            natoms = int(line)
        except ValueError:
            i += 1
            continue
        if i + 1 < len(lines):
            comment = lines[i + 1]
            e = _parse_energy_from_comment(comment)
        else:
            e = None
        energies.append(e)
        i += natoms + 2  # saltar el bloque completo
    return energies


def _parse_energy_from_comment(comment: str) -> Optional[float]:
    """Extrae energía de la línea de comentario de un frame XYZ de CREST."""
    comment = comment.strip()
    # Patrón "Etot= -123.456789"
    m = re.search(r"Etot[=\s]+([-\d.]+)", comment)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    # Número flotante negativo solo (formato CREST clásico)
    m = re.match(r"^\s*([-\d.]+)\s*$", comment)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    # xTB escribe a veces: "-123.456789 Eh"
    m = re.search(r"([-\d.]+)\s*Eh", comment)
    if m:
        try:
            return float(m.group(1))
        except ValueError:
            pass
    return None


# ─────────────────────────────────────────────
#  PARSER TERMOQUÍMICA (NUMFREQ ORCA 6)
# ─────────────────────────────────────────────

def rc_parse_thermochem(freq_out: Path, T: float = T_STD) -> dict:
    """
    Extrae las correcciones termoquímicas del bloque THERMOCHEMISTRY de ORCA 6.

    ORCA 6 reporta un bloque completo a la temperatura especificada en el input
    (por defecto 298.15 K). Este parser extrae todas las contribuciones
    relevantes para el ciclo termodinámico compuesto:

        G_total_i = E_SP(TZVPD/CPCM) + G_corr(SVP/CPCM)

    donde G_corr = H_corr + ZPE - T·S_total se obtiene como:
        G_corr = G_total_freq - E_elec_freq
    siendo G_total_freq la 'Total Gibbs free energy' del NumFreq
    y E_elec_freq la energía electrónica del mismo NumFreq (def2-SVP).

    Retorna dict con las siguientes claves (todas en Hartree salvo indicación):
        E_elec_freq    : energía electrónica del NumFreq (referencia del ciclo)
        ZPE            : corrección de punto cero
        H_corr         : corrección entálpica total (ZPE + contribuciones térmicas)
        TS_total       : término T·S a T_parsed, tal cual lo imprime ORCA
                         ('Final entropy term'), en Eh
        S_total        : entropía pura derivada, S = TS_total / T_parsed (Eh/K)
        G_corr         : H_corr - TS_total  (corrección de Gibbs a T_parsed)
        G_total_freq   : E_elec_freq + G_corr (Gibbs total a nivel SVP)
        T_parsed       : temperatura a la que ORCA calculó la termoquímica (K)
        n_imag         : número de frecuencias imaginarias detectadas
        imag_freqs_cm1 : lista de valores imaginarios (cm⁻¹, negativos)
        ok             : True si todos los campos críticos fueron encontrados

    Si el archivo no existe o el bloque no se puede parsear, retorna un dict
    con ok=False y None en todos los campos numéricos.

    Notas sobre el formato ORCA 6:
    ─────────────────────────────
    El bloque de termoquímica en ORCA 6 (NumFreq + CPCM, con Quasi-RRHO) tiene
    esta estructura:

        THERMOCHEMISTRY AT 298.15K
        Zero point energy                   ...      0.123456 Eh
        Total thermal correction                     0.000789 Eh     (ZPE ya incluido)
        Total Enthalpy                      ...    -500.123456 Eh
        Final entropy term                  ...      0.045678 Eh   48.42 kcal/mol
        Final Gibbs free energy             ...    -500.234567 Eh

    ⚠️ IMPORTANTE: 'Final entropy term' de ORCA 6 NO es la entropía S en
    Eh/K. Es el término T·S ya evaluado a la temperatura del bloque
    (T_ref = 298.15 K por defecto), expresado en Eh. Se puede confirmar
    dividiendo el valor en kcal/mol que ORCA imprime al lado: por ejemplo
    0.07716168 Eh × 627.509 kcal/mol/Eh ≈ 48.42 kcal/mol, que es exactamente
    el valor que ORCA muestra en esa misma línea — si fuera S puro (Eh/K),
    al multiplicar por T=298.15 quedaría ~23 Eh, muy por encima del rango
    físico razonable para esta corrección.

    Por lo tanto:
        G = H - T·S  →  G_corr = H_corr - TS_ref
    donde H_corr = Total Enthalpy - E_elec, TS_ref = 'Final entropy term'
    (ya es T·S a T_ref, NO multiplicar por T de nuevo).

    Si se necesita S puro en Eh/K (p.ej. para recalcular G_corr a otra T),
    se debe derivar como S = TS_ref / T_ref.
    """
    null = {
        "E_elec_freq":    None,
        "ZPE":            None,
        "H_corr":         None,
        "TS_total":       None,
        "S_total":        None,
        "G_corr":         None,
        "G_total_freq":   None,
        "T_parsed":       None,
        "n_imag":         0,
        "imag_freqs_cm1": [],
        "ok":             False,
    }

    if not freq_out.exists():
        return null

    try:
        text = freq_out.read_text(errors="replace")
    except OSError:
        return null

    lines = text.splitlines()

    # ── 1. Frecuencias imaginarias ────────────────────────────────────────
    imag_freqs: list[float] = []
    for line in lines:
        m = re.match(r"\s+\d+:\s+(-\d+\.\d+)\s+cm", line)
        if m:
            imag_freqs.append(float(m.group(1)))

    # ── 2. Energía electrónica del NumFreq (FINAL SINGLE POINT ENERGY) ───
    # Tomamos la PRIMERA ocurrencia: corresponde al SP inicial de la freq,
    # antes de los desplazamientos numéricos.
    E_elec_freq: Optional[float] = None
    for line in lines:
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                E_elec_freq = float(line.split()[-1])
                break
            except ValueError:
                pass

    # ── 3. Bloque THERMOCHEMISTRY ─────────────────────────────────────────
    # Buscar el bloque más reciente (última ocurrencia, por si hay más de uno)
    thermo_start = -1
    for idx, line in enumerate(lines):
        if re.search(r"THERMOCHEMISTRY AT\s+[\d.]+\s*K", line, re.IGNORECASE):
            thermo_start = idx

    if thermo_start < 0:
        warnings.warn(
            f"[THERMOCHEM] No se encontró bloque THERMOCHEMISTRY en {freq_out}. "
            f"¿El cálculo terminó normalmente? ¿Es un NumFreq?"
        )
        result = dict(null)
        result.update({
            "E_elec_freq":    E_elec_freq,
            "n_imag":         len(imag_freqs),
            "imag_freqs_cm1": imag_freqs,
        })
        return result

    # Extraer temperatura del encabezado del bloque
    T_parsed: Optional[float] = None
    m = re.search(r"THERMOCHEMISTRY AT\s+([\d.]+)\s*K",
                  lines[thermo_start], re.IGNORECASE)
    if m:
        try:
            T_parsed = float(m.group(1))
        except ValueError:
            pass

    # Parsear líneas del bloque hasta el siguiente separador o fin de archivo
    ZPE:          Optional[float] = None
    H_total:      Optional[float] = None   # Total Enthalpy (Eh)  = E_elec + H_corr
    TS_total:     Optional[float] = None   # Final entropy term = T·S a T_parsed (Eh)
    G_total_freq: Optional[float] = None   # Total/Final Gibbs free energy (Eh)

    # Patrones robustos para ORCA 6 (espaciado variable, unidades variables)
    _pat_zpe    = re.compile(r"Zero point energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_h      = re.compile(r"Total [Ee]nthalpy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_s      = re.compile(r"Final entropy term\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    # ORCA 6.1.1 con Quasi-RRHO (Quasi RRHO ... True) imprime esta línea como
    # "Final Gibbs free energy", no "Total Gibbs free energy". Sin Quasi-RRHO
    # (o en otras versiones/anotaciones) puede aparecer la otra forma.
    # Aceptamos ambas para que el parser sea robusto a los dos modos.
    _pat_g      = re.compile(r"(?:Total|Final) Gibbs free energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)

    for line in lines[thermo_start:]:
        if ZPE is None:
            m = _pat_zpe.search(line)
            if m:
                try: ZPE = float(m.group(1))
                except ValueError: pass
        if H_total is None:
            m = _pat_h.search(line)
            if m:
                try: H_total = float(m.group(1))
                except ValueError: pass
        if TS_total is None:
            m = _pat_s.search(line)
            if m:
                try: TS_total = float(m.group(1))
                except ValueError: pass
        if G_total_freq is None:
            m = _pat_g.search(line)
            if m:
                try: G_total_freq = float(m.group(1))
                except ValueError: pass

    # ── 4. Derivar correcciones ───────────────────────────────────────────
    H_corr:  Optional[float] = None
    G_corr:  Optional[float] = None
    S_total: Optional[float] = None   # entropía pura derivada (Eh/K)

    if E_elec_freq is not None and H_total is not None:
        H_corr = H_total - E_elec_freq

    if TS_total is not None and T_parsed is not None and T_parsed != 0:
        # TS_total ('Final entropy term') es T·S evaluado a T_parsed (Eh).
        # Despejamos S puro para poder reescalar a la T solicitada por el
        # caller (parámetro T), que puede diferir de T_parsed.
        S_total = TS_total / T_parsed

    if H_corr is not None and S_total is not None:
        # G_corr a la temperatura T solicitada (no necesariamente T_parsed).
        # Si T == T_parsed, esto reproduce exactamente G_corr = H_corr - TS_total.
        G_corr = H_corr - T * S_total

    # Verificación de consistencia (solo válida si T == T_parsed, ya que
    # G_total_freq de ORCA está fijo a T_parsed):
    if (G_corr is not None and G_total_freq is not None
            and E_elec_freq is not None and T_parsed is not None
            and abs(T - T_parsed) < 1e-6):
        G_reconstructed = E_elec_freq + G_corr
        discrepancy = abs(G_reconstructed - G_total_freq)
        if discrepancy > 1e-5:  # > 6 kcal/mol: algo muy mal
            warnings.warn(
                f"[THERMOCHEM] Discrepancia G reconstruida vs ORCA: "
                f"{discrepancy*HARTREE_TO_KCAL:.4f} kcal/mol en {freq_out}. "
                f"Revisar el output manualmente."
            )

    ok = all(x is not None for x in [ZPE, H_corr, TS_total, G_corr, G_total_freq])

    return {
        "E_elec_freq":    E_elec_freq,
        "ZPE":            ZPE,
        "H_corr":         H_corr,
        "TS_total":       TS_total,          # Eh, T·S a T_parsed (tal cual ORCA)
        "S_total":        S_total,           # Eh/K, derivado = TS_total / T_parsed
        "G_corr":         G_corr,            # Eh, a la temperatura T solicitada
        "G_total_freq":   G_total_freq,      # Eh
        "T_parsed":       T_parsed,
        "n_imag":         len(imag_freqs),
        "imag_freqs_cm1": imag_freqs,
        "ok":             ok,
    }


def parse_gas_phase_sp(sp_gas_out: Path) -> Optional[float]:
    """
    Extrae la energía de un SP en fase gas (sin CPCM).
    Usado en el ciclo de pKa para calcular ΔG_solv de cada especie.
    Retorna la última 'FINAL SINGLE POINT ENERGY' en Hartree, o None.
    """
    return rc_parse_sp_energy(sp_gas_out)


# ─────────────────────────────────────────────
#  ANÁLISIS CREST → DFT
# ─────────────────────────────────────────────

def analyze_crest_to_dft(mol_dir: Path, df_conf: pd.DataFrame) -> pd.DataFrame:
    """
    Compara el orden conformacional de CREST (GFN2-xTB) con el orden DFT post-opt.

    Añade al DataFrame de confórmeros:
      - E_xtb_hartree    : energía xTB del confórmero en CREST (Hartree)
      - dE_xtb_kcal      : ΔE(xTB) respecto al mínimo xTB (kcal/mol)
      - rank_crest        : rango CREST (1 = más estable por xTB)
      - rank_dft          : rango DFT (1 = más estable por SP)
      - rank_change       : rank_dft - rank_crest (0 = no hubo inversión)
      - dE_reranking_kcal : ΔE_DFT - ΔE_xTB (cuánto cambió la brecha energética)

    Retorna df_conf enriquecido.
    """
    conformers_xyz = mol_dir / "01_crest" / "crest_conformers.xyz"
    xtb_energies = parse_crest_energies(conformers_xyz)

    if not xtb_energies:
        warnings.warn("[CREST→DFT] No se pudieron leer energías xTB de CREST.")
        return df_conf

    # Mapear conf_id → índice CREST (conf_001 = índice 0, etc.)
    xtb_map: dict[str, Optional[float]] = {}
    for conf_id in df_conf["conformer"]:
        try:
            idx = int(conf_id.split("_")[1]) - 1   # "conf_003" → 2
        except (IndexError, ValueError):
            continue
        xtb_map[conf_id] = xtb_energies[idx] if idx < len(xtb_energies) else None

    df = df_conf.copy()
    df["E_xtb_hartree"] = df["conformer"].map(xtb_map)

    valid_xtb = df["E_xtb_hartree"].notna()
    valid_dft = df["E_hartree"].notna()

    if valid_xtb.any():
        e_xtb_min = df.loc[valid_xtb, "E_xtb_hartree"].min()
        df["dE_xtb_kcal"] = (df["E_xtb_hartree"] - e_xtb_min) * HARTREE_TO_KCAL

    # Rankings (1 = más estable)
    if valid_xtb.any():
        df.loc[valid_xtb, "rank_crest"] = (
            df.loc[valid_xtb, "E_xtb_hartree"]
            .rank(method="min")
            .astype("Int64")
        )
    if valid_dft.any():
        df.loc[valid_dft, "rank_dft"] = (
            df.loc[valid_dft, "E_hartree"]
            .rank(method="min")
            .astype("Int64")
        )

    if "rank_crest" in df.columns and "rank_dft" in df.columns:
        df["rank_change"] = df["rank_dft"] - df["rank_crest"]

    # Cuánto cambió la brecha energética relativa al pasar de xTB a DFT
    if "dE_xtb_kcal" in df.columns and "dE_kcal" in df.columns:
        df["dE_reranking_kcal"] = df["dE_kcal"] - df["dE_xtb_kcal"]

    return df


def print_crest_dft_analysis(df: pd.DataFrame):
    """Imprime tabla de análisis de reordenamiento CREST→DFT."""
    needed = ["conformer", "rank_crest", "rank_dft", "rank_change",
              "dE_xtb_kcal", "dE_kcal", "dE_reranking_kcal"]
    cols = [c for c in needed if c in df.columns]
    if len(cols) < 3:
        print("  (energías xTB no disponibles — omitiendo análisis CREST→DFT)")
        return

    sub = df[cols].copy()
    inversions = df["rank_change"].abs() if "rank_change" in df.columns else pd.Series(dtype=float)
    n_inv = (inversions > 0).sum() if not inversions.empty else 0

    print(f"\n  ── Reordenamiento conformacional CREST (GFN2-xTB) → DFT ──")
    print(f"  Inversiones de rango detectadas: {n_inv}/{len(df)}")

    def fmt(x):
        if pd.isna(x):
            return "    —   "
        return f"{x:+8.3f}" if isinstance(x, float) else f"{x:>8}"

    # Cabecera
    header_map = {
        "conformer":         "conf",
        "rank_crest":        "rank_xTB",
        "rank_dft":          "rank_DFT",
        "rank_change":       "Δrank",
        "dE_xtb_kcal":       "ΔE_xTB(kcal)",
        "dE_kcal":           "ΔE_DFT(kcal)",
        "dE_reranking_kcal": "ΔΔE(kcal)",
    }
    print("\n  " + "  ".join(f"{header_map.get(c,c):>13}" for c in cols))
    print("  " + "─" * (15 * len(cols)))
    for _, row in sub.iterrows():
        vals = []
        for c in cols:
            v = row[c]
            if c == "conformer":
                vals.append(f"{v:>13}")
            elif pd.isna(v):
                vals.append(f"{'—':>13}")
            elif isinstance(v, (int, np.integer)):
                flag = " ◀" if c == "rank_change" and abs(int(v)) > 0 else ""
                vals.append(f"{int(v):>12}{flag}")
            else:
                vals.append(f"{float(v):>13.4f}")
        print("  " + "  ".join(vals))

    if "dE_reranking_kcal" in df.columns:
        max_shift = df["dE_reranking_kcal"].abs().max()
        print(f"\n  Máximo desplazamiento energético xTB→DFT: {max_shift:.3f} kcal/mol")
        if max_shift > 2.0:
            print("  ⚠  La función de energía xTB difiere significativamente de DFT")
            print("     en este sistema. Recomendable ampliar el muestreo CREST.")


# ─────────────────────────────────────────────
#  CICLO TERMODINÁMICO COMPUESTO + BOLTZMANN
# ─────────────────────────────────────────────

def compute_G_total(df_conf: pd.DataFrame) -> pd.DataFrame:
    """
    Aplica el ciclo termodinámico compuesto para calcular G_total de cada
    confórmero, usando las temperaturas definidas en TEMPS_K.

    Ciclo implementado (estándar en literatura de MAAs, Ho & Coote 2010):

        G_total_i(T) = E_SP(TZVPD/CPCM) + G_corr_i(T, SVP/CPCM)

    donde G_corr_i(T) = H_corr_i - T·S_total_i
    y H_corr, S_total provienen del NumFreq a def2-SVP/CPCM.

    IMPORTANTE — temperatura del NumFreq:
    El NumFreq de ORCA se corre a 298.15 K por defecto. Para las otras
    temperaturas (250 K, 274.15 K, 310 K) se recalcula G_corr(T) usando
    la aproximación rígido-rotor/oscilador-armónico (RRHO):

        G_corr(T) = H_corr(T_ref) - T · S_total(T_ref)

    Esto es una aproximación válida cuando |T - T_ref| < 60 K (Grimme 2012),
    adecuada para nuestro rango 250–310 K. Para rigor máximo se debería
    relanzar el NumFreq a cada T, pero el costo computacional no lo justifica
    en este sistema.

    Columnas añadidas al DataFrame:
        G_corr_298K       : corrección de Gibbs a 298.15 K (Hartree)
        G_total_{T}K      : energía libre de Gibbs compuesta a cada T (Hartree)
        dG_{T}K_kcal      : ΔG relativo al mínimo a cada T (kcal/mol)
        H_corr            : corrección entálpica (Hartree, independiente de T)
        S_total_EhK       : entropía total (Eh/K)
        ZPE               : energía de punto cero (Hartree)
        thermo_ok         : True si el NumFreq fue parseado correctamente
        thermo_warning    : descripción de problemas encontrados
    """
    df = df_conf.copy()

    # Columnas nuevas con valores por defecto
    df["G_corr_298K"]     = np.nan
    df["H_corr"]          = np.nan
    df["S_total_EhK"]     = np.nan
    df["ZPE"]             = np.nan
    df["thermo_ok"]       = False
    df["thermo_warning"]  = ""

    for label in TEMP_LABELS:
        df[f"G_total_{label}"] = np.nan
        df[f"dG_{label}_kcal"] = np.nan

    for idx, row in df.iterrows():
        tc = row.get("_thermochem", None)
        if tc is None or not isinstance(tc, dict):
            df.at[idx, "thermo_warning"] = "sin datos termoquímicos"
            continue

        if not tc.get("ok", False):
            issues = []
            for field in ("ZPE", "H_corr", "S_total", "G_corr"):
                if tc.get(field) is None:
                    issues.append(field)
            df.at[idx, "thermo_warning"] = f"campos faltantes: {', '.join(issues)}"
            # Intentar uso parcial si al menos H_corr y S_total están
            if tc.get("H_corr") is None or tc.get("S_total") is None:
                continue

        e_sp = row.get("E_hartree")
        if e_sp is None or pd.isna(e_sp):
            df.at[idx, "thermo_warning"] = (
                df.at[idx, "thermo_warning"] + " | E_SP faltante"
            ).lstrip(" | ")
            continue

        H_corr  = tc["H_corr"]
        S_total = tc["S_total"]   # Eh/K
        ZPE     = tc.get("ZPE")

        df.at[idx, "H_corr"]       = H_corr
        df.at[idx, "S_total_EhK"]  = S_total
        df.at[idx, "ZPE"]          = ZPE
        df.at[idx, "thermo_ok"]    = True

        # G_corr a cada T  (aproximación RRHO, válida ±60 K respecto a T_ref)
        for T, label in zip(TEMPS_K, TEMP_LABELS):
            G_corr_T = H_corr - T * S_total
            G_total  = e_sp + G_corr_T
            df.at[idx, f"G_total_{label}"] = G_total

        # G_corr de referencia a 298.15 K
        df.at[idx, "G_corr_298K"] = H_corr - T_STD * S_total

    # ΔG relativo al mínimo para cada temperatura
    n_ok = df["thermo_ok"].sum()
    for T, label in zip(TEMPS_K, TEMP_LABELS):
        col_g   = f"G_total_{label}"
        col_dg  = f"dG_{label}_kcal"
        valid_g = df["thermo_ok"] & df[col_g].notna()
        if valid_g.any():
            g_min = df.loc[valid_g, col_g].min()
            df.loc[valid_g, col_dg] = (df.loc[valid_g, col_g] - g_min) * HARTREE_TO_KCAL

    # ΔE electrónica pura (de SP) — mantenida para comparación
    valid_e = df["E_hartree"].notna()
    if valid_e.any():
        e_min = df.loc[valid_e, "E_hartree"].min()
        df["dE_kcal"] = (df["E_hartree"] - e_min) * HARTREE_TO_KCAL

    # Advertencia global si hay confórmeros sin termoquímica
    n_missing = (~df["thermo_ok"]).sum()
    if n_missing > 0:
        warnings.warn(
            f"[THERMOCHEM] {n_missing}/{len(df)} confórmeros sin correcciones "
            f"térmicas válidas. El Boltzmann de esos confs se omitirá. "
            f"Detalles en columna 'thermo_warning'."
        )

    return df


def compute_boltzmann_multiT(df_conf: pd.DataFrame) -> pd.DataFrame:
    """
    Calcula pesos de Boltzmann a todas las temperaturas en TEMPS_K.

    Usa ΔG = G_total(T) - G_min(T) cuando las correcciones térmicas están
    disponibles (ciclo compuesto completo). Si thermo_ok es False para algún
    confórmero, ese confórmero queda excluido del ensemble y su peso es NaN.

    Fallback: si NINGÚN confórmero tiene termoquímica, usa ΔE electrónica
    pura con advertencia explícita en df.attrs["boltzmann_mode"].

    Métricas adicionales calculadas para cada T:
        N_eff(T)      = exp(-Σ p_i · ln(p_i))  (número efectivo de confs)
        S_mix(T)/R    = -Σ p_i · ln(p_i)       (entropía de mezcla en unidades R)

    Columna de conveniencia:
        boltz_298K    = alias de boltz_298K (label "298K"), para compatibilidad
                        con el resto del código.
    """
    df = df_conf.copy()

    # Determinar si el ciclo completo está disponible
    has_thermo = "thermo_ok" in df.columns and df["thermo_ok"].any()

    if not has_thermo:
        warnings.warn(
            "[BOLTZMANN] No hay correcciones térmicas disponibles. "
            "Usando solo energía electrónica SP (ΔE). "
            "Los pesos Boltzmann NO incluyen ZPE ni entropía conformacional. "
            "Reportar como 'Boltzmann weighting based on electronic energies'."
        )
        df.attrs["boltzmann_mode"] = "electronic_only"
        valid = df["E_hartree"].notna()
        e_min = df.loc[valid, "E_hartree"].min()
        df.loc[valid, "dE_kcal"] = (df.loc[valid, "E_hartree"] - e_min) * HARTREE_TO_KCAL

        for T, label in zip(TEMPS_K, TEMP_LABELS):
            col = f"boltz_{label}"
            raw = np.exp(-df.loc[valid, "dE_kcal"] / (KB_KCAL * T))
            df.loc[valid, col] = raw / raw.sum()
    else:
        df.attrs["boltzmann_mode"] = "gibbs_composite"
        for T, label in zip(TEMPS_K, TEMP_LABELS):
            col_g  = f"G_total_{label}"
            col_dg = f"dG_{label}_kcal"
            col_b  = f"boltz_{label}"

            valid = df["thermo_ok"] & df[col_g].notna()
            if not valid.any():
                warnings.warn(f"[BOLTZMANN] Sin G_total válido a {label}.")
                continue

            dG = df.loc[valid, col_dg]
            raw = np.exp(-dG / (KB_KCAL * T))
            df.loc[valid, col_b] = raw / raw.sum()

    # ── Métricas del ensemble ─────────────────────────────────────────────
    for T, label in zip(TEMPS_K, TEMP_LABELS):
        col_b = f"boltz_{label}"
        if col_b not in df.columns:
            continue
        valid_b = df[col_b].notna()
        if not valid_b.any():
            continue
        p = df.loc[valid_b, col_b].clip(lower=1e-15)
        s_mix  = -(p * np.log(p)).sum()     # S_mix / R (adimensional)
        n_eff  = np.exp(s_mix)
        df.attrs[f"S_mix_{label}"] = float(s_mix)
        df.attrs[f"N_eff_{label}"] = float(n_eff)
        df.attrs[f"n_above_1pct_{label}"] = int((df[col_b] > 0.01).sum())

    # ── Alias de conveniencia para 298K ──────────────────────────────────
    # La etiqueta "298K" corresponde al índice 2 en TEMP_LABELS
    lbl_298 = "298K"
    col_298 = f"boltz_{lbl_298}"
    if col_298 in df.columns:
        df["boltz_298K"] = df[col_298]
        df.attrs["S_mix_298K"] = df.attrs.get(f"S_mix_{lbl_298}")
        df.attrs["N_eff_298K"] = df.attrs.get(f"N_eff_{lbl_298}")
        df.attrs["n_above_1pct_298K"] = df.attrs.get(f"n_above_1pct_{lbl_298}")

    return df


def print_boltzmann_analysis(df: pd.DataFrame):
    """Imprime tabla de pesos Boltzmann multi-temperatura con modo termodinámico."""
    mode = df.attrs.get("boltzmann_mode", "unknown")
    mode_str = {
        "gibbs_composite":  "ΔG compuesto [E_SP(TZVPD) + G_corr(SVP/CPCM)]",
        "electronic_only":  "ΔE electrónica pura — SIN ZPE/entropía ⚠",
        "unknown":          "modo desconocido",
    }.get(mode, mode)

    print(f"\n  ── Pesos de Boltzmann por temperatura ──")
    print(f"  Modo: {mode_str}")

    avail_labels = [lbl for lbl in TEMP_LABELS if f"boltz_{lbl}" in df.columns]
    if not avail_labels:
        print("  (sin columnas boltz_* disponibles)")
        return

    # Encabezado con temperatura en K
    header = f"\n  {'conf':>10}  "
    header += "  ".join(f"{lbl:>8}" for lbl in avail_labels)
    if mode == "gibbs_composite":
        header += "   ΔG_298(kcal)"
    else:
        header += "   ΔE(kcal)"
    print(header)
    print("  " + "─" * (14 + 10 * len(avail_labels) + 16))

    for _, row in df.iterrows():
        thermo_flag = ""
        if "thermo_ok" in df.columns and not row.get("thermo_ok", True):
            thermo_flag = " [no G_corr]"

        vals = [f"  {row['conformer']:>10}"]
        for lbl in avail_labels:
            col = f"boltz_{lbl}"
            v   = row.get(col, np.nan)
            if pd.isna(v):
                vals.append(f"  {'—':>8}")
            else:
                marker = " ★" if v > 0.50 else (" ·" if v > 0.10 else "  ")
                vals.append(f"  {v:6.4f}{marker}")

        # ΔG o ΔE según modo
        dg_col = "dG_298K_kcal" if mode == "gibbs_composite" else "dE_kcal"
        dg_val = row.get(dg_col, np.nan)
        if pd.isna(dg_val):
            vals.append(f"   {'—':>12}")
        else:
            vals.append(f"   {dg_val:+12.4f}")

        print("".join(vals) + thermo_flag)

    # Métricas del ensemble
    print()
    for lbl in avail_labels:
        s   = df.attrs.get(f"S_mix_{lbl}")
        n   = df.attrs.get(f"N_eff_{lbl}")
        n1  = df.attrs.get(f"n_above_1pct_{lbl}")
        if s is not None:
            print(f"  {lbl:>5}:  S_mix/R = {s:.4f}  │  N_eff = {n:.2f}  │  "
                  f"confs >1% = {n1}/{len(df)}")

    # Sensibilidad a temperatura: ¿cambia el conf dominante entre extremos?
    b_low  = f"boltz_{TEMP_LABELS[0]}"    # 250K
    b_high = f"boltz_{TEMP_LABELS[-1]}"   # 310K
    if b_low in df.columns and b_high in df.columns:
        valid_low  = df[b_low].notna()
        valid_high = df[b_high].notna()
        if valid_low.any() and valid_high.any():
            dom_low  = df.loc[df[b_low].idxmax(),  "conformer"]
            dom_high = df.loc[df[b_high].idxmax(), "conformer"]
            if dom_low != dom_high:
                print(f"\n  ⚠  Confórmero dominante cambia con T:")
                print(f"     {TEMP_LABELS[0]} → {dom_low}   │   "
                      f"{TEMP_LABELS[-1]} → {dom_high}")
            else:
                print(f"\n  Confórmero dominante estable entre "
                      f"{TEMP_LABELS[0]}–{TEMP_LABELS[-1]}: {dom_low}")


# ─────────────────────────────────────────────
#  CICLO TERMODINÁMICO DE pKa
# ─────────────────────────────────────────────

def rc_compute_pka(
    mol_acid: str,
    mol_base: str,
    T: float = T_STD,
) -> dict:
    """
    Calcula el pKa teórico de la reacción AH(aq) → A⁻(aq) + H⁺(aq) usando
    el ciclo termodinámico de Born (ciclo termoquímico de dos-ramas).

    ════════════════════════════════════════════════════════════════════════
    CICLO DE BORN (esquema de dos niveles)
    ════════════════════════════════════════════════════════════════════════

      AH(gas) ──────────────────────────────────────► A⁻(gas) + H⁺(gas)
         │                  ΔG_gas                          │
     ΔG_solv(AH)                              ΔG_solv(A⁻) + ΔG_solv(H⁺)
         │                                                  │
         ▼                                                  ▼
      AH(aq)  ──────────────────────────────────────► A⁻(aq) + H⁺(aq)
                              ΔG_aq

    ΔG_gas   = G_total(A⁻, gas) + G(H⁺, gas) - G_total(AH, gas)
    ΔG_solv  = E_SP(CPCM) - E_SP(gas)   para cada especie
    ΔG_aq    = ΔG_gas + ΔΔG_solv
    pKa      = ΔG_aq / (RT · ln10)

    ════════════════════════════════════════════════════════════════════════
    NIVEL DE TEORÍA
    ════════════════════════════════════════════════════════════════════════

    G_total en fase gas se aproxima como:
        G_gas(X) = E_SP_gas(X) + G_corr(X)
    donde G_corr viene del NumFreq (que usa CPCM). Esta es la aproximación
    de Ho & Coote (Theor. Chem. Acc. 2010): la corrección térmica es casi
    independiente del solvente para moléculas rígidas.

    ΔG_solv(X) = E_SP_aq(X, CPCM) - E_SP_gas(X)
    usando el mismo funcional/base en ambos (ωB97X-D3/def2-TZVPD).

    ════════════════════════════════════════════════════════════════════════
    CONSTANTES DE REFERENCIA (todas a 298.15 K)
    ════════════════════════════════════════════════════════════════════════

    G(H⁺, gas)    = -6.28 kcal/mol    [Tissandier et al. 1998]
    ΔG_solv(H⁺)   = -265.9 kcal/mol   [Kelly, Cramer & Truhlar 2006]
    δG_conc       = +1.894 kcal/mol    [corrección 1 atm → 1 M, RT·ln(24.5)]

    ════════════════════════════════════════════════════════════════════════
    PARÁMETROS
    ════════════════════════════════════════════════════════════════════════

    mol_acid : nombre del directorio en runs/ para AH (forma protonada)
    mol_base : nombre del directorio en runs/ para A⁻ (forma desprotonada)
    T        : temperatura de cálculo (default 298.15 K)

    ════════════════════════════════════════════════════════════════════════
    ARCHIVOS REQUERIDOS
    ════════════════════════════════════════════════════════════════════════

    Para cada mol (acid y base):
      · runs/<mol>/04_sp/<best_conf>/sp.out          — SP acuoso (CPCM)
      · runs/<mol>/06_sp_gas/<best_conf>/sp_gas.out  — SP fase gas (sin CPCM)
      · runs/<mol>/03_numfreq/<best_conf>/freq.out   — correcciones térmicas

    El confórmero representativo es el de menor G_total a 298 K.
    Si no existe 06_sp_gas, el ciclo de pKa no puede calcularse.

    ════════════════════════════════════════════════════════════════════════
    RETORNA
    ════════════════════════════════════════════════════════════════════════

    dict con:
        pKa              : valor teórico (float o None si falta algún dato)
        dG_aq_kcal       : ΔG_aq en kcal/mol
        dG_gas_kcal      : ΔG_gas en kcal/mol
        dG_solv_acid     : ΔG_solv(AH) en kcal/mol
        dG_solv_base     : ΔG_solv(A⁻) en kcal/mol
        dG_solv_hplus    : ΔG_solv(H⁺) constante de literatura
        G_total_acid_gas : G_total(AH) en fase gas (Hartree)
        G_total_base_gas : G_total(A⁻) en fase gas (Hartree)
        conf_acid        : confórmero representativo de AH
        conf_base        : confórmero representativo de A⁻
        T                : temperatura usada (K)
        ok               : True si el pKa se calculó correctamente
        warnings         : lista de advertencias
    """
    result: dict = {
        "pKa":              None,
        "dG_aq_kcal":       None,
        "dG_gas_kcal":      None,
        "dG_solv_acid":     None,
        "dG_solv_base":     None,
        "dG_solv_hplus":    DG_SOLV_HPLUS,
        "G_total_acid_gas": None,
        "G_total_base_gas": None,
        "conf_acid":        None,
        "conf_base":        None,
        "T":                T,
        "ok":               False,
        "warnings":         [],
    }

    def _warn(msg: str):
        result["warnings"].append(msg)
        warnings.warn(f"[pKa] {msg}")

    # ── Cargar resultados de ambas formas ─────────────────────────────────
    try:
        df_acid, _ = collect_results(mol_acid)
        df_base, _ = collect_results(mol_base)
    except Exception as e:
        _warn(f"Error al cargar resultados: {e}")
        return result

    if df_acid.empty or df_base.empty:
        _warn("DataFrames vacíos para ácido o base.")
        return result

    # ── Confórmero representativo: mínimo G_total a 298.15 K ─────────────
    lbl_298 = "298K"
    col_g298 = f"G_total_{lbl_298}"

    def _best_conf(df: pd.DataFrame, mol: str) -> Optional[str]:
        if col_g298 in df.columns and df["thermo_ok"].any():
            valid = df["thermo_ok"] & df[col_g298].notna()
            if valid.any():
                return df.loc[df.loc[valid, col_g298].idxmin(), "conformer"]
        # Fallback: mínima E_SP
        if "E_hartree" in df.columns and df["E_hartree"].notna().any():
            _warn(f"[{mol}] Usando E_SP (sin G_corr) para seleccionar confórmero.")
            return df.loc[df["E_hartree"].idxmin(), "conformer"]
        return None

    conf_acid = _best_conf(df_acid, mol_acid)
    conf_base = _best_conf(df_base, mol_base)
    result["conf_acid"] = conf_acid
    result["conf_base"] = conf_base

    if conf_acid is None or conf_base is None:
        _warn("No se pudo determinar el confórmero representativo.")
        return result

    # ── Energías SP acuosas (CPCM) ───────────────────────────────────────
    acid_sp_aq_file = RC_WORK_ROOT / mol_acid / "04_sp" / conf_acid / "sp.out"
    base_sp_aq_file = RC_WORK_ROOT / mol_base / "04_sp" / conf_base / "sp.out"

    E_acid_aq = rc_parse_sp_energy(acid_sp_aq_file)
    E_base_aq = rc_parse_sp_energy(base_sp_aq_file)

    if E_acid_aq is None:
        _warn(f"No se encontró SP acuoso de AH en {acid_sp_aq_file}")
    if E_base_aq is None:
        _warn(f"No se encontró SP acuoso de A⁻ en {base_sp_aq_file}")

    # ── Energías SP en fase gas ───────────────────────────────────────────
    acid_sp_gas_file = (RC_WORK_ROOT / mol_acid / PKA_GAS_SP_SUBDIR
                        / conf_acid / "sp_gas.out")
    base_sp_gas_file = (RC_WORK_ROOT / mol_base / PKA_GAS_SP_SUBDIR
                        / conf_base / "sp_gas.out")

    E_acid_gas = parse_gas_phase_sp(acid_sp_gas_file)
    E_base_gas = parse_gas_phase_sp(base_sp_gas_file)

    if E_acid_gas is None:
        _warn(
            f"SP en fase gas de AH no encontrado en {acid_sp_gas_file}. "
            f"Corre: orca sp_gas.inp > sp_gas.out (sin CPCM, mismo funcional/base)."
        )
    if E_base_gas is None:
        _warn(
            f"SP en fase gas de A⁻ no encontrado en {base_sp_gas_file}. "
            f"Corre: orca sp_gas.inp > sp_gas.out (sin CPCM, mismo funcional/base)."
        )

    if E_acid_gas is None or E_base_gas is None:
        return result  # no hay suficiente para continuar

    # ── Correcciones térmicas (G_corr) ────────────────────────────────────
    acid_freq_file = RC_WORK_ROOT / mol_acid / "03_numfreq" / conf_acid / "freq.out"
    base_freq_file = RC_WORK_ROOT / mol_base / "03_numfreq" / conf_base / "freq.out"

    tc_acid = rc_parse_thermochem(acid_freq_file, T=T)
    tc_base = rc_parse_thermochem(base_freq_file, T=T)

    G_corr_acid: Optional[float] = None
    G_corr_base: Optional[float] = None

    if tc_acid.get("ok"):
        G_corr_acid = tc_acid["G_corr"]
    else:
        _warn(
            f"Correcciones térmicas de AH no disponibles o incompletas "
            f"({acid_freq_file}). Se usará G_gas ≈ E_gas (sin ZPE ni entropía). "
            f"El pKa resultante tendrá mayor incertidumbre."
        )

    if tc_base.get("ok"):
        G_corr_base = tc_base["G_corr"]
    else:
        _warn(
            f"Correcciones térmicas de A⁻ no disponibles o incompletas "
            f"({base_freq_file}). Se usará G_gas ≈ E_gas."
        )

    # ── G_total en fase gas para cada especie ─────────────────────────────
    #   G_gas = E_SP_gas + G_corr(RRHO, de NumFreq con CPCM)
    #   La aproximación Ho & Coote justifica usar G_corr del cálculo CPCM
    #   incluso para la fase gas (la corrección térmica es casi la misma).
    G_acid_gas = E_acid_gas + (G_corr_acid if G_corr_acid is not None else 0.0)
    G_base_gas = E_base_gas + (G_corr_base if G_corr_base is not None else 0.0)

    result["G_total_acid_gas"] = G_acid_gas
    result["G_total_base_gas"] = G_base_gas

    # ── ΔG_gas (en kcal/mol) ─────────────────────────────────────────────
    #   ΔG_gas = G(A⁻, gas) + G(H⁺, gas) - G(AH, gas)
    dG_gas_hartree = G_base_gas - G_acid_gas   # en Hartree (H⁺ se suma aparte)
    dG_gas_kcal    = dG_gas_hartree * HARTREE_TO_KCAL + G_HPLUS_GAS
    result["dG_gas_kcal"] = dG_gas_kcal

    # ── ΔG_solv para cada especie ─────────────────────────────────────────
    #   ΔG_solv(X) = E_SP_aq(X, CPCM) - E_SP_gas(X)
    #   En kcal/mol. Incluye corrección de concentración estándar δG_conc.
    if E_acid_aq is not None:
        dG_solv_acid = (E_acid_aq - E_acid_gas) * HARTREE_TO_KCAL + DG_CONC_CORR
        result["dG_solv_acid"] = dG_solv_acid
    else:
        _warn("No se puede calcular ΔG_solv(AH) sin E_SP acuoso.")
        return result

    if E_base_aq is not None:
        dG_solv_base = (E_base_aq - E_base_gas) * HARTREE_TO_KCAL + DG_CONC_CORR
        result["dG_solv_base"] = dG_solv_base
    else:
        _warn("No se puede calcular ΔG_solv(A⁻) sin E_SP acuoso.")
        return result

    # ── ΔG_aq (ciclo de Born completo) ───────────────────────────────────
    #   ΔG_aq = ΔG_gas + ΔG_solv(A⁻) + ΔG_solv(H⁺) - ΔG_solv(AH)
    dG_aq_kcal = (
        dG_gas_kcal
        + dG_solv_base
        + DG_SOLV_HPLUS
        - dG_solv_acid
    )
    result["dG_aq_kcal"] = dG_aq_kcal

    # ── pKa ──────────────────────────────────────────────────────────────
    pKa = dG_aq_kcal / (R_KCAL * T * LN10)
    result["pKa"]  = pKa
    result["ok"]   = True

    return result


def print_pka_analysis(pka_result: dict, mol_acid: str, mol_base: str):
    """Imprime el análisis completo del ciclo termodinámico de pKa."""
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  ANÁLISIS DE pKa TEÓRICO")
    print(f"  AH = {mol_acid}   │   A⁻ = {mol_base}")
    print(f"  T  = {pka_result['T']:.2f} K")
    print(sep)

    if not pka_result["ok"]:
        print("\n  ✗  El ciclo de pKa no pudo completarse.")
        for w in pka_result["warnings"]:
            print(f"     ⚠  {w}")
        return

    print(f"\n  Confórmero representativo AH : {pka_result['conf_acid']}")
    print(f"  Confórmero representativo A⁻ : {pka_result['conf_base']}")

    print(f"\n  ── Contribuciones al ciclo de Born ──")
    print(f"  {'ΔG_gas':>30} = {pka_result['dG_gas_kcal']:+10.3f} kcal/mol")
    print(f"  {'ΔG_solv(AH)':>30} = {pka_result['dG_solv_acid']:+10.3f} kcal/mol")
    print(f"  {'ΔG_solv(A⁻)':>30} = {pka_result['dG_solv_base']:+10.3f} kcal/mol")
    print(f"  {'ΔG_solv(H⁺) [Kelly 2006]':>30} = {pka_result['dG_solv_hplus']:+10.3f} kcal/mol")
    print(f"  {'─'*48}")
    print(f"  {'ΔG_aq':>30} = {pka_result['dG_aq_kcal']:+10.3f} kcal/mol")
    print(f"\n  pKa (teórico, {pka_result['T']:.2f} K) = {pka_result['pKa']:.2f}")
    print(f"\n  Nota: incertidumbre típica ±1–2 unidades de pKa para ciclos")
    print(f"  Born/CPCM a nivel ωB97X-D3/def2-TZVPD (Ho & Coote 2010).")

    if pka_result["warnings"]:
        print(f"\n  Advertencias:")
        for w in pka_result["warnings"]:
            print(f"    ⚠  {w}")

    print(f"\n{sep}\n")


# ─────────────────────────────────────────────
#  COLECTOR PRINCIPAL
# ─────────────────────────────────────────────

def collect_results(mol_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Devuelve (df_conf, df_tddft) con todos los campos enriquecidos:
      df_conf  — una fila por confórmero
      df_tddft — una fila por excitación TDDFT de cada confórmero
    """
    mol_dir = RC_WORK_ROOT / mol_name
    if not mol_dir.exists():
        raise FileNotFoundError(f"No se encontró {mol_dir}")

    conf_rows: list[dict] = []
    tddft_rows: list[dict] = []

    sp_dirs = sorted((mol_dir / "04_sp").glob("conf_*")) \
              if (mol_dir / "04_sp").exists() else []

    for sp_dir in sp_dirs:
        conf_id  = sp_dir.name
        sp_out   = sp_dir / "sp.out"
        e_sp     = rc_parse_sp_energy(sp_out)

        # Energía DFT post-opt (de la optimización, no del SP)
        opt_out  = mol_dir / "02_opt" / conf_id / "opt.out"
        e_opt    = parse_opt_energy(opt_out)

        # Frecuencias imaginarias
        freq_out   = mol_dir / "03_numfreq" / conf_id / "freq.out"
        imag_freqs = parse_imaginary_freqs(freq_out)

        conf_rows.append({
            "conformer":      conf_id,
            "E_hartree":      e_sp,       # SP de alta precisión (def2-TZVPD)
            "E_opt_hartree":  e_opt,      # energía al final de la Opt (def2-SVP)
            "has_imaginary":  len(imag_freqs) > 0,
            "n_imaginary":    len(imag_freqs),
            "imag_freqs_cm1": str(imag_freqs) if imag_freqs else "—",
        })

        # TDDFT
        tddft_out = mol_dir / "05_tddft" / conf_id / "tddft.out"
        for exc in parse_tddft_excitations(tddft_out):
            tddft_rows.append({"conformer": conf_id, **exc})

    df_conf = pd.DataFrame(conf_rows)

    if df_conf.empty:
        print(f"⚠  No se encontraron resultados SP para '{mol_name}'")
        return df_conf, pd.DataFrame(tddft_rows)

    # ── Correcciones térmicas ORCA (NumFreq) por confórmero ────────────────
    # compute_boltzmann_multiT() necesita la columna "thermo_ok" (creada por
    # compute_G_total) para usar el ciclo de Gibbs compuesto; sin este paso
    # caía siempre al modo "electronic_only" (ΔE pura) aunque el NumFreq
    # se hubiera parseado bien.
    _thermochem_cache: list[Optional[dict]] = []
    for _, row in df_conf.iterrows():
        conf_id  = row["conformer"]
        freq_out = mol_dir / "03_numfreq" / conf_id / "freq.out"
        tc = rc_parse_thermochem(freq_out, T=T_STD)
        _thermochem_cache.append(tc if tc.get("ok") else None)
    df_conf["_thermochem"] = _thermochem_cache
    df_conf = compute_G_total(df_conf)

    # Boltzmann multi-temperatura
    df_conf = compute_boltzmann_multiT(df_conf)

    # Análisis CREST→DFT
    df_conf = analyze_crest_to_dft(mol_dir, df_conf)

    df_tddft = pd.DataFrame(tddft_rows)

    # λ_max Boltzmann-promediado (S1, 298 K) — adjunto como atributo del DF
    if not df_tddft.empty and "boltz_298K" in df_conf.columns:
        s1 = df_tddft[df_tddft["state"] == 1].copy()
        s1 = s1.merge(df_conf[["conformer", "boltz_298K"]], on="conformer", how="left")
        w = s1["boltz_298K"].fillna(0)
        if w.sum() > 0:
            df_conf.attrs["lambda_max_boltz_298K"] = float(
                (s1["wavelength_nm"] * w).sum() / w.sum()
            )
            df_conf.attrs["f_max_boltz_298K"] = float(
                (s1["fosc"] * w).sum() / w.sum()
            )

    return df_conf, df_tddft


# ─────────────────────────────────────────────
#  IMPRESIÓN DE RESULTADOS
# ─────────────────────────────────────────────

def print_summary(df_conf: pd.DataFrame, df_tddft: pd.DataFrame, mol_name: str):
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  RESULTADOS: {mol_name}")
    print(sep)

    if df_conf.empty:
        print("  Sin datos de confórmeros.")
        return

    # ── Tabla de confórmeros base ──────────────────────────────────────────
    print("\n  ── Confórmeros (SP def2-TZVPD / CPCM) ──")
    base_cols = ["conformer", "E_hartree", "dE_kcal", "boltz_298K",
                 "has_imaginary", "n_imaginary"]
    base_cols = [c for c in base_cols if c in df_conf.columns]
    print(df_conf[base_cols].to_string(
        index=False, float_format=lambda x: f"{x:.6f}"
    ))

    lmax = df_conf.attrs.get("lambda_max_boltz_298K")
    if lmax:
        fmax = df_conf.attrs.get("f_max_boltz_298K", float("nan"))
        print(f"\n  λ_max Boltzmann-promediado (298 K, S1): {lmax:.1f} nm  "
              f"(f_avg = {fmax:.5f})")

    # ── Análisis CREST→DFT ────────────────────────────────────────────────
    print_crest_dft_analysis(df_conf)

    # ── Boltzmann multi-temperatura ────────────────────────────────────────
    print_boltzmann_analysis(df_conf)

    # ── Excitaciones TDDFT ────────────────────────────────────────────────
    if not df_tddft.empty:
        print("\n  ── Excitaciones S1 por confórmero ──")
        s1 = df_tddft[df_tddft["state"] == 1].copy()
        if not s1.empty:
            s1 = s1.merge(
                df_conf[["conformer", "boltz_298K", "dE_kcal"]],
                on="conformer", how="left"
            )
            disp_cols = ["conformer", "wavelength_nm", "E_eV", "fosc",
                         "boltz_298K", "dE_kcal"]
            disp_cols = [c for c in disp_cols if c in s1.columns]
            print(s1[disp_cols].to_string(
                index=False, float_format=lambda x: f"{x:.4f}"
            ))

        # Top 5 excitaciones por fosc ponderada
        if "boltz_298K" in df_conf.columns:
            top = df_tddft.merge(
                df_conf[["conformer", "boltz_298K"]], on="conformer", how="left"
            ).copy()
            top["f_weighted"] = top["fosc"] * top["boltz_298K"].fillna(0)
            top5 = top.nlargest(5, "f_weighted")[
                ["conformer", "state", "wavelength_nm", "E_eV", "fosc", "f_weighted"]
            ]
            print("\n  ── Top 5 excitaciones por intensidad ponderada (298 K) ──")
            print(top5.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    print(f"\n{sep}\n")


# ─────────────────────────────────────────────
#  COMPARACIÓN ENTRE ESTADOS DE PROTONACIÓN
# ─────────────────────────────────────────────

def compare_protonation_states(
    ref_mol: Optional[str] = None,
) -> pd.DataFrame:
    """
    Recopila resultados de todos los estados de protonación disponibles en runs/
    y construye una tabla comparativa.

    La tabla incluye por mol:
      - Número de confórmeros con SP completado
      - Confórmero de mínima energía
      - Distribución de peso Boltzmann (confs > 1%, > 5%, > 20%)
      - N_eff a 298 K y 250 K
      - λ_max Boltzmann-promediado (S1, 298 K)
      - f_avg (298 K)
      - E_min (Hartree) — para comparar entre estados si mismo nivel de teoría

    Si ref_mol se especifica, calcula también ΔE relativa entre estados
    (solo significativa si tienen la misma carga y nivel de teoría).

    Retorna DataFrame de comparación.
    """
    mols = sorted(p.name for p in RC_WORK_ROOT.iterdir()
                  if p.is_dir() and (p / "04_sp").exists())

    if not mols:
        print("No se encontraron resultados SP en ninguna molécula bajo runs/")
        return pd.DataFrame()

    rows = []
    dfs_conf: dict[str, pd.DataFrame] = {}
    dfs_tddft: dict[str, pd.DataFrame] = {}

    for mol in mols:
        try:
            df_c, df_t = collect_results(mol)
        except Exception as e:
            print(f"  ⚠  {mol}: error al leer resultados — {e}")
            continue

        dfs_conf[mol]  = df_c
        dfs_tddft[mol] = df_t

        if df_c.empty:
            continue

        valid = df_c["E_hartree"].notna()
        n_confs = valid.sum()
        e_min   = df_c.loc[valid, "E_hartree"].min() if n_confs > 0 else None
        best    = df_c.loc[df_c["E_hartree"].idxmin(), "conformer"] \
                  if n_confs > 0 else "—"

        # Distribución conformacional (298 K)
        col298 = "boltz_298K"
        n_1pct  = int((df_c[col298] > 0.01).sum()) if col298 in df_c.columns else None
        n_5pct  = int((df_c[col298] > 0.05).sum()) if col298 in df_c.columns else None
        n_20pct = int((df_c[col298] > 0.20).sum()) if col298 in df_c.columns else None
        n_eff   = df_c.attrs.get("N_eff_298K")

        # N_eff a 250 K
        col250 = "boltz_250K"
        if col250 in df_c.columns:
            p250 = df_c.loc[valid, col250].clip(lower=1e-15)
            s250 = -(p250 * np.log(p250)).sum()
            n_eff_250 = np.exp(s250)
        else:
            n_eff_250 = None

        lmax = df_c.attrs.get("lambda_max_boltz_298K")
        favg = df_c.attrs.get("f_max_boltz_298K")

        # Inversiones CREST→DFT
        n_inv = None
        if "rank_change" in df_c.columns:
            n_inv = int((df_c["rank_change"].abs() > 0).sum())

        rows.append({
            "molecule":       mol,
            "n_confs_sp":     n_confs,
            "best_conformer": best,
            "E_min_hartree":  e_min,
            "N_eff_298K":     n_eff,
            "N_eff_250K":     n_eff_250,
            "n_confs_>1%":    n_1pct,
            "n_confs_>5%":    n_5pct,
            "n_confs_>20%":   n_20pct,
            "lambda_max_nm":  lmax,
            "f_avg_298K":     favg,
            "n_rank_inversions": n_inv,
        })

    df_cmp = pd.DataFrame(rows)

    if ref_mol and ref_mol in df_cmp["molecule"].values:
        e_ref = df_cmp.loc[df_cmp["molecule"] == ref_mol, "E_min_hartree"].values[0]
        if e_ref is not None:
            df_cmp["dE_vs_ref_kcal"] = (df_cmp["E_min_hartree"] - e_ref) * HARTREE_TO_KCAL
            df_cmp["dE_vs_ref_kJ"]   = df_cmp["dE_vs_ref_kcal"] * 4.184

    return df_cmp, dfs_conf, dfs_tddft


def print_comparison(df_cmp: pd.DataFrame, ref_mol: Optional[str] = None):
    sep = "─" * 80
    print(f"\n{sep}")
    print(f"  COMPARACIÓN DE ESTADOS DE PROTONACIÓN  "
          f"({'ref: ' + ref_mol if ref_mol else 'sin referencia'})")
    print(sep)

    if df_cmp.empty:
        print("  Sin datos para comparar.")
        return

    # Tabla principal
    core_cols = ["molecule", "n_confs_sp", "N_eff_298K", "N_eff_250K",
                 "n_confs_>1%", "lambda_max_nm", "f_avg_298K",
                 "n_rank_inversions"]
    if "dE_vs_ref_kcal" in df_cmp.columns:
        core_cols.insert(4, "dE_vs_ref_kcal")
    core_cols = [c for c in core_cols if c in df_cmp.columns]

    print("\n  ── Resumen por estado ──")
    print(df_cmp[core_cols].to_string(
        index=False, float_format=lambda x: f"{x:.3f}"
    ))

    # Destacar quién tiene mayor peso UV (λ_max más alto)
    if "lambda_max_nm" in df_cmp.columns and df_cmp["lambda_max_nm"].notna().any():
        idx_max = df_cmp["lambda_max_nm"].idxmax()
        mol_max = df_cmp.loc[idx_max, "molecule"]
        lmax    = df_cmp.loc[idx_max, "lambda_max_nm"]
        print(f"\n  Estado con λ_max más alto: {mol_max}  ({lmax:.1f} nm)")

    # N_eff: cuál tiene ensemble más desordenado (más conformacionalmente flexible)
    if "N_eff_298K" in df_cmp.columns and df_cmp["N_eff_298K"].notna().any():
        idx_neff = df_cmp["N_eff_298K"].idxmax()
        mol_neff = df_cmp.loc[idx_neff, "molecule"]
        neff_v   = df_cmp.loc[idx_neff, "N_eff_298K"]
        print(f"  Estado más conformacionalmente flexible (298 K): {mol_neff}  "
              f"(N_eff = {neff_v:.2f})")

    if "dE_vs_ref_kcal" in df_cmp.columns:
        print(f"\n  ΔE relativas (mínimo de cada estado vs {ref_mol}):")
        print("  NOTA: ΔE entre estados de protonación con distinta carga no son")
        print("        termodinámicamente directamente comparables sin correcciones")
        print("        de solvatación absoluta y energía del protón.")

    print(f"\n{sep}\n")


# ─────────────────────────────────────────────
#  ESPECTROS UV-Vis
# ─────────────────────────────────────────────

def _build_spectrum(df_tddft: pd.DataFrame, df_conf: pd.DataFrame,
                    sigma: float, x: np.ndarray,
                    temp_label: str = "298K") -> np.ndarray:
    """Construye el espectro Gaussiano-suavizado Boltzmann-promediado."""
    col = f"boltz_{temp_label}"
    if col not in df_conf.columns:
        col = "boltz_298K"

    y = np.zeros_like(x)
    if df_tddft.empty:
        return y

    df = df_tddft.merge(df_conf[["conformer", col]], on="conformer", how="left")
    df[col] = df[col].fillna(1.0 / max(len(df_tddft), 1))

    for _, row in df.iterrows():
        gauss = np.exp(-0.5 * ((x - row["wavelength_nm"]) / sigma) ** 2)
        y += row["fosc"] * row[col] * gauss
    return y


def plot_spectrum(df_tddft: pd.DataFrame, df_conf: pd.DataFrame, mol_name: str,
                  sigma: float = 20.0, x_range: tuple = (200, 500)):
    """Espectro UV-Vis individual: Boltzmann-promediado a múltiples temperaturas."""
    try:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
    except ImportError:
        print("matplotlib no disponible — saltando plot")
        return

    if df_tddft.empty:
        print("Sin datos TDDFT para graficar.")
        return

    x = np.linspace(x_range[0], x_range[1], 1500)
    colors = ["#0077bb", "#33bbee", "#ee7733", "#cc3311"]

    fig, ax = plt.subplots(figsize=(9, 4.5))

    for lbl, color in zip(TEMP_LABELS, colors):
        y = _build_spectrum(df_tddft, df_conf, sigma, x, temp_label=lbl)
        ax.plot(x, y, lw=1.8, label=f"T = {lbl}", color=color)

    # Palillos de excitaciones individuales (ponderados, 298 K)
    col298 = "boltz_298K"
    df_sticks = df_tddft.merge(df_conf[["conformer", col298]], on="conformer", how="left")
    df_sticks[col298] = df_sticks[col298].fillna(0)
    ax.vlines(df_sticks["wavelength_nm"], 0,
              df_sticks["fosc"] * df_sticks[col298],
              colors="gray", alpha=0.35, lw=0.8)

    ax.set_xlabel("Longitud de onda (nm)", fontsize=12)
    ax.set_ylabel("Intensidad (u.a.)", fontsize=12)
    ax.set_title(f"UV-Vis  —  {mol_name}  (σ = {sigma} nm)", fontsize=13)
    ax.legend(fontsize=10)
    plt.tight_layout()

    out_png = RC_WORK_ROOT / mol_name / "uvvis_spectrum.png"
    plt.savefig(out_png, dpi=150)
    print(f"  Espectro guardado en: {out_png}")
    plt.show()


def plot_comparison(dfs_conf: dict, dfs_tddft: dict,
                    sigma: float = 20.0, x_range: tuple = (200, 500)):
    """
    Superpone espectros UV-Vis (298 K) de todos los estados de protonación,
    normalizados al máximo individual para facilitar comparación de forma.
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib no disponible — saltando plot")
        return

    x = np.linspace(x_range[0], x_range[1], 1500)

    # Paleta colorblind-friendly (Paul Tol)
    palette = [
        "#4477AA", "#EE6677", "#228833", "#CCBB44",
        "#66CCEE", "#AA3377", "#BBBBBB"
    ]

    fig, axes = plt.subplots(2, 1, figsize=(10, 8),
                             gridspec_kw={"height_ratios": [3, 1]})
    ax_spec, ax_bar = axes

    mol_names = sorted(dfs_conf.keys())
    lambda_maxs = []

    for i, mol in enumerate(mol_names):
        df_t = dfs_tddft.get(mol, pd.DataFrame())
        df_c = dfs_conf.get(mol, pd.DataFrame())
        if df_t.empty or df_c.empty:
            continue

        y = _build_spectrum(df_t, df_c, sigma, x, temp_label="298K")
        y_norm = y / y.max() if y.max() > 0 else y

        color = palette[i % len(palette)]
        ax_spec.plot(x, y_norm, lw=2, label=mol.replace("_", " "), color=color)

        # Marcar λ_max de S1 (Boltzmann-promediado, 298 K) — NO el argmax del
        # espectro combinado, que mezcla todos los estados (S1, S2, S3...)
        # ponderados por fosc y puede quedar dominado por un estado superior
        # de mayor fuerza de oscilador. Usamos el mismo valor que va al xlsx
        # (df_conf.attrs["lambda_max_boltz_298K"], calculado solo con S1) para
        # que la figura y la tabla exportada sean consistentes.
        lmax = df_c.attrs.get("lambda_max_boltz_298K")
        if lmax is None:
            # Fallback: recalcular filtrando a S1 si el atributo no está disponible
            s1 = df_t[df_t["state"] == 1].copy()
            if not s1.empty:
                s1 = s1.merge(df_c[["conformer", "boltz_298K"]], on="conformer", how="left")
                w = s1["boltz_298K"].fillna(0)
                lmax = float((s1["wavelength_nm"] * w).sum() / w.sum()) if w.sum() > 0 else None

        if lmax is not None:
            lambda_maxs.append((mol, lmax, color))
            ax_spec.axvline(lmax, color=color, lw=0.7, ls="--", alpha=0.6)

    ax_spec.set_xlabel("Longitud de onda (nm)", fontsize=12)
    ax_spec.set_ylabel("Intensidad normalizada", fontsize=12)
    ax_spec.set_title(f"Comparación UV-Vis — Estados de protonación  "
                      f"(σ = {sigma} nm, T = 298 K; curva = todos los estados, "
                      f"líneas punteadas = λ_max de S1)", fontsize=13)
    ax_spec.legend(fontsize=9, loc="upper right")

    # Gráfico de barras de λ_max
    mols_bar  = [m for m, _, _ in lambda_maxs]
    lmax_bar  = [l for _, l, _ in lambda_maxs]
    cols_bar  = [c for _, _, c in lambda_maxs]
    y_pos = range(len(mols_bar))
    ax_bar.barh(y_pos, lmax_bar, color=cols_bar, alpha=0.85)
    ax_bar.set_yticks(list(y_pos))
    ax_bar.set_yticklabels([m.replace("_", " ") for m in mols_bar], fontsize=9)
    ax_bar.set_xlabel("λ_max (nm)", fontsize=11)
    ax_bar.set_title("λ_max Boltzmann-promediado (S1, 298 K)", fontsize=11)

    plt.tight_layout()
    out_png = RC_WORK_ROOT / "comparison_uvvis.png"
    plt.savefig(out_png, dpi=150)
    print(f"  Espectro comparativo guardado en: {out_png}")
    plt.show()


# ═══════════════════════════════════════════════════════════════════════════════
#  APÉNDICE GOODVIBES — corrección cuasi-armónica Grimme mRRHO
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Fundamento:
#  ORCA calcula G usando el tratamiento rígido-rotor/oscilador-armónico (RRHO)
#  estándar, que sobreestima la contribución entrópica de modos de baja
#  frecuencia (< ~100 cm⁻¹), produciendo G demasiado negativo.
#  GoodVibes aplica la corrección cuasi-RRHO de Grimme (J. Chem. Theory Comput.
#  2012, 8, 4165): por debajo del cutoff f₀ (default 100 cm⁻¹), la entropía
#  vibracional se sustituye por una interpolación entre el oscilador armónico
#  y un rotor libre, suavizando el error sistemático.
#
#  Lo que se compara en este apéndice:
#
#    G_ORCA(i) = E_SP(TZVPD) + G_corr_ORCA(RRHO, SVP/CPCM)
#    G_GV(i)   = E_SP(TZVPD) + G_corr_GV(mRRHO Grimme, SVP/CPCM)
#
#    ΔG_corr(i)   = G_corr_GV - G_corr_ORCA           [Eh]  — diferencia pura
#    ΔΔG_rel(i,T) = ΔG_GV_rel(i,T) - ΔG_ORCA_rel(i,T) [kcal/mol] — impacto en ranking
#    Δboltz(i,T)  = boltz_GV(i,T) - boltz_ORCA(i,T)    — impacto en Boltzmann
#
#  Flujo:
#    run_goodvibes_on_mol()  → invoca `python -m goodvibes` como subprocess,
#                              guarda JSON en runs/<mol>/<conf>/goodvibes/,
#                              retorna dict conf_id → GV_thermochem
#    compute_gv_appendix()  → construye DataFrame de comparación
#    print_gv_appendix()    → tabla formateada al estilo del script
#    export_gv_appendix()   → guarda CSV/XLSX del apéndice en runs/<mol>/
#
# ═══════════════════════════════════════════════════════════════════════════════


# ── Parámetros GoodVibes (modificables) ─────────────────────────────────────
GV_FREQ_CUTOFF = 100.0    # cm⁻¹ — cutoff para corrección mRRHO Grimme
GV_SUBDIR      = "goodvibes"  # subcarpeta creada dentro de cada conf_xxx

# Entorno conda donde está instalado GoodVibes.
# Si se deja en None, el script lo detecta automáticamente buscando en los
# entornos conocidos de miniforge3/miniconda3/anaconda3.


def _rc_resolve_gv_python() -> str:
    """
    Devuelve la ruta absoluta al intérprete Python que tiene goodvibes instalado.

    Estrategia (en orden de prioridad):
      1. Si GV_CONDA_ENV está definido, buscar ese entorno en los prefijos
         estándar de miniforge3, miniconda3 y anaconda3 en /home/kmilozeta.
      2. Intentar `python -m goodvibes --version` con el ejecutable actual
         (sys.executable) por si el entorno ya lo incluye.
      3. Buscar `goodvibes` en el PATH como script de entrada.
      4. Fallback: devolver sys.executable (fallará con el error habitual).
    """
    import shutil

    # ── 1. Entorno conda explícito ────────────────────────────────────────
    if GV_CONDA_ENV:
        home = Path.home()
        for base in ["miniforge3", "miniconda3", "anaconda3", "conda"]:
            candidate = home / base / "envs" / GV_CONDA_ENV / "bin" / "python"
            if candidate.exists():
                # Verificar que goodvibes esté realmente disponible allí
                try:
                    r = subprocess.run(
                        [str(candidate), "-c", "import goodvibes"],
                        capture_output=True, timeout=10
                    )
                    if r.returncode == 0:
                        return str(candidate)
                except Exception:
                    pass

    # ── 2. Python actual ──────────────────────────────────────────────────
    try:
        r = subprocess.run(
            [sys.executable, "-c", "import goodvibes"],
            capture_output=True, timeout=10
        )
        if r.returncode == 0:
            return sys.executable
    except Exception:
        pass

    # ── 3. goodvibes en PATH (instalado como script de consola) ──────────
    gv_script = shutil.which("goodvibes")
    if gv_script:
        # Derivar el Python del mismo directorio bin/
        py_candidate = Path(gv_script).parent / "python"
        if py_candidate.exists():
            return str(py_candidate)

    # ── 4. Fallback ───────────────────────────────────────────────────────
    return sys.executable


def _probe_gv_capabilities(python_exe: str) -> dict:
    """
    Ejecuta `python -m goodvibes --help` para detectar qué flags están
    disponibles en la versión instalada.

    Retorna un dict con:
        supports_json_format  : bool — existe --output_format json (v4+)
        supports_csv          : bool — existe --csv (v3.x)
        version               : str  — número de versión o "?"
    """
    caps = {"supports_json_format": False, "supports_csv": False, "version": "?"}
    try:
        r = subprocess.run(
            [python_exe, "-m", "goodvibes", "--help"],
            capture_output=True, text=True, timeout=15
        )
        help_text = r.stdout + r.stderr

        # Detectar flags de salida disponibles
        if "--output_format" in help_text or "output_format" in help_text:
            caps["supports_json_format"] = True
        if "--csv" in help_text:
            caps["supports_csv"] = True

        # Intentar versión
        r2 = subprocess.run(
            [python_exe, "-m", "goodvibes", "--version"],
            capture_output=True, text=True, timeout=10
        )
        for line in (r2.stdout + r2.stderr).splitlines():
            m = re.search(r"([\d]+\.[\d]+\.[\d]+|v[\d.]+)", line)
            if m:
                caps["version"] = m.group(1).lstrip("v")
                break
    except Exception:
        pass
    return caps


def _parse_goodvibes_dat(dat_path: Path) -> dict:
    """
    Parsea el archivo de texto tabulado que GoodVibes genera siempre
    (Goodvibes_output.dat o similar) cuando no se solicita salida JSON.

    GoodVibes imprime una tabla con el formato (columnas separadas por espacios):
        Structure  E_SPC  ZPE  H  T.S  T.qh-S  G(T)  qh-G(T)
    o bien (versión 3.x, con nombres ligeramente distintos):
        Structure  E  ZPE  H  T.S  T.qh-S  G(T)  qh-G(T)

    Las unidades son Hartree para energías y Hartree/K para S si se reporta.
    Las líneas de datos comienzan con el nombre del archivo (sin extensión).

    Retorna el mismo formato que _parse_goodvibes_json():
        H_GV, S_GV (Eh/K), TS_GV (Eh), G_GV (Eh), ZPE_GV (Eh),
        E_elec_GV (Eh), T_GV (K), G_corr_GV (Eh), ok, warning
    """
    null = {
        "H_GV": None, "S_GV": None, "TS_GV": None, "G_GV": None,
        "ZPE_GV": None, "E_elec_GV": None, "T_GV": None,
        "G_corr_GV": None, "ok": False, "warning": "dat no parseado",
    }

    if not dat_path.exists():
        null["warning"] = f"Archivo .dat no encontrado: {dat_path}"
        return null

    try:
        text = dat_path.read_text(errors="replace")
    except OSError as exc:
        null["warning"] = f"No se pudo leer {dat_path}: {exc}"
        return null

    lines = text.splitlines()

    # ── Detectar temperatura del encabezado ──────────────────────────────
    T_GV: Optional[float] = None
    for line in lines:
        m = re.search(r"Temp\s*=\s*([\d.]+)\s*K", line, re.I)
        if not m:
            m = re.search(r"at\s+([\d.]+)\s*K", line, re.I)
        if m:
            try:
                T_GV = float(m.group(1))
                break
            except ValueError:
                pass
    if T_GV is None:
        T_GV = 298.15  # asumir condición estándar si no se detecta

    # ── Detectar línea de encabezado de columnas ─────────────────────────
    # GoodVibes imprime una línea con "Structure" y los nombres de columnas.
    # Buscamos esa línea y la siguiente que sea numérica (datos).
    header_idx = -1
    for i, line in enumerate(lines):
        if re.search(r"\bStructure\b", line, re.I):
            header_idx = i
            break

    if header_idx < 0:
        null["warning"] = f"No se encontró encabezado 'Structure' en {dat_path}"
        return null

    header = lines[header_idx].split()

    # Normalizar nombres de columna a minúsculas para búsqueda robusta
    def _col(names: list[str]) -> Optional[int]:
        """Devuelve el índice (0-based) de la primera columna que coincida."""
        for name in names:
            nl = name.lower()
            for j, h in enumerate(header):
                if h.lower() == nl:
                    return j
        return None

    # Columnas posibles según versión de GoodVibes
    idx_e    = _col(["E_SPC", "E_SP", "E", "E_DFT"])
    idx_zpe  = _col(["ZPE"])
    idx_h    = _col(["H"])
    idx_ts   = _col(["T.S", "TS", "T*S"])       # T·S a T_GV (Eh)
    idx_tqs  = _col(["T.qh-S", "T.QH-S"])       # T·qhS (corrección GoodVibes)
    idx_g    = _col(["G(T)", "G"])               # G RRHO (sin corrección quasi)
    idx_qhg  = _col(["qh-G(T)", "qh-G", "QH-G(T)"])  # G quasi-RRHO (el que queremos)

    # Buscar primera línea de datos después del encabezado
    # NOTA: GoodVibes v4 antepone un símbolo de estado (p.ej. "o", "x", "?")
    # ANTES del nombre de la estructura, produciendo una columna extra al inicio
    # que NO existe en el encabezado. Ejemplo real:
    #   encabezado : Structure   E   ZPE   H   qh-H   T.S   T.qh-S   G(T)   qh-G(T)
    #   datos      : o   freq   -1220.967288   0.403359   ...
    # Los índices derivados del encabezado (base 0) deben desplazarse +1 para
    # las líneas de datos cuando se detecta este símbolo inicial.
    E_elec: Optional[float] = None
    ZPE:    Optional[float] = None
    H_GV:   Optional[float] = None
    TS_GV:  Optional[float] = None   # T·qhS (Eh) — con corrección cuasi-RRHO
    G_GV:   Optional[float] = None   # qh-G(T) absoluto (Eh)

    for line in lines[header_idx + 1:]:
        parts = line.split()
        if len(parts) < 3:
            continue

        # Detectar offset: si parts[0] es un símbolo de 1 carácter (o/x/?)
        # y parts[1] no es numérico (es el nombre del archivo), hay offset=1.
        # Si parts[0] no es numérico directamente (es el nombre), offset=0.
        offset = 0
        if len(parts[0]) == 1 and not parts[0].lstrip('-').replace('.','',1).isdigit():
            # parts[0] es símbolo de estado → parts[1] es nombre, datos desde parts[2]
            try:
                float(parts[2])
                offset = 1
            except (ValueError, IndexError):
                continue
        else:
            # parts[0] es nombre → datos desde parts[1]
            try:
                float(parts[1])
                offset = 0
            except (ValueError, IndexError):
                continue

        def _get_col(idx: Optional[int]) -> Optional[float]:
            if idx is None:
                return None
            real_idx = idx + offset
            if real_idx >= len(parts):
                return None
            try:
                return float(parts[real_idx])
            except ValueError:
                return None

        E_elec = _get_col(idx_e)
        ZPE    = _get_col(idx_zpe)
        H_GV   = _get_col(idx_h)
        TS_GV  = _get_col(idx_tqs) or _get_col(idx_ts)  # preferir quasi-RRHO
        G_GV   = _get_col(idx_qhg) or _get_col(idx_g)   # preferir quasi-RRHO
        break  # solo necesitamos la primera fila de datos

    if G_GV is None:
        null["warning"] = (
            f"No se pudo extraer G(T) del .dat ({dat_path}). "
            f"Columnas detectadas: {header}"
        )
        return null

    S_GV: Optional[float] = None
    if TS_GV is not None and T_GV > 0:
        S_GV = TS_GV / T_GV  # Eh/K

    G_corr_GV: Optional[float] = None
    if G_GV is not None and E_elec is not None:
        G_corr_GV = G_GV - E_elec

    return {
        "H_GV":      H_GV,
        "S_GV":      S_GV,        # Eh/K
        "TS_GV":     TS_GV,       # Eh (T·qhS)
        "G_GV":      G_GV,        # Eh absoluto (qh-G(T))
        "ZPE_GV":    ZPE,         # Eh
        "E_elec_GV": E_elec,      # Eh
        "T_GV":      T_GV,        # K
        "G_corr_GV": G_corr_GV,   # Eh (puede ser None si no hay E_elec)
        "ok":        True,
        "warning":   "",
    }


def run_goodvibes_on_mol(mol_name: str) -> dict[str, dict]:
    """
    Invoca GoodVibes sobre todos los freq.out de 03_numfreq/<conf_*>/freq.out
    para la molécula mol_name.

    Flujo para cada confórmero:
      1. Resuelve el Python correcto (entorno sciencebitch o auto-detección).
      2. Detecta las capacidades de la versión instalada (--output_format json
         vs parseo del .dat de salida clásico).
      3. Crea runs/<mol>/<conf>/goodvibes/ y ejecuta GoodVibes desde allí.
      4. Parsea JSON (v4+) o el .dat tabulado (v3.x) según disponibilidad.

    Retorna dict:
        { "conf_001": {"H_GV": float, "S_GV": float, "G_GV": float,
                       "ZPE_GV": float, "T_GV": float, "ok": bool,
                       "json_path": Path, "warning": str}, ... }

    Si GoodVibes no está disponible o falla para un confórmero, "ok" = False
    y "warning" describe el error — el resto del script no se interrumpe.
    """
    mol_dir  = RC_WORK_ROOT / mol_name
    freq_dir = mol_dir / "03_numfreq"

    if not freq_dir.exists():
        warnings.warn(f"[GoodVibes] {freq_dir} no existe — omitiendo {mol_name}.")
        return {}

    conf_dirs = sorted(freq_dir.glob("conf_*"))
    if not conf_dirs:
        warnings.warn(f"[GoodVibes] Sin subdirectorios conf_* en {freq_dir}.")
        return {}

    # ── Resolver Python y capacidades UNA SOLA VEZ para todos los conf ────
    python_exe = _rc_resolve_gv_python()
    caps       = _probe_gv_capabilities(python_exe)

    # Verificar que goodvibes sea invocable con el ejecutable resuelto
    try:
        probe = subprocess.run(
            [python_exe, "-c", "import goodvibes"],
            capture_output=True, timeout=10
        )
        if probe.returncode != 0:
            warnings.warn(
                f"[GoodVibes] No se puede importar goodvibes con {python_exe}.\n"
                f"  Instalar con: conda activate sciencebitch && pip install goodvibes\n"
                f"  stderr: {probe.stderr.strip()[:200]}"
            )
            return {c.name: {
                "H_GV": None, "S_GV": None, "G_GV": None,
                "ZPE_GV": None, "T_GV": None, "ok": False,
                "json_path": None,
                "warning": f"goodvibes no disponible en {python_exe}",
            } for c in conf_dirs}
    except Exception as exc:
        warnings.warn(f"[GoodVibes] Error al verificar goodvibes: {exc}")
        return {}

    print(f"  [GoodVibes v{caps['version']}] Python: {python_exe}")
    use_json = caps["supports_json_format"]
    print(f"  [GoodVibes] Modo salida: {'JSON (--output_format json)' if use_json else '.dat tabulado (clásico)'}")

    results: dict[str, dict] = {}

    for conf_dir in conf_dirs:
        conf_id  = conf_dir.name
        freq_out = conf_dir / "freq.out"

        null_entry = {
            "H_GV":    None, "S_GV":   None, "G_GV":  None,
            "ZPE_GV":  None, "T_GV":   None,
            "ok":      False, "json_path": None, "warning": "",
        }

        if not freq_out.exists():
            null_entry["warning"] = f"freq.out no encontrado en {conf_dir}"
            results[conf_id] = null_entry
            continue

        # Crear subcarpeta goodvibes/ dentro del directorio del confórmero
        gv_dir = conf_dir / GV_SUBDIR
        gv_dir.mkdir(exist_ok=True)

        # ── Construir comando según capacidades de la versión ─────────────
        # GoodVibes siempre necesita la ruta al .out; corremos desde gv_dir
        # para que los archivos de salida queden ahí con nombres predecibles.
        cmd_base = [
            python_exe, "-m", "goodvibes",
            str(freq_out.resolve()),
            "-q",
            "--qs", GV_QS,
            "-f", str(GV_FREQ_CUTOFF),
            "-t", "298.15",
        ]

        if use_json:
            # v4+: --output_format json → escribe <stem>_goodvibes.json en cwd
            json_path = gv_dir / f"{freq_out.stem}_goodvibes.json"
            cmd = cmd_base + ["--output_format", "json"]
        else:
            # v3.x: sin flag de formato → escribe Goodvibes_output.dat en cwd
            json_path = None
            cmd = cmd_base  # salida .dat clásica

        try:
            proc = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=120,
                cwd=str(gv_dir),
            )
        except subprocess.TimeoutExpired:
            null_entry["warning"] = "GoodVibes timeout (>120 s)"
            results[conf_id] = null_entry
            continue
        except FileNotFoundError:
            null_entry["warning"] = (
                f"Ejecutable no encontrado: {python_exe}\n"
                "  Instalar con: conda activate sciencebitch && pip install goodvibes"
            )
            results[conf_id] = null_entry
            continue

        if proc.returncode != 0:
            null_entry["warning"] = (
                f"GoodVibes retornó código {proc.returncode}. "
                f"stderr: {proc.stderr[:300]}"
            )
            results[conf_id] = null_entry
            continue

        # ── Parsear la salida según el formato producido ───────────────────
        if use_json and json_path is not None:
            if not json_path.exists():
                # v4+ con --output_format json pero JSON no generado →
                # intentar fallback con el .dat que puede haber escrito de todas formas
                dat_candidates = list(gv_dir.glob("*.dat")) + list(gv_dir.glob("Goodvibes_output*"))
                if dat_candidates:
                    null_entry["warning"] = (
                        f"JSON no generado por GoodVibes v4 en {json_path}; "
                        f"usando .dat de fallback."
                    )
                    gv_data = _parse_goodvibes_dat(dat_candidates[0])
                    if gv_data["ok"]:
                        gv_data["json_path"] = dat_candidates[0]
                        results[conf_id] = gv_data
                    else:
                        null_entry["warning"] += f" {gv_data['warning']}"
                        results[conf_id] = null_entry
                else:
                    null_entry["warning"] = (
                        f"JSON no generado por GoodVibes en {json_path}. "
                        f"stdout: {proc.stdout[:300]}"
                    )
                    results[conf_id] = null_entry
                continue

            try:
                gv_data = _parse_goodvibes_json(json_path, freq_out.name)
            except Exception as exc:
                null_entry["warning"] = f"Error parseando JSON de GoodVibes: {exc}"
                results[conf_id] = null_entry
                continue

            gv_data["json_path"] = json_path
            results[conf_id] = gv_data

        else:
            # Parsear .dat tabulado (v3.x o cuando JSON no está disponible)
            dat_candidates = (
                list(gv_dir.glob("Goodvibes_output*"))
                + list(gv_dir.glob("*.dat"))
            )
            if not dat_candidates:
                null_entry["warning"] = (
                    f"GoodVibes no generó ningún archivo de salida en {gv_dir}. "
                    f"stdout: {proc.stdout[:400]}"
                )
                results[conf_id] = null_entry
                continue

            dat_path = dat_candidates[0]
            try:
                gv_data = _parse_goodvibes_dat(dat_path)
            except Exception as exc:
                null_entry["warning"] = f"Error parseando .dat de GoodVibes: {exc}"
                results[conf_id] = null_entry
                continue

            if not gv_data["ok"]:
                null_entry["warning"] = gv_data.get("warning", "Parseo .dat falló")
                results[conf_id] = null_entry
                continue

            gv_data["json_path"] = dat_path  # reaprovechamos el campo para la ruta
            results[conf_id] = gv_data

    return results


def _parse_goodvibes_json(json_path: Path, filename: str) -> dict:
    """
    Parsea el JSON producido por GoodVibes v4.x y extrae H, S, G (mRRHO), ZPE.

    El JSON de GoodVibes v4 tiene la forma:
        {
          "timestamp": "...",
          "version": "4.x",
          "molecules": [
            {
              "filename": "freq.out",
              "t":  298.15,
              "zpe": <float Eh>,
              "h":  <float Eh>,    ← entalpía total (E_elec + H_corr_GV)
              "ts": <float Eh>,    ← T·S (en Eh)
              "g":  <float Eh>,    ← G total (E_elec + G_corr_GV)
              ...
            }
          ]
        }

    Nota: GoodVibes reporta H y G como energías absolutas que incluyen E_elec.
    Para obtener G_corr_GV (compatible con el ciclo compuesto del script):
        G_corr_GV = G_GV - E_elec_GV
    donde E_elec_GV es la energía electrónica que GoodVibes lee del .out
    (campo "sp" o "e_sp" en el JSON, según versión).
    """
    with open(json_path) as f:
        data = json.load(f)

    # GoodVibes v4 agrupa bajo "molecules" (lista) o directamente como dict
    molecules = data.get("molecules") or data.get("data") or []

    # Buscar la entrada que corresponde a nuestro freq.out
    entry = None
    stem  = Path(filename).stem  # "freq"

    for mol in molecules:
        fn = mol.get("filename", mol.get("file", ""))
        if Path(fn).stem == stem or fn == filename:
            entry = mol
            break

    # Fallback: si solo hay una entrada, usarla directamente
    if entry is None and len(molecules) == 1:
        entry = molecules[0]

    if entry is None:
        raise ValueError(
            f"No se encontró entrada para '{filename}' en {json_path}. "
            f"Claves disponibles: {[m.get('filename','?') for m in molecules]}"
        )

    # Extraer campos — GoodVibes v4 usa minúsculas
    def _get(*keys):
        for k in keys:
            v = entry.get(k)
            if v is not None:
                try:
                    return float(v)
                except (TypeError, ValueError):
                    pass
        return None

    T_GV   = _get("t", "temp", "temperature")
    ZPE_GV = _get("zpe", "ZPE")
    H_GV   = _get("h", "H", "enthalpy")           # absoluto (Eh)
    TS_GV  = _get("ts", "tds", "t*s")             # T·S en Eh
    G_GV   = _get("g", "G", "gibbs", "qh_g")      # absoluto (Eh), con mRRHO
    E_elec = _get("sp", "e_sp", "e_elec",
                  "electronic_energy", "scf_energy")

    # Algunas versiones reportan S directamente (en Eh/K), otras T·S
    S_GV = None
    if TS_GV is not None and T_GV is not None and T_GV > 0:
        S_GV = TS_GV / T_GV   # Eh/K
    else:
        S_GV = _get("s", "S", "entropy")

    # G_corr_GV = G_total_GV - E_elec
    # Si GoodVibes no reporta E_elec en el JSON, lo dejamos None;
    # el apéndice usará G_GV directamente para ΔG_rel (diferencias).
    G_corr_GV = None
    if G_GV is not None and E_elec is not None:
        G_corr_GV = G_GV - E_elec

    ok = G_GV is not None

    return {
        "H_GV":      H_GV,
        "S_GV":      S_GV,       # Eh/K
        "TS_GV":     TS_GV,      # Eh
        "G_GV":      G_GV,       # Eh absoluto (E_elec + G_corr_GV)
        "G_corr_GV": G_corr_GV,  # Eh — corrección pura (None si no hay E_elec)
        "ZPE_GV":    ZPE_GV,     # Eh
        "E_elec_GV": E_elec,     # Eh — energía electrónica leída por GoodVibes
        "T_GV":      T_GV,       # K
        "ok":        ok,
        "warning":   "" if ok else "G_GV no encontrado en JSON",
    }


def compute_gv_appendix(
    mol_name: str,
    df_conf: pd.DataFrame,
    gv_results: dict[str, dict],
) -> pd.DataFrame:
    """
    Construye el DataFrame del apéndice GoodVibes, comparando G_ORCA vs G_GV
    para cada confórmero a todas las temperaturas.

    Requiere que df_conf tenga las columnas producidas por compute_G_total()
    (G_total_{label}, dG_{label}_kcal, boltz_{label}, H_corr, S_total_EhK).
    Si compute_G_total() no fue llamado previamente (el flujo normal de
    collect_results no lo llama), este módulo lo invoca internamente.

    Columnas del DataFrame resultante:
        conformer          — id del confórmero
        E_SP_hartree       — energía SP ORCA (def2-TZVPD, Eh)
        — ORCA (RRHO) —
        G_corr_ORCA_Eh     — corrección G de ORCA a 298 K (Eh)
        H_corr_ORCA_Eh     — corrección H de ORCA (Eh)
        S_ORCA_EhK         — entropía ORCA (Eh/K)
        G_total_ORCA_Eh    — G compuesto ORCA a 298 K (Eh)
        dG_ORCA_kcal       — ΔG relativo ORCA a 298 K (kcal/mol)
        boltz_ORCA_298K    — peso Boltzmann ORCA a 298 K
        — GoodVibes (mRRHO Grimme) —
        G_GV_abs_Eh        — G absoluto GoodVibes (Eh, incluye E_elec)
        G_corr_GV_Eh       — corrección G de GV (Eh, None si sin E_elec en JSON)
        S_GV_EhK           — entropía quasi-RRHO GV (Eh/K)
        ZPE_GV_Eh          — ZPE de GoodVibes (Eh)
        G_total_GV_Eh      — G compuesto GV a 298 K = E_SP + G_corr_GV (Eh)
        dG_GV_kcal         — ΔG relativo GV a 298 K (kcal/mol)
        boltz_GV_298K      — peso Boltzmann GV a 298 K
        — Diferencias —
        dG_corr_Eh         — G_corr_GV - G_corr_ORCA (Eh)
        dG_corr_kcal       — ídem en kcal/mol
        dS_EhK             — S_GV - S_ORCA (Eh/K)  ← origen del Δ
        ddG_rel_kcal       — ΔG_GV_rel - ΔG_ORCA_rel a 298 K (kcal/mol)
        d_boltz_298K       — boltz_GV - boltz_ORCA a 298 K
        rank_ORCA          — ranking conformacional ORCA (1=más estable)
        rank_GV            — ranking conformacional GV
        rank_change_GV     — rank_GV - rank_ORCA (0 = sin inversión)
        gv_ok              — True si GoodVibes procesó este confórmero
        gv_warning         — mensaje de advertencia si gv_ok=False
    """
    mol_dir = RC_WORK_ROOT / mol_name

    # ── Asegurar correcciones térmicas ORCA en df_conf ────────────────────
    # Si compute_G_total ya fue llamado, G_total_298K existirá.
    # Si no (flujo normal), lo recalculamos internamente aquí.
    if "G_total_298K" not in df_conf.columns or "H_corr" not in df_conf.columns:
        # Cargar termoquímica ORCA para cada confórmero
        df_work = df_conf.copy()
        _thermochem_cache: list[Optional[dict]] = []
        for _, row in df_work.iterrows():
            conf_id  = row["conformer"]
            freq_out = mol_dir / "03_numfreq" / conf_id / "freq.out"
            tc = rc_parse_thermochem(freq_out, T=T_STD)
            _thermochem_cache.append(tc if tc.get("ok") else None)
        df_work["_thermochem"] = _thermochem_cache
        df_work = compute_G_total(df_work)
        df_work = compute_boltzmann_multiT(df_work)
    else:
        df_work = df_conf.copy()

    rows: list[dict] = []

    for _, row in df_work.iterrows():
        conf_id  = row["conformer"]
        gv       = gv_results.get(conf_id, {})
        gv_ok    = gv.get("ok", False)

        e_sp = row.get("E_hartree")

        # ── ORCA RRHO ──────────────────────────────────────────────────────
        g_corr_orca = row.get("G_corr_298K")     # Eh
        h_corr_orca = row.get("H_corr")           # Eh
        s_orca      = row.get("S_total_EhK")      # Eh/K
        g_total_orca = row.get("G_total_298K")    # Eh
        dg_orca      = row.get("dG_298K_kcal")    # kcal/mol
        boltz_orca   = row.get("boltz_298K")

        # ── GoodVibes mRRHO ────────────────────────────────────────────────
        g_gv_abs   = gv.get("G_GV")              # Eh absoluto
        g_corr_gv  = gv.get("G_corr_GV")         # Eh (puede ser None)
        s_gv       = gv.get("S_GV")              # Eh/K
        zpe_gv     = gv.get("ZPE_GV")            # Eh
        e_elec_gv  = gv.get("E_elec_GV")         # Eh

        # G compuesto GV = E_SP(TZVPD) + G_corr_GV
        # Si GoodVibes no reporta E_elec separado, usamos G_GV directamente
        # como proxy (válido para diferencias relativas entre confórmeros).
        g_total_gv: Optional[float] = None
        if g_corr_gv is not None and e_sp is not None:
            g_total_gv = e_sp + g_corr_gv          # ciclo compuesto correcto
        elif g_gv_abs is not None and e_elec_gv is not None and e_sp is not None:
            # Reconstruir G_corr a partir del absoluto
            g_corr_gv  = g_gv_abs - e_elec_gv
            g_total_gv = e_sp + g_corr_gv
        elif g_gv_abs is not None:
            # Último recurso: usar G absoluto de GV para los ΔG (pérdida del SP TZVPD)
            g_total_gv = g_gv_abs

        # ── Diferencias ────────────────────────────────────────────────────
        dg_corr_eh   = None
        dg_corr_kcal = None
        ds_ehk       = None

        if g_corr_gv is not None and g_corr_orca is not None:
            dg_corr_eh   = g_corr_gv - g_corr_orca
            dg_corr_kcal = dg_corr_eh * HARTREE_TO_KCAL
        if s_gv is not None and s_orca is not None:
            ds_ehk = s_gv - s_orca

        rows.append({
            "conformer":        conf_id,
            "E_SP_hartree":     e_sp,
            # ORCA
            "G_corr_ORCA_Eh":   g_corr_orca,
            "H_corr_ORCA_Eh":   h_corr_orca,
            "S_ORCA_EhK":       s_orca,
            "G_total_ORCA_Eh":  g_total_orca,
            "dG_ORCA_kcal":     dg_orca,
            "boltz_ORCA_298K":  boltz_orca,
            # GoodVibes
            "G_GV_abs_Eh":      g_gv_abs,
            "G_corr_GV_Eh":     g_corr_gv,
            "S_GV_EhK":         s_gv,
            "ZPE_GV_Eh":        zpe_gv,
            "G_total_GV_Eh":    g_total_gv,
            # (dG_GV_kcal y boltz_GV_298K se calculan abajo tras tener todos los g_total)
            "gv_ok":            gv_ok,
            "gv_warning":       gv.get("warning", ""),
            # Diferencias parciales (las que no dependen del mínimo del ensemble)
            "dG_corr_Eh":       dg_corr_eh,
            "dG_corr_kcal":     dg_corr_kcal,
            "dS_EhK":           ds_ehk,
        })

    df_gv = pd.DataFrame(rows)

    if df_gv.empty:
        return df_gv

    # ── ΔG relativo y Boltzmann para GV (necesita el mínimo del ensemble) ──
    valid_gv = df_gv["G_total_GV_Eh"].notna() & df_gv["gv_ok"]

    if valid_gv.any():
        g_min_gv = df_gv.loc[valid_gv, "G_total_GV_Eh"].min()
        df_gv["dG_GV_kcal"] = np.nan
        df_gv.loc[valid_gv, "dG_GV_kcal"] = (
            (df_gv.loc[valid_gv, "G_total_GV_Eh"] - g_min_gv) * HARTREE_TO_KCAL
        )

        # Boltzmann GV a 298 K
        dg_gv = df_gv.loc[valid_gv, "dG_GV_kcal"]
        raw   = np.exp(-dg_gv / (KB_KCAL * T_STD))
        df_gv["boltz_GV_298K"] = np.nan
        df_gv.loc[valid_gv, "boltz_GV_298K"] = raw / raw.sum()
    else:
        df_gv["dG_GV_kcal"]    = np.nan
        df_gv["boltz_GV_298K"] = np.nan

    # ── ΔΔG relativo y Δboltz ─────────────────────────────────────────────
    df_gv["ddG_rel_kcal"] = df_gv["dG_GV_kcal"] - df_gv["dG_ORCA_kcal"]
    df_gv["d_boltz_298K"] = df_gv["boltz_GV_298K"] - df_gv["boltz_ORCA_298K"]

    # ── Rankings ──────────────────────────────────────────────────────────
    valid_orca = df_gv["G_total_ORCA_Eh"].notna()
    if valid_orca.any():
        df_gv.loc[valid_orca, "rank_ORCA"] = (
            df_gv.loc[valid_orca, "G_total_ORCA_Eh"]
            .rank(method="min").astype("Int64")
        )
    if valid_gv.any():
        df_gv.loc[valid_gv, "rank_GV"] = (
            df_gv.loc[valid_gv, "G_total_GV_Eh"]
            .rank(method="min").astype("Int64")
        )
    if "rank_ORCA" in df_gv.columns and "rank_GV" in df_gv.columns:
        df_gv["rank_change_GV"] = df_gv["rank_GV"] - df_gv["rank_ORCA"]

    return df_gv


def print_gv_appendix(df_gv: pd.DataFrame, mol_name: str):
    """
    Imprime el apéndice GoodVibes con el mismo estilo que el resto del script.
    Muestra tres secciones:
      1. Correcciones absolutas ORCA vs GV (G_corr, S, ZPE)
      2. Impacto en ΔG relativo y pesos Boltzmann a 298 K
      3. Cambios de ranking conformacional (si los hay)
    """
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  APÉNDICE GOODVIBES — Corrección cuasi-armónica mRRHO (Grimme)")
    print(f"  Molécula : {mol_name}")
    print(f"  Método   : -q --qs grimme -f {GV_FREQ_CUTOFF:.0f} cm⁻¹ (v{_gv_version()})")
    print(f"  Ciclo    : G_GV = E_SP(TZVPD/CPCM) + G_corr_GV(mRRHO, SVP/CPCM)")
    print(sep)

    if df_gv.empty:
        print("  Sin datos de GoodVibes.")
        return

    # Confórmeros sin GV
    failed = df_gv[~df_gv["gv_ok"]]
    if not failed.empty:
        print(f"\n  ⚠  {len(failed)} confórmero(s) sin datos GoodVibes:")
        for _, r in failed.iterrows():
            print(f"     {r['conformer']:>10}  →  {r['gv_warning']}")

    ok = df_gv[df_gv["gv_ok"]].copy()
    if ok.empty:
        print("  Sin confórmeros válidos de GoodVibes para comparar.")
        return

    # ── Sección 1: Correcciones absolutas ────────────────────────────────
    print(f"\n  ── 1. Correcciones absolutas (Eh) ── "
          f"[ORCA=RRHO  |  GV=mRRHO Grimme, f₀={GV_FREQ_CUTOFF:.0f} cm⁻¹]")
    print(f"\n  {'conf':>10}  {'G_corr_ORCA':>14}  {'G_corr_GV':>12}  "
          f"{'ΔG_corr(Eh)':>13}  {'ΔG_corr(kcal)':>14}  "
          f"{'S_ORCA(Eh/K)':>13}  {'S_GV(Eh/K)':>11}  {'ΔS(Eh/K)':>10}")
    print("  " + "─" * 100)

    for _, r in ok.iterrows():
        def _f(v, fmt=".6f"):
            return f"{v:{fmt}}" if pd.notna(v) else "      —     "

        dg_sign = "+" if (r.get("dG_corr_Eh") or 0) >= 0 else ""
        ds_sign = "+" if (r.get("dS_EhK") or 0) >= 0 else ""
        print(
            f"  {r['conformer']:>10}  "
            f"{_f(r.get('G_corr_ORCA_Eh')):>14}  "
            f"{_f(r.get('G_corr_GV_Eh')):>12}  "
            f"{dg_sign}{_f(r.get('dG_corr_Eh'), '.6f'):>12}  "
            f"{dg_sign}{_f(r.get('dG_corr_kcal'), '.4f'):>13}  "
            f"{_f(r.get('S_ORCA_EhK')):>13}  "
            f"{_f(r.get('S_GV_EhK')):>11}  "
            f"{ds_sign}{_f(r.get('dS_EhK'), '.6f'):>9}"
        )

    # Estadística del Δ(G_corr)
    dg_vals = ok["dG_corr_kcal"].dropna()
    if not dg_vals.empty:
        print(f"\n  Δ(G_corr) medio  : {dg_vals.mean():+.4f} kcal/mol")
        print(f"  Δ(G_corr) rango  : [{dg_vals.min():+.4f}, {dg_vals.max():+.4f}] kcal/mol")
        print(f"  Interpretación   : GV{'aumenta' if dg_vals.mean()>0 else 'reduce'} G_corr "
              f"{'(mRRHO menos negativo → G más alto)' if dg_vals.mean()>0 else '(mRRHO más negativo → G más bajo)'}")

    # ── Sección 2: Impacto en ΔG relativo y Boltzmann ────────────────────
    print(f"\n  ── 2. Impacto en ΔG relativo y pesos Boltzmann (T = 298.15 K) ──")
    print(f"\n  {'conf':>10}  {'ΔG_ORCA':>9}  {'ΔG_GV':>8}  {'ΔΔG_rel':>9}  "
          f"{'boltz_ORCA':>11}  {'boltz_GV':>10}  {'Δboltz':>8}  {'dominante':>9}")
    print("  " + "─" * 85)

    for _, r in ok.sort_values("dG_GV_kcal").iterrows():
        dg_o   = r.get("dG_ORCA_kcal")
        dg_gv  = r.get("dG_GV_kcal")
        ddg    = r.get("ddG_rel_kcal")
        bo     = r.get("boltz_ORCA_298K")
        bg     = r.get("boltz_GV_298K")
        db     = r.get("d_boltz_298K")

        def _fv(v, fmt):
            return f"{v:{fmt}}" if pd.notna(v) else "   —   "

        dom_gv = " ★" if pd.notna(bg) and bg > 0.50 else ("  ·" if pd.notna(bg) and bg > 0.10 else "   ")
        print(
            f"  {r['conformer']:>10}  "
            f"{_fv(dg_o,  '+9.3f')}  "
            f"{_fv(dg_gv, '+8.3f')}  "
            f"{_fv(ddg,   '+9.4f')}  "
            f"{_fv(bo,    '11.4f')}  "
            f"{_fv(bg,    '10.4f')}  "
            f"{_fv(db,    '+8.4f')}"
            f"{dom_gv}"
        )

    # Resumen del impacto sobre Boltzmann
    ddg_vals = ok["ddG_rel_kcal"].dropna()
    db_vals  = ok["d_boltz_298K"].dropna()
    if not ddg_vals.empty:
        print(f"\n  |ΔΔG_rel| máximo : {ddg_vals.abs().max():.4f} kcal/mol")
        if ddg_vals.abs().max() < 0.1:
            print("  → Impacto negligible en los ΔG relativos (< 0.1 kcal/mol)")
        elif ddg_vals.abs().max() < 0.5:
            print("  → Impacto pequeño pero no nulo en los ΔG relativos (0.1–0.5 kcal/mol)")
        else:
            print("  → Impacto significativo en los ΔG relativos (> 0.5 kcal/mol) ⚠")
    if not db_vals.empty:
        print(f"  |Δboltz| máximo  : {db_vals.abs().max():.4f} "
              f"({'cambio relevante en ensemble' if db_vals.abs().max() > 0.05 else 'cambio menor en ensemble'})")

    # ── Sección 3: Rankings ───────────────────────────────────────────────
    if "rank_change_GV" in df_gv.columns:
        inversions = df_gv["rank_change_GV"].dropna()
        n_inv = (inversions.abs() > 0).sum()
        print(f"\n  ── 3. Ranking conformacional: ORCA vs GoodVibes ──")
        print(f"  Inversiones detectadas: {n_inv}/{len(ok)}")

        if n_inv > 0:
            print(f"\n  {'conf':>10}  {'rank_ORCA':>10}  {'rank_GV':>8}  {'Δrank':>7}")
            print("  " + "─" * 42)
            for _, r in df_gv[df_gv["rank_change_GV"].abs() > 0].iterrows():
                rch = r.get("rank_change_GV", 0)
                flag = " ◀ inversión"
                print(
                    f"  {r['conformer']:>10}  "
                    f"{str(r.get('rank_ORCA','—')):>10}  "
                    f"{str(r.get('rank_GV','—')):>8}  "
                    f"{int(rch):>+7}{flag}"
                )
        else:
            print("  → El ranking conformacional es idéntico con ORCA-RRHO y GV-mRRHO.")

    print(f"\n  Archivos GV guardados en:")
    for conf_id in df_gv["conformer"]:
        p = RC_WORK_ROOT / mol_name / "03_numfreq" / conf_id / GV_SUBDIR
        if p.exists():
            # El nombre del JSON depende del stem del .out (normalmente freq_goodvibes.json)
            json_files = list(p.glob("*_goodvibes.json"))
            fname = json_files[0].name if json_files else "*.json"
            print(f"    {p}/{fname}")

    print(f"\n{sep}\n")


def _gv_version() -> str:
    """Devuelve la versión de GoodVibes instalada, o '?' si no se puede determinar."""
    try:
        python_exe = _rc_resolve_gv_python()
        result = subprocess.run(
            [python_exe, "-m", "goodvibes", "--version"],
            capture_output=True, text=True, timeout=10
        )
        for line in (result.stdout + result.stderr).splitlines():
            m = re.search(r"([\d]+\.[\d]+\.[\d]+|v[\d.]+)", line)
            if m:
                return m.group(1).lstrip("v")
    except Exception:
        pass
    return "?"


def export_gv_appendix(df_gv: pd.DataFrame, mol_name: str, fmt: str):
    """
    Guarda el DataFrame del apéndice GoodVibes en runs/<mol>/.
    Formatos: "csv", "xlsx" o "both".
    Si fmt incluye xlsx, lo agrega como hoja extra del Excel existente si ya
    existe, o crea un archivo dedicado goodvibes_appendix.xlsx.
    """
    if df_gv.empty:
        return

    mol_dir   = RC_WORK_ROOT / mol_name
    base_name = f"{mol_name}_goodvibes_appendix"

    # Columnas a exportar (en orden lógico)
    export_cols = [
        "conformer", "E_SP_hartree",
        "G_corr_ORCA_Eh", "G_corr_GV_Eh", "dG_corr_Eh", "dG_corr_kcal",
        "S_ORCA_EhK", "S_GV_EhK", "dS_EhK", "ZPE_GV_Eh",
        "G_total_ORCA_Eh", "G_total_GV_Eh",
        "dG_ORCA_kcal", "dG_GV_kcal", "ddG_rel_kcal",
        "boltz_ORCA_298K", "boltz_GV_298K", "d_boltz_298K",
        "rank_ORCA", "rank_GV", "rank_change_GV",
        "gv_ok", "gv_warning",
    ]
    export_cols = [c for c in export_cols if c in df_gv.columns]
    df_out = df_gv[export_cols].copy()

    if fmt in ("csv", "both"):
        csv_path = mol_dir / f"{base_name}.csv"
        df_out.to_csv(csv_path, index=False)
        print(f"  GoodVibes appendix CSV → {csv_path}")

    if fmt in ("xlsx", "both"):
        # Intentar agregar hoja al Excel existente
        xlsx_main = mol_dir / f"{mol_name}_results.xlsx"
        xlsx_gv   = mol_dir / f"{base_name}.xlsx"

        target = xlsx_main if xlsx_main.exists() else xlsx_gv
        mode   = "a" if xlsx_main.exists() else "w"
        kw     = {"if_sheet_exists": "replace"} if mode == "a" else {}

        try:
            with pd.ExcelWriter(target, engine="openpyxl", mode=mode, **kw) as writer:
                df_out.to_excel(writer, sheet_name="GoodVibes_mRRHO", index=False)
            print(f"  GoodVibes appendix XLSX → {target} (hoja: GoodVibes_mRRHO)")
        except Exception as exc:
            # Si falla al agregar al Excel existente, escribir archivo separado
            warnings.warn(f"[GoodVibes export] No se pudo agregar a {target}: {exc}. "
                          f"Escribiendo {xlsx_gv}.")
            with pd.ExcelWriter(xlsx_gv, engine="openpyxl") as writer:
                df_out.to_excel(writer, sheet_name="GoodVibes_mRRHO", index=False)
            print(f"  GoodVibes appendix XLSX → {xlsx_gv}")


# ═══════════════════════════════════════════════════════════════════════════════
#  MÓDULO MULTIWFN/ESP — potencial electrostático de superficie molecular
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Fundamento:
#  El ESP mapeado sobre la superficie de isodensidad electrónica (0.001 u.a.,
#  criterio Bader/Politzer) identifica sitios susceptibles a ataque
#  electrofílico/nucleofílico -- relevante para racionalizar reactividad y
#  patrones de solvatación de las MAAs. Multiwfn calcula extremos (mínimos y
#  máximos locales), área positiva/negativa, promedio, varianza e índices de
#  polaridad molecular (MPI, Π, ν) a partir de la función de onda (.molden).
#
#  Problema resuelto por este módulo:
#  Multiwfn es 100% interactivo por menú numérico via stdin. Dentro del
#  submenú de análisis cuantitativo de superficie existe un comando de
#  VISUALIZACIÓN en consola ("ver mínimos y máximos") que activa un modo de
#  paginación y espera un Enter adicional antes de devolver el control al
#  menú. Si el resto de comandos de post-procesado se manda por stdin sin
#  contar exactamente ese Enter extra, la secuencia se desincroniza y
#  Multiwfn se queda "colgado" (o traga comandos posteriores como si fueran
#  la respuesta a la paginación).
#
#  Estrategia adoptada (sin depender de contar Enters a ciegas):
#    1. Nunca se invoca el comando interactivo de "ver" extremos en consola.
#       En su lugar se usa el comando de EXPORTAR esos mismos datos a un
#       archivo -- es una operación de escritura a disco, siempre retorna
#       limpiamente al menú sin paginación.
#    2. El batch completo de números de menú se arma como un solo string
#       (_build_multiwfn_esp_batch) y se manda de una sola vez via
#       subprocess.run(input=...), nunca interactuando línea a línea con el
#       proceso vivo: el pipe del SO buferiza el stdin completo sin importar
#       el timing exacto en que Multiwfn lo vaya leyendo.
#    3. subprocess.run(..., timeout=MULTIWFN_TIMEOUT) acota el peor caso: si
#       de todas formas algo queda esperando input, el proceso se mata solo
#       en vez de trabar el pipeline completo.
#    4. La verificación de éxito es POR CONTENIDO del log (se buscan las
#       cadenas que Multiwfn imprime al completar el análisis cuantitativo:
#       "Minimal value" / "Maximal value" / "Global surface..."), nunca por
#       conteo de líneas o de prompts respondidos.
#
#  Flujo por confórmero (runs/<mol>/05_tddft/conf_XXX/):
#    _find_esp_wavefunction()    → localiza el .gbw de esa etapa
#    _generate_molden_input()    → orca_2mkl <base> -molden
#    _build_multiwfn_esp_batch() → arma el batch de comandos de menú
#    run_multiwfn_esp_on_conf()  → corre Multiwfn con ese batch (subprocess.run
#                                   con input=..., timeout=..., sin interacción viva)
#                                   y parsea el log + localiza el archivo exportado
#    run_multiwfn_esp_on_mol()   → itera todos los conf_* de 05_tddft/
#    compute_esp_summary_df()    → arma el DataFrame resumen (+ Boltzmann si aplica)
#    print_esp_summary()         → tabla formateada al estilo del script
#    export_esp_results()        → CSV/XLSX (hoja "ESP_Multiwfn")
#
# ═══════════════════════════════════════════════════════════════════════════════

# ── Secuencia de menú Multiwfn para el análisis ESP ─────────────────────────
# Verificar contra la versión instalada (2026.6.2 en este equipo, ver skill
# compquim-kmilozeta) corriendo Multiwfn interactivamente una vez sobre un
# .molden.input de prueba. Si los números difieren, solo hay que tocar estas
# constantes -- el resto del harness (batch único + timeout + verificación
# por contenido) no depende de los valores exactos.
MULTIWFN_ESP_MENU_MAIN           = "12"  # Quantitative analysis of molecular surface
MULTIWFN_ESP_MENU_START          = "0"   # Start analysis now! (ESP es la función mapeada default)
MULTIWFN_ESP_MENU_EXPORT_EXTREMA = "2"   # Exportar extremos de superficie a archivo
                                          # (deliberadamente NO la opción de "ver"/mostrar
                                          #  en consola -- esa es la que se queda pegada)
MULTIWFN_ESP_MENU_BACK            = "q"   # volver un nivel de menú
MULTIWFN_ESP_MENU_QUIT           = "q"   # salir de Multiwfn limpiamente

# NOTA IMPORTANTE (bug encontrado 2026-07-16): en el menú de post-análisis de
# superficie, la opción numerada "0" NO es un "volver" genérico -- antes de
# correr el análisis "0" = "Start analysis now!", pero UNA VEZ que el análisis
# ya corrió, ese mismo "0" cambia de significado a "Enter the graphical
# interface" (abre una ventana GUI de visualización que solo se cierra con
# clic manual en "Return" dentro de la ventana; ver sobereva.com/196, sección
# 2: "此时会看到后处理菜单。选选项0进入图形界面"). Bajo subprocess.run(input=...)
# sin display interactivo, Multiwfn se queda esperando esa interacción para
# siempre -- este era el punto exacto donde el pipeline se colgaba (justo
# después de exportar el .pdb con la opción "2", que sí es no-interactiva y
# termina normalmente). La convención correcta de Multiwfn para retroceder de
# nivel de menú es "q", nunca un número -- los números son siempre opciones
# de menú específicas del contexto.


def _build_multiwfn_esp_batch() -> str:
    """
    Arma el batch completo de comandos de menú de Multiwfn para el análisis
    cuantitativo de superficie ESP como UN SOLO string separado por saltos
    de línea, listo para pasarse entero a subprocess.run(input=...).

    Al mandarse de una sola vez (en vez de interactuar con el proceso vivo),
    el timing exacto de cada prompt deja de importar; lo único que importa es
    que la SECUENCIA de respuestas sea correcta. Se evita explícitamente
    cualquier comando de "ver"/"mostrar en consola" (ver docstring del
    módulo) usando en su lugar el comando de EXPORTAR a archivo, y se evita
    reusar "0" para "volver" (ver nota arriba) usando "q" en su lugar.
    """
    commands = [
        MULTIWFN_ESP_MENU_MAIN,              # 1. entrar a análisis cuantitativo de superficie
        MULTIWFN_ESP_MENU_START,              # 2. iniciar análisis (imprime resumen: min/max/área/...)
        MULTIWFN_ESP_MENU_EXPORT_EXTREMA,     # 3. exportar extremos a archivo (NO "ver" en consola)
        MULTIWFN_ESP_MENU_BACK,               # 4. volver al menú de función 12 (NO "0" -- eso abre la GUI)
        MULTIWFN_ESP_MENU_BACK,               # 5. volver al menú principal de Multiwfn
        MULTIWFN_ESP_MENU_QUIT,               # 6. salir de Multiwfn
    ]
    return "\n".join(commands) + "\n"


def _find_esp_wavefunction(conf_dir: Path) -> Optional[Path]:
    """
    Localiza el .gbw a convertir para el análisis ESP dentro de conf_dir
    (runs/<mol>/05_tddft/conf_XXX/). Se prioriza "tddft.gbw" (nombre estándar
    del job de esa etapa); si no existe, se cae al .gbw más reciente por mtime
    presente en el directorio.
    """
    preferred = conf_dir / "tddft.gbw"
    if preferred.exists():
        return preferred
    candidates = sorted(conf_dir.glob("*.gbw"), key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0] if candidates else None


def _generate_molden_input(gbw_path: Path, esp_dir: Path) -> tuple[Optional[Path], str]:
    """
    Corre `orca_2mkl <base> -molden` (cwd = directorio del .gbw) para generar
    <base>.molden.input, y lo mueve a esp_dir para no ensuciar 05_tddft/.

    Retorna (path_al_molden_input | None, mensaje_de_advertencia).
    """
    base = gbw_path.stem  # p.ej. "tddft" desde "tddft.gbw"
    try:
        proc = subprocess.run(
            [ORCA_2MKL_BIN, base, "-molden"],
            capture_output=True, text=True, timeout=ORCA_2MKL_TIMEOUT,
            cwd=str(gbw_path.parent),
        )
    except subprocess.TimeoutExpired:
        return None, f"orca_2mkl timeout (>{ORCA_2MKL_TIMEOUT}s) sobre {gbw_path.name}"
    except FileNotFoundError:
        return None, f"orca_2mkl no encontrado en {ORCA_2MKL_BIN}"

    molden_src = gbw_path.parent / f"{base}.molden.input"
    if proc.returncode != 0 or not molden_src.exists():
        return None, (f"orca_2mkl retornó código {proc.returncode} o no generó "
                       f"{molden_src.name}. stderr: {proc.stderr[:300]}")

    esp_dir.mkdir(exist_ok=True)
    molden_dst = esp_dir / molden_src.name
    try:
        molden_src.replace(molden_dst)
    except OSError as exc:
        return None, f"No se pudo mover {molden_src} -> {molden_dst}: {exc}"

    return molden_dst, ""


def _find_extrema_export(esp_dir: Path) -> Optional[Path]:
    """Ubica el archivo de extremos de superficie exportado por Multiwfn (si existe)."""
    for pattern in ("*extrem*", "*CPs*", "surfanalysis*.txt", "surfanalysis*.pdb"):
        hits = sorted(esp_dir.glob(pattern))
        if hits:
            return hits[0]
    return None


def _parse_multiwfn_esp_log(log_text: str) -> dict:
    """
    Extrae del log de Multiwfn el bloque de resumen cuantitativo de
    superficie ESP que se imprime automáticamente tras "Start analysis now!".

    Tolerante a variaciones menores de formato entre versiones (regex laxas);
    cualquier campo no encontrado queda en None sin abortar el parseo de los
    demás -- mismo criterio defensivo que el resto de parsers del script.
    """
    result = {
        "ESP_min_kcal": None, "ESP_max_kcal": None,
        "area_total_ang2": None, "area_pos_ang2": None, "area_neg_ang2": None,
        "avg_kcal": None, "variance_kcal2": None,
        "MPI_kcal": None, "nu": None, "pi_kcal2": None,
        "ok": False, "warning": "",
    }
    if not log_text or ("Minimal value" not in log_text and "Global surface" not in log_text):
        result["warning"] = ("Bloque de resumen ESP no detectado en el log de Multiwfn "
                              "(revisar log completo y ajustar menú si la versión difiere)")
        return result

    def _grab(pattern: str) -> Optional[float]:
        m = re.search(pattern, log_text)
        if not m:
            return None
        try:
            return float(m.group(1))
        except (ValueError, IndexError):
            return None

    # Multiwfn suele reportar en u.a. y, entre paréntesis, en kcal/mol —
    # priorizamos kcal/mol por consistencia con el resto del script.
    result["ESP_min_kcal"] = (
        _grab(r"Minimal value[^\n]*?\(\s*([-\d.]+)\s*kcal/mol\s*\)")
        or _grab(r"Global surface minimum[^\n]*?\(\s*([-\d.]+)\s*kcal/mol\s*\)")
    )
    result["ESP_max_kcal"] = (
        _grab(r"Maximal value[^\n]*?\(\s*([-\d.]+)\s*kcal/mol\s*\)")
        or _grab(r"Global surface maximum[^\n]*?\(\s*([-\d.]+)\s*kcal/mol\s*\)")
    )
    result["area_total_ang2"] = _grab(r"Overall surface area[^\n]*?\(\s*([\d.]+)\s*Angstrom")
    result["area_pos_ang2"]   = _grab(r"Positive surface area[^\n]*?\(\s*([\d.]+)\s*Angstrom")
    result["area_neg_ang2"]   = _grab(r"Negative surface area[^\n]*?\(\s*([\d.]+)\s*Angstrom")
    result["avg_kcal"]        = _grab(r"Overall average value[^\n]*?\(\s*([-\d.]+)\s*kcal/mol\s*\)")
    result["variance_kcal2"]  = _grab(r"Overall variance[^\n]*?\(\s*([-\d.]+)\s*kcal")
    result["nu"]              = _grab(r"[Nn]u[^\n]*?[:=]\s*([-\d.]+)")
    result["MPI_kcal"]        = _grab(r"MPI[^\n]*?[:=]\s*([-\d.]+)")
    result["pi_kcal2"]        = _grab(r"[Ii]nternal charge separation[^\n]*?[:=]\s*([-\d.]+)")

    # Éxito mínimo: al menos el par min/max de superficie (el dato original
    # que se quería obtener sin que el proceso se quedara pegado).
    result["ok"] = result["ESP_min_kcal"] is not None or result["ESP_max_kcal"] is not None
    if not result["ok"]:
        result["warning"] = ("No se pudo extraer ESP_min/ESP_max del log -- revisar formato "
                              "de salida de esta versión de Multiwfn y ajustar los patrones "
                              "regex en _parse_multiwfn_esp_log().")
    return result


def run_multiwfn_esp_on_conf(conf_dir: Path) -> dict:
    """
    Ejecuta el flujo ESP/Multiwfn completo para un único confórmero:
      1. Localiza el .gbw de la etapa TDDFT.
      2. orca_2mkl -> .molden.input (en conf_dir/esp_multiwfn/)
      3. Multiwfn con el batch pre-armado, vía subprocess.run(input=...)
      4. Parsea el log buscando el bloque de resumen cuantitativo.

    No interrumpe el pipeline si algo falla: retorna dict con "ok": False y
    "warning" describiendo el problema, igual que run_goodvibes_on_mol().
    """
    null_entry = {
        "ESP_min_kcal": None, "ESP_max_kcal": None,
        "area_total_ang2": None, "area_pos_ang2": None, "area_neg_ang2": None,
        "avg_kcal": None, "variance_kcal2": None,
        "MPI_kcal": None, "nu": None, "pi_kcal2": None,
        "ok": False, "log_path": None, "extrema_path": None, "warning": "",
    }

    gbw_path = _find_esp_wavefunction(conf_dir)
    if gbw_path is None:
        null_entry["warning"] = f"No se encontró .gbw en {conf_dir}"
        return null_entry

    esp_dir = conf_dir / ESP_SUBDIR
    molden_path, warn = _generate_molden_input(gbw_path, esp_dir)
    if molden_path is None:
        null_entry["warning"] = warn
        return null_entry

    batch = _build_multiwfn_esp_batch()

    # OMP_NUM_THREADS controla la paralelización interna de Multiwfn para el
    # análisis de superficie ESP. MKL/OpenBLAS se fijan a 1 para evitar
    # sobre-suscripción de núcleos (mismo criterio que el resto del pipeline
    # de compquim-kmilozeta).
    multiwfn_env = os.environ.copy()
    multiwfn_env["OMP_NUM_THREADS"] = str(MULTIWFN_NTHREADS)
    multiwfn_env["MKL_NUM_THREADS"] = "1"
    multiwfn_env["OPENBLAS_NUM_THREADS"] = "1"

    try:
        proc = subprocess.run(
            [MULTIWFN_BIN, str(molden_path)],
            input=batch,
            capture_output=True, text=True,
            timeout=MULTIWFN_TIMEOUT,
            cwd=str(esp_dir),
            env=multiwfn_env,
        )
    except subprocess.TimeoutExpired as exc:
        # Aunque haya vencido el timeout, guardamos lo que alcanzó a
        # producirse en stdout para poder diagnosticar el punto de bloqueo.
        log_path = esp_dir / "multiwfn_esp.log"
        partial = exc.stdout if isinstance(exc.stdout, str) else ""
        if partial:
            log_path.write_text(partial, errors="replace")
        null_entry["warning"] = (
            f"Multiwfn timeout (>{MULTIWFN_TIMEOUT}s) sobre {molden_path.name}. "
            f"{'Log parcial en ' + str(log_path) if partial else 'Sin salida capturada.'}"
        )
        return null_entry
    except FileNotFoundError:
        null_entry["warning"] = f"Multiwfn no encontrado en {MULTIWFN_BIN}"
        return null_entry

    log_path = esp_dir / "multiwfn_esp.log"
    log_path.write_text(proc.stdout + "\n\n[stderr]\n" + proc.stderr, errors="replace")

    parsed = _parse_multiwfn_esp_log(proc.stdout)
    parsed["log_path"] = log_path
    parsed["extrema_path"] = _find_extrema_export(esp_dir) if parsed["ok"] else None
    if not parsed["ok"] and not parsed["warning"]:
        parsed["warning"] = f"Ver log completo en {log_path}"
    return parsed


def run_multiwfn_esp_on_mol(mol_name: str) -> dict[str, dict]:
    """
    Corre el análisis ESP/Multiwfn sobre todos los confórmeros de
    runs/<mol>/05_tddft/conf_*/.

    Retorna dict conf_id -> resultado (ver run_multiwfn_esp_on_conf). Un
    confórmero que falla no interrumpe el resto, igual que en GoodVibes.
    """
    mol_dir   = RC_WORK_ROOT / mol_name
    tddft_dir = mol_dir / ESP_SOURCE_SUBDIR

    if not tddft_dir.exists():
        warnings.warn(f"[Multiwfn/ESP] {tddft_dir} no existe — omitiendo {mol_name}.")
        return {}

    conf_dirs = sorted(tddft_dir.glob("conf_*"))
    if not conf_dirs:
        warnings.warn(f"[Multiwfn/ESP] Sin subdirectorios conf_* en {tddft_dir}.")
        return {}

    if not Path(MULTIWFN_BIN).exists():
        warnings.warn(f"[Multiwfn/ESP] No se encontró el binario en {MULTIWFN_BIN}.")
        return {c.name: {"ok": False, "warning": f"Multiwfn no disponible en {MULTIWFN_BIN}"}
                for c in conf_dirs}
    if not Path(ORCA_2MKL_BIN).exists():
        warnings.warn(f"[Multiwfn/ESP] No se encontró orca_2mkl en {ORCA_2MKL_BIN}.")
        return {c.name: {"ok": False, "warning": f"orca_2mkl no disponible en {ORCA_2MKL_BIN}"}
                for c in conf_dirs}

    print(f"  [Multiwfn/ESP] {len(conf_dirs)} confórmero(s) en {tddft_dir}")

    results: dict[str, dict] = {}
    for conf_dir in conf_dirs:
        conf_id = conf_dir.name
        print(f"    · {conf_id} ... ", end="", flush=True)
        res = run_multiwfn_esp_on_conf(conf_dir)
        results[conf_id] = res
        print("OK" if res["ok"] else f"OMITIDO ({res['warning'][:80]})")

    return results


def compute_esp_summary_df(df_conf: "pd.DataFrame",
                            esp_results: dict[str, dict]) -> "pd.DataFrame":
    """
    Arma el DataFrame resumen ESP/Multiwfn. Si df_conf (de collect_results)
    trae columnas boltz_*, se hace merge por conformer para poder cruzar
    extremos ESP con población Boltzmann de cada confórmero.
    """
    rows = []
    for conf_id, res in esp_results.items():
        row = {"conformer": conf_id}
        row.update({k: v for k, v in res.items()
                     if k not in ("log_path", "extrema_path")})
        row["log_path"]     = str(res.get("log_path") or "")
        row["extrema_path"] = str(res.get("extrema_path") or "")
        rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        return df

    if df_conf is not None and not df_conf.empty and "conformer" in df_conf.columns:
        boltz_cols = [c for c in df_conf.columns if c.startswith("boltz_")]
        if boltz_cols:
            df = df.merge(df_conf[["conformer"] + boltz_cols], on="conformer", how="left")

    return df


def print_esp_summary(df_esp: "pd.DataFrame", mol_name: str):
    print(f"\n{chr(9472)*74}")
    print(f"  ESP DE SUPERFICIE (Multiwfn) -- {mol_name}")
    print(f"{chr(9472)*74}")

    if df_esp.empty:
        print("  Sin resultados ESP disponibles (ver advertencias arriba).")
        return

    n_ok = int(df_esp["ok"].sum()) if "ok" in df_esp.columns else 0
    print(f"  Confórmeros con análisis exitoso: {n_ok}/{len(df_esp)}")

    cols = [c for c in ("conformer", "ESP_min_kcal", "ESP_max_kcal",
                         "area_total_ang2", "avg_kcal", "MPI_kcal", "ok")
            if c in df_esp.columns]
    print(df_esp[cols].to_string(index=False))


def export_esp_results(df_esp: "pd.DataFrame", mol_name: str, fmt: str):
    if df_esp.empty:
        return
    mol_dir = RC_WORK_ROOT / mol_name
    if fmt in ("csv", "both"):
        out = mol_dir / "esp_multiwfn_summary.csv"
        df_esp.to_csv(out, index=False)
        print(f"  ESP CSV -> {out}")
    if fmt in ("xlsx", "both"):
        xlsx_path = mol_dir / f"{mol_name}_results.xlsx"
        mode = "a" if xlsx_path.exists() else "w"
        kw   = {"if_sheet_exists": "replace"} if mode == "a" else {}
        try:
            with pd.ExcelWriter(xlsx_path, engine="openpyxl", mode=mode, **kw) as w:
                df_esp.to_excel(w, sheet_name="ESP_Multiwfn", index=False)
            print(f"  ESP -> {xlsx_path} (hoja: ESP_Multiwfn)")
        except Exception as exc:
            warnings.warn(f"[ESP export] {exc}")


# ═══════════════════════════════════════════════════════════════════════════════
#  MÓDULO ADC(2) — extracción de λ_max y comparación con TD-DFT
# ═══════════════════════════════════════════════════════════════════════════════
#
#  Fundamento:
#  ADC(2) (Algebraic Diagrammatic Construction, segundo orden) es un método
#  ab initio para estados excitados más preciso que TD-DFT pero también más
#  costoso. La comparación sistemática entre λ_max(TD-DFT) y λ_max(ADC2)
#  permite evaluar si el error del funcional es sistemático (desplazamiento
#  constante) y si puede corregirse empíricamente aplicando un Δλ = λ_TDDFT - λ_ADC2
#  uniforme al espectro UV-Vis calculado con TD-DFT.
#
#  Estructura de archivos esperada:
#    runs/<mol>/07_adc2/conf_<id>/adc2.out
#
#  Por cada conformero se parsea UN root (nroots 1), extrayendo:
#    - E_eV     : energía de excitación vertical (eV)
#    - wl_nm    : longitud de onda (nm)
#    - fosc     : fuerza de oscilador (del bloque LEFT-RIGHT, el más equilibrado)
#
#  Luego se calcula λ_max ADC2 Boltzmann-promediada al igual que en TDDFT y
#  se compara estado a estado con λ_max TDDFT.
#
# ═══════════════════════════════════════════════════════════════════════════════

ADC2_SUBDIR = "07_adc2"   # subcarpeta dentro de runs/<mol>/


def parse_adc2_excitations(adc2_out: Path) -> list[dict]:
    """
    Extrae la(s) excitación(es) del bloque ADC(2) de un output ORCA 6.

    Busca el bloque "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
    del espectro LEFT-RIGHT (el más equilibrado entre los tres que imprime ORCA
    para ADC2: RIGHT, LEFT y LEFT-RIGHT). Si no encuentra el bloque LEFT-RIGHT,
    usa el bloque RIGHT como fallback.

    Formato de la línea de datos (ORCA 6):
        0-1A  ->  1-1A    3.464792   27945.4   357.8   0.416153266   ...

    Retorna lista de dicts con:
        state      : int  — número de root (1, 2, …)
        E_eV       : float — energía de excitación (eV)
        E_cm1      : float — energía en cm⁻¹
        wl_nm      : float — longitud de onda (nm)
        fosc       : float — fuerza de oscilador (LEFT-RIGHT si disponible)
        fosc_right : float — fuerza de oscilador RIGHT (siempre disponible)
        fosc_left  : float — fuerza de oscilador LEFT (si disponible)
        source     : str  — qué bloque se usó como fosc principal
    """
    if not adc2_out.exists():
        return []

    try:
        text = adc2_out.read_text(errors="replace")
    except OSError:
        return []

    lines = text.splitlines()

    # Patrón para líneas de excitación del bloque de absorción:
    # "  0-1A  ->  1-1A    3.464792   27945.4   357.8   0.416153266  ..."
    pat_line = re.compile(
        r"^\s*(\d+)-\S+\s*->\s*(\d+)-\S+\s+"   # initial_root -> final_root
        r"([\d.]+)\s+"                           # E_eV
        r"([\d.]+)\s+"                           # E_cm-1
        r"([\d.]+)\s+"                           # wl_nm
        r"([\d.]+)"                              # fosc
    )

    # Buscar los tres bloques de absorción que imprime ADC2.
    # ORCA 6 los imprime en este orden: RIGHT, LEFT, LEFT-RIGHT.
    # Cada bloque está precedido por su título y la cabecera de columnas.
    block_tags = {
        "right":      "SPECTRUM FOR RIGHT TRANSITION MOMENTS",
        "left":       "SPECTRUM FOR LEFT TRANSITION MOMENTS",
        "left_right": "SPECTRUM FOR LEFT-RIGHT TRANSITION MOMENTS",
    }
    # La línea de datos de absorción viene tras el separador "---" y la cabecera;
    # marcamos el inicio del bloque con "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC"
    abs_header = "ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"

    # Recolectar excitaciones por tipo de bloque
    collected: dict[str, list[dict]] = {"right": [], "left": [], "left_right": []}

    current_block: Optional[str] = None
    in_abs_block   = False

    for line in lines:
        # Detectar a qué bloque de momentos de transición entramos
        for key, tag in block_tags.items():
            if tag in line:
                current_block = key
                in_abs_block  = False
                break

        # Dentro del bloque correcto, esperar el encabezado de absorción
        if current_block and not in_abs_block and abs_header in line:
            in_abs_block = True
            continue

        if in_abs_block and current_block:
            m = pat_line.match(line)
            if m:
                root_initial = int(m.group(1))
                root_final   = int(m.group(2))
                # El estado excitado es root_final (convenio ORCA: "0->1" = S1)
                state_n = root_final if root_final > 0 else 1
                collected[current_block].append({
                    "state":  state_n,
                    "E_eV":   float(m.group(3)),
                    "E_cm1":  float(m.group(4)),
                    "wl_nm":  float(m.group(5)),
                    "fosc":   float(m.group(6)),
                })
            # Salir del bloque si encontramos una línea vacía seguida de "---"
            # o el comienzo de otro bloque (se detectará en el bucle externo)

    # Construir lista final cruzando los tres bloques
    # Usamos LEFT-RIGHT como fosc principal; si no existe, RIGHT.
    primary   = collected["left_right"] or collected["right"]
    right_map = {e["state"]: e["fosc"] for e in collected["right"]}
    left_map  = {e["state"]: e["fosc"] for e in collected["left"]}

    result = []
    for exc in primary:
        s = exc["state"]
        result.append({
            "state":      s,
            "E_eV":       exc["E_eV"],
            "E_cm1":      exc["E_cm1"],
            "wl_nm":      exc["wl_nm"],
            "fosc":       exc["fosc"],
            "fosc_right": right_map.get(s, float("nan")),
            "fosc_left":  left_map.get(s, float("nan")),
            "source":     "left_right" if collected["left_right"] else "right",
        })
    return result


def collect_adc2_results(mol_name: str) -> pd.DataFrame:
    """
    Recopila todas las excitaciones ADC(2) de los confórmeros de una molécula.

    Estructura esperada:
        runs/<mol>/07_adc2/conf_<id>/adc2.out

    Retorna DataFrame con columnas:
        conformer, state, E_eV, E_cm1, wl_nm, fosc, fosc_right, fosc_left, source
    O DataFrame vacío si no hay archivos.
    """
    mol_dir   = RC_WORK_ROOT / mol_name
    adc2_dir  = mol_dir / ADC2_SUBDIR
    if not adc2_dir.exists():
        return pd.DataFrame()

    rows = []
    for conf_dir in sorted(adc2_dir.glob("conf_*")):
        if not conf_dir.is_dir():
            continue
        conf_id  = conf_dir.name
        adc2_out = conf_dir / "adc2.out"
        for exc in parse_adc2_excitations(adc2_out):
            rows.append({"conformer": conf_id, **exc})

    return pd.DataFrame(rows) if rows else pd.DataFrame()


def compute_adc2_lambda_max(
    df_adc2: pd.DataFrame,
    df_conf: pd.DataFrame,
    temp_label: str = "298K",
) -> Optional[float]:
    """
    Calcula λ_max ADC(2) Boltzmann-promediado (igual que en TDDFT).

    Usa la columna boltz_{temp_label} de df_conf.
    Sólo considera el estado S1 (state == 1).

    Retorna λ_max en nm o None si no hay datos suficientes.
    """
    if df_adc2.empty or df_conf.empty:
        return None

    col_b = f"boltz_{temp_label}"
    if col_b not in df_conf.columns:
        return None

    s1 = df_adc2[df_adc2["state"] == 1].copy()
    if s1.empty:
        return None

    s1 = s1.merge(df_conf[["conformer", col_b]], on="conformer", how="left")
    w  = s1[col_b].fillna(0)
    if w.sum() <= 0:
        return None

    return float((s1["wl_nm"] * w).sum() / w.sum())


# ─────────────────────────────────────────────
#  COMPARACIÓN ADC2 vs TDDFT (ΔTDDFT-ADC2)
# ─────────────────────────────────────────────

def compare_adc2_tddft(
    mols: Optional[list[str]] = None,
    temp_label: str = "298K",
) -> pd.DataFrame:
    """
    Compara λ_max(TD-DFT) vs λ_max(ADC2) entre estados de protonación.

    Para cada molécula disponible en runs/ que tenga datos TD-DFT Y ADC2,
    calcula:
        lambda_tddft_nm   : λ_max Boltzmann-promediado TD-DFT (nm)
        lambda_adc2_nm    : λ_max Boltzmann-promediado ADC2 (nm)
        delta_nm          : λ_TDDFT - λ_ADC2 (nm)
                            > 0 → TDDFT sobreestima λ (redshift)
                            < 0 → TDDFT subestima λ (blueshift)
        delta_eV          : conversión a eV  (Δ positivo = TDDFT demasiado bajo en energía)
        n_confs_adc2      : número de confórmeros con ADC2 completado
        n_confs_tddft     : número de confórmeros con TDDFT disponible

    También imprime un análisis de sistematicidad:
        - ¿Es Δ constante entre estados? (std_delta → corrección sistemática viable)
        - Corrección empírica recomendada: shift = mean(Δ) ± std(Δ)

    Parámetros
    ----------
    mols        : lista de nombres de mol (directorios en runs/). None → todos.
    temp_label  : temperatura para el Boltzmann (default "298K").

    Retorna DataFrame con una fila por molécula.
    """
    if mols is None:
        mols = sorted(p.name for p in RC_WORK_ROOT.iterdir()
                      if p.is_dir() and (p / "04_sp").exists())

    rows = []
    for mol in mols:
        # ── Cargar TDDFT ───────────────────────────────────────────────────
        try:
            df_conf, df_tddft = collect_results(mol)
        except Exception as e:
            warnings.warn(f"[ADC2 vs TDDFT] {mol}: error al cargar TDDFT — {e}")
            df_conf  = pd.DataFrame()
            df_tddft = pd.DataFrame()

        lmax_tddft: Optional[float] = None
        n_confs_tddft = 0
        if not df_tddft.empty and not df_conf.empty:
            lmax_tddft    = df_conf.attrs.get("lambda_max_boltz_298K") \
                            if temp_label == "298K" else \
                            _boltz_weighted_lmax(df_tddft, df_conf, temp_label)
            n_confs_tddft = int(df_tddft["conformer"].nunique())

        # ── Cargar ADC2 ────────────────────────────────────────────────────
        df_adc2      = collect_adc2_results(mol)
        lmax_adc2: Optional[float] = None
        n_confs_adc2 = 0
        if not df_adc2.empty and not df_conf.empty:
            lmax_adc2    = compute_adc2_lambda_max(df_adc2, df_conf, temp_label)
            n_confs_adc2 = int(df_adc2["conformer"].nunique())

        # ── Δ = TDDFT − ADC2 ──────────────────────────────────────────────
        delta_nm: Optional[float] = None
        delta_eV: Optional[float] = None
        if lmax_tddft is not None and lmax_adc2 is not None and lmax_adc2 > 0:
            delta_nm = lmax_tddft - lmax_adc2
            # Conversión eV: E = hc/λ → ΔE ≈ -hc·Δλ/λ² (diferencia de energías)
            # Calculamos directamente como diferencia de energías en eV:
            ev_tddft = 1239.84193 / lmax_tddft   # eV
            ev_adc2  = 1239.84193 / lmax_adc2    # eV
            delta_eV = ev_tddft - ev_adc2         # >0 → TDDFT subestima E (redshift)

        rows.append({
            "molecule":        mol,
            "lambda_tddft_nm": lmax_tddft,
            "lambda_adc2_nm":  lmax_adc2,
            "delta_nm":        delta_nm,        # TDDFT − ADC2 en nm
            "delta_eV":        delta_eV,        # TDDFT − ADC2 en eV
            "n_confs_tddft":   n_confs_tddft,
            "n_confs_adc2":    n_confs_adc2,
        })

    return pd.DataFrame(rows)


def _boltz_weighted_lmax(
    df_tddft: pd.DataFrame,
    df_conf:  pd.DataFrame,
    temp_label: str,
) -> Optional[float]:
    """Helper: λ_max TDDFT Boltzmann-promediado a temperatura arbitraria."""
    col_b = f"boltz_{temp_label}"
    if col_b not in df_conf.columns:
        return None
    s1 = df_tddft[df_tddft["state"] == 1].copy()
    if s1.empty:
        return None
    s1 = s1.merge(df_conf[["conformer", col_b]], on="conformer", how="left")
    w  = s1[col_b].fillna(0)
    return float((s1["wavelength_nm"] * w).sum() / w.sum()) if w.sum() > 0 else None


def print_adc2_summary(df_adc2: pd.DataFrame, mol_name: str,
                       df_conf: pd.DataFrame):
    """Imprime tabla de excitaciones ADC2 de una molécula."""
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  ADC(2) EXCITACIONES — {mol_name}")
    print(sep)

    if df_adc2.empty:
        print(f"  Sin datos ADC2 para '{mol_name}'.")
        print(f"  Ruta esperada: runs/{mol_name}/{ADC2_SUBDIR}/conf_*/adc2.out")
        print(f"{sep}\n")
        return

    # S1 por confórmero con pesos Boltzmann
    col_b = "boltz_298K"
    s1 = df_adc2[df_adc2["state"] == 1].copy()
    if not s1.empty and not df_conf.empty and col_b in df_conf.columns:
        s1 = s1.merge(df_conf[["conformer", col_b]], on="conformer", how="left")
        disp_cols = ["conformer", "wl_nm", "E_eV", "fosc", col_b]
        disp_cols = [c for c in disp_cols if c in s1.columns]
        print("\n  ── S1 por confórmero (ADC2, ponderado 298 K) ──")
        print(s1[disp_cols].to_string(index=False,
              float_format=lambda x: f"{x:.4f}"))

    lmax = compute_adc2_lambda_max(df_adc2, df_conf, "298K")
    if lmax:
        print(f"\n  λ_max ADC2 Boltzmann-promediado (298 K, S1): {lmax:.1f} nm")

    n = df_adc2["conformer"].nunique()
    print(f"  Confórmeros ADC2 completados: {n}")
    print(f"{sep}\n")


def print_adc2_tddft_comparison(df_cmp: pd.DataFrame):
    """
    Imprime la tabla comparativa ADC2 vs TDDFT con análisis de sistematicidad.
    """
    sep = "═" * 80
    print(f"\n{sep}")
    print(f"  COMPARACIÓN ADC(2) vs TD-DFT — λ_max Boltzmann-promediado (298 K)")
    print(sep)

    if df_cmp.empty:
        print("  Sin datos para comparar.")
        return

    # ── Tabla por estado de protonación ─────────────────────────────────────
    print(f"\n  {'Molécula':<25} {'λ_TDDFT(nm)':>12} {'λ_ADC2(nm)':>12} "
          f"{'Δ(nm)':>10} {'Δ(eV)':>10} {'N_ADC2':>7}")
    print(f"  {'─'*25} {'─'*12} {'─'*12} {'─'*10} {'─'*10} {'─'*7}")

    valid_rows = []
    for _, row in df_cmp.iterrows():
        tddft = row.get("lambda_tddft_nm")
        adc2  = row.get("lambda_adc2_nm")
        delta = row.get("delta_nm")
        devv  = row.get("delta_eV")
        n     = row.get("n_confs_adc2", 0)

        flag = ""
        if tddft is None or adc2 is None:
            flag = "  ⚠ datos incompletos"

        print(f"  {row['molecule']:<25} "
              f"{'—':>12}" if tddft is None else
              f"  {row['molecule']:<25} {tddft:>12.1f}", end="")
        print(f" {adc2:>12.1f}" if adc2 is not None else f" {'—':>12}", end="")
        print(f" {delta:>+10.1f}" if delta is not None else f" {'—':>10}", end="")
        print(f" {devv:>+10.4f}" if devv is not None else f" {'—':>10}", end="")
        print(f" {int(n):>7}{flag}")

        if delta is not None:
            valid_rows.append(row)

    if not valid_rows:
        print(f"\n  Sin pares TDDFT+ADC2 completos para análisis de sistematicidad.")
        print(f"{sep}\n")
        return

    # ── Análisis de sistematicidad ───────────────────────────────────────────
    deltas_nm = [r["delta_nm"] for r in valid_rows]
    deltas_eV = [r["delta_eV"] for r in valid_rows]

    mean_nm = float(np.mean(deltas_nm))
    std_nm  = float(np.std(deltas_nm, ddof=1)) if len(deltas_nm) > 1 else 0.0
    mean_eV = float(np.mean(deltas_eV))
    std_eV  = float(np.std(deltas_eV, ddof=1)) if len(deltas_eV) > 1 else 0.0

    print(f"\n  ── Sistematicidad del desplazamiento TDDFT − ADC2 ──")
    print(f"  N pares válidos       : {len(deltas_nm)}")
    print(f"  Δ medio (nm)          : {mean_nm:+.2f} ± {std_nm:.2f}")
    print(f"  Δ medio (eV)          : {mean_eV:+.4f} ± {std_eV:.4f}")

    # Evaluación de sistematicidad: si std < |mean|/3 el desplazamiento es sistemático
    systematic = (std_nm < abs(mean_nm) / 3.0) if mean_nm != 0 else (std_nm < 2.0)

    if systematic and len(deltas_nm) > 1:
        direction = "redshift" if mean_nm > 0 else "blueshift"
        print(f"\n  ✓  Desplazamiento SISTEMÁTICO ({direction}):")
        print(f"     TD-DFT predice λ_max {abs(mean_nm):.1f} nm más {direction.replace('shift','ed')}"
              f" que ADC2 (CV = {std_nm/abs(mean_nm)*100:.1f}%)")
        print(f"\n  ► Corrección empírica recomendada:")
        print(f"     λ_corr(nm) = λ_TDDFT(nm) − ({mean_nm:+.1f} nm)")
        print(f"     E_corr(eV) = E_TDDFT(eV) − ({mean_eV:+.4f} eV)")
        print(f"     Incertidumbre: ±{std_nm:.1f} nm  (±{std_eV:.4f} eV)")
    elif len(deltas_nm) == 1:
        print(f"\n  ⚠  Solo 1 par: imposible evaluar sistematicidad.")
        print(f"     Desplazamiento observado: {mean_nm:+.1f} nm ({mean_eV:+.4f} eV)")
    else:
        print(f"\n  ⚠  Desplazamiento NO sistemático (std/|mean| = "
              f"{std_nm/abs(mean_nm)*100:.1f}% > 33%).")
        print(f"     La corrección empírica tiene alta incertidumbre (±{std_nm:.1f} nm).")
        print(f"     Considerar revisar geometrías o ampliar ensemble conformacional.")

    # ADC2 inter-estados: comparación relativa
    if len(valid_rows) > 1:
        lmax_adc2_vals = [r["lambda_adc2_nm"] for r in valid_rows]
        mols_valid     = [r["molecule"]        for r in valid_rows]
        print(f"\n  ── Comparación λ_max ADC2 entre estados de protonación ──")
        adc2_min  = min(lmax_adc2_vals)
        adc2_max  = max(lmax_adc2_vals)
        mol_min   = mols_valid[lmax_adc2_vals.index(adc2_min)]
        mol_max   = mols_valid[lmax_adc2_vals.index(adc2_max)]
        print(f"  Rango ADC2: {adc2_min:.1f} – {adc2_max:.1f} nm")
        print(f"  Estado con λ_max ADC2 menor: {mol_min} ({adc2_min:.1f} nm)")
        print(f"  Estado con λ_max ADC2 mayor: {mol_max} ({adc2_max:.1f} nm)")

    print(f"\n{sep}\n")


def export_adc2_comparison(df_cmp: pd.DataFrame, fmt: str):
    """Exporta la tabla comparativa ADC2 vs TDDFT a CSV y/o XLSX."""
    if df_cmp.empty:
        return
    out = RC_WORK_ROOT / "comparison_adc2_tddft"
    if fmt in ("csv", "both"):
        df_cmp.to_csv(out.with_suffix(".csv"), index=False)
        print(f"  CSV ADC2 vs TDDFT → {out.with_suffix('.csv')}")
    if fmt in ("xlsx", "both"):
        df_cmp.to_excel(out.with_suffix(".xlsx"), index=False)
        print(f"  XLSX ADC2 vs TDDFT → {out.with_suffix('.xlsx')}")


def plot_adc2_tddft_comparison(df_cmp: pd.DataFrame,
                                dfs_conf: dict, dfs_tddft: dict,
                                sigma: float = 20.0,
                                x_range: tuple = (200, 500)):
    """
    Grafica la comparación ADC2 vs TDDFT mostrando para cada estado:
      - Espectro TDDFT (gaussiano, normalizado)
      - Línea vertical λ_max TDDFT
      - Línea vertical λ_max ADC2
      - Flechas mostrando el desplazamiento Δ

    Panel inferior: barra de Δ(nm) por estado de protonación.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.patches as mpatches
    except ImportError:
        print("matplotlib no disponible — saltando plot ADC2")
        return

    if df_cmp.empty:
        return

    # Filtrar solo filas con datos completos
    valid = df_cmp.dropna(subset=["lambda_tddft_nm", "lambda_adc2_nm"])
    if valid.empty:
        print("  Sin pares TDDFT+ADC2 para graficar.")
        return

    n_mols = len(valid)
    palette = [
        "#4477AA", "#EE6677", "#228833", "#CCBB44",
        "#66CCEE", "#AA3377", "#BBBBBB"
    ]

    fig, axes = plt.subplots(
        2, 1, figsize=(11, 9),
        gridspec_kw={"height_ratios": [3, 1.2]}
    )
    ax_spec, ax_delta = axes

    x = np.linspace(x_range[0], x_range[1], 1500)

    for i, (_, row) in enumerate(valid.iterrows()):
        mol   = row["molecule"]
        color = palette[i % len(palette)]

        df_t = dfs_tddft.get(mol, pd.DataFrame())
        df_c = dfs_conf.get(mol, pd.DataFrame())

        # Espectro TDDFT
        if not df_t.empty and not df_c.empty:
            y = _build_spectrum(df_t, df_c, sigma, x, temp_label="298K")
            y_norm = y / y.max() if y.max() > 0 else y
            ax_spec.plot(x, y_norm, lw=1.8, alpha=0.7,
                         label=mol.replace("_", " "), color=color)

        # λ_max TDDFT
        lmax_td = row["lambda_tddft_nm"]
        ax_spec.axvline(lmax_td, color=color, lw=1.5, ls="--", alpha=0.85,
                        label=f"λ_TDDFT {mol.replace('_', ' ')} = {lmax_td:.0f} nm")

        # λ_max ADC2
        lmax_a2 = row["lambda_adc2_nm"]
        ax_spec.axvline(lmax_a2, color=color, lw=2.0, ls=":",
                        label=f"λ_ADC2  {mol.replace('_', ' ')} = {lmax_a2:.0f} nm")

        # Flecha TDDFT → ADC2
        y_arrow = 0.85 - i * 0.12
        if 0 < y_arrow < 1:
            ax_spec.annotate(
                "", xy=(lmax_a2, y_arrow), xytext=(lmax_td, y_arrow),
                arrowprops=dict(arrowstyle="<->", color=color, lw=1.4)
            )
            mid = (lmax_td + lmax_a2) / 2
            ax_spec.text(mid, y_arrow + 0.03,
                         f"Δ={row['delta_nm']:+.1f}nm",
                         ha="center", va="bottom", fontsize=8, color=color)

    ax_spec.set_xlabel("Longitud de onda (nm)", fontsize=12)
    ax_spec.set_ylabel("Intensidad TD-DFT (norm.)", fontsize=12)
    ax_spec.set_title(
        f"TD-DFT vs ADC(2) — λ_max por estado de protonación  "
        f"(σ = {sigma} nm, T = 298 K)", fontsize=12
    )
    # Leyenda compacta: solo las líneas de espectro
    handles = [
        plt.Line2D([0], [0], color=palette[i % len(palette)], lw=1.8, alpha=0.7,
                   label=row["molecule"].replace("_", " "))
        for i, (_, row) in enumerate(valid.iterrows())
    ]
    handles += [
        plt.Line2D([0], [0], color="gray", lw=1.5, ls="--", label="λ_max TD-DFT"),
        plt.Line2D([0], [0], color="gray", lw=2.0, ls=":", label="λ_max ADC(2)"),
    ]
    ax_spec.legend(handles=handles, fontsize=8, loc="upper right")

    # Panel inferior: Δ(nm) por estado
    mols_bar  = valid["molecule"].tolist()
    delta_bar = valid["delta_nm"].tolist()
    colors_bar = [palette[i % len(palette)] for i in range(len(mols_bar))]
    y_pos = range(len(mols_bar))

    bars = ax_delta.barh(list(y_pos), delta_bar, color=colors_bar, alpha=0.85)
    ax_delta.axvline(0, color="black", lw=0.8, ls="-")

    # Añadir valor numérico en cada barra
    for bar, val in zip(bars, delta_bar):
        xpos = val + 0.3 if val >= 0 else val - 0.3
        ha   = "left" if val >= 0 else "right"
        ax_delta.text(xpos, bar.get_y() + bar.get_height() / 2,
                      f"{val:+.1f} nm", va="center", ha=ha, fontsize=9)

    ax_delta.set_yticks(list(y_pos))
    ax_delta.set_yticklabels([m.replace("_", " ") for m in mols_bar], fontsize=9)
    ax_delta.set_xlabel("Δλ = λ_TDDFT − λ_ADC2  (nm)", fontsize=11)
    ax_delta.set_title("Desplazamiento TD-DFT vs ADC(2)  "
                       "[positivo = TDDFT redshifted]", fontsize=10)

    plt.tight_layout()
    out_png = RC_WORK_ROOT / "comparison_adc2_tddft.png"
    plt.savefig(out_png, dpi=150)
    print(f"  Plot ADC2 vs TDDFT guardado en: {out_png}")
    plt.show()


# ─────────────────────────────────────────────
#  EXPORTACIÓN
# ─────────────────────────────────────────────

def rc_export_results(df_conf: pd.DataFrame, df_tddft: pd.DataFrame,
                   mol_name: str, fmt: str):
    mol_dir = RC_WORK_ROOT / mol_name

    # Añadir columnas de Boltzmann para todas las temperaturas a df_tddft
    if not df_tddft.empty and not df_conf.empty:
        boltz_cols = [c for c in df_conf.columns if c.startswith("boltz_")]
        df_tddft = df_tddft.merge(
            df_conf[["conformer"] + boltz_cols + ["dE_kcal"]],
            on="conformer", how="left"
        )

    if fmt in ("csv", "both"):
        df_conf.to_csv(mol_dir / "conformers.csv", index=False)
        df_tddft.to_csv(mol_dir / "tddft_excitations.csv", index=False)
        print(f"  CSVs guardados en {mol_dir}/")

    if fmt in ("xlsx", "both"):
        xlsx_path = mol_dir / f"{mol_name}_results.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            df_conf.to_excel(writer, sheet_name="Conformers", index=False)
            df_tddft.to_excel(writer, sheet_name="TDDFT", index=False)

            # Hoja resumen CREST→DFT
            crest_cols = [c for c in df_conf.columns
                          if c in ("conformer", "rank_crest", "rank_dft",
                                   "rank_change", "dE_xtb_kcal",
                                   "dE_kcal", "dE_reranking_kcal")]
            if len(crest_cols) > 2:
                df_conf[crest_cols].to_excel(
                    writer, sheet_name="CREST_vs_DFT", index=False
                )

            # Hoja multi-T Boltzmann
            boltz_cols = ["conformer"] + [c for c in df_conf.columns
                                           if c.startswith("boltz_")]
            if len(boltz_cols) > 1:
                df_conf[boltz_cols].to_excel(
                    writer, sheet_name="Boltzmann_multiT", index=False
                )

        print(f"  Excel guardado en: {xlsx_path}")


def export_adc2_mol(df_adc2: "pd.DataFrame", mol_name: str, fmt: str):
    """Exporta excitaciones ADC2 de una molecula al CSV/XLSX del estado."""
    if df_adc2.empty:
        return
    mol_dir = RC_WORK_ROOT / mol_name
    if fmt in ("csv", "both"):
        out = mol_dir / "adc2_excitations.csv"
        df_adc2.to_csv(out, index=False)
        print(f"  ADC2 CSV -> {out}")
    if fmt in ("xlsx", "both"):
        xlsx_path = mol_dir / f"{mol_name}_results.xlsx"
        mode = "a" if xlsx_path.exists() else "w"
        kw   = {"if_sheet_exists": "replace"} if mode == "a" else {}
        try:
            with pd.ExcelWriter(xlsx_path, engine="openpyxl", mode=mode, **kw) as w:
                df_adc2.to_excel(w, sheet_name="ADC2", index=False)
            print(f"  ADC2 -> {xlsx_path} (hoja: ADC2)")
        except Exception as exc:
            warnings.warn(f"[ADC2 export] {exc}")


def export_comparison(df_cmp: pd.DataFrame, fmt: str):
    if df_cmp.empty:
        return
    out = RC_WORK_ROOT / "comparison_protonation_states"
    if fmt in ("csv", "both"):
        df_cmp.to_csv(out.with_suffix(".csv"), index=False)
        print(f"  CSV comparativo guardado en: {out.with_suffix('.csv')}")
    if fmt in ("xlsx", "both"):
        df_cmp.to_excel(out.with_suffix(".xlsx"), index=False)
        print(f"  Excel comparativo guardado en: {out.with_suffix('.xlsx')}")


# ─────────────────────────────────────────────
#  PUNTO DE ENTRADA
# ─────────────────────────────────────────────

def _main_collect():
    parser = argparse.ArgumentParser(
        description="Recolecta y analiza resultados ORCA: energías, "
                    "Boltzmann multi-T, reordenamiento CREST→DFT, "
                    "y comparación entre estados de protonación.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  python results_collector.py --mol MGLYVAL --state D
  python results_collector.py --mol MGLYVAL --state D --export xlsx --plot
  python results_collector.py --mol MGLYVAL --state D --goodvibes
  python results_collector.py --mol MGLYVAL --state D --goodvibes --export csv
  python results_collector.py --mol MGLYVAL --state D --esp --export xlsx
  python results_collector.py --mol MGLYVAL --compare-states
  python results_collector.py --mol MGLYVAL --compare-states --plot --ref N
  python results_collector.py --mol MGLYVAL --state D --project-root /home/kmilozeta/PAPER1
        """,
    )

    parser.add_argument("--mol", metavar="MOL", required=True,
                        help="nombre de la molecula (carpeta bajo "
                             "PROJECT_ROOT/<MOL>/runs/), p.ej. MGLYVAL")
    parser.add_argument("--project-root", metavar="PATH", default=None,
                        help=f"raiz del proyecto donde vive <MOL>/runs/ "
                             f"(default: {PROJECT_ROOT})")

    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--state",   metavar="NAME",
                      help="analizar un solo estado de protonacion "
                           "(nombre del directorio en <MOL>/runs/), p.ej. D")
    mode.add_argument("--compare-states", action="store_true",
                      help="comparar todos los estados de protonacion "
                           "disponibles en <MOL>/runs/")

    parser.add_argument("--ref",    metavar="NAME", default=None,
                        help="estado de referencia para dE relativa (solo con --compare-states)")
    parser.add_argument("--export", choices=["csv", "xlsx", "both"], default=None,
                        help="exportar tablas")
    parser.add_argument("--plot",   action="store_true",
                        help="generar espectro(s) UV-Vis")
    parser.add_argument("--goodvibes", action="store_true",
                        help="correr GoodVibes (mRRHO Grimme) sobre cada "
                             "freq.out y comparar contra el G_corr de ORCA "
                             "(solo con --state). Imprime y exporta el apendice "
                             "de comparacion a <MOL>/runs/<state>/.")
    parser.add_argument("--adc2",   action="store_true",
                        help="extraer lambda_max ADC(2) y comparar con TD-DFT. "
                             "Con --state: analiza ese estado. "
                             "Con --compare-states: compara todos los estados y "
                             "evalua la sistematicidad del desplazamiento "
                             "TDDFT-ADC2.")
    parser.add_argument("--esp",    action="store_true",
                        help="correr analisis ESP de superficie molecular via "
                             "Multiwfn (orca_2mkl -> molden -> Multiwfn) sobre "
                             "cada conf_* de 05_tddft/ (solo con --state). "
                             "Imprime y exporta el resumen ESP_min/ESP_max/MPI "
                             "a <MOL>/runs/<state>/.")
    parser.add_argument("--sigma",  type=float, default=20.0,
                        help="anchura Gaussiana en nm para los plots (default: 20)")
    parser.add_argument("--xmin",   type=float, default=200.0,
                        help="limite inferior del eje lambda en nm (default: 200)")
    parser.add_argument("--xmax",   type=float, default=500.0,
                        help="limite superior del eje lambda en nm (default: 500)")
    args = parser.parse_args()

    # ── Resolver WORK_ROOT = PROJECT_ROOT / <MOL> / "runs" ──────────────────
    global RC_WORK_ROOT
    project_root = Path(args.project_root) if args.project_root else PROJECT_ROOT
    RC_WORK_ROOT = project_root / args.mol / "runs"

    if not RC_WORK_ROOT.exists():
        parser.error(
            f"No se encontro la carpeta de runs para la molecula '{args.mol}': "
            f"{RC_WORK_ROOT}\n"
            f"Verifica --mol y/o --project-root (default: {PROJECT_ROOT})."
        )

    x_range = (args.xmin, args.xmax)

    if args.state:
        df_conf, df_tddft = collect_results(args.state)
        print_summary(df_conf, df_tddft, args.state)

        if args.adc2:
            df_adc2 = collect_adc2_results(args.state)
            print_adc2_summary(df_adc2, args.state, df_conf)
            if args.export:
                export_adc2_mol(df_adc2, args.state, args.export)
            if not df_adc2.empty:
                df_cmp_adc2 = compare_adc2_tddft(mols=[args.state])
                print_adc2_tddft_comparison(df_cmp_adc2)
                if args.export:
                    export_adc2_comparison(df_cmp_adc2, args.export)

        if args.export:
            rc_export_results(df_conf, df_tddft, args.state, args.export)
        if args.plot:
            plot_spectrum(df_tddft, df_conf, args.state,
                          sigma=args.sigma, x_range=x_range)
        if args.goodvibes:
            print(f"\n{chr(9472)*74}")
            print(f"  GOODVIBES (mRRHO Grimme) vs ORCA (RRHO) -- {args.mol}/{args.state}")
            print(f"{chr(9472)*74}")
            gv_results = run_goodvibes_on_mol(args.state)
            df_gv = compute_gv_appendix(args.state, df_conf, gv_results)
            print(df_gv.to_string(index=False))
            if args.export:
                export_gv_appendix(df_gv, args.state, args.export)

        if args.esp:
            esp_results = run_multiwfn_esp_on_mol(args.state)
            df_esp = compute_esp_summary_df(df_conf, esp_results)
            print_esp_summary(df_esp, args.state)
            if args.export:
                export_esp_results(df_esp, args.state, args.export)

    elif args.compare_states:
        result = compare_protonation_states(ref_mol=args.ref)
        if isinstance(result, tuple):
            df_cmp, dfs_conf, dfs_tddft = result
        else:
            df_cmp = result
            dfs_conf = dfs_tddft = {}

        print_comparison(df_cmp, ref_mol=args.ref)

        for mol in sorted(dfs_conf.keys()):
            print_summary(dfs_conf[mol], dfs_tddft.get(mol, pd.DataFrame()), mol)

        if args.adc2:
            mols_available = sorted(dfs_conf.keys()) or None
            df_cmp_adc2 = compare_adc2_tddft(mols=mols_available)
            print_adc2_tddft_comparison(df_cmp_adc2)
            for mol in sorted(dfs_conf.keys()):
                df_c   = dfs_conf.get(mol, pd.DataFrame())
                df_adc2 = collect_adc2_results(mol)
                print_adc2_summary(df_adc2, mol, df_c)
                if args.export and not df_adc2.empty:
                    export_adc2_mol(df_adc2, mol, args.export)
            if args.export:
                export_adc2_comparison(df_cmp_adc2, args.export)
            if args.plot:
                plot_adc2_tddft_comparison(
                    df_cmp_adc2, dfs_conf, dfs_tddft,
                    sigma=args.sigma, x_range=x_range
                )

        if args.export:
            export_comparison(df_cmp, args.export)
            for mol in sorted(dfs_conf.keys()):
                if not dfs_conf[mol].empty:
                    rc_export_results(dfs_conf[mol],
                                   dfs_tddft.get(mol, pd.DataFrame()),
                                   mol, args.export)
        if args.plot and not args.adc2:
            plot_comparison(dfs_conf, dfs_tddft,
                            sigma=args.sigma, x_range=x_range)


# ══════════════════════════════════════════════════════════════════════════════
#  §  PROTONATION
#     Ciclo termodinámico de pKa, cadena de 5 estados y Monte Carlo
# ══════════════════════════════════════════════════════════════════════════════

PC_WORK_ROOT        = Path.cwd() / "runs"  # default; sobreescribible con --work-root

# Temperaturas de análisis (contexto FAE antártico)

# G(H⁺, aq) — DOS CONVENCIONES DISTINTAS, NO INTERCAMBIABLES:
#
# 1) Para el CICLO DE BORN (gas + ΔG_solv separados):
#      G(H⁺,gas) = -6.28 kcal/mol      [Tissandier et al. 1998]
#      ΔG_solv(H⁺) = -265.9 kcal/mol   [Kelly, Cramer & Truhlar 2006]
#      G(H⁺,aq) implícito = -272.18 kcal/mol
#    Esta convención está calibrada para combinarse con E_SP en FASE GAS
#    pura (sin CPCM). Úsala solo si G_gas_acid/G_gas_base vienen de SP sin
#    solvente implícito.
#
# 2) Para el MÉTODO DIRECTO (solo energías ya en CPCM, sin SP fase gas):
#    Las energías E_SP_CPCM(X) ya incluyen el efecto de solvatación de la
#    molécula vía el modelo de continuo. Usar G(H⁺,aq) de Tissandier/Kelly
#    aquí sobreestima la estabilización del protón en ~8-10 kcal/mol, porque
#    esa constante fue derivada para cerrar un ciclo gas+solv explícito, no
#    para emparejarse con energías que ya son "solución implícita".
#    El valor correcto para esquemas de fase condensada / esquema directo
#    es el de Pliego & Riveros, J. Phys. Chem. A 2002, 106, 7434:
#      G(H⁺,aq) = -264.61 kcal/mol  (escala molar, 298.15 K, calibrado con
#                                     ciclos de fase condensada, consistente
#                                     con Zhan & Dixon 2001 ≈ -262.4 kcal/mol)

# G(H⁺, aq) PARA EL CICLO DE BORN (gas + solv explícitos, NO usar con el método directo)
G_HPLUS_AQ_BORN   = G_HPLUS_GAS + DG_SOLV_HPLUS   # kcal/mol  ≈ -272.18

# G(H⁺, aq) PARA EL MÉTODO DIRECTO (energías ya en CPCM) — Pliego & Riveros 2002
G_HPLUS_AQ_DIRECT = -264.61   # kcal/mol   [Pliego & Riveros, JPCA 2002, 106, 7434]

# Alias retro-compatible: NO usar en código nuevo, mantenido solo por si algo
# externo importa G_HPLUS_AQ directamente.
G_HPLUS_AQ = G_HPLUS_AQ_DIRECT

# Subdirectorios dentro de cada mol en runs/
SUBDIR_SP      = "04_sp"
SUBDIR_FREQ    = "03_numfreq"
SUBDIR_SP_GAS  = "04b_sp_gas"  # SPs sin CPCM para el ciclo de pKa
SUBDIR_SP_SMD  = "04c_sp_smd"  # SPs con SMD(water) en vez de CPCM(water) [--spsmd]

# ─────────────────────────────────────────────
#  GOODVIBES  (corrección cuasi-armónica opcional)
# ─────────────────────────────────────────────
# G_corr estándar (parse_thermochem) usa el oscilador armónico puro para
# TODOS los modos vibracionales, incluidos los rotores internos de baja
# frecuencia (<100 cm⁻¹) — típicos en M-Gly:Val (cadenas laterales Gly/Val,
# grupos OH/NH libres). Ahí la aproximación armónica sobreestima la entropía
# y por ende puede desviar bastante el G_corr. GoodVibes (Grimme, Head-Gordon)
# corrige esto con un tratamiento cuasi-armónico (mRRHO) para esos modos.
# Estos valores son sobreescritos por --gv-cutoff/--gv-qh/--gv-qs en main().
GV_CUTOFF = 100.0      # cm⁻¹, frequency cutoff (default de GoodVibes)
GV_QH     = False      # aplicar también la corrección cuasi-armónica de
                        # Head-Gordon a la entalpía (si no, H se deja armónica)

# Entorno conda donde está instalado GoodVibes. Igual que en results_collector.py:
# si se deja en None, se salta el paso 1 de _resolve_gv_python() y se prueba
# directamente con sys.executable / PATH.

# Caché del intérprete resuelto (evita repetir la búsqueda en disco + los
# subprocess de verificación en cada llamada a run_goodvibes(), que se invoca
# una vez por confórmero — potencialmente decenas de veces por corrida).
_GV_PYTHON_CACHE: Optional[str] = None

# pH de referencia oceánico antártico (WAP superficial, verano)
PH_ANTARCTIC   = 8.1
PH_STANDARD    = 7.0


# ═════════════════════════════════════════════════════════════════════════
#  CORRECCIÓN EMPÍRICA  (calibración externa: Gly / Val / Asp / Glu / AP)
# ═════════════════════════════════════════════════════════════════════════
#
# El ciclo de Born/CPCM tiene un error sistemático que NO es uniforme entre
# tipos de grupo (Ho & Coote 2010). Por eso la corrección se calibra POR TIPO
# DE TRANSICIÓN, no por posición en la cadena. El "slot" (pKa0/pKa1/pKa2) es
# solo la posición de la transición dentro de la cadena A→[R/D]→B→C de la
# molécula que se corre; el mismo tipo químico cae en slots distintos según
# la molécula (el NH₃⁺ de Gly/Val es el slot pKa1; el de Glu es el pKa2).
#
#   tipo           transición                       carga    slot natural  calibradores
#   carboxyl       1er COOH → COO⁻  (A → R/D)       +1 → 0   pKa0          Gly, Val, Asp(α), Glu(α)
#   carboxyl_2     2º  COOH → COO⁻  (R/D → B)        0 → −1  pKa1          Asp(β), Glu(γ)
#   amino          NH₃⁺ → NH₂, sin carga vecina      0 → −1  pKa1          Gly, Val
#   amino_dianion  NH₃⁺ → NH₂ con 2×COO⁻ (B → C)    −1 → −2  pKa2          Glu
#   pyridinium     N–H⁺ aromático → N (par APA→APR)  +1 → 0  pKa2 / pair   AP (4-aminopiridina)
#
# Correcciones respecto de la versión anterior (verificadas contra la química
# de los .xyz y contra la tabla de McMurry/Carey, no contra este archivo):
#   · Asp NO calibra 'amino': su paso R/D→B (AD→AB) es el β-COOH, es decir
#     'carboxyl_2'. (El comentario previo lo trataba como NH₃⁺.)
#   · 'carboxyl' (+1→0) y 'carboxyl_2' (0→−1) son clases de carga distintas y
#     no se mezclan: el 2º carboxilo se desprotona con un COO⁻ ya presente.
#   · pKa exp. de Glu (α-NH₃⁺) = 9.67 (antes 9.60).
#   · RMSE con n=1 es indefinido (antes se reportaba 0.000).
#   · modo 'auto': pendiente libre solo si n≥3 Y el rango de pKa calculado
#     cubre ≥1.0 unidad; si no, desplazamiento puro (ver fit_empirical_correction).
#
# Los pKa calculados vienen de las corridas de calibración. Para poblar/
# actualizar esta tabla sin editar el script:
#     maaqc.py calibrate --base <dir> --out calibration.json
#     maaqc.py pka ... --calibration-file calibration.json --empirical-correction ...
# Entradas con calc_*=None están "pendientes de calcular": el ajuste las ignora.
#
# IMPORTANTE — alcance: la corrección vale para el MISMO tipo de grupo al
# MISMO nivel de teoría, método (born/direct) y T. Verifica que la transición
# del MAA realmente sea de ese tipo antes de aplicarla; si no, repórtalo como
# limitación. Para pyridinium (AP, +1→0, sin carga vecina) transferido a un
# B→C (−1→−2) vale lo dicho abajo para 'amino_dianion': el apantallamiento por
# agua explícita de dos COO⁻ no lo captura un N–H⁺ aislado.
#
# 'amino_dianion' — CASO APARTE, NO es una variante de 'amino':
#   El tipo 'amino' de arriba (Gly/Val, n=2) calibra la desprotonación de un
#   α-NH₃⁺ que NO tiene ninguna otra carga formal cercana en la molécula
#   (la especie de partida, R/D, es neutra). Esa calibración no transfiere a
#   una desprotonación de NH₃⁺ que ocurre con DOS grupos carboxilato ya
#   formados en la misma molécula (especie de partida con carga -1, producto
#   con carga -2) — la electrostática de campo cercano y sobre todo el
#   apantallamiento de primera esfera de agua que le falta al continuo CPCM
#   son cualitativamente distintos en ese entorno.
#
#   Esto se verificó corriendo el ciclo completo de ácido glutámico (GA/GR/
#   GD/GB/GC), que tiene exactamente esa topología (GB→GC: NH₃⁺ flanqueado
#   por 2×COO⁻ → NH₂ + 2×COO⁻, igual que el paso B→C de M-Gly:Val). El
#   pKa₂ calculado con el método 'born' (17.89 armónico / 17.89 goodvibes,
#   ambas pasadas dan prácticamente el mismo valor) queda ~8 unidades por
#   encima del pKa₂ experimental de Glu (9.67, McMurry/Carey), un error MUY superior al ±1-2 típico del método y
#   a la corrección 'amino' de Gly/Val (que solo cubre ~2.4-2.5 unidades).
#   Esto descarta que el problema sea específico de M-Gly:Val: es un sesgo
#   sistemático del ciclo Born/CPCM para transiciones NH₃⁺→NH₂ con doble
#   carga negativa vecina, reproducido de forma casi idéntica en dos
#   moléculas independientes (M-Gly:Val cruda: pKa₂=18.02; Glu cruda:
#   pKa₂=17.89 — diferencia de solo 0.13 unidades entre sí, ambas ~8
#   unidades por encima de lo esperado). La hipótesis física más plausible
#   es que CPCM, al ser un continuo dieléctrico isotrópico, no resuelve el
#   apantallamiento explícito de agua de primera esfera alrededor de los dos
#   carboxilatos, que en realidad amortigua mucho del costo electrostático
#   de desprotonar el amonio vecino — un efecto de solvatación específica
#   que un continuo no puede capturar por construcción.
#
#   n=1 (un solo punto, Glu): no permite estimar RMSE de calibración: el
#   ajuste por desplazamiento puro (m=1, b=exp−calc) es exacto por
#   definición con un único punto. Se reporta igualmente porque un análogo
#   topológico exacto con n=1 es preferible a extrapolar 'amino' (n=2, sin
#   ninguna carga vecina) fuera de su dominio de calibración — pero el
#   resultado corregido debe presentarse con esta salvedad explícita en el
#   paper, no como un número con la misma confianza que 'carboxyl' (n=3).
_REF_AA = "McMurry Table 26-1 / Carey & Giuliano (consistente con CRC)"
_REF_AP = "4-aminopiridina, 9.17 a 25 °C (dato del proyecto)"

EMPIRICAL_CALIBRATION: dict[str, list[dict]] = {
    "carboxyl": [
        {"mol": "Glycine",  "calc_harmonic": 2.50, "calc_goodvibes": 2.45,
         "exp": 2.34, "exp_ref": _REF_AA},
        {"mol": "Valine",   "calc_harmonic": 2.38, "calc_goodvibes": 2.37,
         "exp": 2.32, "exp_ref": _REF_AA},
        {"mol": "Aspartic", "calc_harmonic": 2.09, "calc_goodvibes": 2.01,
         "exp": 1.88, "exp_ref": _REF_AA + " (α-COOH)"},
        {"mol": "Glutamic", "calc_harmonic": None, "calc_goodvibes": None,
         "exp": 2.19, "exp_ref": _REF_AA + " (α-COOH)"},
    ],
    "carboxyl_2": [
        {"mol": "Aspartic", "calc_harmonic": None, "calc_goodvibes": None,
         "exp": 3.65, "exp_ref": _REF_AA + " (β-COOH)"},
        {"mol": "Glutamic", "calc_harmonic": None, "calc_goodvibes": None,
         "exp": 4.25, "exp_ref": _REF_AA + " (γ-COOH)"},
    ],
    "amino": [
        {"mol": "Glycine",  "calc_harmonic": 11.95, "calc_goodvibes": 11.91,
         "exp": 9.60, "exp_ref": _REF_AA},
        {"mol": "Valine",   "calc_harmonic": 12.18, "calc_goodvibes": 12.19,
         "exp": 9.62, "exp_ref": _REF_AA},
    ],
    "amino_dianion": [
        {"mol": "Glutamic", "calc_harmonic": 17.89, "calc_goodvibes": 17.89,
         "exp": 9.67, "exp_ref": _REF_AA + " (α-NH₃⁺ de Glu)"},
    ],
    "pyridinium": [
        {"mol": "AP", "calc_harmonic": None, "calc_goodvibes": None,
         "exp": 9.17, "exp_ref": _REF_AP},
    ],
}

# Metadatos por tipo: slot(s) donde tiene sentido aplicarlo y clase de carga.
# Solo se usan para validar/avisar en --empirical-correction; no cambian números.
EMPIRICAL_TYPE_META: dict[str, dict] = {
    "carboxyl":      {"slots": ("pKa0",),          "dq": "+1→0",  "grupo": "1er COOH (A→R/D)"},
    "carboxyl_2":    {"slots": ("pKa1",),          "dq": "0→−1",  "grupo": "2º COOH (R/D→B)"},
    "amino":         {"slots": ("pKa1",),          "dq": "0→−1",  "grupo": "NH₃⁺ sin carga vecina"},
    "amino_dianion": {"slots": ("pKa2",),          "dq": "−1→−2", "grupo": "NH₃⁺ con 2×COO⁻"},
    "pyridinium":    {"slots": ("pKa2", "pair"),   "dq": "+1→0",  "grupo": "N–H⁺ aromático → N"},
}
EMPIRICAL_VALID_SLOTS = ("pKa0", "pKa1", "pKa2", "pair")

# Metadatos de la última calibración cargada con --calibration-file (si hay).
CALIBRATION_META: dict = {}

def _ols_fit(x: list[float], y: list[float]) -> dict:
    """
    Regresión lineal simple y = m·x + b por mínimos cuadrados (sin
    dependencia de scipy — solo numpy). Retorna pendiente, intercepto, R² y
    RMSE de los residuos. ok=False si hay <2 puntos o varianza nula en x.
    """
    xa = np.asarray(x, dtype=float)
    ya = np.asarray(y, dtype=float)
    n = len(xa)
    if n < 2:
        return {"ok": False, "n": n}
    mean_x, mean_y = xa.mean(), ya.mean()
    sxx = float(np.sum((xa - mean_x) ** 2))
    if sxx == 0:
        return {"ok": False, "n": n}
    sxy = float(np.sum((xa - mean_x) * (ya - mean_y)))
    m = sxy / sxx
    b = mean_y - m * mean_x
    y_pred = m * xa + b
    resid = ya - y_pred
    ss_res = float(np.sum(resid ** 2))
    ss_tot = float(np.sum((ya - mean_y) ** 2))
    r2 = (1 - ss_res / ss_tot) if ss_tot > 0 else float("nan")
    rmse = float(np.sqrt(ss_res / n))
    return {"ok": True, "n": n, "slope": float(m), "intercept": float(b),
            "r2": float(r2), "rmse": rmse, "residuals": resid.tolist()}


def fit_empirical_correction(group_type: str, gcorr_source: str = "goodvibes",
                             mode: str = "auto") -> dict:
    """
    Ajusta la corrección empírica pKa_exp ≈ m·pKa_calc + b para un tipo de
    transición (claves de EMPIRICAL_CALIBRATION), a partir de los puntos con
    pKa calculado disponible para ese gcorr_source (los None se ignoran).

    mode:
      'auto'   (default) — regresión lineal libre SOLO si hay ≥3 puntos Y el
               rango de pKa calculado cubre ≥1.0 unidad (si no, la pendiente
               no es estimable y extrapolaría mal a otros pKa). En cualquier
               otro caso, desplazamiento puro (m=1, b=sesgo medio).
      'linear' — fuerza regresión lineal (≥2 puntos; con 2 el ajuste es exacto).
      'offset' — fuerza desplazamiento puro con ≥1 punto.

    Con n=1 el RMSE es indefinido (NaN): el ajuste es exacto por construcción.
    Con n≥2 en modo offset se reporta además loo_rmse (leave-one-out).
    """
    calib = EMPIRICAL_CALIBRATION.get(group_type)
    if not calib:
        return {"ok": False,
                "reason": f"tipo de grupo desconocido: '{group_type}' "
                          f"(válidos: {list(EMPIRICAL_CALIBRATION)})"}

    calc_key = "calc_goodvibes" if gcorr_source == "goodvibes" else "calc_harmonic"
    pts = [(pt[calc_key], pt["exp"], pt["mol"]) for pt in calib
           if pt.get(calc_key) is not None and pt.get("exp") is not None]
    n = len(pts)
    if n < 1:
        pend = [pt["mol"] for pt in calib]
        return {"ok": False,
                "reason": f"'{group_type}' sin pKa calculado para "
                          f"gcorr_source='{gcorr_source}' (pendientes: "
                          f"{', '.join(pend)}). Corre `maaqc.py calibrate` y "
                          f"pasa --calibration-file."}

    x = [p[0] for p in pts]
    y = [p[1] for p in pts]
    mols = [p[2] for p in pts]
    warn: list[str] = []
    x_range = max(x) - min(x)

    if mode == "auto":
        use_mode = "linear" if (n >= 3 and x_range >= 1.0) else "offset"
        if n >= 3 and use_mode == "offset":
            warn.append(f"n={n} pero el rango de pKa calculado es solo "
                        f"{x_range:.2f} unidades (<1.0): la pendiente no es "
                        f"estimable; se usa desplazamiento puro.")
    else:
        use_mode = mode

    if use_mode == "linear" and n < 2:
        warn.append("Se pidió regresión lineal pero solo hay 1 punto de "
                    "calibración disponible; se usa desplazamiento puro.")
        use_mode = "offset"

    loo_rmse = float("nan")
    if use_mode == "linear":
        fit = _ols_fit(x, y)
        slope, intercept, r2, rmse = fit["slope"], fit["intercept"], fit["r2"], fit["rmse"]
        if n == 2:
            warn.append("Regresión lineal con solo 2 puntos: el ajuste es EXACTO "
                        "(0 grados de libertad), R² no informativo y pendiente "
                        "no verificada. Considera mode='offset'.")
    else:  # 'offset'
        bias = [yi - xi for xi, yi in zip(x, y)]
        intercept = float(np.mean(bias))
        slope = 1.0
        r2 = float("nan")
        if n >= 2:
            resid = [yi - (xi + intercept) for xi, yi in zip(x, y)]
            rmse = float(np.sqrt(np.mean(np.square(resid))))
            loo = []
            for k in range(n):
                b_k = float(np.mean([bias[j] for j in range(n) if j != k]))
                loo.append(y[k] - (x[k] + b_k))
            loo_rmse = float(np.sqrt(np.mean(np.square(loo))))
        else:
            rmse = float("nan")

    if n < 3:
        warn.append(f"Calibración con solo n={n} punto(s) para '{group_type}'. "
                    f"La incertidumbre de esta corrección es grande y se suma "
                    f"(no reemplaza) a la incertidumbre típica ±1-2 unidades "
                    f"de pKa del método Born/CPCM.")
    if n == 1:
        warn.append("n=1: el ajuste es exacto por construcción; el RMSE de "
                    "calibración es indefinido y no hay validación cruzada.")

    return {"ok": True, "group_type": group_type, "gcorr_source": gcorr_source,
            "mode": use_mode, "n": n, "molecules": mols,
            "slope": slope, "intercept": intercept, "r2": r2, "rmse": rmse,
            "loo_rmse": loo_rmse,
            "calib_x": x, "calib_y": y, "warnings": warn}


def apply_empirical_correction(pka_value: Optional[float], group_type: str,
                               gcorr_source: str = "goodvibes",
                               mode: str = "auto") -> dict:
    """
    Aplica fit_empirical_correction() a un valor de pKa crudo (Born/CPCM).
    Retorna dict con pKa_raw, pKa_corrected, delta, y los metadatos del
    ajuste; ok=False (con 'reason') si no se puede aplicar (pKa None o sin
    calibración disponible para ese tipo/gcorr_source).
    """
    if pka_value is None:
        return {"ok": False, "reason": "pKa crudo no disponible (None)",
                "group_type": group_type}
    fit = fit_empirical_correction(group_type, gcorr_source=gcorr_source, mode=mode)
    if not fit["ok"]:
        return {"ok": False, "reason": fit.get("reason", "ajuste no disponible"),
                "group_type": group_type}
    corrected = fit["slope"] * pka_value + fit["intercept"]
    out = {"ok": True, "pKa_raw": pka_value, "pKa_corrected": corrected,
           "delta": corrected - pka_value}
    out.update({k: fit[k] for k in
                ("group_type", "gcorr_source", "mode", "n", "molecules",
                 "slope", "intercept", "r2", "rmse", "loo_rmse", "warnings")})
    return out


def _print_empirical_correction(corr: Optional[dict], label: str):
    """Imprime el bloque de corrección empírica para un pKa, si aplica."""
    if not corr:
        return
    sep = "  " + "─" * 70
    if not corr.get("ok"):
        print(f"\n{sep}")
        print(f"  CORRECCIÓN EMPÍRICA  —  {label}  (tipo: {corr.get('group_type', '?')})")
        print(f"    ⚠  NO aplicada: {corr.get('reason', 'sin detalle')}")
        print(sep)
        return
    print(f"\n{sep}")
    print(f"  CORRECCIÓN EMPÍRICA  —  {label}  (tipo: {corr['group_type']})")
    print(sep)
    mode_txt = {
        "linear": "regresión lineal  (pKa_exp ≈ m·pKa_calc + b)",
        "offset": "desplazamiento puro  (pKa_exp ≈ pKa_calc + b, m=1)",
    }[corr["mode"]]
    print(f"  Calibración: {', '.join(corr['molecules'])}  "
          f"(n={corr['n']}, gcorr_source={corr['gcorr_source']})")
    print(f"  Ajuste: {mode_txt}")
    r2_txt = f"   R² = {corr['r2']:.4f}" if not np.isnan(corr["r2"]) else ""
    print(f"    m = {corr['slope']:.4f}   b = {corr['intercept']:+.4f}{r2_txt}")
    if np.isnan(corr["rmse"]):
        print("    RMSE (residuos de calibración) = n/a  (n=1: ajuste exacto)")
    else:
        print(f"    RMSE (residuos de calibración) = {corr['rmse']:.3f} unidades de pKa")
    if not np.isnan(corr.get("loo_rmse", float("nan"))):
        print(f"    RMSE leave-one-out             = {corr['loo_rmse']:.3f} unidades de pKa")
    for w in corr["warnings"]:
        print(f"    ⚠  {w}")
    print(f"\n    pKa (Born/CPCM, crudo)         = {corr['pKa_raw']:.2f}")
    print(f"    pKa (corregido empíricamente)  = {corr['pKa_corrected']:.2f}"
          f"   (Δ = {corr['delta']:+.2f})")
    print(sep)


def _parse_empirical_correction_spec(spec: Optional[str]) -> dict[str, str]:
    """
    Parsea el argumento --empirical-correction, formato
    'pKa0=carboxyl,pKa1=carboxyl_2,pKa2=pyridinium' (o 'pair=pyridinium' con
    --pair). Retorna {} si spec es None/vacío. Aborta con SystemExit si el
    formato, el slot o el tipo de grupo no son válidos. Si el tipo se aplica
    a un slot distinto de su(s) slot(s) natural(es) (EMPIRICAL_TYPE_META) solo
    se avisa por stderr: puede ser intencional (p. ej. pyridinium en pKa2).
    """
    if not spec:
        return {}
    out: dict[str, str] = {}
    for part in spec.split(","):
        part = part.strip()
        if not part:
            continue
        if "=" not in part:
            raise SystemExit(
                f"--empirical-correction: formato inválido en '{part}' "
                f"(se espera 'pKaN=tipo', p. ej. 'pKa0=carboxyl')")
        key, val = part.split("=", 1)
        key, val = key.strip(), val.strip().lower()
        if key not in EMPIRICAL_VALID_SLOTS:
            raise SystemExit(
                f"--empirical-correction: slot desconocido '{key}' en "
                f"'{part}'. Slots válidos: {list(EMPIRICAL_VALID_SLOTS)}")
        if val not in EMPIRICAL_CALIBRATION:
            raise SystemExit(
                f"--empirical-correction: tipo desconocido '{val}' en "
                f"'{part}'. Tipos válidos: {list(EMPIRICAL_CALIBRATION)}")
        meta = EMPIRICAL_TYPE_META.get(val)
        if meta and key not in meta["slots"]:
            print(f"AVISO: --empirical-correction: el tipo '{val}' "
                  f"({meta['grupo']}, {meta['dq']}) se calibró para el/los slot(s) "
                  f"{list(meta['slots'])}, pero lo estás aplicando a '{key}'. "
                  f"Verifica que la química de esa transición sea la misma.",
                  file=sys.stderr)
        out[key] = val
    return out


def load_calibration_file(path: Path) -> dict:
    """
    Carga un calibration.json (generado por `maaqc.py calibrate`) y REEMPLAZA
    por tipo las entradas de EMPIRICAL_CALIBRATION con las del archivo. Los
    tipos que el archivo no trae se dejan como están. Devuelve el bloque 'meta'.
    """
    import json
    with open(path, encoding="utf-8") as fh:
        blob = json.load(fh)
    cal = blob.get("EMPIRICAL_CALIBRATION")
    if not isinstance(cal, dict) or not cal:
        raise SystemExit(f"--calibration-file: '{path}' no contiene la clave "
                         f"'EMPIRICAL_CALIBRATION'.")
    for t, pts in cal.items():
        EMPIRICAL_CALIBRATION[t] = pts
    CALIBRATION_META.clear()
    CALIBRATION_META.update(blob.get("meta", {}))
    return CALIBRATION_META


# ─────────────────────────────────────────────
#  PARSERS  (extraídos de results_collector.py,
#            sin dependencias cruzadas)
# ─────────────────────────────────────────────

def parse_sp_energy(out_file: Path) -> Optional[float]:
    """
    Extrae la última 'FINAL SINGLE POINT ENERGY' del .out de ORCA.
    Retorna energía en Hartree, o None si no se encuentra.
    """
    if not out_file.exists():
        return None
    energy = None
    for line in out_file.read_text(errors="replace").splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                energy = float(line.split()[-1])
            except ValueError:
                pass
    return energy


def parse_theory_level(out_file: Path) -> Optional[str]:
    """
    Extrae la línea de keywords ORCA ('! funcional base CPCM ...') del .out,
    tal como aparece ecoizada al inicio del cálculo (línea '|  1> ! ...').

    Se usa solo para diagnóstico de consistencia entre estados — NO se usa
    para ningún cálculo numérico. Retorna None si no se encuentra.
    """
    if not out_file.exists():
        return None
    pat = re.compile(r"^\s*\|\s*1>\s*!\s*(.+)$")
    for line in out_file.read_text(errors="replace").splitlines():
        m = pat.match(line)
        if m:
            return m.group(1).strip()
    return None


def parse_charge_mult_natoms(out_file: Path) -> dict:
    """
    Extrae carga total, multiplicidad y número de átomos de un .out de ORCA.

    Se usa para validar que los estados de protonación pasados por
    --protonated/--neutral/--anion/--dianion/--zwitterion realmente difieren
    entre sí en la forma esperada (Δq=-1/ΔN_átomos=-1 para una desprotonación,
    Δq=0/ΔN_átomos=0 para un par de tautómeros) — ver
    _validate_state_charges(). Sin este chequeo, un mapeo de carpetas
    desalineado (p. ej. --anion apuntando al estado equivocado, o un
    best_conf que no corresponde a la geometría que se cree) se procesa sin
    quejarse y produce un pKa numéricamente válido pero físicamente absurdo.

    Prioriza la línea de eco del input ('|  N> * xyz Q M'), que es
    exactamente lo que se envió a ORCA. Si no aparece (p. ej. geometría
    leída de un archivo externo en vez de un bloque xyz inline), cae a los
    patrones de la sección de configuración SCF ('Total Charge ... Charge
    ... N', 'Multiplicity ... Mult ... N').

    n_atoms se cuenta a partir del primer bloque
    'CARTESIAN COORDINATES (ANGSTROEM)' del archivo.

    Retorna {"charge": int|None, "mult": int|None, "n_atoms": int|None}.
    """
    null = {"charge": None, "mult": None, "n_atoms": None}
    if not out_file.exists():
        return null
    try:
        lines = out_file.read_text(errors="replace").splitlines()
    except OSError:
        return null

    charge: Optional[int] = None
    mult: Optional[int] = None
    n_atoms: Optional[int] = None

    # 1) Eco del input: "|  N> * xyz Q M"
    pat_echo = re.compile(r"^\s*\|\s*\d+>\s*\*\s*xyz\s+(-?\d+)\s+(\d+)", re.IGNORECASE)
    for line in lines:
        m = pat_echo.match(line)
        if m:
            charge, mult = int(m.group(1)), int(m.group(2))
            break

    # 2) Fallback: bloque de configuración SCF ("Total Charge ... N")
    if charge is None:
        pat_q = re.compile(r"Total Charge\s+Charge\s+\.+\s*(-?\d+)", re.IGNORECASE)
        for line in lines:
            m = pat_q.search(line)
            if m:
                charge = int(m.group(1))
                break
    if mult is None:
        pat_m = re.compile(r"Multiplicity\s+Mult\s+\.+\s*(\d+)", re.IGNORECASE)
        for line in lines:
            m = pat_m.search(line)
            if m:
                mult = int(m.group(1))
                break

    # 3) Número de átomos: primer bloque CARTESIAN COORDINATES (ANGSTROEM)
    for idx, line in enumerate(lines):
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line.upper():
            j = idx + 2  # saltar encabezado y línea separadora '---...'
            count = 0
            while j < len(lines) and lines[j].strip():
                count += 1
                j += 1
            n_atoms = count if count > 0 else None
            break

    return {"charge": charge, "mult": mult, "n_atoms": n_atoms}


def _validate_state_charges(
    data_a: dict,
    data_b: dict,
    dq_expected: int,
    dn_expected: int,
    label_a: str,
    label_b: str,
) -> list[str]:
    """
    Verifica que data_b difiera de data_a exactamente en dq_expected unidades
    de carga y dn_expected átomos. Ejemplos:
      - compute_pka(acid, base):              dq_expected=-1, dn_expected=-1
      - compute_tautomer_equilibrium(R, D):    dq_expected=0,  dn_expected=0

    No aborta nada por sí sola: retorna una lista de mensajes. Los mensajes
    que empiezan con '❌' son fatales — el llamador (compute_pka /
    compute_tautomer_equilibrium) los trata como error y NO calcula el
    pKa/equilibrio, en vez de seguir de largo con un resultado numérico sin
    sentido físico. Los mensajes sin '❌' son informativos (p. ej. cuando no
    se pudo parsear la carga/N° de átomos de alguno de los dos estados) y
    solo se agregan a las advertencias.
    """
    msgs: list[str] = []
    qa, qb = data_a.get("charge"), data_b.get("charge")
    na, nb = data_a.get("n_atoms"), data_b.get("n_atoms")

    if qa is None or qb is None:
        msgs.append(
            f"Could not verify charge between '{label_a}' and '{label_b}' "
            f"(charge not determined from the ORCA .out for at least one of "
            f"the two states). Manually verify that the expected charge "
            f"relation (Δq={dq_expected:+d}) is correct."
        )
    elif (qb - qa) != dq_expected:
        msgs.append(
            f"❌ INCONSISTENT CHARGE: '{label_a}' (q={qa:+d}) → '{label_b}' "
            f"(q={qb:+d}) implies Δq={qb - qa:+d}, but Δq={dq_expected:+d} "
            f"was expected. Check whether the folder mapping "
            f"(--protonated/--neutral/--anion/--dianion/--zwitterion) "
            f"corresponds to the state you think it is."
        )

    if na is None or nb is None:
        msgs.append(
            f"Could not verify atom count between '{label_a}' and '{label_b}' "
            f"(atom count not determined from the ORCA .out for at least one "
            f"of the two states). Manually verify that the expected relation "
            f"(ΔN={dn_expected:+d}) is correct."
        )
    elif (nb - na) != dn_expected:
        msgs.append(
            f"❌ INCONSISTENT ATOM COUNT: '{label_a}' (N={na}) → "
            f"'{label_b}' (N={nb}) implies ΔN={nb - na:+d}, but ΔN={dn_expected:+d} "
            f"was expected. Check whether the folder mapping "
            f"corresponds to the state you think it is, or whether the "
            f"representative conformer (best_conf) corresponds to the "
            f"expected geometry."
        )

    return msgs


def parse_thermochem(freq_out: Path, T: float = T_STD) -> dict:
    """
    Extrae correcciones termoquímicas del bloque THERMOCHEMISTRY de ORCA 6.

    Lógica idéntica a results_collector.py:
      - 'Final entropy term' es T·S a T_parsed (en Eh), NO es S puro.
      - S_total = TS_total / T_parsed   (Eh/K)
      - G_corr(T) = H_corr − T · S_total

    Retorna dict con claves:
        E_elec_freq, ZPE, H_corr, TS_total, S_total,
        G_corr, G_total_freq, T_parsed, n_imag, ok
    Todos los valores en Hartree (salvo T_parsed en K y n_imag entero).
    """
    null = {
        "E_elec_freq": None, "ZPE": None,       "H_corr": None,
        "TS_total":    None, "S_total": None,    "G_corr": None,
        "G_total_freq": None, "T_parsed": None,
        "n_imag": 0,         "ok": False,
    }

    if not freq_out.exists():
        return null

    try:
        lines = freq_out.read_text(errors="replace").splitlines()
    except OSError:
        return null

    # Frecuencias imaginarias
    n_imag = sum(
        1 for l in lines
        if re.match(r"\s+\d+:\s+(-\d+\.\d+)\s+cm", l)
    )

    # Energía electrónica (primera ocurrencia = SP del NumFreq)
    E_elec_freq: Optional[float] = None
    for line in lines:
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                E_elec_freq = float(line.split()[-1])
                break
            except ValueError:
                pass

    # Bloque THERMOCHEMISTRY (última ocurrencia)
    thermo_start = -1
    for idx, line in enumerate(lines):
        if "THERMOCHEMISTRY AT" in line.upper():
            thermo_start = idx

    if thermo_start < 0:
        result = dict(null)
        result["E_elec_freq"] = E_elec_freq
        result["n_imag"] = n_imag
        return result

    # Temperatura del bloque
    T_parsed: Optional[float] = None
    m = re.search(r"THERMOCHEMISTRY AT\s+([\d.]+)\s*K",
                  lines[thermo_start], re.IGNORECASE)
    if m:
        try:
            T_parsed = float(m.group(1))
        except ValueError:
            pass

    # Patrones ORCA 6
    _pat_zpe = re.compile(r"Zero point energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_h   = re.compile(r"Total [Ee]nthalpy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_s   = re.compile(r"Final entropy term\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I)
    _pat_g   = re.compile(
        r"(?:Total|Final) Gibbs free energy\s+\.\.\.\s+([-\d.]+)\s*Eh", re.I
    )

    ZPE = H_total = TS_total = G_total_freq = None
    for line in lines[thermo_start:]:
        if ZPE is None:
            m = _pat_zpe.search(line)
            if m:
                try: ZPE = float(m.group(1))
                except ValueError: pass
        if H_total is None:
            m = _pat_h.search(line)
            if m:
                try: H_total = float(m.group(1))
                except ValueError: pass
        if TS_total is None:
            m = _pat_s.search(line)
            if m:
                try: TS_total = float(m.group(1))
                except ValueError: pass
        if G_total_freq is None:
            m = _pat_g.search(line)
            if m:
                try: G_total_freq = float(m.group(1))
                except ValueError: pass

    # Derivar correcciones
    H_corr  = (H_total - E_elec_freq) if (H_total is not None and E_elec_freq is not None) else None
    S_total = (TS_total / T_parsed)   if (TS_total is not None and T_parsed and T_parsed != 0) else None
    G_corr  = (H_corr - T * S_total)  if (H_corr is not None and S_total is not None) else None

    ok = all(x is not None for x in [ZPE, H_corr, TS_total, G_corr, G_total_freq])

    return {
        "E_elec_freq":  E_elec_freq,
        "ZPE":          ZPE,
        "H_corr":       H_corr,
        "TS_total":     TS_total,
        "S_total":      S_total,       # Eh/K
        "G_corr":       G_corr,        # Eh, a temperatura T
        "G_total_freq": G_total_freq,  # Eh
        "T_parsed":     T_parsed,
        "n_imag":       n_imag,
        "ok":           ok,
    }


def _resolve_gv_python(force: bool = False) -> str:
    """
    Devuelve la ruta absoluta al intérprete Python que tiene goodvibes instalado.

    Idéntico en espíritu a _resolve_gv_python() de results_collector.py (misma
    estrategia en el mismo orden de prioridad), para que ambos scripts se
    comporten igual sea cual sea el shell/entorno desde el que se invoquen:

      1. Si GV_CONDA_ENV está definido, buscar ese entorno en los prefijos
         estándar de miniforge3, miniconda3, anaconda3 y conda en $HOME.
      2. Intentar `python -c "import goodvibes"` con el ejecutable actual
         (sys.executable), por si el entorno activo ya lo incluye.
      3. Buscar `goodvibes` en el PATH como script de consola y derivar el
         Python de su mismo directorio bin/.
      4. Fallback: sys.executable (fallará más adelante con un warning claro
         en vez de un `shutil.which` silencioso).

    El resultado se cachea en _GV_PYTHON_CACHE porque run_goodvibes() se llama
    una vez por confórmero — resolver esto en disco + subprocess en cada
    llamada sería lento e innecesario.
    """
    global _GV_PYTHON_CACHE
    if _GV_PYTHON_CACHE is not None and not force:
        return _GV_PYTHON_CACHE

    if GV_CONDA_ENV:
        home = Path.home()
        for base in ["miniforge3", "miniconda3", "anaconda3", "conda"]:
            candidate = home / base / "envs" / GV_CONDA_ENV / "bin" / "python"
            if candidate.exists():
                try:
                    r = subprocess.run(
                        [str(candidate), "-c", "import goodvibes"],
                        capture_output=True, timeout=10,
                    )
                    if r.returncode == 0:
                        _GV_PYTHON_CACHE = str(candidate)
                        return _GV_PYTHON_CACHE
                except (subprocess.TimeoutExpired, OSError):
                    pass

    try:
        r = subprocess.run(
            [sys.executable, "-c", "import goodvibes"],
            capture_output=True, timeout=10,
        )
        if r.returncode == 0:
            _GV_PYTHON_CACHE = sys.executable
            return _GV_PYTHON_CACHE
    except (subprocess.TimeoutExpired, OSError):
        pass

    gv_script = shutil.which("goodvibes")
    if gv_script:
        py_candidate = Path(gv_script).parent / "python"
        if py_candidate.exists():
            _GV_PYTHON_CACHE = str(py_candidate)
            return _GV_PYTHON_CACHE

    _GV_PYTHON_CACHE = sys.executable
    return _GV_PYTHON_CACHE


def run_goodvibes(
    freq_out: Path,
    T: float = T_STD,
    cutoff: float = GV_CUTOFF,
    qs: str = GV_QS,
    qh: bool = GV_QH,
) -> dict:
    """
    Corre GoodVibes UNA vez sobre un freq.out de ORCA y extrae de la MISMA
    tabla dos juegos de valores paralelos:

      • "armónico puro"     → columnas T.S / G(T)       (nunca amortiguadas
                               por el cutoff — así las define GoodVibes; T.S
                               es el RRHO estándar, sin mRRHO)
      • "cuasi-armónico"    → columnas T.qh-S / qh-G(T)  (mRRHO de Grimme,
                               amortiguado con `cutoff`, default 100 cm⁻¹)

    Esto es clave para que la comparación "armónico vs. GoodVibes" del
    script sea limpia: antes, la pasada "armónica" salía de parse_thermochem()
    (bloque THERMOCHEMISTRY de ORCA, que ya viene con `Quasi RRHO = True` y
    su propio cutoff interno ~1 cm⁻¹) mientras que la pasada "GoodVibes"
    salía de un cálculo aparte — dos fuentes de datos distintas, no
    comparables limpiamente. Ahora ambas pasadas vienen del mismo cálculo de
    frecuencias, la misma detección de simetría, etc.; la única diferencia
    real entre ellas es si se aplica o no el amortiguamiento de Grimme.

    Invocación (idéntica a la de terminal):

        python -m goodvibes.GoodVibes <freq.out> \\
            --qs grimme --freq_cutoff 100 -t 298.15

    Se parsea la tabla de stdout (no `--json`/`--output_format`, cuyo nombre
    de archivo de salida varía entre versiones de GoodVibes) — más simple,
    algo más lento pero robusto a esos cambios de flags.

    Retorna dict con:
        scf_energy, ZPE, n_imag, ok, warning
        G_corr_harmonic, H_corr_harmonic, S_total_harmonic, G_total_harmonic
        G_corr_qh,       H_corr_qh,       S_total_qh,       G_total_qh
    más alias retro-compatibles G_corr/H_corr/S_total/qh_gibbs_free_energy
    (= los *_qh, que es lo que antes usaba la pasada "goodvibes").
    """
    null = {
        "G_corr": None, "H_corr": None, "S_total": None,
        "scf_energy": None, "ZPE": None, "qh_gibbs_free_energy": None,
        "G_corr_harmonic": None, "H_corr_harmonic": None,
        "S_total_harmonic": None, "G_total_harmonic": None,
        "G_corr_qh": None, "H_corr_qh": None,
        "S_total_qh": None, "G_total_qh": None,
        "n_imag": 0, "ok": False, "warning": None,
    }

    if not freq_out.exists():
        return {**null, "warning": f"No existe {freq_out}"}

    python_exe = _resolve_gv_python()

    try:
        r = subprocess.run(
            [python_exe, "-c", "import goodvibes"],
            capture_output=True, timeout=10,
        )
        if r.returncode != 0:
            return {**null, "warning": (
                f"No se puede importar goodvibes con {python_exe}. "
                f"Instalar con: conda activate {GV_CONDA_ENV or '<env>'} && "
                f"pip install goodvibes --break-system-packages"
            )}
    except (subprocess.TimeoutExpired, OSError) as e:
        return {**null, "warning": f"No se pudo verificar goodvibes con {python_exe}: {e}"}

    cmd = [
        python_exe, "-m", "goodvibes.GoodVibes", str(freq_out),
        "--qs", qs, "--freq_cutoff", str(cutoff), "-t", str(T),
    ]
    if qh:
        cmd.append("--qh")

    try:
        with tempfile.TemporaryDirectory() as tmp:
            proc = subprocess.run(
                cmd, cwd=tmp, capture_output=True, text=True, timeout=180,
            )
    except (subprocess.TimeoutExpired, OSError) as e:
        return {**null, "warning": (
            f"Fallo al ejecutar goodvibes ({python_exe} -m goodvibes.GoodVibes): {e}"
        )}

    if proc.returncode != 0:
        stderr_tail = (proc.stderr or "").strip()[-400:]
        return {**null, "warning": (
            f"goodvibes ({python_exe}) terminó con returncode={proc.returncode}. "
            f"{stderr_tail}"
        )}

    stdout = proc.stdout or ""
    parsed = _parse_goodvibes_table(stdout, freq_out.stem)
    if parsed is None:
        stdout_tail = stdout.strip()[-600:]
        return {**null, "warning": (
            f"No se pudo parsear la tabla de salida de goodvibes para "
            f"{freq_out.name}. Salida (cola):\n{stdout_tail}"
        )}

    scf_energy = parsed.get("E")
    ZPE        = parsed.get("ZPE")
    G_qh       = parsed.get("qh-G(T)")
    G_harm     = parsed.get("G(T)")
    if scf_energy is None or G_qh is None or G_harm is None:
        return {**null, "warning": (
            "goodvibes no reportó columnas E / G(T) / qh-G(T) en la tabla "
            "(freq.out posiblemente incompleto o mal formado)."
        )}

    # Entalpía armónica pura: siempre la columna H estándar.
    H_harm = parsed.get("H")
    H_corr_harmonic = (H_harm - scf_energy) if H_harm is not None else None

    # Entalpía cuasi-armónica: columna qh-H (Head-Gordon) si --qh está
    # activo; si no, GoodVibes no amortigua la entalpía y
    # H_corr_qh == H_corr_harmonic.
    H_qh = parsed.get("qh-H") if "qh-H" in parsed else H_harm
    H_corr_qh = (H_qh - scf_energy) if H_qh is not None else None

    # T.S / T.qh-S ya vienen multiplicadas por T; se recupera S dividiendo por T.
    T_S    = parsed.get("T.S")
    T_qh_S = parsed.get("T.qh-S")
    S_total_harmonic = (T_S / T) if T_S is not None else None
    S_total_qh       = (T_qh_S / T) if T_qh_S is not None else None

    n_imag = 0
    m = re.search(r"(\d+)\s+imaginary frequenc", stdout, re.IGNORECASE)
    if m:
        n_imag = int(m.group(1))

    G_corr_harmonic = G_harm - scf_energy
    G_corr_qh       = G_qh - scf_energy

    return {
        "scf_energy":         scf_energy,
        "ZPE":                ZPE,
        "n_imag":             n_imag,
        "ok":                 True,
        "warning":            None,
        "G_corr_harmonic":    G_corr_harmonic,
        "H_corr_harmonic":    H_corr_harmonic,
        "S_total_harmonic":   S_total_harmonic,
        "G_total_harmonic":   G_harm,
        "G_corr_qh":          G_corr_qh,
        "H_corr_qh":          H_corr_qh,
        "S_total_qh":         S_total_qh,
        "G_total_qh":         G_qh,
        # alias retro-compatibles (= *_qh, comportamiento previo de "goodvibes")
        "G_corr":               G_corr_qh,
        "H_corr":               H_corr_qh,
        "S_total":              S_total_qh,
        "qh_gibbs_free_energy": G_qh,
    }


def _parse_goodvibes_table(stdout: str, structure_stem: str) -> Optional[dict]:
    """
    Parsea la tabla de resultados que GoodVibes imprime por stdout, del
    estilo:

          Structure              E        ZPE              H        T.S     T.qh-S           G(T)        qh-G(T)
        ────────────────────────────────────────────────────────────────────────────────────────────────────────
             freq        -1220.971714   0.402267   -1220.542612   0.085076   0.079189   -1220.627688   -1220.621801
        ────────────────────────────────────────────────────────────────────────────────────────────────────────

    Devuelve un dict {nombre_columna: valor_float} para la fila cuya
    columna "Structure" coincide con `structure_stem` (nombre del archivo
    sin extensión; en la convención del pipeline siempre es "freq", ya que
    el archivo es conf_NNN/freq.out), o None si la tabla no aparece o no se
    pudo parsear.

    NOTA: GoodVibes antepone un marcador de estado de una sola columna
    (aparece como un bullet/checkmark, p. ej. una "o" suelta) ANTES del
    nombre de estructura en cada fila de datos:

        o     freq        -1220.967288   0.403345   ...

    Si se asumiera que la estructura es siempre el primer token, esa "o"
    correría el conteo de columnas en uno y la fila se descartaría por
    "no calzar" con el header — silenciosamente, sin fallar con un error
    claro. Por eso acá se busca `structure_stem` en la posición que sea
    dentro de los tokens, y se toman los valores inmediatamente después.
    """
    lines = stdout.splitlines()
    header_tokens = None
    header_idx = None
    for i, line in enumerate(lines):
        if "Structure" in line and "G(T)" in line:
            header_tokens = re.split(r"\s{2,}", line.strip())
            header_idx = i
            break
    if header_tokens is None:
        return None

    n_values = len(header_tokens) - 1  # columnas numéricas tras "Structure"

    for line in lines[header_idx + 1:]:
        stripped = line.strip()
        if not stripped or set(stripped) <= {"─", "-", "_"}:
            continue
        tokens = stripped.split()
        if structure_stem not in tokens:
            continue
        idx = tokens.index(structure_stem)
        value_tokens = tokens[idx + 1:]
        if len(value_tokens) != n_values:
            continue
        values: dict = {"Structure": structure_stem}
        for key, val in zip(header_tokens[1:], value_tokens):
            try:
                values[key] = float(val)
            except ValueError:
                values[key] = val
        return values

    return None


def get_thermochem(
    freq_out: Path,
    T: float = T_STD,
    source: str = "harmonic",
    gv_cutoff: float = GV_CUTOFF,
    gv_qs: str = GV_QS,
    gv_qh: bool = GV_QH,
) -> dict:
    """
    Wrapper único para obtener G_corr, con dos fuentes intercambiables que
    ahora comparten la MISMA corrida de GoodVibes sobre el mismo freq.out
    (ver run_goodvibes()) — solo difieren en qué columnas de esa tabla usan:

      source="harmonic"  → columnas T.S / G(T) de GoodVibes
                            (RRHO puro, NUNCA amortiguado por el cutoff)
      source="goodvibes" → columnas T.qh-S / qh-G(T) de GoodVibes
                            (mRRHO de Grimme, amortiguado con gv_cutoff)

    Antes, "harmonic" salía de parse_thermochem() (bloque THERMOCHEMISTRY
    de ORCA, que ya trae `Quasi RRHO = True` con su propio cutoff interno
    ~1 cm⁻¹) mientras que "goodvibes" salía de un cálculo aparte — dos
    fuentes de datos distintas, no comparables. parse_thermochem() se deja
    definida arriba solo como utilidad de referencia/inspección manual, ya
    no la usa esta función.

    Retorna siempre el shape de parse_thermochem() (E_elec_freq, ZPE,
    H_corr, TS_total, S_total, G_corr, G_total_freq, T_parsed, n_imag, ok),
    más "source" y "warning" para trazabilidad. Esto permite que
    collect_mol_data(), _state_G_total_direct() y compute_pka() usen la
    misma interfaz sin importar la fuente.
    """
    if source not in ("harmonic", "goodvibes"):
        raise ValueError(f"source debe ser 'harmonic' o 'goodvibes', no {source!r}")

    gv = run_goodvibes(freq_out, T=T, cutoff=gv_cutoff, qs=gv_qs, qh=gv_qh)

    if source == "harmonic":
        S_total  = gv.get("S_total_harmonic")
        TS_total = (T * S_total) if S_total is not None else None
        return {
            "E_elec_freq":  gv.get("scf_energy"),
            "ZPE":          gv.get("ZPE"),
            "H_corr":       gv.get("H_corr_harmonic"),
            "TS_total":     TS_total,
            "S_total":      S_total,
            "G_corr":       gv.get("G_corr_harmonic"),
            "G_total_freq": gv.get("G_total_harmonic"),
            "T_parsed":     T,
            "n_imag":       gv.get("n_imag", 0),
            "ok":           gv.get("ok", False),
            "source":       "harmonic",
            "warning":      gv.get("warning"),
        }

    # source == "goodvibes"
    S_total  = gv.get("S_total_qh")
    TS_total = (T * S_total) if S_total is not None else None
    return {
        "E_elec_freq":  gv.get("scf_energy"),
        "ZPE":          gv.get("ZPE"),
        "H_corr":       gv.get("H_corr_qh"),
        "TS_total":     TS_total,
        "S_total":      S_total,
        "G_corr":       gv.get("G_corr_qh"),
        "G_total_freq": gv.get("G_total_qh"),
        "T_parsed":     T,
        "n_imag":       gv.get("n_imag", 0),
        "ok":           gv.get("ok", False),
        "source":       "goodvibes",
        "warning":      gv.get("warning"),
    }


# ─────────────────────────────────────────────
#  COLECTOR POR MOL  (versión ligera, sin TDDFT
#  ni CREST reranking — solo lo necesario para
#  el ciclo de pKa)
# ─────────────────────────────────────────────

def collect_mol_data(mol_name: str, verbose: bool = False,
                      gcorr_source: str = "harmonic",
                      use_smd: bool = False) -> dict:
    """
    Recolecta los datos termodinámicos de todos los confórmeros de una molécula.

    Para cada confórmero con SP acuoso:
      - E_sp_aq   : energía SP "acuosa" (Hartree). Por default, CPCM(water)
                    de 04_sp/conf_NNN/sp.out. Si use_smd=True, en su lugar se
                    usa SMD(water) de 04c_sp_smd/conf_NNN/sp_smd.out (falta
                    silenciosamente a NaN si ese confórmero no tiene SP_SMD
                    calculado — igual que E_sp_gas con SP_GAS).
      - E_sp_gas  : energía SP sin CPCM (Hartree, de 04b_sp_gas — opcional)
      - G_corr(T) : corrección de Gibbs de NumFreq, re-evaluada a cada T (Hartree)
      - G_total(T): E_sp_aq + G_corr(T)    ciclo compuesto Ho & Coote

    Retorna dict con:
        mol_name        : str
        conformers      : list[dict]  — datos por confórmero
        df              : pd.DataFrame
        best_conf_298K  : str — id del confórmero con menor G_total a 298.15 K
        G_boltz(T)      : dict {label: G_boltzmann en Hartree} — ensemble efectivo
        ok              : bool
        warnings        : list[str]
        aq_source       : str — "CPCM" o "SMD", según use_smd
    """
    mol_dir = PC_WORK_ROOT / mol_name
    warns: list[str] = []

    if not mol_dir.exists():
        return {"mol_name": mol_name, "ok": False,
                "warnings": [f"Directory not found: {mol_dir}"],
                "conformers": [], "df": pd.DataFrame()}

    sp_dirs = sorted((mol_dir / SUBDIR_SP).glob("conf_*")) \
              if (mol_dir / SUBDIR_SP).exists() else []

    if not sp_dirs:
        return {"mol_name": mol_name, "ok": False,
                "warnings": [f"No SP conformers in {mol_dir / SUBDIR_SP}"],
                "conformers": [], "df": pd.DataFrame()}

    aq_source = "SMD" if use_smd else "CPCM"
    rows: list[dict] = []
    n_missing_smd = 0

    for sp_dir in sp_dirs:
        conf_id = sp_dir.name

        # SP acuoso — CPCM (default) o SMD (si use_smd=True). El nombre de
        # columna sigue siendo E_sp_aq en ambos casos: el resto del pipeline
        # (compute_pka, _state_G_total_direct, compute_tautomer_equilibrium)
        # no necesita saber cuál de los dos modelos de solvente se usó, solo
        # que "E_sp_aq" es la energía acuosa vigente para esta corrida.
        if use_smd:
            sp_aq_file = mol_dir / SUBDIR_SP_SMD / conf_id / "sp_smd.out"
            e_sp_aq    = parse_sp_energy(sp_aq_file)
            if e_sp_aq is None:
                n_missing_smd += 1
        else:
            sp_aq_file  = sp_dir / "sp.out"
            e_sp_aq     = parse_sp_energy(sp_aq_file)

        # SP fase gas (sin CPCM) — para el ciclo de pKa
        sp_gas_file = mol_dir / SUBDIR_SP_GAS / conf_id / "sp_gas.out"
        e_sp_gas    = parse_sp_energy(sp_gas_file)

        # Correcciones térmicas — leídas a T_STD, S_total re-escalado.
        # gcorr_source="harmonic" → columnas T.S/G(T) de GoodVibes (RRHO
        # puro, sin amortiguar); "goodvibes" → T.qh-S/qh-G(T) (mRRHO
        # cuasi-armónico). Ambas vienen de la misma corrida de GoodVibes
        # (ver get_thermochem()/run_goodvibes()), así que un fallo de
        # GoodVibes afecta a las dos fuentes por igual.
        freq_file = mol_dir / SUBDIR_FREQ / conf_id / "freq.out"
        tc = get_thermochem(freq_file, T=T_STD, source=gcorr_source,
                            gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)
        if not tc["ok"] and tc.get("warning"):
            warns.append(f"[{mol_name}/{conf_id}] GoodVibes ({gcorr_source}): {tc['warning']}")

        # Carga/multiplicidad/N° átomos — para validación cruzada entre
        # estados (ver _validate_state_charges). Se intenta primero desde el
        # SP acuoso (CPCM, siempre disponible como referencia estructural) y,
        # si falla, desde el freq.out del mismo confórmero.
        cmn = parse_charge_mult_natoms(sp_dir / "sp.out")
        if cmn["charge"] is None and cmn["n_atoms"] is None:
            cmn = parse_charge_mult_natoms(freq_file)

        row: dict = {
            "conformer":   conf_id,
            "E_sp_aq":     e_sp_aq,
            "E_sp_gas":    e_sp_gas,
            "thermo_ok":   tc["ok"],
            "n_imag":      tc["n_imag"],
            "H_corr":      tc["H_corr"],
            "S_total_EhK": tc["S_total"],
            "ZPE":         tc["ZPE"],
            "charge":      cmn["charge"],
            "mult":        cmn["mult"],
            "n_atoms":     cmn["n_atoms"],
        }

        # G_total compuesto a cada temperatura
        for T, lbl in zip(TEMPS_K, TEMP_LABELS):
            if e_sp_aq is not None and tc["ok"]:
                G_corr_T  = tc["H_corr"] - T * tc["S_total"]
                row[f"G_total_{lbl}"] = e_sp_aq + G_corr_T
            else:
                row[f"G_total_{lbl}"] = None

        rows.append(row)

        if verbose:
            gas_str = f"{e_sp_gas:.8f}" if e_sp_gas is not None else "—"
            aq_str  = f"{e_sp_aq:.8f}"  if e_sp_aq is not None  else "—"
            print(f"    {conf_id}: E_aq[{aq_source}]={aq_str}  E_gas={gas_str}  "
                  f"thermo={'✓' if tc['ok'] else '✗'}  n_imag={tc['n_imag']}")

    df = pd.DataFrame(rows)

    # Pesos de Boltzmann conformacional a cada T
    G_boltz: dict[str, Optional[float]] = {}
    for T, lbl in zip(TEMPS_K, TEMP_LABELS):
        col = f"G_total_{lbl}"
        valid = df[col].notna() & df["thermo_ok"]
        if valid.any():
            G_vals = df.loc[valid, col].values.astype(float)
            G_min  = G_vals.min()
            dG     = G_vals - G_min                           # en Hartree
            dG_kcal= dG * HARTREE_TO_KCAL
            weights = np.exp(-dG_kcal / (R_KCAL * T))
            weights /= weights.sum()
            df.loc[valid, f"boltz_{lbl}"] = weights
            df.loc[valid, f"dG_{lbl}_kcal"] = dG_kcal

            # G_Boltzmann del ensemble (free energy efectivo):
            # -RT·ln(Σ exp(-G_i/RT)) = G_min - RT·ln(Σ exp(-ΔG_i/RT))
            ln_Z   = np.log(np.sum(np.exp(-dG_kcal / (R_KCAL * T))))
            G_boltz[lbl] = G_min - (R_KCAL * T * ln_Z) / HARTREE_TO_KCAL
        else:
            G_boltz[lbl] = None
            if not df["thermo_ok"].any():
                warns.append(
                    f"[{mol_name}] No valid thermal corrections. "
                    f"Conformational Boltzmann distribution not available."
                )

    # Confórmero representativo (mínimo G_total a 298.15 K)
    col_298 = "G_total_298K"
    best_conf = None
    if col_298 in df.columns:
        valid = df[col_298].notna() & df["thermo_ok"]
        if valid.any():
            best_conf = df.loc[df.loc[valid, col_298].idxmin(), "conformer"]

    # Advertencia sobre frecuencias imaginarias
    n_imag_any = (df["n_imag"] > 0).sum()
    if n_imag_any > 0:
        warns.append(
            f"[{mol_name}] {n_imag_any} conformer(s) with imaginary frequencies. "
            f"Check geometries before using them in the pKa cycle."
        )

    # Warning if gas-phase SP is missing
    n_missing_gas = df["E_sp_gas"].isna().sum()
    if n_missing_gas > 0:
        warns.append(
            f"[{mol_name}] {n_missing_gas}/{len(df)} conformers without gas-phase SP "
            f"(expected at runs/{mol_name}/{SUBDIR_SP_GAS}/conf_NNN/sp_gas.out). "
            f"The pKa cycle requires the gas-phase SP of the representative conformer."
        )

    # Warning if SMD SP is missing (only relevant when use_smd=True — in
    # that case E_sp_aq itself comes from SUBDIR_SP_SMD, so a missing file
    # here means the conformer has NO aqueous energy at all, not just a
    # secondary/optional one like E_sp_gas).
    if use_smd and n_missing_smd > 0:
        warns.append(
            f"[{mol_name}] {n_missing_smd}/{len(df)} conformers without SMD SP "
            f"(expected at runs/{mol_name}/{SUBDIR_SP_SMD}/conf_NNN/sp_smd.out). "
            f"Run the pipeline with --spsmd for this molecule before using "
            f"--spsmd here. Those conformers will have E_sp_aq=NaN and be "
            f"excluded from the Boltzmann ensemble and from best_conf."
        )

    # ── Carga / multiplicidad / N° átomos a nivel de estado ────────────────
    # Todos los confórmeros de un mismo estado deben compartir carga, mult.
    # y N° de átomos (son el mismo sistema químico, solo geometrías
    # distintas). Si no es así, algo está mal en runs/<mol_name>/ (p. ej.
    # confórmeros de dos estados de protonación mezclados en la misma
    # carpeta) y NO se puede confiar en ningún valor único para la
    # validación cruzada con estados vecinos (_validate_state_charges).
    charges_seen  = sorted(set(df["charge"].dropna().astype(int)))
    mults_seen    = sorted(set(df["mult"].dropna().astype(int)))
    natoms_seen   = sorted(set(df["n_atoms"].dropna().astype(int)))

    state_charge  = charges_seen[0] if len(charges_seen) == 1 else None
    state_mult    = mults_seen[0]   if len(mults_seen)   == 1 else None
    state_n_atoms = natoms_seen[0]  if len(natoms_seen)  == 1 else None

    if len(charges_seen) > 1:
        warns.append(
            f"[{mol_name}] Inconsistent charge across conformers: "
            f"{charges_seen}. Cannot validate against neighboring states; "
            f"check runs/{mol_name}/{SUBDIR_SP}/conf_*/sp.out."
        )
    if len(natoms_seen) > 1:
        warns.append(
            f"[{mol_name}] Inconsistent atom count across conformers: "
            f"{natoms_seen}. Cannot validate against neighboring states; "
            f"check runs/{mol_name}/{SUBDIR_SP}/conf_*/sp.out."
        )
    if not charges_seen or not natoms_seen:
        warns.append(
            f"[{mol_name}] Could not determine charge and/or atom count "
            f"from any .out file. Cross-validation of charge/atom count "
            f"against neighboring states will be disabled for this "
            f"state (only an informational warning will be issued, not "
            f"an error)."
        )

    return {
        "mol_name":     mol_name,
        "conformers":   rows,
        "df":           df,
        "best_conf":    best_conf,
        "G_boltz":      G_boltz,      # Hartree, ensemble efectivo por T
        "charge":       state_charge,
        "mult":         state_mult,
        "n_atoms":      state_n_atoms,
        "ok":           len(rows) > 0,
        "warnings":     warns,
        "gcorr_source": gcorr_source,
        "aq_source":    aq_source,    # "CPCM" o "SMD" — de dónde vino E_sp_aq
        "use_smd":      use_smd,
    }


# ─────────────────────────────────────────────
#  ENERGÍA LIBRE REPRESENTATIVA DE UN ESTADO (método directo)
# ─────────────────────────────────────────────

def _state_G_total_direct(data: dict, T: float = T_STD,
                           gcorr_source: str = "harmonic",
                           use_ensemble: bool = False) -> dict:
    """
    G_total_CPCM del estado, para el método directo. Funciona con dos tipos
    de entrada:

      1. Un dict "normal" salido de collect_mol_data(): por default
         (use_ensemble=False) se relee freq.out del confórmero representativo
         (best_conf) y se recalcula G_corr a T (E_SP_aq(best_conf) + G_corr).

         Si use_ensemble=True, en vez de eso se usa data["G_boltz"][T_label]
         — el G_Boltzmann del ENSEMBLE conformacional completo, ya calculado
         por collect_mol_data() (línea ~1005-1008) — evitando descartar la
         población de confórmeros que no son el mínimo. Si no hay una entrada
         válida de G_boltz para la T pedida, cae automáticamente al
         comportamiento de best_conf y lo deja registrado en 'warn' (no falla
         silenciosamente).

      2. Un dict "macroestado" (data["is_macrostate"] = True), construido por
         compute_tautomer_equilibrium(): ya trae G_total_direct_Eh precalculado
         (combinación de Boltzmann entre neutro y zwitterión) y se usa tal cual,
         sin volver a tocar archivos. use_ensemble no aplica aquí (el macro-
         estado ya es, por construcción, una combinación de Boltzmann).

    Retorna dict: G_total_Eh, conf, theory_sp, theory_freq, ok, warn (str|None).
    """
    if data.get("is_macrostate"):
        return {
            "G_total_Eh":  data["G_total_direct_Eh"],
            "conf":        data.get("conf_label", "macrostate"),
            "theory_sp":   data.get("theory_label"),
            "theory_freq": data.get("theory_label"),
            "ok":          True,
            "warn":        None,
        }

    mol  = data["mol_name"]
    conf = data.get("best_conf")
    if conf is None:
        return {"ok": False, "warn": f"No representative conformer for {mol}."}

    use_smd = data.get("use_smd", False)
    freq_file = PC_WORK_ROOT / mol / SUBDIR_FREQ / conf / "freq.out"
    if use_smd:
        theory_sp = parse_theory_level(PC_WORK_ROOT / mol / SUBDIR_SP_SMD / conf / "sp_smd.out")
    else:
        theory_sp = parse_theory_level(PC_WORK_ROOT / mol / SUBDIR_SP / conf / "sp.out")
    theory_freq = parse_theory_level(freq_file)

    warn_fallback = None
    if use_ensemble:
        T_label = next(
            (lbl for T_i, lbl in zip(TEMPS_K, TEMP_LABELS) if np.isclose(T_i, T, atol=1e-3)),
            None,
        )
        G_ens = data.get("G_boltz", {}).get(T_label) if T_label is not None else None
        if G_ens is not None:
            return {
                "G_total_Eh":  G_ens,
                "conf":        f"ensemble Boltzmann ({conf}=best_conf)",
                "theory_sp":   theory_sp,
                "theory_freq": theory_freq,
                "ok":          True,
                "warn":        None,
            }
        warn_fallback = (
            f"use_ensemble=True but no valid G_boltz at T={T:.2f} K "
            f"for {mol} (T_label={T_label!r}); using the representative "
            f"conformer (best_conf={conf}) instead."
        )

    df = data["df"]
    row = df[df["conformer"] == conf].iloc[0]
    E_aq = row["E_sp_aq"]
    if E_aq is None or pd.isna(E_aq):
        return {"ok": False, "warn": f"No aqueous SP for {mol}/{conf}."}

    tc = get_thermochem(freq_file, T=T, source=gcorr_source,
                        gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)
    G_corr = tc["G_corr"] if tc["ok"] else None
    warn = None
    if not tc["ok"]:
        extra = f" ({tc['warning']})" if tc.get("warning") else ""
        warn = f"G_corr ({gcorr_source}) not available for {mol}/{conf}. Using G≈E_SP instead (no ZPE/entropy).{extra}"

    if warn_fallback:
        warn = f"{warn_fallback} {warn}" if warn else warn_fallback

    G_total = E_aq + (G_corr if G_corr is not None else 0.0)

    return {
        "G_total_Eh":  G_total,
        "conf":        conf,
        "theory_sp":   theory_sp,
        "theory_freq": theory_freq,
        "ok":          True,
        "warn":        warn,
    }


# ─────────────────────────────────────────────
#  EQUILIBRIO DE TAUTÓMEROS  (neutro ⇌ zwitterión, sin H⁺)
# ─────────────────────────────────────────────

def compute_tautomer_equilibrium(
    data_R: dict,
    data_D: dict,
    T: float = T_STD,
    verbose: bool = False,
    gcorr_source: str = "harmonic",
) -> dict:
    """
    Compara el neutro (R) y el zwitterión (D) como TAUTÓMEROS: mismo q=0,
    mismo número de átomos (solo se desplaza un H interno, p. ej. de COOH a
    N). A diferencia del ciclo de pKa, esto NO involucra H⁺ de bulk — es una
    comparación directa de G_total_CPCM entre confórmeros representativos.

      ΔG_taut = G_total(D) − G_total(R)     [kcal/mol]
      K_taut  = [D]/[R] = exp(−ΔG_taut / RT)
      f_R     = 1 / (1 + K_taut)
      f_D     = K_taut / (1 + K_taut)

    También construye el "macroestado neutro" combinado (R+D) vía Boltzmann,
    que se usa como base/ácido en compute_pka() para los pasos con el catión
    y con el anión:

      G_macro = G_min − RT·ln( Σ exp(−(G_i − G_min)/RT) )   i ∈ {R, D}

    Nota: f_R y f_D son independientes del pH (es un equilibrio intramolecular),
    así que se calculan una sola vez aquí y se reusan en toda la distribución.
    """
    mol_R = data_R["mol_name"]
    mol_D = data_D["mol_name"]
    warns: list[str] = []
    result: dict = {
        "mol_R": mol_R, "mol_D": mol_D, "T": T, "ok": False, "warnings": warns,
        "conf_R": None, "conf_D": None, "G_R_Eh": None, "G_D_Eh": None,
        "dG_taut_kcal": None, "K_taut": None, "f_R": None, "f_D": None,
        "G_macro_Eh": None, "theory_R": None, "theory_D": None,
    }

    def _warn(msg: str):
        warns.append(msg)
        warnings.warn(f"[tautómero {mol_R}⇌{mol_D}] {msg}")

    # ── Validación de carga / N° de átomos entre R y D ──────────────────────
    # R y D son tautómeros: mismo q, mismo N° de átomos (solo se desplaza un
    # H interno). Si difieren, no son tautómeros del mismo par ácido-base y
    # el "equilibrio" que se calcularía no tendría sentido físico.
    charge_msgs = _validate_state_charges(
        data_R, data_D, dq_expected=0, dn_expected=0,
        label_a=mol_R, label_b=mol_D,
    )
    for m in charge_msgs:
        _warn(m)
    if any(m.startswith("❌") for m in charge_msgs):
        _warn(
            f"Aborting tautomer equilibrium ({mol_R}⇌{mol_D}) due to "
            f"charge/atom-count inconsistency."
        )
        return result

    if data_R.get("use_smd", False) != data_D.get("use_smd", False):
        _warn(
            f"INCONSISTENT SOLVATION MODEL between tautomers: "
            f"{mol_R} uses aqueous SP from {data_R.get('aq_source', 'CPCM')}, "
            f"{mol_D} uses {data_D.get('aq_source', 'CPCM')}. Both tautomers "
            f"must use the same solvent model (both --spsmd or neither)."
        )

    gR = _state_G_total_direct(data_R, T=T, gcorr_source=gcorr_source)
    gD = _state_G_total_direct(data_D, T=T, gcorr_source=gcorr_source)

    if not gR.get("ok"):
        _warn(gR.get("warn") or f"Could not obtain G_total for {mol_R}.")
        return result
    if not gD.get("ok"):
        _warn(gD.get("warn") or f"Could not obtain G_total for {mol_D}.")
        return result
    if gR.get("warn"):
        warns.append(gR["warn"])
    if gD.get("warn"):
        warns.append(gD["warn"])

    if gR["theory_sp"] and gD["theory_sp"] and gR["theory_sp"] != gD["theory_sp"]:
        _warn(
            f"INCONSISTENT THEORY LEVEL in the aqueous SP between "
            f"{mol_R} ('{gR['theory_sp']}') and {mol_D} ('{gD['theory_sp']}'). "
            f"The tautomer comparison is not valid if basis set/functional/solvent model differ."
        )

    G_R = gR["G_total_Eh"]
    G_D = gD["G_total_Eh"]
    result["conf_R"], result["conf_D"] = gR["conf"], gD["conf"]
    result["G_R_Eh"], result["G_D_Eh"] = G_R, G_D
    result["theory_R"], result["theory_D"] = gR["theory_sp"], gD["theory_sp"]

    dG_taut_kcal = (G_D - G_R) * HARTREE_TO_KCAL
    K_taut = np.exp(-dG_taut_kcal / (R_KCAL * T))
    f_R = 1.0 / (1.0 + K_taut)
    f_D = K_taut / (1.0 + K_taut)

    # Macroestado combinado (Boltzmann de dos "confórmeros" R y D)
    G_min = min(G_R, G_D)
    dG_R_kcal = (G_R - G_min) * HARTREE_TO_KCAL
    dG_D_kcal = (G_D - G_min) * HARTREE_TO_KCAL
    lnZ = np.log(np.exp(-dG_R_kcal / (R_KCAL * T)) + np.exp(-dG_D_kcal / (R_KCAL * T)))
    G_macro = G_min - (R_KCAL * T * lnZ) / HARTREE_TO_KCAL

    result.update({
        "dG_taut_kcal": dG_taut_kcal, "K_taut": K_taut,
        "f_R": f_R, "f_D": f_D, "G_macro_Eh": G_macro, "ok": True,
        "aq_source": data_R.get("aq_source", "CPCM"),
    })

    if verbose:
        _print_tautomer_cycle(result)

    return result


def _print_tautomer_cycle(r: dict):
    """Imprime el equilibrio de tautómeros neutro⇌zwitterión."""
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  TAUTOMER EQUILIBRIUM (no H⁺) — {r['mol_R']} ⇌ {r['mol_D']}")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. R (neutral)     : {r['conf_R']}")
    print(f"  Conf. D (zwitterion)  : {r['conf_D']}")
    print(f"  SP level (R) : {r['theory_R'] or '—'}")
    print(f"  SP level (D) : {r['theory_D'] or '—'}")
    print()
    aq_src_taut = r.get("aq_source", "CPCM")
    print(f"  {'G_total_' + aq_src_taut + '(R)  [Eh]':<40} {r['G_R_Eh']:>14.8f}")
    print(f"  {'G_total_' + aq_src_taut + '(D)  [Eh]':<40} {r['G_D_Eh']:>14.8f}")
    print(f"  {'─'*40}  {'─'*14}")
    print(f"  {'ΔG_taut = G(D) − G(R)  [kcal/mol]':<40} {r['dG_taut_kcal']:>+14.3f}")
    print(f"  {'K_taut = [D]/[R]':<40} {r['K_taut']:>14.4g}")
    print()
    print(f"  f(R, neutral)     = {r['f_R']*100:6.2f} %")
    print(f"  f(D, zwitterion)  = {r['f_D']*100:6.2f} %")
    dominant = "R (neutral)" if r["f_R"] >= r["f_D"] else "D (zwitterion)"
    print(f"  Dominant tautomer: {dominant}")
    print(f"\n  Neutral macrostate G (Boltzmann R+D) = {r['G_macro_Eh']:.8f} Eh")
    print(f"  (this combined value is what enters the pKa cycle as the")
    print(f"   'neutral/zwitterion' state versus the cation and the anion)")
    if r["warnings"]:
        print(f"\n  Warnings:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


# ─────────────────────────────────────────────
#  CICLO TERMODINÁMICO DE pKa
# ─────────────────────────────────────────────

def compute_pka(
    data_acid: dict,
    data_base: dict,
    T: float = T_STD,
    verbose: bool = False,
    method: str = "auto",
    gcorr_source: str = "harmonic",
) -> dict:
    """
    Calcula el pKa de AH → A⁻ + H⁺.

    Dos métodos disponibles:

    1. "born"   — Ciclo de Born completo (Ho & Coote 2010). Requiere SP en
                  fase gas (04b_sp_gas) además del SP acuoso (CPCM):
                    ΔG_solv(X) = E_SP(CPCM) − E_SP(gas) + δG_conc
                    ΔG_aq = ΔG_gas + ΔG_solv(A⁻) + ΔG_solv(H⁺) − ΔG_solv(AH)
                  Usa G(H⁺,aq) implícito de Tissandier 1998 + Kelly 2006
                  (−272.18 kcal/mol), calibrado para cerrar el ciclo gas+solv
                  explícito.

    2. "direct" — Método directo solo-CPCM. No requiere SP fase gas; usa
                  directamente G_total ya calculado con CPCM (E_SP_CPCM + G_corr)
                  y la energía libre del protón de Pliego & Riveros 2002
                  (G_HPLUS_AQ_DIRECT = −264.61 kcal/mol):
                    ΔG_aq = G_total_CPCM(A⁻) − G_total_CPCM(AH) + G(H⁺, aq)
                  IMPORTANTE: esta constante NO es la misma que la de "born".
                  Las energías E_SP_CPCM ya incluyen el efecto de solvatación
                  de la molécula vía el modelo de continuo; usar la constante
                  de Tissandier/Kelly aquí sobreestima la estabilización del
                  protón en ~8-10 kcal/mol porque esa constante está calibrada
                  para emparejarse con energías en fase gas pura, no con
                  energías ya en solución implícita.

    method="auto" (default): usa "born" si hay SP fase gas disponibles para
    ambos confórmeros representativos; si no, cae automáticamente a "direct"
    y lo deja registrado en result["method"] / result["warnings"].

    Usa el confórmero representativo (mínimo G_total a 298 K) de cada estado.
    Retorna dict con todas las contribuciones al ciclo y el pKa final.
    """
    mol_acid = data_acid["mol_name"]
    mol_base = data_base["mol_name"]
    warns: list[str] = []
    result: dict = {
        "mol_acid": mol_acid, "mol_base": mol_base,
        "T": T, "pKa": None, "ok": False, "warnings": warns,
        "method": None,
        # contribuciones (ciclo de Born)
        "dG_gas_kcal": None, "dG_solv_acid_kcal": None,
        "dG_solv_base_kcal": None, "dG_solv_hplus_kcal": DG_SOLV_HPLUS,
        "dG_aq_kcal": None,
        # contribuciones (método directo)
        "G_total_acid_Eh": None, "G_total_base_Eh": None,
        "G_hplus_aq_kcal": None,
        # datos intermedios
        "conf_acid": None, "conf_base": None,
        "G_gas_acid_Eh": None, "G_gas_base_Eh": None,
        "E_sp_aq_acid_Eh": None, "E_sp_aq_base_Eh": None,
        "E_sp_gas_acid_Eh": None, "E_sp_gas_base_Eh": None,
        "G_corr_acid_Eh": None, "G_corr_base_Eh": None,
        "theory_sp_acid": None, "theory_sp_base": None,
        "theory_freq_acid": None, "theory_freq_base": None,
    }

    def _warn(msg: str):
        warns.append(msg)
        warnings.warn(f"[pKa {mol_acid}→{mol_base}] {msg}")

    if method not in ("auto", "born", "direct"):
        raise ValueError(f"method debe ser 'auto', 'born' o 'direct', no {method!r}")

    # ── Validación de carga / N° de átomos entre ácido y base ──────────────
    # Un pKa modela AH → A⁻ + H⁺: la base SIEMPRE debe tener Δq=-1 y
    # ΔN_átomos=-1 respecto al ácido. Este chequeo corre ANTES de cualquier
    # cálculo, para atrapar de inmediato un mapeo de carpetas desalineado
    # (p. ej. --anion apuntando al estado equivocado, o un best_conf que no
    # corresponde a la geometría esperada) en vez de dejar que el script
    # calcule un pKa numéricamente válido pero físicamente sin sentido.
    charge_msgs = _validate_state_charges(
        data_acid, data_base, dq_expected=-1, dn_expected=-1,
        label_a=mol_acid, label_b=mol_base,
    )
    for m in charge_msgs:
        _warn(m)
    if any(m.startswith("❌") for m in charge_msgs):
        _warn(
            f"Aborting pKa calculation ({mol_acid}→{mol_base}) due to "
            f"charge/atom-count inconsistency. Check the state mapping "
            f"before continuing."
        )
        return result

    # ── Rama simplificada: uno de los dos lados es un macroestado combinado ─
    # (neutro+zwitterión, construido por compute_tautomer_equilibrium).
    #
    # method='born' SÍ se puede resolver aquí, y de forma EXACTA, sin ningún
    # SP de fase gas para la mezcla. La razón es la misma identidad algebraica
    # ya documentada arriba de compute_pka_montecarlo() para el caso de una
    # sola especie: con las definiciones que usa este script,
    #
    #   G_gas(X)    = E_gas(X) + G_corr(X)
    #   ΔG_solv(X)  = E_aq(X) − E_gas(X) + δG_conc
    #
    # TODOS los términos en E_gas se cancelan al construir ΔG_aq del ciclo de
    # Born, y lo que queda es exactamente la resta de energías CPCM (idéntica
    # a la del método directo) más la constante G(H⁺,aq) correspondiente:
    #
    #   ΔG_aq(born)   = G_total_CPCM(base) − G_total_CPCM(acid) + G_HPLUS_AQ_BORN
    #   ΔG_aq(direct) = G_total_CPCM(base) − G_total_CPCM(acid) + G_HPLUS_AQ_DIRECT
    #
    # Esa cancelación se extiende sin cambios al macroestado combinado: si se
    # define el "G_gas" de la mezcla con la MISMA entropía de mezcla (mismos
    # pesos f_R/f_D) que ya usa G_macro_Eh en compute_tautomer_equilibrium, el
    # término de entropía de mezcla aparece igual en G_macro_aq y en
    # G_macro_gas y se cancela exactamente igual que G_corr en el caso de una
    # sola especie. No hace falta calcular ni siquiera construir ese
    # G_macro_gas explícitamente — el resultado final es el mismo de arriba,
    # con G_total_CPCM(macro) = G_macro_Eh ya calculado por
    # compute_tautomer_equilibrium(). Por eso tampoco hacía falta el SP de
    # fase gas para el lado que NO es macroestado: esos términos también se
    # cancelan ahí.
    if data_acid.get("is_macrostate") or data_base.get("is_macrostate"):
        # use_ensemble=True: para el lado que NO es macroestado (A o B), usa
        # el G_Boltzmann del ensemble conformacional completo en vez de solo
        # el best_conf — ver nota #3 del review. Para el lado que SÍ es
        # macroestado, _state_G_total_direct ignora use_ensemble (ya entra
        # por la rama is_macrostate de arriba).
        use_smd_acid = data_acid.get("use_smd", False)
        use_smd_base = data_base.get("use_smd", False)
        if use_smd_acid != use_smd_base:
            _warn(
                f"INCONSISTENT SOLVATION MODEL between states: "
                f"{mol_acid} uses aqueous SP from "
                f"{data_acid.get('aq_source', 'CPCM')}, {mol_base} uses "
                f"{data_base.get('aq_source', 'CPCM')}. Both --spsmd or "
                f"neither — the energy difference is not physically "
                f"meaningful otherwise."
            )

        g_acid = _state_G_total_direct(data_acid, T=T, gcorr_source=gcorr_source,
                                       use_ensemble=True)
        g_base = _state_G_total_direct(data_base, T=T, gcorr_source=gcorr_source,
                                       use_ensemble=True)

        result["conf_acid"] = g_acid.get("conf")
        result["conf_base"] = g_base.get("conf")
        result["theory_sp_acid"] = g_acid.get("theory_sp")
        result["theory_sp_base"] = g_base.get("theory_sp")
        result["aq_source_acid"] = data_acid.get("aq_source", "CPCM")
        result["aq_source_base"] = data_base.get("aq_source", "CPCM")

        if not g_acid.get("ok"):
            _warn(g_acid.get("warn") or f"Could not obtain G_total for {mol_acid}.")
            return result
        if not g_base.get("ok"):
            _warn(g_base.get("warn") or f"Could not obtain G_total for {mol_base}.")
            return result
        if g_acid.get("warn"):
            warns.append(g_acid["warn"])
        if g_base.get("warn"):
            warns.append(g_base["warn"])

        if (g_acid["theory_sp"] and g_base["theory_sp"]
                and g_acid["theory_sp"] != g_base["theory_sp"]):
            _warn(
                f"INCONSISTENT THEORY LEVEL between {mol_acid} "
                f"('{g_acid['theory_sp']}') and {mol_base} ('{g_base['theory_sp']}')."
            )

        G_total_acid = g_acid["G_total_Eh"]
        G_total_base = g_base["G_total_Eh"]
        result["G_total_acid_Eh"] = G_total_acid
        result["G_total_base_Eh"] = G_total_base

        if method == "born":
            use_born_macro = True
        elif method == "direct":
            use_born_macro = False
        else:  # "auto"
            # Exacto sin SP de fase gas, por la derivación algebraica de
            # arriba: el ciclo de Born se reduce a la misma resta CPCM que
            # 'direct', solo con la constante G_HPLUS_AQ_BORN en vez de
            # G_HPLUS_AQ_DIRECT. No hay razón para que 'auto' degrade a
            # 'direct' aquí — a diferencia de la rama simple, donde 'auto'
            # SÍ depende de si existen los SP de fase gas.
            use_born_macro = True

        if use_born_macro:
            result["G_hplus_aq_kcal"] = G_HPLUS_AQ_BORN
            result["method"] = "born"
            if method == "auto":
                _warn(
                    f"method='auto' resolved to 'born' for the combined "
                    f"neutral+zwitterion macrostate "
                    f"({mol_acid if data_acid.get('is_macrostate') else mol_base}): "
                    f"no gas-phase SP is needed here (see derivation note "
                    f"above), so 'born' is used exactly, same as it would be "
                    f"for a single species with gas-phase SPs available."
                )
            else:
                _warn(
                    f"method='born' with a combined neutral+zwitterion macrostate "
                    f"({mol_acid if data_acid.get('is_macrostate') else mol_base}): "
                    f"no gas-phase SP is used (or needed) for the mixture. ΔG_aq "
                    f"is obtained via the algebraic reduction of the Born cycle "
                    f"to G_total_CPCM(base) − G_total_CPCM(acid) + G_HPLUS_AQ_BORN "
                    f"— the same reduction that already applies to a single "
                    f"species in this script, since the E_gas terms cancel "
                    f"identically (see the derivation note above "
                    f"compute_pka_montecarlo())."
                )
        else:
            result["G_hplus_aq_kcal"] = G_HPLUS_AQ_DIRECT
            result["method"] = "direct"
            if method == "direct":
                _warn(
                    f"method='direct' explicitly requested for the combined "
                    f"neutral+zwitterion macrostate "
                    f"({mol_acid if data_acid.get('is_macrostate') else mol_base})."
                )

        dG_aq_kcal = (
            (G_total_base - G_total_acid) * HARTREE_TO_KCAL
            + result["G_hplus_aq_kcal"]
        )
        result["dG_aq_kcal"] = dG_aq_kcal
        result["pKa"] = dG_aq_kcal / (R_KCAL * T * LN10)
        result["ok"] = True

        if verbose:
            if result["method"] == "born":
                _print_born_cycle_macrostate(result)
            else:
                _print_direct_cycle(result)

        return result

    # ── Obtener confórmero representativo de cada estado ──────────────────
    conf_acid = data_acid.get("best_conf")
    conf_base = data_base.get("best_conf")
    result["conf_acid"] = conf_acid
    result["conf_base"] = conf_base

    if conf_acid is None:
        _warn(f"No representative conformer for {mol_acid}.")
        return result
    if conf_base is None:
        _warn(f"No representative conformer for {mol_base}.")
        return result

    df_acid = data_acid["df"]
    df_base = data_base["df"]

    row_acid = df_acid[df_acid["conformer"] == conf_acid].iloc[0]
    row_base = df_base[df_base["conformer"] == conf_base].iloc[0]

    # ── ¿Alguno de los dos estados usa SMD en vez de CPCM? ─────────────────
    use_smd_acid = data_acid.get("use_smd", False)
    use_smd_base = data_base.get("use_smd", False)
    aq_source_acid = data_acid.get("aq_source", "CPCM")
    aq_source_base = data_base.get("aq_source", "CPCM")

    if use_smd_acid != use_smd_base:
        _warn(
            f"INCONSISTENT SOLVATION MODEL between states: "
            f"{mol_acid} uses aqueous SP from {aq_source_acid}, "
            f"{mol_base} uses {aq_source_base}. Both states of a pKa pair "
            f"must use the same solvent model (both --spsmd or neither) — "
            f"the energy difference is not physically meaningful otherwise."
        )

    # ── Verificación de consistencia de nivel de teoría ────────────────────
    # Compara las líneas de keywords ORCA del SP acuoso y del freq.out entre
    # acid y base. Un desajuste aquí (distinta base, funcional, o modelo de
    # solvente/grid) invalida la resta de energías de ~1000+ Eh para sacar
    # una diferencia de ~0.4 Eh, independientemente de qué convención de
    # G(H⁺) se use después.
    sp_acid_file = (
        PC_WORK_ROOT / mol_acid / SUBDIR_SP_SMD / conf_acid / "sp_smd.out"
        if use_smd_acid else
        PC_WORK_ROOT / mol_acid / SUBDIR_SP / conf_acid / "sp.out"
    )
    sp_base_file = (
        PC_WORK_ROOT / mol_base / SUBDIR_SP_SMD / conf_base / "sp_smd.out"
        if use_smd_base else
        PC_WORK_ROOT / mol_base / SUBDIR_SP / conf_base / "sp.out"
    )
    freq_acid_lvl = parse_theory_level(PC_WORK_ROOT / mol_acid / SUBDIR_FREQ / conf_acid / "freq.out")
    freq_base_lvl = parse_theory_level(PC_WORK_ROOT / mol_base / SUBDIR_FREQ / conf_base / "freq.out")
    sp_acid_lvl   = parse_theory_level(sp_acid_file)
    sp_base_lvl   = parse_theory_level(sp_base_file)

    result["theory_sp_acid"]   = sp_acid_lvl
    result["theory_sp_base"]   = sp_base_lvl
    result["theory_freq_acid"] = freq_acid_lvl
    result["theory_freq_base"] = freq_base_lvl
    result["aq_source_acid"]   = aq_source_acid
    result["aq_source_base"]   = aq_source_base

    if sp_acid_lvl and sp_base_lvl and sp_acid_lvl != sp_base_lvl:
        _warn(
            f"INCONSISTENT THEORY LEVEL in the aqueous SP between "
            f"{mol_acid} ('{sp_acid_lvl}') and {mol_base} ('{sp_base_lvl}'). "
            f"The pKa cycle is NOT valid if the basis set/functional/solvent "
            f"model differ between states — check before reporting."
        )
    if freq_acid_lvl and freq_base_lvl and freq_acid_lvl != freq_base_lvl:
        _warn(
            f"INCONSISTENT THEORY LEVEL in NumFreq between "
            f"{mol_acid} ('{freq_acid_lvl}') and {mol_base} ('{freq_base_lvl}')."
        )

    # ── Energías SP acuosas (CPCM o SMD, según use_smd) ────────────────────
    E_aq_acid = row_acid["E_sp_aq"]
    E_aq_base = row_base["E_sp_aq"]
    result["E_sp_aq_acid_Eh"] = E_aq_acid
    result["E_sp_aq_base_Eh"] = E_aq_base

    if E_aq_acid is None or pd.isna(E_aq_acid):
        _warn(f"No aqueous SP ({aq_source_acid}) for {mol_acid}/{conf_acid}.")
        return result
    if E_aq_base is None or pd.isna(E_aq_base):
        _warn(f"No aqueous SP ({aq_source_base}) for {mol_base}/{conf_base}.")
        return result

    # ── Energías SP fase gas (pueden faltar) ───────────────────────────────
    E_gas_acid = row_acid["E_sp_gas"]
    E_gas_base = row_base["E_sp_gas"]
    result["E_sp_gas_acid_Eh"] = E_gas_acid
    result["E_sp_gas_base_Eh"] = E_gas_base

    have_gas = not (
        E_gas_acid is None or pd.isna(E_gas_acid) or
        E_gas_base is None or pd.isna(E_gas_base)
    )

    # ── Decidir método ──────────────────────────────────────────────────────
    if method == "born" and not have_gas:
        _warn(
            f"method='born' requested but gas-phase SPs are missing for "
            f"{mol_acid}/{conf_acid} and/or {mol_base}/{conf_base}. "
            f"Create runs/{mol_acid}/{SUBDIR_SP_GAS}/{conf_acid}/sp_gas.inp "
            f"and the equivalent for {mol_base} (same functional/basis, no CPCM)."
        )
        return result

    use_born = have_gas and method in ("auto", "born")

    if use_born and (use_smd_acid or use_smd_base):
        _warn(
            f"method='born' with --spsmd: the Born-cycle constants "
            f"(G_HPLUS_GAS, DG_SOLV_HPLUS, DG_CONC_CORR) were calibrated "
            f"against CPCM(water) aqueous energies, not SMD(water). "
            f"ΔG_solv(X) = E_SP(SMD) − E_SP(gas) + δG_conc mixes solvent "
            f"models against a cycle closed for CPCM — treat this pKa as "
            f"uncalibrated/exploratory. Consider method='direct' with "
            f"--spsmd instead, which only compares SMD energies to each "
            f"other and does not depend on the Born-cycle gas-phase constants."
        )

    if method == "auto" and not have_gas:
        _warn(
            f"No gas-phase SP for {mol_acid}/{conf_acid} and/or {mol_base}/{conf_base}. "
            f"Using the direct method ({aq_source_acid} only) instead: "
            f"ΔG_aq = G_total_{aq_source_base}(base) − G_total_{aq_source_acid}(acid) + G(H⁺,aq). "
            f"For the full Born cycle, create the gas-phase SPs in "
            f"{SUBDIR_SP_GAS}/conf_NNN/sp_gas.out."
        )

    # ── Correcciones térmicas G_corr a temperatura T (necesarias en ambos métodos) ─
    freq_acid = PC_WORK_ROOT / mol_acid / SUBDIR_FREQ / conf_acid / "freq.out"
    freq_base = PC_WORK_ROOT / mol_base / SUBDIR_FREQ / conf_base / "freq.out"
    tc_acid = get_thermochem(freq_acid, T=T, source=gcorr_source,
                             gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)
    tc_base = get_thermochem(freq_base, T=T, source=gcorr_source,
                             gv_cutoff=GV_CUTOFF, gv_qs=GV_QS, gv_qh=GV_QH)

    G_corr_acid: Optional[float] = None
    G_corr_base: Optional[float] = None

    if tc_acid["ok"]:
        G_corr_acid = tc_acid["G_corr"]
    else:
        extra = f" ({tc_acid['warning']})" if tc_acid.get("warning") else ""
        _warn(
            f"G_corr ({gcorr_source}) not available for {mol_acid}/{conf_acid}. "
            f"Using G ≈ E_SP instead (no ZPE/entropy; pKa has higher uncertainty).{extra}"
        )
    if tc_base["ok"]:
        G_corr_base = tc_base["G_corr"]
    else:
        extra = f" ({tc_base['warning']})" if tc_base.get("warning") else ""
        _warn(
            f"G_corr ({gcorr_source}) not available for {mol_base}/{conf_base}. "
            f"Using G ≈ E_SP.{extra}"
        )

    result["G_corr_acid_Eh"] = G_corr_acid
    result["G_corr_base_Eh"] = G_corr_base

    if use_born:
        # ── G_total en fase gas ───────────────────────────────────────────
        # G_gas(X) = E_SP_gas(X) + G_corr(X)
        # Aproximación Ho & Coote: G_corr del cálculo CPCM es válido en gas.
        G_gas_acid = E_gas_acid + (G_corr_acid if G_corr_acid is not None else 0.0)
        G_gas_base = E_gas_base + (G_corr_base if G_corr_base is not None else 0.0)
        result["G_gas_acid_Eh"] = G_gas_acid
        result["G_gas_base_Eh"] = G_gas_base

        # ── ΔG_gas ─────────────────────────────────────────────────────────
        # ΔG_gas = G(A⁻, gas) + G(H⁺, gas) − G(AH, gas)
        dG_gas_kcal = (G_gas_base - G_gas_acid) * HARTREE_TO_KCAL + G_HPLUS_GAS
        result["dG_gas_kcal"] = dG_gas_kcal

        # ── ΔG_solv por especie ───────────────────────────────────────────
        # ΔG_solv(X) = E_SP(aq_source) − E_SP(gas) + δG_conc
        # aq_source es CPCM salvo que use_smd=True para ese estado (ver
        # warning arriba: mezclar SMD con las constantes de Born calibradas
        # para CPCM deja el resultado como uncalibrado/exploratorio).
        dG_solv_acid = (E_aq_acid - E_gas_acid) * HARTREE_TO_KCAL + DG_CONC_CORR
        dG_solv_base = (E_aq_base - E_gas_base) * HARTREE_TO_KCAL + DG_CONC_CORR
        result["dG_solv_acid_kcal"] = dG_solv_acid
        result["dG_solv_base_kcal"] = dG_solv_base

        # ── ΔG_aq (ciclo de Born completo) ────────────────────────────────
        # ΔG_aq = ΔG_gas + ΔG_solv(A⁻) + ΔG_solv(H⁺) − ΔG_solv(AH)
        dG_aq_kcal = (
            dG_gas_kcal
            + dG_solv_base
            + DG_SOLV_HPLUS
            - dG_solv_acid
        )
        result["dG_aq_kcal"] = dG_aq_kcal
        result["method"] = "born"

    else:
        # ── Método directo: solo energías acuosas (CPCM o SMD) ────────────
        # G_total(X) = E_SP(aq_source) + G_corr(X)
        G_total_acid = E_aq_acid + (G_corr_acid if G_corr_acid is not None else 0.0)
        G_total_base = E_aq_base + (G_corr_base if G_corr_base is not None else 0.0)
        result["G_total_acid_Eh"] = G_total_acid
        result["G_total_base_Eh"] = G_total_base
        result["G_hplus_aq_kcal"] = G_HPLUS_AQ_DIRECT

        # ΔG_aq = G_total(A⁻) − G_total(AH) + G(H⁺, aq)   [aq_source: CPCM o SMD]
        dG_aq_kcal = (G_total_base - G_total_acid) * HARTREE_TO_KCAL + G_HPLUS_AQ_DIRECT
        result["dG_aq_kcal"] = dG_aq_kcal
        result["method"] = "direct"

    # ── pKa ──────────────────────────────────────────────────────────────
    pKa = result["dG_aq_kcal"] / (R_KCAL * T * LN10)
    result["pKa"] = pKa
    result["ok"]  = True

    if verbose:
        if result["method"] == "born":
            _print_born_cycle(result)
        else:
            _print_direct_cycle(result)

    return result


# ─────────────────────────────────────────────
#  MONTE CARLO — PROPAGACIÓN DE INCERTIDUMBRE EN pKa (FASE 1)
# ─────────────────────────────────────────────
#
# NOTA DE DERIVACIÓN (por qué el ruido se aplica solo a E_aq y G_corr, no a
# E_gas por separado, en el método "born"):
#
#   dG_gas       = (G_gas_base − G_gas_acid)·H2K + G_HPLUS_GAS
#   dG_solv_X    = (E_aq_X − E_gas_X)·H2K + DG_CONC_CORR
#   dG_aq(born)  = dG_gas + dG_solv_base + DG_SOLV_HPLUS − dG_solv_acid
#
#   Sustituyendo G_gas_X = E_gas_X + G_corr_X y expandiendo, TODOS los
#   términos en E_gas_acid/E_gas_base se cancelan algebraicamente (y también
#   DG_CONC_CORR, que aparece una vez con cada signo). Lo que queda es
#   idéntico en forma al método directo:
#
#   dG_aq(born) = (E_aq_base − E_aq_acid + G_corr_base − G_corr_acid)·H2K
#                 + (G_HPLUS_GAS + DG_SOLV_HPLUS)      # = G_HPLUS_AQ_BORN
#
#   dG_aq(direct) = (E_aq_base − E_aq_acid + G_corr_base − G_corr_acid)·H2K
#                   + G_HPLUS_AQ_DIRECT
#
#   O sea: dado cómo compute_pka() arma el ciclo de Born hoy, "born" y
#   "direct" tienen EXACTAMENTE la misma sensibilidad al ruido en E_aq/G_corr;
#   solo difieren en qué constante de G(H⁺,aq) se suma. Por eso el modelo de
#   ruido de abajo es el mismo para ambos métodos, y solo cambia el punto de
#   partida (base_result["dG_aq_kcal"]). Si algún día cambia la fórmula de
#   compute_pka() (p. ej. DG_CONC_CORR deja de cancelarse), esta nota y el
#   test test_born_and_direct_same_noise_sensitivity() dejan de ser válidos
#   y hay que revisar esta función.
#
# FUENTES DE INCERTIDUMBRE MODELADAS (ver plan, sección 1):
#   - sigma_elec_kcal  : energía electrónica SP (funcional/base), por estado,
#                        independiente entre ácido y base.
#   - sigma_gcorr_kcal : corrección térmica G_corr (ZPE+TS), por estado,
#                        independiente entre ácido y base.
#   - sigma_hplus_kcal : G(H⁺,aq) — constante compartida por TODA la cadena
#                        de pKa's. Por eso admite shared_hplus_draw: al
#                        correr varios pKa's de la misma cadena (pKa0, pKa1,
#                        pKa2) hay que pasarles el MISMO array de muestras
#                        para esta constante, o la incertidumbre de
#                        cantidades derivadas (Δ entre pKa's, fracciones de
#                        especie) queda sobreestimada.

def compute_pka_montecarlo(
    data_acid: dict,
    data_base: dict,
    T: float = T_STD,
    method: str = "auto",
    gcorr_source: str = "harmonic",
    n_samples: int = 10_000,
    sigma_elec_kcal: float = 0.5,
    sigma_gcorr_kcal: float = 0.4,
    sigma_hplus_kcal: float = 2.5,
    shared_hplus_draw: Optional[np.ndarray] = None,
    rng: Optional[np.random.Generator] = None,
) -> dict:
    """
    Propaga incertidumbre en el pKa de AH → A⁻ + H⁺ vía Monte Carlo,
    reusando compute_pka() como fuente única de verdad para el valor
    determinístico (mismo mapeo de estados, misma validación de carga/N°
    de átomos, mismo método born/direct/auto).

    No vuelve a leer archivos ni reimplementa el ciclo termodinámico: toma
    dG_aq_kcal ya calculado por compute_pka() como centro de la distribución
    y le suma perturbaciones gaussianas sobre los términos que compute_pka()
    ya expone en su dict de resultado (E_sp_aq_acid/base_Eh — si están
    disponibles — o, si no, sobre un sigma combinado, ver el caso de
    macroestados más abajo).

    Parameters
    ----------
    data_acid, data_base : dict
        Mismos dicts que recibiría compute_pka() (salida de
        collect_mol_data() o de compute_tautomer_equilibrium()).
    n_samples : int
        Número de muestras Monte Carlo.
    sigma_elec_kcal, sigma_gcorr_kcal : float
        Desviación estándar asumida para la energía electrónica SP y la
        corrección térmica, respectivamente, por estado (kcal/mol).
    sigma_hplus_kcal : float
        Desviación estándar asumida para G(H⁺,aq) — dispersión típica entre
        valores de literatura (Pliego-Riveros, Kelly-Cramer-Truhlar,
        Zhan-Dixon); ver plan, sección 1.
    shared_hplus_draw : np.ndarray, opcional
        Si se provee (shape (n_samples,)), se usa como la perturbación de
        G(H⁺,aq) en vez de samplear una nueva — para correlacionar el mismo
        error sistemático entre varios pKa's de una misma cadena (ver
        run_chain_montecarlo(), fase 2).
    rng : np.random.Generator, opcional
        Generador de números aleatorios; si no se provee, usa uno sin
        semilla fija (no reproducible entre corridas).

    Returns
    -------
    dict con:
        ok               : bool
        base_result      : dict — salida de compute_pka() (determinístico)
        method_used      : 'born' | 'direct' (resuelto, nunca 'auto')
        pka_samples      : np.ndarray, shape (n_samples,)
        pka_mean         : float
        pka_std          : float (ddof=1)
        pka_median       : float
        pka_ci95         : tuple (lo, hi) — percentiles 2.5/97.5
        hplus_draw       : np.ndarray — la perturbación de G(H⁺,aq) usada
                            (para reusar como shared_hplus_draw en otro pKa)
        warnings         : list[str]
    """
    base_result = compute_pka(
        data_acid, data_base, T=T, method=method,
        gcorr_source=gcorr_source, verbose=False,
    )

    out: dict = {
        "ok": False,
        "base_result": base_result,
        "method_used": base_result.get("method"),
        "n_samples": n_samples,
        "pka_samples": None,
        "pka_mean": None, "pka_std": None, "pka_median": None,
        "pka_ci95": None,
        "hplus_draw": None,
        "sigma_elec_kcal": sigma_elec_kcal,
        "sigma_gcorr_kcal": sigma_gcorr_kcal,
        "sigma_hplus_kcal": sigma_hplus_kcal,
        "warnings": list(base_result.get("warnings", [])),
    }

    if not base_result.get("ok"):
        out["warnings"].append(
            "compute_pka() determinístico falló (ok=False); no se puede "
            "construir una distribución Monte Carlo encima de un ciclo "
            "inválido. Revisa los warnings de base_result."
        )
        return out

    if n_samples < 2:
        out["warnings"].append(f"n_samples={n_samples} demasiado bajo; se requieren ≥2.")
        return out

    if shared_hplus_draw is not None and len(shared_hplus_draw) != n_samples:
        raise ValueError(
            f"shared_hplus_draw tiene {len(shared_hplus_draw)} muestras, "
            f"pero n_samples={n_samples}. Deben coincidir para poder sumarlo "
            f"término a término con el resto de las perturbaciones."
        )

    if rng is None:
        rng = np.random.default_rng()

    d_hplus = (
        shared_hplus_draw if shared_hplus_draw is not None
        else rng.normal(0.0, sigma_hplus_kcal, size=n_samples)
    )

    base_dG_aq = base_result["dG_aq_kcal"]

    # ── ¿Tenemos E_sp_aq por separado (caso normal) o solo G_total combinado
    #    (caso macroestado neutro+zwitterión, ver compute_tautomer_equilibrium)? ──
    have_separate_terms = (
        base_result.get("E_sp_aq_acid_Eh") is not None and
        base_result.get("E_sp_aq_base_Eh") is not None
    )

    if have_separate_terms:
        # Ruido independiente en energía electrónica y corrección térmica,
        # por estado (ácido y base son cálculos distintos → errores no
        # correlacionados entre sí). Ver nota de derivación arriba: esto
        # es válido tanto para method_used='born' como 'direct'.
        d_elec_acid  = rng.normal(0.0, sigma_elec_kcal,  size=n_samples)
        d_elec_base  = rng.normal(0.0, sigma_elec_kcal,  size=n_samples)
        d_gcorr_acid = rng.normal(0.0, sigma_gcorr_kcal, size=n_samples)
        d_gcorr_base = rng.normal(0.0, sigma_gcorr_kcal, size=n_samples)

        dG_aq_samples = (
            base_dG_aq
            + (d_elec_base - d_elec_acid)
            + (d_gcorr_base - d_gcorr_acid)
            + d_hplus
        )
    else:
        # Macroestado combinado (o, por alguna razón, sin desglose): no hay
        # forma de separar contribución electrónica de térmica dentro de
        # G_total_Eh ya mezclado con pesos de Boltzmann internos. Se usa un
        # sigma combinado en cuadratura como aproximación conservadora y se
        # deja registrado en warnings — esto es una simplificación de Fase 1,
        # no un tratamiento riguroso del error del ensemble.
        out["warnings"].append(
            "Sin desglose E_sp_aq/G_corr por separado (probablemente un "
            "macroestado neutro+zwitterión); se usa un sigma combinado "
            "sqrt(sigma_elec² + sigma_gcorr²) por estado en vez de separar "
            "las dos fuentes de error."
        )
        sigma_total = float(np.hypot(sigma_elec_kcal, sigma_gcorr_kcal))
        d_total_acid = rng.normal(0.0, sigma_total, size=n_samples)
        d_total_base = rng.normal(0.0, sigma_total, size=n_samples)

        dG_aq_samples = (
            base_dG_aq
            + (d_total_base - d_total_acid)
            + d_hplus
        )

    pka_samples = dG_aq_samples / (R_KCAL * T * LN10)
    ci_lo, ci_hi = np.percentile(pka_samples, [2.5, 97.5])

    out.update({
        "ok": True,
        "pka_samples": pka_samples,
        "pka_mean":   float(np.mean(pka_samples)),
        "pka_std":    float(np.std(pka_samples, ddof=1)),
        "pka_median": float(np.median(pka_samples)),
        "pka_ci95":   (float(ci_lo), float(ci_hi)),
        "hplus_draw": d_hplus,
    })
    return out


def _print_born_cycle(r: dict):
    """Imprime las contribuciones detalladas del ciclo de Born."""
    aq_acid_bc = r.get("aq_source_acid", "CPCM")
    aq_base_bc = r.get("aq_source_base", "CPCM")
    aq_label_bc = aq_acid_bc if aq_acid_bc == aq_base_bc else f"{aq_acid_bc}/{aq_base_bc}"
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  BORN CYCLE — {r['mol_acid']} → {r['mol_base']} + H⁺")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. AH  : {r['conf_acid']}")
    print(f"  Conf. A⁻  : {r['conf_base']}")
    print(f"  SP level (AH)  : {r['theory_sp_acid'] or '—'}")
    print(f"  SP level (A⁻)  : {r['theory_sp_base'] or '—'}")
    if aq_label_bc != "CPCM":
        print(f"  Aqueous SP source: {aq_label_bc}  (Born-cycle constants "
              f"calibrated for CPCM — see warnings below)")
    print()
    print(f"  {'Term':<40} {'Value':>12}")
    print(f"  {'─'*40}  {'─'*12}")

    def fmt(v): return f"{v:+12.3f}" if v is not None else "       N/A"

    print(f"  {'G_gas(AH)  [Eh]':<40} {r['G_gas_acid_Eh']:>12.8f}")
    print(f"  {'G_gas(A⁻)  [Eh]':<40} {r['G_gas_base_Eh']:>12.8f}")
    print(f"  {'G(H⁺, gas)  [kcal/mol]':<40} {G_HPLUS_GAS:>12.3f}")
    print()
    print(f"  {'ΔG_gas  [kcal/mol]':<40} {fmt(r['dG_gas_kcal'])}")
    print(f"  {'ΔG_solv(AH)  [kcal/mol]':<40} {fmt(r['dG_solv_acid_kcal'])}")
    print(f"  {'ΔG_solv(A⁻)  [kcal/mol]':<40} {fmt(r['dG_solv_base_kcal'])}")
    print(f"  {'ΔG_solv(H⁺) [Kelly 2006]  [kcal/mol]':<40} {DG_SOLV_HPLUS:>12.3f}")
    print(f"  {'─'*40}  {'─'*12}")
    print(f"  {'ΔG_aq  [kcal/mol]':<40} {fmt(r['dG_aq_kcal'])}")
    print(f"  {'ΔG_aq  [kJ/mol]':<40} {r['dG_aq_kcal']*4.184:>12.3f}" if r['dG_aq_kcal'] else "")
    print()
    pka_str = f"{r['pKa']:.2f}" if r['pKa'] is not None else "N/A"
    print(f"  pKa = ΔG_aq / (RT·ln10)  =  {pka_str}")
    print(f"  (typical uncertainty ±1–2 units; Born/{aq_label_bc}/TZVPD)")
    if r["warnings"]:
        print(f"\n  Warnings:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


def _print_direct_cycle(r: dict):
    """Imprime las contribuciones del método directo (solo energía acuosa, sin fase gas)."""
    aq_acid = r.get("aq_source_acid", "CPCM")
    aq_base = r.get("aq_source_base", "CPCM")
    aq_label = aq_acid if aq_acid == aq_base else f"{aq_acid}/{aq_base}"
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  DIRECT METHOD ({aq_label} only) — {r['mol_acid']} → {r['mol_base']} + H⁺")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. AH  : {r['conf_acid']}")
    print(f"  Conf. A⁻  : {r['conf_base']}")
    print(f"  SP level (AH)  : {r['theory_sp_acid'] or '—'}")
    print(f"  SP level (A⁻)  : {r['theory_sp_base'] or '—'}")
    print()
    print(f"  {f'G_total_{aq_acid}(AH)  [Eh]':<40} {r['G_total_acid_Eh']:>12.8f}")
    print(f"  {f'G_total_{aq_base}(A⁻)  [Eh]':<40} {r['G_total_base_Eh']:>12.8f}")
    print(f"  {'G(H⁺, aq)  [kcal/mol]':<40} {r['G_hplus_aq_kcal']:>12.3f}"
          f"   [Pliego & Riveros 2002, JPCA 106, 7434]")
    print(f"  {'─'*40}  {'─'*12}")
    dG = r["dG_aq_kcal"]
    print(f"  {'ΔG_aq  [kcal/mol]':<40} {dG:+12.3f}" if dG is not None else "")
    print(f"  {'ΔG_aq  [kJ/mol]':<40} {dG*4.184:>12.3f}" if dG is not None else "")
    print()
    pka_str = f"{r['pKa']:.2f}" if r['pKa'] is not None else "N/A"
    print(f"  pKa = ΔG_aq / (RT·ln10)  =  {pka_str}")
    print(f"  (direct method: valid if the solvation error of the theory")
    print(f"  level is consistent across species; typical uncertainty ±1–2")
    print(f"  units, similar to the full Born cycle)")
    if r["warnings"]:
        print(f"\n  Warnings:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


def _print_born_cycle_macrostate(r: dict):
    """
    Imprime el ciclo de Born cuando uno de los dos lados es un macroestado
    combinado neutro+zwitterión. No hay G_gas/ΔG_solv por especie que mostrar
    (no se calculan ni se necesitan — ver la nota en compute_pka()): en vez
    de eso, se muestra explícitamente la reducción algebraica a la resta de
    energías CPCM con la constante G_HPLUS_AQ_BORN.
    """
    aq_acid_ms = r.get("aq_source_acid", "CPCM")
    aq_base_ms = r.get("aq_source_base", "CPCM")
    aq_label_ms = aq_acid_ms if aq_acid_ms == aq_base_ms else f"{aq_acid_ms}/{aq_base_ms}"
    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  BORN CYCLE (macrostate, algebraic reduction) — "
          f"{r['mol_acid']} → {r['mol_base']} + H⁺")
    print(f"  T = {r['T']:.2f} K")
    print(f"  {sep}")
    print(f"  Conf. AH  : {r['conf_acid']}")
    print(f"  Conf. A⁻  : {r['conf_base']}")
    print(f"  SP level (AH)  : {r['theory_sp_acid'] or '—'}")
    print(f"  SP level (A⁻)  : {r['theory_sp_base'] or '—'}")
    print()
    print(f"  No explicit gas-phase SP is used for the mixture — the Born")
    print(f"  cycle reduces algebraically to the {aq_label_ms}-only energy difference")
    print(f"  (same E_gas-cancellation already present for a single species")
    print(f"  in this script; see the derivation note above")
    print(f"  compute_pka_montecarlo()).")
    print()
    print(f"  {'G_total_' + aq_acid_ms + '(AH)  [Eh]':<40} {r['G_total_acid_Eh']:>12.8f}")
    print(f"  {'G_total_' + aq_base_ms + '(A⁻)  [Eh]':<40} {r['G_total_base_Eh']:>12.8f}")
    print(f"  {'G(H⁺, aq)  [kcal/mol]':<40} {r['G_hplus_aq_kcal']:>12.3f}"
          f"   [Tissandier 1998 + Kelly 2006]")
    print(f"  {'─'*40}  {'─'*12}")
    dG = r["dG_aq_kcal"]
    print(f"  {'ΔG_aq  [kcal/mol]':<40} {dG:+12.3f}" if dG is not None else "")
    print(f"  {'ΔG_aq  [kJ/mol]':<40} {dG*4.184:>12.3f}" if dG is not None else "")
    print()
    pka_str = f"{r['pKa']:.2f}" if r['pKa'] is not None else "N/A"
    print(f"  pKa = ΔG_aq / (RT·ln10)  =  {pka_str}")
    print(f"  (typical uncertainty ±1–2 units; Born/{aq_label_ms}/TZVPD)")
    if r["warnings"]:
        print(f"\n  Warnings:")
        for w in r["warnings"]:
            print(f"    ⚠  {w}")
    print(f"  {sep}\n")


# ─────────────────────────────────────────────
#  DISTRIBUCIÓN DE BOLTZMANN ENTRE ESTADOS
# ─────────────────────────────────────────────

def boltzmann_protonation_states(
    pKa1: Optional[float],
    pKa2: Optional[float],
    pH_values: list[float],
    T: float = T_STD,
) -> pd.DataFrame:
    """
    Calcula la distribución de especies entre estados de protonación
    a múltiples valores de pH usando Henderson-Hasselbalch encadenado.

    Para un diácido con pKa1 (AH₂ → AH⁻) y pKa2 (AH⁻ → A²⁻):

      f(AH₂) = 1 / D
      f(AH⁻) = 10^(pH − pKa1) / D
      f(A²⁻) = 10^(2·pH − pKa1 − pKa2) / D
      D = 1 + 10^(pH − pKa1) + 10^(2·pH − pKa1 − pKa2)

    Si solo hay pKa1 (sin pKa2 / sin dianión):
      f(AH)  = 1 / (1 + 10^(pH − pKa1))
      f(A⁻)  = 10^(pH − pKa1) / (1 + 10^(pH − pKa1))

    Retorna DataFrame con columnas:
      pH, f_neutral, f_anion, f_dianion (si pKa2 disponible)
    """
    rows = []
    for pH in pH_values:
        if pKa1 is not None and pKa2 is not None:
            a0 = 1.0
            a1 = 10 ** (pH - pKa1)
            a2 = 10 ** (2 * pH - pKa1 - pKa2)
            D  = a0 + a1 + a2
            rows.append({
                "pH":        pH,
                "f_neutral": a0 / D,
                "f_anion":   a1 / D,
                "f_dianion": a2 / D,
            })
        elif pKa1 is not None:
            a0 = 1.0
            a1 = 10 ** (pH - pKa1)
            D  = a0 + a1
            rows.append({
                "pH":        pH,
                "f_neutral": a0 / D,
                "f_anion":   a1 / D,
            })
        else:
            rows.append({"pH": pH})

    return pd.DataFrame(rows)


def boltzmann_5state(
    pKa0: Optional[float],
    pKa1: Optional[float],
    pKa2: Optional[float],
    f_R: Optional[float],
    f_D: Optional[float],
    pH_values: list[float],
) -> pd.DataFrame:
    """
    Distribución de Henderson-Hasselbalch encadenada para la cadena completa:

        A⁺ (catión)  ⇌  [R (neutro) ⇌ D (zwitterión)]  ⇌  B⁻ (anión)  ⇌  C²⁻ (dianión)
              pKa0                  (f_R / f_D, sin pH)        pKa1          pKa2

    pKa0 = A → neutro/zwitterión + H⁺
    pKa1 = neutro/zwitterión → B + H⁺
    pKa2 = B → C + H⁺

    f_R, f_D: fracción de R y D DENTRO del macroestado neutro (constante con
    el pH, viene de compute_tautomer_equilibrium). Si no se proveen, se asume
    f_R = 1, f_D = 0 (todo el macroestado se reporta como "neutro").

    Cualquier pKa que sea None se trata como "ese estado no existe" (se omite
    su rama de la cadena, igual que boltzmann_protonation_states).

    Retorna DataFrame con columnas: pH, f_cation, f_R, f_D, f_neutral_macro,
    f_anion, f_dianion (las que apliquen).
    """
    if f_R is None or f_D is None:
        f_R, f_D = 1.0, 0.0

    def _pow10(exp):
        # clip para evitar overflow cuando un pKa está muy lejos del rango de pH
        # (con pKa's físicamente razonables esto nunca se activa)
        return 10.0 ** np.clip(exp, -300.0, 300.0)

    rows = []
    for pH in pH_values:
        r0 = _pow10(pKa0 - pH) if pKa0 is not None else 0.0
        r1 = _pow10(pH - pKa1) if pKa1 is not None else 0.0
        r2 = (_pow10(2 * pH - pKa1 - pKa2)
              if (pKa1 is not None and pKa2 is not None) else 0.0)

        D = 1.0 + r0 + r1 + r2
        f_neutral_macro = 1.0 / D
        row = {
            "pH":               pH,
            "f_neutral_macro":  f_neutral_macro,
            "f_R":              f_neutral_macro * f_R,
            "f_D":              f_neutral_macro * f_D,
        }
        if pKa0 is not None:
            row["f_cation"] = r0 / D
        if pKa1 is not None:
            row["f_anion"] = r1 / D
        if pKa1 is not None and pKa2 is not None:
            row["f_dianion"] = r2 / D
        rows.append(row)

    return pd.DataFrame(rows)


# ─────────────────────────────────────────────
#  MONTE CARLO — CADENA COMPLETA DE 5 ESTADOS (FASE 2)
# ─────────────────────────────────────────────

def _boltzmann_5state_vectorized(
    pKa0_samples: Optional[np.ndarray],
    pKa1_samples: Optional[np.ndarray],
    pKa2_samples: Optional[np.ndarray],
    f_R: Optional[float],
    f_D: Optional[float],
    pH_values: list[float],
) -> pd.DataFrame:
    """
    Versión vectorizada de boltzmann_5state() sobre un eje adicional de
    muestras Monte Carlo. En vez de llamar boltzmann_5state() una vez por
    cada una de las N muestras (con N=10 000 y ~140 valores de pH serían
    ~1.4M evaluaciones en loop de Python puro), arma arrays 2D
    (n_samples, n_pH) con broadcasting de numpy y reduce a percentiles por
    columna de pH al final.

    Reproduce exactamente la misma fórmula que boltzmann_5state() (mismo
    r0/r1/r2/D) — si esa función cambia, esta hay que actualizarla en
    paralelo; no hay una única fuente de verdad aquí porque una versión
    escalar con loop sobre N=10⁴ muestras sería demasiado lenta para uso
    interactivo.

    Parameters
    ----------
    pKa0_samples, pKa1_samples, pKa2_samples : np.ndarray o None
        Arrays de muestras MC (mismo n_samples cada uno, generados con el
        mismo shared_hplus_draw si vienen de la misma cadena — ver
        run_chain_montecarlo()). None si ese pKa no existe en la cadena
        (mismo criterio que boltzmann_5state: se omite esa rama).
    f_R, f_D : float, opcional
        Igual que en boltzmann_5state(); default 1.0/0.0 si no se proveen.
    pH_values : list[float]

    Returns
    -------
    pd.DataFrame con columna 'pH' y, para cada fracción que aplique
    (f_neutral_macro, f_R, f_D, f_cation, f_anion, f_dianion), tres
    columnas: '{frac}_lo' (percentil 2.5), '{frac}_median', '{frac}_hi'
    (percentil 97.5). DataFrame vacío si no hay ninguna muestra disponible.
    """
    if f_R is None or f_D is None:
        f_R, f_D = 1.0, 0.0

    have_pka0 = pKa0_samples is not None
    have_pka1 = pKa1_samples is not None
    have_pka2 = have_pka1 and (pKa2_samples is not None)

    if not have_pka0 and not have_pka1:
        return pd.DataFrame()

    n_samples = (
        len(pKa0_samples) if have_pka0 else
        len(pKa1_samples)
    )
    for name, arr in [("pKa0", pKa0_samples), ("pKa1", pKa1_samples), ("pKa2", pKa2_samples)]:
        if arr is not None and len(arr) != n_samples:
            raise ValueError(
                f"{name}_samples tiene {len(arr)} muestras, se esperaban "
                f"{n_samples} (todas las ramas de la cadena deben compartir "
                f"el mismo n_samples)."
            )

    pH_arr = np.asarray(pH_values, dtype=float)             # shape (M,)

    def _pow10_2d(exp_2d):
        return 10.0 ** np.clip(exp_2d, -300.0, 300.0)

    zeros_2d = np.zeros((n_samples, len(pH_arr)))

    r0 = (_pow10_2d(pKa0_samples[:, None] - pH_arr[None, :])
          if have_pka0 else zeros_2d)
    r1 = (_pow10_2d(pH_arr[None, :] - pKa1_samples[:, None])
          if have_pka1 else zeros_2d)
    r2 = (_pow10_2d(2 * pH_arr[None, :] - pKa1_samples[:, None] - pKa2_samples[:, None])
          if have_pka2 else zeros_2d)

    D = 1.0 + r0 + r1 + r2
    f_neutral_macro = 1.0 / D
    f_R_arr = f_neutral_macro * f_R
    f_D_arr = f_neutral_macro * f_D

    def _summ(arr_2d):
        lo, med, hi = np.percentile(arr_2d, [2.5, 50.0, 97.5], axis=0)
        return lo, med, hi

    cols: dict = {"pH": pH_arr}

    for label, arr_2d in [
        ("f_neutral_macro", f_neutral_macro),
        ("f_R", f_R_arr),
        ("f_D", f_D_arr),
    ]:
        lo, med, hi = _summ(arr_2d)
        cols[f"{label}_lo"], cols[f"{label}_median"], cols[f"{label}_hi"] = lo, med, hi

    if have_pka0:
        lo, med, hi = _summ(r0 / D)
        cols["f_cation_lo"], cols["f_cation_median"], cols["f_cation_hi"] = lo, med, hi
    if have_pka1:
        lo, med, hi = _summ(r1 / D)
        cols["f_anion_lo"], cols["f_anion_median"], cols["f_anion_hi"] = lo, med, hi
    if have_pka2:
        lo, med, hi = _summ(r2 / D)
        cols["f_dianion_lo"], cols["f_dianion_median"], cols["f_dianion_hi"] = lo, med, hi

    return pd.DataFrame(cols)


def run_chain_montecarlo(
    data_A: Optional[dict],
    data_neutral_macro: Optional[dict],
    data_B: Optional[dict],
    data_C: Optional[dict],
    f_R: Optional[float],
    f_D: Optional[float],
    pH_values: list[float],
    T: float = T_STD,
    method: str = "auto",
    gcorr_source: str = "harmonic",
    n_samples: int = 10_000,
    sigma_elec_kcal: float = 0.5,
    sigma_gcorr_kcal: float = 0.4,
    sigma_hplus_kcal: float = 2.5,
    seed: Optional[int] = None,
) -> dict:
    """
    Propaga incertidumbre por Monte Carlo a través de TODA la cadena de 5
    estados (A⁺ ⇌ [R⇌D] ⇌ B⁻ ⇌ C²⁻), preservando la correlación física de
    G(H⁺,aq) entre pKa0/pKa1/pKa2, y devuelve bandas de confianza para las
    fracciones de especie vs. pH.

    Mismos pares acid/base que usa _run_extended() para la cadena
    determinística:
        pKa0 = compute_pka(data_A, data_neutral_macro, ...)
        pKa1 = compute_pka(data_neutral_macro, data_B,  ...)
        pKa2 = compute_pka(data_B, data_C, ...)

    Cualquier par con un miembro en None se omite (mismo criterio que
    _run_extended: pKa0 requiere data_A, pKa2 requiere data_B y data_C).

    Punto clave de diseño: la perturbación de G(H⁺,aq) (shared_hplus_draw)
    se genera UNA sola vez acá y se pasa idéntica a las tres llamadas de
    compute_pka_montecarlo. Esto hace que, en las N muestras generadas,
    pKa0/pKa1/pKa2 se muevan de forma correlacionada cuando la constante del
    protón se desplaza — que es lo físicamente correcto, porque es la MISMA
    constante en los tres pasos de desprotonación. El ruido de energía
    electrónica/térmica (sigma_elec_kcal, sigma_gcorr_kcal) en cambio usa
    generadores independientes por tramo (rng.spawn), porque A, R/D, B y C
    son cálculos ORCA distintos con errores no correlacionados entre sí.

    Returns
    -------
    dict con:
        ok             : bool
        mc_pka0/1/2    : dict o None — salida completa de
                         compute_pka_montecarlo() para cada tramo presente
        df_dist_bands  : pd.DataFrame — salida de
                         _boltzmann_5state_vectorized() (columnas
                         {frac}_lo/_median/_hi por cada fracción aplicable)
        warnings       : list[str]
    """
    rng = np.random.default_rng(seed)

    # Constante compartida por toda la cadena — un solo draw, reusado en
    # las tres llamadas de compute_pka_montecarlo de abajo.
    shared_hplus_draw = rng.normal(0.0, sigma_hplus_kcal, size=n_samples)

    # Sub-generadores independientes para el ruido de E_elec/G_corr de cada
    # tramo (no deben compartirse entre tramos: cada uno involucra cálculos
    # ORCA distintos).
    rng_pka0, rng_pka1, rng_pka2 = rng.spawn(3)

    out: dict = {
        "ok": False,
        "n_samples": n_samples,
        "mc_pka0": None, "mc_pka1": None, "mc_pka2": None,
        "df_dist_bands": None,
        "warnings": [],
    }

    mc_pka0 = mc_pka1 = mc_pka2 = None

    if data_A is not None and data_neutral_macro is not None:
        mc_pka0 = compute_pka_montecarlo(
            data_A, data_neutral_macro, T=T, method=method, gcorr_source=gcorr_source,
            n_samples=n_samples,
            sigma_elec_kcal=sigma_elec_kcal, sigma_gcorr_kcal=sigma_gcorr_kcal,
            sigma_hplus_kcal=sigma_hplus_kcal,
            shared_hplus_draw=shared_hplus_draw, rng=rng_pka0,
        )
        out["warnings"].extend(mc_pka0["warnings"])

    if data_neutral_macro is not None and data_B is not None:
        mc_pka1 = compute_pka_montecarlo(
            data_neutral_macro, data_B, T=T, method=method, gcorr_source=gcorr_source,
            n_samples=n_samples,
            sigma_elec_kcal=sigma_elec_kcal, sigma_gcorr_kcal=sigma_gcorr_kcal,
            sigma_hplus_kcal=sigma_hplus_kcal,
            shared_hplus_draw=shared_hplus_draw, rng=rng_pka1,
        )
        out["warnings"].extend(mc_pka1["warnings"])

    if data_B is not None and data_C is not None:
        mc_pka2 = compute_pka_montecarlo(
            data_B, data_C, T=T, method=method, gcorr_source=gcorr_source,
            n_samples=n_samples,
            sigma_elec_kcal=sigma_elec_kcal, sigma_gcorr_kcal=sigma_gcorr_kcal,
            sigma_hplus_kcal=sigma_hplus_kcal,
            shared_hplus_draw=shared_hplus_draw, rng=rng_pka2,
        )
        out["warnings"].extend(mc_pka2["warnings"])

    out["mc_pka0"], out["mc_pka1"], out["mc_pka2"] = mc_pka0, mc_pka1, mc_pka2

    if (mc_pka0 is None or not mc_pka0["ok"]) and (mc_pka1 is None or not mc_pka1["ok"]):
        out["warnings"].append(
            "No hay suficientes tramos válidos para construir la cadena MC "
            "(se necesita al menos A+neutro/zwitterión, o neutro/zwitterión+B, "
            "con compute_pka() determinístico exitoso)."
        )
        return out

    pKa0_samples = mc_pka0["pka_samples"] if (mc_pka0 and mc_pka0["ok"]) else None
    pKa1_samples = mc_pka1["pka_samples"] if (mc_pka1 and mc_pka1["ok"]) else None
    pKa2_samples = mc_pka2["pka_samples"] if (mc_pka2 and mc_pka2["ok"]) else None

    df_bands = _boltzmann_5state_vectorized(
        pKa0_samples, pKa1_samples, pKa2_samples, f_R, f_D, pH_values,
    )

    out["ok"] = True
    out["df_dist_bands"] = df_bands
    return out


def print_mc_summary(chain_mc: Optional[dict], pka_labels: list) -> None:
    """
    Imprime tabla comparativa pKa determinístico vs. Monte Carlo (media±std,
    CI95) — análoga a print_pka_summary(), pero para la salida de
    run_chain_montecarlo().

    pka_labels: lista de tuplas (key, resultado_dict_o_None, etiqueta_str),
    p. ej. [("pKa0", pka0_result, "pKa₀ (A→R/D)"), ("pKa1", pka1_result,
    "pKa₁ (R/D→B)"), ("pKa2", pka2_result, "pKa₂ (B→C)")] — key debe
    coincidir con mc_pka0/mc_pka1/mc_pka2 dentro de chain_mc.
    """
    sep = "═" * 74
    print(f"\n{sep}")
    print("  PROPAGACIÓN DE INCERTIDUMBRE — MONTE CARLO")
    print(sep)

    if chain_mc is None or not chain_mc.get("ok"):
        print("\n  ⚠  El Monte Carlo no pudo ejecutarse.")
        for w in (chain_mc or {}).get("warnings", []):
            print(f"    ⚠  {w}")
        print(f"{sep}\n")
        return

    mc_map = {
        "pKa0": chain_mc.get("mc_pka0"),
        "pKa1": chain_mc.get("mc_pka1"),
        "pKa2": chain_mc.get("mc_pka2"),
    }

    print(f"\n  N = {chain_mc['n_samples']:,} muestras  "
          f"(G(H⁺,aq) compartido entre pKa0/pKa1/pKa2 — ver nota abajo)")
    print(f"\n  {'':<22} {'determinístico':>14} {'media MC':>10} "
          f"{'± std':>8} {'CI 95%':>18}")
    print(f"  {'-'*22} {'-'*14} {'-'*10} {'-'*8} {'-'*18}")

    any_row = False
    for key, det, label in pka_labels:
        mc = mc_map.get(key)
        if det is None and mc is None:
            continue
        any_row = True
        det_v = f"{det['pKa']:.2f}" if (det and det.get("ok")) else "—"
        if mc is not None and mc.get("ok"):
            mean_v = f"{mc['pka_mean']:.2f}"
            std_v  = f"±{mc['pka_std']:.2f}"
            lo, hi = mc["pka_ci95"]
            ci_v   = f"[{lo:.2f}, {hi:.2f}]"
        else:
            mean_v = std_v = ci_v = "—"
        print(f"  {label:<22} {det_v:>14} {mean_v:>10} {std_v:>8} {ci_v:>18}")

    if not any_row:
        print("  (sin pKa's disponibles para Monte Carlo)")

    sample_mc = next((mc for mc in mc_map.values() if mc is not None and mc.get("ok")), None)
    if sample_mc is not None:
        print(f"\n  σ asumidos: E_elec=±{sample_mc['sigma_elec_kcal']:.2f}  "
              f"G_corr=±{sample_mc['sigma_gcorr_kcal']:.2f}  "
              f"G(H⁺,aq)=±{sample_mc['sigma_hplus_kcal']:.2f} kcal/mol")

    if chain_mc["warnings"]:
        print(f"\n  ⚠  Advertencias:")
        for w in chain_mc["warnings"]:
            print(f"     ⚠  {w}")

    print(f"\n  Nota: G(H⁺,aq) es la misma constante física en los tres pasos")
    print(f"  de desprotonación, así que su perturbación se comparte entre")
    print(f"  pKa0/pKa1/pKa2 en cada muestra MC (no se samplea independiente")
    print(f"  por tramo) — esto correlaciona el sesgo sistemático entre pKa's")
    print(f"  y evita inflar artificialmente la incertidumbre de cantidades")
    print(f"  derivadas como Δ(pKa2−pKa1) o la fracción de especie a un pH")
    print(f"  dado. Las bandas de confianza de la distribución de especies")
    print(f"  quedan en df_dist_bands (percentiles 2.5/50/97.5).")
    print(f"{sep}\n")


# ─────────────────────────────────────────────
#  IMPRESIÓN DE RESULTADOS
# ─────────────────────────────────────────────

def print_mol_summary(data: dict, verbose: bool = False):
    """Imprime resumen del ensemble conformacional de una molécula."""
    mol = data["mol_name"]
    df  = data["df"]
    sep = "─" * 74

    print(f"\n{sep}")
    print(f"  ESTADO: {mol}")
    print(sep)

    if df.empty:
        print("  ⚠  Sin datos.")
        return

    q, mult, n_at = data.get("charge"), data.get("mult"), data.get("n_atoms")
    q_str  = f"{q:+d}"  if q  is not None else "?"
    m_str  = f"{mult}"  if mult is not None else "?"
    n_str  = f"{n_at}"  if n_at is not None else "?"
    print(f"  Carga: {q_str}   Multiplicidad: {m_str}   N° átomos: {n_str}")

    # Tabla de confórmeros
    cols_base = ["conformer", "E_sp_aq", "thermo_ok", "n_imag"]
    cols_g    = [c for c in df.columns if c.startswith("G_total_")]
    cols_b    = [c for c in df.columns if c.startswith("boltz_")]
    cols_dg   = [c for c in df.columns if c.startswith("dG_")]

    print(f"\n  Confórmeros SP (acuoso, {data.get('aq_source', 'CPCM')}):")
    disp = df[cols_base + cols_g[:2] + cols_b[:2] + cols_dg[:2]].copy()
    print(disp.to_string(index=False, float_format=lambda x: f"{x:.5f}"))

    print(f"\n  Confórmero representativo (min G_total @ 298K): {data['best_conf']}")
    print(f"\n  G_Boltzmann del ensemble (energía libre efectiva, Hartree):")
    for lbl, G in data["G_boltz"].items():
        if G is not None:
            print(f"    {lbl}: {G:.8f} Eh")
        else:
            print(f"    {lbl}: N/A")

    for w in data["warnings"]:
        print(f"\n  ⚠  {w}")


def _print_single_pka(pka: dict, label: str):
    """Imprime un único resultado de pKa con su etiqueta (pKa₀, pKa₁, pKa₂, ...)."""
    if pka is None:
        return
    if pka["ok"]:
        print(f"\n  {label}  ({pka['mol_acid']} → {pka['mol_base']} + H⁺)"
              f"   [método: {pka['method']}]")
        if pka["method"] == "born" and pka.get("dG_gas_kcal") is not None:
            print(f"  ΔG_gas        = {pka['dG_gas_kcal']:+.3f} kcal/mol")
            print(f"  ΔG_solv(ácido) = {pka['dG_solv_acid_kcal']:+.3f} kcal/mol")
            print(f"  ΔG_solv(base)  = {pka['dG_solv_base_kcal']:+.3f} kcal/mol")
            print(f"  ΔG_solv(H⁺)   = {DG_SOLV_HPLUS:+.3f} kcal/mol  [Kelly 2006]")
        elif pka["method"] == "born":
            # Macroestado: ciclo de Born resuelto vía la reducción algebraica
            # (sin SP de fase gas para la mezcla) — ver compute_pka().
            print(f"  G(H⁺, aq)     = {pka['G_hplus_aq_kcal']:+.3f} kcal/mol "
                  f"[Tissandier 1998 + Kelly 2006; macroestado, reducción algebraica]")
        else:
            print(f"  G(H⁺, aq)     = {pka['G_hplus_aq_kcal']:+.3f} kcal/mol "
                  f"[Pliego & Riveros 2002]")
        print(f"  ΔG_aq         = {pka['dG_aq_kcal']:+.3f} kcal/mol")
        print(f"\n  ━━►  {label}  =  {pka['pKa']:.2f}  (T = {pka['T']:.2f} K)")
        if pka["warnings"]:
            print(f"\n  ⚠  Advertencias (revisar antes de reportar):")
            for w in pka["warnings"]:
                print(f"     ⚠  {w}")
    else:
        print(f"\n  {label}: NO calculado.")
        for w in pka["warnings"]:
            print(f"    ⚠  {w}")


def print_pka_summary(pka_list: list):
    """
    Imprime resumen de todos los pKa's calculados.

    pka_list: lista de tuplas (resultado_dict_o_None, etiqueta_str), en orden,
    p. ej. [(pka0, "pKa₀ (A→R/D)"), (pka1, "pKa₁ (R/D→B)"), (pka2, "pKa₂ (B→C)")]
    """
    sep = "═" * 74
    print(f"\n{sep}")
    print(f"  RESULTADOS FINALES — CICLOS DE pKa")
    print(sep)

    for pka, label in pka_list:
        if pka is not None:
            _print_single_pka(pka, label)

    print(f"\n  Nota: incertidumbre típica ±1–2 unidades de pKa para métodos")
    print(f"  Born/CPCM a nivel ωB97X-D3/def2-TZVPD (Ho & Coote 2010). El método")
    print(f"  'direct' usa solo energías CPCM (sin SP fase gas) con G(H⁺,aq) de")
    print(f"  Pliego & Riveros 2002 — esta constante NO es intercambiable con la")
    print(f"  de Tissandier/Kelly usada en 'born', que está calibrada para SP en")
    print(f"  fase gas pura. Cada método usa su propia referencia del protón.")
    print(f"  Los pasos que involucran el macroestado neutro/zwitterión combinado")
    print(f"  resuelven method='auto' a 'born' (exacto, sin SP fase gas — ver")
    print(f"  advertencias arriba). Con --method direct explícito usan 'direct'.")
    print(f"{sep}\n")


def _cross_validate_ph_scale(
    data_acid: dict,
    data_base: dict,
    result: dict,
    T: float = T_STD,
    gcorr_source: str = "harmonic",
) -> dict:
    """
    Diagnóstico de consistencia de escala de energía libre del protón.

    Problema que resuelve: con method='auto', pKa₁ (R/D→B, macroestado) y
    pKa₂ (B→C, simple) pueden resolver a métodos distintos. pKa₁ resuelve a
    'born' siempre que involucre el macroestado (exacto ahí, sin necesitar
    SP fase gas — ver rama is_macrostate en compute_pka), mientras que pKa₂
    depende de si existen SP fase gas para B y C. 'direct' y 'born' usan
    constantes G(H⁺,aq) distintas y NO intercambiables (Pliego & Riveros
    −264.61 kcal/mol vs. Tissandier/Kelly −272.18 kcal/mol) — sumar un pKa
    'direct' con un pKa 'born' en la misma cadena mezcla dos escalas de
    referencia del protón. Bajo method='auto' esto ahora se trata como error
    duro en el caller (ver chequeo tras _cross_validate_ph_scale más abajo),
    no solo como diagnóstico.

    Esta función NO reemplaza el resultado ya calculado (`result`); solo
    lo audita. Si `result['method'] == 'born'`, recalcula el mismo par
    ácido/base forzando method='direct' (reusando compute_pka, que ya
    implementa ambas rutas) y reporta el delta entre ambos pKa como
    diagnóstico. Si `result['method']` ya es 'direct', no hay nada que
    cross-validar (ambos pasos ya comparten la misma escala) y se retorna
    ok=False sin advertencia.
    """
    diag: dict = {
        "ok": False, "applicable": False,
        "pKa_born": None, "pKa_direct": None,
        "delta_pKa": None, "delta_dG_kcal": None,
        "warnings": [],
    }

    if result is None or not result.get("ok"):
        return diag

    if result.get("method") != "born":
        # Ya está en 'direct': mismo pie que pKa₁/pKa₀, nada que cruzar.
        return diag

    diag["applicable"] = True
    diag["pKa_born"] = result["pKa"]

    direct_result = compute_pka(
        data_acid, data_base, T=T, verbose=False,
        method="direct", gcorr_source=gcorr_source,
    )

    if not direct_result.get("ok"):
        diag["warnings"].append(
            "No se pudo recalcular con method='direct' para cross-validar "
            "contra el resultado 'born' (revisar warnings del sub-cálculo)."
        )
        return diag

    diag["ok"] = True
    diag["pKa_direct"] = direct_result["pKa"]
    diag["delta_pKa"] = result["pKa"] - direct_result["pKa"]
    diag["delta_dG_kcal"] = result["dG_aq_kcal"] - direct_result["dG_aq_kcal"]
    return diag


def _print_ph_scale_diagnostic(diag: dict, label: str = "pKa₂ (B→C)"):
    """Imprime el resultado de _cross_validate_ph_scale(), si aplica."""
    if not diag.get("applicable"):
        return

    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  DIAGNÓSTICO — consistencia de escala del protón ({label})")
    print(f"  {sep}")

    if not diag["ok"]:
        print(f"  ⚠  No se pudo completar el cross-check:")
        for w in diag["warnings"]:
            print(f"     {w}")
        print(f"  {sep}\n")
        return

    print(f"  {'method=\"born\"  -> pKa':<40} {diag['pKa_born']:>10.2f}")
    print(f"  {'method=\"direct\" (forzado) -> pKa':<40} {diag['pKa_direct']:>10.2f}")
    print(f"  {'Delta pKa (born - direct)':<40} {diag['delta_pKa']:>+10.2f}")
    print(f"  {'Delta dG_aq (born - direct)  [kcal/mol]':<40} {diag['delta_dG_kcal']:>+10.2f}")
    print()
    if abs(diag["delta_pKa"]) > 2.0:
        print(f"  \u26a0  Delta pKa > 2 unidades: born y direct discrepan mas de lo")
        print(f"     tipicamente esperable entre metodos. Si este resultado se")
        print(f"     va a combinar con pKa1 (que siempre usa 'direct'), evaluar")
        print(f"     reportar tambien la version 'direct' de este paso para que")
        print(f"     toda la cadena comparta la misma referencia de G(H+,aq).")
    else:
        print(f"  Delta dentro de lo esperable entre convenciones de G(H+,aq);")
        print(f"  no invalida el resultado, pero recordar que 'direct' es la")
        print(f"  escala que comparte el resto de la cadena (pKa0/pKa1).")
    print(f"  {sep}\n")


def _pka_temperature_scan(
    steps: list,
    gcorr_source: str = "harmonic",
    method: str = "auto",
) -> dict:
    """
    Escanea pKa en las 4 temperaturas de TEMPS_K/TEMP_LABELS (mismas T a las
    que collect_mol_data() ya tabula G_boltz por ensemble conformacional).

    Nota #8 del review: T estaba hardcodeada a T_STD=298.15 K en toda la
    orquestación (_run_classic/_run_extended), pese a que el pipeline ya
    calcula G a 4 temperaturas para el Boltzmann conformacional. Esta
    función reusa esa infraestructura vía compute_pka(T=...) para reportar
    cómo cambia cada pKa con T, sin fijar una sola T "de fábrica" — relevante
    para MAAs en agua de mar antártica (~274 K), no solo a 298.15 K estándar.

    steps: lista de tuplas (label, data_acid, data_base). Si data_acid o
    data_base es None, ese step se omite (igual que en el flujo principal).

    Retorna {label: {T_label: pka_result_dict}}.
    """
    scan: dict = {}
    for label, data_acid, data_base in steps:
        if data_acid is None or data_base is None:
            continue
        scan[label] = {}
        for T, lbl in zip(TEMPS_K, TEMP_LABELS):
            scan[label][lbl] = compute_pka(
                data_acid, data_base, T=T, verbose=False,
                method=method, gcorr_source=gcorr_source,
            )
    return scan


def _print_pka_temperature_scan(scan: dict, reported_T: float):
    """Imprime tabla pKa vs. T para cada step de _pka_temperature_scan()."""
    if not scan:
        return

    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  ESCANEO DE TEMPERATURA — pKa(T)")
    print(f"  {sep}")
    print(f"  {'Paso':<18}" + "".join(f"{lbl:>12}" for lbl in TEMP_LABELS))
    print(f"  {'─'*18}" + "".join(f"{'─'*12}" for _ in TEMP_LABELS))

    for label, by_T in scan.items():
        row = f"  {label:<18}"
        for T, lbl in zip(TEMPS_K, TEMP_LABELS):
            r = by_T.get(lbl)
            marker = "*" if np.isclose(T, reported_T, atol=1e-3) else " "
            if r and r.get("ok"):
                row += f"{r['pKa']:>11.2f}{marker}"
            else:
                row += f"{'N/A':>11}{marker}"
        print(row)

    print(f"\n  * = T reportada en el resto del output (--temp, default "
          f"{T_STD:.2f} K).")
    print(f"  Nota: el pKa incluye T tanto en G_corr (vía GoodVibes/NumFreq)")
    print(f"  como en RT·ln10 al final del ciclo; el efecto no es solo el")
    print(f"  denominador. No confundir con --goodvibes (harmonic vs.")
    print(f"  quasi-harmonic), que es un eje ortogonal a esta T.")
    print(f"  {sep}\n")


def _extended_temperature_scan(
    data_A: Optional[dict], data_R: dict, data_D: dict,
    data_B: Optional[dict], data_C: Optional[dict],
    gcorr_source: str = "harmonic",
    method: str = "auto",
) -> dict:
    """
    Análogo a _pka_temperature_scan() pero para la cadena extendida de 5
    estados, donde pKa₀/pKa₁ dependen del macroestado neutro/zwitterión
    (data_neutral_macro), que a su vez está "horneado" a una única T dentro
    de taut_result["G_macro_Eh"]. No alcanza con pasarle T distinta a
    compute_pka() para esos dos pasos — hay que reconstruir el equilibrio
    de tautómeros y el macroestado en cada T del escaneo.

    Retorna {label: {T_label: pka_result_dict}}, mismo formato que
    _pka_temperature_scan(), más "f_D" (fracción de zwitterión por T) para
    contexto, ya que también depende de T.
    """
    scan: dict = {"f_D": {}}

    for T, lbl in zip(TEMPS_K, TEMP_LABELS):
        taut_T = compute_tautomer_equilibrium(data_R, data_D, T=T, verbose=False,
                                              gcorr_source=gcorr_source)
        if not taut_T["ok"]:
            gR_T = _state_G_total_direct(data_R, T=T, gcorr_source=gcorr_source)
            if not gR_T.get("ok"):
                scan["f_D"][lbl] = None
                continue
            taut_T = {
                "f_R": 1.0, "f_D": 0.0,
                "G_macro_Eh": gR_T["G_total_Eh"],
                "theory_R": gR_T["theory_sp"], "theory_D": None,
            }

        scan["f_D"][lbl] = taut_T["f_D"]

        macro_T = {
            "mol_name":          "neutro/zwitterión",
            "is_macrostate":     True,
            "G_total_direct_Eh": taut_T["G_macro_Eh"],
            "conf_label":        f"macroestado @ {lbl}",
            "theory_label":      taut_T["theory_R"] or taut_T["theory_D"],
            "charge":            data_R.get("charge"),
            "n_atoms":           data_R.get("n_atoms"),
            "use_smd":           data_R.get("use_smd", False),
            "aq_source":         data_R.get("aq_source", "CPCM"),
        }

        if data_A is not None:
            scan.setdefault("pKa0", {})[lbl] = compute_pka(
                data_A, macro_T, T=T, verbose=False,
                method=method, gcorr_source=gcorr_source,
            )
        if data_B is not None:
            scan.setdefault("pKa1", {})[lbl] = compute_pka(
                macro_T, data_B, T=T, verbose=False,
                method=method, gcorr_source=gcorr_source,
            )
        if data_B is not None and data_C is not None:
            scan.setdefault("pKa2", {})[lbl] = compute_pka(
                data_B, data_C, T=T, verbose=False,
                method=method, gcorr_source=gcorr_source,
            )

    return scan


def _print_extended_temperature_scan(scan: dict, reported_T: float):
    """Imprime tabla pKa(T) + f_D(T) para la cadena extendida de 5 estados."""
    pka_labels = [k for k in ("pKa0", "pKa1", "pKa2") if k in scan]
    if not pka_labels:
        return

    sep = "─" * 66
    print(f"\n  {sep}")
    print(f"  ESCANEO DE TEMPERATURA — pKa(T) y f_D(T)  (cadena extendida)")
    print(f"  {sep}")
    print(f"  {'Paso':<18}" + "".join(f"{lbl:>12}" for lbl in TEMP_LABELS))
    print(f"  {'─'*18}" + "".join(f"{'─'*12}" for _ in TEMP_LABELS))

    row = f"  {'f_D (zwitt.)':<18}"
    for lbl in TEMP_LABELS:
        fD = scan["f_D"].get(lbl)
        row += f"{fD*100:>10.1f}% " if fD is not None else f"{'N/A':>11} "
    print(row)

    for label in pka_labels:
        row = f"  {label:<18}"
        for T, lbl in zip(TEMPS_K, TEMP_LABELS):
            r = scan[label].get(lbl)
            marker = "*" if np.isclose(T, reported_T, atol=1e-3) else " "
            if r and r.get("ok"):
                row += f"{r['pKa']:>11.2f}{marker}"
            else:
                row += f"{'N/A':>11}{marker}"
        print(row)

    print(f"\n  * = T reportada en el resto del output (--temp, default "
          f"{T_STD:.2f} K).")
    print(f"  f_D depende de T (equilibrio R⇌D) y se re-cierra el macroestado")
    print(f"  neutro/zwitterión en cada T — no es solo re-evaluar RT·ln10 con")
    print(f"  un G_macro fijo. pKa₂ (B→C) no usa el macroestado y puede diferir")
    print(f"  en método ('born' vs 'direct') de pKa₀/pKa₁; ver diagnóstico de")
    print(f"  cross-validación de escala del protón más arriba si aplica.")
    print(f"  {sep}\n")


def print_distribution_table(
    df_dist: pd.DataFrame,
    pKa1: Optional[float],
    pKa2: Optional[float],
    ph_highlight: Optional[float] = None,
):
    """Imprime tabla de distribución de especies por pH."""
    sep = "─" * 74
    print(f"\n{sep}")
    print(f"  DISTRIBUCIÓN DE ESTADOS DE PROTONACIÓN  "
          f"(pKa₁ = {pKa1:.2f}" + (f", pKa₂ = {pKa2:.2f}" if pKa2 else "") + ")")
    print(sep)

    has_dianion = "f_dianion" in df_dist.columns

    # Cabecera
    if has_dianion:
        print(f"\n  {'pH':>6}  {'f(AH₂)':>10}  {'f(AH⁻)':>10}  {'f(A²⁻)':>10}  {'Dominante'}")
        print(f"  {'─'*6}  {'─'*10}  {'─'*10}  {'─'*10}  {'─'*12}")
    else:
        print(f"\n  {'pH':>6}  {'f(AH)':>10}  {'f(A⁻)':>10}  {'Dominante'}")
        print(f"  {'─'*6}  {'─'*10}  {'─'*10}  {'─'*12}")

    for _, row in df_dist.iterrows():
        pH   = row["pH"]
        fn   = row.get("f_neutral", 0.0)
        fa   = row.get("f_anion",   0.0)
        fd   = row.get("f_dianion", 0.0) if has_dianion else None

        fracs = {"neutral": fn, "anion": fa}
        if has_dianion:
            fracs["dianion"] = fd
        dominant = max(fracs, key=fracs.get)

        marker = " ◀" if (ph_highlight is not None and abs(pH - ph_highlight) < 0.05) else ""

        if has_dianion:
            print(f"  {pH:6.2f}  {fn:10.4f}  {fa:10.4f}  {fd:10.4f}  "
                  f"{dominant}{marker}")
        else:
            print(f"  {pH:6.2f}  {fn:10.4f}  {fa:10.4f}  {dominant}{marker}")

    # Puntos de cambio de especie dominante
    print(f"\n  Punto isoelectrónico (dominancia igual):")
    if pKa1:
        print(f"    AH₂ = AH⁻  →  pH = pKa₁ = {pKa1:.2f}")
    if pKa2:
        print(f"    AH⁻ = A²⁻  →  pH = pKa₂ = {pKa2:.2f}")
    if pKa1 and pKa2:
        pH_PI = (pKa1 + pKa2) / 2
        print(f"    Punto isoeléctrico (pI) ≈ (pKa₁ + pKa₂)/2 = {pH_PI:.2f}")

    print(sep)


def print_distribution_table_5(
    df_dist: pd.DataFrame,
    pKa0: Optional[float],
    pKa1: Optional[float],
    pKa2: Optional[float],
    ph_highlight: Optional[float] = None,
):
    """Imprime tabla de distribución de 5 especies (catión, R, D, anión, dianión)."""
    sep = "─" * 92
    print(f"\n{sep}")
    hdr = "  DISTRIBUCIÓN DE ESTADOS DE PROTONACIÓN  ("
    parts = []
    if pKa0 is not None: parts.append(f"pKa₀ = {pKa0:.2f}")
    if pKa1 is not None: parts.append(f"pKa₁ = {pKa1:.2f}")
    if pKa2 is not None: parts.append(f"pKa₂ = {pKa2:.2f}")
    print(hdr + ", ".join(parts) + ")")
    print(sep)

    has_cation  = "f_cation"  in df_dist.columns
    has_anion   = "f_anion"   in df_dist.columns
    has_dianion = "f_dianion" in df_dist.columns

    cols   = ["pH"]
    labels = [f"{'pH':>6}"]
    if has_cation:  labels.append(f"{'f(A⁺)':>9}")
    labels.append(f"{'f(R)':>9}")
    labels.append(f"{'f(D)':>9}")
    if has_anion:   labels.append(f"{'f(B⁻)':>9}")
    if has_dianion: labels.append(f"{'f(C²⁻)':>9}")
    labels.append("Dominante")
    print(f"\n  {'  '.join(labels)}")
    print(f"  {'─'*(len(''.join(labels)) + 2*(len(labels)-1))}")

    for _, row in df_dist.iterrows():
        pH = row["pH"]
        fracs = {"R": row.get("f_R", 0.0), "D": row.get("f_D", 0.0)}
        vals = [f"{pH:6.2f}"]
        if has_cation:
            fracs["A⁺"] = row.get("f_cation", 0.0)
            vals.append(f"{row.get('f_cation', 0.0):9.4f}")
        vals.append(f"{row.get('f_R', 0.0):9.4f}")
        vals.append(f"{row.get('f_D', 0.0):9.4f}")
        if has_anion:
            fracs["B⁻"] = row.get("f_anion", 0.0)
            vals.append(f"{row.get('f_anion', 0.0):9.4f}")
        if has_dianion:
            fracs["C²⁻"] = row.get("f_dianion", 0.0)
            vals.append(f"{row.get('f_dianion', 0.0):9.4f}")
        dominant = max(fracs, key=fracs.get)
        marker = " ◀" if (ph_highlight is not None and abs(pH - ph_highlight) < 0.05) else ""
        print(f"  {'  '.join(vals)}  {dominant}{marker}")

    print(f"\n  Nota: f(R) y f(D) son las fracciones DENTRO del macroestado neutro")
    print(f"  (equilibrio de tautómeros, independiente del pH); f(R)+f(D) = f(neutro/zwitterión total).")
    if pKa0 is not None and pKa1 is not None:
        pH_PI = (pKa0 + pKa1) / 2
        print(f"\n  pI aproximado ≈ (pKa₀ + pKa₁)/2 = {pH_PI:.2f}")
    print(sep)


# ─────────────────────────────────────────────
#  VISUALIZACIÓN (matplotlib)
# ─────────────────────────────────────────────

def plot_distribution(
    df_dist: pd.DataFrame,
    pKa1: Optional[float],
    pKa2: Optional[float],
    output_path: Optional[Path] = None,
    ph_highlight: Optional[float] = None,
    chain_mc: Optional[dict] = None,
):
    """
    Diagrama de distribución de especies (fracciones molares vs pH).
    Guarda como PNG si se provee output_path; si no, muestra en pantalla.

    chain_mc : dict de run_chain_montecarlo(), opcional (--montecarlo). En
        _run_classic() se construye con data_A=None, f_R=1.0, f_D=0.0, así
        que la columna 'f_neutral_macro' de df_dist_bands es equivalente a
        'f_neutral' acá (f_R=1 → f_R_macro == f_neutral_macro).
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: conda install matplotlib")
        return

    has_dianion = "f_dianion" in df_dist.columns
    pH_vals = df_dist["pH"].values

    df_bands = None
    if chain_mc is not None and chain_mc.get("ok"):
        cand = chain_mc.get("df_dist_bands")
        if cand is not None and not cand.empty and len(cand) == len(pH_vals):
            df_bands = cand
        elif cand is not None and not cand.empty:
            print("  ⚠  df_dist_bands tiene distinto largo que df_dist "
                  "(pH_all no coincide) — se omiten las bandas MC en el gráfico.")

    fig, ax = plt.subplots(figsize=(8, 5))

    colors = ["#2196F3", "#4CAF50", "#F44336"]  # azul, verde, rojo
    labels = ["AH₂ (q=0)", "AH⁻ (q=−1)", "A²⁻ (q=−2)"]
    cols   = ["f_neutral", "f_anion"] + (["f_dianion"] if has_dianion else [])
    # mapeo a columnas de df_dist_bands (run_chain_montecarlo siempre usa
    # nomenclatura de 5 estados internamente)
    band_cols = {"f_neutral": "f_neutral_macro", "f_anion": "f_anion", "f_dianion": "f_dianion"}

    for col, label, color in zip(cols, labels, colors):
        if col in df_dist.columns:
            bcol = band_cols.get(col)
            if df_bands is not None and bcol and f"{bcol}_lo" in df_bands.columns and f"{bcol}_hi" in df_bands.columns:
                ax.fill_between(pH_vals,
                                 df_bands[f"{bcol}_lo"].values * 100,
                                 df_bands[f"{bcol}_hi"].values * 100,
                                 color=color, alpha=0.15, linewidth=0, zorder=1)
            ax.plot(pH_vals, df_dist[col].values * 100, label=label,
                    color=color, linewidth=2.0, zorder=2)

    if df_bands is not None:
        ax.plot([], [], color="gray", alpha=0.3, linewidth=8,
                label="banda MC 95% (percentil 2.5–97.5)")

    # Líneas verticales en pKa
    if pKa1:
        ax.axvline(pKa1, color="#2196F3", linestyle="--", alpha=0.6,
                   label=f"pKa₁ = {pKa1:.2f}")
    if pKa2:
        ax.axvline(pKa2, color="#4CAF50", linestyle="--", alpha=0.6,
                   label=f"pKa₂ = {pKa2:.2f}")

    # pH oceánico antártico
    if ph_highlight:
        ax.axvline(ph_highlight, color="gray", linestyle=":", alpha=0.8,
                   label=f"pH = {ph_highlight} (referencia)")
        ax.axvspan(ph_highlight - 0.1, ph_highlight + 0.1,
                   alpha=0.08, color="gray")

    ax.set_xlabel("pH", fontsize=12)
    ax.set_ylabel("Fracción molar (%)", fontsize=12)
    ax.set_title("Distribución de estados de protonación\n"
                 "(Henderson-Hasselbalch, ciclo de Born)", fontsize=11)
    ax.set_ylim(-2, 102)
    ax.set_xlim(pH_vals[0], pH_vals[-1])
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g%%"))
    ax.legend(fontsize=9, loc="center right")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Gráfico guardado en: {output_path}")
    else:
        plt.show()
    plt.close(fig)


def plot_distribution_5(
    df_dist: pd.DataFrame,
    pKa0: Optional[float],
    pKa1: Optional[float],
    pKa2: Optional[float],
    output_path: Optional[Path] = None,
    ph_highlight: Optional[float] = None,
    chain_mc: Optional[dict] = None,
):
    """
    Diagrama de distribución de 5 especies (catión, R, D, anión, dianión) vs pH.

    chain_mc : dict de run_chain_montecarlo(), opcional (--montecarlo). Si se
        provee y trae df_dist_bands no vacío, se dibuja una banda sombreada
        de percentil 2.5–97.5 alrededor de cada curva (misma columna de
        pH que df_dist, ya que ambas se generan con el mismo pH_all — ver
        _run_extended()/_run_classic()).
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.ticker as ticker
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: conda install matplotlib")
        return

    pH_vals = df_dist["pH"].values

    df_bands = None
    if chain_mc is not None and chain_mc.get("ok"):
        cand = chain_mc.get("df_dist_bands")
        if cand is not None and not cand.empty and len(cand) == len(pH_vals):
            df_bands = cand
        elif cand is not None and not cand.empty:
            print("  ⚠  df_dist_bands tiene distinto largo que df_dist "
                  "(pH_all no coincide) — se omiten las bandas MC en el gráfico.")

    fig, ax = plt.subplots(figsize=(9, 5.5))

    series = []
    if "f_cation" in df_dist.columns:
        series.append(("f_cation", "A⁺ catión (q=+1)", "#9C27B0", "-"))
    series.append(("f_R", "R neutro (q=0)",       "#2196F3", "-"))
    series.append(("f_D", "D zwitterión (q=0)",   "#00BCD4", "--"))
    if "f_anion" in df_dist.columns:
        series.append(("f_anion", "B⁻ anión (q=−1)",   "#4CAF50", "-"))
    if "f_dianion" in df_dist.columns:
        series.append(("f_dianion", "C²⁻ dianión (q=−2)", "#F44336", "-"))

    for col, label, color, ls in series:
        if col in df_dist.columns:
            if df_bands is not None and f"{col}_lo" in df_bands.columns and f"{col}_hi" in df_bands.columns:
                ax.fill_between(pH_vals,
                                 df_bands[f"{col}_lo"].values * 100,
                                 df_bands[f"{col}_hi"].values * 100,
                                 color=color, alpha=0.15, linewidth=0, zorder=1)
            ax.plot(pH_vals, df_dist[col].values * 100, label=label,
                    color=color, linewidth=2.0, linestyle=ls, zorder=2)

    if df_bands is not None:
        ax.plot([], [], color="gray", alpha=0.3, linewidth=8,
                label="banda MC 95% (percentil 2.5–97.5)")

    if pKa0 is not None:
        ax.axvline(pKa0, color="#9C27B0", linestyle=":", alpha=0.6, label=f"pKa₀ = {pKa0:.2f}")
    if pKa1 is not None:
        ax.axvline(pKa1, color="#4CAF50", linestyle=":", alpha=0.6, label=f"pKa₁ = {pKa1:.2f}")
    if pKa2 is not None:
        ax.axvline(pKa2, color="#F44336", linestyle=":", alpha=0.6, label=f"pKa₂ = {pKa2:.2f}")

    if ph_highlight is not None:
        ax.axvline(ph_highlight, color="gray", linestyle=":", alpha=0.8,
                   label=f"pH = {ph_highlight} (referencia)")
        ax.axvspan(ph_highlight - 0.1, ph_highlight + 0.1, alpha=0.08, color="gray")

    ax.set_xlabel("pH", fontsize=12)
    ax.set_ylabel("Fracción molar (%)", fontsize=12)
    ax.set_title("Distribución de estados de protonación\n"
                 "catión ⇌ [neutro ⇌ zwitterión] ⇌ anión ⇌ dianión", fontsize=11)
    ax.set_ylim(-2, 102)
    ax.set_xlim(pH_vals[0], pH_vals[-1])
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g%%"))
    ax.legend(fontsize=8, loc="center right", ncol=1)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Gráfico guardado en: {output_path}")
    else:
        plt.show()
    plt.close(fig)


def plot_conformer_boltzmann(data_list: list[dict], output_path: Optional[Path] = None):
    """
    Gráfico de barras apiladas: distribución conformacional (298 K) por estado.
    Muestra los confórmeros con mayor peso en cada estado de protonación.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
    except ImportError:
        return

    fig, axes = plt.subplots(1, len(data_list), figsize=(5 * len(data_list), 4),
                             sharey=True)
    if len(data_list) == 1:
        axes = [axes]

    for ax, data in zip(axes, data_list):
        df = data["df"]
        if df.empty or "boltz_298K" not in df.columns:
            ax.set_title(data["mol_name"])
            ax.text(0.5, 0.5, "Sin datos", ha="center", va="center",
                    transform=ax.transAxes)
            continue

        # Filtrar confs con boltz > 1%
        col_b = "boltz_298K"
        sub = df[df[col_b] > 0.01].sort_values(col_b, ascending=False)

        colors = cm.viridis(np.linspace(0.2, 0.85, len(sub)))
        bars = ax.bar(range(len(sub)), sub[col_b].values * 100,
                      color=colors, width=0.7)
        ax.set_xticks(range(len(sub)))
        ax.set_xticklabels(sub["conformer"].tolist(), rotation=45, ha="right",
                           fontsize=7)
        ax.set_title(data["mol_name"], fontsize=9)
        ax.set_ylabel("Boltzmann @ 298 K (%)" if ax == axes[0] else "")
        ax.set_ylim(0, 105)

        # Anotar porcentaje sobre cada barra
        for bar, pct in zip(bars, sub[col_b].values * 100):
            if pct > 5:
                ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                        f"{pct:.1f}%", ha="center", fontsize=7)

    fig.suptitle("Distribución conformacional (Boltzmann, 298 K)", fontsize=11)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Gráfico conformacional guardado en: {output_path}")
    else:
        plt.show()
    plt.close(fig)


def plot_pka_histograms(chain_mc: Optional[dict], pka_labels: list,
                        output_path: Optional[Path] = None):
    """
    Histograma de las muestras Monte Carlo de cada pKa presente en la
    cadena (chain_mc['mc_pka0']/['mc_pka1']/['mc_pka2']), con el pKa
    determinístico (compute_pka(), línea vertical) y la media/CI95 de la
    distribución MC superpuestos. Complementa las bandas de confianza de
    plot_distribution_5()/plot_distribution() con la vista "por pKa" que
    pide la Fase 4 del plan (histograma de cada pKa individual).

    pka_labels : misma estructura que print_mc_summary() — lista de
        tuplas (key, base_result, label) con key en {'pKa0','pKa1','pKa2'}
        y base_result el resultado determinístico de compute_pka() (para
        la línea vertical de referencia).
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("  ⚠  matplotlib no disponible. Instalar con: conda install matplotlib")
        return

    if chain_mc is None or not chain_mc.get("ok"):
        return

    key_to_mc = {"pKa0": "mc_pka0", "pKa1": "mc_pka1", "pKa2": "mc_pka2"}
    panels = []
    for key, base_result, label in pka_labels:
        mc = chain_mc.get(key_to_mc.get(key))
        if mc is not None and mc.get("ok") and mc.get("pka_samples") is not None:
            panels.append((key, base_result, label, mc))

    if not panels:
        return

    fig, axes = plt.subplots(1, len(panels), figsize=(4.5 * len(panels), 4),
                             sharey=False)
    if len(panels) == 1:
        axes = [axes]

    for ax, (key, base_result, label, mc) in zip(axes, panels):
        samples = mc["pka_samples"]
        ax.hist(samples, bins=60, color="#2196F3", alpha=0.75, edgecolor="none")

        pka_det = base_result["pKa"] if (base_result and base_result.get("ok")) else None
        if pka_det is not None:
            ax.axvline(pka_det, color="black", linestyle="-", linewidth=1.5,
                       label=f"determinístico = {pka_det:.2f}")

        ax.axvline(mc["pka_mean"], color="#F44336", linestyle="--", linewidth=1.5,
                   label=f"media MC = {mc['pka_mean']:.2f}")
        lo95, hi95 = mc["pka_ci95"]
        ax.axvspan(lo95, hi95, color="#F44336", alpha=0.08,
                   label=f"CI95 = [{lo95:.2f}, {hi95:.2f}]")

        ax.set_title(label, fontsize=10)
        ax.set_xlabel("pKa")
        ax.set_ylabel("Frecuencia" if ax == axes[0] else "")
        ax.legend(fontsize=7, loc="upper right")
        ax.grid(True, alpha=0.25)

    fig.suptitle(f"Distribución Monte Carlo de pKa (N={chain_mc['n_samples']:,} muestras)",
                fontsize=11)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=150, bbox_inches="tight")
        print(f"  Histogramas de pKa (MC) guardados en: {output_path}")
    else:
        plt.show()
    plt.close(fig)


# ─────────────────────────────────────────────
#  EXPORTACIÓN
# ─────────────────────────────────────────────

def export_results(
    pka_list: list,
    df_dist: pd.DataFrame,
    data_list: list[dict],
    fmt: str,
    taut_result: Optional[dict] = None,
    chain_mc: Optional[dict] = None,
):
    """
    Exporta resultados del ciclo de pKa y distribución de especies.

    pka_list: lista de tuplas (resultado_dict_o_None, etiqueta_str), p. ej.
      [(pka0, "pKa0"), (pka1, "pKa1"), (pka2, "pKa2")]
    taut_result: dict de compute_tautomer_equilibrium(), opcional (se agrega
      como fila/hoja adicional si se provee).
    chain_mc: dict de run_chain_montecarlo(), opcional (--montecarlo). Si se
      provee y ok=True, agrega las muestras crudas de pKa0/pKa1/pKa2
      (mc_samples.csv / hoja MC_samples) y las bandas de confianza de la
      distribución de especies (mc_distribution_bands.csv / hoja
      MC_distribution_bands) — útil para reproducir el análisis de
      incertidumbre sin recorrer todo el pipeline.
    """
    out_dir = PC_WORK_ROOT
    out_dir.mkdir(parents=True, exist_ok=True)

    # ── DataFrame del ciclo de Born / directo ────────────────────────────
    cycle_rows = []
    for pka_r, label in pka_list:
        if pka_r is None:
            continue
        cycle_rows.append({
            "reaction":          label,
            "mol_acid":          pka_r.get("mol_acid"),
            "mol_base":          pka_r.get("mol_base"),
            "method":            pka_r.get("method"),
            "aq_source_acid":    pka_r.get("aq_source_acid", "CPCM"),
            "aq_source_base":    pka_r.get("aq_source_base", "CPCM"),
            "conf_acid":         pka_r.get("conf_acid"),
            "conf_base":         pka_r.get("conf_base"),
            "T_K":               pka_r.get("T"),
            "dG_gas_kcal":       pka_r.get("dG_gas_kcal"),
            "dG_solv_acid_kcal": pka_r.get("dG_solv_acid_kcal"),
            "dG_solv_base_kcal": pka_r.get("dG_solv_base_kcal"),
            "dG_solv_Hplus_kcal": DG_SOLV_HPLUS if pka_r.get("method") == "born" else None,
            "G_hplus_aq_kcal":   pka_r.get("G_hplus_aq_kcal"),
            "dG_aq_kcal":        pka_r.get("dG_aq_kcal"),
            "dG_aq_kJ":          (pka_r.get("dG_aq_kcal") or 0) * 4.184,
            "pKa":               pka_r.get("pKa"),
            "ok":                pka_r.get("ok"),
            # Corrección empírica (--empirical-correction), si se aplicó a
            # este pKa; columnas quedan en None/NaN si no se pidió corrección.
            "pKa_raw":           pka_r.get("pKa_raw"),
            "pKa_corrected":     pka_r.get("pKa_corrected"),
            "correction_type":   (pka_r.get("correction") or {}).get("group_type"),
            "correction_mode":   (pka_r.get("correction") or {}).get("mode"),
            "correction_n_calib": (pka_r.get("correction") or {}).get("n"),
            "correction_slope":  (pka_r.get("correction") or {}).get("slope"),
            "correction_intercept": (pka_r.get("correction") or {}).get("intercept"),
        })

    df_cycle = pd.DataFrame(cycle_rows)

    # ── DataFrame del equilibrio de tautómeros (opcional) ────────────────
    df_taut = None
    if taut_result is not None and taut_result.get("ok"):
        df_taut = pd.DataFrame([{
            "mol_R":         taut_result.get("mol_R"),
            "mol_D":         taut_result.get("mol_D"),
            "aq_source":     taut_result.get("aq_source", "CPCM"),
            "conf_R":        taut_result.get("conf_R"),
            "conf_D":        taut_result.get("conf_D"),
            "T_K":           taut_result.get("T"),
            "G_R_Eh":        taut_result.get("G_R_Eh"),
            "G_D_Eh":        taut_result.get("G_D_Eh"),
            "dG_taut_kcal":  taut_result.get("dG_taut_kcal"),
            "K_taut":        taut_result.get("K_taut"),
            "f_R":           taut_result.get("f_R"),
            "f_D":           taut_result.get("f_D"),
            "G_macro_Eh":    taut_result.get("G_macro_Eh"),
        }])

    # ── DataFrame de distribución conformacional ─────────────────────────
    conf_rows = []
    for d in data_list:
        df_c = d["df"]
        for _, row in df_c.iterrows():
            r = {"mol": d["mol_name"], "aq_source": d.get("aq_source", "CPCM")}
            r.update(row.to_dict())
            conf_rows.append(r)
    df_confs = pd.DataFrame(conf_rows)

    # ── DataFrame de muestras crudas Monte Carlo (opcional) ───────────────
    df_mc_samples = None
    df_mc_bands = None
    if chain_mc is not None and chain_mc.get("ok"):
        mc_cols = {}
        for key, colname in [("mc_pka0", "pKa0"), ("mc_pka1", "pKa1"), ("mc_pka2", "pKa2")]:
            mc = chain_mc.get(key)
            if mc is not None and mc.get("ok"):
                mc_cols[colname] = mc["pka_samples"]
        if mc_cols:
            df_mc_samples = pd.DataFrame(mc_cols)
        if chain_mc.get("df_dist_bands") is not None and not chain_mc["df_dist_bands"].empty:
            df_mc_bands = chain_mc["df_dist_bands"]

    if fmt in ("csv", "both"):
        df_cycle.to_csv(out_dir / "pka_cycle.csv", index=False)
        df_dist.to_csv(out_dir / "protonation_distribution.csv", index=False)
        df_confs.to_csv(out_dir / "conformer_boltzmann.csv", index=False)
        if df_taut is not None:
            df_taut.to_csv(out_dir / "tautomer_equilibrium.csv", index=False)
        if df_mc_samples is not None:
            df_mc_samples.to_csv(out_dir / "mc_samples.csv", index=False)
        if df_mc_bands is not None:
            df_mc_bands.to_csv(out_dir / "mc_distribution_bands.csv", index=False)
        print(f"  CSVs guardados en: {out_dir}/")

    if fmt in ("xlsx", "both"):
        xlsx_path = out_dir / "protonation_cycle_results.xlsx"
        with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
            if not df_cycle.empty:
                df_cycle.to_excel(writer, sheet_name="pKa_cycle", index=False)
            if df_taut is not None:
                df_taut.to_excel(writer, sheet_name="Tautomer_equilibrium", index=False)
            df_dist.to_excel(writer, sheet_name="Species_distribution", index=False)
            if not df_confs.empty:
                # columnas relevantes
                keep = ["mol", "conformer", "E_sp_aq", "thermo_ok", "n_imag",
                        "G_total_298K", "boltz_298K", "dG_298K_kcal"]
                keep = [c for c in keep if c in df_confs.columns]
                df_confs[keep].to_excel(writer, sheet_name="Conformer_Boltzmann", index=False)
            if df_mc_samples is not None:
                df_mc_samples.to_excel(writer, sheet_name="MC_samples", index=False)
            if df_mc_bands is not None:
                df_mc_bands.to_excel(writer, sheet_name="MC_distribution_bands", index=False)
        print(f"  Excel guardado en: {xlsx_path}")


# ─────────────────────────────────────────────
#  PUNTO DE ENTRADA
# ─────────────────────────────────────────────

class _Tee:
    """
    Duplica cada escritura hacia un stream original (stdout/stderr de
    consola) y, además, hacia un archivo .log en disco. Se usa flush()
    inmediato en cada write() para que el .log quede completo incluso si
    el proceso muere a medias (ORCA se cae, Ctrl+C, etc.) en vez de
    perderse en el buffer.
    """

    def __init__(self, primary, logfile):
        self._primary = primary
        self._logfile = logfile

    def write(self, data):
        self._primary.write(data)
        self._logfile.write(data)
        self._logfile.flush()
        return len(data)

    def flush(self):
        self._primary.flush()
        self._logfile.flush()

    def isatty(self):
        return self._primary.isatty()


def _setup_run_logger(work_root: Path, argv: list[str]):
    """
    Prepara un .log COMPLETO y EXCLUSIVO para esta corrida.

    Ubicación: <carpeta_de_mol>/logs/, es decir un nivel arriba de runs/
    (la carpeta de grupo determinada por --mol / --work-root), para que
    los logs de todas las corridas sobre una misma molécula queden
    juntos y separados de los de otra molécula.

    Nombre: run_<timestamp con microsegundos>_pid<PID>.log — así, sin
    importar cuántas veces o en qué modos/configuraciones se corra el
    script (incluso si se lanzan dos corridas casi al mismo tiempo),
    cada invocación obtiene su propio archivo y ninguna sobreescribe a
    otra.

    Redirige sys.stdout y sys.stderr a través de _Tee, de modo que TODO
    lo que ya se imprime en el script (incluyendo el detalle verbose del
    ciclo de Born, que ahora es el comportamiento por default) queda
    también volcado al .log, sin tener que tocar cada print() del
    código.

    Devuelve (log_path, log_file, orig_stdout, orig_stderr) para que el
    llamador pueda cerrar/restaurar todo correctamente en un finally.
    """
    group_dir = work_root.parent  # carpeta de la molécula (p. ej. <cwd>/MGLYVAL/)
    log_dir = group_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)

    ts = _dt_datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    log_path = log_dir / f"run_{ts}_pid{os.getpid()}.log"

    log_file = open(log_path, "w", encoding="utf-8")
    log_file.write("═" * 74 + "\n")
    log_file.write("  protonation_cycle.py — log completo de corrida\n")
    log_file.write(f"  fecha    : {_dt_datetime.now().isoformat(timespec='seconds')}\n")
    log_file.write(f"  cwd      : {Path.cwd()}\n")
    log_file.write(f"  work_root: {work_root}\n")
    log_file.write(f"  comando  : {' '.join(argv)}\n")
    log_file.write("═" * 74 + "\n\n")
    log_file.flush()

    orig_stdout, orig_stderr = sys.stdout, sys.stderr
    sys.stdout = _Tee(orig_stdout, log_file)
    sys.stderr = _Tee(orig_stderr, log_file)

    return log_path, log_file, orig_stdout, orig_stderr


def _main_pka():
    # Debe ir ANTES de cualquier lectura de estas variables en esta función
    # (p. ej. como default= en add_argument más abajo). Python exige que la
    # declaración `global` preceda todo uso del nombre dentro del mismo
    # scope, o produce SyntaxError: "name used prior to global declaration".
    global PC_WORK_ROOT, GV_CUTOFF, GV_QS, GV_QH

    parser = argparse.ArgumentParser(
        description="Ciclo termodinámico de pKa y distribución de Boltzmann "
                    "entre estados de protonación (Born cycle, Ho & Coote 2010).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Ejemplos:
  # Los tres estados clásicos (pKa1 + pKa2 + distribución):
  python protonation_cycle.py \\
      --mol MGLYVAL --neutral R --anion B --dianion C

  # Solo pKa1:
  python protonation_cycle.py --mol MGLYVAL --neutral R --anion B

  # Cadena completa de 5 estados con catión y zwitterión (tautómero de R):
  python protonation_cycle.py \\
      --mol MGLYVAL \\
      --neutral R --anion B --dianion C --protonated A --zwitterion D \\
      --ph 7 --method direct --export xlsx --plot --verbose

  # Con exportación, gráficos y pH de referencia antártico:
  python protonation_cycle.py \\
      --mol SHINORINE --neutral R --anion B --dianion C \\
      --ph 8.1 --export xlsx --plot --verbose
        """,
    )

    parser.add_argument("--neutral",    metavar="DIR",
                        help="Carpeta en runs/ para el neutro R (q=0), p. ej. 'R'")
    parser.add_argument("--anion",      metavar="DIR",
                        help="Carpeta en runs/ para el anión B (q=−1), p. ej. 'B'")
    parser.add_argument("--dianion",    metavar="DIR",
                        help="Carpeta en runs/ para el dianión C (q=−2), p. ej. 'C'")
    parser.add_argument("--protonated", metavar="DIR",
                        help="Carpeta en runs/ para el catión protonado A (q=+1), "
                             "p. ej. 'A'. OPCIONAL: solo tiene efecto junto con "
                             "--zwitterion. Si se omite, no se calcula pKa₀ y la "
                             "cadena arranca en el macroestado neutro/zwitterión.")
    parser.add_argument("--zwitterion", metavar="DIR",
                        help="Carpeta en runs/ para el zwitterión D — TAUTÓMERO de "
                             "--neutral (mismo q=0, mismos átomos), p. ej. 'D'. "
                             "Activa la cadena extendida; requiere también "
                             "--neutral. --protonated es opcional dentro de esta "
                             "cadena.")
    parser.add_argument("--mol", metavar="MOL", default=None,
                        help="Nombre de la carpeta de la molécula dentro del "
                             "directorio actual, p. ej. 'MGLYVAL' o 'SHINORINE'. "
                             "Se usa para construir la ruta por default: "
                             "<cwd>/<MOL>/runs. Ignorado si se pasa --work-root.")
    parser.add_argument("--work-root", metavar="DIR", default=None,
                        help="Ruta explícita a la carpeta runs/ que contiene los "
                             "estados (default: <cwd>/<MOL>/runs usando --mol, "
                             "p. ej. ~/PAPER1/MGLYVAL/runs). Si se pasa, tiene "
                             "prioridad sobre --mol.")
    parser.add_argument("--ph",       type=float, default=PH_ANTARCTIC,
                        help=f"pH de referencia para la distribución "
                             f"(default: {PH_ANTARCTIC}, agua de mar antártica)")
    parser.add_argument("--ph-range", nargs=2, type=float, default=[0.0, 14.0],
                        metavar=("pH_MIN", "pH_MAX"),
                        help="Rango de pH para el diagrama (default: 0 14)")
    parser.add_argument("--ph-step",  type=float, default=0.1,
                        help="Paso del eje de pH (default: 0.1)")
    parser.add_argument("--method",   choices=["auto", "born", "direct"], default="auto",
                        help="Método para ΔG_aq: 'born' (ciclo completo, requiere "
                             "SP fase gas), 'direct' (solo energía acuosa, sin SP "
                             "fase gas), 'auto' usa born si hay SP gas y si no cae "
                             "a direct (default: auto). Los pasos que involucran el "
                             "macroestado neutro/zwitterión (con --protonated/"
                             "--zwitterion) siempre usan 'direct'.")
    parser.add_argument("--spsmd", action="store_true",
                        help="Usa las energías SP con SMD(water) "
                             "(04c_sp_smd/conf_NNN/sp_smd.out, calculadas con "
                             "pipeline.py --spsmd) como energía SP acuosa en vez de "
                             "las CPCM(water) normales (04_sp/conf_NNN/sp.out) para "
                             "TODOS los estados de esta corrida (neutro, anión, "
                             "dianión, catión, zwitterión). No mezcla — si un "
                             "confórmero no tiene SP_SMD calculado queda con "
                             "E_sp_aq=NaN y se excluye del ensemble/best_conf "
                             "(mismo comportamiento que un SP fase gas faltante). "
                             "Con method='born' se emite una advertencia: las "
                             "constantes del ciclo de Born fueron calibradas contra "
                             "CPCM, no SMD; se recomienda combinar con "
                             "--method direct.")
    parser.add_argument("--temp",     type=float, default=T_STD, choices=TEMPS_K,
                        help=f"Temperatura T (K) para el pKa reportado (default: "
                             f"{T_STD} K, estándar). Restringido a {TEMPS_K} porque "
                             f"G_boltz (usado en los pasos con macroestado "
                             f"neutro/zwitterión, --protonated/--zwitterion) solo "
                             f"está tabulado a esas 4 T; pedir otra T degrada "
                             f"silenciosamente esos pasos a un solo confórmero en "
                             f"vez del ensemble Boltzmann. Útil para reportar el "
                             f"pKa a T de agua de mar antártica (274.15 K) en vez "
                             f"del estándar de 298.15 K. Ver también el escaneo de "
                             f"las 4 T que se imprime siempre como diagnóstico.")
    parser.add_argument("--export",   choices=["csv", "xlsx", "both"], default=None,
                        help="Exportar resultados")
    parser.add_argument("--plot",     action="store_true",
                        help="Generar gráficos de distribución y Boltzmann conformacional")
    parser.add_argument("--verbose",  action=argparse.BooleanOptionalAction, default=True,
                        help="Mostrar contribuciones detalladas del ciclo de Born "
                             "(default: activado). Usa --no-verbose para el modo "
                             "resumido. El .log de la corrida (ver --mol) siempre "
                             "guarda TODO lo impreso, independientemente de este flag.")
    parser.add_argument("--goodvibes", action="store_true",
                        help="Corre el ciclo COMPLETO dos veces: una con G_corr "
                             "armónico estándar de NumFreq (RRHO, tal cual lo "
                             "reporta ORCA) y otra con la corrección cuasi-armónica "
                             "de GoodVibes (mRRHO de Grimme para los modos de baja "
                             "frecuencia), y al final compara cuánto cambian los "
                             "pKa's y las fracciones de especie. Requiere el "
                             "binario 'goodvibes' instalado "
                             "(pip install goodvibes --break-system-packages).")
    parser.add_argument("--gv-cutoff", type=float, default=GV_CUTOFF, metavar="CM-1",
                        help=f"Frequency cutoff de GoodVibes en cm⁻¹ "
                             f"(default: {GV_CUTOFF:.0f})")
    parser.add_argument("--gv-qs", choices=["grimme", "truhlar"], default=GV_QS,
                        help=f"Método de interpolación de entropía cuasi-armónica "
                             f"de GoodVibes (default: {GV_QS})")
    parser.add_argument("--gv-qh", action="store_true",
                        help="Además de la entropía, aplica también la corrección "
                             "cuasi-armónica de Head-Gordon a la entalpía "
                             "(default: solo la entropía, la más relevante)")
    parser.add_argument("--montecarlo", action="store_true",
                        help="Propaga incertidumbre por Monte Carlo a través de "
                             "toda la cadena de pKa's (ver --mc-*), preservando "
                             "la correlación física de G(H⁺,aq) entre pKa0/pKa1/"
                             "pKa2. Imprime media±std/CI95 junto al valor "
                             "determinístico y, si además se pasa --export, "
                             "agrega las muestras crudas y las bandas de "
                             "confianza de la distribución de especies.")
    parser.add_argument("--mc-samples", type=int, default=10_000, metavar="N",
                        help="Número de muestras Monte Carlo (default: 10000)")
    parser.add_argument("--mc-sigma-elec", type=float, default=0.5, metavar="KCAL",
                        help="σ asumida para la energía electrónica SP por "
                             "estado, en kcal/mol (default: 0.5, típico para "
                             "ωB97X-D3/def2-TZVPD)")
    parser.add_argument("--mc-sigma-gcorr", type=float, default=0.4, metavar="KCAL",
                        help="σ asumida para la corrección térmica G_corr por "
                             "estado, en kcal/mol (default: 0.4; mayor si hay "
                             "n_imag>0 o gcorr_source='harmonic' en vez de mRRHO)")
    parser.add_argument("--mc-sigma-hplus", type=float, default=2.5, metavar="KCAL",
                        help="σ asumida para G(H⁺,aq) — dispersión típica entre "
                             "valores de literatura (Pliego-Riveros −264.61, "
                             "Kelly-Cramer-Truhlar −265.9, Zhan-Dixon ≈−262.4; "
                             "default: 2.5). Esta es sistemática y se comparte "
                             "entre los tres pKa's de la cadena, no se muestrea "
                             "independiente por tramo.")
    parser.add_argument("--mc-seed", type=int, default=None, metavar="N",
                        help="Semilla para reproducibilidad del Monte Carlo "
                             "(default: sin fijar)")
    parser.add_argument("--pair", nargs=2, metavar=("ACIDO", "BASE"), default=None,
                        help="Modo par: calcula UN pKa entre dos carpetas de "
                             "runs/ (ácido → base + H⁺), sin cadena ni "
                             "macroestado. Pensado para AP: --pair APA APR. "
                             "Incompatible con --neutral/--anion/--dianion/"
                             "--protonated/--zwitterion. La corrección se pide "
                             "con --empirical-correction pair=pyridinium.")
    parser.add_argument("--calibration-file", metavar="JSON", default=None,
                        help="calibration.json generado por `maaqc.py calibrate`: "
                             "reemplaza, por tipo, la tabla EMPIRICAL_CALIBRATION "
                             "embebida. Avisa si --method o --temp difieren de los "
                             "usados al calibrar.")
    parser.add_argument("--empirical-correction", metavar="pKaN=TIPO[,pKaN=TIPO...]",
                        default=None,
help="Aplica la corrección empírica calibrada (ver "
                             "EMPIRICAL_CALIBRATION) a los slots indicados. "
                             "Formato: 'pKa0=carboxyl,pKa1=carboxyl_2,"
                             "pKa2=pyridinium' (con --pair: 'pair=pyridinium'). "
                             "Tipos: 'carboxyl' [1er COOH, +1→0], 'carboxyl_2' "
                             "[2º COOH, 0→−1], 'amino' [NH₃⁺ sin carga vecina, "
                             "0→−1], 'amino_dianion' [NH₃⁺ con 2×COO⁻, −1→−2], "
                             "'pyridinium' [N–H⁺ aromático, +1→0, AP]. El pKa "
                             "corregido reemplaza al crudo en la distribución, el "
                             "pI y los gráficos; el crudo se sigue imprimiendo y "
                             "exportando. Sin este flag no hay corrección.")
    parser.add_argument("--empirical-correction-mode",
                        choices=["auto", "linear", "offset"], default="auto",
                        help="Modo de ajuste de --empirical-correction (ver "
                             "fit_empirical_correction()). 'auto' (default): "
                             "regresión lineal si hay ≥3 puntos de "
                             "calibración para ese tipo, desplazamiento puro "
                             "si hay exactamente 2 (recomendado — con 2 "
                             "puntos la pendiente de una regresión no está "
                             "verificada).")
    args = parser.parse_args()

    if args.calibration_file:
        _meta = load_calibration_file(Path(args.calibration_file).expanduser())
        _cm, _ct = _meta.get("method_resolved"), _meta.get("T")
        if _cm and args.method not in ("auto", _cm):
            print(f"AVISO: la calibración se hizo con method='{_cm}' y estás "
                  f"corriendo --method {args.method}: usan referencias de "
                  f"G(H⁺,aq) distintas; la corrección no es transferible.",
                  file=sys.stderr)
        if _ct is not None and abs(float(_ct) - float(args.temp)) > 1e-6:
            print(f"AVISO: la calibración es a T={_ct} K y corres a "
                  f"T={args.temp} K; la corrección empírica se ajustó solo a "
                  f"la T de calibración.", file=sys.stderr)

    args.empirical_correction_spec = _parse_empirical_correction_spec(
        args.empirical_correction)

    GV_CUTOFF = args.gv_cutoff
    GV_QS     = args.gv_qs
    GV_QH     = args.gv_qh
    if args.work_root:
        PC_WORK_ROOT = Path(args.work_root).expanduser().resolve()
    elif args.mol:
        PC_WORK_ROOT = (Path.cwd() / args.mol / "runs").resolve()
    else:
        parser.error(
            "Debes indicar --mol MOL (p. ej. --mol MGLYVAL) o, alternativamente, "
            "--work-root /ruta/a/runs. Ya no se asume ./runs en el directorio actual, "
            "porque ahora conviven varias moléculas (MGLYVAL, SHINORINE, etc.) bajo "
            "el mismo PAPER1/."
        )
    if not PC_WORK_ROOT.exists():
        parser.error(
            f"No existe la carpeta runs/: {PC_WORK_ROOT}\n"
            f"Verifica que la molécula '{args.mol}' exista bajo {Path.cwd()} "
            f"(se esperaba {Path.cwd() / (args.mol or '<MOL>') / 'runs'}), "
            f"o usa --work-root /ruta/a/runs para indicar la correcta."
        )

    # Modo par: excluyente con la cadena
    if args.pair:
        _chain_flags = [n for n, v in (("--neutral", args.neutral), ("--anion", args.anion),
                                       ("--dianion", args.dianion),
                                       ("--protonated", args.protonated),
                                       ("--zwitterion", args.zwitterion)) if v is not None]
        if _chain_flags:
            parser.error("--pair es incompatible con " + ", ".join(_chain_flags))

    # Verificar que se pasó al menos un par de estados
    if not args.pair and args.neutral is None and args.anion is None:
        parser.error("Debes proporcionar al menos --neutral y --anion (o usar --pair).")

    # El catión protonado A⁺ (--protonated) es OPCIONAL dentro de la cadena
    # extendida: si no se pasa, simplemente se omite pKa₀ y la distribución
    # arranca en el macroestado neutro/zwitterión. Lo único que activa el
    # modo extendido (con tautómero R⇌D) es --zwitterion.
    use_extended = args.zwitterion is not None
    if use_extended:
        missing = [name for name, val in
                   [("--zwitterion", args.zwitterion), ("--neutral", args.neutral)]
                   if val is None]
        if missing:
            parser.error(
                "La cadena extendida (zwitterión) requiere --zwitterion y "
                "--neutral simultáneamente. Falta(n): " + ", ".join(missing) + ". "
                "--protonated (catión A⁺, q=+1) es opcional: si se omite, "
                "simplemente no se calcula pKa₀ y la cadena arranca en el "
                "macroestado neutro/zwitterión."
            )
    elif args.protonated is not None and not args.pair:
        parser.error(
            "--protonated (catión A⁺) solo tiene sentido junto con "
            "--zwitterion (cadena extendida); para la cadena clásica de 3 "
            "estados no se usa. Agrega --zwitterion o quita --protonated."
        )

    log_path, log_file, orig_stdout, orig_stderr = _setup_run_logger(PC_WORK_ROOT, sys.argv)

    try:
        print("\n" + "═" * 74)
        print("  PROTONATION CYCLE  —  FAE-MAA / kmilozeta")
        print("═" * 74)
        if args.mol:
            print(f"  molécula : {args.mol}")
        print(f"  runs/ : {PC_WORK_ROOT}")
        print(f"  log   : {log_path}")

        if args.pair:
            _run_pair(args, "harmonic")
            if args.goodvibes:
                print("\n" + "═" * 74 + "\n  SEGUNDA PASADA: gcorr_source=goodvibes\n" + "═" * 74)
                _run_pair(args, "goodvibes")
        elif args.goodvibes:
            _run_comparison(args, use_extended)
        elif use_extended:
            _run_extended(args)
        else:
            _run_classic(args)
    finally:
        # Restaurar stdout/stderr originales y cerrar el .log de esta
        # corrida ANTES de salir, sin importar si terminó bien o con error
        # (así una excepción o un Ctrl+C también quedan registrados).
        sys.stdout, sys.stderr = orig_stdout, orig_stderr
        print(f"\n  (log completo de esta corrida guardado en: {log_path})\n")
        log_file.close()


def _run_pair(args, gcorr_source: str = "harmonic"):
    """
    Modo par: UN solo pKa HA → A⁻ + H⁺ entre dos carpetas (p. ej. AP:
    --pair APA APR). No construye cadena ni macroestado. La corrección
    empírica se pide con --empirical-correction pair=<tipo>.
    """
    acid, base = args.pair
    print(f"\n▸ Cargando ácido: {acid}")
    data_a = collect_mol_data(acid, verbose=args.verbose,
                              gcorr_source=gcorr_source, use_smd=args.spsmd)
    print_mol_summary(data_a, verbose=args.verbose)
    print(f"\n▸ Cargando base: {base}")
    data_b = collect_mol_data(base, verbose=args.verbose,
                              gcorr_source=gcorr_source, use_smd=args.spsmd)
    print_mol_summary(data_b, verbose=args.verbose)

    print(f"\n▸ Calculando pKa  ({acid} → {base} + H⁺)  @ T={args.temp:.2f} K")
    res = compute_pka(data_a, data_b, T=args.temp, verbose=args.verbose,
                      method=args.method, gcorr_source=gcorr_source)
    label = f"pKa ({acid}→{base})"
    print_pka_summary([(res, label)])

    ec_spec = getattr(args, "empirical_correction_spec", {}) or {}
    if "pair" in ec_spec:
        raw = res["pKa"] if (res and res["ok"]) else None
        corr = apply_empirical_correction(raw, ec_spec["pair"],
                                          gcorr_source=gcorr_source,
                                          mode=getattr(args, "empirical_correction_mode", "auto"))
        _print_empirical_correction(corr, label)
    return {"pka": res}


def _run_classic(args, gcorr_source: str = "harmonic"):
    """
    Cadena original de 3 estados: neutro → anión → dianión.

    gcorr_source controla de dónde sale G_corr para TODOS los estados:
      "harmonic"  → RRHO puro de NumFreq (default, comportamiento original)
      "goodvibes" → mRRHO cuasi-armónico de GoodVibes (ver --goodvibes)

    Retorna un dict con los resultados clave (pKa1, pKa2, df_dist, ...) para
    que _run_comparison() pueda diferenciar entre las dos pasadas.
    """
    data_neutral = data_anion = data_dianion = None

    if args.neutral:
        print(f"\n▸ Cargando estado neutro: {args.neutral}")
        data_neutral = collect_mol_data(args.neutral, verbose=args.verbose,
                                        gcorr_source=gcorr_source, use_smd=args.spsmd)
        print_mol_summary(data_neutral, verbose=args.verbose)

    if args.anion:
        print(f"\n▸ Cargando anión: {args.anion}")
        data_anion = collect_mol_data(args.anion, verbose=args.verbose,
                                      gcorr_source=gcorr_source, use_smd=args.spsmd)
        print_mol_summary(data_anion, verbose=args.verbose)

    if args.dianion:
        print(f"\n▸ Cargando dianión: {args.dianion}")
        data_dianion = collect_mol_data(args.dianion, verbose=args.verbose,
                                        gcorr_source=gcorr_source, use_smd=args.spsmd)
        print_mol_summary(data_dianion, verbose=args.verbose)

    pka1_result = pka2_result = None

    if data_neutral and data_anion:
        print(f"\n▸ Calculando pKa₁  ({args.neutral} → {args.anion} + H⁺)  @ T={args.temp:.2f} K")
        pka1_result = compute_pka(data_neutral, data_anion,
                                  T=args.temp, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if data_anion and data_dianion:
        print(f"\n▸ Calculando pKa₂  ({args.anion} → {args.dianion} + H⁺)  @ T={args.temp:.2f} K")
        pka2_result = compute_pka(data_anion, data_dianion,
                                  T=args.temp, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if pka1_result is None and pka2_result is None:
        print("\n  ⚠  No se pudo calcular ningún pKa. "
              "Revisa los directorios y los SP fase gas.")
        sys.exit(1)

    temp_scan = _pka_temperature_scan(
        [("pKa1", data_neutral, data_anion), ("pKa2", data_anion, data_dianion)],
        gcorr_source=gcorr_source, method=args.method,
    )
    _print_pka_temperature_scan(temp_scan, reported_T=args.temp)

    print_pka_summary([(pka1_result, "pKa₁"), (pka2_result, "pKa₂")])

    pKa1_raw = pka1_result["pKa"] if (pka1_result and pka1_result["ok"]) else None
    pKa2_raw = pka2_result["pKa"] if (pka2_result and pka2_result["ok"]) else None

    # ── Corrección empírica (--empirical-correction) ───────────────────────
    ec_spec = getattr(args, "empirical_correction_spec", {}) or {}
    ec_mode = getattr(args, "empirical_correction_mode", "auto")
    corr1 = corr2 = None
    if ec_spec:
        if "pKa1" in ec_spec:
            corr1 = apply_empirical_correction(pKa1_raw, ec_spec["pKa1"],
                                               gcorr_source=gcorr_source, mode=ec_mode)
        if "pKa2" in ec_spec:
            corr2 = apply_empirical_correction(pKa2_raw, ec_spec["pKa2"],
                                               gcorr_source=gcorr_source, mode=ec_mode)
        _print_empirical_correction(corr1, "pKa₁")
        _print_empirical_correction(corr2, "pKa₂")
        if any(c and c.get("ok") for c in (corr1, corr2)):
            print(f"\n  Nota: la distribución de especies, el pI y los "
                  f"gráficos de aquí en adelante usan los pKa's CORREGIDOS "
                  f"EMPÍRICAMENTE donde se pidió corrección (ver arriba). "
                  f"El escaneo de temperatura (arriba) sigue mostrando "
                  f"valores CRUDOS Born/CPCM sin corregir: la calibración "
                  f"empírica solo se ajustó a T=298.15 K. Si además usas "
                  f"--montecarlo, esa propagación de incertidumbre corre "
                  f"sobre el ciclo de Born crudo (no conoce la corrección "
                  f"empírica); súmale el RMSE de calibración reportado "
                  f"arriba como una fuente de error sistemático aparte.")
        for c, pka_result in ((corr1, pka1_result), (corr2, pka2_result)):
            if c and c.get("ok") and pka_result is not None:
                pka_result["pKa_raw"] = c["pKa_raw"]
                pka_result["pKa_corrected"] = c["pKa_corrected"]
                pka_result["correction"] = c

    pKa1 = corr1["pKa_corrected"] if (corr1 and corr1.get("ok")) else pKa1_raw
    pKa2 = corr2["pKa_corrected"] if (corr2 and corr2.get("ok")) else pKa2_raw

    if pKa1 is not None:
        pH_range = np.arange(args.ph_range[0], args.ph_range[1] + args.ph_step / 2,
                             args.ph_step)
        pH_highlights = [args.ph, PH_STANDARD]
        if pKa1: pH_highlights.append(round(pKa1, 1))
        if pKa2: pH_highlights.append(round(pKa2, 1))
        pH_all = sorted(set(list(pH_range) + pH_highlights))
        pH_all = [p for p in pH_all if args.ph_range[0] <= p <= args.ph_range[1]]

        df_dist = boltzmann_protonation_states(pKa1, pKa2, pH_all)

        pH_display = np.arange(args.ph_range[0], args.ph_range[1] + 0.5, 0.5)
        pH_display_list = sorted(set(
            [round(p, 1) for p in pH_display] + pH_highlights
        ))
        pH_display_list = [p for p in pH_display_list
                           if args.ph_range[0] <= p <= args.ph_range[1]]
        df_display = boltzmann_protonation_states(pKa1, pKa2, pH_display_list)
        print_distribution_table(df_display, pKa1, pKa2, ph_highlight=args.ph)

        print(f"\n  ── Fracciones a pH de referencia ──")
        for pH_ref, label in [
            (args.ph,     f"pH {args.ph} (agua de mar antártica)"),
            (PH_STANDARD, f"pH {PH_STANDARD} (estándar)"),
        ]:
            row = boltzmann_protonation_states(pKa1, pKa2, [pH_ref])
            if not row.empty:
                fn = row["f_neutral"].iloc[0] * 100 if "f_neutral" in row else None
                fa = row["f_anion"].iloc[0]   * 100 if "f_anion"   in row else None
                fd = row["f_dianion"].iloc[0] * 100 if "f_dianion" in row else None
                parts = []
                if fn is not None: parts.append(f"f(AH₂)={fn:.1f}%")
                if fa is not None: parts.append(f"f(AH⁻)={fa:.1f}%")
                if fd is not None: parts.append(f"f(A²⁻)={fd:.1f}%")
                print(f"    {label}: {', '.join(parts)}")

        # ── Monte Carlo (cadena de 3 estados: no hay catión A ni macroestado,
        #    se reusa run_chain_montecarlo() con data_A=None, f_R=1/f_D=0) ──
        chain_mc = None
        if args.montecarlo:
            print(f"\n▸ Propagando incertidumbre Monte Carlo "
                  f"(N={args.mc_samples:,} muestras)...")
            chain_mc = run_chain_montecarlo(
                None, data_neutral, data_anion, data_dianion,
                f_R=1.0, f_D=0.0, pH_values=pH_all,
                T=args.temp, method=args.method, gcorr_source=gcorr_source,
                n_samples=args.mc_samples,
                sigma_elec_kcal=args.mc_sigma_elec,
                sigma_gcorr_kcal=args.mc_sigma_gcorr,
                sigma_hplus_kcal=args.mc_sigma_hplus,
                seed=args.mc_seed,
            )
            print_mc_summary(chain_mc, [
                ("pKa1", pka1_result, "pKa₁"),
                ("pKa2", pka2_result, "pKa₂"),
            ])
    else:
        df_dist = pd.DataFrame()
        chain_mc = None

    data_list = [d for d in [data_neutral, data_anion, data_dianion] if d is not None]

    if args.export and not df_dist.empty:
        export_results([(pka1_result, "pKa1"), (pka2_result, "pKa2")],
                       df_dist, data_list, args.export, chain_mc=chain_mc)

    if args.plot and pKa1 is not None:
        out_dir = PC_WORK_ROOT
        plot_distribution(df_dist, pKa1, pKa2,
                          output_path=out_dir / "species_distribution.png",
                          ph_highlight=args.ph, chain_mc=chain_mc)
        plot_conformer_boltzmann(data_list,
                                 output_path=out_dir / "conformer_boltzmann.png")
        if args.montecarlo:
            plot_pka_histograms(chain_mc, [
                ("pKa1", pka1_result, "pKa₁"),
                ("pKa2", pka2_result, "pKa₂"),
            ], output_path=out_dir / "pka_mc_histograms.png")

    print("\nTerminado.\n")

    return {
        "gcorr_source": gcorr_source,
        "pKa1": pKa1, "pKa2": pKa2,
        "pKa1_raw": pKa1_raw, "pKa2_raw": pKa2_raw,
        "correction": {"pKa1": corr1, "pKa2": corr2},
        "df_dist": df_dist,
        "data_list": data_list,
        "chain_mc": chain_mc,
    }


def _build_neutral_macrostate(data_R: dict, data_D: dict, neutral: str,
                              zwitterion: str, T: float, verbose: bool,
                              gcorr_source: str = "harmonic") -> tuple[dict, dict]:
    """
    Construye el macroestado neutro/zwitterión (R ⇌ D, tautómeros) usado como
    estado intermedio de la cadena A → [R/D] → B → C. Retorna
    (taut_result, data_neutral_macro). Extraído tal cual de _run_extended()
    para que `calibrate` use exactamente la misma construcción. Si ni siquiera
    R tiene G_total, termina con sys.exit(1) como antes.
    """
    # ── Equilibrio de tautómeros R ⇌ D (sin H⁺) ────────────────────────────
    print(f"\n▸ Calculando equilibrio de tautómeros  ({neutral} ⇌ {zwitterion})  @ T={T:.2f} K")
    taut_result = compute_tautomer_equilibrium(data_R, data_D, T=T, verbose=verbose,
                                               gcorr_source=gcorr_source)

    if not taut_result["ok"]:
        # No se pudo cerrar el equilibrio R⇌D (típicamente porque D no tiene
        # confórmero representativo). Esto NO debe abortar todo el script:
        # se cae a usar R solo como estado "neutro" (f_R=100%, f_D=0%) y se
        # sigue con la cadena de pKa para los estados que sí tienen datos.
        print(f"\n  ⚠  No se pudo calcular el equilibrio de tautómeros "
              f"({neutral}⇌{zwitterion}). Se continúa usando solo "
              f"'{neutral}' como estado neutro "
              f"(f_{neutral}=100%, f_{zwitterion}=0%).")
        for w in taut_result["warnings"]:
            print(f"    ⚠  {w}")

        gR_fallback = _state_G_total_direct(data_R, T=T, gcorr_source=gcorr_source)
        if not gR_fallback.get("ok"):
            print(f"\n  ⚠  Tampoco se pudo obtener G_total para "
                  f"'{neutral}'. Abortando.")
            sys.exit(1)

        taut_result = {
            "mol_R": neutral, "mol_D": zwitterion, "T": T,
            "ok": True,
            "warnings": taut_result["warnings"] + [
                f"Macroestado neutro construido solo con '{neutral}' "
                f"(sin datos válidos para '{zwitterion}')."
            ],
            "conf_R": gR_fallback["conf"], "conf_D": None,
            "G_R_Eh": gR_fallback["G_total_Eh"], "G_D_Eh": None,
            "dG_taut_kcal": None, "K_taut": None,
            "f_R": 1.0, "f_D": 0.0,
            "G_macro_Eh": gR_fallback["G_total_Eh"],
            "theory_R": gR_fallback["theory_sp"], "theory_D": None,
        }
        print(f"    → G_total({neutral}, confórmero {gR_fallback['conf']}) "
              f"= {gR_fallback['G_total_Eh']:.8f} Eh  usado como G_macro")
    elif not verbose:
        _print_tautomer_cycle(taut_result)

    if taut_result["conf_D"] is not None:
        conf_label = (f"{neutral}:{taut_result['conf_R']} + "
                      f"{zwitterion}:{taut_result['conf_D']} (Boltzmann)")
        mol_label  = f"{neutral}/{zwitterion}"
    else:
        conf_label = f"{neutral}:{taut_result['conf_R']} (sin {zwitterion})"
        mol_label  = neutral

    data_neutral_macro = {
        "mol_name":          mol_label,
        "is_macrostate":     True,
        "G_total_direct_Eh": taut_result["G_macro_Eh"],
        "conf_label":        conf_label,
        "theory_label":      taut_result["theory_R"] or taut_result["theory_D"],
        # R y D son tautómeros (mismo q, mismo N° átomos, ya validado en
        # compute_tautomer_equilibrium), así que cualquiera de los dos sirve
        # como representante del macroestado para la validación cruzada de
        # compute_pka() frente al catión (A) y al anión (B).
        "charge":            data_R.get("charge"),
        "n_atoms":           data_R.get("n_atoms"),
        "use_smd":           data_R.get("use_smd", False),
        "aq_source":         data_R.get("aq_source", "CPCM"),
    }

    return taut_result, data_neutral_macro


def _run_extended(args, gcorr_source: str = "harmonic"):
    """
    Cadena de 5 estados:
      A⁺ (catión) ⇌ [R (neutro) ⇌ D (zwitterión)] ⇌ B⁻ (anión) ⇌ C²⁻ (dianión)

    gcorr_source: ver docstring de _run_classic(). Retorna dict con los
    resultados clave (pKa0/1/2, f_R, f_D, df_dist, ...) para _run_comparison().
    """
    data_A = None
    if args.protonated:
        print(f"\n▸ Cargando catión protonado: {args.protonated}")
        data_A = collect_mol_data(args.protonated, verbose=args.verbose,
                                  gcorr_source=gcorr_source, use_smd=args.spsmd)
        print_mol_summary(data_A, verbose=args.verbose)
    else:
        print("\n▸ Catión protonado A⁺ no proporcionado (--protonated omitido): "
              "se omite pKa₀, la cadena arranca en el macroestado neutro/zwitterión.")

    print(f"\n▸ Cargando estado neutro: {args.neutral}")
    data_R = collect_mol_data(args.neutral, verbose=args.verbose,
                              gcorr_source=gcorr_source, use_smd=args.spsmd)
    print_mol_summary(data_R, verbose=args.verbose)

    print(f"\n▸ Cargando zwitterión (tautómero de {args.neutral}): {args.zwitterion}")
    data_D = collect_mol_data(args.zwitterion, verbose=args.verbose,
                              gcorr_source=gcorr_source, use_smd=args.spsmd)
    print_mol_summary(data_D, verbose=args.verbose)

    data_B = data_C = None
    if args.anion:
        print(f"\n▸ Cargando anión: {args.anion}")
        data_B = collect_mol_data(args.anion, verbose=args.verbose,
                                  gcorr_source=gcorr_source, use_smd=args.spsmd)
        print_mol_summary(data_B, verbose=args.verbose)

    if args.dianion:
        print(f"\n▸ Cargando dianión: {args.dianion}")
        data_C = collect_mol_data(args.dianion, verbose=args.verbose,
                                  gcorr_source=gcorr_source, use_smd=args.spsmd)
        print_mol_summary(data_C, verbose=args.verbose)

    taut_result, data_neutral_macro = _build_neutral_macrostate(
        data_R, data_D, args.neutral, args.zwitterion, args.temp,
        verbose=args.verbose, gcorr_source=gcorr_source)

    # ── Cadena de pKa: A→[R/D], [R/D]→B, B→C ───────────────────────────────
    pka0_result = pka1_result = pka2_result = None

    if data_A:
        print(f"\n▸ Calculando pKa₀  ({args.protonated} → neutro/zwitterión + H⁺)  @ T={args.temp:.2f} K")
        pka0_result = compute_pka(data_A, data_neutral_macro,
                                  T=args.temp, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if data_B:
        print(f"\n▸ Calculando pKa₁  (neutro/zwitterión → {args.anion} + H⁺)  @ T={args.temp:.2f} K")
        pka1_result = compute_pka(data_neutral_macro, data_B,
                                  T=args.temp, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

    if data_B and data_C:
        print(f"\n▸ Calculando pKa₂  ({args.anion} → {args.dianion} + H⁺)  @ T={args.temp:.2f} K")
        pka2_result = compute_pka(data_B, data_C,
                                  T=args.temp, verbose=args.verbose, method=args.method,
                                  gcorr_source=gcorr_source)

        # Nota #2 del review: pKa1 siempre usa 'direct' (macroestado), pero
        # pKa2 puede resolver a 'born' con method='auto'. Ambos métodos usan
        # constantes G(H+,aq) distintas y no intercambiables — este chequeo
        # cuantifica la discrepancia cuando aplica, sin alterar pka2_result.
        ph_scale_diag = _cross_validate_ph_scale(
            data_B, data_C, pka2_result, T=args.temp, gcorr_source=gcorr_source,
        )
        if ph_scale_diag.get("applicable"):
            _print_ph_scale_diagnostic(ph_scale_diag, label="pKa₂ (B→C)")

        # Con method='auto', pKa1 (macroestado) y pKa2 (simple) pueden ahora
        # resolver a métodos distintos ('born' vs 'direct') dependiendo de si
        # existen SP de fase gas para B y C. Antes esto quedaba escondido
        # porque pKa1 SIEMPRE caía a 'direct' (bug ya corregido arriba); con
        # el fix, mezclar 'born' (G_HPLUS_AQ_BORN) y 'direct'
        # (G_HPLUS_AQ_DIRECT) en la misma cadena de pKa's suma dos escalas de
        # G(H+,aq) no intercambiables sin que nada lo señale. Se aborta en
        # vez de solo advertir.
        if (args.method == "auto" and pka1_result is not None
                and pka1_result.get("ok") and pka2_result.get("ok")
                and pka1_result.get("method") != pka2_result.get("method")):
            raise SystemExit(
                f"\n❌ INCONSISTENT G(H+,aq) REFERENCE within the same pKa "
                f"chain under method='auto':\n"
                f"   pKa₁ ([R/D]→{args.anion}) resolved to "
                f"'{pka1_result.get('method')}'\n"
                f"   pKa₂ ({args.anion}→{args.dianion}) resolved to "
                f"'{pka2_result.get('method')}'\n"
                f"   These use different, non-interchangeable G(H+,aq) "
                f"constants (Pliego-Riveros {G_HPLUS_AQ_DIRECT:.2f} kcal/mol "
                f"vs. Tissandier/Kelly {G_HPLUS_AQ_BORN:.2f} kcal/mol).\n"
                f"   Fix: either provide gas-phase SPs for all species in the "
                f"chain (so 'auto' resolves to 'born' everywhere), or force "
                f"a single method explicitly with --method direct or "
                f"--method born."
            )

    # Nota #8 del review: reconstruye tautómero+macroestado y recalcula
    # pKa0/1/2 en las 4 T de TEMPS_K (no solo re-evalúa RT·ln10 a T fija).
    ext_temp_scan = _extended_temperature_scan(
        data_A, data_R, data_D, data_B, data_C,
        gcorr_source=gcorr_source, method=args.method,
    )
    _print_extended_temperature_scan(ext_temp_scan, reported_T=args.temp)

    print_pka_summary([
        (pka0_result, "pKa₀ (A→R/D)"),
        (pka1_result, "pKa₁ (R/D→B)"),
        (pka2_result, "pKa₂ (B→C)"),
    ])

    pKa0_raw = pka0_result["pKa"] if (pka0_result and pka0_result["ok"]) else None
    pKa1_raw = pka1_result["pKa"] if (pka1_result and pka1_result["ok"]) else None
    pKa2_raw = pka2_result["pKa"] if (pka2_result and pka2_result["ok"]) else None

    # ── Corrección empírica (--empirical-correction) ───────────────────────
    ec_spec = getattr(args, "empirical_correction_spec", {}) or {}
    ec_mode = getattr(args, "empirical_correction_mode", "auto")
    corr0 = corr1 = corr2 = None
    if ec_spec:
        if "pKa0" in ec_spec:
            corr0 = apply_empirical_correction(pKa0_raw, ec_spec["pKa0"],
                                               gcorr_source=gcorr_source, mode=ec_mode)
        if "pKa1" in ec_spec:
            corr1 = apply_empirical_correction(pKa1_raw, ec_spec["pKa1"],
                                               gcorr_source=gcorr_source, mode=ec_mode)
        if "pKa2" in ec_spec:
            corr2 = apply_empirical_correction(pKa2_raw, ec_spec["pKa2"],
                                               gcorr_source=gcorr_source, mode=ec_mode)
        _print_empirical_correction(corr0, "pKa₀ (A→R/D)")
        _print_empirical_correction(corr1, "pKa₁ (R/D→B)")
        _print_empirical_correction(corr2, "pKa₂ (B→C)")
        if any(c and c.get("ok") for c in (corr0, corr1, corr2)):
            print(f"\n  Nota: la distribución de especies, el pI y los "
                  f"gráficos de aquí en adelante usan los pKa's CORREGIDOS "
                  f"EMPÍRICAMENTE donde se pidió corrección (ver arriba). "
                  f"El escaneo de temperatura (arriba) sigue mostrando "
                  f"valores CRUDOS Born/CPCM sin corregir: la calibración "
                  f"empírica solo se ajustó a T=298.15 K. Si además usas "
                  f"--montecarlo, esa propagación de incertidumbre corre "
                  f"sobre el ciclo de Born crudo (no conoce la corrección "
                  f"empírica); súmale el RMSE de calibración reportado "
                  f"arriba como una fuente de error sistemático aparte.")
        for c, pka_result in ((corr0, pka0_result), (corr1, pka1_result), (corr2, pka2_result)):
            if c and c.get("ok") and pka_result is not None:
                pka_result["pKa_raw"] = c["pKa_raw"]
                pka_result["pKa_corrected"] = c["pKa_corrected"]
                pka_result["correction"] = c

    pKa0 = corr0["pKa_corrected"] if (corr0 and corr0.get("ok")) else pKa0_raw
    pKa1 = corr1["pKa_corrected"] if (corr1 and corr1.get("ok")) else pKa1_raw
    pKa2 = corr2["pKa_corrected"] if (corr2 and corr2.get("ok")) else pKa2_raw
    f_R, f_D = taut_result["f_R"], taut_result["f_D"]

    if pKa0 is None and pKa1 is None:
        print("\n  ⚠  No se pudo calcular ningún pKa de la cadena principal. "
              "Revisa los directorios.")
        sys.exit(1)

    # ── Distribución de especies ──────────────────────────────────────────
    pH_range = np.arange(args.ph_range[0], args.ph_range[1] + args.ph_step / 2, args.ph_step)
    pH_highlights = [args.ph, PH_STANDARD]
    for pk in (pKa0, pKa1, pKa2):
        if pk is not None:
            pH_highlights.append(round(pk, 1))
    pH_all = sorted(set(list(pH_range) + pH_highlights))
    pH_all = [p for p in pH_all if args.ph_range[0] <= p <= args.ph_range[1]]

    df_dist = boltzmann_5state(pKa0, pKa1, pKa2, f_R, f_D, pH_all)

    pH_display = np.arange(args.ph_range[0], args.ph_range[1] + 0.5, 0.5)
    pH_display_list = sorted(set([round(p, 1) for p in pH_display] + pH_highlights))
    pH_display_list = [p for p in pH_display_list if args.ph_range[0] <= p <= args.ph_range[1]]
    df_display = boltzmann_5state(pKa0, pKa1, pKa2, f_R, f_D, pH_display_list)
    print_distribution_table_5(df_display, pKa0, pKa1, pKa2, ph_highlight=args.ph)

    print(f"\n  ── Fracciones a pH de referencia ──")
    for pH_ref, label in [
        (args.ph,     f"pH {args.ph} (agua de mar antártica)"),
        (PH_STANDARD, f"pH {PH_STANDARD} (estándar)"),
    ]:
        row = boltzmann_5state(pKa0, pKa1, pKa2, f_R, f_D, [pH_ref])
        if not row.empty:
            r0 = row.iloc[0]
            parts = []
            if "f_cation" in row: parts.append(f"f(A⁺)={r0['f_cation']*100:.1f}%")
            parts.append(f"f(R)={r0['f_R']*100:.1f}%")
            parts.append(f"f(D)={r0['f_D']*100:.1f}%")
            if "f_anion" in row: parts.append(f"f(B⁻)={r0['f_anion']*100:.1f}%")
            if "f_dianion" in row: parts.append(f"f(C²⁻)={r0['f_dianion']*100:.1f}%")
            print(f"    {label}: {', '.join(parts)}")

    # ── Monte Carlo (cadena completa de 5 estados, G(H⁺,aq) correlacionado
    #    entre pKa0/pKa1/pKa2 — ver run_chain_montecarlo()) ─────────────────
    chain_mc = None
    if args.montecarlo:
        print(f"\n▸ Propagando incertidumbre Monte Carlo por toda la cadena "
              f"(N={args.mc_samples:,} muestras)...")
        chain_mc = run_chain_montecarlo(
            data_A, data_neutral_macro, data_B, data_C,
            f_R=f_R, f_D=f_D, pH_values=pH_all,
            T=args.temp, method=args.method, gcorr_source=gcorr_source,
            n_samples=args.mc_samples,
            sigma_elec_kcal=args.mc_sigma_elec,
            sigma_gcorr_kcal=args.mc_sigma_gcorr,
            sigma_hplus_kcal=args.mc_sigma_hplus,
            seed=args.mc_seed,
        )
        print_mc_summary(chain_mc, [
            ("pKa0", pka0_result, "pKa₀ (A→R/D)"),
            ("pKa1", pka1_result, "pKa₁ (R/D→B)"),
            ("pKa2", pka2_result, "pKa₂ (B→C)"),
        ])

    # ── Exportar ─────────────────────────────────────────────────────────
    data_list = [d for d in [data_A, data_R, data_D, data_B, data_C] if d is not None]

    if args.export and not df_dist.empty:
        export_results(
            [(pka0_result, "pKa0"), (pka1_result, "pKa1"), (pka2_result, "pKa2")],
            df_dist, data_list, args.export, taut_result=taut_result,
            chain_mc=chain_mc,
        )

    # ── Gráficos ─────────────────────────────────────────────────────────
    if args.plot:
        out_dir = PC_WORK_ROOT
        plot_distribution_5(
            df_dist, pKa0, pKa1, pKa2,
            output_path=out_dir / "species_distribution.png",
            ph_highlight=args.ph,
            chain_mc=chain_mc,
        )
        plot_conformer_boltzmann(
            data_list, output_path=out_dir / "conformer_boltzmann.png",
        )
        if args.montecarlo:
            plot_pka_histograms(chain_mc, [
                ("pKa0", pka0_result, "pKa₀ (A→R/D)"),
                ("pKa1", pka1_result, "pKa₁ (R/D→B)"),
                ("pKa2", pka2_result, "pKa₂ (B→C)"),
            ], output_path=out_dir / "pka_mc_histograms.png")

    print("\nTerminado.\n")

    return {
        "gcorr_source": gcorr_source,
        "pKa0": pKa0, "pKa1": pKa1, "pKa2": pKa2,
        "pKa0_raw": pKa0_raw, "pKa1_raw": pKa1_raw, "pKa2_raw": pKa2_raw,
        "correction": {"pKa0": corr0, "pKa1": corr1, "pKa2": corr2},
        "f_R": f_R, "f_D": f_D,
        "df_dist": df_dist,
        "data_list": data_list,
        "chain_mc": chain_mc,
    }


def _run_comparison(args, use_extended: bool):
    """
    Corre el ciclo COMPLETO dos veces con distinta fuente de G_corr —
    armónico estándar (NumFreq/RRHO) vs. cuasi-armónico (GoodVibes/mRRHO) —
    y compara al final cuánto cambian los pKa's y la distribución de especies.

    Activado con --goodvibes. Usa _run_classic() o _run_extended() según
    corresponda, sin duplicar lógica: solo cambia gcorr_source entre pasadas.
    """
    gv_python = _resolve_gv_python()
    try:
        check = subprocess.run(
            [gv_python, "-c", "import goodvibes"],
            capture_output=True, timeout=10,
        )
    except (subprocess.TimeoutExpired, OSError) as e:
        print(f"\n  ⚠  No se pudo verificar goodvibes con {gv_python}: {e}\n")
        sys.exit(1)
    if check.returncode != 0:
        print(f"\n  ⚠  No se puede importar goodvibes con {gv_python}.")
        print(f"     Instalar con: conda activate {GV_CONDA_ENV or '<env>'} && "
              f"pip install goodvibes --break-system-packages")
        print("     (o, en un entorno conda: conda install -c conda-forge goodvibes)\n")
        sys.exit(1)

    runner = _run_extended if use_extended else _run_classic

    print("\n" + "═" * 74)
    print("  PASADA 1/2  —  G_corr ARMÓNICO ESTÁNDAR (NumFreq, RRHO)")
    print("═" * 74)
    res_harm = runner(args, gcorr_source="harmonic")

    qh_tag = ", +qh (Head-Gordon)" if args.gv_qh else ""
    print("\n" + "═" * 74)
    print(f"  PASADA 2/2  —  G_corr CUASI-ARMÓNICO (GoodVibes, "
          f"qs={args.gv_qs}, cutoff={args.gv_cutoff:.0f} cm⁻¹{qh_tag})")
    print("═" * 74)
    res_gv = runner(args, gcorr_source="goodvibes")

    _print_goodvibes_comparison(res_harm, res_gv, extended=use_extended, args=args)


def _print_goodvibes_comparison(res_harm: dict, res_gv: dict,
                                 extended: bool, args) -> None:
    """
    Imprime una tabla comparativa entre la pasada armónica (NumFreq) y la
    cuasi-armónica (GoodVibes): Δ en cada pKa y Δ en las fracciones de
    especie a los pH de referencia (--ph y pH 7.0 estándar).
    """
    sep = "═" * 74
    print(f"\n{sep}")
    print("  COMPARACIÓN  —  G_corr ARMÓNICO (NumFreq) vs. CUASI-ARMÓNICO (GoodVibes)")
    print(sep)

    if res_harm is None or res_gv is None:
        print("\n  ⚠  Una de las dos pasadas no produjo resultados; no se puede comparar.\n")
        return

    pka_keys = ([("pKa0", "pKa₀ (A→R/D)")] if extended else []) + [
        ("pKa1", "pKa₁"), ("pKa2", "pKa₂"),
    ]

    print(f"\n  {'':<20} {'armónico':>12} {'GoodVibes':>12} {'Δ (GV−arm)':>12}")
    print(f"  {'-'*20} {'-'*12} {'-'*12} {'-'*12}")
    any_pka = False
    for key, label in pka_keys:
        v1, v2 = res_harm.get(key), res_gv.get(key)
        if v1 is None and v2 is None:
            continue
        any_pka = True
        v1s = f"{v1:.2f}" if v1 is not None else "—"
        v2s = f"{v2:.2f}" if v2 is not None else "—"
        ds  = f"{(v2 - v1):+.2f}" if (v1 is not None and v2 is not None) else "—"
        print(f"  {label:<20} {v1s:>12} {v2s:>12} {ds:>12}")
    if not any_pka:
        print("  (sin pKa's calculados en ninguna de las dos pasadas)")

    # ── Fracciones de especie a los pH de referencia ────────────────────
    def _row_at_ph(df: Optional[pd.DataFrame], ph: float):
        if df is None or df.empty or "pH" not in df.columns:
            return None
        sub = df[np.isclose(df["pH"].astype(float), ph, atol=1e-6)]
        return sub.iloc[0] if not sub.empty else None

    for ph_ref, ph_label in [
        (args.ph,     f"pH {args.ph} (referencia)"),
        (PH_STANDARD, f"pH {PH_STANDARD} (estándar)"),
    ]:
        row_h = _row_at_ph(res_harm.get("df_dist"), ph_ref)
        row_g = _row_at_ph(res_gv.get("df_dist"), ph_ref)
        if row_h is None or row_g is None:
            continue
        frac_cols = [c for c in row_h.index if c.startswith("f_")]
        if not frac_cols:
            continue
        print(f"\n  Fracciones de especie a {ph_label}:")
        print(f"  {'especie':<14} {'armónico':>10} {'GoodVibes':>10} {'Δ (pp)':>10}")
        for c in frac_cols:
            fh = float(row_h[c]) * 100 if pd.notna(row_h[c]) else None
            fg = float(row_g[c]) * 100 if (c in row_g and pd.notna(row_g[c])) else None
            if fh is None and fg is None:
                continue
            fhs = f"{fh:.2f}%" if fh is not None else "—"
            fgs = f"{fg:.2f}%" if fg is not None else "—"
            dps = f"{(fg - fh):+.2f}" if (fh is not None and fg is not None) else "—"
            print(f"  {c:<14} {fhs:>10} {fgs:>10} {dps:>10}")

    print(f"\n  Nota: la aproximación armónica pura sobreestima la entropía de los")
    print(f"  modos de baja frecuencia (rotores internos, <{args.gv_cutoff:.0f} cm⁻¹) —")
    print(f"  típicos en las cadenas laterales Gly/Val y grupos OH/NH libres de")
    print(f"  M-Gly:Val. Diferencias de pKa >0.3-0.5 unidades, o de fracción de")
    print(f"  especie >5 puntos porcentuales entre ambas columnas, son razón para")
    print(f"  preferir el valor de GoodVibes (mRRHO) en el paper en vez del RRHO puro.")
    print(f"{sep}\n")


# ══════════════════════════════════════════════════════════════════════════════
#  §  CUBES
#     Extracción de cubos HOMO/LUMO y NTO hole/particle vía orca_plot
# ══════════════════════════════════════════════════════════════════════════════

_cubes_logger = logging.getLogger("extract_tddft_cubes")

# --------------------------------------------------------------------------- #
# Estructuras de datos
# --------------------------------------------------------------------------- #


@dataclass
class NTOPair:
    hole: int
    particle: int
    coef: float


@dataclass
class ExcitedState:
    state_idx: int
    ntos: list[NTOPair] = _dc_field(default_factory=list)

    @property
    def participation_ratio(self) -> Optional[float]:
        """PR = 1 / sum(lambda_i^2). lambda_i normalizado a partir de coef."""
        if not self.ntos:
            return None
        total = sum(p.coef for p in self.ntos)
        if total <= 0:
            return None
        lambdas = [p.coef / total for p in self.ntos]
        denom = sum(l ** 2 for l in lambdas)
        if denom <= 0:
            return None
        return 1.0 / denom


@dataclass
class ConformerJob:
    conf_dir: Path
    molecule: str
    run_label: str
    conf_id: str
    homo_idx: Optional[int] = None
    lumo_idx: Optional[int] = None
    is_open_shell: bool = False
    states: dict[int, ExcitedState] = _dc_field(default_factory=dict)


@dataclass
class ReportRow:
    molecule: str
    run_label: str
    conf_id: str
    kind: str  # "ground_state" | "state1", "state2", ...
    homo_or_hole: Optional[int]
    lumo_or_particle: Optional[int]
    coef: Optional[float]
    participation_ratio: Optional[float]
    fosc: Optional[float]
    cube_files: str
    status: str
    error: str = ""


# --------------------------------------------------------------------------- #
# Descubrimiento de carpetas
# --------------------------------------------------------------------------- #


def discover_conformers(root: Path) -> list[Path]:
    """Encuentra carpetas conf_* bajo */runs/*/05_tddft/ con tddft.out + tddft.gbw."""
    pattern = "*/runs/*/05_tddft/conf_*"
    candidates = sorted(root.glob(pattern))
    valid = []
    for conf_dir in candidates:
        if not conf_dir.is_dir():
            continue
        out_file = conf_dir / "tddft.out"
        gbw_file = conf_dir / "tddft.gbw"
        if out_file.exists() and gbw_file.exists():
            valid.append(conf_dir)
        else:
            missing = [
                name
                for name, f in (("tddft.out", out_file), ("tddft.gbw", gbw_file))
                if not f.exists()
            ]
            _cubes_logger.warning(
                "Saltando %s: faltan %s", conf_dir, ", ".join(missing)
            )
    return valid


def parse_conformer_path(conf_dir: Path, root: Path) -> tuple[str, str, str]:
    """Extrae (molecula, run_label, conf_id) de PAPER1/MOL/runs/B/05_tddft/conf_6."""
    rel = conf_dir.relative_to(root)
    parts = rel.parts
    # parts: (MOL, "runs", RUN_LABEL, "05_tddft", "conf_X")
    molecule = parts[0]
    run_label = parts[2] if len(parts) > 2 else "unknown"
    conf_id = conf_dir.name
    return molecule, run_label, conf_id


# --------------------------------------------------------------------------- #
# Parseo de tddft.out
# --------------------------------------------------------------------------- #

ORBITAL_ENERGIES_HEADER = re.compile(r"ORBITAL ENERGIES")
# Filas tipo:   6   2.0000    -0.522130      -14.2079
ORBITAL_ROW = re.compile(
    r"^\s*(\d+)\s+([\-\d.]+)\s+([\-\d.]+)\s+([\-\d.]+)\s*$"
)
NTO_STATE_HEADER = re.compile(r"NATURAL TRANSITION ORBITALS FOR STATE\s+(\d+)")
# Filas tipo:   54a -> 55a  :  n= 0.523456
NTO_ROW = re.compile(
    r"^\s*(\d+)([ab])\s*->\s*(\d+)([ab])\s*:\s*n\s*=\s*([\-\d.]+)"
)

ABSORPTION_HEADER = re.compile(
    r"ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS"
)
# Fila tipo:  0-1A  ->  1-1A    4.203068   33900.0   295.0   0.752016216 ...
ABSORPTION_ROW = re.compile(
    r"^\s*\d+-\S+\s*->\s*(\d+)-\S+\s+([\-\d.]+)\s+([\-\d.]+)\s+([\-\d.]+)\s+([\-\d.]+)"
)


def parse_oscillator_strengths(out_text: str) -> dict[int, float]:
    """Extrae fosc por estado desde la primera tabla 'ABSORPTION SPECTRUM VIA
    TRANSITION ELECTRIC DIPOLE MOMENTS' (dipolo eléctrico, no velocidad ni CD).
    Devuelve {estado: fosc}."""
    lines = out_text.splitlines()
    header_positions = [i for i, l in enumerate(lines) if ABSORPTION_HEADER.search(l)]
    if not header_positions:
        return {}

    start = header_positions[0]
    fosc_by_state: dict[int, float] = {}
    dash_count = 0
    for line in lines[start:]:
        if "-----" in line:
            dash_count += 1
            if dash_count >= 3:
                break  # tercera línea de guiones = fin de la tabla de datos
            continue
        if dash_count < 2:
            continue  # todavía en encabezado (nombre, unidades)
        m = ABSORPTION_ROW.match(line)
        if m:
            state_idx = int(m.group(1))
            fosc = float(m.group(5))
            fosc_by_state[state_idx] = fosc
    return fosc_by_state


def select_dominant_state(
    fosc_by_state: dict[int, float],
    fosc_min: float,
) -> Optional[int]:
    """Selecciona el único estado 'dominante' (el de deltamax): el de mayor
    fosc, solo si ese máximo supera fosc_min. Si no hay tabla de fosc,
    devuelve None con significado especial de 'sin filtrar' (fallback inseguro
    ya no aplica aquí: sin tabla no hay forma de saber cuál es el dominante,
    así que no se extrae ningún NTO y se advierte)."""
    if not fosc_by_state:
        return None
    brightest = max(fosc_by_state, key=fosc_by_state.get)
    if fosc_by_state[brightest] <= fosc_min:
        return None
    return brightest


def parse_ground_state_frontier_orbitals(
    out_text: str,
) -> tuple[Optional[int], Optional[int], bool]:
    """Devuelve (homo_idx, lumo_idx, is_open_shell) desde el último bloque
    ORBITAL ENERGIES en tddft.out."""
    lines = out_text.splitlines()
    header_positions = [
        i for i, l in enumerate(lines) if ORBITAL_ENERGIES_HEADER.search(l)
    ]
    if not header_positions:
        return None, None, False

    # Tomar el bloque a partir de la última ocurrencia
    start = header_positions[-1]
    block = lines[start:]

    # Heurística open-shell: ORCA imprime "SPIN UP ORBITALS" / "SPIN DOWN ORBITALS"
    is_open_shell = any(
        "SPIN UP ORBITALS" in l or "SPIN DOWN ORBITALS" in l for l in block
    )

    homo_idx: Optional[int] = None
    lumo_idx: Optional[int] = None
    prev_occ_idx: Optional[int] = None

    for line in block:
        m = ORBITAL_ROW.match(line)
        if not m:
            # Si ya encontramos HOMO/LUMO y llegamos a una sección nueva, cortar
            if homo_idx is not None and lumo_idx is not None:
                break
            continue
        idx = int(m.group(1))
        occ = float(m.group(2))
        if occ > 0:
            prev_occ_idx = idx
        elif occ == 0 and prev_occ_idx is not None and lumo_idx is None:
            homo_idx = prev_occ_idx
            lumo_idx = idx
            if not is_open_shell:
                break

    return homo_idx, lumo_idx, is_open_shell


def parse_nto_states(out_text: str) -> dict[int, ExcitedState]:
    """Extrae, para cada estado excitado, TODOS los pares NTO reportados
    (con su coeficiente), ordenados por coeficiente descendente."""
    lines = out_text.splitlines()
    states: dict[int, ExcitedState] = {}

    current_state: Optional[int] = None
    for line in lines:
        header_m = NTO_STATE_HEADER.search(line)
        if header_m:
            current_state = int(header_m.group(1))
            states[current_state] = ExcitedState(state_idx=current_state)
            continue
        if current_state is None:
            continue
        row_m = NTO_ROW.match(line)
        if row_m:
            hole = int(row_m.group(1))
            particle = int(row_m.group(3))
            coef = float(row_m.group(5))
            states[current_state].ntos.append(
                NTOPair(hole=hole, particle=particle, coef=coef)
            )

    # Ordenar cada lista por coeficiente descendente (por si el output no
    # viniera ya ordenado)
    for state in states.values():
        state.ntos.sort(key=lambda p: p.coef, reverse=True)

    return states


def select_ntos_for_state(
    state: ExcitedState,
    threshold: float,
    min_coef: float,
    max_pairs: Optional[int],
) -> list[NTOPair]:
    """Selecciona los pares NTO a extraer: el dominante siempre, más los que
    cumplan coef >= threshold * coef_dominante Y coef >= min_coef."""
    if not state.ntos:
        return []
    dominant = state.ntos[0]
    selected = [dominant]
    for pair in state.ntos[1:]:
        if pair.coef >= min_coef and pair.coef >= threshold * dominant.coef:
            selected.append(pair)
    if max_pairs is not None:
        selected = selected[:max_pairs]
    return selected


def match_ntos_to_existing_files(
    conf_dir: Path, states: dict[int, ExcitedState]
) -> dict[int, ExcitedState]:
    """Cruza los estados parseados con los archivos tddft.sN.nto que
    realmente existen, para no inventar estados no calculados."""
    existing_state_ids = set()
    for f in conf_dir.glob("tddft.s*.nto"):
        m = re.match(r"tddft\.s(\d+)\.nto", f.name)
        if m:
            existing_state_ids.add(int(m.group(1)))

    filtered = {
        idx: state for idx, state in states.items() if idx in existing_state_ids
    }
    missing = set(states.keys()) - existing_state_ids
    if missing:
        _cubes_logger.warning(
            "%s: estados en tddft.out sin archivo .nto correspondiente: %s",
            conf_dir,
            sorted(missing),
        )
    return filtered


# --------------------------------------------------------------------------- #
# Construcción y ejecución de orca_plot
# --------------------------------------------------------------------------- #


def build_orca_plot_stdin(orbital_indices: list[int], grid: int) -> str:
    """Arma el bloque de comandos interactivo de orca_plot para una lista de
    índices de orbitales, reutilizando la sesión."""
    commands = ["5", "7", "4", str(grid)]
    for idx in orbital_indices:
        commands += ["1", "1", "2", str(idx), "3", "0", "11"]
    commands.append("12")
    return "\n".join(commands) + "\n"


def run_orca_plot(
    orca_plot_bin: str,
    target_file: Path,
    orbital_indices: list[int],
    grid: int,
    timeout: int,
) -> tuple[bool, str]:
    """Ejecuta orca_plot en modo interactivo vía stdin sobre target_file.
    Devuelve (exito, mensaje)."""
    stdin_text = build_orca_plot_stdin(orbital_indices, grid)
    try:
        result = subprocess.run(
            [orca_plot_bin, target_file.name, "-i"],
            input=stdin_text,
            cwd=target_file.parent,
            capture_output=True,
            text=True,
            timeout=timeout,
        )
    except subprocess.TimeoutExpired:
        return False, f"timeout tras {timeout}s"
    except FileNotFoundError:
        return False, f"binario orca_plot no encontrado: {orca_plot_bin}"

    if result.returncode != 0:
        return False, f"returncode={result.returncode}: {result.stderr[-500:]}"
    if "PLOTTING FINISHED" not in result.stdout:
        return False, "no se encontró 'PLOTTING FINISHED' en stdout"
    return True, "ok"


def find_generated_cubes(target_file: Path, orbital_indices: list[int]) -> list[Path]:
    """orca_plot nombra los cubos como basename.mo{N}a.cube (o .b. si beta).
    Los busca en la carpeta del target_file."""
    found = []
    basename = target_file.stem  # ej. "tddft" o "tddft.s6"
    for idx in orbital_indices:
        for spin in ("a", "b"):
            candidate = target_file.parent / f"{basename}.mo{idx}{spin}.cube"
            if candidate.exists():
                found.append(candidate)
    return found


# --------------------------------------------------------------------------- #
# Procesamiento de un confórmero
# --------------------------------------------------------------------------- #


def process_conformer(
    conf_dir: Path,
    root: Path,
    args: argparse.Namespace,
) -> list[ReportRow]:
    rows: list[ReportRow] = []
    molecule, run_label, conf_id = parse_conformer_path(conf_dir, root)

    out_text = (conf_dir / "tddft.out").read_text(errors="replace")

    homo_idx, lumo_idx, is_open_shell = parse_ground_state_frontier_orbitals(out_text)
    if is_open_shell:
        _cubes_logger.warning(
            "%s: sistema open-shell detectado, asumiendo alpha por defecto "
            "(revisar manualmente si esto no es correcto)",
            conf_dir,
        )

    cubes_dir = conf_dir / "cubes"
    if args.dry_run:
        cubes_dir_display = str(cubes_dir)
    else:
        cubes_dir.mkdir(exist_ok=True)
        cubes_dir_display = str(cubes_dir)

    # --- Estado fundamental: HOMO/LUMO ------------------------------------
    if homo_idx is None or lumo_idx is None:
        rows.append(
            ReportRow(
                molecule=molecule,
                run_label=run_label,
                conf_id=conf_id,
                kind="ground_state",
                homo_or_hole=homo_idx,
                lumo_or_particle=lumo_idx,
                coef=None,
                participation_ratio=None,
                fosc=None,
                cube_files="",
                status="error",
                error="no se pudo identificar HOMO/LUMO en ORBITAL ENERGIES",
            )
        )
    else:
        gbw_target = conf_dir / "tddft.gbw"
        cube_targets = [
            conf_dir / f"tddft.mo{homo_idx}a.cube",
            conf_dir / f"tddft.mo{lumo_idx}a.cube",
        ]
        already_exist = all(c.exists() for c in cube_targets)
        final_names = [cubes_dir / "HOMO_mo{}a.cube".format(homo_idx),
                        cubes_dir / "LUMO_mo{}a.cube".format(lumo_idx)]
        final_exist = all(c.exists() for c in final_names)

        if final_exist and not args.force:
            status, error = "skipped_exists", ""
            cube_list = final_names
        elif args.dry_run:
            status, error = "dry_run", ""
            cube_list = final_names
        else:
            ok, msg = run_orca_plot(
                args.orca_plot, gbw_target, [homo_idx, lumo_idx], args.grid, args.timeout
            )
            if not ok:
                status, error = "error", msg
                cube_list = []
            else:
                generated = find_generated_cubes(gbw_target, [homo_idx, lumo_idx])
                cube_list = []
                labels = {homo_idx: "HOMO", lumo_idx: "LUMO"}
                for g in generated:
                    m = re.search(r"\.mo(\d+)([ab])\.cube$", g.name)
                    if not m:
                        continue
                    idx = int(m.group(1))
                    label = labels.get(idx, f"mo{idx}")
                    dest = cubes_dir / f"{label}_mo{idx}{m.group(2)}.cube"
                    shutil.move(str(g), str(dest))
                    cube_list.append(dest)
                status, error = ("ok", "") if cube_list else ("error", "cubos no encontrados tras orca_plot")

        rows.append(
            ReportRow(
                molecule=molecule,
                run_label=run_label,
                conf_id=conf_id,
                kind="ground_state",
                homo_or_hole=homo_idx,
                lumo_or_particle=lumo_idx,
                coef=None,
                participation_ratio=None,
                fosc=None,
                cube_files=";".join(str(c.name) for c in cube_list),
                status=status,
                error=error,
            )
        )

    # --- Estados excitados: pares NTO --------------------------------------
    all_states = parse_nto_states(out_text)
    states = match_ntos_to_existing_files(conf_dir, all_states)

    fosc_by_state = parse_oscillator_strengths(out_text)
    dominant_state = select_dominant_state(fosc_by_state, args.fosc_min)

    if not fosc_by_state:
        _cubes_logger.warning(
            "%s: no se encontró tabla de fosc, no se puede identificar el "
            "estado dominante (deltamax) — se omite la extracción de NTO "
            "para este confórmero", conf_dir,
        )
    elif dominant_state is None:
        brightest = max(fosc_by_state, key=fosc_by_state.get)
        _cubes_logger.warning(
            "%s: ningún estado supera fosc_min=%.3f (máximo: estado %d, "
            "fosc=%.4f) — se omite la extracción de NTO para este confórmero",
            conf_dir, args.fosc_min, brightest, fosc_by_state[brightest],
        )

    for state_idx in sorted(states):
        fosc = fosc_by_state.get(state_idx)

        if state_idx != dominant_state:
            rows.append(
                ReportRow(
                    molecule=molecule,
                    run_label=run_label,
                    conf_id=conf_id,
                    kind=f"state{state_idx}",
                    homo_or_hole=None,
                    lumo_or_particle=None,
                    coef=None,
                    participation_ratio=states[state_idx].participation_ratio,
                    fosc=fosc,
                    cube_files="",
                    status="skipped_not_dominant",
                    error="",
                )
            )
            continue

        state = states[state_idx]
        selected = select_ntos_for_state(
            state, args.nto_threshold, args.nto_min_coef, args.nto_max_pairs
        )
        pr = state.participation_ratio

        if not selected:
            rows.append(
                ReportRow(
                    molecule=molecule,
                    run_label=run_label,
                    conf_id=conf_id,
                    kind=f"state{state_idx}",
                    homo_or_hole=None,
                    lumo_or_particle=None,
                    coef=None,
                    participation_ratio=pr,
                    fosc=fosc,
                    cube_files="",
                    status="error",
                    error="sin pares NTO parseados",
                )
            )
            continue

        nto_target = conf_dir / f"tddft.s{state_idx}.nto"
        orbital_indices = []
        for pair in selected:
            orbital_indices += [pair.hole, pair.particle]

        final_names = []
        for i, pair in enumerate(selected, start=1):
            suffix = "" if i == 1 else f"_alt{i}"
            final_names.append(
                cubes_dir / f"state{state_idx}_hole{suffix}_mo{pair.hole}a.cube"
            )
            final_names.append(
                cubes_dir / f"state{state_idx}_particle{suffix}_mo{pair.particle}a.cube"
            )
        final_exist = all(c.exists() for c in final_names)

        if final_exist and not args.force:
            status, error = "skipped_exists", ""
            cube_list = final_names
        elif args.dry_run:
            status, error = "dry_run", ""
            cube_list = final_names
        else:
            ok, msg = run_orca_plot(
                args.orca_plot, nto_target, orbital_indices, args.grid, args.timeout
            )
            if not ok:
                status, error = "error", msg
                cube_list = []
            else:
                generated = find_generated_cubes(nto_target, orbital_indices)
                cube_list = []
                idx_to_role = {}
                for i, pair in enumerate(selected, start=1):
                    suffix = "" if i == 1 else f"_alt{i}"
                    idx_to_role[pair.hole] = f"hole{suffix}"
                    idx_to_role[pair.particle] = f"particle{suffix}"
                for g in generated:
                    m = re.search(r"\.mo(\d+)([ab])\.cube$", g.name)
                    if not m:
                        continue
                    idx = int(m.group(1))
                    role = idx_to_role.get(idx, f"mo{idx}")
                    dest = cubes_dir / f"state{state_idx}_{role}_mo{idx}{m.group(2)}.cube"
                    shutil.move(str(g), str(dest))
                    cube_list.append(dest)
                status, error = ("ok", "") if cube_list else ("error", "cubos no encontrados tras orca_plot")

        rows.append(
            ReportRow(
                molecule=molecule,
                run_label=run_label,
                conf_id=conf_id,
                kind=f"state{state_idx}",
                homo_or_hole=selected[0].hole,
                lumo_or_particle=selected[0].particle,
                coef=selected[0].coef,
                participation_ratio=pr,
                fosc=fosc,
                cube_files=";".join(str(c.name) for c in cube_list),
                status=status,
                error=error,
            )
        )

        caracter = "mixto" if len(selected) > 1 else "puro"
        _cubes_logger.info(
            "%s estado %d (dominante, fosc=%.4f): carácter %s (PR=%.2f), "
            "%d par(es) NTO extraídos",
            conf_dir, state_idx, fosc or -1, caracter, pr or -1, len(selected),
        )

    return rows


# --------------------------------------------------------------------------- #
# Reporte final
# --------------------------------------------------------------------------- #


def write_report(rows: list[ReportRow], out_path: Path) -> None:
    fieldnames = [
        "molecule", "run_label", "conf_id", "kind",
        "homo_or_hole", "lumo_or_particle", "coef", "participation_ratio",
        "fosc", "cube_files", "status", "error",
    ]
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "molecule": row.molecule,
                    "run_label": row.run_label,
                    "conf_id": row.conf_id,
                    "kind": row.kind,
                    "homo_or_hole": row.homo_or_hole,
                    "lumo_or_particle": row.lumo_or_particle,
                    "coef": row.coef,
                    "participation_ratio": (
                        f"{row.participation_ratio:.3f}"
                        if row.participation_ratio is not None
                        else ""
                    ),
                    "fosc": f"{row.fosc:.6f}" if row.fosc is not None else "",
                    "cube_files": row.cube_files,
                    "status": row.status,
                    "error": row.error,
                }
            )


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def _cubes_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Extrae cubos HOMO/LUMO y NTO hole/particle desde runs TDDFT."
    )
    parser.add_argument("--root", type=Path, default=Path("PAPER1"),
                         help="raíz de búsqueda (default: PAPER1)")
    parser.add_argument("--grid", type=int, default=80,
                         help="resolución de grilla para orca_plot (default: 80)")
    parser.add_argument(
        "--orca-plot",
        default="/home/kmilozeta/orca_6_1_1_avx2/orca_plot",
        help="ruta al binario orca_plot",
    )
    parser.add_argument(
        "--nto-threshold", type=float, default=0.15,
        help="fracción relativa al coeficiente dominante para incluir un par "
             "NTO adicional (default: 0.15, calibrado contra conf_5/B: separa "
             "limpiamente PR<=1.37 (puro) de PR>=1.39 (mixto) en ese conformero; "
             "ej. dominante=0.64 -> se incluyen pares con coef >= 0.096)",
    )
    parser.add_argument(
        "--nto-min-coef", type=float, default=0.05,
        help="piso absoluto de coeficiente para considerar un par NTO real, "
             "evita ruido numérico (default: 0.05)",
    )
    parser.add_argument(
        "--nto-max-pairs", type=int, default=2,
        help="tope duro de pares NTO por estado, aunque el umbral sugiera "
             "más (default: 2, evita nomenclatura de carpeta compleja)",
    )
    parser.add_argument(
        "--fosc-min", type=float, default=0.1,
        help="fuerza de oscilador mínima para considerar un estado como el "
             "'dominante' (deltamax) de un confórmero (default: 0.1). Solo se "
             "extraen NTO del único estado con mayor fosc, y solo si ese "
             "máximo supera este umbral; el resto queda como "
             "status=skipped_not_dominant en el reporte. Si el máximo no "
             "supera el umbral, o si no se encuentra la tabla ABSORPTION "
             "SPECTRUM en el .out, no se extrae ningún NTO para ese confórmero.",
    )
    parser.add_argument("--force", action="store_true",
                         help="regenerar aunque ya existan cubos")
    parser.add_argument("--dry-run", action="store_true",
                         help="solo lista qué haría, sin correr orca_plot")
    parser.add_argument("--timeout", type=int, default=600,
                         help="timeout por llamada a orca_plot en segundos (default: 600)")
    parser.add_argument(
        "--report", type=Path, default=None,
        help="ruta del CSV de salida (default: <root>/tddft_cubes_report.csv)",
    )
    parser.add_argument("-v", "--verbose", action="store_true",
                         help="logging DEBUG")
    return parser


def _main_cubes(argv: Optional[list[str]] = None) -> int:
    parser = _cubes_arg_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(levelname)s: %(message)s",
    )

    root = args.root.resolve()
    if not root.exists():
        _cubes_logger.error("La raíz de búsqueda no existe: %s", root)
        return 1

    report_path = args.report or (root / "tddft_cubes_report.csv")

    conformers = discover_conformers(root)
    if not conformers:
        _cubes_logger.warning("No se encontraron confórmeros TDDFT bajo %s", root)
        return 0

    _cubes_logger.info("Encontrados %d confórmeros TDDFT válidos", len(conformers))

    all_rows: list[ReportRow] = []
    n_errors = 0
    for conf_dir in conformers:
        try:
            rows = process_conformer(conf_dir, root, args)
            all_rows.extend(rows)
            n_errors += sum(1 for r in rows if r.status == "error")
        except Exception as exc:  # noqa: BLE001 - no queremos matar el batch
            _cubes_logger.exception("Fallo procesando %s", conf_dir)
            all_rows.append(
                ReportRow(
                    molecule="", run_label="", conf_id=conf_dir.name,
                    kind="conformer", homo_or_hole=None, lumo_or_particle=None,
                    coef=None, participation_ratio=None, fosc=None, cube_files="",
                    status="error", error=str(exc),
                )
            )
            n_errors += 1

    if not args.dry_run:
        write_report(all_rows, report_path)
        _cubes_logger.info("Reporte escrito en %s", report_path)
    else:
        _cubes_logger.info("[dry-run] no se escribió reporte; %d filas planificadas", len(all_rows))

    _cubes_logger.info(
        "Terminado: %d confórmeros, %d filas, %d errores",
        len(conformers), len(all_rows), n_errors,
    )
    return 1 if n_errors else 0


# ══════════════════════════════════════════════════════════════════════════════
#  §  CALIBRACIÓN SISTEMÁTICA DE pKa   (subcomando `calibrate`)
#     Plan de transiciones + cálculo de pKa crudos + ajuste por tipo
# ══════════════════════════════════════════════════════════════════════════════
#
# Estados (inferidos de la conectividad de los .xyz, no de los encabezados):
#   A  NH₃⁺, todos COOH (+1)      D  NH₃⁺, α-COO⁻, resto COOH (0, zwitterión)
#   R  NH₂, todos COOH (0)        B  Asp/Glu: NH₃⁺ + todos COO⁻ (−1)
#                                    Gly/Val: NH₂ + COO⁻ (−1)
#   C  Glu: NH₂ + 2×COO⁻ (−2)     B2 (AB2, GB2): NH₂, α-COO⁻, resto COOH (−1)
#                                    — tautómero de B, NO entra en la cadena.
# Cadena por molécula:  A →[R/D]→ B (→ C)
#   pKa0 = A→R/D     pKa1 = R/D→B     pKa2 = B→C
# AP (4-aminopiridina, verificada por conectividad) es un par aislado APA→APR.
#
# 'chain' = {rol: carpeta dentro de <base>/<MOL>/runs/}. Si una carpeta se
# llama distinto (p. ej. glycine en vez de GR) usa --override GLYCINE.neutral=glycine.
_CALIB_REF_AA = "McMurry Table 26-1 / Carey & Giuliano"
CALIBRATION_SETS: dict[str, dict] = {
    "ASPARTIC": {
        "label": "Aspartic",
        "chain": {"protonated": "AA", "neutral": "AR", "zwitterion": "AD", "anion": "AB"},
        "transitions": [
            {"slot": "pKa0", "type": "carboxyl",   "step": "A→R/D", "group": "α-COOH", "exp": 1.88},
            {"slot": "pKa1", "type": "carboxyl_2", "step": "R/D→B", "group": "β-COOH", "exp": 3.65},
        ],
    },
    "GLUTAMIC": {
        "label": "Glutamic",
        "chain": {"protonated": "GA", "neutral": "GR", "zwitterion": "GD",
                  "anion": "GB", "dianion": "GC"},
        "transitions": [
            {"slot": "pKa0", "type": "carboxyl",      "step": "A→R/D", "group": "α-COOH",  "exp": 2.19},
            {"slot": "pKa1", "type": "carboxyl_2",    "step": "R/D→B", "group": "γ-COOH",  "exp": 4.25},
            {"slot": "pKa2", "type": "amino_dianion", "step": "B→C",   "group": "NH₃⁺ (2×COO⁻)", "exp": 9.67},
        ],
    },
    "GLYCINE": {
        "label": "Glycine",
        "chain": {"protonated": "GA", "neutral": "GR", "zwitterion": "GD", "anion": "GB"},
        "transitions": [
            {"slot": "pKa0", "type": "carboxyl", "step": "A→R/D", "group": "COOH", "exp": 2.34},
            {"slot": "pKa1", "type": "amino",    "step": "R/D→B", "group": "NH₃⁺",  "exp": 9.60},
        ],
    },
    "VALINE": {
        "label": "Valine",
        "chain": {"protonated": "VA", "neutral": "VR", "zwitterion": "VD", "anion": "VB"},
        "transitions": [
            {"slot": "pKa0", "type": "carboxyl", "step": "A→R/D", "group": "COOH", "exp": 2.32},
            {"slot": "pKa1", "type": "amino",    "step": "R/D→B", "group": "NH₃⁺",  "exp": 9.62},
        ],
    },
    "AP": {
        "label": "AP",
        "pair": ("APA", "APR"),
        "transitions": [
            {"slot": "pair", "type": "pyridinium", "step": "APA→APR",
             "group": "N–H⁺ (4-aminopiridinio)", "exp": 9.17},
        ],
    },
}
_CALIB_ROLE_ORDER = ("protonated", "neutral", "zwitterion", "anion", "dianion")


def _calib_files(spec: dict, slot: str) -> tuple[str, str]:
    """Archivos/carpetas (ácido, base) que participan en una transición."""
    if "pair" in spec:
        return spec["pair"][0], spec["pair"][1]
    ch = spec["chain"]
    macro = f"{ch['neutral']}/{ch['zwitterion']}"
    if slot == "pKa0":
        return ch["protonated"], macro
    if slot == "pKa1":
        return macro, ch["anion"]
    if slot == "pKa2":
        return ch["anion"], ch["dianion"]
    raise KeyError(slot)


def _calib_apply_overrides(sets: dict, overrides: list[str]) -> dict:
    import copy
    sets = copy.deepcopy(sets)
    for item in overrides or []:
        try:
            lhs, folder = item.split("=", 1)
            mol, role = lhs.split(".", 1)
        except ValueError:
            raise SystemExit(f"--override: formato inválido '{item}' "
                             f"(se espera MOL.rol=carpeta, p. ej. GLYCINE.neutral=glycine)")
        mol, role = mol.upper(), role.lower()
        if mol not in sets or "chain" not in sets[mol] or role not in sets[mol]["chain"]:
            raise SystemExit(f"--override: '{mol}.{role}' no existe en CALIBRATION_SETS")
        sets[mol]["chain"][role] = folder.strip()
    return sets


def _calib_plan_rows(sets: dict, only: Optional[list[str]]) -> list[dict]:
    rows = []
    for mol, spec in sets.items():
        if only and mol not in only:
            continue
        for tr in spec["transitions"]:
            acid, base = _calib_files(spec, tr["slot"])
            rows.append({**tr, "mol": mol, "label": spec["label"],
                         "acid": acid, "base": base})
    return rows


def _calib_print_plan(rows: list[dict]) -> None:
    sep = "═" * 100
    print(f"\n{sep}\n  PLAN DE CALIBRACIÓN — transiciones, archivos y pKa experimental\n{sep}")
    hdr = (f"  {'Molécula':<9s} {'slot':<5s} {'tipo':<14s} {'paso':<8s} "
           f"{'grupo':<22s} {'ácido':<10s} {'base':<10s} {'exp':>5s}")
    print(hdr)
    print("  " + "─" * 96)
    last = None
    for r in rows:
        mol = r["mol"] if r["mol"] != last else ""
        last = r["mol"]
        print(f"  {mol:<9s} {r['slot']:<5s} {r['type']:<14s} {r['step']:<8s} "
              f"{r['group']:<22s} {r['acid']:<10s} {r['base']:<10s} {r['exp']:>5.2f}")
    print(sep)


def _calib_check_dirs(base: Path, sets: dict, only: Optional[list[str]]) -> list[str]:
    """Lista de problemas de disposición de carpetas (vacía = todo presente)."""
    probs: list[str] = []
    for mol, spec in sets.items():
        if only and mol not in only:
            continue
        runs = base / mol / "runs"
        if not runs.is_dir():
            probs.append(f"{mol}: no existe {runs}")
            continue
        names = list(spec["pair"]) if "pair" in spec else \
            [spec["chain"][r] for r in _CALIB_ROLE_ORDER if r in spec["chain"]]
        for nm in names:
            d = runs / nm
            if not d.is_dir():
                probs.append(f"{mol}: falta la carpeta de estado '{nm}' ({d})")
            elif not (d / SUBDIR_SP).is_dir():
                probs.append(f"{mol}: '{nm}' no tiene {SUBDIR_SP}/ (¿SP sin terminar?)")
    return probs


def _calib_failed(msg: str) -> dict:
    return {"ok": False, "pKa": None, "method": None, "warnings": [msg]}


def _calib_compute_mol(mol: str, spec: dict, base: Path, T: float,
                       method: str, gcorr_source: str) -> dict[str, dict]:
    """
    pKa crudos de todas las transiciones de UNA molécula. Usa exactamente
    las mismas piezas que `pka` (collect_mol_data, _build_neutral_macrostate,
    compute_pka), sin imprimir el ciclo. Retorna {slot: resultado_compute_pka}.
    """
    global PC_WORK_ROOT
    PC_WORK_ROOT = (base / mol / "runs").resolve()
    wanted = [t["slot"] for t in spec["transitions"]]
    out: dict[str, dict] = {}

    def load(name: str) -> dict:
        return collect_mol_data(name, verbose=False, gcorr_source=gcorr_source,
                                use_smd=False)

    if "pair" in spec:
        a, b = spec["pair"]
        da, db = load(a), load(b)
        if not (da.get("ok") and db.get("ok")):
            out["pair"] = _calib_failed(
                f"datos incompletos: {a} ok={da.get('ok')}, {b} ok={db.get('ok')}")
        else:
            out["pair"] = compute_pka(da, db, T=T, verbose=False, method=method,
                                      gcorr_source=gcorr_source)
        return out

    ch = spec["chain"]
    D = {role: load(name) for role, name in ch.items()}
    bad = [f"{ch[r]}" for r, d in D.items() if not d.get("ok")]

    macro = None
    if D["neutral"].get("ok") and D["zwitterion"].get("ok"):
        try:
            _taut, macro = _build_neutral_macrostate(
                D["neutral"], D["zwitterion"], ch["neutral"], ch["zwitterion"],
                T, verbose=False, gcorr_source=gcorr_source)
        except SystemExit:
            macro = None

    def need(*roles) -> Optional[str]:
        miss = [ch[r] for r in roles if r not in D or not D[r].get("ok")]
        return f"estados sin datos válidos: {', '.join(miss)}" if miss else None

    if "pKa0" in wanted:
        e = need("protonated") or (None if macro else "macroestado R/D no construible")
        out["pKa0"] = _calib_failed(e) if e else compute_pka(
            D["protonated"], macro, T=T, verbose=False, method=method,
            gcorr_source=gcorr_source)
    if "pKa1" in wanted:
        e = need("anion") or (None if macro else "macroestado R/D no construible")
        out["pKa1"] = _calib_failed(e) if e else compute_pka(
            macro, D["anion"], T=T, verbose=False, method=method,
            gcorr_source=gcorr_source)
    if "pKa2" in wanted:
        e = need("anion", "dianion")
        out["pKa2"] = _calib_failed(e) if e else compute_pka(
            D["anion"], D["dianion"], T=T, verbose=False, method=method,
            gcorr_source=gcorr_source)
    return out


def _calib_print_fits(types: list[str]) -> None:
    sep = "─" * 100
    print(f"\n{'═' * 100}\n  AJUSTE POR TIPO  (modo auto; offset salvo rango≥1.0 y n≥3)\n{'═' * 100}")
    for gsrc in ("harmonic", "goodvibes"):
        shown = False
        for t in types:
            fit = fit_empirical_correction(t, gcorr_source=gsrc, mode="auto")
            if not fit.get("ok"):
                continue
            if not shown:
                print(f"\n  gcorr_source = {gsrc}")
                print(f"  {'tipo':<14s} {'n':>2s} {'modo':<7s} {'m':>7s} {'b':>8s} "
                      f"{'RMSE':>7s} {'LOO':>7s}  calibradores")
                print("  " + sep)
                shown = True
            rm = "n/a" if np.isnan(fit["rmse"]) else f"{fit['rmse']:.3f}"
            lo = "n/a" if np.isnan(fit["loo_rmse"]) else f"{fit['loo_rmse']:.3f}"
            print(f"  {t:<14s} {fit['n']:>2d} {fit['mode']:<7s} {fit['slope']:>7.3f} "
                  f"{fit['intercept']:>+8.3f} {rm:>7s} {lo:>7s}  {', '.join(fit['molecules'])}")
            for w in fit["warnings"]:
                print(f"      ⚠  {w}")


def _main_calibrate():
    import copy, csv, datetime, json

    ap = argparse.ArgumentParser(
        prog="maaqc.py calibrate",
        description="Calibración sistemática del ciclo de pKa: imprime el plan de "
                    "transiciones (archivos y pKa exp.), calcula los pKa crudos "
                    "de Asp/Glu/Gly/Val/AP con la misma cadena que `pka`, ajusta "
                    "una corrección por tipo y escribe calibration.json "
                    "(consumible con `pka --calibration-file`).",
        epilog="""
Disposición esperada:  <base>/<MOL>/runs/<estado>/   (como en `pka --mol`)
  ASPARTIC/runs/{AA,AR,AD,AB}   GLUTAMIC/runs/{GA,GR,GD,GB,GC}
  GLYCINE/runs/{GA,GR,GD,GB}    VALINE/runs/{VA,VR,VD,VB}    AP/runs/{APA,APR}

Ejemplos:
  maaqc.py calibrate --list --base ~/PAPER1/calib
  maaqc.py calibrate --base ~/PAPER1/calib --method born --gcorr both
  maaqc.py calibrate --base ~/PAPER1/calib --override GLYCINE.neutral=glycine
""",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", default=".", metavar="DIR",
                    help="carpeta que contiene ASPARTIC/, GLUTAMIC/, GLYCINE/, "
                         "VALINE/ y AP/ (default: cwd)")
    ap.add_argument("--mols", nargs="+", choices=list(CALIBRATION_SETS), default=None,
                    help="subconjunto de moléculas (default: todas)")
    ap.add_argument("--list", action="store_true",
                    help="solo imprime el plan y verifica carpetas; no calcula")
    ap.add_argument("--method", choices=["auto", "born", "direct"], default="auto",
                    help="método de compute_pka (default auto). Usa EL MISMO que "
                         "usarás en M-Gly:Val; se guarda en el JSON.")
    ap.add_argument("--temp", type=float, default=T_STD, choices=TEMPS_K,
                    help=f"T de la calibración (default {T_STD}); los pKa exp. son a 25 °C")
    ap.add_argument("--gcorr", choices=["harmonic", "goodvibes", "both"],
                    default="harmonic", help="fuente de G_corr (default harmonic)")
    ap.add_argument("--override", nargs="*", default=[], metavar="MOL.rol=CARPETA",
                    help="renombra una carpeta de estado, p. ej. GLYCINE.neutral=glycine")
    ap.add_argument("--out", default="calibration.json", metavar="JSON")
    ap.add_argument("--csv", default=None, metavar="CSV",
                    help="además, tabla plana de transiciones en CSV")
    args = ap.parse_args()

    sets = _calib_apply_overrides(CALIBRATION_SETS, args.override)
    rows = _calib_plan_rows(sets, args.mols)
    _calib_print_plan(rows)

    base = Path(args.base).expanduser().resolve()
    probs = _calib_check_dirs(base, sets, args.mols)
    if probs:
        print("\n  ⚠  Problemas en la disposición de carpetas:")
        for p_ in probs:
            print(f"     · {p_}")
    else:
        print(f"\n  ✔ Todas las carpetas de estado existen bajo {base}")
    if args.list:
        return
    if probs:
        print("\n  Se calcularán solo las transiciones cuyos estados estén completos.")

    gsources = ["harmonic", "goodvibes"] if args.gcorr == "both" else [args.gcorr]
    # resultados[(mol, slot)][gsrc] = resultado compute_pka
    results: dict[tuple[str, str], dict[str, dict]] = {}
    for mol, spec in sets.items():
        if args.mols and mol not in args.mols:
            continue
        for gsrc in gsources:
            print(f"\n▸ {mol}  (gcorr_source={gsrc}, method={args.method}, T={args.temp:.2f} K)")
            try:
                res = _calib_compute_mol(mol, spec, base, args.temp, args.method, gsrc)
            except SystemExit as exc:
                res = {t["slot"]: _calib_failed(f"abortado (SystemExit {exc.code})")
                       for t in spec["transitions"]}
            for slot, r in res.items():
                results.setdefault((mol, slot), {})[gsrc] = r

    # ── tabla de resultados ────────────────────────────────────────────────
    methods: set[str] = set()
    out_rows: list[dict] = []
    sep = "═" * 100
    print(f"\n{sep}\n  pKa CRUDOS vs EXPERIMENTAL\n{sep}")
    print(f"  {'Molécula':<9s} {'slot':<5s} {'tipo':<14s} {'paso':<8s} "
          f"{'gcorr':<10s} {'método':<7s} {'calc':>7s} {'exp':>6s} {'calc−exp':>9s}")
    print("  " + "─" * 96)
    for r in rows:
        for gsrc in gsources:
            res = results.get((r["mol"], r["slot"]), {}).get(gsrc)
            ok = bool(res and res.get("ok"))
            calc = res["pKa"] if ok else None
            meth = res.get("method") if res else None
            if ok and meth:
                methods.add(meth)
            err = (calc - r["exp"]) if calc is not None else None
            cs = f"{calc:7.2f}" if calc is not None else "    n/a"
            es = f"{err:+9.2f}" if err is not None else "      n/a"
            print(f"  {r['mol']:<9s} {r['slot']:<5s} {r['type']:<14s} {r['step']:<8s} "
                  f"{gsrc:<10s} {str(meth or '-'):<7s} {cs} {r['exp']:6.2f} {es}")
            if not ok and res:
                for w in res.get("warnings", [])[:2]:
                    print(f"      ⚠  {w}")
            out_rows.append({
                "mol": r["mol"], "slot": r["slot"], "type": r["type"], "step": r["step"],
                "group": r["group"], "acid_files": r["acid"], "base_files": r["base"],
                "gcorr_source": gsrc, "method": meth, "pKa_calc": calc,
                "pKa_exp": r["exp"], "error": err, "ok": ok})

    if len(methods) > 1:
        print(f"\n  ⚠  MÉTODOS MEZCLADOS en la calibración: {sorted(methods)}. "
              f"G(H⁺,aq) difiere entre 'born' y 'direct': fuerza --method para "
              f"que todas las transiciones usen el mismo.")
    method_resolved = (next(iter(methods)) if len(methods) == 1
                       else ("mixed" if methods else args.method))

    # ── EMPIRICAL_CALIBRATION actualizada ──────────────────────────────────
    cal = copy.deepcopy(EMPIRICAL_CALIBRATION)
    for r in rows:
        pts = cal.setdefault(r["type"], [])
        pt = next((p_ for p_ in pts if p_["mol"] == r["label"]), None)
        if pt is None:
            pt = {"mol": r["label"], "calc_harmonic": None, "calc_goodvibes": None,
                  "exp": r["exp"], "exp_ref": _CALIB_REF_AA}
            pts.append(pt)
        pt["exp"] = r["exp"]
        for gsrc in gsources:
            res = results.get((r["mol"], r["slot"]), {}).get(gsrc)
            if res and res.get("ok"):
                pt[f"calc_{gsrc}"] = round(float(res["pKa"]), 3)
    # Instala en memoria para que los ajustes de abajo usen estos valores.
    EMPIRICAL_CALIBRATION.clear()
    EMPIRICAL_CALIBRATION.update(cal)
    _calib_print_fits(list(cal))

    blob = {
        "meta": {"created": datetime.datetime.now().isoformat(timespec="seconds"),
                 "T": args.temp, "method_requested": args.method,
                 "method_resolved": method_resolved, "gcorr": gsources,
                 "base": str(base), "script": "maaqc.py calibrate",
                 "note": "pKa exp. a 25 °C; ver exp_ref por punto"},
        "points": out_rows,
        "EMPIRICAL_CALIBRATION": cal,
    }
    outp = Path(args.out).expanduser()
    outp.write_text(json.dumps(blob, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n  → {outp.resolve()}")
    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(out_rows[0]))
            w.writeheader()
            w.writerows(out_rows)
        print(f"  → {Path(args.csv).resolve()}")
    print("\n  Uso en M-Gly:Val:\n"
          f"    maaqc.py pka --mol MGLYVAL --protonated A --neutral R --zwitterion D \\\n"
          f"        --anion B --dianion C --method {method_resolved if method_resolved != 'mixed' else '<born|direct>'} \\\n"
          f"        --calibration-file {outp} \\\n"
          f"        --empirical-correction pKa0=carboxyl,pKa1=carboxyl_2,pKa2=pyridinium \\\n"
          f"        --empirical-correction-mode offset\n")


# ══════════════════════════════════════════════════════════════════════════════
#  §  CLI UNIFICADO
#     Despacha a los main() originales conservando sus flags intactos
# ══════════════════════════════════════════════════════════════════════════════

_SUBCOMMANDS = {
    "run":     ("etapas de cómputo CREST→CENSO→OPT→NumFreq→SP→TDDFT→ADC2",
                "ex pipeline.py"),
    "status":  ("estado de todos los runs (atajo de `run --status`)",
                "ex pipeline.py --status"),
    "collect": ("análisis de resultados: Boltzmann multi-T, GoodVibes, ESP, ADC2",
                "ex results_collector.py"),
    "pka":     ("ciclo termodinámico de pKa, cadena de 5 estados y Monte Carlo",
                "ex protonation_cycle.py"),
    "cubes":   ("cubos HOMO/LUMO y NTO hole/particle vía orca_plot",
                "ex extract_tddft_cubes.py"),
    "calibrate": ("calibración sistemática de pKa (Asp/Glu/Gly/Val/AP) → calibration.json",
                  "nuevo"),
}

_USAGE = """\
maaqc.py — MAA Quantum Chemistry, herramienta única

uso: maaqc.py <subcomando> [opciones...]

subcomandos:
""" + "\n".join(
    f"  {name:<9s} {desc}\n{'':11s}({origin})"
    for name, (desc, origin) in _SUBCOMMANDS.items()
) + """

Cada subcomando conserva exactamente los flags de su script de origen:

  maaqc.py run --mol molecules/mglyval.xyz --charge -1 --cores 8
  maaqc.py run --mol molecules/mglyval.xyz --start opt --skipnumfreq
  maaqc.py status
  maaqc.py collect --mol MGLYVAL --state D --export xlsx --plot
  maaqc.py collect --mol MGLYVAL --compare-states --adc2
  maaqc.py pka --mol MGLYVAL --neutral R --anion B --dianion C --plot
  maaqc.py pka --mol MGLYVAL --neutral R --anion B --dianion C \\
               --protonated A --zwitterion D --montecarlo
  maaqc.py cubes --root PAPER1 --grid 80 --dry-run
  maaqc.py calibrate --list --base PAPER1/calib
  maaqc.py calibrate --base PAPER1/calib --method born --gcorr both
  maaqc.py pka --mol AP --pair APA APR --method born
  maaqc.py pka --mol MGLYVAL --protonated A --neutral R --zwitterion D --anion B \\
               --dianion C --calibration-file calibration.json \\
               --empirical-correction pKa0=carboxyl,pKa1=carboxyl_2,pKa2=pyridinium

Ayuda de cada uno:  maaqc.py <subcomando> --help
"""


def main(argv: "list[str] | None" = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)

    if not argv or argv[0] in ("-h", "--help", "help"):
        print(_USAGE)
        return 0

    sub, rest = argv[0], argv[1:]

    if sub not in _SUBCOMMANDS:
        print(f"maaqc.py: subcomando desconocido: {sub!r}\n", file=sys.stderr)
        print(_USAGE, file=sys.stderr)
        return 2

    # Los main() originales leen sys.argv directamente (pipeline lo necesita
    # para partir la línea de comandos por --mol repetido), así que se
    # reescribe prog + args y se restaura al salir.
    saved_argv = sys.argv
    try:
        if sub == "status":
            sys.argv = ["maaqc.py status", "--status", *rest]
            _main_run()
            return 0
        if sub == "run":
            sys.argv = ["maaqc.py run", *rest]
            _main_run()
            return 0
        if sub == "collect":
            sys.argv = ["maaqc.py collect", *rest]
            _main_collect()
            return 0
        if sub == "pka":
            sys.argv = ["maaqc.py pka", *rest]
            _main_pka()
            return 0
        if sub == "cubes":
            sys.argv = ["maaqc.py cubes", *rest]
            return _main_cubes(rest)
        if sub == "calibrate":
            sys.argv = ["maaqc.py calibrate", *rest]
            _main_calibrate()
            return 0
    except SystemExit as exc:
        return exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
    finally:
        sys.argv = saved_argv
    return 0


if __name__ == "__main__":
    sys.exit(main())
