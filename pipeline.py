#!/usr/bin/env python3
"""
pipeline.py — Pipeline ORCA automatizado para cálculos conformacionales
Flujo: CREST → CENSO (prescreen) → Opt → NumFreq → SP → TDDFT → ADC2

El pipeline siempre empieza en CREST. CENSO pre-screening es obligatorio
entre CREST y OPT. El número de confórmeros que pasan a OPT se determina
automáticamente a partir del panorama conformacional post-CENSO (mín 3, máx 20).

Uso (modo normal):
    python pipeline.py --mol molecules/mol_A.xyz --charge -1
    python pipeline.py --mol molecules/mol_A.xyz --charge 0 --start opt
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --start sp
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --skipnumfreq   # OPT → SP → TDDFT directo
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --nroot         # SP+NRoots1 top-5 por E_OPT
    python pipeline.py --mol molecules/mol_A.xyz --charge -1 --nroot --nroot-top 3  # top-3
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

Uso (modo expansión conformacional, --expand-crest):
    Diseñado para casos como el estado B donde la inversión energética GFN2→DFT
    exige explorar más confórmeros más allá de los primeros N ya calculados.
    El flujo es: extrae confs CREST[skip:skip+window] → OPT en paralelo (lotes
    de 3 × 4 cores) → filtro ΔE(DFT) < threshold kcal/mol → NUMFREQ → SP → TDDFT.

    python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest
    python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest \\
        --crest-skip 3 --crest-window 17 --dft-threshold 3.0 --parallel 3
    python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest \\
        --start numfreq   # salta expand-OPT, retoma desde NUMFREQ
"""

import argparse
import datetime
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

ORCA_BIN   = Path("/home/kmilozeta/orca_6_1_1_avx2/orca")
CREST_BIN  = Path("crest")          # asumido en PATH; ajusta si hace falta
CENSO_BIN  = Path("censo")          # asumido en PATH (env sciencebitch); ajusta si hace falta
WORK_ROOT  = Path("runs")           # directorio de salida
TMPL_DIR   = Path("templates")
MULT       = 1                      # multiplicidad siempre 1
CREST_TIMEOUT_HOURS = 24.0          # timeout CREST (horas)
CENSO_TIMEOUT_HOURS = 24.0          # timeout CENSO pre-screening (horas)

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

# OMP_ENV se construye bajo demanda en build_omp_env() para que siempre
# refleje el valor actualizado de NPROCS (nunca queda obsoleto si el
# usuario pasa --nprocs desde CLI).
def build_omp_env(nprocs: int) -> dict:
    return {
        **os.environ,
        "OMP_NUM_THREADS": str(nprocs),
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        # ── Fix para CPU híbrida P-core/E-core + taskset ──────────────
        # Cuando taskset restringe el proceso a un subconjunto de cores,
        # el mpirun interno de ORCA (orca_leanscf_mpi en la fase LEANSCF)
        # invoca hwloc para detectar topología, y hwloc lee /sys/.../cluster
        # directamente — sin filtrar por el cpuset heredado de taskset.
        # En esta CPU híbrida eso produce cpusets de cluster L1d que no son
        # subconjunto coherente del cpuset restringido ("intersection
        # without inclusion"), y aunque hwloc dice que ignora el error y
        # continúa, el reparto de trabajo entre ranks MPI en LEANSCF queda
        # corrupto y aborta justo al iniciar la formación de Fock.
        # Forzar el binding por núcleo en vez de por topología (que es lo
        # que realmente queremos al usar taskset) evita que hwloc intente
        # razonar sobre clusters L1d/L2 inconsistentes.
        "OMPI_MCA_hwloc_base_binding_policy": "none",
        "OMPI_MCA_rmaps_base_mapping_policy": "core",
        "OMPI_MCA_hwloc_base_use_hwthreads_as_cpus": "1",
    }

# Alias global mantenido para compatibilidad con run_orca(env=None)
OMP_ENV: dict = build_omp_env(NPROCS)

STEPS = ["crest", "opt", "numfreq", "sp", "tddft", "adc2"]

# ─────────────────────────────────────────────
#  GESTIÓN DE SEÑALES / CTRL-C
# ─────────────────────────────────────────────

_active_procs: list[subprocess.Popen] = []   # procesos ORCA activos
_procs_lock   = threading.Lock()

def _kill_active_procs() -> None:
    """Mata todos los subprocesos ORCA/taskset registrados en _active_procs."""
    with _procs_lock:
        procs = list(_active_procs)
    for p in procs:
        try:
            p.terminate()
        except OSError:
            pass
    # Espera breve y SIGKILL si siguen vivos
    for p in procs:
        try:
            p.wait(timeout=3)
        except subprocess.TimeoutExpired:
            try:
                p.kill()
            except OSError:
                pass

def _sigint_handler(signum, frame) -> None:  # noqa: ANN001
    print("\n[pipeline] Ctrl-C recibido — terminando procesos ORCA activos…",
          file=sys.stderr)
    _kill_active_procs()
    sys.exit(130)   # código estándar para Ctrl-C (128 + SIGINT)

signal.signal(signal.SIGINT,  _sigint_handler)
signal.signal(signal.SIGTERM, _sigint_handler)

# ─────────────────────────────────────────────
#  CONSTANTES FÍSICAS / MODO EXPAND
# ─────────────────────────────────────────────

HARTREE_TO_KCAL     = 627.5094740631   # factor de conversión
EXPAND_BATCH_SIZE   = 3                # OPTs en paralelo (default; override con --total-cores)
EXPAND_DFT_THRESHOLD = 3.0             # kcal/mol (default filtro post-OPT)


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
        MAX_OPT_NPROCS_PER_JOB — cores/job en OPT (modo normal y expand)
        OPT_BATCH_SIZE         — jobs paralelos en OPT modo normal
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
            batch_size   : jobs simultáneos en modo expand (default 1 = modo normal).
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
        Variante para modo expand: recibe total_cores explícito (de --cores)
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

        if self.batch_size > 1:
            expand_info = (
                f"\n  Expand batch    : {self.batch_size} jobs × {self.nprocs_opt} cores "
                f"= {self.batch_size * self.nprocs_opt} cores/lote"
                f"\n                    RAM/lote ≈ {self.batch_size * ram_opt:.1f} GB"
            )
        else:
            expand_info = ""

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
            f"→ {ram_adc2:.1f} GB{expand_info}"
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
        # Modo normal: siempre OPT_BATCH_SIZE jobs × MAX_OPT_NPROCS_PER_JOB cores.
        # Modo expand (p>1): batch_size viene de --parallel; mismo límite por job.
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
        # Modo normal: siempre SP_BATCH_SIZE jobs × MAX_SP_NPROCS cores.
        # Modo expand (p>1): mismo esquema que OPT.
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
        log.info(f"  ▶  lanzando ORCA: {inp_file.name}  (cores {cpu_str})")
    else:
        log.info(f"  ▶  lanzando ORCA: {inp_file.name}")

    t0       = time.time()
    timeout  = timeout_hours * 3600
    proc_env = env if env is not None else build_omp_env(NPROCS)

    try:
        with open(out_file, "w") as fout:
            proc = subprocess.Popen(
                cmd,
                stdout=fout,
                stderr=subprocess.STDOUT,
                cwd=inp_file.parent,
                env=proc_env,
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
                log: logging.Logger) -> Optional[Path]:
    """
    Etapa 01_crest: corre CREST sobre la molécula de entrada.
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

    log.info(f"[CREST] lanzando sobre {mol_xyz.name}  (carga={charge})")
    cmd = [
        str(CREST_BIN),
        str(input_xyz.name),
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
                cmd,
                stdout=flog,
                stderr=subprocess.STDOUT,
                cwd=crest_dir,
                env=build_omp_env(NPROCS_CREST),
                timeout=CREST_TIMEOUT_HOURS * 3600,
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

    # CREST puede escribir el archivo con distintos nombres
    for candidate in ["crest_conformers.xyz", "crest_best.xyz"]:
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
        log.info(f"[CREST] {n} confórmeros generados en {elapsed:.1f} min")
        return conformers_xyz
    else:
        log.error(f"[CREST] no se encontró conformers.xyz "
                  f"(returncode={returncode}) — revisa {crest_dir / 'crest.log'}")
        return None


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
        mol_dir        : directorio raíz del run (runs/<mol_name>/).
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

    output_xyz = censo_dir / "censo_ensemble_sorted.xyz"

    # Skip si ya existe y no está vacío
    if output_xyz.exists() and output_xyz.stat().st_size > 0:
        log.info(f"[CENSO]  ↩  pre-screening ya completado — saltando ({output_xyz.name})")
        return output_xyz

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

    cmd = [
        str(CENSO_BIN),
        "--prescreening",
        "-i", local_input.name,
        "-c", str(charge),
        "-u", "0",                  # MULT global siempre 1 → singlete, capa cerrada
        "--solvent", "water",
        # prog/func/basis se leen de ~/.censo2rc ([prescreening] prog=xtb, func=gfn2, smgsolv=alpb)
        # --prog fue eliminado en CENSO 2.x; ya no se pasa por CLI
        "--maxcores", str(NPROCS_CREST),
        "-T", f"{temperature:.2f}",             # temperatura configurable (default 274.15 K = 1°C Antártico)
    ]

    log.info(f"[CENSO]  ▶  lanzando pre-screening sobre {conformers_xyz.name} "
             f"(carga={charge}, T={temperature:.2f} K)")
    t0 = time.time()
    try:
        with open(censo_dir / "censo.log", "w") as flog:
            proc = subprocess.run(
                cmd, stdout=flog, stderr=subprocess.STDOUT,
                cwd=censo_dir, env=build_omp_env(NPROCS_CREST),
                timeout=CENSO_TIMEOUT_HOURS * 3600,
            )
        returncode = proc.returncode
    except subprocess.TimeoutExpired:
        log.error(f"[CENSO]  ✗  TIMEOUT ({CENSO_TIMEOUT_HOURS} h) — "
                  f"usando ensemble CREST sin filtrar")
        return conformers_xyz
    except FileNotFoundError:
        log.error(f"[CENSO]  ✗  no se encontró el ejecutable ({CENSO_BIN}) — "
                  f"¿'conda activate sciencebitch'? Usando ensemble CREST sin filtrar")
        return conformers_xyz
    except Exception as e:
        log.error(f"[CENSO]  ✗  error inesperado: {e} — usando ensemble CREST sin filtrar")
        return conformers_xyz

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

    # ── Búsqueda robusta de output: nombres fijos + glob de respaldo ─────────
    # CENSO 2/3 puede escribir en censo_dir directamente o en subdirectorios
    # (CENSO/, prescreening/, etc.) con nombres que varían entre versiones.
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


def stage_opt(conformers_xyz: Path, mol_dir: Path,
              charge: int, nconf: int, log: logging.Logger,
              conf_ids: Optional[list[str]] = None) -> list[Path]:
    """
    Etapa 02_opt: optimización DFT de cada confórmero.
    `nconf` viene de select_nconf_auto() o select_nconf_boltzmann().

    Args:
        conf_ids : lista explícita de CONF# a optimizar (modo --boltzmann-threshold).
                   Si se especifica, se ignora nconf y se seleccionan los frames
                   correspondientes del ensemble. Si es None, se toman los primeros
                   nconf frames (comportamiento clásico).

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
    else:
        # Modo clásico: primeros nconf frames
        log.info(f"[OPT] optimizando {nconf} confórmeros (seleccionados automáticamente)")
        try:
            conformers = read_conformers(conformers_xyz, nconf)
        except RuntimeError as e:
            log.error(f"[OPT] no se pudo leer {conformers_xyz}: {e}")
            return []

    if not conformers:
        log.error("[OPT] el archivo de confórmeros está vacío")
        return []

    log.info(f"[OPT] {len(conformers)} confórmeros disponibles — "
             f"batch {OPT_BATCH_SIZE}×{NPROCS} cores")
    template   = get_template("opt")
    successful = []

    tasks_all = [
        {
            "conf_id":  f"conf_{i:03d}",
            "conf_dir": mol_dir / "02_opt" / f"conf_{i:03d}",
            "atoms":    atoms,
            "charge":   charge,
            "template": template,
            "nprocs":   NPROCS,
            "maxcore":  MAXCORE,
            "log":      log,
        }
        for i, atoms in enumerate(conformers, start=1)
    ]

    n_lotes = -(-len(tasks_all) // OPT_BATCH_SIZE)   # ceiling division
    _check_cpu_budget(OPT_BATCH_SIZE, NPROCS, log)

    for b_start in range(0, len(tasks_all), OPT_BATCH_SIZE):
        batch    = tasks_all[b_start : b_start + OPT_BATCH_SIZE]
        lote_num = b_start // OPT_BATCH_SIZE + 1
        ids      = [t["conf_id"] for t in batch]

        cpu_lists = assign_cpu_lists(len(batch), NPROCS)
        for t, cl in zip(batch, cpu_lists):
            t["cpu_list"] = cl

        log.info(f"[OPT] ── Lote {lote_num}/{n_lotes}: {ids}  "
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
                sys.exit(130)

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
    batch_size>0  → valor explícito (usado por --expand-crest con --parallel N).
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

    for b_start in range(0, len(tasks_all), _batch):
        batch    = tasks_all[b_start : b_start + _batch]
        lote_num = b_start // _batch + 1
        ids      = [t["conf_id"] for t in batch]

        cpu_lists = assign_cpu_lists(len(batch), _nprocs)
        for t, cl in zip(batch, cpu_lists):
            t["cpu_list"] = cl

        log.info(f"[SP] ── Lote {lote_num}/{n_lotes}: {ids}  "
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
                sys.exit(130)

    log.info(f"[SP] {len(successful)}/{len(freq_dirs)} completados")
    return successful


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
    worker_env = build_omp_env(nprocs)
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

    n_lotes = -(-len(tasks_all) // TDDFT_BATCH_SIZE)
    _check_cpu_budget(TDDFT_BATCH_SIZE, NPROCS_TDDFT, log)

    for b_start in range(0, len(tasks_all), TDDFT_BATCH_SIZE):
        batch    = tasks_all[b_start : b_start + TDDFT_BATCH_SIZE]
        lote_num = b_start // TDDFT_BATCH_SIZE + 1
        ids      = [t["conf_id"] for t in batch]

        cpu_lists = assign_cpu_lists(len(batch), NPROCS_TDDFT)
        for t, cl in zip(batch, cpu_lists):
            t["cpu_list"] = cl

        log.info(f"[TDDFT] ── Lote {lote_num}/{n_lotes}: {ids}  "
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
                sys.exit(130)

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
#  ETAPA ADC2  (séptima etapa, obligatoria)
# ─────────────────────────────────────────────

def stage_adc2(mol_dir: Path, charge: int, log: logging.Logger) -> Optional[Path]:
    """
    Etapa 07_adc2: validación ADC2 (NRoots 1) sobre un único confórmero —
    el dominante, definido como el de menor energía SP (FINAL SINGLE POINT
    ENERGY en 04_sp/conf_*/sp.out).

    La geometría usada es la misma que TDDFT: sp_opt.xyz copiada en
    05_tddft/<conf_id>/ (idéntica a opt_opt.xyz de 02_opt/numfreq, pero se
    reutiliza la copia ya validada en esa etapa para evitar resolver de nuevo
    el archivo de origen).

    Requiere que 04_sp/ y 05_tddft/<conf_id>/ ya existan para el confórmero
    elegido. Si 05_tddft no tiene la copia (p.ej. corrida con --start previo
    a tddft incompleta), aborta esa etapa con error explicado.

    Devuelve el directorio 07_adc2/<conf_id>/ si el cálculo terminó bien,
    o None si falló o no se pudo seleccionar confórmero.
    """
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

    # ── Localizar geometría: misma que usó TDDFT (sp_opt.xyz) ────
    tddft_dir = mol_dir / "05_tddft" / conf_id
    geom_src  = tddft_dir / "sp_opt.xyz"
    if not geom_src.exists():
        log.error(
            f"[ADC2] No se encontró {geom_src} — el confórmero dominante "
            f"({conf_id}) necesita haber pasado por TDDFT primero (la geometría "
            f"de ADC2 se toma de ahí). Corre la etapa TDDFT para este confórmero."
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
        return adc2_dir

    log.error(f"[ADC2] {conf_id} falló.")
    return None


# ─────────────────────────────────────────────
#  EXPANSIÓN CONFORMACIONAL (--expand-crest)
# ─────────────────────────────────────────────

def find_opt_minimum(mol_dir: Path, log: logging.Logger) -> Optional[float]:
    """
    Escanea runs/<mol>/02_opt/conf_*/opt.out para encontrar la energía DFT
    mínima entre los confórmeros ya optimizados.

    Se usa como referencia E_ref en el filtro de expansión; si no hay opts
    previos exitosos devuelve None (el mínimo se determinará sobre el lote nuevo).
    """
    opt_root = mol_dir / "02_opt"
    if not opt_root.exists():
        return None

    min_e  : Optional[float] = None
    min_id : Optional[str]   = None

    for conf_dir in sorted(opt_root.glob("conf_*")):
        out_file = conf_dir / "opt.out"
        if not check_success(out_file):
            continue
        e = extract_sp_energy(out_file)
        if e is not None and (min_e is None or e < min_e):
            min_e  = e
            min_id = conf_dir.name

    if min_e is not None:
        log.info(f"[EXPAND] Mínimo DFT de opts previos: {min_id}  "
                 f"E = {min_e:.10f} Ha")
    else:
        log.info("[EXPAND] Sin opts previos con energía válida — "
                 "se usará el mínimo dinámico del lote actual")
    return min_e


def _opt_worker(task: dict) -> dict:
    """
    Worker para un único OPT, diseñado para ejecutarse en un hilo de
    ThreadPoolExecutor.  ORCA corre como subproceso externo (libera el GIL),
    por lo que el paralelismo real es efectivo.

    task debe contener:
        conf_id, conf_dir, atoms, charge, template, nprocs, maxcore, log
    Opcional:
        cpu_list: list[int] de CPUs lógicas exclusivas para este worker,
                  asignado por el orquestador del lote (stage_expand_opt)
                  vía assign_cpu_lists() para evitar que OpenMPI distribuya
                  mal los procesos entre workers concurrentes.

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
    worker_env = build_omp_env(nprocs)
    success = run_orca(inp_file, out_file, log, env=worker_env, cpu_list=cpu_list)
    elapsed = (time.time() - t0) / 60

    energy: Optional[float] = None
    if success:
        if resolve_opt_xyz(conf_dir) is None:
            log.warning(f"[EXPAND] {conf_id}: ORCA OK pero sin geometría optimizada")
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
    worker_env = build_omp_env(nprocs)
    success    = run_orca(inp_file, out_file, log, env=worker_env, cpu_list=cpu_list)
    elapsed    = (time.time() - t0) / 60

    return {
        "conf_id": conf_id,
        "sp_dir":  sp_dir if success else None,
        "success": success,
        "elapsed": elapsed,
    }


def assign_cpu_lists(batch_size: int, nprocs_per_job: int) -> list[list[int]]:
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
      - == 2: preferentemente 1 P-core completo (2 hilos SMT); si no
        alcanzan los P-cores, 2 E-cores individuales (no comparten core
        físico entre sí, pero tampoco se parte ningún core entre jobs).
      - == 1: 1 E-core (preferido) o, si no alcanzan, 1 hilo SMT de P-core.
      - otro valor (ej. 6 para NumFreq): reparto secuencial sobre el pool
        completo en el mismo orden de prioridad.

    Devuelve una lista de listas de CPUs, una por worker, en el mismo
    orden que `batch`. Si no hay suficientes CPUs físicas para satisfacer
    batch_size workers, la lista devuelta es más corta que batch_size
    (el caller debe manejarlo o confiar en _check_cpu_budget para advertir
    antes de llegar aquí).
    """
    P_CORE_BLOCKS   = [[0, 1], [2, 3]]                 # CORE 0, CORE 1 (SMT)
    E_CORE_SINGLES  = [4, 5, 6, 7, 8, 9, 10, 11]        # CORE 2-9
    LP_CORE_SINGLES = [12, 13]                          # CORE 10-11 (lentos)

    lists: list[list[int]] = []

    if nprocs_per_job == 2:
        for block in P_CORE_BLOCKS:
            if len(lists) >= batch_size:
                break
            lists.append(list(block))

        e_pool = list(E_CORE_SINGLES)
        while len(lists) < batch_size and len(e_pool) >= 2:
            lists.append([e_pool.pop(0), e_pool.pop(0)])

        lp_pool = list(LP_CORE_SINGLES)
        while len(lists) < batch_size and len(lp_pool) >= 2:
            lists.append([lp_pool.pop(0), lp_pool.pop(0)])

    elif nprocs_per_job == 1:
        pool = list(E_CORE_SINGLES) + list(LP_CORE_SINGLES) + \
               [c for block in P_CORE_BLOCKS for c in block]
        for cpu in pool:
            if len(lists) >= batch_size:
                break
            lists.append([cpu])

    else:
        full_pool = [c for block in P_CORE_BLOCKS for c in block] + \
                    E_CORE_SINGLES + LP_CORE_SINGLES
        idx = 0
        for _ in range(batch_size):
            chunk = full_pool[idx: idx + nprocs_per_job]
            if len(chunk) < nprocs_per_job:
                break
            lists.append(chunk)
            idx += nprocs_per_job

    return lists


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


def stage_expand_opt(
    conformers_xyz : Path,
    mol_dir        : Path,
    charge         : int,
    skip           : int,
    window         : int,
    dft_threshold  : float,
    batch_size     : int,
    nprocs_per_job : int,
    maxcore_per_job: int,
    log            : logging.Logger,
) -> list[Path]:
    """
    Expansión conformacional desde CREST con OPT en paralelo y filtro ΔE DFT.

    Flujo detallado:
      1. Lee CREST conformers en el rango [skip : skip+window] (índices 0-based).
      2. Optimiza en lotes de `batch_size` con ThreadPoolExecutor.
         Cada OPT usa `nprocs_per_job` cores  →  batch_size × nprocs_per_job cores
         en paralelo. Los lotes son secuenciales; dentro de cada lote los OPTs
         corren simultáneamente.
      3. Determina E_ref = mínimo DFT global (opts previos ∪ lote actual).
      4. Descarta confórmeros con ΔE = E - E_ref > dft_threshold kcal/mol.
      5. Escribe trazabilidad en mol_dir/expand_discarded.log.

    Args:
        conformers_xyz  : ruta a 01_crest/crest_conformers.xyz
        mol_dir         : directorio raíz runs/<mol>/
        charge          : carga formal
        skip            : confórmeros CREST ya procesados (se omiten)
        window          : cuántos confórmeros adicionales intentar
        dft_threshold   : umbral ΔE en kcal/mol (recomendado 3.0)
        batch_size      : OPTs simultáneos por lote
        nprocs_per_job  : cores por job OPT (viene de --nprocs o calculado por --cores)
        maxcore_per_job : MB/core por job OPT
        log             : logger del pipeline principal

    Returns:
        Lista de directorios 02_opt/conf_XXX/ que pasan el filtro,
        listos para continuar en stage_numfreq.
    """
    log.info("=" * 60)
    log.info(f"[EXPAND] skip={skip}  window={window}  "
             f"threshold={dft_threshold} kcal/mol  "
             f"batch={batch_size} × {nprocs_per_job} cores = "
             f"{batch_size * nprocs_per_job} cores/lote")

    # ── Leer el multi-XYZ de CREST ──────────────────────────────
    try:
        all_confs = read(str(conformers_xyz), index=":")
        if not isinstance(all_confs, list):
            all_confs = [all_confs]
    except Exception as exc:
        log.error(f"[EXPAND] No se pudo leer {conformers_xyz}: {exc}")
        return []

    total = len(all_confs)
    log.info(f"[EXPAND] Total confórmeros en CREST: {total}")

    end_idx   = min(skip + window, total)
    new_confs = all_confs[skip:end_idx]

    if not new_confs:
        log.warning(f"[EXPAND] Rango [{skip}:{end_idx}] vacío — nada que procesar")
        return []

    log.info(f"[EXPAND] Procesando CREST[{skip}:{end_idx}] → "
             f"conf_{skip+1:03d} … conf_{end_idx:03d}  "
             f"({len(new_confs)} estructuras)")

    # ── Referencia energética de opts previos ────────────────────
    e_ref = find_opt_minimum(mol_dir, log)

    # ── Armar lista de tareas ────────────────────────────────────
    template  = get_template("opt")
    _check_cpu_budget(batch_size, nprocs_per_job, log)
    tasks_all = [
        {
            "conf_id":  f"conf_{skip + rel + 1:03d}",
            "conf_dir": mol_dir / "02_opt" / f"conf_{skip + rel + 1:03d}",
            "atoms":    atoms,
            "charge":   charge,
            "template": template,
            "nprocs":   nprocs_per_job,
            "maxcore":  maxcore_per_job,
            "log":      log,
        }
        for rel, atoms in enumerate(new_confs)
    ]

    n_lotes = -(-len(tasks_all) // batch_size)   # ceiling division
    log.info(f"[EXPAND] {len(tasks_all)} OPTs en {n_lotes} lote(s) de ≤{batch_size} "
             f"(CPUs disponibles: {os.cpu_count()})")
    log.info(f"[EXPAND] Cores por OPT: {nprocs_per_job}  "
             f"Cores en paralelo por lote: {batch_size * nprocs_per_job}  "
             f"RAM por lote ≈ {batch_size * nprocs_per_job * maxcore_per_job / 1024:.1f} GB")

    # ── Lotes paralelos ──────────────────────────────────────────
    opt_results: list[dict] = []

    for b_start in range(0, len(tasks_all), batch_size):
        batch    = tasks_all[b_start : b_start + batch_size]
        lote_num = b_start // batch_size + 1
        ids      = [t["conf_id"] for t in batch]

        # Cores exclusivos por worker dentro de ESTE lote (no solapan entre
        # sí; se reasignan en cada lote ya que los lotes son secuenciales).
        cpu_lists = assign_cpu_lists(len(batch), nprocs_per_job)
        for t, cl in zip(batch, cpu_lists):
            t["cpu_list"] = cl

        log.info(f"[EXPAND] ── Lote {lote_num}/{n_lotes}: {ids}  "
                 f"(cores: {[','.join(map(str, cl)) for cl in cpu_lists]})")

        with ThreadPoolExecutor(max_workers=len(batch)) as pool:
            futures = {pool.submit(_opt_worker, t): t["conf_id"] for t in batch}
            try:
                for fut in as_completed(futures):
                    res = fut.result()
                    opt_results.append(res)
                    tag = "✓" if res["success"] else "✗"
                    e_s = f"{res['energy']:.10f} Ha" if res["energy"] is not None else "—"
                    log.info(f"[EXPAND]   {tag} {res['conf_id']}  "
                             f"E = {e_s}  ({res['elapsed']:.1f} min)")
            except KeyboardInterrupt:
                log.warning("[EXPAND] Ctrl-C — cancelando lote, esperando workers…")
                pool.shutdown(wait=False, cancel_futures=True)
                _kill_active_procs()
                sys.exit(130)

    # ── Actualizar E_ref con el mínimo del lote actual ───────────
    ok_new = [r for r in opt_results if r["success"] and r["energy"] is not None]

    if ok_new:
        e_min_new = min(r["energy"] for r in ok_new)
        if e_ref is None or e_min_new < e_ref:
            e_ref = e_min_new
            log.info(f"[EXPAND] Nuevo mínimo DFT global (lote actual): "
                     f"{e_ref:.10f} Ha")

    if e_ref is None:
        log.error("[EXPAND] No se pudo determinar E_ref — imposible aplicar filtro")
        return []

    # ── Filtro ΔE < dft_threshold ────────────────────────────────
    log.info(f"[EXPAND] Filtro ΔE < {dft_threshold} kcal/mol  "
             f"(E_ref = {e_ref:.10f} Ha)")

    passing  : list[Path] = []
    discarded: list[str]  = []

    for r in sorted(ok_new, key=lambda x: x["energy"]):  # type: ignore[return-value]
        dE = (r["energy"] - e_ref) * HARTREE_TO_KCAL     # type: ignore[operator]
        if dE <= dft_threshold:
            log.info(f"[EXPAND]   ✓  {r['conf_id']}  ΔE = {dE:+.4f} kcal/mol → PASA")
            passing.append(r["conf_dir"])
        else:
            log.info(f"[EXPAND]   ✗  {r['conf_id']}  ΔE = {dE:+.4f} kcal/mol → DESCARTADO")
            discarded.append(r["conf_id"])

    failed = [r["conf_id"] for r in opt_results if not r["success"]]

    log.info(f"[EXPAND] ─────────────────────────────────────────────")
    log.info(f"[EXPAND]   OPTs intentados : {len(opt_results)}")
    log.info(f"[EXPAND]   Fallidos OPT     : {len(failed)}  "
             f"{'  '.join(failed) if failed else '—'}")
    log.info(f"[EXPAND]   Pasan filtro     : {len(passing)}")
    log.info(f"[EXPAND]   Descartados      : {len(discarded)}  "
             f"{'  '.join(discarded) if discarded else '—'}")
    log.info(f"[EXPAND] ─────────────────────────────────────────────")

    # ── Trazabilidad en disco (texto legible) ─────────────────────
    discard_log = mol_dir / "expand_discarded.log"
    e_ref_str = f"{e_ref:.10f} Ha" if e_ref is not None else "N/A"
    with open(discard_log, "a") as flog:
        flog.write(
            f"\n# {datetime.datetime.now().isoformat()}"
            f"  threshold={dft_threshold} kcal/mol"
            f"  E_ref={e_ref_str}\n"
        )
        for r in sorted(opt_results, key=lambda x: x["conf_id"]):
            if r["energy"] is not None and e_ref is not None:
                dE      = (r["energy"] - e_ref) * HARTREE_TO_KCAL
                verdict = "PASA" if dE <= dft_threshold else "DESCARTADO"
            else:
                dE      = float("nan")
                verdict = "FALLO_OPT"
            flog.write(
                f"  {r['conf_id']}  E={r['energy']}  "
                f"dE={dE:.4f} kcal/mol  {verdict}\n"
            )

    # ── Trazabilidad en disco (JSON estructurado, fuente de verdad) ─
    # Esto es lo que lee la recuperación de --start posteriores a opt,
    # para no tener que adivinar ni re-derivar qué confs pasaron el
    # filtro ΔE a partir de qué carpetas quedaron físicamente en 02_opt/.
    filter_result_path = mol_dir / "expand_filter_result.json"
    filter_result = {
        "timestamp":     datetime.datetime.now().isoformat(),
        "skip":          skip,
        "window":        window,
        "dft_threshold": dft_threshold,
        "e_ref_ha":      e_ref,
        "passing":       sorted(p.name for p in passing),
        "discarded":     sorted(discarded),
        "failed_opt":    sorted(failed),
    }
    with open(filter_result_path, "w") as fjson:
        json.dump(filter_result, fjson, indent=2)
    log.info(f"[EXPAND] Resultado del filtro guardado en {filter_result_path}")

    return passing


# ─────────────────────────────────────────────
#  STATUS GLOBAL
# ─────────────────────────────────────────────

def print_status():
    """Imprime estado de todos los cálculos en runs/."""
    runs = WORK_ROOT
    if not runs.exists():
        print("No se encontró el directorio 'runs/'. Aún no hay cálculos.")
        return

    print(f"\n{'─'*60}")
    print(f"  Estado de runs en: {runs.resolve()}")
    print(f"{'─'*60}")

    for mol_dir in sorted(runs.iterdir()):
        if not mol_dir.is_dir():
            continue
        print(f"\n  📁 {mol_dir.name}")
        for step_num in ["01_crest", "02_opt", "03_numfreq", "04_sp", "05_tddft", "06_sp_nroot", "07_adc2"]:
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

def main():
    global NPROCS_CREST, NPROCS, MAXCORE, NPROCS_NUMFREQ, MAXCORE_NUMFREQ, \
           NPROCS_TDDFT, MAXCORE_TDDFT, NPROCS_ADC2, MAXCORE_ADC2, \
           OPT_BATCH_SIZE, SP_BATCH_SIZE, TDDFT_BATCH_SIZE, OMP_ENV
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

  # Ver estado de todos los runs
  python pipeline.py --status

Ejemplos (modo expansión conformacional):
  # Expande B: 12 cores totales → 3 OPTs simultáneos × 4 cores cada uno
  python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest --cores 12

  # Control explícito: 2 OPTs simultáneos × 6 cores c/u (= 12 cores)
  python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest \\
      --cores 12 --parallel 2

  # Ventana y threshold personalizados
  python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest \\
      --cores 12 --crest-skip 3 --crest-window 17 --dft-threshold 3.0

  # Retomar desde NUMFREQ (expand-OPT ya hecho)
  python pipeline.py --mol molecules/mol_B.xyz --charge -1 --expand-crest \\
      --cores 8 --start numfreq
        """,
    )
    parser.add_argument("--mol",    type=Path, help="ruta al .xyz de entrada")
    parser.add_argument("--charge", type=int, default=0,
                        help="carga formal de la molécula (default: 0)")
    parser.add_argument("--start",  choices=STEPS, default="crest",
                        help="etapa desde la que iniciar (default: crest)")
    parser.add_argument("--status", action="store_true",
                        help="mostrar estado de todos los runs y salir")
    parser.add_argument(
        "--skipnumfreq", action="store_true",
        help=(
            "Salta NUMFREQ y pasa directo de OPT a SP/TDDFT (modo normal, "
            "no aplica a --expand-crest). Usa la geometría de 02_opt/ "
            "directamente como input de SP, sin correcciones térmicas de NumFreq."
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
            "automáticamente:\n"
            "  • Modo normal  → todos los cores van a cada job secuencial\n"
            "                   (OPT / NumFreq / SP / TDDFT)\n"
            "  • Modo expand  → se reparten entre --parallel jobs en paralelo.\n"
            "                   Con --parallel P: nprocs_por_job = cores // P\n"
            "                   (el resto de cores se ignora para no saturar RAM).\n"
            "  maxcore se recalcula como (RAM_total_MB × 0.85) // cores,\n"
            "  con RAM_total_MB = 28000 por defecto (ajustar en MAX_RAM_MB).\n"
            "  Ejemplo: --cores 12 --parallel 3  →  4 cores × 3 OPTs paralelos."
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

    # ── Modo expansión conformacional ──────────────────────────────────────
    parser.add_argument(
        "--expand-crest", action="store_true",
        help=(
            "Modo expansión: extrae confórmeros adicionales de CREST "
            "(más allá de los ya procesados), los optimiza en paralelo "
            "y continúa los que pasan el filtro ΔE por NUMFREQ → SP → TDDFT. "
            "Requiere que CREST ya haya completado (01_crest/ debe existir)."
        ),
    )
    parser.add_argument(
        "--crest-skip", type=int, default=3, metavar="N",
        help="confórmeros CREST ya procesados a omitir (default: 3)",
    )
    parser.add_argument(
        "--crest-window", type=int, default=17, metavar="N",
        help="confórmeros adicionales a intentar (default: 17, i.e. conf_4…20)",
    )
    parser.add_argument(
        "--dft-threshold", type=float, default=EXPAND_DFT_THRESHOLD, metavar="X",
        help=f"umbral ΔE(DFT) en kcal/mol para el filtro post-OPT "
             f"(default: {EXPAND_DFT_THRESHOLD})",
    )
    parser.add_argument(
        "--parallel", type=int, default=EXPAND_BATCH_SIZE, metavar="N",
        help=(
            f"Número de cálculos simultáneos por lote en modo expand (default: {EXPAND_BATCH_SIZE}). "
            f"Con --cores C, cada proceso usará C // N cores."
        ),
    )

    # ── Modo SP + NRoots 1 ─────────────────────────────────────────────────
    parser.add_argument(
        "--nroot", action="store_true",
        help=(
            "Modo SP+NRoots 1: calcula un SP con TD-DFT mínimo (NRoots 1) "
            "sobre los mejores confórmeros según energía DFT de OPT (02_opt/). "
            "Los inputs se escriben en 06_sp_nroot/<conf_id>/sp_nroot.inp. "
            "Requiere que OPT haya completado. Incompatible con --expand-crest."
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
            "--trange es opcional (default: 274 298 24). "
            "Compatible con modo normal y --expand-crest."
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

    args = parser.parse_args()

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

    # Construir lista de temperaturas
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

    if args.skipnumfreq and args.expand_crest:
        sys.exit(
            "ERROR: --skipnumfreq no aplica en modo --expand-crest "
            "(usa --start numfreq/sp/tddft para saltar etapas ahí)."
        )
    if args.skipnumfreq and args.start in ("numfreq", "sp", "tddft"):
        sys.exit(
            f"ERROR: --skipnumfreq es incompatible con --start {args.start}. "
            f"--skipnumfreq ya implica saltar NUMFREQ desde OPT; usa solo "
            f"--start opt (o --start crest) junto con --skipnumfreq."
        )
    if args.nroot and args.expand_crest:
        sys.exit(
            "ERROR: --nroot es incompatible con --expand-crest. "
            "Corre --nroot como paso independiente sobre un run ya completado."
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
            batch_size  = args.parallel if args.expand_crest else 1,
        )
    else:
        rp = ResourceProfile.detect(
            batch_size = args.parallel if args.expand_crest else 1,
        )

    # Extraer valores del perfil autodetectado
    nprocs_crest          = rp.nprocs_crest
    nprocs                = rp.nprocs_opt
    maxcore               = rp.maxcore_opt
    nprocs_numfreq        = rp.nprocs_numfreq
    maxcore_nf            = rp.maxcore_numfreq
    nprocs_per_expand_job = rp.nprocs_opt     # alias para log posterior

    # Sobreescritura por flags manuales (mayor precedencia)
    if args.nprocs is not None:
        nprocs = args.nprocs
        nprocs_per_expand_job = args.nprocs
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

    if args.status:
        print_status()
        return

    if args.mol is None:
        parser.error("--mol es obligatorio salvo que uses --status")

    mol_xyz = args.mol.resolve()
    if not mol_xyz.exists():
        sys.exit(f"ERROR: no se encontró {mol_xyz}")

    if not mol_xyz.suffix.lower() == ".xyz":
        sys.exit(f"ERROR: el archivo de entrada debe ser un .xyz: {mol_xyz}")

    # Validar binarios antes de empezar
    if not ORCA_BIN.exists():
        sys.exit(f"ERROR: no se encontró el ejecutable de ORCA en: {ORCA_BIN}\n"
                 f"       Ajusta ORCA_BIN al inicio del script.")

    mol_name = mol_xyz.stem
    mol_dir  = WORK_ROOT / mol_name
    mol_dir.mkdir(parents=True, exist_ok=True)

    log = setup_logging(mol_dir)
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
    if args.expand_crest:
        total_expand_cores = args.parallel * nprocs_per_expand_job
        total_expand_ram   = total_expand_cores * maxcore / 1024
        log.info(f"Expand batch           : {args.parallel} × {nprocs_per_expand_job} = "
                 f"{total_expand_cores} cores en paralelo / RAM estimada ≈ {total_expand_ram:.1f} GB")
    log.info("=" * 60)

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
    #  MODO EXPANSIÓN CONFORMACIONAL  (--expand-crest)
    # ════════════════════════════════════════════════════════════
    if args.expand_crest:
        log.info("MODO: Expansión conformacional (--expand-crest)")
        log.info(f"skip={args.crest_skip}  window={args.crest_window}  "
                 f"threshold={args.dft_threshold} kcal/mol  "
                 f"paralelo={args.parallel}")

        conformers_xyz = mol_dir / "01_crest" / "crest_conformers.xyz"
        if not conformers_xyz.exists():
            sys.exit(
                f"ERROR: no existe {conformers_xyz}\n"
                f"       --expand-crest requiere que CREST haya completado "
                f"previamente (corre primero sin --expand-crest)."
            )

        # ── expand-OPT: lotes paralelos + filtro ────────────────
        if start_idx <= STEPS.index("opt"):
            # Obtenemos el nombre de la etapa actual para decidir qué función llamar
            current_stage = STEPS[start_idx]

            # Si estamos en expand-OPT, usamos la lógica de batching
            if current_stage == "opt":
                output_dirs = stage_expand_opt(
                    conformers_xyz=conformers_xyz,
                    mol_dir=mol_dir,
                    charge=args.charge,
                    skip=args.crest_skip,
                    window=args.crest_window,
                    dft_threshold=args.dft_threshold,
                    batch_size=args.parallel,
                    nprocs_per_job=nprocs_per_expand_job,
                    maxcore_per_job=maxcore,
                    log=log,
                )

                if not output_dirs:
                    log.error(f"[EXPAND-OPT] Ningún confórmero nuevo pasó el filtro — abortando.")
                    sys.exit(1)
                expand_opt_dirs = output_dirs

            else:
                # Recuperar SOLO los confs que pasaron el filtro ΔE.
                # Fuente de verdad: expand_filter_result.json, escrito por
                # stage_expand_opt() en la corrida que hizo el OPT+filtro.
                # NO se reconstruye por índice de carpeta — eso ignoraría el
                # filtro ΔE y arrastraría confs ya DESCARTADOs (bug histórico).
                opt_root = mol_dir / "02_opt"
                if not opt_root.exists():
                    sys.exit(
                        f"ERROR: no existe {opt_root}\n"
                        f"       Con --start {args.start} se esperan OPTs previos "
                        f"en 02_opt/."
                    )

                filter_result_path = mol_dir / "expand_filter_result.json"
                if not filter_result_path.exists():
                    sys.exit(
                        f"ERROR: no existe {filter_result_path}\n"
                        f"       --expand-crest --start {args.start} necesita el "
                        f"resultado del filtro ΔE de una corrida previa de "
                        f"stage_expand_opt (etapa OPT). Sin ese archivo no hay "
                        f"forma segura de saber qué confórmeros pasaron el "
                        f"filtro vs. cuáles fueron DESCARTADOs.\n"
                        f"       Corre primero con --start opt (o sin --start) "
                        f"para regenerar {filter_result_path.name}, o, si lo "
                        f"perdiste pero conservas el log, reconstrúyelo a mano "
                        f"a partir de expand_discarded.log."
                    )

                with open(filter_result_path) as fjson:
                    filter_result = json.load(fjson)

                passing_names = filter_result.get("passing", [])
                if not passing_names:
                    sys.exit(
                        f"[EXPAND] {filter_result_path} no registra ningún "
                        f"confórmero PASA — nada que hacer."
                    )

                expand_opt_dirs = sorted(
                    opt_root / name for name in passing_names
                    if (opt_root / name).is_dir()
                )
                missing = [n for n in passing_names if not (opt_root / n).is_dir()]
                if missing:
                    log.warning(
                        f"[EXPAND] {len(missing)} conf(s) marcados PASA en "
                        f"{filter_result_path.name} pero sin carpeta en 02_opt/: "
                        f"{'  '.join(missing)}"
                    )
                if not expand_opt_dirs:
                    sys.exit(
                        f"[EXPAND] Ningún directorio de 02_opt/ coincide con los "
                        f"confs PASA registrados en {filter_result_path.name}."
                    )

                log.info(
                    f"[EXPAND] OPT saltado — {len(expand_opt_dirs)} dir(s) "
                    f"recuperados de 02_opt/ según filtro ΔE registrado "
                    f"({filter_result_path.name}): "
                    f"{'  '.join(d.name for d in expand_opt_dirs)}"
                )
        else:
            # Saltamos OPT directamente (start_idx > opt): cargar resultados del filtro
            opt_root = mol_dir / "02_opt"
            if not opt_root.exists():
                sys.exit(
                    f"ERROR: no existe {opt_root}\n"
                    f"       Con --start {args.start} se esperan OPTs previos en 02_opt/."
                )

            filter_result_path = mol_dir / "expand_filter_result.json"
            if not filter_result_path.exists():
                sys.exit(
                    f"ERROR: no existe {filter_result_path}\n"
                    f"       --expand-crest --start {args.start} necesita el resultado "
                    f"del filtro ΔE de una corrida previa. Corre primero con --start opt."
                )

            with open(filter_result_path) as fjson:
                filter_result = json.load(fjson)

            passing_names = filter_result.get("passing", [])
            if not passing_names:
                sys.exit(
                    f"[EXPAND] {filter_result_path} no registra ningún confórmero PASA."
                )

            expand_opt_dirs = sorted(
                opt_root / name for name in passing_names
                if (opt_root / name).is_dir()
            )
            if not expand_opt_dirs:
                sys.exit(
                    f"[EXPAND] Ningún directorio de 02_opt/ coincide con los confs PASA."
                )

            log.info(
                f"[EXPAND] OPT saltado — {len(expand_opt_dirs)} dir(s) recuperados "
                f"según filtro ΔE: {'  '.join(d.name for d in expand_opt_dirs)}"
            )

        # ── NUMFREQ (solo nuevos confs) ──────────────────────────
        if start_idx <= STEPS.index("numfreq"):
            freq_dirs = stage_numfreq(expand_opt_dirs, mol_dir, args.charge, log)
            if not freq_dirs:
                log.error("[EXPAND] Ningún NUMFREQ completado. Abortando.")
                sys.exit(1)
        else:
            opt_names = {d.name for d in expand_opt_dirs}
            freq_root = mol_dir / "03_numfreq"
            freq_dirs = (
                sorted(d for d in freq_root.glob("conf_*") if d.name in opt_names)
                if freq_root.exists() else []
            )
            log.info(f"[EXPAND] NUMFREQ saltado — "
                     f"{len(freq_dirs)} dir(s) recuperados de 03_numfreq/")

        # ── SP ───────────────────────────────────────────────────
        if start_idx <= STEPS.index("sp"):
            sp_dirs = stage_sp(freq_dirs, mol_dir, args.charge, log,
                               batch_size=args.parallel,
                               nprocs_per_job=nprocs_per_expand_job,
                               maxcore_per_job=maxcore)
            if not sp_dirs:
                log.error("[EXPAND] Ningún SP completado. Abortando.")
                sys.exit(1)
        else:
            freq_names = {d.name for d in freq_dirs}
            sp_root = mol_dir / "04_sp"
            sp_dirs = (
                sorted(d for d in sp_root.glob("conf_*") if d.name in freq_names)
                if sp_root.exists() else []
            )
            log.info(f"[EXPAND] SP saltado — "
                     f"{len(sp_dirs)} dir(s) recuperados de 04_sp/")

        # ── TDDFT ────────────────────────────────────────────────
        if start_idx <= STEPS.index("tddft"):
            tddft_dirs = stage_tddft(sp_dirs, mol_dir, args.charge, log)
            if not tddft_dirs:
                log.warning("[EXPAND] Ningún TDDFT completado.")
            else:
                log.info(f"[EXPAND] {len(tddft_dirs)} espectros TD-DFT calculados ✓")
        else:
            log.info("[EXPAND] TDDFT no incluido en esta ejecución")

        # ── ADC2 (obligatoria, confórmero dominante por E_SP) ────
        if start_idx <= STEPS.index("adc2"):
            adc2_dir = stage_adc2(mol_dir, args.charge, log)
            if adc2_dir is None:
                log.warning("[EXPAND] No se pudo completar la validación ADC2.")
            else:
                log.info(f"[EXPAND] ADC2 completado ✓ → {adc2_dir}")
        else:
            log.info("[EXPAND] ADC2 no incluido en esta ejecución")

        log.info("=" * 60)
        log.info("PIPELINE (EXPAND) FINALIZADO")
        log.info(f"Resultados en: {mol_dir.resolve()}")
        log.info("=" * 60)
        return   # no continuar al flujo normal

    # ════════════════════════════════════════════════════════════
    #  MODO NORMAL
    # ════════════════════════════════════════════════════════════

    # ── Etapa 1: CREST ──────────────────────────────────────────
    # CREST siempre corre cuando start_idx == 0 (default).
    # Si se usa --start opt/numfreq/sp/tddft/adc2, se asume que
    # 01_crest/crest_conformers.xyz ya existe.
    conformers_xyz = mol_dir / "01_crest" / "crest_conformers.xyz"

    if start_idx == 0:
        conformers_xyz = stage_crest(mol_xyz, mol_dir, args.charge, log)
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
        log.info("[CREST] saltado (--start indica etapa posterior)")

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
                charge=args.charge,
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
                conformers_xyz, mol_dir, args.charge, log,
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

        opt_dirs = stage_opt(
            conformers_xyz, mol_dir, args.charge, nconf, log,
            conf_ids=boltz_conf_ids if boltz_conf_ids else None,
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
        freq_dirs = stage_numfreq(opt_dirs, mol_dir, args.charge, log)
        if not freq_dirs:
            log.error("Ningún cálculo de frecuencias completado. Abortando.")
            sys.exit(1)
    else:
        freq_dirs = sorted((mol_dir / "03_numfreq").glob("conf_*")) \
                    if (mol_dir / "03_numfreq").exists() else []
        log.info(f"[NUMFREQ] saltado — encontrados {len(freq_dirs)} directorios previos")

    # ── Etapa 4: SP ─────────────────────────────────────────────
    if start_idx <= STEPS.index("sp"):
        sp_dirs = stage_sp(freq_dirs, mol_dir, args.charge, log)
        if not sp_dirs:
            log.error("Ningún SP completado. Abortando.")
            sys.exit(1)
    else:
        sp_dirs = sorted((mol_dir / "04_sp").glob("conf_*")) \
                  if (mol_dir / "04_sp").exists() else []
        log.info(f"[SP] saltado — encontrados {len(sp_dirs)} directorios previos")

    # ── Etapa 5: TDDFT ──────────────────────────────────────────
    if start_idx <= STEPS.index("tddft"):
        tddft_dirs = stage_tddft(sp_dirs, mol_dir, args.charge, log)
        if not tddft_dirs:
            log.warning("Ningún TDDFT completado.")
        else:
            log.info(f"[TDDFT] {len(tddft_dirs)} espectros calculados ✓")
    else:
        log.info("[TDDFT] no incluido en esta ejecución")

    # ── Etapa 6: ADC2 (obligatoria, confórmero dominante por E_SP) ──
    if start_idx <= STEPS.index("adc2"):
        adc2_dir = stage_adc2(mol_dir, args.charge, log)
        if adc2_dir is None:
            log.warning("[ADC2] No se pudo completar la validación ADC2.")
        else:
            log.info(f"[ADC2] validación completada ✓ → {adc2_dir}")
    else:
        log.info("[ADC2] no incluido en esta ejecución")

    log.info("=" * 60)
    log.info("PIPELINE FINALIZADO")
    log.info(f"Resultados en: {mol_dir.resolve()}")
    log.info("=" * 60)


if __name__ == "__main__":
    main()
