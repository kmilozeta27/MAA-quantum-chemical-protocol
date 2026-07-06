"""ResourceProfile — detección topology-aware de CPU (P-core/E-core) y
distribución de cores/memoria entre stages.

`build_omp_env()` y `assign_cpu_lists()` ya NO viven acá: fueron migradas
1:1 desde pipeline.py directamente a `maaqc.stages._orca_shared`, que es
la fuente única de verdad y de donde las importan opt.py/sp.py/tddft.py/
crest.py/_censo_shared.py. La verificación de presupuesto de cores
también quedó ahí, como `check_cpu_budget()` (pública, sin underscore).

`ResourceProfile` (esta clase) es la última pieza pendiente de Issue B —
migrada 1:1 desde pipeline.py:414-650. Mismos nombres, mismas constantes
de política, mismo `_compute()`. A propósito queda AISLADA: no está
cableada a `config/schema.py` (campo `cores_total`) ni a `cli.py`
(no existe `--cores` todavía). Cablearla implica tocar opt.py/sp.py/
tddft.py/crest.py/numfreq.py a la vez para que consuman el profile en
vez de los valores estáticos de `RecipeConfig.resources` — eso queda
para un issue de integración aparte, una vez que esta clase tenga tests
propios y kmilozeta decida explícitamente activar el cableado.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field


@dataclass
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

    Nota de migración: en pipeline.py legacy esta era una clase "normal"
    (no dataclass) con __init__ explícito que llamaba a self._compute()
    al final. Acá se preserva ese mismo flujo (post_init dispara
    _compute()) en vez de convertir los campos derivados en
    dataclass fields, porque son ~14 atributos calculados y no tiene
    sentido declararlos como parámetros del constructor -- seguirían
    siendo sobreescritos por _compute() en cada instancia.
    """

    # ── Parámetros de política (ajustar aquí si cambia el hardware) ──────
    SAFETY_FACTOR: float = field(default=0.85, init=False, repr=False)
    MAX_CREST_NPROCS: int = field(default=12, init=False, repr=False)
    MAX_NUMFREQ_NPROCS: int = field(default=12, init=False, repr=False)
    MAX_OPT_NPROCS_PER_JOB: int = field(default=4, init=False, repr=False)
    OPT_BATCH_SIZE: int = field(default=3, init=False, repr=False)
    MAX_SP_NPROCS: int = field(default=4, init=False, repr=False)
    SP_BATCH_SIZE: int = field(default=3, init=False, repr=False)
    TDDFT_NPROCS_PER_JOB: int = field(default=6, init=False, repr=False)
    TDDFT_BATCH_SIZE: int = field(default=2, init=False, repr=False)
    RAM_PENALTY_TDDFT: float = field(default=1.4, init=False, repr=False)
    NPROCS_ADC2_MAX: int = field(default=4, init=False, repr=False)
    MAXCORE_ADC2_MIN: int = field(default=5000, init=False, repr=False)

    # ── Parámetros de instancia (pipeline.py:483-497, __init__) ──────────
    total_cores: int = 4
    ram_mb: int = 28_000
    batch_size: int = 1

    def __post_init__(self) -> None:
        """Equivalente al __init__ legacy: guarda los parámetros crudos
        y dispara _compute() para poblar los atributos derivados
        (nprocs_*/maxcore_* por stage)."""
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
            except ImportError:
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
        # iMTD-GC independiente → escala casi lineal hasta saturar BW de
        # memoria. Se usa MAX_CREST_NPROCS para priorizar throughput sin
        # desperdiciar cores que beneficiarían a ORCA (que escala peor).
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
