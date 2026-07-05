"""ResourceProfile — detección topology-aware de CPU (P-core/E-core) y
distribución de cores/memoria entre stages.

Migrar desde pipeline.py:
  - class ResourceProfile        líneas 414-650
      .detect()                  líneas 503-531
      .detect_with_batch()       líneas 532-554
      .summary()                 líneas 555-585
      ._compute()                líneas 586-650
  - build_omp_env()              líneas 197-233 (incluye el workaround de
                                  OMPI_MCA_hwloc_base_binding_policy para la
                                  CPU híbrida del ThinkPad T14, ver comentario
                                  extenso en el original — NO simplificar,
                                  ese fix es producto de debugging real)
  - assign_cpu_lists()           líneas 3259-3351
  - _check_cpu_budget()          líneas 3352-3370

Nota: esta es de las piezas más autocontenidas del legacy (poca
dependencia cruzada con parsing de ORCA/CENSO) — buena candidata para
migrar temprano, en paralelo a stages/opt.py, sin bloquear el resto.
"""

from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass
class ResourceProfile:
    """Ver pipeline.py líneas 414-650 para los campos reales (topología
    P-core/E-core, nprocs/maxcore derivados, batch sizes)."""

    # TODO: portar campos de __init__ (líneas 483-502) y _compute().

    @classmethod
    def detect(cls, batch_size: int = 1) -> "ResourceProfile":
        raise NotImplementedError("Migrar desde pipeline.py:ResourceProfile.detect (503-531)")

    @classmethod
    def detect_with_batch(cls, total_cores: int, batch_size: int, ram_mb: int = 0) -> "ResourceProfile":
        raise NotImplementedError("Migrar desde pipeline.py:ResourceProfile.detect_with_batch (532-554)")


def build_omp_env(nprocs: int) -> dict:
    """Construye el entorno OMP/MKL/OpenBLAS + workarounds OMPI_MCA para
    la CPU híbrida. Migrar tal cual desde pipeline.py líneas 197-233,
    incluyendo el comentario que explica el fix de hwloc/taskset — es
    contexto necesario para no revertir el fix por accidente."""
    raise NotImplementedError("Migrar desde pipeline.py:build_omp_env (líneas 197-233)")


def assign_cpu_lists(batch_size: int, nprocs_per_job: int) -> list[list[int]]:
    raise NotImplementedError("Migrar desde pipeline.py:assign_cpu_lists (líneas 3259-3351)")
