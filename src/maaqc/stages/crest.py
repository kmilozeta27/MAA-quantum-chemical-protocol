"""CrestStage — búsqueda conformacional GFN2-xTB.

Migrado desde pipeline.py:stage_crest() (líneas 1152-1233 en el legacy).
No tiene stages previos (requires = []).

⚠️ NOTA DE FIDELIDAD: este módulo se escribió sin tener a la vista el
`pipeline.py` legacy real (no estaba disponible al momento de portar) —
se reconstruyó a partir de: (1) el resumen de comportamiento dado en el
issue de migración, (2) los patrones ya establecidos en `opt.py` y
`_censo_shared.py` (uso de `safe_copy`/`build_omp_env`/`USE_NOHUP`/
`_apply_nohup`/`resolve_bin`/`kill_active_procs`), y (3) convenciones
de CREST 3.x ya documentadas en la skill de química computacional.
Antes de mergear, comparar contra pipeline.py líneas 1152-1233 línea a
línea — en particular el orden exacto de flags y qué pasa con stdout vs.
un log file dedicado (acá asumí patrón `censo.log`-like, i.e. un
`crest.log` con stdout+stderr combinados, igual que hace
`run_censo_prescreen` en `_censo_shared.py`).

Puntos a preservar al migrar (recordatorio del stub original):
  - Uso de CREST_BIN resuelto vía _resolve_bin() -> system/platform_detect.py.
  - Timeout CREST_TIMEOUT_HOURS (24h default).
  - Energías embebidas en el multi-XYZ de salida se parsean después con
    io/parsing.py::parse_crest_energies (NO reimplementar el parseo acá).

Preguntas abiertas (marcadas también en el issue, sin resolver todavía):
  1. Origen del .xyz de entrada: acá se asume convención `workdir /
     "input.xyz"`, escrita por el orchestrator antes de invocar
     CrestStage (simétrico a cómo OptStage busca su ensemble "hacia
     atrás" en el workdir en vez de recibir una ruta como parámetro).
     Confirmar con kmilo cómo el orchestrator resuelve el input inicial
     de la cadena para los distintos `--entry` (smiles/xyz/...).
  2. Flag `-T`: se asume que es *threads* (uso estándar de CREST 3.x,
     análogo a `--maxcores` en CENSO), NO temperatura — de ahí que se
     mapea a `config.resources.nprocs_crest`, no a `config.censo.temperature`.
     Si el legacy lo usaba para temperatura, esto está mal y hay que
     corregirlo contra pipeline.py real.
"""

from __future__ import annotations

import logging
import subprocess
import time
from pathlib import Path
from typing import Optional

from maaqc.config.schema import RecipeConfig
from maaqc.system._orca_shared import (
    USE_NOHUP,
    _apply_nohup,
    build_omp_env,
    kill_active_procs,
    safe_copy,
)
from maaqc.system.platform_detect import resolve_bin
from maaqc.stages.base import Stage, StageResult

log = logging.getLogger("maaqc.stages.crest")


class CrestStage(Stage):
    name = "crest"
    requires: list[str] = []

    STAGE_DIR = "01_crest"
    INPUT_FILENAME = "input.xyz"          # ← supuesto, ver nota de arriba
    LOCAL_INPUT_FILENAME = "crest_input.xyz"
    OUTPUT_FILENAME = "crest_conformers.xyz"
    # CREST puede nombrar la salida de cualquiera de estas dos formas
    # según versión/flags (--nci, ensembles chicos colapsan a un único
    # mínimo, etc.) — mismo fallback que el legacy (líneas 1211-1217).
    _RESULT_CANDIDATES = ["crest_conformers.xyz", "crest_best.xyz"]

    # ── contrato Stage ─────────────────────────────────────────────────

    def is_complete(self, workdir: Path) -> bool:
        out = workdir / self.STAGE_DIR / self.OUTPUT_FILENAME
        return out.exists() and out.stat().st_size > 0

    def inputs(self, workdir: Path) -> list[Path]:
        return [workdir / self.INPUT_FILENAME]

    def outputs(self, workdir: Path) -> list[Path]:
        return [workdir / self.STAGE_DIR / self.OUTPUT_FILENAME]

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        if config.charge is None:
            return StageResult(ok=False, message="config.charge no está definido")

        mol_xyz = workdir / self.INPUT_FILENAME
        if not mol_xyz.exists() or mol_xyz.stat().st_size == 0:
            return StageResult(
                ok=False, message=f"no se encontró geometría de entrada en {mol_xyz}",
            )

        crest_dir = workdir / self.STAGE_DIR
        crest_dir.mkdir(parents=True, exist_ok=True)
        output_xyz = crest_dir / self.OUTPUT_FILENAME

        if self.is_complete(workdir):
            log.info(f"[CREST] {output_xyz} ya existe — skip")
            return StageResult(ok=True, message="ya completado (skip)", data={"output": str(output_xyz)})

        local_input = crest_dir / self.LOCAL_INPUT_FILENAME
        if not safe_copy(mol_xyz, local_input):
            return StageResult(ok=False, message=f"no se pudo copiar {mol_xyz.name} a {crest_dir}")

        charge = config.charge
        # multiplicity 1 -> 0 electrones desapareados (uhf); si algún día
        # se corren radicales, esto necesita generalizarse a
        # (multiplicity - 1) en vez de hardcodear 0.
        uhf = config.multiplicity - 1
        nprocs_crest = config.resources.nprocs_crest
        timeout_hours = config.crest.timeout_hours
        alpb_solvent = config.crest.alpb_solvent

        crest_bin = resolve_bin("crest")
        cmd = [
            str(crest_bin), local_input.name,
            "--gfn2",
            "--alpb", alpb_solvent,
            "-T", str(nprocs_crest),
            "--chrg", str(charge),
            "--uhf", str(uhf),
        ]

        log.info(
            f"[CREST]  ▶  lanzando búsqueda conformacional sobre {mol_xyz.name} "
            f"(carga={charge}, uhf={uhf}, alpb={alpb_solvent}, T={nprocs_crest})"
        )
        t0 = time.time()
        log_path = crest_dir / "crest.log"
        try:
            with open(log_path, "w") as flog:
                proc = subprocess.run(
                    _apply_nohup(cmd), stdout=flog, stderr=subprocess.STDOUT,
                    cwd=crest_dir, env=build_omp_env(nprocs_crest),
                    timeout=timeout_hours * 3600, start_new_session=USE_NOHUP,
                )
            returncode = proc.returncode
        except subprocess.TimeoutExpired:
            log.error(f"[CREST]  ✗  TIMEOUT ({timeout_hours} h) sobre {mol_xyz.name}")
            kill_active_procs()
            return StageResult(ok=False, message=f"timeout tras {timeout_hours}h")
        except FileNotFoundError:
            log.error(f"[CREST]  ✗  no se encontró el ejecutable ({crest_bin})")
            return StageResult(ok=False, message=f"ejecutable CREST no encontrado ({crest_bin})")
        except KeyboardInterrupt:
            log.warning("[CREST]  Ctrl-C — cancelando, matando procesos huérfanos…")
            kill_active_procs()
            raise
        except Exception as e:
            log.error(f"[CREST]  ✗  error inesperado: {e}")
            return StageResult(ok=False, message=f"error inesperado: {e}")

        elapsed = (time.time() - t0) / 60

        # Búsqueda del resultado bajo los nombres posibles (líneas 1211-
        # 1217 del legacy) — normalizamos siempre al nombre canónico
        # crest_conformers.xyz para que OptStage/CensoStage no tengan
        # que saber de esta ambigüedad.
        for candidate_name in self._RESULT_CANDIDATES:
            candidate = crest_dir / candidate_name
            if candidate.exists() and candidate.stat().st_size > 0:
                if candidate != output_xyz:
                    if not safe_copy(candidate, output_xyz):
                        return StageResult(
                            ok=False,
                            message=f"no se pudo normalizar {candidate.name} → {output_xyz.name}",
                        )
                message = (
                    f"CREST terminó en {elapsed:.1f} min "
                    f"(returncode={returncode}, fuente={candidate.name})"
                )
                log.info(f"[CREST]  ✓  {message}")
                return StageResult(ok=True, message=message, data={"output": str(output_xyz)})

        log.error(
            f"[CREST]  ✗  no se encontró ensemble de salida reconocible "
            f"(returncode={returncode}) en {crest_dir} — revisa {log_path.name} a mano."
        )
        return StageResult(
            ok=False,
            message=f"sin ensemble de salida (returncode={returncode}), revisar {log_path}",
        )
