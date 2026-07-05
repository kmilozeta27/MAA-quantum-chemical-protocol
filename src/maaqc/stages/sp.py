"""stages/sp.py — Etapa 04_sp: single point de alta precisión.

Migrado desde pipeline.py (`stage_sp` + `_sp_worker`, líneas ~2672 y
~3201). Replica el patrón de "dependencia que crece" documentado en
dag.py: `SpStage.requires` se queda estático en `["opt"]` — NumFreq es
un input OPCIONAL cuya necesidad depende de `config.skip_numfreq`, así
que esa decisión se toma acá adentro (en `_source_dir`/`run`), no en
`Stage.requires` ni en `orchestrator/dag.py` (ver docstring de dag.py,
sección sobre `sp.py con skip_numfreq`).

Geometría de entrada: SP toma `opt_opt.xyz` de 03_numfreq/<conf_id>/ (si
`skip_numfreq=False`, comportamiento estándar — la NumFreq no cambia la
geometría, solo la recalcula desde la misma) o directamente de
02_opt/<conf_id>/ (si `skip_numfreq=True`), exactamente como en
pipeline.py:3782-3788 (`freq_dirs = opt_dirs` bajo `--skipnumfreq`).
"""

from __future__ import annotations

import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from maaqc.config.schema import RecipeConfig
from maaqc.system import _orca_shared as shared
from maaqc.stages.base import Stage, StageResult

log = logging.getLogger("maaqc.stages.sp")


class SpStage(Stage):
    name = "sp"
    requires: list[str] = ["opt"]

    STAGE_DIR = "04_sp"
    OPT_DIR = "02_opt"
    NUMFREQ_DIR = "03_numfreq"

    # ── localización de la geometría de entrada ────────────────────────

    def _numfreq_conf_dirs(self, workdir: Path) -> list[Path]:
        numfreq_dir = workdir / self.NUMFREQ_DIR
        if not numfreq_dir.is_dir():
            return []
        return sorted(d for d in numfreq_dir.glob("conf_*") if d.is_dir())

    def _opt_conf_dirs(self, workdir: Path) -> list[Path]:
        opt_dir = workdir / self.OPT_DIR
        if not opt_dir.is_dir():
            return []
        return sorted(d for d in opt_dir.glob("conf_*") if d.is_dir())

    def _source_dirs_for_validation(self, workdir: Path) -> list[Path]:
        """Usado por `inputs()`, que no recibe `config`: prioriza NumFreq
        si ya existe en disco, si no cae a OPT — solo para que dag.py
        pueda validar existencia antes de correr. La decisión REAL
        (según `config.skip_numfreq`) se hace en `run()`."""
        numfreq_dirs = self._numfreq_conf_dirs(workdir)
        if numfreq_dirs:
            return [workdir / self.NUMFREQ_DIR]
        return [workdir / self.OPT_DIR]

    # ── contrato Stage ─────────────────────────────────────────────────

    def is_complete(self, workdir: Path) -> bool:
        stage_dir = workdir / self.STAGE_DIR
        if not stage_dir.is_dir():
            return False
        return any(
            shared.check_success(conf_dir / "sp.out")
            for conf_dir in stage_dir.glob("conf_*")
            if conf_dir.is_dir()
        )

    def inputs(self, workdir: Path) -> list[Path]:
        return self._source_dirs_for_validation(workdir)

    def outputs(self, workdir: Path) -> list[Path]:
        return [workdir / self.STAGE_DIR]

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        if config.charge is None:
            return StageResult(ok=False, message="config.charge no está definido")

        if config.skip_numfreq:
            source_dirs = self._opt_conf_dirs(workdir)
            source_label = self.OPT_DIR
        else:
            source_dirs = self._numfreq_conf_dirs(workdir)
            source_label = self.NUMFREQ_DIR

        if not source_dirs:
            hint = (
                f"corre primero OPT (o desactiva skip_numfreq incorrecto)"
                if config.skip_numfreq else
                f"corre primero NumFreq, o pon skip_numfreq=true para usar {self.OPT_DIR}/ directamente"
            )
            return StageResult(ok=False, message=f"no hay confórmeros en {source_label}/ — {hint}")

        nprocs = config.resources.nprocs
        maxcore = config.resources.maxcore
        batch_size = config.resources.sp_batch_size
        template = shared.get_template("sp")
        stage_dir = workdir / self.STAGE_DIR

        log.info(
            f"[SP] {len(source_dirs)} confórmeros desde {source_label}/ — "
            f"batch {batch_size}x{nprocs} cores"
        )

        tasks_all = [
            {
                "conf_id": d.name,
                "source_dir": d,
                "sp_dir": stage_dir / d.name,
                "charge": config.charge,
                "mult": config.multiplicity,
                "template": template,
                "nprocs": nprocs,
                "maxcore": maxcore,
            }
            for d in source_dirs
        ]

        shared.check_cpu_budget(batch_size, nprocs, log)
        successful: list[Path] = []

        for b_start in range(0, len(tasks_all), batch_size):
            batch = tasks_all[b_start: b_start + batch_size]
            lote_num = b_start // batch_size + 1
            n_lotes = -(-len(tasks_all) // batch_size)
            ids = [t["conf_id"] for t in batch]

            cpu_lists = shared.assign_cpu_lists(len(batch), nprocs)
            for t, cl in zip(batch, cpu_lists):
                t["cpu_list"] = cl

            log.info(f"[SP] ── Lote {lote_num}/{n_lotes}: {ids}")

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
                    shared.kill_active_procs()
                    raise

        ok = len(successful) > 0
        message = f"{len(successful)}/{len(source_dirs)} completados"
        log.info(f"[SP] {message}")
        return StageResult(
            ok=ok, message=message,
            data={"successful_dirs": [str(p) for p in successful]},
        )


def _sp_worker(task: dict) -> dict:
    """Worker de un único SP. Ported 1:1 desde pipeline.py:3201."""
    conf_id: str = task["conf_id"]
    source_dir: Path = task["source_dir"]
    sp_dir: Path = task["sp_dir"]
    cpu_list: Optional[list[int]] = task.get("cpu_list")

    opt_xyz = shared.resolve_opt_xyz(source_dir)
    if opt_xyz is None:
        log.warning(f"[SP] {conf_id}: no se encontró opt_opt.xyz ni opt.xyz — saltando")
        return {"conf_id": conf_id, "sp_dir": None, "success": False, "elapsed": 0.0}

    sp_dir.mkdir(parents=True, exist_ok=True)
    out_file = sp_dir / "sp.out"
    inp_file = sp_dir / "sp.inp"

    if not shared.check_success(out_file):
        if not shared.safe_copy(opt_xyz, sp_dir / "opt_opt.xyz"):
            log.warning(f"[SP] {conf_id}: error copiando opt_opt.xyz — saltando")
            return {"conf_id": conf_id, "sp_dir": None, "success": False, "elapsed": 0.0}
        inp_text = task["template"].format(
            CHARGE=task["charge"], MULT=task["mult"],
            NPROCS=task["nprocs"], MAXCORE=task["maxcore"],
        )
        inp_file.write_text(inp_text)

    t0 = time.time()
    worker_env = shared.build_omp_env(task["nprocs"])
    success = shared.run_orca(inp_file, out_file, log, env=worker_env, cpu_list=cpu_list)
    elapsed = (time.time() - t0) / 60

    return {
        "conf_id": conf_id,
        "sp_dir": sp_dir if success else None,
        "success": success,
        "elapsed": elapsed,
    }
