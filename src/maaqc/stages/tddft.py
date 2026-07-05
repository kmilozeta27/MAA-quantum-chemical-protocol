"""TddftStage — TD-DFT sobre los confórmeros que completaron SP.

Migrado desde pipeline.py:
  - stage_tddft()       líneas 2814-2880
  - _tddft_worker()     líneas 2748-2813

Corre después de SP salvo que config.entry sea posterior en la cadena
(ej. --entry sp-only detiene antes de llegar acá).

Geometría/orbitales de entrada: reusa `sp.gbw` (guess de orbitales del
SP) y la geometría resuelta por `resolve_opt_xyz(sp_dir)` — exactamente
como el legacy, copiados a `sp.gbw`/`sp_opt.xyz` dentro de 05_tddft/<id>/.

⚠️ CAMBIO RESPECTO AL LEGACY (confirmar con kmilozeta antes de correr algo
citable): TEMPLATE_TDDFT ya no hardcodea "wB97X-D3"/"NRoots 15" — usa
`config.tddft.functional`/`config.tddft.states`, que ya existían en el
schema pero nunca se conectaban a nada real. El default de `states` en
el schema es 10, no el 15 del legacy — un recipe con `tddft.states: 15`
restaura el comportamiento legacy exacto.
"""

from __future__ import annotations

import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

from maaqc.config.schema import RecipeConfig
from maaqc.stages.base import Stage, StageResult
from maaqc.system import _orca_shared as shared

log = logging.getLogger("maaqc.stages.tddft")


class TddftStage(Stage):
    name = "tddft"
    requires = ["sp"]

    STAGE_DIR = "05_tddft"
    SP_DIR = "04_sp"

    # ── contrato Stage ─────────────────────────────────────────────────

    def is_complete(self, workdir: Path) -> bool:
        stage_dir = workdir / self.STAGE_DIR
        if not stage_dir.is_dir():
            return False
        return any(
            shared.check_success(conf_dir / "tddft.out")
            for conf_dir in stage_dir.glob("conf_*")
            if conf_dir.is_dir()
        )

    def inputs(self, workdir: Path) -> list[Path]:
        return [workdir / self.SP_DIR]

    def outputs(self, workdir: Path) -> list[Path]:
        return [workdir / self.STAGE_DIR]

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        if config.charge is None:
            return StageResult(ok=False, message="config.charge no está definido")

        sp_dir = workdir / self.SP_DIR
        sp_conf_dirs = sorted(d for d in sp_dir.glob("conf_*") if d.is_dir()) \
            if sp_dir.is_dir() else []
        if not sp_conf_dirs:
            return StageResult(ok=False, message=f"no hay confórmeros en {sp_dir}/ — corre primero SP")

        nprocs = config.resources.nprocs_tddft
        maxcore = config.resources.maxcore_tddft
        batch_size = config.resources.tddft_batch_size
        template = shared.get_template("tddft")
        stage_dir = workdir / self.STAGE_DIR

        log.info(
            f"[TDDFT] {len(sp_conf_dirs)} confórmeros — batch {batch_size}x{nprocs} cores "
            f"({config.tddft.functional}, NRoots {config.tddft.states})"
        )

        tasks_all = [
            {
                "conf_id": d.name,
                "sp_dir": d,
                "tddft_dir": stage_dir / d.name,
                "charge": config.charge,
                "mult": config.multiplicity,
                "template": template,
                "functional": config.tddft.functional,
                "nroots": config.tddft.states,
                "nprocs": nprocs,
                "maxcore": maxcore,
            }
            for d in sp_conf_dirs
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

            log.info(f"[TDDFT] ── Lote {lote_num}/{n_lotes}: {ids}")

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
                    shared.kill_active_procs()
                    raise

        ok = len(successful) > 0
        message = f"{len(successful)}/{len(sp_conf_dirs)} completados"
        log.info(f"[TDDFT] {message}")
        return StageResult(
            ok=ok, message=message,
            data={"successful_dirs": [str(p) for p in successful]},
        )


def _tddft_worker(task: dict) -> dict:
    """Worker de un único TDDFT. Ported 1:1 desde pipeline.py:2748-2813,
    adaptado para tomar functional/nroots de la task en vez de constantes
    de módulo."""
    conf_id: str = task["conf_id"]
    sp_dir: Path = task["sp_dir"]
    tddft_dir: Path = task["tddft_dir"]
    cpu_list: Optional[list[int]] = task.get("cpu_list")

    gbw_file = sp_dir / "sp.gbw"
    opt_xyz = shared.resolve_opt_xyz(sp_dir)

    if not gbw_file.exists():
        return {"conf_id": conf_id, "tddft_dir": None, "success": False, "elapsed": 0.0,
                "message": "no se encontró sp.gbw"}
    if opt_xyz is None:
        return {"conf_id": conf_id, "tddft_dir": None, "success": False, "elapsed": 0.0,
                "message": "no se encontró opt_opt.xyz ni opt.xyz"}

    tddft_dir.mkdir(parents=True, exist_ok=True)
    out_file = tddft_dir / "tddft.out"
    inp_file = tddft_dir / "tddft.inp"

    if not shared.check_success(out_file):
        if not shared.safe_copy(opt_xyz, tddft_dir / "sp_opt.xyz"):
            return {"conf_id": conf_id, "tddft_dir": None, "success": False, "elapsed": 0.0,
                    "message": "error copiando sp_opt.xyz"}
        if not shared.safe_copy(gbw_file, tddft_dir / "sp.gbw"):
            return {"conf_id": conf_id, "tddft_dir": None, "success": False, "elapsed": 0.0,
                    "message": "error copiando sp.gbw"}
        inp_text = task["template"].format(
            CHARGE=task["charge"], MULT=task["mult"],
            NPROCS=task["nprocs"], MAXCORE=task["maxcore"],
            FUNCTIONAL=task["functional"], NROOTS=task["nroots"],
            GBW_PATH="sp.gbw",
        )
        inp_file.write_text(inp_text)

    t0 = time.time()
    worker_env = shared.build_omp_env(task["nprocs"])
    success = shared.run_orca(inp_file, out_file, log, env=worker_env, cpu_list=cpu_list)
    elapsed = (time.time() - t0) / 60

    return {
        "conf_id": conf_id,
        "tddft_dir": tddft_dir if success else None,
        "success": success,
        "elapsed": elapsed,
    }
