"""stages/opt.py — Etapa 02_opt: optimización DFT de cada confórmero.

Migrado desde pipeline.py (`stage_opt` + `_opt_worker`, líneas ~2462 y
~3136). Cambio de diseño deliberado respecto al legacy, a confirmar con
kmilozeta:

    pipeline.py legacy mezclaba dentro de `stage_opt` dos
    responsabilidades: (1) SELECCIONAR cuántos/cuáles confórmeros optimizar
    (nconf de select_nconf_auto/select_nconf_boltzmann, o --boltzmann-threshold
    con conf_ids explícitos) y (2) OPTIMIZARLOS. Como `Stage.run(config,
    workdir)` no recibe ese `nconf`/`conf_ids` como parámetro aparte,
    `OptStage` asume que la selección YA ocurrió upstream (en CrestStage o
    en la futura CensoStage, que debería escribir el ensemble ya recortado)
    y simplemente optimiza TODOS los frames del archivo de confórmeros que
    encuentra en `_source_xyz()`. Si la lógica de Boltzmann-threshold /
    select_nconf_auto todavía no vive en ningún stage migrado, esto
    optimizaría de más (todo el ensemble CREST) hasta que CensoStage (o
    quien corresponda) quede migrado — avisar si esto no calza con el
    plan real.
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

log = logging.getLogger("maaqc.stages.opt")


class OptStage(Stage):
    name = "opt"
    requires: list[str] = ["crest"]

    STAGE_DIR = "02_opt"
    # Prioridad de fuente del ensemble: si ya existe un prescreen de CENSO,
    # se usa ese (ensemble recortado); si no, el crudo de CREST.
    _CANDIDATE_SOURCE_DIRS = ["01b_censo_prescreen", "01_crest"]
    _ENSEMBLE_FILENAME = "crest_conformers.xyz"

    # ── localización del ensemble de entrada ──────────────────────────

    def _source_xyz(self, workdir: Path) -> Optional[Path]:
        for dirname in self._CANDIDATE_SOURCE_DIRS:
            candidate = workdir / dirname / self._ENSEMBLE_FILENAME
            if candidate.exists() and candidate.stat().st_size > 0:
                return candidate
        return None

    # ── contrato Stage ─────────────────────────────────────────────────

    def is_complete(self, workdir: Path) -> bool:
        """Coarse: al menos un conf_* con opt.out exitoso. Igual criterio
        que `print_status()` en pipeline.py legacy — completitud "fina"
        (¿optimizamos TODOS los confórmeros esperados?) queda a cargo de
        `state.py::needs_rerun` comparando config_hash/fingerprint, no de
        este método."""
        stage_dir = workdir / self.STAGE_DIR
        if not stage_dir.is_dir():
            return False
        return any(
            shared.check_success(conf_dir / "opt.out")
            for conf_dir in stage_dir.glob("conf_*")
            if conf_dir.is_dir()
        )

    def inputs(self, workdir: Path) -> list[Path]:
        source = self._source_xyz(workdir)
        # Si no hay fuente todavía, apunta al candidato de mayor prioridad
        # para que el mensaje de "falta input" de dag.py sea específico.
        return [source] if source else [
            workdir / self._CANDIDATE_SOURCE_DIRS[-1] / self._ENSEMBLE_FILENAME
        ]

    def outputs(self, workdir: Path) -> list[Path]:
        return [workdir / self.STAGE_DIR]

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        if config.charge is None:
            return StageResult(ok=False, message="config.charge no está definido")

        source_xyz = self._source_xyz(workdir)
        if source_xyz is None:
            return StageResult(
                ok=False,
                message=f"no se encontró ensemble de confórmeros en {self._CANDIDATE_SOURCE_DIRS}",
            )

        try:
            conformers = shared.read_all_conformers(source_xyz)
        except RuntimeError as e:
            return StageResult(ok=False, message=str(e))
        if not conformers:
            return StageResult(ok=False, message=f"{source_xyz} está vacío")

        conf_id_map = shared.extract_conf_ids_from_xyz(source_xyz)
        crest_labels = [
            conf_id_map.get(i, f"CONF{i + 1}") for i in range(len(conformers))
        ]

        nprocs = config.resources.nprocs
        maxcore = config.resources.maxcore
        batch_size = config.resources.opt_batch_size
        template = shared.get_template("opt")
        stage_dir = workdir / self.STAGE_DIR

        log.info(
            f"[OPT] {len(conformers)} confórmeros desde {source_xyz} — "
            f"batch {batch_size}x{nprocs} cores"
        )

        tasks_all = [
            {
                "conf_id": f"conf_{lbl[4:]}",
                "conf_dir": stage_dir / f"conf_{lbl[4:]}",
                "atoms": atoms,
                "charge": config.charge,
                "mult": config.multiplicity,
                "template": template,
                "nprocs": nprocs,
                "maxcore": maxcore,
            }
            for lbl, atoms in zip(crest_labels, conformers)
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

            log.info(f"[OPT] ── Lote {lote_num}/{n_lotes}: {ids}")

            with ThreadPoolExecutor(max_workers=len(batch)) as pool:
                futures = {pool.submit(_opt_worker, t): t["conf_id"] for t in batch}
                try:
                    for fut in as_completed(futures):
                        res = fut.result()
                        tag = "✓" if res["success"] else "✗"
                        log.info(f"[OPT]   {tag} {res['conf_id']}  ({res['elapsed']:.1f} min)")
                        if res["success"]:
                            successful.append(res["conf_dir"])
                        else:
                            log.warning(f"[OPT] {res['conf_id']} falló — excluido")
                except KeyboardInterrupt:
                    log.warning("[OPT] Ctrl-C — cancelando lote, esperando workers…")
                    pool.shutdown(wait=False, cancel_futures=True)
                    shared.kill_active_procs()
                    raise

        ok = len(successful) > 0
        message = f"{len(successful)}/{len(conformers)} confórmeros optimizados"
        log.info(f"[OPT] {message}")
        return StageResult(
            ok=ok, message=message,
            data={"successful_dirs": [str(p) for p in successful]},
        )


def _opt_worker(task: dict) -> dict:
    """Worker de un único OPT (ThreadPoolExecutor — ORCA libera el GIL al
    correr como subproceso, el paralelismo real es efectivo). Ported 1:1
    desde pipeline.py:3136."""
    conf_id: str = task["conf_id"]
    conf_dir: Path = task["conf_dir"]
    atoms = task["atoms"]
    cpu_list: Optional[list[int]] = task.get("cpu_list")

    inp_file = conf_dir / "opt.inp"
    out_file = conf_dir / "opt.out"
    conf_dir.mkdir(parents=True, exist_ok=True)

    if not shared.check_success(out_file):
        coords = shared.atoms_to_xyz_block(atoms)
        inp_text = task["template"].format(
            CHARGE=task["charge"], MULT=task["mult"],
            NPROCS=task["nprocs"], MAXCORE=task["maxcore"],
            COORDS=coords,
        )
        inp_file.write_text(inp_text)

    t0 = time.time()
    worker_env = shared.build_omp_env(task["nprocs"])
    success = shared.run_orca(inp_file, out_file, log, env=worker_env, cpu_list=cpu_list)
    elapsed = (time.time() - t0) / 60

    if success and shared.resolve_opt_xyz(conf_dir) is None:
        log.warning(f"[OPT] {conf_id}: ORCA OK pero sin geometría optimizada")
        success = False

    return {
        "conf_id": conf_id,
        "conf_dir": conf_dir,
        "success": success,
        "elapsed": elapsed,
    }
