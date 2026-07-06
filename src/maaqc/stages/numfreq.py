"""NumFreqStage — frecuencias vibracionales numéricas (termoquímica).

Migrado desde pipeline.py:
  - stage_numfreq()          líneas 2601-2671
  - has_numfreq_restart()    líneas 903-912 (ahora vive en system/_orca_shared.py)

Nota: config.skip_numfreq (= --skipnumfreq legacy) hace que el
orchestrator salte este stage completo (OPT -> SP -> TDDFT directo),
no que NumFreqStage se ejecute como no-op. Esa decisión vive en
orchestrator/dag.py, no acá.

Cambio de diseño respecto al legacy (mismo patrón que OptStage):
`stage_numfreq(opt_dirs, mol_dir, charge, log)` recibía `opt_dirs` ya
resuelto por el caller (la lista de `successful` que devolvía
`stage_opt`). Como `Stage.run(config, workdir)` no recibe ese output
intermedio de OptStage, `NumFreqStage` descubre sus propios confórmeros
de entrada re-escaneando `02_opt/conf_*` y filtrando por
`check_success(opt.out)` — igual que `OptStage._source_xyz()` re-escanea
`01_crest`/`01b_censo_prescreen` en vez de recibir el ensemble ya
resuelto por CrestStage/CensoStage.

Discrepancia a confirmar con kmilozeta: `_opt_worker` (OptStage) arma un
`env` explícito vía `shared.build_omp_env(nprocs)` para cada ORCA, pero
`stage_numfreq` legacy llama `run_orca(inp_file, out_file, log)` SIN pasar
`env` (usa el default de `run_orca`, `build_omp_env(4)`, sin importar
`NPROCS_NUMFREQ`). Acá se portó 1:1 ese comportamiento del legacy
(sin `env` propio) — si era un bug no detectado y no una decisión
deliberada, avisar para agregar `env=shared.build_omp_env(nprocs)`.
"""

from __future__ import annotations

import logging
from pathlib import Path

from maaqc.config.schema import RecipeConfig
from maaqc.system import _orca_shared as shared
from maaqc.stages.base import Stage, StageResult

log = logging.getLogger("maaqc.stages.numfreq")


class NumFreqStage(Stage):
    name = "numfreq"
    requires = ["opt"]

    STAGE_DIR = "03_numfreq"
    _OPT_DIRNAME = "02_opt"

    # ── localización de los confórmeros OPT exitosos ──────────────────

    def _opt_conf_dirs(self, workdir: Path) -> list[Path]:
        opt_dir = workdir / self._OPT_DIRNAME
        if not opt_dir.is_dir():
            return []
        return sorted(
            conf_dir
            for conf_dir in opt_dir.glob("conf_*")
            if conf_dir.is_dir() and shared.check_success(conf_dir / "opt.out")
        )

    # ── contrato Stage ─────────────────────────────────────────────────

    def is_complete(self, workdir: Path) -> bool:
        """Coarse: al menos un conf_* con freq.out exitoso. Igual criterio
        que OptStage.is_complete — completitud "fina" (¿corrimos NumFreq
        para TODOS los confórmeros OPT esperados?) queda a cargo de
        state.py::needs_rerun, no de este método."""
        stage_dir = workdir / self.STAGE_DIR
        if not stage_dir.is_dir():
            return False
        return any(
            shared.check_success(conf_dir / "freq.out")
            for conf_dir in stage_dir.glob("conf_*")
            if conf_dir.is_dir()
        )

    def inputs(self, workdir: Path) -> list[Path]:
        opt_dir = workdir / self._OPT_DIRNAME
        # Igual que OptStage.inputs: si el directorio de OPT no existe
        # todavía, apuntar ahí igual para que el mensaje de "falta input"
        # de dag.py sea específico.
        return [opt_dir]

    def outputs(self, workdir: Path) -> list[Path]:
        return [workdir / self.STAGE_DIR]

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        if config.charge is None:
            return StageResult(ok=False, message="config.charge no está definido")

        opt_conf_dirs = self._opt_conf_dirs(workdir)
        if not opt_conf_dirs:
            return StageResult(
                ok=False,
                message=f"no se encontraron confórmeros OPT exitosos en "
                        f"{workdir / self._OPT_DIRNAME}",
            )

        charge = config.charge
        mult = config.multiplicity
        nprocs = config.resources.nprocs_numfreq
        maxcore = config.resources.maxcore_numfreq

        template_fresh = shared.get_template("numfreq")
        template_restart = shared.get_template("numfreq_restart")
        stage_dir = workdir / self.STAGE_DIR
        successful: list[Path] = []

        log.info(f"[NUMFREQ] calculando frecuencias para {len(opt_conf_dirs)} confórmeros")

        for opt_dir in opt_conf_dirs:
            conf_id = opt_dir.name
            opt_xyz = shared.resolve_opt_xyz(opt_dir)

            if opt_xyz is None:
                log.warning(f"[NUMFREQ] {conf_id}: no se encontró opt_opt.xyz — saltando")
                continue

            freq_dir = stage_dir / conf_id
            freq_dir.mkdir(parents=True, exist_ok=True)

            out_file = freq_dir / "freq.out"
            inp_file = freq_dir / "freq.inp"

            if not shared.check_success(out_file):
                if not shared.safe_copy(opt_xyz, freq_dir / "opt_opt.xyz"):
                    log.warning(f"[NUMFREQ] {conf_id}: error copiando opt_opt.xyz — saltando")
                    continue

                if shared.has_numfreq_restart(freq_dir):
                    n_done = len(list(freq_dir.glob("freq_D*.res.Gradients")))
                    log.info(
                        f"[NUMFREQ] {conf_id}: restart detectado "
                        f"({n_done} gradientes guardados) — activando Restart true"
                    )
                    inp_text = template_restart.format(
                        CHARGE=charge, MULT=mult,
                        NPROCS_NUMFREQ=nprocs, MAXCORE_NUMFREQ=maxcore,
                    )
                else:
                    inp_text = template_fresh.format(
                        CHARGE=charge, MULT=mult,
                        NPROCS_NUMFREQ=nprocs, MAXCORE_NUMFREQ=maxcore,
                    )
                inp_file.write_text(inp_text)

            ok = shared.run_orca(inp_file, out_file, log)
            if not ok:
                log.warning(f"[NUMFREQ] {conf_id} falló — excluido")
                continue

            if shared.check_no_imaginary(out_file):
                log.info(f"[NUMFREQ] {conf_id}: sin frecuencias imaginarias ✓")
            else:
                log.warning(
                    f"[NUMFREQ] {conf_id}: frecuencias imaginarias detectadas "
                    f"— incluido de todos modos (revisar manualmente)"
                )
            successful.append(freq_dir)

        ok = len(successful) > 0
        message = f"{len(successful)}/{len(opt_conf_dirs)} confórmeros completados"
        log.info(f"[NUMFREQ] {message}")
        return StageResult(
            ok=ok, message=message,
            data={"successful_dirs": [str(p) for p in successful]},
        )
