"""stages/censo.py — Etapa 01b_censo_prescreen: pre-screening CENSO
(obligatorio entre CREST y OPT).

Migrado desde pipeline.py, cableando `CensoStage` a la lógica ya
implementada en `system/_censo_shared.py` (run_censo_prescreen,
run_censo_multitemp, select_nconf_auto, select_nconf_boltzmann,
extract_conformers_by_id) — no hay lógica nueva acá, solo el stage
delgado que la orquesta, siguiendo el mismo patrón que OptStage/SpStage
envolviendo a `_orca_shared`.

Detalle importante que no estaba explícito en el stub original:
`run_censo_prescreen` deja el ensemble *ordenado* (`censo_ensemble_sorted.xyz`)
pero sin recortar a `nconf`. `OptStage._source_xyz()` busca
`01b_censo_prescreen/crest_conformers.xyz` — no `censo_ensemble_sorted.xyz`.
Sin un paso intermedio que aplique select_nconf_auto/boltzmann y escriba
el resultado con el nombre que OptStage espera, CENSO corría pero OPT
nunca lo veía y caía silenciosamente al ensemble crudo de 01_crest/. Este
stage hace ese último paso (`_select_final_ensemble`) explícito.

Supuestos sobre `config.censo` (no tengo `config/schema.py` a la vista —
avisar si los nombres reales difieren):
    config.censo.multitemp            bool
    config.censo.temperature          float   (K, modo single)
    config.censo.t_ini/t_fin/n_temps  float/float/int  (modo multitemp)
    config.censo.timeout_hours        float
    config.censo.selection            "auto" | "boltzmann"
    config.censo.boltzmann_threshold  float   (%, default 5.0)
    config.censo.ncon_min/ncon_max    int/int
    config.resources.nprocs_crest     int
    config.multiplicity               int  (unpaired = multiplicity - 1)

is_complete() sigue siendo un chequeo "grueso" (¿existe un
crest_conformers.xyz utilizable?), igual que en OptStage/SpStage —
delegar en orchestrator.state.rebuild_censo_state() queda pendiente de
que ese módulo exista; por ahora reimplementa el mismo criterio de
"screening parseable" que ya usa run_censo_prescreen internamente, en
vez de confiar en la sola presencia del archivo.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

from maaqc.config.schema import RecipeConfig
from maaqc.stages.base import Stage, StageResult
from maaqc.system import _censo_shared as shared

log = logging.getLogger("maaqc.stages.censo")


class CensoStage(Stage):
    name = "censo"
    requires: list[str] = ["crest"]

    STAGE_DIR = "01b_censo_prescreen"
    CREST_DIR = "01_crest"
    _ENSEMBLE_FILENAME = "crest_conformers.xyz"  # nombre que espera OptStage
    _SCREENING_RELPATHS = (
        "1_SCREENING.out",
        "CENSO/1_SCREENING.out",
        "prescreening/1_SCREENING.out",
    )

    # ── localización del ensemble de entrada (crudo de CREST) ──────────

    def _crest_xyz(self, workdir: Path) -> Optional[Path]:
        candidate = workdir / self.CREST_DIR / self._ENSEMBLE_FILENAME
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate
        return None

    def _find_screening_out(self, censo_dir: Path) -> Optional[Path]:
        for rel in self._SCREENING_RELPATHS:
            p = censo_dir / rel
            if shared.parse_censo_screening_out(p):
                return p
        return None

    # ── selección post-CENSO (Boltzmann o auto) + escritura final ──────

    def _select_final_ensemble(
        self, censo_dir: Path, output_xyz: Path, config: RecipeConfig, target: Path,
    ) -> tuple[bool, str]:
        """A partir del ensemble ordenado que devolvió CENSO en `censo_dir`,
        selecciona nconf confórmeros y escribe el resultado en `target`
        (crest_conformers.xyz) para que OptStage lo encuentre."""
        censo_cfg = config.censo
        screening_out = self._find_screening_out(censo_dir)
        if screening_out is None:
            return False, "CENSO no produjo un 1_SCREENING.out parseable — nada que seleccionar"

        all_ids = shared.parse_censo_screening_out(screening_out) or []
        selection_mode = getattr(censo_cfg, "selection", "auto")

        conf_ids: list[str] = []
        if selection_mode == "boltzmann":
            _, conf_ids = shared.select_nconf_boltzmann(
                conformers_xyz=output_xyz, screening_out=screening_out, log=log,
                threshold=censo_cfg.boltzmann_threshold, nconf_min=censo_cfg.ncon_min,
            )
        if not conf_ids:
            # selection_mode == "auto", o boltzmann no encontró pesos
            # (parse_boltzmann_weights vacío) — cae en select_nconf_auto.
            nconf = shared.select_nconf_auto(
                output_xyz, log,
                nconf_min=censo_cfg.ncon_min, nconf_max=censo_cfg.ncon_max,
            )
            conf_ids = all_ids[:nconf]

        if not conf_ids:
            return False, f"selección vacía a partir de {screening_out.name}"

        if not shared.extract_conformers_by_id(output_xyz, conf_ids, target, log):
            return False, f"extract_conformers_by_id falló escribiendo {target.name}"

        return True, f"{len(conf_ids)} confórmero(s) seleccionados → {target.name} ({censo_dir.name})"

    # ── contrato Stage ───────────────────────────────────────────────────

    def is_complete(self, workdir: Path) -> bool:
        stage_dir = workdir / self.STAGE_DIR
        if not stage_dir.is_dir():
            return False

        summary = stage_dir / "multitemp_summary.json"
        if summary.exists():
            temp_dirs = [d for d in stage_dir.iterdir() if d.is_dir() and d.name.startswith("T")]
            if not temp_dirs:
                return False
            return all(
                (d / self._ENSEMBLE_FILENAME).exists()
                and (d / self._ENSEMBLE_FILENAME).stat().st_size > 0
                for d in temp_dirs
            )

        target = stage_dir / self._ENSEMBLE_FILENAME
        return target.exists() and target.stat().st_size > 0

    def inputs(self, workdir: Path) -> list[Path]:
        source = self._crest_xyz(workdir)
        return [source] if source else [workdir / self.CREST_DIR / self._ENSEMBLE_FILENAME]

    def outputs(self, workdir: Path) -> list[Path]:
        return [workdir / self.STAGE_DIR]

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        if config.charge is None:
            return StageResult(ok=False, message="config.charge no está definido")

        source_xyz = self._crest_xyz(workdir)
        if source_xyz is None:
            return StageResult(
                ok=False,
                message=f"no se encontró ensemble en {self.CREST_DIR}/{self._ENSEMBLE_FILENAME}",
            )

        unpaired = config.multiplicity - 1
        censo_cfg = config.censo
        stage_dir = workdir / self.STAGE_DIR
        stage_dir.mkdir(parents=True, exist_ok=True)
        nprocs_crest = config.resources.nprocs_crest
        timeout_hours = censo_cfg.timeout_hours

        log.info(f"[CENSO] pre-screening desde {source_xyz} (carga={config.charge})")

        # ── modo multitemp: una corrida de CENSO por temperatura ────────
        if censo_cfg.multitemp:
            temperatures = shared.trange_to_temperatures(
                censo_cfg.t_ini, censo_cfg.t_fin, censo_cfg.n_temps,
            )
            raw_results = shared.run_censo_multitemp(
                conformers_xyz=source_xyz, base_censo_dir=stage_dir,
                charge=config.charge, unpaired=unpaired,
                temperatures=temperatures, nprocs_crest=nprocs_crest,
                timeout_hours=timeout_hours, log=log,
            )

            n_ok = 0
            details: dict[str, str] = {}
            for temp, out_xyz in raw_results.items():
                subdir = stage_dir / f"T{temp:.2f}K"
                target = subdir / self._ENSEMBLE_FILENAME
                if out_xyz == source_xyz:
                    details[f"{temp:.2f}"] = "CENSO falló (fail-open) — sin ensemble filtrado"
                    continue
                ok, msg = self._select_final_ensemble(subdir, out_xyz, config, target)
                details[f"{temp:.2f}"] = msg
                if ok:
                    n_ok += 1

            message = f"[multitemp] {n_ok}/{len(temperatures)} temperaturas con ensemble seleccionado"
            log.info(f"[CENSO] {message}")
            return StageResult(ok=n_ok > 0, message=message, data=details)

        # ── modo single-T ────────────────────────────────────────────────
        output_xyz = shared.run_censo_prescreen(
            conformers_xyz=source_xyz, censo_dir=stage_dir,
            charge=config.charge, unpaired=unpaired,
            temperature=censo_cfg.temperature, nprocs_crest=nprocs_crest,
            timeout_hours=timeout_hours, log=log,
        )

        if output_xyz == source_xyz:
            # Fail-open: run_censo_prescreen ya logueó la causa exacta
            # (CENSO ausente, timeout, sin ensemble reconocible).
            message = "CENSO no produjo un ensemble filtrado — OPT usará 01_crest/ sin pasar por CENSO"
            log.warning(f"[CENSO] {message}")
            return StageResult(ok=False, message=message)

        target = stage_dir / self._ENSEMBLE_FILENAME
        ok, msg = self._select_final_ensemble(stage_dir, output_xyz, config, target)
        log.info(f"[CENSO] {msg}")
        return StageResult(ok=ok, message=msg, data={"censo_output": str(output_xyz)})
