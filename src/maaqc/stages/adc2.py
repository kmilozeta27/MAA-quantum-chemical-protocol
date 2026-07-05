"""Adc2Stage — ADC2 sobre el confórmero de menor energía SP. OPCIONAL.

Migrar desde pipeline.py:
  - stage_adc2()          líneas 2972-3135
  - extract_adc2_lambda() líneas 929-998
  - load_exp_lambda()     líneas 270-351 (usa EXP_LAMBDA_FALLBACK / molecules.yml
                                           -> ahora vive en config.adc2)

No corre por defecto: requiere config.adc2.mandatory = True (equivalente
a --adc2 explícito en el legacy) O config.entry == "sp-only" (equivalente
a --start adc2 standalone).

Caso especial --entry sp-only (antes "--start adc2"): NO retoma la
cadena normal de requires. Solo exige que exista 04_sp/ (.out +
opt_opt.xyz por confórmero) — ignora crest/censo/opt/numfreq/tddft
incluso si no existen en disco. Esto es intencional: el flujo real de
Kmilo respalda todo en disco externo y borra los archivos pesados antes
de correr ADC2 (~400GB de disco requeridos), dejando solo 04_sp/
disponible localmente.

TODO de diseño: esto rompe el contrato "requires = lista estática de
stages previos que deben estar is_complete()" definido en base.py. Con
--entry sp-only, Adc2Stage.requires debería evaluarse como [] (no
validar crest/censo/opt/numfreq/tddft) en vez de la cadena completa.
Ver la misma discusión pendiente en sp.py sobre requires dinámico según
config — resolver ambos casos con el mismo mecanismo en
orchestrator/dag.py, no con un if especial solo acá.
"""

from __future__ import annotations

from pathlib import Path

from maaqc.config.schema import RecipeConfig
from maaqc.stages.base import Stage, StageResult


class Adc2Stage(Stage):
    name = "adc2"
    requires = ["tddft"]  # TODO: [] cuando config.entry == "sp-only" (ver docstring)

    def is_complete(self, workdir: Path) -> bool:
        raise NotImplementedError("Migrar desde pipeline.py:stage_adc2 (líneas 2972-3135)")

    def inputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def outputs(self, workdir: Path) -> list[Path]:
        raise NotImplementedError

    def run(self, config: RecipeConfig, workdir: Path) -> StageResult:
        raise NotImplementedError
