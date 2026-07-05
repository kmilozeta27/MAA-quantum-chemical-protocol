"""orchestrator/dag.py — Resolución de stages por --entry (Opción 2).

Contexto (ver adc2.py y la discusión de decisión (a)):

  sp.py con skip_numfreq: dependencia que CRECE (opt siempre requerido,
  a veces +numfreq) — se resuelve dentro de SpStage.is_complete(), sin
  tocar este archivo ni Stage.

  adc2.py con --entry sp-only: dependencia que A VECES NO EXISTE (0
  requisitos en vez de ["tddft"]) — reemplaza el viejo `--start adc2`
  standalone de pipeline.py (líneas 3641-3669). Esto SÍ necesita algo
  fuera de Stage, porque el orchestrator debe decidir qué validar antes
  de instanciar nada.

Solución elegida: separar dos preguntas que antes vivían mezcladas en
`requires`:

  1. "¿el stage previo terminó bien?"      → solo se chequea si ese stage
     previo forma parte de la cadena resuelta para este --entry.
  2. "¿existen los archivos que necesito?" → SIEMPRE se chequea, vía
     stage.inputs(workdir), sea cual sea el --entry.

`Stage.requires` NO cambia de tipo (sigue siendo list[str] estático en
los 7 stages). Todo el mecanismo vive acá.
"""

from __future__ import annotations

from pathlib import Path

from maaqc.stages.base import Stage

STAGE_ORDER: list[str] = [
    "crest", "censo", "opt", "numfreq", "sp", "tddft", "adc2",
]

# Entries que NO retoman la cadena normal, sino que corren un único stage
# standalone y solo validan sus inputs() en disco. Reemplaza los viejos
# `--start <stage>` "especiales" de pipeline.py (hoy solo `--start adc2`,
# ver pipeline.py:3641-3669 y el docstring del módulo, líneas 60-68).
#
# Si en el futuro aparece otro --entry con este mismo patrón de "saltar
# cadena completa" (la skill menciona `calibrated --reference shinorine`
# como candidato), agregarlo acá — no como un nuevo `if entry == "..."`
# suelto en el orchestrator.
STANDALONE_ENTRIES: dict[str, str] = {
    "sp-only": "adc2",  # antes: --start adc2
}


def resolve_stages_for_entry(entry: str) -> list[str]:
    """Nombres de stage que el orchestrator debe instanciar/ejecutar.

    - "smiles" / "xyz" / "crest-ensemble" / "calibrated": cadena completa.
      El punto real de retoma dentro de la cadena lo decide is_complete()
      de cada stage (resume robusto vía state.json), no este mapeo.
    - "sp-only" (o cualquier otro en STANDALONE_ENTRIES): un único stage,
      el resto de STAGE_ORDER queda fuera del DAG resuelto para este run.
    """
    if entry in STANDALONE_ENTRIES:
        return [STANDALONE_ENTRIES[entry]]
    return list(STAGE_ORDER)


def is_standalone_entry(entry: str) -> bool:
    return entry in STANDALONE_ENTRIES


def validate_stage_ready(
    stage: Stage,
    workdir: Path,
    entry: str,
    stages_by_name: dict[str, Stage],
) -> tuple[bool, str | None]:
    """Decide si `stage` puede correr, sin tocar el atributo estático
    `stage.requires`.

    Siempre valida `stage.inputs(workdir)` (pregunta 2). Si `entry` es
    standalone para este stage, NO valida la cadena `stage.requires`
    (pregunta 1) — exactamente lo que --entry sp-only necesita para
    Adc2Stage sin volver dinámico Stage.requires. En cadena normal, sí
    valida is_complete() de cada requisito, igual que antes.
    """
    missing = [p for p in stage.inputs(workdir) if not p.exists()]
    if missing:
        faltantes = ", ".join(str(p) for p in missing)
        return False, f"Faltan inputs para '{stage.name}': {faltantes}"

    if is_standalone_entry(entry):
        return True, None

    for prev_name in stage.requires:
        prev_stage = stages_by_name.get(prev_name)
        if prev_stage is None:
            return False, f"Stage previo desconocido en el registro: '{prev_name}'"
        if not prev_stage.is_complete(workdir):
            return False, (
                f"'{stage.name}' requiere que '{prev_name}' esté completo "
                f"y no lo está."
            )
    return True, None
