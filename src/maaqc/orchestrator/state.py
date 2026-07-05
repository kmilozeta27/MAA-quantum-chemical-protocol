"""orchestrator/state.py — Resume robusto vía state.json explícito.

Generaliza el parche puntual que ya se hizo para recovery de
1_SCREENING.out (matching por energía DFT en vez de tokens de string):
en vez de que cada stage decida "¿ya corrí?" mirando patrones frágiles
dentro de un .out, el orchestrator guarda un registro explícito por
molécula con:

  1. un hash de la porción de `config` relevante para ESE stage (para
     detectar "config cambió, hay que re-correr" incluso si los archivos
     de output todavía existen en disco — ej. cambiaste el funcional de
     TDDFT pero 05_tddft/ quedó de la corrida anterior); y
  2. un fingerprint de los outputs() del stage en el momento en que
     terminó (para detectar "algo tocó/borró estos archivos por fuera
     del orchestrator" — ej. Kmilo corrió un .out a mano o restauró un
     backup parcial).

`is_complete()` de cada Stage sigue siendo la fuente de verdad sobre
"¿el .out dice ORCA TERMINATED NORMALLY?" — state.json no la reemplaza,
la complementa: un stage puede seguir siendo is_complete()==True y aun
así necesitar re-correr si config_hash cambió o el fingerprint no calza.

No requiere tocar Stage ni dag.py: state.py es una capa aparte que el
orchestrator (dag.py / la futura CLI `maaqc run`) consulta antes de
decidir si saltarse un stage.
"""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

STATE_FILENAME = "state.json"
SCHEMA_VERSION = 1


# ─────────────────────────────────────────────
#  MODELO DE ESTADO
# ─────────────────────────────────────────────

@dataclass
class StageRecord:
    completed: bool
    config_hash: str
    output_fingerprint: dict[str, str]
    timestamp: str
    message: str = ""


@dataclass
class RunState:
    schema_version: int = SCHEMA_VERSION
    stages: dict[str, StageRecord] = field(default_factory=dict)

    def to_json(self) -> str:
        payload = {
            "schema_version": self.schema_version,
            "stages": {name: asdict(rec) for name, rec in self.stages.items()},
        }
        return json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n"

    @classmethod
    def from_json(cls, text: str) -> "RunState":
        payload = json.loads(text)
        stages = {
            name: StageRecord(**rec)
            for name, rec in payload.get("stages", {}).items()
        }
        return cls(schema_version=payload.get("schema_version", SCHEMA_VERSION), stages=stages)


# ─────────────────────────────────────────────
#  IO DE state.json
# ─────────────────────────────────────────────

def state_path(workdir: Path) -> Path:
    return workdir / STATE_FILENAME


def load_state(workdir: Path) -> RunState:
    path = state_path(workdir)
    if not path.exists():
        return RunState()
    try:
        return RunState.from_json(path.read_text())
    except (OSError, json.JSONDecodeError, TypeError):
        # state.json corrupto/incompatible: no abortar el pipeline por
        # esto — se trata como "sin estado previo" y se reconstruye desde
        # cero al terminar el próximo stage.
        return RunState()


def save_state(workdir: Path, state: RunState) -> None:
    state_path(workdir).write_text(state.to_json())


# ─────────────────────────────────────────────
#  HASH DE CONFIG (por stage)
# ─────────────────────────────────────────────

def _to_jsonable(value: Any) -> Any:
    """Convierte objetos de config (pydantic/dataclass/dict/lo que sea) a
    algo serializable de forma determinista, sin asumir una implementación
    concreta de RecipeConfig."""
    if hasattr(value, "model_dump"):       # pydantic v2
        return value.model_dump()
    if hasattr(value, "dict"):             # pydantic v1
        return value.dict()
    if hasattr(value, "__dict__"):
        return {k: v for k, v in vars(value).items() if not k.startswith("_")}
    return value


def compute_config_hash(config: Any, stage_name: str) -> str:
    """Hash estable de los campos de `config` relevantes para `stage_name`.

    Incluye:
      - el bloque de config específico del stage (ej. config.adc2), que en
        el YAML de recipes vive bajo una clave con el mismo nombre que el
        stage (ver presets/publication.yml en la skill: crest:, tddft:,
        adc2:, ...).
      - `config.charge`, porque afecta a TODOS los stages (multiplicidad/
        carga cambian el .inp aunque el resto del bloque no cambie).

    Si `config` no tiene el bloque del stage (ej. stage sin sección propia
    en el YAML), hashea solo charge — no es un error, solo un hash menos
    granular para ese stage.
    """
    stage_block = getattr(config, stage_name, None)
    charge = getattr(config, "charge", None)
    payload = {
        "charge": charge,
        "stage": stage_name,
        "config": _to_jsonable(stage_block) if stage_block is not None else None,
    }
    canonical = json.dumps(payload, sort_keys=True, default=str)
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


# ─────────────────────────────────────────────
#  FINGERPRINT DE OUTPUTS
# ─────────────────────────────────────────────

def _fingerprint_file(path: Path) -> str:
    """sha256 del contenido. Los outputs de stage relevantes para el
    fingerprint (adc2_results.json, crest_conformers.xyz, etc.) son chicos;
    para directorios completos (04_sp/ con .out de ORCA de GBs) se usa
    _fingerprint_dir en su lugar, sin leer contenido."""
    h = hashlib.sha256()
    try:
        with path.open("rb") as fh:
            for chunk in iter(lambda: fh.read(1 << 20), b""):
                h.update(chunk)
    except OSError:
        return "unreadable"
    return h.hexdigest()


def _fingerprint_dir(path: Path) -> str:
    """Fingerprint liviano de un directorio: nombre + tamaño de cada
    archivo, recursivo, SIN leer contenido (los .out de ORCA/CREST pueden
    ser grandes y esto se corre en cada resume). Suficiente para detectar
    "algo se borró/agregó/cambió de tamaño" sin el costo de un hash de
    contenido completo."""
    h = hashlib.sha256()
    try:
        entries = sorted(
            (p.relative_to(path).as_posix(), p.stat().st_size)
            for p in path.rglob("*") if p.is_file()
        )
    except OSError:
        return "unreadable"
    for rel, size in entries:
        h.update(f"{rel}:{size}".encode("utf-8"))
    return h.hexdigest()


def compute_output_fingerprint(paths: list[Path]) -> dict[str, str]:
    """Fingerprint por path de outputs()/inputs(). Directorios usan
    _fingerprint_dir (liviano); archivos usan _fingerprint_file (hash de
    contenido, son chicos: JSON de resultados, xyz de confórmeros, etc.)."""
    fp: dict[str, str] = {}
    for p in paths:
        if not p.exists():
            fp[str(p)] = "missing"
        elif p.is_dir():
            fp[str(p)] = _fingerprint_dir(p)
        else:
            fp[str(p)] = _fingerprint_file(p)
    return fp


# ─────────────────────────────────────────────
#  API DE RESUME
# ─────────────────────────────────────────────

def record_stage_result(
    workdir: Path, stage_name: str, config: Any, output_paths: list[Path],
    message: str = "",
) -> None:
    """Registra en state.json que `stage_name` terminó exitosamente con
    esta `config`. Llamar justo después de un StageResult(ok=True)."""
    state = load_state(workdir)
    state.stages[stage_name] = StageRecord(
        completed=True,
        config_hash=compute_config_hash(config, stage_name),
        output_fingerprint=compute_output_fingerprint(output_paths),
        timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"),
        message=message,
    )
    save_state(workdir, state)


def needs_rerun(
    workdir: Path, stage_name: str, config: Any, output_paths: list[Path],
) -> tuple[bool, str]:
    """True + motivo si `stage_name` necesita (re)correr, aunque
    stage.is_complete(workdir) diga True. Motivos posibles:

      - "sin registro previo"        : nunca se registró como completo.
      - "config cambió"              : config_hash no calza con el actual.
      - "outputs modificados en disco": el fingerprint no calza (archivos
                                         borrados/agregados/con otro tamaño
                                         desde el último resume).

    Si no necesita re-correr, devuelve (False, "").
    """
    state = load_state(workdir)
    record = state.stages.get(stage_name)

    if record is None or not record.completed:
        return True, "sin registro previo en state.json"

    current_hash = compute_config_hash(config, stage_name)
    if current_hash != record.config_hash:
        return True, "config cambió desde la última corrida completa"

    current_fp = compute_output_fingerprint(output_paths)
    if current_fp != record.output_fingerprint:
        return True, "outputs modificados en disco desde la última corrida completa"

    return False, ""


def resolve_resume_index(
    workdir: Path, stage_names: list[str], stages_by_name: dict, config: Any,
) -> int:
    """Índice (0-based) dentro de `stage_names` donde debe RETOMAR la
    ejecución. Generaliza el viejo --start: escanea la cadena en orden y
    se detiene en el primer stage que (a) no está is_complete() según su
    propia lógica, o (b) needs_rerun() según state.json (config cambió u
    outputs tocados por fuera del orchestrator).

    Devuelve len(stage_names) si TODA la cadena ya está al día (nada que
    correr).
    """
    for idx, name in enumerate(stage_names):
        stage = stages_by_name[name]
        if not stage.is_complete(workdir):
            return idx
        rerun, _reason = needs_rerun(workdir, name, config, stage.outputs(workdir))
        if rerun:
            return idx
    return len(stage_names)
