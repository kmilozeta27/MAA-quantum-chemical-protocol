"""Schema de validación (pydantic) para los recipes YAML de maaqc.

Los campos y sus defaults están tomados 1:1 de las constantes globales
y flags de argparse de pipeline.py (legacy), para que migrar un recipe
existente sea mecánico:

    NPROCS_CREST/NPROCS/NPROCS_NUMFREQ/NPROCS_TDDFT/NPROCS_ADC2
    MAXCORE/MAXCORE_NUMFREQ/MAXCORE_TDDFT/MAXCORE_ADC2
    OPT_BATCH_SIZE/SP_BATCH_SIZE/TDDFT_BATCH_SIZE
        -> pipeline.py líneas 180-192

    --temperature / --multitemp / --trange / --boltzmann-threshold
        -> pipeline.py líneas 4040-4082 (bloque "Modo multi-temperatura")

    EXP_LAMBDA_FALLBACK / molecules.yml
        -> pipeline.py líneas 222-259 (usado por stage_adc2)

TODO: una vez migrados los stages reales, cada Stage.run() debe recibir
un StageConfig (no el RecipeConfig completo) — ver contrato en
stages/base.py. Falta decidir si StageConfig se deriva de RecipeConfig
por sección (config.opt, config.sp, ...) o si cada stage define su
propio sub-modelo pydantic que RecipeConfig anida.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal, Optional

import yaml
from pydantic import BaseModel, Field, model_validator


class ResourceConfig(BaseModel):
    """Cores/memoria por etapa. Defaults = constantes actuales de pipeline.py."""

    nprocs_crest: int = 12
    nprocs: int = 4                 # OPT / SP / TDDFT
    nprocs_numfreq: int = 8
    nprocs_tddft: int = 6
    nprocs_adc2: int = 4

    maxcore: int = 3500              # MB/core, OPT / SP / TDDFT
    maxcore_numfreq: int = 2800
    maxcore_tddft: int = 2000
    maxcore_adc2: int = 5500

    opt_batch_size: int = 3
    sp_batch_size: int = 3
    tddft_batch_size: int = 2

    # Reemplaza --cores: si se define, redistribuye nprocs/maxcore
    # automáticamente vía system/affinity.py::ResourceProfile.detect().
    cores_total: Optional[int] = None
    max_ram_mb: int = 28000


class CensoConfig(BaseModel):
    """Config de CENSO prescreen, incluye modo multi-temperatura."""

    temperature: float = 298.15
    multitemp: bool = False
    # (T_ini, T_fin, N_puntos) — default histórico de --trange.
    trange: tuple[float, float, int] = (274.0, 298.0, 24)
    boltzmann_threshold: Optional[float] = None  # None -> select_nconf_auto
    ncon_min: int = 3
    ncon_max: int = 20
    timeout_hours: float = 24.0

    @model_validator(mode="after")
    def _validate_trange(self) -> "CensoConfig":
        t_ini, t_fin, n = self.trange
        if n < 1:
            raise ValueError("trange: N_puntos debe ser >= 1")
        if n == 1 and t_ini != t_fin:
            # Comportamiento legacy: N=1 usa solo T_INI.
            object.__setattr__(self, "trange", (t_ini, t_ini, 1))
        return self


class Adc2Config(BaseModel):
    mandatory: bool = False          # equivalente a --adc2 (opt-in explícito)
    basis: str = "def2-TZVP"
    exp_lambda_fallback: dict[str, float] = Field(default_factory=dict)
    molecules_yml: Optional[Path] = None  # override de _MOLECULES_YML


class TddftConfig(BaseModel):
    functional: str = "wb97x-d3"
    states: int = 10


class NrootConfig(BaseModel):
    enabled: bool = False
    top_n: int = 5


class RecipeConfig(BaseModel):
    """Recipe completo. Un preset YAML (quick/standard/publication) puebla esto."""

    extends: Optional[str] = None    # otro preset YAML como base (herencia simple)

    entry: Literal["smiles", "xyz", "crest-ensemble", "sp-only", "calibrated"] = "xyz"
    reference: Optional[str] = None  # ej. "shinorine" para --entry calibrated

    # Etiqueta del estado de protonación (ej. "A", "B", "R" — ver
    # EXP_LAMBDA_FALLBACK en stages/adc2.py). Antes, adc2.py leía este
    # campo vía getattr(config, "molecule_name", None) sin que existiera
    # en el schema: el override nunca podía tener efecto, siempre caía a
    # workdir.name. Ahora es un campo real; sigue siendo opcional porque
    # workdir.name es un fallback razonable cuando coincide con el label.
    molecule_name: Optional[str] = None

    skip_numfreq: bool = False       # reemplaza --skipnumfreq
    charge: Optional[int] = None
    multiplicity: int = 1            # MULT, siempre 1 en pipeline.py legacy

    resources: ResourceConfig = Field(default_factory=ResourceConfig)
    censo: CensoConfig = Field(default_factory=CensoConfig)
    tddft: TddftConfig = Field(default_factory=TddftConfig)
    adc2: Adc2Config = Field(default_factory=Adc2Config)
    nroot: NrootConfig = Field(default_factory=NrootConfig)


def _read_yaml_mapping(path: Path) -> dict[str, Any]:
    """Lee un YAML y exige que la raíz sea un mapping (no lista/escalar)."""
    try:
        raw = path.read_text()
    except OSError as e:
        raise FileNotFoundError(f"no se pudo leer el recipe: {path} ({e})") from e

    data = yaml.safe_load(raw)
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(
            f"{path}: la raíz de un recipe debe ser un mapping YAML "
            f"(clave: valor), no {type(data).__name__}"
        )
    return data


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Merge recursivo: `override` gana en cada clave, pero si ambos lados
    tienen un dict en la misma clave (ej. `censo:`, `tddft:`), se mergean
    campo a campo en vez de que `override` reemplace el bloque entero.
    Así `publication.yml` puede pisar solo `censo.multitemp` sin tener que
    repetir el resto de las claves de `censo` que ya trae `standard.yml`.

    La clave `extends` nunca sobrevive al merge: ya se usó para localizar
    al padre, y el `RecipeConfig` resultante siempre queda "aplanado".
    """
    merged = dict(base)
    for key, value in override.items():
        if key == "extends":
            continue
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _load_recipe_dict(path: Path, _chain: tuple[Path, ...] = ()) -> dict[str, Any]:
    """Resuelve `extends` recursivamente y devuelve el dict ya mergeado
    (sin validar todavía contra `RecipeConfig` — eso lo hace `load_recipe`
    una sola vez, al final de toda la cadena)."""
    path = path.resolve()
    if path in _chain:
        ciclo = " -> ".join(str(p) for p in (*_chain, path))
        raise ValueError(f"'extends' circular detectado: {ciclo}")

    own = _read_yaml_mapping(path)
    parent_name = own.get("extends")
    if parent_name is None:
        return own

    parent_path = (path.parent / parent_name).resolve()
    if not parent_path.exists():
        raise FileNotFoundError(
            f"{path}: extends '{parent_name}' no existe (buscado en {parent_path})"
        )

    parent_dict = _load_recipe_dict(parent_path, _chain=(*_chain, path))
    return _deep_merge(parent_dict, own)


def load_recipe(path: Path) -> RecipeConfig:
    """Carga y valida un recipe YAML, resolviendo `extends` recursivamente.

    - Sin `extends` (ej. `quick.yml`, `standard.yml`): se parsea y se
      valida tal cual contra `RecipeConfig`.
    - Con `extends: <otro.yml>` (ej. `publication.yml` -> `standard.yml`):
      se resuelve el padre primero (recursivo, soporta cadenas de más de
      un nivel), y se hace deep-merge campo a campo — el hijo solo
      necesita declarar lo que cambia respecto al padre.
    - `extends` se busca como ruta relativa al directorio del propio
      archivo (así los tres presets en `config/presets/` se referencian
      entre sí por nombre de archivo, sin rutas absolutas).
    - Herencia circular (A extends B extends A) levanta `ValueError` con
      la cadena completa, en vez de recursión infinita.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"recipe no encontrado: {path}")

    merged = _load_recipe_dict(path)
    merged.pop("extends", None)
    return RecipeConfig(**merged)
