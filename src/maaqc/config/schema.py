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
from typing import Literal, Optional

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

    skip_numfreq: bool = False       # reemplaza --skipnumfreq
    charge: Optional[int] = None
    multiplicity: int = 1            # MULT, siempre 1 en pipeline.py legacy

    resources: ResourceConfig = Field(default_factory=ResourceConfig)
    censo: CensoConfig = Field(default_factory=CensoConfig)
    tddft: TddftConfig = Field(default_factory=TddftConfig)
    adc2: Adc2Config = Field(default_factory=Adc2Config)
    nroot: NrootConfig = Field(default_factory=NrootConfig)


def load_recipe(path: Path) -> RecipeConfig:
    """Carga y valida un recipe YAML, resolviendo `extends` recursivamente.

    TODO: implementar merge de `extends` (deep-merge de dicts antes de
    validar con pydantic) una vez que exista más de un preset real que
    lo use. Por ahora los tres presets en config/presets/ son standalone.
    """
    raise NotImplementedError("load_recipe: pendiente merge de `extends` + parseo YAML")
