# maaqc — MAA Quantum Characterization

Paquete instalable para el pipeline de caracterización cuántica de MAAs
(mycosporine-like amino acids) usado en PAPER1 (Marine Drugs):

```
CREST -> CENSO (prescreen) -> OPT -> NumFreq -> SP -> TDDFT -> [ADC2]
```

Reemplaza los scripts monolíticos `pipeline.py` / `results_collector.py` /
`protonation_cycle.py` por un paquete con la misma "sensación" que
CREST/CENSO/ORCA: un único CLI de subcomandos, recipes declarativos en
YAML, y una interfaz `Stage` común para que el resume sea robusto sin
depender de flags como `--start adc2` o `--skipnumfreq`.

> **Estado:** en migración activa desde el pipeline legacy (tag
> `v0-legacy` en git). Ver `docs/` para el estado real de cada stage y
> `CHANGELOG.md` para el progreso de la migración.

## Instalación

```bash
conda env create -f environment.yml
conda activate maaqc
pip install -e ".[dev]"
maaqc doctor   # verifica ORCA / CREST / CENSO en PATH o env conda
```

## Uso básico

```bash
maaqc run       molecule.xyz --config recipe_publication.yml
maaqc resume    --from sp
maaqc collect   --states A,B,C,D,R
maaqc report    --format html
maaqc doctor
```

Ver `docs/config_reference.md` para cada keyword de configuración y
`docs/stages.md` para inputs/outputs de cada stage.

## Citar este software

Ver `CITATION.cff`. Sigue la convención de citar CREST/CENSO/ORCA con
paper + versión (pendiente: DOI de Zenodo tras el primer release).

## Estructura del repo

Ver `docs/stages.md` y el árbol en `src/maaqc/` — sigue un src-layout
estándar con `stages/`, `orchestrator/`, `analysis/`, `system/`, `io/`
y `config/`.
