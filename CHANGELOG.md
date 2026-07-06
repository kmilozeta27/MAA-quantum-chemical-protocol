# Changelog

Formato basado en [Keep a Changelog](https://keepachangelog.com/es/1.0.0/).

## [Unreleased]

### Fixed
- `io/xyz.py` ya no tiene stubs (`NotImplementedError`): las 3 funciones
  de I/O de XYZ (`atoms_to_xyz_block`, `read_all_conformers`,
  `extract_conf_ids_from_xyz`) estaban duplicadas en
  `system/_orca_shared.py`. Ahora `io/xyz.py` es la fuente única de
  verdad y `_orca_shared.py` importa desde ahí (sin cambios para
  `opt.py`/`sp.py`/`crest.py`/`tddft.py`/`adc2.py`, que siguen usando
  `shared.atoms_to_xyz_block(...)` etc. vía re-export). Se descartó
  formalmente la firma legacy `read_conformers(xyz_file, nconf)` — ver
  `docs/stages.md` y el docstring de `io/xyz.py`. Agregado
  `tests/test_io/test_xyz.py` (usa `tests/fixtures/glycine.xyz`).

### Added
- Esqueleto inicial del paquete `maaqc` (src-layout, `pyproject.toml`,
  contrato `Stage` en `stages/base.py`).

### Migration status (legacy -> maaqc)
Ver `docs/stages.md` para detalle por stage. Estado actual (según
`maaqc doctor`): 5 de 7 stages migrados con lógica real —
`crest`, `opt`, `sp`, `tddft`, `adc2`. Quedan como stubs
(`NotImplementedError`): `censo`, `numfreq`. `pipeline.py` v0-legacy
sigue siendo la fuente de verdad para correr química real en esas dos
etapas hasta que tengan su contraparte migrada y verificada (outputs
idénticos, ver sección "Orden de migración" en
`docs/troubleshooting.md`).

## [0.1.0.dev0] - TODO
- Tag `v0-legacy` en git sobre `pipeline.py` / `results_collector.py` /
  `protonation_cycle.py` / `recover_prescreening.py` antes de empezar
  la migración, para no perder el estado que ya produjo resultados
  usados en PAPER1.
