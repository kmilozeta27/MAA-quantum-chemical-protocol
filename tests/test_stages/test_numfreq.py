"""Tests de NumFreqStage.

Mismo patrón que test_opt.py: chequeo de contrato (name/requires) que no
necesita ORCA ni fixture real.

TODO (bloqueado en tests/fixtures/glycine.xyz, ver README de fixtures):
test de integración real que corra `NumFreqStage.run()` sobre una
glicina ya optimizada (fixture con `02_opt/conf_001/opt_opt.xyz` +
`opt.out` con "ORCA TERMINATED NORMALLY") con ORCA real o mockeado, y
compare outputs contra una corrida de pipeline.py legacy sobre el mismo
input — no solo contra valores hardcodeados. También cubrir el camino de
restart: fixture con `freq_D00001.res.Gradients` presente para verificar
que se arma el input con `TEMPLATE_NUMFREQ_RESTART` en vez de
`TEMPLATE_NUMFREQ`.
"""

from maaqc.stages.numfreq import NumFreqStage


def test_numfreq_stage_contract():
    stage = NumFreqStage()
    assert stage.name == "numfreq"
    assert stage.requires == ["opt"]


def test_numfreq_opt_conf_dirs_empty_when_no_opt_dir(tmp_path):
    stage = NumFreqStage()
    assert stage._opt_conf_dirs(tmp_path) == []


def test_numfreq_is_complete_false_when_no_stage_dir(tmp_path):
    stage = NumFreqStage()
    assert stage.is_complete(tmp_path) is False
