"""Tests de OptStage. Skipeados hasta que stages/opt.py tenga lógica real.

Cuando se migre, correr esto contra tests/fixtures/ (molécula chica,
glicina, NO MGLYVAL) y comparar outputs contra una corrida del
pipeline.py legacy sobre el mismo input, no solo contra valores
hardcodeados en el test.
"""

import pytest

from maaqc.stages.opt import OptStage


@pytest.mark.skip(reason="OptStage.run() todavía no está migrado (ver stages/opt.py)")
def test_opt_stage_completes_on_glycine_fixture(tmp_path):
    stage = OptStage()
    assert stage.name == "opt"
    assert stage.requires == ["censo"]
