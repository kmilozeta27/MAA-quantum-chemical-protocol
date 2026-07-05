"""Tests de OptStage.

stages/opt.py YA tiene lógica real (migrado) — este archivo estaba
marcado @skip con un mensaje que decía lo contrario y además esperaba
`requires == ["censo"]` (no calzaba con el código real: `["crest"]`).
Se corrigió lo esperado y se quitó el skip para el chequeo de contrato,
que no necesita ORCA ni fixture.

TODO (bloqueado en tests/fixtures/glycine.xyz, ver README de fixtures):
test de integración real que corra `OptStage.run()` contra la glicina
de test con ORCA real o mockeado, y compare outputs contra una corrida
de pipeline.py legacy sobre el mismo input — no solo contra valores
hardcodeados.
"""

from maaqc.stages.opt import OptStage


def test_opt_stage_contract():
    stage = OptStage()
    assert stage.name == "opt"
    assert stage.requires == ["crest"]
