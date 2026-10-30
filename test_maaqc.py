"""
Tests de maaqc.py — `pytest test_maaqc.py`

Cubren las piezas donde un error cuesta horas de cómputo o un pKa
equivocado: detección de éxito/convergencia de ORCA, manifiesto y
reanudación de jobs, reserva de CPUs, constantes del protón y cachés.
Nunca escriben en los runs/ reales: los outputs de ORCA se copian a
directorios temporales.
"""
from __future__ import annotations

import json
import logging
import math
import os
import shutil
from pathlib import Path

import pytest

import maaqc as m

ROOT = Path(__file__).resolve().parent
REAL_OPT   = ROOT / "AP/runs/APR/02_opt/conf_1/opt.out"
REAL_FREQ  = ROOT / "AP/runs/APR/03_numfreq/conf_1/freq.out"
REAL_ADC2  = ROOT / "MGLYVAL/runs/B/07_adc2/conf_006/adc2.out"

needs = lambda p: pytest.mark.skipif(not p.exists(), reason=f"falta {p} (sparse checkout)")
LOG = logging.getLogger("test_maaqc")


# ── carga por sufijo ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("name,charge", [
    ("SA", 1), ("SB", -1), ("SC", -2), ("SD", 0), ("MR", 0),
    ("AB2", -1), ("mol_a", 1), ("MGLYVAL", None), ("", None), ("123", None),
])
def test_infer_charge_from_name(name, charge):
    assert m.infer_charge_from_name(name) == charge


# ── éxito / convergencia / imaginarias ──────────────────────────────────────

@needs(REAL_OPT)
def test_check_success_and_convergence_real(tmp_path):
    assert m.check_success(REAL_OPT)
    assert m.check_opt_converged(REAL_OPT)
    trunc = tmp_path / "opt.out"
    trunc.write_bytes(REAL_OPT.read_bytes()[:-5000])
    assert not m.check_success(trunc)
    assert not m.check_opt_converged(trunc)


def test_opt_maxiter_is_not_converged(tmp_path):
    out = tmp_path / "opt.out"
    out.write_text("...\nThe optimization did not converge but reached the maximum "
                   "number of optimization cycles.\n****ORCA TERMINATED NORMALLY****\n")
    assert m.check_success(out)
    assert not m.check_opt_converged(out)


def test_check_success_missing(tmp_path):
    assert not m.check_success(tmp_path / "nope.out")


@needs(REAL_FREQ)
def test_check_no_imaginary_real():
    assert m.check_no_imaginary(REAL_FREQ) in (True, False)


def test_check_no_imaginary_synthetic(tmp_path):
    block = ("VIBRATIONAL FREQUENCIES\n-----\n"
             "   0:         0.00 cm**-1\n"
             "   6:      -123.45 cm**-1 ***imaginary mode***\n"
             "   7:       250.00 cm**-1\nNORMAL MODES\n")
    bad = tmp_path / "bad.out"
    bad.write_text(block)
    good = tmp_path / "good.out"
    good.write_text(block.replace("-123.45", "123.45"))
    assert not m.check_no_imaginary(bad)
    assert m.check_no_imaginary(good)
    assert not m.check_no_imaginary(tmp_path / "missing.out")   # antes: True


@needs(REAL_ADC2)
def test_extract_adc2_lambda_real():
    res = m.extract_adc2_lambda(REAL_ADC2)
    assert res is not None
    lam, ev, fosc = res
    assert 150 < lam < 600 and 2 < ev < 9 and fosc >= 0
    assert math.isclose(lam, 1239.84198 / ev, rel_tol=5e-3)


# ── constantes del protón dependientes de T ──────────────────────────────────

def test_proton_constants_exact_at_Tstd():
    assert m.g_hplus_gas(m.T_STD) == pytest.approx(m.G_HPLUS_GAS, abs=1e-12)
    assert m.dg_conc_corr(m.T_STD) == pytest.approx(m.DG_CONC_CORR, abs=1e-12)
    assert m.g_hplus_aq(m.T_STD) == pytest.approx(m.G_HPLUS_AQ_BORN, abs=1e-12)


def test_proton_constants_temperature_dependence():
    # δG_conc(T) = RT ln(R'T): valor analítico independiente del anclaje
    T = 250.0
    analytic = m.R_KCAL * T * math.log(m.R_LATM * T)
    anchored_offset = m.DG_CONC_CORR - m.R_KCAL * m.T_STD * math.log(m.R_LATM * m.T_STD)
    assert m.dg_conc_corr(T) == pytest.approx(analytic + anchored_offset, abs=1e-12)
    # H⁺ gas: dG/dT = −S(T) < 0, así que G sube al bajar T
    assert m.g_hplus_gas(250.0) > m.g_hplus_gas(274.15) > m.g_hplus_gas(298.15) > m.g_hplus_gas(310.0)
    assert m.g_hplus_gas(250.0) == pytest.approx(-5.05, abs=0.01)


def test_born_equals_direct_algebraically():
    """El ciclo de Born se reduce a la resta acuosa + G(H⁺,aq) a cualquier T."""
    H2K = m.HARTREE_TO_KCAL
    Eaq_a, Eaq_b, Egas_a, Egas_b, Gc_a, Gc_b = -500.1, -499.6, -500.0, -499.4, 0.21, 0.20
    for T in m.TEMPS_K:
        conc = m.dg_conc_corr(T)
        dG_gas = ((Egas_b + Gc_b) - (Egas_a + Gc_a)) * H2K + m.g_hplus_gas(T)
        born = (dG_gas + (Eaq_b - Egas_b) * H2K + conc + m.DG_SOLV_HPLUS + conc
                - ((Eaq_a - Egas_a) * H2K + conc))
        direct = ((Eaq_b + Gc_b) - (Eaq_a + Gc_a)) * H2K + m.g_hplus_aq(T)
        assert born == pytest.approx(direct, abs=1e-9)


# ── huella del job / manifiesto / archivo ────────────────────────────────────

SP_TMPL = m.TEMPLATE_SP


def _sp(charge=0, nprocs=4, maxcore=3000):
    return SP_TMPL.format(CHARGE=charge, MULT=1, NPROCS=nprocs, MAXCORE=maxcore)


def test_spec_ignores_resources_but_not_chemistry():
    assert m.job_spec_sha1(_sp(nprocs=4, maxcore=3000)) == m.job_spec_sha1(_sp(nprocs=8, maxcore=1000))
    assert m.job_spec_sha1(_sp(charge=0)) != m.job_spec_sha1(_sp(charge=-1))


def test_numfreq_restart_text_same_spec():
    fresh = m.TEMPLATE_NUMFREQ.format(CHARGE=-1, MULT=1, NPROCS_NUMFREQ=8, MAXCORE_NUMFREQ=2000)
    rst = m._numfreq_restart_text(fresh)
    assert "Restart true" in rst and rst.index("%freq") < rst.index("* xyzfile")
    assert m.job_spec_sha1(rst) == m.job_spec_sha1(fresh)


def test_tddft_guess_moread_does_not_change_spec():
    new = m.TEMPLATE_TDDFT.format(CHARGE=0, MULT=1, NPROCS=6, MAXCORE=2000, GBW_PATH="sp.gbw")
    assert "Guess   MORead" in new
    old = new.replace("  Guess   MORead\n", "")
    assert m.job_spec_sha1(new) == m.job_spec_sha1(old)


def _xyz(tmp, name="opt_opt.xyz", dx=0.0):
    p = tmp / name
    p.write_text(f"2\nE -1.0\nH 0.0 0.0 0.0\nH 0.0 0.0 {0.74 + dx:.4f}\n")
    return p


def _finish(job_dir, stem):
    (job_dir / f"{stem}.out").write_text("...\n****ORCA TERMINATED NORMALLY****\n")


def test_prepare_orca_job_lifecycle(tmp_path, monkeypatch, caplog):
    src = _xyz(tmp_path)
    job = tmp_path / "job"
    inp = m.prepare_orca_job(job, "sp.inp", _sp(0), LOG, "SP", xyz_src=src)
    assert inp.exists() and (job / "opt_opt.xyz").exists()
    man = m.read_manifest(job)
    assert man["status"] == "prepared" and man["spec_sha1"]
    _finish(job, "sp")

    # misma huella (solo cambian recursos) → se reutiliza, sin aviso
    caplog.clear()
    m.prepare_orca_job(job, "sp.inp", _sp(0, nprocs=8), LOG, "SP", xyz_src=src)
    assert "parámetros distintos" not in caplog.text

    # cambia la carga → aviso, se conserva el resultado
    monkeypatch.setattr(m, "RERUN_STALE", False)
    m.prepare_orca_job(job, "sp.inp", _sp(-1), LOG, "SP", xyz_src=src)
    assert "parámetros distintos" in caplog.text
    assert (job / "sp.out").exists()
    assert "stale_detected" in m.read_manifest(job)

    # con --rerun-stale → se archiva y se reescribe el input
    monkeypatch.setattr(m, "RERUN_STALE", True)
    m.prepare_orca_job(job, "sp.inp", _sp(-1), LOG, "SP", xyz_src=src)
    assert not (job / "sp.out").exists()
    archived = list((job / "_archive").glob("*_stale/sp.out"))
    assert archived, "el .out obsoleto debe quedar archivado, no borrado"
    assert "* xyzfile -1 1" in (job / "sp.inp").read_text()


def test_prepare_detects_changed_source_geometry(tmp_path, caplog, monkeypatch):
    monkeypatch.setattr(m, "RERUN_STALE", False)
    job = tmp_path / "job"
    m.prepare_orca_job(job, "sp.inp", _sp(0), LOG, "SP", xyz_src=_xyz(tmp_path))
    _finish(job, "sp")
    caplog.clear()
    other = _xyz(tmp_path, "other.xyz", dx=0.05)
    m.prepare_orca_job(job, "sp.inp", _sp(0), LOG, "SP", xyz_src=other)
    assert "parámetros distintos" in caplog.text


def test_prepare_backfills_legacy_job(tmp_path):
    """Un job anterior al manifiesto con el mismo input se acepta y se registra."""
    src = _xyz(tmp_path)
    job = tmp_path / "legacy"
    job.mkdir()
    shutil.copy(src, job / "opt_opt.xyz")
    (job / "sp.inp").write_text(_sp(0, nprocs=12))
    _finish(job, "sp")
    m.prepare_orca_job(job, "sp.inp", _sp(0), LOG, "SP", xyz_src=src)
    man = m.read_manifest(job)
    assert man.get("backfilled") and man["status"] == "ok"
    assert (job / "sp.inp").read_text() == _sp(0, nprocs=12)   # no se tocó


def test_archive_files_never_overwrites(tmp_path):
    f = tmp_path / "a.out"
    f.write_text("1")
    d1 = m.archive_files(tmp_path, [f], "failed")
    f.write_text("2")
    d2 = m.archive_files(tmp_path, [f], "failed")
    assert d1 != d2 and (d1 / "a.out").read_text() == "1" and (d2 / "a.out").read_text() == "2"


def test_atomic_write(tmp_path):
    p = tmp_path / "x.json"
    m.atomic_write(p, '{"a": 1}')
    assert json.loads(p.read_text()) == {"a": 1}
    assert not [q for q in tmp_path.iterdir() if q.name.endswith(".tmp")]


# ── reanudación de OPT ──────────────────────────────────────────────────────

def _opt_inp(coords="H 0 0 0\nH 0 0 0.9"):
    return m.TEMPLATE_OPT.format(CHARGE=0, MULT=1, NPROCS=4, MAXCORE=3000, COORDS=coords)


def test_opt_restart_uses_partial_geometry(tmp_path):
    canonical = _opt_inp()
    m.update_manifest(tmp_path, spec_sha1=m.job_spec_sha1(canonical))
    (tmp_path / "opt.xyz").write_text("2\ncycle 7\nH 0.0 0.0 0.0\nH 0.0 0.0 0.7412\n")
    (tmp_path / "opt.gbw").write_bytes(b"\0" * 16)
    txt = m._opt_restart_input(tmp_path, canonical, 2, LOG, "conf_1")
    assert "H 0.0 0.0 0.7412" in txt and "H 0 0 0.9" not in txt
    assert "! MORead" in txt and (tmp_path / "opt_restart.gbw").exists()
    assert txt.rstrip().endswith("*")


def test_opt_restart_refuses_other_job_or_wrong_natoms(tmp_path):
    (tmp_path / "opt.xyz").write_text("2\nc\nH 0 0 0\nH 0 0 0.74\n")
    m.update_manifest(tmp_path, spec_sha1="otro-job")
    assert m._opt_restart_input(tmp_path, _opt_inp(), 2, LOG, "c") is None
    m.update_manifest(tmp_path, spec_sha1=m.job_spec_sha1(_opt_inp()))
    assert m._opt_restart_input(tmp_path, _opt_inp(), 3, LOG, "c") is None


def test_resolve_opt_xyz_refreshes_stale_copy(tmp_path):
    canon = tmp_path / "opt_opt.xyz"
    canon.write_text("viejo")
    os.utime(canon, (1, 1))
    (tmp_path / "opt.xyz").write_text("nuevo")
    assert m.resolve_opt_xyz(tmp_path) == canon
    assert canon.read_text() == "nuevo"


# ── ensemble CONF# ──────────────────────────────────────────────────────────

def test_index_frames_by_conf_id_plain_and_extxyz(tmp_path):
    plain = tmp_path / "plain.xyz"
    plain.write_text("2\nCONF7 -1.0\nH 0 0 0\nH 0 0 0.74\n2\nCONF3 -0.9\nH 0 0 0\nH 0 0 0.75\n")
    idx = m.index_frames_by_conf_id(plain)
    assert set(idx) == {"CONF7", "CONF3"}
    assert idx["CONF3"].positions[1][2] == pytest.approx(0.75)
    ext = tmp_path / "ext.xyz"
    ext.write_text('2\nCONF5=T energy=-1.0\nH 0 0 0\nH 0 0 0.74\n')
    assert set(m.index_frames_by_conf_id(ext)) == {"CONF5"}


# ── CPUs: topología y lock entre procesos ───────────────────────────────────

def test_topology_covers_affinity():
    pus = set(os.sched_getaffinity(0))
    assert set(m.PU_TO_CORE) >= pus
    assert set(m.usable_pus()).isdisjoint(m.LP_CORE_SINGLES)


@pytest.mark.skipif(len(os.sched_getaffinity(0)) < 8,
                    reason="necesita ≥ 8 CPUs lógicas (p.ej. runner de CI con 4)")
def test_cpu_lock_is_transactional(tmp_path, monkeypatch):
    lock = tmp_path / "lock.json"
    monkeypatch.setattr(m, "CPU_LOCK_FILE", lock)
    # PID 1 vivo (init), PID absurdo muerto → se poda
    lock.write_text(json.dumps({"1": [4, 5], "999999999": [6, 7]}))
    lists = m._acquire_cpu_lists(1, 2, LOG, "T")
    data = json.loads(lock.read_text())
    assert "999999999" not in data
    assert lists and not (set(lists[0]) & {4, 5})
    assert data[str(os.getpid())] == sorted(lists[0])
    m._cpu_lock_release()
    assert str(os.getpid()) not in json.loads(lock.read_text())


# ── GoodVibes: caché ────────────────────────────────────────────────────────

def test_goodvibes_cache(tmp_path, monkeypatch):
    freq = tmp_path / "freq.out"
    freq.write_text("fake freq")
    calls = []

    def fake(freq_out, T, cutoff, qs, qh, python_exe, null):
        calls.append(T)
        return {**null, "ok": True, "G_corr": 0.2 + T * 1e-6, "scf_energy": -1.0}

    monkeypatch.setattr(m, "_run_goodvibes_uncached", fake)
    monkeypatch.setattr(m, "_gv_importable", lambda exe: True)
    monkeypatch.setattr(m, "_gv_version_str", lambda exe: "test")
    m._GV_MEM_CACHE.clear()
    a = m.run_goodvibes(freq, T=298.15)
    b = m.run_goodvibes(freq, T=298.15)
    assert a == b and calls == [298.15]
    m._GV_MEM_CACHE.clear()                     # "otro proceso": lee de disco
    m.run_goodvibes(freq, T=298.15)
    assert calls == [298.15]
    assert (tmp_path / "freq.out.gv_cache.json").exists()
    freq.write_text("freq recalculado")         # cambia el .out → invalida
    m.run_goodvibes(freq, T=298.15)
    assert len(calls) == 2


# ── logging ──────────────────────────────────────────────────────────────────

def test_setup_logging_distinct_per_path(tmp_path):
    a = m.setup_logging(tmp_path / "G1" / "runs" / "A")
    b = m.setup_logging(tmp_path / "G2" / "runs" / "A")
    assert a is not b
    b.info("solo-en-G2")
    for h in a.handlers + b.handlers:
        h.flush()
    assert "solo-en-G2" in (tmp_path / "G2/runs/A/pipeline.log").read_text()
    assert "solo-en-G2" not in (tmp_path / "G1/runs/A/pipeline.log").read_text()
    for lg in (a, b):
        for h in list(lg.handlers):
            h.close()
            lg.removeHandler(h)
