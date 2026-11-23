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
    # adc2.out real de MGLYVAL/B/07_adc2/conf_006 fue cortado por timeout
    # (sin TERMINATED NORMALLY): el extractor no debe fiarse de él.
    assert not m.check_success(REAL_ADC2)
    assert m.extract_adc2_lambda(REAL_ADC2) is None
    assert m.parse_adc2_excitations(REAL_ADC2) == []


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


# ═════════════════════════════════════════════════════════════════════════════
#  P0-2 · P0-3 · P0-5 · P0-7  (fixtures sintéticos en tmp_path)
# ═════════════════════════════════════════════════════════════════════════════

import numpy as np

_ABS_BLOCK = (
    "-----\n ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS\n-----\n"
    "     Transition      Energy     Energy  Wavelength fosc(D2)\n-----\n"
    "  0-1A  ->  1-1A    {ev:.6f}   28000.0   {nm:.1f}   0.400000000   4.0  0.1  0.1  0.1\n"
    "-----\n ABSORPTION SPECTRUM VIA TRANSITION VELOCITY DIPOLE MOMENTS\n"
)


def _mk_sp(path: Path, e: float):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f"FINAL SINGLE POINT ENERGY      {e:.8f}\n****ORCA TERMINATED NORMALLY****\n")


# ── P0-2: ids numéricos de conformero ────────────────────────────────────────

def test_conf_num_numeric_ids():
    assert m._conf_num("conf_5") == m._conf_num("conf_005") == 5
    assert m._conf_num("conf_6_PBE0") is None
    assert m._conf_num("conf6") is None and m._conf_num("sel") is None


def test_collect_tddft_includes_all_conformer_dir_names(tmp_path, monkeypatch):
    """conf_5/conf_6 (hechos a mano) deben cruzarse con conf_005/conf_006 de
    04_sp y con sus pesos de Boltzmann; conf_6_PBE0 no es parte del ensemble."""
    mol = tmp_path / "X"
    energies = {"conf_001": -100.000, "conf_004": -99.999,
                "conf_005": -99.998, "conf_006": -100.002}   # 006 = dominante
    for cid, e in energies.items():
        _mk_sp(mol / "04_sp" / cid / "sp.out", e)
    _mk_sp(mol / "04_sp" / "conf_6_PBE0" / "sp.out", -100.5)   # no debe entrar
    for cid, tdir, ev in (("conf_001", "conf_001", 3.5), ("conf_004", "conf_004", 3.6),
                          ("conf_005", "conf_5", 3.7), ("conf_006", "conf_6", 3.4),
                          (None, "conf_6_PBE0", 4.2)):
        d = mol / "05_tddft" / tdir
        d.mkdir(parents=True)
        (d / "tddft.out").write_text(_ABS_BLOCK.format(ev=ev, nm=1239.84198 / ev))
    monkeypatch.setattr(m, "RC_WORK_ROOT", tmp_path)
    with pytest.warns(UserWarning):          # sin NumFreq: Boltzmann solo electrónico
        df_conf, df_tddft = m.collect_results("X")
    assert sorted(df_conf["conformer"]) == sorted(energies)
    assert "conf_6_PBE0" not in set(df_conf["conformer"])
    s1 = df_tddft[df_tddft["state"] == 1]
    assert sorted(s1["conformer"]) == sorted(energies)          # 4/4, no 2/4
    assert not s1["E_eV"].eq(4.2).any()
    w = df_conf.set_index("conformer")["boltz_298K"]
    assert w.sum() == pytest.approx(1.0)
    assert w.idxmax() == "conf_006"
    assert df_conf.attrs["lambda_max_boltz_298K"] > 0   # el cruce con pesos no pierde conf_5/6


def test_ensemble_conf_dirs_duplicates_prefer_file(tmp_path):
    (tmp_path / "conf_6").mkdir()
    d = tmp_path / "conf_006"
    d.mkdir()
    (d / "tddft.out").write_text("x")
    with pytest.warns(UserWarning):
        got = m._ensemble_conf_dirs(tmp_path, "tddft.out")
    assert [g.name for g in got] == ["conf_006"]


# ── P0-3: CENSO ausente aborta ───────────────────────────────────────────────

def _patch_censo_missing(monkeypatch):
    def boom(*a, **k):
        raise FileNotFoundError("censo")
    monkeypatch.setattr(m, "_acquire_cpu_pool", lambda *a, **k: [0, 1])
    monkeypatch.setattr(m, "_reconcile_censo_budget", lambda *a, **k: (2, 1))
    monkeypatch.setattr(m, "_cpu_lock_release", lambda: None)
    monkeypatch.setattr(m, "_run_censo_locked", boom)


def test_censo_missing_aborts(tmp_path, monkeypatch):
    xyz = tmp_path / "crest_conformers.xyz"
    xyz.write_text("1\nE=-1.0\nH 0 0 0\n")
    _patch_censo_missing(monkeypatch)
    monkeypatch.setattr(m, "ALLOW_NO_CENSO", False)
    with pytest.raises(m.CensoUnavailableError):
        m.stage_censo_prescreen(xyz, tmp_path / "run", 0, LOG)
    # con el flag explícito se recupera el comportamiento antiguo (con aviso)
    monkeypatch.setattr(m, "ALLOW_NO_CENSO", True)
    assert m.stage_censo_prescreen(xyz, tmp_path / "run", 0, LOG) == xyz


# ── P0-5: calibración obsoleta fuera del camino por defecto ──────────────────

def _fake_calib(**kw):
    return {"amino": [
        {"mol": "A", "calc_harmonic": 10.0, "calc_goodvibes": 10.0, "exp": 9.0, "exp_ref": "t"},
        {"mol": "B", "calc_harmonic": 12.0, "calc_goodvibes": 12.0, "exp": 9.5, "exp_ref": "t"},
    ]}


def test_calibration_file_overrides_builtin(tmp_path, monkeypatch):
    # prueba la mecánica del ajuste con n=2 (bajo CALIB_MIN_N): override explícito
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_below_min", True)
    import copy
    monkeypatch.setattr(m, "EMPIRICAL_CALIBRATION", copy.deepcopy(m.EMPIRICAL_CALIBRATION))
    before = m.fit_empirical_correction("amino", mode="offset")["intercept"]
    f = tmp_path / "cal.json"
    f.write_text(json.dumps({"meta": {"T": 298.15}, "EMPIRICAL_CALIBRATION": {
        "amino": [{"mol": "Z1", "calc_goodvibes": 20.0, "calc_harmonic": 20.0, "exp": 10.0},
                  {"mol": "Z2", "calc_goodvibes": 22.0, "calc_harmonic": 22.0, "exp": 12.0}]}}))
    meta = m.load_calibration_file(f)
    assert meta["T"] == 298.15
    fit = m.fit_empirical_correction("amino", mode="offset")
    assert fit["molecules"] == ["Z1", "Z2"] and fit["intercept"] == pytest.approx(-10.0)
    assert fit["intercept"] != pytest.approx(before)
    # los tipos que el archivo no trae se conservan
    assert "carboxyl" in m.EMPIRICAL_CALIBRATION


def test_empirical_correction_requires_calibration_file(tmp_path):
    with pytest.raises(SystemExit, match="calibration"):
        m.resolve_calibration_source(None, "pKa0=carboxyl", search_roots=[tmp_path])
    # archivo explícito gana; sin --empirical-correction no hace falta nada
    assert m.resolve_calibration_source("x.json", "pKa0=carboxyl") == Path("x.json")
    assert m.resolve_calibration_source(None, None, search_roots=[tmp_path]) is None
    # escape explícito
    assert m.resolve_calibration_source(None, "pKa0=carboxyl", allow_builtin=True,
                                        search_roots=[tmp_path]) is None
    # calibration_v4.json por defecto, si existe
    (tmp_path / "calibration").mkdir()
    (tmp_path / "calibration" / "calibration_v4.json").write_text("{}")
    got = m.resolve_calibration_source(None, "pKa0=carboxyl", search_roots=[tmp_path])
    assert got == tmp_path / "calibration" / "calibration_v4.json"


def test_refuses_cross_class_without_force(monkeypatch):
    # prueba la mecánica del ajuste con n=2 (bajo CALIB_MIN_N): override explícito
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_below_min", True)
    import copy
    monkeypatch.setattr(m, "EMPIRICAL_CALIBRATION", copy.deepcopy(_fake_calib()))
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "force_class", False)
    # 'amino' es 0→−1; el tramo a corregir es B→C (−1→−2)
    res = m.apply_empirical_correction(15.0, "amino", mode="offset",
                                       pka_result={"q_acid": -1, "q_base": -2})
    assert not res["ok"] and "clase de carga" in res["reason"]
    # misma clase → se aplica
    ok = m.apply_empirical_correction(15.0, "amino", mode="offset",
                                      pka_result={"q_acid": 0, "q_base": -1})
    assert ok["ok"]
    # forzado
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "force_class", True)
    forced = m.apply_empirical_correction(15.0, "amino", mode="offset",
                                          pka_result={"q_acid": -1, "q_base": -2})
    assert forced["ok"]


# ── P0-7: ajuste empírico ────────────────────────────────────────────────────

def _calib_xy(monkeypatch, pts):
    monkeypatch.setattr(m, "EMPIRICAL_CALIBRATION", {"t": [
        {"mol": f"M{i}", "calc_goodvibes": x, "calc_harmonic": x, "exp": y}
        for i, (x, y) in enumerate(pts)]})


def test_fit_empirical_offset_n1_rmse_nan(monkeypatch):
    _calib_xy(monkeypatch, [(11.0, 9.0)])
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_n1", False)
    assert not m.fit_empirical_correction("t")["ok"]          # n=1 rechazado por defecto
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_n1", True)
    fit = m.fit_empirical_correction("t")
    assert fit["ok"] and fit["n"] == 1 and fit["mode"] == "offset"
    assert math.isnan(fit["rmse"]) and math.isnan(fit["loo_rmse"])
    assert fit["slope"] == 1.0 and fit["intercept"] == pytest.approx(-2.0)


def test_fit_auto_switches_to_linear_n3_span1(monkeypatch):
    # prueba la mecánica del ajuste con n=2 (bajo CALIB_MIN_N): override explícito
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_below_min", True)
    _calib_xy(monkeypatch, [(2.0, 1.0), (3.0, 2.2), (4.0, 2.9)])    # rango = 2 ≥ 1
    assert m.fit_empirical_correction("t")["mode"] == "linear"
    _calib_xy(monkeypatch, [(2.0, 1.0), (2.4, 1.5), (2.8, 2.1)])    # rango = 0.8 < 1
    fit = m.fit_empirical_correction("t")
    assert fit["mode"] == "offset" and fit["slope"] == 1.0
    assert any("rango" in w for w in fit["warnings"])
    _calib_xy(monkeypatch, [(2.0, 1.0), (3.0, 2.2)])                # n=2 → offset
    assert m.fit_empirical_correction("t")["mode"] == "offset"


def test_loo_rmse_two_points(monkeypatch):
    # prueba la mecánica del ajuste con n=2 (bajo CALIB_MIN_N): override explícito
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_below_min", True)
    _calib_xy(monkeypatch, [(10.0, 8.0), (12.0, 9.0)])   # sesgos −2 y −3
    fit = m.fit_empirical_correction("t", mode="offset")
    assert fit["intercept"] == pytest.approx(-2.5)
    assert fit["loo_rmse"] == pytest.approx(1.0)         # |b0 − b1|
    assert fit["rmse"] == pytest.approx(0.5)


# ── P0-7: especiación ────────────────────────────────────────────────────────

def test_speciation_fractions_sum_to_one_and_monotonic():
    phs = [float(x) for x in np.arange(0, 14.5, 0.5)]
    d = m.boltzmann_protonation_states(3.0, 9.0, phs)
    tot = d["f_neutral"] + d["f_anion"] + d["f_dianion"]
    assert np.allclose(tot, 1.0)
    assert d["f_neutral"].is_monotonic_decreasing
    assert d["f_dianion"].is_monotonic_increasing
    d5 = m.boltzmann_5state(2.0, 5.0, 9.0, 0.3, 0.7, phs)
    tot5 = d5["f_cation"] + d5["f_R"] + d5["f_D"] + d5["f_anion"] + d5["f_dianion"]
    assert np.allclose(tot5, 1.0)
    assert d5["f_cation"].is_monotonic_decreasing
    assert d5["f_dianion"].is_monotonic_increasing
    # pH = pKa → mitad/mitad en el par de dos estados
    d1 = m.boltzmann_protonation_states(7.0, None, [7.0])
    assert d1["f_neutral"][0] == pytest.approx(0.5)


# ── P0-7: Monte Carlo ────────────────────────────────────────────────────────

def _fake_compute_pka(monkeypatch):
    monkeypatch.setattr(
        m, "compute_pka",
        lambda *a, **k: {"ok": True, "dG_aq_kcal": 13.0, "method": "direct", "warnings": []})


def test_mc_reproducible_with_seed(monkeypatch):
    _fake_compute_pka(monkeypatch)
    kw = dict(data_A={"a": 1}, data_neutral_macro={"n": 1}, data_B={"b": 1},
              data_C={"c": 1}, f_R=1.0, f_D=0.0, pH_values=[7.0, 8.1], n_samples=500)
    r1 = m.run_chain_montecarlo(seed=7, **kw)
    r2 = m.run_chain_montecarlo(seed=7, **kw)
    r3 = m.run_chain_montecarlo(seed=8, **kw)
    assert r1["ok"]
    for k in ("mc_pka0", "mc_pka1", "mc_pka2"):
        assert np.array_equal(r1[k]["pka_samples"], r2[k]["pka_samples"])
    assert not np.array_equal(r1["mc_pka1"]["pka_samples"], r3["mc_pka1"]["pka_samples"])
    assert r1["df_dist_bands"].equals(r2["df_dist_bands"])


def test_mc_shared_proton_draw_reduces_dpKa_variance(monkeypatch):
    _fake_compute_pka(monkeypatch)
    n = 4000
    rng = np.random.default_rng(1)
    shared = rng.normal(0.0, 2.5, size=n)

    def two(shared_draw):
        a = m.compute_pka_montecarlo({}, {}, n_samples=n, shared_hplus_draw=shared_draw,
                                     rng=np.random.default_rng(11))
        b = m.compute_pka_montecarlo({}, {}, n_samples=n, shared_hplus_draw=shared_draw,
                                     rng=np.random.default_rng(12))
        return a["pka_samples"] - b["pka_samples"]

    var_shared = np.var(two(shared))
    # sin compartir: cada tramo sortea su propio protón (rng distinto)
    a = m.compute_pka_montecarlo({}, {}, n_samples=n, rng=np.random.default_rng(21))
    b = m.compute_pka_montecarlo({}, {}, n_samples=n, rng=np.random.default_rng(22))
    var_indep = np.var(a["pka_samples"] - b["pka_samples"])
    assert var_shared < 0.25 * var_indep   # teórico ≈ 0.12 (σ_H⁺=2.5 vs σ_e,g combinado)


# ── P0-7: Boltzmann multi-T ──────────────────────────────────────────────────

def test_boltzmann_weights_multiT_sum_one():
    import pandas as pd
    e = [-100.000, -99.999, -99.997]
    df = pd.DataFrame({"conformer": ["conf_001", "conf_002", "conf_003"], "E_hartree": e})
    with pytest.warns(UserWarning):
        out = m.compute_boltzmann_multiT(df)             # modo electrónico
    assert out.attrs["boltzmann_mode"] == "electronic_only"
    for lbl in m.TEMP_LABELS:
        assert out[f"boltz_{lbl}"].sum() == pytest.approx(1.0)
    # más frío → más peso en el mínimo
    assert out["boltz_250K"][0] > out["boltz_310K"][0]
    # modo Gibbs compuesto
    df2 = df.copy()
    df2["thermo_ok"] = True
    for T, lbl in zip(m.TEMP_LABELS, m.TEMP_LABELS):
        g = np.array(e) + np.array([0.0, 0.0005, 0.0002])
        df2[f"G_total_{lbl}"] = g
        df2[f"dG_{lbl}_kcal"] = (g - g.min()) * m.HARTREE_TO_KCAL
    out2 = m.compute_boltzmann_multiT(df2)
    assert out2.attrs["boltzmann_mode"] == "gibbs_composite"
    for lbl in m.TEMP_LABELS:
        assert out2[f"boltz_{lbl}"].sum() == pytest.approx(1.0)
    assert out2["boltz_298K"].equals(out2["boltz_298K"])


# ── P0-7: NTO, Multiwfn, extractores ORCA ────────────────────────────────────

def test_nto_selection_dual_criterion():
    st = m.ExcitedState(state_idx=1, ntos=[
        m.NTOPair(10, 11, 0.95), m.NTOPair(9, 12, 0.40),
        m.NTOPair(8, 13, 0.06), m.NTOPair(7, 14, 0.02)])
    # relativo (≥ 0.1·0.95 = 0.095) Y absoluto (≥ min_coef)
    sel = m.select_ntos_for_state(st, threshold=0.1, min_coef=0.05, max_pairs=None)
    assert [p.coef for p in sel] == [0.95, 0.40]             # 0.06 falla el relativo
    sel = m.select_ntos_for_state(st, threshold=0.01, min_coef=0.05, max_pairs=None)
    assert [p.coef for p in sel] == [0.95, 0.40, 0.06]       # 0.02 falla el absoluto
    sel = m.select_ntos_for_state(st, threshold=0.0, min_coef=0.0, max_pairs=2)
    assert len(sel) == 2
    assert m.select_ntos_for_state(m.ExcitedState(2), 0.1, 0.05, None) == []


def test_parse_multiwfn_esp_log():
    log = (
        " Global surface minimum: -0.151 a.u.\n"
        " Minimal value:   -0.15127 a.u. (  -94.92 kcal/mol )\n"
        " Maximal value:    0.12790 a.u. (   80.26 kcal/mol )\n"
        " Overall surface area:      412.5 Bohr^2  (  115.5 Angstrom^2)\n"
        " Positive surface area:     200.0 Bohr^2  (   56.0 Angstrom^2)\n"
        " Negative surface area:     212.5 Bohr^2  (   59.5 Angstrom^2)\n"
        " Overall average value:   0.01000 a.u. (    6.28 kcal/mol )\n"
    )
    r = m._parse_multiwfn_esp_log(log)
    assert r["ok"]
    assert r["ESP_min_kcal"] == pytest.approx(-94.92)
    assert r["ESP_max_kcal"] == pytest.approx(80.26)
    assert r["area_total_ang2"] == pytest.approx(115.5)
    assert r["avg_kcal"] == pytest.approx(6.28)
    bad = m._parse_multiwfn_esp_log("nada útil")
    assert not bad["ok"] and bad["warning"]


REAL_STEOM = ROOT / "MGLYVAL/runs/B/07_adc2/conf_6/CCSD_CPCM.out"
REAL_MATZVP = ROOT / "MGLYVAL/runs/B/07_adc2/conf_6/CCSD_maTZVP.out"


@needs(REAL_STEOM)
def test_extract_steom_vs_cis_spectrum():
    """El output STEOM trae primero el espectro CIS de referencia (4.943 eV);
    el extractor debe quedarse con el bloque RIGHT de STEOM (3.549 eV)."""
    res = m.extract_adc2_lambda(REAL_STEOM)
    assert res is not None
    lam, ev, fosc = res
    assert ev == pytest.approx(3.549, abs=1e-3)
    assert not math.isclose(ev, 4.943, abs_tol=0.1)
    assert lam == pytest.approx(349.4, abs=0.1)


@needs(REAL_MATZVP)
def test_extract_steom_matzvp_failed_returns_none():
    assert not m.check_success(REAL_MATZVP)
    assert m.extract_adc2_lambda(REAL_MATZVP) is None


@needs(REAL_ADC2)
def test_adc2_requires_normal_termination(tmp_path):
    # real: cortado por timeout → None / []
    assert m.extract_adc2_lambda(REAL_ADC2) is None
    assert m.parse_adc2_excitations(REAL_ADC2) == []
    # sintético: mismo bloque RIGHT con y sin terminación normal
    body = ("SPECTRUM FOR RIGHT TRANSITION MOMENTS\n" + _ABS_BLOCK.format(ev=3.46, nm=357.8))
    cut = tmp_path / "cut.out"
    cut.write_text(body)
    done = tmp_path / "done.out"
    done.write_text(body + "\n   ****ORCA TERMINATED NORMALLY****\n")
    assert m.extract_adc2_lambda(cut) is None
    assert m.extract_adc2_lambda(done)[1] == pytest.approx(3.46)
    assert m.parse_adc2_excitations(cut) == []
    assert len(m.parse_adc2_excitations(done)) == 1


# ── P1-11 / P1-12: calibración, N microestados, MC calibrado ─────────────────

def test_pka_exp_single_source():
    # Asp: IUPAC Reliable (1.99/3.90/9.84), no McMurry (1.88/3.65/9.60)
    assert m.PKA_EXP["Aspartic"] == {"carboxyl": 1.99, "carboxyl_2": 3.90, "amino": 9.84}
    asp = {p["mol"]: p["exp"] for p in m.EMPIRICAL_CALIBRATION["carboxyl"]}
    assert asp["Aspartic"] == m.PKA_EXP["Aspartic"]["carboxyl"]
    assert m.EMPIRICAL_CALIBRATION["carboxyl_2"][0]["exp"] == m.PKA_EXP["Aspartic"]["carboxyl_2"]
    tr = {(mol, t["type"]): t["exp"] for mol, spec in m.CALIBRATION_SETS.items()
          for t in spec["transitions"]}
    assert tr[("ASPARTIC", "carboxyl")] == 1.99 and tr[("ASPARTIC", "carboxyl_2")] == 3.90
    # toda 'exp' de CALIBRATION_SETS coincide con PKA_EXP
    for mol, spec in m.CALIBRATION_SETS.items():
        for t in spec["transitions"]:
            assert t["exp"] == m.PKA_EXP[spec["label"]][t["type"]]


def test_loo_linear(monkeypatch):
    # prueba la mecánica del ajuste con n=2 (bajo CALIB_MIN_N): override explícito
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_below_min", True)
    pts = [(2.0, 1.0), (3.0, 2.3), (4.0, 2.8), (6.0, 5.1)]
    _calib_xy(monkeypatch, pts)
    fit = m.fit_empirical_correction("t", mode="linear")
    assert fit["mode"] == "linear"
    # LOO manual con numpy
    res = []
    for k in range(len(pts)):
        xs = [p[0] for i, p in enumerate(pts) if i != k]
        ys = [p[1] for i, p in enumerate(pts) if i != k]
        sl, ic = np.polyfit(xs, ys, 1)
        res.append(pts[k][1] - (sl * pts[k][0] + ic))
    assert fit["loo_rmse"] == pytest.approx(float(np.sqrt(np.mean(np.square(res)))))
    assert fit["loo_rmse"] > fit["rmse"]
    # n=2 lineal: sin grados de libertad para el LOO
    _calib_xy(monkeypatch, pts[:2])
    assert math.isnan(m.fit_empirical_correction("t", mode="linear")["loo_rmse"])


def test_aicc_selects_offset_vs_linear(monkeypatch):
    # sesgo constante + ruido pequeño, n=6, rango ≥3 → offset gana
    rng = np.random.default_rng(0)
    xs = [2, 3.5, 5, 7, 9, 11]
    _calib_xy(monkeypatch, [(x, x - 1.5 + rng.normal(0, 0.05)) for x in xs])
    fit = m.fit_empirical_correction("t", mode="aicc")
    assert fit["mode"] == "offset" and fit["aicc_offset"] < fit["aicc_linear"] + 2
    # pendiente real 0.7 → lineal gana
    _calib_xy(monkeypatch, [(x, 0.7 * x + 0.5 + rng.normal(0, 0.05)) for x in xs])
    assert m.fit_empirical_correction("t", mode="aicc")["mode"] == "linear"
    # n<5 → offset siempre, AICc lineal indefinido
    _calib_xy(monkeypatch, [(2, 1), (5, 3), (8, 4.5), (11, 7)])
    f4 = m.fit_empirical_correction("t", mode="aicc")
    assert f4["mode"] == "offset" and math.isnan(f4["aicc_linear"])


def _micro(name, G, q=0, n=10):
    return {"mol_name": name, "is_macrostate": True, "G_total_direct_Eh": G,
            "charge": q, "n_atoms": n, "use_smd": False}


def test_macrostate_n_micro_degeneracy(tmp_path):
    T = m.T_STD
    kT = m.R_KCAL * T / m.HARTREE_TO_KCAL      # Eh
    # N=3 microestados degenerados en G, g = 1, 1, 2 → G_macro = G − kT·ln 4
    G0 = -100.0
    micro = [_micro("R", G0), _micro("D", G0), _micro("Dp", G0)]
    res = m.compute_micro_equilibrium(micro, degeneracy={"Dp": 2}, T=T)
    assert res["ok"]
    assert res["G_macro_Eh"] == pytest.approx(G0 - kT * math.log(4.0))
    assert res["fractions"]["Dp"] == pytest.approx(0.5)
    assert sum(res["fractions"].values()) == pytest.approx(1.0)
    # sin degeneración y N=2 coincide con el macroestado histórico R/D
    dG = 1.0 / m.HARTREE_TO_KCAL
    two = m.compute_micro_equilibrium([_micro("R", G0), _micro("D", G0 + dG)], T=T)
    pair = m.compute_tautomer_equilibrium(_micro("R", G0), _micro("D", G0 + dG), T=T)
    assert two["G_macro_Eh"] == pytest.approx(pair["G_macro_Eh"])
    # carga inconsistente → ok=False, sin fallback silencioso
    bad = m.compute_micro_equilibrium([_micro("R", G0), _micro("X", G0, q=-1)], T=T)
    assert not bad["ok"]
    # _build_neutral_macrostate con extra_micro + degeneración
    taut, macro = m._build_neutral_macrostate(
        _micro("R", G0), _micro("D", G0), "R", "D", T, verbose=False,
        extra_micro=[_micro("Dp", G0)], degeneracy={"Dp": 2})
    assert macro["G_total_direct_Eh"] == pytest.approx(G0 - kT * math.log(4.0))
    assert macro["mol_name"] == "R/D/Dp"
    assert taut["f_R"] == pytest.approx(0.5)      # f_R/f_D siguen siendo los del par R/D


def test_mc_calibrated_no_hplus(monkeypatch):
    _fake_compute_pka(monkeypatch)                # dG=13 kcal/mol → pKa crudo ≈ 9.53
    pts = [(10.0, 8.0), (12.0, 9.0), (11.0, 8.6)]
    _calib_xy(monkeypatch, pts)
    fit = m.fit_empirical_correction("t", mode="offset")
    fit["group_type"] = "t"
    n = 20000
    kw = dict(n_samples=n)
    draw = m._calibration_intercept_draw(fit, n, np.random.default_rng(3))
    # σ_H⁺ enorme y shared_hplus_draw no tienen efecto en modo calibrado
    a = m.compute_pka_montecarlo({}, {}, sigma_hplus_kcal=50.0, calibration_fit=fit,
                                 calib_intercept_draw=draw, rng=np.random.default_rng(1), **kw)
    b = m.compute_pka_montecarlo({}, {}, sigma_hplus_kcal=0.0,
                                 shared_hplus_draw=np.full(n, 99.0), calibration_fit=fit,
                                 calib_intercept_draw=draw, rng=np.random.default_rng(2), **kw)
    assert a["ok"] and a["hplus_draw"] is None
    assert np.array_equal(a["pka_samples"], b["pka_samples"])
    pka_raw = 13.0 / (m.R_KCAL * m.T_STD * m.LN10)
    assert a["pka_mean"] == pytest.approx(pka_raw + fit["intercept"], abs=0.05)
    nn = fit["n"]
    expected = fit["loo_rmse"] * math.sqrt(1 + 1 / nn)
    # t con 2 g.l. tiene varianza infinita: la mediana del |desvío| es la medida robusta
    mad = float(np.median(np.abs(a["pka_samples"] - a["pka_median"])))
    assert 0.4 * expected < mad < 1.6 * expected
    # run_chain_montecarlo: un solo sorteo por clase, compartido entre tramos
    r = m.run_chain_montecarlo(
        data_A={"a": 1}, data_neutral_macro={"n": 1}, data_B={"b": 1}, data_C=None,
        f_R=1.0, f_D=0.0, pH_values=[7.0], n_samples=300, seed=5,
        calibrations={"pKa0": fit, "pKa1": fit})
    assert r["ok"] and r["mc_pka0"]["calibrated"] and r["mc_pka1"]["calibrated"]
    assert np.array_equal(r["mc_pka0"]["pka_samples"], r["mc_pka1"]["pka_samples"])
    # n=1 → sin error estimable: cae al MC crudo con aviso
    n1 = dict(fit, n=1, loo_rmse=float("nan"))
    r1 = m.run_chain_montecarlo(
        data_A={"a": 1}, data_neutral_macro={"n": 1}, data_B={"b": 1}, data_C=None,
        f_R=1.0, f_D=0.0, pH_values=[7.0], n_samples=300, seed=5,
        calibrations={"pKa0": n1})
    assert not r1["mc_pka0"].get("calibrated") and any(
        f"n<{m.CALIB_MIN_N}" in w for w in r1["warnings"])


def test_double_difference_class(monkeypatch):
    # prueba la mecánica del ajuste con n=2 (bajo CALIB_MIN_N): override explícito
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_below_min", True)
    def pts(mol_pts):
        return [{"mol": f"M{i}", "calc_goodvibes": x, "calc_harmonic": x, "exp": y}
                for i, (x, y) in enumerate(mol_pts)]
    monkeypatch.setattr(m, "EMPIRICAL_CALIBRATION", {
        "ammonium_dianion": pts([(17.0, 9.0), (18.0, 10.0)]),     # sesgo −8 → Δ=−8
        "amidinium":        pts([(10.0, 9.0), (12.0, 11.0)]),     # Δ=−1
        "ammonium":         pts([(12.0, 10.0), (13.0, 11.0)]),    # Δ=−2
    })
    assert m.EMPIRICAL_TYPE_META["double_difference"]["dq"] == "−1→−2"
    for t in ("amidinium", "ammonium", "ammonium_dianion", "pyridinium_anion"):
        assert t in m.EMPIRICAL_TYPE_META
    assert m._dq_of_type("amidinium") == (1, 0) and m._dq_of_type("ammonium_dianion") == (-1, -2)
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "force_class", False)
    res = m.apply_double_difference_correction(
        16.0, mode="offset", pka_result={"q_acid": -1, "q_base": -2})
    # pKa_corr = 16 + (−8) + (−1 − (−2)) = 9
    assert res["ok"] and res["pKa_corrected"] == pytest.approx(9.0)
    assert res["components"] == pytest.approx({"ammonium_dianion": -8.0, "amidinium": -1.0, "ammonium": -2.0})
    # LOO de dos puntos con sesgos iguales → 0 en los tres
    assert res["sigma"] == pytest.approx(0.0)
    # clase de carga incorrecta → rechazada (no pasa por --force-charge-class)
    bad = m.apply_double_difference_correction(
        16.0, pka_result={"q_acid": 0, "q_base": -1})
    assert not bad["ok"] and "−1→−2" in bad["reason"]
    # componente faltante → error explícito
    monkeypatch.setattr(m, "EMPIRICAL_CALIBRATION", {"amidinium": pts([(10.0, 9.0), (12.0, 11.0)])})
    miss = m.apply_double_difference_correction(16.0)
    assert not miss["ok"] and "ammonium_dianion" in miss["reason"]


# ── pka-spec-v2: mínimo de calibradores, inversión de imaginarias, peor caso ──

def test_calib_min_n_rejects_n2_unless_override(monkeypatch):
    assert m.CALIB_MIN_N == 3
    _calib_xy(monkeypatch, [(2.0, 1.0), (3.0, 2.2)])
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_n1", False)
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_below_min", False)
    fit = m.fit_empirical_correction("t")
    assert not fit["ok"] and "CALIB_MIN_N" in fit["reason"]
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_below_min", True)
    fit = m.fit_empirical_correction("t")
    assert fit["ok"] and any("sensibilidad" in w for w in fit["warnings"])
    # el MC calibrado no usa una clase bajo el mínimo aunque el ajuste pase
    assert m._calibration_intercept_draw(fit, 10, np.random.default_rng(0)) is None


def test_calib_min_n_accepts_n3(monkeypatch):
    _calib_xy(monkeypatch, [(2.0, 1.0), (2.2, 1.1), (2.4, 1.5)])
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_below_min", False)
    fit = m.fit_empirical_correction("t", mode="offset")
    assert fit["ok"] and fit["n"] == 3
    assert not any("sensibilidad" in w for w in fit["warnings"])


def test_gv_invert_args(monkeypatch):
    m._GV_INVERT_FLAG.clear()
    monkeypatch.setitem(m._GV_INVERT_FLAG, "py", "--invert")
    assert m._gv_invert_args("py") == ["--invert", "20"]
    monkeypatch.setattr(m, "GV_INVERT_IFREQ", None)
    assert m._gv_invert_args("py") == []


def test_worst_case_m3_test():
    # pKa1 bajo: B domina a pH 8.1 si pKa2 está alto
    wc = m.worst_case_m3_test(2.0, 3.5, 15.0, 0.0, 1.0, bias=4.0, pH=8.1)
    assert wc["robust"] and wc["f_B_shifted"] >= 0.95
    # b* = sesgo con el que f(B) cae exactamente al 95 %
    pk_crit = 15.0 - wc["bias_star"]
    f_at = m.boltzmann_5state(2.0, 3.5, pk_crit, 0.0, 1.0, [8.1]).iloc[0]["f_anion"]
    assert f_at == pytest.approx(0.95, abs=1e-6)
    # un sesgo mayor que b* rompe la robustez
    wc2 = m.worst_case_m3_test(2.0, 3.5, 15.0, 0.0, 1.0, bias=wc["bias_star"] + 0.5, pH=8.1)
    assert wc2["robust"] is False
    # sin calibrador −1→−2 se reporta solo b*
    wc3 = m.worst_case_m3_test(2.0, 3.5, 15.0, 0.0, 1.0, bias=None, pH=8.1)
    assert wc3["robust"] is None and wc3["bias_star"] is not None


def test_worst_case_bias_uses_minus1_minus2_classes(monkeypatch):
    monkeypatch.setattr(m, "EMPIRICAL_CALIBRATION", {
        "amino_dianion": [{"mol": "Glu", "calc_goodvibes": 15.7, "calc_harmonic": 15.6, "exp": 9.7}],
        "ammonium_dianion": [{"mol": "IDA", "calc_goodvibes": 14.0, "calc_harmonic": None, "exp": 9.3}],
        "amino": [{"mol": "Gly", "calc_goodvibes": 30.0, "calc_harmonic": 30.0, "exp": 9.6}],
    })
    b = m.worst_case_bias_m3("goodvibes")
    assert b["ok"] and b["mol"] == "Glu" and b["bias"] == pytest.approx(6.0) and b["n"] == 2


# ── --censo-refine: lectura de 3_REFINEMENT ────────────────────────────────

def _censo_out(path: Path, rows: list[tuple[str, float, float]]) -> None:
    """Tabla con el formato de N_PARTE.out de CENSO 3 (Gtot, ΔGtot, peso)."""
    lines = ["CONF#  E (DFT)  ΔGsolv  GmRRHO  Gtot  ΔGtot  Boltzmann weight"]
    for name, gtot, w in rows:
        lines.append(f"{name}  -1.0  0.000000  0.01  {gtot:.6f}  0.00  {w:.2f}")
    path.write_text("\n".join(lines) + "\n")


def _censo_xyz(path: Path, names: list[str], x: float) -> None:
    path.write_text("".join(f"1\n{n}\nH {x} 0.0 0.0\n" for n in names))


def test_find_censo_result_prefers_refinement_geometries(tmp_path, monkeypatch):
    crest = tmp_path / "crest_conformers.xyz"
    _censo_xyz(crest, ["CONF1", "CONF2", "CONF3"], 0.0)
    _censo_out(tmp_path / "1_SCREENING.out", [("CONF1", -2.0, 60), ("CONF2", -1.9, 40)])
    _censo_xyz(tmp_path / "1_SCREENING.xyz", ["CONF1", "CONF2"], 0.0)

    monkeypatch.setattr(m, "CENSO_REFINE", True)
    out, ids, geoms = m.find_censo_result(tmp_path, crest)
    assert out.name == "1_SCREENING.out" and geoms[0] == crest
    assert not m.censo_result_complete(out)

    _censo_out(tmp_path / "3_REFINEMENT.out", [("CONF2", -3.0, 80), ("CONF1", -2.9, 20)])
    _censo_xyz(tmp_path / "3_REFINEMENT.xyz", ["CONF2", "CONF1"], 9.0)
    out, ids, geoms = m.find_censo_result(tmp_path, crest)
    assert out.name == "3_REFINEMENT.out" and ids == ["CONF2", "CONF1"]
    # nunca las geometrías GFN2 de CREST para confs ya optimizados con DFT
    assert geoms == [tmp_path / "3_REFINEMENT.xyz"]
    assert m.censo_result_complete(out)


def test_censo_cli_parts(monkeypatch, tmp_path):
    monkeypatch.setattr(m, "CENSO_RC", None)
    monkeypatch.setattr(m, "CENSO_REFINE", False)
    assert m.censo_cli_parts() == ["--prescreening", "--screening"]
    monkeypatch.setattr(m, "CENSO_REFINE", True)
    parts = m.censo_cli_parts(include_prescreening=False)
    assert parts[:4] == ["--screening", "--optimization", "--refinement", "--ignore-failed"]
    assert parts[-2:] == ["--inprc", str(m.CENSO_REFINE_RC)]
    assert m.censo_timeout_hours() == m.CENSO_TIMEOUT_HOURS_REFINE


def test_stage_censo_reuses_existing_refinement(tmp_path, monkeypatch):
    mol_dir = tmp_path / "run"
    censo_dir = mol_dir / "01b_censo_prescreen"
    censo_dir.mkdir(parents=True)
    crest = tmp_path / "crest_conformers.xyz"
    _censo_xyz(crest, ["CONF1", "CONF2"], 0.0)
    _censo_out(censo_dir / "3_REFINEMENT.out", [("CONF2", -3.0, 100)])
    _censo_xyz(censo_dir / "3_REFINEMENT.xyz", ["CONF2"], 9.0)
    monkeypatch.setattr(m, "CENSO_REFINE", True)
    monkeypatch.setattr(m, "_run_censo_locked",
                        lambda *a, **k: pytest.fail("no debía relanzar CENSO"))
    out = m.stage_censo_prescreen(crest, mol_dir, 0, LOG)
    assert out.name == "censo_ensemble_sorted.xyz"
    text = out.read_text()
    assert "CONF2" in text and "9.0" in text and "CONF1" not in text
    assert m._find_screening_out(mol_dir) == censo_dir / "3_REFINEMENT.out"


def test_cli_flags_pka_calibrate(capsys, monkeypatch):
    # el parser de `pka` expone --mc-calibrated, --mc-seed, --mc-samples y modo aicc
    m.main(["pka", "--help"])        # main captura el SystemExit de argparse
    h = capsys.readouterr().out
    assert all(f in h for f in ("--mc-calibrated", "--mc-seed", "--mc-samples", "aicc"))
    # `calibrate --fit-mode aicc` llega a la tabla de ajustes
    m.main(["calibrate", "--help"])
    assert "--fit-mode" in capsys.readouterr().out
    _calib_xy(monkeypatch, [(2, 1), (5, 3), (8, 4.5), (11, 7)])
    monkeypatch.setitem(m.EMPIRICAL_OPTS, "allow_below_min", True)
    m._calib_print_fits(["t"], "aicc")
    assert "modo aicc" in capsys.readouterr().out


def test_calibrate_no_stale_builtin_points(monkeypatch, tmp_path):
    import copy, json
    monkeypatch.setattr(m, "EMPIRICAL_CALIBRATION", copy.deepcopy(m.EMPIRICAL_CALIBRATION))
    monkeypatch.setattr(m, "_calib_compute_mol", lambda mol, spec, base, T, meth, gsrc, use_smd=False:
                        {t["slot"]: {"ok": True, "pKa": 3.5, "method": "born"}
                         for t in spec["transitions"]})
    out = tmp_path / "c.json"
    m.main(["calibrate", "--base", str(tmp_path), "--mols", "VALINE", "--gcorr", "goodvibes",
            "--out", str(out)])
    cal = json.loads(out.read_text())["EMPIRICAL_CALIBRATION"]
    gly = next(p for p in cal["carboxyl"] if p["mol"] == "Glycine")
    val = next(p for p in cal["carboxyl"] if p["mol"] == "Valine")
    assert gly["calc_goodvibes"] is None and gly["calc_harmonic"] is None
    assert val["calc_goodvibes"] == 3.5 and val["calc_harmonic"] is None


def test_fit_harmonic_none_not_ok(monkeypatch):
    monkeypatch.setattr(m, "EMPIRICAL_CALIBRATION", {"t": [
        {"mol": "M", "calc_goodvibes": 3.0, "calc_harmonic": None, "exp": 2.0}]})
    assert not m.fit_empirical_correction("t", gcorr_source="harmonic")["ok"]


def test_mc_linear_slope_uncertainty(monkeypatch):
    # fit lineal sintético (n=6): la std del MC crece con |x0 − x̄| según la fórmula
    xs = [2.0, 3.0, 4.5, 6.0, 8.0, 10.0]
    ys = [1.0, 2.4, 3.0, 4.9, 6.1, 7.2]
    _calib_xy(monkeypatch, list(zip(xs, ys)))
    fit = m.fit_empirical_correction("t", mode="linear")
    fit["group_type"] = "t"
    n, xb = fit["n"], float(np.mean(xs))
    sxx = float(np.sum((np.array(xs) - xb) ** 2))
    N = 40000
    draw = m._calibration_intercept_draw(fit, N, np.random.default_rng(11))
    stds = {}
    for dG in (13.0, 40.0):            # pKa crudo ≈ 9.5 (cerca de x̄=5.6) y ≈ 29 (lejos)
        monkeypatch.setattr(m, "compute_pka", lambda *a, dG=dG, **k:
                            {"ok": True, "dG_aq_kcal": dG, "method": "direct", "warnings": []})
        r = m.compute_pka_montecarlo({}, {}, calibration_fit=fit, calib_intercept_draw=draw,
                                     rng=np.random.default_rng(1), n_samples=N)
        x0 = dG / (m.R_KCAL * m.T_STD * m.LN10)
        scale = fit["loo_rmse"] * math.sqrt(1 + 1 / n + (x0 - xb) ** 2 / sxx)
        mad = float(np.median(np.abs(r["pka_samples"] - r["pka_median"])))
        # mediana del |t_{n-2}| · escala (robusto a la cola pesada)
        t_med = float(np.median(np.abs(np.random.default_rng(5).standard_t(n - 2, size=400000))))
        assert mad == pytest.approx(t_med * scale, rel=0.08)
        stds[dG] = mad
    assert stds[40.0] > 2 * stds[13.0]


def test_cli_extra_micro_degeneracy_flags(capsys):
    m.main(["pka", "--help"])
    h = capsys.readouterr().out
    assert "--extra-micro" in h and "--degeneracy" in h
    members = ["GR", "GD", "GD2"]
    assert m._parse_degeneracy(["GD=2", "GR=0.5"], members) == {"GD": 2.0, "GR": 0.5}
    for bad in (["GD"], ["GD=x"], ["GD=0"], ["GD=-1"], ["ZZ=2"]):
        with pytest.raises(SystemExit):
            m._parse_degeneracy(bad, members)
