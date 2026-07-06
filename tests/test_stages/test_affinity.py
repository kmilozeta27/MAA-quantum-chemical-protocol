"""tests/test_system/test_affinity.py — ResourceProfile (Issue B).

Mockea os.cpu_count() y la RAM (psutil o /proc/meminfo) para simular el
T14 Gen5 (12 cores lógicas, ~32 GB RAM) sin depender del hardware real
de CI. No usa fixtures de molécula (glicina/MGLYVAL): ResourceProfile no
toca química, solo aritmética de asignación de recursos.
"""

from __future__ import annotations

import builtins

import pytest

from maaqc.system.affinity import ResourceProfile


# ─────────────────────────────────────────────
#  _compute() — construcción directa (sin detect())
# ─────────────────────────────────────────────

class TestCompute:
    def test_t14_defaults_batch_size_1(self):
        """Perfil de referencia del T14 Gen5: 12 cores lógicas, 32 GB RAM,
        batch_size=1 (modo normal, sin --parallel)."""
        profile = ResourceProfile(total_cores=12, ram_mb=32_000, batch_size=1)

        assert profile.nprocs_crest == 12  # min(12, MAX_CREST_NPROCS=12)
        assert profile.nprocs_opt == 4     # 12 // OPT_BATCH_SIZE(3) = 4, cap MAX_OPT=4
        assert profile.nprocs_sp == 4
        assert profile.nprocs_tddft == 6   # 12 // TDDFT_BATCH_SIZE(2) = 6, cap 6
        assert profile.nprocs_numfreq == 12
        assert profile.nprocs_adc2 == 4    # cap NPROCS_ADC2_MAX, no 12

    def test_maxcore_never_below_500(self):
        """max(500, ...) es un piso duro en las 5 ramas de _compute();
        con muy poca RAM ninguna debería caer por debajo."""
        profile = ResourceProfile(total_cores=4, ram_mb=100, batch_size=1)

        assert profile.maxcore_opt >= 500
        assert profile.maxcore_sp >= 500
        assert profile.maxcore_numfreq >= 500
        assert profile.maxcore_tddft >= 500

    def test_adc2_maxcore_floor_independent_of_ram(self):
        """MAXCORE_ADC2_MIN (5000) es un piso propio de ADC2, distinto
        del piso genérico de 500 — con RAM abundante y pocos cores para
        ADC2, ram/cores puede superarlo; con RAM escasa, el piso manda."""
        profile = ResourceProfile(total_cores=4, ram_mb=1_000, batch_size=1)
        assert profile.maxcore_adc2 == profile.MAXCORE_ADC2_MIN

    def test_tddft_ram_penalty_applied(self):
        """maxcore_tddft debe ser menor que el cálculo "crudo" sin
        RAM_PENALTY_TDDFT, en la misma proporción (1.4x)."""
        profile = ResourceProfile(total_cores=12, ram_mb=32_000, batch_size=1)
        ram = int(32_000 * profile.SAFETY_FACTOR)
        active_tddft_cores = profile.nprocs_tddft * profile.TDDFT_BATCH_SIZE
        raw = max(500, ram // active_tddft_cores)
        expected = max(500, int(raw / profile.RAM_PENALTY_TDDFT))
        assert profile.maxcore_tddft == expected

    def test_batch_size_greater_than_1_overrides_opt_sp_split(self):
        """Con batch_size > 1 (--parallel), nprocs_opt/sp se derivan de
        c // p en vez de c // OPT_BATCH_SIZE fijo."""
        profile = ResourceProfile(total_cores=12, ram_mb=32_000, batch_size=4)
        assert profile.nprocs_opt == min(max(1, 12 // 4), profile.MAX_OPT_NPROCS_PER_JOB)
        assert profile.nprocs_sp == min(max(1, 12 // 4), profile.MAX_SP_NPROCS)

    def test_low_core_count_does_not_zero_out(self):
        """max(1, ...) evita nprocs=0 en sistemas con pocas cores lógicas
        (ej. CI con 2 cores) incluso con batch_size alto."""
        profile = ResourceProfile(total_cores=2, ram_mb=8_000, batch_size=8)
        assert profile.nprocs_opt >= 1
        assert profile.nprocs_sp >= 1


# ─────────────────────────────────────────────
#  detect() — mock de os.cpu_count() y psutil/meminfo
# ─────────────────────────────────────────────

class TestDetect:
    def test_detect_uses_cpu_count_and_psutil(self, monkeypatch):
        monkeypatch.setattr("os.cpu_count", lambda: 12)

        class _FakeMem:
            total = 32 * 1024 * 1024 * 1024  # 32 GB en bytes

        class _FakePsutil:
            @staticmethod
            def virtual_memory():
                return _FakeMem()

        monkeypatch.setitem(__import__("sys").modules, "psutil", _FakePsutil())

        profile = ResourceProfile.detect()
        assert profile.total_cores == 12
        assert profile.ram_mb == 32 * 1024

    def test_detect_falls_back_to_proc_meminfo_without_psutil(self, monkeypatch, tmp_path):
        monkeypatch.setattr("os.cpu_count", lambda: 8)
        monkeypatch.setitem(__import__("sys").modules, "psutil", None)

        real_import = builtins.__import__

        def _fake_import(name, *args, **kwargs):
            if name == "psutil":
                raise ImportError("no psutil en este entorno de CI")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", _fake_import)

        meminfo = tmp_path / "meminfo"
        meminfo.write_text("MemTotal:       16000000 kB\nMemFree:  1000 kB\n")

        real_open = builtins.open

        def _fake_open(path, *args, **kwargs):
            if path == "/proc/meminfo":
                return real_open(meminfo, *args, **kwargs)
            return real_open(path, *args, **kwargs)

        monkeypatch.setattr(builtins, "open", _fake_open)

        profile = ResourceProfile.detect()
        assert profile.total_cores == 8
        assert profile.ram_mb == 16000000 // 1024

    def test_detect_cpu_count_none_defaults_to_4(self, monkeypatch):
        monkeypatch.setattr("os.cpu_count", lambda: None)

        class _FakeMem:
            total = 8 * 1024 * 1024 * 1024

        class _FakePsutil:
            @staticmethod
            def virtual_memory():
                return _FakeMem()

        monkeypatch.setitem(__import__("sys").modules, "psutil", _FakePsutil())

        profile = ResourceProfile.detect()
        assert profile.total_cores == 4


# ─────────────────────────────────────────────
#  detect_with_batch() — total_cores/batch_size explícitos (--cores/--parallel)
# ─────────────────────────────────────────────

class TestDetectWithBatch:
    def test_explicit_ram_mb_skips_psutil_and_meminfo(self):
        """Cuando ram_mb > 0 se pasa explícito, no debería tocar
        psutil ni /proc/meminfo en absoluto."""
        profile = ResourceProfile.detect_with_batch(total_cores=6, batch_size=2, ram_mb=12_000)
        assert profile.total_cores == 6
        assert profile.ram_mb == 12_000
        assert profile.batch_size == 2

    def test_zero_ram_mb_triggers_autodetect(self, monkeypatch):
        class _FakeMem:
            total = 24 * 1024 * 1024 * 1024

        class _FakePsutil:
            @staticmethod
            def virtual_memory():
                return _FakeMem()

        monkeypatch.setitem(__import__("sys").modules, "psutil", _FakePsutil())

        profile = ResourceProfile.detect_with_batch(total_cores=6, batch_size=1, ram_mb=0)
        assert profile.ram_mb == 24 * 1024


# ─────────────────────────────────────────────
#  summary() — solo verifica que no reviente y contenga los campos clave
# ─────────────────────────────────────────────

def test_summary_contains_all_stage_sections():
    profile = ResourceProfile(total_cores=12, ram_mb=32_000, batch_size=1)
    text = profile.summary()
    for label in ("CREST/CENSO", "OPT", "SP", "TDDFT", "NumFreq", "ADC2"):
        assert label in text
