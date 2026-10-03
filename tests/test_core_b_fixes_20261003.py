"""Fixes from the 2026-10-03 core/codes audit that have no natural home in another test file."""

from __future__ import annotations

import pytest

from adit import lang
from adit.config import Profile, config_path, default_config, save_config
from adit.spec import Cp2kMethod, EspressoMethod, LammpsMethod, Runtime, Task, VaspMethod, XtbMethod
from tests.conftest import cfg_for, water_spec


@pytest.fixture(autouse=True)
def _restore_language():
    before = lang.LANGUAGE
    yield
    lang.set_language(before)


def _japanese(text: str) -> bool:
    return any("\u3040" <= ch <= "\u30ff" or "\u4e00" <= ch <= "\u9fff" for ch in text)


def test_convert_cli_uses_the_configured_language(tmp_path, capsys, monkeypatch):
    from adit.convert import main

    monkeypatch.delenv("ADIT_LANG", raising=False)
    missing, out = str(tmp_path / "missing.xyz"), str(tmp_path / "out.xyz")
    save_config(default_config().model_copy(update={"language": "en"}), config_path())
    lang.set_language("ja")
    assert main(["structure", missing, out]) == 1
    err = capsys.readouterr().err
    assert err.strip() and not _japanese(err)
    monkeypatch.setenv("ADIT_LANG", "ja")
    assert main(["structure", missing, out]) == 1
    assert _japanese(capsys.readouterr().err)
    monkeypatch.delenv("ADIT_LANG")
    config_path().write_text("language = [broken", encoding="utf-8")
    lang.set_language("ja")
    assert main(["structure", missing, out]) == 1
    assert _japanese(capsys.readouterr().err)


def test_documented_settings_are_english_in_english():
    from adit.docvalues import documented_settings, format_settings

    settings = documented_settings(code="vasp", task="vibrations")
    assert settings and all(s.source_en for s in settings)
    lang.set_language("en")
    text = format_settings(settings)
    assert "retrieved 2026-09-12" in text and "VASP tutorial" in text and not _japanese(text)
    lang.set_language("ja")
    assert "2026-09-12 取得" in format_settings(settings)


@pytest.mark.parametrize("method,expected", [
    (VaspMethod(), "mpirun -np 6 vasp_std"), (EspressoMethod(), "mpirun -np 6 pw.x"),
    (Cp2kMethod(), "mpirun -np 6 cp2k.psmp"), (LammpsMethod(), "mpirun -np 6 lmp")])
def test_default_mpi_commands_count_processes_over_all_nodes(method, expected):
    import adit.codes.cp2k  # noqa: F401
    import adit.codes.espresso  # noqa: F401
    import adit.codes.lammps  # noqa: F401
    import adit.codes.vasp  # noqa: F401
    from adit.codes.base import GENERATORS

    spec = water_spec(method=method, runtime=Runtime(profile="local", nodes=2, ncpus=3, mpiprocs=3, omp_threads=1))
    generator = GENERATORS[method.code]
    assert expected in generator.run_command(spec, Profile(kind="direct"))
    one_node = spec.model_copy(update={"runtime": Runtime(profile="local", nodes=1, ncpus=3, mpiprocs=3, omp_threads=1)})
    assert expected.replace("-np 6", "-np 3") in generator.run_command(one_node, Profile(kind="direct"))
    custom = Profile(kind="direct", commands={method.code: "mpirun -np {mpiprocs} -x OMP_NUM_THREADS={omp_threads} x_{binary}"})
    assert "mpirun -np 3 -x OMP_NUM_THREADS=1 x_" in generator.run_command(spec, custom)


def test_sella_run_command_placeholders(tmp_path, sk_root):
    from adit.ts_setup import TsError, write_sella

    cfg = cfg_for(sk_root)
    spec = water_spec(method=XtbMethod(gfn="2"), task=Task(type="geometry_optimization", max_steps=50),
                      runtime=Runtime(profile="local", nodes=2, ncpus=3, mpiprocs=3, omp_threads=1))
    cfg.profiles["local"].commands["python"] = "env OMP_NUM_THREADS=${omp} python3"
    with pytest.raises(TsError, match=r"\{mpiprocs\} \{ntasks\} \{omp_threads\} \{binary\}"):
        write_sella(spec, cfg, tmp_path / "bad")
    assert not (tmp_path / "bad").exists()
    cfg.profiles["local"].commands["python"] = "mpirun -np {ntasks} python3"
    write_sella(spec, cfg, tmp_path / "ts")
    assert "mpirun -np 6 python3 run_sella.py" in (tmp_path / "ts" / "submit.sh").read_text(encoding="utf-8")
