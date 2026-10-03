
import re

import pytest

from adit import lang
from adit.project import build_project
from adit.spec import (AtomsData, CalculationSpec, DftbMethod, KPoints, MDSettings, OrcaMethod,
                        Runtime, Structure, Task, XtbMethod)
from adit.structure import from_bulk
from tests.conftest import cfg_for, make_fake_skset, water_spec
from tests.test_espresso import make_fake_upf, si_spec as qe_si_spec
from tests.test_vasp import h2o_spec as vasp_h2o_spec, si_spec as vasp_si_spec

_JA = re.compile(r"[぀-ヿ一-鿿]")

MARKERS = {
    "dftbplus": {"geometry_optimization": "geom.out.gen", "molecular_dynamics": "geo_end.xyz",
                 "vibrations": "hessian.out", "band_structure": "bands/band.out"},
    "xtb": {"geometry_optimization": "xtbopt.xyz", "molecular_dynamics": "xtb.trj", "vibrations": "vibspectrum"},
    "vasp": {"geometry_optimization": "CONTCAR       最適化後", "molecular_dynamics": "XDATCAR",
             "vibrations": "dynamical matrix", "band_structure": "bands/EIGENVAL"},
    "espresso": {"geometry_optimization": "Begin final coordinates", "molecular_dynamics": "ATOMIC_POSITIONS",
                 "vibrations": "dynmat.out", "band_structure": "bands/output.log"},
    "orca": {"geometry_optimization": "orca.xyz", "molecular_dynamics": "trajectory.xyz", "vibrations": "VIBRATIONAL FREQUENCIES"},
}


@pytest.fixture
def cfg_all(sk_root, tmp_path):
    make_fake_skset(sk_root, "fake-si", ["Si"])
    cfg = cfg_for(sk_root)
    cfg.pseudo_root = str(make_fake_upf(tmp_path / "pseudo"))
    return cfg


def _spec(code: str, task: str) -> CalculationSpec:
    t = Task(type=task, md=MDSettings(ensemble="NVE") if code == "vasp" else MDSettings())
    local = Runtime(profile="local")
    if code == "dftbplus":
        if task == "band_structure":
            a = from_bulk("Si")
            return CalculationSpec(structure=Structure(source="bulk", source_ref="Si", atoms=AtomsData.from_ase(a)),
                                   method=DftbMethod(sk_set="fake-si"), kpoints=KPoints(mode="mesh", mesh=(4, 4, 4)),
                                   task=t, runtime=local)
        return water_spec(task=t)
    if code == "xtb":
        return water_spec(method=XtbMethod(), task=t)
    if code == "orca":
        return water_spec(method=OrcaMethod(), task=t)
    if code == "vasp":
        return (vasp_si_spec if task == "band_structure" else vasp_h2o_spec)(task=t, runtime=local)
    return qe_si_spec(task=t)


def _results_section(readme: str) -> str:
    return readme.split("== 結果の見方 ==", 1)[1].split("== 引用 ==", 1)[0]


CASES = [(code, task) for code, m in MARKERS.items() for task in [*m, "single_point"]]


@pytest.mark.parametrize("code,task", CASES)
def test_outputs_only_for_the_chosen_task(cfg_all, code, task):
    r = build_project(_spec(code, task), cfg_all).texts["README.txt"]
    sec = _results_section(r)
    for other, marker in MARKERS[code].items():
        if other == task:
            assert marker in sec, (code, task, marker)
        else:
            assert marker not in sec, (code, task, marker)
    assert "output.log" in sec


def test_local_readme_has_both_sections(cfg_all, tmp_path):
    out = tmp_path / "water opt"
    r = build_project(water_spec(), cfg_all, output_dir=out).texts["README.txt"]
    here = r.split("== この PC で実行する ==", 1)[1].split("== 研究室のクラスタで実行したいとき ==", 1)
    assert len(here) == 2
    local, cluster = here
    assert "cd '<この計算ディレクトリのパス>'" in local
    assert str(out.resolve()) not in r
    assert "command -v dftb+" in local and "PATH (コマンドを探す場所の一覧)" in local and "bash submit.sh" in local
    for s in ("docs/SETTINGS.md (環境設定) の「クラスタで実行する場合」", "[profiles.remote]", 'kind = "pbs"', "<クラスタのホスト名>", "scp -r", "rsync -av", "ssh ", "qsub submit.sh", "sbatch submit.sh",
              "qstat -u $USER", "squeue -u $USER", "cd <クラスタでの作業ディレクトリ>/'water opt'"):
        assert s in cluster, s
    assert "scp -r '<この計算ディレクトリのパス>'" in cluster
    assert "ジョブスクリプト (" in r


def test_english_readme_uses_current_directory_placeholder(cfg_all, tmp_path, monkeypatch):
    from adit import lang

    monkeypatch.setattr(lang, "LANGUAGE", "en")
    out = tmp_path / "moved calculation"
    readme = build_project(water_spec(), cfg_all, output_dir=out).texts["README.txt"]
    assert "cd '<path to this calculation directory>'" in readme
    assert "scp -r '<path to this calculation directory>'" in readme
    assert "If you copied it to another PC" in readme
    assert str(out.resolve()) not in readme


@pytest.mark.parametrize("language", ["ja", "en"])
@pytest.mark.parametrize("selection", ["explicit", "environment", "current", "legacy", "both", "first_run"])
def test_cli_readme_identifies_the_settings_file_used(tmp_path, monkeypatch, language, selection):
    from adit.cli import main
    from adit.config import default_config, save_config

    monkeypatch.chdir(tmp_path)
    monkeypatch.delenv("ADIT_CONFIG", raising=False)
    monkeypatch.delenv("QCGUI_CONFIG", raising=False)
    monkeypatch.setenv("ADIT_LANG", language)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "settings"))
    monkeypatch.setenv("APPDATA", str(tmp_path / "settings"))
    monkeypatch.setattr(lang, "LANGUAGE", language)
    current = tmp_path / "settings" / "adit" / "cluster.toml"
    legacy = tmp_path / "settings" / "qcgui" / "cluster.toml"
    custom = tmp_path / "custom" / "cluster.toml"
    expected = current
    args = []
    if selection in {"explicit", "environment", "legacy", "both"}:
        save_config(default_config(), legacy)
    if selection in {"current", "both"}:
        save_config(default_config(), current)
    if selection in {"explicit", "environment"}:
        expected = custom
        save_config(default_config(), custom)
        if selection == "explicit":
            args = ["--config", "custom/cluster.toml"]  # Relative CLI paths must remain unambiguous in the README.
        else:
            monkeypatch.setenv("ADIT_CONFIG", str(custom))
    elif selection == "legacy":
        expected = legacy
    spec = tmp_path / "source.json"
    water_spec(method=XtbMethod()).save(spec)
    output = tmp_path / "output"
    assert main([str(spec), str(output), *args]) == 0
    readme = (output / "README.txt").read_text(encoding="utf-8")
    assert str(expected) in readme
    for unused in {current, legacy, custom} - {expected}:
        assert str(unused) not in readme


def test_readme_keeps_loaded_path_when_config_environment_changes(tmp_path, monkeypatch):
    from adit.config import default_config, load_config, save_config

    original = tmp_path / "original.toml"
    save_config(default_config(), original)
    cfg = load_config(original)
    monkeypatch.setenv("ADIT_CONFIG", str(tmp_path / "other.toml"))
    readme = build_project(water_spec(method=XtbMethod()), cfg).texts["README.txt"]
    assert str(original) in readme
    assert str(tmp_path / "other.toml") not in readme


def test_vibrations_readme_identifies_reference_structure(cfg_all):
    r = build_project(_spec("dftbplus", "vibrations"), cfg_all).texts["README.txt"]
    assert "== 振動解析に使う構造 ==" in r
    assert "入力した構造をそのまま基準" in r
    assert "adit-gen --continue-from <最適化のディレクトリ>" in r


def test_english_vibrations_readme_identifies_reference_structure(cfg_all, english):
    r = build_project(_spec("dftbplus", "vibrations"), cfg_all).texts["README.txt"]
    assert "== Structure used for the vibrational analysis ==" in r
    assert "uses the input structure as-is" in r
    assert "adit-gen --continue-from <optimization directory>" in r


def test_xtb_md_readme_records_mass_constraints_and_scc(cfg_all):
    r = build_project(_spec("xtb", "molecular_dynamics"), cfg_all).texts["README.txt"]
    assert "hmass = 4 u" in r and "shake = 2" in r and "sccacc = 2" in r


@pytest.mark.parametrize("profile,submit,other", [("cluster", "qsub submit.sh", "sbatch"), ("slurm", "sbatch submit.sh", "qsub")])
def test_cluster_profile_readme_is_not_duplicated(cfg_all, profile, submit, other):
    spec = water_spec(runtime=Runtime(profile=profile, ncpus=8, omp_threads=8, walltime="01:00:00", job_name="w"))
    r = build_project(spec, cfg_all).texts["README.txt"]
    assert r.count(submit) == 1 and other not in r
    assert "bash submit.sh" not in r and "[profiles.remote]" not in r
    assert "local" in r.split("== この PC で実行する ==", 1)[1].split("==", 1)[0]
    assert "dftbplus/25.1" not in r and "#PBS -q" not in r and "--partition" not in r


@pytest.fixture
def english():
    lang.set_language("en")
    yield
    lang.set_language("ja")


@pytest.mark.parametrize("code", list(MARKERS))
def test_english_readme_for_every_code(cfg_all, english, code):
    r = build_project(_spec(code, "geometry_optimization"), cfg_all).texts["README.txt"]
    assert not _JA.search(r), [l for l in r.splitlines() if _JA.search(l)]
    assert "== Run on this PC ==" in r and "== Running on your group's cluster ==" in r and "PATH (the list of places" in r
