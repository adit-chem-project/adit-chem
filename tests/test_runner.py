import os
import subprocess
import sys
from pathlib import Path

import pytest

from adit.runner import RunTarget, block_reason, executable_of, run_and_wait, target_from_dir

REPO = Path(__file__).resolve().parents[1]


def _cfg(sk_root, *, enable_run=True):
    from tests.conftest import cfg_for

    cfg = cfg_for(sk_root)
    cfg.enable_run = enable_run
    return cfg


def _generated(tmp_path, sk_root, profile="local"):
    from tests.conftest import cfg_for, make_fake_skset
    from adit.project import write_project
    from adit.samples import copy_sample
    from adit.spec import CalculationSpec

    make_fake_skset(sk_root, "mio-1-1", ["H", "O"])
    spec = CalculationSpec.load(copy_sample("water_generated", tmp_path / f"spec_{profile}.json"))
    spec.runtime.profile = profile
    spec.runtime.ncpus = spec.runtime.mpiprocs * spec.runtime.omp_threads   # a cluster profile checks that the cores fit
    cfg = cfg_for(sk_root)
    out = tmp_path / f"run_{profile}"
    write_project(spec, cfg, out)
    return out


def test_the_executable_is_taken_from_the_command():
    assert executable_of("dftb+ > output.log 2>&1") == "dftb+"
    assert executable_of("cd x && mpirun -np 4 pw.x -in pw.in") == "pw.x"
    assert executable_of("srun -n 8 /opt/vasp/bin/vasp_std") == "/opt/vasp/bin/vasp_std"
    assert executable_of("") == ""


def test_running_is_refused_when_the_setting_is_off(tmp_path, sk_root):
    out = _generated(tmp_path, sk_root)
    target = target_from_dir(out, _cfg(sk_root))
    reason = block_reason(_cfg(sk_root, enable_run=False), target)
    assert "enable_run" in reason
    assert "投入はしません" in reason or "never submits" in reason


def test_cluster_output_is_never_run(tmp_path, sk_root):
    out = _generated(tmp_path, sk_root, profile="cluster")
    target = target_from_dir(out, _cfg(sk_root))
    assert target.kind != "direct"
    reason = block_reason(_cfg(sk_root), target)
    assert "クラスタ用" in reason or "Cluster" in reason


def test_a_missing_executable_says_where_to_look(tmp_path, sk_root):
    out = _generated(tmp_path, sk_root)
    target = target_from_dir(out, _cfg(sk_root))
    target.exe = "no_such_program_12345"
    reason = block_reason(_cfg(sk_root), target)
    assert (os.name == "nt") or "PATH" in reason


def test_a_directory_without_submit_sh_is_refused(tmp_path, sk_root):
    target = RunTarget(tmp_path, "direct", "bash")
    assert "submit.sh" in block_reason(_cfg(sk_root), target)


def test_nothing_generated_yet(sk_root):
    reason = block_reason(_cfg(sk_root), None)
    assert "生成" in reason or "generate" in reason


@pytest.mark.skipif(os.name == "nt", reason="Windows では submit.sh を実行しない (block_reason で止める)")
def test_it_runs_and_reports_the_exit_code(tmp_path):
    d = tmp_path / "fake"
    d.mkdir()
    (d / "submit.sh").write_text("echo hello\nexit 3\n", encoding="utf-8")
    lines = []
    code = run_and_wait(RunTarget(d, "direct", "echo"), on_line=lines.append)
    assert code == 3 and lines == ["hello"]
    assert (d / "run.log").read_text(encoding="utf-8").strip() == "hello"


@pytest.mark.skipif(os.name == "nt", reason="Windows では実行しない")
def test_the_cli_can_run_what_it_generated(tmp_path, sk_root):
    from adit.config import save_config

    from tests.conftest import make_fake_skset

    make_fake_skset(sk_root, "mio-1-1", ["H", "O"])
    config = tmp_path / "cluster.toml"
    save_config(_cfg(sk_root), config)
    spec = tmp_path / "spec.json"
    env = {**os.environ, "ADIT_CONFIG": str(config), "PATH": os.environ.get("PATH", "")}
    assert subprocess.run([sys.executable, "-m", "adit.cli", "--sample", "water_generated", str(spec)],
                          capture_output=True, text=True, cwd=REPO, timeout=300, env=env).returncode == 0
    r = subprocess.run([sys.executable, "-m", "adit.cli", str(spec), str(tmp_path / "out"), "--run",
                        "--config", str(config)], capture_output=True, text=True, cwd=REPO, timeout=600, env=env)
    if "が見つかりません" in r.stderr or "not found" in r.stderr:
        pytest.skip("dftb+ が PATH にない")
    assert r.returncode in (0, 4), r.stderr
    assert "実行します" in r.stdout or "running" in r.stdout
    assert (tmp_path / "out" / "run.log").is_file()


def test_the_cli_refuses_with_the_same_words_as_the_screens(tmp_path, sk_root):
    from adit.config import save_config

    from tests.conftest import make_fake_skset

    make_fake_skset(sk_root, "mio-1-1", ["H", "O"])
    config = tmp_path / "cluster.toml"
    save_config(_cfg(sk_root, enable_run=False), config)
    spec = tmp_path / "spec.json"
    env = {**os.environ, "ADIT_CONFIG": str(config)}
    subprocess.run([sys.executable, "-m", "adit.cli", "--sample", "water_generated", str(spec)],
                   capture_output=True, text=True, cwd=REPO, timeout=300, env=env)
    r = subprocess.run([sys.executable, "-m", "adit.cli", str(spec), str(tmp_path / "out"), "--run",
                        "--config", str(config)], capture_output=True, text=True, cwd=REPO, timeout=300, env=env)
    assert r.returncode == 3
    assert "enable_run" in r.stderr and ("投入はしません" in r.stderr or "never submits" in r.stderr)


def test_only_the_command_line_runs_the_generated_input():
    # The desktop window does not run calculations: people run them in the terminal, only the CLI keeps --run.
    root = REPO / "src" / "adit"
    assert "from adit.runner import" in (root / "cli.py").read_text(encoding="utf-8")
    for path in ("gui/main_window.py", "web/server.py"):
        text = (root / path).read_text(encoding="utf-8")
        assert "block_reason" not in text and "runner import start" not in text, path
        assert "run_and_wait" not in text, path
