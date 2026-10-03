import importlib.util
import os
import stat
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parent.parent
SPEC = REPO / "packaging" / "adit.spec"
LAUNCH = REPO / "packaging" / "launch.py"


def test_spec_bundles_the_scripts_copied_at_run_time_and_the_figure_backends():
    text = SPEC.read_text(encoding="utf-8")
    assert '"handoff.py"' in text and '"msd_worker.py"' in text
    for backend in ("backend_svg", "backend_pdf", "backend_ps"):
        assert f"matplotlib.backends.{backend}" in text
    assert 'collect_data_files("winpty"' in text and 'collect_dynamic_libs("winpty")' in text and 'platform == "win32"' in text
    assert "Path(SPECPATH).parent" in text and "os.getcwd()" not in text
    assert (REPO / "src" / "adit" / "handoff.py").is_file() and (REPO / "src" / "adit" / "analysis" / "msd_worker.py").is_file()


def test_stages_read_handoff_through_importlib_resources(tmp_path):
    from adit import handoff as hf
    from adit.stages import _handoff_source

    assert _handoff_source() == Path(hf.__file__).read_text(encoding="utf-8")


def _launch():
    spec = importlib.util.spec_from_file_location("adit_launch_under_test", LAUNCH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.skipif(os.name == "nt", reason="uses a shell script as the fake console executable")
def test_windowed_launcher_hands_commands_to_the_console_executable(tmp_path, monkeypatch):
    mod = _launch()
    exe = tmp_path / "ADIT"
    exe.write_text("", encoding="utf-8")
    cli = tmp_path / "adit-cli"
    cli.write_text('#!/bin/sh\n[ "$1" = gen ] && [ "$2" = --version ] && exit 7\nexit 9\n', encoding="utf-8")
    cli.chmod(cli.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setattr(sys, "executable", str(exe))
    monkeypatch.setattr(sys, "argv", ["ADIT", "gen", "--version"])
    assert mod._console_executable() == str(cli)
    monkeypatch.setattr(sys, "stdout", None)
    assert mod.main() == 7
    cli.unlink()
    assert mod._console_executable() is None
