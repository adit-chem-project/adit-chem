
from __future__ import annotations

import io
import sys
from pathlib import Path

import pytest

from adit import compat, config



def test_readme_path_windows_drive_becomes_wsl_form():
    assert compat.readme_shell_path(r"C:\Users\taro\adit_runs\run_1", windows=True) == "/mnt/c/Users/taro/adit_runs/run_1"


def test_readme_path_windows_spaces_and_japanese_are_quoted():
    got = compat.readme_shell_path(r"D:\研究 データ\水 1", windows=True)
    assert got == "'/mnt/d/研究 データ/水 1'"


def test_readme_path_windows_unc_uses_slashes():
    got = compat.readme_shell_path(r"\\server\share\run", windows=True)
    assert "\\" not in got and got.endswith("/run")


def test_readme_path_posix_unchanged():
    assert compat.readme_shell_path("/home/a b/run", windows=False) == "'/home/a b/run'"



@pytest.mark.parametrize("sent, name", [
    (r"C:\Users\taro\Desktop\water.xyz", "water.xyz"),
    ("water.xyz", "water.xyz"),
    ("../../etc/passwd", "passwd"),
    ("..", "upload"),
    ("水 分子.xyz", "水 分子.xyz"),
])
def test_upload_basename(sent, name):
    assert compat.upload_basename(sent) == name



def test_printable_stdio_on_cp932(monkeypatch):
    buf = io.BytesIO()
    fake = io.TextIOWrapper(buf, encoding="cp932")
    monkeypatch.setattr(sys, "stdout", fake)
    compat.ensure_printable_stdio()
    print("振動数 1600 cm⁻¹、距離 1.0 Å")
    fake.flush()
    out = buf.getvalue().decode("cp932")
    assert "振動数" in out and "\\u207b" in out


def test_printable_stdio_leaves_utf8(monkeypatch):
    fake = io.TextIOWrapper(io.BytesIO(), encoding="utf-8")
    monkeypatch.setattr(sys, "stdout", fake)
    compat.ensure_printable_stdio()
    assert fake.errors == "strict"



def test_config_with_bom_loads(tmp_path):
    p = tmp_path / "cluster.toml"
    p.write_bytes(b"\xef\xbb\xbf" + 'sk_root = "C:/Users/taro/slakos"\n[profiles.local]\nkind = "direct"\n'.encode())
    cfg = config.load_config(p)
    assert cfg.sk_root == "C:/Users/taro/slakos"
    config.set_top_level_value(p, "theme", "dark")
    assert config.load_config(p).theme == "dark"


def test_config_in_cp932_is_config_error(tmp_path):
    p = tmp_path / "cluster.toml"
    p.write_bytes('# 研究室の設定\nsk_root = "/home/太郎/slakos"\n'.encode("cp932"))
    with pytest.raises(config.ConfigError, match="UTF-8"):
        config.load_config(p)


def test_config_path_per_os(monkeypatch, tmp_path):
    monkeypatch.delenv("ADIT_CONFIG", raising=False)
    monkeypatch.delenv("QCGUI_CONFIG", raising=False)
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "xdg"))
    monkeypatch.setenv("APPDATA", str(tmp_path / "appdata"))
    import types
    monkeypatch.setattr(config, "os", types.SimpleNamespace(name="posix", environ=config.os.environ))
    assert config.config_path() == tmp_path / "xdg" / "adit" / "cluster.toml"
    monkeypatch.setattr(config, "os", types.SimpleNamespace(name="nt", environ=config.os.environ))
    assert config.config_path() == tmp_path / "appdata" / "adit" / "cluster.toml"



def test_analysis_writes_are_utf8_even_if_locale_is_cp932(monkeypatch, tmp_path):
    from adit.analysis import report

    src = Path(__file__).resolve().parents[1] / "examples" / "water_generated"
    if not (src / "detailed.out").is_file():
        pytest.skip("基準の実行結果が無い")
    import shutil
    run = tmp_path / "run"
    shutil.copytree(src, run)
    orig_write = Path.write_text

    def cp932_write(self, data, encoding=None, errors=None, newline=None):
        return orig_write(self, data, encoding=encoding or "cp932", errors=errors, newline=newline)

    monkeypatch.setattr(Path, "write_text", cp932_write)
    monkeypatch.setattr(report.AnalysisResult, "summary_text", lambda self: "距離 1.0 Å、1600 cm⁻¹")
    report.run_analysis(run, report.AnalysisOptions())
    assert (run / "analysis" / "summary.txt").read_text(encoding="utf-8").startswith("距離 1.0 Å")


def test_report_source_has_no_default_encoding_writes():
    import ast

    repo = Path(__file__).resolve().parents[1]
    bad = []
    for p in sorted((repo / "src" / "adit").rglob("*.py")) + sorted((repo / "tests").glob("*.py")):
        for node in ast.walk(ast.parse(p.read_text(encoding="utf-8"))):
            if not isinstance(node, ast.Call):
                continue
            f = node.func
            name = f.attr if isinstance(f, ast.Attribute) else f.id if isinstance(f, ast.Name) else ""
            if isinstance(f, ast.Attribute) and isinstance(f.value, ast.Name) and f.value.id == "webbrowser":
                continue
            if name not in ("read_text", "write_text", "open"):
                continue
            if any(k.arg == "encoding" for k in node.keywords):
                continue
            mode = None
            if name == "open" and node.args:
                mode = node.args[1] if (isinstance(f, ast.Name) and len(node.args) > 1) else (
                    node.args[0] if isinstance(f, ast.Attribute) else None)
            if isinstance(mode, ast.Constant) and "b" in str(mode.value):
                continue
            src = ast.get_source_segment(p.read_text(encoding="utf-8"), node) or ""
            if "/proc/version" in src:
                continue
            bad.append(f"{p.relative_to(repo)}:{node.lineno}: {src[:100]}")
    assert bad == []



def test_custom_frame_only_on_linux(monkeypatch):
    pytest.importorskip("PySide6")
    from adit.gui import titlebar

    monkeypatch.delenv("ADIT_FRAME", raising=False)
    monkeypatch.delenv("QCGUI_FRAME", raising=False)
    for plat, want in (("linux", True), ("darwin", False), ("win32", False)):
        monkeypatch.setattr(titlebar.sys, "platform", plat)
        assert titlebar.use_custom_frame("auto") is want, plat


def test_xcb_switch_only_on_linux(monkeypatch):
    pytest.importorskip("PySide6")
    from adit.gui import app

    monkeypatch.delenv("QT_QPA_PLATFORM", raising=False)
    monkeypatch.setenv("DISPLAY", ":0")
    for plat in ("darwin", "win32"):
        monkeypatch.setattr(app.sys, "platform", plat)
        app._prefer_xcb_on_wsl()
        assert "QT_QPA_PLATFORM" not in __import__("os").environ, plat


