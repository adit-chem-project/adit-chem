import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
EX = REPO / "examples"
MD = EX / "dftb_md_water_generated"
VIB = EX / "xtb_vib_water_generated"


def _run(args, **kw):
    return subprocess.run([sys.executable, "-m", "adit.analysis.cli", *args],
                          capture_output=True, text=True, cwd=REPO, timeout=600, **kw)


def test_the_fit_range_the_summary_suggests_actually_works(tmp_path):
    for form in (["--msd-fit", "10,45"], ["--msd-fit", "10", "45"]):
        r = _run([str(MD), "-o", str(tmp_path / "o"), "--msd", *form])
        assert r.returncode == 0, r.stderr
        assert "拡散係数" in r.stdout or "diffusion" in r.stdout


def test_a_bad_fit_range_is_explained_in_japanese(tmp_path):
    r = _run([str(MD), "-o", str(tmp_path / "o"), "--msd", "--msd-fit", "45,10"])
    assert r.returncode != 0
    assert "小さく" in r.stderr or "smaller" in r.stderr


def test_the_font_warning_is_not_printed(tmp_path):
    r = _run([str(MD), "-o", str(tmp_path / "o")])
    assert "findfont" not in r.stdout + r.stderr


def test_the_diffusion_coefficient_comes_with_the_size_of_the_run(tmp_path):
    out = tmp_path / "o"
    r = _run([str(MD), "-o", str(out), "--msd"])
    assert r.returncode == 0, r.stderr
    assert "規模" in r.stdout
    summary = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    scale = summary["tables"]["msd_scale"]
    assert scale["n_atoms"] == 3 and scale["n_frames"] == 21 and scale["periodic"] is False


def test_frequencies_are_written_as_numbers(tmp_path):
    out = tmp_path / "o"
    assert _run([str(VIB), "-o", str(out)]).returncode == 0
    rows = (out / "frequencies.csv").read_text(encoding="utf-8").splitlines()
    assert rows[0] == "frequency_cm-1,ir_intensity,raman_activity"
    values = [r.split(",") for r in rows[1:]]
    assert any(float(v[0]) > 1600 and float(v[1]) > 0 for v in values)


def test_the_scaling_factor_is_applied_and_recorded(tmp_path):
    out = tmp_path / "o"
    r = _run([str(VIB), "-o", str(out), "--freq-scale", "0.5"])
    assert r.returncode == 0, r.stderr
    assert "補正係数 0.5" in r.stdout
    freqs = json.loads((out / "summary.json").read_text(encoding="utf-8"))["tables"]["frequencies"]
    assert max(freqs) == pytest.approx(1656.71 * 0.5, rel=1e-6)


def test_a_measured_spectrum_can_be_overlaid(tmp_path):
    measured = tmp_path / "m.txt"
    measured.write_text("# cm-1 intensity\n800 0.2\n1600 1.0\n", encoding="utf-8")
    r = _run([str(VIB), "-o", str(tmp_path / "o"), "--spectrum-measured", str(measured)])
    assert r.returncode == 0, r.stderr
    assert "測定したスペクトル" in r.stdout
    bad = _run([str(VIB), "-o", str(tmp_path / "o2"), "--spectrum-measured", str(tmp_path / "none.txt")])
    assert "測定データのファイルがありません" in bad.stdout


def test_a_directory_with_only_an_output_can_be_analyzed(tmp_path):
    import gzip

    d = tmp_path / "only_out"
    d.mkdir()
    (d / "output.log").write_bytes(gzip.decompress((REPO / "tests" / "data" / "gaussian_dvb_raman.out.gz").read_bytes()))
    assert _run([str(d), "-o", str(tmp_path / "o")]).returncode != 0
    r = _run([str(d), "-o", str(tmp_path / "o"), "--code", "gaussian"])
    assert r.returncode == 0, r.stderr
    assert "ラマン" in r.stdout
    bad = _run([str(d), "-o", str(tmp_path / "o"), "--code", "gauss"])
    assert "知らないコードです" in bad.stdout + bad.stderr


def test_help_is_short_and_help_all_is_complete():
    short = _run(["--help"]).stdout
    full = _run(["--help-all"]).stdout
    assert len(short.splitlines()) < 60
    assert len(full.splitlines()) > 150
    assert "--help-all" in short and "--rdf" in short
    assert "--steinhardt" not in short and "--steinhardt" in full
    assert ("よく使う指定" in short or "Only the common options" in short)
    assert "よく使う指定" not in full and "Only the common options" not in full


def test_the_file_count_message_shows_what_is_in_subdirectories(tmp_path, sk_root, monkeypatch):
    from adit.config import save_config

    from tests.conftest import cfg_for

    from tests.conftest import make_fake_skset

    make_fake_skset(sk_root, "mio-1-1", ["H", "O"])
    config = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), config)
    env = {**dict(__import__("os").environ), "ADIT_CONFIG": str(config)}
    spec = tmp_path / "s.json"
    assert subprocess.run([sys.executable, "-m", "adit.cli", "--sample", "water_generated", str(spec)],
                          capture_output=True, text=True, cwd=REPO, timeout=300, env=env).returncode == 0
    r = subprocess.run([sys.executable, "-m", "adit.cli", str(spec), str(tmp_path / "run")],
                       capture_output=True, text=True, cwd=REPO, timeout=300, env=env)
    assert r.returncode == 0, r.stderr + r.stdout
    assert "skf/ の中に" in r.stdout


def test_the_readme_shows_how_to_start():
    text = (REPO / "README.md").read_text(encoding="utf-8")
    assert "## インストール" in text and "## できること" in text
    assert "[チュートリアル](docs/USAGE.md)" in text          # the steps live in the tutorial
    assert text.count("![") >= 4                              # shown with pictures
    assert len(text.splitlines()) < 150
    usage = (REPO / "docs" / "USAGE.md").read_text(encoding="utf-8")
    for command in ("--list-samples", "--sample water_generated", "adit-analyze"):
        assert command in usage, command
    assert usage.count("![") >= 10, "チュートリアルは画面の画像で説明する"
    install = REPO / "docs" / "INSTALL.md"
    assert install.is_file() and "実行ファイルを使う" in install.read_text(encoding="utf-8")


def test_the_web_form_says_the_detailed_fields_are_optional():
    html = (REPO / "src" / "adit" / "web" / "templates" / "analysis.html").read_text(encoding="utf-8")
    assert "すべて任意です" in html
    index = (REPO / "src" / "adit" / "web" / "templates" / "index.html").read_text(encoding="utf-8")
    assert "ADIT は計算を投入しません" in index


def test_the_disabled_run_button_says_how_to_enable_it():
    from adit import runner

    src = Path(runner.__file__).read_text(encoding="utf-8")
    assert "enable_run = true" in src and "投入はしません" in src
    server_src = Path(__file__).resolve().parents[1] / "src" / "adit" / "web" / "server.py"
    assert "from adit.runner import" in server_src.read_text(encoding="utf-8")


def test_the_figure_style_can_be_changed(tmp_path):
    from adit.analysis.plotstyle import PlotStyleError, from_text

    style = from_text(colors="black,#d62728", tick_direction="in", grid="none",
                      spines="left-bottom", line_width="1.6", font_size="9", dpi="300")
    style.apply()
    from matplotlib import rcParams

    assert rcParams["xtick.direction"] == "in" and rcParams["ytick.direction"] == "in"
    assert rcParams["axes.grid"] is False
    assert rcParams["axes.spines.top"] is False and rcParams["axes.spines.left"] is True
    assert rcParams["lines.linewidth"] == 1.6 and rcParams["savefig.dpi"] == 300
    assert "線の色" in style.summary() and "目盛り" in style.summary()
    for bad in ({"colors": "not-a-color"}, {"tick_direction": "sideways"},
                {"grid": "diagonal"}, {"spines": "round"}, {"line_width": "太め"}):
        with pytest.raises(PlotStyleError):
            from_text(**bad).apply() if "colors" in bad or "tick" in str(bad) else from_text(**bad).apply()
    from matplotlib import rcdefaults
    rcdefaults()


def test_the_figure_style_reaches_the_figures(tmp_path):
    plain = tmp_path / "plain"
    styled = tmp_path / "styled"
    assert _run([str(MD), "-o", str(plain)]).returncode == 0
    r = _run([str(MD), "-o", str(styled), "--plot-colors", "black", "--plot-grid", "none",
              "--plot-ticks", "in", "--plot-spines", "left-bottom"])
    assert r.returncode == 0, r.stderr
    assert (plain / "temperature.png").read_bytes() != (styled / "temperature.png").read_bytes()


def test_the_figure_style_is_offered_in_both_screens():
    from adit.gui import analysis_fields as AF

    for key in ("plot_colors", "plot_ticks", "plot_grid", "plot_spines", "plot_line_width",
                "plot_font_size", "plot_dpi"):
        assert key in AF.LABELS
    html = (REPO / "src" / "adit" / "web" / "templates" / "analysis.html").read_text(encoding="utf-8")
    assert '"plot_colors"' in html and '"plot_grid"' in html
