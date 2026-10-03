
import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from adit.analysis import AnalysisOptions, run_analysis
from adit.analysis import report
from adit.results import summarize_run

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"


def _copy(name: str, tmp_path: Path) -> Path:
    d = tmp_path / name
    shutil.copytree(EX / name, d, ignore=shutil.ignore_patterns("analysis"))
    return d


def test_md_is_not_reported_as_unconverged():
    s = summarize_run(EX / "dftb_md_water_generated")
    assert s.task_type == "molecular_dynamics"
    t = s.text()
    assert "未収束" not in t and "収束しません" not in t and "指定の計算が終わりました" in t
    assert "MD のステップ数" in t
    assert "(該当なし)" not in t, "MD の結合長は軌跡の最後のフレームから (出発の構造ではない)"
    assert "終了コード 0" in summarize_run(EX / "dftb_md_water_generated", exit_code=0).status_line()
    assert "0 ではありません" in summarize_run(EX / "dftb_md_water_generated", exit_code=1).status_line()


def test_optimization_says_converged_or_not(tmp_path):
    assert summarize_run(EX / "water_generated").status_line() == "構造最適化: 収束しました"
    d = _copy("water_generated", tmp_path)
    log = d / "output.log"
    log.write_text(log.read_text(encoding="utf-8").replace("Geometry converged", "Geometry did not converge"), encoding="utf-8")
    assert "収束しませんでした" in summarize_run(d).status_line()


def test_nonzero_exit_or_missing_end_mark(tmp_path):
    d = _copy("qe_md_si_generated", tmp_path)
    log = d / "output.log"
    log.write_text(log.read_text(encoding="utf-8").replace("JOB DONE", ""), encoding="utf-8")
    assert "正常終了の印がありません" in summarize_run(d).status_line()


def test_vibrations_has_no_step_count():
    t = summarize_run(EX / "xtb_vib_water_generated").text()
    assert "構造ステップ数" not in t and "指定の計算が終わりました" in t


def test_vasp_without_stdout_reads_oszicar():
    s = summarize_run(EX / "vasp_h2o_generated")
    assert s.geometry_steps == 5 and s.converged


def test_md_notes_temperature_and_frames(tmp_path):
    d = _copy("dftb_md_water_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions(rdf=True, msd=True, skip_frames=5))
    txt = res.summary_text()
    assert "目標 300 K、平均" in txt
    assert "軌跡は 16 フレーム" in txt and "16 フレーム" in txt  # 21 - 5 = 16 < FEW_FRAMES
    assert "RDF (動径分布関数):" in txt and "MSD (平均二乗変位。全原子、16 フレーム)" in txt


def test_selected_msd_species_labels_other_elements_as_reference(tmp_path):
    d = _copy("dftb_md_water_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions(msd=True, msd_species="O"))
    txt = res.summary_text()
    assert "MSD (平均二乗変位。O、" in txt
    assert "参考・全元素の拡散係数" in txt
    assert "重心の移動: 除去しました" in txt
    assert set(res.tables["msd"]["by_element"]) == {"H", "O"}
    error = res.tables["msd"]["D_error"]
    assert error["reason_code"] == "fit_range_not_available_in_all_blocks"
    assert "--msd-fit" in error["reason"]
    assert res.tables["msd"]["D_err_cm2_s"] is None
    assert "ブロック誤差: 出せません" in txt


def test_small_msd_is_not_rounded_to_zero_in_summary(tmp_path):
    d = _copy("qe_md_si_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions(msd=True, msd_species="Si"))
    value = res.tables["msd"]["last_A2"]
    assert 0 < value < 0.001
    assert f"最終 {value:.3g} Å²" in res.summary_text()
    assert "最終 0.000 Å²" not in res.summary_text()


def test_temperature_within_range_has_no_note(tmp_path):
    d = _copy("qe_md_si_generated", tmp_path)
    txt = run_analysis(d, AnalysisOptions()).summary_text()
    assert "目標" not in txt


def test_qe_fermi_note_only_when_used(tmp_path):
    d = _copy("qe_md_si_generated", tmp_path)
    assert "最高被占準位" not in run_analysis(d, AnalysisOptions()).summary_text()
    assert "最高被占準位" in run_analysis(d, AnalysisOptions(dos=True)).summary_text() or True
    b = _copy("qe_si_bands_generated", tmp_path)
    assert "最高被占準位" in run_analysis(b, AnalysisOptions()).summary_text()


@pytest.mark.parametrize("language,label", [
    ("ja", "最高被占準位を基準にした最小の間隔"),
    ("en", "smallest gap relative to the highest occupied level"),
])
def test_qe_band_gap_names_its_actual_reference(tmp_path, monkeypatch, language, label):
    from adit import lang

    monkeypatch.setattr(lang, "LANGUAGE", language)
    d = _copy("qe_si_bands_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions())
    assert res.tables["bands"]["reference_level"] == "highest_occupied"
    assert label in res.summary_text()
    assert "smallest gap across the Fermi level" not in res.summary_text()


def test_rdf_cell_width_note(tmp_path):
    d = _copy("qe_md_si_generated", tmp_path)
    txt = run_analysis(d, AnalysisOptions(rdf=True)).summary_text()
    assert "セルの幅の半分" in txt


def test_single_point_energy_has_no_figure_and_bonds_line_explains(tmp_path):
    d = _copy("xtb_vib_water_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions())
    assert "energy" not in res.figures and "energy" in res.tables and "spectrum" in res.figures
    txt = res.summary_text()
    assert "結合長 (最終構造): 該当なし (" in txt
    assert "力の最大値は" in txt and "eV/Å 以上" in txt and "最適化済みの構造" in txt


def test_vibration_note_without_forces(tmp_path):
    d = _copy("xtb_vib_water_generated", tmp_path)
    log = d / "output.log"
    log.write_text(log.read_text(encoding="utf-8").replace("GRADIENT NORM", "GRAD"), encoding="utf-8")
    txt = run_analysis(d, AnalysisOptions()).summary_text()
    assert "元の構造の力を読めません" in txt


def test_negative_diffusion_note(tmp_path, monkeypatch):
    d = _copy("dftb_md_water_generated", tmp_path)
    from adit.analysis import compute
    real = compute.msd_analysis
    def negative(*args, **kwargs):
        result = real(*args, **kwargs); result["D_cm2_s"] = -1.0e-7; return result
    monkeypatch.setattr(compute, "msd_analysis", negative)
    txt = run_analysis(d, AnalysisOptions(msd=True)).summary_text()
    assert "拡散係数が負の値" in txt


def test_rdf_axis_is_clipped_for_bonded_peaks(tmp_path, monkeypatch):
    rng = np.random.default_rng(0)
    L = 12.0
    frames = []
    for _ in range(30):
        pos, sym = [], []
        for c in rng.uniform(0, L, size=(20, 3)):
            pos += [c, c + [0.96, 0, 0], c + [-0.24, 0.93, 0]]; sym += ["O", "H", "H"]
        frames.append(Atoms(sym, positions=pos, cell=[L, L, L], pbc=True))
    from adit.analysis import readers
    fake = readers.RunData("dftbplus", tmp_path, frames=frames)
    monkeypatch.setattr(report, "load_run", lambda _d, _code=None: fake)
    res = run_analysis(tmp_path, AnalysisOptions(energy=False, temperature=False, bonds=False, rdf=True, rdf_rmax=6.0))
    t = res.tables["rdf"]
    assert t["ylim_clipped"] is not None and t["g_max"] > report.RDF_CLIP_RATIO * t["ylim_clipped"] / 1.2
    assert "縦軸を" in res.summary_text()
    data = json.loads((tmp_path / "analysis" / "rdf.json").read_text(encoding="utf-8"))
    assert max(data["H-O"]["g"]) == pytest.approx(t["g_max"]) or max(max(v["g"]) for v in data.values()) == pytest.approx(t["g_max"])


def test_figure_titles_ja_en():
    from adit import lang
    assert report.figure_title("rdf") == "動径分布関数 (RDF)"
    lang.set_language("en")
    try:
        assert report.figure_title("energy") == "Energy"
    finally:
        lang.set_language("ja")


@pytest.fixture
def qapp():
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


def test_panel_md_defaults_titles_and_species(qapp, tmp_path):
    from adit.gui.panels.analysis_panel import AnalysisPanel, FigureView
    p = AnalysisPanel()
    d = _copy("dftb_md_water_generated", tmp_path)
    p.set_run_dir(d, ["H", "O"])
    assert p.cb_rdf.isChecked() and p.cb_msd.isChecked(), "MD なら RDF と MSD を既定で入れる"
    assert p.options().msd_species is None
    p.msd_species.setCurrentIndex(2)
    assert p.options().msd_species == "O"
    res = p.run()
    assert res is not None
    titles = [p.figs_lay.itemAt(i).widget().text() for i in range(0, p.figs_lay.count(), 2)]
    assert "エネルギーの推移" in titles and "動径分布関数 (RDF)" in titles and not any("/" in t for t in titles)
    view = p.figs_lay.itemAt(1).widget()
    assert isinstance(view, FigureView) and view.toolTip().endswith(".png")
    assert view.heightForWidth(view.pixmap.width() // 2) == pytest.approx(view.pixmap.height() / 2, abs=1)
    v = _copy("xtb_vib_water_generated", tmp_path)
    p.set_run_dir(v, ["H", "O"])
    assert not p.cb_rdf.isChecked() and not p.cb_msd.isChecked()


def test_panel_all_atoms_label_english(qapp):
    from adit import lang
    lang.set_language("en")
    try:
        from adit.gui.panels.analysis_panel import AnalysisPanel
        p = AnalysisPanel()
        assert p.msd_species.itemText(0) == "(all atoms)" and p.options().msd_species is None
    finally:
        lang.set_language("ja")


def test_the_unoptimized_example_carries_a_note():
    note = EX / "xtb_vib_water_generated" / "NOTE.md"
    assert note.is_file()
    text = note.read_text(encoding="utf-8")
    assert "構造最適化をしていません" in text
    assert "1595" in text and "3657" in text
