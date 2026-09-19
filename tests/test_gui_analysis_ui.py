
import importlib.util
import json
import os
import shutil
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QDialog  # noqa: E402

from adit.analysis import AnalysisResult  # noqa: E402
from adit.analysis import trajectory as trj  # noqa: E402
from adit.analysis.report import figure_title  # noqa: E402
from adit.gui import analysis_fields as AF  # noqa: E402
from adit.gui.panels.analysis_panel import AnalysisPanel  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture(autouse=True)
def _delete_windows_after_test():
    yield
    from PySide6.QtCore import QCoreApplication, QEvent
    a = QApplication.instance()
    if a is not None:
        for w in a.topLevelWidgets():
            if w.parentWidget() is None:
                w.close(); w.deleteLater()
        QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


def _copy(name: str, dest: Path) -> Path:
    d = dest / name
    shutil.copytree(EX / name, d, ignore=shutil.ignore_patterns("analysis"))
    return d


def _panel(d: Path, elements=None) -> AnalysisPanel:
    p = AnalysisPanel()
    p.set_run_dir(d, elements or [])
    return p


def _view(p: AnalysisPanel, key: str):
    views = {v.section.key: v for v in p.section_views()}
    assert key in views, list(views)
    return views[key]


def _fig_titles(p: AnalysisPanel) -> list[str]:
    return [p.figs_lay.itemAt(i).widget().text() for i in range(0, p.figs_lay.count(), 2)]


def test_md_tables_charges_statistics_and_diffusion(app, tmp_path):
    d = _copy("dftb_md_water_generated", tmp_path)
    p = _panel(d, ["O", "H"])
    assert p.cb_rdf.isChecked() and p.cb_msd.isChecked()
    assert not p.more.is_expanded()
    p.stride.setValue(2); p.msd_fit_from.setText("10"); p.msd_fit_to.setText("50")
    res = p.run()
    assert res is not None, p.summary.toPlainText()
    assert res.tables["trajectory"]["stride"] == 2 and res.tables["msd"]["fit_range_user"]
    ch = _view(p, "charges0")
    notes = " ".join(w.text() for w in ch.notes)
    assert "Mulliken" in notes and "detailed.out:" in notes
    assert [ch.table.item(i, 0).text() for i in range(3)] == ["O1", "H2", "H3"]
    tr = _view(p, "trajectory")
    cells = {tr.table.item(i, 0).text(): tr.table.item(i, 1).text() for i in range(tr.table.rowCount())}
    assert cells["間引き (N フレームおき)"] == "2" and cells["使ったフレーム数"] == "11" and cells["読んだファイル"] == "geo_end.xyz"
    msd = _view(p, "msd")
    assert [msd.table.item(i, 0).text() for i in range(msd.table.rowCount())] == ["全原子", "H", "O"]
    st = _view(p, "timeseries_stats")
    assert st.table.rowCount() >= 1 and "ブロック平均" in st.table.horizontalHeaderItem(6).text()
    titles = _fig_titles(p)
    for name in ("coordination", "blocking", "msd", "rdf"):
        assert figure_title(name) in titles, (name, titles)


def test_zdensity_pressure_and_spacegroup_for_lammps(app, tmp_path):
    d = _copy("lammps_cu_nvt_generated", tmp_path)
    p = _panel(d)
    p.cb_zdens.setChecked(True); p.zdens_bin.setValue(0.5)
    res = p.run()
    assert res is not None, p.summary.toPlainText()
    z = res.tables["zdensity"]
    assert z["bin_A"] == pytest.approx(z["axis_length_A"] / round(z["axis_length_A"] / 0.5))
    titles = _fig_titles(p)
    assert figure_title("zdensity") in titles and figure_title("pressure") in titles
    if importlib.util.find_spec("spglib"):
        sg = _view(p, "spacegroup")
        assert sg.table.rowCount() == 3
    p.symprec.setText("0.01")
    res = p.run()
    if importlib.util.find_spec("spglib"):
        assert [r["symprec_A"] for r in res.tables["spacegroup"]["results"]] == [0.01]


def test_thermo_has_no_defaults_and_shows_the_reason(app, tmp_path):
    d = _copy("xtb_vib_water_generated", tmp_path)
    p = _panel(d)
    assert p.th_model.currentData() == "" and not p.th_temps.text() and not p.th_pressure.text() and p.th_geometry.currentData() == ""
    res = p.run()
    assert "thermo_ase" not in res.tables
    code = _view(p, "thermochemistry")
    assert any(code.table.item(i, 3).text().startswith("output.log:") for i in range(code.table.rowCount()))
    el = _view(p, "electronic")
    assert "6.3383" in el.table.item(0, 1).text() and el.table.item(0, 2).text().startswith("xtbout.json:")
    p.th_model.setCurrentIndex(p.th_model.findData("ideal_gas"))
    res = p.run()
    th = _view(p, "thermo_ase")
    assert not res.tables["thermo_ase"]["computed"] and th.table is None
    assert "圧力" in " ".join(w.text() for w in th.notes)
    p.th_temps.setText("298.15, 400"); p.th_pressure.setText("100000"); p.th_sigma.setText("2")
    p.th_geometry.setCurrentIndex(p.th_geometry.findData("nonlinear")); p.th_spin.setText("0")
    res = p.run()
    t = res.tables["thermo_ase"]
    assert t["computed"] and [r["T_K"] for r in t["rows"]] == [298.15, 400.0]
    th = _view(p, "thermo_ase")
    assert th.table.rowCount() == 2 and th.table.horizontalHeaderItem(0).text() == "T [K]"
    assert th.table.item(0, 1).text() == "100000"
    assert figure_title("spectrum") in _fig_titles(p)


def test_unreadable_fields_name_the_field(app, tmp_path):
    d = _copy("xtb_vib_water_generated", tmp_path)
    p = _panel(d)
    p.th_model.setCurrentIndex(p.th_model.findData("harmonic")); p.th_exclude.setText("abc")
    assert p.run() is None and p.summary.toPlainText().startswith("解析の条件を読めません: 除く低い振動の本数")
    p.th_exclude.setText(""); p.uv_shape.setCurrentIndex(p.uv_shape.findData("gauss"))
    assert p.run() is None and "UV-Vis の半値全幅" in p.summary.toPlainText()
    p.uv_shape.setCurrentIndex(0); p.msd_fit_from.setText("30"); p.msd_fit_to.setText("10")
    assert p.run() is None and "MSD の当てはめ範囲" in p.summary.toPlainText()


def test_too_large_trajectory_offers_the_stride(app, tmp_path, monkeypatch):
    d = _copy("dftb_md_water_generated", tmp_path)
    monkeypatch.setattr(trj.Trajectory, "file_size", lambda self: 50 * 2**30)
    p = _panel(d)
    p.cb_rdf.setChecked(False); p.cb_msd.setChecked(True)
    assert p.run() is None
    assert p.too_large.isVisibleTo(p) and "--stride" in p.too_large_text.text()
    n = p._suggested_stride
    assert n and n > 1 and p.btn_use_stride.text() == f"間引きを {n} にする"
    assert p.too_large_text.text() in p.summary.toPlainText()
    p.btn_use_stride.click()
    assert p.stride.value() == n and p.more.is_expanded() and not p.too_large.isVisibleTo(p)
    assert p.run() is not None


def test_export_shows_readme_and_opens_the_folder(app, tmp_path, monkeypatch):
    import adit.gui.panels.analysis_panel as ap
    d = _copy("cp2k_h2o_md_generated", tmp_path)
    p = _panel(d)
    assert not p.export_box.isVisibleTo(p)
    p.cb_unwrap.setChecked(True)
    opened = []
    monkeypatch.setattr(ap, "open_folder", lambda path: opened.append(str(path)) or True)
    res = p.export()
    assert res is not None, p.summary.toPlainText()
    exp = Path(res.tables["export"]["dir"])
    assert res.tables["export"]["unwrap_requested"] and not res.tables["export"]["unwrap_molecules"]
    assert (exp / "trajectory.extxyz").is_file()
    assert p.export_box.isVisibleTo(p) and p.export_readme.toPlainText() == (exp / "export_README.txt").read_text(encoding="utf-8")
    assert str(exp) in p.export_path.text()
    p.btn_open_export.click()
    assert opened == [str(exp)] and p.last_open_ok
    p.run()
    assert not p.export_box.isVisibleTo(p)


def test_open_folder_per_os(app, tmp_path, monkeypatch):
    from adit.gui import analysis_views as av
    calls = []
    monkeypatch.setattr(av, "_is_wsl", lambda: False)
    monkeypatch.setattr(av.QDesktopServices, "openUrl", staticmethod(lambda url: calls.append(url.toLocalFile()) or True))
    assert av.open_folder(tmp_path) and [Path(c) for c in calls] == [tmp_path]
    assert not av.open_folder(tmp_path / "absent")
    monkeypatch.setattr(av, "_is_wsl", lambda: True)
    monkeypatch.setattr(av.shutil, "which", lambda name: "/mnt/c/WINDOWS/explorer.exe" if name == "explorer.exe" else None)

    class Done:
        stdout = "\\\\wsl.localhost\\Ubuntu\\tmp\\x\n"

    monkeypatch.setattr(av.subprocess, "run", lambda *a, **k: Done())
    popen = []
    monkeypatch.setattr(av.subprocess, "Popen", lambda args: popen.append(args))
    assert av.open_folder(tmp_path) and popen == [["explorer.exe", "\\\\wsl.localhost\\Ubuntu\\tmp\\x"]]


def test_compare_dialog_and_tables(app, tmp_path):
    from adit.gui.compare_dialog import CompareDialog
    for n in ("xtb_water_generated", "xtb_vib_water_generated", "water_generated"):
        _copy(n, tmp_path)
    dlg = CompareDialog(str(tmp_path))
    dlg.set_rows([("same", "1", "xtb_water_generated"), ("", "-1", "xtb_vib_water_generated"), ("unbal", "1", "water_generated")])
    dlg._try_accept()
    assert dlg.result() == QDialog.DialogCode.Accepted
    rx = dlg.reactions()
    assert [r.name for r in rx] == ["same", "unbal"] and rx[0].terms == [(1.0, "xtb_water_generated"), (-1.0, "xtb_vib_water_generated")]
    p = AnalysisPanel()
    cres = p.run_compare(dlg.base_dir(), rx)
    assert cres is not None
    views = p.section_views()
    assert [v.section.key for v in views] == ["compare_reactions", "compare_runs", "compare_conditions"]
    reac = views[0]
    assert reac.table.rowCount() == 2 and reac.table.item(0, 10).text() == "釣り合う"
    assert "H:+2 O:+1" in reac.table.item(1, 10).text()
    assert "thermo.csv" in reac.table.item(0, 6).text() and reac.table.item(0, 7).text() == "-"
    assert float(reac.table.item(0, 2).text()) == pytest.approx(cres.reactions[0]["delta_e_ev"], abs=1e-6)
    cond = views[2]
    assert [cond.table.horizontalHeaderItem(j).text() for j in range(cond.table.columnCount())] == \
        ["項目", "xtb_water_generated", "xtb_vib_water_generated", "water_generated"]
    assert "task.type" in [cond.table.item(i, 0).text() for i in range(cond.table.rowCount())]
    assert figure_title("compare_energy") in _fig_titles(p)
    assert "比べた計算: 3 個" in p.summary.toPlainText()
    bad = CompareDialog(str(tmp_path))
    bad.set_rows([("x", "abc", "water_generated")]); bad._try_accept()
    assert bad.result() != QDialog.DialogCode.Accepted and "係数" in bad.error.text()
    (tmp_path / "compare.json").write_text(json.dumps({"reactions": [{"name": "same", "terms": [
        {"dir": "xtb_water_generated", "nu": 1}, {"dir": "xtb_vib_water_generated", "nu": -1}]}]}), encoding="utf-8")
    j = CompareDialog(str(tmp_path))
    assert j.btn_load.isEnabled() and j.load_compare_json()
    assert j.rows()[:2] == [("same", "1", "xtb_water_generated"), ("", "-1", "xtb_vib_water_generated")]


@pytest.mark.parametrize("lang", ["ja", "en"])
def test_analysis_tab_fits_at_1366(app, sk_root, tmp_path, lang):
    from adit.config import default_config
    from adit.gui.i18n import set_language, translate_widgets
    from adit.gui.main_window import MainWindow
    from adit.gui.style import apply_base_style
    apply_base_style(app, "light")
    set_language(lang)
    try:
        d = _copy("dftb_md_water_generated", tmp_path)
        win = MainWindow(default_config(sk_root=str(sk_root)), tmp_path / "cluster.toml"); translate_widgets(win)
        win.resize(1366, 768); win.show(); win.set_mode(win.MODE_ANALYSIS)
        p = win.analysis
        p.set_run_dir(d, ["O", "H"]); p.more.set_expanded(True); assert p.run() is not None
        for _ in range(6):
            app.processEvents()
        need, have = p.scroll.widget().minimumSizeHint().width(), p.scroll.viewport().width()
        assert need <= have, (lang, need, have)
    finally:
        set_language("ja")


def test_long_tables_show_the_head_and_the_file(tmp_path):
    res = AnalysisResult(code="x", run_dir=str(tmp_path))
    res.tables["charges"] = [{"definition": "Mulliken", "values": [0.01 * i for i in range(30)], "source": "detailed.out:13", "note": "", "sum": 4.35}]
    sec = AF.result_sections(res)[0]
    assert len(sec.rows) == AF.ROW_LIMIT and sec.total_rows == 30
    assert "全 30 行" in sec.more_text() and str(tmp_path / "analysis" / "summary.json") in sec.more_text()


def test_fields_roundtrip():
    f = {"rdf": "on", "rmax": "6", "stride": "3", "msd_fit_from": "100", "msd_fit_to": "900", "zdens": "on", "zdens_bin": "0.25",
         "stats": "", "memory_mb": "512", "symprec": "1e-4, 0.01", "th_model": "harmonic", "th_temps": "300", "th_exclude": "6",
         "th_imag": "ignore", "uv_shape": "lorentz", "uv_fwhm": "0.2", "export_unwrap": "on"}
    o = AF.options_from_fields(f)
    assert (o.stride, o.msd_fit_fs, o.zdens, o.zdens_bin_ang, o.stats, o.memory_budget_mb, o.symprecs) == (3, (100.0, 900.0), True, 0.25, False, 512.0, (1e-4, 0.01))
    assert o.thermo.model == "harmonic" and o.thermo.exclude_lowest == 6 and o.thermo.pressure_pa is None and o.uvvis_broadening == ("lorentz", 0.2)
    back = AF.fields_from_options(o)
    assert AF.options_from_fields(back) == o
    assert AF.details_changed(back) and not AF.details_changed(AF.fields_from_options(AF.options_from_fields({"stats": "on"})))


@pytest.mark.parametrize("language,range_text,formula_text", [
    ("ja", "既定: 最大ずれ時間の 10〜50 %", "MSD = 4 D t + c"),
    ("en", "default: 10–50% of the maximum lag time", "MSD = 4 D t + c"),
])
def test_msd_screen_text_matches_fit_default_and_dimension(tmp_path, monkeypatch, language, range_text, formula_text):
    from adit import lang
    from adit.gui.help import help_for

    monkeypatch.setattr(lang, "LANGUAGE", language)
    res = AnalysisResult(code="test", run_dir=str(tmp_path))
    res.tables["msd"] = {"species": "O", "D_cm2_s": 1e-6, "last_A2": 0.002,
                         "fit_range_fs": [10, 50], "fit_range_user": False,
                         "dimension": 2, "formula": "MSD = 4 D t + c", "by_element": {}}
    section = next(s for s in AF.result_sections(res) if s.key == "msd")
    assert range_text in " ".join(section.notes)
    assert formula_text in " ".join(section.notes)
    assert "後半" not in " ".join(section.notes) and "second half" not in " ".join(section.notes)
    assert "10" in AF.ph("msd_fit")
    h = help_for("MSD の当てはめ範囲 [fs]")
    assert h and "10" in (h.ja if language == "ja" else h.en)


def test_labels_have_english_and_help(app):
    from adit.gui.help import help_for
    from adit.gui.i18n import _EN, set_language, translate_widgets
    rows = ("stride", "msd_fit", "memory_mb", "zdens", "zdens_bin", "stats", "pdos", "symprec", "th_model", "th_temps", "th_pressure", "th_sigma",
            "th_geometry", "th_spin", "th_imag", "th_exclude", "th_qh", "th_tau", "uv_shape", "uv_fwhm", "export_unwrap", "compare_base")
    for key in rows:
        ja = AF.LABELS[key][0]
        h = help_for(ja)
        assert h is not None and h.ja and h.en, key
        assert h.required == (key == "compare_base"), key
    for ja, en in AF.LABELS.values():
        assert ja in _EN, ja
    set_language("en")
    try:
        p = AnalysisPanel(); translate_widgets(p)
        assert p.more.button.text() == "More options" and p.btn_export.text() == "Export for TRAVIS, VMD and OVITO"
        assert p.th_model.itemText(0) == "(not computed)" and p.cb_unwrap.text() == "Make molecules whole"
        assert p.btn_compare.text() == "Compare runs…"
    finally:
        set_language("ja")
