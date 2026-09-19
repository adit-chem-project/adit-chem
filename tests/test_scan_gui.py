
import json
import os
import shutil
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QLabel, QMessageBox  # noqa: E402

from tests.conftest import pbs_profile  # noqa: E402
from adit.analysis.report import figure_title  # noqa: E402
from adit.config import default_config  # noqa: E402
from adit.gui.main_window import MainWindow  # noqa: E402
from adit.gui.panels.analysis_panel import AnalysisPanel  # noqa: E402
from adit.gui.scan_dialog import OTHER, choices_for  # noqa: E402

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def boxes(monkeypatch):
    seen = []
    for name in ("information", "critical", "warning"):
        monkeypatch.setattr(QMessageBox, name, staticmethod(lambda *a, _n=name, **k: seen.append((_n, a[2] if len(a) > 2 else ""))))
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))
    return seen


def make_window(sk_root, tmp_path) -> MainWindow:
    cfg = default_config(sk_root=str(sk_root))
    cfg.profiles["cluster"] = pbs_profile()
    win = MainWindow(cfg, tmp_path / "cluster.toml")
    win.method.sk_set.setCurrentText("fake-1-0")
    win.runtime.outdir.setText(str(tmp_path / "out"))
    win.refresh_preview()
    return win


def paths(code, periodic):
    return [c.path for c in choices_for(code, periodic)]


def test_choices_follow_code_and_periodicity():
    assert paths("espresso", True) == ["method.ecutwfc", "kpoints.mesh", "kpoints.density", "scale", OTHER]
    assert paths("vasp", True)[0] == "method.encut" and paths("vasp", False) == ["method.encut", OTHER]
    assert paths("dftbplus", True) == ["kpoints.mesh", "kpoints.density", "scale", OTHER]
    for code in ("dftbplus", "xtb", "orca"):
        assert paths(code, False) == [OTHER]
    for c in choices_for("espresso", True):
        assert c.example and c.purpose


def test_density_choice_switches_the_kpoint_mode():
    from ase.build import bulk
    from tests.conftest import water_spec
    from adit.scan import apply_value
    from adit.spec import AtomsData, KPoints, Structure
    spec = water_spec().model_copy(update={"structure": Structure(source="bulk", source_ref="Si", atoms=AtomsData.from_ase(bulk("Si", a=5.43))),
                                           "kpoints": KPoints(mode="mesh", mesh=(2, 2, 2))})
    new = apply_value(spec, "kpoints.density", "3")
    assert new.kpoints.mode == "density" and new.kpoints.density == 3.0


def test_dialog_in_window_uses_current_code_and_structure(app, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    dlg = win.scan_dialog()
    assert [dlg.item.itemData(i) for i in range(dlg.item.count())] == [OTHER]
    assert dlg.folder.text() == str(tmp_path / "out") + "_scan" and dlg.values.text() == ""
    assert not dlg.path.isHidden()
    win.structure.set_source("bulk"); win._on_context(); win.refresh_preview()
    win.method.code.setCurrentIndex(win.method.code.findData("espresso"))
    dlg = win.scan_dialog()
    assert [dlg.item.itemData(i) for i in range(dlg.item.count())] == paths("espresso", True)
    dlg.item.setCurrentIndex(dlg.item.findData("scale"))
    assert dlg.path.isHidden() and "5" in dlg.purpose.text() and "1.00" in dlg.values.placeholderText()
    assert dlg.scan_text().startswith("scale=")


def test_generate_writes_each_value_and_fills_the_analysis_tab(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    dlg = win.scan_dialog()
    dlg.path.setText("method.max_scc_iterations"); dlg.values.setText("50, 100, 200")
    dlg.generate()
    out = tmp_path / "out_scan"
    assert [d.name for d in dlg.dirs] == ["max_scc_iterations_50", "max_scc_iterations_100", "max_scc_iterations_200"]
    assert all((d / "dftb_in.hsd").is_file() for d in dlg.dirs) and (out / "scan.json").is_file()
    assert boxes[-1][0] == "information" and "3" in boxes[-1][1] and "submit.sh" in boxes[-1][1]
    win.scan_written(dlg.out_dir)
    assert win.analysis.run_dir.text() == str(out)
    assert not win.analysis.scan_note.isHidden() and not win.analysis.cb_rdf.isEnabled()
    win.analysis.set_run_dir(tmp_path / "out")
    assert win.analysis.scan_note.isHidden() and win.analysis.cb_energy.isEnabled()


def test_generate_error_is_shown_in_words_and_writes_nothing(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    dlg = win.scan_dialog()
    dlg.path.setText("method.no_such"); dlg.values.setText("1, 2")
    dlg.generate()
    assert boxes[-1][0] == "critical" and "method.no_such" in boxes[-1][1] and "Error" not in boxes[-1][1]
    assert dlg.out_dir is None and not any((tmp_path / "out_scan").glob("*_1"))
    dlg.path.setText("method.max_scc_iterations"); dlg.values.setText("50")
    dlg.generate()
    assert boxes[-1][0] == "critical" and "2" in boxes[-1][1]
    dlg.path.setText("")
    dlg.generate()
    assert boxes[-1][0] == "warning"


def test_analysis_tab_shows_table_and_figures_for_a_scan(app, tmp_path):
    out = tmp_path / "scan"
    for name in ("x_1", "x_2"):
        shutil.copytree(REPO / "examples" / "dftb_tio2_generated", out / name)
    (out / "scan.json").write_text(json.dumps({"path": "method.x", "values": ["1", "2"], "dirs": ["x_1", "x_2"]}), encoding="utf-8")
    panel = AnalysisPanel()
    panel.set_run_dir(out)
    res = panel.run()
    assert res is not None and len(res.rows) == 2
    text = panel.summary.toPlainText()
    assert "method.x" in text and "-425.07" not in text and "\t" not in text
    t = panel.scan_table
    assert not t.isHidden() and t.rowCount() == 2 and t.columnCount() == 6
    heads = [t.horizontalHeaderItem(j).text() for j in range(6)]
    assert heads[0] == "値" and heads[2] == "最後の値との差\n[meV/原子]" and heads[4] == "圧力\n[GPa]"
    assert t.item(0, 1).text() == res.rows[0]["energy_ev"] and t.item(0, 4).text() == res.rows[0]["pressure_gpa"]
    assert t.editTriggers() == t.EditTrigger.NoEditTriggers
    assert panel.summary.toPlainText() and res.summary_text().count("\n") > text.count("\n")
    from adit.gui.panels.analysis_panel import FigureHeader

    titles = [panel.figs_lay.itemAt(i).widget().text() for i in range(panel.figs_lay.count())
              if isinstance(panel.figs_lay.itemAt(i).widget(), (QLabel, FigureHeader))]
    assert titles == [figure_title("scan_energy")] and titles[0] != "scan_energy"
    assert (out / "scan_energy.png").stat().st_size > 1000
    plain = tmp_path / "plain"; shutil.copytree(REPO / "examples" / "dftb_tio2_generated", plain)
    panel.set_run_dir(plain); panel.run()
    assert t.isHidden()


def test_scan_overwrite_keeps_a_pruned_backup(app, boxes, sk_root, tmp_path, monkeypatch):
    from adit.project import BACKUP_DIR, restore_backup

    asked = []
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: (asked.append(a[2]), QMessageBox.StandardButton.Yes)[1]))
    win = make_window(sk_root, tmp_path)
    dlg = win.scan_dialog()
    dlg.path.setText("method.max_scc_iterations"); dlg.values.setText("50, 100")
    dlg.generate()
    out = tmp_path / "out_scan"
    assert not asked and dlg.backup_dir is None
    (out / "max_scc_iterations_50" / "dftb_in.hsd").write_text("old\n", encoding="utf-8")
    (out / "extra.txt").write_text("keep\n", encoding="utf-8")
    dlg = win.scan_dialog()
    dlg.path.setText("method.max_scc_iterations"); dlg.values.setText("50, 100")
    dlg.generate()
    assert len(asked) == 1 and "同じ名前のファイルを上書きします" in asked[0] and BACKUP_DIR in asked[0]
    d = dlg.backup_dir
    assert d is not None and d.parent == out / BACKUP_DIR
    saved = sorted(p.relative_to(d).as_posix() for p in d.rglob("*") if p.is_file())
    assert "max_scc_iterations_50/dftb_in.hsd" in saved and "extra.txt" not in saved
    win.offer_undo(d)
    assert not win.undo_link.isHidden()
    restore_backup(d)
    assert (out / "max_scc_iterations_50" / "dftb_in.hsd").read_text(encoding="utf-8") == "old\n"
