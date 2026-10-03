"""Fixes from the desktop-GUI audit of 2026-10-03: measurements, undo, the changed-only filter, the workspace, i18n."""

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtCore import QCoreApplication, QEvent, Qt  # noqa: E402
from PySide6.QtWidgets import QApplication, QComboBox, QFileDialog, QFormLayout, QMessageBox  # noqa: E402

from adit.gui.panels.method_panel import CODES  # noqa: E402
from adit.measure import distance, measure  # noqa: E402
from adit.spec import AtomsData, KPoints, Structure  # noqa: E402
from tests.conftest import water_spec  # noqa: E402
from tests.test_gui import make_window  # noqa: E402

JA = re.compile(r"[぀-ヿ㐀-鿿]")
SHELL = ["cmd.exe"] if os.name == "nt" else ["/bin/sh"]


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def quiet(monkeypatch):
    monkeypatch.setattr(QMessageBox, "information", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda *a, **k: None))
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))


@pytest.fixture
def english():
    from adit.gui.i18n import set_language
    set_language("en")
    yield
    set_language("ja")


def _flush_deletes() -> None:
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


class FakeSession:
    """Stands in for a ShellSession: records what is written, reports the state the test sets."""

    rows, cols = 24, 80
    cursor = (0, 0)
    history_above = history_below = 0

    def __init__(self, alive: bool = True, busy: bool = False, rows: int = 24, cols: int = 80):
        self.alive, self._busy, self.rows, self.cols = alive, busy, rows, cols
        self.written: list[str] = []

    def busy(self) -> bool:
        return self._busy

    def write(self, text) -> None:
        self.written.append(text)

    def lines(self):
        from adit.gui.terminal_backend import Cell
        return [[Cell("a", "default", "default", False, False) for _ in range(self.cols)] for _ in range(self.rows)]

    def drain(self) -> bool:
        return False

    def resize(self, rows: int, cols: int) -> None:
        pass

    def close(self) -> None:
        pass

    def text(self) -> str:
        return ""

    def mouse_wanted(self) -> bool:
        return False

    def mouse_motion_wanted(self) -> bool:
        return False


# ---- H1: the minimum image follows the periodic directions only ----------------

def test_measurements_do_not_wrap_across_the_vacuum_of_a_slab(app):
    from adit.gui.viewer3d import Viewer3D

    cell = np.diag([10.0, 10.0, 22.0])
    pos = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 12.0]])
    assert distance(pos, 0, 1, cell, pbc=(True, True, False)) == pytest.approx(12.0)
    assert distance(pos, 0, 1, cell, pbc=True) == pytest.approx(10.0)              # fully periodic: the image is closer
    assert measure(pos, [0, 1], cell, pbc=[True, True, False])["value"] == pytest.approx(12.0)
    slab = Atoms("CC", positions=pos, cell=cell, pbc=(True, True, False))
    v = Viewer3D(); v.set_atoms(slab); v.set_selection([0, 1])
    assert v.measurement()["value"] == pytest.approx(12.0)
    assert slab.get_distance(0, 1, mic=True) == pytest.approx(v.measurement()["value"])
    # an angle across the vacuum is measured on the real atoms as well
    three = Atoms("CCC", positions=[(0, 0, 0), (0, 0, 12.0), (3.0, 0, 12.0)], cell=cell, pbc=(True, True, False))
    v.set_atoms(three); v.set_selection([0, 1, 2])
    assert v.measurement()["value"] == pytest.approx(three.get_angle(0, 1, 2, mic=True))
    # in-plane wrapping still works on the same slab
    v.set_atoms(Atoms("CC", positions=[(0.5, 5, 5), (9.5, 5, 5)], cell=cell, pbc=(True, True, False))); v.set_selection([0, 1])
    assert v.measurement()["value"] == pytest.approx(1.0)


# ---- M1: undo inside the preview debounce --------------------------------------

def test_undo_right_after_an_edit_keeps_the_edit_in_the_history(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.task.max_steps.setValue(77); win._timer.stop(); win.refresh_preview()
    n = len(win._history)
    win.task.max_steps.setValue(88)
    assert win._timer.isActive()                                        # the preview has not caught up yet
    win.go_back()
    assert win.task.max_steps.value() == 77 and len(win._history) == n + 1
    assert win.act_forward.isEnabled()
    win.go_forward()
    assert win.task.max_steps.value() == 88
    win.task.max_steps.setValue(99)
    win.go_forward()                                                    # nothing to redo: the pending edit is only recorded
    assert win.task.max_steps.value() == 99 and not win.act_forward.isEnabled()
    win.close()


# ---- M2: the changed-only filter survives an inline error row --------------------

def _visible_anchors(form) -> set:
    out = set()
    for r in range(form.rowCount()):
        if not form.isRowVisible(r):
            continue
        for role in (QFormLayout.ItemRole.LabelRole, QFormLayout.ItemRole.SpanningRole, QFormLayout.ItemRole.FieldRole):
            item = form.itemAt(r, role)
            if item is not None and item.widget() is not None:
                out.add(item.widget()); break
    return out


def test_filter_off_shows_every_row_after_an_inline_error_was_inserted(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    win.set_mode(win.MODE_SETTINGS)
    win.method.code.setCurrentIndex(list(CODES).index("lammps")); win._timer.stop(); win.refresh_preview()
    lp = win.method.lammps; form = lp.layout()
    good = tmp_path / "ok.data"
    good.write_text("LAMMPS data\n\n3 atoms\n2 atom types\n\n0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\nMasses\n\n1 15.999\n2 1.008\n\n"
                    "Atoms\n\n1 1 0 0 0\n2 2 1 0 0\n3 2 0 1 0\n", encoding="utf-8")
    lp.units.setCurrentIndex(1); lp.pair_style.setText("lj/cut 10.0"); lp.pair_coeff.setPlainText("* * 0.1 3.0")
    lp.data_file.setText(str(good)); lp.type_elements.setText("O H")
    win._timer.stop(); win.refresh_preview()
    assert win._errors == [], win._errors
    before = _visible_anchors(form)
    rows = form.rowCount()
    win.filter_changed.setChecked(True)
    assert len(_visible_anchors(form)) < len(before)
    lp.data_file.setText(str(tmp_path / "missing.data")); win._timer.stop(); win.refresh_preview()
    assert any("data_file" in loc for loc in win._error_locations), win._error_locations
    assert form.rowCount() == rows + 1                                  # the inline error row was inserted
    win.filter_changed.setChecked(False)
    after = _visible_anchors(form)
    assert before <= after, [w.text() if hasattr(w, "text") else w for w in before - after]
    # the only rows still hidden are inline error lines of errors that were cleared earlier
    for r in range(form.rowCount()):
        if not form.isRowVisible(r):
            field = form.itemAt(r, QFormLayout.ItemRole.FieldRole)
            assert field is not None and field.widget() is not None and field.widget().objectName() == "field_error", r
    win.close()


# ---- M3 / L4: cd into the terminal --------------------------------------------

@pytest.mark.skipif(os.name == "nt", reason="the shell on Windows is different")
def test_cd_is_not_typed_into_a_busy_or_dead_terminal(app, tmp_path):
    from adit.gui.panels.workspace_panel import WorkspacePanel

    panel = WorkspacePanel(root=tmp_path)
    term = panel.terminal.current; real = term.session
    msgs = []; panel.status.connect(msgs.append)
    try:
        fake = FakeSession(busy=True); term.session = fake
        assert panel.cd_to(tmp_path) is False and fake.written == [] and msgs and "cd" in msgs[-1]
        assert panel.terminal.send("x") is False
        fake._busy = False
        assert panel.cd_to(tmp_path) is True
        assert fake.written == [f"cd {shlex.quote(str(tmp_path))}\r"]
        fake.alive = False
        assert panel.cd_to(tmp_path) is False and len(fake.written) == 1 and len(msgs) == 2
        term.session = None
        assert panel.cd_to(tmp_path) is False
    finally:
        term.session = real
        panel.close_session()


def test_cd_is_quoted_for_the_windows_shells(monkeypatch):
    from adit.gui.panels.workspace_panel import WorkspacePanel

    from pathlib import PurePosixPath, PureWindowsPath   # Path("/tmp/x") would turn into \\tmp\\x on Windows

    monkeypatch.setattr(os, "name", "nt")
    assert WorkspacePanel._quoted(PureWindowsPath("C:\\My Runs\\x")) == '"C:\\My Runs\\x"'
    monkeypatch.setattr(os, "name", "posix")
    assert WorkspacePanel._quoted(PurePosixPath("/tmp/my runs/x")) == "'/tmp/my runs/x'"


# ---- M4: the editor notices changes on disk ------------------------------------

def _touch_later(path: Path) -> None:
    st = path.stat()
    os.utime(path, ns=(st.st_atime_ns + 2_000_000_000, st.st_mtime_ns + 2_000_000_000))


@pytest.mark.skipif(os.name == "nt", reason="the shell on Windows is different")
def test_saving_over_a_file_changed_on_disk_asks_first(app, tmp_path, monkeypatch):
    from adit.gui.panels.workspace_panel import WorkspacePanel

    panel = WorkspacePanel(root=tmp_path)
    f = tmp_path / "submit.sh"; f.write_text("one\n", encoding="utf-8")
    asked = []
    try:
        assert panel.open_file(f) == ""
        panel.editor.appendPlainText("edit")
        f.write_text("two\n", encoding="utf-8"); _touch_later(f)       # another program rewrote it
        assert panel.editor.changed_on_disk()
        monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: asked.append(a[2]) or QMessageBox.StandardButton.No))
        assert panel.save() != "" and f.read_text(encoding="utf-8") == "two\n" and len(asked) == 1
        assert "submit.sh" in asked[0]
        monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: asked.append(a[2]) or QMessageBox.StandardButton.Yes))
        assert panel.save() == "" and "edit" in f.read_text(encoding="utf-8") and len(asked) == 2
        assert not panel.editor.changed_on_disk()
        panel.editor.appendPlainText("more")
        assert panel.save() == "" and len(asked) == 2                   # no change on disk: no question
    finally:
        panel.close_session()


@pytest.mark.skipif(os.name == "nt", reason="the shell on Windows is different")
def test_the_open_file_is_reloaded_after_it_was_regenerated(app, tmp_path, monkeypatch):
    from adit.gui.panels.workspace_panel import WorkspacePanel

    panel = WorkspacePanel(root=tmp_path)
    f = tmp_path / "submit.sh"; f.write_text("one\n", encoding="utf-8")
    msgs = []; panel.status.connect(msgs.append)
    try:
        assert panel.open_file(f) == ""
        assert panel.files_written([tmp_path / "other.txt"]) == ""
        f.write_text("regen\n", encoding="utf-8"); _touch_later(f)
        assert panel.files_written([f]) == "reloaded"
        assert panel.editor.toPlainText() == "regen\n" and not panel.editor.dirty and msgs
        panel.editor.appendPlainText("mine")
        f.write_text("regen2\n", encoding="utf-8"); _touch_later(f)
        monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.No))
        assert panel.files_written([str(f)]) == "kept" and "mine" in panel.editor.toPlainText() and panel.editor.dirty
        monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))
        assert panel.files_written([f]) == "reloaded" and panel.editor.toPlainText() == "regen2\n" and not panel.editor.dirty
    finally:
        panel.close_session()


def test_generate_reloads_the_open_file_and_cds_the_terminal(app, quiet, sk_root, tmp_path):
    from adit.gui.panels.workspace_panel import WorkspacePanel

    win = make_window(sk_root, tmp_path)
    win.method.sk_set.setCurrentText("fake-1-0"); win.refresh_preview()
    win.generate()
    out = tmp_path / "out"
    assert win.workspace.open_file(out / "submit.sh") == ""
    text = win.workspace.editor.toPlainText()
    (out / "submit.sh").write_text("stale\n", encoding="utf-8"); _touch_later(out / "submit.sh")
    sent = []
    win.workspace.terminal.current.session = FakeSession()
    win.workspace.terminal.current.session.write = sent.append
    win.generate()                                                      # overwrite confirmed by the quiet fixture
    assert win.workspace.editor.toPlainText() == text and not win.workspace.editor.dirty
    assert sent == [f"cd {WorkspacePanel._quoted(out)}\r"]                 # double quotes on Windows
    win.workspace.terminal.current.session = None
    win.close()


# ---- M5 / M6: English mode -----------------------------------------------------

def _texts(root) -> list[tuple[str, str]]:
    from PySide6.QtWidgets import QAbstractButton, QGroupBox, QLabel, QLineEdit, QPlainTextEdit, QTabBar, QTabWidget, QToolButton, QWidget

    out = []
    for w in [root, *root.findChildren(QWidget)]:
        if isinstance(w, QGroupBox):
            out.append(("title", w.title()))
        if isinstance(w, (QLabel, QAbstractButton)):
            out.append(("text", w.text()))
        if isinstance(w, QComboBox):
            out += [("item", w.itemText(i)) for i in range(w.count())]
        if isinstance(w, (QLineEdit, QPlainTextEdit)):
            out.append(("placeholder", w.placeholderText()))
        if isinstance(w, (QTabWidget, QTabBar)):
            out += [("tab", w.tabText(i)) for i in range(w.count())]
        if isinstance(w, QToolButton) and w.defaultAction() is not None:
            out += [("action", w.defaultAction().text()), ("iconText", w.defaultAction().iconText())]
        out.append(("tooltip", w.toolTip()))
    return out


PROPER_NOUNS = {"日本語"}       # the name of the language on its own button


def _japanese_leaks(root) -> list[tuple[str, str]]:
    return sorted({(k, t) for k, t in _texts(root) if t and JA.search(t) and t not in PROPER_NOUNS})


def test_the_analysis_combo_boxes_use_their_english_texts(app, english):
    from adit.gui.panels.analysis_panel import AnalysisPanel

    panel = AnalysisPanel()
    items = [c.itemText(i) for c in panel.findChildren(QComboBox) for i in range(c.count())]
    assert items and not [t for t in items if JA.search(t)]
    assert panel.msd_axes.itemText(0) == "3D (xyz, MSD = 6Dt)" and panel.plot_grid.itemText(0) == "(default)"
    assert panel.btn_report.text() == "Build a report…"


def test_a_translated_main_window_shows_no_japanese(app, quiet, english, sk_root, tmp_path):
    from adit.gui.i18n import translate_widgets

    win = make_window(sk_root, tmp_path); translate_widgets(win)
    assert [b.text() for b in win.mode_bar.buttons] == ["Structure", "Calculation settings", "Analysis", "Workspace"]
    assert win.structure_view.rotation.itemText(0) == "Along x+y+z"
    assert any(it.title == "Workspace" for it in win.palette_items())
    assert _japanese_leaks(win) == []
    # texts written by slots after the translation: molecule -> periodic -> molecule
    win.structure.set_source("bulk"); win.structure._rebuild(); win._on_context(); win._timer.stop(); win.refresh_preview()
    assert win.method.dftb.solv_file.placeholderText() == "not available for periodic systems (molecules only)"
    win.structure.set_source("preset"); win.structure._rebuild(); win._on_context(); win._timer.stop(); win.refresh_preview()
    assert win.method.dftb.solv_file.placeholderText() == "Empty = no solvent (a file like param_gbsa_<solvent>.txt)"
    assert win.kpoints.info.text() == "Not used for molecules"
    for i in range(win.method.code.count()):
        win.method.code.setCurrentIndex(i); win._timer.stop(); win.refresh_preview()
    win.set_mode(win.MODE_ANALYSIS); win.analysis.more.set_expanded(True)
    assert _japanese_leaks(win) == []
    win.close()


def test_file_dialog_titles_and_filters_are_language_pairs():
    root = Path(__file__).resolve().parents[1] / "src" / "adit" / "gui" / "panels"
    for name in ("gromacs_panel.py", "lammps_panel.py"):
        text = (root / name).read_text(encoding="utf-8")
        assert "すべて (*)\"" not in text.replace("すべて (*)\", \"", ""), name   # every filter string is wrapped in L(ja, en)
        for line in text.splitlines():
            if "file_row(" in line and "def file_row" not in line:
                assert "L(" in line, line


def test_the_sk_set_info_is_english_in_english_mode(app, english, sk_root):
    from adit.gui.panels.method_panel import DftbMethodPanel

    panel = DftbMethodPanel(str(sk_root))
    panel.sk_set.setCurrentText("fake-1-0")
    assert "Hubbard" not in panel.sk_info.text() or "derivatives" in panel.sk_info.text()
    assert not JA.search(panel.sk_info.text()), panel.sk_info.text()


def test_i18n_set_language_drives_the_message_language():
    from adit import lang
    from adit.gui import i18n

    i18n.set_language("en")
    try:
        assert lang.LANGUAGE == "en" and i18n.LANGUAGE == "en"
    finally:
        i18n.set_language("ja")
    assert lang.LANGUAGE == "ja" and i18n.LANGUAGE == "ja"


# ---- M7 / L1: the terminal widget ---------------------------------------------

@pytest.mark.skipif(os.name == "nt", reason="the shell on Windows is different")
def test_copying_a_selection_after_the_terminal_shrank_does_not_fail(app, tmp_path):
    from adit.gui.terminal import TerminalWidget

    term = TerminalWidget(cwd=tmp_path, command=SHELL)
    real = term.session
    try:
        term.session = FakeSession(rows=6, cols=10)
        term._sel_from, term._sel_to = (30, 5), (40, 100)              # made when the terminal was tall and wide
        assert term.selected_text() == "a" * 5
        term._sel_from, term._sel_to = (2, 3), (3, 50)
        assert term.selected_text() == "a" * 7 + "\n" + "a" * 10
        term.copy()
        assert app.clipboard().text() == "a" * 7 + "\n" + "a" * 10
    finally:
        term.session = real
        term.close_session()


@pytest.mark.skipif(os.name == "nt", reason="the shell on Windows is different")
def test_ctrl_slash_reaches_the_shell_and_ctrl_f1_the_window(app, tmp_path):
    from PySide6.QtGui import QAction
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QMainWindow

    from adit.gui.terminal import TerminalWidget

    win = QMainWindow()
    term = TerminalWidget(cwd=tmp_path, command=SHELL)
    win.setCentralWidget(term)
    fired = []
    for name in ("Ctrl+/", "Ctrl+F1", "Ctrl+K"):
        act = QAction(name, win); act.setShortcut(name); act.triggered.connect(lambda _c=False, n=name: fired.append(n))
        win.addAction(act)
    sent = []
    term.send = sent.append
    try:
        win.show(); term.setFocus(); app.processEvents()
        QTest.keyClick(term, Qt.Key.Key_Slash, Qt.KeyboardModifier.ControlModifier); app.processEvents()
        assert sent == ["\x1f"] and fired == []                          # readline's undo
        QTest.keyClick(term, Qt.Key.Key_K, Qt.KeyboardModifier.ControlModifier); app.processEvents()
        assert sent == ["\x1f", "\x0b"] and fired == []
        QTest.keyClick(term, Qt.Key.Key_F1, Qt.KeyboardModifier.ControlModifier); app.processEvents()
        assert fired == ["Ctrl+F1"] and len(sent) == 2                   # function keys stay with the window
    finally:
        term.close_session()
        win.close()


# ---- M8: one analysis at a time -------------------------------------------------

def test_a_second_analysis_cannot_start_while_one_runs(app, tmp_path, monkeypatch):
    from adit.gui.panels.analysis_panel import AnalysisPanel

    panel = AnalysisPanel()
    panel.run_dir.setText(str(tmp_path))
    started = []

    def fake_start(fn, *args, **kwargs):
        started.append(fn); panel.job._running = True
        return 1

    monkeypatch.setattr(panel.job, "start", fake_start)
    assert panel.run_async() is True
    assert not panel.btn_run.isEnabled() and not panel.btn_export.isEnabled() and not panel.empty.button.isEnabled()
    assert panel.run_async() is False and panel.run_async(export=True) is False and len(started) == 1
    panel.cancel()
    assert panel.btn_run.isEnabled() and panel.empty.button.isEnabled()
    assert panel.run_async() is True and len(started) == 2
    panel.cancel()


# ---- L2: the isosurface level ---------------------------------------------------

def test_isosurface_level_commas(app):
    from adit.analysis.isosurface import IsosurfaceError
    from adit.gui.isosurface_panel import IsosurfacePanel

    panel = IsosurfacePanel()
    panel.level.setText("0,05")
    assert panel.parse_level() == pytest.approx(0.05) and "0.05" in panel._level_note
    panel.level.setText("0.05")
    assert panel.parse_level() == pytest.approx(0.05) and panel._level_note == ""
    for bad in ("1,000.5", "0.1,0.2", "1,2,3"):
        panel.level.setText(bad)
        with pytest.raises(IsosurfaceError):
            panel.parse_level()


# ---- L3: the k-point shift -------------------------------------------------------

def test_a_loaded_kpoint_shift_is_shown_as_it_is(app):
    from adit.gui.panels.kpoints_panel import KPointsPanel

    panel = KPointsPanel()
    panel.set_kpoints(KPoints(mode="mesh", mesh=(2, 2, 2), shift=(0.5, 0.5, 0.0)))
    assert panel.shift.currentText() == "0.5 0.5 0" and panel.kpoints().shift == (0.5, 0.5, 0.0)
    panel.shift.setCurrentText("0")
    assert panel.kpoints().shift == (0.0, 0.0, 0.0)
    panel.shift.setCurrentText("0.5 0.5 0")
    assert panel.kpoints().shift == (0.5, 0.5, 0.0)
    panel.set_kpoints(KPoints(mode="mesh", mesh=(2, 2, 2), shift=(0.25, 0.25, 0.25)))
    assert panel.shift.currentText() == "0.25" and panel.kpoints().shift == (0.25, 0.25, 0.25)
    assert panel.shift.findText("0.5 0.5 0") < 0                        # the previous extra item is gone
    panel.set_kpoints(KPoints(mode="mesh", mesh=(2, 2, 2), shift=(0.5, 0.5, 0.5)))
    assert panel.shift.currentText() == "0.5" and panel.shift.count() == 2 and panel.kpoints().shift == (0.5, 0.5, 0.5)


# ---- L5: the README action --------------------------------------------------------

def test_readme_falls_back_to_the_project_url(monkeypatch):
    import importlib.metadata as md

    from adit.gui.main_window import project_url, readme_location

    assert readme_location().endswith("README.md") or readme_location().startswith("http")
    assert project_url().startswith("https://github.com/") and project_url() != "https://github.com/"
    monkeypatch.setattr(md, "metadata", lambda name: (_ for _ in ()).throw(md.PackageNotFoundError(name)))
    assert project_url() == ""


# ---- L6 / L9: widgets that are freed ---------------------------------------------------

@pytest.mark.skipif(os.name == "nt", reason="the shell on Windows is different")
def test_a_closed_terminal_tab_is_deleted(app, tmp_path, monkeypatch):
    import shiboken6

    from adit.gui.terminal_pane import TerminalTabs

    monkeypatch.setenv("ADIT_SHELL", SHELL[0])
    tabs = TerminalTabs(cwd=tmp_path)
    try:
        view = tabs.add_tab()
        assert tabs.tabs.count() == 2
        tabs.close_tab(1); _flush_deletes()
        assert not shiboken6.isValid(view) and tabs.tabs.count() == 1
    finally:
        tabs.close_session()


def test_the_command_palette_is_deleted_when_closed(app, quiet, sk_root, tmp_path):
    import shiboken6

    win = make_window(sk_root, tmp_path)
    dlg = win.open_palette(); dlg.reject(); _flush_deletes()
    assert not shiboken6.isValid(dlg)
    win.close()


# ---- L7 / testsci M5: Open in ASE GUI ----------------------------------------------------

def _two_atoms() -> Structure:
    return Structure(source="file", source_ref="x", atoms=AtomsData(symbols=["Ar", "Ar"], positions=[(0, 0, 0), (3, 0, 0)]))


def test_the_view_panel_keeps_one_temporary_directory_and_removes_it(app):
    from adit.gui.panels.structure_view_panel import StructureViewPanel

    panel = StructureViewPanel()
    assert panel._tmp is None
    d = panel.temp_dir()
    assert d.is_dir() and panel.temp_dir() == d
    panel.cleanup()
    assert not d.exists() and panel._tmp is None


def test_ase_gui_is_not_started_through_the_frozen_executable(app, monkeypatch):
    from adit.gui.panels import structure_view_panel as svp

    calls = []
    monkeypatch.setattr(subprocess, "Popen", lambda cmd, *a, **k: calls.append(list(cmd)))
    assert svp.ase_gui_command() == [sys.executable, "-m", "ase", "gui"]
    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(shutil, "which", lambda name: None)
    assert svp.ase_gui_command() is None
    panel = svp.StructureViewPanel()
    panel.set_structure(_two_atoms())
    assert not panel.asegui.isEnabled() and "ASE" in panel.asegui.toolTip()
    panel._open_ase_gui()
    assert calls == []
    monkeypatch.setattr(shutil, "which", lambda name: "/opt/bin/ase" if name == "ase" else None)
    panel = svp.StructureViewPanel()
    panel.set_structure(_two_atoms())
    assert panel.asegui.isEnabled()
    panel._open_ase_gui()
    assert len(calls) == 1 and calls[0][:2] == ["/opt/bin/ase", "gui"] and calls[0][2].endswith("structure.xyz")
    assert Path(calls[0][2]).is_file()
    panel.cleanup()
    assert not Path(calls[0][2]).exists()


def test_the_structure_panel_button_follows_the_same_rule(app, monkeypatch):
    from adit.gui.panels.structure_panel import StructurePanel

    monkeypatch.setattr(sys, "frozen", True, raising=False)
    monkeypatch.setattr(shutil, "which", lambda name: None)
    panel = StructurePanel()
    assert not panel.asegui.isEnabled() and panel._ase_cmd is None
    calls = []
    monkeypatch.setattr(subprocess, "Popen", lambda cmd, *a, **k: calls.append(list(cmd)))
    panel._open_ase_gui()
    assert calls == []


def test_closing_the_window_removes_the_temporary_directories(app, quiet, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    a, b = win.structure_view.temp_dir(), win.structure.temp_dir()
    assert a.is_dir() and b.is_dir() and a != b
    assert win.close()
    assert not a.exists() and not b.exists()


# ---- S1 / S2: errors and progress ---------------------------------------------------------

def test_every_adit_error_in_the_preview_reaches_the_error_list(app, quiet, sk_root, tmp_path, monkeypatch):
    from adit.structure import StructureError

    win = make_window(sk_root, tmp_path)
    assert win.btn_generate.isEnabled()

    def boom():
        raise StructureError("no such atom")

    monkeypatch.setattr(win, "build", boom)
    win.refresh_preview()
    assert not win.btn_generate.isEnabled() and win.error_badge_action.isVisible()
    assert "no such atom" in win.preview.status.text() and win._errors == ["no such atom"]

    def disk():
        raise OSError("disk gone")

    monkeypatch.setattr(win, "build", disk)
    win.refresh_preview()
    assert not win.btn_generate.isEnabled() and "disk gone" in win.preview.status.text()
    monkeypatch.undo()
    win.refresh_preview()
    assert win.btn_generate.isEnabled()
    win.close()


def test_progress_from_a_superseded_job_is_dropped_in_the_worker(app):
    from adit import progress as reports
    from adit.gui.progress import Job

    job = Job()
    raw, shown = [], []
    job._report.connect(lambda *a: raw.append(a))
    job.progress.connect(lambda *a: shown.append(a))

    def work():
        reports.report(1, 2, "step")
        return "done"

    job._token = 5
    job._work(4, work, (), {})                                          # an older run: its reports are not even emitted
    app.processEvents()
    assert raw == [] and shown == []
    job._work(5, work, (), {})
    app.processEvents()
    assert raw == [(5, 1, 2, "step")] and shown == [(1, 2, "step")]
    assert reports.report(1, 1) is None                                  # the reporter was cleared afterwards


# ---- unknown keys in an opened spec.json -----------------------------------------------------------

def test_opening_a_spec_with_unknown_keys_warns_and_still_loads(app, quiet, sk_root, tmp_path, monkeypatch):
    proj = tmp_path / "proj"; proj.mkdir()
    data = json.loads(water_spec().model_dump_json())
    data["task"]["max_stepz"] = 5
    data["bogus"] = 1
    (proj / "spec.json").write_text(json.dumps(data), encoding="utf-8")
    win = make_window(sk_root, tmp_path)
    warned = []
    monkeypatch.setattr(QFileDialog, "getOpenFileName", staticmethod(lambda *a, **k: (str(proj / "spec.json"), "")))
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda *a, **k: warned.append(a[2])))
    win.open_spec()
    assert len(warned) == 1 and "task.max_stepz" in warned[0] and "bogus" in warned[0]
    assert win.runtime.outdir.text() == str(proj)
    assert win.current_spec().task.max_steps == 100                     # the file's valid values were applied
    warned.clear()
    (proj / "spec.json").write_text(json.dumps(json.loads(water_spec().model_dump_json())), encoding="utf-8")
    win.open_spec()
    assert warned == []
    win.close()
