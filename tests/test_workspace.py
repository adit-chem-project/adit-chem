import os
import time
from pathlib import Path

import pytest

pytest.importorskip("pyte")
if os.name != "nt":
    pytest.importorskip("ptyprocess")

from adit.gui.terminal_backend import ShellSession, available, default_shell

SHELL = ["cmd.exe"] if os.name == "nt" else ["/bin/sh"]
ECHO = "echo adit-terminal" + ("\r\n" if os.name == "nt" else "\n")


def _wait_for(session: ShellSession, needle: str, seconds: float = 8.0) -> bool:
    deadline = time.time() + seconds
    while time.time() < deadline:
        session.drain()
        if needle in session.text():
            return True
        time.sleep(0.05)
    return False


def test_a_shell_runs_and_its_output_reaches_the_screen(tmp_path):
    assert available()[0]
    session = ShellSession(cwd=tmp_path, command=SHELL, rows=12, cols=60)
    try:
        session.write(ECHO)
        assert _wait_for(session, "adit-terminal"), session.text()
        assert session.alive
    finally:
        session.close()


def test_the_screen_keeps_the_size_it_was_given(tmp_path):
    session = ShellSession(cwd=tmp_path, command=SHELL, rows=10, cols=40)
    try:
        assert len(session.lines()) == 10 and len(session.lines()[0]) == 40
        session.resize(20, 100)
        assert (session.rows, session.cols) == (20, 100)
        assert len(session.lines()) == 20 and len(session.lines()[0]) == 100
    finally:
        session.close()


@pytest.mark.skipif(os.name == "nt", reason="cmd.exe の出力の書式が違う")
def test_wide_characters_take_two_columns(tmp_path):
    session = ShellSession(cwd=tmp_path, command=SHELL, rows=8, cols=40)
    try:
        session.write("echo 日本語\n")
        assert _wait_for(session, "日本語"), session.text()
        row = next(r for r in session.lines() if any(c.text == "日" for c in r))
        i = next(i for i, c in enumerate(row) if c.text == "日")
        assert row[i + 1].text == ""      # second cell of a wide character is empty
    finally:
        session.close()


def test_the_default_shell_is_a_real_program():
    argv = default_shell()
    assert argv and (Path(argv[0]).name or argv[0])


@pytest.mark.skipif(os.name == "nt", reason="Windows は別のシェル")
def test_the_editor_opens_saves_and_refuses_binaries(tmp_path):
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from adit.gui.panels.workspace_panel import WorkspacePanel

    QApplication.instance() or QApplication([])

    (tmp_path / "in.hsd").write_text("Driver = {}\n", encoding="utf-8")
    (tmp_path / "blob.bin").write_bytes(bytes(range(256)))
    panel = WorkspacePanel(root=tmp_path)
    try:
        assert panel.open_file(tmp_path / "in.hsd") == ""
        panel.editor.insertPlainText("# note\n")
        assert panel.editor.dirty and panel.save() == ""
        assert (tmp_path / "in.hsd").read_text(encoding="utf-8").startswith("# note")
        assert panel.open_file(tmp_path / "blob.bin") != ""      # returns the reason
    finally:
        panel.close_session()


def test_the_mouse_is_reported_when_a_program_asks_for_it(tmp_path):
    session = ShellSession(cwd=tmp_path, command=SHELL, rows=10, cols=40)
    try:
        assert not session.mouse_wanted()
        session.screen.set_mode(1000, private=True)
        assert session.mouse_wanted()
        assert session.mouse_report(0, 4, 2, True) == "\x1b[M \x25\x23".replace("\x25", chr(32 + 5)).replace("\x23", chr(32 + 3))
        session.screen.set_mode(1006, private=True)      # SGR form
        assert session.mouse_report(0, 4, 2, True) == "\x1b[<0;5;3M"
        assert session.mouse_report(0, 4, 2, False) == "\x1b[<0;5;3m"
    finally:
        session.close()


def test_the_history_can_be_scrolled_back(tmp_path):
    session = ShellSession(cwd=tmp_path, command=SHELL, rows=6, cols=40)
    try:
        for i in range(1, 13):           # more lines than the screen, whatever the shell
            session.write(f"echo line-{i}" + ("\r\n" if os.name == "nt" else "\n"))
        assert _wait_for(session, "line-12"), session.text()
        assert session.history_above > 0                  # lines scrolled away
        before = session.text()
        session.scroll_pages(-1)
        assert session.text() != before                   # scrolled back
        session.scroll_pages(1)
    finally:
        session.close()


@pytest.mark.skipif(os.name == "nt", reason="Windows は別のシェル")
def test_several_terminals_can_be_open_at_once(tmp_path):
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from adit.gui.terminal_pane import TerminalTabs

    QApplication.instance() or QApplication([])
    tabs = TerminalTabs(cwd=tmp_path)
    try:
        assert tabs.tabs.count() == 1
        tabs.add_tab()
        _wait_idle(tabs.current.session)      # a login shell still running its profile counts as busy
        assert tabs.tabs.count() == 2 and tabs.current is not None
        tabs.close_tab(1)
        assert tabs.tabs.count() == 1                     # the last tab stays
        _wait_idle(tabs.current.session)
        tabs.close_tab(0)
        assert tabs.tabs.count() == 1
    finally:
        tabs.close_session()


def _wait_idle(session, seconds: float = 8.0) -> None:
    deadline = time.time() + seconds
    while time.time() < deadline and session is not None and session.busy():
        time.sleep(0.05)


def _wait_busy(session, seconds: float = 8.0) -> bool:
    deadline = time.time() + seconds
    while time.time() < deadline:
        if session.busy():
            return True
        time.sleep(0.05)
    return False


@pytest.mark.skipif(os.name == "nt", reason="Windows は別のシェル")
def test_closing_a_tab_asks_while_a_program_is_running(tmp_path, monkeypatch):
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication, QMessageBox

    from adit.gui.terminal_pane import TerminalTabs

    QApplication.instance() or QApplication([])
    monkeypatch.setenv("ADIT_SHELL", "/bin/sh")
    tabs = TerminalTabs(cwd=tmp_path)
    try:
        session = tabs.current.session
        _wait_idle(session)
        assert not session.busy()
        session.write("sleep 30\n")
        assert _wait_busy(session)
        answers = []
        monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: answers.append(a[2]) or QMessageBox.StandardButton.No))
        tabs.close_tab(0)
        assert len(answers) == 1 and tabs.current.session is session     # answered No: still there
        monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))
        tabs.close_tab(0)
        assert not session.alive and tabs.current.session is not session
    finally:
        tabs.close_session()


@pytest.mark.skipif(os.name == "nt", reason="printf の書式が違う")
def test_bytes_that_are_not_utf8_do_not_stop_the_terminal(tmp_path):
    session = ShellSession(cwd=tmp_path, command=["/bin/sh", "-c", "printf 'before\\n\\377\\376\\nafter\\n'; sleep 1"], rows=6, cols=40)
    try:
        assert _wait_for(session, "after"), session.text()      # the reader survived the bad bytes
    finally:
        session.close()


@pytest.mark.skipif(os.name == "nt", reason="Windows は別のシェル")
def test_control_keys_and_tab_reach_the_shell_not_the_window(tmp_path):
    pytest.importorskip("PySide6")
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QAction
    from PySide6.QtTest import QTest
    from PySide6.QtWidgets import QApplication, QLineEdit, QMainWindow, QVBoxLayout, QWidget

    from adit.gui.terminal import TerminalWidget

    app = QApplication.instance() or QApplication([])
    win = QMainWindow()
    box = QWidget(); lay = QVBoxLayout(box)
    other = QLineEdit(); term = TerminalWidget(cwd=tmp_path, command=SHELL)
    lay.addWidget(other); lay.addWidget(term)
    win.setCentralWidget(box)
    fired = []
    for name in ("Ctrl+Z", "Ctrl+D", "Ctrl+Q", "Ctrl+S", "Ctrl+O", "Ctrl+G"):
        act = QAction(name, win); act.setShortcut(name); act.triggered.connect(lambda _c=False, n=name: fired.append(n))
        win.addAction(act)
    sent = []
    term.send = sent.append
    try:
        win.show(); term.setFocus(); app.processEvents()
        for key in (Qt.Key.Key_Z, Qt.Key.Key_D, Qt.Key.Key_Q, Qt.Key.Key_S, Qt.Key.Key_O, Qt.Key.Key_G):
            QTest.keyClick(term, key, Qt.KeyboardModifier.ControlModifier); app.processEvents()
        assert fired == []
        assert sent == ["\x1a", "\x04", "\x11", "\x13", "\x0f", "\x07"]
        sent.clear()
        QTest.keyClick(term, Qt.Key.Key_Tab); app.processEvents()
        QTest.keyClick(term, Qt.Key.Key_Backtab, Qt.KeyboardModifier.ShiftModifier); app.processEvents()
        assert sent == ["\t", "\x1b[Z"] and app.focusWidget() is term     # completion, not focus change
        sent.clear()
        QTest.keyClick(term, Qt.Key.Key_C, Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.ShiftModifier)
        assert sent == []                                                  # Ctrl+Shift+C stays a copy
    finally:
        term.close_session()
        win.close()


@pytest.mark.skipif(os.name == "nt", reason="Windows は別のシェル")
def test_pasted_line_breaks_become_single_enters(tmp_path):
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from adit.gui.terminal import TerminalWidget

    app = QApplication.instance() or QApplication([])
    term = TerminalWidget(cwd=tmp_path, command=SHELL)
    sent = []
    term.send = sent.append
    try:
        app.clipboard().setText("echo a\r\necho b\n")
        term.paste()
        assert sent == ["echo a\recho b\r"]
    finally:
        term.close_session()


@pytest.mark.skipif(os.name == "nt", reason="Windows は別のシェル")
def test_the_editor_keeps_crlf_and_does_not_rewrite_other_encodings(tmp_path):
    pytest.importorskip("PySide6")
    from PySide6.QtGui import QTextCursor
    from PySide6.QtWidgets import QApplication

    from adit.gui.panels.workspace_panel import WorkspacePanel

    QApplication.instance() or QApplication([])
    crlf = tmp_path / "run.bat"; crlf.write_bytes(b"echo one\r\necho two\r\n")
    latin = tmp_path / "old.log"; latin.write_bytes(b"caf\xe9\n")
    panel = WorkspacePanel(root=tmp_path)
    try:
        assert panel.open_file(crlf) == ""
        assert panel.editor.toPlainText() == "echo one\necho two\n"
        panel.editor.moveCursor(QTextCursor.MoveOperation.End); panel.editor.insertPlainText("echo three\n")
        assert panel.save() == ""
        assert crlf.read_bytes() == b"echo one\r\necho two\r\necho three\r\n"
        why = panel.open_file(latin)
        assert "UTF-8" in why and panel.editor.isReadOnly() and panel.editor.path is None
        assert panel.editor.save() != "" and latin.read_bytes() == b"caf\xe9\n"
        assert panel.open_file(crlf) == "" and not panel.editor.isReadOnly()
    finally:
        panel.close_session()


@pytest.mark.skipif(os.name == "nt", reason="Windows は別のシェル")
def test_rename_refuses_an_existing_name_and_follows_the_open_file(tmp_path, monkeypatch):
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication, QMessageBox

    from adit.gui.panels.workspace_panel import WorkspacePanel

    QApplication.instance() or QApplication([])
    warned = []
    monkeypatch.setattr(QMessageBox, "warning", staticmethod(lambda *a, **k: warned.append(a[2])))
    a = tmp_path / "a.txt"; a.write_text("AAA", encoding="utf-8")
    b = tmp_path / "b.txt"; b.write_text("BBB", encoding="utf-8")
    panel = WorkspacePanel(root=tmp_path)
    try:
        assert panel.open_file(a) == ""
        assert panel.rename(a, "b.txt") != "" and len(warned) == 1
        assert b.read_text(encoding="utf-8") == "BBB" and a.exists()
        assert panel.rename(a, "c.txt") == ""
        assert panel.editor.path == tmp_path / "c.txt" and "c.txt" in panel.file_label.text()
        panel.editor.appendPlainText("x")
        assert panel.save() == "" and (tmp_path / "c.txt").read_text(encoding="utf-8").startswith("AAA")
    finally:
        panel.close_session()


@pytest.mark.skipif(os.name == "nt", reason="Windows は別のシェル")
def test_the_restart_button_follows_the_shell_of_every_tab(tmp_path, monkeypatch):
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from adit.gui.panels.workspace_panel import WorkspacePanel

    app = QApplication.instance() or QApplication([])
    monkeypatch.setenv("ADIT_SHELL", "/bin/sh")
    panel = WorkspacePanel(root=tmp_path)
    try:
        view = panel.terminal.add_tab()
        (tmp_path / "sub").mkdir()
        panel.set_root(tmp_path / "sub")
        assert panel.terminal.cwd == tmp_path / "sub"           # the next tab opens where the tree is
        view.terminal.session.write("exit\n")
        deadline = time.time() + 8
        while time.time() < deadline and view.terminal.session is not None and view.terminal.session.alive:
            app.processEvents(); time.sleep(0.05)
        app.processEvents(); time.sleep(0.1); app.processEvents()
        assert not panel.btn_restart.isHidden()                  # second tab, and still noticed
        panel.terminal.tabs.setCurrentIndex(0)
        assert panel.btn_restart.isHidden()
        panel.terminal.tabs.setCurrentIndex(1)
        assert not panel.btn_restart.isHidden()
        panel._restart_terminal()
        assert view.terminal.session is not None and view.terminal.session.alive
    finally:
        panel.close_session()


@pytest.mark.skipif(os.name == "nt", reason="Windows は別のシェル")
def test_the_font_size_changes_the_number_of_columns(tmp_path):
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from adit.gui.terminal import TerminalWidget

    QApplication.instance() or QApplication([])
    term = TerminalWidget(cwd=tmp_path, command=SHELL)
    try:
        term.resize(800, 400)
        term.set_font_size(9)
        small = term.cols()
        term.set_font_size(18)
        assert term.cols() < small                        # bigger font, fewer columns
        assert term.session is not None and term.session.cols == term.cols()
    finally:
        term.close_session()


@pytest.mark.skipif(os.name == "nt", reason="Windows は別のシェル")
def test_the_editor_and_terminal_can_be_rearranged_and_folded(tmp_path):
    pytest.importorskip("PySide6")
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QApplication

    from adit.gui.panels.workspace_panel import WorkspacePanel

    QApplication.instance() or QApplication([])
    panel = WorkspacePanel(root=tmp_path)
    try:
        assert panel.right.orientation() == Qt.Orientation.Vertical
        assert panel.right.widget(0) is panel.edit_box
        panel.layout_choice.setCurrentIndex(panel.layout_choice.findData("h_rev"))
        assert panel.right.orientation() == Qt.Orientation.Horizontal
        assert panel.right.widget(0) is panel.term_box          # swapped
        panel.btn_fold_editor.setChecked(True)                   # no window shown, so use isHidden
        assert panel.editor.isHidden() and not panel.terminal.isHidden()
        panel.btn_fold_terminal.setChecked(True)                 # never fold both
        assert not (panel.editor.isHidden() and panel.terminal.isHidden())
    finally:
        panel.close_session()


def test_the_tree_marks_the_kind_of_each_file(tmp_path):
    pytest.importorskip("PySide6")
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QApplication

    from adit.gui.panels.workspace_panel import FileIcons

    QApplication.instance() or QApplication([])
    (tmp_path / "dftb_in.hsd").write_text("x", encoding="utf-8")
    (tmp_path / "notes.zzz").write_text("x", encoding="utf-8")
    model = FileIcons()
    model.setRootPath(str(tmp_path))
    icon = model.data(model.index(str(tmp_path / "dftb_in.hsd")), Qt.ItemDataRole.DecorationRole)
    assert icon is not None and not icon.isNull()
    assert FileIcons.KIND[".hsd"][0] == "in"                     # marked as an input file
