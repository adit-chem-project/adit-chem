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
        assert tabs.tabs.count() == 2 and tabs.current is not None
        tabs.close_tab(1)
        assert tabs.tabs.count() == 1                     # the last tab stays
        tabs.close_tab(0)
        assert tabs.tabs.count() == 1
    finally:
        tabs.close_session()


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
