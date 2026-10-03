"""Several terminals in tabs, each with a scrollbar."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QHBoxLayout, QMessageBox, QPushButton, QScrollBar, QTabWidget, QVBoxLayout, QWidget

from adit.gui.terminal import TerminalWidget
from adit.lang import L


class TerminalView(QWidget):
    """One terminal and its scrollbar."""

    def __init__(self, cwd=None, dark: bool = False, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.terminal = TerminalWidget(cwd=cwd, dark=dark)
        self.bar = QScrollBar(Qt.Orientation.Vertical)
        self.bar.setRange(0, 0)
        self.bar.valueChanged.connect(self._on_bar)
        self.terminal.scrolled.connect(self._sync_bar)
        self._syncing = False
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(0)
        lay.addWidget(self.terminal, 1); lay.addWidget(self.bar)
        self.terminal.installEventFilter(self)

    def eventFilter(self, obj, event) -> bool:  # noqa: N802 (Qt)
        from PySide6.QtCore import QEvent

        if obj is self.terminal and event.type() in (QEvent.Type.Paint, QEvent.Type.Resize):
            self._sync_bar()
        return super().eventFilter(obj, event)

    def _sync_bar(self) -> None:
        session = self.terminal.session
        if session is None or self._syncing:
            return
        above, below = session.history_above, session.history_below
        self._syncing = True
        self.bar.setRange(0, above + below)
        self.bar.setPageStep(max(1, session.rows))
        self.bar.setValue(above)
        self._syncing = False

    def _on_bar(self, value: int) -> None:
        session = self.terminal.session
        if session is None or self._syncing:
            return
        steps = value - session.history_above
        if steps:
            session.scroll_pages(steps // max(1, int(session.rows * 0.5)) or (1 if steps > 0 else -1))
            self.terminal.update()
            self._sync_bar()

    def close_session(self) -> None:
        self.terminal.close_session()


class TerminalTabs(QWidget):
    """Terminals in tabs: one per remote host, or one per calculation."""

    tab_added = Signal(object)

    def __init__(self, cwd=None, dark: bool = False, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.cwd = Path(cwd).expanduser() if cwd else Path.home()
        self.dark = dark
        self.tabs = QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.setDocumentMode(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        self.btn_new = QPushButton("+")
        self.btn_new.setFixedWidth(28)
        self.btn_new.setToolTip(L("ターミナルをもう 1 つ開く", "Open another terminal"))
        self.btn_new.clicked.connect(lambda: self.add_tab())
        self.tabs.setCornerWidget(self.btn_new, Qt.Corner.TopRightCorner)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(0)
        lay.addWidget(self.tabs)
        self.add_tab()

    # ---- tabs ----
    def add_tab(self, cwd=None) -> TerminalView:
        view = TerminalView(cwd=cwd or self.cwd, dark=self.dark)
        index = self.tabs.addTab(view, str(self.tabs.count() + 1))
        self.tabs.setCurrentIndex(index)
        self.tab_added.emit(view)
        return view

    def close_tab(self, index: int) -> None:
        view = self.tabs.widget(index)
        if isinstance(view, TerminalView):
            session = view.terminal.session
            if session is not None and session.busy():
                ans = QMessageBox.question(self, L("計算が走っています", "A program is running"),
                                           L("このターミナルでは計算が走っています。閉じますか?",
                                             "A program is still running in this terminal. Close it?"))
                if ans != QMessageBox.StandardButton.Yes:
                    return
            view.close_session()
        self.tabs.removeTab(index)
        if view is not None:
            view.deleteLater()      # removeTab only detaches it; the screen history would stay in memory
        if self.tabs.count() == 0:
            self.add_tab()

    @property
    def current(self) -> TerminalWidget | None:
        view = self.tabs.currentWidget()
        return view.terminal if isinstance(view, TerminalView) else None

    # ---- what the workspace uses ----
    def why_blocked(self) -> str:
        # "" when a command may be typed into the current terminal's prompt; otherwise the reason, for the user.
        term = self.current
        session = term.session if term is not None else None
        if session is None or not session.alive:
            return L("シェルが動いていません (「シェルを起動し直す」を押してください)",
                     "the shell is not running (press \"Restart the shell\")")
        if session.busy():
            return L("ターミナルでプログラムが動いています (終わるのを待つか、Ctrl+C で止めてください)",
                     "a program is running in the terminal (wait for it, or stop it with Ctrl+C)")
        return ""

    def send(self, text: str) -> bool:
        if self.why_blocked():
            return False
        self.current.send(text)
        return True

    def set_dark(self, dark: bool) -> None:
        self.dark = dark
        for i in range(self.tabs.count()):
            view = self.tabs.widget(i)
            if isinstance(view, TerminalView):
                view.terminal.dark = dark
                view.terminal.update()

    def close_session(self) -> None:
        for i in range(self.tabs.count()):
            view = self.tabs.widget(i)
            if isinstance(view, TerminalView):
                view.close_session()
