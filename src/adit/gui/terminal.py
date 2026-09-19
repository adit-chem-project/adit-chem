"""Terminal widget: draws the screen of a shell session and sends keystrokes to it."""

from __future__ import annotations

from PySide6.QtCore import QEvent, QRect, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QFontDatabase, QFontMetricsF, QKeyEvent, QPainter
from PySide6.QtWidgets import QApplication, QWidget

from adit.gui.terminal_backend import ShellSession, TerminalError, available
from adit.lang import L

POLL_MS = 30
PADDING = 6

# xterm's 8 colours plus the default; pyte gives a name or 6 hex digits
NAMED = {"black": "#2e3436", "red": "#cc0000", "green": "#4e9a06", "brown": "#c4a000", "yellow": "#c4a000",
         "blue": "#3465a4", "magenta": "#75507b", "cyan": "#06989a", "white": "#d3d7cf"}
BRIGHT = {"black": "#555753", "red": "#ef2929", "green": "#8ae234", "brown": "#fce94f", "yellow": "#fce94f",
          "blue": "#729fcf", "magenta": "#ad7fa8", "cyan": "#34e2e2", "white": "#eeeeec"}


def _color(name: str, default: str, bold: bool = False) -> QColor:
    if not name or name == "default":
        return QColor(default)
    table = BRIGHT if bold else NAMED
    if name in table:
        return QColor(table[name])
    if len(name) == 6:
        return QColor("#" + name)
    return QColor(default)


class TerminalWidget(QWidget):
    """A shell in a widget. Left click focuses it; keys go straight to the shell."""

    finished = Signal()
    scrolled = Signal()

    def __init__(self, cwd=None, command=None, dark: bool = False, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, True)
        self.dark = dark
        self.session: ShellSession | None = None
        self.error = ""
        self._sel_from: tuple[int, int] | None = None   # (row, column)
        self._sel_to: tuple[int, int] | None = None
        self._selecting = False
        self.setMouseTracking(True)
        font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        font.setPointSizeF(max(9.0, font.pointSizeF()))
        font.setStyleHint(QFont.StyleHint.Monospace)
        self.setFont(font)
        self._measure()
        ok, missing = available()
        if not ok:
            self.error = L(f"ターミナルを開けません ({missing} が入っていません)",
                           f"cannot open a terminal ({missing} is not installed)")
        else:
            try:
                self.session = ShellSession(cwd=cwd, command=command, rows=self.rows(), cols=self.cols())
            except TerminalError as ex:
                self.error = L(f"シェルを起動できません: {ex}", f"cannot start a shell: {ex}")
        self._timer = QTimer(self)
        self._timer.timeout.connect(self._poll)
        self._timer.start(POLL_MS)

    def _measure(self) -> None:
        self._metrics = QFontMetricsF(self.font())
        self._cw = max(1.0, self._metrics.horizontalAdvance("M"))
        self._ch = max(1.0, self._metrics.height())

    def set_font_size(self, points: float) -> None:
        """Change the font size (Ctrl+= / Ctrl+- / Ctrl+0)."""
        font = self.font()
        font.setPointSizeF(max(6.0, min(36.0, points)))
        self.setFont(font)
        self._measure()
        if self.session is not None:
            self.session.resize(self.rows(), self.cols())
        self.update()

    def zoom(self, step: float) -> None:
        self.set_font_size(self.font().pointSizeF() + step)

    # ---- size in characters ----
    def cols(self) -> int:
        return max(20, int((self.width() - 2 * PADDING) / self._cw))

    def rows(self) -> int:
        return max(4, int((self.height() - 2 * PADDING) / self._ch))

    def resizeEvent(self, event) -> None:  # noqa: N802 (Qt)
        super().resizeEvent(event)
        if self.session is not None:
            self.session.resize(self.rows(), self.cols())
            self.update()

    # ---- running ----
    def _poll(self) -> None:
        if self.session is None:
            return
        if self.session.drain():
            self.update()
        if not self.session.alive:
            self._timer.stop()
            self.finished.emit()
            self.update()

    def send(self, text: str) -> None:
        """Type text into the shell (used by the buttons that paste a command)."""
        if self.session is not None:
            self.session.write(text)

    def close_session(self) -> None:
        self._timer.stop()
        if self.session is not None:
            self.session.close()
            self.session = None

    # ---- drawing ----
    def paintEvent(self, _event) -> None:  # noqa: N802 (Qt)
        p = QPainter(self)
        back = QColor("#1e1e1e") if self.dark else QColor("#ffffff")
        fore = "#d7d7d7" if self.dark else "#1e1e1e"
        p.fillRect(self.rect(), back)
        p.setFont(self.font())
        if self.session is None:
            p.setPen(QColor("#b00020"))
            p.drawText(self.rect().adjusted(PADDING, PADDING, -PADDING, -PADDING),
                       int(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft | Qt.TextFlag.TextWordWrap),
                       self.error)
            p.end()
            return
        ascent = self._metrics.ascent()
        for y, row in enumerate(self.session.lines()):
            x = 0
            while x < len(row):
                cell = row[x]
                if cell.text == "":            # second cell of a wide character
                    x += 1
                    continue
                run, start = cell.text, x
                x += 1
                while x < len(row) and row[x].fg == cell.fg and row[x].bg == cell.bg and row[x].bold == cell.bold \
                        and row[x].reverse == cell.reverse:
                    run += row[x].text         # keep the empty second cell of a wide character
                    x += 1
                fg = _color(cell.fg, fore, cell.bold)
                bg = _color(cell.bg, back.name())
                if cell.reverse:
                    fg, bg = bg, fg
                left = PADDING + start * self._cw
                top = PADDING + y * self._ch
                if bg != back:
                    p.fillRect(QRect(int(left), int(top), int(len(run) * self._cw) + 1, int(self._ch) + 1), bg)
                if run.strip():
                    p.setPen(fg)
                    f = p.font()
                    f.setBold(cell.bold)
                    p.setFont(f)
                    p.drawText(int(left), int(top + ascent), run)
        if self._sel_from is not None and self._sel_to is not None and self._sel_from != self._sel_to:
            (r1, c1), (r2, c2) = sorted([self._sel_from, self._sel_to])
            for row in range(r1, r2 + 1):
                start = c1 if row == r1 else 0
                end = c2 if row == r2 else self.session.cols - 1
                rect = QRect(int(PADDING + start * self._cw), int(PADDING + row * self._ch),
                             int((end - start + 1) * self._cw), int(self._ch) + 1)
                p.fillRect(rect, QColor(60, 120, 220, 70))
        cy, cx = self.session.cursor
        if self.hasFocus():
            p.fillRect(QRect(int(PADDING + cx * self._cw), int(PADDING + cy * self._ch),
                             max(2, int(self._cw)), int(self._ch)), QColor(fore))
        p.end()

    # ---- keyboard ----
    KEYS = {Qt.Key.Key_Return: "\r", Qt.Key.Key_Enter: "\r", Qt.Key.Key_Backspace: "\x7f",
            Qt.Key.Key_Tab: "\t", Qt.Key.Key_Backtab: "\x1b[Z", Qt.Key.Key_Escape: "\x1b", Qt.Key.Key_Delete: "\x1b[3~",
            Qt.Key.Key_Up: "\x1b[A", Qt.Key.Key_Down: "\x1b[B", Qt.Key.Key_Right: "\x1b[C", Qt.Key.Key_Left: "\x1b[D",
            Qt.Key.Key_Home: "\x1b[H", Qt.Key.Key_End: "\x1b[F",
            Qt.Key.Key_PageUp: "\x1b[5~", Qt.Key.Key_PageDown: "\x1b[6~"}

    @staticmethod
    def _wants_key(event: QKeyEvent) -> bool:
        key, mods = event.key(), event.modifiers()
        if key in (Qt.Key.Key_Tab, Qt.Key.Key_Backtab):
            return True
        if not mods & Qt.KeyboardModifier.ControlModifier or not Qt.Key.Key_A <= key <= Qt.Key.Key_Z:
            return False
        return not (mods & Qt.KeyboardModifier.ShiftModifier and key in (Qt.Key.Key_C, Qt.Key.Key_V))

    def event(self, event) -> bool:  # noqa: N802 (Qt)
        # Ctrl+letter and Tab belong to the shell; without this the window's
        # QActions (Ctrl+Z, Ctrl+D, ...) and the focus chain would take them first.
        if self.session is not None and event.type() in (QEvent.Type.ShortcutOverride, QEvent.Type.KeyPress) \
                and self._wants_key(event):
            if event.type() == QEvent.Type.KeyPress:
                self.keyPressEvent(event)
            event.accept()
            return True
        return super().event(event)

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802 (Qt)
        if self.session is None:
            return
        mods, key = event.modifiers(), event.key()
        if mods & Qt.KeyboardModifier.ControlModifier and mods & Qt.KeyboardModifier.ShiftModifier and key == Qt.Key.Key_C:
            self.copy()                                               # the selection, or the whole screen
            return
        if mods & Qt.KeyboardModifier.ControlModifier and mods & Qt.KeyboardModifier.ShiftModifier and key == Qt.Key.Key_V:
            self.paste()
            return
        if mods & Qt.KeyboardModifier.ControlModifier and key in (Qt.Key.Key_Plus, Qt.Key.Key_Equal):
            self.zoom(1.0)
            return
        if mods & Qt.KeyboardModifier.ControlModifier and key == Qt.Key.Key_Minus:
            self.zoom(-1.0)
            return
        if mods & Qt.KeyboardModifier.ControlModifier and key == Qt.Key.Key_0:
            self.set_font_size(QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont).pointSizeF())
            return
        if key in self.KEYS:
            self.send(self.KEYS[key])
            return
        if mods & Qt.KeyboardModifier.ControlModifier and Qt.Key.Key_A <= key <= Qt.Key.Key_Z:
            self.send(chr(key - Qt.Key.Key_A + 1))                    # Ctrl+C sends \x03
            return
        if event.text():
            self.send(event.text())

    # ---- mouse ----
    def _cell_at(self, pos) -> tuple[int, int]:
        row = int((pos.y() - PADDING) / self._ch)
        col = int((pos.x() - PADDING) / self._cw)
        rows = self.session.rows if self.session else 1
        cols = self.session.cols if self.session else 1
        return max(0, min(rows - 1, row)), max(0, min(cols - 1, col))

    def _button_code(self, button) -> int:
        return {Qt.MouseButton.LeftButton: 0, Qt.MouseButton.MiddleButton: 1, Qt.MouseButton.RightButton: 2}.get(button, 0)

    def mousePressEvent(self, event) -> None:  # noqa: N802 (Qt)
        if self.session is None:
            return
        self.setFocus()
        row, col = self._cell_at(event.position())
        if self.session.mouse_wanted():
            self.send(self.session.mouse_report(self._button_code(event.button()), col, row, True))
            return
        if event.button() == Qt.MouseButton.LeftButton:
            self._sel_from = self._sel_to = (row, col)
            self._selecting = True
            self.update()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802 (Qt)
        if self.session is None:
            return
        row, col = self._cell_at(event.position())
        if self.session.mouse_wanted():
            if self.session.mouse_motion_wanted() and event.buttons():
                self.send(self.session.mouse_report(32 + self._button_code(event.buttons()), col, row, True))
            return
        if self._selecting:
            self._sel_to = (row, col)
            self.update()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802 (Qt)
        if self.session is None:
            return
        row, col = self._cell_at(event.position())
        if self.session.mouse_wanted():
            self.send(self.session.mouse_report(self._button_code(event.button()), col, row, False))
            return
        self._selecting = False

    def selected_text(self) -> str:
        """The selected text, or "" when nothing is selected."""
        if self.session is None or self._sel_from is None or self._sel_to is None:
            return ""
        (r1, c1), (r2, c2) = sorted([self._sel_from, self._sel_to])
        lines = self.session.lines()
        out = []
        for row in range(r1, r2 + 1):
            start = c1 if row == r1 else 0
            end = c2 if row == r2 else len(lines[row]) - 1
            out.append("".join(cell.text for cell in lines[row][start:end + 1]).rstrip())
        return "\n".join(out)

    def copy(self) -> None:
        text = self.selected_text() or (self.session.text() if self.session else "")
        if text:
            QApplication.clipboard().setText(text)

    def paste(self) -> None:
        # A pasted line break must reach the shell as a single Enter (CRLF would give two).
        self.send(QApplication.clipboard().text().replace("\r\n", "\n").replace("\n", "\r"))

    def contextMenuEvent(self, event) -> None:  # noqa: N802 (Qt)
        if self.session is not None and self.session.mouse_wanted():
            return                                  # the program handles right clicks itself
        from PySide6.QtWidgets import QMenu

        menu = QMenu(self)
        act_copy = menu.addAction(L("コピー", "Copy"))
        act_copy.setEnabled(bool(self.selected_text()))
        act_all = menu.addAction(L("画面全体をコピー", "Copy the whole screen"))
        act_paste = menu.addAction(L("貼り付け", "Paste"))
        chosen = menu.exec(event.globalPos())
        if chosen is act_copy:
            QApplication.clipboard().setText(self.selected_text())
        elif chosen is act_all and self.session is not None:
            QApplication.clipboard().setText(self.session.text())
        elif chosen is act_paste:
            self.paste()

    def wheelEvent(self, event) -> None:  # noqa: N802 (Qt)
        """Zoom with Ctrl, report to the program, or scroll the history."""
        if self.session is None:
            return
        steps = event.angleDelta().y()
        if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
            self.zoom(1.0 if steps > 0 else -1.0)
            return
        if self.session.mouse_wanted():
            row, col = self._cell_at(event.position())
            self.send(self.session.mouse_report(64 if steps > 0 else 65, col, row, True))
            return
        self.session.scroll_pages(-1 if steps > 0 else 1)
        self.scrolled.emit()
        self.update()

    def restart(self, cwd=None) -> None:
        """Start the shell again after it exited."""
        self.close_session()
        try:
            self.session = ShellSession(cwd=cwd, rows=self.rows(), cols=self.cols())
            self.error = ""
        except TerminalError as ex:
            self.error = L(f"シェルを起動できません: {ex}", f"cannot start a shell: {ex}")
        self._timer.start(POLL_MS)
        self.update()

    def focusInEvent(self, event) -> None:  # noqa: N802 (Qt)
        super().focusInEvent(event)
        self.update()

    def focusOutEvent(self, event) -> None:  # noqa: N802 (Qt)
        super().focusOutEvent(event)
        self.update()
