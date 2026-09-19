"""Shell session behind the terminal widget: a pseudo terminal plus a screen model."""
# No Qt here, so the behaviour can be tested without a display. POSIX uses ptyprocess,
# Windows uses pywinpty (ConPTY); pyte turns the byte stream into a screen of characters.

from __future__ import annotations

import os
import queue
import shutil
import sys
import threading
from dataclasses import dataclass

DEFAULT_ROWS, DEFAULT_COLS = 24, 80
READ_CHUNK = 65536
SCROLLBACK = 5000


def default_shell() -> list[str]:
    """The shell to start. ADIT_SHELL overrides it (a command line, for example "bash --norc")."""
    override = os.environ.get("ADIT_SHELL", "").strip()
    if override:
        import shlex

        return shlex.split(override, posix=os.name != "nt")
    if os.name == "nt":
        for name in ("pwsh.exe", "powershell.exe", "cmd.exe"):
            found = shutil.which(name)
            if found:
                return [found]
        return ["cmd.exe"]
    shell = os.environ.get("SHELL") or shutil.which("bash") or "/bin/sh"
    return [shell, "-l"]


class TerminalError(RuntimeError):
    pass


@dataclass
class Cell:
    text: str          # "" is the second cell of the wide character before it
    fg: str
    bg: str
    bold: bool
    reverse: bool


class ShellSession:
    """One running shell. Feed it keystrokes, read the screen."""

    def __init__(self, cwd: str | os.PathLike | None = None, command: list[str] | None = None,
                 rows: int = DEFAULT_ROWS, cols: int = DEFAULT_COLS) -> None:
        import pyte

        self.rows, self.cols = rows, cols
        self.screen = pyte.HistoryScreen(cols, rows, history=SCROLLBACK, ratio=0.5)
        self.screen.set_mode(pyte.modes.LNM)
        self.stream = pyte.ByteStream(self.screen)
        self._out: queue.Queue[bytes] = queue.Queue()
        self._alive = True
        self._proc = self._spawn(list(command or default_shell()), str(cwd) if cwd else None)
        self._reader = threading.Thread(target=self._read_loop, daemon=True)
        self._reader.start()

    # ---- process ----
    def _spawn(self, argv: list[str], cwd: str | None):
        env = dict(os.environ, TERM="xterm-256color", COLUMNS=str(self.cols), LINES=str(self.rows))
        env.pop("LINES", None) if os.name == "nt" else None
        try:
            if os.name == "nt":
                import winpty

                proc = winpty.PtyProcess.spawn(argv, cwd=cwd, env=env, dimensions=(self.rows, self.cols))
            else:
                from ptyprocess import PtyProcessUnicode

                proc = PtyProcessUnicode.spawn(argv, cwd=cwd, env=env, dimensions=(self.rows, self.cols))
        except (ImportError, OSError, FileNotFoundError) as ex:
            raise TerminalError(str(ex)) from ex
        return proc

    def _read_loop(self) -> None:
        while self._alive:
            try:
                data = self._proc.read(READ_CHUNK)
            except EOFError:
                break
            except OSError:
                break
            if not data:
                break
            self._out.put(data.encode("utf-8", "replace") if isinstance(data, str) else data)
        self._alive = False
        self._out.put(b"")

    # ---- input and output ----
    def write(self, data: str | bytes) -> None:
        if not self._alive:
            return
        text = data.decode("utf-8", "replace") if isinstance(data, bytes) else data
        try:
            self._proc.write(text)
        except (OSError, EOFError):
            self._alive = False

    def drain(self) -> bool:
        """Move whatever the shell printed into the screen. True when something changed."""
        changed = False
        while True:
            try:
                data = self._out.get_nowait()
            except queue.Empty:
                break
            if data:
                self.stream.feed(data)
            changed = True
        return changed

    def resize(self, rows: int, cols: int) -> None:
        rows, cols = max(1, rows), max(1, cols)
        if (rows, cols) == (self.rows, self.cols):
            return
        self.rows, self.cols = rows, cols
        self.screen.resize(rows, cols)
        try:
            self._proc.setwinsize(rows, cols)
        except (AttributeError, OSError):
            try:
                self._proc.set_size(rows, cols)      # pywinpty
            except (AttributeError, OSError):
                pass

    @property
    def alive(self) -> bool:
        return self._alive and self._proc.isalive()

    def close(self) -> None:
        self._alive = False
        try:
            self._proc.terminate(force=True)
        except (OSError, AttributeError):
            pass

    # ---- reading the screen ----
    def text(self) -> str:
        return "\n".join(self.screen.display).rstrip()

    def lines(self) -> list[list[Cell]]:
        out = []
        for y in range(self.rows):
            row = self.screen.buffer[y]
            out.append([Cell(row[x].data, row[x].fg, row[x].bg, row[x].bold, row[x].reverse)
                        for x in range(self.cols)])
        return out

    @property
    def cursor(self) -> tuple[int, int]:
        return self.screen.cursor.y, self.screen.cursor.x

    # ---- history (what scrolled off the top) ----
    @property
    def history_above(self) -> int:
        return len(self.screen.history.top)

    @property
    def history_below(self) -> int:
        return len(self.screen.history.bottom)

    def scroll_pages(self, pages: int) -> None:
        """Negative: go back into the history. Positive: come back down."""
        for _ in range(abs(pages)):
            (self.screen.prev_page if pages < 0 else self.screen.next_page)()

    # ---- mouse reporting (set by programs such as htop or vim) ----
    MOUSE_MODES = (1000, 1002, 1003)
    SGR_MODE = 1006

    def mouse_wanted(self) -> bool:
        return any((m << 5) in self.screen.mode for m in self.MOUSE_MODES)

    def mouse_motion_wanted(self) -> bool:
        return any((m << 5) in self.screen.mode for m in (1002, 1003))

    def mouse_report(self, button: int, col: int, row: int, pressed: bool = True) -> str:
        """The escape sequence for one mouse event (columns and rows are 0-based)."""
        col, row = max(0, col) + 1, max(0, row) + 1
        if (self.SGR_MODE << 5) in self.screen.mode:
            return f"\x1b[<{button};{col};{row}{'M' if pressed else 'm'}"
        code = button if pressed else 3
        if col > 223 or row > 223:
            return ""
        return "\x1b[M" + chr(32 + code) + chr(32 + col) + chr(32 + row)


def available() -> tuple[bool, str]:
    """Whether a terminal can run here, and why not when it cannot."""
    try:
        import pyte  # noqa: F401
    except ImportError:
        return False, "pyte"
    if os.name == "nt":
        try:
            import winpty  # noqa: F401
        except ImportError:
            return False, "pywinpty"
    else:
        try:
            import ptyprocess  # noqa: F401
        except ImportError:
            return False, "ptyprocess"
    if os.name == "nt" and sys.getwindowsversion().build < 17763:   # ConPTY needs Windows 10 1809
        return False, "ConPTY"
    return True, ""
