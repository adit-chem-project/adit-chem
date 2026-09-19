
from __future__ import annotations

import threading
import time
from typing import Callable

from PySide6.QtCore import QObject, Qt, QTimer, Signal
from PySide6.QtWidgets import QHBoxLayout, QLabel, QProgressBar, QPushButton, QSizePolicy, QVBoxLayout, QWidget

from adit import progress as reports
from adit.lang import L

SHOW_AFTER = 1.0       # seconds before anything is shown
PERCENT_AFTER = 10.0   # seconds before a known fraction is shown as a percentage


class Job(QObject):
    # Runs one function on a worker thread; a cancelled or superseded run has its result dropped.

    finished = Signal(object, object)    # result, exception or None
    progress = Signal(int, int, str)
    _done = Signal(int, object, object)
    _report = Signal(int, int, int, str)

    def __init__(self, parent: QObject | None = None):
        super().__init__(parent)
        self._token = 0
        self._running = False
        self._done.connect(self._on_done, Qt.ConnectionType.QueuedConnection)
        self._report.connect(self._on_report, Qt.ConnectionType.QueuedConnection)

    def is_running(self) -> bool:
        return self._running

    def start(self, fn: Callable, *args, **kwargs) -> int:
        self._token += 1
        token = self._token
        self._running = True
        threading.Thread(target=self._work, args=(token, fn, args, kwargs), daemon=True).start()
        return token

    def cancel(self) -> None:
        if self._running:
            self._token += 1
            self._running = False

    def _work(self, token: int, fn: Callable, args, kwargs) -> None:
        reports.set_reporter(lambda d, t, w="": self._emit(self._report, token, d, t, w))
        try:
            result, error = fn(*args, **kwargs), None
        except Exception as ex:  # reported to the caller as the error
            result, error = None, ex
        finally:
            reports.set_reporter(None)
        self._emit(self._done, token, result, error)

    @staticmethod
    def _emit(signal, *args) -> None:
        try:
            signal.emit(*args)
        except RuntimeError:      # the owner was deleted while the work was still running
            pass

    def _on_done(self, token: int, result, error) -> None:
        if token != self._token:
            return
        self._running = False
        self.finished.emit(result, error)

    def _on_report(self, token: int, done: int, total: int, what: str) -> None:
        if token == self._token:
            self.progress.emit(done, total, what)


class ProgressStrip(QWidget):
    # Nothing under 1 s, a looping bar with Cancel from 1 s, a percentage from 10 s when the fraction is known.

    cancel_requested = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("progress_strip")
        self.clock = time.monotonic
        self.label = QLabel(""); self.label.setObjectName("hint")
        self.bar = QProgressBar(); self.bar.setTextVisible(False); self.bar.setRange(0, 0)
        self.bar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.btn_cancel = QPushButton(L("中止", "Cancel"))
        self.btn_cancel.clicked.connect(self.cancel_requested.emit)
        row = QHBoxLayout(); row.setContentsMargins(0, 0, 0, 0); row.setSpacing(8)
        row.addWidget(self.bar, 1); row.addWidget(self.btn_cancel)
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(4)
        lay.addWidget(self.label); lay.addLayout(row)
        self._timer = QTimer(self); self._timer.setInterval(200); self._timer.timeout.connect(self._tick)
        self._active = False
        self._t0 = 0.0
        self._what = ""
        self._done: int | None = None
        self._total: int | None = None
        self._cancellable = True
        self.hide()

    def begin(self, what: str = "", *, cancellable: bool = True) -> None:
        self._active, self._t0, self._what = True, self.clock(), what
        self._done = self._total = None
        self._cancellable = cancellable
        self.btn_cancel.setVisible(cancellable)
        self.hide()
        self._timer.start()

    def update(self, done: int, total: int, what: str = "") -> None:
        if not self._active:
            return
        self._done, self._total = done, total if total > 0 else None
        if what:
            self._what = what
        if not self.isHidden():
            self._paint()

    def end(self) -> None:
        self._active = False
        self._timer.stop()
        self.hide()

    def phase(self) -> str:
        if not self._active:
            return "idle"
        elapsed = self.clock() - self._t0
        if elapsed < SHOW_AFTER:
            return "hidden"
        if elapsed >= PERCENT_AFTER and self._total:
            return "percent"
        return "busy"

    def _tick(self) -> None:
        if self.phase() == "hidden":
            return
        self._paint()
        if self.isHidden():
            self.show()

    def _paint(self) -> None:
        phase = self.phase()
        if phase == "percent":
            self.bar.setRange(0, 100)
            self.bar.setValue(int(round(100.0 * self._done / self._total)))
            self.label.setText(self._text_with_estimate())
        else:
            self.bar.setRange(0, 0)
            self.label.setText(self._text())

    def _text(self) -> str:
        count = f"  {self._done} / {self._total}" if self._done is not None and self._total else ""
        return (self._what or L("処理しています…", "Working…")) + count

    def _text_with_estimate(self) -> str:
        elapsed = self.clock() - self._t0
        pct = int(round(100.0 * self._done / self._total))
        text = f"{self._what or L('処理しています…', 'Working…')}  {pct} %"
        if 0 < self._done < self._total:
            left = elapsed * (self._total - self._done) / self._done
            text += L(f" (残り約 {max(1, round(left))} 秒)", f" (about {max(1, round(left))} s left)")
        return text


class EmptyState(QWidget):
    # Heading + one line + the next action, centred in an area that has nothing to show yet.

    clicked = Signal()

    def __init__(self, title: str, line: str, button: str = "", parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("empty_state")
        self.title = QLabel(title); self.title.setObjectName("title"); self.title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.line = QLabel(line); self.line.setObjectName("hint"); self.line.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.line.setWordWrap(True)
        self.button = QPushButton(button); self.button.setObjectName("primary")
        self.button.setVisible(bool(button)); self.button.clicked.connect(self.clicked.emit)
        self.button.setCursor(Qt.CursorShape.PointingHandCursor)
        lay = QVBoxLayout(self); lay.setContentsMargins(24, 12, 24, 12); lay.setSpacing(8)
        lay.addStretch(1); lay.addWidget(self.title); lay.addWidget(self.line)
        lay.addSpacing(6); lay.addWidget(self.button, 0, Qt.AlignmentFlag.AlignHCenter); lay.addStretch(1)

    def set_texts(self, title: str | None = None, line: str | None = None, button: str | None = None) -> None:
        if title is not None:
            self.title.setText(title)
        if line is not None:
            self.line.setText(line)
        if button is not None:
            self.button.setText(button); self.button.setVisible(bool(button))


__all__ = ["Job", "ProgressStrip", "EmptyState", "SHOW_AFTER", "PERCENT_AFTER"]
