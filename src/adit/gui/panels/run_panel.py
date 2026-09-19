
from __future__ import annotations

import re
import time
from pathlib import Path

from PySide6.QtCore import QTimer
from PySide6.QtWidgets import QLabel, QPlainTextEdit, QProgressBar, QVBoxLayout, QWidget

from adit.lang import L

STEP_RE = re.compile(r"Geometry step:\s*(\d+)")


def last_step(text: str) -> int | None:
    found = STEP_RE.findall(text)
    return int(found[-1]) if found else None


def progress_text(step: int, total: int | None, elapsed_s: float) -> str:
    done = step + 1
    rate = elapsed_s / done if done > 0 else 0.0
    head = L(f"ステップ {done}", f"step {done}") + (f" / {total} ({100 * done / total:.0f} %)" if total else "")
    speed = L(f"、1 ステップ {rate:.2g} 秒", f", measured {rate:.2g} s/step") if elapsed_s > 0 else ""
    left = ""
    if total and elapsed_s > 0 and done < total:
        rest = rate * (total - done)
        left = L(f"、残り約 {rest / 60:.0f} 分", f", about {rest / 60:.0f} min left") if rest >= 90 else L(f"、残り約 {rest:.0f} 秒", f", about {rest:.0f} s left")
    return head + speed + left


class RunPanel(QWidget):
    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.status = QLabel(L("まだ実行していません", "Not run yet")); self.status.setObjectName("title")
        self.where = QLabel(""); self.where.setObjectName("hint"); self.where.setWordWrap(True)
        self.progress = QProgressBar(); self.progress.setTextVisible(False); self.progress.hide()
        self.progress_label = QLabel(""); self.progress_label.setObjectName("hint"); self.progress_label.hide()
        self.log = QPlainTextEdit(); self.log.setReadOnly(True); self.log.setStyleSheet("font-family: 'SF Mono', Menlo, monospace;")
        lay = QVBoxLayout(self); lay.setContentsMargins(12, 12, 12, 12); lay.setSpacing(10)
        lay.addWidget(self.status); lay.addWidget(self.where); lay.addWidget(self.progress); lay.addWidget(self.progress_label); lay.addWidget(self.log, 1)
        self._dir: Path | None = None
        self._total: int | None = None
        self._t0 = 0.0
        self._timer = QTimer(self); self._timer.setInterval(1000); self._timer.timeout.connect(self.refresh)

    def start(self, run_dir: Path, total_steps: int | None = None) -> None:
        self._dir = Path(run_dir); self._total = total_steps if total_steps and total_steps > 0 else None
        self._t0 = time.monotonic()
        self.status.setText(L("実行中", "Running")); self.where.setText(str(self._dir)); self.log.setPlainText("")
        self.progress.setRange(0, self._total or 0); self.progress.setValue(0)
        self.progress.show(); self.progress_label.setText(L("開始しました", "started")); self.progress_label.show()
        self._timer.start(); self.refresh()

    def finish(self, code: int, summary: str) -> None:
        self._timer.stop(); self.refresh()
        self.status.setText(L(f"終了 (終了コード {code})", f"Finished (exit code {code})"))
        self.progress.hide()
        p = self._dir / "output.log" if self._dir else None
        step = last_step(p.read_text(encoding="utf-8", errors="replace")) if p and p.is_file() else None
        took = time.monotonic() - self._t0
        took_s = L(f"{took / 60:.1f} 分", f"{took / 60:.1f} min") if took >= 90 else L(f"{took:.0f} 秒", f"{took:.0f} s")
        self.progress_label.setText(L(f"{step + 1} ステップで終了しました (かかった時間 {took_s})", f"finished after {step + 1} steps ({took_s})")
                                    if step is not None else L(f"終了しました (かかった時間 {took_s})", f"finished ({took_s})"))
        self.log.appendPlainText("\n" + summary)

    def refresh(self) -> None:
        if self._dir is None:
            return
        p = self._dir / "output.log"
        if p.is_file():
            text = p.read_text(encoding="utf-8", errors="replace")
            self.log.setPlainText("\n".join(text.splitlines()[-200:]))
            self.log.verticalScrollBar().setValue(self.log.verticalScrollBar().maximum())
            step = last_step(text)
            if step is not None:
                if self._total:
                    self.progress.setValue(min(step + 1, self._total))
                self.progress_label.setText(progress_text(step, self._total, time.monotonic() - self._t0))
