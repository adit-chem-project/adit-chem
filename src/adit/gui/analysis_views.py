
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (QAbstractItemView, QHeaderView, QLabel, QTableWidget, QTableWidgetItem, QToolButton, QVBoxLayout, QWidget)

MAX_COLUMN_PX = 360


class Collapsible(QWidget):

    def __init__(self, title: str, expanded: bool = False, parent: QWidget | None = None):
        super().__init__(parent)
        self.button = QToolButton()
        self.button.setObjectName("collapsible")
        self.button.setText(title)
        self.button.setCheckable(True)
        self.button.setAutoRaise(True)
        self.button.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self.button.setStyleSheet("font-weight: 600;")
        self.body = QWidget()
        self.body_layout = QVBoxLayout(self.body)
        self.body_layout.setContentsMargins(18, 4, 0, 4)
        self.body_layout.setSpacing(6)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(2)
        lay.addWidget(self.button)
        lay.addWidget(self.body)
        self.button.toggled.connect(self.set_expanded)
        self.set_expanded(expanded)

    def set_expanded(self, on: bool) -> None:
        if self.button.isChecked() != on:
            self.button.setChecked(on)
            return
        self.button.setArrowType(Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow)
        self.body.setVisible(on)

    def is_expanded(self) -> bool:
        return self.button.isChecked()


def hint_label(text: str) -> QLabel:
    w = QLabel(text)
    w.setObjectName("hint")
    w.setWordWrap(True)
    w.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
    return w


class SectionView(Collapsible):

    def __init__(self, sec, parent: QWidget | None = None):
        super().__init__(sec.title, expanded=True, parent=parent)
        self.section = sec
        self.table: QTableWidget | None = None
        if sec.columns and sec.rows:
            t = QTableWidget(len(sec.rows), len(sec.columns))
            t.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            t.verticalHeader().setVisible(False)
            t.setHorizontalHeaderLabels(sec.columns)
            right = Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            for i, row in enumerate(sec.rows):
                for j, v in enumerate(row):
                    item = QTableWidgetItem(v)
                    item.setToolTip(v)
                    if j < len(sec.numeric) and sec.numeric[j]:
                        item.setTextAlignment(right)
                    t.setItem(i, j, item)
            t.resizeColumnsToContents()
            for j in range(t.columnCount()):
                if t.columnWidth(j) > MAX_COLUMN_PX:
                    t.setColumnWidth(j, MAX_COLUMN_PX)
            t.horizontalHeader().setStretchLastSection(not (sec.numeric and sec.numeric[-1]))
            t.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
            t.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            h = (t.horizontalHeader().sizeHint().height() + sum(t.rowHeight(i) for i in range(t.rowCount())) + 2 * t.frameWidth()
                 + t.horizontalScrollBar().sizeHint().height() + 2)
            t.setFixedHeight(h)
            self.table = t
            self.body_layout.addWidget(t)
        self.notes = [hint_label(n) for n in [*sec.notes, sec.more_text()] if n]
        for w in self.notes:
            self.body_layout.addWidget(w)


def _is_wsl() -> bool:
    if not sys.platform.startswith("linux"):
        return False
    try:
        return "microsoft" in Path("/proc/version").read_text().lower()
    except OSError:
        return False


def open_folder(path: str | Path) -> bool:
    p = Path(path).expanduser()
    if not p.is_dir():
        return False
    if _is_wsl() and shutil.which("explorer.exe"):
        try:
            win = subprocess.run(["wslpath", "-w", str(p)], capture_output=True, text=True, errors="replace", timeout=10, check=False).stdout.strip()
            subprocess.Popen(["explorer.exe", win or str(p)])
            return True
        except (OSError, subprocess.SubprocessError):
            pass
    return bool(QDesktopServices.openUrl(QUrl.fromLocalFile(str(p))))


__all__ = ["Collapsible", "SectionView", "hint_label", "open_folder"]
