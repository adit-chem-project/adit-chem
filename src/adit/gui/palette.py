
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from PySide6.QtCore import QEvent, Qt
from PySide6.QtGui import QAction, QIcon, QKeySequence
from PySide6.QtWidgets import (QAbstractItemView, QDialog, QLabel, QLineEdit, QListWidget, QListWidgetItem,
                               QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from adit.lang import L

KIND_ORDER = {"action": 0, "field": 1, "preset": 2}


@dataclass
class Item:
    kind: str
    title: str
    detail: str
    run: Callable[[], None]
    keys: tuple[str, ...] = ()
    icon: QIcon | None = None
    enabled: bool = True
    score: int = field(default=0, compare=False)


def shortcut_text(action: QAction) -> str:
    return action.shortcut().toString(QKeySequence.SequenceFormat.NativeText) if not action.shortcut().isEmpty() else ""


def with_shortcut(tip: str, action: QAction) -> str:
    # "元に戻す (Ctrl+Z)": the shortcut is appended once, after the tooltip (or the name when there is none).
    key = shortcut_text(action)
    # Qt derives a default tooltip from the text without its "…"; use the text itself so it still translates.
    base = tip if tip and tip != action.text().rstrip("…") else action.text()
    if not key or f"({key})" in base:
        return base
    return f"{base} ({key})"


def split_shortcut(tip: str, action: QAction) -> tuple[str, str]:
    key = shortcut_text(action)
    suffix = f" ({key})"
    if key and tip.endswith(suffix):
        return tip[:-len(suffix)], suffix
    return tip, ""


def _match(needle: str, hay: str) -> int:
    # 0 = starts with, 1 = a word starts with, 2 = contained, -1 = no match.
    if not hay:
        return -1
    if hay.startswith(needle):
        return 0
    pos = hay.find(needle)
    if pos < 0:
        return -1
    return 1 if hay[pos - 1] in " (/[・" else 2


def search(items: list[Item], query: str, limit: int = 40) -> list[Item]:
    words = query.casefold().split()
    out: list[Item] = []
    for it in items:
        hays = [h.casefold() for h in (it.title, *it.keys) if h]
        score = 0
        for w in words:
            best = min((m for m in (_match(w, h) for h in hays) if m >= 0), default=-1)
            if best < 0:
                break
            score += best
        else:
            it.score = score
            out.append(it)
    out.sort(key=lambda it: (it.score, KIND_ORDER.get(it.kind, 9), it.title))
    return out[:limit]


class CommandPalette(QDialog):

    def __init__(self, items: list[Item], parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("palette")
        self.setWindowFlags(Qt.WindowType.Dialog | Qt.WindowType.FramelessWindowHint)
        self.setModal(True)
        self.items = items
        self.ran: Item | None = None
        self.search = QLineEdit(); self.search.setObjectName("palette_search")
        self.search.setPlaceholderText(L("操作・欄・プリセットを名前で検索 (例: 温度, generate, 水)",
                                         "Search commands, fields and presets by name (e.g. temperature, 生成, water)"))
        self.search.setClearButtonEnabled(True)
        self.list = QListWidget(); self.list.setObjectName("palette_list")
        self.list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.list.setUniformItemSizes(True); self.list.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.hint = QLabel(L("↑↓ で選ぶ、Enter で実行、Esc で閉じる", "↑↓ to choose, Enter to run, Esc to close")); self.hint.setObjectName("hint")
        lay = QVBoxLayout(self); lay.setContentsMargins(12, 12, 12, 10); lay.setSpacing(8)
        lay.addWidget(self.search); lay.addWidget(self.list, 1); lay.addWidget(self.hint)
        self.setMinimumWidth(560); self.resize(620, 420)
        self.search.textChanged.connect(self.refresh)
        self.search.installEventFilter(self)
        self.list.itemActivated.connect(lambda *_: self.run_current())
        self.list.itemClicked.connect(lambda *_: self.run_current())
        self.refresh()

    def refresh(self) -> None:
        self.list.clear()
        for it in search(self.items, self.search.text()):
            text = it.title if not it.detail else f"{it.title}    {it.detail}"
            row = QListWidgetItem(it.icon or QIcon(), text)
            row.setData(Qt.ItemDataRole.UserRole, it)
            if not it.enabled:
                row.setFlags(row.flags() & ~Qt.ItemFlag.ItemIsEnabled)
            self.list.addItem(row)
        self._select_first_enabled()

    def _select_first_enabled(self) -> None:
        for i in range(self.list.count()):
            if self.list.item(i).flags() & Qt.ItemFlag.ItemIsEnabled:
                self.list.setCurrentRow(i); return

    def current_item(self) -> Item | None:
        row = self.list.currentItem()
        return row.data(Qt.ItemDataRole.UserRole) if row is not None else None

    def shown(self) -> list[Item]:
        return [self.list.item(i).data(Qt.ItemDataRole.UserRole) for i in range(self.list.count())]

    def run_current(self) -> bool:
        it = self.current_item()
        if it is None or not it.enabled:
            return False
        self.ran = it
        self.accept()
        it.run()
        return True

    def _move(self, delta: int) -> None:
        n = self.list.count()
        if n == 0:
            return
        row = (self.list.currentRow() + delta) % n
        self.list.setCurrentRow(row)

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if obj is self.search and event.type() == QEvent.Type.KeyPress:
            key = event.key()
            if key in (Qt.Key.Key_Down, Qt.Key.Key_PageDown):
                self._move(1); return True
            if key in (Qt.Key.Key_Up, Qt.Key.Key_PageUp):
                self._move(-1); return True
            if key in (Qt.Key.Key_Return, Qt.Key.Key_Enter):
                self.run_current(); return True
            if key == Qt.Key.Key_Escape:
                self.reject(); return True
        return super().eventFilter(obj, event)

    def place_over(self, window: QWidget) -> None:
        g = window.frameGeometry()
        self.move(g.left() + (g.width() - self.width()) // 2, g.top() + min(120, g.height() // 6))


class ShortcutsDialog(QDialog):

    def __init__(self, rows: list[tuple[str, str]], parent: QWidget | None = None):
        super().__init__(parent)
        self.setWindowTitle(L("キーボードショートカット", "Keyboard shortcuts"))
        self.table = QTableWidget(len(rows), 2)
        self.table.setHorizontalHeaderLabels([L("キー", "Key"), L("操作", "Command")])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.table.setShowGrid(False)
        for r, (key, what) in enumerate(rows):
            k = QTableWidgetItem(key); k.setFont(self._mono())
            self.table.setItem(r, 0, k); self.table.setItem(r, 1, QTableWidgetItem(what))
        self.table.resizeColumnsToContents(); self.table.horizontalHeader().setStretchLastSection(True)
        note = QLabel(L("Ctrl+K で、操作・欄・プリセットを名前で探せます", "Ctrl+K searches commands, fields and presets by name"))
        note.setObjectName("hint")
        lay = QVBoxLayout(self); lay.setContentsMargins(14, 14, 14, 14); lay.setSpacing(8)
        lay.addWidget(self.table, 1); lay.addWidget(note)
        self.resize(520, min(560, 80 + 30 * (len(rows) + 1)))

    @staticmethod
    def _mono():
        from PySide6.QtGui import QFont
        f = QFont(); f.setFamily("monospace"); return f

    def rows(self) -> list[tuple[str, str]]:
        return [(self.table.item(r, 0).text(), self.table.item(r, 1).text()) for r in range(self.table.rowCount())]


__all__ = ["Item", "CommandPalette", "ShortcutsDialog", "search", "shortcut_text", "with_shortcut", "split_shortcut"]
