
from __future__ import annotations

import re

from PySide6.QtCore import QEvent, QObject, Qt
from PySide6.QtGui import QFont, QFontMetrics, QValidator
from PySide6.QtWidgets import QApplication, QComboBox, QDoubleSpinBox, QFormLayout, QLabel, QWidget

from adit.gui.style import LABEL_GAP, LABEL_WIDTH, NARROW_FIELD

_SCI = re.compile(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$")
_PARTIAL = re.compile(r"^[+-]?(\d*\.?\d*)([eE][+-]?\d*)?$")


class SciDoubleSpinBox(QDoubleSpinBox):

    def __init__(self, lo: float, hi: float, value: float, step: float, parent: QWidget | None = None):
        super().__init__(parent)
        self.setDecimals(15)
        self.setRange(lo, hi)
        self.setSingleStep(step)
        self.setValue(value)
        self.setKeyboardTracking(False)

    def textFromValue(self, v: float) -> str:  # noqa: N802 
        return f"{v:g}"

    def valueFromText(self, text: str) -> float:  # noqa: N802
        try:
            return float(text.strip())
        except ValueError:
            return self.value()

    def validate(self, text: str, pos: int):  # noqa: N802
        t = text.strip()
        if _SCI.match(t):
            return QValidator.State.Acceptable, text, pos
        if _PARTIAL.match(t):
            return QValidator.State.Intermediate, text, pos
        return QValidator.State.Invalid, text, pos

    def stepBy(self, steps: int) -> None:  # noqa: N802
        v = self.value()
        if v > 0:
            self.setValue(v * (10.0 ** steps))
        else:
            super().stepBy(steps)


PILL_GAP = 8
PILL_HEIGHT = 18
PILL_POINT_SIZE = 8.5
BAR_WIDTH = 2       # the VS Code-style mark on a field whose value differs from the code default
BAR_GAP = 6


class FieldLabel(QLabel):
    """Form label whose text stays the help / translation key; a "required" pill sits at its right edge."""

    def __init__(self, text: str, parent: QWidget | None = None):
        super().__init__(text, parent)
        self.pill: QLabel | None = None
        self.changed_kind = ""      # "", "user" (typed by the user) or "template" (value from a group template)

    def mark_required(self) -> None:
        from adit.lang import L

        self.setObjectName("required")
        pill = QLabel(L("必須", "required"), self)
        pill.setObjectName("required_pill")
        pill.setAlignment(Qt.AlignmentFlag.AlignCenter)
        font = QFont(pill.font()); font.setPointSizeF(PILL_POINT_SIZE)
        pill.setFixedSize(QFontMetrics(font).horizontalAdvance(pill.text()) + 14, PILL_HEIGHT)
        self.pill = pill
        self._apply_margins()

    def set_changed(self, kind: str) -> None:
        if kind == self.changed_kind:
            return
        self.changed_kind = kind
        self._apply_margins()
        self.update()

    def _apply_margins(self) -> None:
        left = BAR_WIDTH + BAR_GAP if self.changed_kind else 0
        right = self.pill.width() + PILL_GAP if self.pill is not None else 0     # keeps the text clear of the pill
        self.setContentsMargins(left, 0, right, 0)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        if self.pill is not None:
            x = min(self.contentsMargins().left() + self.fontMetrics().horizontalAdvance(self.text()) + PILL_GAP,
                    self.width() - self.pill.width())
            self.pill.move(x, (self.height() - self.pill.height()) // 2)

    def paintEvent(self, event) -> None:  # noqa: N802
        super().paintEvent(event)
        if not self.changed_kind:
            return
        from PySide6.QtGui import QColor, QPainter

        from adit.gui.style import current_theme, tokens_for

        t = tokens_for(current_theme())
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor(t.accent if self.changed_kind == "user" else t.muted))
        p.drawRoundedRect(0, 3, BAR_WIDTH, max(4, self.height() - 6), 1, 1)
        p.end()


def label(text: str, *, required: bool | None = None, help_text: str | None = None) -> QLabel:
    from adit.gui.help import help_for

    w = FieldLabel(text)
    w.setProperty("adit_key", text)
    w.setMinimumWidth(LABEL_WIDTH)
    w.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
    h = help_for(text)
    req = required if required is not None else (h.required if h else False)
    tip = help_text if help_text is not None else (h.text() if h else "")
    if req:
        w.mark_required()
    if tip:
        w.setToolTip(tip); w.setStatusTip(tip)
    return w


def add_row(form: QFormLayout, text: str, field, *, required: bool | None = None, help_text: str | None = None) -> None:
    form.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)
    form.setHorizontalSpacing(LABEL_GAP)
    lab = label(text, required=required, help_text=help_text)
    if lab.toolTip() and isinstance(field, QWidget):
        field.setToolTip(lab.toolTip()); field.setStatusTip(lab.toolTip())
    form.addRow(lab, field)


def file_row(edit, title: str, filters: str = "", *, append: bool = False) -> QWidget:
    from PySide6.QtWidgets import QFileDialog, QHBoxLayout, QPushButton

    btn = QPushButton("参照…")
    row = QWidget(); row.setObjectName("rowbox")
    h = QHBoxLayout(row); h.setContentsMargins(0, 0, 0, 0); h.setSpacing(8)
    h.addWidget(edit, 1); h.addWidget(btn, 0, Qt.AlignmentFlag.AlignTop)

    def choose() -> None:
        path, _ = QFileDialog.getOpenFileName(row, title, "", filters or "*")
        if not path:
            return
        if append:
            cur = edit.toPlainText().rstrip("\n")
            edit.setPlainText((cur + "\n" if cur else "") + path)
        else:
            edit.setText(path)

    btn.clicked.connect(choose)
    row.browse = btn
    return row


def narrow(w: QWidget) -> QWidget:
    """Numeric field: grows to the shared width in a form column, may shrink when the column is narrow."""
    w.setMaximumWidth(NARROW_FIELD)
    return w


def unit_row(field: QWidget, unit: str) -> QWidget:
    """A fixed-width numeric field followed by its unit in gray."""
    from PySide6.QtWidgets import QHBoxLayout

    field.setFixedWidth(NARROW_FIELD)
    row = QWidget(); row.setObjectName("rowbox")
    h = QHBoxLayout(row); h.setContentsMargins(0, 0, 0, 0); h.setSpacing(6)
    u = QLabel(unit); u.setObjectName("unit")
    h.addWidget(field); h.addWidget(u); h.addStretch(1)
    return row


POPUP_ROWS = 12


def limit_combo(combo: QComboBox, rows: int = POPUP_ROWS) -> None:
    combo.setMaxVisibleItems(rows)
    combo.view().setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)


class _ComboPopupLimiter(QObject):

    def eventFilter(self, obj, event) -> bool:  # noqa: N802 
        if event.type() == QEvent.Type.Polish and isinstance(obj, QComboBox):
            limit_combo(obj)
        return False


_LIMITER: _ComboPopupLimiter | None = None


def limit_combo_popups(root: QWidget | None = None) -> None:
    global _LIMITER
    app = QApplication.instance()
    if app is not None and _LIMITER is None:
        _LIMITER = _ComboPopupLimiter(app)
        app.installEventFilter(_LIMITER)
    if root is not None:
        for c in root.findChildren(QComboBox):
            limit_combo(c)


def confirm_overwrite(parent: QWidget, out) -> tuple[bool | None, object]:
    """Ask before writing into a non-empty directory. Returns (overwrite, backup): overwrite is None when the user declined.

    The generated file names are not known in advance here, so the backup starts as a snapshot and is pruned to
    the overwritten files by backup.finish() after the write.
    """
    from pathlib import Path

    from PySide6.QtWidgets import QMessageBox

    from adit.lang import L
    from adit.project import Backup, has_files, overwrite_plan

    out = Path(out).expanduser()
    if not has_files(out):
        return False, None
    keep = Backup(out)
    fits = keep.can_snapshot()
    text = overwrite_plan(out, None).message(keep.dir if fits else None, too_large=not fits)
    if QMessageBox.question(parent, L("上書きの確認", "Overwrite?"), text) != QMessageBox.StandardButton.Yes:
        return None, None
    if not fits:
        return True, None
    keep.snapshot()
    return True, keep
