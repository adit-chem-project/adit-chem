
from __future__ import annotations

import re

from PySide6.QtCore import QEvent, QObject, Qt
from PySide6.QtGui import QFont, QFontMetrics, QValidator
from PySide6.QtWidgets import QApplication, QComboBox, QDoubleSpinBox, QFormLayout, QLabel, QWidget

from adit.gui.style import LABEL_WIDTH, NARROW_FIELD

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


class FieldLabel(QLabel):
    """Form label whose text stays the help / translation key; a "required" pill sits at its right edge."""

    def __init__(self, text: str, parent: QWidget | None = None):
        super().__init__(text, parent)
        self.pill: QLabel | None = None

    def mark_required(self) -> None:
        from adit.lang import L

        self.setObjectName("required")
        pill = QLabel(L("必須", "required"), self)
        pill.setObjectName("required_pill")
        pill.setAlignment(Qt.AlignmentFlag.AlignCenter)
        font = QFont(pill.font()); font.setPointSizeF(PILL_POINT_SIZE)
        pill.setFixedSize(QFontMetrics(font).horizontalAdvance(pill.text()) + 14, PILL_HEIGHT)
        self.setContentsMargins(0, 0, pill.width() + PILL_GAP, 0)   # keeps the text clear of the pill
        self.pill = pill

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        if self.pill is not None:
            x = min(self.fontMetrics().horizontalAdvance(self.text()) + PILL_GAP, self.width() - self.pill.width())
            self.pill.move(x, (self.height() - self.pill.height()) // 2)


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
    w.setMaximumWidth(NARROW_FIELD)
    return w


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
