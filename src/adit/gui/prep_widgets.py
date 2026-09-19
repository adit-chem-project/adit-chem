
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QComboBox, QGridLayout, QHBoxLayout, QLabel, QLineEdit, QPushButton, QVBoxLayout, QWidget

from adit.gui.style import NARROW_FIELD
from adit.lang import L
from adit.textparse import short_number
from adit.web.codefields import SHELLS, element_values_from_rows, hubbard_from_rows, hubbard_rows


def _element_combo(elements: list[str], current: str = "") -> QComboBox:
    c = QComboBox()
    c.addItem(L("(元素)", "(element)"), "")
    for e in elements:
        c.addItem(e, e)
    if current and current not in elements:
        c.addItem(L(f"{current} (構造に無い)", f"{current} (not in the structure)"), current)
    c.setCurrentIndex(max(0, c.findData(current)))
    return c


class HubbardTable(QWidget):

    changed = Signal()

    def __init__(self, parent: QWidget | None = None, *, with_j: bool = True):
        super().__init__(parent)
        self.setObjectName("rowbox")
        self.elements: list[str] = []
        self.with_j = with_j
        self.rows: list[tuple[QComboBox, QComboBox, QLineEdit, QLineEdit]] = []
        self.grid = QGridLayout(); self.grid.setContentsMargins(0, 0, 0, 0); self.grid.setHorizontalSpacing(8); self.grid.setVerticalSpacing(4)
        self.btn_add = QPushButton(L("行を追加", "Add row")); self.btn_add.setObjectName("link")
        self.btn_add.setToolTip(L("DFT+U を掛ける元素を 1 つ足します。元素はいまの構造にあるものだけを選べます",
                                  "Add one element to apply DFT+U to; only elements in the current structure can be chosen"))
        self.empty = QLabel(L("DFT+U を使わないなら空のままにします", "Leave empty to not use DFT+U")); self.empty.setObjectName("hint")
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(4)
        lay.addLayout(self.grid); lay.addWidget(self.empty)
        bl = QHBoxLayout(); bl.setContentsMargins(0, 0, 0, 0); bl.addWidget(self.btn_add); bl.addStretch(1); lay.addLayout(bl)
        self.btn_add.clicked.connect(lambda: self.add_row())
        self._relayout()

    def add_row(self, el: str = "", orb: str = "", u: str = "", j: str = "") -> None:
        ec = _element_combo(self.elements, el)
        oc = QComboBox(); oc.setEditable(True); oc.addItems([""] + SHELLS); oc.setCurrentText(orb)
        oc.lineEdit().setPlaceholderText(L("殻 (例 3d)", "shell (e.g. 3d)")); oc.setMinimumContentsLength(4)
        ue = QLineEdit(u); ue.setPlaceholderText(L("U [eV] (必須)", "U [eV] (required)"))
        je = QLineEdit(j); je.setPlaceholderText(L("J [eV] (空欄 = 0)", "J [eV] (empty = 0)"))
        for w in (ue, je):
            w.setFixedWidth(NARROW_FIELD); w.textChanged.connect(self._emit)
        ec.currentIndexChanged.connect(self._emit); oc.currentTextChanged.connect(self._emit)
        self.rows.append((ec, oc, ue, je))
        self._relayout(); self._emit()

    def remove_row(self, i: int) -> None:
        if 0 <= i < len(self.rows):
            for w in self.rows.pop(i):
                w.hide(); w.setParent(None); w.deleteLater()
            self._relayout(); self._emit()

    def _relayout(self) -> None:
        while self.grid.count():
            item = self.grid.takeAt(0)
            w = item.widget()
            if w is not None and getattr(w, "_hub_header", False):
                w.hide(); w.setParent(None); w.deleteLater()
        if self.rows:
            heads = [L("元素", "Element"), L("殻", "Shell"), "U [eV]"] + (["J [eV]"] if self.with_j else [])
            for c, t in enumerate(heads):
                h = QLabel(t); h.setObjectName("hint"); h._hub_header = True; self.grid.addWidget(h, 0, c)
        for r, (ec, oc, ue, je) in enumerate(self.rows, start=1):
            self.grid.addWidget(ec, r, 0); self.grid.addWidget(oc, r, 1); self.grid.addWidget(ue, r, 2)
            if self.with_j:
                self.grid.addWidget(je, r, 3)
            else:
                je.hide()
            rm = QPushButton("×"); rm.setObjectName("link"); rm._hub_header = True
            rm.setToolTip(L("この行を消します", "Remove this row")); rm.setAccessibleName(L("この行を消す", "Remove this row"))
            rm.clicked.connect(lambda _c=False, row=r - 1: self.remove_row(row))
            self.grid.addWidget(rm, r, 4)
        self.grid.setColumnStretch(1, 1)
        self.empty.setVisible(not self.rows)

    def set_elements(self, elements: list[str]) -> None:
        if list(elements) == self.elements:
            return
        self.elements = list(elements)
        for k, (ec, *_rest) in enumerate(self.rows):
            cur = ec.currentData() or ""
            ec.blockSignals(True); ec.clear()
            ec.addItem(L("(元素)", "(element)"), "")
            for e in self.elements:
                ec.addItem(e, e)
            if cur and cur not in self.elements:
                ec.addItem(L(f"{cur} (構造に無い)", f"{cur} (not in the structure)"), cur)
            ec.setCurrentIndex(max(0, ec.findData(cur))); ec.blockSignals(False)
        self._emit()

    def row_texts(self) -> list[tuple[str, str, str, str]]:
        # J is kept even when the column is hidden, so a spec with J survives the round trip.
        return [(ec.currentData() or "", oc.currentText(), ue.text(), je.text()) for ec, oc, ue, je in self.rows]

    def hubbard(self) -> dict:
        return hubbard_from_rows(self.row_texts())

    def set_hubbard(self, hubbard: dict) -> None:
        while self.rows:
            for w in self.rows.pop():
                w.hide(); w.setParent(None); w.deleteLater()
        for el, orb, u, j in hubbard_rows(hubbard):
            self.add_row(el, orb, u, j)
        self._relayout(); self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()


class ElementValues(QWidget):

    changed = Signal()

    def __init__(self, what: str, placeholder: str, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("rowbox")
        self.what = what
        self.placeholder = placeholder
        self.elements: list[str] = []
        self.edits: dict[str, QLineEdit] = {}
        self._kept: dict[str, str] = {}
        self.grid = QGridLayout(self); self.grid.setContentsMargins(0, 0, 0, 0); self.grid.setHorizontalSpacing(8); self.grid.setVerticalSpacing(4)
        self.empty = QLabel(L("構造を作ると、元素ごとの欄が出ます", "Fields appear per element once a structure exists")); self.empty.setObjectName("hint")
        self._rebuild()

    def _rebuild(self) -> None:
        for e, w in self.edits.items():
            self._kept[e] = w.text()
        while self.grid.count():
            item = self.grid.takeAt(0)
            if item.widget() is not None and item.widget() is not self.empty:
                w = item.widget(); w.hide(); w.setParent(None); w.deleteLater()
        self.edits = {}
        cols = 3
        for i, e in enumerate(self.elements):
            r, c = divmod(i, cols)
            lab = QLabel(e); lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            ed = QLineEdit(self._kept.get(e, "")); ed.setPlaceholderText(self.placeholder); ed.setFixedWidth(NARROW_FIELD)
            ed.setAccessibleName(f"{self.what} {e}")
            ed.textChanged.connect(self._emit)
            self.grid.addWidget(lab, r, 2 * c); self.grid.addWidget(ed, r, 2 * c + 1)
            self.edits[e] = ed
        self.grid.setColumnStretch(2 * cols, 1)
        if not self.elements:
            self.grid.addWidget(self.empty, 0, 0, 1, 2 * cols); self.empty.show()
        else:
            self.empty.hide()

    def set_elements(self, elements: list[str]) -> None:
        if list(elements) != self.elements:
            self.elements = list(elements)
            self._rebuild(); self._emit()

    def values(self) -> dict[str, float]:
        return element_values_from_rows([(e, w.text()) for e, w in self.edits.items()], self.what)

    def set_values(self, values: dict[str, float]) -> None:
        self._kept = {e: short_number(v) for e, v in values.items()}
        for e, w in self.edits.items():
            w.blockSignals(True); w.setText(self._kept.get(e, "")); w.blockSignals(False)
        extra = [e for e in values if e not in self.elements]
        if extra:
            self.elements = self.elements + extra
            self._rebuild()
        self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()


class DocValueBox(QWidget):

    use_requested = Signal(str, str)

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("rowbox")
        self.lines: list = []
        self.grid = QGridLayout(self); self.grid.setContentsMargins(0, 0, 0, 0)
        self.grid.setHorizontalSpacing(8); self.grid.setVerticalSpacing(2)
        self.setVisible(False)

    def set_lines(self, lines: list) -> None:
        from adit.web.prep23 import DOC_MAX_LINES, doc_more

        if lines == self.lines:
            return
        self.lines = list(lines)
        while self.grid.count():
            item = self.grid.takeAt(0)
            if item.widget() is not None:
                w = item.widget(); w.hide(); w.setParent(None); w.deleteLater()
        for r, line in enumerate(self.lines[:DOC_MAX_LINES]):
            lab = QLabel(line.text); lab.setObjectName("hint"); lab.setWordWrap(True); lab.setToolTip(line.tip)
            self.grid.addWidget(lab, r, 0)
            if line.field:
                btn = QPushButton(L("入れる", "Use")); btn.setObjectName("link")
                btn.setToolTip(L(f"この値 ({line.value}) を欄に入れます (ADIT が選んだ値ではありません)",
                                 f"Puts this value ({line.value}) in the field (it is not chosen by ADIT)"))
                btn.clicked.connect(lambda _c=False, f=line.field, v=line.value: self.use_requested.emit(f, v))
                self.grid.addWidget(btn, r, 1)
        more = doc_more(self.lines)
        if more:
            lab = QLabel(more); lab.setObjectName("hint"); self.grid.addWidget(lab, DOC_MAX_LINES, 0)
        self.grid.setColumnStretch(0, 1)
        self.setVisible(bool(self.lines))


__all__ = ["HubbardTable", "ElementValues", "DocValueBox", "SHELLS"]
