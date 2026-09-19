
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QFormLayout, QLabel, QSpinBox, QWidget

from adit.codes.xtb import solvent_names
from adit.gui.style import ROW_SPACING
from adit.gui.widgets import SciDoubleSpinBox, add_row, label, narrow
from adit.lang import L
from adit.spec import XtbMethod

GFN = {"2": "GFN2-xTB", "1": "GFN1-xTB", "0": "GFN0-xTB", "ff": "GFN-FF (力場)"}
OPT_LEVELS = ["crude", "sloppy", "loose", "normal", "tight", "verytight", "extreme"]
SOLVATION = {"none": "なし", "alpb": "ALPB", "gbsa": "GBSA"}


def fill_solvents(combo: QComboBox, names: list[str], current: str) -> None:
    combo.blockSignals(True)
    combo.clear()
    combo.addItem(L("(選んでください)", "(choose one)") if names else L("(この組み合わせでは候補がありません)", "(no solvents for this combination)"), "")
    for n in names:
        combo.addItem(n, n)
    if current and current not in names:
        combo.addItem(current, current)
    combo.setCurrentIndex(max(0, combo.findData(current)))
    combo.blockSignals(False)


class XtbMethodPanel(QWidget):
    changed = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.gfn = QComboBox()
        for k, v in GFN.items():
            self.gfn.addItem(v, k)
        self.accuracy = narrow(SciDoubleSpinBox(1e-4, 1000, 1.0, 0.1))
        self.etemp = narrow(SciDoubleSpinBox(0, 1e5, 300.0, 50))
        self.max_iter = narrow(QSpinBox()); self.max_iter.setRange(1, 100000); self.max_iter.setValue(250)
        self.opt_level = QComboBox(); self.opt_level.addItems(OPT_LEVELS); self.opt_level.setCurrentText("normal")
        self.solvation = QComboBox()
        for k, v in SOLVATION.items():
            self.solvation.addItem(v, k)
        self.solvent = QComboBox()
        self.md_hmass = narrow(SciDoubleSpinBox(0.01, 1000, 4.0, 1))
        self.md_shake = QComboBox()
        for value, text in ((0, L("固定しない", "Off")), (1, L("X–H 結合", "X–H bonds")), (2, L("すべての結合", "All bonds"))):
            self.md_shake.addItem(text, value)
        self.md_sccacc = narrow(SciDoubleSpinBox(1e-6, 1e6, 2.0, 0.5))
        info = QLabel("パラメータファイルは不要です (計算手法に内蔵)。分子系のみ。電荷とスピン多重度は構造パネルの値をコマンドライン引数で渡します")
        info.setObjectName("hint"); info.setWordWrap(True)
        form = QFormLayout(self); form.setVerticalSpacing(ROW_SPACING); self._form = form
        add_row(form, "計算手法 (--gfn)", self.gfn)
        add_row(form, "精度 (--acc)", self.accuracy)
        add_row(form, "電子温度 [K] (--etemp)", self.etemp)
        add_row(form, "SCC の最大反復回数", self.max_iter)
        add_row(form, "最適化の収束レベル (--opt)", self.opt_level)
        add_row(form, "溶媒モデル (--alpb / --gbsa)", self.solvation)
        add_row(form, "溶媒", self.solvent)
        add_row(form, L("MD の水素質量 [u] ($md hmass)", "Hydrogen mass in MD [u] ($md hmass)"), self.md_hmass)
        add_row(form, L("MD の結合固定 ($md shake)", "Bond constraints in MD ($md shake)"), self.md_shake)
        add_row(form, L("MD の SCC 精度 ($md sccacc)", "SCC accuracy in MD ($md sccacc)"), self.md_sccacc)
        form.addRow(label(""), info)
        for w in (self.opt_level,):
            w.currentTextChanged.connect(self._emit)
        for w in (self.gfn, self.solvation):
            w.currentIndexChanged.connect(self._refill)
        self.solvent.currentIndexChanged.connect(self._emit)
        for w in (self.accuracy, self.etemp, self.max_iter, self.md_hmass, self.md_sccacc):
            w.valueChanged.connect(self._emit)
        self.md_shake.currentIndexChanged.connect(self._emit)
        self._refill()

    def _refill(self, *_, keep: str | None = None) -> None:
        model = self.solvation.currentData()
        cur = keep if keep is not None else (self.solvent.currentData() or "")
        fill_solvents(self.solvent, solvent_names(model, self.gfn.currentData()) if model != "none" else [], cur if model != "none" else cur)
        self._form.setRowVisible(self.solvent, model != "none" or bool(cur))
        self._emit()

    def method(self) -> XtbMethod:
        return XtbMethod(gfn=self.gfn.currentData(), accuracy=self.accuracy.value(), etemp=self.etemp.value(),
                         max_iterations=self.max_iter.value(), opt_level=self.opt_level.currentText(),
                         solvation=self.solvation.currentData(), solvent=self.solvent.currentData() or "",
                         md_hmass=self.md_hmass.value(), md_shake=self.md_shake.currentData(), md_sccacc=self.md_sccacc.value())

    def set_method(self, m: XtbMethod) -> None:
        self.gfn.blockSignals(True); self.solvation.blockSignals(True)
        self.gfn.setCurrentIndex(list(GFN).index(m.gfn)); self.solvation.setCurrentIndex(max(0, self.solvation.findData(m.solvation)))
        self.gfn.blockSignals(False); self.solvation.blockSignals(False)
        self.accuracy.setValue(m.accuracy); self.etemp.setValue(m.etemp)
        self.max_iter.setValue(m.max_iterations); self.opt_level.setCurrentText(m.opt_level)
        self.md_hmass.setValue(m.md_hmass); self.md_shake.setCurrentIndex(self.md_shake.findData(m.md_shake)); self.md_sccacc.setValue(m.md_sccacc)
        self._refill(keep=m.solvent)

    def _emit(self, *_) -> None:
        self.changed.emit()
