
from __future__ import annotations

import re

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QGridLayout, QHBoxLayout, QLabel,
                               QLineEdit, QListWidget, QPushButton, QSizePolicy, QSpinBox, QStackedWidget, QVBoxLayout, QWidget)

from adit.builder.model import (Adsorb, Box, Fix, MoleculeRef, Remove, Selection, Slab, SolventLayer, Solvate, Substitute, Supercell,
                               Vacuum, interface_steps, op_label)
from adit.gui import icons
from adit.gui.i18n import translate_widgets
from adit.gui.panels.mixture_editor import KINDS, MixtureEditor, RefCell
from adit.gui.style import ROW_SPACING
from adit.gui.widgets import SciDoubleSpinBox, add_row, label, limit_combo, narrow
from adit.lang import L
from adit.textparse import short_number
from adit.mixture import Component
from adit.structure import pretty_formula

from ase.data import chemical_symbols

ELEMENTS = list(chemical_symbols[1:104])

ELECTROLYTE_EXAMPLE = [Component(kind="smiles", ref="O=C1OCCO1", count=32, label="C3H4O3"),
                       Component(kind="smiles", ref="[Li+]", count=1, charge=1, label="Li+"),
                       Component(kind="smiles", ref="F[P-](F)(F)(F)(F)F", count=1, charge=-1, label="PF6-")]

ADD_ORDER = ["slab", "solvent_layer", "fix", "supercell", "remove", "substitute", "adsorb", "vacuum", "box", "solvate"]
INTERFACE = "interface"


def default_step(op: str):
    if op == "supercell":
        return Supercell(repeat=(2, 2, 1))
    if op == "slab":
        return Slab(miller=(1, 0, 0))
    if op == "vacuum":
        return Vacuum()
    if op == "box":
        return Box()
    if op == "adsorb":
        return Adsorb(molecule=MoleculeRef(kind="preset", ref="H2O"), xy=(0.0, 0.0))
    if op == "remove":
        return Remove(count=1)
    if op == "substitute":
        return Substitute(count=1, to="Al")
    if op == "solvent_layer":
        return SolventLayer(components=[Component(ref="H2O", count=32, label="H2O")])
    if op == "solvate":
        return Solvate(components=[Component(ref="H2O", count=0, label="H2O")])
    if op == "fix":
        return Fix(bottom_layers=1)
    raise KeyError(op)


class Num(SciDoubleSpinBox):

    def stepBy(self, steps: int) -> None:  # noqa: N802 
        QDoubleSpinBox.stepBy(self, steps)


def num(lo: float, hi: float, value: float, step: float = 0.5) -> Num:
    return narrow(Num(lo, hi, value, step))


def spin(lo: int, hi: int, value: int, width: int = 72) -> QSpinBox:
    w = QSpinBox(); w.setRange(lo, hi); w.setValue(value); w.setMaximumWidth(width)
    return w


def hrow(*widgets, stretch: bool = True) -> QWidget:
    w = QWidget(); w.setObjectName("rowbox"); h = QHBoxLayout(w); h.setContentsMargins(0, 0, 0, 0)
    for x in widgets:
        h.addWidget(QLabel(x) if isinstance(x, str) else x)
    if stretch:
        h.addStretch()
    return w


def hint(text: str) -> QLabel:
    lab = QLabel(text); lab.setObjectName("hint"); lab.setWordWrap(True)
    return lab


def _float_or_none(text: str, what: str) -> float | None:
    t = text.strip()
    if not t:
        return None
    try:
        return float(t)
    except ValueError as ex:
        raise ValueError(L(f"{what} を数値として読めません: {t!r}", f"{what} is not a number: {t!r}")) from ex


def _fmt(x: float | None) -> str:
    return "" if x is None else short_number(x)


class StepEditor(QWidget):
    changed = Signal()
    draw_requested = Signal(object, int)

    op = ""

    def __init__(self, step, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("rowbox")
        self._model = step
        self.form = QFormLayout(self); self.form.setContentsMargins(0, 0, 0, 0); self.form.setVerticalSpacing(ROW_SPACING)
        self.build()
        self.load(step)
        self.connect_all()

    def build(self) -> None: ...
    def load(self, step) -> None: ...
    def values(self) -> dict: return {}
    def summary(self, step) -> str: return ""

    def step(self):
        data = self._model.model_dump()
        data.update(self.values())
        return type(self._model).model_validate(data)

    def set_step(self, step) -> None:
        self._model = step
        self.blockSignals(True)
        try:
            self.load(step)
        finally:
            self.blockSignals(False)

    def connect_all(self) -> None:
        def own(w: QWidget) -> bool:
            p = w.parentWidget()
            while p is not None and p is not self:
                if isinstance(p, (MixtureEditor, RefCell)):
                    return False
                p = p.parentWidget()
            return True

        for w in self.findChildren(QSpinBox) + self.findChildren(QDoubleSpinBox):
            if own(w):
                w.valueChanged.connect(self._emit)
        for w in self.findChildren(QComboBox):
            if own(w):
                (w.currentIndexChanged if not w.isEditable() else w.currentTextChanged).connect(self._emit)
        for w in self.findChildren(QLineEdit):
            if own(w) and not isinstance(w.parentWidget(), (QSpinBox, QDoubleSpinBox, QComboBox)):
                w.editingFinished.connect(self._emit)
        for w in self.findChildren(QCheckBox):
            if own(w):
                w.toggled.connect(self._emit)

    def _emit(self, *_) -> None:
        if not self.signalsBlocked():
            self.changed.emit()


class SelectionRows:

    def __init__(self, form: QFormLayout):
        self.elements = QLineEdit(); self.elements.setPlaceholderText(L("例: Li (空欄なら全元素)", "e.g. Li (empty = all elements)"))
        self.z_min = QLineEdit(); self.z_max = QLineEdit()
        for w, ja, en in ((self.z_min, "下限 (空欄なら無し)", "min (empty = none)"), (self.z_max, "上限 (空欄なら無し)", "max (empty = none)")):
            w.setPlaceholderText(L(ja, en)); w.setMaximumWidth(150)
        add_row(form, "元素 (条件)", self.elements)
        add_row(form, "z の範囲 [Å]", hrow(self.z_min, "〜", self.z_max))

    def load(self, sel: Selection) -> None:
        self.elements.setText(" ".join(sel.elements)); self.z_min.setText(_fmt(sel.z_min)); self.z_max.setText(_fmt(sel.z_max))

    def value(self, old: Selection) -> dict:
        els = [e for e in re.split(r"[,\s、]+", self.elements.text().strip()) if e]
        return {"elements": els, "indices": list(old.indices), "z_min": _float_or_none(self.z_min.text(), L("z の下限", "z min")),
                "z_max": _float_or_none(self.z_max.text(), L("z の上限", "z max"))}


def describe_selection(sel: Selection) -> str:
    parts = []
    if sel.elements:
        parts.append(" ".join(sel.elements))
    if sel.z_min is not None or sel.z_max is not None:
        parts.append(f"z {_fmt(sel.z_min)}〜{_fmt(sel.z_max)} Å")
    if sel.indices:
        parts.append(L(f"番号 {len(sel.indices)} 個", f"{len(sel.indices)} indices"))
    return ", ".join(parts)


class SupercellEditor(StepEditor):
    op = "supercell"

    def build(self):
        self.mode = QComboBox(); self.mode.addItem(L("各軸の繰り返し", "Repeat along each axis"), "repeat"); self.mode.addItem(L("変換行列 (3×3)", "Transformation matrix (3×3)"), "matrix")
        limit_combo(self.mode)
        self.rep = [spin(1, 50, 1, 64) for _ in range(3)]
        self.rep_row = hrow(self.rep[0], "×", self.rep[1], "×", self.rep[2])
        self.mat = [[spin(-20, 20, int(i == j), 64) for j in range(3)] for i in range(3)]
        self.mat_box = QWidget(); self.mat_box.setObjectName("rowbox"); g = QGridLayout(self.mat_box); g.setContentsMargins(0, 0, 0, 0)
        for i in range(3):
            for j in range(3):
                g.addWidget(self.mat[i][j], i, j)
        add_row(self.form, "超格子の取り方", self.mode)
        add_row(self.form, "繰り返し (a × b × c)", self.rep_row)
        add_row(self.form, "変換行列", self.mat_box)
        self.fit_hint = hint(L("溶液の成分が周期境界で重ならないよう、a・b 方向の繰り返し数を必要に応じて増やします。表示値は繰り返し数の下限です",
                               "The repeats along a and b are increased as needed so that each solution component clears its periodic images. The displayed values are minimum repeats."))
        self.form.addRow(label(""), self.fit_hint)
        self.form.addRow(label(""), hint(L("変換行列の各行が新しい格子ベクトル (元の a, b, c の整数倍の和) です", "each matrix row is a new lattice vector (an integer combination of a, b, c)")))
        self.mode.currentIndexChanged.connect(self._show)

    def _show(self, *_):
        m = self.mode.currentData() == "matrix"
        self.form.setRowVisible(self.rep_row, not m); self.form.setRowVisible(self.mat_box, m)

    def load(self, s: Supercell):
        self.mode.setCurrentIndex(1 if s.matrix is not None else 0)
        for w, v in zip(self.rep, s.repeat or (1, 1, 1)):
            w.setValue(int(v))
        m = s.matrix or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        for i in range(3):
            for j in range(3):
                self.mat[i][j].setValue(int(m[i][j]))
        self.fit_hint.setVisible(bool(s.fit_components))
        self._show()

    def values(self):
        if self.mode.currentData() == "matrix":
            return {"repeat": None, "matrix": [[self.mat[i][j].value() for j in range(3)] for i in range(3)]}
        return {"repeat": tuple(w.value() for w in self.rep), "matrix": None}

    def summary(self, s: Supercell):
        repeat = " × ".join(str(x) for x in s.repeat) if s.repeat else L("変換行列", "matrix")
        return L(f"溶液に合わせて自動調整 (下限 {repeat})", f"auto-sized for solution (minimum {repeat})") if s.fit_components else repeat


class SlabEditor(StepEditor):
    op = "slab"
    terminations_needed = Signal()

    def build(self):
        self.hkl = [spin(-9, 9, 0, 60) for _ in range(3)]
        self.layers = spin(1, 50, 3)
        self.vacuum = num(0.0, 200.0, 10.0, 1.0)
        self.term = QComboBox(); limit_combo(self.term)
        self.term.setToolTip(L("切った面の 1 周期の中の原子面を上から数えた番号と、その面の組成。どの面をいちばん上にするかを選びます (組成は変わりません)",
                               "the atomic planes of one period, counted from the top, with their composition; choose which one is the top surface (the composition does not change)"))
        add_row(self.form, "ミラー指数 (h k l)", hrow(*self.hkl))
        add_row(self.form, "層の数 (周期)", self.layers)
        add_row(self.form, "真空 (片側) [Å]", self.vacuum)
        add_row(self.form, "終端 (最上面)", self.term)
        self.form.addRow(label(""), hint(L("ミラー指数は、いまのセルの格子ベクトルに対する指数です。立方晶の指数で切るなら、バルクは「立方晶セル」にしてください",
                                           "Miller indices refer to the current cell vectors. For cubic indices, tick \"Cubic cell\" for the bulk")))
        self._term_n = 0

    def connect_all(self):
        super().connect_all()
        for w in self.hkl:
            w.valueChanged.connect(lambda *_: self.terminations_needed.emit())

    def miller(self) -> tuple[int, int, int]:
        return tuple(w.value() for w in self.hkl)

    def set_terminations(self, names: list[str] | None) -> None:
        cur = self.term.currentData()
        cur = self._model.termination if cur is None else cur
        self.term.blockSignals(True)
        self.term.clear()
        if names:
            for i, f in enumerate(names):
                self.term.addItem(L(f"{i}: 最上面が {pretty_formula(f)}", f"{i}: {pretty_formula(f)} on top"), i)
        else:
            for i in range(max(6, int(cur) + 1)):
                self.term.addItem(L(f"{i} (組成は「作る」のあとに分かります)", f"{i} (composition known after Build)") if i == 0 else str(i), i)
        k = self.term.findData(int(cur))
        if k < 0:
            self.term.addItem(str(cur), int(cur)); k = self.term.count() - 1
        self.term.setCurrentIndex(k)
        self.term.blockSignals(False)

    def load(self, s: Slab):
        for w, v in zip(self.hkl, s.miller):
            w.setValue(int(v))
        self.layers.setValue(s.layers); self.vacuum.setValue(s.vacuum)
        if self.term.count() == 0 or self.term.findData(s.termination) < 0:
            self.set_terminations(None)
        self.term.setCurrentIndex(max(0, self.term.findData(s.termination)))

    def values(self):
        t = self.term.currentData()
        return {"miller": self.miller(), "layers": self.layers.value(), "vacuum": self.vacuum.value(), "termination": int(t or 0)}

    def summary(self, s: Slab):
        return L(f"({' '.join(map(str, s.miller))})、{s.layers} 層", f"({' '.join(map(str, s.miller))}), {s.layers} layer{'s' if s.layers != 1 else ''}")


class VacuumEditor(StepEditor):
    op = "vacuum"

    def build(self):
        self.axis = QComboBox()
        for i, name in enumerate(("a", "b", L("c (スラブなら面に垂直)", "c (normal to a slab)"))):
            self.axis.addItem(name, i)
        limit_combo(self.axis)
        self.thickness = num(0.0, 500.0, 10.0, 1.0)
        add_row(self.form, "軸", self.axis)
        add_row(self.form, "周期の像との隙間 [Å]", self.thickness)

    def load(self, s: Vacuum):
        self.axis.setCurrentIndex(s.axis); self.thickness.setValue(s.thickness)

    def values(self):
        return {"axis": self.axis.currentData(), "thickness": self.thickness.value()}

    def summary(self, s: Vacuum):
        return f"{'abc'[s.axis]}, {s.thickness:g} Å"


class BoxEditor(StepEditor):
    op = "box"

    def build(self):
        self.padding = num(0.0, 200.0, 5.0, 1.0)
        self.maxm = spin(1, 64, 8)
        add_row(self.form, "余白 [Å]", self.padding)
        add_row(self.form, "倍数の上限", self.maxm)
        self.form.addRow(label(""), hint(L("分子やナノ粒子は、端から余白ずつ離した直方体の箱に入れます。周期の構造は、直交するセルに取り直します (元のセルの「倍数の上限」倍まで)",
                                           "a molecule or nanoparticle is put in a rectangular box with this padding on every side; a periodic structure is re-cut into an orthogonal cell (up to the multiple limit)")))

    def load(self, s: Box):
        self.padding.setValue(s.padding); self.maxm.setValue(s.max_multiple)

    def values(self):
        return {"padding": self.padding.value(), "max_multiple": self.maxm.value()}

    def summary(self, s: Box):
        return L(f"余白 {s.padding:g} Å", f"padding {s.padding:g} Å")


class _PickEditor(StepEditor):
    def build(self):
        self.sel = SelectionRows(self.form)
        self.mode = QComboBox(); self.mode.addItem(L("個数", "Count"), "count"); self.mode.addItem(L("割合 (0〜1)", "Fraction (0-1)"), "fraction")
        limit_combo(self.mode)
        self.count = spin(1, 100000, 1, 96)
        self.fraction = num(0.0, 1.0, 0.1, 0.05)
        self.seed = spin(0, 999999, 0, 96)
        add_row(self.form, "選ぶ数", hrow(self.mode, self.count, self.fraction))
        add_row(self.form, "乱数の種", self.seed)
        self.mode.currentIndexChanged.connect(self._show)

    def _show(self, *_):
        f = self.mode.currentData() == "fraction"
        self.count.setVisible(not f); self.fraction.setVisible(f)

    def load(self, s):
        self.sel.load(s.where)
        self.mode.setCurrentIndex(1 if s.fraction is not None else 0)
        if s.count is not None:
            self.count.setValue(s.count)
        if s.fraction is not None:
            self.fraction.setValue(s.fraction)
        self.seed.setValue(s.seed)
        self._show()

    def values(self):
        f = self.mode.currentData() == "fraction"
        return {"where": self.sel.value(self._model.where), "count": None if f else self.count.value(),
                "fraction": self.fraction.value() if f else None, "seed": self.seed.value()}

    def _amount(self, s) -> str:
        what = describe_selection(s.where) or L("全原子", "all atoms")
        n = f"{s.count}" if s.count is not None else f"{s.fraction:g}"
        return L(f"{what} から {n}", f"{n} of {what}")


class RemoveEditor(_PickEditor):
    op = "remove"

    def summary(self, s: Remove):
        return self._amount(s)


class SubstituteEditor(_PickEditor):
    op = "substitute"

    def build(self):
        super().build()
        self.to = QComboBox(); self.to.addItems(ELEMENTS); limit_combo(self.to)
        add_row(self.form, "置き換え先の元素", self.to)

    def load(self, s: Substitute):
        super().load(s); self.to.setCurrentText(s.to)

    def values(self):
        return {**super().values(), "to": self.to.currentText()}

    def summary(self, s: Substitute):
        return f"{self._amount(s)} → {s.to}"


class AdsorbEditor(StepEditor):
    op = "adsorb"

    def build(self):
        self.kind = QComboBox()
        for k, v in KINDS.items():
            self.kind.addItem(icons.source_icon(k), {"preset": L("プリセット", "Preset"), "smiles": "SMILES", "file": L("ファイル", "File")}[k], k)
        limit_combo(self.kind)
        self.ref = RefCell("preset", "H2O")
        self.ref.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.place = QComboBox()
        self.place.addItem(L("吸着サイトの名前 (ASE の面)", "Site name (ASE surface)"), "site")
        self.place.addItem(L("原子の真上", "Above an atom"), "above_atom")
        self.place.addItem(L("xy 座標", "xy position"), "xy")
        limit_combo(self.place)
        self.site = QComboBox(); self.site.setEditable(True); self.site.addItems(["ontop", "bridge", "fcc", "hcp", "hollow", "shortbridge", "longbridge"])
        self.site.setToolTip(L("土台が「スラブ」(ASE の面) のときだけ使えます。面ごとの吸着サイトの名前は ASE と同じ", "only for a Slab base (ASE surface); the site names are ASE's"))
        self.atom = spin(1, 1000000, 1, 110)
        self.x = num(-1000.0, 1000.0, 0.0, 0.5); self.y = num(-1000.0, 1000.0, 0.0, 0.5)
        self.height = num(0.0, 50.0, 2.0, 0.1)
        self.down = spin(1, 10000, 1)
        add_row(self.form, "分子", hrow(self.kind, self.ref, stretch=False))
        add_row(self.form, "置き場所", self.place)
        add_row(self.form, "吸着サイトの名前", self.site)
        add_row(self.form, "原子の番号 (1 始まり)", self.atom)
        add_row(self.form, "xy 座標 [Å]", hrow(self.x, self.y))
        add_row(self.form, "高さ [Å]", self.height)
        add_row(self.form, "下に向ける原子 (1 始まり)", self.down)
        self.place.currentIndexChanged.connect(self._show)
        self.kind.currentIndexChanged.connect(lambda *_: (self.ref.set_kind(self.kind.currentData()), self._emit()))
        self.ref.changed.connect(self._emit)
        self.ref.draw_requested.connect(lambda: self.draw_requested.emit(self, 0))

    def _show(self, *_):
        p = self.place.currentData()
        self.form.setRowVisible(self.site, p == "site"); self.form.setRowVisible(self.atom, p == "above_atom")
        self.form.setRowVisible(self.x.parentWidget(), p == "xy")

    def row_smiles(self, _row: int) -> str:
        return self.ref.text().strip() if self.ref.kind() == "smiles" else ""

    def set_row_smiles(self, _row: int, smiles: str) -> None:
        self.kind.setCurrentIndex(list(KINDS).index("smiles")); self.ref.set_kind("smiles"); self.ref.setText(smiles); self._emit()

    def load(self, s: Adsorb):
        self.kind.blockSignals(True); self.kind.setCurrentIndex(list(KINDS).index(s.molecule.kind)); self.kind.blockSignals(False)
        self.ref.set_kind(s.molecule.kind); self.ref.setText(s.molecule.ref)
        place = "site" if s.site is not None else "above_atom" if s.above_atom is not None else "xy"
        self.place.setCurrentIndex(self.place.findData(place))
        if s.site is not None:
            self.site.setCurrentText(s.site)
        if s.above_atom is not None:
            self.atom.setValue(s.above_atom + 1)
        if s.xy is not None:
            self.x.setValue(s.xy[0]); self.y.setValue(s.xy[1])
        self.height.setValue(s.height); self.down.setValue(s.down_atom + 1)
        self._show()

    def values(self):
        p = self.place.currentData()
        return {"molecule": {"kind": self.ref.kind(), "ref": self.ref.text().strip()},
                "site": self.site.currentText().strip() if p == "site" else None,
                "above_atom": self.atom.value() - 1 if p == "above_atom" else None,
                "xy": (self.x.value(), self.y.value()) if p == "xy" else None,
                "height": self.height.value(), "down_atom": self.down.value() - 1}

    def summary(self, s: Adsorb):
        where = s.site if s.site is not None else L(f"原子 {s.above_atom + 1} の上", f"above atom {s.above_atom + 1}") if s.above_atom is not None else f"xy ({s.xy[0]:g}, {s.xy[1]:g})"
        mol = pretty_formula(s.molecule.ref) if s.molecule.kind == "preset" else s.molecule.ref
        return f"{mol}, {where}, {s.height:g} Å"


class _ComponentsEditor(StepEditor):

    auto_count = False

    def build(self):
        self.mix = MixtureEditor(cell_row=False, example=False, auto_count=self.auto_count)
        self.form.addRow(self.mix)
        self.mix.changed.connect(self._emit)
        self.mix.draw_requested.connect(lambda row: self.draw_requested.emit(self.mix, row))

    def set_draw_enabled(self, on: bool) -> None:
        self.mix.set_draw_enabled(on)

    def load_components(self, s) -> None:
        from adit.mixture import MixtureSpec
        self.mix.blockSignals(True)
        try:
            self.mix.set_spec(MixtureSpec(components=list(s.components), min_distance=s.min_distance, seed=s.seed))
        finally:
            self.mix.blockSignals(False)

    def component_values(self) -> dict:
        m = self.mix.spec()
        return {"components": [c.model_dump() for c in m.components], "min_distance": m.min_distance, "seed": m.seed}

    def comps_summary(self, s) -> str:
        return ", ".join(f"{c.name()} × {c.count if c.count else L('自動', 'auto')}" for c in s.components)


class SolventLayerEditor(_ComponentsEditor):
    op = "solvent_layer"
    area_provider = None

    def build(self):
        super().build()
        self.mode = QComboBox(); self.mode.addItem(L("密度から [g/cm³]", "From density [g/cm³]"), "density"); self.mode.addItem(L("厚みを指定 [Å]", "Thickness [Å]"), "thickness")
        limit_combo(self.mode)
        self.density = num(0.01, 30.0, 1.0, 0.05)
        self.density.setToolTip(L("成分の総質量をこの密度で割った体積が、溶液の帯になります (スラブの断面は変えずに厚みが決まります)",
                                  "the band volume is the total mass of the components divided by this density (the slab cross-section is kept; the thickness follows)"))
        self.thickness = num(0.1, 1000.0, 20.0, 1.0)
        self.thickness.setToolTip(L("分子の中心を置く帯の厚み。詰めすぎると入り切らないので、はじめは密度から決める方が確実です",
                                    "thickness of the band where molecule centers are placed; too thin and the molecules will not fit, so starting from a density is safer"))
        self.gap = num(0.0, 50.0, 2.0, 0.1)
        self.vacuum = num(0.0, 500.0, 0.0, 1.0)
        self.conc = num(0.0, 50.0, 1.0, 0.1)
        self.btn_conc = QPushButton(L("選んだ行に入れる", "Fill selected row"))
        self.btn_conc.setToolTip(L("溶液の帯の体積と、この濃度 [mol/L] から個数を出し、表で選んだ行 (イオンなど) に入れます。塩なら陽イオンと陰イオンの行それぞれに",
                                   "computes a count from the band volume and this concentration [mol/L] and puts it in the selected table row (for a salt, do it for the cation and the anion rows)"))
        add_row(self.form, "厚みの決め方", hrow(self.mode, self.density, self.thickness))
        add_row(self.form, "隙間 [Å]", self.gap)
        add_row(self.form, "真空 [Å] (0 なら両側が溶液)", self.vacuum)
        add_row(self.form, "濃度から個数", hrow(self.conc, "mol/L", self.btn_conc))
        self.mode.currentIndexChanged.connect(self._show)
        self.btn_conc.clicked.connect(self.fill_count_from_concentration)

    def _show(self, *_):
        d = self.mode.currentData() == "density"
        self.density.setVisible(d); self.thickness.setVisible(not d)

    def load(self, s: SolventLayer):
        self.load_components(s)
        self.mode.setCurrentIndex(0 if s.thickness is None else 1)
        self.density.setValue(s.density_g_cm3)
        if s.thickness is not None:
            self.thickness.setValue(s.thickness)
        self.gap.setValue(s.gap); self.vacuum.setValue(s.vacuum)
        self._show()

    def values(self):
        d = self.mode.currentData() == "density"
        return {**self.component_values(), "density_g_cm3": self.density.value(), "thickness": None if d else self.thickness.value(),
                "gap": self.gap.value(), "vacuum": self.vacuum.value()}

    def summary(self, s: SolventLayer):
        return self.comps_summary(s)

    def band_volume(self, counts: list[int]) -> float | None:
        from adit.mixture import AMU_TO_G, _one_molecule
        m = self.mix.spec()
        if self.mode.currentData() == "density":
            mass = sum(float(_one_molecule(c).get_masses().sum()) * n for c, n in zip(m.components, counts)) * AMU_TO_G
            return mass / self.density.value() * 1e24
        area = self.area_provider() if callable(self.area_provider) else None
        return None if not area else area * self.thickness.value()

    def fill_count_from_concentration(self) -> str:
        from adit.mixture import salt_count
        r = self.mix.current_row()
        m = self.mix.spec()
        if not 0 <= r < len(m.components):
            msg = L("先に表で、個数を決める行 (イオンなど) を選んでください", "select the table row (such as an ion) whose count you want first")
            self.mix.info.setText(msg); return msg
        counts = [c.count for c in m.components]
        try:
            n = counts[r]
            for _ in range(30):
                v = self.band_volume(counts)
                if v is None:
                    msg = L("断面積がまだ分かりません。先に「作る」を押すか、厚みを密度から決めてください",
                            "the cross-section is not known yet; press Build first, or choose the thickness from a density")
                    self.mix.info.setText(msg); return msg
                new, _c = salt_count(self.conc.value(), v)
                if new == n:
                    break
                n = counts[r] = new
            if n < 1:
                msg = L(f"帯の体積 {v:.0f} Å³ では {self.conc.value():g} mol/L は 1 個未満です。成分を増やすか濃度を上げてください",
                        f"{self.conc.value():g} mol/L is less than one ion in a band of {v:.0f} Å³; add more components or raise the concentration")
                self.mix.info.setText(msg); return msg
            actual = salt_count(self.conc.value(), v)[1]
        except Exception as ex:
            msg = str(ex); self.mix.info.setText(msg); return msg
        self.mix.table.cellWidget(r, 2).setValue(n)
        name = pretty_formula(m.components[r].name())
        msg = L(f"帯の体積 {v:.0f} Å³ に {name} を {n} 個 (実際の濃度 {actual:.2f} mol/L)", f"{n} × {name} in a band of {v:.0f} Å³ (actual {actual:.2f} mol/L)")
        if self.mode.currentData() == "density":
            msg += L("。ほかの行の個数を変えたら、もう一度押してください", ". If you change another row's count, press again")
        self.mix.info.setText(msg)
        return msg


class SolvateEditor(_ComponentsEditor):
    op = "solvate"
    auto_count = True

    def build(self):
        super().build()
        self.padding = num(0.0, 200.0, 8.0, 1.0)
        self.density = num(0.01, 30.0, 1.0, 0.05)
        add_row(self.form, "溶質の周りの余白 [Å]", self.padding)
        add_row(self.form, "密度 [g/cm³]", self.density)
        self.form.addRow(label(""), hint(L("いまの構造 (溶質) は動かしません。個数を「自動」にした成分 1 つを、箱全体がこの密度になるまで入れます",
                                           "the current structure (solute) is not moved; one component set to \"auto\" is added until the whole box reaches this density")))

    def load(self, s: Solvate):
        self.load_components(s); self.padding.setValue(s.padding); self.density.setValue(s.density_g_cm3)

    def values(self):
        return {**self.component_values(), "padding": self.padding.value(), "density_g_cm3": self.density.value()}

    def summary(self, s: Solvate):
        return self.comps_summary(s)


class FixEditor(StepEditor):
    op = "fix"

    def build(self):
        self.layers = spin(0, 200, 1)
        self.layers.setSpecialValueText(L("使わない", "not used"))
        add_row(self.form, "下から数えた層", self.layers)
        self.sel = SelectionRows(self.form)
        self.form.addRow(label(""), hint(L("条件はすべて「かつ」。何も指定しないと全原子を固定します。固定の手順があるときは、共通設定の「固定原子」の欄は使いません",
                                           "all conditions must hold; with none, every atom is fixed. When there is a fix step, the \"Fixed atoms\" field below is not used")))

    def load(self, s: Fix):
        self.layers.setValue(s.bottom_layers); self.sel.load(s.where)

    def values(self):
        return {"bottom_layers": self.layers.value(), "where": self.sel.value(self._model.where)}

    def summary(self, s: Fix):
        parts = [L(f"下から {s.bottom_layers} 層", f"bottom {s.bottom_layers} layer{'s' if s.bottom_layers != 1 else ''}")] if s.bottom_layers else []
        d = describe_selection(s.where)
        if d:
            parts.append(d)
        return ", ".join(parts) or L("全原子", "all atoms")


EDITORS = {c.op: c for c in (SupercellEditor, SlabEditor, VacuumEditor, BoxEditor, AdsorbEditor, RemoveEditor, SubstituteEditor,
                             SolventLayerEditor, SolvateEditor, FixEditor)}

_STEP_NO = re.compile(r"(?:手順|step) (\d+) \(")


def failed_step(message: str) -> int | None:
    m = _STEP_NO.search(message or "")
    k = int(m.group(1)) if m else 0
    return k if k >= 1 else None


class RecipeEditor(QWidget):

    changed = Signal()
    build_requested = Signal()
    cancel_requested = Signal()
    draw_requested = Signal(object, int)
    terminations_needed = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("rowbox")
        self.editors: list[StepEditor] = []
        self.base_is_slab = lambda: False
        self.area_before = lambda index: None
        self._draw_enabled = True

        self.list = QListWidget(); self.list.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.list.setMinimumHeight(64); self.list.setMaximumHeight(150)
        self.add_kind = QComboBox()
        for op in ADD_ORDER:
            self.add_kind.addItem(op_label(op), op)
        self.add_kind.insertSeparator(self.add_kind.count())
        self.add_kind.addItem(icons.icon("slab"), L("電極と電解質の界面 (面 + 断面の自動調整 + 溶液 + 固定)", "Electrode–electrolyte interface (surface cut + auto-sized cross-section + solution + fixed layer)"), INTERFACE)
        limit_combo(self.add_kind)
        self.add_kind.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.add_kind.setMinimumContentsLength(12)
        self.add_kind.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.add_kind.view().setMinimumWidth(self.add_kind.view().sizeHintForColumn(0) + 32)
        self.add_kind.setToolTip(L("足す手順の種類。「電極と電解質の界面」は、面で切る・断面を溶液に合わせて広げる・溶液の層を置く・下の層を固定する、の各手順をまとめて追加します。土台がすでにスラブなら、面で切る手順は省きます",
                                   "the kind of step to add. \"Electrode–electrolyte interface\" adds steps to cut a slab, widen its cross-section to fit the solution, add the electrolyte layer, and fix the bottom layer. The slab step is skipped when the base is already a slab"))
        self.btn_add = QPushButton(L("追加", "Add"))
        self.btn_up = QPushButton(L("上へ", "Up")); self.btn_down = QPushButton(L("下へ", "Down")); self.btn_del = QPushButton(L("削除", "Remove"))
        for b in (self.btn_up, self.btn_down, self.btn_del):
            b.setObjectName("cell_button")
        self.empty = hint(L("手順を足すと、土台の構造を面で切る・溶液を足す・原子を抜くなどの加工ができます。手順があるときは「作る」を押して構造を作ります",
                            "add steps to cut the base structure along a plane, add a solution layer, remove atoms and so on. With steps, press Build to make the structure"))
        self.stack = QStackedWidget()
        self.step_title = QLabel(""); self.step_title.setObjectName("subtitle")
        self.btn_build = QPushButton(icons.icon("generate", 16), L("作る", "Build")); self.btn_build.setObjectName("primary")
        self.btn_build.setToolTip(L("組み立て手順を順に実行して構造を作ります (溶液の詰め込みやポリマーは数秒〜数十秒かかります)",
                                    "runs the steps in order to make the structure (packing a solution or a polymer takes seconds to tens of seconds)"))
        from adit.gui.progress import ProgressStrip
        self.progress = ProgressStrip(); self.btn_cancel = self.progress.btn_cancel
        self.status = QLabel(""); self.status.setObjectName("hint"); self.status.setWordWrap(True)
        self.log = QLabel(""); self.log.setObjectName("hint"); self.log.setWordWrap(True)
        self.log.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        add_row_w = hrow(self.add_kind, self.btn_add, stretch=False)
        edit_row = hrow(self.btn_up, self.btn_down, self.btn_del)
        self.build_row = hrow(self.btn_build)
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(6)
        for w in (add_row_w, self.list, edit_row, self.empty, self.step_title, self.stack, self.build_row, self.progress, self.status, self.log):
            lay.addWidget(w)
        self.btn_add.clicked.connect(self._on_add)
        self.btn_up.clicked.connect(lambda: self.move(-1)); self.btn_down.clicked.connect(lambda: self.move(1))
        self.btn_del.clicked.connect(self.remove_current)
        self.list.currentRowChanged.connect(self._on_row)
        self.btn_build.clicked.connect(self.build_requested.emit); self.progress.cancel_requested.connect(self.cancel_requested.emit)
        self._refresh()

    def count(self) -> int:
        return len(self.editors)

    def steps(self) -> list:
        out = []
        for k, ed in enumerate(self.editors, start=1):
            try:
                out.append(ed.step())
            except ValueError as ex:
                msg = str(ex)
                errs = getattr(ex, "errors", None)
                if callable(errs):
                    e = errs()[0] if errs() else {}
                    msg = str(e.get("msg", msg)).removeprefix("Value error, ")
                raise ValueError(L(f"手順 {k} ({op_label(ed.op)}): {msg}", f"step {k} ({op_label(ed.op)}): {msg}")) from ex
        return out

    def has_op(self, op: str) -> bool:
        return any(ed.op == op for ed in self.editors)

    def add_step(self, step, *, select: bool = True, emit: bool = True) -> StepEditor:
        ed = EDITORS[step.op](step)
        translate_widgets(ed)
        if isinstance(ed, _ComponentsEditor):
            ed.set_draw_enabled(self._draw_enabled)
        if isinstance(ed, SolventLayerEditor):
            ed.area_provider = lambda ed=ed: self.area_before(self.editors.index(ed))
        if isinstance(ed, AdsorbEditor):
            ed.ref.btn_draw.setEnabled(self._draw_enabled)
        ed.changed.connect(lambda ed=ed: self._on_editor_changed(ed))
        ed.draw_requested.connect(self.draw_requested.emit)
        if isinstance(ed, SlabEditor):
            ed.terminations_needed.connect(self.terminations_needed.emit)
        self.editors.append(ed); self.stack.addWidget(ed)
        self.list.addItem("")
        self._refresh()
        if select:
            self.list.setCurrentRow(len(self.editors) - 1)
        if emit:
            self.changed.emit()
        return ed

    def add_interface(self) -> None:
        for step in interface_steps(self.base_is_slab(), ELECTROLYTE_EXAMPLE):
            self.add_step(step, select=False, emit=False)
        self.list.setCurrentRow(self.editors.index(next(e for e in self.editors[::-1] if e.op == "solvent_layer")))
        self.status.setText(L("溶液の成分は例です。個数と密度 (または厚み) を決めてから「作る」を押してください。スラブの断面は、成分が周期境界で重ならない広さまで自動的に広がります。イオンの個数は「濃度から個数」でも決められます",
                              "the solution components are examples. Set their counts and the density (or thickness), then press Build. The slab cross-section is enlarged automatically until each component clears its periodic images. You can also set ion counts with \"Count from concentration\""))
        self.changed.emit()

    def set_steps(self, steps: list) -> None:
        self.clear()
        for s in steps:
            self.add_step(s, select=False, emit=False)
        if self.editors:
            self.list.setCurrentRow(0)
        self._refresh()

    def clear(self) -> None:
        for ed in self.editors:
            self.stack.removeWidget(ed); ed.deleteLater()
        self.editors.clear(); self.list.clear(); self._refresh()

    def move(self, delta: int) -> None:
        r = self.list.currentRow(); t = r + delta
        if not (0 <= r < len(self.editors) and 0 <= t < len(self.editors)):
            return
        self.editors[r], self.editors[t] = self.editors[t], self.editors[r]
        self._refresh(); self.list.setCurrentRow(t); self.changed.emit()

    def remove_current(self) -> None:
        r = self.list.currentRow()
        if not 0 <= r < len(self.editors):
            return
        ed = self.editors.pop(r); self.stack.removeWidget(ed); ed.deleteLater()
        self.list.takeItem(r)
        self._refresh()
        if self.editors:
            self.list.setCurrentRow(min(r, len(self.editors) - 1))
        self.changed.emit()

    def select_step(self, k: int) -> None:
        if 1 <= k <= len(self.editors):
            self.list.setCurrentRow(k - 1)

    def set_draw_enabled(self, on: bool) -> None:
        self._draw_enabled = on
        for ed in self.editors:
            if isinstance(ed, _ComponentsEditor):
                ed.set_draw_enabled(on)
            if isinstance(ed, AdsorbEditor):
                ed.ref.btn_draw.setEnabled(on)

    def set_buildable(self, on: bool) -> None:
        self.build_row.setVisible(on)

    def set_building(self, on: bool) -> None:
        self.btn_build.setEnabled(not on)
        if on:
            self.status.setObjectName("hint"); self._restyle(self.status); self.status.setText("")
            self.progress.begin(L("作っています…", "Building…"))
        else:
            self.progress.end()

    def show_stale(self, text: str) -> None:
        self.status.setObjectName("hint"); self._restyle(self.status); self.status.setText(text)

    def show_result(self, lines: list[str], error: str = "") -> None:
        self.set_building(False)
        self.log.setText("\n".join(lines))
        self.status.setObjectName("status_ng" if error else "status_ok"); self._restyle(self.status)
        if error:
            self.status.setText(error)
            k = failed_step(error)
            if k is not None:
                self.select_step(k)
        else:
            self.status.setText(L("できました", "Done"))

    def show_cancelled(self) -> None:
        self.set_building(False)
        self.status.setObjectName("hint"); self._restyle(self.status)
        self.status.setText(L("中止しました (途中の計算は裏で終わるまで続き、結果は捨てます)", "cancelled (the running build finishes in the background and its result is discarded)"))

    @staticmethod
    def _restyle(w: QWidget) -> None:
        w.style().unpolish(w); w.style().polish(w)

    def _on_add(self) -> None:
        op = self.add_kind.currentData()
        if op == INTERFACE:
            self.add_interface()
        elif op:
            self.add_step(default_step(op))

    def _on_editor_changed(self, ed: StepEditor) -> None:
        self._refresh_row(self.editors.index(ed)); self.changed.emit()

    def _on_row(self, r: int) -> None:
        if 0 <= r < len(self.editors):
            for i, ed in enumerate(self.editors):
                pol = QSizePolicy.Policy.Preferred if i == r else QSizePolicy.Policy.Ignored
                ed.setSizePolicy(pol, pol)
            self.stack.setCurrentWidget(self.editors[r]); self.stack.updateGeometry()
            self.step_title.setText(L(f"手順 {r + 1}: {op_label(self.editors[r].op)}", f"Step {r + 1}: {op_label(self.editors[r].op)}"))

    def _refresh_row(self, r: int) -> None:
        ed = self.editors[r]
        try:
            s = ed.summary(ed.step())
        except ValueError:
            s = L("欄に読めない値があります", "a field cannot be read")
        item = self.list.item(r)
        if item is not None:
            item.setText(f"{r + 1}. {op_label(ed.op)}" + (f"  —  {s}" if s else ""))

    def _refresh(self) -> None:
        while self.list.count() < len(self.editors):
            self.list.addItem("")
        for r in range(len(self.editors)):
            self._refresh_row(r)
        for i, ed in enumerate(self.editors):
            if self.stack.indexOf(ed) != i:
                self.stack.removeWidget(ed); self.stack.insertWidget(i, ed)
        has = bool(self.editors)
        for w in (self.list, self.btn_up, self.btn_down, self.btn_del, self.stack, self.step_title):
            w.setVisible(has)
        self.btn_up.parentWidget().setVisible(has)
        self.empty.setVisible(not has)
        r = self.list.currentRow()
        if has:
            self._on_row(r if 0 <= r < len(self.editors) else 0)
