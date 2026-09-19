
from __future__ import annotations

from PySide6.QtCore import QEvent, Signal
from PySide6.QtWidgets import (QAbstractItemView, QComboBox, QFileDialog, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QPushButton, QSizePolicy,
                               QSpinBox, QStackedWidget, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from adit.gui import icons
from adit.gui.i18n import tr
from adit.gui.widgets import SciDoubleSpinBox, limit_combo, narrow
from adit.lang import L
from adit.mixture import Component, MixtureSpec, concentration_mol_per_l
from adit.structure import preset_names, pretty_formula

KINDS = {"preset": "プリセット", "smiles": "SMILES", "file": "ファイル"}
COL_KIND, COL_REF, COL_COUNT, COL_CHARGE, COL_LABEL = range(5)


class RefCell(QWidget):

    changed = Signal()
    draw_requested = Signal()

    def __init__(self, kind: str, ref: str, parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("rowbox")
        self.preset = QComboBox()
        for name in preset_names():
            self.preset.addItem(pretty_formula(name), name)
        limit_combo(self.preset)
        self.preset.view().setMinimumWidth(self.preset.view().sizeHintForColumn(0) + 32)
        self.preset.setToolTip(tr("ASE に収録された分子 (G2 集など) から選びます"))
        self.smiles = QLineEdit(); self.smiles.setPlaceholderText(tr("例 [Na+]、CCO"))
        self.smiles.setToolTip(tr("SMILES を入力して Enter。「Draw」で描くこともできます"))
        self.btn_draw = QPushButton("Draw"); self.btn_draw.setToolTip(tr("この行の分子を描いて SMILES にします (RDKit が要ります)"))
        self.file = QLineEdit(); self.file.setPlaceholderText(tr("構造ファイルのパス"))
        self.btn_browse = QPushButton(tr("参照…"))
        for b in (self.btn_draw, self.btn_browse):
            b.setObjectName("cell_button")
        for e in (self.smiles, self.file):
            e.setMinimumWidth(36)
        self.stack = QStackedWidget()
        self.stack.addWidget(self.preset)
        self.stack.addWidget(self._row(self.smiles, self.btn_draw))
        self.stack.addWidget(self._row(self.file, self.btn_browse))
        lay = QHBoxLayout(self); lay.setContentsMargins(2, 2, 2, 2); lay.addWidget(self.stack)
        self.preset.setCurrentIndex(max(0, self.preset.findData("H2O")))
        self.set_kind(kind); self.setText(ref)
        self.preset.activated.connect(lambda *_: self.changed.emit())
        self.smiles.editingFinished.connect(self._edited)
        self.file.editingFinished.connect(self._edited)
        self.btn_draw.clicked.connect(self.draw_requested.emit)
        self.btn_browse.clicked.connect(self._browse)
        self._last = self.text()

    @staticmethod
    def _row(edit: QLineEdit, button: QPushButton) -> QWidget:
        w = QWidget(); w.setObjectName("rowbox"); h = QHBoxLayout(w); h.setContentsMargins(0, 0, 0, 0); h.setSpacing(4)
        button.ensurePolished(); button.setMinimumWidth(button.sizeHint().width())
        h.addWidget(edit, 1); h.addWidget(button)
        return w

    def kind(self) -> str:
        return list(KINDS)[self.stack.currentIndex()]

    def set_kind(self, kind: str) -> None:
        self.stack.setCurrentIndex(list(KINDS).index(kind) if kind in KINDS else 0)
        self._last = self.text()

    def text(self) -> str:
        k = self.kind()
        if k == "preset":
            return self.preset.currentData() or ""
        return (self.smiles if k == "smiles" else self.file).text()

    def setText(self, ref: str) -> None:  # noqa: N802 
        k = self.kind()
        if k == "preset":
            i = self.preset.findData(ref)
            if i < 0 and ref:
                self.preset.addItem(pretty_formula(ref), ref); i = self.preset.count() - 1
            if i >= 0:
                self.preset.setCurrentIndex(i)
        else:
            (self.smiles if k == "smiles" else self.file).setText(ref)
        self._last = self.text()

    def _edited(self) -> None:
        if self.text() != self._last:
            self._last = self.text(); self.changed.emit()

    def _browse(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, tr("構造ファイル"), "",
                                                   L("構造 (*.xyz *.cif *.pdb *.gen *.vasp POSCAR CONTCAR);;すべて (*)",
                                                     "Structures (*.xyz *.cif *.pdb *.gen *.vasp POSCAR CONTCAR);;All files (*)"))
        if path:
            self.file.setText(path); self._edited()


class MixtureEditor(QWidget):
    changed = Signal()
    draw_requested = Signal(int)

    def __init__(self, parent: QWidget | None = None, *, cell_row: bool = True, example: bool = True, auto_count: bool = False):
        super().__init__(parent)
        self.setObjectName("rowbox")
        self._draw_enabled = True
        self._auto_count = auto_count
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels([tr(t) for t in ("成分の種類", "指定", "個数", "電荷", "表示名")])
        self.table.horizontalHeaderItem(COL_REF).setToolTip(tr("プリセットは一覧から選び、SMILES は入力するか「Draw」で描き、ファイルは「参照…」で選びます"))
        self.table.horizontalHeaderItem(COL_CHARGE).setToolTip(L("SMILES に電荷を書いた成分 ([Na+] など) は、形式電荷から自動で入ります",
                                                                  "for SMILES with a charge (such as [Na+]) this is filled from the formal charge"))
        hdr = self.table.horizontalHeader(); hdr.setSectionResizeMode(COL_REF, QHeaderView.ResizeMode.Stretch); hdr.setMinimumSectionSize(36)
        self.table.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        for c in (COL_KIND, COL_COUNT, COL_CHARGE, COL_LABEL):
            hdr.setSectionResizeMode(c, QHeaderView.ResizeMode.ResizeToContents)
        self.table.verticalHeader().setVisible(False); self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.verticalHeader().setDefaultSectionSize(40)
        self.table.setMinimumHeight(130); self.table.setMaximumHeight(230)
        self.btn_add = QPushButton("成分を追加"); self.btn_del = QPushButton("削除"); self.btn_draw = QPushButton("Draw")
        self.btn_del.setToolTip(tr("選んだ行 (入力欄を触った行) を消します"))
        self.btn_draw.setToolTip(tr("分子を描いて、選んだ行 (無ければ新しい行) に SMILES で入れます (RDKit が要ります)"))
        self.box_mode = QComboBox()
        self.box_mode.addItem(icons.icon("box_density"), "密度から自動 [g/cm³]", "density")
        self.box_mode.addItem(icons.icon("box_edge"), "一辺を指定 [Å]", "edge")
        self.density = narrow(SciDoubleSpinBox(0.01, 30.0, 1.0, 0.05))
        self.edge = narrow(SciDoubleSpinBox(3.0, 500.0, 15.0, 1.0))
        self.min_dist = narrow(SciDoubleSpinBox(0.5, 10.0, 2.0, 0.1))
        self.seed = narrow(QSpinBox()); self.seed.setRange(0, 999999)
        self.info = QLabel(""); self.info.setObjectName("hint"); self.info.setWordWrap(True)
        row1 = QHBoxLayout(); row1.addWidget(self.btn_add); row1.addWidget(self.btn_del); row1.addWidget(self.btn_draw); row1.addStretch()
        self.cell_row = QWidget(); self.cell_row.setObjectName("rowbox")
        row2 = QHBoxLayout(self.cell_row); row2.setContentsMargins(0, 0, 0, 0)
        row2.addWidget(QLabel("セル")); row2.addWidget(self.box_mode); row2.addWidget(self.density); row2.addWidget(self.edge); row2.addStretch()
        self.cell_row.setVisible(cell_row)
        row3 = QHBoxLayout(); row3.addWidget(QLabel("分子間の最短距離 [Å]")); row3.addWidget(self.min_dist); row3.addStretch()
        row4 = QHBoxLayout(); row4.addWidget(QLabel("乱数の種")); row4.addWidget(self.seed); row4.addStretch()
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(6)
        lay.addWidget(self.table); lay.addLayout(row1); lay.addWidget(self.cell_row); lay.addLayout(row3); lay.addLayout(row4); lay.addWidget(self.info)
        self.btn_add.clicked.connect(lambda: self.add_component(Component(ref="H2O", count=1)))
        self.btn_del.clicked.connect(self._delete)
        self.btn_draw.clicked.connect(lambda: self.draw_requested.emit(self.current_row()))
        self.box_mode.currentIndexChanged.connect(self._on_mode)
        for w in (self.density, self.edge, self.min_dist):
            w.valueChanged.connect(self._emit)
        self.seed.valueChanged.connect(self._emit)
        self.table.cellChanged.connect(self._on_cell)
        self._on_mode()
        if example:
            self.add_component(Component(ref="H2O", count=32, label="H2O"))

    def add_component(self, c: Component) -> None:
        r = self.table.rowCount(); self.table.blockSignals(True); self.table.insertRow(r)
        kind = QComboBox()
        for k, v in KINDS.items():
            kind.addItem(icons.source_icon(k), tr(v), k)
        kind.setCurrentIndex(list(KINDS).index(c.kind)); kind.setMaximumWidth(126)
        self.table.setCellWidget(r, COL_KIND, kind)
        ref = RefCell(c.kind, c.ref); ref.btn_draw.setEnabled(self._draw_enabled)
        self.table.setCellWidget(r, COL_REF, ref)
        count = QSpinBox(); count.setRange(0 if self._auto_count else 1, 100000)
        if self._auto_count:
            count.setSpecialValueText(L("自動", "auto")); count.setMaximumWidth(80)
            count.setToolTip(L("0 (自動) にした成分 1 つを、箱全体が指定の密度になる個数だけ入れます", "a component set to 0 (auto) is added until the whole box reaches the given density"))
        count.setValue(c.count); count.setMaximumWidth(64); count.valueChanged.connect(self._emit); self.table.setCellWidget(r, COL_COUNT, count)
        charge = QSpinBox(); charge.setRange(-9, 9); charge.setValue(c.charge); charge.setMaximumWidth(54); charge.valueChanged.connect(self._emit); self.table.setCellWidget(r, COL_CHARGE, charge)
        self.table.setItem(r, COL_LABEL, QTableWidgetItem(c.label))
        kind.currentIndexChanged.connect(lambda *_, w=kind: self._on_kind(self._row_of(w)))
        ref.changed.connect(lambda w=ref: self._on_cell(self._row_of(w), COL_REF))
        ref.draw_requested.connect(lambda w=ref: self.draw_requested.emit(self._row_of(w)))
        for w in (kind, ref, count, charge, *ref.findChildren(QWidget), *count.findChildren(QWidget), *charge.findChildren(QWidget)):
            w.installEventFilter(self)
        self.table.blockSignals(False); self._emit()

    def _row_of(self, w: QWidget) -> int:
        for r in range(self.table.rowCount()):
            for c in (COL_KIND, COL_REF, COL_COUNT, COL_CHARGE):
                cw = self.table.cellWidget(r, c)
                if cw is not None and (cw is w or cw.isAncestorOf(w)):
                    return r
        return -1

    def eventFilter(self, obj, event) -> bool:  # noqa: N802 
        if event.type() in (QEvent.Type.FocusIn, QEvent.Type.MouseButtonPress) and isinstance(obj, QWidget):
            r = self._row_of(obj)
            if r >= 0 and r != self.table.currentRow():
                self.table.selectRow(r)
        return False

    def ref_cell(self, r: int) -> RefCell:
        return self.table.cellWidget(r, COL_REF)

    def set_ref(self, r: int, ref: str) -> None:
        self.ref_cell(r).setText(ref); self._on_cell(r, COL_REF)

    def set_draw_enabled(self, enabled: bool) -> None:
        self._draw_enabled = enabled
        self.btn_draw.setEnabled(enabled)
        for r in range(self.table.rowCount()):
            self.ref_cell(r).btn_draw.setEnabled(enabled)

    def _delete(self) -> None:
        r = self.table.currentRow()
        if r >= 0:
            self.table.removeRow(r); self._emit()

    def current_row(self) -> int:
        return self.table.currentRow()

    def row_smiles(self, r: int) -> str:
        if 0 <= r < self.table.rowCount() and self.ref_cell(r).kind() == "smiles":
            return self.ref_cell(r).text().strip()
        return ""

    def set_row_smiles(self, r: int, smiles: str) -> None:
        if r < 0:
            self.add_component(Component(kind="smiles", ref=smiles, count=1)); r = self.table.rowCount() - 1
            self._on_cell(r, COL_REF); return
        self.table.cellWidget(r, COL_KIND).setCurrentIndex(list(KINDS).index("smiles"))
        self.set_ref(r, smiles)

    def _on_kind(self, r: int) -> None:
        if r < 0:
            return
        self.ref_cell(r).set_kind(self.table.cellWidget(r, COL_KIND).currentData())
        self._on_cell(r, COL_REF)

    def _on_mode(self, *_) -> None:
        dens = self.box_mode.currentData() == "density"
        self.density.setVisible(dens); self.edge.setVisible(not dens); self._emit()

    def spec(self) -> MixtureSpec:
        comps = []
        for r in range(self.table.rowCount()):
            ref = self.ref_cell(r).text().strip()
            label = (self.table.item(r, COL_LABEL).text() if self.table.item(r, COL_LABEL) else "").strip()
            comps.append(Component(kind=self.table.cellWidget(r, COL_KIND).currentData(), ref=ref, count=self.table.cellWidget(r, COL_COUNT).value(),
                                   charge=self.table.cellWidget(r, COL_CHARGE).value(), label=label))
        edge = self.edge.value() if self.box_mode.currentData() == "edge" else 0.0
        return MixtureSpec(components=comps, box_a=edge, density_g_cm3=self.density.value(), min_distance=self.min_dist.value(), seed=self.seed.value())

    def set_spec(self, m: MixtureSpec) -> None:
        self.table.blockSignals(True)
        while self.table.rowCount():
            self.table.removeRow(0)
        self.table.blockSignals(False)
        for c in m.components:
            self.add_component(c)
        self.box_mode.setCurrentIndex(1 if m.box_a > 0 else 0)
        if m.box_a > 0:
            self.edge.setValue(m.box_a)
        self.density.setValue(m.density_g_cm3); self.min_dist.setValue(m.min_distance); self.seed.setValue(m.seed)
        self._emit()

    def show_result(self, box_a: float | None, error: str = "") -> None:
        m = self.spec()
        if error or box_a is None:
            self.info.setText(error); return
        parts = [f"{pretty_formula(c.name())}: {concentration_mol_per_l(c.count, box_a):.2f} mol/L" for c in m.components]
        self.info.setText(L(f"セル {box_a:.2f} Å、全電荷 {m.total_charge():+d}、濃度: ", f"Cell {box_a:.2f} Å, total charge {m.total_charge():+d}, concentrations: ") + ", ".join(parts))

    def _on_cell(self, r: int, c: int) -> None:
        if c == COL_REF and 0 <= r < self.table.rowCount() and self.ref_cell(r) is not None and self.ref_cell(r).kind() == "smiles":
            q = smiles_formal_charge(self.ref_cell(r).text().strip())
            box = self.table.cellWidget(r, COL_CHARGE)
            if q is not None and box is not None and box.value() != q:
                box.blockSignals(True); box.setValue(q); box.blockSignals(False)
        self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()


def smiles_formal_charge(smiles: str) -> int | None:
    if not smiles:
        return None
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        mol = Chem.MolFromSmiles(smiles)
        return None if mol is None else int(Chem.GetFormalCharge(mol))
    except Exception:
        return None
