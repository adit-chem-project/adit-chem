
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QCheckBox, QComboBox, QFormLayout, QGridLayout, QLabel, QPlainTextEdit, QSpinBox, QWidget

from adit.config import Config
from adit.gui.prep_widgets import ElementValues, HubbardTable
from adit.gui.style import ROW_SPACING
from adit.gui.widgets import SciDoubleSpinBox, add_row, label, narrow
from adit.lang import L
from adit.spec import Cp2kMethod
from adit.web.codefields import cp2k_candidates, cp2k_files, parse_sections, sections_text

DISPERSION = {"none": "none", "d3": "D3", "d3bj": "D3(BJ)"}
POISSON = ["", "MT", "WAVELET", "ANALYTIC", "MULTIPOLE"]


def _compact(c: QComboBox) -> None:
    c.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
    c.setMinimumContentsLength(10)


class Cp2kMethodPanel(QWidget):
    changed = Signal()

    def __init__(self, cfg: Config, parent: QWidget | None = None):
        super().__init__(parent)
        self.cfg = cfg
        self.elements: list[str] = []
        self.basis_widgets: dict[str, QComboBox] = {}
        self.potential_widgets: dict[str, QComboBox] = {}

        self.xc = QComboBox(); self.xc.setEditable(True); self.xc.addItems(["", "PBE", "BLYP", "PBE0", "B3LYP", "TPSS"])
        self.xc.lineEdit().setPlaceholderText("例 PBE、BLYP、PBE0 (CP2K の XC_FUNCTIONAL の名前)")
        self.basis_file = QComboBox(); self.basis_file.setEditable(True); _compact(self.basis_file)
        self.potential_file = QComboBox(); self.potential_file.setEditable(True); _compact(self.potential_file)
        self.data_info = QLabel(""); self.data_info.setObjectName("hint"); self.data_info.setWordWrap(True)
        self.kind_box = QWidget(); self.kind_box.setObjectName("rowbox")
        self.kind_grid = QGridLayout(self.kind_box); self.kind_grid.setContentsMargins(0, 0, 0, 0)
        self.kind_grid.setHorizontalSpacing(8); self.kind_grid.setVerticalSpacing(4)
        self.dispersion = QComboBox()
        for k, v in DISPERSION.items():
            self.dispersion.addItem(v, k)
        self.cutoff = narrow(SciDoubleSpinBox(0, 100000, 0.0, 50))
        self.rel_cutoff = narrow(SciDoubleSpinBox(0, 100000, 0.0, 10))
        self.eps_scf = narrow(SciDoubleSpinBox(1e-14, 1.0, 1e-5, 1e-6))
        self.max_scf = narrow(QSpinBox()); self.max_scf.setRange(1, 100000); self.max_scf.setValue(50)
        self.uks = QCheckBox("UKS")
        self.poisson = QComboBox()
        for v in POISSON:
            self.poisson.addItem(v or "(3 方向とも周期なら不要)", v)
        self.box = narrow(SciDoubleSpinBox(0, 10000, 0.0, 1))
        self.mag = ElementValues("MAGNETIZATION", L("空欄 = 指定しない", "empty = not set"))
        self.hubbard = HubbardTable()
        self.plus_u = QComboBox(); self.plus_u.addItems(["MULLIKEN", "LOWDIN", "MULLIKEN_CHARGES"])
        self.sccs = narrow(SciDoubleSpinBox(0, 10000, 0.0, 1))
        self.ot = QCheckBox(L("軌道変換法 (OT) を使う", "Use orbital transformation (OT)")); self.ot_minimizer = QComboBox(); self.ot_minimizer.addItems(["", "SD", "CG", "DIIS", "BROYDEN", "LBFGS"])
        self.ot_preconditioner = QComboBox(); self.ot_preconditioner.addItems(["", "FULL_ALL", "FULL_SINGLE_INVERSE", "FULL_SINGLE", "FULL_KINETIC", "FULL_S_INVERSE", "NONE"])
        self.ot_extra = QPlainTextEdit(); self.ot_extra.setMaximumHeight(60)
        self.surface_dipole = QCheckBox(L("有効にする", "Enable")); self.surf_dip_dir = QComboBox(); self.surf_dip_dir.addItems(["", "X", "Y", "Z"])
        self.extra = QPlainTextEdit(); self.extra.setMaximumHeight(90)
        self.extra.setPlaceholderText("[セクションのパス] の行のあとに、その節の末尾に足す行 (例\n[FORCE_EVAL/DFT/SCF]\nSCF_GUESS ATOMIC)")
        note = QLabel("基底関数と擬ポテンシャルは CP2K の data ディレクトリのファイルから読み、使う項目だけを生成したファイルに写します。"
                      "data ディレクトリの場所は環境設定 (cp2k_data) に書きます")
        note.setObjectName("hint"); note.setWordWrap(True)

        form = QFormLayout(self); form.setVerticalSpacing(ROW_SPACING)
        add_row(form, "汎関数", self.xc)
        add_row(form, "基底関数のファイル", self.basis_file)
        add_row(form, "擬ポテンシャルのファイル", self.potential_file)
        form.addRow(label(""), self.data_info)
        form.addRow(label("元素ごとの基底と擬ポテンシャル"))
        form.addRow(self.kind_box)
        add_row(form, "分散補正", self.dispersion)
        add_row(form, "カットオフ [Ry]", self.cutoff)
        add_row(form, "相対カットオフ [Ry]", self.rel_cutoff)
        add_row(form, "SCF の収束の閾値 (EPS_SCF)", self.eps_scf)
        add_row(form, "SCF の反復の上限 (MAX_SCF)", self.max_scf)
        add_row(form, "スピン分極", self.uks)
        add_row(form, "元素ごとの MAGNETIZATION", self.mag)
        form.addRow(label("DFT+U"))
        form.addRow(self.hubbard)
        add_row(form, "PLUS_U_METHOD", self.plus_u)
        add_row(form, "ポアソン方程式の解き方", self.poisson)
        add_row(form, "分子の箱の一辺 [Å]", self.box)
        add_row(form, "溶媒の比誘電率 (SCCS)", self.sccs)
        add_row(form, L("軌道変換法 (OT)", "Orbital transformation (OT)"), self.ot)
        add_row(form, "OT MINIMIZER", self.ot_minimizer); add_row(form, "OT PRECONDITIONER", self.ot_preconditioner)
        add_row(form, L("OT の追加行", "Additional OT lines"), self.ot_extra)
        add_row(form, "SURFACE_DIPOLE_CORRECTION", self.surface_dipole); add_row(form, "SURF_DIP_DIR", self.surf_dip_dir)
        add_row(form, "追加の行 (節ごと)", self.extra)
        form.addRow(label(""), note)

        self.reload_data()
        self.basis_file.currentTextChanged.connect(self._rebuild_rows)
        self.potential_file.currentTextChanged.connect(self._rebuild_rows)
        for w in (self.xc, self.dispersion, self.poisson, self.plus_u, self.ot_minimizer, self.ot_preconditioner, self.surf_dip_dir):
            w.currentTextChanged.connect(self._emit)
        for w in (self.cutoff, self.rel_cutoff, self.eps_scf, self.max_scf, self.box, self.sccs):
            w.valueChanged.connect(self._emit)
        self.mag.changed.connect(self._emit); self.hubbard.changed.connect(self._emit)
        self.uks.toggled.connect(self._emit); self.ot.toggled.connect(self._emit); self.surface_dipole.toggled.connect(self._emit)
        self.ot_extra.textChanged.connect(self._emit)
        self.extra.textChanged.connect(self._emit)

    def reload_data(self, cfg: Config | None = None) -> None:
        if cfg is not None:
            self.cfg = cfg
        basis, pots = cp2k_files(self.cfg.cp2k_data)
        for combo, names, default in ((self.basis_file, basis, "BASIS_MOLOPT"), (self.potential_file, pots, "GTH_POTENTIALS")):
            cur = combo.currentText() or default
            combo.blockSignals(True); combo.clear(); combo.addItems(names or [default])
            if cur not in names and names:
                combo.addItem(cur)
            combo.setCurrentText(cur); combo.blockSignals(False)
        self._rebuild_rows()

    def set_context(self, elements: list[str], profile_name: str, cfg: Config | None = None) -> None:
        if cfg is not None:
            self.cfg = cfg
        self.mag.set_elements(elements); self.hubbard.set_elements(elements)
        if elements != self.elements:
            self.elements = list(elements)
            self._rebuild_rows()

    def _rebuild_rows(self, *_) -> None:
        old_b = {e: w.currentText() for e, w in self.basis_widgets.items()}
        old_p = {e: w.currentText() for e, w in self.potential_widgets.items()}
        while self.kind_grid.count():
            item = self.kind_grid.takeAt(0)
            if item.widget():
                w = item.widget(); w.hide(); w.setParent(None); w.deleteLater()
        self.basis_widgets, self.potential_widgets = {}, {}
        root, cands = cp2k_candidates(self.cfg.cp2k_data, self.basis_file.currentText().strip(), self.potential_file.currentText().strip(), self.elements)
        if root:
            shown = root.replace("/", "/​")
            self.data_info.setText(L(f"data ディレクトリ: {shown}", f"data directory: {shown}"))
        else:
            self.data_info.setText(L("この PC に CP2K の data ディレクトリが見つかりません。元素ごとの名前を入れてください (名前は確かめずに書きます)。"
                                     "場所は環境設定の cp2k_data に書けます",
                                     "No CP2K data directory on this PC. Enter the names per element (they are written unchecked). "
                                     "The location can be set as cp2k_data in the settings"))
        by_el = {e: (b, p) for e, b, p in cands}
        if self.elements:
            for c, t in ((1, L("基底関数", "Basis set")), (2, L("擬ポテンシャル", "Pseudopotential"))):
                h = QLabel(t); h.setObjectName("hint"); self.kind_grid.addWidget(h, 0, c)
        for i, e in enumerate(self.elements, start=1):
            b_names, p_names = by_el.get(e, ([], []))
            combos = []
            for names, old, store in ((b_names, old_b, self.basis_widgets), (p_names, old_p, self.potential_widgets)):
                c = QComboBox(); c.setEditable(True); _compact(c)
                c.addItems([""] + names)
                ph = (L(f"自動: {names[0]}", f"auto: {names[0]}") if len(names) == 1 else
                      L("選んでください", "choose one") if names else
                      L("ファイルにこの元素がありません", "no entry for this element") if root else L("名前を入力", "enter a name"))
                c.lineEdit().setPlaceholderText(ph)
                c.setCurrentText(old.get(e, ""))
                c.currentTextChanged.connect(self._emit)
                store[e] = c; combos.append(c)
            self.kind_grid.addWidget(QLabel(e), i, 0); self.kind_grid.addWidget(combos[0], i, 1); self.kind_grid.addWidget(combos[1], i, 2)
        self.kind_grid.setColumnStretch(1, 1); self.kind_grid.setColumnStretch(2, 1)
        self._emit()

    def method(self) -> Cp2kMethod:
        basis = {e: w.currentText().strip() for e, w in self.basis_widgets.items() if w.currentText().strip()}
        pots = {e: w.currentText().strip() for e, w in self.potential_widgets.items() if w.currentText().strip()}
        return Cp2kMethod(basis_file=self.basis_file.currentText().strip() or "BASIS_MOLOPT",
                          potential_file=self.potential_file.currentText().strip() or "GTH_POTENTIALS",
                          basis=basis, potential=pots, xc=self.xc.currentText().strip(), dispersion=self.dispersion.currentData(),
                          cutoff_ry=self.cutoff.value(), rel_cutoff_ry=self.rel_cutoff.value(), eps_scf=self.eps_scf.value(),
                          max_scf=self.max_scf.value(), uks=self.uks.isChecked(), poisson_solver=self.poisson.currentData(),
                          isolated_box_ang=self.box.value(), extra_sections=parse_sections(self.extra.toPlainText()),
                          magnetization_by_element=self.mag.values(), hubbard=self.hubbard.hubbard(), plus_u_method=self.plus_u.currentText(),
                          sccs_relative_permittivity=self.sccs.value(), ot=self.ot.isChecked(), ot_minimizer=self.ot_minimizer.currentText(),
                          ot_preconditioner=self.ot_preconditioner.currentText(), ot_extra=self.ot_extra.toPlainText().strip(),
                          surface_dipole_correction=self.surface_dipole.isChecked(), surf_dip_dir=self.surf_dip_dir.currentText())

    def set_method(self, m: Cp2kMethod) -> None:
        self.basis_file.setCurrentText(m.basis_file); self.potential_file.setCurrentText(m.potential_file)
        self.xc.setCurrentText(m.xc); self.dispersion.setCurrentIndex(max(0, self.dispersion.findData(m.dispersion)))
        self.cutoff.setValue(m.cutoff_ry); self.rel_cutoff.setValue(m.rel_cutoff_ry); self.eps_scf.setValue(m.eps_scf)
        self.max_scf.setValue(m.max_scf); self.uks.setChecked(m.uks)
        self.mag.set_values(m.magnetization_by_element); self.hubbard.set_hubbard(m.hubbard)
        self.plus_u.setCurrentText(m.plus_u_method); self.sccs.setValue(m.sccs_relative_permittivity)
        self.ot.setChecked(m.ot); self.ot_minimizer.setCurrentText(m.ot_minimizer); self.ot_preconditioner.setCurrentText(m.ot_preconditioner)
        self.ot_extra.setPlainText(m.ot_extra); self.surface_dipole.setChecked(m.surface_dipole_correction); self.surf_dip_dir.setCurrentText(m.surf_dip_dir)
        self.poisson.setCurrentIndex(max(0, self.poisson.findData(m.poisson_solver))); self.box.setValue(m.isolated_box_ang)
        self.extra.setPlainText(sections_text(m.extra_sections))
        self._rebuild_rows()
        for table, widgets in ((m.basis, self.basis_widgets), (m.potential, self.potential_widgets)):
            for e, w in widgets.items():
                w.setCurrentText(table.get(e, ""))
        self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()
