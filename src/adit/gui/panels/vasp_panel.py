
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import (QCheckBox, QComboBox, QFormLayout, QGridLayout, QLabel, QLineEdit, QPlainTextEdit, QSpinBox, QWidget)

from adit.textparse import parse_extra_incar, short_number
from adit.lang import L
from adit.codes.potcar import PotcarError, PotcarLibrary
from adit.codes.potcar_names import MP_POTCAR_NAMES
from adit.codes.vasp import PP_ENV
from adit.config import Config
from adit.gui.prep_widgets import ElementValues, HubbardTable
from adit.gui.style import ROW_SPACING
from adit.gui.widgets import SciDoubleSpinBox, add_row, label, narrow
from adit.spec import VaspMethod


def _pick_text(combo: QComboBox, text: str) -> None:
    # A value that is not in the list is added rather than silently replaced by the current one.
    if combo.findText(text) < 0:
        combo.addItem(text)
    combo.setCurrentText(text)


class VaspMethodPanel(QWidget):
    changed = Signal()

    def __init__(self, cfg: Config, parent: QWidget | None = None):
        super().__init__(parent)
        self.cfg = cfg
        self.profile_name = cfg.default_profile
        self.elements: list[str] = []
        self.lib: PotcarLibrary | None = None
        self.potcar_widgets: dict[str, QComboBox] = {}

        self.potcar_set = QLineEdit("potpaw_PBE")
        self.lib_info = QLabel(""); self.lib_info.setObjectName("hint"); self.lib_info.setWordWrap(True)
        self.potcar_grid = QGridLayout(); self.potcar_grid.setHorizontalSpacing(8); self.potcar_grid.setVerticalSpacing(4)
        self.binary = QComboBox(); self.binary.addItems(["std", "gam", "ncl"])
        self.ibrion = QComboBox()
        for v, t in ((2, L("2: 共役勾配", "2: conjugate gradient")), (1, L("1: RMM-DIIS (準ニュートン)", "1: RMM-DIIS (quasi-Newton)")), (3, L("3: 減衰 MD", "3: damped MD"))):
            self.ibrion.addItem(t, v)
        self.encut = narrow(SciDoubleSpinBox(0, 5000, 0.0, 10))
        self.prec = QComboBox(); self.prec.addItems(["Normal", "Accurate", "Low", "Medium", "High", "Single"])
        self.algo = QComboBox(); self.algo.addItems(["Normal", "Fast", "VeryFast", "All", "Damped", "Conjugate"])
        self.ediff = narrow(SciDoubleSpinBox(1e-12, 1.0, 1e-4, 1e-5))
        self.nelm = narrow(QSpinBox()); self.nelm.setRange(1, 10000); self.nelm.setValue(60)
        self.nelmin = narrow(QSpinBox()); self.nelmin.setRange(0, 10000)
        self.ismear = narrow(QSpinBox()); self.ismear.setRange(-5, 5)
        self.sigma = narrow(SciDoubleSpinBox(0, 10, 0.1, 0.01))
        self.ispin = QComboBox(); self.ispin.addItems(["1", "2"])
        self.magmom = QLineEdit(); self.magmom.setPlaceholderText("空欄なら指定しない。原子ごとに空白区切り (例 1 1 -1)")
        self.mag = ElementValues("MAGMOM", L("空欄 = 0", "empty = 0"))
        self.hubbard = HubbardTable()
        self.ldau_type = QComboBox()
        for v, t in ((2, L("2: Dudarev (U − J だけが効く)", "2: Dudarev (only U − J matters)")), (1, L("1: Liechtenstein", "1: Liechtenstein")),
                     (4, L("4: 1 と同じで LSDA の交換分裂なし", "4: as 1, without LSDA exchange splitting"))):
            self.ldau_type.addItem(t, v)
        self.ivdw = QComboBox(); self.ivdw.addItems(["none", "10", "11", "12", "20", "4"])
        self.lreal = QComboBox(); self.lreal.addItems(["Auto", ".FALSE.", "On"])
        self.lasph = QCheckBox(L("有効にする", "Enable"))
        self.lmaxmix = narrow(QSpinBox()); self.lmaxmix.setRange(0, 100)
        self.nbands = narrow(QSpinBox()); self.nbands.setRange(0, 10_000_000)
        self.isym = narrow(QLineEdit()); self.isym.setPlaceholderText(L("空欄なら VASP の既定", "Empty = VASP default"))
        self.idipol = QComboBox()
        for value, text in ((0, L("指定しない", "Not set")), (1, "1: a"), (2, "2: b"), (3, "3: c"), (4, "4: molecule")):
            self.idipol.addItem(text, value)
        self.ldipol = QCheckBox(L("ポテンシャルと力も補正", "Correct potential and forces"))
        self.dipol = QLineEdit(); self.dipol.setPlaceholderText(L("格子座標の 3 成分 (例 0.5 0.5 0.5)", "Three direct coordinates (e.g. 0.5 0.5 0.5)"))
        self.kpoints_centering = QComboBox()
        self.kpoints_centering.addItem(L("Monkhorst–Pack", "Monkhorst-Pack"), "monkhorst-pack")
        self.kpoints_centering.addItem(L("Γ 中心", "Gamma-centered"), "gamma")
        self.extra = QPlainTextEdit(); self.extra.setPlaceholderText("画面にない INCAR のキーを 1 行に 1 つ、KEY = value の形で (例 NCORE = 4)")
        self.extra.setMaximumHeight(90)

        form = QFormLayout(self)
        form.setVerticalSpacing(ROW_SPACING)
        add_row(form, "POTCAR のセット", self.potcar_set)
        form.addRow(label(""), self.lib_info)
        add_row(form, "元素ごとの POTCAR 名", self.potcar_grid)
        add_row(form, "実行ファイルの種類", self.binary)
        add_row(form, "IBRION (構造最適化の方法)", self.ibrion)
        add_row(form, "ENCUT [eV] (0 = 指定しない)", self.encut)
        add_row(form, "PREC", self.prec)
        add_row(form, "ALGO", self.algo)
        add_row(form, "EDIFF [eV]", self.ediff)
        add_row(form, "NELM", self.nelm)
        add_row(form, L("NELMIN (0 = 指定しない)", "NELMIN (0 = omit)"), self.nelmin)
        add_row(form, "ISMEAR", self.ismear)
        add_row(form, "SIGMA [eV]", self.sigma)
        add_row(form, "ISPIN", self.ispin)
        add_row(form, "MAGMOM", self.magmom)
        add_row(form, "元素ごとの初期磁気モーメント [μB]", self.mag)
        form.addRow(label("DFT+U"))
        form.addRow(self.hubbard)
        add_row(form, "LDAUTYPE", self.ldau_type)
        add_row(form, "IVDW", self.ivdw)
        add_row(form, "LREAL", self.lreal)
        add_row(form, "LASPH", self.lasph)
        add_row(form, L("LMAXMIX (0 = 指定しない)", "LMAXMIX (0 = omit)"), self.lmaxmix)
        add_row(form, L("NBANDS (0 = 指定しない)", "NBANDS (0 = omit)"), self.nbands)
        add_row(form, L("ISYM (空欄 = 指定しない)", "ISYM (empty = omit)"), self.isym)
        add_row(form, "IDIPOL", self.idipol); add_row(form, "LDIPOL", self.ldipol); add_row(form, "DIPOL", self.dipol)
        add_row(form, L("k 点メッシュの中心", "k-point mesh centering"), self.kpoints_centering)
        add_row(form, "追加の INCAR 設定", self.extra)

        self.potcar_set.editingFinished.connect(self._reload_library)
        for w in (self.binary, self.prec, self.algo, self.ispin, self.ivdw, self.lreal, self.ibrion, self.kpoints_centering):
            w.currentTextChanged.connect(self._emit)
        for w in (self.encut, self.ediff, self.nelm, self.nelmin, self.ismear, self.sigma, self.lmaxmix, self.nbands):
            w.valueChanged.connect(self._emit)
        self.magmom.textChanged.connect(self._emit); self.isym.textChanged.connect(self._emit); self.dipol.textChanged.connect(self._emit)
        self.lasph.toggled.connect(self._emit); self.ldipol.toggled.connect(self._emit); self.idipol.currentIndexChanged.connect(self._emit)
        self.extra.textChanged.connect(self._emit)
        self.mag.changed.connect(self._emit); self.hubbard.changed.connect(self._emit)
        self.ldau_type.currentIndexChanged.connect(self._emit)
        self._reload_library()

    def set_context(self, elements: list[str], profile_name: str, cfg: Config | None = None) -> None:
        if cfg is not None:
            self.cfg = cfg
        changed = elements != self.elements or profile_name != self.profile_name
        self.elements, self.profile_name = list(elements), profile_name
        self.mag.set_elements(self.elements); self.hubbard.set_elements(self.elements)
        if changed:
            self._reload_library()

    def _reload_library(self) -> None:
        prof = self.cfg.profiles.get(self.profile_name)
        root = prof.env.get(PP_ENV) if prof else None
        self.lib = PotcarLibrary.open_if_present(root, self.potcar_set.text().strip() or "potpaw_PBE")
        if self.lib is None:
            self.lib_info.setText(L(f"この PC には POTCAR のライブラリが見つかりません (プロファイル {self.profile_name} の env.{PP_ENV} を確認してください)。"
                                    "POTCAR 名は未確認のまま生成し、実行時に make_potcar.sh が確認します。",
                                    f"POTCAR library not on this machine (env.{PP_ENV} of profile {self.profile_name}). "
                                    "Names are written unverified; make_potcar.sh checks them at run time"))
        else:
            self.lib_info.setText(L(f"POTCAR のライブラリ: {self.lib.set_dir} (ヘッダを確認済みです)", f"POTCAR library: {self.lib.set_dir} (headers checked)"))
        self._rebuild_potcar_rows()

    def _rebuild_potcar_rows(self) -> None:
        old = {e: w.currentText() for e, w in self.potcar_widgets.items()}
        while self.potcar_grid.count():
            item = self.potcar_grid.takeAt(0)
            if item.widget():
                w = item.widget(); w.hide(); w.setParent(None); w.deleteLater()
        self.potcar_widgets = {}
        for i, e in enumerate(self.elements):
            combo = QComboBox(); combo.setEditable(True)
            names = self.lib.names_for(e) if self.lib else []
            default = old.get(e) or MP_POTCAR_NAMES.get(e, e)
            combo.addItems(names or [default])
            if default not in names and names:
                combo.addItem(default)
            combo.setCurrentText(default)
            info = ""
            if self.lib and self.lib.has(combo.currentText()):
                try:
                    h = self.lib.header(combo.currentText())
                    info = f"{h.titel}  ZVAL {h.zval:g}  ENMAX {h.enmax:g}"
                except PotcarError as ex:
                    info = str(ex)
            lab = QLabel(info); lab.setObjectName("hint")
            self.potcar_grid.addWidget(QLabel(e), i, 0); self.potcar_grid.addWidget(combo, i, 1); self.potcar_grid.addWidget(lab, i, 2)
            combo.currentTextChanged.connect(self._emit)
            self.potcar_widgets[e] = combo
        self.potcar_grid.setColumnStretch(2, 1)
        self._emit()

    def method(self) -> VaspMethod:
        magmom = [float(x) for x in self.magmom.text().split()] if self.magmom.text().strip() else None
        ivdw = None if self.ivdw.currentText() == "none" else int(self.ivdw.currentText())
        extra = parse_extra_incar(self.extra.toPlainText())
        potcar = {e: w.currentText().strip() for e, w in self.potcar_widgets.items() if w.currentText().strip()}
        return VaspMethod(potcar_set=self.potcar_set.text().strip() or "potpaw_PBE", potcar=potcar, binary=self.binary.currentText(),
                          encut=self.encut.value(), ediff=self.ediff.value(), nelm=self.nelm.value(), ismear=self.ismear.value(),
                          sigma=self.sigma.value(), ispin=int(self.ispin.currentText()), magmom=magmom, ivdw=ivdw, ibrion=self.ibrion.currentData(),
                          algo=self.algo.currentText(), prec=self.prec.currentText(), lreal=self.lreal.currentText(), extra_incar=extra,
                          kpoints_centering=self.kpoints_centering.currentData(),
                          nelmin=self.nelmin.value(), lasph=self.lasph.isChecked(), lmaxmix=self.lmaxmix.value(), nbands=self.nbands.value(),
                          isym=int(self.isym.text()) if self.isym.text().strip() else None, idipol=self.idipol.currentData(),
                          ldipol=self.ldipol.isChecked(), dipol=self.dipol.text().strip(),
                          magmom_by_element=self.mag.values(), hubbard=self.hubbard.hubbard(), ldau_type=self.ldau_type.currentData())

    def set_method(self, m: VaspMethod) -> None:
        self.potcar_set.setText(m.potcar_set); self.binary.setCurrentText(m.binary); self.encut.setValue(m.encut)
        if self.ibrion.findData(m.ibrion) < 0:
            self.ibrion.addItem(str(m.ibrion), m.ibrion)
        self.ibrion.setCurrentIndex(self.ibrion.findData(m.ibrion))
        _pick_text(self.prec, m.prec); _pick_text(self.algo, m.algo); self.ediff.setValue(m.ediff); self.nelm.setValue(m.nelm)
        self.nelmin.setValue(m.nelmin); self.lasph.setChecked(m.lasph); self.lmaxmix.setValue(m.lmaxmix); self.nbands.setValue(m.nbands)
        self.isym.setText("" if m.isym is None else str(m.isym)); self.idipol.setCurrentIndex(self.idipol.findData(m.idipol))
        self.ldipol.setChecked(m.ldipol); self.dipol.setText(m.dipol)
        self.ismear.setValue(m.ismear); self.sigma.setValue(m.sigma); self.ispin.setCurrentText(str(m.ispin))
        self.magmom.setText(" ".join(short_number(x) for x in m.magmom) if m.magmom else "")
        _pick_text(self.ivdw, "none" if m.ivdw is None else str(m.ivdw)); _pick_text(self.lreal, str(m.lreal))
        self.kpoints_centering.setCurrentIndex(max(0, self.kpoints_centering.findData(m.kpoints_centering)))
        self.extra.setPlainText("\n".join(f"{k} = {'.TRUE.' if v is True else '.FALSE.' if v is False else v}" for k, v in m.extra_incar.items()))
        self.mag.set_values(m.magmom_by_element); self.hubbard.set_hubbard(m.hubbard)
        self.ldau_type.setCurrentIndex(max(0, self.ldau_type.findData(m.ldau_type)))
        self._reload_library()
        for e, name in m.potcar.items():
            if e in self.potcar_widgets:
                self.potcar_widgets[e].setCurrentText(name)
        self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()
