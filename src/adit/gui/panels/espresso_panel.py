
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QComboBox, QFormLayout, QGridLayout, QLabel, QLineEdit, QPlainTextEdit, QPushButton,
                               QSpinBox, QVBoxLayout, QWidget)

from adit.textparse import parse_extra_namelist as parse_extra
from adit.lang import L
from adit.codes.upf import UpfError, UpfLibrary
from adit.config import Config
from adit.gui.prep_widgets import ElementValues, HubbardTable
from adit.gui.style import ROW_SPACING
from adit.gui.widgets import SciDoubleSpinBox, add_row, label, narrow
from adit.spec import EspressoMethod

SSSP_URL = "https://www.materialscloud.org/discover/sssp"


class EspressoMethodPanel(QWidget):
    changed = Signal()
    root_chosen = Signal(str)

    def __init__(self, cfg: Config, parent: QWidget | None = None):
        super().__init__(parent)
        self.cfg = cfg
        self.elements: list[str] = []
        self.lib: UpfLibrary | None = None
        self.pseudo_widgets: dict[str, QComboBox] = {}
        self.pseudo_set = QComboBox(); self.pseudo_set.setEditable(True)
        self.lib_info = QLabel(""); self.lib_info.setObjectName("hint"); self.lib_info.setWordWrap(True)
        self.pseudo_grid = QGridLayout(); self.pseudo_grid.setHorizontalSpacing(8); self.pseudo_grid.setVerticalSpacing(4)
        self.ecutwfc = narrow(SciDoubleSpinBox(0, 10000, 0.0, 5))
        self.ecutrho = narrow(SciDoubleSpinBox(0, 100000, 0.0, 20))
        self.conv_thr = narrow(SciDoubleSpinBox(1e-14, 1, 1e-6, 1e-7))
        self.maxstep = narrow(QSpinBox()); self.maxstep.setRange(1, 10000); self.maxstep.setValue(100)
        self.mixing = narrow(SciDoubleSpinBox(0.01, 1.0, 0.7, 0.1))
        self.occupations = QComboBox(); self.occupations.addItems(["fixed", "smearing", "tetrahedra"])
        self.smearing = QComboBox(); self.smearing.addItems(["gaussian", "methfessel-paxton", "marzari-vanderbilt", "fermi-dirac"])
        self.degauss = narrow(SciDoubleSpinBox(0, 10, 0.0, 0.005))
        self.nspin = QComboBox(); self.nspin.addItems(["1", "2"])
        self.input_dft = QLineEdit(); self.input_dft.setPlaceholderText("空欄なら UPF の汎関数を使用 (例 PBE, PBEsol)")
        self.assume_isolated = QComboBox()
        for value, text in (("", L("指定しない", "Not set")), ("makov-payne", "makov-payne"), ("martyna-tuckerman", "martyna-tuckerman"), ("esm", "esm"), ("2D", "2D")):
            self.assume_isolated.addItem(text, value)
        self.mag = ElementValues("starting_magnetization", L("空欄 = 指定しない", "empty = not set"))
        self.hubbard = HubbardTable(with_j=False)
        self._dipole: dict = {}                 # dipole_* fields have no widgets; kept for the round trip
        self.hub_proj = QComboBox(); self.hub_proj.addItems(["atomic", "ortho-atomic", "norm-atomic", "wf", "pseudo"])
        self.extra = QPlainTextEdit(); self.extra.setMaximumHeight(80)
        self.extra.setPlaceholderText("画面にない変数を 1 行に 1 つ、名前空間.変数 = 値 の形で (例 system.nbnd = 20)")

        self.lib_info.setOpenExternalLinks(True)
        self.root_browse = QPushButton("フォルダを選ぶ…"); self.root_browse.setObjectName("link")
        self.root_browse.setToolTip("UPF のセット (SSSP など) のフォルダを並べた場所を選びます。選んだ場所は環境設定 (pseudo_root) に保存されます")
        self.root_browse.clicked.connect(self._browse_root)

        form = QFormLayout(self); form.setVerticalSpacing(ROW_SPACING)
        add_row(form, "擬ポテンシャルのセット", self.pseudo_set)
        info_row = QVBoxLayout(); info_row.setSpacing(2); info_row.addWidget(self.lib_info); info_row.addWidget(self.root_browse, 0, Qt.AlignmentFlag.AlignLeft)
        form.addRow(label(""), info_row)
        add_row(form, "元素ごとの UPF ファイル", self.pseudo_grid)
        add_row(form, "ecutwfc [Ry]", self.ecutwfc)
        add_row(form, "ecutrho [Ry] (0 = 指定しない)", self.ecutrho)
        add_row(form, "conv_thr [Ry]", self.conv_thr)
        add_row(form, "electron_maxstep", self.maxstep)
        add_row(form, "mixing_beta", self.mixing)
        add_row(form, "occupations", self.occupations)
        add_row(form, "smearing", self.smearing)
        add_row(form, "degauss [Ry]", self.degauss)
        add_row(form, "nspin", self.nspin)
        add_row(form, "starting_magnetization (元素ごと)", self.mag)
        form.addRow(label("DFT+U"))
        form.addRow(self.hubbard)
        add_row(form, "HUBBARD の射影", self.hub_proj)
        add_row(form, "input_dft", self.input_dft)
        add_row(form, "assume_isolated", self.assume_isolated)
        add_row(form, "追加の変数 (名前空間.変数 = 値)", self.extra)

        self.reload_sets()
        self.pseudo_set.currentTextChanged.connect(self._reload_library)
        for w in (self.occupations, self.smearing, self.nspin, self.assume_isolated):
            w.currentTextChanged.connect(self._emit)
        for w in (self.ecutwfc, self.ecutrho, self.conv_thr, self.maxstep, self.mixing, self.degauss):
            w.valueChanged.connect(self._emit)
        self.input_dft.textChanged.connect(self._emit)
        self.extra.textChanged.connect(self._emit)
        self.mag.changed.connect(self._emit); self.hubbard.changed.connect(self._emit)
        self.hub_proj.currentTextChanged.connect(self._emit)
        self.occupations.currentTextChanged.connect(self._on_occ)
        self._on_occ()

    def reload_sets(self, cfg: Config | None = None) -> None:
        if cfg is not None:
            self.cfg = cfg
        root = Path(self.cfg.pseudo_root).expanduser() if self.cfg.pseudo_root else None
        sets = sorted(d.name for d in root.iterdir() if d.is_dir() and any(p.suffix.lower() == ".upf" for p in d.iterdir())) if root and root.is_dir() else []
        cur = self.pseudo_set.currentText()
        self.pseudo_set.clear(); self.pseudo_set.addItems(sets)
        if cur in sets:
            self.pseudo_set.setCurrentText(cur)
        self._reload_library()

    def set_context(self, elements: list[str], profile_name: str, cfg: Config | None = None) -> None:
        if cfg is not None:
            self.cfg = cfg
        self.mag.set_elements(elements); self.hubbard.set_elements(elements)
        if elements != self.elements:
            self.elements = list(elements)
            self._rebuild_rows()

    def _reload_library(self, *_) -> None:
        self.lib = UpfLibrary.open_if_present(self.cfg.pseudo_root or None, self.pseudo_set.currentText().strip() or "-")
        if self.lib is None:
            self.lib_info.setText(L(f"擬ポテンシャル (UPF) が見つかりません。<a href=\"{SSSP_URL}\">SSSP</a> などから入手し、"
                                    "「置き場所/セット名/*.UPF」の形に置いて、「フォルダを選ぶ…」で置き場所を指定してください"
                                    + (f" (いまの置き場所: {self.cfg.pseudo_root})" if self.cfg.pseudo_root else ""),
                                    f"No pseudopotentials (UPF) found. Download e.g. <a href=\"{SSSP_URL}\">SSSP</a>, arrange as "
                                    "\"folder/set name/*.UPF\", and pick the folder with \"Choose folder…\""
                                    + (f" (current: {self.cfg.pseudo_root})" if self.cfg.pseudo_root else "")))
        else:
            docs = ", ".join(self.lib.doc_files()) or L("なし", "none")
            self.lib_info.setText(L(f"ライブラリ: {self.lib.set_dir}   付属の文書: {docs}", f"Library: {self.lib.set_dir}   documents: {docs}"))
        self._rebuild_rows()

    def _browse_root(self) -> None:
        from adit.gui.panels.method_panel import choose_root

        root = choose_root(self, L("擬ポテンシャル (UPF) の置き場所", "Folder of pseudopotential (UPF) sets"), ".upf")
        if root:
            self.root_chosen.emit(root)

    def _rebuild_rows(self) -> None:
        old = {e: w.currentText() for e, w in self.pseudo_widgets.items()}
        while self.pseudo_grid.count():
            item = self.pseudo_grid.takeAt(0)
            if item.widget():
                w = item.widget(); w.hide(); w.setParent(None); w.deleteLater()
        self.pseudo_widgets = {}
        for i, e in enumerate(self.elements):
            combo = QComboBox(); combo.setEditable(True)
            names = self.lib.files_for(e) if self.lib else []
            combo.addItems(names)
            if old.get(e):
                combo.setCurrentText(old[e])
            info = "ライブラリに無い" if self.lib and not names else ""
            if self.lib and combo.currentText():
                try:
                    h = self.lib.header(combo.currentText())
                    info = f"{h.functional} {h.pseudo_type}  Z_val {h.z_valence:g}" + (f"  wfc_cutoff {h.wfc_cutoff:g} Ry" if h.wfc_cutoff else "")
                except UpfError as ex:
                    info = str(ex)
            lab = QLabel(info); lab.setObjectName("hint")
            self.pseudo_grid.addWidget(QLabel(e), i, 0); self.pseudo_grid.addWidget(combo, i, 1); self.pseudo_grid.addWidget(lab, i, 2)
            combo.currentTextChanged.connect(self._emit)
            self.pseudo_widgets[e] = combo
        self.pseudo_grid.setColumnStretch(2, 1)
        self._emit()

    def _on_occ(self, *_) -> None:
        on = self.occupations.currentText() == "smearing"
        self.smearing.setEnabled(on); self.degauss.setEnabled(on)
        self._emit()

    def method(self) -> EspressoMethod:
        pseudo = {e: w.currentText().strip() for e, w in self.pseudo_widgets.items() if w.currentText().strip()}
        return EspressoMethod(pseudo_set=self.pseudo_set.currentText().strip(), pseudo=pseudo, ecutwfc=self.ecutwfc.value(),
                              ecutrho=self.ecutrho.value(), conv_thr=self.conv_thr.value(), electron_maxstep=self.maxstep.value(),
                              mixing_beta=self.mixing.value(), occupations=self.occupations.currentText(), smearing=self.smearing.currentText(),
                              degauss=self.degauss.value(), nspin=int(self.nspin.currentText()), input_dft=self.input_dft.text().strip(),
                              extra=parse_extra(self.extra.toPlainText()), starting_magnetization=self.mag.values(),
                              hubbard=self.hubbard.hubbard(), hubbard_projector=self.hub_proj.currentText(),
                              assume_isolated=self.assume_isolated.currentData(), **self._dipole)

    def set_method(self, m: EspressoMethod) -> None:
        self._dipole = {k: getattr(m, k) for k in ("dipole_correction", "dipole_direction", "dipole_maxpos",
                                                    "dipole_decrease", "dipole_amplitude")}
        self.pseudo_set.setCurrentText(m.pseudo_set); self.ecutwfc.setValue(m.ecutwfc); self.ecutrho.setValue(m.ecutrho)
        self.conv_thr.setValue(m.conv_thr); self.maxstep.setValue(m.electron_maxstep); self.mixing.setValue(m.mixing_beta)
        self.occupations.setCurrentText(m.occupations); self.smearing.setCurrentText(m.smearing); self.degauss.setValue(m.degauss)
        self.nspin.setCurrentText(str(m.nspin)); self.input_dft.setText(m.input_dft)
        self.assume_isolated.setCurrentIndex(max(0, self.assume_isolated.findData(m.assume_isolated)))
        self.extra.setPlainText("\n".join(f"{ns}.{k} = {v}" for ns, d in m.extra.items() for k, v in d.items()))
        self.mag.set_values(m.starting_magnetization); self.hubbard.set_hubbard(m.hubbard); self.hub_proj.setCurrentText(m.hubbard_projector)
        self._reload_library()
        for e, name in m.pseudo.items():
            if e in self.pseudo_widgets:
                self.pseudo_widgets[e].setCurrentText(name)
        self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()
