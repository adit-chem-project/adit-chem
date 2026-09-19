
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QPushButton, QVBoxLayout
from PySide6.QtWidgets import (QCheckBox, QComboBox, QFormLayout, QGridLayout, QGroupBox, QLabel, QLineEdit, QSpinBox,
                               QVBoxLayout, QWidget)

from adit.lang import L
from adit.codes.sk_sets import SKSet, discover_sets
from adit.gui import icons
from adit.gui.style import ROW_SPACING
from adit.gui.widgets import SciDoubleSpinBox, add_row, file_row, label, narrow
from adit.config import Config
from adit.gui.panels.vasp_panel import VaspMethodPanel
from adit.gui.panels.espresso_panel import EspressoMethodPanel
from adit.gui.panels.orca_panel import OrcaMethodPanel
from adit.gui.panels.xtb_panel import XtbMethodPanel
from adit.gui.panels.cp2k_panel import Cp2kMethodPanel
from adit.gui.panels.lammps_panel import LammpsMethodPanel
from adit.gui.panels.gromacs_panel import GromacsMethodPanel
from adit.gui.panels.mlip_panel import MlipMethodPanel
from adit.spec import DftbMethod

DISPERSION = ["none", "dftd3", "lennard-jones"]


SK_DOWNLOAD_URL = "https://dftb.org/parameters/download.html"


def choose_root(parent: QWidget, title: str, marker_suffix: str) -> str:
    from PySide6.QtWidgets import QFileDialog

    d = QFileDialog.getExistingDirectory(parent, title)
    if not d:
        return ""
    p = Path(d)
    if any(f.suffix.lower() == marker_suffix for f in p.iterdir() if f.is_file()):
        p = p.parent
    return str(p)


class DftbMethodPanel(QWidget):
    changed = Signal()
    root_chosen = Signal(str)

    def __init__(self, sk_root: str, parent: QWidget | None = None):
        super().__init__(parent)
        self.sets: dict[str, SKSet] = {}
        self.unshown: list[str] = []
        self._seed = DftbMethod.model_fields["seed"].default   # no widget; kept for the round trip
        self.sk_set = QComboBox()
        self.sk_info = QLabel("")
        self.sk_info.setWordWrap(True)
        self.scc = QCheckBox("SCC (自己無撞着電荷)"); self.scc.setChecked(True)
        self.scc_tol = narrow(SciDoubleSpinBox(1e-12, 1.0, 1e-5, 1e-6))
        self.max_scc = narrow(QSpinBox()); self.max_scc.setRange(1, 100000); self.max_scc.setValue(100)
        self.third = QCheckBox("DFTB3 (ThirdOrderFull)")
        self.dispersion = QComboBox(); self.dispersion.addItems(DISPERSION)
        self.d3 = {k: narrow(SciDoubleSpinBox(-100, 100, 0.0, 0.01)) for k in ("s6", "s8", "a1", "a2")}
        self.temperature = narrow(SciDoubleSpinBox(0, 1e6, 0.0, 10))
        self.solv_file = QLineEdit(); self.solv_file.setPlaceholderText("空欄なら溶媒なし (param_gbsa_<溶媒>.txt の形のファイル)")
        self._solv_ph = self.solv_file.placeholderText()
        self.solv_row = file_row(self.solv_file, L("溶媒のパラメータファイル (GBSA)", "Solvation parameter file (GBSA)"))
        for w in (self.scc, self.third):
            w.setProperty("adit_key", w.text())

        form = QFormLayout(self)
        form.setVerticalSpacing(ROW_SPACING)
        add_row(form, "Slater-Koster パラメータ", self.sk_set)
        self.sk_info.setObjectName("hint"); self.sk_info.setOpenExternalLinks(True)
        self.sk_browse = QPushButton("フォルダを選ぶ…"); self.sk_browse.setObjectName("link")
        self.sk_browse.setToolTip("Slater-Koster パラメータのセット (mio-1-1 など) を入れたフォルダを選びます。選んだ場所は環境設定 (sk_root) に保存されます")
        self.sk_browse.clicked.connect(self._browse_root)
        info_row = QVBoxLayout(); info_row.setSpacing(2); info_row.addWidget(self.sk_info); info_row.addWidget(self.sk_browse, 0, Qt.AlignmentFlag.AlignLeft)
        form.addRow(label(""), info_row)
        form.addRow(self.scc)
        add_row(form, "SCC の収束判定 (SccTolerance)", self.scc_tol)
        add_row(form, "SCC の反復上限 (MaxSccIterations)", self.max_scc)
        form.addRow(self.third)
        add_row(form, "分散力補正", self.dispersion)
        grid = QGridLayout(); grid.setHorizontalSpacing(8); grid.setVerticalSpacing(4)
        for i, (k, w) in enumerate(self.d3.items()):
            r, c = divmod(i, 2)
            grid.addWidget(QLabel(k), r, 2 * c); grid.addWidget(w, r, 2 * c + 1)
        grid.setColumnStretch(4, 1)
        add_row(form, "DFT-D3 係数 (BJ)", grid)
        add_row(form, "電子温度 [K]", self.temperature)
        add_row(form, "溶媒のパラメータファイル (GBSA)", self.solv_row)

        self.reload_sets(sk_root)
        self.sk_set.currentTextChanged.connect(self._on_set)
        self.solv_file.textChanged.connect(self._emit)
        for w in (self.scc, self.third):
            w.toggled.connect(self._emit)
        for w in (self.scc_tol, self.max_scc, self.temperature, *self.d3.values()):
            w.valueChanged.connect(self._emit)
        self.dispersion.currentTextChanged.connect(self._on_dispersion)
        self._on_dispersion(); self._on_set()

    def reload_sets(self, sk_root: str) -> None:
        self.sets = discover_sets(Path(sk_root).expanduser()) if sk_root else {}
        cur = self.sk_set.currentText()
        self.sk_set.clear(); self.sk_set.addItems(sorted(self.sets))
        if cur in self.sets:
            self.sk_set.setCurrentText(cur)
        if not self.sets:
            where = L(f"{sk_root} にセットがありません。", f"no parameter set in {sk_root}. ") if sk_root else L("置き場所がまだ指定されていません。", "the folder is not set yet. ")
            self.sk_info.setText(L(f"{where}DFTB+ の公式サイトの <a href=\"{SK_DOWNLOAD_URL}\">Parameters</a> から "
                                   "mio などを入手して展開し、「フォルダを選ぶ…」で指定してください",
                                   f"{where}Download e.g. mio from <a href=\"{SK_DOWNLOAD_URL}\">Parameters</a> "
                                   "on the DFTB+ website, unpack it, and pick the folder with \"Choose folder…\""))

    def _browse_root(self) -> None:
        root = choose_root(self, L("Slater-Koster パラメータの置き場所", "Folder of Slater-Koster parameter sets"), ".skf")
        if root:
            self.root_chosen.emit(root)

    def current_set(self) -> SKSet | None:
        return self.sets.get(self.sk_set.currentText())

    def method(self) -> DftbMethod:
        d3 = {k: w.value() for k, w in self.d3.items()} if self.dispersion.currentText() == "dftd3" else None
        return DftbMethod(sk_set=self.sk_set.currentText(), scc=self.scc.isChecked(), scc_tolerance=self.scc_tol.value(),
                      max_scc_iterations=self.max_scc.value(), third_order=self.third.isChecked(),
                      dispersion=self.dispersion.currentText(), d3_params=d3, filling_temperature=self.temperature.value(),
                      solvation_param_file=self.solv_file.text().strip(), seed=self._seed)

    def set_periodic(self, periodic: bool) -> None:
        on = not periodic or bool(self.solv_file.text().strip())
        self.solv_file.setEnabled(on); self.solv_row.browse.setEnabled(on)
        self.solv_file.setPlaceholderText(self._solv_ph if not periodic else L("周期系では使えません (分子のときだけ)", "not available for periodic systems (molecules only)"))

    def set_method(self, m: DftbMethod) -> None:
        self._seed = m.seed
        self.unshown = []
        if m.sk_set and self.sk_set.findText(m.sk_set) < 0:
            self.unshown.append(L(f"Slater-Koster パラメータ ({m.sk_set} → {self.sk_set.currentText() or '空欄'})",
                                  f"Slater-Koster parameters ({m.sk_set} → {self.sk_set.currentText() or 'empty'})"))
        self.sk_set.setCurrentText(m.sk_set); self.scc.setChecked(m.scc); self.scc_tol.setValue(m.scc_tolerance)
        self.max_scc.setValue(m.max_scc_iterations); self.third.setChecked(m.third_order)
        self.dispersion.setCurrentText(m.dispersion); self.temperature.setValue(m.filling_temperature)
        self.solv_file.setText(m.solvation_param_file)
        for k, w in self.d3.items():
            w.setValue((m.d3_params or {}).get(k, 0.0))
        self._emit()

    def _on_set(self, *_) -> None:
        s = self.current_set()
        if s is not None:
            extras = [n for n, ok in (("LICENSE/README", bool(s.doc_files()) and len(s.doc_files()) == 2),
                                      ("spinw.txt", s.spin_constants() is not None),
                                      ("Hubbard 微分 (README)", s.hubbard_derivs() is not None)) if ok]
            self.sk_info.setText(L(f"対応元素: {' '.join(s.elements)}   付属ファイル: {', '.join(extras) or 'なし'}", f"Elements: {' '.join(s.elements)}   included files: {', '.join(extras) or 'none'}"))
        self._emit()

    def _on_dispersion(self, *_) -> None:
        on = self.dispersion.currentText() == "dftd3"
        for w in self.d3.values():
            w.setEnabled(on)
        self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()


CODES = {"dftbplus": "DFTB+", "vasp": "VASP", "xtb": "xtb (GFN-xTB)", "espresso": "Quantum ESPRESSO (pw.x)", "orca": "ORCA",
         "cp2k": "CP2K", "lammps": "LAMMPS", "gromacs": "GROMACS", "mlip": "機械学習ポテンシャル (MACE・CHGNet)"}


class _CurrentPageStack(QWidget):

    def __init__(self, parent=None):
        super().__init__(parent)
        self._pages: list[QWidget] = []
        self._current = -1
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(0)

    def addWidget(self, w: QWidget) -> int:
        self._pages.append(w); self.layout().addWidget(w)
        if self._current < 0:
            self._current = 0
        w.setVisible(len(self._pages) - 1 == self._current)
        return len(self._pages) - 1

    def count(self) -> int:
        return len(self._pages)

    def widget(self, i: int) -> QWidget:
        return self._pages[i]

    def currentIndex(self) -> int:
        return self._current

    def currentWidget(self) -> QWidget | None:
        return self._pages[self._current] if 0 <= self._current < len(self._pages) else None

    def setCurrentIndex(self, i: int) -> None:
        self._current = i
        for k, w in enumerate(self._pages):
            w.setVisible(k == i)


class MethodPanel(QGroupBox):

    changed = Signal()

    def __init__(self, cfg: Config, parent: QWidget | None = None):
        super().__init__("計算手法", parent)
        self.setObjectName("method")
        self.unshown: list[str] = []
        self.code = QComboBox()
        self.code.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.code.setMinimumContentsLength(24)
        for k, v in CODES.items():
            self.code.addItem(icons.code_icon(k), v, k)
        self.dftb = DftbMethodPanel(cfg.sk_root)
        self.vasp = VaspMethodPanel(cfg)
        self.xtb = XtbMethodPanel()
        self.espresso = EspressoMethodPanel(cfg)
        self.orca = OrcaMethodPanel()
        self.cp2k = Cp2kMethodPanel(cfg)
        self.lammps = LammpsMethodPanel()
        self.gromacs = GromacsMethodPanel()
        self.mlip = MlipMethodPanel()
        self.stack = _CurrentPageStack()
        for w in self._panels().values():
            self.stack.addWidget(w)
        form = QFormLayout(); form.setVerticalSpacing(ROW_SPACING)
        add_row(form, "計算コード", self.code)
        lay = QVBoxLayout(self); lay.addLayout(form); lay.addWidget(self.stack)
        self.code.currentIndexChanged.connect(self._on_code)
        for w in self._panels().values():
            w.changed.connect(self.changed)
        self._attach_doc_boxes()

    def _panels(self) -> dict:
        return {"dftbplus": self.dftb, "vasp": self.vasp, "xtb": self.xtb, "espresso": self.espresso, "orca": self.orca,
                "cp2k": self.cp2k, "lammps": self.lammps, "gromacs": self.gromacs, "mlip": self.mlip}

    def _attach_doc_boxes(self) -> None:
        from adit.gui.prep_widgets import DocValueBox

        self.doc_boxes: dict[str, DocValueBox] = {}
        places = [("ecutwfc", self.espresso, self.espresso.ecutwfc, False), ("ecutrho", self.espresso, self.espresso.ecutrho, False),
                  ("encut", self.vasp, self.vasp.encut, False), ("cp2k_kind", self.cp2k, self.cp2k.kind_box, False),
                  ("vasp_incar", self.vasp, self.vasp.dipol, False),
                  ("d3", self.dftb, self.dftb.temperature, True), ("dftb_third", self.dftb, self.dftb.third, False),
                  ("dftb_sk", self.dftb, self.dftb.sk_set, False)]
        for key, panel, anchor, before in places:
            form = panel.layout()
            if not isinstance(form, QFormLayout):
                continue
            box = DocValueBox()
            box.use_requested.connect(self._use_doc_value)
            row, _role = form.getWidgetPosition(anchor)
            if row < 0:
                form.addRow(box)
            else:
                form.insertRow(row if before else row + 1, box)
            self.doc_boxes[key] = box

    def set_doc_values(self, groups: dict) -> None:
        for key, box in self.doc_boxes.items():
            box.set_lines(groups.get(key, []))

    def _use_doc_value(self, field: str, value: str) -> None:
        setters = {"ecutwfc": self.espresso.ecutwfc, "ecutrho": self.espresso.ecutrho, "encut": self.vasp.encut,
                   "prec": self.vasp.prec, "ediff": self.vasp.ediff, "nelmin": self.vasp.nelmin, "lreal": self.vasp.lreal,
                   **{f"d3_{k}": w for k, w in self.dftb.d3.items()}}
        w = setters.get(field)
        if w is None:
            return
        if isinstance(w, QComboBox):
            w.setCurrentText(value); return
        try: w.setValue(float(value))
        except ValueError: return

    @property
    def sk_set(self):
        return self.dftb.sk_set

    @property
    def sets(self):
        return self.dftb.sets

    def reload_sets(self, sk_root: str, cfg: Config | None = None) -> None:
        self.dftb.reload_sets(sk_root)
        if cfg is not None:
            self.espresso.reload_sets(cfg)
            self.cp2k.reload_data(cfg)

    def current_code(self) -> str:
        return self.code.currentData()

    def set_context(self, elements: list[str], profile_name: str, cfg: Config | None = None) -> None:
        self.vasp.set_context(elements, profile_name, cfg)
        self.espresso.set_context(elements, profile_name, cfg)
        self.cp2k.set_context(elements, profile_name, cfg)

    def set_periodic(self, periodic: bool) -> None:
        self.dftb.set_periodic(periodic)

    def method(self):
        return self._panels()[self.current_code()].method()

    def set_method(self, m) -> None:
        self.code.setCurrentIndex(list(CODES).index(m.code))
        panel = self._panels()[m.code]
        panel.set_method(m)
        self.unshown = list(getattr(panel, "unshown", []))

    def _on_code(self, *_) -> None:
        self.stack.setCurrentIndex(self.code.currentIndex())
        self.changed.emit()
