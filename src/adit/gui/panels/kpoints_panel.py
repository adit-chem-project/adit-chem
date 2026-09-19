
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QFormLayout, QGroupBox, QHBoxLayout, QLabel, QSpinBox, QWidget

from adit.gui import icons
from adit.gui.style import ROW_SPACING
from adit.gui.widgets import SciDoubleSpinBox, add_row, narrow
from adit.spec import KPoints

MODES = {"gamma": "Γ 点のみ", "mesh": "メッシュを指定", "density": "密度から自動"}


class KPointsPanel(QGroupBox):
    changed = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__("k 点 (周期系)", parent)
        self.setObjectName("kpoints")
        self.mode = QComboBox()
        for k, v in MODES.items():
            self.mode.addItem(icons.icon(f"kpoints_{k}"), v, k)
        self.mesh = [QSpinBox() for _ in range(3)]
        for w in self.mesh:
            w.setRange(1, 999); w.setValue(1); w.setFixedWidth(64)
        self.shift = QComboBox(); self.shift.addItems(["0", "0.5"])
        self.density = narrow(SciDoubleSpinBox(0.0, 1000, 0.0, 0.5))
        self.info = QLabel(""); self.info.setObjectName("hint")
        self._shift: tuple[float, float, float] | None = None   # the loaded 3 components; the field shows one value
        self._shift_shown = ""

        form = QFormLayout(self); form.setVerticalSpacing(ROW_SPACING); self._form = form
        add_row(form, "サンプリング方法", self.mode)
        self.mesh_row = QWidget(); self.mesh_row.setObjectName("rowbox"); row = QHBoxLayout(self.mesh_row); row.setContentsMargins(0, 0, 0, 0)
        for w in self.mesh:
            row.addWidget(w)
        row.addWidget(QLabel("シフト")); row.addWidget(self.shift); row.addStretch()
        add_row(form, "メッシュ (n1 n2 n3)", self.mesh_row)
        add_row(form, "k 点密度 [点/Å⁻¹]", self.density)
        form.addRow(self.info)
        self.mode.currentIndexChanged.connect(self._on_mode)
        for w in (*self.mesh, self.density):
            w.valueChanged.connect(self._emit)
        self.shift.currentTextChanged.connect(self._emit)
        self._on_mode()

    def set_periodic(self, periodic: bool) -> None:
        self.setEnabled(periodic)
        self.info.setText("" if periodic else "分子系 (非周期) では使いません")

    def kpoints(self) -> KPoints:
        s = float(self.shift.currentText())
        shift = self._shift if self._shift is not None and self.shift.currentText() == self._shift_shown else (s, s, s)
        return KPoints(mode=self.mode.currentData(), mesh=tuple(w.value() for w in self.mesh), shift=shift, density=self.density.value())

    def set_kpoints(self, kp: KPoints) -> None:
        self.mode.setCurrentIndex(list(MODES).index(kp.mode))
        for w, v in zip(self.mesh, kp.mesh):
            w.setValue(v)
        self._shift, self._shift_shown = tuple(kp.shift), "0.5" if kp.shift[0] == 0.5 else "0"
        self.shift.setCurrentText(self._shift_shown); self.density.setValue(kp.density)
        self._emit()

    def _on_mode(self, *_) -> None:
        m = self.mode.currentData()
        self._form.setRowVisible(self.mesh_row, m == "mesh")
        self._form.setRowVisible(self.density, m == "density")
        self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()
