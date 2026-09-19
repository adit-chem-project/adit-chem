
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QCheckBox, QComboBox, QFormLayout, QLabel, QLineEdit, QSpinBox, QWidget

from adit.gui.style import ROW_SPACING
from adit.gui.widgets import add_row, label, narrow
from adit.lang import L
from adit.spec import MlipMethod

FAMILIES = [("", L("(選んでください)", "(choose one)")), ("mace_mp", "MACE-MP (mace_mp)"), ("mace_off", "MACE-OFF (mace_off)"),
            ("chgnet", "CHGNet (chgnet)")]
DTYPES = [("", L("(パッケージの既定)", "(package default)")), ("float32", "float32"), ("float64", "float64")]


class MlipMethodPanel(QWidget):
    changed = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.family = QComboBox()
        for k, v in FAMILIES:
            self.family.addItem(v, k)
        self.model = QLineEdit()
        self.model.setPlaceholderText("空欄ならパッケージの既定のモデル (例 small、medium-mpa-0、またはモデルのファイルのパス)")
        self.device = QComboBox(); self.device.setEditable(True)
        self.device.addItems(["", "cpu", "cuda", "mps"])
        self.device.lineEdit().setPlaceholderText("空欄ならパッケージの既定 (例 cpu、cuda)")
        self.dtype = QComboBox()
        for k, v in DTYPES:
            self.dtype.addItem(v, k)
        self.dispersion = QCheckBox("D3 を足す (mace_mp の dispersion)")
        self.dispersion.setProperty("adit_key", self.dispersion.text())
        self.seed = narrow(QSpinBox()); self.seed.setRange(0, 2 ** 31 - 1); self.seed.setValue(12345)
        self.pip = QLabel(""); self.pip.setObjectName("hint"); self.pip.setWordWrap(True); self.pip.setTextInteractionFlags(
            self.pip.textInteractionFlags())
        info = QLabel("分子でも周期系でも使えます。全電荷とスピン多重度は使いません (0 と 1 のまま)。k 点はありません。"
                      "計算は生成した run_mlip.py (ASE) が行い、ADIT 自身はモデルのパッケージを使いません")
        info.setObjectName("hint"); info.setWordWrap(True)

        form = QFormLayout(self); form.setVerticalSpacing(ROW_SPACING)
        add_row(form, "機械学習ポテンシャルの種類", self.family)
        add_row(form, "モデル", self.model)
        add_row(form, "計算に使うデバイス (device)", self.device)
        add_row(form, "数値の精度 (dtype)", self.dtype)
        add_row(form, "分散補正", self.dispersion)
        add_row(form, "乱数の種", self.seed)
        form.addRow(label(""), self.pip)
        form.addRow(label(""), info)

        self.family.currentIndexChanged.connect(self._on_family)
        self.model.textChanged.connect(self._emit)
        self.device.currentTextChanged.connect(self._emit)
        self.dtype.currentIndexChanged.connect(self._emit)
        self.dispersion.toggled.connect(self._emit)
        self.seed.valueChanged.connect(self._emit)
        self._on_family()

    def _on_family(self, *_) -> None:
        from adit.web.prep23 import mlip_pip_line

        fam = self.family.currentData() or ""
        self.pip.setText(mlip_pip_line(fam))
        self.dispersion.setEnabled(fam == "mace_mp" or self.dispersion.isChecked())
        self.dtype.setEnabled(fam != "chgnet" or bool(self.dtype.currentData()))
        self.model.setEnabled(fam != "chgnet" or bool(self.model.text().strip()))
        self._emit()

    def method(self) -> MlipMethod:
        return MlipMethod(model_family=self.family.currentData() or "", model=self.model.text().strip(),
                          device=self.device.currentText().strip(), dtype=self.dtype.currentData() or "",
                          dispersion=self.dispersion.isChecked(), seed=self.seed.value())

    def set_method(self, m: MlipMethod) -> None:
        self.family.setCurrentIndex(max(0, self.family.findData(m.model_family)))
        self.model.setText(m.model)
        self.device.setCurrentText(m.device)
        self.dtype.setCurrentIndex(max(0, self.dtype.findData(m.dtype)))
        self.dispersion.setChecked(m.dispersion)
        self.seed.setValue(m.seed)
        self._on_family()

    def _emit(self, *_) -> None:
        self.changed.emit()
