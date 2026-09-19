
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QFormLayout, QGroupBox, QLabel, QLineEdit, QSpinBox, QWidget

from adit.gui import icons
from adit.gui.style import ROW_SPACING
from adit.gui.widgets import SciDoubleSpinBox, add_row, narrow
from adit.lang import L
from adit.spec import BandSettings, MDSettings, Task

TYPES = {"single_point": "一点計算", "geometry_optimization": "構造最適化", "molecular_dynamics": "分子動力学", "vibrations": "振動解析", "band_structure": "バンド計算"}
OPTIMIZERS = ["Rational", "LBFGS", "FIRE", "SteepestDescent"]
RELAX_CELL = {"no": "固定", "shape_and_volume": "形と体積", "volume_only": "体積だけ"}
ENSEMBLES = ["NVT", "NVE", "NPT"]
THERMOSTATS = {"berendsen": "Berendsen", "andersen": "Andersen", "nose_hoover": "Nosé-Hoover", "langevin": "Langevin", "csvr": "CSVR (速度再スケール)"}


class TaskPanel(QGroupBox):
    changed = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__("計算の種類", parent)
        self.setObjectName("task")
        self.type = QComboBox()
        for k, v in TYPES.items():
            self.type.addItem(icons.task_icon(k), v, k)
        self.optimizer = QComboBox(); self.optimizer.addItems(OPTIMIZERS)
        self.max_steps = narrow(QSpinBox()); self.max_steps.setRange(-1, 1000000); self.max_steps.setValue(200)
        self.force_tol = narrow(SciDoubleSpinBox(1e-8, 10, Task().force_tolerance_ev_per_ang, 1e-3))
        self.relax_cell = QComboBox()
        for k, v in RELAX_CELL.items():
            self.relax_cell.addItem(v, k)
        self.ensemble = QComboBox(); self.ensemble.addItems(ENSEMBLES)
        self.thermostat = QComboBox()
        for k, v in THERMOSTATS.items():
            self.thermostat.addItem(v, k)
        self.temperature = narrow(SciDoubleSpinBox(-1e100, 1e100, 300.0, 50))
        self.timestep = narrow(SciDoubleSpinBox(-1e100, 1e100, 1.0, 0.5))
        self.md_steps = narrow(QSpinBox()); self.md_steps.setRange(-2_000_000_000, 2_000_000_000); self.md_steps.setValue(1000)
        self.dump = narrow(QSpinBox()); self.dump.setRange(1, 1000000); self.dump.setValue(10)
        self.coupling = narrow(SciDoubleSpinBox(1e-3, 1e6, 100.0, 50))
        self.pressure = narrow(SciDoubleSpinBox(0, 1e7, 1.0, 1))
        self.barostat_time = narrow(SciDoubleSpinBox(1e-3, 1e7, 1000.0, 100))
        self.band_path = QLineEdit(); self.band_path.setPlaceholderText("空欄なら格子の標準経路 (例 GXWKGLUWLK,UX)")
        self.band_npoints = narrow(QSpinBox()); self.band_npoints.setRange(2, 10000)
        self.band_empty = narrow(QSpinBox()); self.band_empty.setRange(0, 1000); self.band_empty.setValue(4)

        form = QFormLayout(self); form.setVerticalSpacing(ROW_SPACING); self._form = form
        add_row(form, "種類", self.type)
        self.vibration_note = QLabel(L("入力した構造をそのまま使います。最適化の出力を使うには「前の計算の続き」を選びます。",
                                      "Uses the input structure as-is. Choose Continue a previous calculation to use an optimization result."))
        self.vibration_note.setObjectName("hint"); self.vibration_note.setWordWrap(True)
        form.addRow(self.vibration_note)
        add_row(form, "最適化アルゴリズム", self.optimizer)
        add_row(form, "最大ステップ数 (MaxSteps)", self.max_steps)
        add_row(form, "力の収束判定 [eV/Å]", self.force_tol)
        add_row(form, "セルの緩和 (周期系)", self.relax_cell)
        add_row(form, "アンサンブル", self.ensemble)
        add_row(form, "熱浴", self.thermostat)
        add_row(form, "温度 [K]", self.temperature)
        add_row(form, "時間刻み [fs]", self.timestep)
        add_row(form, "MD ステップ数", self.md_steps)
        add_row(form, "軌跡の出力間隔 [ステップ]", self.dump)
        add_row(form, "熱浴の緩和時間 [fs]", self.coupling)
        add_row(form, "圧力 [bar] (NPT)", self.pressure)
        add_row(form, "圧力浴の緩和時間 [fs] (NPT)", self.barostat_time)
        add_row(form, "k 点の経路 (バンド)", self.band_path)
        add_row(form, "経路上の k 点数", self.band_npoints)
        add_row(form, "空のバンド数 (pw.x)", self.band_empty)
        self._periodic = False
        self.type.currentIndexChanged.connect(self._on_type)
        self.ensemble.currentTextChanged.connect(self._on_type)
        for w in (self.optimizer, self.relax_cell, self.thermostat):
            w.currentIndexChanged.connect(self._emit)
        for w in (self.max_steps, self.force_tol, self.temperature, self.timestep, self.md_steps, self.dump, self.coupling, self.pressure, self.barostat_time,
                  self.band_npoints, self.band_empty):
            w.valueChanged.connect(self._emit)
        self.band_path.textChanged.connect(self._emit)
        self._on_type()

    def task(self) -> Task:
        md = MDSettings(ensemble=self.ensemble.currentText(), thermostat=self.thermostat.currentData(), temperature_k=self.temperature.value(),
                        timestep_fs=self.timestep.value(), steps=self.md_steps.value(), dump_interval=self.dump.value(),
                        coupling_time_fs=self.coupling.value(), pressure_bar=self.pressure.value(), barostat_time_fs=self.barostat_time.value())
        bands = BandSettings(path=self.band_path.text().strip(), npoints=self.band_npoints.value(), empty_bands=self.band_empty.value())
        return Task(type=self.type.currentData(), optimizer=self.optimizer.currentText(), max_steps=self.max_steps.value(),
                    force_tolerance_ev_per_ang=self.force_tol.value(), relax_cell=self.relax_cell.currentData(), md=md, bands=bands)

    def set_periodic(self, periodic: bool) -> None:
        self._periodic = periodic
        self._on_type()

    def set_task(self, t: Task) -> None:
        self.type.setCurrentIndex(list(TYPES).index(t.type)); self.optimizer.setCurrentText(t.optimizer)
        self.max_steps.setValue(t.max_steps); self.force_tol.setValue(t.force_tolerance_ev_per_ang)
        self.relax_cell.setCurrentIndex(list(RELAX_CELL).index(t.relax_cell))
        md = t.md
        self.ensemble.setCurrentText(md.ensemble); self.thermostat.setCurrentIndex(list(THERMOSTATS).index(md.thermostat))
        self.temperature.setValue(md.temperature_k); self.timestep.setValue(md.timestep_fs); self.md_steps.setValue(md.steps)
        self.dump.setValue(md.dump_interval); self.coupling.setValue(md.coupling_time_fs); self.pressure.setValue(md.pressure_bar)
        self.barostat_time.setValue(md.barostat_time_fs)
        self.band_path.setText(t.bands.path); self.band_npoints.setValue(t.bands.npoints); self.band_empty.setValue(t.bands.empty_bands)
        self._emit()

    def _on_type(self, *_) -> None:
        kind = self.type.currentData()
        opt, md = kind == "geometry_optimization", kind == "molecular_dynamics"
        for w in (self.optimizer, self.max_steps, self.force_tol):
            self._form.setRowVisible(w, opt)
        self._form.setRowVisible(self.relax_cell, opt); self.relax_cell.setEnabled(self._periodic)
        ens = self.ensemble.currentText()
        for w in (self.ensemble, self.temperature, self.timestep, self.md_steps, self.dump):
            self._form.setRowVisible(w, md)
        for w in (self.thermostat, self.coupling):
            self._form.setRowVisible(w, md and ens != "NVE")
        for w in (self.pressure, self.barostat_time):
            self._form.setRowVisible(w, md and ens == "NPT")
        band = kind == "band_structure"
        self.vibration_note.setVisible(kind == "vibrations")
        for w in (self.band_path, self.band_npoints, self.band_empty):
            self._form.setRowVisible(w, band)
        self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()
