
from __future__ import annotations

from datetime import datetime
from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QPushButton, QSpinBox, QWidget

from adit.lang import L
from adit.config import Config
from adit.gui.style import ROW_SPACING
from adit.gui.style import NARROW_FIELD
from adit.gui import icons
from adit.gui.widgets import add_row, label, narrow
from adit.spec import Runtime

_LOCAL_JA = "この PC で bash submit.sh を実行 (計算ソフトがインストール済みで PATH にあること)"
_LOCAL_EN = "run bash submit.sh on this PC (the program must be installed and on PATH)"


class _DefaultDescriptions(dict):
    def get(self, key, default=None):
        if key in ("この機械で bash submit.sh (実行ファイルが PATH にあること)", _LOCAL_JA, _LOCAL_EN):
            return L(_LOCAL_JA, _LOCAL_EN)
        return default


DEFAULT_DESCRIPTIONS = _DefaultDescriptions()


class RuntimePanel(QGroupBox):
    changed = Signal()

    def __init__(self, cfg: Config, parent: QWidget | None = None):
        super().__init__("実行環境とリソース", parent)
        self.setObjectName("runtime")
        self.cfg = cfg
        self.unshown: list[str] = []
        self.profile = QComboBox()
        self.profile_info = QLabel("")
        self.nodes = narrow(QSpinBox()); self.nodes.setRange(1, 1000)
        self.ncpus = narrow(QSpinBox()); self.ncpus.setRange(1, 4096); self.ncpus.setValue(8)
        self.mpiprocs = narrow(QSpinBox()); self.mpiprocs.setRange(1, 4096)
        self.omp = narrow(QSpinBox()); self.omp.setRange(1, 4096); self.omp.setValue(8)
        self.walltime = narrow(QLineEdit("01:00:00"))
        self.job_name = QLineEdit("adit"); self.job_name.setMaximumWidth(2 * NARROW_FIELD)
        base = Path.home() / "adit_runs"
        base.mkdir(exist_ok=True)
        self.outdir = QLineEdit(str(base / datetime.now().strftime("run_%Y%m%d_%H%M%S")))
        self.browse = QPushButton("参照…")

        form = QFormLayout(self); self._form = form
        form.setVerticalSpacing(ROW_SPACING)
        add_row(form, "プロファイル", self.profile)
        self.profile_info.setObjectName("hint"); self.profile_info.setWordWrap(True)
        form.addRow(label(""), self.profile_info)
        add_row(form, "ノード数", self.nodes)
        add_row(form, "ノードあたりのコア数", self.ncpus)
        add_row(form, "ノードあたりの MPI プロセス数", self.mpiprocs)
        add_row(form, "OpenMP スレッド数", self.omp)
        add_row(form, "制限時間 (HH:MM:SS)", self.walltime)
        add_row(form, "ジョブ名", self.job_name)
        row = QHBoxLayout(); row.addWidget(self.outdir); row.addWidget(self.browse)
        add_row(form, "出力ディレクトリ", row)

        self.reload_profiles(cfg)
        self.profile.currentTextChanged.connect(self._on_profile)
        for w in (self.nodes, self.ncpus, self.mpiprocs, self.omp):
            w.valueChanged.connect(self._emit)
        for w in (self.walltime, self.job_name, self.outdir):
            w.textChanged.connect(self._emit)
        self.browse.clicked.connect(self._browse)
        self._on_profile()

    def reload_profiles(self, cfg: Config) -> None:
        self.cfg = cfg
        self.profile.clear()
        for name in sorted(cfg.profiles):
            self.profile.addItem(icons.icon(f"profile_{cfg.profiles[name].kind}"), name)
        if cfg.default_profile in cfg.profiles:
            self.profile.setCurrentText(cfg.default_profile)

    def runtime(self) -> Runtime:
        return Runtime(profile=self.profile.currentText(), nodes=self.nodes.value(), ncpus=self.ncpus.value(),
                       mpiprocs=self.mpiprocs.value(), omp_threads=self.omp.value(),
                       walltime=self.walltime.text().strip(), job_name=self.job_name.text().strip())

    def output_dir(self) -> str:
        return self.outdir.text().strip()

    def set_runtime(self, r: Runtime) -> None:
        self.unshown = []
        if self.profile.findText(r.profile) < 0:
            self.unshown.append(L(f"プロファイル ({r.profile} → {self.profile.currentText()})",
                                  f"profile ({r.profile} → {self.profile.currentText()})"))
        self.profile.setCurrentText(r.profile); self.nodes.setValue(r.nodes); self.ncpus.setValue(r.ncpus)
        self.mpiprocs.setValue(r.mpiprocs); self.omp.setValue(r.omp_threads)
        self.walltime.setText(r.walltime); self.job_name.setText(r.job_name)
        self._emit()

    def _on_profile(self, *_) -> None:
        p = self.cfg.profiles.get(self.profile.currentText())
        if p is None:
            self.profile_info.setText(L("プロファイルがありません (cluster.toml を確認してください)", "no profile (check cluster.toml)"))
        else:
            how = f"{p.submit} submit.sh"
            desc = DEFAULT_DESCRIPTIONS.get(p.description or "", p.description or "")
            self.profile_info.setText(L(f"{p.kind}: {desc}  実行は手動で: {how}", f"{p.kind}: {desc}  run manually: {how}"))
            for w in (self.nodes, self.walltime, self.job_name):
                self._form.setRowVisible(w, p.kind != "direct")
        self._emit()

    def _browse(self) -> None:
        d = QFileDialog.getExistingDirectory(self, L("出力先の親ディレクトリ", "Parent directory of the output"), str(Path(self.outdir.text()).parent))
        if d:
            self.outdir.setText(str(Path(d) / Path(self.outdir.text()).name))

    def _emit(self, *_) -> None:
        self.changed.emit()
