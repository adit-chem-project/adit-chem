
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QCheckBox, QComboBox, QFormLayout, QLabel, QLineEdit, QPlainTextEdit, QSpinBox, QWidget

from adit.codes.orca import solvent_names
from adit.gui.panels.xtb_panel import fill_solvents
from adit.gui.style import ROW_SPACING
from adit.gui.widgets import add_row, label, narrow
from adit.spec import OrcaMethod

SOLVATION = {"none": "なし", "cpcm": "CPCM", "smd": "SMD"}


class OrcaMethodPanel(QWidget):
    changed = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self._goat = False
        self._docker_guest_file = ""
        self._docker_assume_neutral_singlet = False
        self.method_kw = QComboBox(); self.method_kw.setEditable(True)
        self.method_kw.addItems(["HF", "PBE", "B3LYP", "PBE0", "TPSS", "wB97X-D3", "MP2", "DLPNO-CCSD(T)"])
        self.basis = QComboBox(); self.basis.setEditable(True)
        self.basis.addItems(["def2-SVP", "def2-TZVP", "def2-TZVPP", "def2-QZVPP", "cc-pVDZ", "cc-pVTZ", "ma-def2-SVP"])
        self.scf_conv = QComboBox(); self.scf_conv.addItems(["NormalSCF", "TightSCF", "LooseSCF"])
        self.scf_maxiter = narrow(QSpinBox()); self.scf_maxiter.setRange(1, 10000); self.scf_maxiter.setValue(125)
        self.maxcore = narrow(QSpinBox()); self.maxcore.setRange(0, 1000000); self.maxcore.setSingleStep(500)
        self.solvation = QComboBox()
        for k, v in SOLVATION.items():
            self.solvation.addItem(v, k)
        self.solvent = QComboBox()
        self.extra_kw = QLineEdit(); self.extra_kw.setPlaceholderText("! 行に追加するキーワード (例 D3BJ RIJCOSX)")
        self.extra_blocks = QPlainTextEdit(); self.extra_blocks.setMaximumHeight(80)
        self.extra_blocks.setPlaceholderText("%ブロックをそのまま追加 (例 %basis newGTO ... end)")
        self.ts_search = QCheckBox("! OptTS にする (計算の種類は構造最適化)")
        self.ts_calc_hess = QCheckBox("%geom Calc_Hess true"); self.ts_calc_hess.setChecked(True)
        self.ts_recalc = narrow(QSpinBox()); self.ts_recalc.setRange(0, 1000)
        self.ts_freq = QCheckBox("! Freq を足す"); self.ts_freq.setChecked(True)
        self.irc = QCheckBox("! Freq IRC にする (計算の種類は一点計算)")
        self.irc_max_iter = narrow(QSpinBox()); self.irc_max_iter.setRange(0, 100000)
        self.irc_direction = QComboBox()
        for k, v in (("both", "both"), ("forward", "forward"), ("backward", "backward"), ("down", "down")):
            self.irc_direction.addItem(v, k)
        for w in (self.ts_search, self.ts_calc_hess, self.ts_freq, self.irc):
            w.setProperty("adit_key", w.text())
        info = QLabel("分子系のみ。電荷とスピン多重度は構造パネルの値を使います。並列数は実行環境の MPI プロセス数 (%pal nprocs)。ORCA 本体は登録して入手し、実行パスは環境設定に書きます")
        info.setObjectName("hint"); info.setWordWrap(True)
        form = QFormLayout(self); form.setVerticalSpacing(ROW_SPACING); self._form = form
        add_row(form, "計算手法 (! 行)", self.method_kw)
        add_row(form, "基底関数", self.basis)
        add_row(form, "SCF の収束判定", self.scf_conv)
        add_row(form, "SCF の最大反復回数", self.scf_maxiter)
        add_row(form, "%maxcore [MB] (0 = 指定しない)", self.maxcore)
        add_row(form, "溶媒モデル (CPCM / SMD)", self.solvation)
        add_row(form, "溶媒", self.solvent)
        add_row(form, "追加のキーワード (! 行)", self.extra_kw)
        add_row(form, "追加の %ブロック", self.extra_blocks)
        add_row(form, "遷移状態の探索 (OptTS)", self.ts_search)
        add_row(form, "最初にヘシアンを計算 (Calc_Hess)", self.ts_calc_hess)
        add_row(form, "ヘシアンを計算し直す間隔 (Recalc_Hess)", self.ts_recalc)
        add_row(form, "最後に振動数を計算 (Freq)", self.ts_freq)
        add_row(form, "反応座標をたどる (IRC)", self.irc)
        add_row(form, "IRC の反復の上限", self.irc_max_iter)
        add_row(form, "IRC の向き", self.irc_direction)
        self.ts_note = QLabel("遷移状態の探索は計算の種類が構造最適化のとき、IRC は一点計算のときに書かれます。"
                             "両方を続けて実行するなら「実行」→「遷移状態と IRC」で段階に分けて生成します")
        self.ts_note.setObjectName("hint"); self.ts_note.setWordWrap(True)
        form.addRow(label(""), self.ts_note)
        form.addRow(label(""), info)
        for w in (self.method_kw, self.basis, self.scf_conv):
            w.currentTextChanged.connect(self._emit)
        for w in (self.scf_maxiter, self.maxcore):
            w.valueChanged.connect(self._emit)
        self.solvation.currentIndexChanged.connect(self._refill)
        self.solvent.currentIndexChanged.connect(self._emit)
        self.extra_kw.textChanged.connect(self._emit); self.extra_blocks.textChanged.connect(self._emit)
        for w in (self.ts_search, self.ts_calc_hess, self.ts_freq, self.irc):
            w.toggled.connect(self._on_ts)
        for w in (self.ts_recalc, self.irc_max_iter):
            w.valueChanged.connect(self._emit)
        self.irc_direction.currentIndexChanged.connect(self._emit)
        self._refill(); self._on_ts()

    def _on_ts(self, *_) -> None:
        for w in (self.ts_calc_hess, self.ts_recalc, self.ts_freq):
            self._form.setRowVisible(w, self.ts_search.isChecked())
        for w in (self.irc_max_iter, self.irc_direction):
            self._form.setRowVisible(w, self.irc.isChecked())
        self._emit()

    def _refill(self, *_, keep: str | None = None) -> None:
        model = self.solvation.currentData()
        cur = keep if keep is not None else (self.solvent.currentData() or "")
        fill_solvents(self.solvent, solvent_names(model) if model != "none" else [], cur)
        self._form.setRowVisible(self.solvent, model != "none" or bool(cur))
        self._emit()

    def method(self) -> OrcaMethod:
        return OrcaMethod(method=self.method_kw.currentText().strip(), basis=self.basis.currentText().strip(),
                          scf_convergence=self.scf_conv.currentText(), scf_maxiter=self.scf_maxiter.value(), maxcore_mb=self.maxcore.value(),
                          extra_keywords=self.extra_kw.text(), extra_blocks=self.extra_blocks.toPlainText(),
                          solvation=self.solvation.currentData(), solvent=self.solvent.currentData() or "",
                          ts_search=self.ts_search.isChecked(), ts_calc_hess=self.ts_calc_hess.isChecked(),
                          ts_recalc_hess=self.ts_recalc.value(), ts_freq=self.ts_freq.isChecked(),
                          irc=self.irc.isChecked(), irc_max_iter=self.irc_max_iter.value(),
                          irc_direction=self.irc_direction.currentData(),
                          goat=self._goat, docker_guest_file=self._docker_guest_file,
                          docker_assume_neutral_singlet=self._docker_assume_neutral_singlet)

    def set_method(self, m: OrcaMethod) -> None:
        self._goat = m.goat
        self._docker_guest_file = m.docker_guest_file
        self._docker_assume_neutral_singlet = m.docker_assume_neutral_singlet
        self.method_kw.setCurrentText(m.method); self.basis.setCurrentText(m.basis); self.scf_conv.setCurrentText(m.scf_convergence)
        self.scf_maxiter.setValue(m.scf_maxiter); self.maxcore.setValue(m.maxcore_mb)
        self.extra_kw.setText(m.extra_keywords); self.extra_blocks.setPlainText(m.extra_blocks)
        for w, v in ((self.ts_search, m.ts_search), (self.ts_calc_hess, m.ts_calc_hess), (self.ts_freq, m.ts_freq), (self.irc, m.irc)):
            w.blockSignals(True); w.setChecked(v); w.blockSignals(False)
        self.ts_recalc.setValue(m.ts_recalc_hess); self.irc_max_iter.setValue(m.irc_max_iter)
        self.irc_direction.setCurrentIndex(max(0, self.irc_direction.findData(m.irc_direction)))
        self._on_ts()
        self.solvation.blockSignals(True); self.solvation.setCurrentIndex(max(0, self.solvation.findData(m.solvation))); self.solvation.blockSignals(False)
        self._refill(keep=m.solvent)

    def _emit(self, *_) -> None:
        self.changed.emit()
