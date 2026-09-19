
from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import QComboBox, QDoubleSpinBox, QHBoxLayout, QLabel, QPushButton, QSlider, QVBoxLayout, QWidget

from adit.gui.flow_layout import FlowLayout
from adit.gui.viewer3d import Viewer3D
from adit.lang import L

PHASES_PER_PERIOD = 24
TICK_MS = 60


class SeriesCanvas(QWidget):
    """A small matplotlib figure with a vertical marker that follows the frame; clicking jumps to that x."""

    clicked = Signal(float)

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
        from matplotlib.figure import Figure

        from adit.analysis.report import use_cjk_font

        use_cjk_font()
        self.figure = Figure(figsize=(6.0, 2.3), dpi=90)
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.canvas.setMinimumHeight(150)
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.addWidget(self.canvas)
        self._axes = []
        self._markers = []
        self.canvas.mpl_connect("button_press_event", self._on_press)

    def set_series(self, x, series: list[tuple[str, np.ndarray]], xlabel: str, stems: bool = False) -> None:
        self.figure.clear(); self._axes = []; self._markers = []
        series = [(lab, np.asarray(y, dtype=float)) for lab, y in series if len(y) == len(x)]
        if not series:
            self.canvas.draw_idle(); return
        n = len(series)
        for k, (lab, y) in enumerate(series):
            ax = self.figure.add_subplot(n, 1, k + 1, sharex=self._axes[0] if self._axes else None)
            if stems:
                ax.vlines(x, 0, y, color="#404048", lw=1.0)
            else:
                ax.plot(x, y, color="#404048", lw=1.0)
            ax.set_ylabel(lab, fontsize=8)
            ax.tick_params(labelsize=8)
            ax.grid(True, lw=0.4, alpha=0.5)
            self._markers.append(ax.axvline(x[0], color="#0A7AFF", lw=1.4))
            self._axes.append(ax)
        self._axes[-1].set_xlabel(xlabel, fontsize=8)
        for ax in self._axes[:-1]:
            ax.tick_params(labelbottom=False)
        self.figure.subplots_adjust(left=0.13, right=0.98, top=0.95, bottom=0.28 if n == 1 else 0.16, hspace=0.15)
        self.canvas.draw_idle()

    def set_marker(self, xv: float) -> None:
        for m in self._markers:
            m.set_xdata([xv, xv])
        self.canvas.draw_idle()

    def _on_press(self, ev) -> None:
        if ev.inaxes is None or ev.xdata is None:
            return
        self.clicked.emit(float(ev.xdata))


class PlaybackPanel(QWidget):
    """Plays a trajectory (MD frames, optimization steps) or a vibrational mode in a 3D view, with a linked series plot."""

    frameChanged = Signal(int)

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.viewer = Viewer3D(); self.viewer.setMinimumHeight(260)
        self.what = QComboBox()
        self.what.addItem(L("軌跡", "Trajectory"), "frames"); self.what.addItem(L("振動モード", "Vibrational mode"), "modes")
        self.mode = QComboBox()
        self.amplitude = QDoubleSpinBox(); self.amplitude.setRange(0.01, 5.0); self.amplitude.setSingleStep(0.1); self.amplitude.setValue(0.5)
        self.amplitude.setDecimals(2); self.amplitude.setSuffix(" Å")
        self.btn_first = QPushButton("|<"); self.btn_prev = QPushButton("<"); self.btn_play = QPushButton(L("再生", "Play"))
        self.btn_next = QPushButton(">"); self.btn_last = QPushButton(">|")
        for b in (self.btn_first, self.btn_prev, self.btn_next, self.btn_last):
            b.setFixedWidth(40)
        self.btn_play.setFixedWidth(64)
        self.btn_first.setToolTip(L("最初のフレーム", "First frame")); self.btn_last.setToolTip(L("最後のフレーム", "Last frame"))
        self.btn_prev.setToolTip(L("1 つ前 (コマ送り)", "Previous frame")); self.btn_next.setToolTip(L("1 つ先 (コマ送り)", "Next frame"))
        self.btn_play.setToolTip(L("再生 / 停止", "Play / pause"))
        self.slider = QSlider(Qt.Orientation.Horizontal); self.slider.setRange(0, 0)
        self.position = QLabel(""); self.position.setObjectName("hint")
        self.note = QLabel(""); self.note.setObjectName("hint"); self.note.setWordWrap(True)
        self.series = SeriesCanvas()
        self.lbl_what = QLabel(L("再生するもの", "Play")); self.lbl_mode = QLabel(L("モード", "Mode")); self.lbl_amp = QLabel(L("振幅 (いちばん動く原子)", "Amplitude (largest atom motion)"))
        top = FlowLayout()
        for w in (self.lbl_what, self.what, self.lbl_mode, self.mode, self.lbl_amp, self.amplitude):
            top.addWidget(w)
        row = QHBoxLayout(); row.setSpacing(4)
        for b in (self.btn_first, self.btn_prev, self.btn_play, self.btn_next, self.btn_last):
            row.addWidget(b)
        row.addWidget(self.slider, 1); row.addWidget(self.position)
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(6)
        lay.addLayout(top); lay.addWidget(self.viewer, 1); lay.addLayout(row); lay.addWidget(self.series); lay.addWidget(self.note)
        self._timer = QTimer(self); self._timer.setInterval(TICK_MS); self._timer.timeout.connect(self._tick)
        self._positions: np.ndarray | None = None       # (n_frames, n_atoms, 3)
        self._x: np.ndarray | None = None               # x value per frame for the series plot
        self._x_label = ""
        self._times_fs: np.ndarray | None = None
        self._modes = None
        self._series_ratio = 1.0
        self.frame = 0
        self.btn_first.clicked.connect(lambda: self.go(0)); self.btn_last.clicked.connect(lambda: self.go(self.count() - 1))
        self.btn_prev.clicked.connect(lambda: self.go(self.frame - 1)); self.btn_next.clicked.connect(lambda: self.go(self.frame + 1))
        self.btn_play.clicked.connect(self.toggle)
        self.slider.valueChanged.connect(self._on_slider)
        self.what.currentIndexChanged.connect(self._on_what)
        self.mode.currentIndexChanged.connect(lambda *_: self._show_modes())
        self.amplitude.valueChanged.connect(lambda *_: self._apply(self.frame))
        self.series.clicked.connect(self._on_series_click)
        self.clear()

    # ---- loading ----
    def clear(self) -> None:
        self.pause()
        self._positions = None; self._x = None; self._times_fs = None; self._modes = None
        self.viewer.set_atoms(None); self.slider.setRange(0, 0); self.position.setText("")
        self.what.setVisible(False); self.lbl_what.setVisible(False)
        self._set_modes_visible(False)
        self.note.setText("")

    def set_dark(self, dark: bool) -> None:
        self.viewer.dark = dark; self.viewer.update()

    def load_run(self, run_dir: Path | str, *, skip: int = 0, stride: int = 1, memory_mb: float | None = None) -> bool:
        """Read a run directory: frames for playback and, when available, vibrational modes. Returns True if anything can be played."""
        from adit.analysis.modes import read_modes, unsupported_reason
        from adit.analysis.readers import detect_code, load_run
        from adit.analysis.trajectory import MEMORY_BUDGET_MB, Trajectory, TrajectoryTooLarge, check_budget

        self.clear()
        d = Path(run_dir)
        notes: list[str] = []
        try:
            data = load_run(d)
        except Exception as ex:
            self.note.setText(L(f"再生できません: {ex}", f"cannot play: {ex}")); return False
        n_total = len(data.frames)
        have_frames = False
        if n_total >= 2:
            stride = max(1, int(stride)); skip = int(skip) if n_total > skip else 0
            sel = data.frames[skip::stride] if (skip or stride > 1) else data.frames
            nat = data.frames.first_natoms() if isinstance(data.frames, Trajectory) else len(data.frames[0])
            try:
                check_budget(sel, nat, memory_mb or MEMORY_BUDGET_MB, what=L("再生", "playback"))
            except TrajectoryTooLarge as ex:
                stride = int(ex.suggested_stride)
                sel = data.frames[skip::stride]
                notes.append(L(f"軌跡が大きいので {stride} フレームに 1 回に間引いて読みました", f"the trajectory is large, so every {stride}th frame was read"))
            frames = list(sel)
            if frames:
                pos = np.array([f.get_positions() for f in frames], dtype=float)
                idx = np.arange(skip, skip + stride * len(frames), stride)[: len(frames)]
                times = None
                if data.times_fs and len(data.times_fs) == n_total:
                    times = np.asarray(data.times_fs, dtype=float)[idx]
                elif data.frame_dt_fs:
                    times = idx * float(data.frame_dt_fs)
                md = bool(data.temperatures_k)
                series = []
                e = np.asarray(data.energies_ev, dtype=float) if data.energies_ev else np.zeros(0)
                if len(e) == n_total:
                    series.append((L("エネルギー [eV]", "energy [eV]"), e[idx]))
                    if md and data.temperatures_k and len(data.temperatures_k) == n_total:
                        series.append((L("温度 [K]", "temperature [K]"), np.asarray(data.temperatures_k, dtype=float)[idx]))
                    x = times if times is not None else idx.astype(float)
                    self._series_ratio = 1.0
                elif len(e) >= 2:
                    # the energy list and the frames have different counts: place the marker proportionally
                    x = np.linspace(0, len(e) - 1, len(frames))
                    series.append((L("エネルギー [eV]", "energy [eV]"), np.interp(x, np.arange(len(e)), e)))
                    notes.append(L(f"エネルギーの点数 ({len(e)}) とフレーム数 ({n_total}) が違うので、縦線の位置は比例で合わせています",
                                   f"the energy list ({len(e)} points) and the frames ({n_total}) differ in count, so the marker is placed proportionally"))
                    times = None
                else:
                    x = idx.astype(float)
                xlabel = L("時刻 [fs]", "time [fs]") if times is not None and x is times else (L("ステップ", "step") if not md else L("フレーム", "frame"))
                self.set_frames(frames[0], pos, x, xlabel, series, times, source=data.frame_source or "")
                have_frames = True
        else:
            notes.append(L("フレームが 1 つしかないので軌跡は再生できません", "only one frame, so there is no trajectory to play"))
        modes = None
        try:
            modes = read_modes(d, data.code)
        except Exception as ex:
            notes.append(L(f"振動モードを読めません: {ex}", f"cannot read the vibrational modes: {ex}"))
        if modes is not None and len(modes):
            self.set_modes(modes)
            notes += modes.notes
        elif data.frequencies_cm1:
            notes.append(unsupported_reason(detect_code(d)))
        both = have_frames and modes is not None
        self.what.setVisible(both); self.lbl_what.setVisible(both)
        if not have_frames and modes is None:
            self.note.setText("\n".join(notes)); return False
        self.what.blockSignals(True); self.what.setCurrentIndex(0 if have_frames else 1); self.what.blockSignals(False)
        self._on_what()
        self.note.setText("\n".join(notes))
        return True

    def set_frames(self, atoms0, positions: np.ndarray, x, xlabel: str, series: list[tuple[str, np.ndarray]] | None = None,
                   times_fs=None, source: str = "") -> None:
        self._positions = np.asarray(positions, dtype=float)
        self._x = np.asarray(x, dtype=float); self._x_label = xlabel
        self._times_fs = None if times_fs is None else np.asarray(times_fs, dtype=float)
        self._frames_atoms = atoms0.copy()
        self._frames_series = list(series or [])
        self._frames_source = source

    def set_modes(self, modes) -> None:
        self._modes = modes
        self.mode.blockSignals(True); self.mode.clear()
        for k in range(len(modes)):
            self.mode.addItem(modes.label(k), k)
        self.mode.blockSignals(False)

    def _set_modes_visible(self, on: bool) -> None:
        for w in (self.lbl_mode, self.mode, self.lbl_amp, self.amplitude):
            w.setVisible(on)

    def playing_modes(self) -> bool:
        return self.what.currentData() == "modes"

    def _on_what(self, *_) -> None:
        self.pause()
        if self.playing_modes() and self._modes is not None:
            self._set_modes_visible(True)
            self._show_modes()
        elif self._positions is not None:
            self._set_modes_visible(False)
            self.viewer.set_atoms(self._frames_atoms.copy())
            self.slider.blockSignals(True); self.slider.setRange(0, len(self._positions) - 1); self.slider.blockSignals(False)
            self.series.set_series(self._x, self._frames_series, self._x_label)
            self.go(0)

    def _show_modes(self) -> None:
        self.pause()
        m = self._modes
        if m is None:
            return
        self.viewer.set_atoms(m.atoms.copy())
        self.slider.blockSignals(True); self.slider.setRange(0, PHASES_PER_PERIOD - 1); self.slider.blockSignals(False)
        f = np.asarray(m.frequencies_cm1, dtype=float)
        self.series.set_series(f, [(L("モード", "modes"), np.ones_like(f))], L("振動数 [cm⁻¹]", "frequency [cm⁻¹]"), stems=True)
        self.go(0)

    # ---- playing ----
    def count(self) -> int:
        if self.playing_modes() and self._modes is not None:
            return PHASES_PER_PERIOD
        return 0 if self._positions is None else len(self._positions)

    def go(self, i: int) -> None:
        n = self.count()
        if n == 0:
            return
        i = int(i) % n
        self.frame = i
        self.slider.blockSignals(True); self.slider.setValue(i); self.slider.blockSignals(False)
        self._apply(i)

    def _apply(self, i: int) -> None:
        if self.playing_modes() and self._modes is not None:
            k = int(self.mode.currentData() or 0)
            phase = 2 * np.pi * i / PHASES_PER_PERIOD
            self.viewer.set_positions(self._modes.displaced(k, self.amplitude.value(), phase))
            self.series.set_marker(self._modes.frequencies_cm1[k])
            self.position.setText(L(f"モード {k + 1}: {self._modes.frequencies_cm1[k]:.1f} cm⁻¹、位相 {i + 1}/{PHASES_PER_PERIOD}",
                                    f"mode {k + 1}: {self._modes.frequencies_cm1[k]:.1f} cm⁻¹, phase {i + 1}/{PHASES_PER_PERIOD}"))
        elif self._positions is not None:
            self.viewer.set_positions(self._positions[i])
            self.series.set_marker(float(self._x[i]))
            t = "" if self._times_fs is None else L(f"、t = {self._times_fs[i]:.1f} fs", f", t = {self._times_fs[i]:.1f} fs")
            self.position.setText(L(f"フレーム {i + 1} / {len(self._positions)}{t}", f"frame {i + 1} / {len(self._positions)}{t}"))
        self.frameChanged.emit(i)

    def _on_slider(self, v: int) -> None:
        self.go(v)

    def _on_series_click(self, xv: float) -> None:
        if self.playing_modes() and self._modes is not None:
            k = int(np.argmin(np.abs(np.asarray(self._modes.frequencies_cm1) - xv)))
            self.mode.setCurrentIndex(k)
        elif self._x is not None and len(self._x):
            self.go(int(np.argmin(np.abs(self._x - xv))))

    def toggle(self) -> None:
        if self._timer.isActive():
            self.pause()
        else:
            self.play()

    def play(self) -> None:
        if self.count() == 0:
            return
        self._timer.start(); self.btn_play.setText(L("停止", "Pause"))

    def pause(self) -> None:
        self._timer.stop(); self.btn_play.setText(L("再生", "Play"))

    def is_playing(self) -> bool:
        return self._timer.isActive()

    def _tick(self) -> None:
        self.go(self.frame + 1)


__all__ = ["PlaybackPanel", "SeriesCanvas", "PHASES_PER_PERIOD"]
