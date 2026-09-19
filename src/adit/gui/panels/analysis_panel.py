
from __future__ import annotations

import shutil
from pathlib import Path

from PySide6.QtCore import QSize, Qt, Signal
from PySide6.QtGui import QPainter, QPixmap
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog, QFileDialog, QFormLayout, QFrame, QGroupBox, QHBoxLayout,
                               QHeaderView, QLabel, QLineEdit, QPlainTextEdit, QPushButton, QScrollArea, QSizePolicy, QSpinBox,
                               QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget)

from adit.analysis import AnalysisOptions, AnalysisResult, run_analysis
from adit.analysis.report import figure_title
from adit.gui import analysis_fields as AF
from adit.gui.analysis_views import Collapsible, SectionView, hint_label, open_folder
from adit.gui.help import help_for
from adit.gui.isosurface_panel import IsosurfacePanel
from adit.gui.playback import PlaybackPanel
from adit.lang import L
from adit.gui.style import ROW_SPACING
from adit.gui.widgets import SciDoubleSpinBox, add_row, narrow


class FigureHeader(QWidget):

    def __init__(self, title: str, path: str, parent: QWidget | None = None):
        super().__init__(parent)
        from adit.gui.copy_save import copy_button, save_button

        self.label = QLabel(title); self.label.setToolTip(path)
        self.label.setStyleSheet("font-weight: 600;")
        self.label.setWordWrap(True)
        self.label.setMinimumWidth(160)
        self.label.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.copy = copy_button(self); self.save = save_button(self)
        row = QHBoxLayout(self); row.setContentsMargins(0, 0, 0, 0); row.setSpacing(6)
        row.addWidget(self.label, 1); row.addWidget(self.copy); row.addWidget(self.save)

    def text(self) -> str:
        return self.label.text()


class FigureView(QWidget):

    def __init__(self, path: str, parent: QWidget | None = None):
        super().__init__(parent)
        self.path = path
        self.pixmap = QPixmap(path)
        self.setToolTip(path)
        sp = QSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred); sp.setHeightForWidth(True)
        self.setSizePolicy(sp)

    def _scaled_size(self, w: int) -> QSize:
        pw, ph = max(1, self.pixmap.width()), max(1, self.pixmap.height())
        w = min(w, pw)
        return QSize(w, round(ph * w / pw))

    def hasHeightForWidth(self) -> bool:
        return True

    def heightForWidth(self, w: int) -> int:
        return self._scaled_size(w).height()

    def sizeHint(self) -> QSize:
        return self.pixmap.size()

    def minimumSizeHint(self) -> QSize:
        return QSize(50, 30)

    def paintEvent(self, _ev) -> None:
        if self.pixmap.isNull():
            return
        s = self._scaled_size(self.width())
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform)
        p.drawPixmap(0, 0, s.width(), s.height(), self.pixmap)
        p.end()

    def resizeEvent(self, ev) -> None:
        super().resizeEvent(ev)
        if self.height() != self.heightForWidth(self.width()):
            self.setFixedHeight(self.heightForWidth(self.width()))


class SummaryLabel(QLabel):

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setWordWrap(True)
        self.setTextFormat(Qt.TextFormat.PlainText)
        self.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.setMinimumWidth(240)
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Minimum)

    def setPlainText(self, text: str) -> None:
        self.setText(text)

    def toPlainText(self) -> str:
        return self.text()


def _clear(lay) -> None:
    while lay.count():
        item = lay.takeAt(0)
        if item.widget():
            w = item.widget(); w.hide(); w.setParent(None); w.deleteLater()


def _combo(items) -> QComboBox:
    c = QComboBox()
    c.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon); c.setMinimumContentsLength(16)
    for v, ja, _en in items:
        c.addItem(ja, v)
    return c


def _check(key: str, checked: bool = False) -> QCheckBox:
    cb = QCheckBox(AF.LABELS[key][0]); cb.setChecked(checked)
    h = help_for(AF.LABELS[key][0])
    if h:
        cb.setToolTip(h.text())
    return cb


class AnalysisPanel(QWidget):
    analyzed = Signal(object)  # AnalysisResult

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.run_dir = QLineEdit(); self.run_dir.setPlaceholderText("計算結果のディレクトリ (生成すると自動で入ります)")
        self.browse = QPushButton("参照…")
        self.cb_energy = QCheckBox("エネルギーの推移"); self.cb_energy.setChecked(True)
        self.cb_temp = QCheckBox("温度の推移 (MD)"); self.cb_temp.setChecked(True)
        self.cb_bonds = QCheckBox("結合長 (最終構造)"); self.cb_bonds.setChecked(True)
        self.cb_rdf = QCheckBox("動径分布関数 (RDF)"); self.rdf_rmax = narrow(SciDoubleSpinBox(1.0, 50.0, 8.0, 1.0))
        self.cb_msd = QCheckBox("平均二乗変位 (MSD) と拡散係数"); self.msd_species = QComboBox(); self.msd_species.setEditable(True)
        self.msd_species.addItem(L("(全原子)", "(all atoms)"), "all")
        self.cb_dos = QCheckBox("状態密度 (DOS)"); self.dos_sigma = narrow(SciDoubleSpinBox(0.001, 5.0, 0.1, 0.05))
        self.cb_vib = QCheckBox("振動数とスペクトル"); self.cb_vib.setChecked(True)
        self.skip = narrow(QSpinBox()); self.skip.setRange(0, 10000000)
        self.btn_run = QPushButton("解析を実行"); self.btn_run.setObjectName("primary")
        self.btn_compare = QPushButton(AF.LABELS["compare"][0])
        from adit.gui import report_fields as RF

        self.btn_report = QPushButton(RF.LABELS["report"][0] + "…")
        self.btn_export = QPushButton(AF.LABELS["export"][0])
        self.cb_unwrap = _check("export_unwrap")
        self.summary = SummaryLabel()
        self.summary.setStyleSheet("font-family: monospace;")

        self.stride = narrow(QSpinBox()); self.stride.setRange(1, 10_000_000); self.stride.setValue(1)
        self.msd_fit_from = narrow(QLineEdit()); self.msd_fit_to = narrow(QLineEdit())
        self.cb_zdens = _check("zdens"); self.zdens_bin = narrow(SciDoubleSpinBox(0.01, 10.0, 0.2, 0.05))
        self.cb_stats = _check("stats", True)
        self.memory_mb = narrow(QLineEdit()); self.memory_mb.setPlaceholderText(AF.PLACEHOLDERS["memory_mb"][0])
        self.cb_pdos = _check("pdos")
        self.symprec = QLineEdit(); self.symprec.setPlaceholderText(AF.PLACEHOLDERS["symprec"][0])
        self.th_model = _combo(AF.THERMO_MODELS)
        self.th_temps = QLineEdit(); self.th_temps.setPlaceholderText(AF.PLACEHOLDERS["th_temps"][0])
        self.th_pressure = narrow(QLineEdit()); self.th_sigma = narrow(QLineEdit()); self.th_geometry = _combo(AF.GEOMETRIES)
        self.th_spin = narrow(QLineEdit()); self.th_imag = _combo(AF.IMAGINARY); self.th_exclude = narrow(QLineEdit())
        self.th_qh = narrow(QLineEdit()); self.th_tau = narrow(QLineEdit())
        self.uv_shape = _combo(AF.UV_SHAPES); self.uv_fwhm = narrow(QLineEdit())
        self.msd_axes = _combo(AF.MSD_AXES); self.msd_blocks = narrow(QLineEdit())
        self.msd_blocks.setPlaceholderText(AF.PLACEHOLDERS["msd_blocks"][0])
        self.cb_msd_keep_drift = _check("msd_keep_drift")
        self.cb_msd_per_atom = _check("msd_per_atom")
        self.cb_vanhove = _check("vanhove")
        self.vanhove_taus = narrow(QLineEdit()); self.vanhove_taus.setPlaceholderText(AF.PLACEHOLDERS["vanhove_taus"][0])
        self.vanhove_displacement = _combo(AF.VANHOVE_DISPLACEMENTS)
        self.select = QLineEdit(); self.select.setPlaceholderText(AF.PLACEHOLDERS["select"][0])
        self.rdf_pairs = QLineEdit(); self.rdf_pairs.setPlaceholderText(AF.PLACEHOLDERS["rdf_pairs"][0])
        self.zdens_axis = _combo(AF.AXES)
        self.coordination = narrow(QLineEdit()); self.coordination.setPlaceholderText(AF.PLACEHOLDERS["coordination"][0])
        self.centrosymmetry = narrow(QLineEdit()); self.centrosymmetry.setPlaceholderText(AF.PLACEHOLDERS["centrosymmetry"][0])
        self.steinhardt = narrow(QLineEdit()); self.steinhardt.setPlaceholderText(AF.PLACEHOLDERS["steinhardt"][0])
        self.clusters = narrow(QLineEdit()); self.clusters.setPlaceholderText(AF.PLACEHOLDERS["clusters"][0])
        self.adf = narrow(QLineEdit()); self.adf.setPlaceholderText(AF.PLACEHOLDERS["adf"][0])
        self.adf_cutoff = narrow(QLineEdit())
        self.cb_sq = _check("sq")
        self.hbond = narrow(QLineEdit()); self.hbond.setPlaceholderText(AF.PLACEHOLDERS["hbond"][0])
        self.cb_hbond_lifetime = _check("hbond_lifetime")
        self.hbond_cdf = narrow(QLineEdit()); self.hbond_cdf.setPlaceholderText(AF.PLACEHOLDERS["hbond_cdf"][0])
        self.cb_rg = _check("rg")
        self.density_grid = narrow(QLineEdit()); self.density_grid.setPlaceholderText(AF.PLACEHOLDERS["density_grid"][0])
        self.cb_voronoi = _check("voronoi")
        self.voronoi_face = narrow(QLineEdit()); self.voronoi_face.setPlaceholderText(AF.PLACEHOLDERS["voronoi_face"][0])
        self.sasa = narrow(QLineEdit()); self.sasa.setPlaceholderText(AF.PLACEHOLDERS["sasa"][0])
        self.distances = QLineEdit(); self.distances.setPlaceholderText(AF.PLACEHOLDERS["distances"][0])
        self.angles = QLineEdit(); self.angles.setPlaceholderText(AF.PLACEHOLDERS["angles"][0])
        self.dihedrals = QLineEdit(); self.dihedrals.setPlaceholderText(AF.PLACEHOLDERS["dihedrals"][0])
        self.rmsd_reference = narrow(QLineEdit()); self.rmsd_reference.setPlaceholderText(AF.PLACEHOLDERS["rmsd_reference"][0])
        self.cb_rmsf = _check("rmsf")
        self.cb_vacf = _check("vacf")
        self.conductivity_charge = narrow(QLineEdit()); self.conductivity_charge.setPlaceholderText(AF.PLACEHOLDERS["conductivity_charge"][0])
        self.conductivity_temperature = narrow(QLineEdit()); self.conductivity_temperature.setPlaceholderText(AF.PLACEHOLDERS["conductivity_temperature"][0])
        self.displacement = narrow(QLineEdit()); self.displacement.setPlaceholderText(AF.PLACEHOLDERS["displacement"][0])
        self.strain = narrow(QLineEdit()); self.strain.setPlaceholderText(AF.PLACEHOLDERS["strain"][0])
        self.pca = narrow(QLineEdit()); self.pca.setPlaceholderText(AF.PLACEHOLDERS["pca"][0])
        self.cluster = narrow(QLineEdit()); self.cluster.setPlaceholderText(AF.PLACEHOLDERS["cluster"][0])
        self.fes = narrow(QLineEdit()); self.fes.setPlaceholderText(AF.PLACEHOLDERS["fes"][0])
        self.fes_bins = narrow(QLineEdit()); self.fes_bins.setPlaceholderText(AF.PLACEHOLDERS["fes_bins"][0])
        self.conformer_temperature = narrow(QLineEdit()); self.conformer_temperature.setPlaceholderText(AF.PLACEHOLDERS["conformer_temperature"][0])
        self.fes_unit = _combo(AF.FES_UNITS)
        self.bands_window = narrow(QLineEdit()); self.bands_window.setPlaceholderText(AF.PLACEHOLDERS["bands_window"][0])
        self.effective_mass_points = narrow(QLineEdit())
        self.effective_mass_points.setPlaceholderText(AF.PLACEHOLDERS["effective_mass_points"][0])
        self.bader = QLineEdit(); self.bader.setPlaceholderText(AF.PLACEHOLDERS["bader"][0])
        self.bader_valence = narrow(QLineEdit()); self.bader_valence.setPlaceholderText(AF.PLACEHOLDERS["bader_valence"][0])
        self.xrd = narrow(QLineEdit()); self.xrd.setPlaceholderText(AF.PLACEHOLDERS["xrd"][0])
        self.xrd_range = narrow(QLineEdit()); self.xrd_range.setPlaceholderText(AF.PLACEHOLDERS["xrd_range"][0])
        self.xrd_measured = QLineEdit(); self.xrd_measured.setPlaceholderText(AF.PLACEHOLDERS["xrd_measured"][0])
        self.cb_viscosity = _check("viscosity")
        self.plane_average = _combo(AF.PLANE_AXES)
        self.cube_unit = _combo(AF.CUBE_UNITS)
        self.cb_work_function = _check("work_function")
        self.heavy_limit = narrow(QLineEdit()); self.heavy_limit.setPlaceholderText(AF.PLACEHOLDERS["heavy_limit"][0])
        self.code = narrow(QLineEdit()); self.code.setPlaceholderText(AF.PLACEHOLDERS["code"][0])
        self.freq_scale = narrow(QLineEdit()); self.freq_scale.setPlaceholderText(AF.PLACEHOLDERS["freq_scale"][0])
        self.spectrum_measured = QLineEdit()
        self.spectrum_measured.setPlaceholderText(AF.PLACEHOLDERS["spectrum_measured"][0])
        self.plot_colors = QLineEdit(); self.plot_colors.setPlaceholderText(AF.PLACEHOLDERS["plot_colors"][0])
        self.plot_ticks = _combo(AF.TICK_DIRECTIONS)
        self.plot_grid = _combo(AF.GRID_CHOICES)
        self.plot_spines = _combo(AF.SPINE_CHOICES)
        self.plot_line_width = narrow(QLineEdit()); self.plot_line_width.setPlaceholderText(AF.PLACEHOLDERS["plot_line_width"][0])
        self.plot_font_size = narrow(QLineEdit()); self.plot_font_size.setPlaceholderText(AF.PLACEHOLDERS["plot_font_size"][0])
        self.plot_dpi = narrow(QLineEdit()); self.plot_dpi.setPlaceholderText(AF.PLACEHOLDERS["plot_dpi"][0])
        self.figure_format = narrow(QLineEdit())
        self.figure_format.setPlaceholderText(AF.PLACEHOLDERS["figure_format"][0])
        self.more = Collapsible(AF.LABELS["more"][0])
        mform = QFormLayout(); mform.setVerticalSpacing(ROW_SPACING)
        mform.addRow(hint_label(L("ここの欄はすべて任意です。空のままでも解析は走ります "
                                  "(エネルギー・温度・結合長・振動数・バンドは常に出ます)。必要な解析の欄だけ埋めてください。",
                                  "Every field here is optional: the analysis runs with all of them empty "
                                  "(energy, temperature, bond lengths, frequencies and bands are always produced). "
                                  "Fill in only what you need.")))
        mform.addRow(self._subhead(L("軌跡", "Trajectory")))
        add_row(mform, AF.LABELS["stride"][0], self.stride)
        row = QHBoxLayout(); row.addWidget(self.msd_fit_from); row.addWidget(QLabel("〜")); row.addWidget(self.msd_fit_to); row.addStretch()
        add_row(mform, AF.LABELS["msd_fit"][0], row)
        mform.addRow(hint_label(AF.PLACEHOLDERS["msd_fit"][0]))
        add_row(mform, AF.LABELS["msd_axes"][0], self.msd_axes)
        add_row(mform, AF.LABELS["msd_blocks"][0], self.msd_blocks)
        mform.addRow(self.cb_msd_keep_drift)
        mform.addRow(self.cb_msd_per_atom)
        mform.addRow(self.cb_vanhove)
        add_row(mform, AF.LABELS["vanhove_taus"][0], self.vanhove_taus)
        add_row(mform, AF.LABELS["vanhove_displacement"][0], self.vanhove_displacement)
        mform.addRow(hint_label(L("変位の分布は重い解析です。見積もりが 60 秒を超えると、その場では計算せず、実行用のファイル "
                                  "(msd_worker.py と msd_run.sh) を計算のディレクトリに置きます。それを実行したあと、もう一度「解析を実行」すると図になります",
                                  "The displacement distribution is heavy: if the estimate exceeds 60 s, files to run it (msd_worker.py, msd_run.sh) "
                                  "are written to the run directory instead; run them and press Run analysis again to get the figures")))
        add_row(mform, AF.LABELS["memory_mb"][0], self.memory_mb)
        mform.addRow(self.cb_zdens)
        add_row(mform, AF.LABELS["zdens_bin"][0], self.zdens_bin)
        mform.addRow(self.cb_stats)
        mform.addRow(self._subhead(L("熱化学 (振動数から ASE で計算。欄に既定値は入れていません。足りない欄があれば、計算しない理由が要約に出ます)",
                                     "Thermochemistry (computed with ASE from the frequencies; no default values are filled in; if something is missing, "
                                     "the summary says why nothing was computed)")))
        for key, w in (("th_model", self.th_model), ("th_temps", self.th_temps), ("th_pressure", self.th_pressure), ("th_sigma", self.th_sigma),
                       ("th_geometry", self.th_geometry), ("th_spin", self.th_spin), ("th_imag", self.th_imag), ("th_exclude", self.th_exclude),
                       ("th_qh", self.th_qh), ("th_tau", self.th_tau)):
            add_row(mform, AF.LABELS[key][0], w)
        mform.addRow(self._subhead(L("スペクトルと構造", "Spectra and structure")))
        add_row(mform, AF.LABELS["freq_scale"][0], self.freq_scale)
        add_row(mform, AF.LABELS["spectrum_measured"][0], self.spectrum_measured)
        add_row(mform, AF.LABELS["code"][0], self.code)
        add_row(mform, AF.LABELS["uv_shape"][0], self.uv_shape)
        add_row(mform, AF.LABELS["uv_fwhm"][0], self.uv_fwhm)
        mform.addRow(self.cb_pdos)
        add_row(mform, AF.LABELS["symprec"][0], self.symprec)
        mform.addRow(self._subhead(L("原子の選び方と近傍の構造", "Atom selection and local structure")))
        add_row(mform, AF.LABELS["select"][0], self.select)
        add_row(mform, AF.LABELS["rdf_pairs"][0], self.rdf_pairs)
        add_row(mform, AF.LABELS["zdens_axis"][0], self.zdens_axis)
        for key, widget in (("coordination", self.coordination), ("centrosymmetry", self.centrosymmetry),
                            ("steinhardt", self.steinhardt), ("clusters", self.clusters),
                            ("adf", self.adf), ("adf_cutoff", self.adf_cutoff), ("hbond", self.hbond), ("hbond_cdf", self.hbond_cdf),
                            ("density_grid", self.density_grid), ("voronoi_face", self.voronoi_face),
                            ("sasa", self.sasa)):
            add_row(mform, AF.LABELS[key][0], widget)
        mform.addRow(self.cb_sq); mform.addRow(self.cb_hbond_lifetime); mform.addRow(self.cb_rg); mform.addRow(self.cb_voronoi)
        mform.addRow(self._subhead(L("時系列と分布", "Time series and distributions")))
        for key, widget in (("distances", self.distances), ("angles", self.angles), ("dihedrals", self.dihedrals),
                            ("rmsd_reference", self.rmsd_reference), ("conductivity_charge", self.conductivity_charge),
                            ("conductivity_temperature", self.conductivity_temperature),
                            ("displacement", self.displacement), ("strain", self.strain),
                            ("pca", self.pca), ("cluster", self.cluster),
                            ("fes", self.fes), ("fes_bins", self.fes_bins), ("fes_unit", self.fes_unit),
                            ("conformer_temperature", self.conformer_temperature)):
            add_row(mform, AF.LABELS[key][0], widget)
        mform.addRow(self.cb_rmsf); mform.addRow(self.cb_vacf); mform.addRow(self.cb_viscosity)
        mform.addRow(self._subhead(L("電子・回折・体積データ", "Electronic structure, diffraction and volumetric data")))
        for key, widget in (("bands_window", self.bands_window), ("effective_mass_points", self.effective_mass_points),
                            ("bader", self.bader), ("bader_valence", self.bader_valence),
                            ("xrd", self.xrd), ("xrd_range", self.xrd_range), ("xrd_measured", self.xrd_measured),
                            ("plane_average", self.plane_average), ("cube_unit", self.cube_unit),
                            ("heavy_limit", self.heavy_limit)):
            add_row(mform, AF.LABELS[key][0], widget)
        mform.addRow(self._subhead(L("図の見た目", "Figure style")))
        for key, widget in (("plot_colors", self.plot_colors), ("plot_ticks", self.plot_ticks),
                            ("plot_grid", self.plot_grid), ("plot_spines", self.plot_spines),
                            ("plot_line_width", self.plot_line_width), ("plot_font_size", self.plot_font_size),
                            ("plot_dpi", self.plot_dpi), ("figure_format", self.figure_format)):
            add_row(mform, AF.LABELS[key][0], widget)
        mform.addRow(hint_label(L("数値は変わりません。論文誌や研究室の作法に合わせるための指定です。",
                                  "These change only the appearance; the numbers stay the same.")))
        mform.addRow(self.cb_work_function)
        mform.addRow(hint_label(L("射影バンド (projwfc.x の filproj)・光学 (epsilon.x)・VASP の PROCAR は、"
                                  "計算結果のディレクトリにあれば自動で読みます。",
                                  "Projected bands (filproj from projwfc.x), optics (epsilon.x) and a VASP PROCAR "
                                  "are read automatically when they are in the run directory.")))
        self.more.body_layout.addLayout(mform)

        box = QGroupBox("解析"); box.setObjectName("analysis")
        form = QFormLayout(box); form.setVerticalSpacing(ROW_SPACING)
        row = QHBoxLayout(); row.addWidget(self.run_dir); row.addWidget(self.browse)
        add_row(form, "計算結果のディレクトリ", row)
        form.addRow(self.cb_energy); form.addRow(self.cb_temp); form.addRow(self.cb_bonds)
        form.addRow(self.cb_rdf); add_row(form, "r の最大値 [Å]", self.rdf_rmax)
        form.addRow(self.cb_msd); add_row(form, "元素", self.msd_species)
        form.addRow(self.cb_dos); add_row(form, "ガウス幅 [eV]", self.dos_sigma)
        form.addRow(self.cb_vib)
        add_row(form, "平衡化として捨てるフレーム数", self.skip)
        form.addRow(self.more)
        self.scan_note = QLabel(L("1 つの条件だけを変えた一連の計算です。「解析を実行」で、値ごとのエネルギー・力・圧力を表と図にします (上のチェック項目は使いません)。",
                                  "This is a parameter scan. \"Run analysis\" tabulates and plots the energy, forces and pressure for each value (the options above are not used)."))
        self.scan_note.setObjectName("hint"); self.scan_note.setWordWrap(True); self.scan_note.hide()
        form.addRow(self.scan_note)
        self._options = [self.cb_energy, self.cb_temp, self.cb_bonds, self.cb_rdf, self.rdf_rmax, self.cb_msd, self.msd_species,
                         self.cb_dos, self.dos_sigma, self.cb_vib, self.skip, self.more, self.btn_export, self.cb_unwrap,
                         self.msd_axes, self.msd_blocks, self.cb_msd_keep_drift, self.cb_msd_per_atom,
                         self.cb_vanhove, self.vanhove_taus, self.vanhove_displacement]
        row = QHBoxLayout(); row.addWidget(self.btn_run); row.addWidget(self.btn_compare)
        row.addWidget(self.btn_report); row.addStretch(); form.addRow(row)
        row = QHBoxLayout(); row.addWidget(self.btn_export); row.addStretch(); form.addRow(row)
        form.addRow(self.cb_unwrap)
        for f in (form, mform):
            f.setRowWrapPolicy(QFormLayout.RowWrapPolicy.WrapLongRows)

        self.too_large = QFrame(); self.too_large.setObjectName("too_large")
        tl = QVBoxLayout(self.too_large); tl.setContentsMargins(0, 0, 0, 0)
        self.too_large_text = QLabel(); self.too_large_text.setObjectName("error"); self.too_large_text.setWordWrap(True)
        self.too_large_text.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.btn_use_stride = QPushButton()
        tl.addWidget(self.too_large_text); tl.addWidget(self.btn_use_stride, 0, Qt.AlignmentFlag.AlignLeft)
        self.too_large.hide()
        self._suggested_stride: int | None = None

        self.scan_table = QTableWidget(); self.scan_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.scan_table.verticalHeader().setVisible(False); self.scan_table.hide()

        self.export_box = Collapsible(L("書き出したファイル (export_README.txt)", "Exported files (export_README.txt)"), expanded=True)
        self.export_path = hint_label("")
        self.btn_open_export = QPushButton(AF.LABELS["open_export"][0])
        self.export_readme = QPlainTextEdit(); self.export_readme.setReadOnly(True); self.export_readme.setFixedHeight(220)
        self.export_readme.setStyleSheet("font-family: monospace;")
        self.export_box.body_layout.addWidget(self.export_path)
        self.export_box.body_layout.addWidget(self.btn_open_export, 0, Qt.AlignmentFlag.AlignLeft)
        self.export_box.body_layout.addWidget(self.export_readme)
        self.export_box.hide()
        self._export_dir = ""
        self.last_open_ok: bool | None = None

        self.play_box = Collapsible(L("3D で再生 (軌跡・最適化のステップ・振動モード)", "Play in 3D (trajectory, optimization steps, vibrational modes)"),
                                    expanded=True)
        self.playback = PlaybackPanel()
        self.play_box.body_layout.addWidget(self.playback)
        self.play_box.hide()
        self.iso_box = Collapsible(L("等値面 (軌道・電子密度・静電ポテンシャル)", "Isosurface (orbitals, electron density, electrostatic potential)"),
                                   expanded=True)
        self.iso = IsosurfacePanel()
        self.iso_box.body_layout.addWidget(self.iso)
        self.iso_box.hide()

        self.sections = QWidget(); self.sections_lay = QVBoxLayout(self.sections)
        self.sections_lay.setContentsMargins(0, 0, 0, 0); self.sections_lay.setSpacing(8)
        self.figs = QWidget(); self.figs_lay = QVBoxLayout(self.figs); self.figs_lay.setAlignment(Qt.AlignmentFlag.AlignTop)
        self.figs_lay.setContentsMargins(0, 0, 0, 0)

        content = QWidget()
        lay = QVBoxLayout(content); lay.setContentsMargins(12, 12, 12, 12); lay.setSpacing(10)
        for w in (box, self.too_large, self.scan_table, self.summary, self.play_box, self.iso_box, self.export_box, self.sections, self.figs):
            lay.addWidget(w)
        lay.addStretch(1)
        self.scroll = QScrollArea(); self.scroll.setWidget(content); self.scroll.setWidgetResizable(True)
        self.scroll.setFrameShape(QFrame.Shape.NoFrame)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        outer = QVBoxLayout(self); outer.setContentsMargins(0, 0, 0, 0); outer.addWidget(self.scroll)

        self._compare_rows: list[tuple[str, str, str]] = []
        self.browse.clicked.connect(self._browse)
        self.btn_run.clicked.connect(self.run)
        self.btn_export.clicked.connect(self.export)
        self.btn_compare.clicked.connect(self._open_compare)
        self.btn_report.clicked.connect(self._open_report)
        self.btn_use_stride.clicked.connect(self.use_suggested_stride)
        self.btn_open_export.clicked.connect(self.open_export_folder)
        self.run_dir.editingFinished.connect(self._on_dir_edited)

    @staticmethod
    def _subhead(text: str) -> QLabel:
        w = QLabel(text); w.setWordWrap(True); w.setStyleSheet("font-weight: 600; margin-top: 6px;")
        return w

    def set_run_dir(self, path: Path | str, elements: list[str] | None = None) -> None:
        self.run_dir.setText(str(path))
        all_atoms = self.msd_species.currentIndex() == 0 and self.msd_species.currentData() == "all"
        cur = self.msd_species.currentText()
        self.msd_species.clear(); self.msd_species.addItem(L("(全原子)", "(all atoms)"), "all"); self.msd_species.addItems(elements or [])
        if cur and not all_atoms:
            self.msd_species.setCurrentText(cur)
        self._apply_defaults(str(path))
        self._refresh_iso(str(path))

    def _on_dir_edited(self) -> None:
        self._apply_defaults(self.run_dir.text().strip())
        self._refresh_iso(self.run_dir.text().strip())

    def _refresh_iso(self, d: str) -> None:
        try:
            n = self.iso.set_run_dir(d if d and Path(d).is_dir() else None)
        except Exception:
            n = 0
        self.iso_box.setVisible(n > 0)

    def _apply_defaults(self, d: str) -> None:
        scan = self.is_scan_dir(d)
        self.scan_note.setVisible(scan)
        for w in self._options:
            w.setEnabled(not scan)
        if not d or not Path(d).is_dir() or scan:
            return
        p = Path(d)
        try:
            from adit.project import load_project
            md = load_project(p).task.type == "molecular_dynamics"
        except Exception:
            md = any((p / n).is_file() for n in ("md.out", "geo_end.xyz", "xtb.trj", "XDATCAR"))
        self.cb_rdf.setChecked(md); self.cb_msd.setChecked(md)

    @staticmethod
    def is_scan_dir(d: str | Path) -> bool:
        from adit.scan import SCAN_FILE
        return bool(str(d).strip()) and (Path(d).expanduser() / SCAN_FILE).is_file()

    def _show_scan_table(self, rows: list[dict], reference: str = "last") -> None:
        from adit.scan import diff_label
        have = [r for r in rows if r["energy_ev"]]
        per_atom = bool(have) and all(r["diff_per_atom_mev"] for r in have)
        diff_key = "diff_per_atom_mev" if per_atom else "diff_from_last_mev"
        ref_ja, ref_en = diff_label(reference)
        cols = [("value", L("値", "Value")), ("energy_ev", L("エネルギー\n[eV]", "Energy\n[eV]")),
                (diff_key, L(f"{ref_ja}\n[meV/原子]", f"{ref_en}\n[meV/atom]") if per_atom
                 else L(f"{ref_ja}\n[meV]", f"{ref_en}\n[meV]")),
                ("fmax_ev_ang", L("力の最大値\n[eV/Å]", "Max force\n[eV/Å]")), ("pressure_gpa", L("圧力\n[GPa]", "Pressure\n[GPa]")),
                ("note", L("注", "Note"))]
        t = self.scan_table
        t.clear(); t.setColumnCount(len(cols)); t.setRowCount(len(rows))
        t.setHorizontalHeaderLabels([c[1] for c in cols])
        right = Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        for i, r in enumerate(rows):
            for j, (key, _) in enumerate(cols):
                item = QTableWidgetItem(r.get(key) or ("" if key == "note" else "-"))
                if key not in ("note", "value"):
                    item.setTextAlignment(right)
                t.setItem(i, j, item)
        t.resizeColumnsToContents()
        t.horizontalHeader().setSectionResizeMode(len(cols) - 1, QHeaderView.ResizeMode.Stretch)
        h = (t.horizontalHeader().sizeHint().height() + sum(t.rowHeight(i) for i in range(t.rowCount())) + 2 * t.frameWidth()
             + t.horizontalScrollBar().sizeHint().height() + 2)
        t.setFixedHeight(min(h, 260))
        t.show()

    def _show_figures(self, figures: dict) -> None:
        from adit.gui.copy_save import copy_button, save_button

        _clear(self.figs_lay)
        for name, path in figures.items():
            wrap = FigureHeader(figure_title(name), path)
            wrap.copy.clicked.connect(lambda _c=False, p=path: self._copy_figure(p))
            wrap.save.clicked.connect(lambda _c=False, p=path, n=name: self._save_figure(p, n))
            self.figs_lay.addWidget(wrap); self.figs_lay.addWidget(FigureView(path))

    @staticmethod
    def figure_siblings(path: str) -> dict[str, Path]:
        base = Path(path)
        return {suffix: base.with_suffix("." + suffix) for suffix in ("svg", "pdf", "eps")
                if base.with_suffix("." + suffix).is_file()}

    def _copy_figure(self, path: str) -> None:
        from adit.gui.copy_save import copy_file

        copy_file(Path(path))
        self.summary.setText(L(f"図をクリップボードにコピーしました ({Path(path).name})。Word や PowerPoint に貼れます",
                               f"the figure is on the clipboard ({Path(path).name}); paste it into Word or PowerPoint"))

    def _save_figure(self, path: str, name: str) -> None:
        from PySide6.QtWidgets import QFileDialog, QMessageBox

        others = self.figure_siblings(path)
        filters = [L("PNG 画像 (*.png)", "PNG image (*.png)")]
        for suffix in others:
            filters.insert(0, L(f"{suffix.upper()} (拡大しても粗くならない) (*.{suffix})",
                                f"{suffix.upper()}, stays sharp when enlarged (*.{suffix})"))
        first = next(iter(others), "png")
        target, _ = QFileDialog.getSaveFileName(self, L("図を保存…", "Save figure…"),
                                                f"{name}.{first}", ";;".join(filters))
        if not target:
            return
        suffix = Path(target).suffix.lstrip(".").lower()
        source = others.get(suffix, Path(path))
        try:
            shutil.copyfile(source, target)
        except OSError as ex:
            QMessageBox.warning(self, L("保存できません", "Cannot save"), str(ex))
            return
        note = "" if suffix in others else L(
            "  (ベクタ形式で保存するには、詳しい条件の「図の追加の形式」に svg と入れて、もう一度解析してください)",
            "  (to save a vector version, set \"Extra figure formats\" to svg and run the analysis again)")
        self.summary.setText(L(f"図を保存しました: {target}", f"saved the figure: {target}") + note)

    def _show_sections(self, sections: list) -> None:
        _clear(self.sections_lay)
        for sec in sections:
            self.sections_lay.addWidget(SectionView(sec))

    def section_views(self) -> list[SectionView]:
        return [self.sections_lay.itemAt(i).widget() for i in range(self.sections_lay.count())]

    def _msd_species(self) -> str | None:
        i = self.msd_species.currentIndex()
        text = self.msd_species.currentText().strip()
        if (i == 0 and self.msd_species.itemData(0) == "all" and text == self.msd_species.itemText(0)) or not text:
            return None
        return text

    def fields(self) -> dict[str, str]:
        def on(cb: QCheckBox) -> str:
            return "on" if cb.isChecked() else ""

        def data(c: QComboBox) -> str:
            return str(c.currentData() or "")

        return {"rdf": on(self.cb_rdf), "rmax": repr(self.rdf_rmax.value()), "msd": on(self.cb_msd), "msd_species": self._msd_species() or "",
                "dos": on(self.cb_dos), "sigma": repr(self.dos_sigma.value()), "skip": str(self.skip.value()), "stride": str(self.stride.value()),
                "msd_fit_from": self.msd_fit_from.text(), "msd_fit_to": self.msd_fit_to.text(), "zdens": on(self.cb_zdens),
                "msd_axes": data(self.msd_axes), "msd_blocks": self.msd_blocks.text(),
                "msd_keep_drift": on(self.cb_msd_keep_drift), "msd_per_atom": on(self.cb_msd_per_atom),
                "vanhove": on(self.cb_vanhove), "vanhove_taus": self.vanhove_taus.text(),
                "vanhove_displacement": data(self.vanhove_displacement),
                "zdens_bin": repr(self.zdens_bin.value()), "stats": on(self.cb_stats), "memory_mb": self.memory_mb.text(), "pdos": on(self.cb_pdos),
                "symprec": self.symprec.text(), "th_model": data(self.th_model), "th_temps": self.th_temps.text(), "th_pressure": self.th_pressure.text(),
                "th_sigma": self.th_sigma.text(), "th_geometry": data(self.th_geometry), "th_spin": self.th_spin.text(), "th_imag": data(self.th_imag),
                "th_exclude": self.th_exclude.text(), "th_qh": self.th_qh.text(), "th_tau": self.th_tau.text(), "uv_shape": data(self.uv_shape),
                "uv_fwhm": self.uv_fwhm.text(), "export_unwrap": on(self.cb_unwrap),
                "select": self.select.text(), "rdf_pairs": self.rdf_pairs.text(), "zdens_axis": data(self.zdens_axis),
                "coordination": self.coordination.text(), "centrosymmetry": self.centrosymmetry.text(),
                "steinhardt": self.steinhardt.text(), "clusters": self.clusters.text(),
                "adf": self.adf.text(), "adf_cutoff": self.adf_cutoff.text(), "sq": on(self.cb_sq),
                "hbond": self.hbond.text(), "hbond_lifetime": on(self.cb_hbond_lifetime), "hbond_cdf": self.hbond_cdf.text(),
                "rg": on(self.cb_rg), "density_grid": self.density_grid.text(),
                "voronoi": on(self.cb_voronoi), "voronoi_face": self.voronoi_face.text(), "sasa": self.sasa.text(),
                "distances": self.distances.text(), "angles": self.angles.text(), "dihedrals": self.dihedrals.text(),
                "rmsd_reference": self.rmsd_reference.text(), "rmsf": on(self.cb_rmsf), "vacf": on(self.cb_vacf),
                "conductivity_charge": self.conductivity_charge.text(),
                "conductivity_temperature": self.conductivity_temperature.text(),
                "displacement": self.displacement.text(), "strain": self.strain.text(),
                "pca": self.pca.text(), "cluster": self.cluster.text(),
                "fes": self.fes.text(), "fes_bins": self.fes_bins.text(), "fes_unit": data(self.fes_unit),
                "conformer_temperature": self.conformer_temperature.text(),
                "bands_window": self.bands_window.text(), "effective_mass_points": self.effective_mass_points.text(),
                "bader": self.bader.text(), "bader_valence": self.bader_valence.text(),
                "xrd": self.xrd.text(), "xrd_range": self.xrd_range.text(), "xrd_measured": self.xrd_measured.text(),
                "viscosity": on(self.cb_viscosity), "plane_average": data(self.plane_average),
                "cube_unit": data(self.cube_unit), "work_function": on(self.cb_work_function),
                "heavy_limit": self.heavy_limit.text(), "code": self.code.text(),
                "freq_scale": self.freq_scale.text(),
                "spectrum_measured": self.spectrum_measured.text(),
                "plot_colors": self.plot_colors.text(), "plot_ticks": data(self.plot_ticks),
                "plot_grid": data(self.plot_grid), "plot_spines": data(self.plot_spines),
                "plot_line_width": self.plot_line_width.text(), "plot_font_size": self.plot_font_size.text(),
                "plot_dpi": self.plot_dpi.text(), "figure_format": self.figure_format.text()}

    def options(self) -> AnalysisOptions:
        o = AF.options_from_fields(self.fields())
        o.energy, o.temperature, o.bonds, o.vibrations = (self.cb_energy.isChecked(), self.cb_temp.isChecked(), self.cb_bonds.isChecked(),
                                                          self.cb_vib.isChecked())
        return o

    def run(self, export: bool = False):
        from adit.analysis.trajectory import TrajectoryTooLarge
        d = self.run_dir.text().strip()
        if not d or not Path(d).is_dir():
            self.summary.setPlainText(L(f"ディレクトリがありません: {d!r}", f"directory not found: {d!r}"))
            return None
        scan = self.is_scan_dir(d)
        self.scan_table.hide(); self.too_large.hide(); self.play_box.hide(); self.playback.clear()
        try:
            if scan:
                from adit.scan import analyze_scan
                res = analyze_scan(d)
            else:
                opts = self.options()
                opts.export = bool(export)
                res = run_analysis(d, opts)
        except AF.FieldError as ex:
            self.summary.setPlainText(L(f"解析の条件を読めません: {ex}", f"cannot read the analysis options: {ex}"))
            return None
        except TrajectoryTooLarge as ex:
            self._show_too_large(ex)
            self.summary.setPlainText(L(f"解析できません: {ex}", f"cannot analyze: {ex}"))
            return None
        except Exception as ex:
            self.summary.setPlainText(L(f"解析できません: {ex}", f"cannot analyze: {ex}"))
            return None
        if scan:
            self._show_scan_table(res.rows, res.reference)
            self.summary.setPlainText(res.summary_text(include_table=False))
            self._show_sections([]); self.export_box.hide()
        else:
            self.summary.setPlainText(res.summary_text())
            self._show_sections(AF.result_sections(res))
            self._show_export(res)
            self._show_playback(d, opts)
        self._show_figures(res.figures or {})
        self.analyzed.emit(res)
        return res

    def export(self):
        return self.run(export=True)

    def _show_playback(self, run_dir: str, opts) -> None:
        try:
            ok = self.playback.load_run(run_dir, skip=int(opts.skip_frames), stride=int(opts.stride), memory_mb=opts.memory_budget_mb)
        except Exception as ex:
            self.playback.note.setText(L(f"再生できません: {ex}", f"cannot play: {ex}")); ok = False
        self.play_box.setVisible(ok or bool(self.playback.note.text()))

    def _show_export(self, res) -> None:
        d, text = AF.read_export_readme(res)
        self._export_dir = d
        if not d:
            self.export_box.hide(); return
        self.export_path.setText(L(f"書き出し先: {d}", f"written to: {d}"))
        self.export_readme.setPlainText(text)
        self.export_box.show()

    def open_export_folder(self) -> bool:
        ok = bool(self._export_dir) and open_folder(self._export_dir)
        self.last_open_ok = ok
        if not ok:
            self.export_path.setText(L(f"フォルダを開けませんでした。場所: {self._export_dir}", f"could not open the folder; it is at: {self._export_dir}"))
        return ok

    def _show_too_large(self, ex) -> None:
        from adit.gui.style import current_theme, tokens_for
        n = getattr(ex, "suggested_stride", None)
        self.too_large_text.setText(str(ex))
        self.too_large_text.setStyleSheet(f"color: {tokens_for(current_theme()).ng}; font-weight: 600;")
        self._suggested_stride = int(n) if n else None
        self.btn_use_stride.setVisible(bool(n))
        if n:
            self.btn_use_stride.setText(L(f"間引きを {n} にする", f"Set the stride to {n}"))
        self.too_large.show()

    def use_suggested_stride(self) -> None:
        n = self._suggested_stride
        if not n:
            return
        self.stride.setValue(n)
        self.more.set_expanded(True)
        self.too_large.hide()
        self.summary.setPlainText(L(f"間引きの欄に {n} を入れました。「解析を実行」をもう一度押してください。",
                                    f"The stride is now {n}. Press \"Run analysis\" again."))

    def _open_report(self) -> None:
        from adit.gui.report_dialog import ReportDialog

        d = self.run_dir.text().strip()
        ReportDialog([d] if d else [], self).exec()

    def _open_compare(self) -> None:
        from adit.gui.compare_dialog import CompareDialog
        d = self.run_dir.text().strip()
        base = str(Path(d).expanduser().parent) if d and Path(d).expanduser().is_dir() else d
        dlg = CompareDialog(base, self)
        if self._compare_rows:
            dlg.set_rows(self._compare_rows)
        else:
            dlg.load_compare_json()
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._compare_rows = dlg.rows()
            self.run_compare(dlg.base_dir(), dlg.reactions())

    def run_compare(self, base: str | Path, reactions):
        from adit.analysis.compare import CompareError, analyze_compare
        self.too_large.hide(); self.scan_table.hide(); self.export_box.hide(); self.play_box.hide(); self.playback.clear()
        try:
            cres = analyze_compare(base, reactions)
        except (CompareError, ValueError, OSError) as ex:
            self.summary.setPlainText(L(f"比べられません: {ex}", f"cannot compare: {ex}"))
            return None
        self.summary.setPlainText(cres.summary_text())
        self._show_sections(AF.compare_sections(cres))
        self._show_figures(cres.figures or {})
        return cres

    def _browse(self) -> None:
        d = QFileDialog.getExistingDirectory(self, L("計算結果のディレクトリ", "Run directory"), self.run_dir.text() or str(Path.home()))
        if d:
            self.run_dir.setText(d)
