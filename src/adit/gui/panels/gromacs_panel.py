
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QComboBox, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit, QPushButton, QSpinBox,
                               QVBoxLayout, QWidget)

from adit.gui.style import ROW_SPACING
from adit.gui.widgets import SciDoubleSpinBox, add_row, file_row, label, narrow
from adit.lang import L
from adit.spec import GromacsMethod
from adit.web.codefields import mdp_text, parse_mdp

COULOMB = ["Cut-off", "PME", "Reaction-Field"]
CONSTRAINTS = ["none", "h-bonds", "all-bonds", "h-angles", "all-angles"]
PCOUPL = {"": "(NPT のときに選びます)", "C-rescale": "C-rescale", "Parrinello-Rahman": "Parrinello-Rahman"}


class GromacsMethodPanel(QWidget):
    changed = Signal()
    use_as_structure = Signal(str)

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.topology = QLineEdit(); self.topology.setPlaceholderText("CHARMM-GUI、acpype、pdb2gmx などで作った .top")
        self.structure_file = QLineEdit(); self.structure_file.setPlaceholderText(".gro か .pdb (トポロジーと同じ原子の並び)")
        self.load_structure = QPushButton("構造の群にも読み込む"); self.load_structure.setObjectName("link")
        self.load_structure.setToolTip("構造の群の作り方を「ファイル」にして、このファイルを読み込みます (原子数の確認と 3D 表示に使います)")
        self.load_structure.clicked.connect(lambda: self.structure_file.text().strip() and self.use_as_structure.emit(self.structure_file.text().strip()))
        self.coulomb = QComboBox(); self.coulomb.addItems(COULOMB)
        self.rcoulomb = narrow(SciDoubleSpinBox(0.0, 100.0, 1.0, 0.1))
        self.rvdw = narrow(SciDoubleSpinBox(0.0, 100.0, 1.0, 0.1))
        self.constraints = QComboBox(); self.constraints.addItems(CONSTRAINTS)
        self.pcoupl = QComboBox()
        for k, v in PCOUPL.items():
            self.pcoupl.addItem(v, k)
        self.compressibility = narrow(SciDoubleSpinBox(0.0, 1.0, 0.0, 1e-5))
        self.define = QLineEdit(); self.define.setPlaceholderText("例 -DPOSRES (位置の拘束)")
        self.checkpoint = QLineEdit(); self.checkpoint.setPlaceholderText("前の段階 (NVT → NPT → 本計算) の .cpt。空欄なら初速を作ります")
        self.gen_seed = narrow(QSpinBox()); self.gen_seed.setRange(0, 2**31 - 1); self.gen_seed.setValue(12345)
        self.extra = QPlainTextEdit(); self.extra.setMaximumHeight(70)
        self.extra.setPlaceholderText("画面にない mdp の項目を 1 行に 1 つ、名前 = 値 の形で (例 nstcalcenergy = 100)")
        note = QLabel("構造の正本は上の構造のファイル (.gro / .pdb) です。構造の群にも同じファイルを読み込んでください (原子数の確認と表示に使います)。"
                      "全電荷とスピン多重度はトポロジーが決めるので使いません (0 と 1 のまま)。k 点はありません")
        note.setObjectName("hint"); note.setWordWrap(True)
        topology_help = QLabel(L(
            'トポロジーの形式は <a href="https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html">GROMACS 公式文書</a>を参照してください。pdb2gmx は GROMACS のコマンド、CHARMM-GUI はウェブサービスです。研究室の力場がある場合は、その .top と付属の .itp を管理者や先輩に確認してください。',
            'See the <a href="https://manual.gromacs.org/current/reference-manual/topologies/topology-file-formats.html">official GROMACS documentation</a> for topology formats. pdb2gmx is a GROMACS command; CHARMM-GUI is a web service. If your group maintains a force field, ask its administrator or an experienced group member for the .top and accompanying .itp files.'))
        topology_help.setObjectName("hint"); topology_help.setWordWrap(True); topology_help.setOpenExternalLinks(True)

        cut = QWidget(); cut.setObjectName("rowbox"); ch = QHBoxLayout(cut); ch.setContentsMargins(0, 0, 0, 0); ch.setSpacing(8)
        ch.addWidget(QLabel("rcoulomb")); ch.addWidget(self.rcoulomb); ch.addWidget(QLabel("rvdw")); ch.addWidget(self.rvdw); ch.addStretch(1)
        conf = QWidget(); conf.setObjectName("rowbox"); cv = QVBoxLayout(conf); cv.setContentsMargins(0, 0, 0, 0); cv.setSpacing(2)
        cv.addWidget(file_row(self.structure_file, L("構造のファイル", "Structure file"),
                              L("GROMACS (*.gro *.pdb);;すべて (*)", "GROMACS (*.gro *.pdb);;All files (*)")))
        cv.addWidget(self.load_structure, 0, Qt.AlignmentFlag.AlignLeft)

        form = QFormLayout(self); form.setVerticalSpacing(ROW_SPACING)
        add_row(form, "トポロジー (.top)", file_row(self.topology, L("トポロジー", "Topology"), L("GROMACS (*.top);;すべて (*)", "GROMACS (*.top);;All files (*)")))
        form.addRow(label(""), topology_help)
        add_row(form, "構造のファイル (.gro / .pdb)", conf)
        add_row(form, "静電相互作用 (coulombtype)", self.coulomb)
        add_row(form, "カットオフ [nm]", cut)
        add_row(form, "拘束 (constraints)", self.constraints)
        add_row(form, "圧力浴 (pcoupl)", self.pcoupl)
        add_row(form, "等温圧縮率 [1/bar]", self.compressibility)
        add_row(form, "define", self.define)
        add_row(form, "前の段階の .cpt", file_row(self.checkpoint, L("前の段階の .cpt", "Previous stage .cpt"),
                                                  L("GROMACS (*.cpt);;すべて (*)", "GROMACS (*.cpt);;All files (*)")))
        add_row(form, "初速の乱数の種 (gen-seed)", self.gen_seed)
        add_row(form, "追加の mdp", self.extra)
        form.addRow(label(""), note)

        for w in (self.coulomb, self.constraints, self.pcoupl):
            w.currentTextChanged.connect(self._emit)
        for w in (self.rcoulomb, self.rvdw, self.compressibility, self.gen_seed):
            w.valueChanged.connect(self._emit)
        for w in (self.topology, self.structure_file, self.define, self.checkpoint):
            w.textChanged.connect(self._emit)
        self.extra.textChanged.connect(self._emit)

    def method(self) -> GromacsMethod:
        return GromacsMethod(topology_file=self.topology.text().strip(), structure_file=self.structure_file.text().strip(),
                             coulombtype=self.coulomb.currentText(), rcoulomb_nm=self.rcoulomb.value(), rvdw_nm=self.rvdw.value(),
                             constraints=self.constraints.currentText(), pcoupl=self.pcoupl.currentData(),
                             compressibility_per_bar=self.compressibility.value(), define=self.define.text().strip(),
                             checkpoint_file=self.checkpoint.text().strip(), gen_seed=self.gen_seed.value(), extra_mdp=parse_mdp(self.extra.toPlainText()))

    def set_method(self, m: GromacsMethod) -> None:
        self.topology.setText(m.topology_file); self.structure_file.setText(m.structure_file); self.coulomb.setCurrentText(m.coulombtype)
        self.rcoulomb.setValue(m.rcoulomb_nm); self.rvdw.setValue(m.rvdw_nm); self.constraints.setCurrentText(m.constraints)
        self.pcoupl.setCurrentIndex(max(0, self.pcoupl.findData(m.pcoupl))); self.compressibility.setValue(m.compressibility_per_bar)
        self.define.setText(m.define); self.checkpoint.setText(m.checkpoint_file); self.gen_seed.setValue(m.gen_seed)
        self.extra.setPlainText(mdp_text(m.extra_mdp))
        self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()
