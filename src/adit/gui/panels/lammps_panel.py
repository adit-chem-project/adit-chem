
from __future__ import annotations

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QFormLayout, QLabel, QLineEdit, QPlainTextEdit, QSpinBox, QWidget

from adit.gui.style import ROW_SPACING
from adit.gui.widgets import add_row, file_row, label, narrow
from adit.lang import L
from adit.spec import LammpsMethod
from adit.web.codefields import parse_lines, parse_words

UNITS = {"": "(選んでください)", "metal": "metal (eV, Å, ps, bar)", "real": "real (kcal/mol, Å, fs, atm)"}
ATOM_STYLES = ["atomic", "charge", "full", "molecular", "bond", "angle"]


class LammpsMethodPanel(QWidget):
    changed = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.units = QComboBox()
        for k, v in UNITS.items():
            self.units.addItem(v, k)
        self.atom_style = QComboBox(); self.atom_style.setEditable(True); self.atom_style.addItems(ATOM_STYLES)
        self.data_file = QLineEdit(); self.data_file.setPlaceholderText("空欄なら構造から data.lammps を書きます (atomic か charge のとき)")
        self.type_elements = QLineEdit(); self.type_elements.setPlaceholderText("型番号 1, 2, … の元素を順に空白区切りで (例 O H)。空欄なら構造の元素の順")
        self.pair_style = QLineEdit(); self.pair_style.setPlaceholderText("例 eam、eam/alloy、reaxff NULL、mace no_domain_decomposition")
        self.pair_coeff = QPlainTextEdit(); self.pair_coeff.setMaximumHeight(70)
        self.pair_coeff.setPlaceholderText("pair_coeff の行 (行頭の pair_coeff は省けます。例 * * Cu_u3.eam)")
        self.potential_files = QPlainTextEdit(); self.potential_files.setMaximumHeight(70)
        self.potential_files.setPlaceholderText("生成先へ写すファイル (力場・モデル) を 1 行に 1 つ。入力の中ではファイル名だけで書きます")
        self.style_commands = QPlainTextEdit(); self.style_commands.setMaximumHeight(60)
        self.style_commands.setPlaceholderText("例 bond_style harmonic / special_bonds lj/coul 0 0 0.5")
        self.extra_commands = QPlainTextEdit(); self.extra_commands.setMaximumHeight(60)
        self.extra_commands.setPlaceholderText("例 kspace_style pppm 1e-4 / neigh_modify every 1")
        self.seed = narrow(QSpinBox()); self.seed.setRange(1, 2**31 - 1); self.seed.setValue(12345)
        self._thermo_pressure_tensor = False    # no widget; kept for the round trip
        note = QLabel("相互作用は外部のもの (力場のファイル、data ファイル、機械学習ポテンシャルのモデル) をそのまま使います。"
                      "全電荷とスピン多重度は使いません (0 と 1 のまま)。k 点はありません")
        note.setObjectName("hint"); note.setWordWrap(True)
        force_field_help = QLabel(L(
            'pair_style と力場ファイルは同じ出典の組を使います。利用できる形式は <a href="https://docs.lammps.org/Howto_force_fields.html">LAMMPS 公式の力場案内</a>と <a href="https://docs.lammps.org/pair_style.html">pair_style 一覧</a>で確認できます。研究室で使う力場が決まっている場合は、対応する data ファイル・pair_style・pair_coeff を管理者や先輩に確認してください。',
            'Use pair_style and force-field files from the same source. See the official LAMMPS <a href="https://docs.lammps.org/Howto_force_fields.html">force-field guide</a> and <a href="https://docs.lammps.org/pair_style.html">pair_style list</a>. If your group has a designated force field, ask its administrator or an experienced group member for the matching data file, pair_style and pair_coeff.'))
        force_field_help.setObjectName("hint"); force_field_help.setWordWrap(True); force_field_help.setOpenExternalLinks(True)

        form = QFormLayout(self); form.setVerticalSpacing(ROW_SPACING)
        add_row(form, "単位系 (units)", self.units)
        add_row(form, "原子の形式 (atom_style)", self.atom_style)
        add_row(form, "data ファイル", file_row(self.data_file, "data ファイル", "LAMMPS data (*.data *.lmp *.lammps data.*);;すべて (*)"))
        add_row(form, "型番号の元素", self.type_elements)
        add_row(form, "pair_style", self.pair_style)
        form.addRow(label(""), force_field_help)
        add_row(form, "pair_coeff", self.pair_coeff)
        add_row(form, "写すファイル", file_row(self.potential_files, "写すファイル", "", append=True))
        add_row(form, "read_data の前の行", self.style_commands)
        add_row(form, "pair_coeff の後の行", self.extra_commands)
        add_row(form, "乱数の種", self.seed)
        form.addRow(label(""), note)

        for w in (self.units, self.atom_style):
            w.currentTextChanged.connect(self._emit)
        for w in (self.data_file, self.type_elements, self.pair_style):
            w.textChanged.connect(self._emit)
        for w in (self.pair_coeff, self.potential_files, self.style_commands, self.extra_commands):
            w.textChanged.connect(self._emit)
        self.seed.valueChanged.connect(self._emit)

    def method(self) -> LammpsMethod:
        return LammpsMethod(units=self.units.currentData(), atom_style=self.atom_style.currentText().strip() or "atomic",
                            data_file=self.data_file.text().strip(), type_elements=parse_words(self.type_elements.text()),
                            pair_style=self.pair_style.text().strip(), pair_coeff=self.pair_coeff.toPlainText().strip(),
                            potential_files=parse_lines(self.potential_files.toPlainText()),
                            style_commands=self.style_commands.toPlainText().strip(), extra_commands=self.extra_commands.toPlainText().strip(),
                            seed=self.seed.value(), thermo_pressure_tensor=self._thermo_pressure_tensor)

    def set_method(self, m: LammpsMethod) -> None:
        self._thermo_pressure_tensor = m.thermo_pressure_tensor
        self.units.setCurrentIndex(max(0, self.units.findData(m.units))); self.atom_style.setCurrentText(m.atom_style)
        self.data_file.setText(m.data_file); self.type_elements.setText(" ".join(m.type_elements)); self.pair_style.setText(m.pair_style)
        self.pair_coeff.setPlainText(m.pair_coeff); self.potential_files.setPlainText("\n".join(m.potential_files))
        self.style_commands.setPlainText(m.style_commands); self.extra_commands.setPlainText(m.extra_commands); self.seed.setValue(m.seed)
        self._emit()

    def _emit(self, *_) -> None:
        self.changed.emit()
