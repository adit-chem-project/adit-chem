
from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtWidgets import (QComboBox, QDialog, QFileDialog, QFormLayout, QHBoxLayout, QLabel, QLineEdit,
                               QMessageBox, QPushButton, QStackedWidget, QVBoxLayout, QWidget)

from adit.gui.style import GROUP_SPACING, PANEL_MARGIN, ROW_SPACING
from adit.gui.widgets import add_row
from adit.lang import L


def _hint(text: str) -> QLabel:
    w = QLabel(text); w.setObjectName("hint"); w.setWordWrap(True); return w


class ConvertDialog(QDialog):
    def __init__(self, make_spec, cfg, parent=None):
        super().__init__(parent)
        self.make_spec, self.cfg = make_spec, cfg
        self.setWindowTitle(L("構造・計算コードの変換", "Convert a structure, or retarget it to another code"))
        self.kind = QComboBox()
        self.kind.addItem(L("構造ファイルの形式を変換", "Convert a structure file"), "structure")
        self.kind.addItem(L("別の計算コード用に変換", "Convert for another calculation code"), "calculation")
        self.stack = QStackedWidget(); self.stack.addWidget(self._structure_page()); self.stack.addWidget(self._calculation_page())
        self.status = _hint("")
        cancel = QPushButton(L("閉じる", "Close")); cancel.clicked.connect(self.reject)
        self.run = QPushButton(L("変換", "Convert")); self.run.setObjectName("primary"); self.run.clicked.connect(self.convert)
        buttons = QHBoxLayout(); buttons.addStretch(); buttons.addWidget(cancel); buttons.addWidget(self.run)
        form = QFormLayout(); form.setVerticalSpacing(ROW_SPACING); add_row(form, L("変換するもの", "What to convert"), self.kind)
        lay = QVBoxLayout(self); lay.setContentsMargins(PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN); lay.setSpacing(GROUP_SPACING)
        lay.addLayout(form); lay.addWidget(self.stack); lay.addWidget(self.status); lay.addLayout(buttons)
        self.kind.currentIndexChanged.connect(self.stack.setCurrentIndex)
        self.kind.currentIndexChanged.connect(self._update_run_state)
        self._update_run_state()
        self.setMinimumWidth(720)

    def _update_run_state(self) -> None:
        self.run.setEnabled(self.kind.currentData() == "structure" or self.template.count() > 0)

    def _path_row(self, edit: QLineEdit, folder: bool, save: bool = False) -> QWidget:
        box = QWidget(); row = QHBoxLayout(box); row.setContentsMargins(0, 0, 0, 0); row.addWidget(edit, 1)
        button = QPushButton(L("参照…", "Browse…")); row.addWidget(button)
        def choose():
            if folder:
                p = QFileDialog.getExistingDirectory(self, L("出力ディレクトリ", "Output directory"), edit.text() or str(Path.home()))
            elif save:
                p, _ = QFileDialog.getSaveFileName(self, L("変換後の構造", "Converted structure"), edit.text() or str(Path.home()))
            else:
                p, _ = QFileDialog.getOpenFileName(self, L("変換元の構造", "Source structure"), edit.text() or str(Path.home()))
            if p: edit.setText(p)
        button.clicked.connect(choose); return box

    def _structure_page(self) -> QWidget:
        page = QWidget(); form = QFormLayout(page); form.setVerticalSpacing(ROW_SPACING)
        self.struct_source, self.struct_output, self.cell = QLineEdit(), QLineEdit(), QLineEdit()
        self.input_format, self.output_format = QLineEdit(), QLineEdit()
        self.cell.setPlaceholderText(L("例: 20 または 20,20,30", "e.g. 20 or 20,20,30"))
        add_row(form, L("変換元", "Source"), self._path_row(self.struct_source, False))
        add_row(form, L("変換後", "Output"), self._path_row(self.struct_output, False, True))
        add_row(form, L("入力形式 (任意)", "Input format (optional)"), self.input_format)
        add_row(form, L("出力形式 (任意)", "Output format (optional)"), self.output_format)
        add_row(form, L("セル [Å]", "Cell [Å]"), self.cell)
        form.addRow(_hint(L("形式はファイル名から判定します。XYZ などセルを持たない構造を POSCAR へ変換するときだけ、セルを指定してください。ADIT はセルの大きさを決めません。",
                              "Formats are inferred from filenames. Specify a cell only when converting a cell-free format such as XYZ to POSCAR. ADIT does not choose the cell size.")))
        return page

    def _calculation_page(self) -> QWidget:
        from adit.templates import list_templates, template_dirs
        page = QWidget(); form = QFormLayout(page); form.setVerticalSpacing(ROW_SPACING)
        self.template = QComboBox()
        for t in list_templates(self.cfg):
            if not t.error: self.template.addItem(f"{t.name} — {t.code}" + (f" — {t.comment}" if t.comment else ""), t.name)
        self.calc_output = QLineEdit()
        add_row(form, L("変換先の雛形", "Target template"), self.template)
        add_row(form, L("出力ディレクトリ", "Output directory"), self._path_row(self.calc_output, True))
        form.addRow(_hint(L("構造・固定原子・速度・計算の種類・MD の条件は現在の画面から保ちます。計算手法、力場、カットオフ、擬ポテンシャル、k 点は変換せず、変換先の雛形から取ります。",
                              "The structure, fixed atoms, velocities, task, and MD settings are preserved from the current screen. The method, force field, cutoffs, pseudopotentials, and k-points are not translated; they come from the target template.")))
        if not self.template.count():
            form.addRow(_hint(L(f"変換先の雛形がありません。先に「雛形として保存」で作ってください。置き場所: {', '.join(map(str, template_dirs(self.cfg)))}",
                                  f"No target templates are available. Create one with Save as template first. Folders: {', '.join(map(str, template_dirs(self.cfg)))}")))
        return page

    @staticmethod
    def _cell(text: str):
        if not text.strip(): return None
        try: values = [float(x.strip()) for x in text.split(",")]
        except ValueError as ex: raise ValueError(L("セルは A または A,B,C の形で入力してください", "enter the cell as A or A,B,C")) from ex
        if len(values) == 1: values *= 3
        if len(values) != 3: raise ValueError(L("セルは A または A,B,C の形で入力してください", "enter the cell as A or A,B,C"))
        return tuple(values)

    def convert(self) -> None:
        try:
            if self.kind.currentData() == "structure":
                from adit.convert import convert_structure
                p = convert_structure(self.struct_source.text().strip(), self.struct_output.text().strip(),
                                      input_format=self.input_format.text().strip() or None,
                                      output_format=self.output_format.text().strip() or None,
                                      cell=self._cell(self.cell.text()))
                msg = L(f"構造を書きました: {p}", f"Wrote the structure: {p}")
            else:
                from adit.convert import REPORT_FILE, retarget_spec
                from adit.project import write_project
                from adit.templates import load_template
                source = self.make_spec(); target = load_template(self.template.currentData(), source, self.cfg)
                converted, report = retarget_spec(source, target)
                written = write_project(converted, self.cfg, self.calc_output.text().strip(),
                                        extra_readme=[L("== 計算コード間の変換 ==", "== Conversion between calculation codes =="),
                                                      f"  {source.method.code} -> {converted.method.code}", f"  {report['rule']}",
                                                      L(f"  詳細: {REPORT_FILE}", f"  Details: {REPORT_FILE}"), ""],
                                        extra_texts={REPORT_FILE: json.dumps(report, ensure_ascii=False, indent=2) + "\n"})
                msg = L(f"{converted.method.code} の入力を生成しました ({len(written)} ファイル): {self.calc_output.text().strip()}",
                        f"Generated {converted.method.code} input ({len(written)} files): {self.calc_output.text().strip()}")
        except Exception as ex:
            QMessageBox.critical(self, L("変換できません", "Cannot convert"), str(ex)); return
        self.status.setText(msg)
