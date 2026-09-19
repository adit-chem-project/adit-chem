
from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (QAbstractItemView, QDialog, QDialogButtonBox, QFileDialog, QFormLayout, QHBoxLayout, QHeaderView, QLabel,
                               QLineEdit, QPushButton, QTableWidget, QTableWidgetItem, QVBoxLayout)

from adit.gui import analysis_fields as AF
from adit.gui.widgets import add_row
from adit.lang import L

START_ROWS = 4


class CompareDialog(QDialog):
    def __init__(self, base: str = "", parent=None):
        super().__init__(parent)
        self.setWindowTitle(L("組にして比べる", "Compare runs"))
        self.setMinimumWidth(640)
        self.base = QLineEdit(base)
        self.browse = QPushButton("参照…")
        row = QHBoxLayout(); row.addWidget(self.base, 1); row.addWidget(self.browse)
        form = QFormLayout()
        add_row(form, AF.LABELS["compare_base"][0], row)
        self.hint = QLabel(L("係数 ν は生成物を正、反応物を負にします (ΔE = ΣνE)。ディレクトリは基準のディレクトリからの相対パスか絶対パス。"
                             "反応の名前が空の行は、上の行と同じ反応に入ります。",
                             "The coefficient ν is positive for products and negative for reactants (ΔE = ΣνE). Directories are relative to the base "
                             "directory or absolute. A row with an empty reaction name belongs to the reaction above it."))
        self.hint.setObjectName("hint"); self.hint.setWordWrap(True)
        self.table = QTableWidget(START_ROWS, 3)
        self.table.setHorizontalHeaderLabels([AF.lab("rx_name"), AF.lab("rx_nu"), AF.lab("rx_dir")])
        self.table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self.table.verticalHeader().setVisible(False)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.btn_add = QPushButton(L("行を追加", "Add row")); self.btn_dir = QPushButton(L("ディレクトリを選んで行を追加…", "Add a row from a folder…"))
        self.btn_del = QPushButton(L("選んだ行を削除", "Remove selected rows")); self.btn_load = QPushButton(L("compare.json を読む", "Load compare.json"))
        tools = QHBoxLayout()
        for b in (self.btn_add, self.btn_dir, self.btn_del, self.btn_load):
            tools.addWidget(b)
        tools.addStretch()
        self.error = QLabel(""); self.error.setObjectName("error"); self.error.setWordWrap(True); self.error.hide()
        self.buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setText(AF.LABELS["compare_run"][0])
        self.buttons.button(QDialogButtonBox.StandardButton.Cancel).setText(L("キャンセル", "Cancel"))
        lay = QVBoxLayout(self); lay.setContentsMargins(16, 16, 16, 16); lay.setSpacing(10)
        lay.addLayout(form); lay.addWidget(self.hint); lay.addWidget(self.table, 1); lay.addLayout(tools); lay.addWidget(self.error); lay.addWidget(self.buttons)
        self.browse.clicked.connect(self._browse_base)
        self.btn_add.clicked.connect(self.add_row)
        self.btn_dir.clicked.connect(self._add_dir)
        self.btn_del.clicked.connect(self._delete_rows)
        self.btn_load.clicked.connect(self.load_compare_json)
        self.base.textChanged.connect(self._update_load)
        self.buttons.accepted.connect(self._try_accept)
        self.buttons.rejected.connect(self.reject)
        self._update_load()
        self._reactions = None
        from adit.gui.i18n import translate_widgets
        translate_widgets(self)

    def add_row(self, values: tuple[str, str, str] = ("", "", "")) -> None:
        r = self.table.rowCount()
        self.table.insertRow(r)
        for c, v in enumerate(values):
            self.table.setItem(r, c, QTableWidgetItem(v))

    def rows(self) -> list[tuple[str, str, str]]:
        out = []
        for r in range(self.table.rowCount()):
            out.append(tuple((self.table.item(r, c).text() if self.table.item(r, c) else "") for c in range(3)))
        return out

    def set_rows(self, rows: list[tuple[str, str, str]]) -> None:
        self.table.setRowCount(0)
        for v in rows:
            self.add_row(v)
        for _ in range(max(0, START_ROWS - len(rows))):
            self.add_row()

    def _delete_rows(self) -> None:
        for r in sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True):
            self.table.removeRow(r)

    def _add_dir(self) -> None:
        d = QFileDialog.getExistingDirectory(self, AF.lab("rx_dir"), self.base.text() or str(Path.home()))
        if not d:
            return
        base = Path(self.base.text().strip()).expanduser()
        try:
            d = str(Path(d).relative_to(base)) if self.base.text().strip() else d
        except ValueError:
            pass
        for r in range(self.table.rowCount()):
            if not any(self.table.item(r, c) and self.table.item(r, c).text().strip() for c in range(3)):
                self.table.setItem(r, 2, QTableWidgetItem(d)); return
        self.add_row(("", "", d))

    def _browse_base(self) -> None:
        d = QFileDialog.getExistingDirectory(self, AF.lab("compare_base"), self.base.text() or str(Path.home()))
        if d:
            self.base.setText(d)

    def _compare_json(self) -> Path | None:
        from adit.analysis.compare import COMPARE_FILE
        t = self.base.text().strip()
        p = Path(t).expanduser() / COMPARE_FILE if t else None
        return p if p is not None and p.is_file() else None

    def _update_load(self) -> None:
        self.btn_load.setEnabled(self._compare_json() is not None)

    def load_compare_json(self) -> bool:
        from adit.analysis.compare import CompareError, load_compare_file
        p = self._compare_json()
        if p is None:
            return False
        try:
            self.set_rows(AF.rows_from_reactions(load_compare_file(p)))
        except (CompareError, ValueError, OSError, KeyError) as ex:
            self._show_error(L(f"compare.json を読めません: {ex}", f"cannot read compare.json: {ex}")); return False
        self.error.hide()
        return True

    def _show_error(self, text: str) -> None:
        from adit.gui.style import current_theme, tokens_for
        self.error.setStyleSheet(f"color: {tokens_for(current_theme()).ng};")
        self.error.setText(text); self.error.show()

    def _try_accept(self) -> None:
        from adit.analysis.compare import CompareError
        base = self.base.text().strip()
        if not base or not Path(base).expanduser().is_dir():
            self._show_error(L(f"{AF.lab('compare_base')}: ディレクトリがありません ({base!r})", f"{AF.lab('compare_base')}: directory not found ({base!r})")); return
        try:
            self._reactions = AF.reactions_from_rows(self.rows())
        except (AF.FieldError, CompareError) as ex:
            self._show_error(str(ex)); return
        self.accept()

    def reactions(self):
        return self._reactions

    def base_dir(self) -> str:
        return self.base.text().strip()
