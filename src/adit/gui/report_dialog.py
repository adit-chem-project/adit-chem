
from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout,
                               QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit, QPushButton, QVBoxLayout)

from adit.gui import report_fields as R
from adit.gui.help import help_for
from adit.gui.i18n import translate_widgets
from adit.gui.style import PANEL_MARGIN, ROW_SPACING
from adit.gui.widgets import add_row
from adit.lang import L


class ReportDialog(QDialog):

    def __init__(self, run_dirs: list[str] | None = None, parent=None):
        super().__init__(parent)
        self.setWindowTitle(R.lab("report"))
        self.dirs = QPlainTextEdit("\n".join(run_dirs or []))
        self.dirs.setPlaceholderText(L("1 行に 1 つ。複数入れると 1 枚の表にまとめます",
                                       "one per line; several directories are combined into one table"))
        self.dirs.setFixedHeight(80)
        self.browse = QPushButton(L("参照…", "Browse…"))
        self.lang = QComboBox()
        for value, ja, en in R.LANGUAGES:
            self.lang.addItem(L(ja, en), value)
        self.methods = QLineEdit(); self.conditions = QLineEdit()
        self.results = QLineEdit(); self.bundle = QLineEdit()
        for widget, hint in ((self.methods, "methods.md"), (self.conditions, "conditions.csv"),
                             (self.results, "results.csv"), (self.bundle, "pack.zip")):
            widget.setPlaceholderText(L(f"空欄なら作りません (例 {hint})", f"empty = not created (e.g. {hint})"))
        self.check = QCheckBox(R.lab("rep_check"))
        self.check.setToolTip(_help(R.LABELS["rep_check"][0]))
        self.out = QPlainTextEdit(); self.out.setReadOnly(True); self.out.setFixedHeight(160)
        self.out.setStyleSheet("font-family: monospace;")

        form = QFormLayout(); form.setVerticalSpacing(ROW_SPACING)
        row = QHBoxLayout(); row.addWidget(self.dirs); row.addWidget(self.browse)
        add_row(form, R.LABELS["rep_dirs"][0], row, required=True, help_text=_help(R.LABELS["rep_dirs"][0]))
        for key, widget in (("rep_lang", self.lang), ("rep_methods", self.methods),
                            ("rep_conditions", self.conditions), ("rep_results", self.results),
                            ("rep_bundle", self.bundle)):
            add_row(form, R.LABELS[key][0], widget, help_text=_help(R.LABELS[key][0]))
        form.addRow(self.check)
        self.note = QLabel(L("コマンドの adit-report と同じものを作ります。判定はしません "
                             "(ハッシュの照合だけは「一致 / 不一致 / 記録なし」を出します)。",
                             "This produces the same output as the adit-report command. Nothing is judged "
                             "(only the fingerprint check reports match / differ / no record)."))
        self.note.setObjectName("hint"); self.note.setWordWrap(True)
        buttons = QDialogButtonBox()
        self.run_button = buttons.addButton(R.lab("rep_run"), QDialogButtonBox.ButtonRole.ActionRole)
        buttons.addButton(QDialogButtonBox.StandardButton.Close)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN)
        layout.addLayout(form); layout.addWidget(self.note); layout.addWidget(self.out); layout.addWidget(buttons)
        self.browse.clicked.connect(self._add_dir)
        self.run_button.clicked.connect(self.build)
        buttons.rejected.connect(self.reject)
        translate_widgets(self)

    def _add_dir(self) -> None:
        path = QFileDialog.getExistingDirectory(self, R.lab("rep_dirs"))
        if path:
            text = self.dirs.toPlainText().strip()
            self.dirs.setPlainText((text + "\n" + path).strip())

    def fields(self) -> dict:
        return {"rep_dirs": self.dirs.toPlainText(), "rep_lang": str(self.lang.currentData() or ""),
                "rep_methods": self.methods.text(), "rep_conditions": self.conditions.text(),
                "rep_results": self.results.text(), "rep_bundle": self.bundle.text(),
                "rep_check": "on" if self.check.isChecked() else ""}

    def build(self) -> R.ReportOutcome | None:
        from adit.report import ReportError

        try:
            outcome = R.build(R.request_from_fields(self.fields()))
        except (R.ReportFieldError, ReportError, OSError) as ex:
            self.out.setPlainText(str(ex))
            return None
        self.out.setPlainText("\n".join([*outcome.lines, "", outcome.methods_text]).strip())
        return outcome


def _help(label: str) -> str:
    h = help_for(label)
    return h.text() if h else ""
