
from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtWidgets import (QComboBox, QDialog, QFileDialog, QFormLayout, QHBoxLayout,
                               QLabel, QLineEdit, QPlainTextEdit, QPushButton,
                               QStackedWidget, QVBoxLayout, QWidget)

from adit.gui.style import GROUP_SPACING, PANEL_MARGIN, ROW_SPACING
from adit.gui.widgets import add_row
from adit.lang import L


class InputAuditDialog(QDialog):
    """Three read/review workflows; only import creates a new review directory."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(L("入力と結果の点検", "Inspect inputs and results"))
        self.setMinimumSize(760, 610)

        self.kind = QComboBox()
        self.kind.addItem(L("既存入力の読み込み", "Import existing input"), "import")
        self.kind.addItem(L("生成入力の照合", "Verify generated input"), "verify")
        self.kind.addItem(L("複数の計算の前提の突き合わせ", "Compare prerequisites across runs"), "audit")
        selector = QFormLayout(); selector.setVerticalSpacing(ROW_SPACING)
        add_row(selector, L("点検するもの", "Inspect"), self.kind)

        self.pages = QStackedWidget()
        self.pages.addWidget(self._import_page())
        self.pages.addWidget(self._verify_page())
        self.pages.addWidget(self._audit_page())
        self.kind.currentIndexChanged.connect(self._change_kind)

        self.summary = QLabel("")
        self.summary.setObjectName("hint"); self.summary.setWordWrap(True)
        self.details = QPlainTextEdit(); self.details.setReadOnly(True)
        self.details.setPlaceholderText(L("点検結果がここに表示されます。", "Inspection results appear here."))
        self.inspect_button = QPushButton(L("点検する", "Inspect"))
        self.inspect_button.setObjectName("primary")
        self.inspect_button.clicked.connect(self.inspect)
        close_button = QPushButton(L("閉じる", "Close")); close_button.clicked.connect(self.reject)
        buttons = QHBoxLayout(); buttons.addStretch(); buttons.addWidget(close_button); buttons.addWidget(self.inspect_button)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN)
        layout.setSpacing(GROUP_SPACING)
        layout.addLayout(selector); layout.addWidget(self.pages)
        layout.addWidget(self.summary); layout.addWidget(self.details, 1); layout.addLayout(buttons)

    def _path_row(self, edit: QLineEdit, *, on_choose=None) -> QWidget:
        box = QWidget(); row = QHBoxLayout(box)
        row.setContentsMargins(0, 0, 0, 0); row.addWidget(edit, 1)
        button = QPushButton(L("参照…", "Browse…")); row.addWidget(button)

        def choose() -> None:
            selected = QFileDialog.getExistingDirectory(
                self, L("ディレクトリを選択", "Choose a directory"),
                edit.text() or str(Path.home()))
            if selected:
                edit.setText(selected)
                if on_choose is not None:
                    on_choose(selected)

        button.clicked.connect(choose)
        return box

    @staticmethod
    def _hint(ja: str, en: str) -> QLabel:
        label = QLabel(L(ja, en)); label.setObjectName("hint"); label.setWordWrap(True)
        return label

    def _import_page(self) -> QWidget:
        page = QWidget(); form = QFormLayout(page); form.setVerticalSpacing(ROW_SPACING)
        self.import_source = QLineEdit()
        self.import_review = QLineEdit()
        self.import_code = QComboBox()
        for label, value in ((L("自動判定", "Auto-detect"), None), ("VASP", "vasp"),
                             ("Quantum ESPRESSO", "qe"), ("LAMMPS", "lammps")):
            self.import_code.addItem(label, value)
        add_row(form, L("既存入力の場所", "Existing input directory"),
                self._path_row(self.import_source, on_choose=self._suggest_review))
        add_row(form, L("コード", "Code"), self.import_code)
        add_row(form, L("新しい点検先", "New review directory"), self.import_review)
        form.addRow(self._hint(
            "点検先には import_report.json と、読めた場合だけ draft_spec.json を作ります。既存入力は変更しません。下書きは再生成の同一性や実行可能性を保証しません。",
            "Creates import_report.json and, only when parsing succeeds, draft_spec.json in a new review directory. The source is unchanged. A draft does not establish equivalent regeneration or executability."))
        return page

    def _suggest_review(self, selected: str) -> None:
        if not self.import_review.text().strip():
            source = Path(selected)
            self.import_review.setText(str(source.parent / (source.name + "_review")))

    def _verify_page(self) -> QWidget:
        page = QWidget(); form = QFormLayout(page); form.setVerticalSpacing(ROW_SPACING)
        self.verify_dir = QLineEdit()
        add_row(form, L("生成ディレクトリ", "Generated directory"), self._path_row(self.verify_dir))
        form.addRow(self._hint(
            "spec.json と、実際に生成した入力ファイルを読み比べます。読み戻せる欄だけが対象で、計算全体が同じか、そのまま走るかは判定しません。ファイルは書き換えません。",
            "Checks only explicitly mapped fields in spec.json against the generated native inputs. Read-only; this does not establish whole-calculation equivalence or executability."))
        return page

    def _audit_page(self) -> QWidget:
        page = QWidget(); form = QFormLayout(page); form.setVerticalSpacing(ROW_SPACING)
        self.audit_first = QLineEdit(); self.audit_second = QLineEdit()
        self.audit_kind = QComboBox()
        self.audit_kind.addItem(L("MSD の前提条件", "MSD prerequisites"), "msd")
        self.audit_kind.addItem(L("エネルギーの前提条件", "Energy prerequisites"), "energy")
        add_row(form, L("基準の計算", "Reference run"), self._path_row(self.audit_first))
        add_row(form, L("比べる計算", "Other run"), self._path_row(self.audit_second))
        add_row(form, L("点検の種類", "Audit kind"), self.audit_kind)
        form.addRow(self._hint(
            "出力を読み取って機械的な前提条件を点検します。図やファイルは作らず、異なる計算法の結果を科学的に直接比較してよいかは判定しません。",
            "Reads outputs to check mechanical prerequisites. It creates no figures or files and does not judge whether results from different methods are scientifically comparable."))
        return page

    def _change_kind(self, index: int) -> None:
        self.pages.setCurrentIndex(index)
        self.summary.clear(); self.details.clear()

    def inspect(self) -> None:
        """Call the same importer, verifier and run audit used by the CLI."""
        try:
            kind = self.kind.currentData()
            if kind == "import":
                from adit.native_workflow import write_import_review

                result, output = write_import_review(
                    self.import_source.text().strip(), self.import_review.text().strip(),
                    code=self.import_code.currentData())
                self.summary.setText(L(
                    f"点検記録: {output / 'import_report.json'}。" +
                    ("下書きの計算設定 (draft_spec.json) も作りました。既定値と外部パラメータを確認してください。" if result.spec is not None else
                     "読み取れない条件があるため下書きの計算設定 (draft_spec.json) は作っていません。"),
                    f"Review report: {output / 'import_report.json'}. " +
                    ("A draft spec was also created; review defaults and external parameters." if result.spec is not None else
                     "No draft spec was created because some settings could not be imported.")))
                details = result.report()
            elif kind == "verify":
                from adit.native_verify import verify_generated_bundle
                from adit.native_workflow import verification_passed

                source = self.verify_dir.text().strip()
                if not source:
                    raise ValueError(L("生成ディレクトリを指定してください", "Specify a generated directory"))
                result = verify_generated_bundle(source)
                passed = verification_passed(result)
                self.summary.setText(L(
                    f"対応項目のみの点検: {'問題なし' if passed else '問題あり'}。一致 {len(result.preserved)} 件、不一致 {len(result.mismatched)} 件、未確認 {len(result.unverifiable)} 件。未確認は検証済みを意味しません。",
                    f"Mapped-field check only: {'passed' if passed else 'failed'}. {len(result.preserved)} matched, {len(result.mismatched)} mismatched, {len(result.unverifiable)} unverified. Unverified does not mean verified."))
                details = result.report()
            else:
                from adit.analysis.comparability import audit_run_dirs

                first, second = self.audit_first.text().strip(), self.audit_second.text().strip()
                if not first or not second:
                    raise ValueError(L("比べる二つの計算ディレクトリを指定してください",
                                       "Specify both run directories to audit"))
                result = audit_run_dirs(
                    [first, second],
                    require_md=self.audit_kind.currentData() == "msd")
                self.summary.setText(L("機械的な点検結果です。化学的に比べてよいかは判断していません。",
                                       "Mechanical check only; chemical comparability is not judged."))
                self.details.setPlainText(result.summary_text())
                return
            self.details.setPlainText(json.dumps(details, ensure_ascii=False, indent=2))
        except Exception as exc:
            self.summary.setText(L(f"点検できません: {exc}", f"Cannot inspect: {exc}"))
            self.details.clear()
