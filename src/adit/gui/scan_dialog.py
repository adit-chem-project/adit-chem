
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from PySide6.QtWidgets import (QDialog, QFileDialog, QFormLayout, QHBoxLayout, QComboBox, QLabel, QLineEdit,
                               QMessageBox, QPushButton, QVBoxLayout, QWidget)

from adit.lang import L
from adit.gui.style import GROUP_SPACING, PANEL_MARGIN, ROW_SPACING
from adit.gui.widgets import add_row, confirm_overwrite

from adit.web.scan_choices import OTHER, Choice, choices_for  # noqa: E402,F401


class ScanDialog(QDialog):

    def __init__(self, code: str, periodic: bool, default_dir: str, make_spec: Callable, cfg, parent: QWidget | None = None):
        super().__init__(parent)
        self.setWindowTitle(L("1 つの条件を変えて一括生成", "Parameter scan"))
        self.make_spec, self.cfg = make_spec, cfg
        self.out_dir: Path | None = None
        self.backup_dir: Path | None = None
        self.dirs: list[Path] = []
        self.choices = choices_for(code, periodic)

        self.item = QComboBox()
        for c in self.choices:
            self.item.addItem(c.label, c.path)
        self.purpose = QLabel(); self.purpose.setObjectName("hint"); self.purpose.setWordWrap(True)
        self.path = QLineEdit(); self.path.setPlaceholderText(L("例 method.scc_tolerance", "e.g. method.scc_tolerance"))
        self.values = QLineEdit()
        self.folder = QLineEdit(default_dir)
        self.browse = QPushButton(L("参照…", "Browse…"))
        self.cost = QLabel(L("値の数だけ計算が増えます。大きな系は、この PC ではなくクラスタで実行してください。",
                             "Each value is a separate calculation. Run large systems on a cluster rather than on this PC."))
        self.cost.setObjectName("hint"); self.cost.setWordWrap(True)
        self.btn_cancel = QPushButton(L("キャンセル", "Cancel"))
        self.btn_ok = QPushButton(L("生成", "Generate")); self.btn_ok.setObjectName("primary"); self.btn_ok.setDefault(True)

        self.form = QFormLayout(); self.form.setVerticalSpacing(ROW_SPACING)
        add_row(self.form, L("変える項目", "Parameter"), self.item,
                help_text=L("1 つだけ値を変え、残りの設定はいまの画面のまま、値ごとに入力を作ります。",
                            "Only this setting changes; everything else is taken from the current settings."))
        self.form.addRow(self.purpose)
        add_row(self.form, L("項目の場所", "Setting path"), self.path,
                help_text=L("spec.json の中の場所をドットでつないで書きます (例 method.scc_tolerance)。",
                            "The dotted path of the setting inside spec.json (e.g. method.scc_tolerance)."))
        add_row(self.form, L("値", "Values"), self.values,
                help_text=L("カンマで区切って 2 つ以上。値ごとにディレクトリを 1 つ作ります。",
                            "Two or more, separated by commas. One directory is created per value."))
        row = QHBoxLayout(); row.addWidget(self.folder, 1); row.addWidget(self.browse)
        add_row(self.form, L("保存先", "Output directory"), row,
                help_text=L("この中に値ごとのディレクトリ (例 ecutwfc_30) を作ります。",
                            "One directory per value (e.g. ecutwfc_30) is created inside it."))
        from adit.gui.progress import EmptyState
        self.empty = EmptyState(L("値ごとの入力はまだありません", "No inputs yet"),
                                L("値をカンマで区切って書き、「生成」を押すと、保存先の中に値ごとのディレクトリを作ります",
                                  "Type the values separated by commas and press Generate; one directory per value is made in the output directory"),
                                L("例の値を入れる", "Use the example values"))
        self.empty.clicked.connect(self.use_example)
        buttons = QHBoxLayout(); buttons.addStretch(); buttons.addWidget(self.btn_cancel); buttons.addWidget(self.btn_ok)
        lay = QVBoxLayout(self); lay.setContentsMargins(PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN); lay.setSpacing(GROUP_SPACING)
        lay.addLayout(self.form); lay.addWidget(self.empty, 1); lay.addWidget(self.cost); lay.addLayout(buttons)
        self.setMinimumWidth(680)

        self.item.currentIndexChanged.connect(self._on_item)
        self.browse.clicked.connect(self._browse)
        self.btn_cancel.clicked.connect(self.reject)
        self.btn_ok.clicked.connect(self.generate)
        self.values.textChanged.connect(self._on_values)
        self._on_item()

    def _on_values(self, *_) -> None:
        n = len([v for v in self.values.text().split(",") if v.strip()])
        if n == 0:
            self.empty.set_texts(L("値ごとの入力はまだありません", "No inputs yet"),
                                 L("値をカンマで区切って書き、「生成」を押すと、保存先の中に値ごとのディレクトリを作ります",
                                   "Type the values separated by commas and press Generate; one directory per value is made in the output directory"),
                                 L("例の値を入れる", "Use the example values"))
        else:
            self.empty.set_texts(L(f"{n} 個の値 → {n} 個のディレクトリ", f"{n} values → {n} directories"),
                                 L("「生成」を押すと、保存先の中に値ごとのディレクトリを作ります",
                                   "Press Generate to make one directory per value in the output directory"), "")

    def use_example(self) -> None:
        self.values.setText(self.current_choice().example)

    def current_choice(self) -> Choice:
        return self.choices[max(0, self.item.currentIndex())]

    def _on_item(self, *_) -> None:
        c = self.current_choice()
        self.purpose.setText(c.purpose)
        self.values.setPlaceholderText(c.example)
        self.form.setRowVisible(self.path, c.path == OTHER)

    def _browse(self) -> None:
        start = self.folder.text().strip() or str(Path.home())
        d = QFileDialog.getExistingDirectory(self, L("保存先", "Output directory"), str(Path(start).parent))
        if d:
            self.folder.setText(d)

    def scan_text(self) -> str:
        c = self.current_choice()
        path = self.path.text().strip() if c.path == OTHER else c.path
        return f"{path}={self.values.text()}"

    def generate(self) -> None:
        from pydantic import ValidationError as PydanticError
        from adit.config import ConfigError
        from adit.project import ProjectError
        from adit.scan import ScanError, parse_scan, write_scan

        c = self.current_choice()
        if c.path == OTHER and not self.path.text().strip():
            QMessageBox.warning(self, L("生成できません", "Cannot generate"),
                                L("項目の場所を書いてください (例 method.scc_tolerance)。", "Type the setting path (e.g. method.scc_tolerance)."))
            return
        if not self.folder.text().strip():
            QMessageBox.warning(self, L("生成できません", "Cannot generate"),
                                L("保存先を指定してください。", "Choose an output directory."))
            return
        out = Path(self.folder.text().strip()).expanduser()
        try:
            scan = parse_scan(self.scan_text())
            spec = self.make_spec()
            overwrite, keep = confirm_overwrite(self, out)
            if overwrite is None:
                return
            dirs = write_scan(spec, self.cfg, out, scan, overwrite=overwrite)
            self.backup_dir = keep.finish() if keep is not None else None
        except (ScanError, ProjectError, ConfigError, ValueError, OSError) as ex:
            if isinstance(ex, PydanticError):
                from adit.validate_types import friendly_pydantic
                ex = ProjectError(friendly_pydantic(ex))
            QMessageBox.critical(self, L("生成できません", "Cannot generate"), str(ex))
            return
        self.out_dir, self.dirs = out, dirs
        names = "\n".join(f"  {d.name}" for d in dirs)
        QMessageBox.information(self, L("生成しました", "Generated"), L(
            f"{out} の中に {len(dirs)} 個のディレクトリを作りました。\n{names}\n\n"
            "次にすること\n"
            "1. それぞれの計算を実行します。この PC なら各ディレクトリで bash submit.sh、"
            "クラスタならジョブとして投入します (手順は各ディレクトリの README.txt にあります)。\n"
            "2. すべて終わったら、解析タブで「解析を実行」を押します。計算結果のディレクトリ欄には、この保存先を入れておきました。",
            f"Created {len(dirs)} directories in {out}:\n{names}\n\n"
            "Next steps\n"
            "1. Run each calculation: on this PC, run bash submit.sh in each directory; on a cluster, submit each one as a job "
            "(see README.txt in each directory).\n"
            "2. When they have all finished, click \"Run analysis\" in the Analysis tab. The Run directory field is already filled in."))
        self.accept()
