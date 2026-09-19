
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import QLabel, QPlainTextEdit, QTabWidget, QVBoxLayout, QWidget

from adit.lang import L
from adit.gui.i18n import tr
from adit.project import ProjectFiles

ORDER = ["dftb_in.hsd", "geometry.gen", "INCAR", "POSCAR", "KPOINTS", "potcar.spec", "make_potcar.sh", "struct.xyz", "xtb.inp", "pw.in", "orca.inp", "analyze.py",
         "submit.sh", "spec.json", "README.txt"]


class PreviewPanel(QWidget):
    fix_requested = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        from PySide6.QtWidgets import QHBoxLayout, QPushButton
        from adit.gui.progress import EmptyState
        self.empty = EmptyState(L("生成できません", "Cannot generate"), "", L("欄へ移動", "Go to the field"))
        self.empty.clicked.connect(self.fix_requested.emit); self.empty.hide()
        self.status = QLabel("")
        self.detail = QLabel(""); self.detail.setObjectName("hint"); self.detail.setWordWrap(True)
        self.origin = QLabel(""); self.origin.setObjectName("hint"); self.origin.setWordWrap(True)
        self.btn_clear_origin = QPushButton(L("外す", "Clear")); self.btn_clear_origin.setObjectName("link")
        self.btn_clear_origin.setToolTip(L("前の計算の続きと雛形の控えを外します (画面の値はそのまま)", "Drop the continuation and template records (the fields keep their values)"))
        self.origin_row = QWidget(); self.origin_row.setObjectName("rowbox")
        orl = QHBoxLayout(self.origin_row); orl.setContentsMargins(0, 0, 0, 0); orl.setSpacing(6)
        orl.addWidget(self.origin, 1); orl.addWidget(self.btn_clear_origin, 0, Qt.AlignmentFlag.AlignTop)
        self.origin_row.hide()
        self.tabs = QTabWidget()
        self.editors: dict[str, QPlainTextEdit] = {}
        self._set_tabs(["dftb_in.hsd", "geometry.gen", "submit.sh", "spec.json", "README.txt"])
        self.prov = QLabel(""); self.prov.setObjectName("hint"); self.prov.setWordWrap(True)
        self.btn_prov = QPushButton(L("SHA-256 を表示", "Show SHA-256")); self.btn_prov.setObjectName("link"); self.btn_prov.setCheckable(True)
        self.prov_text = QPlainTextEdit(); self.prov_text.setReadOnly(True); self.prov_text.setMaximumHeight(96)
        self.prov_text.setStyleSheet("font-family: monospace; font-size: 9pt;"); self.prov_text.hide()
        self.prov_row = QWidget(); self.prov_row.setObjectName("rowbox")
        prl = QHBoxLayout(self.prov_row); prl.setContentsMargins(0, 0, 0, 0); prl.setSpacing(6)
        prl.addWidget(self.prov, 1); prl.addWidget(self.btn_prov, 0, Qt.AlignmentFlag.AlignTop)
        self.prov_row.hide()
        self.btn_prov.toggled.connect(self._toggle_prov)
        lay = QVBoxLayout(self); lay.setContentsMargins(12, 12, 12, 12); lay.setSpacing(10)
        lay.addWidget(self.status); lay.addWidget(self.detail); lay.addWidget(self.origin_row); lay.addWidget(self.tabs, 1)
        lay.addWidget(self.empty, 1)
        lay.addWidget(self.prov_row); lay.addWidget(self.prov_text)

    def _toggle_prov(self, on: bool) -> None:
        self.prov_text.setVisible(on and bool(self.prov_text.toPlainText()))
        self.btn_prov.setText(L("SHA-256 を隠す", "Hide SHA-256") if on else L("SHA-256 を表示", "Show SHA-256"))

    def set_origin(self, text: str) -> None:
        self.origin.setText(text); self.origin_row.setVisible(bool(text))

    def set_provenance(self, prov: dict | None) -> None:
        from adit.web.codefields import provenance_summary
        head, lines = provenance_summary(prov)
        self.prov.setText(head); self.prov.setToolTip("\n".join(lines))
        self.prov_text.setPlainText("\n".join(lines))
        self.btn_prov.setVisible(bool(lines)); self.prov_row.setVisible(bool(head))
        self._toggle_prov(self.btn_prov.isChecked())

    def _set_tabs(self, names: list[str]) -> None:
        if list(self.editors) == names:
            return
        cur = self.tabs.tabText(self.tabs.currentIndex()) if self.tabs.count() else ""
        while self.tabs.count():
            self.tabs.removeTab(0)
        self.editors = {}
        for name in names:
            ed = QPlainTextEdit(); ed.setReadOnly(True); ed.setStyleSheet("font-family: monospace;")
            self.editors[name] = ed; self.tabs.addTab(ed, name)
        if cur in self.editors:
            self.tabs.setCurrentWidget(self.editors[cur])

    def show_files(self, files: ProjectFiles) -> None:
        names = [n for n in ORDER if n in files.texts] + [n for n in files.texts if n not in ORDER]
        self._set_tabs(names)
        self.tabs.show(); self.empty.hide()
        for name, ed in self.editors.items():
            ed.setPlainText(files.texts.get(name, ""))
        copies = ", ".join(sorted(files.copies))
        self.status.setText(tr("生成できます"))
        self.status.setObjectName("status_ok"); self.status.style().polish(self.status)
        self.detail.setText(L(f"コピーされるファイル: {copies}", f"Files to copy: {copies}"))
        import json
        try:
            prov = json.loads(files.texts.get("spec.json") or "{}").get("provenance")
        except ValueError:
            prov = None
        self.set_provenance(prov)

    def show_errors(self, message: str, *, can_jump: bool = False) -> None:
        self.set_provenance(None)
        error_tab = L("エラー", "Errors")
        self._set_tabs([error_tab])
        self.editors[error_tab].setPlainText(message)
        first, _, rest = message.partition("\n")
        self.status.setText(first)
        self.status.setObjectName("status_ng"); self.status.style().polish(self.status)
        self.detail.setText(rest)
        items = [x.strip() for x in rest.splitlines() if x.strip()]
        line = first if not items else items[0] + (L(f" (ほか {len(items) - 1} 件)", f" (+{len(items) - 1} more)") if len(items) > 1 else "")
        self.empty.set_texts(line=line, button=L("欄へ移動", "Go to the field") if can_jump else "")
        self.tabs.hide(); self.empty.show()
