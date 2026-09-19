
from __future__ import annotations

import os
import shutil
import shlex
from pathlib import Path

from PySide6.QtCore import QEvent, QProcess, Qt, QTimer
from PySide6.QtGui import QAction, QActionGroup, QIcon
from PySide6.QtWidgets import (QBoxLayout, QFileDialog, QHBoxLayout, QMainWindow, QMessageBox, QPushButton, QScrollArea,
                               QLabel, QMenuBar, QSizePolicy, QSplitter, QStackedWidget, QTabWidget, QToolBar,
                               QVBoxLayout,
                               QWidget)

from adit.config import Config, ConfigError, load_config
from adit.project import ProjectError, ProjectFiles, build_project, load_project, write_project
from adit.gui import icons
from adit.gui.i18n import tr
from adit.lang import L
from adit.results import summarize_run
from adit.spec import CalculationSpec
from adit.gui.panels.analysis_panel import AnalysisPanel
from adit.gui.panels.kpoints_panel import KPointsPanel
from adit.gui.panels.method_panel import CODES, MethodPanel
from adit.gui.panels.preview_panel import PreviewPanel
from adit.gui.panels.runtime_panel import RuntimePanel
from adit.gui.panels.structure_view_panel import StructureViewPanel
from adit.gui.panels.structure_panel import StructurePanel
from adit.gui.panels.task_panel import TaskPanel
from adit.gui.style import GROUP_SPACING, PANEL_MARGIN
from adit.gui.ribbon import ModeBar
from adit.gui.widgets import limit_combo_popups
from adit.web.codefields import PrepOrigin, template_marks


class MainWindow(QMainWindow):
    def __init__(self, cfg: Config, cfg_path: Path):
        super().__init__()
        self.cfg, self.cfg_path = cfg, cfg_path
        self.setWindowTitle("ADIT - Atomistic Design and Interpretation Toolkit")
        available = self.screen().availableGeometry() if self.screen() is not None else None
        if available is not None:
            self.resize(min(1800, int(available.width() * 0.95)), min(1000, int(available.height() * 0.95)))
            self.move(available.left() + max(0, (available.width() - self.width()) // 2),
                      available.top() + max(0, (available.height() - self.height()) // 2))
        else:
            self.resize(1800, 1000)
        self.last_written: Path | None = None
        self.last_written_kind: str | None = None
        self.last_written_exe: str | None = None
        self._proc: QProcess | None = None
        self.origin = PrepOrigin()
        self._marked: set = set()

        self.structure = StructurePanel()
        self.method = MethodPanel(cfg)
        self.kpoints = KPointsPanel()
        self.task = TaskPanel()
        self.runtime = RuntimePanel(cfg)
        self.preview = PreviewPanel()
        self.analysis = AnalysisPanel()
        self.structure_view = StructureViewPanel()
        self.files_pane = self.preview
        self.right_tabs = QTabWidget(); self.right_tabs.setDocumentMode(True)
        self.right_tabs.addTab(self.structure_view, icons.icon("tab_structure"), "構造")
        self.right_tabs.addTab(self.files_pane, icons.icon("tab_files"), "生成ファイル")
        self.right_tabs.setCurrentWidget(self.structure_view)
        from adit.gui.style import current_theme
        self.structure_view.set_dark(current_theme() == "dark")
        self._dark = current_theme() == "dark"

        limit_combo_popups()
        self.act_open = QAction(icons.icon("open", 32), "計算設定 (spec.json) を開く…", self); self.act_open.setShortcut("Ctrl+O"); self.act_open.setIconText("開く")
        self.act_reload = QAction(icons.icon("reload", 32), "環境設定を再読み込み", self); self.act_reload.setIconText("再読み込み")
        self.act_back = QAction(icons.icon("undo", 32), "元に戻す", self); self.act_back.setToolTip("1 つ前の設定に戻す"); self.act_back.setShortcut("Ctrl+Z")
        self.act_forward = QAction(icons.icon("redo", 32), "やり直す", self); self.act_forward.setToolTip("やり直す"); self.act_forward.setShortcut("Ctrl+Shift+Z")
        self.act_settings = QAction(icons.icon("settings", 32), "環境設定…", self); self.act_settings.setToolTip("環境設定ファイル (cluster.toml) を編集"); self.act_settings.setIconText("環境設定")
        self._history: list[str] = []; self._hist_pos = -1; self._applying = False
        tb = QToolBar("操作"); tb.setMovable(False); tb.setObjectName("ribbon_actions"); self.action_bar = tb
        self.act_lang = QAction("English / 日本語", self)
        self.act_lang.setToolTip("表示言語を切り替えます (環境設定の language に保存。再起動後に反映)")
        self.btn_generate = QPushButton(icons.icon("generate", 16), "生成"); self.btn_generate.setObjectName("primary"); self.btn_generate.setDefault(True)
        self.run_hint = QLabel(""); self.run_hint.setObjectName("hint")
        home = str(Path.home())
        shown = str(cfg_path).replace(home, "~", 1) if str(cfg_path).startswith(home) else str(cfg_path)
        self.cfg_label = QLabel(L(f"環境設定: {shown}", f"preferences: {shown}")); self.cfg_label.setToolTip(str(cfg_path)); self.cfg_label.setObjectName("hint")
        self.gen_hint = QPushButton(""); self.gen_hint.setObjectName("gen_hint"); self.gen_hint.setFlat(True)
        self.gen_hint.setMaximumWidth(560); self.gen_hint.setCursor(Qt.CursorShape.PointingHandCursor)
        self.gen_hint.clicked.connect(self._on_gen_hint_clicked)
        spacer = QWidget(); spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        tb.addWidget(spacer)
        self.gen_hint_action = tb.addWidget(self.gen_hint)
        self.gen_hint_action.setVisible(False)
        tb.addWidget(self.btn_generate)
        self.statusBar().addWidget(self.run_hint)
        self.legend = QLabel(L("青字は必須です。最初から入っている値は、説明に出典がない限り ADIT が置いた値です。ラベルにカーソルを合わせると説明が表示されます",
                               "Blue fields are required. A pre-filled value was supplied by ADIT unless its explanation names another source. Hover a label for details.")); self.legend.setObjectName("hint")
        self.statusBar().addPermanentWidget(self.cfg_label)
        self.statusBar().addPermanentWidget(self.legend)
        self.act_back.triggered.connect(self.go_back); self.act_forward.triggered.connect(self.go_forward)
        self.act_settings.triggered.connect(self.open_settings)
        self._update_history_buttons()
        self._build_ribbon()

        left = QWidget(); self._left = left
        ll = QBoxLayout(QBoxLayout.Direction.LeftToRight, left); self._left_layout = ll
        ll.setSpacing(GROUP_SPACING); ll.setContentsMargins(PANEL_MARGIN, 8, 8, PANEL_MARGIN)
        self._col1w, self._col2w = QWidget(), QWidget()
        col1 = QVBoxLayout(self._col1w); col2 = QVBoxLayout(self._col2w)
        for cw, col in ((self._col1w, col1), (self._col2w, col2)):
            cw.setObjectName("rowbox")
            col.setContentsMargins(0, 0, 0, 0); col.setSpacing(GROUP_SPACING)
        for w in (self.task, self.kpoints, self.runtime):
            col1.addWidget(w)
        col2.addWidget(self.method)
        col1.addStretch(); col2.addStretch()
        ll.addWidget(self._col1w, 1); ll.addWidget(self._col2w, 1)
        scroll = QScrollArea(); scroll.setWidget(left); scroll.setWidgetResizable(True); self._left_scroll = scroll
        scroll.viewport().installEventFilter(self); left.installEventFilter(self)

        structure_page = QWidget()
        sl = QVBoxLayout(structure_page); sl.setContentsMargins(PANEL_MARGIN, 8, 8, PANEL_MARGIN); sl.setSpacing(GROUP_SPACING)
        sl.addWidget(self.structure); sl.addStretch()
        structure_scroll = QScrollArea(); structure_scroll.setWidget(structure_page); structure_scroll.setWidgetResizable(True)

        analysis_scroll = QScrollArea(); analysis_scroll.setWidget(self.analysis); analysis_scroll.setWidgetResizable(True)

        from adit.gui.panels.workspace_panel import WorkspacePanel

        start = Path(self.runtime.output_dir() or "~").expanduser().parent
        self.workspace = WorkspacePanel(root=start if start.is_dir() else Path.home(), dark=self._dark)
        self.main_stack = QStackedWidget()
        for w in (structure_scroll, scroll, analysis_scroll, self.workspace):
            self.main_stack.addWidget(w)
        self.mode_bar = ModeBar([("tab_structure", "構造"), ("settings", "計算条件"), ("tab_analysis", "解析"),
                                 ("tab_workspace", "ワークスペース")])
        self.mode_bar.changed.connect(self.set_mode)

        split = QSplitter(); split.addWidget(self.main_stack); split.addWidget(self.right_tabs); split.setSizes([1200, 580])
        split.setHandleWidth(GROUP_SPACING)
        root = QWidget(); rl = QVBoxLayout(root); rl.setContentsMargins(PANEL_MARGIN, 2, PANEL_MARGIN, 0); rl.setSpacing(4)
        rl.addWidget(self.mode_bar); rl.addWidget(split, 1); rl.addWidget(self.action_bar)
        self.setCentralWidget(root)

        self._timer = QTimer(self); self._timer.setSingleShot(True); self._timer.setInterval(250)
        self._timer.timeout.connect(self.refresh_preview)
        for p in (self.structure, self.method, self.kpoints, self.task, self.runtime):
            p.changed.connect(self._timer.start)
        self.structure.changed.connect(self._on_context)
        self.runtime.changed.connect(self._on_context)
        self._on_context()
        self.method.dftb.root_chosen.connect(lambda p: self._set_root("sk_root", p))
        self.method.espresso.root_chosen.connect(lambda p: self._set_root("pseudo_root", p))
        self.method.code.currentIndexChanged.connect(self._on_context)
        self.method.gromacs.use_as_structure.connect(self._use_structure_file)
        self.act_open.triggered.connect(self.open_spec)
        self.preview.btn_clear_origin.clicked.connect(self.clear_origin)
        self.preview.fix_requested.connect(self._on_gen_hint_clicked)
        self.act_reload.triggered.connect(self.reload_config)
        self.btn_generate.clicked.connect(self.generate)
        self.act_lang.triggered.connect(self._toggle_language)
        self.refresh_preview()
        limit_combo_popups(self)

    def _fit_gen_hint(self) -> None:
        width = self.width()
        self.gen_hint.setMaximumWidth(560 if width >= 1600 else 380 if width >= 1400 else 260)

    def resizeEvent(self, event):  # noqa: N802 
        super().resizeEvent(event)
        self._fit_gen_hint()
        text = getattr(self, "_gen_hint_full", "")
        if text and self.gen_hint_action.isVisible():
            fm = self.gen_hint.fontMetrics()
            self.gen_hint.setText(fm.elidedText(text, Qt.TextElideMode.ElideRight, self.gen_hint.maximumWidth() - 16))

    def _find_field(self, location: str):
        from PySide6.QtWidgets import QFormLayout, QLabel

        from adit.validate_types import place_key

        mode = self.MODE_STRUCTURE if location.startswith("structure") else self.MODE_SETTINGS
        key = place_key(location)
        page = self.main_stack.widget(mode)
        labels = [w for w in page.findChildren(QLabel) if w.property("adit_key")]
        exact = [w for w in labels if w.property("adit_key") == key]
        starts = [w for w in labels if str(w.property("adit_key")).startswith(key)]
        holds = [w for w in labels if key in str(w.property("adit_key"))]
        found = (exact or starts or holds)
        if not found:
            return mode, None, None
        return mode, found[0], self._field_of(found[0])

    @staticmethod
    def _field_of(label):
        from PySide6.QtWidgets import QFormLayout

        parent = label.parentWidget()
        if parent is None:
            return None
        forms = [lay for lay in (parent.layout(), *parent.findChildren(QFormLayout)) if isinstance(lay, QFormLayout)]
        for layout in forms:
            for row in range(layout.rowCount()):
                item = layout.itemAt(row, QFormLayout.ItemRole.LabelRole)
                if item is not None and item.widget() is label:
                    field_item = layout.itemAt(row, QFormLayout.ItemRole.FieldRole)
                    if field_item is None:
                        return None
                    return field_item.widget() or (field_item.layout() and field_item.layout().itemAt(0).widget())
        return None

    def _reveal(self, mode: int, label, field) -> bool:
        from PySide6.QtWidgets import QScrollArea

        self.set_mode(mode)
        target = field or label
        if target is None:
            return False
        parent = target.parentWidget()
        while parent is not None and not isinstance(parent, QScrollArea):
            parent = parent.parentWidget()
        if isinstance(parent, QScrollArea):
            parent.ensureWidgetVisible(target, 60, 60)
        if field is not None:
            field.setFocus(Qt.FocusReason.OtherFocusReason)
        return True

    def _mode_of(self, w) -> int:
        for i in range(self.main_stack.count()):
            if self.main_stack.widget(i).isAncestorOf(w):
                return i
        return -1

    def field_entries(self) -> list[tuple[object, int, str]]:
        # (label or checkbox, mode, group title) for every field that is shown right now
        from PySide6.QtWidgets import QCheckBox, QGroupBox

        out = []
        for mode in range(self.main_stack.count()):
            page = self.main_stack.widget(mode)
            for w in [*page.findChildren(QLabel), *page.findChildren(QCheckBox)]:
                if not w.property("adit_key") or not w.isVisibleTo(page):
                    continue
                box = w.parentWidget()
                while box is not None and not isinstance(box, QGroupBox):
                    box = box.parentWidget()
                out.append((w, mode, box.title() if box is not None else ""))
        return out

    def jump_to_label(self, label) -> bool:
        from PySide6.QtWidgets import QCheckBox

        mode = self._mode_of(label)
        if mode < 0:
            return False
        field = label if isinstance(label, QCheckBox) else self._field_of(label)
        return self._reveal(mode, label, field)

    def clear_error_marks(self) -> None:
        for w in getattr(self, "_error_marked", []):
            w.setProperty("adit_error", False)
            w.style().unpolish(w); w.style().polish(w)
        self._error_marked = []

    def focus_error(self, index: int = 0) -> bool:
        locations = getattr(self, "_error_locations", [])
        if not locations:
            return False
        index %= len(locations)
        self._error_index = index
        mode, label, field = self._find_field(locations[index])
        self.clear_error_marks()
        marks = [w for w in (label, field) if w is not None]
        for w in marks:
            w.setProperty("adit_error", True)
            w.style().unpolish(w); w.style().polish(w)
        self._error_marked = marks
        return self._reveal(mode, label, field)

    def _on_gen_hint_clicked(self) -> None:
        if not self.focus_error(getattr(self, "_error_index", -1) + 1):
            self.right_tabs.setCurrentWidget(self.files_pane)

    MODE_STRUCTURE, MODE_SETTINGS, MODE_ANALYSIS, MODE_WORKSPACE = 0, 1, 2, 3

    def set_mode(self, index: int) -> None:
        if not 0 <= index < self.main_stack.count():
            return
        self.main_stack.setCurrentIndex(index)
        self.right_tabs.setVisible(index != self.MODE_WORKSPACE)   # the workspace uses the whole width
        self.mode_bar.set_current(index)
        if hasattr(self, "act_save"):
            # In the workspace Ctrl+S saves the edited file; two shortcuts on one key would cancel each other.
            self.act_save.setShortcut("" if index == self.MODE_WORKSPACE else "Ctrl+S")
        if hasattr(self, "_mode_actions"):
            self._mode_actions[index].setChecked(True)

    def mode(self) -> int:
        return self.main_stack.currentIndex()

    def left_columns(self) -> int:
        return 2 if self._left_layout.direction() == QBoxLayout.Direction.LeftToRight else 1

    def _fit_left_columns(self) -> None:
        m = self._left_layout.contentsMargins()
        need = m.left() + m.right() + self._col1w.minimumSizeHint().width() + self._left_layout.spacing() + self._col2w.minimumSizeHint().width()
        want = QBoxLayout.Direction.LeftToRight if self._left_scroll.viewport().width() >= need else QBoxLayout.Direction.TopToBottom
        if self._left_layout.direction() != want:
            self._left_layout.setDirection(want)
            self._left_layout.invalidate(); self._left_layout.activate(); self._left.updateGeometry()

    def eventFilter(self, obj, event) -> bool:  # noqa: N802 
        scroll = getattr(self, "_left_scroll", None)
        if scroll is not None and event.type() in (QEvent.Type.Resize, QEvent.Type.LayoutRequest) and obj in (scroll.viewport(), self._left):
            self._fit_left_columns()
        return super().eventFilter(obj, event)

    def _on_context(self, *_) -> None:
        st = self.structure.structure()
        periodic = bool(st and st.periodic)
        self.structure_view.set_structure(st, self.structure.error())
        self.kpoints.set_periodic(periodic)
        self.kpoints.setVisible(periodic and self._uses_kpoints())
        self.task.set_periodic(periodic)
        self.method.set_periodic(periodic)
        elements = []
        if st is not None:
            for s in st.atoms.symbols:
                if s not in elements:
                    elements.append(s)
        self.method.set_context(elements, self.runtime.runtime().profile, self.cfg)
        self._fit_left_columns()

    def _uses_kpoints(self) -> bool:
        from adit.codes import GENERATORS
        gen = GENERATORS.get(self.method.current_code())
        return gen is None or gen.uses_kpoints

    def _use_structure_file(self, path: str) -> None:
        self.structure.set_source("file")
        self.structure._file_restored = None
        self.structure.file.setText(path)
        self.structure._rebuild()

    def current_spec(self) -> CalculationSpec:
        st = self.structure.structure()
        if st is None:
            raise ProjectError(L("構造: ", "Structure: ") + self.structure.error())
        kp = self.kpoints.kpoints() if st.periodic and self._uses_kpoints() else None
        spec = CalculationSpec(structure=st, method=self.method.method(), kpoints=kp, task=self.task.task(), runtime=self.runtime.runtime())
        return self.origin.apply(spec)

    def build(self) -> ProjectFiles:
        return build_project(self.current_spec(), self.cfg, output_dir=self.runtime.output_dir() or None)

    def _show_origin(self) -> None:
        try:
            spec = self.current_spec()
        except Exception:
            spec = None
        self.preview.set_origin(self.origin.describe(spec) if self.origin.active else "")
        self._apply_template_marks(template_marks(spec) if spec is not None else {})

    def _apply_template_marks(self, marks: dict[str, str]) -> None:
        import html
        from PySide6.QtWidgets import QCheckBox
        for w in [*self._left.findChildren(QLabel), *self._left.findChildren(QCheckBox)]:
            key = w.property("adit_key")
            if not key:
                continue
            name = marks.get(key)
            if name is None:
                if w in self._marked:
                    w.setText(tr(key)); self._marked.discard(w)
                continue
            note = L(f"雛形 {name} の値", f"value from template {name}")
            if isinstance(w, QLabel):
                w.setText(f"{html.escape(tr(key))}<br><span style='font-size:9pt; color:#8e8e93;'>{html.escape(note)}</span>")
            else:
                w.setText(f"{tr(key)}  ({note})")
            self._marked.add(w)

    def clear_origin(self) -> None:
        self.origin = PrepOrigin()
        self.refresh_preview()

    def refresh_preview(self) -> None:
        try:
            if not self.runtime.output_dir().strip():
                raise ProjectError(L("出力ディレクトリ: 生成したファイルを置く場所を指定してください", "Output directory: choose where to put the generated files"))
            self.preview.show_files(self.build())
            self.btn_generate.setEnabled(True); self.act_generate.setEnabled(True)
            self.gen_hint_action.setVisible(False)
            self._error_locations = []; self._error_index = -1
            self.clear_error_marks()
        except (ProjectError, ConfigError, ValueError, TypeError) as ex:
            from pydantic import ValidationError as PydanticError
            if isinstance(ex, PydanticError):
                from adit.validate_types import friendly_pydantic
                ex = ProjectError(friendly_pydantic(ex))
            self.btn_generate.setEnabled(False); self.act_generate.setEnabled(False)
            self._show_gen_hint(ex)
            self.preview.show_errors(str(ex), can_jump=bool(self._error_locations))
        self._show_origin()
        self._show_doc_values()
        self._on_structure_or_errors_changed()
        self._record_history()
        if self.last_written is not None:
            self.last_written = None
            self.run_hint.clear()

    def _show_doc_values(self) -> None:
        from adit.web.prep23 import doc_lines

        try:
            spec = self.current_spec()
        except Exception:
            self.method.set_doc_values({}); self._doc_key = None
            return
        key = (spec.method.model_dump_json(), tuple(spec.elements), self.cfg.sk_root, self.cfg.pseudo_root,
               self.cfg.cp2k_data, spec.runtime.profile)
        if key == getattr(self, "_doc_key", None):
            return
        self._doc_key = key
        try:
            groups = doc_lines(spec, self.cfg)
        except (OSError, ValueError):
            groups = {}
        self.method.set_doc_values(groups)

    def _show_gen_hint(self, ex: Exception) -> None:
        errs = getattr(ex, "errors", None) or (ex.args[1] if len(ex.args) > 1 and isinstance(ex.args[1], list) else None)
        first = str(errs[0]) if errs else (str(ex).splitlines() or [""])[0]
        more = L(f" (ほか {len(errs) - 1} 件)", f" (+{len(errs) - 1} more)") if errs and len(errs) > 1 else ""
        full = L("生成できません: ", "cannot generate: ") + first + more
        fm = self.gen_hint.fontMetrics()
        self._fit_gen_hint()
        self._gen_hint_full = full
        self.gen_hint.setText(fm.elidedText(full, Qt.TextElideMode.ElideRight, self.gen_hint.maximumWidth() - 16))
        self._error_locations = [e.location for e in errs] if errs else []
        self._error_index = -1
        count = L(f" (全 {len(errs)} 件。押すたびに次へ)", f" ({len(errs)} in total; click again for the next one)") if errs and len(errs) > 1 else ""
        self.gen_hint.setToolTip(L("押すと、その欄へ移動して赤い枠で囲みます", "click to jump to the field and outline it in red") + count + "\n\n" + str(ex))
        self.gen_hint_action.setVisible(True)

    def _on_structure_or_errors_changed(self) -> None:
        pass

    def _set_root(self, key: str, path: str) -> None:
        setattr(self.cfg, key, path); self._save_cfg()
        self.reload_config()
        self.statusBar().showMessage(L(f"{key} を {path} にしました (環境設定に保存しました)", f"{key} set to {path} (saved in the preferences)"))


    def _toggle_language(self) -> None:
        self.cfg.language = "en" if self.cfg.language != "en" else "ja"
        try:
            from adit.config import save_config
            save_config(self.cfg, self.cfg_path)
            self.statusBar().showMessage(tr("保存しました。表示言語は次回の起動から: ") + self.cfg.language)
        except OSError as ex:
            self.statusBar().showMessage(L(f"設定を保存できません: {ex}", f"cannot save the settings: {ex}"))



    @staticmethod
    def _executable_of(run_command: str) -> str:
        tail = run_command.split("&&")[-1].strip()
        words = tail.split()
        if not words:
            return ""
        if words[0] in ("mpirun", "mpiexec", "srun"):
            rest = [w for w in words[1:] if not w.startswith("-") and not w.isdigit()]
            return rest[0] if rest else ""
        return words[0]


    def generate(self) -> None:
        if not self._flush():
            return
        out = Path(self.runtime.output_dir()).expanduser()
        try:
            spec = self.current_spec()
            overwrite = False
            if out.exists() and any(out.iterdir()):
                ans = QMessageBox.question(self, tr("上書きの確認"), L(f"{out} は空ではありません。中のファイルを上書きしますか?", f"{out} is not empty. Overwrite the files inside?"))
                if ans != QMessageBox.StandardButton.Yes:
                    return
                overwrite = True
            written = write_project(spec, self.cfg, out, overwrite=overwrite)
        except (ProjectError, ConfigError, ValueError, OSError) as ex:
            QMessageBox.critical(self, tr("生成できません"), str(ex))
            return
        self.last_written = out
        self.analysis.set_run_dir(out, spec.elements)
        self.workspace.set_root(out)          # point the tree and the terminal at the directory just written
        self.workspace.terminal.send(f"cd {shlex.quote(str(out))}\r")
        self.last_written_kind = self.cfg.profiles[spec.runtime.profile].kind
        from adit.codes import GENERATORS
        self.last_written_exe = self._executable_of(GENERATORS[spec.method.code].run_command(spec, self.cfg.profiles[spec.runtime.profile]))
        t = spec.task
        self.last_written_steps = t.md.steps if t.type == "molecular_dynamics" else (t.max_steps if t.type == "geometry_optimization" and t.max_steps > 0 else None)
        kind = self.cfg.profiles[spec.runtime.profile].kind
        step = (L("ワークスペースのターミナルで  bash submit.sh  を実行します。",
                  "In the workspace terminal, run  bash submit.sh")
                if kind == "direct" else
                L("transfer_and_submit.sh のコマンドで、クラスタへ送って投入します。",
                  "Use the commands in transfer_and_submit.sh to send it to the cluster and submit it."))
        self.run_hint.setText(L(f"次: {step}", f"Next: {step}"))
        self.statusBar().showMessage(L(f"{len(written)} ファイルを {out} に書きました", f"wrote {len(written)} files to {out}"))
        QMessageBox.information(self, tr("生成しました"), f"{out}\n\n{step}")

    def scan_dialog(self):
        from adit.gui.scan_dialog import ScanDialog
        st = self.structure.structure()
        base = self.runtime.output_dir()
        return ScanDialog(self.method.current_code(), bool(st and st.periodic), (base.rstrip("/") + "_scan") if base else "",
                          self.current_spec, self.cfg, self)

    def open_scan_dialog(self) -> None:
        if self._timer.isActive():
            self._timer.stop(); self.refresh_preview()
        dlg = self.scan_dialog()
        if dlg.exec() and dlg.out_dir is not None:
            self.scan_written(dlg.out_dir)

    def scan_written(self, out: Path) -> None:
        self.analysis.set_run_dir(out)
        self.statusBar().showMessage(L(f"値ごとの入力を {out} に作りました", f"wrote the scan inputs to {out}"))

    def _flush(self) -> bool:
        # The preview may still be pending (250 ms debounce): settle it before acting on the button state.
        if self._timer.isActive():
            self._timer.stop(); self.refresh_preview()
        return self.btn_generate.isEnabled()

    def continue_dialog(self):
        from adit.gui.prep_dialogs import ContinueDialog
        return ContinueDialog(self.current_spec, str(self.last_written or ""), self)

    def open_continue_dialog(self) -> None:
        self._flush()
        dlg = self.continue_dialog()
        if dlg.exec() and dlg.spec is not None:
            self.apply_continuation(dlg.spec)
            QMessageBox.information(self, L("前の計算の続きを作りました", "Continuation created"), dlg.summary() + "\n\n" + L(
                "出力ディレクトリは、前の計算とは別の場所にしてください (同じ場所に書くと前の出力を上書きします)。",
                "Use an output directory other than the previous run's (writing there would overwrite its output)."))

    def apply_continuation(self, spec: CalculationSpec) -> None:
        self.apply_spec(spec)
        prev = Path((spec.meta.continued_from or {}).get("dir", ""))
        if prev.name and Path(self.runtime.output_dir()).expanduser().resolve() == prev.resolve():
            self.runtime.outdir.setText(str(prev) + "_cont")
        self.statusBar().showMessage(L(f"前の計算 {prev} の続きを読み込みました", f"loaded a continuation of {prev}"))

    def stages_dialog(self):
        from adit.gui.prep_dialogs import StagesDialog
        base = self.runtime.output_dir()
        return StagesDialog((base.rstrip("/") + "_stages") if base else "", self.current_spec, self.cfg, self)

    def open_stages_dialog(self) -> None:
        self._flush()
        dlg = self.stages_dialog()
        if dlg.exec() and dlg.out_dir is not None:
            if dlg.dirs:
                self.analysis.set_run_dir(dlg.dirs[0])
            self.statusBar().showMessage(L(f"段階に分けた入力を {dlg.out_dir} に作りました", f"wrote the staged inputs to {dlg.out_dir}"))

    def batch_dialog(self, kind: str):
        from adit.gui.prep23_dialogs import DIALOGS
        base = self.runtime.output_dir()
        suffix = {"compare": "_set", "conformers": "_conf", "neb": "_neb", "phonons": "_phonon", "elastic": "_elastic", "ts": "_ts"}[kind]
        default = (base.rstrip("/") + suffix) if base else ""
        if kind == "ts":
            return DIALOGS[kind](self.current_spec, self.cfg, default, self.method.current_code(), self)
        return DIALOGS[kind](self.current_spec, self.cfg, default, self)

    def open_batch_dialog(self, kind: str) -> None:
        self._flush()
        dlg = self.batch_dialog(kind)
        if dlg.exec() and dlg.result is not None:
            self.batch_written(dlg.result)

    def batch_written(self, res) -> None:
        if res.compare_base and res.dirs:
            self.analysis.set_run_dir(res.dirs[0])
        elif res.analysis_dir:
            self.analysis.set_run_dir(Path(res.analysis_dir))
        self.statusBar().showMessage(L(f"{res.out} に {len(res.dirs)} 個のディレクトリを作りました",
                                       f"wrote {len(res.dirs)} directories in {res.out}"))

    def template_load_dialog(self):
        from adit.gui.prep_dialogs import TemplateLoadDialog
        return TemplateLoadDialog(self.cfg, self)

    def open_template_load_dialog(self) -> None:
        self._flush()
        dlg = self.template_load_dialog()
        if dlg.exec() and dlg.chosen:
            self.load_template(dlg.chosen)

    def load_template(self, name: str) -> bool:
        from adit.templates import TemplateError, load_template
        st = self.structure.structure()
        if st is None:
            QMessageBox.critical(self, L("読み込めません", "Cannot load"), L("先に構造を作ってください (雛形は構造を持たないため)。",
                                                                          "Make a structure first (templates have no structure).")); return False
        try:
            spec = load_template(name, st, self.cfg)
        except TemplateError as ex:
            QMessageBox.critical(self, L("読み込めません", "Cannot load"), str(ex)); return False
        keep = self.origin
        self.apply_spec(spec)
        if keep.continued_from and keep.same_structure(st):
            self.origin.continued_from, self.origin.handoff, self.origin.velocities = keep.continued_from, keep.handoff, keep.velocities
            self.origin.code, self.origin.source_ref, self.origin.symbols = keep.code, keep.source_ref, keep.symbols
            self.refresh_preview()
        self.statusBar().showMessage(L(f"雛形 {name} を読み込みました (構造はそのまま)", f"loaded the template {name} (structure kept)"))
        return True

    def template_save_dialog(self):
        from adit.gui.prep_dialogs import TemplateSaveDialog
        return TemplateSaveDialog(self.cfg, self.current_spec, self)

    def open_template_save_dialog(self) -> None:
        self._flush()
        dlg = self.template_save_dialog()
        if dlg.exec() and dlg.path is not None:
            self.statusBar().showMessage(L(f"雛形を保存しました: {dlg.path}", f"saved the template: {dlg.path}"))

    def open_convert_dialog(self) -> None:
        self._flush()
        from adit.gui.convert_dialog import ConvertDialog
        ConvertDialog(self.current_spec, self.cfg, self).exec()

    def open_input_audit_dialog(self) -> None:
        from adit.gui.input_audit_dialog import InputAuditDialog
        InputAuditDialog(self).exec()

    def menuBar(self) -> QMenuBar:  # noqa: N802 
        if self._menubar_stub is None:
            self._menubar_stub = QMenuBar(self); self._menubar_stub.hide()
        return self._menubar_stub

    def _choice_actions(self, items, current: str, slot, icon_prefix: str) -> list[tuple[str, QAction]]:
        group = QActionGroup(self); group.setExclusive(True)
        made = []
        for code, text, short in items:
            a = QAction(icons.icon(f"{icon_prefix}{code}", 32) if icon_prefix else QIcon(), text, self)
            a.setIconText(short); a.setCheckable(True); a.setChecked(current == code); group.addAction(a)
            a.triggered.connect(lambda _c=False, code=code: slot(code)); made.append((code, a))
        return made

    def _build_ribbon(self) -> None:
        from adit.gui import titlebar
        from adit.gui.ribbon import QuickAccessBar, Ribbon
        self._menubar_stub = None
        self.act_save = QAction(icons.icon("save", 32), "計算設定 (spec.json) を保存…", self); self.act_save.setShortcut("Ctrl+S"); self.act_save.setIconText("保存"); self.act_save.triggered.connect(self.save_spec)
        self.act_quit = QAction(icons.icon("quit", 32), "終了", self); self.act_quit.setShortcut("Ctrl+Q"); self.act_quit.triggered.connect(self.close)
        self.act_draw = QAction(icons.icon("draw", 32), "Draw", self); self.act_draw.setShortcut("Ctrl+D"); self.act_draw.triggered.connect(self.structure._draw_smiles)
        self.act_draw.setToolTip("Draw: 分子を描いて SMILES にします (RDKit が要ります)"); self.act_draw.setEnabled(self.structure.draw_button.isEnabled())
        self.act_load_structure = QAction(icons.icon("file", 32), "構造ファイルを読み込む…", self); self.act_load_structure.setIconText("構造ファイル"); self.act_load_structure.triggered.connect(self._insert_structure_file)
        self.act_add_component = QAction(icons.icon("mixture", 32), "溶液に成分を追加", self); self.act_add_component.setIconText("成分を追加"); self.act_add_component.triggered.connect(self._insert_component)
        self._tab_actions = []
        tab_group = QActionGroup(self); tab_group.setExclusive(True)
        for i in range(self.right_tabs.count()):
            a = QAction(icons.icon(("tab_structure", "tab_files")[i], 32), self.right_tabs.tabText(i), self)
            a.setCheckable(True); a.setChecked(i == self.right_tabs.currentIndex()); tab_group.addAction(a)
            a.triggered.connect(lambda _c=False, i=i: self.right_tabs.setCurrentIndex(i)); self._tab_actions.append(a)
        self.right_tabs.currentChanged.connect(lambda i: 0 <= i < len(self._tab_actions) and self._tab_actions[i].setChecked(True))
        self._mode_actions = []
        mode_group = QActionGroup(self); mode_group.setExclusive(True)
        for i, (icon_name, text) in enumerate((("tab_structure", "構造"), ("settings", "計算条件"), ("tab_analysis", "解析"),
                                               ("tab_workspace", "ワークスペース"))):
            a = QAction(icons.icon(icon_name, 32), text, self)
            a.setCheckable(True); a.setChecked(i == 0); mode_group.addAction(a)
            a.triggered.connect(lambda _c=False, i=i: self.set_mode(i)); self._mode_actions.append(a)
        self.act_generate = QAction(icons.icon("generate", 32), "生成", self); self.act_generate.setShortcut("Ctrl+G"); self.act_generate.triggered.connect(self.generate)
        self.act_scan = QAction(icons.icon("scan", 32), L("1 つの条件を変えて一括生成…", "Generate a parameter scan…"), self)
        self.act_scan.setIconText(L("一括生成", "Parameter scan"))
        self.act_scan.setToolTip(L("カットオフや k 点など 1 つの値だけを変えた入力を、値ごとのディレクトリにまとめて作ります (収束の確認・格子定数の探索)",
                                   "Create one input per value of a single setting, such as the cutoff or k-points (convergence tests, lattice constants)"))
        self.act_scan.triggered.connect(self.open_scan_dialog)
        self.act_stages = QAction(icons.icon("stages", 32), L("段階に分けて生成…", "Generate a staged calculation…"), self)
        self.act_stages.setIconText(L("段階に分けて生成", "Stages"))
        self.act_stages.setToolTip(L("最小化 → NVT → NPT → 本計算のように、段階ごとのディレクトリと、段階を順に実行する submit.sh を作ります",
                                     "Create one directory per stage (such as minimization, NVT, NPT, production) and a submit.sh that runs them in order"))
        self.act_stages.triggered.connect(self.open_stages_dialog)
        from adit.web import prep23 as P
        self.act_batch: dict[str, QAction] = {}
        for kind, icon_name, short_ja, short_en in (("compare", "batch_compare", "組にして比べる", "Compare a set"),
                                                    ("conformers", "batch_conformers", "配座の候補", "Conformers"),
                                                    ("neb", "batch_neb", "反応経路 (NEB)", "Reaction path"),
                                                    ("phonons", "batch_phonons", "フォノン", "Phonons"),
                                                    ("elastic", "batch_elastic", "弾性定数", "Elastic constants"),
                                                    ("ts", "batch_ts", "遷移状態と IRC", "TS and IRC")):
            a = QAction(icons.icon(icon_name, 32), P.kind_title(kind) + "…", self)
            a.setIconText(L(short_ja, short_en)); a.setToolTip(P.kind_note(kind))
            a.triggered.connect(lambda _c=False, k=kind: self.open_batch_dialog(k))
            self.act_batch[kind] = a
        self.act_continue = QAction(icons.icon("continue", 32), L("前の計算の続きを作る…", "Continue a previous calculation…"), self)
        self.act_continue.setIconText(L("続きを作る", "Continue"))
        self.act_continue.setToolTip(L("ADIT が生成して実行したディレクトリを選び、その最終構造 (MD なら速度も) から続きの計算を画面に読み込みます",
                                       "Pick a directory generated by ADIT and already run; its final structure (and MD velocities) is loaded as a new calculation"))
        self.act_continue.triggered.connect(self.open_continue_dialog)
        self.act_template_load = QAction(icons.icon("template_load", 32), L("研究室の雛形を読み込む…", "Load a group template…"), self)
        self.act_template_load.setIconText(L("雛形を読み込む", "Load template"))
        self.act_template_load.setToolTip(L("雛形 (構造を除いた計算の条件) を一覧から選んで読み込みます。構造は今のままです",
                                            "Load a template (calculation conditions without a structure) from the list; the structure is kept"))
        self.act_template_load.triggered.connect(self.open_template_load_dialog)
        self.act_template_save = QAction(icons.icon("template_save", 32), L("今の設定を雛形として保存…", "Save the current settings as a template…"), self)
        self.act_template_save.setIconText(L("雛形として保存", "Save as template"))
        self.act_template_save.setToolTip(L("いまの画面の設定から構造を除いたものを、研究室の雛形として保存します (置き場所は環境設定の templates_dir)",
                                            "Save the current settings without the structure as a group template (folder: templates_dir in the preferences)"))
        self.act_template_save.triggered.connect(self.open_template_save_dialog)
        self.act_convert = QAction(icons.icon("file", 32), L("構造・計算コードを変換…", "Convert structures and calculation codes…"), self)
        self.act_convert.setIconText(L("変換", "Convert"))
        self.act_convert.setToolTip(L("構造ファイルの形式を変換するか、共通条件を保って別の計算コードの入力を生成します",
                                      "Convert a structure format, or generate input for another code while preserving shared settings"))
        self.act_convert.triggered.connect(self.open_convert_dialog)
        self.act_input_audit = QAction(icons.icon("file", 32), L("既存入力・生成入力・計算結果を点検…", "Inspect existing inputs, generated inputs, and runs…"), self)
        self.act_input_audit.setIconText(L("入力と結果の点検", "Inspect inputs/results"))
        self.act_input_audit.setToolTip(L("既存入力の読み込み、生成した入力の読み戻し照合、コードをまたぐ複数の計算の前提の突き合わせができます",
                                         "Import existing inputs, verify mapped generated fields, or audit mechanical prerequisites across runs"))
        self.act_input_audit.triggered.connect(self.open_input_audit_dialog)
        self._lang_actions = self._choice_actions((("ja", "日本語", "日本語"), ("en", "English", "English")), self.cfg.language, self._set_language, "")
        for _code, a in self._lang_actions:
            a.setIcon(icons.icon("language", 32))
        self._theme_actions = self._choice_actions((("auto", "システムに従う", "システムに従う"), ("light", "ライト", "ライト"), ("dark", "ダーク", "ダーク")),
                                                   self.cfg.theme, self._set_theme, "theme_")
        self._frame_actions = self._choice_actions((("auto", "自動 (Linux では ADIT が描く)", "自動"), ("custom", "ADIT が描く (ボタンにカーソルで色が点く)", "ADIT が描く"),
                                                    ("native", "OS に任せる", "OS の枠")), self.cfg.window_frame, self._set_frame, "frame_")
        self.act_readme = QAction(icons.icon("help", 32), "README を開く", self); self.act_readme.setIconText("README"); self.act_readme.triggered.connect(self._open_readme)
        self.act_about = QAction(icons.icon("about", 32), "ADIT について", self); self.act_about.triggered.connect(self._about)
        self.act_palette = QAction(icons.icon("search", 32), L("コマンドパレット…", "Command palette…"), self)
        self.act_palette.setIconText(L("コマンド", "Commands")); self.act_palette.setShortcut("Ctrl+K")
        self.act_palette.setToolTip(L("操作・欄・プリセットを名前で探して実行します", "Find a command, a field or a preset by name and run it"))
        self.act_palette.triggered.connect(self.open_palette)
        self.act_shortcuts = QAction(icons.icon("keyboard", 32), L("キーボードショートカット一覧", "Keyboard shortcuts"), self)
        self.act_shortcuts.setIconText(L("ショートカット", "Shortcuts")); self.act_shortcuts.setShortcuts(["Ctrl+/", "?"])
        self.act_shortcuts.triggered.connect(self.open_shortcuts)

        rb = self.ribbon = Ribbon(collapsed=bool(self.cfg.ribbon_collapsed))
        rb.collapsed_changed.connect(self._remember_ribbon)
        p = rb.add_page("ファイル")
        g = p.add_group("計算設定"); g.add_large(self.act_open); g.add_large(self.act_save)
        g = p.add_group(L("前の計算と雛形", "Previous runs and templates")); g.add_large(self.act_continue)
        g.add_large(self.act_template_load); g.add_large(self.act_template_save)
        g = p.add_group(L("変換と点検", "Convert and inspect")); g.add_large(self.act_convert); g.add_large(self.act_input_audit)
        g = p.add_group("アプリ"); g.add_large(self.act_quit)
        p = rb.add_page("挿入")
        g = p.add_group("構造"); g.add_large(self.act_draw); g.add_large(self.act_load_structure); g.add_large(self.act_add_component)
        p = rb.add_page("表示")
        g = p.add_group("メインの面")
        for a in self._mode_actions:
            g.add_large(a)
        g = p.add_group("右側のタブ")
        for a in self._tab_actions:
            g.add_large(a)
        p = rb.add_page("実行")
        g = p.add_group("入力"); g.add_large(self.act_generate); g.add_large(self.act_scan); g.add_large(self.act_stages)
        g = p.add_group("まとめて作る")
        g.add_small([self.act_batch[k] for k in ("compare", "conformers", "neb")])
        g.add_small([self.act_batch[k] for k in ("phonons", "elastic", "ts")])
        p = rb.add_page("設定")
        g = p.add_group("環境設定"); g.add_large(self.act_settings); g.add_small([self.act_reload])
        g = p.add_group("言語 (再起動後に反映)"); g.add_small([a for _c, a in self._lang_actions])
        g = p.add_group("テーマ (再起動後に反映)"); g.add_small([a for _c, a in self._theme_actions])
        g = p.add_group("ウィンドウの枠 (再起動後に反映)"); g.add_small([a for _c, a in self._frame_actions])
        p = rb.add_page("ヘルプ")
        g = p.add_group(L("探す", "Find")); g.add_large(self.act_palette); g.add_large(self.act_shortcuts)
        g = p.add_group("ヘルプ"); g.add_large(self.act_readme); g.add_large(self.act_about)
        for a in (self.act_back, self.act_forward, *rb.all_actions()):
            if not a.shortcut().isEmpty():
                self.addAction(a)

        self.quick_access = QuickAccessBar([self.act_back, self.act_forward])
        self.top_area = QWidget(); self.top_area.setObjectName("top_area")
        tl = QVBoxLayout(self.top_area); tl.setContentsMargins(0, 0, 0, 0); tl.setSpacing(0)
        self.titlebar = None
        if titlebar.use_custom_frame(self.cfg.window_frame):
            self.titlebar = titlebar.install(self, self.quick_access); tl.addWidget(self.titlebar)
        else:
            rb.set_leading(self.quick_access)
        tl.addWidget(rb)
        self.setMenuWidget(self.top_area)
        self._apply_shortcut_tips()

    def _remember_ribbon(self, collapsed: bool) -> None:
        if bool(self.cfg.ribbon_collapsed) != collapsed:
            self.cfg.ribbon_collapsed = collapsed; self._save_cfg()

    def _apply_shortcut_tips(self) -> None:
        from adit.gui.palette import with_shortcut

        for a in self.findChildren(QAction):
            if not a.shortcut().isEmpty() and a.text():
                a.setToolTip(with_shortcut(a.toolTip(), a))

    def shortcut_rows(self) -> list[tuple[str, str]]:
        from adit.gui.palette import shortcut_text

        rows, seen = [], set()
        for a in self.findChildren(QAction):
            key = shortcut_text(a)
            if key and a.text() and key not in seen:
                seen.add(key); rows.append((key, a.text().rstrip("…")))
        rows.append(("Ctrl+S", L("ワークスペース: 開いているファイルを保存", "Workspace: save the open file")))
        return sorted(rows, key=lambda r: (len(r[0].split("+")), r[0]))

    def palette_items(self) -> list:
        from adit.gui.i18n import _EN
        from adit.gui.palette import Item, shortcut_text
        from adit.structure import preset_search_text, pretty_formula

        items = []
        names = [tr(t) for t in ("構造", "計算条件", "解析", "ワークスペース")]
        for a in [self.act_back, self.act_forward, *self.ribbon.all_actions()]:
            if not a.text() or a is self.act_palette:
                continue
            en = _EN.get(a.text(), "")
            items.append(Item("action", a.text().rstrip("…"), shortcut_text(a), a.trigger, (a.iconText(), en, _EN.get(a.iconText(), "")),
                              a.icon(), a.isEnabled()))
        for w, mode, group in self.field_entries():
            key = str(w.property("adit_key"))
            where = names[mode] + (f" › {tr(group)}" if group else "")
            items.append(Item("field", tr(key), where, lambda w=w: self.jump_to_label(w), (key, _EN.get(key, ""), group, tr(group))))
        for name in getattr(self.structure, "_preset_names", []):
            items.append(Item("preset", pretty_formula(name), L("プリセットの分子", "Preset molecule"), lambda n=name: self.use_preset(n),
                              (preset_search_text(name),), icons.icon("preset", 16)))
        return items

    def use_preset(self, name: str) -> None:
        self.set_mode(self.MODE_STRUCTURE)
        self.structure.set_source("preset")
        if hasattr(self.structure, "preset_search"):
            self.structure.preset_search.setText("")
        self.structure.preset.setCurrentIndex(max(0, self.structure.preset.findData(name)))
        self.structure.preset.setFocus(Qt.FocusReason.OtherFocusReason)

    def open_palette(self):
        from adit.gui.i18n import translate_widgets
        from adit.gui.palette import CommandPalette

        dlg = CommandPalette(self.palette_items(), self)
        translate_widgets(dlg)
        dlg.place_over(self)
        dlg.show(); dlg.search.setFocus(Qt.FocusReason.OtherFocusReason)
        return dlg

    def open_shortcuts(self):
        from adit.gui.i18n import translate_widgets
        from adit.gui.palette import ShortcutsDialog

        dlg = ShortcutsDialog(self.shortcut_rows(), self)
        translate_widgets(dlg)
        dlg.show()
        return dlg

    def _save_cfg(self) -> None:
        try:
            from adit.config import save_config
            save_config(self.cfg, self.cfg_path)
        except OSError as ex:
            self.statusBar().showMessage(L(f"設定を保存できません: {ex}", f"cannot save the settings: {ex}"))

    def _set_language(self, code: str) -> None:
        self.cfg.language = code; self._save_cfg()
        self.statusBar().showMessage(L(f"言語を {code} にしました。次回の起動から反映されます", f"language set to {code}; takes effect after restart"))

    def _set_theme(self, code: str) -> None:
        self.cfg.theme = code; self._save_cfg()
        for c, a in self._theme_actions:
            a.setChecked(c == code)
        self.statusBar().showMessage(L(f"テーマを {code} にしました。次回の起動から反映されます", f"theme set to {code}; takes effect after restart"))

    def _set_frame(self, code: str) -> None:
        self.cfg.window_frame = code; self._save_cfg()
        for c, a in self._frame_actions:
            a.setChecked(c == code)
        self.statusBar().showMessage(L(f"ウィンドウの枠を {code} にしました。次回の起動から反映されます", f"window frame set to {code}; takes effect after restart"))

    def save_spec(self) -> None:
        try:
            spec = self.current_spec()
        except Exception as ex:
            QMessageBox.critical(self, tr("生成できません"), str(ex)); return
        path, _ = QFileDialog.getSaveFileName(self, "spec.json", "spec.json", "JSON (*.json)")
        if path:
            try:
                Path(path).write_text(spec.model_dump_json(indent=2), encoding="utf-8")
            except OSError as ex:
                QMessageBox.critical(self, L("保存できません", "Cannot save"), str(ex)); return
            self.statusBar().showMessage(L(f"{path} に保存しました", f"saved to {path}"))

    def _insert_structure_file(self) -> None:
        self.structure.set_source("file"); self.structure._browse()

    def _insert_component(self) -> None:
        self.structure.set_source("mixture"); self.structure.mixture.btn_add.click()

    def _open_readme(self) -> None:
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices
        for cand in (Path(__file__).resolve().parents[3] / "README.md", Path.cwd() / "README.md"):
            if cand.is_file():
                QDesktopServices.openUrl(QUrl.fromLocalFile(str(cand))); return
        QDesktopServices.openUrl(QUrl("https://github.com/"))

    def _about(self) -> None:
        from adit import __version__
        QMessageBox.about(self, tr("ADIT について"), f"ADIT {__version__}\nAtomistic Design and Interpretation Toolkit\n\n"
                          + L("計算化学の入力の準備、実行、解析を 1 つの画面で行うためのツール。MIT ライセンス。", "Prepare, run and analyze computational-chemistry calculations in one place. MIT license."))

    def open_spec(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "spec.json", "", "spec (spec.json);;JSON (*.json)")
        if not path:
            return
        try:
            spec = load_project(Path(path).parent) if Path(path).name == "spec.json" else CalculationSpec.load(path)
        except Exception as ex:
            QMessageBox.critical(self, tr("読めません"), str(ex))
            return
        if spec.method.code not in CODES:
            QMessageBox.critical(self, tr("読めません"),
                                 L(f"この計算コード ({spec.method.code}) は GUI では編集できません。CLI で spec.json を使ってください。",
                                   f"This code ({spec.method.code}) cannot be edited in the GUI. Use spec.json with the CLI."))
            return
        self.apply_spec(spec)
        self.runtime.outdir.setText(str(Path(path).parent))

    def apply_spec(self, spec: CalculationSpec, *, origin: bool = True) -> None:
        if origin:
            self.origin = PrepOrigin(spec)
        from adit.structure import StructureError
        self._applying = True
        try:
            try:
                self.structure.set_structure(spec.structure)
            except StructureError as ex:
                QMessageBox.critical(self, tr("読めません"), str(ex)); return
            self._on_context()
            self.method.set_method(spec.method)
            from adit.spec import KPoints
            self.kpoints.set_kpoints(spec.kpoints if spec.kpoints is not None else KPoints())
            self.task.set_task(spec.task)
            self.runtime.set_runtime(spec.runtime)
        finally:
            self._applying = False
        self.refresh_preview()
        unshown = [*self.method.unshown, *self.runtime.unshown]
        if unshown:
            self.statusBar().showMessage(L("次の欄は画面で表せないので変わりました: " + "、".join(unshown),
                                           "These fields cannot be shown here and were changed: " + ", ".join(unshown)))

    def _record_history(self) -> None:
        if self._applying:
            return
        try:
            snap = self.current_spec().model_dump_json(exclude={"meta"})
        except Exception:
            return
        if self._hist_pos >= 0 and self._history[self._hist_pos] == snap:
            return
        del self._history[self._hist_pos + 1:]
        self._history.append(snap); self._hist_pos = len(self._history) - 1
        del self._history[:-50]; self._hist_pos = len(self._history) - 1
        self._update_history_buttons()

    def _update_history_buttons(self) -> None:
        self.act_back.setEnabled(self._hist_pos > 0)
        self.act_forward.setEnabled(0 <= self._hist_pos < len(self._history) - 1)

    def _goto_history(self, pos: int) -> None:
        if not (0 <= pos < len(self._history)):
            return
        self._hist_pos = pos
        self._applying = True
        try:
            self.apply_spec(CalculationSpec.model_validate_json(self._history[pos]), origin=False)
        finally:
            self._applying = False
        self._update_history_buttons()

    def go_back(self) -> None:
        self._goto_history(self._hist_pos - 1)

    def go_forward(self) -> None:
        self._goto_history(self._hist_pos + 1)

    def open_settings(self) -> None:
        from adit.gui.settings_dialog import SettingsDialog

        dlg = SettingsDialog(self.cfg_path, self)
        if dlg.exec():
            self.reload_config()

    def reload_config(self) -> None:
        try:
            self.cfg = load_config(self.cfg_path)
        except ConfigError as ex:
            QMessageBox.critical(self, tr("設定を読めません"), str(ex)); return
        self.method.reload_sets(self.cfg.sk_root, self.cfg)
        self.runtime.reload_profiles(self.cfg)
        self._on_context()                       # the code panels keep their own cfg (POTCAR root etc.)
        self._sync_menu_checks()
        self.refresh_preview()

    def closeEvent(self, event) -> None:  # noqa: N802 (Qt)
        if self.workspace.editor.dirty and not self.workspace._ask_discard():
            event.ignore()
            return
        self.workspace.editor.discard()          # do not ask again when Qt closes the window a second time
        self.workspace.close_session()
        super().closeEvent(event)

    def _sync_menu_checks(self) -> None:
        for items, value in ((self._lang_actions, self.cfg.language), (self._theme_actions, self.cfg.theme), (self._frame_actions, self.cfg.window_frame)):
            for code, a in items:
                a.setChecked(code == value)


