"""Workspace tab: a file tree, a text editor and a terminal on the same screen."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QDir, QModelIndex, QRect, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontDatabase, QIcon, QKeySequence, QPainter, QPixmap, QShortcut
from PySide6.QtWidgets import (QComboBox, QFileSystemModel, QHBoxLayout, QInputDialog, QLabel, QLineEdit, QMenu, QMessageBox,
                               QPlainTextEdit, QPushButton, QSplitter, QTreeView, QVBoxLayout, QWidget)

from adit.gui.style import GROUP_SPACING, PANEL_MARGIN
from adit.gui.terminal_pane import TerminalTabs
from adit.lang import L

MAX_EDIT_BYTES = 2_000_000
TEXT_SUFFIXES = {".txt", ".md", ".json", ".toml", ".yaml", ".yml", ".sh", ".py", ".hsd", ".gen", ".xyz", ".extxyz",
                 ".in", ".inp", ".out", ".log", ".dat", ".csv", ".cfg", ".conf", ".ini", ".gjf", ".com", ".nw",
                 ".control", ".mdp", ".top", ".itp", ".gro", ".pdb", ".cif", ".POSCAR", ".KPOINTS", ".INCAR"}
NAME_ONLY = {"INCAR", "POSCAR", "CONTCAR", "KPOINTS", "POTCAR", "OUTCAR", "README", "LICENSE", "Makefile"}


def _looks_like_text(path: Path) -> bool:
    if path.name in NAME_ONLY or path.suffix.lower() in TEXT_SUFFIXES:
        return True
    try:
        chunk = path.open("rb").read(4096)
    except OSError:
        return False
    return b"\0" not in chunk


class Editor(QPlainTextEdit):
    """Plain text editor for one file at a time."""

    dirty_changed = Signal(bool)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        font.setStyleHint(QFont.StyleHint.Monospace)
        self.setFont(font)
        self.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.path: Path | None = None
        self._dirty = False
        self._newline = "\n"
        self.textChanged.connect(self._on_changed)

    def _on_changed(self) -> None:
        if self.path is not None and not self._dirty:
            self._dirty = True
            self.dirty_changed.emit(True)

    @property
    def dirty(self) -> bool:
        return self._dirty

    def open_file(self, path: Path) -> str:
        """Show the file. Returns an empty string, or the reason it cannot be shown."""
        try:
            size = path.stat().st_size
        except OSError as ex:
            return str(ex)
        if size > MAX_EDIT_BYTES:
            return L(f"大きすぎて開けません ({size // 1024} KB)。ターミナルで開いてください",
                     f"too large to open ({size // 1024} KB); open it in the terminal")
        if not _looks_like_text(path):
            return L("テキストではないので開けません", "not a text file")
        try:
            raw = path.read_bytes()
        except OSError as ex:
            return str(ex)
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError:
            # Show it, but never write it back: saving would replace the unknown bytes.
            self.path = None
            self.setPlainText(raw.decode("utf-8", errors="replace"))
            self.setReadOnly(True)
            self._dirty = False
            self.dirty_changed.emit(False)
            return L("文字コードが UTF-8 でないので編集できません", "not UTF-8, so it cannot be edited here")
        self._newline = "\r\n" if "\r\n" in text else "\n"
        self.path = None                      # textChanged during loading is not an edit
        self.setReadOnly(False)
        self.setPlainText(text.replace("\r\n", "\n"))
        self.path, self._dirty = path, False
        self.dirty_changed.emit(False)
        return ""

    def discard(self) -> None:
        self.path = None
        self.clear()
        self._dirty = False
        self.dirty_changed.emit(False)

    def save(self) -> str:
        if self.path is None:
            return L("開いているファイルがありません", "no file is open")
        try:
            with open(self.path, "w", encoding="utf-8", newline="") as f:
                f.write(self.toPlainText().replace("\n", self._newline))
        except OSError as ex:
            return str(ex)
        self._dirty = False
        self.dirty_changed.emit(False)
        return ""


class FileIcons(QFileSystemModel):
    """ファイルの種類が見て分かるように、拡張子ごとにしるしを付ける。"""

    KIND = {
        ".hsd": ("in", "#2f7ae5"), ".in": ("in", "#2f7ae5"), ".inp": ("in", "#2f7ae5"), ".gjf": ("in", "#2f7ae5"),
        ".nw": ("in", "#2f7ae5"), ".dat": ("in", "#2f7ae5"), ".mdp": ("in", "#2f7ae5"), ".conf": ("in", "#2f7ae5"),
        ".log": ("out", "#12876f"), ".out": ("out", "#12876f"), ".tag": ("out", "#12876f"),
        ".xyz": ("st", "#d2691e"), ".gen": ("st", "#d2691e"), ".cif": ("st", "#d2691e"), ".pdb": ("st", "#d2691e"),
        ".gro": ("st", "#d2691e"), ".extxyz": ("st", "#d2691e"),
        ".json": ("cfg", "#7b4fd0"), ".toml": ("cfg", "#7b4fd0"), ".yaml": ("cfg", "#7b4fd0"),
        ".sh": ("run", "#b8860b"), ".py": ("py", "#3776ab"), ".j2": ("cfg", "#7b4fd0"),
        ".png": ("img", "#c2185b"), ".svg": ("img", "#c2185b"), ".csv": ("tbl", "#12876f"), ".md": ("doc", "#555f6d"),
    }
    NAMED = {"INCAR": ("in", "#2f7ae5"), "POSCAR": ("st", "#d2691e"), "CONTCAR": ("st", "#d2691e"),
             "KPOINTS": ("in", "#2f7ae5"), "POTCAR": ("in", "#2f7ae5"), "OUTCAR": ("out", "#12876f"),
             "README.txt": ("doc", "#555f6d"), "submit.sh": ("run", "#b8860b")}

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if role == Qt.ItemDataRole.DecorationRole and index.column() == 0:
            path = Path(self.filePath(index))
            if path.is_dir():
                return super().data(index, role)
            kind = self.NAMED.get(path.name) or self.KIND.get(path.suffix.lower())
            if kind is not None:
                return self._badge(*kind)
        return super().data(index, role)

    _CACHE: dict[tuple[str, str], QIcon] = {}

    @classmethod
    def _badge(cls, text: str, color: str) -> QIcon:
        got = cls._CACHE.get((text, color))
        if got is not None:
            return got
        size = 16
        pix = QPixmap(size, size)
        pix.fill(Qt.GlobalColor.transparent)
        p = QPainter(pix)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor(color))
        p.drawRoundedRect(1, 2, size - 2, size - 4, 3, 3)
        font = QFont(); font.setPointSizeF(5.5); font.setBold(True)
        p.setFont(font); p.setPen(QColor("#ffffff"))
        p.drawText(QRect(1, 2, size - 2, size - 4), int(Qt.AlignmentFlag.AlignCenter), text)
        p.end()
        icon = QIcon(pix)
        cls._CACHE[(text, color)] = icon
        return icon


class WorkspacePanel(QWidget):
    """File tree on the left, editor above the terminal on the right."""

    def __init__(self, root: Path | str | None = None, dark: bool = False, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.root = Path(root).expanduser() if root else Path.home()

        self.model = FileIcons(self)
        self.model.setRootPath(str(self.root))
        self.model.setFilter(QDir.Filter.AllEntries | QDir.Filter.NoDotAndDotDot | QDir.Filter.Hidden)
        self.tree = QTreeView()
        self.tree.setModel(self.model)
        self.tree.setRootIndex(self.model.index(str(self.root)))
        for column in (1, 2, 3):
            self.tree.hideColumn(column)
        self.tree.setHeaderHidden(True)
        self.tree.setIndentation(16)
        self.tree.setAnimated(False)
        self.tree.setRootIsDecorated(True)
        self.tree.setExpandsOnDoubleClick(False)
        self.tree.setStyleSheet("QTreeView { show-decoration-selected: 1; }")
        self.tree.setUniformRowHeights(True)
        self.model.directoryLoaded.connect(self._expand_new)
        self.tree.setDragDropMode(QTreeView.DragDropMode.InternalMove)   # ドラッグで移動できる
        self.tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._menu)
        self.tree.doubleClicked.connect(self._open_index)

        self.path_label = QLabel(str(self.root))
        self.path_label.setObjectName("hint")
        self.btn_up = QPushButton(L("上へ", "Up"))
        self.btn_up.clicked.connect(lambda: self.set_root(self.root.parent))

        self.editor = Editor()
        self.file_label = QLabel(L("ファイルを選ぶと、ここで編集できます", "pick a file to edit it here"))
        self.file_label.setObjectName("hint")
        self.btn_save = QPushButton(L("保存", "Save"))
        self.btn_save.setEnabled(False)
        self.btn_save.clicked.connect(self.save)
        self.editor.dirty_changed.connect(self._on_dirty)

        self.layout_choice = QComboBox()
        for value, text in (("v", L("エディタが上", "Editor on top")), ("v_rev", L("ターミナルが上", "Terminal on top")),
                            ("h", L("エディタが左", "Editor on the left")), ("h_rev", L("ターミナルが左", "Terminal on the left"))):
            self.layout_choice.addItem(text, value)
        self.layout_choice.setToolTip(L("エディタとターミナルの並べ方", "How the editor and the terminal are arranged"))
        self.layout_choice.currentIndexChanged.connect(lambda *_: self.apply_layout())
        self.btn_fold_editor = QPushButton(L("エディタを畳む", "Collapse the editor"))
        self.btn_fold_editor.setCheckable(True)
        self.btn_fold_editor.toggled.connect(lambda *_: self.apply_layout())
        self.btn_fold_terminal = QPushButton(L("ターミナルを畳む", "Collapse the terminal"))
        self.btn_fold_terminal.setCheckable(True)
        self.btn_fold_terminal.toggled.connect(lambda *_: self.apply_layout())

        self.terminal = TerminalTabs(cwd=self.root, dark=dark)
        self.btn_here = QPushButton(L("ここへ移動 (cd)", "cd here"))
        self.btn_here.clicked.connect(lambda: self.terminal.send(f"cd {self._quoted(self.root)}\r"))
        self.btn_restart = QPushButton(L("シェルを起動し直す", "Restart the shell"))
        self.btn_restart.setVisible(False)
        self.btn_restart.clicked.connect(self._restart_terminal)
        self._watch_terminal()

        save = QShortcut(QKeySequence.StandardKey.Save, self)      # Ctrl+S
        save.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        save.activated.connect(self.save)

        left = QWidget()
        ll = QVBoxLayout(left); ll.setContentsMargins(0, 0, 0, 0); ll.setSpacing(4)
        top = QHBoxLayout(); top.addWidget(self.btn_up); top.addWidget(self.path_label, 1)
        ll.addLayout(top); ll.addWidget(self.tree, 1)

        self.edit_box = QWidget()
        el = QVBoxLayout(self.edit_box); el.setContentsMargins(0, 0, 0, 0); el.setSpacing(4)
        head = QHBoxLayout()
        head.addWidget(self.file_label, 1); head.addWidget(self.btn_fold_editor); head.addWidget(self.btn_save)
        el.addLayout(head); el.addWidget(self.editor, 1)

        self.term_box = QWidget()
        tl = QVBoxLayout(self.term_box); tl.setContentsMargins(0, 0, 0, 0); tl.setSpacing(4)
        thead = QHBoxLayout()
        label = QLabel(L("ターミナル", "Terminal")); label.setObjectName("section")
        thead.addWidget(label); thead.addStretch()
        thead.addWidget(self.layout_choice); thead.addWidget(self.btn_fold_terminal)
        thead.addWidget(self.btn_restart); thead.addWidget(self.btn_here)
        tl.addLayout(thead); tl.addWidget(self.terminal, 1)

        self.right = QSplitter(Qt.Orientation.Vertical)
        self.right.addWidget(self.edit_box); self.right.addWidget(self.term_box)

        self.split = QSplitter(Qt.Orientation.Horizontal)
        self.split.addWidget(left); self.split.addWidget(self.right)
        self.split.setStretchFactor(0, 0); self.split.setStretchFactor(1, 1)
        left.setMinimumWidth(150)
        self.split.setSizes([230, 1070])            # ツリーは細く、エディタとターミナルを広く
        lay = QVBoxLayout(self)
        lay.setContentsMargins(PANEL_MARGIN, 8, PANEL_MARGIN, PANEL_MARGIN); lay.setSpacing(GROUP_SPACING)
        lay.addWidget(self.split)
        self.apply_layout()

    # ---- layout ----
    def apply_layout(self) -> None:
        """並べ方 (上下・左右・入れ替え) と、畳む指定を反映する。"""
        choice = self.layout_choice.currentData() or "v"
        vertical = choice.startswith("v")
        reverse = choice.endswith("_rev")
        self.right.setOrientation(Qt.Orientation.Vertical if vertical else Qt.Orientation.Horizontal)
        first, second = (self.term_box, self.edit_box) if reverse else (self.edit_box, self.term_box)
        if self.right.widget(0) is not first:
            self.right.insertWidget(0, first)
            self.right.insertWidget(1, second)
        fold_editor, fold_terminal = self.btn_fold_editor.isChecked(), self.btn_fold_terminal.isChecked()
        if fold_editor and fold_terminal:           # 両方は畳めない (どちらかは残す)
            self.btn_fold_terminal.setChecked(False)
            fold_terminal = False
        self.editor.setVisible(not fold_editor)
        self.terminal.setVisible(not fold_terminal)
        self.btn_fold_editor.setText(L("エディタを開く", "Expand the editor") if fold_editor
                                     else L("エディタを畳む", "Collapse the editor"))
        self.btn_fold_terminal.setText(L("ターミナルを開く", "Expand the terminal") if fold_terminal
                                       else L("ターミナルを畳む", "Collapse the terminal"))
        total = max(400, (self.right.height() if vertical else self.right.width()) or 800)
        if fold_editor:
            self.right.setSizes([40, total - 40] if not reverse else [total - 40, 40])
        elif fold_terminal:
            self.right.setSizes([total - 40, 40] if not reverse else [40, total - 40])
        else:
            self.right.setSizes([total // 2, total - total // 2])

    def _expand_new(self, path: str) -> None:
        """読み込めたフォルダを、根から 3 階層まで開く。"""
        index = self.model.index(path)
        if not index.isValid():
            return
        depth, walk = 0, index
        root = self.model.index(str(self.root))
        while walk.isValid() and walk != root and depth < 8:
            walk = walk.parent(); depth += 1
        if depth < 3:
            self.tree.expand(index)

    # ---- helpers ----
    @staticmethod
    def _quoted(path: Path) -> str:
        import shlex

        return shlex.quote(str(path))

    def _selected(self) -> Path | None:
        index = self.tree.currentIndex()
        return Path(self.model.filePath(index)) if index.isValid() else None

    def set_root(self, path: Path | str) -> None:
        path = Path(path).expanduser()
        if not path.is_dir():
            return
        self.root = path
        self.terminal.cwd = path
        self.model.setRootPath(str(path))
        self.tree.setRootIndex(self.model.index(str(path)))
        self.path_label.setText(str(path))

    def set_dark(self, dark: bool) -> None:
        self.terminal.set_dark(dark)

    def _watch_terminal(self) -> None:
        for i in range(self.terminal.tabs.count()):
            self._watch_view(self.terminal.tabs.widget(i))
        self.terminal.tab_added.connect(self._watch_view)
        self.terminal.tabs.currentChanged.connect(lambda *_: self._show_restart())

    def _watch_view(self, view) -> None:
        view.terminal.finished.connect(self._show_restart)

    def _show_restart(self) -> None:
        current = self.terminal.current
        self.btn_restart.setVisible(current is not None and (current.session is None or not current.session.alive))

    def _restart_terminal(self) -> None:
        current = self.terminal.current
        if current is not None:
            current.restart(self.root)
        self.btn_restart.setVisible(False)

    # ---- editor ----
    def _open_index(self, index: QModelIndex) -> None:
        path = Path(self.model.filePath(index))
        if path.is_dir():
            self.set_root(path)
            return
        self.open_file(path)

    def open_file(self, path: Path) -> str:
        if self.editor.dirty and not self._ask_discard():
            return ""
        why = self.editor.open_file(path)
        if why:
            self.file_label.setText(f"{path.name}: {why}")
            return why
        self.file_label.setText(str(path))
        return ""

    def _ask_discard(self) -> bool:
        answer = QMessageBox.question(self, L("保存していません", "Not saved"),
                                      L("編集した内容を保存していません。捨ててよいですか?",
                                        "The file has unsaved changes. Discard them?"))
        return answer == QMessageBox.StandardButton.Yes

    def _on_dirty(self, dirty: bool) -> None:
        self.btn_save.setEnabled(dirty)
        if self.editor.path is not None:
            self.file_label.setText(("* " if dirty else "") + str(self.editor.path))

    def save(self) -> str:
        why = self.editor.save()
        if why:
            QMessageBox.warning(self, L("保存できません", "Cannot save"), why)
        return why

    # ---- file operations ----
    def _menu(self, point) -> None:
        path = self._selected()
        if path is None:
            return
        menu = QMenu(self)
        act_open = menu.addAction(L("開く", "Open"))
        act_term = menu.addAction(L("ターミナルでここへ移動", "cd here in the terminal"))
        act_rename = menu.addAction(L("名前を変える…", "Rename…"))
        act_new = menu.addAction(L("フォルダを作る…", "New folder…"))
        chosen = menu.exec(self.tree.viewport().mapToGlobal(point))
        if chosen is act_open:
            self._open_index(self.tree.currentIndex())
        elif chosen is act_term:
            self.terminal.send(f"cd {self._quoted(path if path.is_dir() else path.parent)}\r")
        elif chosen is act_rename:
            name, ok = QInputDialog.getText(self, L("名前を変える", "Rename"), L("新しい名前", "New name"),
                                            QLineEdit.EchoMode.Normal, path.name)
            if ok:
                self.rename(path, name)
        elif chosen is act_new:
            name, ok = QInputDialog.getText(self, L("フォルダを作る", "New folder"), L("名前", "Name"))
            if ok and name.strip():
                base = path if path.is_dir() else path.parent
                try:
                    (base / name.strip()).mkdir()
                except OSError as ex:
                    QMessageBox.warning(self, L("作れません", "Cannot create"), str(ex))

    def rename(self, path: Path, name: str) -> str:
        name = name.strip()
        if not name or name == path.name:
            return ""
        target = path.with_name(name)
        if target.exists():
            why = L(f"{target.name} はすでにあります", f"{target.name} already exists")
        else:
            try:
                path.rename(target)
                why = ""
            except OSError as ex:
                why = str(ex)
        if why:
            QMessageBox.warning(self, L("変えられません", "Cannot rename"), why)
            return why
        if self.editor.path == path:
            self.editor.path = target
            self._on_dirty(self.editor.dirty)
        return ""

    def close_session(self) -> None:
        self.terminal.close_session()
