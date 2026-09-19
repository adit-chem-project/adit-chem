
from __future__ import annotations

import re
import tomllib
from pathlib import Path

import tomli_w
from pydantic import ValidationError as PydanticError
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog, QFormLayout, QGroupBox, QHBoxLayout, QHeaderView,
                               QLabel, QLineEdit, QMessageBox, QPlainTextEdit, QPushButton, QScrollArea, QTableWidget, QTableWidgetItem,
                               QTabWidget, QVBoxLayout, QWidget)

from adit.config import Config, ConfigError, env_var, load_config
from adit.lang import L

FIELD_KEYS = ("sk_root", "pseudo_root", "cp2k_data", "templates_dir", "language", "theme", "window_frame", "default_profile")


def _choices() -> dict[str, list[tuple[str, str]]]:
    return {
        "language": [("ja", "日本語"), ("en", "English")],
        "theme": [("auto", L("システムに従う", "Follow the system")), ("light", L("ライト", "Light")), ("dark", L("ダーク", "Dark"))],
        "window_frame": [("auto", L("自動 (Linux では ADIT が描く)", "Automatic (ADIT draws it on Linux)")),
                         ("custom", L("ADIT が描く", "Drawn by ADIT")), ("native", L("OS に任せる", "Use the system frame"))],
    }


def _kind_label(kind: str) -> str:
    return {"direct": L("この PC で実行 (direct)", "This PC (direct)"),
            "pbs": L("クラスタ・PBS (pbs)", "Cluster, PBS (pbs)"),
            "slurm": L("クラスタ・Slurm (slurm)", "Cluster, Slurm (slurm)")}.get(kind, kind)


_TOML_MSG_JA = {
    "Invalid value": "値の書き方が正しくありません。文字列は \"…\" で囲みます",
    "Expected '=' after a key in a key/value pair": "項目の名前のあとに = がありません",
    "Invalid statement": "この行は「名前 = 値」の形になっていません",
    "Unclosed string": "文字列の \" が閉じていません",
    "Illegal character in string": "文字列の中に使えない文字があります。\\ は \\\\ と 2 つ重ねて書きます",
    "Unescaped '\\' in a string": "文字列の中の \\ は \\\\ と 2 つ重ねて書きます",
    "Cannot overwrite a value": "同じ項目が 2 回書かれています",
    "Cannot declare": "同じ見出し ([…]) が 2 回書かれています",
}


def describe_problem(text: str) -> str | None:
    try:
        data = tomllib.loads(text)
    except tomllib.TOMLDecodeError as ex:
        msg = getattr(ex, "msg", str(ex))
        line, col = getattr(ex, "lineno", None), getattr(ex, "colno", None)
        if line is None:
            import re as _re

            m = _re.search(r"at line (\d+), column (\d+)", str(ex))
            if m:
                line, col = int(m.group(1)), int(m.group(2))
                msg = str(ex).split(" (at line", 1)[0]
        ja = next((v for k, v in _TOML_MSG_JA.items() if msg.startswith(k)), "")
        where_ja = f"{line} 行目 {col} 文字目" if line else "どこか"
        where_en = f"line {line}, column {col}" if line else "somewhere"
        return L(f"環境設定ファイルの {where_ja}が TOML の書き方として正しくありません" + (f"。{ja}" if ja else "") + f" ({msg})",
                 f"The settings file is not valid TOML at {where_en} ({msg}).")
    try:
        Config.model_validate(data)
    except PydanticError as ex:
        lines = []
        for e in ex.errors():
            loc = ".".join(str(x) for x in e.get("loc", ()))
            lines.append(f"  {loc or L('全体', 'whole file')}: {e.get('msg', '')}")
        return L("環境設定ファイルの書き方は正しいのですが、次の項目の値が使えません。\n", "The settings file is valid TOML, but these values cannot be used:\n") + "\n".join(lines)
    return None


def _top_region_end(lines: list[str]) -> int:
    return next((i for i, s in enumerate(lines) if s.lstrip().startswith("[")), len(lines))


def set_top_level(text: str, values: dict) -> str:
    try:
        current = tomllib.loads(text)
    except tomllib.TOMLDecodeError:
        current = None
    lines = text.splitlines()
    for key, value in values.items():
        if current is not None and current.get(key) == value:
            continue
        new_line = tomli_w.dumps({key: value}).strip()
        end = _top_region_end(lines)
        # Same shape as config.set_top_level_value: replace the value, keep a trailing comment.
        pat = re.compile(rf"""^(\s*{re.escape(key)}\s*=\s*)("[^"]*"|'[^']*'|[^\s#]+)(.*)$""")
        idx = next((i for i in range(end) if pat.match(lines[i])), None)
        if idx is not None:
            m = pat.match(lines[idx])
            lines[idx] = m.group(1) + new_line.split("=", 1)[1].strip() + m.group(3)
        else:
            at = end
            while at > 0 and not lines[at - 1].strip():
                at -= 1
            lines.insert(at, new_line)
    out = "\n".join(lines)
    return out + "\n" if out and not out.endswith("\n") else out


class SettingsDialog(QDialog):

    def __init__(self, path: Path, parent: QWidget | None = None):
        super().__init__(parent)
        self.path = Path(path)
        self.setWindowTitle(L("設定", "Settings")); self.resize(780, 640)
        self.last_error = ""

        self.tabs = QTabWidget()
        self.tabs.addTab(self._build_basic(), L("基本", "Basic"))
        self.tabs.addTab(self._build_raw(), L("詳細 (TOML を直接編集)", "Advanced (edit the TOML directly)"))
        self._tab = 0
        self.tabs.currentChanged.connect(self._on_tab)

        where = QLabel(L(f"環境設定ファイル: {self.path}", f"Settings file: {self.path}")); where.setObjectName("hint")
        where.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse); where.setWordWrap(True)

        self.btn_cancel = QPushButton(L("キャンセル", "Cancel")); self.btn_cancel.clicked.connect(self.reject)
        self.btn_save = QPushButton(L("保存", "Save")); self.btn_save.setObjectName("primary"); self.btn_save.setDefault(True)
        self.btn_save.clicked.connect(self._save)
        buttons = QHBoxLayout(); buttons.addStretch(1); buttons.addWidget(self.btn_cancel); buttons.addWidget(self.btn_save)

        lay = QVBoxLayout(self); lay.setContentsMargins(16, 16, 16, 16); lay.setSpacing(10)
        lay.addWidget(self.tabs, 1); lay.addWidget(where); lay.addLayout(buttons)

        try:
            self.editor.setPlainText(self.path.read_bytes().decode("utf-8-sig", errors="replace") if self.path.is_file() else "")
        except OSError as ex:
            self.editor.setPlainText("")
            self.last_error = str(ex)
            self.raw_note.setObjectName("status_ng"); self.raw_note.setText(L(f"設定ファイルを読めません: {ex}", f"cannot read the settings file: {ex}"))
        problem = describe_problem(self.editor.toPlainText())
        if problem:
            self.tabs.setCurrentIndex(1); self._tab = 1
            self.raw_note.setText(problem); self.raw_note.setObjectName("status_ng")
        else:
            self._text_to_fields()

    def _desc(self, text: str) -> QLabel:
        w = QLabel(text); w.setObjectName("hint"); w.setWordWrap(True)
        return w

    def _field(self, form: QFormLayout, label: str, widget, desc: str, extra: QLabel | None = None) -> None:
        cell = QWidget(); cell.setObjectName("rowbox")
        box = QVBoxLayout(cell); box.setSpacing(2); box.setContentsMargins(0, 0, 0, 4)
        box.addWidget(widget); box.addWidget(self._desc(desc))
        if extra is not None:
            box.addWidget(extra)
        lab = QLabel(label); lab.setMinimumWidth(170); lab.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        form.addRow(lab, cell)

    def _folder_row(self, marker: str, title: str) -> tuple[QWidget, QLineEdit, QLabel]:
        edit = QLineEdit(); edit.setPlaceholderText(L("まだ選んでいません", "not chosen yet"))
        btn = QPushButton(L("フォルダを選ぶ…", "Choose folder…"))
        btn.clicked.connect(lambda: self._choose(edit, title, marker))
        row = QWidget(); row.setObjectName("rowbox")
        h = QHBoxLayout(row); h.setContentsMargins(0, 0, 0, 0); h.setSpacing(8); h.addWidget(edit, 1); h.addWidget(btn)
        note = QLabel(""); note.setObjectName("status_ng"); note.setWordWrap(True); note.hide()
        edit.textChanged.connect(lambda _t: self._check_folder(edit, note))
        return row, edit, note

    def _build_basic(self) -> QWidget:
        page = QWidget(); outer = QVBoxLayout(page); outer.setContentsMargins(12, 12, 12, 12); outer.setSpacing(14)

        g_files = QGroupBox(L("計算に使うファイルの置き場所", "Where the parameter files are")); f1 = QFormLayout(g_files)
        f1.setVerticalSpacing(10)
        row, self.sk_root, self.sk_note = self._folder_row(".skf", L("Slater-Koster パラメータの置き場所", "Slater-Koster parameter folder"))
        self._field(f1, L("Slater-Koster パラメータ", "Slater-Koster parameters"), row,
                    L("DFTB+ が使うパラメータ (.skf) のセットを入れたフォルダです。中に mio-1-1 などのフォルダが並んでいる親を選びます。",
                      "Folder holding the DFTB+ parameter sets (.skf). Choose the parent folder that contains mio-1-1 and similar."),
                    self.sk_note)
        row, self.pseudo_root, self.pseudo_note = self._folder_row(".upf", L("擬ポテンシャル (UPF) の置き場所", "Pseudopotential (UPF) folder"))
        self._field(f1, L("擬ポテンシャル (UPF)", "Pseudopotentials (UPF)"), row,
                    L("Quantum ESPRESSO が使う擬ポテンシャル (.UPF) を入れたフォルダです。QE を使わなければ空のままで構いません。",
                      "Folder holding the Quantum ESPRESSO pseudopotentials (.UPF). Leave it empty if you do not use QE."),
                    self.pseudo_note)
        row, self.cp2k_data, self.cp2k_note = self._folder_row(".cp2k-data", L("CP2K の data ディレクトリ", "CP2K data directory"))
        self._field(f1, L("CP2K の data ディレクトリ", "CP2K data directory"), row,
                    L("CP2K の基底関数と擬ポテンシャルのファイル (BASIS_MOLOPT、GTH_POTENTIALS など) があるフォルダです。"
                      "空のままなら、環境変数 CP2K_DATA_DIR、PATH にある cp2k の隣の share/cp2k/data の順に探します。",
                      "Folder holding the CP2K basis-set and pseudopotential files (BASIS_MOLOPT, GTH_POTENTIALS, ...). "
                      "If empty, CP2K_DATA_DIR and then share/cp2k/data next to the cp2k on PATH are searched."),
                    self.cp2k_note)
        row, self.templates_dir, self.templates_note = self._folder_row(".adit-no-marker", L("研究室の雛形の置き場所", "Group template folder"))
        self._field(f1, L("研究室の雛形", "Group templates"), row,
                    L("研究室で共有する雛形 (構造を除いた計算の条件、名前.json) を置くフォルダです。空のままなら、"
                      "この環境設定ファイルと同じ場所の templates/ だけを使います。ここに書いたフォルダを templates/ より先に探します。",
                      "Folder holding shared group templates (calculation conditions without a structure, name.json). If empty, only templates/ "
                      "next to this settings file is used; a folder given here is searched before templates/."),
                    self.templates_note)

        g_look = QGroupBox(L("画面", "Appearance")); f2 = QFormLayout(g_look); f2.setVerticalSpacing(10)
        self.combos: dict[str, QComboBox] = {}
        descs = {
            "language": (L("表示する言語", "Language"), L("メニューやボタンの言葉です。", "Language of menus and buttons.")),
            "theme": (L("テーマ", "Theme"), L("明るい画面 (ライト) か暗い画面 (ダーク) か。「システムに従う」なら OS の設定に合わせます。",
                                                   "Light or dark. \"Follow the system\" follows your OS setting.")),
            "window_frame": (L("ウィンドウの枠", "Window frame"), L("タイトルバーとボタンを誰が描くか。表示がおかしいときは「OS に任せる」を試してください。",
                                                                "Who draws the title bar. Try \"Use the system frame\" if the window looks wrong.")),
        }
        for key, items in _choices().items():
            cb = QComboBox()
            for value, text in items:
                cb.addItem(text, value)
            self.combos[key] = cb
            label, desc = descs[key]
            self._field(f2, label, cb, desc)
        later = self._desc(L("言語・色・枠は、ADIT を次に起動したときに変わります。", "Language, theme and frame take effect the next time ADIT starts."))
        overridden = [n for n in ("LANG", "THEME", "FRAME") if env_var(n)]
        if overridden:
            later.setText(later.text() + L(f" 環境変数 {', '.join('ADIT_' + n for n in overridden)} が設定されているので、そちらが優先されます。",
                                           f" The environment variable(s) {', '.join('ADIT_' + n for n in overridden)} take priority."))
        f2.addRow(later)

        g_run = QGroupBox(L("実行", "Running")); v3 = QVBoxLayout(g_run); v3.setSpacing(8)
        f3_box = QWidget(); f3_box.setObjectName("rowbox"); f3 = QFormLayout(f3_box); f3.setContentsMargins(0, 0, 0, 0); f3.setVerticalSpacing(10)
        v3.addWidget(f3_box)
        self.default_profile = QComboBox()
        self._field(f3, L("いつも使う実行先", "Default target"), self.default_profile,
                    L("新しく計算を作るときに最初に選ばれている実行先です。実行先は下の一覧にあります。",
                      "The target selected first for a new calculation. Targets are listed below."))
        self.profiles = QTableWidget(0, 3)
        self.profiles.setHorizontalHeaderLabels([L("名前", "Name"), L("種類", "Kind"), L("説明", "Description")])
        self.profiles.verticalHeader().setVisible(False)
        self.profiles.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.profiles.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        hh = self.profiles.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents); hh.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        v3.addWidget(QLabel(L("実行先 (プロファイル) の一覧", "Execution targets (profiles)")))
        v3.addWidget(self.profiles)
        v3.addWidget(self._desc(L("実行先を足したり変えたりするには「詳細 (TOML を直接編集)」タブを使います。クラスタの書き方の例は README の付録にあります。",
                               "To add or change a target, use the \"Advanced\" tab. Examples for clusters are in the appendix of the README.")))

        for g in (g_files, g_look, g_run):
            outer.addWidget(g)
        outer.addStretch(1)
        scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll.setWidget(page)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        return scroll

    def _build_raw(self) -> QWidget:
        page = QWidget(); lay = QVBoxLayout(page); lay.setContentsMargins(12, 12, 12, 12); lay.setSpacing(8)
        lay.addWidget(self._desc(L("環境設定ファイルの本文です。「基本」タブの欄の値はここに書き込まれています。"
                                   "ここで直した値は「基本」タブに戻ったときに欄へ入ります。# から後ろはメモ (コメント) として残ります。",
                                   "The settings file itself. Values from the Basic tab are already written here; "
                                   "edits made here are loaded into the Basic tab when you switch back. Text after # is kept as a comment.")))
        self.editor = QPlainTextEdit(); self.editor.setStyleSheet("font-family: 'SF Mono', Menlo, 'DejaVu Sans Mono', monospace;")
        self.editor.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        lay.addWidget(self.editor, 1)
        self.raw_note = QLabel(""); self.raw_note.setObjectName("hint"); self.raw_note.setWordWrap(True)
        lay.addWidget(self.raw_note)
        return page

    def _choose(self, edit: QLineEdit, title: str, marker: str) -> None:
        from adit.gui.panels.method_panel import choose_root

        d = choose_root(self, title, marker)
        if d:
            edit.setText(d)

    def _check_folder(self, edit: QLineEdit, note: QLabel) -> None:
        t = edit.text().strip()
        bad = bool(t) and not Path(t).expanduser().is_dir()
        note.setText(L("このフォルダは見つかりません。場所を確かめてください (保存はできます)。",
                       "This folder does not exist. Please check the location (you can still save).") if bad else "")
        note.setVisible(bad)

    def field_values(self) -> dict:
        return {
            "sk_root": self.sk_root.text().strip(),
            "pseudo_root": self.pseudo_root.text().strip(),
            "cp2k_data": self.cp2k_data.text().strip(),
            "templates_dir": self.templates_dir.text().strip(),
            "language": self.combos["language"].currentData(),
            "theme": self.combos["theme"].currentData(),
            "window_frame": self.combos["window_frame"].currentData(),
            "default_profile": self.default_profile.currentData(),
        }

    def _text_to_fields(self) -> None:
        cfg = Config.model_validate(tomllib.loads(self.editor.toPlainText()))
        self.sk_root.setText(cfg.sk_root); self.pseudo_root.setText(cfg.pseudo_root); self.cp2k_data.setText(cfg.cp2k_data)
        self.templates_dir.setText(cfg.templates_dir)
        for key, cb in self.combos.items():
            value = getattr(cfg, key)
            i = cb.findData(value)
            if i < 0:
                cb.addItem(value, value); i = cb.count() - 1
            cb.setCurrentIndex(i)
        self.default_profile.clear()
        for name in cfg.profiles:
            self.default_profile.addItem(name, name)
        if cfg.default_profile not in cfg.profiles:
            self.default_profile.addItem(L(f"{cfg.default_profile} (一覧に無い名前)", f"{cfg.default_profile} (not in the list)"), cfg.default_profile)
        self.default_profile.setCurrentIndex(self.default_profile.findData(cfg.default_profile))
        self.profiles.setRowCount(len(cfg.profiles))
        for r, (name, p) in enumerate(cfg.profiles.items()):
            mark = L(" (いつも使う)", " (default)") if name == cfg.default_profile else ""
            for c, text in enumerate((name + mark, _kind_label(p.kind), p.description)):
                item = QTableWidgetItem(text); item.setToolTip(text); self.profiles.setItem(r, c, item)
        if not cfg.profiles:
            self.profiles.setRowCount(1)
            self.profiles.setItem(0, 0, QTableWidgetItem(L("(実行先がありません)", "(no targets)")))
        self.profiles.resizeRowsToContents()
        rows = sum(self.profiles.rowHeight(r) for r in range(self.profiles.rowCount()))
        self.profiles.setFixedHeight(self.profiles.horizontalHeader().height() + rows + 2 * self.profiles.frameWidth() + 2)

    def _fields_to_text(self) -> None:
        text = self.editor.toPlainText()
        values = self.field_values()
        new = set_top_level(text, values)
        try:
            ok = all(tomllib.loads(new).get(k) == v for k, v in values.items())
        except tomllib.TOMLDecodeError:
            ok = False
        if not ok:
            data = tomllib.loads(text) if describe_problem(text) is None else {}
            data.update(values)
            new = tomli_w.dumps(Config.model_validate(data).model_dump(mode="json"))
        if new != text:
            self.editor.setPlainText(new)

    def _on_tab(self, index: int) -> None:
        if index == self._tab:
            return
        if index == 1:
            self._fields_to_text(); self._tab = 1
            self.raw_note.setObjectName("hint"); self.raw_note.setText(""); self._restyle(self.raw_note)
            return
        problem = describe_problem(self.editor.toPlainText())
        if problem:
            self._show_problem(L("「基本」タブに切り替えられません", "Cannot switch to the Basic tab"),
                               problem + L("\n\n「詳細」タブで直してから切り替えてください。", "\n\nFix it in the Advanced tab first."))
            self.tabs.blockSignals(True); self.tabs.setCurrentIndex(1); self.tabs.blockSignals(False)
            return
        self._text_to_fields(); self._tab = 0

    def _restyle(self, w: QWidget) -> None:
        w.style().unpolish(w); w.style().polish(w)

    def _show_problem(self, title: str, text: str) -> None:
        self.last_error = text
        self.raw_note.setObjectName("status_ng"); self.raw_note.setText(text); self._restyle(self.raw_note)
        QMessageBox.critical(self, title, text)

    def _save(self) -> None:
        if self.tabs.currentIndex() == 0:
            self._fields_to_text()
        text = self.editor.toPlainText()
        problem = describe_problem(text)
        if problem:
            self._show_problem(L("保存できません", "Cannot save"), problem + L("\n\n環境設定ファイルは書き換えていません。", "\n\nThe settings file was not changed."))
            return
        tmp = self.path.with_suffix(".toml.tmp")
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
            tmp.write_text(text, encoding="utf-8")
            load_config(tmp)
        except (ConfigError, OSError) as ex:
            tmp.unlink(missing_ok=True)
            msg = str(ex).replace(str(tmp), L("環境設定ファイル", "the settings file"))
            self._show_problem(L("保存できません", "Cannot save"), msg); return
        try:
            tmp.replace(self.path)
        except OSError as ex:
            tmp.unlink(missing_ok=True)
            self._show_problem(L("保存できません", "Cannot save"), str(ex)); return
        self.accept()
