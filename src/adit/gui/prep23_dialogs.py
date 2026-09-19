
from __future__ import annotations

from pathlib import Path
from typing import Callable

from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDialog, QFileDialog, QFormLayout, QGridLayout, QHBoxLayout,
                               QHeaderView, QLabel, QLineEdit, QMessageBox, QPushButton, QTableWidget, QVBoxLayout, QWidget)

from adit.gui.style import GROUP_SPACING, PANEL_MARGIN, ROW_SPACING
from adit.gui.widgets import add_row, narrow
from adit.lang import L
from adit.web import prep23 as P

RX_ROWS = 3
RX_COLUMNS = ("name", "file", "nu", "charge", "mult")


def _hint(text: str) -> QLabel:
    w = QLabel(text); w.setObjectName("hint"); w.setWordWrap(True)
    return w


def _combo(items: list[tuple[str, str]], current: str = "") -> QComboBox:
    c = QComboBox()
    for k, v in items:
        c.addItem(v, k)
    c.setCurrentIndex(max(0, c.findData(current)))
    return c


def _file_row(edit: QLineEdit, title: str, parent: QWidget) -> QHBoxLayout:
    btn = QPushButton(L("参照…", "Browse…"))

    def choose() -> None:
        path, _ = QFileDialog.getOpenFileName(parent, title, edit.text().strip() or str(Path.home()), "*")
        if path:
            edit.setText(path)
    btn.clicked.connect(choose)
    row = QHBoxLayout(); row.addWidget(edit, 1); row.addWidget(btn)
    edit.browse = btn
    return row


def _folder_row(edit: QLineEdit, title: str, parent: QWidget) -> QHBoxLayout:
    btn = QPushButton(L("参照…", "Browse…"))

    def choose() -> None:
        start = edit.text().strip() or str(Path.home())
        d = QFileDialog.getExistingDirectory(parent, title, str(Path(start).expanduser().parent))
        if d:
            edit.setText(d)
    btn.clicked.connect(choose)
    row = QHBoxLayout(); row.addWidget(edit, 1); row.addWidget(btn)
    edit.browse = btn
    return row


class BatchDialog(QDialog):

    kind = ""

    def __init__(self, make_spec: Callable, cfg, default_dir: str, parent: QWidget | None = None):
        super().__init__(parent)
        self.make_spec, self.cfg = make_spec, cfg
        self.out_dir: Path | None = None
        self.dirs: list[Path] = []
        self.result = None
        self._rows: dict = {}
        self.setWindowTitle(P.kind_title(self.kind))
        self.form = QFormLayout(); self.form.setVerticalSpacing(ROW_SPACING)
        self.note = _hint(P.kind_note(self.kind))
        self.folder = QLineEdit(default_dir)
        self.cost = _hint(L("ディレクトリの数だけ計算が増えます。大きな系は、この PC ではなくクラスタで実行してください。",
                            "Each directory is a separate calculation. Run large systems on a cluster rather than on this PC."))
        self.btn_cancel = QPushButton(L("キャンセル", "Cancel"))
        self.btn_ok = QPushButton(L("生成", "Generate")); self.btn_ok.setObjectName("primary"); self.btn_ok.setDefault(True)
        self.build()
        add_row(self.form, P.LABELS["out"][0], _folder_row(self.folder, P.lab("out"), self),
                help_text=L("この中にディレクトリを作ります。空でないときは上書きの確認をします。",
                            "The directories are created inside it; a non-empty directory asks before overwriting."))
        buttons = QHBoxLayout(); buttons.addStretch(); buttons.addWidget(self.btn_cancel); buttons.addWidget(self.btn_ok)
        lay = QVBoxLayout(self); lay.setContentsMargins(PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN, PANEL_MARGIN); lay.setSpacing(GROUP_SPACING)
        lay.addWidget(self.note); lay.addLayout(self.form); lay.addWidget(self.cost); lay.addLayout(buttons)
        self.setMinimumWidth(760)
        self.btn_cancel.clicked.connect(self.reject)
        self.btn_ok.clicked.connect(self.generate)

    def build(self) -> None:
        raise NotImplementedError

    def values(self) -> dict[str, str]:
        raise NotImplementedError

    def row(self, key: str, widget, *, help_text: str | None = None, required: bool | None = None):
        add_row(self.form, P.LABELS[key][0], widget, required=required, help_text=help_text)
        return widget

    def line(self, key: str, placeholder: str = "", *, browse: str = "", required: bool | None = None) -> QLineEdit:
        edit = QLineEdit()
        if placeholder:
            edit.setPlaceholderText(placeholder)
        field = _file_row(edit, browse, self) if browse else edit
        self.row(key, field, required=required)
        if field is not edit:
            self._rows[edit] = field
        return edit

    def set_row_visible(self, anchor, visible: bool) -> None:
        self.form.setRowVisible(self._rows.get(anchor, anchor), visible)

    def generate(self) -> None:
        from adit.config import ConfigError
        from adit.project import ProjectError

        if not self.folder.text().strip():
            QMessageBox.warning(self, L("生成できません", "Cannot generate"), L("保存先を指定してください。", "Choose an output directory."))
            return
        out = Path(self.folder.text().strip()).expanduser()
        try:
            spec = self.make_spec()
            overwrite = False
            if out.exists() and any(out.iterdir()):
                ans = QMessageBox.question(self, L("上書きの確認", "Overwrite?"),
                                           L(f"{out} は空ではありません。中のファイルを上書きしますか?", f"{out} is not empty. Overwrite the files inside?"))
                if ans != QMessageBox.StandardButton.Yes:
                    return
                overwrite = True
            res = P.run_batch(self.kind, spec, self.cfg, out, self.values(), overwrite=overwrite)
        except ImportError as ex:
            QMessageBox.critical(self, L("生成できません", "Cannot generate"),
                                 L(f"必要なパッケージがありません: {ex}", f"a required package is missing: {ex}"))
            return
        except (ProjectError, ConfigError, ValueError, OSError) as ex:
            from pydantic import ValidationError as PydanticError
            text = str(ex)
            if isinstance(ex, PydanticError):
                from adit.validate_types import friendly_pydantic
                text = friendly_pydantic(ex)
            QMessageBox.critical(self, L("生成できません", "Cannot generate"), text)
            return
        self.out_dir, self.dirs, self.result = out, list(res.dirs), res
        QMessageBox.information(self, L("生成しました", "Generated"), P.result_message(res))
        self.accept()


class CompareSetDialog(BatchDialog):
    kind = "compare"

    def build(self) -> None:
        self.kind_box = self.row("cmp_kind", _combo(P.compare_kind_names()), required=True)
        self.slab = self.line("cmp_slab", L("スラブの構造のファイル (空欄なら画面の構造)", "structure file of the slab (empty = the structure on screen)"),
                              browse=P.lab("cmp_slab"))
        self.slab_fixed = self.line("cmp_slab_fixed", L("1 始まりの番号 (例 1-4,7)。空欄なら固定しません", "1-based numbers (e.g. 1-4,7); empty = none fixed"))
        self.mol = self.line("cmp_mol", L("分子の構造のファイル (空欄なら画面の構造)", "structure file of the molecule (empty = the structure on screen)"),
                             browse=P.lab("cmp_mol"))
        self.box = self.row("cmp_box", _combo(P.box_names(), "as_is"))
        self.ads = self.line("cmp_ads", L("吸着した構造のファイル (空欄なら画面の構造)", "structure file of the adsorbed state (empty = the structure on screen)"),
                             browse=P.lab("cmp_ads"))
        self.ads_fixed = self.line("cmp_ads_fixed", L("1 始まりの番号 (例 1-4,7)", "1-based numbers (e.g. 1-4,7)"))
        self.use_slab = QPushButton(L("スラブに使う", "Use as slab"))
        self.use_mol = QPushButton(L("分子に使う", "Use as molecule"))
        self.use_ads = QPushButton(L("吸着構造に使う", "Use as adsorbed state"))
        self.save_structure = QPushButton(L("extended XYZ で保存…", "Save as extended XYZ…"))
        uses = QGridLayout()
        for i, button in enumerate((self.use_slab, self.use_mol, self.use_ads, self.save_structure)):
            button.setObjectName("link"); uses.addWidget(button, i // 2, i % 2)
        self.structure_tools = self.row("cmp_structure_tools", uses)
        self.use_slab.clicked.connect(lambda: self.slab.clear())
        self.use_mol.clicked.connect(lambda: self.mol.clear())
        self.use_ads.clicked.connect(lambda: self.ads.clear())
        self.save_structure.clicked.connect(self._save_structure)
        self.table = QTableWidget(0, len(RX_COLUMNS))
        self.table.setHorizontalHeaderLabels([L("名前 (ディレクトリ)", "Name (directory)"), L("構造のファイル", "Structure file"),
                                              L("係数 ν", "Coefficient ν"), L("全電荷", "Total charge"), L("スピン多重度", "Spin multiplicity")])
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        self.table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.table.setMinimumHeight(150)
        self.btn_add = QPushButton(L("行を追加", "Add row")); self.btn_add.setObjectName("link")
        self.btn_add.clicked.connect(lambda: self.add_row())
        rx = QWidget(); rx.setObjectName("rowbox"); rxl = QVBoxLayout(rx); rxl.setContentsMargins(0, 0, 0, 0); rxl.setSpacing(4)
        rxl.addWidget(self.table)
        rxl.addWidget(_hint(L("生成物は正の係数、反応物は負の係数。構造のファイルが空欄の行は画面の構造を使います。",
                              "Products have a positive coefficient, reactants a negative one; an empty structure file means the structure on screen.")))
        brow = QHBoxLayout(); brow.addWidget(self.btn_add); brow.addStretch(1); rxl.addLayout(brow)
        self.rx_box = self.row("cmp_rows", rx, required=True)
        for _ in range(RX_ROWS):
            self.add_row()
        self.solvent = self.line("cmp_solvent", '{"solvation": "alpb", "solvent": "water"}', required=True)
        self.gas = self.line("cmp_gas", L('空欄なら画面の条件のまま (例 {"solvation": "none"})', 'empty = the settings on screen (e.g. {"solvation": "none"})'))
        self.set_file = self.line("cmp_set_file", L("組の定義を書いた JSON を使うなら、そのパス (上の欄は使いません)",
                                                    "path of a JSON set definition, if you use one (the fields above are then unused)"),
                                  browse=P.lab("cmp_set_file"))
        self.kind_box.currentIndexChanged.connect(self._on_kind)
        self.set_file.textChanged.connect(self._on_kind)
        self._on_kind()

    def add_row(self, values: tuple[str, ...] = ()) -> None:
        r = self.table.rowCount(); self.table.insertRow(r)
        for c, col in enumerate(RX_COLUMNS):
            w = QLineEdit(values[c] if c < len(values) else "")
            w.setAccessibleName(f"{self.table.horizontalHeaderItem(c).text()} {r + 1}")
            if col == "nu":
                w.setPlaceholderText(L("例 -1", "e.g. -1"))
            elif col == "name":
                w.setPlaceholderText(L("例 h2o", "e.g. h2o"))
            self.table.setCellWidget(r, c, w)

    def _save_structure(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, L("画面の構造を保存", "Save the screen structure"), "structure.extxyz",
                                              L("extended XYZ (*.extxyz *.xyz)", "extended XYZ (*.extxyz *.xyz)"))
        if not path:
            return
        try:
            from adit.convert import write_spec_structure
            written = write_spec_structure(self.make_spec(), path)
        except (ValueError, OSError) as ex:
            QMessageBox.critical(self, L("保存できません", "Cannot save"), str(ex)); return
        QMessageBox.information(self, L("保存しました", "Saved"), L(f"画面の構造を保存しました: {written}", f"Saved the screen structure: {written}"))

    def _on_kind(self, *_) -> None:
        kind = self.kind_box.currentData()
        from_file = bool(self.set_file.text().strip())
        for w, kinds in ((self.slab, ("adsorption",)), (self.slab_fixed, ("adsorption",)), (self.mol, ("adsorption",)),
                         (self.box, ("adsorption",)), (self.ads, ("adsorption",)), (self.ads_fixed, ("adsorption",)),
                         (self.structure_tools, ("adsorption", "reaction")),
                         (self.rx_box, ("reaction",)), (self.solvent, ("solvation",)), (self.gas, ("solvation",))):
            self.set_row_visible(w, kind in kinds and not from_file)
        self.set_row_visible(self.kind_box, not from_file)

    def values(self) -> dict[str, str]:
        v = {"cmp_kind": self.kind_box.currentData(), "cmp_box": self.box.currentData(),
             "cmp_slab": self.slab.text(), "cmp_slab_fixed": self.slab_fixed.text(), "cmp_mol": self.mol.text(),
             "cmp_ads": self.ads.text(), "cmp_ads_fixed": self.ads_fixed.text(),
             "cmp_solvent": self.solvent.text(), "cmp_gas": self.gas.text(), "cmp_set_file": self.set_file.text(),
             "cmp_rx_rows": str(self.table.rowCount())}
        for r in range(self.table.rowCount()):
            for c, col in enumerate(RX_COLUMNS):
                v[f"cmp_rx{r + 1}_{col}"] = self.table.cellWidget(r, c).text()
        return v


class ConformerDialog(BatchDialog):
    kind = "conformers"

    def build(self) -> None:
        self.n = self.line("conf_n", L("例 50 (RDKit の ETKDG で作る数)", "e.g. 50 (embedded with RDKit ETKDG)"), required=True)
        self.n.setMaximumWidth(220)
        self.rmsd = self.line("conf_rmsd", L("これ以下の RMSD の配座は重複として外します (既定値はありません)",
                                             "conformers closer than this RMSD count as duplicates (there is no default)"), required=True)
        self.seed = self.line("conf_seed", "12345")
        self.ff = self.row("conf_ff", _combo(P.ff_names(), "MMFF94"))
        self.iters = self.line("conf_iters", "200")
        self.all_atoms = QCheckBox(L("水素も含めた全原子で RMSD を測る (既定は水素を除く)", "measure the RMSD over all atoms including hydrogens (default: without)"))
        self.row("conf_all", self.all_atoms)
        self.keep = self.line("conf_keep", L("空欄なら残った配座を全部使います", "empty = keep all remaining conformers"))
        self.smiles = self.line("conf_smiles", L("画面の構造が SMILES 由来でないときに書きます", "write this if the structure on screen does not come from SMILES"))
        for w in (self.seed, self.iters, self.keep):
            narrow(w)
        self.note.setText(self.note.text() + "\n" + L("乱数の種を 0 にすると、RDKit 2026.03 では配座がすべて同じ座標になりました (別の種にしてください)。",
                                                      "With seed 0, RDKit 2026.03 embedded identical coordinates for every conformer; use another seed."))

    def values(self) -> dict[str, str]:
        return {"conf_n": self.n.text(), "conf_rmsd": self.rmsd.text(), "conf_seed": self.seed.text(),
                "conf_ff": self.ff.currentData(), "conf_iters": self.iters.text(), "conf_all": "on" if self.all_atoms.isChecked() else "",
                "conf_keep": self.keep.text(), "conf_smiles": self.smiles.text()}


class NebDialog(BatchDialog):
    kind = "neb"

    def build(self) -> None:
        self.start = self.line("neb_start", L("空欄なら画面の構造を始状態にします", "empty = the structure on screen is the initial state"),
                               browse=P.lab("neb_start"))
        self.end = self.line("neb_end", L("終状態の構造のファイル (始状態と同じ原子の並び)", "structure file of the final state (same atom order as the initial one)"),
                             browse=P.lab("neb_end"), required=True)
        self.images = self.line("neb_images", L("例 5 (始状態と終状態を除く数)", "e.g. 5 (not counting the end points)"), required=True)
        narrow(self.images)
        self.mode = self.row("neb_mode", _combo(P.neb_mode_names(), "native"))
        self.interp = self.row("neb_interp", _combo(P.interp_names(), "idpp"))
        self.climb = QCheckBox(L("climbing image を使う (QE と CP2K。VASP 本体にはありません)", "use the climbing image (QE and CP2K; not in plain VASP)"))
        self.row("neb_climb", self.climb)
        self.spring = narrow(self.line("neb_spring", L("空欄なら VASP の既定", "empty = the VASP default")))
        self.kspring = narrow(self.line("neb_kspring", L("空欄なら CP2K の既定", "empty = the CP2K default")))
        self.opt = self.row("neb_opt", _combo(P.qe_opt_names(), ""))

    def values(self) -> dict[str, str]:
        return {"neb_start": self.start.text(), "neb_end": self.end.text(), "neb_images": self.images.text(),
                "neb_mode": self.mode.currentData(), "neb_interp": self.interp.currentData(),
                "neb_climb": "on" if self.climb.isChecked() else "", "neb_spring": self.spring.text(),
                "neb_kspring": self.kspring.text(), "neb_opt": self.opt.currentData() or ""}


class PhononDialog(BatchDialog):
    kind = "phonons"

    def build(self) -> None:
        self.dim = narrow(self.line("ph_dim", L("例 2x2x2", "e.g. 2x2x2"), required=True))
        self.disp = narrow(self.line("ph_disp", L("空欄なら phonopy / ASE の既定", "empty = the phonopy / ASE default")))
        self.backend = self.row("ph_backend", _combo(P.backend_names(), "auto"))
        self.dos_mesh = narrow(self.line("ph_dos_mesh", L("空欄なら状態密度を書きません (例 10x10x10)", "empty = no DOS is written (e.g. 10x10x10)")))
        self.dos_width = narrow(self.line("ph_dos_width", L("ASE で集めるときは必須 (phonopy では使いません)", "required for the ASE backend (unused with phonopy)")))

    def values(self) -> dict[str, str]:
        return {"ph_dim": self.dim.text(), "ph_disp": self.disp.text(), "ph_backend": self.backend.currentData(),
                "ph_dos_mesh": self.dos_mesh.text(), "ph_dos_width": self.dos_width.text()}


class ElasticDialog(BatchDialog):
    kind = "elastic"

    def build(self) -> None:
        self.strains = self.line("el_strains", L("カンマで区切って (例 -0.01,-0.005,0.005,0.01)。歪みなしの e0 は自動で作ります",
                                                 "comma separated (e.g. -0.01,-0.005,0.005,0.01); the unstrained e0 is added automatically"), required=True)
        box = QWidget(); box.setObjectName("rowbox"); row = QHBoxLayout(box); row.setContentsMargins(0, 0, 0, 0)
        self.comps = {}
        for j, name in enumerate(("1 (xx)", "2 (yy)", "3 (zz)", "4 (yz)", "5 (xz)", "6 (xy)"), start=1):
            c = QCheckBox(name); c.setChecked(True); row.addWidget(c); self.comps[j] = c
        row.addStretch(1)
        self.row("el_comps", box)

    def values(self) -> dict[str, str]:
        return {"el_strains": self.strains.text(), "el_comps": ",".join(str(j) for j, c in self.comps.items() if c.isChecked())}


class TsDialog(BatchDialog):
    kind = "ts"

    def __init__(self, make_spec: Callable, cfg, default_dir: str, code: str = "", parent: QWidget | None = None):
        self.code = code
        super().__init__(make_spec, cfg, default_dir, parent)

    def build(self) -> None:
        default = "sella" if self.code in ("xtb", "mlip") else "ts"
        self.mode = self.row("ts_mode", _combo(P.ts_mode_names(), default), required=True)
        self.irc = QCheckBox(L("鞍点を見つけたあと、IRC を前向き・後ろ向きにたどる", "after the saddle point, follow the IRC forward and backward"))
        self.irc.setChecked(True)
        self.row("ts_irc", self.irc)
        self.mode.currentIndexChanged.connect(self._on_mode)
        self._on_mode()

    def _on_mode(self, *_) -> None:
        self.form.setRowVisible(self.irc, self.mode.currentData() == "sella")

    def values(self) -> dict[str, str]:
        return {"ts_mode": self.mode.currentData(), "ts_irc": "on" if self.irc.isChecked() else ""}


DIALOGS = {"compare": CompareSetDialog, "conformers": ConformerDialog, "neb": NebDialog, "phonons": PhononDialog,
           "elastic": ElasticDialog, "ts": TsDialog}

__all__ = ["BatchDialog", "CompareSetDialog", "ConformerDialog", "NebDialog", "PhononDialog", "ElasticDialog", "TsDialog", "DIALOGS"]
