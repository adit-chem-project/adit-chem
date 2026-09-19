
from __future__ import annotations

import subprocess
import sys
import tempfile
import threading
from pathlib import Path

import numpy as np
from ase.io import write
from PySide6.QtCore import Qt, Signal
from ase.data import chemical_symbols
from PySide6.QtWidgets import (QCheckBox, QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
                               QPushButton, QSpinBox, QVBoxLayout, QWidget)

from adit.builder.model import Base, ClusterRef, PolymerRef, Recipe, TwoD
from adit.lang import L
from adit.gui.style import ROW_SPACING
from adit.gui.style import LABEL_WIDTH  # noqa: F401  
from adit.gui import icons
from adit.gui.i18n import tr
from adit.gui.panels.mixture_editor import MixtureEditor
from adit.gui.panels.recipe_editor import Num, RecipeEditor, hint, hrow, spin
from adit.gui.widgets import SciDoubleSpinBox, add_row, label, limit_combo, narrow
from adit.spec import AtomsData, Structure
from adit.textparse import parse_constraints, parse_indices  # noqa: F401  
from adit.structure import (CRYSTAL_STRUCTURES, SURFACE_FUNCTIONS, StructureError, build_structure, default_bulk,  # noqa: F401
                             has_rdkit, preset_matches, preset_names, pretty_formula)

ELEMENTS = [s for s in chemical_symbols[1:104]]


SOURCES = {"preset": "プリセット", "smiles": "SMILES", "file": "ファイル", "bulk": "バルク", "surface": "スラブ", "mixture": "溶液・混合物",
           "2d": "2 次元材料・ナノチューブ", "cluster": "ナノ粒子", "polymer": "ポリマー"}
NEW_BASES = ("2d", "cluster", "polymer")
SLOW_BASES = ("polymer",)
TWOD_KINDS = {"graphene": "グラフェン型 (C₂、BN)", "mx2": "MX₂ 型 (MoS₂ など)", "nanoribbon": "ナノリボン", "nanotube": "ナノチューブ"}
CLUSTER_KINDS = {"icosahedron": "正二十面体", "decahedron": "十面体", "octahedron": "八面体", "wulff": "Wulff 形 (面のエネルギーから)"}


def _compact(model, extra: dict) -> dict:
    d = model.model_dump(mode="json", exclude_defaults=True)
    full = model.model_dump(mode="json")
    for k in extra:
        if k not in d and k in full:
            d[k] = full[k]
    return d


class StructurePanel(QGroupBox):
    changed = Signal()
    built = Signal()
    _build_done = Signal(int, object, object, str)

    def __init__(self, parent: QWidget | None = None):
        super().__init__("構造", parent)
        self.setObjectName("structure")
        self._structure: Structure | None = None
        self._error: str = ""
        self._restoring = False
        self._raw: Structure | None = None
        self._built_ref = ""
        self._token = 0
        self._building = False
        self._pending_ref = ""
        self._auto_charge: int | None = None
        self._extra: dict[str, dict] = {}
        self._file_restored: tuple[str, str] | None = None
        self._sha_cache: dict[tuple, str] = {}
        self._before_cache: dict[str, object] = {}

        self.source = QComboBox()
        self.source.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
        for key, text in SOURCES.items():
            self.source.addItem(icons.source_icon(key), text, key)
        self.preset = QComboBox()
        self._preset_names = preset_names()
        for name in self._preset_names:
            self.preset.addItem(pretty_formula(name), name)
        self.preset.setCurrentIndex(max(0, self.preset.findData("H2O")))
        self.preset_search = QLineEdit(); self.preset_search.setPlaceholderText(L("式・英語名・日本語名で絞り込み", "Filter by formula, English name or Japanese name"))
        self.preset_search.setClearButtonEnabled(True)
        self.smiles = QLineEdit()
        self.smiles.setPlaceholderText("例 CCO")
        self.smiles_go = QPushButton("3D 構造を生成")
        self.file = QLineEdit()
        self.file_browse = QPushButton("参照…")
        self.mixture = MixtureEditor()
        self.smiles_draw = QPushButton("Draw"); self.smiles_draw.setToolTip("Draw: 分子を描いて SMILES にします (RDKit が要ります)")
        self.bulk_el = QComboBox(); self.bulk_el.addItems(ELEMENTS); self.bulk_el.setCurrentText("Si")
        self.bulk_struct = QComboBox(); self.bulk_struct.addItems(["(元素の既定値)"] + CRYSTAL_STRUCTURES)
        self.bulk_a = narrow(SciDoubleSpinBox(0.0, 100.0, 0.0, 0.1)); self.bulk_a.setToolTip("0 なら元素の既定値 (ASE の参照状態)")
        self.bulk_cubic = QCheckBox("立方晶セル")
        self.surf_facet = QComboBox(); self.surf_facet.addItems(SURFACE_FUNCTIONS); self.surf_facet.setCurrentText("fcc111")
        self.surf_el = QComboBox(); self.surf_el.addItems(ELEMENTS); self.surf_el.setCurrentText("Al")
        self.surf_n = [narrow(QSpinBox()) for _ in range(3)]
        for w in self.surf_n:
            w.setMaximumWidth(64)
        for w, v in zip(self.surf_n, (2, 2, 3)):
            w.setRange(1, 50); w.setValue(v); w.setMaximumWidth(60)
        self.surf_vac = narrow(SciDoubleSpinBox(0.0, 100.0, 10.0, 1.0))
        self._build_new_base_widgets()
        self.box = QCheckBox("周期セルに入れる"); self.box.setToolTip("分子を立方体の周期セルに置きます。VASP と pw.x では周期セルが必須です")
        self.box_size = narrow(SciDoubleSpinBox(1.0, 500.0, 15.0, 1.0))
        self._fixed_placeholder = "1 始まりの番号 (例 1-4,7)。軸を指定するなら 7:xy (x, y を固定)。空欄なら全原子を動かします"
        self.fixed = QLineEdit(); self.fixed.setPlaceholderText(self._fixed_placeholder)
        self.charge = narrow(QSpinBox()); self.charge.setRange(-20, 20)
        self.multiplicity = narrow(QSpinBox()); self.multiplicity.setRange(1, 20)
        self.info = QLabel("")
        self.asegui = QPushButton("ASE GUI で開く"); self.asegui.setVisible(False)

        if not has_rdkit():
            for w in (self.smiles, self.smiles_go):
                w.setEnabled(False)
            for key, text in (("smiles", "SMILES (RDKit が未導入のため使えません)"), ("polymer", L("ポリマー (RDKit が未導入のため使えません)", "Polymer (RDKit not installed)"))):
                i = self.source.findData(key)
                self.source.setItemText(i, text)
                item = self.source.model().item(i)
                if item is not None:
                    item.setEnabled(False)

        form = QFormLayout(); self._form = form
        form.setVerticalSpacing(ROW_SPACING)
        self.draw_button = QPushButton("Draw")
        self.draw_button.setToolTip("Draw: 分子を描いて SMILES にします (RDKit が要ります)")
        src_row = QWidget(); src_row.setObjectName("rowbox"); sl = QHBoxLayout(src_row); sl.setContentsMargins(0, 0, 0, 0)
        sl.addWidget(self.source, 1); sl.addWidget(self.draw_button)
        add_row(form, "構造の作り方", src_row)
        self._source_rows: dict[str, list[QWidget]] = {}
        add_row(form, L("プリセットを検索", "Search presets"), self.preset_search)
        add_row(form, "プリセット", self.preset); self._source_rows["preset"] = [self.preset_search, self.preset]
        row_w = QWidget(); row_w.setObjectName("rowbox"); row = QHBoxLayout(row_w); row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self.smiles); row.addWidget(self.smiles_go); row.addWidget(self.smiles_draw)
        add_row(form, "SMILES", row_w); self._source_rows["smiles"] = [row_w]
        row_w = QWidget(); row_w.setObjectName("rowbox"); row = QHBoxLayout(row_w); row.setContentsMargins(0, 0, 0, 0); row.addWidget(self.file); row.addWidget(self.file_browse)
        add_row(form, "ファイル", row_w); self._source_rows["file"] = [row_w]
        row_w = QWidget(); row_w.setObjectName("rowbox"); row = QHBoxLayout(row_w); row.setContentsMargins(0, 0, 0, 0)
        for w in (self.bulk_el, self.bulk_struct):
            row.addWidget(w)
        row.addStretch(); add_row(form, "バルク", row_w)
        row_b2 = QWidget(); row_b2.setObjectName("rowbox"); row = QHBoxLayout(row_b2); row.setContentsMargins(0, 0, 0, 0)
        for w in (QLabel("a [Å]"), self.bulk_a, self.bulk_cubic):
            row.addWidget(w)
        row.addStretch(); form.addRow(label(""), row_b2); self._source_rows["bulk"] = [row_w, row_b2]
        row_w = QWidget(); row_w.setObjectName("rowbox"); row = QHBoxLayout(row_w); row.setContentsMargins(0, 0, 0, 0)
        for w in (self.surf_facet, self.surf_el, self.surf_n[0], QLabel("×"), self.surf_n[1]):
            row.addWidget(w)
        row.addStretch(); add_row(form, "スラブ", row_w)
        row_w2 = QWidget(); row_w2.setObjectName("rowbox"); row = QHBoxLayout(row_w2); row.setContentsMargins(0, 0, 0, 0)
        for w in (QLabel("層"), self.surf_n[2], QLabel("真空層 [Å]"), self.surf_vac):
            row.addWidget(w)
        row.addStretch(); form.addRow(label(""), row_w2); self._source_rows["surface"] = [row_w, row_w2]
        form.addRow(self.mixture); self._source_rows["mixture"] = [self.mixture]
        self._add_new_base_rows(form)
        sub = QLabel("組み立て手順"); sub.setObjectName("subtitle")
        sub.setToolTip(L("土台の構造に加工 (面で切る、溶液の層、原子を抜く・置換、吸着、固定、超格子、真空、直方体、溶媒和) を上から順に行います",
                         "processes the base structure from top to bottom: cut a slab, add a solvent layer, remove or substitute atoms, adsorb, fix, supercell, vacuum, box, solvate"))
        form.addRow(sub)
        self.recipe = RecipeEditor(); form.addRow(self.recipe)
        sub = QLabel("共通設定"); sub.setObjectName("subtitle"); form.addRow(sub)
        row_w = QWidget(); row_w.setObjectName("rowbox"); row = QHBoxLayout(row_w); row.setContentsMargins(0, 0, 0, 0)
        row.addWidget(self.box); row.addWidget(QLabel("一辺 [Å]")); row.addWidget(self.box_size); row.addStretch()
        add_row(form, "周期セルに入れる", row_w); self._box_row = row_w
        self.box.setText("")
        add_row(form, "固定原子", self.fixed)
        add_row(form, "全電荷", self.charge)
        add_row(form, "スピン多重度", self.multiplicity)
        self.info.setObjectName("hint"); self.info.setWordWrap(True)
        form.addRow(label(""), self.info)
        lay = QVBoxLayout(self); lay.addLayout(form)
        self.source.currentIndexChanged.connect(self._show_source_rows)
        self._show_source_rows()
        self.mixture.changed.connect(self._rebuild)
        self.mixture.draw_requested.connect(self._draw_into_mixture)
        self.smiles_draw.clicked.connect(self._draw_smiles); self.draw_button.clicked.connect(self._draw_smiles)
        if not has_rdkit():
            self.smiles_draw.setEnabled(False); self.mixture.set_draw_enabled(False); self.draw_button.setEnabled(False)
            self.recipe.set_draw_enabled(False)

        self.source.currentIndexChanged.connect(self._rebuild)
        for w in (self.bulk_el, self.bulk_struct, self.surf_facet, self.surf_el):
            w.currentTextChanged.connect(self._rebuild)
        for w in (self.bulk_a, self.surf_vac, self.box_size, *self.surf_n):
            w.valueChanged.connect(self._rebuild)
        self.bulk_cubic.toggled.connect(self._rebuild)
        self.box.toggled.connect(self._rebuild)
        self.fixed.editingFinished.connect(self._rebuild)
        self.preset.currentIndexChanged.connect(self._rebuild)
        self.preset_search.textChanged.connect(self._filter_presets)
        self.smiles_go.clicked.connect(self._rebuild)
        self.smiles.returnPressed.connect(self._rebuild)
        from PySide6.QtCore import QTimer
        self._smiles_timer = QTimer(self); self._smiles_timer.setSingleShot(True); self._smiles_timer.setInterval(500)
        self._smiles_timer.timeout.connect(self._rebuild)
        self.smiles.textChanged.connect(lambda *_: self._smiles_timer.start())
        self.file.editingFinished.connect(self._on_file_edited)
        self.file_browse.clicked.connect(self._browse)
        self.charge.valueChanged.connect(self._rebuild)
        self.multiplicity.valueChanged.connect(self._rebuild)
        self.asegui.clicked.connect(self._open_ase_gui)
        self.recipe.base_is_slab = self._base_is_slab
        self.recipe.area_before = self._area_before
        self.recipe.changed.connect(self._rebuild)
        self.recipe.build_requested.connect(self.build_recipe)
        self.recipe.cancel_requested.connect(self.cancel_build)
        self.recipe.draw_requested.connect(self._draw_into)
        self._build_done.connect(self._on_built, Qt.ConnectionType.QueuedConnection)
        self._rebuild()

    def _build_new_base_widgets(self) -> None:
        self.td_kind = QComboBox()
        for k, v in TWOD_KINDS.items():
            self.td_kind.addItem(v, k)
        self.td_formula = QLineEdit(); self.td_formula.setMaximumWidth(150)
        self.td_nx, self.td_ny = spin(1, 50, 2, 64), spin(1, 50, 2, 64)
        self.td_n, self.td_m, self.td_len = spin(0, 60, 6, 64), spin(0, 60, 0, 64), spin(1, 50, 1, 64)
        self.td_n.setToolTip(L("ナノチューブ: カイラル指数 n / ナノリボン: 幅", "nanotube: chiral index n / nanoribbon: width"))
        self.td_m.setToolTip(L("ナノチューブ: カイラル指数 m / ナノリボン: 長さ (周期の方向の繰り返し)", "nanotube: chiral index m / nanoribbon: length (repeats along the periodic direction)"))
        self.td_ribbon = QComboBox(); self.td_ribbon.addItem(L("ジグザグ端", "zigzag edge"), "zigzag"); self.td_ribbon.addItem(L("アームチェア端", "armchair edge"), "armchair")
        self.td_sat = QCheckBox(L("端を H で終端", "H-terminated edges")); self.td_sat.setChecked(True)
        self.td_vac = narrow(Num(0.1, 200.0, 10.0, 1.0))
        self.td_size_w = hrow(L("繰り返し", "repeat"), self.td_nx, "×", self.td_ny, stretch=False)
        self.td_nm_w = hrow("n", self.td_n, "m", self.td_m, stretch=False)
        self.td_len_w = hrow(L("長さ", "length"), self.td_len, stretch=False)
        self.td_rib_w = hrow(self.td_ribbon, self.td_sat, stretch=False)

        self.cl_kind = QComboBox()
        for k, v in CLUSTER_KINDS.items():
            self.cl_kind.addItem(v, k)
        self.cl_el = QComboBox(); self.cl_el.addItems(ELEMENTS); self.cl_el.setCurrentText("Cu")
        self.cl_shells = spin(1, 15, 3, 64)
        self.cl_p, self.cl_q, self.cl_r = spin(1, 20, 2, 60), spin(1, 20, 2, 60), spin(0, 20, 0, 60)
        self.cl_len, self.cl_cut = spin(1, 30, 5, 64), spin(0, 15, 0, 64)
        self.cl_size = spin(1, 50000, 100, 96)
        self.cl_struct = QComboBox(); self.cl_struct.addItems(["fcc", "bcc", "sc"])
        self.cl_a = narrow(Num(0.0, 20.0, 0.0, 0.05)); self.cl_a.setToolTip(L("0 なら元素の既定値 (ASE の参照状態)", "0 = element default (ASE reference state)"))
        self.cl_shells_w = hrow(L("殻の数", "shells"), self.cl_shells, stretch=False)
        self.cl_pqr_w = hrow("p", self.cl_p, "q", self.cl_q, "r", self.cl_r, stretch=False)
        self.cl_oct_w = hrow(L("稜の原子数", "atoms per edge"), self.cl_len, L("頂点を切る層", "corner cut"), self.cl_cut, stretch=False)
        self.cl_wulff_w = hrow(L("原子数の目安", "target atoms"), self.cl_size, self.cl_struct, stretch=False)
        self.cl_pqr_w.setToolTip(L("ASE の Decahedron の p, q, r ((100) 面の数、稜の長さ、頂点の切り込み)", "p, q, r of ASE's Decahedron ((100) facets, edge length, re-entrance)"))

        self.pl_unit = QLineEdit(); self.pl_unit.setPlaceholderText(L("例: *CC* (連結点の * を 2 つ)", "e.g. *CC* (two * connection points)"))
        self.pl_unit.setToolTip(L("繰り返し単位の SMILES。つなぐ位置に * を 2 つ書きます (例 *CC*、*CC(*)c1ccccc1)。末端は H で終端します",
                                  "SMILES of the repeat unit with two * where the units join (e.g. *CC*, *CC(*)c1ccccc1); the chain ends are capped with H"))
        self.pl_n = spin(1, 500, 10, 80); self.pl_seed = spin(0, 999999, 0, 96)

    def _add_new_base_rows(self, form: QFormLayout) -> None:
        r1 = hrow(self.td_kind, self.td_formula)
        r2 = QWidget(); r2.setObjectName("rowbox"); h = QHBoxLayout(r2); h.setContentsMargins(0, 0, 0, 0)
        for w in (self.td_size_w, self.td_nm_w, self.td_len_w, self.td_rib_w):
            h.addWidget(w)
        h.addStretch()
        r3 = hrow(L("真空 (片側) [Å]", "vacuum (each side) [Å]"), self.td_vac)
        add_row(form, "2 次元材料・ナノチューブ", r1); form.addRow(label(""), r2); form.addRow(label(""), r3)
        self._source_rows["2d"] = [r1, r2, r3]
        c1 = hrow(self.cl_kind, self.cl_el)
        c2 = QWidget(); c2.setObjectName("rowbox"); h = QHBoxLayout(c2); h.setContentsMargins(0, 0, 0, 0)
        for w in (self.cl_shells_w, self.cl_pqr_w, self.cl_oct_w, self.cl_wulff_w):
            h.addWidget(w)
        h.addStretch()
        c3 = hrow(L("格子定数 a [Å]", "lattice constant a [Å]"), self.cl_a)
        add_row(form, "ナノ粒子", c1); form.addRow(label(""), c2); form.addRow(label(""), c3)
        self._source_rows["cluster"] = [c1, c2, c3]
        p2 = hrow(L("重合度", "units"), self.pl_n, L("乱数の種", "seed"), self.pl_seed)
        p3 = hint(L("鎖を 1 本作ります (H を含めて 500 原子まで)。3 次元化に数秒〜20 秒ほどかかるので、「作る」を押して作ります",
                    "builds one chain (up to 500 atoms including H). Making the 3D structure takes a few to about 20 seconds, so press Build"))
        add_row(form, "ポリマー", self.pl_unit); form.addRow(label(""), p2); form.addRow(label(""), p3)
        self._source_rows["polymer"] = [self.pl_unit, p2, p3]
        for w in (self.td_kind, self.td_ribbon, self.cl_kind, self.cl_el, self.cl_struct):
            limit_combo(w)
            w.currentIndexChanged.connect(self._rebuild)
        self.td_kind.currentIndexChanged.connect(self._show_new_base_parts); self.cl_kind.currentIndexChanged.connect(self._show_new_base_parts)
        for w in (self.td_nx, self.td_ny, self.td_n, self.td_m, self.td_len, self.td_vac, self.cl_shells, self.cl_p, self.cl_q, self.cl_r,
                  self.cl_len, self.cl_cut, self.cl_size, self.cl_a, self.pl_n, self.pl_seed):
            w.valueChanged.connect(self._rebuild)
        self.td_sat.toggled.connect(self._rebuild)
        self.td_formula.editingFinished.connect(self._rebuild); self.pl_unit.editingFinished.connect(self._rebuild)
        self._show_new_base_parts()

    def _show_new_base_parts(self, *_) -> None:
        k = self.td_kind.currentData()
        self.td_size_w.setVisible(k in ("graphene", "mx2")); self.td_nm_w.setVisible(k in ("nanotube", "nanoribbon"))
        self.td_len_w.setVisible(k == "nanotube"); self.td_rib_w.setVisible(k == "nanoribbon")
        self.td_formula.setVisible(k != "nanoribbon")
        self.td_formula.setPlaceholderText({"graphene": L("C2 (空欄なら C2。BN も可)", "C2 (empty = C2; BN also works)"), "mx2": L("MoS2 (空欄なら MoS2)", "MoS2 (empty = MoS2)"),
                                            "nanotube": L("元素 (空欄なら C)", "element (empty = C)")}.get(k, ""))
        c = self.cl_kind.currentData()
        self.cl_shells_w.setVisible(c == "icosahedron"); self.cl_pqr_w.setVisible(c == "decahedron")
        self.cl_oct_w.setVisible(c == "octahedron"); self.cl_wulff_w.setVisible(c == "wulff")

    def _twod_ref(self) -> dict:
        k = self.td_kind.currentData()
        d = dict(self._extra.get("2d", {})); d["kind"] = k; d["vacuum"] = self.td_vac.value()
        text = self.td_formula.text().strip()
        if k in ("graphene", "mx2"):
            third = (d.get("size") or (1, 1, 1))[2]
            d["formula"] = text; d["size"] = [self.td_nx.value(), self.td_ny.value(), int(third)]
        elif k == "nanotube":
            d["n"], d["m"], d["length"], d["symbol"] = self.td_n.value(), self.td_m.value(), self.td_len.value(), text or "C"
        else:
            d["n"], d["m"], d["ribbon_type"], d["saturated"] = self.td_n.value(), self.td_m.value(), self.td_ribbon.currentData(), self.td_sat.isChecked()
        return _compact(TwoD.model_validate(d), self._extra.get("2d", {}))

    def _cluster_ref(self) -> dict:
        k = self.cl_kind.currentData()
        d = dict(self._extra.get("cluster", {})); d["kind"] = k; d["symbol"] = self.cl_el.currentText()
        d["lattice_constant"] = self.cl_a.value() or None
        if k == "icosahedron":
            d["shells"] = self.cl_shells.value()
        elif k == "decahedron":
            d["p"], d["q"], d["r"] = self.cl_p.value(), self.cl_q.value(), self.cl_r.value()
        elif k == "octahedron":
            d["length"], d["cutoff"] = self.cl_len.value(), self.cl_cut.value()
        else:
            d["size"], d["structure"] = self.cl_size.value(), self.cl_struct.currentText()
        return _compact(ClusterRef.model_validate(d), self._extra.get("cluster", {}))

    def _polymer_ref(self) -> dict:
        unit = self.pl_unit.text().strip()
        if not unit:
            raise StructureError(L("繰り返し単位の SMILES を入力してください (例 *CC*)", "enter the SMILES of the repeat unit (e.g. *CC*)"))
        d = dict(self._extra.get("polymer", {})); d.update(unit=unit, n=self.pl_n.value(), seed=self.pl_seed.value())
        return _compact(PolymerRef.model_validate(d), self._extra.get("polymer", {}))

    def _load_new_base(self, src: str, ref: dict) -> None:
        if src == "2d":
            t = TwoD.model_validate(ref)
            self.td_kind.setCurrentIndex(max(0, self.td_kind.findData(t.kind)))
            self.td_formula.setText(t.symbol if t.kind == "nanotube" and t.symbol != "C" else (t.formula if t.kind != "nanotube" else ""))
            self.td_nx.setValue(t.size[0]); self.td_ny.setValue(t.size[1])
            self.td_n.setValue(t.n); self.td_m.setValue(t.m); self.td_len.setValue(t.length)
            self.td_ribbon.setCurrentIndex(max(0, self.td_ribbon.findData(t.ribbon_type))); self.td_sat.setChecked(t.saturated)
            self.td_vac.setValue(t.vacuum)
        elif src == "cluster":
            c = ClusterRef.model_validate(ref)
            self.cl_kind.setCurrentIndex(max(0, self.cl_kind.findData(c.kind))); self.cl_el.setCurrentText(c.symbol)
            self.cl_shells.setValue(c.shells); self.cl_p.setValue(c.p); self.cl_q.setValue(c.q); self.cl_r.setValue(c.r)
            self.cl_len.setValue(c.length); self.cl_cut.setValue(c.cutoff); self.cl_size.setValue(c.size); self.cl_struct.setCurrentText(c.structure)
            self.cl_a.setValue(c.lattice_constant or 0.0)
        elif src == "polymer":
            p = PolymerRef.model_validate(ref)
            self.pl_unit.setText(p.unit); self.pl_n.setValue(p.n); self.pl_seed.setValue(p.seed)
        self._extra[src] = dict(ref)
        self._show_new_base_parts()

    def _show_source_rows(self, *_) -> None:
        cur = self.current_source()
        for src, widgets in self._source_rows.items():
            for w in widgets:
                self._form.setRowVisible(w, src == cur)
        self._form.setRowVisible(self._box_row, cur != "mixture")

    def structure(self) -> Structure | None:
        return self._structure

    def error(self) -> str:
        return self._error

    def is_building(self) -> bool:
        return self._building

    def set_structure(self, s: Structure) -> None:
        self._token += 1; self._building = False
        rec = None
        if s.source == "recipe":
            try:
                rec = Recipe.from_ref(s.source_ref)
            except StructureError:
                rec = None
        self._restoring = True
        try:
            if rec is not None:
                self._restore_base(rec.base)
                self.recipe.set_steps(list(rec.steps))
            else:
                self.recipe.clear()
                if s.source != "recipe":
                    self._restore_base(Base(source=s.source, ref=s.source_ref))
            has_fix = rec is not None and any(st.op == "fix" for st in rec.steps)
            if has_fix:
                self.fixed.setText("")
            else:
                parts = [str(i + 1) for i in s.fixed_atoms]
                parts += [f"{int(k) + 1}:{''.join(c for c, m in zip('xyz', v) if not m)}" for k, v in s.fixed_axes.items()]
                self.fixed.setText(",".join(parts))
            self.box.setChecked(False)
            self.charge.setValue(s.charge)
            self.multiplicity.setValue(s.multiplicity)
        finally:
            self._restoring = False
        self._raw, self._built_ref = None, ""
        if rec is not None:
            self._raw = s
            self._auto_charge = rec.total_charge()
            try:
                self._built_ref = self.current_recipe().to_ref()
            except (StructureError, ValueError, OSError):
                self._built_ref = ""
            self.recipe.set_buildable(True)
            self.recipe.show_stale(L("spec.json の原子座標を使っています。手順から作り直すなら「作る」を押してください",
                                     "using the coordinates saved in spec.json; press Build to rebuild from the steps"))
        self._sync_fixed_field()
        self._structure, self._error = s, ""
        self._show_source_rows()
        self._show()
        self.changed.emit()

    def _restore_base(self, base: Base) -> None:
        src, ref = base.source, base.ref
        self.set_source(src)
        self._file_restored = None
        if src in NEW_BASES and isinstance(ref, dict):
            try:
                self._load_new_base(src, ref)
            except ValueError:
                pass
            return
        if not isinstance(ref, str):
            return
        if src == "mixture":
            from adit.mixture import MixtureError, MixtureSpec
            try:
                self.mixture.blockSignals(True); self.mixture.set_spec(MixtureSpec.from_ref(ref))
            except MixtureError:
                pass
            finally:
                self.mixture.blockSignals(False)
        elif src == "preset":
            self.preset.setCurrentIndex(max(0, self.preset.findData(ref)))
        elif src == "smiles":
            self.smiles.setText(ref)
        elif src == "file":
            self.file.setText(ref)
            if base.sha256:
                self._file_restored = (ref, base.sha256)
        elif src == "bulk":
            self._load_bulk(ref)
        elif src == "surface":
            self._load_surface(ref)

    def _load_bulk(self, ref: str) -> None:
        tokens = ref.split()
        if not tokens:
            return
        self.bulk_el.setCurrentText(tokens[0])
        rest = tokens[1:]
        self.bulk_cubic.setChecked("cubic" in rest)
        rest = [t for t in rest if t != "cubic"]
        if rest and rest[0] in CRYSTAL_STRUCTURES:
            self.bulk_struct.setCurrentIndex(self.bulk_struct.findText(rest[0])); rest = rest[1:]
        else:
            self.bulk_struct.setCurrentIndex(0)
        try:
            self.bulk_a.setValue(float(rest[0]) if rest else 0.0)
        except ValueError:
            pass

    def _load_surface(self, ref: str) -> None:
        tokens = ref.split()
        if len(tokens) < 3:
            return
        self.surf_facet.setCurrentText(tokens[0]); self.surf_el.setCurrentText(tokens[1])
        try:
            for w, v in zip(self.surf_n, tokens[2].lower().split("x")):
                w.setValue(int(v))
            for t in tokens[3:]:
                k, _, v = t.partition("=")
                if k == "vacuum":
                    self.surf_vac.setValue(float(v))
        except ValueError:
            pass

    def set_source(self, key: str) -> None:
        self.source.setCurrentIndex(max(0, self.source.findData(key)))

    def current_source(self) -> str:
        return self.source.currentData() or "preset"

    def recipe_mode(self) -> bool:
        return self.recipe.count() > 0 or self.current_source() in NEW_BASES

    def needs_build_button(self) -> bool:
        return self.recipe.count() > 0 or self.current_source() in SLOW_BASES

    def current_recipe(self) -> Recipe:
        return Recipe(base=self._base(), steps=self.recipe.steps())

    def _base(self) -> Base:
        src = self.current_source()
        if src == "2d":
            return Base(source="2d", ref=self._twod_ref())
        if src == "cluster":
            return Base(source="cluster", ref=self._cluster_ref())
        if src == "polymer":
            return Base(source="polymer", ref=self._polymer_ref())
        if src == "file":
            path = self.file.text().strip()
            if not path:
                raise StructureError(L("構造ファイルを指定してください", "choose a structure file"))
            if self._file_restored and self._file_restored[0] == path:
                return Base(source="file", ref=path, sha256=self._file_restored[1])
            p = Path(path).expanduser()
            if not p.is_file():
                raise StructureError(L(f"構造ファイル {p} がありません", f"structure file {p} was not found"))
            stat = p.stat(); key = (str(p.resolve()), stat.st_mtime_ns, stat.st_size)
            if key not in self._sha_cache:
                from adit.builder import file_sha256
                self._sha_cache[key] = file_sha256(p)
            return Base(source="file", ref=path, sha256=self._sha_cache[key])
        source, ref = self._source()
        if not ref:
            raise StructureError(self._empty_message(source))
        return Base(source=source, ref=ref)

    @staticmethod
    def _empty_message(source: str) -> str:
        empty = {"file": L("構造ファイルを指定してください", "choose a structure file"),
                 "smiles": L("SMILES を入力するか、「Draw」で描いてください", "enter a SMILES string, or draw the molecule with \"Draw\""),
                 "mixture": L("成分を 1 つ以上追加してください", "add at least one component")}
        return empty.get(source, L("構造が指定されていません", "no structure specified"))

    def _source(self) -> tuple[str, str]:
        kind = self.current_source()
        if kind == "preset":
            return "preset", self.preset.currentData() or ""
        if kind == "smiles":
            return "smiles", self.smiles.text().strip()
        if kind == "bulk":
            el = self.bulk_el.currentText()
            struct = self.bulk_struct.currentText()
            parts = [el]
            if self.bulk_struct.currentIndex() != 0:
                parts.append(struct)
            if self.bulk_a.value() > 0:
                parts.append(f"{self.bulk_a.value():g}")
            if self.bulk_cubic.isChecked():
                parts.append("cubic")
            return "bulk", " ".join(parts)
        if kind == "mixture":
            return "mixture", self.mixture.spec().to_ref()
        if kind == "surface":
            n = "x".join(str(w.value()) for w in self.surf_n)
            return "surface", f"{self.surf_facet.currentText()} {self.surf_el.currentText()} {n} vacuum={self.surf_vac.value():g}"
        return "file", self.file.text().strip()

    def _rebuild(self, *_) -> None:
        if self._restoring:
            return
        self._sync_fixed_field()
        if self.recipe_mode():
            self._rebuild_recipe(); return
        self.recipe.set_buildable(False); self.recipe.show_stale(""); self.recipe.log.setText("")
        self._raw, self._built_ref, self._auto_charge = None, "", None
        source, ref = self._source()
        self.charge.setEnabled(source != "mixture")
        if source == "mixture":
            total = self.mixture.spec().total_charge()
            if self.charge.value() != total:
                self.charge.blockSignals(True); self.charge.setValue(total); self.charge.blockSignals(False)
        if not ref:
            self._structure, self._error = None, self._empty_message(source)
        else:
            try:
                st = build_structure(source, ref, charge=self.charge.value(), multiplicity=self.multiplicity.value())
                self._structure, self._error = self._apply_common(st, use_field=True), ""
            except (StructureError, ValueError) as ex:
                self._structure, self._error = None, str(ex)
        if source == "mixture":
            a = self._structure.atoms.to_ase() if self._structure is not None else None
            self.mixture.show_result(float(a.cell.lengths()[0]) if a is not None else None, self._error)
        self._show()
        self.changed.emit()

    def _apply_common(self, st: Structure, *, use_field: bool) -> Structure:
        atoms = st.atoms.to_ase()
        upd: dict = {"charge": self.charge.value(), "multiplicity": self.multiplicity.value()}
        if self.box.isChecked() and not any(atoms.pbc):
            edge = self.box_size.value()
            atoms.set_cell([edge, edge, edge]); atoms.center(); atoms.pbc = True
            upd["atoms"] = AtomsData.from_ase(atoms)
        if use_field:
            fixed_atoms, fixed_axes = parse_constraints(self.fixed.text(), len(atoms))
            upd.update(fixed_atoms=fixed_atoms, fixed_axes=fixed_axes)
        return st.model_copy(update=upd)

    def _rebuild_recipe(self) -> None:
        self.recipe.set_buildable(True)
        self.charge.setEnabled(True)
        try:
            rec = self.current_recipe()
            ref = rec.to_ref()
        except (StructureError, ValueError, OSError) as ex:
            self._structure, self._error = None, str(ex)
            if not self._building:
                self.recipe.show_stale("")
            self._show(); self.changed.emit(); return
        q = rec.total_charge()
        if q != self._auto_charge:
            self._auto_charge = q
            if self.charge.value() != q:
                self.charge.blockSignals(True); self.charge.setValue(q); self.charge.blockSignals(False)
        self._update_terminations(rec)
        if self._raw is not None and ref == self._built_ref:
            self._finish()
        elif self.needs_build_button():
            msg = (L("組み立て手順を変えました。「作る」を押すと作り直します", "the steps have changed; press Build to rebuild") if self._raw is not None
                   else L("「作る」を押すと、組み立て手順から構造を作ります", "press Build to make the structure from the steps"))
            self._structure, self._error = None, msg
            if not self._building:
                self.recipe.show_stale(msg)
        else:
            from adit.builder import recipe_structure
            try:
                st, logs = recipe_structure(rec, charge=self.charge.value(), multiplicity=self.multiplicity.value())
                self._raw, self._built_ref = st, ref
                self.recipe.show_result([x.line() for x in logs])
                self._finish()
            except StructureError as ex:
                self._raw = None
                self._structure, self._error = None, str(ex)
                self.recipe.show_result([], str(ex))
        self._show()
        self.changed.emit()

    def _finish(self) -> None:
        try:
            self._structure = self._apply_common(self._raw, use_field=not self.recipe.has_op("fix"))
            self._error = ""
        except (StructureError, ValueError) as ex:
            self._structure, self._error = None, str(ex)

    def _sync_fixed_field(self) -> None:
        by_step = self.recipe.has_op("fix")
        self.fixed.setEnabled(not by_step)
        self.fixed.setPlaceholderText(L("組み立て手順の「固定」で選んだ原子を固定します (この欄は使いません)",
                                        "atoms chosen by the fix step are fixed (this field is not used)") if by_step else tr(self._fixed_placeholder))

    def build_recipe(self) -> None:
        try:
            rec = self.current_recipe()
        except (StructureError, ValueError, OSError) as ex:
            self._structure, self._error = None, str(ex)
            self.recipe.show_result([], str(ex)); self._show(); self.changed.emit(); self.built.emit(); return
        self._token += 1
        token, ref = self._token, rec.to_ref()
        self._pending_ref = ref
        self._building = True
        self.recipe.set_building(True)
        threading.Thread(target=self._build_worker, args=(token, ref, self.charge.value(), self.multiplicity.value()), daemon=True).start()

    def _build_worker(self, token: int, ref: str, charge: int, multiplicity: int) -> None:
        from adit.builder import recipe_structure
        try:
            st, logs = recipe_structure(ref, charge=charge, multiplicity=multiplicity)
            self._build_done.emit(token, st, [x.line() for x in logs], "")
        except StructureError as ex:
            self._build_done.emit(token, None, [], str(ex))
        except Exception as ex:
            self._build_done.emit(token, None, [], L(f"作れませんでした ({type(ex).__name__}: {ex})", f"could not build ({type(ex).__name__}: {ex})"))

    def cancel_build(self) -> None:
        if not self._building:
            return
        self._token += 1; self._building = False
        self.recipe.show_cancelled()

    def _on_built(self, token: int, st, lines, error: str) -> None:
        if token != self._token:
            return
        self._building = False
        if error:
            self._raw = None
            self._structure, self._error = None, error
            self.recipe.show_result(list(lines), error)
            self._show(); self.changed.emit()
        else:
            self._raw, self._built_ref = st, self._pending_ref
            self.recipe.show_result(list(lines))
            self._rebuild()
        self.built.emit()

    def _atoms_before(self, index: int, rec: Recipe | None = None):
        try:
            rec = rec or self.current_recipe()
        except (StructureError, ValueError, OSError):
            return None
        steps = rec.steps[:index]
        if rec.base.source in ("mixture", "polymer") or any(s.op in ("solvent_layer", "solvate") for s in steps):
            return None
        part = Recipe(base=rec.base, steps=steps)
        key = part.to_ref()
        if key not in self._before_cache:
            from adit.builder import build_recipe
            try:
                atoms = build_recipe(part)[0]
            except Exception:
                atoms = None
            if len(self._before_cache) > 16:
                self._before_cache.clear()
            self._before_cache[key] = atoms
        return self._before_cache[key]

    def _update_terminations(self, rec: Recipe) -> None:
        from adit.builder import list_terminations
        from adit.gui.panels.recipe_editor import SlabEditor
        for i, ed in enumerate(self.recipe.editors):
            if not isinstance(ed, SlabEditor):
                continue
            atoms = self._atoms_before(i, rec)
            names = None
            if atoms is not None and ed.miller() != (0, 0, 0):
                key = f"terms:{Recipe(base=rec.base, steps=rec.steps[:i]).to_ref()}:{ed.miller()}"
                if key not in self._before_cache:
                    try:
                        self._before_cache[key] = list_terminations(atoms, ed.miller())
                    except Exception:
                        self._before_cache[key] = None
                names = self._before_cache[key]
            ed.set_terminations(names)

    def _area_before(self, index: int) -> float | None:
        atoms = self._atoms_before(index)
        if atoms is None and self._raw is not None:
            atoms = self._raw.atoms.to_ase()
        if atoms is None or not all(atoms.pbc):
            return None
        return float(np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1])))

    def _base_is_slab(self) -> bool:
        src = self.current_source()
        if src in ("surface", "2d"):
            return True
        if src != "file":
            return False
        from adit.builder import slab_cell_problem
        from adit.structure import from_file
        try:
            atoms = from_file(self.file.text().strip())
        except StructureError:
            return False
        if slab_cell_problem(atoms) is not None:
            return False
        f = np.sort(atoms.get_scaled_positions(wrap=True)[:, 2])
        gaps = np.diff(np.concatenate([f, [f[0] + 1.0]]))
        return float(gaps.max() * atoms.cell.lengths()[2]) > 4.0

    def _draw_smiles(self) -> None:
        from adit.gui.sketcher import SketchDialog

        dlg = SketchDialog(self.smiles.text().strip(), self)
        if dlg.exec():
            self.smiles.setText(dlg.smiles()); self.set_source("smiles"); self._rebuild()

    def _draw_into_mixture(self, row: int) -> None:
        self._draw_into(self.mixture, row)

    def _draw_into(self, target, row: int) -> None:
        from adit.gui.sketcher import SketchDialog

        dlg = SketchDialog(target.row_smiles(row), self)
        if dlg.exec():
            target.set_row_smiles(row, dlg.smiles())

    def _show(self) -> None:
        if self._structure is None:
            self.info.setText(L(f"構造がありません: {self._error}", f"no structure: {self._error}"))
            self.asegui.setEnabled(False)
        else:
            a = self._structure.atoms.to_ase()
            per = L("周期系", "periodic") if any(a.pbc) else L("分子 (非周期)", "molecule (non-periodic)")
            nfix = len(self._structure.fixed_atoms)
            fixed = L(f"、固定 {nfix} 原子", f", {nfix} fixed") if self._structure.fixed_atoms else ""
            elems = ' '.join(sorted(set(a.get_chemical_symbols())))
            self.info.setText(L(f"{pretty_formula(a.get_chemical_formula())}、{len(a)} 原子、元素: {elems}、{per}{fixed}",
                                f"{pretty_formula(a.get_chemical_formula())}, {len(a)} atoms, elements: {elems}, {per}{fixed}"))
            self.asegui.setEnabled(True)

    def _filter_presets(self, text: str) -> None:
        from PySide6.QtCore import QSignalBlocker

        needle = text.strip().casefold()
        names = [name for name in self._preset_names if preset_matches(name, needle)]
        current = self.preset.currentData()
        with QSignalBlocker(self.preset):
            self.preset.clear()
            for name in names:
                self.preset.addItem(pretty_formula(name), name)
            wanted = current if current in names else (names[0] if len(names) == 1 else None)
            if wanted is not None:
                self.preset.setCurrentIndex(self.preset.findData(wanted))
        if current != self.preset.currentData():
            self._rebuild()

    def _on_file_edited(self) -> None:
        if self._file_restored and self.file.text().strip() != self._file_restored[0]:
            self._file_restored = None
        self._rebuild()

    def _browse(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, L("構造ファイル", "Structure file"), "",
                                                   L("構造 (*.xyz *.extxyz *.cif *.pdb *.gen *.vasp POSCAR CONTCAR);;すべて (*)",
                                                     "Structures (*.xyz *.extxyz *.cif *.pdb *.gen *.vasp POSCAR CONTCAR);;All files (*)"))
        if path:
            self._file_restored = None
            self.set_source("file")
            self.file.setText(path)
            self._rebuild()

    def _open_ase_gui(self) -> None:
        if self._structure is None:
            return
        tmp = Path(tempfile.mkdtemp(prefix="adit_")) / "structure.xyz"
        write(tmp, self._structure.atoms.to_ase())
        subprocess.Popen([sys.executable, "-m", "ase", "gui", str(tmp)])
