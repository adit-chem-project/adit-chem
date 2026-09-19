
from __future__ import annotations

import tempfile
from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QComboBox, QHBoxLayout, QLabel, QMenu, QPushButton, QVBoxLayout, QWidget

from adit.gui.atom_select import SERIES_FIELD, base_indices
from adit.gui.flow_layout import FlowLayout
from adit.gui.viewer3d import Viewer3D
from adit.lang import L
from adit.spec import Structure
from adit.structure import pretty_formula

ROTATIONS = {"-60x,30y,0z": "x+y+z 方向から", "0x,0y,0z": "z 軸から", "-90x,0y,0z": "y 軸から", "0x,90y,0z": "x 軸から"}


class StructureViewPanel(QWidget):
    send_selection = Signal(str, list)   # destination field key, 1-based indices of the base cell

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self._structure: Structure | None = None
        self._tmp = Path(tempfile.mkdtemp(prefix="adit_view_"))
        self.info = QLabel(""); self.info.setObjectName("hint"); self.info.setWordWrap(True)
        self.rotation = QComboBox()
        for k, v in ROTATIONS.items():
            self.rotation.addItem(v, k)
        self.repeat = QComboBox(); self.repeat.addItems(["1×1×1", "2×2×1", "2×2×2", "3×3×1"])
        self.asegui = QPushButton("ASE GUI で開く"); self.asegui.setObjectName("link")
        from adit.gui.copy_save import copy_button, save_button

        self.btn_copy = copy_button(self)
        self.btn_save = save_button(self)
        self.viewer = Viewer3D()
        self.image = self.viewer
        hint = QLabel("左ドラッグで回転、ホイールで拡大縮小、右ドラッグで移動、ダブルクリックでリセット。原子をクリックで選択、Shift+クリックで追加 (2 個で距離、3 個で角度、4 個で二面角)")
        hint.setObjectName("hint"); hint.setWordWrap(True)
        self.measure = QLabel(""); self.measure.setObjectName("measure"); self.measure.setWordWrap(True)
        self.measure.setStyleSheet("font-weight: 600;")
        self.btn_clear_sel = QPushButton("選択を消す"); self.btn_clear_sel.setObjectName("link")
        self.btn_send = QPushButton("欄に送る"); self.btn_send.setObjectName("link")
        self.send_menu = QMenu(self.btn_send)
        self.act_send_fixed = self.send_menu.addAction("固定原子の欄へ")
        self.act_send_select = self.send_menu.addAction("解析の「原子の選び方」へ")
        self.act_send_series = self.send_menu.addAction("解析の時系列の欄へ (距離・角度・二面角)")
        self.btn_send.setMenu(self.send_menu)
        from adit.gui.help import help_for

        for b, key in ((self.btn_send, "欄に送る"), (self.btn_clear_sel, "選択を消す")):
            h = help_for(key)
            if h is not None:
                b.setToolTip(h.text())
        top = FlowLayout()
        for w in (QLabel("視点"), self.rotation, QLabel("繰り返し数"), self.repeat, self.btn_copy, self.btn_save, self.asegui,
                  self.btn_clear_sel, self.btn_send):
            top.addWidget(w)
        lay = QVBoxLayout(self); lay.setContentsMargins(12, 12, 12, 12); lay.setSpacing(8)
        lay.addLayout(top); lay.addWidget(self.viewer, 1); lay.addWidget(self.measure); lay.addWidget(hint); lay.addWidget(self.info)
        self.viewer.selectionChanged.connect(self._on_selection)
        self.btn_clear_sel.clicked.connect(self.viewer.clear_selection)
        self.act_send_fixed.triggered.connect(lambda: self._send("fixed"))
        self.act_send_select.triggered.connect(lambda: self._send("select"))
        self.act_send_series.triggered.connect(lambda: self._send("series"))
        self._on_selection([])
        self.rotation.currentIndexChanged.connect(lambda *_: self.viewer.set_view(self.rotation.currentData()))
        self.repeat.currentIndexChanged.connect(self._render)
        self.repeat.activated.connect(self._on_user_repeat)
        self._user_repeat = False
        self.asegui.clicked.connect(self._open_ase_gui)
        self.btn_copy.clicked.connect(self._copy_image)
        self.btn_save.clicked.connect(self._save_image)

    def set_dark(self, dark: bool) -> None:
        self.viewer.dark = dark; self.viewer.update()

    def n_base(self) -> int:
        return len(self._structure.atoms.symbols) if self._structure is not None else 0

    def selected_indices(self) -> list[int]:
        return base_indices(self.viewer.selected, self.n_base()) if self.n_base() else []

    def _on_selection(self, sel: list) -> None:
        n = len(sel)
        self.btn_clear_sel.setEnabled(n > 0); self.btn_send.setEnabled(n > 0)
        self.act_send_series.setEnabled(n in SERIES_FIELD)
        text = self.viewer.measurement_text()
        if n == 1:
            i = self.viewer.selected[0]
            text = L(f"選択: {self.viewer._sym[i]}{self.selected_indices()[0]} (もう 1 個を Shift+クリックすると距離)",
                     f"selected: {self.viewer._sym[i]}{self.selected_indices()[0]} (Shift+click another atom for the distance)")
        elif n >= 2 and self.n_base() and len(self.viewer.selected) != len(self.selected_indices()):
            text += L(" (繰り返した像の同じ原子です)", " (the same atom in repeated images)")
        self.measure.setText(text)

    def _send(self, dest: str) -> None:
        idx = self.selected_indices()
        if not idx:
            return
        if dest == "series":
            dest = SERIES_FIELD.get(len(idx), "")
            if not dest:
                return
        self.send_selection.emit(dest, idx)

    def set_structure(self, st: Structure | None, error: str = "") -> None:
        self._structure = st
        self.measure.setText("")
        if st is None:
            self.viewer.set_atoms(None)
            self.info.setText(L(f"構造がありません: {error}", f"no structure: {error}")); self.asegui.setEnabled(False)
            return
        self.asegui.setEnabled(True)
        a = st.atoms.to_ase()
        per = L("周期系", "periodic") if any(a.pbc) else L("分子 (非周期)", "molecule (non-periodic)")
        cell = ""
        if any(a.pbc):
            la, lb, lc = a.cell.lengths(); cell = L(f"、セル {la:.3f} × {lb:.3f} × {lc:.3f} Å", f", cell {la:.3f} × {lb:.3f} × {lc:.3f} Å")
        elems = " ".join(sorted(set(a.get_chemical_symbols())))
        if not self._user_repeat:
            want = "2×2×2" if any(a.pbc) and len(a) < 8 else "1×1×1"
            if self.repeat.currentText() != want:
                self.repeat.blockSignals(True); self.repeat.setCurrentText(want); self.repeat.blockSignals(False)
        shown = ""
        if any(a.pbc) and self.repeat.currentText() != "1×1×1":
            shown = L(f"。表示は基本セルを {self.repeat.currentText()} 並べたもの (計算するのは基本セル {len(a)} 原子)",
                      f". Shown as {self.repeat.currentText()} copies of the cell (the calculation uses the {len(a)}-atom cell)")
        self.info.setText(L(f"{pretty_formula(a.get_chemical_formula())}、{len(a)} 原子、元素: {elems}、{per}{cell}{shown}", f"{pretty_formula(a.get_chemical_formula())}, {len(a)} atoms, elements: {elems}, {per}{cell}{shown}"))
        self._render()

    def _on_user_repeat(self, *_) -> None:
        self._user_repeat = True
        if self._structure is not None:
            self.set_structure(self._structure)

    def _render(self, *_) -> None:
        if self._structure is None:
            return
        a = self._structure.atoms.to_ase()
        rep = tuple(int(x) for x in self.repeat.currentText().split("×"))
        if any(a.pbc) and rep != (1, 1, 1):
            a = a * rep
        self.viewer.set_atoms(a)

    def _open_ase_gui(self) -> None:
        if self._structure is None:
            return
        import subprocess, sys
        from ase.io import write

        path = self._tmp / "structure.xyz"
        write(str(path), self._structure.atoms.to_ase())
        subprocess.Popen([sys.executable, "-m", "ase", "gui", str(path)])

    def _scene_svg(self) -> str | None:
        from adit.export_image import scene_to_svg
        from adit.web.structure3d import scene_from_atoms

        atoms = self.viewer.atoms() if hasattr(self.viewer, "atoms") else getattr(self.viewer, "_atoms", None)
        if atoms is None or len(atoms) == 0:
            return None
        scene = scene_from_atoms(atoms)
        rot = getattr(self.viewer, "_rot", None)
        try:
            return scene_to_svg(scene)
        except ValueError:
            return None

    def _copy_image(self) -> None:
        from adit.gui.copy_save import copy_widget

        copy_widget(self.viewer)

    def _save_image(self) -> None:
        from adit.gui.copy_save import save_dialog

        save_dialog(self, "structure.svg", svg_text=self._scene_svg(), widget=self.viewer)
