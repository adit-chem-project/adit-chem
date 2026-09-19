"""Isosurfaces (orbitals, densities, potentials) of cube / CHGCAR files over the atoms in a 3D view."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (QColorDialog, QComboBox, QDoubleSpinBox, QHBoxLayout, QLabel, QLineEdit, QPushButton, QSlider,
                               QVBoxLayout, QWidget)

from adit.gui.flow_layout import FlowLayout
from adit.gui.viewer3d import Viewer3D
from adit.lang import L

DEFAULT_POSITIVE = QColor(0x2B, 0x6F, 0xE6)   # blue
DEFAULT_NEGATIVE = QColor(0xE6, 0x4B, 0x3C)   # red


def _swatch(button: QPushButton, color: QColor) -> None:
    button.setStyleSheet(f"background: {color.name()}; color: {'#fff' if color.lightness() < 128 else '#000'}; min-width: 56px;")


class IsosurfacePanel(QWidget):

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self._dir: Path | None = None
        self._grid = None
        self._grid_key: tuple | None = None
        self._stats: dict | None = None
        self._pos_color = QColor(DEFAULT_POSITIVE); self._neg_color = QColor(DEFAULT_NEGATIVE)
        self.file = QComboBox(); self.file.setMinimumWidth(220)
        self.level = QLineEdit(); self.level.setPlaceholderText(L("等値 (ファイルの単位)", "level (in the file's unit)")); self.level.setFixedWidth(140)
        self.range = QLabel(""); self.range.setObjectName("hint"); self.range.setWordWrap(True)
        self.range.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.btn_pos = QPushButton(L("正の色", "Positive colour")); self.btn_neg = QPushButton(L("負の色", "Negative colour"))
        _swatch(self.btn_pos, self._pos_color); _swatch(self.btn_neg, self._neg_color)
        self.opacity = QSlider(Qt.Orientation.Horizontal); self.opacity.setRange(10, 100); self.opacity.setValue(55); self.opacity.setFixedWidth(120)
        self.opacity_label = QLabel(L("透明度", "Opacity"))
        self.stride = QDoubleSpinBox(); self.stride.setDecimals(0); self.stride.setRange(1, 64); self.stride.setValue(1); self.stride.setFixedWidth(64)
        self.lbl_stride = QLabel(L("間引き (格子点)", "Stride (grid points)"))
        self.btn_show = QPushButton(L("表示", "Show")); self.btn_show.setObjectName("primary")
        self.btn_clear = QPushButton(L("等値面を消す", "Remove the surface")); self.btn_clear.setObjectName("link")
        self.btn_use_stride = QPushButton(); self.btn_use_stride.hide()
        self.viewer = Viewer3D(); self.viewer.setMinimumHeight(300)
        self.note = QLabel(""); self.note.setObjectName("hint"); self.note.setWordWrap(True)
        self.note.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        top = FlowLayout()
        for w in (QLabel(L("ファイル", "File")), self.file, QLabel(L("等値", "Level")), self.level, self.btn_pos, self.btn_neg,
                  self.opacity_label, self.opacity, self.lbl_stride, self.stride, self.btn_show, self.btn_clear):
            top.addWidget(w)
        row = QHBoxLayout(); row.addWidget(self.btn_use_stride); row.addStretch()
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.setSpacing(6)
        lay.addLayout(top); lay.addWidget(self.range); lay.addLayout(row); lay.addWidget(self.viewer, 1); lay.addWidget(self.note)
        self.file.currentIndexChanged.connect(lambda *_: self._on_file())
        self.btn_pos.clicked.connect(lambda: self._pick_colour(True)); self.btn_neg.clicked.connect(lambda: self._pick_colour(False))
        self.opacity.valueChanged.connect(lambda *_: self._recolour())
        self.btn_show.clicked.connect(self.show_surface)
        self.btn_clear.clicked.connect(self.clear_surface)
        self.btn_use_stride.clicked.connect(self._use_suggested_stride)
        self._suggested: int | None = None
        self._meshes: list = []
        from adit.analysis.isosurface import MAX_TRIANGLES
        from adit.gui.help import help_for

        self.max_triangles = MAX_TRIANGLES

        for w, key in ((self.level, "等値"), (self.stride, "間引き (格子点)"), (self.opacity, "透明度"), (self.btn_show, "等値面を表示")):
            h = help_for(key)
            if h is not None:
                w.setToolTip(h.text())

    def set_dark(self, dark: bool) -> None:
        self.viewer.dark = dark; self.viewer.update()

    # ---- files ----
    def set_run_dir(self, run_dir: Path | str | None) -> int:
        """List the volumetric files of the directory; returns how many there are."""
        from adit.analysis.volumetric import find_files

        self._dir = Path(run_dir) if run_dir else None
        files = find_files(self._dir) if self._dir else []
        self.file.blockSignals(True); self.file.clear()
        for f in files:
            self.file.addItem(f.name, str(f))
        self.file.blockSignals(False)
        self._grid = None; self._grid_key = None; self._stats = None; self._meshes = []
        self.viewer.set_atoms(None); self.viewer.clear_surfaces(); self.note.setText(""); self.range.setText("")
        self.btn_use_stride.hide()
        if files:
            self._on_file()
        return len(files)

    def files(self) -> list[str]:
        return [self.file.itemData(i) for i in range(self.file.count())]

    def _load(self):
        from adit.analysis.volumetric import read_grid

        path = self.file.currentData()
        if not path:
            return None
        p = Path(path)
        try:
            key = (str(p), p.stat().st_mtime_ns)
        except OSError:
            key = (str(p), 0)
        if self._grid is None or key != self._grid_key:
            self._grid = read_grid(p); self._grid_key = key; self._stats = None
        return self._grid

    def _on_file(self) -> None:
        from adit.analysis.isosurface import grid_stats

        self.viewer.clear_surfaces(); self._meshes = []; self.btn_use_stride.hide()
        try:
            grid = self._load()
        except Exception as ex:
            self.range.setText(L(f"読めません: {ex}", f"cannot read: {ex}")); self.viewer.set_atoms(None); return
        if grid is None:
            return
        self._stats = st = grid_stats(grid)
        unit = f" {st['unit']}" if st["unit"] else ""
        shape = "×".join(str(n) for n in st["shape"])
        self.range.setText(L(
            f"値の範囲: 最小 {st['min']:.4g}{unit}、最大 {st['max']:.4g}{unit} (格子 {shape}、{st['n_points']:,} 点"
            f"{'、周期的' if st['periodic'] else ''})。参考: |最大| の 10 % = {st['tenth_of_abs_max']:.4g}{unit} "
            "(この 10 % は単に |最大| の 1/10 で、化学的な根拠はありません。等値は利用者が決めます)",
            f"value range: min {st['min']:.4g}{unit}, max {st['max']:.4g}{unit} (grid {shape}, {st['n_points']:,} points"
            f"{', periodic' if st['periodic'] else ''}). For reference: 10 % of |max| = {st['tenth_of_abs_max']:.4g}{unit} "
            "(that 10 % is just one tenth of |max|, with no chemical basis; the level is your choice)"))
        atoms = grid.atoms.copy()
        if not grid.periodic:
            atoms.pbc = False      # the cube box starts at the grid origin, not at the cell origin ASE draws; show no box
        self.viewer.set_atoms(atoms)

    # ---- colours ----
    def _pick_colour(self, positive: bool) -> None:
        current = self._pos_color if positive else self._neg_color
        c = QColorDialog.getColor(current, self, L("等値面の色", "Isosurface colour"))
        if not c.isValid():
            return
        if positive:
            self._pos_color = c; _swatch(self.btn_pos, c)
        else:
            self._neg_color = c; _swatch(self.btn_neg, c)
        self._recolour()

    def colour(self, level: float) -> QColor:
        c = QColor(self._pos_color if level >= 0 else self._neg_color)
        c.setAlpha(int(round(255 * self.opacity.value() / 100)))
        return c

    def _recolour(self) -> None:
        if self._meshes:
            self.viewer.set_surfaces([(m.triangles(), m.face_normals(), self.colour(m.level)) for m in self._meshes])

    # ---- surfaces ----
    def parse_level(self) -> float:
        from adit.analysis.isosurface import IsosurfaceError

        text = self.level.text().strip().replace("，", ",").replace(",", "")
        if not text:
            raise IsosurfaceError(L("等値を入れてください (上の範囲を見て決めます。既定値はありません)",
                                    "enter a level (choose it from the range above; there is no default)"))
        try:
            return float(text)
        except ValueError as ex:
            raise IsosurfaceError(L(f"等値を数で入れてください: {text!r}", f"the level must be a number: {text!r}")) from ex

    def show_surface(self) -> bool:
        from adit.analysis.isosurface import IsosurfaceError, IsosurfaceTooLarge, isosurface, isosurfaces

        self.btn_use_stride.hide(); self._suggested = None
        try:
            grid = self._load()
            if grid is None:
                self.note.setText(L("ファイルがありません", "no file")); return False
            level = self.parse_level()
            stride = int(self.stride.value())
            limit = self.max_triangles
            meshes = isosurfaces(grid, level, stride, limit) if level > 0 else [isosurface(grid, level, stride, limit)]
        except IsosurfaceTooLarge as ex:
            self._suggested = int(ex.suggested_stride)
            self.note.setText(str(ex))
            self.btn_use_stride.setText(L(f"間引きを {self._suggested} にする", f"Set the stride to {self._suggested}")); self.btn_use_stride.show()
            return False
        except IsosurfaceError as ex:
            self.note.setText(str(ex)); return False
        except Exception as ex:
            self.note.setText(L(f"等値面を作れません: {ex}", f"cannot build the isosurface: {ex}")); return False
        self._meshes = meshes
        self._recolour()
        parts = [L(f"{m.level:+.4g}: {m.n_triangles:,} 枚", f"{m.level:+.4g}: {m.n_triangles:,} triangles") for m in meshes]
        notes = []
        for m in meshes:
            for n in m.notes:
                if n not in notes:
                    notes.append(n)
        empty = [m for m in meshes if m.n_triangles == 0]
        if empty:
            notes.append(L("三角形が 0 枚の等値があります (その値を通る場所が格子にありません)",
                           "a level produced no triangles (no grid values cross it)"))
        self.note.setText(L(f"三角形 {' / '.join(parts)}。", f"Triangles {' / '.join(parts)}.") + (" " + " ".join(notes) if notes else ""))
        return True

    def clear_surface(self) -> None:
        self._meshes = []; self.viewer.clear_surfaces(); self.note.setText("")

    def _use_suggested_stride(self) -> None:
        if self._suggested:
            self.stride.setValue(self._suggested); self.btn_use_stride.hide(); self.show_surface()

    def meshes(self) -> list:
        return list(self._meshes)


__all__ = ["IsosurfacePanel", "DEFAULT_POSITIVE", "DEFAULT_NEGATIVE"]
