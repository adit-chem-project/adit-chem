
from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from PySide6.QtCore import QPointF, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QMouseEvent, QPainter, QPen, QRadialGradient, QWheelEvent
from PySide6.QtWidgets import QWidget

_BOND_FACTOR = 1.2
_UNKNOWN_RGB = (0.6, 0.6, 0.6)


def _radius(z: int) -> float:
    return float(covalent_radii[z]) if 0 <= z < len(covalent_radii) else 1.5


def _rgb(z: int) -> tuple[float, float, float]:
    return tuple(jmol_colors[z]) if 0 <= z < len(jmol_colors) else _UNKNOWN_RGB
_MAX_BOND_ATOMS = 2000


def _rot_x(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def _rot_y(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


class Viewer3D(QWidget):
    changed = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setMinimumSize(200, 200)
        self.setMouseTracking(False)
        self._atoms: Atoms | None = None
        self._pos = np.zeros((0, 3)); self._sym: list[str] = []; self._num: list[int] = []
        self._bonds: list[tuple[int, int]] = []
        self._cell_lines: list[tuple[np.ndarray, np.ndarray]] = []
        self._rot = _rot_x(np.radians(-60)) @ _rot_y(np.radians(30))
        self._zoom = 1.0
        self._pan = np.zeros(2)
        self._last: QPointF | None = None
        self._button = Qt.MouseButton.NoButton
        self.dark = False

    def set_atoms(self, atoms: Atoms | None) -> None:
        self._atoms = atoms
        if atoms is None or len(atoms) == 0:
            self._pos = np.zeros((0, 3)); self._sym = []; self._num = []; self._bonds = []; self._cell_lines = []
            self.update(); return
        pos = atoms.get_positions()
        self._sym = atoms.get_chemical_symbols(); self._num = list(atoms.get_atomic_numbers())
        self._cell_lines = []
        if any(atoms.pbc):
            cell = np.asarray(atoms.cell)
            corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)]) @ cell
            for a, b in ((0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7)):
                self._cell_lines.append((corners[a], corners[b]))
            center = cell.sum(axis=0) / 2
        else:
            center = pos.mean(axis=0)
        self._pos = pos - center
        self._cell_lines = [(a - center, b - center) for a, b in self._cell_lines]
        self._bonds = self._find_bonds(atoms)
        self.update()

    def _find_bonds(self, atoms: Atoms) -> list[tuple[int, int]]:
        n = len(atoms)
        if n < 2 or n > _MAX_BOND_ATOMS:
            return []
        r = np.array([_radius(z) for z in self._num])
        d = atoms.get_all_distances(mic=False)
        cut = _BOND_FACTOR * (r[:, None] + r[None, :])
        i, j = np.where((d < cut) & (d > 0.1))
        return [(int(a), int(b)) for a, b in zip(i, j) if a < b]

    def set_view(self, rotation: str) -> None:
        m = np.eye(3)
        for part in rotation.split(","):
            part = part.strip()
            if not part:
                continue
            ang, axis = float(part[:-1]), part[-1]
            if axis == "x":
                m = _rot_x(np.radians(ang)) @ m
            elif axis == "y":
                m = _rot_y(np.radians(ang)) @ m
            elif axis == "z":
                c, s = np.cos(np.radians(ang)), np.sin(np.radians(ang))
                m = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]) @ m
        self._rot = m; self.update()

    def reset(self) -> None:
        self._zoom = 1.0; self._pan = np.zeros(2); self.update()

    def _scale(self) -> float:
        if len(self._pos) == 0:
            return 20.0
        pts = np.vstack([self._pos] + [np.array([a, b]) for a, b in self._cell_lines]) if self._cell_lines else self._pos
        extent = max(1.0, float(np.abs(pts).max()) + 1.5)
        return 0.45 * min(self.width(), self.height()) / extent * self._zoom

    def paintEvent(self, _ev) -> None:  # noqa: N802 
        p = QPainter(self); p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(QColor(255, 255, 255, 18 if self.dark else 150))
        p.drawRoundedRect(self.rect().adjusted(0, 0, -1, -1), 10, 10)
        if len(self._pos) == 0:
            p.end(); return
        s = self._scale(); cx, cy = self.width() / 2 + self._pan[0], self.height() / 2 + self._pan[1]
        rp = self._pos @ self._rot.T
        xy = rp[:, :2] * s; xy[:, 1] *= -1; xy += (cx, cy)
        depth = rp[:, 2]
        pen = QPen(QColor("#9A9AA0" if self.dark else "#808088")); pen.setStyle(Qt.PenStyle.DashLine); pen.setWidthF(1.2); p.setPen(pen)
        for a, b in self._cell_lines:
            ra, rb = a @ self._rot.T, b @ self._rot.T
            p.drawLine(QPointF(cx + ra[0] * s, cy - ra[1] * s), QPointF(cx + rb[0] * s, cy - rb[1] * s))
        bonds = sorted(self._bonds, key=lambda ij: depth[ij[0]] + depth[ij[1]])
        pen = QPen(QColor("#B0B0B6" if self.dark else "#606068")); pen.setWidthF(max(1.5, 0.12 * s)); pen.setCapStyle(Qt.PenCapStyle.RoundCap); p.setPen(pen)
        for i, j in bonds:
            p.drawLine(QPointF(*xy[i]), QPointF(*xy[j]))
        order = np.argsort(depth)
        zmin, zmax = float(depth.min()), float(depth.max()); zr = max(zmax - zmin, 1e-6)
        p.setPen(QPen(QColor(0, 0, 0, 90), 0.8))
        for i in order:
            z = self._num[i]
            r = max(2.0, _radius(z) * 0.55 * s)
            col = QColor.fromRgbF(*_rgb(z))
            shade = 0.75 + 0.25 * (depth[i] - zmin) / zr
            col = QColor(int(col.red() * shade), int(col.green() * shade), int(col.blue() * shade))
            g = QRadialGradient(QPointF(xy[i][0] - r * 0.35, xy[i][1] - r * 0.35), r * 1.4)
            g.setColorAt(0.0, col.lighter(140)); g.setColorAt(1.0, col.darker(120))
            p.setBrush(QBrush(g)); p.drawEllipse(QPointF(*xy[i]), r, r)
        self._draw_axes(p)
        p.end()

    AXIS_COLORS = ("#D9534F", "#5CB85C", "#4A90D9")

    def _draw_axes(self, p: QPainter) -> None:
        arm = 22.0
        ox, oy = 14 + arm, self.height() - 14 - arm
        font = p.font(); font.setPointSizeF(max(7.0, font.pointSizeF() - 1)); p.setFont(font)
        for k, name in enumerate("xyz"):
            v = self._rot[:, k]
            ex, ey = ox + v[0] * arm, oy - v[1] * arm
            col = QColor(self.AXIS_COLORS[k])
            if self.dark:
                col = col.lighter(125)
            pen = QPen(col); pen.setWidthF(1.6); pen.setCapStyle(Qt.PenCapStyle.RoundCap); p.setPen(pen)
            p.drawLine(QPointF(ox, oy), QPointF(ex, ey))
            p.setBrush(QBrush(col)); p.setPen(Qt.PenStyle.NoPen)
            p.drawEllipse(QPointF(ex, ey), 2.2, 2.2)
            p.setPen(QPen(col))
            p.drawText(QPointF(ox + v[0] * (arm + 9) - 3, oy - v[1] * (arm + 9) + 4), name)

    def mousePressEvent(self, ev: QMouseEvent) -> None:  # noqa: N802
        self._last = ev.position(); self._button = ev.button()

    def mouseMoveEvent(self, ev: QMouseEvent) -> None:  # noqa: N802
        if self._last is None:
            return
        d = ev.position() - self._last; self._last = ev.position()
        if self._button == Qt.MouseButton.LeftButton:
            self._rot = _rot_x(np.radians(d.y() * 0.5)) @ _rot_y(np.radians(d.x() * 0.5)) @ self._rot
        elif self._button == Qt.MouseButton.RightButton:
            self._pan += (d.x(), d.y())
        self.update(); self.changed.emit()

    def mouseReleaseEvent(self, _ev: QMouseEvent) -> None:  # noqa: N802
        self._last = None; self._button = Qt.MouseButton.NoButton

    def mouseDoubleClickEvent(self, _ev: QMouseEvent) -> None:  # noqa: N802
        self.reset()

    def wheelEvent(self, ev: QWheelEvent) -> None:  # noqa: N802
        self._zoom *= 1.1 ** (ev.angleDelta().y() / 120); self._zoom = min(20.0, max(0.1, self._zoom)); self.update()
