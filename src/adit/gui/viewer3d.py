
from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from ase.data.colors import jmol_colors
from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QMouseEvent, QPainter, QPen, QRadialGradient, QWheelEvent
from PySide6.QtWidgets import QWidget

from adit.measure import measure, measure_text

_BOND_FACTOR = 1.2
_UNKNOWN_RGB = (0.6, 0.6, 0.6)


def _radius(z: int) -> float:
    return float(covalent_radii[z]) if 0 <= z < len(covalent_radii) else 1.5


def _rgb(z: int) -> tuple[float, float, float]:
    return tuple(jmol_colors[z]) if 0 <= z < len(jmol_colors) else _UNKNOWN_RGB
_MAX_BOND_ATOMS = 2000
_REBOND_ATOMS = 400          # bonds are recomputed per frame only up to this size
_CLICK_PX = 4.0
MAX_SELECTED = 4


def _rot_x(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def _rot_y(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


class Viewer3D(QWidget):
    changed = Signal()
    selectionChanged = Signal(list)   # 0-based indices of the displayed atoms, in click order

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        self.setMinimumSize(200, 200)
        self.setMouseTracking(False)
        self._atoms: Atoms | None = None
        self._pos = np.zeros((0, 3)); self._sym: list[str] = []; self._num: list[int] = []
        self._bonds: list[tuple[int, int]] = []
        self._cell_lines: list[tuple[np.ndarray, np.ndarray]] = []
        self._center = np.zeros(3)
        self._rot = _rot_x(np.radians(-60)) @ _rot_y(np.radians(30))
        self._zoom = 1.0
        self._pan = np.zeros(2)
        self._last: QPointF | None = None
        self._press: QPointF | None = None
        self._dragged = False
        self._button = Qt.MouseButton.NoButton
        self._xy = np.zeros((0, 2)); self._r = np.zeros(0); self._depth = np.zeros(0)
        self.selected: list[int] = []
        self.dark = False
        self._surfaces: list[tuple[np.ndarray, np.ndarray, QColor]] = []   # (triangles (m,3,3) Å absolute, normals (m,3), colour)
        self._light = np.array([0.3, 0.5, 0.81]); self._light /= np.linalg.norm(self._light)

    def set_atoms(self, atoms: Atoms | None) -> None:
        self._atoms = atoms
        had = bool(self.selected)
        self.selected = []
        if had:
            self.selectionChanged.emit([])
        if atoms is None or len(atoms) == 0:
            self._pos = np.zeros((0, 3)); self._sym = []; self._num = []; self._bonds = []; self._cell_lines = []
            self._xy = np.zeros((0, 2)); self._r = np.zeros(0); self._depth = np.zeros(0)
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
        self._center = center
        self._pos = pos - center
        self._cell_lines = [(a - center, b - center) for a, b in self._cell_lines]
        self._bonds = self._find_bonds(atoms)
        self.update()

    def set_positions(self, pos) -> None:
        """Move the atoms (same count) without resetting the view, the center or the selection."""
        if self._atoms is None or len(self._atoms) != len(pos):
            return
        self._atoms.set_positions(np.asarray(pos, dtype=float))
        self._pos = np.asarray(pos, dtype=float) - self._center
        if len(self._atoms) <= _REBOND_ATOMS:
            self._bonds = self._find_bonds(self._atoms)
        self.update()

    def atoms(self) -> Atoms | None:
        return self._atoms

    # ---- isosurfaces ----
    def set_surfaces(self, surfaces) -> None:
        """surfaces: iterable of (triangles (m,3,3) in Å, absolute coordinates; face normals (m,3); QColor with alpha)."""
        out = []
        for tri, nrm, color in surfaces:
            tri = np.asarray(tri, dtype=float).reshape(-1, 3, 3)
            nrm = np.asarray(nrm, dtype=float).reshape(-1, 3)
            if len(tri) and len(nrm) == len(tri):
                out.append((tri, nrm, QColor(color)))
        self._surfaces = out
        self.update()

    def clear_surfaces(self) -> None:
        if self._surfaces:
            self._surfaces = []
            self.update()

    def n_triangles(self) -> int:
        return sum(len(t) for t, _, _ in self._surfaces)

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

    # ---- selection and measurement ----
    def set_selection(self, indices) -> None:
        n = len(self._pos)
        sel: list[int] = []
        for i in indices:
            i = int(i)
            if 0 <= i < n and i not in sel:
                sel.append(i)
        self.selected = sel[-MAX_SELECTED:]
        self.update(); self.selectionChanged.emit(list(self.selected))

    def clear_selection(self) -> None:
        if self.selected:
            self.set_selection([])

    def measurement(self) -> dict | None:
        if self._atoms is None or len(self.selected) < 2:
            return None
        cell = np.asarray(self._atoms.cell, dtype=float) if any(self._atoms.pbc) else None
        return measure(self._atoms.get_positions(), self.selected, cell)

    def measurement_text(self) -> str:
        return measure_text(self.measurement(), self._sym)

    def pick(self, x: float, y: float) -> int | None:
        """Index of the front-most atom drawn under the point, or None."""
        if len(self._xy) == 0:
            return None
        d = np.hypot(self._xy[:, 0] - x, self._xy[:, 1] - y)
        hit = np.nonzero(d <= self._r + 1.5)[0]
        if len(hit) == 0:
            return None
        return int(hit[np.argmax(self._depth[hit])])

    def _click(self, x: float, y: float, additive: bool) -> None:
        i = self.pick(x, y)
        if i is None:
            if not additive:
                self.clear_selection()
            return
        if additive:
            sel = [k for k in self.selected if k != i] if i in self.selected else self.selected + [i]
        else:
            sel = [i]
        self.set_selection(sel)

    def _scale(self) -> float:
        if len(self._pos) == 0:
            return 20.0
        pts = np.vstack([self._pos] + [np.array([a, b]) for a, b in self._cell_lines]) if self._cell_lines else self._pos
        extent = max(1.0, float(np.abs(pts).max()) + 1.5)
        return 0.45 * min(self.width(), self.height()) / extent * self._zoom

    def _accent(self) -> QColor:
        return QColor("#3A95FF" if self.dark else "#0A7AFF")

    def paintEvent(self, _ev) -> None:  # noqa: N802
        p = QPainter(self); p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(QColor(255, 255, 255, 18 if self.dark else 150))
        p.drawRoundedRect(self.rect().adjusted(0, 0, -1, -1), 10, 10)
        if len(self._pos) == 0:
            self._xy = np.zeros((0, 2)); self._r = np.zeros(0); self._depth = np.zeros(0)
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
        zmin, zmax = float(depth.min()), float(depth.max()); zr = max(zmax - zmin, 1e-6)
        radii = np.zeros(len(self._pos))
        selected = set(self.selected)
        tri_xy, tri_depth, tri_col = self._project_surfaces(s, cx, cy)
        n_atoms = len(depth)
        order = np.argsort(np.concatenate([depth, tri_depth])) if len(tri_depth) else np.argsort(depth)
        for k in order:
            if k >= n_atoms:
                self._draw_triangle(p, tri_xy[k - n_atoms], tri_col[k - n_atoms])
                continue
            i = int(k)
            z = self._num[i]
            r = max(2.0, _radius(z) * 0.55 * s); radii[i] = r
            col = QColor.fromRgbF(*_rgb(z))
            shade = 0.75 + 0.25 * (depth[i] - zmin) / zr
            col = QColor(int(col.red() * shade), int(col.green() * shade), int(col.blue() * shade))
            g = QRadialGradient(QPointF(xy[i][0] - r * 0.35, xy[i][1] - r * 0.35), r * 1.4)
            g.setColorAt(0.0, col.lighter(140)); g.setColorAt(1.0, col.darker(120))
            p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            p.setBrush(QBrush(g)); p.setPen(QPen(QColor(0, 0, 0, 90), 0.8)); p.drawEllipse(QPointF(*xy[i]), r, r)
            if i in selected:
                ring = QPen(self._accent()); ring.setWidthF(2.5); p.setPen(ring); p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawEllipse(QPointF(*xy[i]), r + 3.0, r + 3.0)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        self._xy, self._r, self._depth = xy, radii, depth
        self._draw_measurement(p, xy)
        self._draw_axes(p)
        p.end()

    def _project_surfaces(self, s: float, cx: float, cy: float):
        """Screen triangles of all surfaces: (m,3,2) points, (m,) depth, list of shaded colours."""
        if not self._surfaces:
            return np.zeros((0, 3, 2)), np.zeros(0), []
        xy_all, depth_all, cols = [], [], []
        for tri, nrm, color in self._surfaces:
            rp = (tri - self._center) @ self._rot.T                       # (m,3,3) rotated
            xy = rp[..., :2] * s; xy[..., 1] *= -1; xy[..., 0] += cx; xy[..., 1] += cy
            depth = rp[..., 2].mean(axis=1)
            shade = 0.45 + 0.55 * np.abs((nrm @ self._rot.T) @ self._light)   # two-sided flat shading
            r, g, b, a = color.red(), color.green(), color.blue(), color.alpha()
            cols += [QColor(int(r * f), int(g * f), int(b * f), a) for f in shade]
            xy_all.append(xy); depth_all.append(depth)
        return np.concatenate(xy_all), np.concatenate(depth_all), cols

    @staticmethod
    def _draw_triangle(p: QPainter, xy: np.ndarray, color: QColor) -> None:
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)         # no seams between neighbouring triangles
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(color)
        p.drawConvexPolygon([QPointF(xy[0][0], xy[0][1]), QPointF(xy[1][0], xy[1][1]), QPointF(xy[2][0], xy[2][1])])

    def _draw_measurement(self, p: QPainter, xy: np.ndarray) -> None:
        if not self.selected:
            return
        accent = self._accent()
        pen = QPen(accent); pen.setWidthF(2.0); pen.setStyle(Qt.PenStyle.DashLine); pen.setCapStyle(Qt.PenCapStyle.RoundCap); p.setPen(pen)
        sel = self.selected
        for a, b in zip(sel, sel[1:]):
            p.drawLine(QPointF(*xy[a]), QPointF(*xy[b]))
        font = p.font(); font.setBold(True); p.setFont(font)
        for k, i in enumerate(sel):
            self._label(p, xy[i][0] + self._r[i] + 4, xy[i][1] - self._r[i] - 4, str(k + 1), accent)
        text = self.measurement_text()
        if not text:
            return
        if len(sel) == 2:
            at = (xy[sel[0]] + xy[sel[1]]) / 2
        elif len(sel) == 3:
            at = xy[sel[1]] + (0, 18)
        else:
            at = (xy[sel[1]] + xy[sel[2]]) / 2
        self._label(p, float(at[0]) + 6, float(at[1]) + 6, text, accent)

    def _label(self, p: QPainter, x: float, y: float, text: str, color: QColor) -> None:
        fm = p.fontMetrics()
        w, h = fm.horizontalAdvance(text) + 8, fm.height() + 4
        x = min(max(2.0, x), max(2.0, self.width() - w - 2)); y = min(max(2.0, y), max(2.0, self.height() - h - 2))
        box = QRectF(x, y, w, h)
        p.setPen(Qt.PenStyle.NoPen); p.setBrush(QColor(30, 30, 30, 200) if self.dark else QColor(255, 255, 255, 220))
        p.drawRoundedRect(box, 4, 4)
        p.setPen(QPen(color)); p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawText(box, Qt.AlignmentFlag.AlignCenter, text)

    AXIS_COLORS = ("#D9534F", "#5CB85C", "#4A90D9")

    def _draw_axes(self, p: QPainter) -> None:
        arm = 22.0
        ox, oy = 14 + arm, self.height() - 14 - arm
        font = p.font(); font.setBold(False); font.setPointSizeF(max(7.0, font.pointSizeF() - 1)); p.setFont(font)
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
        self._last = ev.position(); self._press = ev.position(); self._dragged = False; self._button = ev.button()

    def mouseMoveEvent(self, ev: QMouseEvent) -> None:  # noqa: N802
        if self._last is None:
            return
        d = ev.position() - self._last; self._last = ev.position()
        if self._press is not None and (ev.position() - self._press).manhattanLength() > _CLICK_PX:
            self._dragged = True
        if self._button == Qt.MouseButton.LeftButton:
            self._rot = _rot_x(np.radians(d.y() * 0.5)) @ _rot_y(np.radians(d.x() * 0.5)) @ self._rot
        elif self._button == Qt.MouseButton.RightButton:
            self._pan += (d.x(), d.y())
        self.update(); self.changed.emit()

    def mouseReleaseEvent(self, ev: QMouseEvent) -> None:  # noqa: N802
        if self._button == Qt.MouseButton.LeftButton and not self._dragged and self._press is not None:
            self._click(self._press.x(), self._press.y(), bool(ev.modifiers() & Qt.KeyboardModifier.ShiftModifier))
        self._last = None; self._press = None; self._button = Qt.MouseButton.NoButton

    def mouseDoubleClickEvent(self, _ev: QMouseEvent) -> None:  # noqa: N802
        self.reset()

    def wheelEvent(self, ev: QWheelEvent) -> None:  # noqa: N802
        self._zoom *= 1.1 ** (ev.angleDelta().y() / 120); self._zoom = min(20.0, max(0.1, self._zoom)); self.update()
