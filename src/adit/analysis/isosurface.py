"""Isosurfaces of volumetric grids (orbitals, densities, potentials) as triangles, with numpy only.

Marching tetrahedra: every grid cube is split into six tetrahedra around its main diagonal and each
tetrahedron is cut by the level with the case logic of PolygoniseTri in
https://paulbourke.net/geometry/polygonise/ (Bourke, "Polygonising a scalar field"; source1.c).
The cube corner numbering is the one of that page: 0 (0,0,0), 1 (1,0,0), 2 (1,1,0), 3 (0,1,0),
4 (0,0,1), 5 (1,0,1), 6 (1,1,1), 7 (0,1,1). The six tetrahedra all contain the diagonal 0-6, so the
cut of every cube face is the same from both neighbouring cubes and the surface has no cracks.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from adit.errors import AditValueError
from adit.lang import L

MAX_TRIANGLES = 20_000            # QPainter redraw measured at 0.3 s for 22k triangles; the browser canvas allows more
TRIANGLES_PER_ACTIVE_CUBE = 6.0   # measured on Gaussian spheres (40^3 and 96^3 grids); used for the estimate only

CORNERS = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
TETRAHEDRA = ((0, 1, 2, 6), (0, 2, 3, 6), (0, 3, 7, 6), (0, 7, 4, 6), (0, 4, 5, 6), (0, 5, 1, 6))


class IsosurfaceError(AditValueError):
    pass


class IsosurfaceTooLarge(IsosurfaceError):

    def __init__(self, message: str, suggested_stride: int):
        super().__init__(message)
        self.suggested_stride = suggested_stride


@dataclass
class Mesh:
    vertices: np.ndarray          # (3 * n_triangles, 3) in Å; three consecutive rows form one triangle
    normals: np.ndarray           # (3 * n_triangles, 3) unit normals pointing away from the enclosed values
    faces: np.ndarray             # (n_triangles, 3) indices into vertices
    level: float
    stride: int = 1
    notes: list[str] = field(default_factory=list)

    @property
    def n_triangles(self) -> int:
        return int(len(self.faces))

    def triangles(self) -> np.ndarray:
        return self.vertices[self.faces]

    def face_normals(self) -> np.ndarray:
        return self.normals[self.faces].mean(axis=1)


def _tetra_triangles(case: int) -> list[tuple[tuple[int, int], ...]]:
    """Edges (as local corner pairs) of the triangles cut out of one tetrahedron for a 4-bit below/above case."""
    below = [t for t in range(4) if case & (1 << t)]
    above = [t for t in range(4) if t not in below]
    if len(below) in (0, 4):
        return []
    if len(below) in (1, 3):
        lone = below[0] if len(below) == 1 else above[0]
        others = [t for t in range(4) if t != lone]
        return [tuple((lone, o) for o in others)]
    a, b = below
    c, d = above
    return [((a, c), (a, d), (b, d)), ((a, c), (b, d), (b, c))]


TRI_TABLE = {case: _tetra_triangles(case) for case in range(16)}


def _prepare(values: np.ndarray, stride: int, periodic: bool) -> np.ndarray:
    v = np.asarray(values, dtype=float)
    if v.ndim != 3:
        raise IsosurfaceError(L(f"3 次元の格子が要ります (次元 {v.ndim})", f"a 3D grid is required (got {v.ndim} dimensions)"))
    stride = max(1, int(stride))
    if stride > 1:
        v = v[::stride, ::stride, ::stride]
    if periodic:
        v = np.pad(v, ((0, 1), (0, 1), (0, 1)), mode="wrap")
    if min(v.shape) < 2:
        raise IsosurfaceError(L(f"間引き {stride} では各方向 2 点未満になります: {v.shape}", f"stride {stride} leaves fewer than 2 points per axis: {v.shape}"))
    return v


def _active_cubes(v: np.ndarray, level: float) -> np.ndarray:
    """Boolean (nx-1, ny-1, nz-1): cubes whose eight corners are not all on one side of the level."""
    inside = v > level
    any_in = np.zeros(tuple(n - 1 for n in v.shape), dtype=bool)
    all_in = np.ones_like(any_in)
    for dx, dy, dz in CORNERS:
        block = inside[dx:dx + any_in.shape[0], dy:dy + any_in.shape[1], dz:dz + any_in.shape[2]]
        any_in |= block
        all_in &= block
    return any_in & ~all_in


def check_budget(values: np.ndarray, levels, stride: int = 1, periodic: bool = False,
                 max_triangles: int | None = MAX_TRIANGLES) -> dict:
    """Estimate the triangle count before extracting; raise IsosurfaceTooLarge with a stride that should fit."""
    v = _prepare(values, stride, periodic)
    levels = [float(x) for x in (levels if isinstance(levels, (list, tuple, np.ndarray)) else [levels])]
    n_active = sum(int(_active_cubes(v, lv).sum()) for lv in levels)
    est = n_active * TRIANGLES_PER_ACTIVE_CUBE
    info = {"active_cubes": n_active, "estimated_triangles": int(est), "stride": int(max(1, stride)),
            "grid_shape": tuple(int(n) for n in v.shape), "max_triangles": max_triangles}
    if max_triangles is not None and est > max_triangles:
        # triangles scale with the surface area in grid units, i.e. with 1 / stride^2
        factor = math.sqrt(est / max_triangles)
        suggested = max(int(max(1, stride)) + 1, math.ceil(max(1, stride) * factor))
        raise IsosurfaceTooLarge(L(
            f"等値面の三角形が多すぎるので描かずに止めました。見積もり: 約 {int(est):,} 枚 (上限 {max_triangles:,} 枚、"
            f"格子 {'×'.join(str(n) for n in values.shape)})。格子を間引いてください: 例 {suggested} 点に 1 点。"
            "上限そのものは max_triangles で変えられます",
            f"the isosurface has too many triangles, so it was not drawn. Estimate: about {int(est):,} (limit {max_triangles:,}; "
            f"grid {'×'.join(str(n) for n in values.shape)}). Thin the grid, e.g. every {suggested}th point. "
            "The limit itself can be changed with max_triangles"), suggested)
    return info


def _gradient_at(v: np.ndarray, idx: np.ndarray) -> np.ndarray:
    """Central differences of v (in index units) at integer grid points idx (..., 3), one-sided at the borders."""
    shape = np.array(v.shape)
    out = np.empty(idx.shape, dtype=float)
    for ax in range(3):
        lo = idx.copy(); hi = idx.copy()
        lo[..., ax] = np.maximum(idx[..., ax] - 1, 0)
        hi[..., ax] = np.minimum(idx[..., ax] + 1, shape[ax] - 1)
        span = (hi[..., ax] - lo[..., ax]).astype(float)
        span[span == 0] = 1.0
        out[..., ax] = (v[hi[..., 0], hi[..., 1], hi[..., 2]] - v[lo[..., 0], lo[..., 1], lo[..., 2]]) / span
    return out


def isosurface(grid, level: float, stride: int = 1, max_triangles: int | None = MAX_TRIANGLES,
               periodic: bool | None = None) -> Mesh:
    """Triangles of the surface values == level. grid: adit.analysis.volumetric.Grid (values, cell, origin, periodic)."""
    level = float(level)
    if not math.isfinite(level):
        raise IsosurfaceError(L("等値は有限の数で指定してください", "the level must be a finite number"))
    values = np.asarray(grid.values, dtype=float)
    per = bool(grid.periodic if periodic is None else periodic)
    stride = max(1, int(stride))
    check_budget(values, level, stride, per, max_triangles)
    v = _prepare(values, stride, per)
    step = np.asarray(grid.cell, dtype=float) / np.asarray(values.shape, dtype=float)[:, None] * stride   # rows: one index step
    origin = np.asarray(getattr(grid, "origin", np.zeros(3)), dtype=float)
    inv_step = np.linalg.inv(step)
    active = np.argwhere(_active_cubes(v, level))
    tri_pos: list[np.ndarray] = []
    tri_grad: list[np.ndarray] = []
    if len(active):
        corner_idx = active[:, None, :] + CORNERS[None, :, :]                     # (n, 8, 3)
        corner_val = v[corner_idx[..., 0], corner_idx[..., 1], corner_idx[..., 2]]  # (n, 8)
        corner_grad = _gradient_at(v, corner_idx)                                   # (n, 8, 3)
        for tet in TETRAHEDRA:
            t_idx = corner_idx[:, tet, :].astype(float)                             # (n, 4, 3)
            t_val = corner_val[:, tet]
            t_grad = corner_grad[:, tet, :]
            case = np.zeros(len(active), dtype=int)
            for t in range(4):
                case |= (t_val[:, t] < level).astype(int) << t
            for c in range(1, 15):
                rows = np.nonzero(case == c)[0]
                if len(rows) == 0:
                    continue
                for tri in TRI_TABLE[c]:
                    pts = np.empty((len(rows), 3, 3)); grads = np.empty((len(rows), 3, 3))
                    for k, (a, b) in enumerate(tri):
                        va, vb = t_val[rows, a], t_val[rows, b]
                        t = (level - va) / (vb - va)
                        pts[:, k, :] = t_idx[rows, a, :] + t[:, None] * (t_idx[rows, b, :] - t_idx[rows, a, :])
                        grads[:, k, :] = t_grad[rows, a, :] + t[:, None] * (t_grad[rows, b, :] - t_grad[rows, a, :])
                    tri_pos.append(pts); tri_grad.append(grads)
    if tri_pos:
        pos = np.concatenate(tri_pos).reshape(-1, 3)
        grad = np.concatenate(tri_grad).reshape(-1, 3)
    else:
        pos = np.zeros((0, 3)); grad = np.zeros((0, 3))
    verts = origin + pos @ step
    normals = grad @ inv_step.T                       # index-space gradient -> Cartesian gradient
    if level >= 0:
        normals = -normals                            # the enclosed region has the larger values
    lengths = np.linalg.norm(normals, axis=1)
    lengths[lengths == 0] = 1.0
    normals = normals / lengths[:, None]
    faces = np.arange(len(verts)).reshape(-1, 3)
    if len(faces):
        p = verts[faces]
        fn = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
        flip = np.einsum("ij,ij->i", fn, normals[faces].mean(axis=1)) < 0
        faces[flip, 1], faces[flip, 2] = faces[flip, 2], faces[flip, 1]
    mesh = Mesh(vertices=verts, normals=normals, faces=faces, level=level, stride=stride)
    if stride > 1:
        mesh.notes.append(L(f"格子を {stride} 点に 1 点に間引いて描いています", f"drawn from every {stride}th grid point"))
    if per:
        mesh.notes.append(L("周期的な格子として、セルの端で値をつないでいます", "treated as a periodic grid: the values are wrapped at the cell faces"))
    return mesh


def isosurfaces(grid, level: float, stride: int = 1, max_triangles: int | None = MAX_TRIANGLES,
                periodic: bool | None = None) -> list[Mesh]:
    """Surfaces at +level and, when the grid has values below -level, also at -level (orbital lobes). level > 0."""
    level = float(level)
    if not level > 0:
        raise IsosurfaceError(L("等値は正の数で指定してください (負の側は自動で -等値 に描きます)",
                                "give a positive level (the negative side is drawn at -level automatically)"))
    values = np.asarray(grid.values, dtype=float)
    levels = [level] + ([-level] if float(values.min()) < -level else [])
    per = bool(grid.periodic if periodic is None else periodic)
    check_budget(values, levels, stride, per, max_triangles)
    return [isosurface(grid, lv, stride, None, per) for lv in levels]


def grid_stats(grid) -> dict:
    """Min / max / |max| of the values and the size of the grid, for choosing a level (no default is suggested)."""
    values = np.asarray(grid.values, dtype=float)
    vmin, vmax = float(values.min()), float(values.max())
    amax = max(abs(vmin), abs(vmax))
    step = np.asarray(grid.cell, dtype=float) / np.asarray(values.shape, dtype=float)[:, None]
    return {"min": vmin, "max": vmax, "abs_max": amax, "tenth_of_abs_max": 0.1 * amax, "shape": tuple(int(n) for n in values.shape),
            "n_points": int(values.size), "step_ang": [float(np.linalg.norm(s)) for s in step],
            "has_negative": bool(vmin < 0), "unit": getattr(grid, "unit", ""), "kind": getattr(grid, "kind", ""),
            "periodic": bool(getattr(grid, "periodic", False))}


def payload(meshes: list[Mesh], digits: int = 3) -> list[dict]:
    """Flat triangle lists (9 floats per triangle) with one normal per triangle, for the browser and JSON."""
    out = []
    for m in meshes:
        tri = np.round(m.triangles().reshape(-1), digits)
        nrm = np.round(m.face_normals().reshape(-1), 3)
        out.append({"level": m.level, "n_triangles": m.n_triangles, "stride": m.stride, "notes": list(m.notes),
                    "tri": tri.tolist(), "normal": nrm.tolist()})
    return out


__all__ = ["MAX_TRIANGLES", "Mesh", "IsosurfaceError", "IsosurfaceTooLarge", "check_budget", "isosurface", "isosurfaces",
           "grid_stats", "payload"]
