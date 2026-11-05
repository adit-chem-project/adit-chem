"""Voronoi volumes and face counts. Periodic cells are computed with surrounding image atoms."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from adit.lang import L


@dataclass
class VoronoiResult:
    volume: np.ndarray
    faces: np.ndarray
    max_face_area: np.ndarray
    threshold_ang2: float


def _images(positions: np.ndarray, cell: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    shifts = np.array([[i, j, k] for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
    out = [positions]
    owner = [np.arange(len(positions))]
    for s in shifts:
        if not np.any(s):
            continue
        out.append(positions + s @ cell)
        owner.append(np.full(len(positions), -1))
    return np.vstack(out), np.concatenate(owner)


def _polyhedron_volume(points: np.ndarray) -> float:
    from scipy.spatial import ConvexHull

    return float(ConvexHull(points).volume)


def voronoi(positions, cell=None, face_area_threshold_ang2: float = 0.0) -> VoronoiResult:
    """Per-atom Voronoi volume and face count. Pass cell for periodic systems."""
    from scipy.spatial import Voronoi

    pos = np.asarray(positions, dtype=float)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(L("positions は (原子数, 3) にしてください", "positions must have shape (atoms, 3)"))
    if len(pos) < 4:
        raise ValueError(L("Voronoi には 4 原子以上が要ります", "Voronoi analysis needs 4 or more atoms"))
    if cell is not None:
        all_pos, owner = _images(pos, np.asarray(cell, dtype=float))
    else:
        all_pos, owner = pos, np.arange(len(pos))
    vor = Voronoi(all_pos)
    n = len(pos)
    volume = np.full(n, np.nan)
    faces = np.zeros(n, dtype=int)
    max_area = np.zeros(n)
    for i in range(n):
        region = vor.regions[vor.point_region[i]]
        if not region or -1 in region:
            continue
        verts = vor.vertices[region]
        volume[i] = _polyhedron_volume(verts)
        for (p, q), ridge in zip(vor.ridge_points, vor.ridge_vertices):
            if i not in (p, q) or -1 in ridge or len(ridge) < 3:
                continue
            area = _face_area(vor.vertices[ridge])
            if area > face_area_threshold_ang2:
                faces[i] += 1
                max_area[i] = max(max_area[i], area)
    _ = owner
    return VoronoiResult(volume, faces, max_area, face_area_threshold_ang2)


def _face_area(verts: np.ndarray) -> float:
    center = verts.mean(axis=0)
    v = verts - center
    total = 0.0
    for a, b in zip(v, np.roll(v, -1, axis=0)):
        total += np.linalg.norm(np.cross(a, b)) / 2.0
    return float(total)


def summary_lines(vr: VoronoiResult, cell_volume_ang3: float | None = None) -> list[str]:
    ok = np.isfinite(vr.volume)
    lines = [L(f"Voronoi (閉じた多面体 {int(ok.sum())} / {vr.volume.size} 原子)",
               f"Voronoi (closed polyhedra for {int(ok.sum())} / {vr.volume.size} atoms)"),
             L(f"  体積: 平均 {vr.volume[ok].mean():.4f} Å³、最小 {vr.volume[ok].min():.4f}、最大 {vr.volume[ok].max():.4f}",
               f"  volume: mean {vr.volume[ok].mean():.4f} Å³, min {vr.volume[ok].min():.4f}, max {vr.volume[ok].max():.4f}"),
             L(f"  面の数: 平均 {vr.faces[ok].mean():.2f} (面積 {vr.threshold_ang2:g} Å² より大きい面だけ)",
               f"  faces: mean {vr.faces[ok].mean():.2f} (only faces larger than {vr.threshold_ang2:g} Å²)")]
    if cell_volume_ang3 is not None and ok.all():
        lines.append(L(f"  体積の合計 {vr.volume.sum():.4f} Å³ / セル {cell_volume_ang3:.4f} Å³ "
                       f"(比 {vr.volume.sum() / cell_volume_ang3:.6f})",
                       f"  total volume {vr.volume.sum():.4f} Å³ / cell {cell_volume_ang3:.4f} Å³ "
                       f"(ratio {vr.volume.sum() / cell_volume_ang3:.6f})"))
    return lines
