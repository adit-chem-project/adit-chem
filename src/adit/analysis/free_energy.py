"""Free-energy surface A(s) = -kT ln P(s) from a distribution."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from adit.lang import L

KB_EV_K = 8.617333262e-5
KB_KJ_MOL_K = 0.008314462618
KB_KCAL_MOL_K = 0.001987204259

UNITS = {"eV": KB_EV_K, "kJ/mol": KB_KJ_MOL_K, "kcal/mol": KB_KCAL_MOL_K}


@dataclass
class FreeEnergySurface:
    edges: list
    centers: list
    counts: np.ndarray
    free_energy: np.ndarray
    unit: str
    temperature_k: float


def read_colvar(path) -> tuple[list[str], np.ndarray]:
    from pathlib import Path

    names: list[str] = []
    rows = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        s = line.strip()
        if s.startswith("#!"):
            parts = s.split()
            if len(parts) > 2 and parts[1] == "FIELDS":
                names = parts[2:]
            continue
        if not s or s.startswith("#"):
            continue
        rows.append([float(x) for x in s.split()])
    if not rows:
        raise ValueError(L(f"{Path(path).name} に数値の行がありません", f"{Path(path).name} has no numeric rows"))
    data = np.array(rows)
    if names and len(names) != data.shape[1]:
        raise ValueError(L(f"{Path(path).name}: 列の名前 {len(names)} 個と数値の列 {data.shape[1]} 個が合いません",
                           f"{Path(path).name}: {len(names)} column names do not match {data.shape[1]} numeric columns"))
    return names, data


def free_energy_surface(values, temperature_k: float, bins=50, weights=None,
                        unit: str = "kJ/mol", ranges=None) -> FreeEnergySurface:
    """A(s) = -kT ln P(s) for 1D or 2D data. Optional weights are used as counting weights."""
    if unit not in UNITS:
        raise ValueError(L(f"単位は {', '.join(UNITS)} のどれかにしてください", f"the unit must be one of {', '.join(UNITS)}"))
    if temperature_k <= 0:
        raise ValueError(L("温度は正の値にしてください", "the temperature must be positive"))
    v = np.asarray(values, dtype=float)
    if v.ndim == 1:
        v = v[:, None]
    if v.ndim != 2 or v.shape[1] not in (1, 2):
        raise ValueError(L("values は 1 次元か、(点の数, 2) の 2 次元にしてください",
                           "values must be 1-dimensional or 2-dimensional with shape (points, 2)"))
    w = None if weights is None else np.asarray(weights, dtype=float)
    if w is not None and w.shape[0] != v.shape[0]:
        raise ValueError(L(f"重み {w.shape[0]} 個が点の数 {v.shape[0]} 個と合いません",
                           f"{w.shape[0]} weights do not match {v.shape[0]} points"))
    kt = UNITS[unit] * temperature_k
    if v.shape[1] == 1:
        counts, edges = np.histogram(v[:, 0], bins=bins, weights=w, range=ranges)
        edge_list, center_list = [edges], [0.5 * (edges[:-1] + edges[1:])]
    else:
        counts, ex, ey = np.histogram2d(v[:, 0], v[:, 1], bins=bins, weights=w, range=ranges)
        edge_list = [ex, ey]
        center_list = [0.5 * (ex[:-1] + ex[1:]), 0.5 * (ey[:-1] + ey[1:])]
    total = counts.sum()
    if total <= 0:
        raise ValueError(L("その範囲に点がありません", "no points in that range"))
    p = counts / total
    with np.errstate(divide="ignore"):
        a = -kt * np.log(p)
    a -= a[np.isfinite(a)].min()
    return FreeEnergySurface(edge_list, center_list, counts, a, unit, temperature_k)


def summary_lines(fes: FreeEnergySurface) -> list[str]:
    finite = fes.free_energy[np.isfinite(fes.free_energy)]
    empty = int((fes.counts == 0).sum())
    dim = len(fes.centers)
    lines = [L(f"自由エネルギー面 ({dim} 次元、{fes.temperature_k:g} K、単位 {fes.unit})",
               f"free-energy surface ({dim}D, {fes.temperature_k:g} K, unit {fes.unit})"),
             L(f"  最小を 0 にしたときの最大 {finite.max():.3f} {fes.unit}"
               f" (点の入っている区間 {finite.size} 個、空の区間 {empty} 個)",
               f"  maximum {finite.max():.3f} {fes.unit} with the minimum set to 0"
               f" ({finite.size} bins with points, {empty} empty bins)")]
    if dim == 1:
        lo = int(np.argmin(fes.free_energy))
        lines.append(L(f"  最小の位置 {fes.centers[0][lo]:.4f}", f"  position of the minimum {fes.centers[0][lo]:.4f}"))
    lines.append(L("  空の区間は値が出ません。収束したかどうかは判定しません (同じ計算を分けて比べてください)。",
                   "  Empty bins have no value. Convergence is not judged (split the same run into parts and compare them)."))
    return lines


@dataclass
class Distribution2D:
    x_centers: np.ndarray
    y_centers: np.ndarray
    counts: np.ndarray
    density: np.ndarray


def distribution_2d(x, y, bins=50, ranges=None, weights=None) -> Distribution2D:
    xs = np.asarray(x, dtype=float).ravel()
    ys = np.asarray(y, dtype=float).ravel()
    if xs.size != ys.size:
        raise ValueError(L(f"2 つの列の長さが違います: {xs.size} と {ys.size}", f"the two columns differ in length: {xs.size} and {ys.size}"))
    if xs.size == 0:
        raise ValueError(L("点がありません", "no points"))
    w = None if weights is None else np.asarray(weights, dtype=float).ravel()
    if w is not None and w.size != xs.size:
        raise ValueError(L(f"重み {w.size} 個が点の数 {xs.size} 個と合いません", f"{w.size} weights do not match {xs.size} points"))
    counts, ex, ey = np.histogram2d(xs, ys, bins=bins, range=ranges, weights=w)
    area = np.outer(np.diff(ex), np.diff(ey))
    total = counts.sum()
    density = counts / (total * area) if total > 0 else np.zeros_like(counts)
    return Distribution2D(0.5 * (ex[:-1] + ex[1:]), 0.5 * (ey[:-1] + ey[1:]), counts, density)
