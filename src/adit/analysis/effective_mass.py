"""Effective masses from the curvature of band edges."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from adit.lang import L

HBAR_SQ_OVER_ME_EV_ANG2 = 7.6199682


@dataclass
class EffectiveMass:
    m_over_me: float
    curvature_ev_ang2: float  # d^2E/dk^2 [eV·Å²]
    band_index: int
    kpoint_index: int
    n_points: int
    r_squared: float
    k_window_inv_ang: float
    kind: str


def _fit(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    coeff = np.polyfit(x, y, 2)
    resid = y - np.polyval(coeff, x)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - float((resid ** 2).sum()) / ss_tot if ss_tot > 0 else float("nan")
    return 2.0 * float(coeff[0]), r2


def effective_mass(kdist_inv_ang, energies_ev, band_index: int, kpoint_index: int,
                   n_points: int = 5, kind: str = "electron") -> EffectiveMass:
    if n_points < 3:
        raise ValueError(L("n_points は 3 以上にしてください", "n_points must be 3 or more"))
    x = np.asarray(kdist_inv_ang, dtype=float)
    e = np.asarray(energies_ev, dtype=float)
    if e.ndim != 2:
        raise ValueError(L("energies_ev は (k 点, バンド) の 2 次元にしてください", "energies_ev must be 2-dimensional (k-points, bands)"))
    if not (0 <= band_index < e.shape[1]) or not (0 <= kpoint_index < e.shape[0]):
        raise ValueError(L("band_index または kpoint_index が範囲の外です", "band_index or kpoint_index is out of range"))
    half = n_points // 2
    lo, hi = kpoint_index - half, kpoint_index + half + 1
    if lo < 0 or hi > x.size:
        raise ValueError(L(f"端の周りに {n_points} 点が取れません (k 点 {kpoint_index}、全 {x.size} 点)",
                           f"cannot take {n_points} points around the band edge (k-point {kpoint_index}, {x.size} points in all)"))
    xs, ys = x[lo:hi], e[lo:hi, band_index]
    if np.any(np.diff(xs) <= 0):
        raise ValueError(L("この範囲の k 点の距離が単調に増えていません (経路の折れ目をまたいでいます)",
                           "the k-point distance does not increase steadily in this range (it crosses a corner of the path)"))
    curvature, r2 = _fit(xs - xs[half], ys)
    if curvature == 0.0:
        raise ValueError(L("曲率が 0 です (この範囲では平らなので有効質量を出せません)",
                           "the curvature is 0 (the band is flat in this range, so there is no effective mass)"))
    return EffectiveMass(m_over_me=HBAR_SQ_OVER_ME_EV_ANG2 / curvature, curvature_ev_ang2=curvature,
                         band_index=band_index, kpoint_index=kpoint_index, n_points=n_points,
                         r_squared=r2, k_window_inv_ang=float(xs[-1] - xs[0]), kind=kind)


def at_band_edges(kdist_inv_ang, energies_ev, fermi_ev: float | None,
                  n_points: int = 5) -> list[EffectiveMass]:
    """Effective masses at the valence-band maximum and conduction-band minimum. Empty for metals."""
    if fermi_ev is None:
        return []
    e = np.asarray(energies_ev, dtype=float)
    occupied = e <= fermi_ev
    if not occupied.any() or occupied.all():
        return []
    if np.any(occupied.any(axis=0) & ~occupied.all(axis=0)):
        return []
    vb = np.where(occupied, e, -np.inf)
    cb = np.where(~occupied, e, np.inf)
    out: list[EffectiveMass] = []
    for arr, pick, kind in ((vb, np.argmax, "hole"), (cb, np.argmin, "electron")):
        per_k = arr.max(axis=1) if kind == "hole" else arr.min(axis=1)
        k = int(pick(per_k))
        b = int(pick(arr[k]))
        for width in range(n_points, 2, -2):
            try:
                out.append(effective_mass(kdist_inv_ang, e, b, k, width, kind))
                break
            except ValueError:
                continue
    return out


def summary_lines(masses: list[EffectiveMass]) -> list[str]:
    if not masses:
        return []
    lines = [L("有効質量 (経路に沿った向きの値。テンソルではありません)", "effective mass (along the path; not a tensor)")]
    for m in masses:
        who = (L("正孔 (価電子帯の頂上)", "hole (valence-band maximum)") if m.kind == "hole"
               else L("電子 (伝導帯の底)", "electron (conduction-band minimum)"))
        lines.append(L(f"  {who}: m*/m_e = {m.m_over_me:+.3f}  "
                       f"(バンド {m.band_index}、k 点 {m.kpoint_index}、{m.n_points} 点、"
                       f"幅 {m.k_window_inv_ang:.4f} 1/Å、決定係数 {m.r_squared:.4f})",
                       f"  {who}: m*/m_e = {m.m_over_me:+.3f}  "
                       f"(band {m.band_index}, k-point {m.kpoint_index}, {m.n_points} points, "
                       f"width {m.k_window_inv_ang:.4f} 1/Å, R² {m.r_squared:.4f})"))
    lines.append(L("  当てはめの幅と決定係数を見てください。幅が広すぎると 2 次式から外れます。",
                   "  Check the fit width and R². If the width is too large, the band departs from a quadratic."))
    return lines
