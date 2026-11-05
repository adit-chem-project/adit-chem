"""Displacement vectors and atomic strain relative to a reference structure."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from adit.analysis.msd_worker import mic
from adit.lang import L


@dataclass
class LocalStrain:
    deformation: np.ndarray
    strain: np.ndarray
    volumetric: np.ndarray
    shear: np.ndarray
    residual: np.ndarray
    neighbors: np.ndarray


def displacements(reference, current, cell=None) -> np.ndarray:
    """Displacement from the reference, shape (n, 3) in Angstrom. Pass cell for minimum-image distances."""
    r0 = np.asarray(reference, dtype=float)
    r1 = np.asarray(current, dtype=float)
    if r0.shape != r1.shape:
        raise ValueError(L(f"原子数が違います: 基準 {r0.shape[0]} 個、いま {r1.shape[0]} 個",
                           f"the atom counts differ: {r0.shape[0]} in the reference, {r1.shape[0]} now"))
    d = r1 - r0
    return mic(d, np.asarray(cell, dtype=float)) if cell is not None else d


def local_strain(reference, current, cutoff_ang: float, cell=None, cell_current=None) -> LocalStrain:
    r0 = np.asarray(reference, dtype=float)
    r1 = np.asarray(current, dtype=float)
    if r0.shape != r1.shape:
        raise ValueError(L(f"原子数が違います: 基準 {r0.shape[0]} 個、いま {r1.shape[0]} 個",
                           f"the atom counts differ: {r0.shape[0]} in the reference, {r1.shape[0]} now"))
    if cutoff_ang <= 0:
        raise ValueError(L("カットオフは正の値にしてください", "the cutoff must be positive"))
    cell0 = np.asarray(cell, dtype=float) if cell is not None else None
    cell1 = np.asarray(cell_current, dtype=float) if cell_current is not None else cell0
    n = r0.shape[0]
    F = np.zeros((n, 3, 3))
    strain = np.zeros((n, 3, 3))
    vol = np.zeros(n)
    shear = np.zeros(n)
    resid = np.zeros(n)
    counts = np.zeros(n, dtype=int)
    identity = np.eye(3)
    for i in range(n):
        d0 = r0 - r0[i]
        if cell0 is not None:
            d0 = mic(d0, cell0)
        dist = np.linalg.norm(d0, axis=1)
        idx = np.nonzero((dist > 0) & (dist <= cutoff_ang))[0]
        vec0 = d0[idx]
        counts[i] = len(idx)
        if len(idx) < 3:
            F[i] = identity
            continue
        d1 = r1[idx] - r1[i]
        if cell1 is not None:
            d1 = mic(d1, cell1)
        sol, *_ = np.linalg.lstsq(vec0, d1, rcond=None)
        F[i] = sol.T
        resid[i] = float(np.sqrt(((vec0 @ sol - d1) ** 2).sum(axis=1).mean()))
        e = 0.5 * (F[i].T @ F[i] - identity)
        strain[i] = e
        vol[i] = np.trace(e) / 3.0
        dev = e - np.trace(e) / 3.0 * identity
        shear[i] = float(np.sqrt(2.0 / 3.0 * (dev * dev).sum()))
    return LocalStrain(F, strain, vol, shear, resid, counts)


def summary_lines(u: np.ndarray, ls: LocalStrain | None = None) -> list[str]:
    mag = np.linalg.norm(u, axis=1)
    lines = [L(f"変位 (原子 {mag.size} 個): 平均 {mag.mean():.4f} Å、最大 {mag.max():.4f} Å",
               f"displacement ({mag.size} atoms): mean {mag.mean():.4f} Å, max {mag.max():.4f} Å")]
    if ls is not None:
        ok = ls.neighbors >= 3
        lines.append(L(f"  局所ひずみ (近傍 3 個以上の {int(ok.sum())} 原子): "
                       f"体積ひずみ 平均 {ls.volumetric[ok].mean():+.5f}、"
                       f"せん断 平均 {ls.shear[ok].mean():.5f}、最大 {ls.shear[ok].max():.5f}",
                       f"  local strain ({int(ok.sum())} atoms with 3 or more neighbors): "
                       f"volumetric strain mean {ls.volumetric[ok].mean():+.5f}, "
                       f"shear mean {ls.shear[ok].mean():.5f}, max {ls.shear[ok].max():.5f}"))
        lines.append(L(f"  当てはめの残差 平均 {ls.residual[ok].mean():.4f} Å "
                       "(大きいときは、その原子のまわりを一様な変形で表せていません)",
                       f"  mean fit residual {ls.residual[ok].mean():.4f} Å "
                       "(a large value means a uniform deformation does not describe the surroundings of that atom)"))
    lines.append(L("  原子の対応が基準と同じ順であることが前提です (入れ替わると意味がありません)。",
                   "  Atoms must be in the same order as in the reference (the result means nothing if they swap)."))
    return lines
