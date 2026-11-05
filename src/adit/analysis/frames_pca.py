"""Principal component analysis and k-means clustering of trajectory frames."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from adit.analysis.geometry_series import kabsch_rmsd
from adit.lang import L


@dataclass
class PcaResult:
    components: np.ndarray       # (n_components, 3 * n_atoms)
    projections: np.ndarray      # (n_frames, n_components)
    variance: np.ndarray
    explained: np.ndarray
    mean: np.ndarray


@dataclass
class ClusterResult:
    labels: np.ndarray
    sizes: np.ndarray
    centers: np.ndarray
    representative: np.ndarray
    inertia: float


def _aligned(coords: np.ndarray) -> np.ndarray:
    ref = coords[0]
    center = ref.mean(axis=0)
    out = np.empty_like(coords)
    out[0] = ref - center
    for i in range(1, coords.shape[0]):
        _, rotated = kabsch_rmsd(ref, coords[i])
        out[i] = rotated - center
    return out


def pca(coords, n_components: int = 3, align: bool = True) -> PcaResult:
    """Principal components of coordinates with shape (frames, atoms, 3)."""
    x = np.asarray(coords, dtype=float)
    if x.ndim != 3 or x.shape[2] != 3:
        raise ValueError(L("coords は (フレーム, 原子, 3) にしてください", "coords must have shape (frames, atoms, 3)"))
    if x.shape[0] < 3:
        raise ValueError(L("主成分分析には 3 フレーム以上が要ります", "principal component analysis needs 3 or more frames"))
    if align:
        x = _aligned(x)
    flat = x.reshape(x.shape[0], -1)
    mean = flat.mean(axis=0)
    centered = flat - mean
    n_components = min(n_components, min(centered.shape) - 1)
    if n_components < 1:
        raise ValueError(L("成分の数が 1 未満になります (フレームか原子が少なすぎます)",
                           "the number of components falls below 1 (too few frames or atoms)"))
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    var = s ** 2 / (centered.shape[0] - 1)
    return PcaResult(vt[:n_components], (u[:, :n_components] * s[:n_components]),
                     var[:n_components], var[:n_components] / var.sum(), mean)


def kmeans(points, k: int, iterations: int = 100, seed: int = 0) -> ClusterResult:
    x = np.asarray(points, dtype=float)
    if x.ndim != 2:
        raise ValueError(L("points は (点, 次元) にしてください", "points must have shape (points, dimensions)"))
    if not 1 <= k <= x.shape[0]:
        raise ValueError(L(f"群の数 k は 1 以上 {x.shape[0]} 以下にしてください", f"the number of clusters k must be from 1 to {x.shape[0]}"))
    rng = np.random.default_rng(seed)
    centers = x[rng.choice(x.shape[0], size=k, replace=False)].copy()
    labels = np.zeros(x.shape[0], dtype=int)
    for _ in range(iterations):
        d = np.linalg.norm(x[:, None, :] - centers[None, :, :], axis=2)
        new_labels = np.argmin(d, axis=1)
        if np.array_equal(new_labels, labels) and _ > 0:
            break
        labels = new_labels
        for c in range(k):
            if np.any(labels == c):
                centers[c] = x[labels == c].mean(axis=0)
    d = np.linalg.norm(x - centers[labels], axis=1)
    rep = np.array([int(np.nonzero(labels == c)[0][np.argmin(d[labels == c])]) if np.any(labels == c) else -1
                    for c in range(k)])
    sizes = np.array([int((labels == c).sum()) for c in range(k)])
    return ClusterResult(labels, sizes, centers, rep, float((d ** 2).sum()))


def summary_lines(res: PcaResult, cl: ClusterResult | None = None) -> list[str]:
    ratios = ", ".join(f"{v * 100:.1f} %" for v in res.explained)
    lines = [L(f"主成分分析 (成分 {res.explained.size} 個): 寄与率 {ratios} (合わせて {res.explained.sum() * 100:.1f} %)",
               f"principal component analysis ({res.explained.size} components): explained variance {ratios} "
               f"(total {res.explained.sum() * 100:.1f} %)")]
    if cl is not None:
        sizes = ", ".join(str(int(s)) for s in cl.sizes)
        frames = ", ".join(str(int(i) + 1) for i in cl.representative)
        lines.append(L(f"  分類 (k-means、群 {cl.sizes.size} 個): 大きさ {sizes}、中心に近いフレーム {frames}",
                       f"  clustering (k-means, {cl.sizes.size} clusters): sizes {sizes}; frames closest to the centers {frames}"))
        lines.append(L("  群の数は利用者が決めた値です。いくつが正しいかは判定しません。",
                       "  The number of clusters is your choice. How many is right is not judged."))
    return lines
