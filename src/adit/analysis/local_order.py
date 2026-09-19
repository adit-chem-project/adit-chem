
from __future__ import annotations

import numpy as np
from ase import Atoms

from adit.errors import AditValueError
from adit.lang import L


class LocalOrderError(AditValueError):
    pass


def _neighbors(atoms: Atoms, cutoff: float, indices=None):
    from adit.analysis.compute import pairs_within

    if cutoff <= 0:
        raise LocalOrderError(L("カットオフ [Å] は正の値です", "the cutoff (Å) must be positive"))
    i, j, _ = pairs_within(atoms, cutoff)
    pos = atoms.get_positions()
    cell = np.asarray(atoms.cell, dtype=float)
    periodic = bool(np.any(atoms.pbc)) and abs(np.linalg.det(cell)) > 0
    from adit.analysis.compute import _mic_step

    out: dict[int, list] = {k: [] for k in range(len(atoms))}
    if len(i):
        diff = pos[j] - pos[i]
        if periodic:
            diff = _mic_step(diff, cell)
        for a, b, d in zip(i.tolist(), j.tolist(), diff):
            out[a].append((b, d))
            out[b].append((a, -d))
    if indices is not None:
        keep = set(int(x) for x in indices)
        out = {k: v for k, v in out.items() if k in keep}
    return out


def coordination(atoms: Atoms, cutoff: float, indices=None) -> dict:
    nb = _neighbors(atoms, cutoff, indices)
    keys = sorted(nb)
    values = [len(nb[k]) for k in keys]
    syms = atoms.get_chemical_symbols()
    return {"index": [k + 1 for k in keys], "symbol": [syms[k] for k in keys], "coordination": values,
            "cutoff_A": float(cutoff), "mean": float(np.mean(values)) if values else None,
            "note": L("カットオフ以内にある原子の数です。カットオフは利用者が指定した値で、"
                      "ADIT は既定値を持ちません (第一配位圏の切り方は化学的な判断です)。",
                      "the number of atoms within the cutoff; the cutoff is the value you gave. "
                      "ADIT has no default (where to cut the first coordination shell is a chemical judgement).")}


def centrosymmetry(atoms: Atoms, n_neighbors: int = 12, indices=None) -> dict:
    if n_neighbors < 2 or n_neighbors % 2:
        raise LocalOrderError(L("相手の数は 2 以上の偶数です (FCC なら 12、BCC なら 8)",
                                "the number of neighbors must be an even number of at least 2 (12 for FCC, 8 for BCC)"))
    pos = atoms.get_positions()
    cell = np.asarray(atoms.cell, dtype=float)
    periodic = bool(np.any(atoms.pbc)) and abs(np.linalg.det(cell)) > 0
    from adit.analysis.compute import _mic_step

    targets = range(len(atoms)) if indices is None else [int(x) for x in indices]
    values = []
    for a in targets:
        diff = pos - pos[a]
        if periodic:
            diff = _mic_step(diff, cell)
        d = np.linalg.norm(diff, axis=1)
        order = np.argsort(d)
        near = [k for k in order if k != a][:n_neighbors]
        if len(near) < n_neighbors:
            values.append(float("nan"))
            continue
        vectors = diff[near]
        used = set()
        total = 0.0
        for p in range(len(vectors)):
            if p in used:
                continue
            best, best_val = None, None
            for q in range(len(vectors)):
                if q == p or q in used:
                    continue
                val = float(np.sum((vectors[p] + vectors[q]) ** 2))
                if best_val is None or val < best_val:
                    best, best_val = q, val
            if best is None:
                break
            used.update({p, best})
            total += best_val
        values.append(total)
    syms = atoms.get_chemical_symbols()
    return {"index": [a + 1 for a in targets], "symbol": [syms[a] for a in targets],
            "centrosymmetry_A2": values, "n_neighbors": int(n_neighbors),
            "note": L("中心対称性パラメータ [Å²] (Kelchner ら 1998)。完全な中心対称の格子で 0 になります。"
                      "どこからを欠陥・表面と呼ぶかは判定していません (利用者が決めます)。",
                      "centrosymmetry parameter in Å² (Kelchner et al. 1998); zero for a perfectly centrosymmetric "
                      "lattice. No threshold for calling an atom a defect or a surface is applied.")}


def steinhardt(atoms: Atoms, cutoff: float, ls=(4, 6), indices=None) -> dict:
    from scipy.special import sph_harm_y

    nb = _neighbors(atoms, cutoff, indices)
    keys = sorted(nb)
    syms = atoms.get_chemical_symbols()
    out: dict[str, list] = {f"q{l}": [] for l in ls}
    for k in keys:
        vectors = np.array([d for _, d in nb[k]], dtype=float)
        if len(vectors) == 0:
            for l in ls:
                out[f"q{l}"].append(float("nan"))
            continue
        r = np.linalg.norm(vectors, axis=1)
        theta = np.arccos(np.clip(vectors[:, 2] / r, -1.0, 1.0))
        phi = np.arctan2(vectors[:, 1], vectors[:, 0])
        for l in ls:
            total = 0.0
            for m in range(-l, l + 1):
                qlm = np.mean(sph_harm_y(l, m, theta, phi))
                total += abs(qlm) ** 2
            out[f"q{l}"].append(float(np.sqrt(4 * np.pi / (2 * l + 1) * total)))
    return {"index": [k + 1 for k in keys], "symbol": [syms[k] for k in keys], "cutoff_A": float(cutoff),
            **out,
            "note": L("Steinhardt の q_l (1983)。カットオフ以内の相手の向きを球面調和関数で平均した量です。"
                      "値がどの構造に当たるかは判定していません (文献の値と比べるのは利用者です)。",
                      "Steinhardt's q_l (1983), from the directions of the neighbors within the cutoff. "
                      "ADIT does not say which structure a value corresponds to.")}


def clusters(atoms: Atoms, cutoff: float, indices=None) -> dict:
    nb = _neighbors(atoms, cutoff, indices)
    seen: dict[int, int] = {}
    sizes: list[int] = []
    for start in sorted(nb):
        if start in seen:
            continue
        label = len(sizes)
        stack, size = [start], 0
        seen[start] = label
        while stack:
            node = stack.pop()
            size += 1
            for other, _ in nb.get(node, []):
                if other in nb and other not in seen:
                    seen[other] = label
                    stack.append(other)
        sizes.append(size)
    order = np.argsort(sizes)[::-1]
    return {"cutoff_A": float(cutoff), "n_clusters": len(sizes), "sizes": [int(sizes[i]) for i in order],
            "largest": int(max(sizes)) if sizes else 0,
            "cluster_of_atom": {str(k + 1): int(v) for k, v in sorted(seen.items())},
            "note": L("カットオフ以内でつながった原子のかたまりです。カットオフは利用者が指定します。"
                      "かたまりが「分子」「凝集体」かどうかは判定していません。",
                      "groups of atoms connected within the cutoff; the cutoff is yours. "
                      "ADIT does not decide whether a group is a molecule or an aggregate.")}


def angle_distribution(frames, center: str | None = None, cutoff: float = 0.0, nbins: int = 90,
                       outer: str | None = None) -> dict:
    if cutoff <= 0:
        raise LocalOrderError(L("角度の分布にはカットオフ [Å] が要ります (第一配位圏の切り方は利用者が決めます)",
                                "an angle distribution needs a cutoff in Å (where to cut the first shell is your choice)"))
    edges = np.linspace(0.0, 180.0, nbins + 1)
    hist = np.zeros(nbins)
    pairs = 0
    for fr in frames:
        syms = fr.get_chemical_symbols()
        nb = _neighbors(fr, cutoff)
        for k, items in nb.items():
            if center and syms[k] != center:
                continue
            vectors = [(other, d) for other, d in items if not outer or syms[other] == outer]
            for p in range(len(vectors)):
                for q in range(p + 1, len(vectors)):
                    u, v = vectors[p][1], vectors[q][1]
                    cos = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
                    hist[min(int(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))) / 180.0 * nbins), nbins - 1)] += 1
                    pairs += 1
    centers = 0.5 * (edges[1:] + edges[:-1])
    return {"angle_deg": centers.tolist(), "counts": hist.tolist(),
            "normalized": (hist / pairs).tolist() if pairs else hist.tolist(),
            "center_element": center, "outer_element": outer, "cutoff_A": float(cutoff), "n_angles": int(pairs),
            "note": L("カットオフ以内の相手 2 つが中心の原子に対して作る角度の分布です。"
                      "カットオフは利用者の指定で、山の帰属 (四面体・八面体など) は判定していません。",
                      "distribution of angles formed at each central atom by two neighbors within the cutoff; "
                      "the cutoff is yours and no assignment of the peaks is made.")}


def structure_factor(r: np.ndarray, g: np.ndarray, number_density: float, q: np.ndarray | None = None) -> dict:
    r = np.asarray(r, dtype=float)
    g = np.asarray(g, dtype=float)
    if r.size < 4 or r.size != g.size:
        raise LocalOrderError(L("g(r) の点が足りません (r と g は同じ長さで 4 点以上)",
                                "not enough points in g(r) (r and g must have the same length, at least 4)"))
    if number_density <= 0:
        raise LocalOrderError(L("数密度 [Å⁻³] が要ります (正の値)", "a positive number density in Å⁻³ is required"))
    if q is None:
        q = np.linspace(0.3, 12.0, 240)
    q = np.asarray(q, dtype=float)
    dr = float(r[1] - r[0])
    integrand = r ** 2 * (g - 1.0)
    window = np.sinc(r / r[-1])
    s = np.empty_like(q)
    for k, value in enumerate(q):
        s[k] = 1.0 + 4.0 * np.pi * number_density * np.sum(integrand * np.sinc(value * r / np.pi) * window) * dr
    return {"q_1_A": q.tolist(), "s_q": s.tolist(), "number_density_A3": float(number_density),
            "rmax_A": float(r[-1]), "reliable_above_q": float(2 * np.pi / r[-1]),
            "note": L(f"g(r) を {r[-1]:.1f} Å まで積んだ S(q) です (Lorch 窓を掛けて打ち切りのリップルを抑えています)。"
                      f"q が {2 * np.pi / r[-1]:.2f} Å⁻¹ より小さいところは、積む範囲が有限なことの影響を受けます。"
                      "測定との一致は判定していません。",
                      f"S(q) from g(r) integrated to {r[-1]:.1f} Å with a Lorch window to damp truncation ripples; "
                      f"below q = {2 * np.pi / r[-1]:.2f} A^-1 the finite integration range matters. "
                      "No comparison with experiment is judged.")}
