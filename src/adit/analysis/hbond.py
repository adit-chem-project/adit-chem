"""Hydrogen bonds: counts, lifetimes and the distance-angle map. The geometric criteria are supplied by the caller; there are no defaults."""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from ase import Atoms

from adit.analysis.trajectory import MEMORY_BUDGET_MB, Trajectory
from adit.errors import AditValueError
from adit.lang import L

KNOWN_CRITERIA = {
    "gromacs": {"donor_acceptor_A": 3.5, "angle_deg": 150.0,
                "source": "GROMACS gmx hbond (既定: D–A 0.35 nm、H–D–A 30 度)"},
    "vmd": {"donor_acceptor_A": 3.0, "angle_deg": 160.0,
            "source": "VMD HBonds プラグイン (既定: D–A 3.0 Å、角度のずれ 20 度)"},
    "cpptraj": {"donor_acceptor_A": 3.0, "angle_deg": 135.0,
                "source": "AmberTools cpptraj hbond (既定: D–A 3.0 Å、D–H–A 135 度)"},
}
DEFAULT_DONORS = ("N", "O", "F")
CHUNK_BYTES = 32 * 2 ** 20
MAP_BINS = (60, 60)
# definitions follow MDAnalysis HydrogenBondAnalysis.lifetime / lib.correlations.autocorrelation:
# S(tau) = < N(t0, t0 + tau) / N(t0) >_t0 ; intermittency 0 = continuous, unlimited = intermittent
LIFETIME_SOURCE = "https://docs.mdanalysis.org/stable/documentation_pages/analysis/hydrogenbonds.html"


class HydrogenBondError(AditValueError):
    pass


def criteria_note() -> str:
    parts = [f"{v['source']}" for v in KNOWN_CRITERIA.values()]
    return L("よく使われる値: " + " / ".join(parts) + "。ADIT は既定値を持ちません (どれを採るかは利用者が決めます)",
             "commonly used criteria: " + " / ".join(parts) + ". ADIT has no default; the choice is yours")


def _check_criteria(donor_acceptor_A: float, angle_deg: float) -> None:
    if donor_acceptor_A <= 0 or not 0 < angle_deg <= 180:
        raise HydrogenBondError(L("距離 [Å] と角度 [度] の両方を指定してください (角度は 0 より大きく 180 以下)",
                                  "give both a distance in Å and an angle in degrees (0 < angle <= 180)"))


def _geometry(atoms: Atoms, donors, hydrogen: str, max_dh_A: float, rmax_A: float):
    # every D-H...A triplet with D-A <= rmax_A; H is assigned to its nearest heavy atom within max_dh_A
    from adit.analysis.compute import _mic_step

    syms = np.array(atoms.get_chemical_symbols())
    pos = atoms.get_positions()
    cell = np.asarray(atoms.cell, dtype=float)
    periodic = bool(np.any(atoms.pbc)) and abs(np.linalg.det(cell)) > 0
    heavy = np.flatnonzero(np.isin(syms, list(donors)))
    hyd = np.flatnonzero(syms == hydrogen)
    empty = (np.zeros(0, dtype=int),) * 3 + (np.zeros(0), np.zeros(0))
    if not len(heavy) or not len(hyd):
        return empty

    def rel(a_idx, b_idx):
        d = pos[b_idx][None, :, :] - pos[a_idx][:, None, :]
        return _mic_step(d.reshape(-1, 3), cell).reshape(d.shape) if periodic else d

    step = max(1, CHUNK_BYTES // (len(heavy) * 24))
    parts: list[tuple] = []
    for s in range(0, len(hyd), step):
        h = hyd[s:s + step]
        v = rel(h, heavy)
        dist = np.linalg.norm(v, axis=2)
        j = np.argmin(dist, axis=1)
        ok = dist[np.arange(len(h)), j] <= max_dh_A
        if not ok.any():
            continue
        h, j, v = h[ok], j[ok], v[ok]
        donor = heavy[j]
        u = v[np.arange(len(h)), j]
        da = rel(donor, heavy)
        dda = np.linalg.norm(da, axis=2)
        cand = dda <= rmax_A
        cand[np.arange(len(h)), j] = False
        ih, ia = np.nonzero(cand)
        if not len(ih):
            continue
        w, uu = v[ih, ia], u[ih]
        cos = np.einsum("ij,ij->i", uu, w) / (np.linalg.norm(uu, axis=1) * np.linalg.norm(w, axis=1))
        ang = np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))   # D-H...A
        parts.append((donor[ih], h[ih], heavy[ia], dda[ih, ia], ang))
    if not parts:
        return empty
    return tuple(np.concatenate([p[k] for p in parts]) for k in range(5))


def count_frame(atoms: Atoms, donor_acceptor_A: float, angle_deg: float, donors=DEFAULT_DONORS,
                hydrogen: str = "H", max_dh_A: float = 1.3) -> dict:
    _check_criteria(donor_acceptor_A, angle_deg)
    d, h, a, r, ang = _geometry(atoms, donors, hydrogen, max_dh_A, donor_acceptor_A)
    keep = np.flatnonzero(ang >= angle_deg)
    pairs = [{"donor": int(d[k]) + 1, "hydrogen": int(h[k]) + 1, "acceptor": int(a[k]) + 1,
              "distance_A": float(r[k]), "angle_deg": float(ang[k])} for k in keep]
    return {"count": len(pairs), "pairs": pairs, "donors": list(donors)}


def count_series(frames, donor_acceptor_A: float, angle_deg: float, donors=DEFAULT_DONORS) -> dict:
    """Hydrogen-bond count per frame, with the mean and the count per atom."""
    counts = []
    natoms = 0
    for fr in frames:
        got = count_frame(fr, donor_acceptor_A, angle_deg, donors)
        counts.append(got["count"])
        natoms = len(fr)
    values = np.array(counts, dtype=float)
    return {"counts": counts, "mean": float(values.mean()) if len(values) else None,
            "std": float(values.std(ddof=1)) if len(values) > 1 else None,
            "per_atom": float(values.mean() / natoms) if natoms and len(values) else None,
            "criteria": {"donor_acceptor_A": float(donor_acceptor_A), "angle_deg": float(angle_deg),
                         "donors": list(donors)},
            "note": L("D–H···A の距離と角度が、指定した条件を満たす組の数です。"
                      "条件は利用者が入れた値で、ADIT は既定値を持ちません。"
                      "本数が多い・少ないの判定もしていません。",
                      "the number of D-H...A pairs that satisfy the distance and angle you gave; "
                      "adit has no default criteria and makes no judgement about the count.")}


def _too_large(n_frames: int, n_pairs: int, budget_mb: float, what: str) -> HydrogenBondError:
    need = n_frames * n_pairs
    stride = math.ceil(need / (budget_mb * 2 ** 20))
    return HydrogenBondError(L(
        f"{what}: 存在の表が大きすぎるので計算せずに止めました。見積もり {n_frames:,} フレーム × {n_pairs:,} 組 = 約 {need / 2**20:,.0f} MB "
        f"(上限 {budget_mb:,.0f} MB)。間引いてください (例 --stride {stride})。上限は --memory-mb で変えられます",
        f"{what}: the presence table is too large, so it was not computed. Estimate {n_frames:,} frames x {n_pairs:,} pairs = about {need / 2**20:,.0f} MB "
        f"(limit {budget_mb:,.0f} MB). Thin the trajectory (e.g. --stride {stride}); the limit itself can be changed with --memory-mb"))


def _cross_correlation_sum(w: np.ndarray, h: np.ndarray, tau_max: int) -> np.ndarray:
    # sum_j sum_t0 w(t0, j) h(t0 + tau, j) for tau = 0..tau_max, by FFT along time in column chunks
    n = w.shape[0]
    m = 1 << (2 * n - 1).bit_length()
    acc = np.zeros(m // 2 + 1, dtype=complex)
    step = max(1, CHUNK_BYTES // (m * 16))
    for c in range(0, w.shape[1], step):
        wf = np.fft.rfft(w[:, c:c + step], n=m, axis=0)
        hf = np.fft.rfft(h[:, c:c + step].astype(float), n=m, axis=0)
        acc += (np.conj(wf) * hf).sum(axis=1)
    return np.fft.irfft(acc, n=m)[:tau_max + 1]


def _lifetime_of(c: np.ndarray, dt: float) -> dict:
    lag = np.arange(len(c)) * dt
    trapezoid = getattr(np, "trapezoid", None) or getattr(np, "trapz")
    integral = float(trapezoid(c, lag))
    below = np.flatnonzero(c <= math.exp(-1))
    t_e = None
    if len(below):
        k = int(below[0])
        if k == 0:
            t_e = 0.0
        else:
            c0, c1 = float(c[k - 1]), float(c[k])
            t_e = float(lag[k - 1] + (c0 - math.exp(-1)) / (c0 - c1) * dt) if c0 != c1 else float(lag[k])
    return {"integral": integral, "one_over_e": t_e, "reached_one_over_e": t_e is not None, "last_value": float(c[-1])}


def lifetime(frames, donor_acceptor_A: float, angle_deg: float, donors=DEFAULT_DONORS, dt_fs: float | None = None,
             tau_max: int | None = None, budget_mb: float = MEMORY_BUDGET_MB) -> dict:
    """Intermittent and continuous autocorrelation of the hydrogen-bond presence h(t), and the lifetimes from them."""
    _check_criteria(donor_acceptor_A, angle_deg)
    n_est = frames.estimate_len() if isinstance(frames, Trajectory) else len(frames)
    limit = budget_mb * 2 ** 20
    cols: dict[tuple[int, int], int] = {}
    per_frame: list[np.ndarray] = []
    estimate = None
    for fr in frames:
        got = count_frame(fr, donor_acceptor_A, angle_deg, donors)
        # unique bonds are hydrogen-acceptor pairs, as in MDAnalysis
        idx = np.array([cols.setdefault((p["hydrogen"], p["acceptor"]), len(cols)) for p in got["pairs"]], dtype=np.int64)
        per_frame.append(idx)
        if estimate is None:
            # bonds exchange partners over time: allow four times the pairs of the first frame
            estimate = n_est * max(4 * len(idx), 1)
            if estimate > limit:
                raise _too_large(n_est, max(4 * len(idx), 1), budget_mb, L("水素結合の寿命", "hydrogen-bond lifetime"))
        if len(per_frame) * len(cols) > limit:
            raise _too_large(n_est, len(cols), budget_mb, L("水素結合の寿命", "hydrogen-bond lifetime"))
    n, p = len(per_frame), len(cols)
    if n < 2 or p == 0:
        raise HydrogenBondError(L(f"寿命を出せません (フレーム {n} 個、水素結合の組 {p} 個。2 フレーム以上と 1 組以上が要ります)",
                                  f"cannot compute lifetimes ({n} frames, {p} bonded pairs; at least 2 frames and 1 pair are needed)"))
    h = np.zeros((n, p), dtype=bool)
    for t, idx in enumerate(per_frame):
        h[t, idx] = True
    n0 = h.sum(axis=1)
    tau_max = n // 2 if tau_max is None else int(tau_max)
    tau_max = int(min(max(tau_max, 1), n - 1))
    w = h / np.maximum(n0, 1)[:, None]
    valid = np.array([int(np.count_nonzero(n0[:n - tau])) for tau in range(tau_max + 1)], dtype=float)
    c_int = _cross_correlation_sum(w, h, tau_max) / np.maximum(valid, 1.0)
    # continuous: a bond present without interruption from t0 for `run` frames contributes to every tau < run
    run = np.zeros(p, dtype=np.int64)
    contrib = np.zeros(n + 2)
    for t in range(n - 1, -1, -1):
        run = np.where(h[t], run + 1, 0)
        nz = h[t]
        if nz.any():
            contrib += np.bincount(run[nz], weights=w[t, nz], minlength=n + 2)
    survive = np.cumsum(contrib[::-1])[::-1]
    c_cont = survive[1:tau_max + 2] / np.maximum(valid, 1.0)
    dt = float(dt_fs) if dt_fs else 1.0
    unit = "fs" if dt_fs else "frame"
    return {"n_frames": n, "n_pairs": p, "tau_max": tau_max, "dt_fs": dt_fs, "unit": unit,
            "lag": [float(x) for x in np.arange(tau_max + 1) * dt],
            "intermittent": [float(x) for x in c_int], "continuous": [float(x) for x in c_cont],
            "lifetime_intermittent": _lifetime_of(c_int, dt), "lifetime_continuous": _lifetime_of(c_cont, dt),
            "n_origins": [int(v) for v in valid],
            "memory": {"estimated_mb": estimate / 2 ** 20, "actual_mb": n * p / 2 ** 20, "budget_mb": float(budget_mb)},
            "criteria": {"donor_acceptor_A": float(donor_acceptor_A), "angle_deg": float(angle_deg), "donors": list(donors)},
            "source": LIFETIME_SOURCE,
            "definition": L(
                "C(τ) = ⟨ Σ_ij h_ij(t0) h_ij(t0+τ) / Σ_ij h_ij(t0) ⟩ (t0 で平均。h_ij は水素 i と受容体 j の組がその時刻に条件を満たせば 1)。"
                "intermittent は途中で切れて戻っても数え、continuous は t0 から τ まで切れずに続いた組だけを数える "
                "(MDAnalysis の HydrogenBondAnalysis.lifetime の intermittency = ∞ と 0 に当たる)。"
                f"寿命は C(τ) の積分 (τ_max = {tau_max} {unit} で打ち切り) と、C が 1/e を切る τ。指数関数の当てはめはしていません",
                "C(tau) = < sum_ij h_ij(t0) h_ij(t0+tau) / sum_ij h_ij(t0) > averaged over t0, where h_ij is 1 when hydrogen i and acceptor j "
                "satisfy the criteria at that time. Intermittent counts a bond that breaks and re-forms; continuous counts only bonds present "
                "without interruption from t0 to t0+tau (intermittency = infinity and 0 of MDAnalysis HydrogenBondAnalysis.lifetime). "
                f"The lifetime is the integral of C(tau) (truncated at tau_max = {tau_max} {unit}) and the tau where C falls below 1/e; no exponential fit")}


def distance_angle_map(frames, rmax_A: float, donors=DEFAULT_DONORS, bins=MAP_BINS, hydrogen: str = "H", max_dh_A: float = 1.3) -> dict:
    """2D histogram of the D-A distance (0..rmax_A) and the D-H...A angle (0..180 deg) over all D-H...A triplets of every frame."""
    from adit.analysis.free_energy import Distribution2D

    if rmax_A <= 0:
        raise HydrogenBondError(L("距離の上限 [Å] は正の値にしてください", "the upper distance in Å must be positive"))
    ex = np.linspace(0.0, float(rmax_A), int(bins[0]) + 1)
    ey = np.linspace(0.0, 180.0, int(bins[1]) + 1)
    counts = np.zeros((int(bins[0]), int(bins[1])))
    n_frames = n_points = 0
    for fr in frames:
        _, _, _, r, ang = _geometry(fr, donors, hydrogen, max_dh_A, rmax_A)
        if len(r):
            c, _, _ = np.histogram2d(r, ang, bins=[ex, ey])
            counts += c
        n_frames += 1
        n_points += int(len(r))
    area = np.outer(np.diff(ex), np.diff(ey))
    total = counts.sum()
    density = counts / (total * area) if total > 0 else np.zeros_like(counts)
    dist = Distribution2D(0.5 * (ex[:-1] + ex[1:]), 0.5 * (ey[:-1] + ey[1:]), counts, density)
    peak = np.unravel_index(int(np.argmax(counts)), counts.shape) if total > 0 else None
    return {"distribution": dist, "n_frames": n_frames, "n_points": n_points, "rmax_A": float(rmax_A), "bins": [int(b) for b in bins],
            "donors": list(donors), "max_dh_A": float(max_dh_A),
            "peak": None if peak is None else {"distance_A": float(dist.x_centers[peak[0]]), "angle_deg": float(dist.y_centers[peak[1]])},
            "note": L("D–H···A の全部の組 (距離の上限以内、角度は制限なし) の 2 次元ヒストグラム。密度は面積で割って全体を 1 にした値。"
                      "しきい値の候補を見るためのもので、ADIT はしきい値を決めません",
                      "2D histogram of every D-H...A triplet within the upper distance, any angle; the density is normalised to 1 over the area. "
                      "It is meant for choosing the thresholds; ADIT does not choose them")}


def write_lifetime_csv(t: dict, path: Path | str) -> Path:
    path = Path(path)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow([f"lag_{t['unit']}", "c_intermittent", "c_continuous", "n_origins"])
        for lag, ci, cc, n in zip(t["lag"], t["intermittent"], t["continuous"], t["n_origins"]):
            w.writerow([f"{lag:g}", f"{ci:.6f}", f"{cc:.6f}", n])
    return path


def write_map_csv(t: dict, path: Path | str) -> Path:
    path = Path(path)
    d = t["distribution"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["distance_A", "angle_deg", "counts", "density_per_A_deg"])
        for i, x in enumerate(d.x_centers):
            for j, y in enumerate(d.y_centers):
                w.writerow([f"{x:.4f}", f"{y:.2f}", int(d.counts[i, j]), f"{d.density[i, j]:.6g}"])
    return path


__all__ = ["HydrogenBondError", "KNOWN_CRITERIA", "DEFAULT_DONORS", "criteria_note", "count_frame", "count_series",
           "lifetime", "distance_angle_map", "write_lifetime_csv", "write_map_csv"]
