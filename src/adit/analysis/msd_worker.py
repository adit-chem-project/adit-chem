#!/usr/bin/env python3

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path

import numpy as np

SHIFTS = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], dtype=float)
CHUNK_BYTES = 32 * 2 ** 20
NEAR_CAP = 0.9
CAP_FRACTION = 0.01


def mic(vectors: np.ndarray, cell: np.ndarray) -> np.ndarray:
    cell = np.asarray(cell, dtype=float)
    shape = vectors.shape
    flat = vectors.reshape(-1, 3)
    frac = np.linalg.solve(cell.T, flat.T).T
    frac -= np.round(frac)
    cand = (frac[:, None, :] + SHIFTS[None, :, :]) @ cell
    best = np.argmin(np.einsum("nij,nij->ni", cand, cand), axis=1)
    return cand[np.arange(len(flat)), best].reshape(shape)


def cell_half_width(cell) -> float | None:
    if cell is None:
        return None
    cell = np.asarray(cell, dtype=float)
    vol = abs(float(np.linalg.det(cell)))
    if vol <= 0:
        return None
    widths = [vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)]
    return float(min(widths) / 2.0)


def unwrap(wrapped: np.ndarray, cell) -> np.ndarray:
    if cell is None:
        return np.asarray(wrapped, dtype=float)
    steps = mic(np.diff(wrapped, axis=0), cell)
    out = np.empty_like(wrapped, dtype=float)
    out[0] = wrapped[0]
    np.cumsum(steps, axis=0, out=out[1:])
    out[1:] += wrapped[0]
    return out


# ---------------- MSD ----------------
def msd_per_atom(unwrapped: np.ndarray, chunk_bytes: int = CHUNK_BYTES) -> np.ndarray:
    unwrapped = np.asarray(unwrapped, dtype=float)
    frames, natoms, _ = unwrapped.shape
    size = 1 << (2 * frames - 1).bit_length()
    out = np.empty((frames, natoms))
    per_chunk = max(1, int(chunk_bytes // max(1, size * 3 * 16)))
    counts = (frames - np.arange(frames))[:, None]
    for start in range(0, natoms, per_chunk):
        x = unwrapped[:, start:start + per_chunk, :]
        k = x.shape[1]
        sq = np.einsum("tac,tac->ta", x, x)
        total = np.sum(sq, axis=0)
        left = total - np.concatenate([np.zeros((1, k)), np.cumsum(sq, axis=0)[:-1]])
        right = total - np.concatenate([np.zeros((1, k)), np.cumsum(sq[::-1], axis=0)[:-1]])
        f = np.fft.rfft(x, n=size, axis=0)
        corr = np.sum(np.fft.irfft(f * np.conjugate(f), n=size, axis=0)[:frames], axis=2)
        out[:, start:start + k] = (left + right) / counts - 2.0 * corr / counts
    out[0] = 0.0
    return out


def fit_slope_d(times: np.ndarray, msd: np.ndarray, lo: float, hi: float, dim: int = 3) -> float | None:
    sel = (times >= lo) & (times <= hi)
    if int(sel.sum()) < 2:
        return None
    x, y = times[sel], msd[sel]
    n = len(x)
    denom = (x ** 2).sum() - x.sum() ** 2 / n
    if denom <= 0:
        return None
    return float(((x * y).sum() - y.sum() * x.sum() / n) / denom / (2 * dim))


def fit_gaussian_d(centers: np.ndarray, counts: np.ndarray, tau_fs: float, guess: float) -> float:
    from scipy.optimize import least_squares

    if len(centers) < 2 or counts.sum() <= 0 or tau_fs <= 0:
        return float("nan")
    dr = float(centers[1] - centers[0])
    n = float(counts.sum())
    r2 = centers ** 2

    def model(d: float) -> np.ndarray:
        return n * dr * r2 * np.exp(-r2 / (4 * d * tau_fs)) / (2 * np.sqrt(math.pi) * (d * tau_fs) ** 1.5)

    def jac(p):
        d = p[0]
        return (-model(d) * (r2 / (4 * d * d * tau_fs) - 1.5 / d)).reshape(-1, 1)

    start = max(float(guess), 1e-12)
    res = least_squares(lambda p: counts - model(p[0]), [start], jac=jac, bounds=(1e-30, np.inf))
    return float(res.x[0])


def moments(pos: np.ndarray, cell, tau: int, use_mic: bool, chunk_bytes: int = CHUNK_BYTES, cap: float | None = None):
    frames, natoms, _ = pos.shape
    origins = frames - tau
    per_chunk = max(1, int(chunk_bytes // max(1, natoms * 3 * 8 * (27 if use_mic else 1))))
    s2 = s4 = 0.0
    count = near = 0
    parts = []
    for start in range(0, origins, per_chunk):
        stop = min(start + per_chunk, origins)
        step = pos[start + tau:stop + tau] - pos[start:stop]
        if use_mic:
            step = mic(step, cell)
        r2 = np.einsum("tac,tac->ta", step, step)
        s2 += float(r2.sum()); s4 += float((r2 ** 2).sum()); count += r2.size
        if cap:
            near += int((r2 > (NEAR_CAP * cap) ** 2).sum())
        parts.append(np.sqrt(r2).ravel())
    return s2 / count, s4 / count, np.concatenate(parts) if parts else np.zeros(0), (near / count if count else 0.0)


def van_hove_self(pos: np.ndarray, cell, taus, dt_fs: float, *, displacement: str = "mic",
                  chunk_bytes: int = CHUNK_BYTES) -> dict:
    pos = np.asarray(pos, dtype=float)
    use_mic = displacement == "mic" and cell is not None
    half = cell_half_width(cell) if use_mic else None
    times, d_fit, d_direct, alpha2, rms = [], [], [], [], []
    guess = None
    truncated = None
    near_fraction = 0.0
    for tau in taus:
        tau = int(tau)
        if not 0 < tau < pos.shape[0]:
            continue
        tau_fs = tau * dt_fs
        m2, m4, r, near = moments(pos, cell, tau, use_mic, chunk_bytes, half)
        direct = m2 / (6 * tau_fs)
        guess = guess if guess is not None else direct
        counts, edges = np.histogram(r, bins="fd")
        centers = (edges[:-1] + edges[1:]) / 2.0
        times.append(tau_fs)
        d_fit.append(fit_gaussian_d(centers, counts, tau_fs, guess))
        d_direct.append(direct)
        alpha2.append(3.0 * m4 / (5.0 * m2 ** 2) - 1.0 if m2 > 0 else float("nan"))
        rms.append(math.sqrt(m2))
        if use_mic and truncated is None and near > CAP_FRACTION:
            truncated, near_fraction = tau_fs, near
    return {"times_fs": times, "d_fit_A2_fs": d_fit, "d_direct_A2_fs": d_direct, "alpha2": alpha2, "rms_A": rms,
            "displacement": "mic" if use_mic else "unwrapped", "half_width_A": half,
            "truncated_from_fs": truncated, "near_cap_fraction": near_fraction}


def default_taus(frames: int, count: int = 100) -> list[int]:
    top = max(2, frames // 2)
    return sorted({int(v) for v in np.linspace(1, top, min(count, top)) if v >= 1})


def read_cell(path: Path):
    rows = [list(map(float, line.split()[1:4])) for line in Path(path).read_text(encoding="utf-8").splitlines() if "TV" in line]
    if len(rows) != 3:
        raise SystemExit(f"{path} に TV の行が 3 本ありません ({len(rows)} 本)")
    return np.array(rows, dtype=float)


def read_table(path: Path, natoms: int) -> np.ndarray:
    flat = np.loadtxt(path, usecols=(1, 2, 3), dtype=float)
    frames, rest = divmod(len(flat), natoms)
    if rest:
        raise SystemExit(f"行数 {len(flat)} が原子数 {natoms} で割り切れません (最後のフレームが欠けています)")
    return flat.reshape(frames, natoms, 3)


def read_with_ase(path: Path, index: str = ":"):
    from ase.io import read

    frames = read(str(path), index=index)
    frames = frames if isinstance(frames, list) else [frames]
    pos = np.array([f.get_positions() for f in frames], dtype=float)
    cell = np.asarray(frames[-1].cell, dtype=float) if any(frames[-1].pbc) else None
    return pos, frames[0].get_chemical_symbols(), cell


def run(pos: np.ndarray, symbols, cell, dt_fs: float, *, species: str | None = None, taus=None,
        displacement: str = "mic", fit_fraction=(0.1, 0.5), chunk_bytes: int = CHUNK_BYTES) -> dict:
    t0 = time.perf_counter()
    pos = np.asarray(pos, dtype=float)
    symbols = list(symbols) if symbols is not None else [""] * pos.shape[1]
    idx = [i for i, s in enumerate(symbols) if s == species] if species else list(range(pos.shape[1]))
    if not idx:
        raise SystemExit(f"元素 {species} の原子がありません")
    unwrapped = unwrap(pos, cell)[:, idx, :]
    msd = msd_per_atom(unwrapped, chunk_bytes)
    times = np.arange(len(msd), dtype=float) * dt_fs
    lo, hi = times[-1] * fit_fraction[0], times[-1] * fit_fraction[1]
    per_atom = [fit_slope_d(times, msd[:, a], lo, hi) for a in range(msd.shape[1])]
    mean_msd = msd.mean(axis=1)
    t_msd = time.perf_counter() - t0
    vh = van_hove_self(pos[:, idx, :] if displacement == "mic" else unwrapped, cell,
                       taus if taus is not None else default_taus(pos.shape[0]), dt_fs,
                       displacement=displacement, chunk_bytes=chunk_bytes)
    good = [v for v in per_atom if v is not None]
    warnings = []
    if cell is None:
        warnings.append("セルが無いので座標を巻き戻していません (折り返された周期系の軌跡なら --cell でセルを渡してください。"
                        "変位の分布も巻き戻していない座標から求めています)")
    return {
        "frames": int(pos.shape[0]), "atoms": len(idx), "species": species, "dt_fs": dt_fs,
        "cell_A": None if cell is None else np.asarray(cell, dtype=float).tolist(), "warnings": warnings,
        "symbols": [symbols[i] for i in idx], "atom_index": idx,
        "fit_range_fs": [lo, hi], "fit_fraction": list(fit_fraction),
        "msd": {"lag_fs": times.tolist(), "mean_A2": mean_msd.tolist(),
                "per_atom_A2_last": msd[-1].tolist()},
        "d_per_atom_cm2_s": [None if v is None else v * 1e-1 for v in per_atom],
        "d_mean_cm2_s": (sum(good) / len(good) * 1e-1) if good else None,
        "d_from_mean_msd_cm2_s": (lambda v: None if v is None else v * 1e-1)(fit_slope_d(times, mean_msd, lo, hi)),
        "vanhove": vh,
        "seconds": {"msd": t_msd, "total": time.perf_counter() - t0},
        "units": "D は cm²/s、MSD は Å²、時間は fs (Å²/fs → cm²/s は ×0.1)",
        "note": ("時間原点を全部使った MSD と、変位の分布から出した D。どちらも判定はしていません。"
                 "変位の分布の D は、ガウス (フィックの拡散) なら MSD の傾きと一致します。"),
    }


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="重い MSD と変位の分布を計算して JSON に書く (ADIT を import しない単体のファイル)")
    ap.add_argument("trajectory", type=Path)
    ap.add_argument("--natoms", type=int, default=0, help="表として読むときの 1 フレームの原子数 (省くと ASE で読む)")
    ap.add_argument("--cell", type=Path, default=None, help="TV の行があるファイル (表として読むとき、または軌跡にセルが無いときに指定)")
    ap.add_argument("--dt", type=float, required=True, help="フレームの間隔 [fs]")
    ap.add_argument("--species", default=None, help="元素を 1 つに絞る (例 Na)")
    ap.add_argument("--taus", type=int, default=100, help="変位の分布を見る遅れ時間の点数 (既定 100)")
    ap.add_argument("--displacement", choices=("mic", "unwrapped"), default="mic")
    ap.add_argument("--out", type=Path, default=Path("msd_vanhove.json"))
    a = ap.parse_args(argv)

    if a.natoms:
        pos = read_table(a.trajectory, a.natoms)
        symbols = None
        cell = read_cell(a.cell) if a.cell else None
    else:
        pos, symbols, cell = read_with_ase(a.trajectory)
        if a.cell:
            cell = read_cell(a.cell)
    result = run(pos, symbols, cell, a.dt, species=a.species, taus=default_taus(pos.shape[0], a.taus),
                 displacement=a.displacement)
    a.out.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"{result['frames']} フレーム × {result['atoms']} 原子 を {result['seconds']['total']:.1f} 秒で処理しました")
    for w in result["warnings"]:
        print(f"  注意: {w}")
    print(f"  D (原子ごとの平均): {result['d_mean_cm2_s']:.4g} cm²/s" if result["d_mean_cm2_s"] else "  D: 求められません")
    print(f"  書き出し: {a.out}  (図にするのは adit-analyze か画面の「解析」です)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
