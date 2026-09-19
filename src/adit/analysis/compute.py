
from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers, covalent_radii
AMU_G = 1.66053906660e-24  # 1 u [g]


class RDFAccumulator:

    def __init__(self, pairs: list[tuple[str, str]], rmax: float = 8.0, nbins: int = 160):
        self.pairs = list(pairs)
        self.rmax = float(rmax)
        self.edges = np.linspace(0.0, self.rmax, nbins + 1)
        self.shell = 4.0 / 3.0 * np.pi * (self.edges[1:] ** 3 - self.edges[:-1] ** 3)
        self.g = {p: np.zeros(nbins) for p in self.pairs}
        self.n = {p: np.zeros(nbins) for p in self.pairs}
        self.n_rev = {p: np.zeros(nbins) for p in self.pairs}
        self.dens = {p: np.zeros(nbins) for p in self.pairs}
        self.dens_rev = {p: np.zeros(nbins) for p in self.pairs}
        self.frames = {p: 0 for p in self.pairs}
        self.n_frames = 0
        self.vol_sum = 0.0
        self.vol_frames = 0
        self.periodic_frames = 0

    def add(self, fr: Atoms) -> None:
        self.n_frames += 1
        sym = np.array(fr.get_chemical_symbols())
        periodic = bool(any(fr.pbc))
        need = {e for p in self.pairs for e in p}
        mask = np.isin(sym, list(need))
        sub = fr[mask] if not mask.all() else fr
        ssym = sym[mask]
        if len(sub) < 2:
            return
        i, j, d = pairs_within(sub, self.rmax)
        vol = fr.get_volume() if periodic else 4.0 / 3.0 * np.pi * self.rmax**3
        self.vol_sum += float(vol)
        self.vol_frames += 1
        self.periodic_frames += int(periodic)
        si, sj = ssym[i], ssym[j]
        for a, b in self.pairs:
            na, nb = int(np.sum(ssym == a)), int(np.sum(ssym == b))
            if na == 0 or nb == 0 or (a == b and na < 2):
                continue
            if a == b:
                sel = (si == a) & (sj == a)
                hist = 2.0 * np.histogram(d[sel], bins=self.edges)[0]
                npairs = na * (na - 1)
            else:
                sel = ((si == a) & (sj == b)) | ((si == b) & (sj == a))
                hist = np.histogram(d[sel], bins=self.edges)[0].astype(float)
                npairs = na * nb
            self.g[(a, b)] += hist * vol / (npairs * self.shell)
            self.n[(a, b)] += np.cumsum(hist) / na
            self.n_rev[(a, b)] += np.cumsum(hist) / (na if a == b else nb)
            self.dens[(a, b)] += hist / (na * self.shell)
            self.dens_rev[(a, b)] += hist / ((na if a == b else nb) * self.shell)
            self.frames[(a, b)] += 1

    def result(self, pair: tuple[str, str]) -> dict:
        r = 0.5 * (self.edges[1:] + self.edges[:-1])
        k = max(self.frames[pair], 1)
        return {"r": r, "g": self.g[pair] / k, "n": self.n[pair] / k, "n_reverse": self.n_rev[pair] / k, "r_edges_upper": self.edges[1:],
                "density": self.dens[pair] / k, "density_reverse": self.dens_rev[pair] / k, "n_frames": self.frames[pair]}

    def normalization(self) -> dict:
        periodic = self.vol_frames > 0 and self.periodic_frames == self.vol_frames
        mixed = 0 < self.periodic_frames < self.vol_frames
        return {"kind": ("cell_volume" if periodic else "mixed" if mixed else "sphere_of_rmax"),
                "mean_volume_A3": (self.vol_sum / self.vol_frames) if self.vol_frames else None,
                "rmax_A": self.rmax, "depends_on_rmax": not periodic,
                "periodic_frames": self.periodic_frames, "n_frames": self.vol_frames}


def pairs_within(atoms: Atoms, rmax: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    from scipy.spatial import cKDTree

    pos = atoms.get_positions()
    n = len(pos)
    empty = (np.zeros(0, int), np.zeros(0, int), np.zeros(0))
    if n < 2:
        return empty
    pbc = np.asarray(atoms.pbc, dtype=bool)
    if not pbc.any() or atoms.cell.rank < 3:
        pr = cKDTree(pos).query_pairs(rmax, output_type="ndarray")
        i, j = pr[:, 0], pr[:, 1]
        return i, j, np.linalg.norm(pos[i] - pos[j], axis=1)
    cell = np.asarray(atoms.cell, dtype=float)
    frac = np.linalg.solve(cell.T, pos.T).T
    frac[:, pbc] %= 1.0
    p = frac @ cell
    shifts = np.array([(a, b, c) for a in ((-1, 0, 1) if pbc[0] else (0,)) for b in ((-1, 0, 1) if pbc[1] else (0,))
                       for c in ((-1, 0, 1) if pbc[2] else (0,))], dtype=float)
    img = (p[None, :, :] + (shifts @ cell)[:, None, :]).reshape(-1, 3)
    sdm = cKDTree(p).sparse_distance_matrix(cKDTree(img), rmax, output_type="ndarray")
    i = sdm["i"].astype(int)
    j = (sdm["j"] % n).astype(int)
    d = sdm["v"].astype(float)
    keep = j > i
    i, j, d = i[keep], j[keep], d[keep]
    if len(i) == 0:
        return empty
    order = np.lexsort((d, j, i))
    i, j, d = i[order], j[order], d[order]
    first = np.ones(len(i), bool)
    first[1:] = (i[1:] != i[:-1]) | (j[1:] != j[:-1])
    return i[first], j[first], d[first]


def bond_pairs(atoms: Atoms, factor: float = 1.2) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nums = atoms.get_atomic_numbers()
    if len(nums) < 2:
        return np.zeros(0, int), np.zeros(0, int), np.zeros(0)
    rad = covalent_radii[nums]
    i, j, d = pairs_within(atoms, 2.0 * factor * float(rad.max()))
    keep = d <= factor * (rad[i] + rad[j])
    return i[keep], j[keep], d[keep]


def rdf(frames: Iterable[Atoms], pair: tuple[str, str], rmax: float = 8.0, nbins: int = 160) -> tuple[np.ndarray, np.ndarray]:
    acc = RDFAccumulator([pair], rmax, nbins)
    for fr in frames:
        acc.add(fr)
    out = acc.result(pair)
    return out["r"], out["g"]


# ---------------- MSD ----------------
_MIC_SHIFTS = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)], dtype=float)


def _mic_step(step: np.ndarray, cell: np.ndarray) -> np.ndarray:
    frac = np.linalg.solve(cell.T, step.T).T
    frac -= np.round(frac)
    cand = (frac[:, None, :] + _MIC_SHIFTS[None, :, :]) @ cell
    best = np.argmin(np.einsum("nij,nij->ni", cand, cand), axis=1)
    return cand[np.arange(len(frac)), best]


class UnwrapAccumulator:

    def __init__(self, n_frames: int | None = None):
        self.n_frames = n_frames
        self.buf: np.ndarray | None = None
        self.rows: list[np.ndarray] = []
        self.k = 0
        self.prev = None
        self.cur = None
        self.symbols: list[str] = []
        self.max_step_fraction = 0.0
        self.large_step_count = 0

    def add(self, fr: Atoms) -> None:
        p = fr.get_positions()
        if self.prev is None:
            self.cur = p.copy()
            self.symbols = fr.get_chemical_symbols()
            if self.n_frames:
                self.buf = np.empty((self.n_frames, len(p), 3))
        else:
            step = p - self.prev
            if any(fr.pbc):
                step = _mic_step(step, np.asarray(fr.cell, dtype=float))
                cell = np.asarray(fr.cell, dtype=float)
                vol = abs(float(np.linalg.det(cell)))
                widths = [vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)]
                shortest = min(widths)
                if shortest > 0 and len(step):
                    ratio = float(np.linalg.norm(step, axis=1).max() / shortest)
                    self.max_step_fraction = max(self.max_step_fraction, ratio)
                    if ratio >= 0.4:
                        self.large_step_count += 1
            self.cur = self.cur + step
        if self.buf is not None and self.k < len(self.buf):
            self.buf[self.k] = self.cur
        else:
            self.rows.append(self.cur)
        self.k += 1
        self.prev = p

    def result(self) -> tuple[np.ndarray, list[str]]:
        if self.k == 0:
            return np.zeros((0, 0, 3)), []
        if self.buf is not None and not self.rows:
            return self.buf[: self.k], self.symbols
        head = [self.buf[: min(self.k, len(self.buf))]] if self.buf is not None else []
        return np.concatenate(head + [np.stack(self.rows)]) if head else np.stack(self.rows), self.symbols


def unwrapped_positions(frames: Iterable[Atoms], idx: np.ndarray | None = None) -> tuple[np.ndarray, list[str]]:
    acc = UnwrapAccumulator(len(frames) if hasattr(frames, "__len__") else None)
    for fr in frames:
        acc.add(fr)
    pos, syms = acc.result()
    if idx is not None and pos.size:
        pos, syms = pos[:, idx, :], [syms[i] for i in idx]
    return pos, syms


def _autocorr_fft(x: np.ndarray) -> np.ndarray:
    n = x.shape[0]
    f = np.fft.rfft(x, n=2 * n, axis=0)
    res = np.fft.irfft(f * f.conj(), axis=0)[:n]
    return res / (n - np.arange(n))[:, None]


def remove_com_drift(pos: np.ndarray, symbols: list[str]) -> tuple[np.ndarray, dict]:
    from ase.data import atomic_masses

    T = pos.shape[0]
    empty = {"removed": False, "displacement_A": None, "max_displacement_A": None, "n_atoms": int(pos.shape[1]) if pos.size else 0,
             "mass_weighted": True, "reason": "no frames"}
    if T == 0 or pos.shape[1] == 0:
        return pos, empty
    if symbols and len(symbols) == pos.shape[1]:
        w = np.array([atomic_masses[atomic_numbers[s]] for s in symbols], dtype=float)
        weighted = True
    else:
        w = np.ones(pos.shape[1], dtype=float)
        weighted = False
    com = np.einsum("tnc,n->tc", pos, w) / w.sum()  # (T, 3)
    d = np.linalg.norm(com - com[0], axis=1)
    info = {"removed": True, "displacement_A": float(d[-1]), "max_displacement_A": float(d.max()),
            "n_atoms": int(pos.shape[1]), "mass_weighted": bool(weighted)}
    return pos - com[:, None, :], info


AXIS_INDEX = {"x": 0, "y": 1, "z": 2}


def parse_axes(axes: str) -> tuple[int, ...]:
    a = "".join(dict.fromkeys(str(axes).strip().lower().replace(",", "")))
    if not a or any(c not in AXIS_INDEX for c in a):
        raise ValueError(f"MSD の成分は x / y / z の組み合わせで指定します (受け取った値: {axes!r})")
    return tuple(AXIS_INDEX[c] for c in a)


def msd_formula(dim: int) -> str:
    return f"MSD = {2 * dim} D t + c"


def msd_fft(pos: np.ndarray, chunk_bytes: int = 64 * 2**20) -> np.ndarray:
    T, n, nc = pos.shape
    if T == 0 or n == 0 or nc == 0:
        return np.zeros(T)
    total = np.zeros(T)
    per_atom = max(1, int(chunk_bytes // max(1, T * nc * 16 * 3)))
    m = np.arange(T)
    for s in range(0, n, per_atom):
        x = pos[:, s:s + per_atom, :]
        d = np.sum(x**2, axis=2)  # (T, k)
        front = np.vstack([np.zeros((1, d.shape[1])), np.cumsum(d, axis=0)[:-1]])  # Σ_{k<m} D[k]
        back = np.vstack([np.zeros((1, d.shape[1])), np.cumsum(d[::-1], axis=0)[:-1]])  # Σ_{k≥T−m} D[k]
        s1 = (2 * d.sum(axis=0)[None, :] - front - back) / (T - m)[:, None]
        s2 = sum(_autocorr_fft(x[:, :, c]) for c in range(nc))
        total += np.sum(s1 - 2 * s2, axis=1)
    return total / n


DEFAULT_FIT_FRACTION = (0.1, 0.5)


def fit_range_fs(t: np.ndarray, fit_fs: tuple[float, float] | None, frac: tuple[float, float] = DEFAULT_FIT_FRACTION) -> tuple[float, float]:
    if fit_fs is not None:
        return float(fit_fs[0]), float(fit_fs[1])
    if len(t) < 2:
        return float(t[0]) if len(t) else 0.0, float(t[-1]) if len(t) else 0.0
    lo, hi = float(frac[0]) * float(t[-1]), float(frac[1]) * float(t[-1])
    if np.count_nonzero((t >= lo - 1e-9) & (t <= hi + 1e-9)) < 2:
        lo, hi = float(t[1]), float(t[-1])
    return lo, hi


def diffusion_fit(t: np.ndarray, m: np.ndarray, rng: tuple[float, float], dim: int = 3) -> float | None:
    sel = (t >= rng[0] - 1e-9) & (t <= rng[1] + 1e-9)
    if np.count_nonzero(sel) < 2:
        return None
    slope = np.polyfit(t[sel], m[sel], 1)[0]  # Å²/fs
    return float(slope / (2.0 * dim) * 1e-16 / 1e-15)  # Å²/fs → cm^2/s


def loglog_slope(t: np.ndarray, m: np.ndarray, rng: tuple[float, float]) -> float | None:
    sel = (t >= rng[0] - 1e-9) & (t <= rng[1] + 1e-9) & (t > 0) & (m > 0)
    if np.count_nonzero(sel) < 2:
        return None
    return float(np.polyfit(np.log(t[sel]), np.log(m[sel]), 1)[0])


def diffusion_blocks(pos: np.ndarray, dt_fs: float, dim: int, frac: tuple[float, float] = DEFAULT_FIT_FRACTION,
                     n_blocks: int = 5, min_block: int = 4, *, fit_fs: tuple[float, float] | None = None) -> dict:
    from ..lang import L

    T = pos.shape[0]
    rng = fit_range_fs(np.arange(T, dtype=float) * dt_fs, fit_fs, frac)
    k = int(min(n_blocks, T // min_block))
    out = {"d_err_cm2_s": None, "n_blocks": max(0, k), "d_blocks_cm2_s": [], "block_frames": 0,
           "block_frame_counts": [], "blocks": [], "fit_range_fs": list(rng), "n_frames_used": 0,
           "d_blocks_mean_cm2_s": None, "reason": None, "reason_code": None,
           "estimator": "standard_error_of_block_mean", "estimate_for": "mean_of_block_D"}
    if k < 2:
        out.update(n_blocks=0, reason_code="insufficient_blocks", reason=L(
            f"ブロックを 2 つ以上 (1 ブロック {min_block} フレーム以上) 取れません (使ったフレーム {T})",
            f"Cannot form at least two blocks of {min_block} or more frames ({T} frames available)"))
        return out
    requested = k
    while k >= 2:
        length_k, extra_k = divmod(T, k)
        shortest = length_k
        if rng[1] <= (shortest - 1) * dt_fs + 1e-9 and rng[0] >= -1e-9:
            break
        k -= 1
    if k < 2:
        length_k = T // 2
        span = (length_k - 1) * dt_fs
        hint = ""
        if span > rng[0] + 1e-9 and length_k >= min_block:
            hint = L(f" ブロックを 2 つ取るなら {rng[0]:g}〜{span:g} fs までが入ります (--msd-fit {rng[0]:g},{span:g})。"
                     "ただし上の D とは当てはめ範囲が変わります",
                     f" With two blocks, {rng[0]:g}-{span:g} fs would fit (--msd-fit {rng[0]:g},{span:g}), "
                     "but that is a different fit range from the D above")
        out.update(n_blocks=0, reason_code="fit_range_not_available_in_all_blocks", reason=L(
            f"主当てはめ範囲 {rng[0]:g}〜{rng[1]:g} fs を含むブロックを 2 つ取れません (使ったフレーム {T})。",
            f"cannot form two blocks that both cover the main fit range {rng[0]:g}-{rng[1]:g} fs ({T} frames). ") + hint)
        return out
    out["n_blocks"] = k
    if k != requested:
        out["n_blocks_requested"] = requested
        out["adjusted"] = L(
            f"ブロックの数を {requested} から {k} に減らしました (当てはめ範囲 {rng[0]:g}〜{rng[1]:g} fs が全ブロックに収まる最大の数)",
            f"the number of blocks was reduced from {requested} to {k} (the largest number for which every block covers "
            f"the fit range {rng[0]:g}-{rng[1]:g} fs)")
    length, extra = divmod(T, k)
    start = 0
    for b in range(k):
        size = length + int(b < extra)
        tb = np.arange(size, dtype=float) * dt_fs
        sel = (tb >= rng[0] - 1e-9) & (tb <= rng[1] + 1e-9)
        out["blocks"].append({"frame_start": start, "frame_stop": start + size, "n_frames": size,
                              "duration_fs": float(tb[-1]), "fit_range_fs": list(rng),
                              "fit_points": int(sel.sum()), "D_cm2_s": None})
        start += size
    out.update(block_frame_counts=[b["n_frames"] for b in out["blocks"]], n_frames_used=T,
               block_frames=length if extra == 0 else None)
    if any(rng[0] < -1e-9 or rng[1] > b["duration_fs"] + 1e-9 or b["fit_points"] < 2 for b in out["blocks"]):
        out.update(reason_code="fit_range_not_available_in_all_blocks", reason=L(
            f"全ブロックで主当てはめ範囲 {rng[0]:g}〜{rng[1]:g} fs を含む2点以上を取れません。範囲は変更していません",
            f"Not every block covers the main fit range {rng[0]:g}-{rng[1]:g} fs with at least two points; the range was not changed"))
        return out
    for b in out["blocks"]:
        tb = np.arange(b["n_frames"], dtype=float) * dt_fs
        b["D_cm2_s"] = diffusion_fit(tb, msd_fft(pos[b["frame_start"]:b["frame_stop"]]), rng, dim)
    a = np.array([b["D_cm2_s"] for b in out["blocks"]], dtype=float)
    out.update(d_err_cm2_s=float(a.std(ddof=1) / np.sqrt(k)), d_blocks_cm2_s=a.tolist(),
               d_blocks_mean_cm2_s=float(a.mean()))
    return out


def msd(frames, species: str | None, times_fs: list[float] | None, *, symbols: list[str] | None = None, dt_fs: float | None = None,
        fit_fs: tuple[float, float] | None = None, axes: str = "xyz", remove_drift: bool = False,
        drift_info: dict | None = None) -> tuple[np.ndarray, np.ndarray, float | None]:
    out = msd_analysis(frames, species, times_fs, symbols=symbols, dt_fs=dt_fs, fit_fs=fit_fs, axes=axes,
                       remove_drift=remove_drift, error=False)
    if drift_info is not None:
        drift_info.update(out["drift"])
    return out["lag"], out["msd_A2"], out["D_cm2_s"]


def msd_analysis(frames, species: str | None, times_fs: list[float] | None, *, symbols: list[str] | None = None,
                 dt_fs: float | None = None, fit_fs: tuple[float, float] | None = None, axes: str = "xyz",
                 remove_drift: bool = True, error: bool = True, n_blocks: int = 5) -> dict:
    comps = parse_axes(axes)
    dim = len(comps)
    if isinstance(frames, np.ndarray):
        pos, syms = frames, list(symbols or [])
    else:
        pos, syms = unwrapped_positions(frames)
    base = {"species": species, "axes": "".join("xyz"[c] for c in comps), "dimension": dim, "formula": msd_formula(dim),
            "drift": {"removed": False, "displacement_A": None, "max_displacement_A": None, "requested": bool(remove_drift)}}
    if pos.shape[0] < 2:
        return {**base, "lag": np.array([0.0]), "msd_A2": np.array([0.0]), "D_cm2_s": None, "fit_range_fs": None,
                "fit_range_user": fit_fs is not None, "fit_fraction": None if fit_fs else list(DEFAULT_FIT_FRACTION),
                "loglog_slope": None, "error": None, "dt_fs": dt_fs,
                "reason": L("軌跡のフレームが 1 つしかないので MSD を計算できません",
                            "the trajectory has only one frame, so no MSD can be computed")}
    if remove_drift:
        pos, info = remove_com_drift(pos, syms)
        base["drift"] = {**info, "requested": True}
    idx = np.where(np.array(syms) == species)[0] if species else np.arange(pos.shape[1])
    sub = pos[:, idx, :][:, :, comps]
    m = msd_fft(sub)
    dt = dt_fs if dt_fs else (times_fs[1] - times_fs[0] if times_fs and len(times_fs) >= 2 else None)
    t = np.arange(len(m), dtype=float) * (dt if dt else 1.0)
    D = rng = slope = err = None
    if dt and len(m) >= 4:
        rng = fit_range_fs(t, fit_fs)
        D = diffusion_fit(t, m, rng, dim)
        slope = loglog_slope(t, m, rng)
        if error:
            err = diffusion_blocks(sub, dt, dim, n_blocks=n_blocks, fit_fs=rng)
    return {**base, "lag": t, "msd_A2": m, "D_cm2_s": D, "fit_range_fs": list(rng) if rng else None,
            "fit_range_user": fit_fs is not None, "fit_fraction": None if fit_fs else list(DEFAULT_FIT_FRACTION),
            "loglog_slope": slope, "error": err, "dt_fs": dt}


def msd_per_atom(pos: np.ndarray, dt_fs: float, dim: int, fit_fs: tuple[float, float] | None = None,
                 chunk_bytes: int = 64 * 2 ** 20) -> dict:
    T, n, nc = pos.shape
    if T < 4 or n == 0:
        return {"d_cm2_s": [], "msd_A2": np.zeros((max(T, 1), 0)), "fit_range_fs": None}
    m = np.arange(T)
    counts = (T - m)[:, None]
    out = np.empty((T, n))
    per_atom = max(1, int(chunk_bytes // max(1, (1 << (2 * T - 1).bit_length()) * nc * 16)))
    for start in range(0, n, per_atom):
        x = pos[:, start:start + per_atom, :]
        d = np.sum(x ** 2, axis=2)
        front = np.vstack([np.zeros((1, d.shape[1])), np.cumsum(d, axis=0)[:-1]])
        back = np.vstack([np.zeros((1, d.shape[1])), np.cumsum(d[::-1], axis=0)[:-1]])
        s1 = (2 * d.sum(axis=0)[None, :] - front - back) / counts
        s2 = sum(_autocorr_fft(x[:, :, c]) for c in range(nc))
        out[:, start:start + x.shape[1]] = s1 - 2 * s2
    t = np.arange(T, dtype=float) * dt_fs
    rng = fit_range_fs(t, fit_fs)
    values = [diffusion_fit(t, out[:, a], rng, dim) for a in range(n)]
    return {"d_cm2_s": values, "msd_A2": out, "fit_range_fs": list(rng) if rng else None,
            "spread": {"min": min(v for v in values if v is not None) if any(v is not None for v in values) else None,
                       "max": max(v for v in values if v is not None) if any(v is not None for v in values) else None}}


def block_average(x: np.ndarray) -> list[dict]:
    x = np.asarray(x, dtype=float)
    out = []
    size = 1
    while len(x) >= 2:
        n = len(x)
        c0 = float(np.var(x))
        sem = np.sqrt(c0 / (n - 1))
        out.append({"block_size": size, "n_blocks": n, "sem": float(sem), "sem_err": float(sem / np.sqrt(2 * (n - 1)))})
        if n % 2:
            x = x[:-1]
        x = 0.5 * (x[0::2] + x[1::2])
        size *= 2
    return out


def autocorrelation_time(x: np.ndarray, c: float = 5.0) -> dict:
    x = np.asarray(x, dtype=float)
    n = len(x)
    if n < 4 or np.var(x) == 0:
        return {"tau_int": None, "window": None, "sem": None}
    y = x - x.mean()
    f = np.fft.rfft(y, n=2 * n)
    acf = np.fft.irfft(f * f.conj())[:n]
    rho = acf / acf[0]
    tau = 0.5 + np.cumsum(rho[1:])  # tau[M-1] = τ_int(M)
    M = np.arange(1, n)
    ok = np.where(M >= c * tau)[0]
    reached = ok.size > 0
    w = int(ok[0]) if reached else n - 2
    t_int = float(tau[w])
    sem = float(np.std(x) * np.sqrt(2 * t_int / n)) if t_int > 0 else None
    return {"tau_int": t_int, "window": int(M[w]), "sem": sem, "window_reached": bool(reached)}


def series_stats(x: list[float] | np.ndarray, dt_fs: float | None = None) -> dict:
    x = np.asarray(x, dtype=float)
    ac = autocorrelation_time(x)
    return {"n": int(len(x)), "mean": float(x.mean()), "std": float(x.std()), "blocks": block_average(x), "tau_int_samples": ac["tau_int"],
            "tau_int_fs": (ac["tau_int"] * dt_fs if ac["tau_int"] is not None and dt_fs else None), "acf_window": ac["window"],
            "acf_window_reached": ac.get("window_reached"), "sem_acf": ac["sem"], "dt_fs": dt_fs}


def density_g_cm3(atoms: Atoms) -> float | None:
    if not all(atoms.pbc) or atoms.cell.rank < 3:
        return None
    return float(atoms.get_masses().sum() * AMU_G / (atoms.get_volume() * 1e-24))


class ZDensityAccumulator:

    def __init__(self, bin_ang: float = 0.2, axis: int = 2):
        if axis not in (0, 1, 2):
            raise ValueError("axis は 0 (a) / 1 (b) / 2 (c) のどれか")
        self.bin_ang = float(bin_ang)
        self.axis = int(axis)
        self.nbins = None
        self.count: dict[str, np.ndarray] = {}
        self.mass = None
        self.d_sum = 0.0
        self.n_frames = 0
        self.c_perpendicular = True

    def add(self, fr: Atoms) -> None:
        cell = np.asarray(fr.cell, dtype=float)
        other = [i for i in (0, 1, 2) if i != self.axis]
        area_vec = np.cross(cell[other[0]], cell[other[1]])
        area = float(np.linalg.norm(area_vec))
        vol = abs(float(np.linalg.det(cell)))
        d = vol / area
        if self.nbins is None:
            self.nbins = max(1, int(round(d / self.bin_ang)))
            self.mass = np.zeros(self.nbins)
        nhat = area_vec / area
        axis_vec = cell[self.axis]
        if abs(abs(np.dot(axis_vec, nhat)) - np.linalg.norm(axis_vec)) > 1e-6 * np.linalg.norm(axis_vec):
            self.c_perpendicular = False
        s3 = np.linalg.solve(cell.T, fr.get_positions().T).T[:, self.axis] % 1.0
        b = np.minimum((s3 * self.nbins).astype(int), self.nbins - 1)
        slab = vol / self.nbins  # Å³
        sym = np.array(fr.get_chemical_symbols())
        for el in np.unique(sym):
            h = np.bincount(b[sym == el], minlength=self.nbins) / slab
            self.count.setdefault(el, np.zeros(self.nbins))
            self.count[el] += h
        m = np.bincount(b, weights=fr.get_masses(), minlength=self.nbins) / slab
        self.mass += m * AMU_G / 1e-24
        self.d_sum += d
        self.n_frames += 1

    def result(self) -> dict:
        k = max(self.n_frames, 1)
        d = self.d_sum / k
        z = (np.arange(self.nbins) + 0.5) * d / self.nbins
        return {"z_A": z, "number_density_A3": {el: v / k for el, v in sorted(self.count.items())}, "mass_density_g_cm3": self.mass / k,
                "axis_length_A": d, "bin_A": d / self.nbins, "n_frames": self.n_frames, "c_perpendicular": self.c_perpendicular,
                "axis": "abc"[self.axis]}


class MoleculeUnwrapper:

    def __init__(self, first: Atoms, factor: float = 1.2):
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import breadth_first_order, connected_components

        self.n = len(first)
        self.active = bool(any(first.pbc)) and first.cell.rank == 3
        nums = first.get_atomic_numbers()
        i, j, _ = bond_pairs(first, factor)  # i < j
        g = coo_matrix((np.ones(len(i)), (i, j)), shape=(self.n, self.n)).tocsr()
        ncomp, labels = connected_components(g, directed=False)
        self.labels = labels
        self.n_molecules = int(ncomp)
        self.sizes = np.bincount(labels, minlength=ncomp) if self.n else np.array([], int)
        self.bonds = (i, j)
        self.limit = factor * (covalent_radii[nums[i]] + covalent_radii[nums[j]]) if len(i) else np.zeros(0)
        order, parent, depth = [], np.full(self.n, -1), np.zeros(self.n, int)
        for c in range(ncomp):
            root = int(np.where(labels == c)[0][0])
            o, pred = breadth_first_order(g, root, directed=False, return_predecessors=True)
            for a in o[1:]:
                parent[a] = pred[a]; depth[a] = depth[pred[a]] + 1
            order.extend(o.tolist())
        self.parent, self.depth = parent, depth
        self.levels = [np.where(depth == k)[0] for k in range(1, int(depth.max()) + 1)] if self.n else []
        self.network = np.zeros(ncomp, bool)
        self.bonds_in_cell_pairs = None

    def apply(self, fr: Atoms) -> Atoms:
        if not self.active or len(fr) != self.n:
            return fr
        cell = np.asarray(fr.cell, dtype=float)
        orig = fr.get_positions()
        p = orig.copy()
        for lv in self.levels:
            par = self.parent[lv]
            p[lv] = p[par] + _mic_step(orig[lv] - p[par], cell)
        a, b = self.bonds
        if len(a):
            bad = np.linalg.norm(p[a] - p[b], axis=1) > self.limit + 1e-6
            if bad.any():
                self.network[np.unique(self.labels[a[bad]])] = True
        ncomp = self.n_molecules
        center = np.zeros((ncomp, 3))
        np.add.at(center, self.labels, p)
        center /= np.maximum(self.sizes, 1)[:, None]
        shift = -np.floor(np.linalg.solve(cell.T, center.T).T) @ cell
        p += shift[self.labels]
        keep = self.network[self.labels]
        p[keep] = orig[keep]
        out = fr.copy()
        out.set_positions(p)
        return out


def dos(eigs_ev: np.ndarray, weights: np.ndarray, sigma: float = 0.1, npts: int = 800, emin: float | None = None, emax: float | None = None):
    if len(eigs_ev) > 200 and np.all(np.diff(eigs_ev) > 0):
        return eigs_ev, weights
    lo = emin if emin is not None else float(eigs_ev.min()) - 2
    hi = emax if emax is not None else float(eigs_ev.max()) + 2
    x = np.linspace(lo, hi, npts)
    y = np.zeros_like(x)
    for e, w in zip(eigs_ev, weights):
        y += w * np.exp(-0.5 * ((x - e) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
    return x, y


def bonds(atoms: Atoms, factor: float = 1.2) -> list[tuple[str, float]]:
    sym = atoms.get_chemical_symbols()
    i, j, d = bond_pairs(atoms, factor)
    order = np.lexsort((j, i))
    return [(f"{sym[a]}{a + 1}-{sym[b]}{b + 1}", float(x)) for a, b, x in zip(i[order], j[order], d[order])]


def spectrum(freqs_cm1: list[float], intensities: list[float] | None, sigma: float = 20.0, npts: int = 2000):
    f = np.array(freqs_cm1)
    w = np.array(intensities) if intensities else np.ones_like(f)
    keep = f > 1.0
    f, w = f[keep], w[keep]
    x = np.linspace(0, max(4000.0, float(f.max()) + 200) if len(f) else 4000.0, npts)
    y = np.zeros_like(x)
    for fi, wi in zip(f, w):
        y += wi * np.exp(-0.5 * ((x - fi) / sigma) ** 2)
    return x, y
