
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from adit.analysis import msd_worker as W
from adit.errors import AditValueError
from adit.lang import L

mic = W.mic
cell_half_width = W.cell_half_width
fit_gaussian_d = W.fit_gaussian_d
default_taus = W.default_taus
NEAR_CAP, CAP_FRACTION, CHUNK_BYTES = W.NEAR_CAP, W.CAP_FRACTION, W.CHUNK_BYTES


class VanHoveError(AditValueError):
    pass


@dataclass
class VanHoveResult:
    times: np.ndarray
    d_fit: np.ndarray
    d_direct: np.ndarray
    alpha2: np.ndarray
    rms: np.ndarray
    truncated_from: float | None
    half_width_A: float | None
    displacement: str
    near_cap_fraction: float = 0.0
    histograms: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {"times_fs": self.times.tolist(),
                "d_fit_cm2_s": [float(v) * 1e-1 for v in self.d_fit],
                "d_direct_cm2_s": [float(v) * 1e-1 for v in self.d_direct],
                "alpha2": self.alpha2.tolist(), "rms_A": self.rms.tolist(),
                "mean_d_fit_cm2_s": float(np.nanmean(self.d_fit) * 1e-1),
                "mean_d_direct_cm2_s": float(np.mean(self.d_direct) * 1e-1),
                "max_alpha2": float(np.nanmax(self.alpha2)) if len(self.alpha2) else None,
                "displacement": self.displacement, "half_width_A": self.half_width_A,
                "truncated_from_fs": self.truncated_from, "near_cap_fraction": self.near_cap_fraction}


def van_hove_self(pos: np.ndarray, cell, taus, dt_fs: float, *, displacement: str = "mic",
                  keep_histograms=(), chunk_bytes: int = W.CHUNK_BYTES) -> VanHoveResult:
    pos = np.asarray(pos, dtype=float)
    if pos.ndim != 3 or pos.shape[0] < 2:
        raise VanHoveError(L("フレームが 2 つ以上ある (フレーム, 原子, 3) の配列が要ります",
                             "a (frames, atoms, 3) array with at least two frames is required"))
    if displacement not in ("mic", "unwrapped"):
        raise VanHoveError(L("displacement は mic か unwrapped です", "displacement must be 'mic' or 'unwrapped'"))
    if displacement == "mic" and cell is None:
        raise VanHoveError(L("最小像を使うにはセルが要ります", "a cell is required for minimum-image displacements"))
    usable = [int(t) for t in taus if 0 < int(t) < pos.shape[0]]
    if not usable:
        raise VanHoveError(L("使える遅れ時間がありません (フレーム数より小さい値を指定してください)",
                             "no usable lag times (give values smaller than the number of frames)"))
    raw = W.van_hove_self(pos, cell, usable, dt_fs, displacement=displacement, chunk_bytes=chunk_bytes)
    hists = {}
    for tau in keep_histograms:
        if 0 < int(tau) < pos.shape[0]:
            _, _, r, _ = W.moments(pos, cell, int(tau), displacement == "mic" and cell is not None, chunk_bytes)
            counts, edges = np.histogram(r, bins="fd")
            hists[int(tau) * dt_fs] = ((edges[:-1] + edges[1:]) / 2.0, counts)
    return VanHoveResult(times=np.array(raw["times_fs"]), d_fit=np.array(raw["d_fit_A2_fs"]),
                         d_direct=np.array(raw["d_direct_A2_fs"]), alpha2=np.array(raw["alpha2"]),
                         rms=np.array(raw["rms_A"]), truncated_from=raw["truncated_from_fs"],
                         half_width_A=raw["half_width_A"], displacement=raw["displacement"],
                         near_cap_fraction=raw["near_cap_fraction"], histograms=hists)


def notes(result: VanHoveResult) -> list[str]:
    out = [L(f"変位の分布から出した D: 当てはめ {np.nanmean(result.d_fit) * 1e-1:.4g} cm²/s、"
             f"⟨r²⟩/(6τ) {np.mean(result.d_direct) * 1e-1:.4g} cm²/s "
             f"(遅れ時間 {result.times[0]:g}〜{result.times[-1]:g} fs の {len(result.times)} 点の平均)",
             f"D from the displacement distribution: fit {np.nanmean(result.d_fit) * 1e-1:.4g} cm^2/s, "
             f"<r^2>/(6 tau) {np.mean(result.d_direct) * 1e-1:.4g} cm^2/s "
             f"(mean over {len(result.times)} lags from {result.times[0]:g} to {result.times[-1]:g} fs)")]
    top = float(np.nanmax(result.alpha2)) if len(result.alpha2) else float("nan")
    out.append(L(f"非ガウス因子 α₂ の最大 {top:.3g} (0 ならガウス = フィックの拡散。正なら跳びの多い動き)。"
                 "どこからを「非ガウス」と呼ぶかは判定していません",
                 f"largest non-Gaussian parameter alpha2 {top:.3g} (0 means Gaussian/Fickian; positive means jump-like motion). "
                 "No threshold for 'non-Gaussian' is applied"))
    if result.displacement == "mic" and result.truncated_from is not None:
        out.append(L(f"注意: 遅れ時間 {result.truncated_from:g} fs から先は、変位の {result.near_cap_fraction:.1%} が "
                     f"最小像で測れる上限 {result.half_width_A:.2f} Å (セルの最小の幅の半分) の 9 割を超えています。"
                     "最小像の変位は頭打ちになるので、この範囲の D は低めに出ます (--vanhove-displacement unwrapped で巻き戻した変位を使えます)",
                     f"note: from lag {result.truncated_from:g} fs, {result.near_cap_fraction:.1%} of the displacements are beyond "
                     f"90 % of {result.half_width_A:.2f} Å, half the smallest cell width, which is the largest displacement a "
                     "minimum-image convention can represent; D is biased low there (use --vanhove-displacement unwrapped)"))
    return out
