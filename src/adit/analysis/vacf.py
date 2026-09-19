"""Velocity autocorrelation, the vibrational spectrum derived from it, and the Green-Kubo diffusion coefficient."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from adit.errors import AditValueError
from adit.lang import L

CM1_PER_FS1 = 1.0 / 2.99792458e-5
# 1 Å²/fs = 0.1 cm²/s
A2_FS_TO_CM2_S = 0.1


class VacfError(AditValueError):
    pass


@dataclass
class VacfResult:
    times_fs: np.ndarray
    vacf: np.ndarray
    running_d_cm2_s: np.ndarray
    freq_cm1: np.ndarray
    spectrum: np.ndarray
    source: str
    dt_fs: float
    max_lag_fs: float

    @property
    def d_cm2_s(self) -> float:
        return float(self.running_d_cm2_s[-1]) if len(self.running_d_cm2_s) else float("nan")

    def as_dict(self) -> dict:
        return {"dt_fs": self.dt_fs, "max_lag_fs": self.max_lag_fs, "source": self.source,
                "d_green_kubo_cm2_s": self.d_cm2_s,
                "times_fs": self.times_fs.tolist(), "vacf": self.vacf.tolist(),
                "running_d_cm2_s": self.running_d_cm2_s.tolist(),
                "freq_cm1": self.freq_cm1.tolist(), "spectrum": self.spectrum.tolist(),
                "note": L("VACF は時間原点を全部使った平均です。D は Green-Kubo の式 (1/3)∫⟨v(0)·v(t)⟩dt で、"
                          "MSD の傾きから出した D と一致するはずの量です。スペクトルは VACF のフーリエ変換 "
                          "(振動の状態密度)。収束もピークの帰属も判定していません。",
                          "The VACF averages over all time origins. D is the Green-Kubo integral (1/3) int <v(0).v(t)> dt, "
                          "which should agree with the slope of the MSD. The spectrum is the Fourier transform of the VACF "
                          "(vibrational density of states). Neither convergence nor peak assignment is judged.")}


def velocities_from_positions(pos: np.ndarray, dt_fs: float) -> np.ndarray:
    pos = np.asarray(pos, dtype=float)
    if pos.shape[0] < 3:
        raise VacfError(L("速度を差分で作るにはフレームが 3 つ以上要ります",
                          "at least three frames are needed to make velocities by finite differences"))
    v = np.empty_like(pos)
    v[1:-1] = (pos[2:] - pos[:-2]) / (2.0 * dt_fs)
    v[0] = (pos[1] - pos[0]) / dt_fs
    v[-1] = (pos[-1] - pos[-2]) / dt_fs
    return v


def _autocorrelation(v: np.ndarray, chunk_bytes: int = 32 * 2 ** 20) -> np.ndarray:
    frames, natoms, _ = v.shape
    size = 1 << (2 * frames - 1).bit_length()
    total = np.zeros(frames)
    per_chunk = max(1, int(chunk_bytes // max(1, size * 3 * 16)))
    for start in range(0, natoms, per_chunk):
        x = v[:, start:start + per_chunk, :]
        f = np.fft.rfft(x, n=size, axis=0)
        corr = np.fft.irfft(f * np.conjugate(f), n=size, axis=0)[:frames]
        total += np.sum(corr, axis=(1, 2))
    counts = (frames - np.arange(frames)).astype(float)
    return total / counts / natoms


def vacf(velocities: np.ndarray, dt_fs: float, *, max_lag_fraction: float = 0.5,
         source: str = "velocities", window: str = "hann") -> VacfResult:
    """VACF, Green-Kubo D and power spectrum. velocities has shape (frames, atoms, 3) in Angstrom/fs; window is "hann" or "none"."""
    v = np.asarray(velocities, dtype=float)
    if v.ndim != 3 or v.shape[0] < 4:
        raise VacfError(L("(フレーム, 原子, 3) の速度が 4 フレーム以上必要です",
                          "velocities of shape (frames, atoms, 3) with at least four frames are required"))
    if dt_fs <= 0:
        raise VacfError(L("フレームの間隔 [fs] が要ります", "the time between frames (fs) is required"))
    raw = _autocorrelation(v)                       # [Å²/fs²]
    keep = max(4, int(len(raw) * max_lag_fraction))
    raw = raw[:keep]
    times = np.arange(keep, dtype=float) * dt_fs
    running = np.concatenate([[0.0], np.cumsum((raw[1:] + raw[:-1]) / 2.0) * dt_fs]) / 3.0 * A2_FS_TO_CM2_S
    norm = raw / raw[0] if raw[0] != 0 else raw
    w = np.hanning(2 * keep)[keep:] if window == "hann" else np.ones(keep)
    size = 1 << (4 * keep - 1).bit_length()
    spec = np.abs(np.fft.rfft(raw * w, n=size))
    freq = np.fft.rfftfreq(size, d=dt_fs) * CM1_PER_FS1
    return VacfResult(times_fs=times, vacf=norm, running_d_cm2_s=running, freq_cm1=freq, spectrum=spec,
                      source=source, dt_fs=float(dt_fs), max_lag_fs=float(times[-1]))


def notes(result: VacfResult) -> list[str]:
    nyquist = CM1_PER_FS1 / (2.0 * result.dt_fs)
    top = float(result.freq_cm1[int(np.argmax(result.spectrum[1:])) + 1]) if len(result.spectrum) > 1 else float("nan")
    how = (L("軌跡に入っていた速度", "velocities stored in the trajectory") if result.source == "velocities"
           else L("座標の差分で作った速度", "velocities made by finite differences of the positions"))
    return [
        L(f"速度自己相関 (VACF): {how}、遅れ時間 {result.max_lag_fs:g} fs まで。"
          f"Green-Kubo の D = {result.d_cm2_s:.4g} cm²/s (MSD の傾きと突き合わせてください)",
          f"velocity autocorrelation (VACF): {how}, lags up to {result.max_lag_fs:g} fs. "
          f"Green-Kubo D = {result.d_cm2_s:.4g} cm^2/s (compare it with the slope of the MSD)"),
        L(f"振動スペクトル (VACF のフーリエ変換): いちばん強い振動数は {top:.0f} cm⁻¹。"
          f"書き出し間隔 {result.dt_fs:g} fs で見えるのは {nyquist:.0f} cm⁻¹ まで (それより速い振動は出ません)",
          f"vibrational spectrum (Fourier transform of the VACF): the strongest frequency is {top:.0f} cm^-1. "
          f"With a frame spacing of {result.dt_fs:g} fs, frequencies above {nyquist:.0f} cm^-1 cannot appear"),
    ]
