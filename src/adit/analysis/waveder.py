"""Read VASP WAVEDER files (band-to-band dipole matrix elements)."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from adit.lang import L

HBAR_EV_S = 6.582119569e-16
RYDBERG_EV = 13.605693122994
BOHR_ANG = 0.529177210903
EDEPS = 4 * np.pi * 2 * RYDBERG_EV * BOHR_ANG


@dataclass
class Waveder:
    cder: np.ndarray

    @property
    def n_bands(self) -> int:
        return self.cder.shape[0]

    @property
    def n_kpoints(self) -> int:
        return self.cder.shape[2]

    def matrix_element(self, band_i: int, band_j: int, kpoint: int, spin: int = 0) -> np.ndarray:
        return self.cder[band_i, band_j, kpoint, spin, :]


def _int32(stream) -> int:
    raw = stream.read(4)
    if len(raw) < 4:
        raise ValueError(L("WAVEDER の途中でファイルが終わりました", "the WAVEDER file ended early"))
    return int(np.frombuffer(raw, dtype=np.int32)[0])


def _read_record(stream, dtype):
    data = b""
    while True:
        n = _int32(stream)
        data += stream.read(abs(n))
        if abs(n) != abs(_int32(stream)):
            raise ValueError(L("WAVEDER のレコードの長さが前後で合いません (別の書式かもしれません)",
                               "the WAVEDER record lengths do not match (it may be another format)"))
        if n > 0:
            break
    return np.frombuffer(data, dtype=dtype)


def read_waveder(path: Path, dtype: str = "complex64") -> Waveder:
    """Read a WAVEDER file. For gamma-only VASP builds pass dtype as "float64" or "float32"."""
    p = Path(path)
    opener = gzip.open if p.suffix == ".gz" else open
    with opener(p, "rb") as stream:
        header = _read_record(stream, np.int32)
        if header.size < 4:
            raise ValueError(L("WAVEDER の見出し (nbands, nelect, nk, ispin) を読めません",
                               "cannot read the WAVEDER header (nbands, nelect, nk, ispin)"))
        nbands, nelect, nk, ispin = (int(x) for x in header[:4])
        _read_record(stream, np.float64)
        _read_record(stream, np.float64)
        flat = _read_record(stream, np.dtype(dtype))
    expected = 3 * ispin * nk * nelect * nbands
    if flat.size != expected:
        raise ValueError(L(f"WAVEDER の行列要素が {flat.size} 個で、見出しから求めた {expected} 個と違います "
                           "(ガンマ点版なら dtype に float64 か float32 を渡してください)",
                           f"WAVEDER has {flat.size} matrix elements, not the {expected} expected from the header "
                           "(for a gamma-only build pass dtype as float64 or float32)"))
    return Waveder(flat.reshape((3, ispin, nk, nelect, nbands)).T)


def dielectric_imag(wd: Waveder, energies_ev, occupations, weights, volume_ang3: float,
                    grid_ev, broadening_ev: float, spin: int = 0) -> np.ndarray:
    e = np.asarray(energies_ev, dtype=float)          # (nk, nbands)
    f = np.asarray(occupations, dtype=float)          # (nk, nbands)
    w = np.asarray(weights, dtype=float)              # (nk,)
    grid = np.asarray(grid_ev, dtype=float)
    if e.shape != f.shape:
        raise ValueError(L(f"固有値 {e.shape} と占有数 {f.shape} の形が違います",
                           f"the eigenvalues {e.shape} and occupations {f.shape} differ in shape"))
    if w.shape[0] != e.shape[0]:
        raise ValueError(L(f"k 点の重み {w.shape[0]} 個が固有値の k 点 {e.shape[0]} 個と合いません",
                           f"{w.shape[0]} k-point weights do not match {e.shape[0]} k-points in the eigenvalues"))
    if broadening_ev <= 0:
        raise ValueError(L("広がりの幅は正にしてください", "the broadening width must be positive"))
    nk, nb = e.shape
    if wd.cder.shape[2] < nk or wd.cder.shape[0] < nb:
        raise ValueError(L("WAVEDER の大きさが、渡された固有値の大きさに足りません", "WAVEDER is smaller than the eigenvalues passed in"))
    rspin = 3 - wd.cder.shape[3]
    w = w / w.sum()
    out = np.zeros((grid.size, 3))
    norm = 1.0 / (np.sqrt(2 * np.pi) * broadening_ev)
    for k in range(nk):
        de = e[k][None, :] - e[k][:, None]             # (v, c)
        df = f[k][:, None] - f[k][None, :]
        mask = (de > 0) & (df > 1e-8)
        if not mask.any():
            continue
        v_idx, c_idx = np.nonzero(mask)
        me = wd.cder[v_idx, c_idx, k, spin, :]         # (n, 3)
        strength = (df[v_idx, c_idx] * rspin)[:, None] * np.abs(me) ** 2
        gauss = norm * np.exp(-0.5 * ((grid[:, None] - de[v_idx, c_idx][None, :]) / broadening_ev) ** 2)
        out += w[k] * (gauss @ strength)
    return (EDEPS * np.pi / volume_ang3) * out
