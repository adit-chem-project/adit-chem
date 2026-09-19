
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from ase import Atoms

from adit.errors import AditValueError
from adit.lang import L

DEFAULT_RADIATION = "CuKa"
DEFAULT_RANGE = (5.0, 90.0)


class XrdError(AditValueError):
    pass


@dataclass
class Pattern:
    two_theta: np.ndarray
    intensity: np.ndarray
    hkl: list[str]
    d_ang: np.ndarray
    radiation: str
    wavelength_ang: float
    note: str

    def as_dict(self) -> dict:
        return {"radiation": self.radiation, "wavelength_ang": float(self.wavelength_ang),
                "n_peaks": int(len(self.two_theta)),
                "peaks": [{"two_theta_deg": float(t), "intensity": float(i), "hkl": h, "d_ang": float(d)}
                          for t, i, h, d in zip(self.two_theta, self.intensity, self.hkl, self.d_ang)],
                "note": self.note}


def has_pymatgen() -> bool:
    try:
        import pymatgen.analysis.diffraction.xrd  # noqa: F401
    except Exception:
        return False
    return True


def _hkl_text(hkls) -> str:
    out = []
    for item in hkls:
        idx = item["hkl"] if isinstance(item, dict) else item
        out.append("(" + " ".join(str(int(v)) for v in idx) + ")")
    return ", ".join(dict.fromkeys(out))


def powder_pattern(atoms: Atoms, radiation: str = DEFAULT_RADIATION,
                   two_theta_range: tuple[float, float] = DEFAULT_RANGE) -> Pattern:
    if not has_pymatgen():
        raise XrdError(L("pymatgen が入っていないので粉末回折は計算できません (pip install pymatgen)",
                         "pymatgen is not installed, so no powder pattern can be computed (pip install pymatgen)"))
    from pymatgen.analysis.diffraction.xrd import XRDCalculator
    from pymatgen.io.ase import AseAtomsAdaptor

    if not all(atoms.get_pbc()) or atoms.get_volume() <= 0:
        raise XrdError(L("粉末回折は周期系 (セルのある構造) にだけ計算します",
                         "a powder pattern is only computed for a periodic structure (one with a cell)"))
    if radiation not in XRDCalculator.AVAILABLE_RADIATION:
        raise XrdError(L(f"線源が分かりません: {radiation} (使えるのは {', '.join(XRDCalculator.AVAILABLE_RADIATION)})",
                         f"unknown radiation: {radiation} (available: {', '.join(XRDCalculator.AVAILABLE_RADIATION)})"))
    lo, hi = float(two_theta_range[0]), float(two_theta_range[1])
    if not hi > lo:
        raise XrdError(L("2θ の範囲は「下限 < 上限」で書いてください", "the 2-theta range must have lower < upper"))
    calc = XRDCalculator(wavelength=radiation)
    structure = AseAtomsAdaptor.get_structure(atoms)
    try:
        pattern = calc.get_pattern(structure, two_theta_range=(lo, hi))
    except Exception as ex:
        raise XrdError(L(f"粉末回折を計算できません: {ex}", f"cannot compute the powder pattern: {ex}")) from ex
    if not len(pattern.x):
        raise XrdError(L(f"2θ = {lo:g}〜{hi:g} 度にピークがありません", f"no peaks between 2-theta = {lo:g} and {hi:g} degrees"))
    note = L(f"pymatgen の XRDCalculator ({radiation}、λ = {calc.wavelength:.5g} Å、2θ = {lo:g}〜{hi:g} 度)。"
             "強度は最大を 100 にそろえた相対値です。測定との一致は判定していません。"
             "熱振動 (デバイ・ワラー因子) と選択配向は入れていません。指数 hkl は、与えた構造のセルでの値です (基本セルと慣用セルでは指数の付き方が変わります)。",
             f"pymatgen's XRDCalculator ({radiation}, lambda = {calc.wavelength:.5g} Å, 2-theta = {lo:g}-{hi:g} degrees). "
             "Intensities are relative with the maximum set to 100; agreement with a measurement is not assessed. "
             "Thermal (Debye-Waller) factors and preferred orientation are not included. "
             "The hkl indices refer to the cell of the structure as given (a primitive cell gives different indices from the conventional cell).")
    return Pattern(two_theta=np.asarray(pattern.x, dtype=float), intensity=np.asarray(pattern.y, dtype=float),
                   hkl=[_hkl_text(h) for h in pattern.hkls], d_ang=np.asarray(pattern.d_hkls, dtype=float),
                   radiation=radiation, wavelength_ang=float(calc.wavelength), note=note)


def read_measured(path: Path | str) -> tuple[np.ndarray, np.ndarray]:
    p = Path(path).expanduser()
    if not p.is_file():
        raise XrdError(L(f"測定データのファイルがありません: {p}", f"no such measured-pattern file: {p}"))
    xs: list[float] = []
    ys: list[float] = []
    for line in p.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = [v for v in line.replace(",", " ").replace("\t", " ").split() if v]
        if len(parts) < 2:
            continue
        try:
            x, y = float(parts[0]), float(parts[1])
        except ValueError:
            continue
        xs.append(x); ys.append(y)
    if len(xs) < 2:
        raise XrdError(L(f"{p.name} から 2 点以上の (2θ, 強度) を読めません", f"cannot read at least two (2-theta, intensity) pairs from {p.name}"))
    return np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)


def scale_to_100(values: np.ndarray) -> np.ndarray:
    top = float(np.max(values)) if len(values) else 0.0
    return values * (100.0 / top) if top > 0 else np.asarray(values, dtype=float)


def write_csv(path: Path, pattern: Pattern) -> Path:
    import csv

    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["two_theta_deg", "intensity_rel", "d_ang", "hkl"])
        writer.writerows([[f"{t:.4f}", f"{i:.4f}", f"{d:.5f}", h]
                          for t, i, d, h in zip(pattern.two_theta, pattern.intensity, pattern.d_ang, pattern.hkl)])
    return path
