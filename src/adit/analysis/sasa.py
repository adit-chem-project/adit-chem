"""Solvent-accessible surface area by the Shrake-Rupley method."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from adit.lang import L

RADII_SETS: dict[str, tuple[str, dict[str, float]]] = {
    "bondi": ("Bondi, J. Phys. Chem. 68, 441 (1964)",
              {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47, "P": 1.80, "S": 1.80,
               "Cl": 1.75, "Br": 1.85, "I": 1.98, "Na": 2.27, "K": 2.75, "Mg": 1.73, "Ca": 2.31}),
    "ase": ("ASE ase.data.vdw_radii", {}),
}


def radii_description(name: str) -> str:
    if name == "bondi":
        return L("Bondi, J. Phys. Chem. 68, 441 (1964) のファンデルワールス半径",
                 "van der Waals radii from Bondi, J. Phys. Chem. 68, 441 (1964)")
    if name == "ase":
        return L("ASE の ase.data.vdw_radii (Bondi と Alvarez の値)",
                 "ASE ase.data.vdw_radii (values from Bondi and Alvarez)")
    return name


@dataclass
class SasaResult:
    per_atom_ang2: np.ndarray
    radii_ang: np.ndarray
    probe_ang: float
    n_points: int
    radii_source: str

    @property
    def total_ang2(self) -> float:
        return float(self.per_atom_ang2.sum())


def _sphere_points(n: int) -> np.ndarray:
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.column_stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)])


def radii_for(symbols, radii_set: str = "bondi") -> np.ndarray:
    if radii_set == "ase":
        from ase.data import atomic_numbers, vdw_radii

        out = []
        missing = []
        for s in symbols:
            r = vdw_radii[atomic_numbers[s]]
            (out if np.isfinite(r) else missing).append(r if np.isfinite(r) else s)
        if missing:
            raise ValueError(L(f"ASE の表に半径がありません: {', '.join(sorted(set(missing)))}",
                               f"no radius in the ASE table for: {', '.join(sorted(set(missing)))}"))
        return np.array(out, dtype=float)
    if radii_set not in RADII_SETS:
        raise ValueError(L(f"半径の表は {', '.join(RADII_SETS)} のどれかにしてください",
                           f"the radii table must be one of {', '.join(RADII_SETS)}"))
    table = RADII_SETS[radii_set][1]
    missing = sorted({s for s in symbols if s not in table})
    if missing:
        raise ValueError(L(f"{radii_set} の表にない元素です: {', '.join(missing)}。"
                           "半径を自分で渡すか、別の表を選んでください",
                           f"elements not in the {radii_set} table: {', '.join(missing)}. "
                           "Pass the radii yourself or choose another table"))
    return np.array([table[s] for s in symbols], dtype=float)


def sasa(positions, radii_ang, probe_ang: float = 1.4, n_points: int = 960,
         radii_source: str = "") -> SasaResult:
    """Shrake-Rupley SASA. Atom radii and the probe radius are supplied by the caller."""
    pos = np.asarray(positions, dtype=float)
    r = np.asarray(radii_ang, dtype=float)
    if pos.ndim != 2 or pos.shape[1] != 3:
        raise ValueError(L("positions は (原子数, 3) にしてください", "positions must have shape (atoms, 3)"))
    if r.shape != (pos.shape[0],):
        raise ValueError(L(f"半径 {r.size} 個が原子数 {pos.shape[0]} 個と合いません", f"{r.size} radii do not match {pos.shape[0]} atoms"))
    if probe_ang < 0:
        raise ValueError(L("プローブ半径は 0 以上にしてください", "the probe radius must be 0 or more"))
    if n_points < 20:
        raise ValueError(L("点の数は 20 以上にしてください", "the number of points must be 20 or more"))
    sphere = _sphere_points(n_points)
    big = r + probe_ang
    area = np.zeros(pos.shape[0])
    for i in range(pos.shape[0]):
        test = pos[i] + big[i] * sphere
        d = np.linalg.norm(test[:, None, :] - pos[None, :, :], axis=2)
        d[:, i] = np.inf
        exposed = np.all(d >= big[None, :], axis=1)
        area[i] = 4 * np.pi * big[i] ** 2 * exposed.mean()
    return SasaResult(area, r, probe_ang, n_points, radii_source)


def summary_lines(res: SasaResult) -> list[str]:
    return [L(f"溶媒接触表面積 (Shrake–Rupley、点 {res.n_points} 個、プローブ {res.probe_ang:g} Å)",
              f"solvent-accessible surface area (Shrake-Rupley, {res.n_points} points, probe {res.probe_ang:g} Å)"),
            L(f"  合計 {res.total_ang2:.2f} Å²、原子あたり 平均 {res.per_atom_ang2.mean():.3f} Å²、"
              f"表に出ていない原子 {int((res.per_atom_ang2 == 0).sum())} 個",
              f"  total {res.total_ang2:.2f} Å², mean per atom {res.per_atom_ang2.mean():.3f} Å², "
              f"{int((res.per_atom_ang2 == 0).sum())} buried atoms"),
            L(f"  半径の出典: {res.radii_source or '利用者が渡した値'}",
              f"  source of the radii: {res.radii_source or 'values passed by the user'}"),
            L("  半径とプローブ半径で値が変わります。比べるときは同じ値を使ってください。",
              "  The result depends on the radii and the probe radius. Use the same values when comparing.")]
