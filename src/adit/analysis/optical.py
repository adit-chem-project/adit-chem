"""Optical quantities from the dielectric function written by Quantum ESPRESSO epsilon.x."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from adit.lang import L

HBAR_EV_S = 6.582119569e-16
C_CM_S = 2.99792458e10


@dataclass
class Optical:
    energy_ev: np.ndarray     # (n,)
    eps1: np.ndarray          # (n, 3) x y z
    eps2: np.ndarray          # (n, 3)

    def _abs_eps(self) -> np.ndarray:
        return np.hypot(self.eps1, self.eps2)

    def refractive_index(self) -> np.ndarray:
        return np.sqrt(np.maximum((self._abs_eps() + self.eps1) / 2.0, 0.0))

    def extinction(self) -> np.ndarray:
        return np.sqrt(np.maximum((self._abs_eps() - self.eps1) / 2.0, 0.0))

    def absorption_cm(self) -> np.ndarray:
        omega = self.energy_ev[:, None] / HBAR_EV_S
        return 2.0 * omega * self.extinction() / C_CM_S

    def reflectivity(self) -> np.ndarray:
        n, k = self.refractive_index(), self.extinction()
        return ((n - 1.0) ** 2 + k ** 2) / ((n + 1.0) ** 2 + k ** 2)

    def eels(self) -> np.ndarray:
        denom = self.eps1 ** 2 + self.eps2 ** 2
        return np.divide(self.eps2, denom, out=np.zeros_like(denom), where=denom > 0)

    def static_dielectric(self) -> np.ndarray:
        return self.eps1[int(np.argmin(self.energy_ev))]


def _read_columns(path: Path) -> tuple[np.ndarray, np.ndarray]:
    rows = [line for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    if not rows:
        raise ValueError(L(f"{Path(path).name} に数値の行がありません", f"{Path(path).name} has no numeric rows"))
    data = np.array([[float(x) for x in re.split(r"\s+", r.strip())] for r in rows])
    if data.shape[1] < 4:
        raise ValueError(L(f"{Path(path).name} の列が 4 つ (エネルギーと x y z) ありません",
                           f"{Path(path).name} does not have 4 columns (energy and x y z)"))
    return data[:, 0], data[:, 1:4]


def read_epsilon(run_dir: Path, prefix: str | None = None) -> Optical:
    """Read the epsr_<prefix>.dat and epsi_<prefix>.dat files written by epsilon.x."""
    run_dir = Path(run_dir)
    real = sorted(run_dir.glob(f"epsr_{prefix}.dat" if prefix else "epsr_*.dat"))
    imag = sorted(run_dir.glob(f"epsi_{prefix}.dat" if prefix else "epsi_*.dat"))
    if not real or not imag:
        raise FileNotFoundError(L(f"{run_dir} に epsr_*.dat と epsi_*.dat がありません "
                                  "(epsilon.x を実行してください)",
                                  f"{run_dir} has no epsr_*.dat and epsi_*.dat (run epsilon.x)"))
    e1_grid, eps1 = _read_columns(real[0])
    e2_grid, eps2 = _read_columns(imag[0])
    if e1_grid.shape != e2_grid.shape or not np.allclose(e1_grid, e2_grid):
        raise ValueError(L("epsr と epsi のエネルギーの刻みが違います", "the energy grids of epsr and epsi differ"))
    return Optical(e1_grid, eps1, eps2)


def summary_lines(op: Optical) -> list[str]:
    eps0 = op.static_dielectric()
    alpha = op.absorption_cm()
    lines = [L(f"光学 (epsilon.x、{op.energy_ev.size} 点、"
               f"{op.energy_ev.min():.2f}〜{op.energy_ev.max():.2f} eV)",
               f"optical properties (epsilon.x, {op.energy_ev.size} points, "
               f"{op.energy_ev.min():.2f}-{op.energy_ev.max():.2f} eV)"),
             L("  静的誘電率 eps1(E→0): x {:.3f} / y {:.3f} / z {:.3f}",
               "  static dielectric constant eps1(E→0): x {:.3f} / y {:.3f} / z {:.3f}").format(*eps0)]
    for axis, i in (("x", 0), ("y", 1), ("z", 2)):
        peak = int(np.argmax(op.eps2[:, i]))
        lines.append(L(f"  {axis}: eps2 の最大 {op.eps2[peak, i]:.3f} @ {op.energy_ev[peak]:.2f} eV、"
                       f"吸収係数の最大 {alpha[:, i].max():.3e} 1/cm",
                       f"  {axis}: maximum eps2 {op.eps2[peak, i]:.3f} @ {op.energy_ev[peak]:.2f} eV, "
                       f"maximum absorption coefficient {alpha[:, i].max():.3e} 1/cm"))
    lines.append(L("  k 点の数と経験的な広がりで値が変わります。収束は利用者が確かめてください。",
                   "  The values depend on the number of k-points and the empirical broadening. Check convergence yourself."))
    return lines
