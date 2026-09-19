
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from adit.errors import AditValueError
from adit.lang import L


class ElasticModuliError(AditValueError):
    pass


@dataclass
class Moduli:
    k_voigt: float
    k_reuss: float
    g_voigt: float
    g_reuss: float
    stable: bool | None

    @property
    def k_hill(self) -> float:
        return 0.5 * (self.k_voigt + self.k_reuss)

    @property
    def g_hill(self) -> float:
        return 0.5 * (self.g_voigt + self.g_reuss)

    @property
    def young(self) -> float:
        k, g = self.k_hill, self.g_hill
        return 9.0 * k * g / (3.0 * k + g) if (3.0 * k + g) else float("nan")

    @property
    def poisson(self) -> float:
        k, g = self.k_hill, self.g_hill
        return (3.0 * k - 2.0 * g) / (2.0 * (3.0 * k + g)) if (3.0 * k + g) else float("nan")

    def as_dict(self) -> dict:
        return {"bulk_modulus_gpa": {"voigt": self.k_voigt, "reuss": self.k_reuss, "hill": self.k_hill},
                "shear_modulus_gpa": {"voigt": self.g_voigt, "reuss": self.g_reuss, "hill": self.g_hill},
                "young_modulus_gpa": self.young, "poisson_ratio": self.poisson,
                "pugh_ratio": (self.k_hill / self.g_hill) if self.g_hill else None,
                "born_stability": self.stable,
                "note": L("Voigt (歪み一様) と Reuss (応力一様) が上限と下限、Hill はその平均です "
                          "(Hill 1952)。ヤング率とポアソン比は Hill の値から出しています。"
                          "born_stability は C の固有値がすべて正かどうかという機械的な事実で、"
                          "材料の良し悪しの判定ではありません。",
                          "Voigt (uniform strain) and Reuss (uniform stress) bound the result and Hill is their average "
                          "(Hill 1952). Young's modulus and Poisson's ratio come from the Hill values. "
                          "born_stability only states whether all eigenvalues of C are positive; it is not a judgment "
                          "about the material.")}


def from_cij(cij) -> Moduli:
    c = np.asarray(cij, dtype=float)
    if c.shape != (6, 6):
        raise ElasticModuliError(L(f"弾性定数は 6×6 の行列です (いまの形 {c.shape})",
                                   f"the elastic constants must be a 6x6 matrix (got {c.shape})"))
    if not np.all(np.isfinite(c)):
        raise ElasticModuliError(L("弾性定数に数でない値があります", "the elastic constants contain non-finite values"))
    sym = 0.5 * (c + c.T)
    try:
        s = np.linalg.inv(sym)
    except np.linalg.LinAlgError as ex:
        raise ElasticModuliError(L("弾性定数の行列が逆行列を持ちません (Reuss の平均を出せません)",
                                   "the elastic matrix is singular, so the Reuss average cannot be computed")) from ex
    k_v = ((sym[0, 0] + sym[1, 1] + sym[2, 2]) + 2 * (sym[0, 1] + sym[1, 2] + sym[2, 0])) / 9.0
    g_v = ((sym[0, 0] + sym[1, 1] + sym[2, 2]) - (sym[0, 1] + sym[1, 2] + sym[2, 0])
           + 3 * (sym[3, 3] + sym[4, 4] + sym[5, 5])) / 15.0
    k_r_inv = (s[0, 0] + s[1, 1] + s[2, 2]) + 2 * (s[0, 1] + s[1, 2] + s[2, 0])
    g_r_inv = 4 * (s[0, 0] + s[1, 1] + s[2, 2]) - 4 * (s[0, 1] + s[1, 2] + s[2, 0]) + 3 * (s[3, 3] + s[4, 4] + s[5, 5])
    k_r = 1.0 / k_r_inv if k_r_inv else float("nan")
    g_r = 15.0 / g_r_inv if g_r_inv else float("nan")
    stable = bool(np.all(np.linalg.eigvalsh(sym) > 0))
    return Moduli(k_voigt=float(k_v), k_reuss=float(k_r), g_voigt=float(g_v), g_reuss=float(g_r), stable=stable)


def summary_lines(m: Moduli) -> list[str]:
    out = [L(f"弾性率 (Voigt–Reuss–Hill): 体積弾性率 K = {m.k_hill:.1f} GPa "
             f"(Voigt {m.k_voigt:.1f}、Reuss {m.k_reuss:.1f})、剛性率 G = {m.g_hill:.1f} GPa "
             f"(Voigt {m.g_voigt:.1f}、Reuss {m.g_reuss:.1f})",
             f"elastic moduli (Voigt-Reuss-Hill): bulk K = {m.k_hill:.1f} GPa "
             f"(Voigt {m.k_voigt:.1f}, Reuss {m.k_reuss:.1f}), shear G = {m.g_hill:.1f} GPa "
             f"(Voigt {m.g_voigt:.1f}, Reuss {m.g_reuss:.1f})"),
          L(f"ヤング率 E = {m.young:.1f} GPa、ポアソン比 ν = {m.poisson:.3f}、K/G = {m.k_hill / m.g_hill:.2f} "
            "(Hill の値から。良し悪しは判定していません)",
            f"Young's modulus E = {m.young:.1f} GPa, Poisson's ratio nu = {m.poisson:.3f}, K/G = {m.k_hill / m.g_hill:.2f} "
            "(from the Hill values; no judgment is made)")]
    if m.stable is False:
        out.append(L("弾性定数の行列に 0 以下の固有値があります (Born の安定条件を満たしません)。"
                     "計算の設定か構造を確かめてください",
                     "the elastic matrix has a non-positive eigenvalue (it does not satisfy the Born stability criteria); "
                     "check the settings or the structure"))
    return out
