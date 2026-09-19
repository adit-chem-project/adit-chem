"""Count hydrogen bonds. The geometric criteria are supplied by the caller; there are no defaults."""

from __future__ import annotations

import numpy as np
from ase import Atoms

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


class HydrogenBondError(AditValueError):
    pass


def criteria_note() -> str:
    parts = [f"{v['source']}" for v in KNOWN_CRITERIA.values()]
    return L("よく使われる値: " + " / ".join(parts) + "。ADIT は既定値を持ちません (どれを採るかは利用者が決めます)",
             "commonly used criteria: " + " / ".join(parts) + ". ADIT has no default; the choice is yours")


def count_frame(atoms: Atoms, donor_acceptor_A: float, angle_deg: float, donors=DEFAULT_DONORS,
                hydrogen: str = "H", max_dh_A: float = 1.3) -> dict:
    if donor_acceptor_A <= 0 or not 0 < angle_deg <= 180:
        raise HydrogenBondError(L("距離 [Å] と角度 [度] の両方を指定してください (角度は 0 より大きく 180 以下)",
                                  "give both a distance in Å and an angle in degrees (0 < angle <= 180)"))
    from adit.analysis.compute import _mic_step

    syms = np.array(atoms.get_chemical_symbols())
    pos = atoms.get_positions()
    cell = np.asarray(atoms.cell, dtype=float)
    periodic = bool(np.any(atoms.pbc)) and abs(np.linalg.det(cell)) > 0

    def rel(a, b):
        d = pos[b] - pos[a]
        return _mic_step(d.reshape(1, 3), cell).reshape(3) if periodic else d

    heavy = np.flatnonzero(np.isin(syms, list(donors)))
    hydrogens = np.flatnonzero(syms == hydrogen)
    if not len(heavy) or not len(hydrogens):
        return {"count": 0, "pairs": [], "donors": list(donors)}
    bonded: dict[int, int] = {}
    for h in hydrogens:
        best, best_d = None, None
        for d_atom in heavy:
            dist = float(np.linalg.norm(rel(h, d_atom)))
            if dist <= max_dh_A and (best_d is None or dist < best_d):
                best, best_d = int(d_atom), dist
        if best is not None:
            bonded[int(h)] = best
    pairs = []
    for h, donor in bonded.items():
        for acceptor in heavy:
            acceptor = int(acceptor)
            if acceptor == donor:
                continue
            da = float(np.linalg.norm(rel(donor, acceptor)))
            if da > donor_acceptor_A:
                continue
            u, v = rel(h, donor), rel(h, acceptor)
            cos = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
            angle = float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))   # D–H···A
            if angle >= angle_deg:
                pairs.append({"donor": donor + 1, "hydrogen": h + 1, "acceptor": acceptor + 1,
                              "distance_A": da, "angle_deg": angle})
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
                      "ADIT has no default criteria and makes no judgement about the count.")}
