"""Distance, angle and dihedral between chosen atoms, with the minimum image in periodic cells."""

from __future__ import annotations

import numpy as np

from adit.lang import L


def mic_vector(a, b, cell=None) -> np.ndarray:
    """Vector from a to b; the shortest periodic image when a cell is given."""
    d = np.asarray(b, dtype=float) - np.asarray(a, dtype=float)
    if cell is None:
        return d
    c = np.asarray(cell, dtype=float)
    if abs(np.linalg.det(c)) < 1e-12:
        return d
    frac = np.linalg.solve(c.T, d)
    frac -= np.round(frac)
    best, best_len = None, None
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for k in (-1, 0, 1):
                v = (frac + (i, j, k)) @ c
                n = float(v @ v)
                if best_len is None or n < best_len:
                    best, best_len = v, n
    return best


def distance(pos, i: int, j: int, cell=None) -> float:
    return float(np.linalg.norm(mic_vector(pos[i], pos[j], cell)))


def angle(pos, i: int, j: int, k: int, cell=None) -> float:
    """Angle i-j-k in degrees, at atom j."""
    u, v = mic_vector(pos[j], pos[i], cell), mic_vector(pos[j], pos[k], cell)
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu < 1e-12 or nv < 1e-12:
        return float("nan")
    cos = float(np.dot(u, v) / (nu * nv))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def dihedral(pos, i: int, j: int, k: int, l: int, cell=None) -> float:
    """Dihedral i-j-k-l in degrees, in (-180, 180]."""
    b1, b2, b3 = mic_vector(pos[i], pos[j], cell), mic_vector(pos[j], pos[k], cell), mic_vector(pos[k], pos[l], cell)
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    nb2 = np.linalg.norm(b2)
    if nb2 < 1e-12 or np.linalg.norm(n1) < 1e-12 or np.linalg.norm(n2) < 1e-12:
        return float("nan")
    m = np.cross(n1, b2 / nb2)
    return float(np.degrees(np.arctan2(np.dot(m, n2), np.dot(n1, n2))))


def measure(pos, indices, cell=None) -> dict | None:
    """Measure 2 (distance), 3 (angle) or 4 (dihedral) atoms; None for other counts."""
    idx = [int(i) for i in indices]
    if len(idx) == 2:
        return {"kind": "distance", "value": distance(pos, *idx, cell=cell), "unit": "Å", "indices": idx}
    if len(idx) == 3:
        return {"kind": "angle", "value": angle(pos, *idx, cell=cell), "unit": "°", "indices": idx}
    if len(idx) == 4:
        return {"kind": "dihedral", "value": dihedral(pos, *idx, cell=cell), "unit": "°", "indices": idx}
    return None


def measure_text(m: dict | None, symbols=None) -> str:
    """One line for the screen: what was measured, between which atoms (1-based), and the value."""
    if m is None:
        return ""
    names = [f"{symbols[i]}{i + 1}" if symbols is not None else str(i + 1) for i in m["indices"]]
    kind = {"distance": L("距離", "distance"), "angle": L("角度", "angle"), "dihedral": L("二面角", "dihedral")}[m["kind"]]
    digits = 3 if m["kind"] == "distance" else 1
    return f"{kind} {'-'.join(names)}: {m['value']:.{digits}f} {m['unit']}"


__all__ = ["mic_vector", "distance", "angle", "dihedral", "measure", "measure_text"]
