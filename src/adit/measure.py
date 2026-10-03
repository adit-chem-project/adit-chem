"""Distance, angle and dihedral between chosen atoms, with the minimum image in periodic cells."""

from __future__ import annotations

import numpy as np

from adit.lang import L


def _periodic_axes(pbc) -> tuple[bool, bool, bool]:
    if pbc is None:
        return (True, True, True)
    flags = np.asarray(pbc, dtype=bool).reshape(-1)
    if flags.size == 1:
        return (bool(flags[0]),) * 3
    return tuple(bool(x) for x in flags[:3])


def mic_vector(a, b, cell=None, pbc=None) -> np.ndarray:
    """Vector from a to b; the shortest image along the periodic directions when a cell is given."""
    d = np.asarray(b, dtype=float) - np.asarray(a, dtype=float)
    if cell is None:
        return d
    axes = _periodic_axes(pbc)
    if not any(axes):
        return d
    c = np.asarray(cell, dtype=float)
    if abs(np.linalg.det(c)) < 1e-12:
        return d
    frac = np.linalg.solve(c.T, d)
    frac -= np.round(frac) * np.asarray(axes, dtype=float)
    # Only the periodic directions may be shifted by a lattice vector (a slab keeps its vacuum).
    shifts = [(-1, 0, 1) if on else (0,) for on in axes]
    best, best_len = None, None
    for i in shifts[0]:
        for j in shifts[1]:
            for k in shifts[2]:
                v = (frac + (i, j, k)) @ c
                n = float(v @ v)
                if best_len is None or n < best_len:
                    best, best_len = v, n
    return best


def distance(pos, i: int, j: int, cell=None, pbc=None) -> float:
    return float(np.linalg.norm(mic_vector(pos[i], pos[j], cell, pbc)))


def angle(pos, i: int, j: int, k: int, cell=None, pbc=None) -> float:
    """Angle i-j-k in degrees, at atom j."""
    u, v = mic_vector(pos[j], pos[i], cell, pbc), mic_vector(pos[j], pos[k], cell, pbc)
    nu, nv = np.linalg.norm(u), np.linalg.norm(v)
    if nu < 1e-12 or nv < 1e-12:
        return float("nan")
    cos = float(np.dot(u, v) / (nu * nv))
    return float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0))))


def dihedral(pos, i: int, j: int, k: int, l: int, cell=None, pbc=None) -> float:
    """Dihedral i-j-k-l in degrees, in (-180, 180]."""
    b1, b2, b3 = (mic_vector(pos[i], pos[j], cell, pbc), mic_vector(pos[j], pos[k], cell, pbc),
                  mic_vector(pos[k], pos[l], cell, pbc))
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    nb2 = np.linalg.norm(b2)
    if nb2 < 1e-12 or np.linalg.norm(n1) < 1e-12 or np.linalg.norm(n2) < 1e-12:
        return float("nan")
    m = np.cross(n1, b2 / nb2)
    return float(np.degrees(np.arctan2(np.dot(m, n2), np.dot(n1, n2))))


def measure(pos, indices, cell=None, pbc=None) -> dict | None:
    """Measure 2 (distance), 3 (angle) or 4 (dihedral) atoms; None for other counts."""
    idx = [int(i) for i in indices]
    if len(idx) == 2:
        return {"kind": "distance", "value": distance(pos, *idx, cell=cell, pbc=pbc), "unit": "Å", "indices": idx}
    if len(idx) == 3:
        return {"kind": "angle", "value": angle(pos, *idx, cell=cell, pbc=pbc), "unit": "°", "indices": idx}
    if len(idx) == 4:
        return {"kind": "dihedral", "value": dihedral(pos, *idx, cell=cell, pbc=pbc), "unit": "°", "indices": idx}
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
