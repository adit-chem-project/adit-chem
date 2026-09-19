
from __future__ import annotations

import numpy as np

from adit.errors import AditValueError
from adit.lang import L


class GeometryError(AditValueError):
    pass


def _relative(a: np.ndarray, b: np.ndarray, cell) -> np.ndarray:
    diff = b - a
    if cell is None:
        return diff
    from adit.analysis.compute import _mic_step

    return _mic_step(np.atleast_2d(diff), np.asarray(cell, dtype=float)).reshape(diff.shape)


def distance_series(frames, i: int, j: int, cell=None) -> np.ndarray:
    out = []
    for fr in frames:
        pos = fr.get_positions()
        c = np.asarray(fr.cell, dtype=float) if any(fr.pbc) else cell
        out.append(float(np.linalg.norm(_relative(pos[i], pos[j], c))))
    return np.array(out)


def angle_series(frames, i: int, j: int, k: int, cell=None) -> np.ndarray:
    out = []
    for fr in frames:
        pos = fr.get_positions()
        c = np.asarray(fr.cell, dtype=float) if any(fr.pbc) else cell
        u, v = _relative(pos[j], pos[i], c), _relative(pos[j], pos[k], c)
        cos = float(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)))
        out.append(float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))))
    return np.array(out)


def dihedral_series(frames, i: int, j: int, k: int, l: int, cell=None) -> np.ndarray:
    out = []
    for fr in frames:
        pos = fr.get_positions()
        c = np.asarray(fr.cell, dtype=float) if any(fr.pbc) else cell
        b1, b2, b3 = _relative(pos[i], pos[j], c), _relative(pos[j], pos[k], c), _relative(pos[k], pos[l], c)
        n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
        m = np.cross(n1, b2 / np.linalg.norm(b2))
        out.append(float(np.degrees(np.arctan2(np.dot(m, n2), np.dot(n1, n2)))))
    return np.array(out)


def kabsch_rmsd(reference: np.ndarray, moving: np.ndarray, *, superpose: bool = True) -> tuple[float, np.ndarray]:
    ref = np.asarray(reference, dtype=float)
    mov = np.asarray(moving, dtype=float)
    if ref.shape != mov.shape or ref.ndim != 2 or ref.shape[1] != 3:
        raise GeometryError(L("同じ形 (原子数, 3) の座標が 2 つ要ります", "two coordinate arrays of the same shape (atoms, 3) are required"))
    if not superpose:
        return float(np.sqrt(np.mean(np.sum((mov - ref) ** 2, axis=1)))), mov
    p, q = ref - ref.mean(axis=0), mov - mov.mean(axis=0)
    u, _, vt = np.linalg.svd(q.T @ p)
    d = np.sign(np.linalg.det(u @ vt))
    rot = u @ np.diag([1.0, 1.0, d]) @ vt
    fitted = q @ rot
    return float(np.sqrt(np.mean(np.sum((fitted - p) ** 2, axis=1)))), fitted + ref.mean(axis=0)


def rmsd_series(frames, indices=None, *, reference=0, superpose: bool = True) -> np.ndarray:
    frames = list(frames)
    if not frames:
        raise GeometryError(L("フレームがありません", "there are no frames"))
    idx = np.arange(len(frames[0])) if indices is None else np.asarray(indices, dtype=int)
    ref = frames[reference].get_positions()[idx]
    return np.array([kabsch_rmsd(ref, fr.get_positions()[idx], superpose=superpose)[0] for fr in frames])


def rmsf(frames, indices=None, *, superpose: bool = True) -> dict:
    frames = list(frames)
    if len(frames) < 2:
        raise GeometryError(L("RMSF にはフレームが 2 つ以上要ります", "RMSF needs at least two frames"))
    idx = np.arange(len(frames[0])) if indices is None else np.asarray(indices, dtype=int)
    ref = frames[0].get_positions()[idx]
    aligned = []
    for fr in frames:
        pos = fr.get_positions()[idx]
        aligned.append(kabsch_rmsd(ref, pos, superpose=superpose)[1] if superpose else pos)
    stack = np.array(aligned)
    mean = stack.mean(axis=0)
    values = np.sqrt(np.mean(np.sum((stack - mean) ** 2, axis=2), axis=0))
    syms = np.array(frames[0].get_chemical_symbols())[idx]
    return {"index": (idx + 1).tolist(), "symbol": syms.tolist(), "rmsf_A": values.tolist(),
            "superposed": bool(superpose), "n_frames": len(frames),
            "note": L("平均構造からの揺らぎです。重ね合わせ (Kabsch) をしてから測っています。"
                      "大きい・小さいの判定はしていません。",
                      "fluctuation around the mean structure, measured after a Kabsch superposition; "
                      "no judgment of large or small is made.")}


def radius_of_gyration(frames, indices=None, *, mass_weighted: bool = True) -> dict:
    frames = list(frames)
    if not frames:
        raise GeometryError(L("フレームがありません", "there are no frames"))
    idx = np.arange(len(frames[0])) if indices is None else np.asarray(indices, dtype=int)
    values, principal = [], []
    for fr in frames:
        pos = fr.get_positions()[idx]
        w = fr.get_masses()[idx] if mass_weighted else np.ones(len(idx))
        center = (w[:, None] * pos).sum(axis=0) / w.sum()
        d = pos - center
        values.append(float(np.sqrt((w * np.sum(d ** 2, axis=1)).sum() / w.sum())))
        tensor = (w[:, None, None] * np.einsum("ni,nj->nij", d, d)).sum(axis=0) / w.sum()
        principal.append(np.sqrt(np.sort(np.linalg.eigvalsh(tensor))[::-1]).tolist())
    return {"rg_A": values, "principal_A": principal, "mass_weighted": bool(mass_weighted),
            "n_atoms": int(len(idx)),
            "note": L("慣性半径 Rg [Å] です (質量で重み付け)。principal_A は慣性テンソルの主値の平方根で、"
                      "大きい順に 3 つ (形の偏りを見るときに使います)。広がりの良し悪しは判定していません。",
                      "radius of gyration in Å (mass weighted); principal_A holds the square roots of the "
                      "eigenvalues of the gyration tensor, largest first. No judgment about the size is made.")}


def parse_atom_list(text: str, expected: int) -> list[int]:
    parts = [p for p in text.replace(" ", "").split(",") if p]
    if len(parts) != expected:
        raise GeometryError(L(f"原子を {expected} 個、カンマ区切りで指定してください (例 1,2{',3' if expected > 2 else ''})",
                              f"give {expected} atoms separated by commas (e.g. 1,2{',3' if expected > 2 else ''})"))
    out = []
    for p in parts:
        try:
            value = int(p)
        except ValueError as ex:
            raise GeometryError(L(f"原子の番号を読めません: {p!r}", f"cannot read the atom index {p!r}")) from ex
        if value < 1:
            raise GeometryError(L("原子の番号は 1 から始まります", "atom indices start at 1"))
        out.append(value - 1)
    return out
