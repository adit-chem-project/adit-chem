
from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.data import atomic_numbers

from adit.builder.bases import limit
from adit.builder.model import (Adsorb, Box, Fix, MoleculeRef, RecipeError, Remove, Selection, Slab, SolventLayer, Solvate,
                               Substitute, Supercell, Vacuum)
from adit.lang import L
from adit.mixture import AMU_TO_G, Component, MixtureError, pack_failure_message, pack_molecules

_TOL = 1e-6


def prepared(atoms: Atoms) -> Atoms:
    a = atoms.copy()
    if "molecule" not in a.arrays:
        a.set_array("molecule", np.full(len(a), -1, dtype=int))
    if "fixed" not in a.arrays:
        fixed = np.zeros(len(a), dtype=bool)
        from ase.constraints import FixAtoms
        for c in a.constraints:
            if isinstance(c, FixAtoms):
                fixed[c.index] = True
        a.set_array("fixed", fixed)
    a.set_constraint()
    return a


def concat(a: Atoms, b: Atoms) -> Atoms:
    a, b = prepared(a), prepared(b)
    mb = b.get_array("molecule").copy()
    start = int(a.get_array("molecule").max(initial=-1)) + 1
    mb[mb >= 0] += start
    b.set_array("molecule", mb)
    out = a.copy()
    out.extend(b)
    return out


def _periodic(atoms: Atoms) -> bool:
    return bool(all(atoms.pbc)) and abs(np.linalg.det(np.asarray(atoms.cell))) > 1e-6


def _need_periodic(atoms: Atoms, what: str) -> None:
    if not _periodic(atoms):
        raise RecipeError(L(f"{what}は 3 方向とも周期の構造 (結晶・スラブ・溶液の箱) にだけ使えます。分子やクラスターなら、先に box で直方体の箱に入れてください",
                            f"{what} needs a structure periodic in all three directions (crystal, slab, liquid box). For a molecule or cluster, put it in a box first"))


def slab_cell_problem(atoms: Atoms) -> str | None:
    if not _periodic(atoms):
        return L("3 方向とも周期のセルではありません", "the cell is not periodic in all three directions")
    a, b, c = np.asarray(atoms.cell)
    for v, name in ((a, "a"), (b, "b")):
        cos = abs(float(np.dot(v, c))) / (np.linalg.norm(v) * np.linalg.norm(c))
        if cos > 1e-4:
            ang = float(np.degrees(np.arccos(min(1.0, cos))))
            return L(f"c が {name} に垂直ではありません ({name} と c のなす角 {ang:.2f}° または {180 - ang:.2f}°)",
                     f"c is not perpendicular to {name} (angle {ang:.2f} or {180 - ang:.2f} deg)")
    return None


def standardized(atoms: Atoms) -> Atoms:
    a = atoms.copy()
    rcell, _q = a.cell.standard_form()
    if not np.allclose(np.asarray(rcell), np.asarray(a.cell), atol=1e-8):
        a.set_cell(rcell, scale_atoms=True)
    return a


def unwrap_axis(atoms: Atoms, axis: int) -> Atoms:
    a = atoms.copy()
    f = a.get_scaled_positions(wrap=False)
    f[:, axis] %= 1.0
    s = np.sort(f[:, axis])
    gaps = np.diff(np.concatenate([s, [s[0] + 1.0]]))
    start = s[(int(np.argmax(gaps)) + 1) % len(s)]
    f[:, axis] = (f[:, axis] - start) % 1.0
    f[:, axis][f[:, axis] > 1.0 - 1e-9] = 0.0
    a.set_scaled_positions(f)
    return a


def make_molecules_whole(atoms: Atoms) -> Atoms:
    if "molecule" not in atoms.arrays or not _periodic(atoms):
        return atoms
    mol = atoms.get_array("molecule")
    if not (mol >= 0).any():
        return atoms
    a = atoms.copy()
    cell = np.asarray(a.cell); inv = np.linalg.inv(cell)
    pos = a.get_positions()
    ids, first = np.unique(mol[mol >= 0], return_index=True)
    ref_index = dict(zip(ids.tolist(), np.flatnonzero(mol >= 0)[first].tolist()))
    sel = np.flatnonzero(mol >= 0)
    ref = np.array([ref_index[int(m)] for m in mol[sel]])
    d = (pos[sel] - pos[ref]) @ inv
    d -= np.round(d)
    pos[sel] = pos[ref] + d @ cell
    a.set_positions(pos)
    return a


def split_molecules(atoms: Atoms) -> list[int]:
    if "molecule" not in atoms.arrays or not _periodic(atoms):
        return []
    mol = atoms.get_array("molecule")
    sel = np.flatnonzero(mol >= 0)
    if not len(sel):
        return []
    whole = make_molecules_whole(atoms).get_positions()
    bad = np.linalg.norm(whole[sel] - atoms.get_positions()[sel], axis=1) > 1e-6
    return sorted(set(mol[sel][bad].tolist()))


def planes(z: np.ndarray, tol: float) -> list[np.ndarray]:
    order = np.argsort(z, kind="stable")
    out = [[order[0]]] if len(order) else []
    for p, q in zip(order[:-1], order[1:]):
        if z[q] - z[p] <= tol:
            out[-1].append(q)
        else:
            out.append([q])
    return [np.array(g) for g in out]


def select(atoms: Atoms, sel: Selection) -> np.ndarray:
    n = len(atoms)
    m = np.ones(n, dtype=bool)
    if sel.elements:
        bad = [e for e in sel.elements if e not in atomic_numbers]
        if bad:
            raise RecipeError(L(f"元素記号として解釈できません: {bad}", f"not chemical symbols: {bad}"))
        m &= np.isin(np.array(atoms.get_chemical_symbols()), sel.elements)
    if sel.indices:
        bad = [i for i in sel.indices if not 0 <= i < n]
        if bad:
            raise RecipeError(L(f"原子の番号が範囲外です (0〜{n - 1}): {bad}", f"atom index out of range (0..{n - 1}): {bad}"))
        m &= np.isin(np.arange(n), sel.indices)
    z = atoms.positions[:, 2]
    if sel.z_min is not None:
        m &= z >= sel.z_min
    if sel.z_max is not None:
        m &= z <= sel.z_max
    return m


def describe(sel: Selection) -> str:
    parts = []
    if sel.elements:
        parts.append(L(f"元素 {', '.join(sel.elements)}", f"elements {', '.join(sel.elements)}"))
    if sel.indices:
        parts.append(L(f"番号 {sel.indices[:8]}{'…' if len(sel.indices) > 8 else ''}", f"indices {sel.indices[:8]}"))
    if sel.z_min is not None or sel.z_max is not None:
        lo = "" if sel.z_min is None else f"{sel.z_min:g}"
        hi = "" if sel.z_max is None else f"{sel.z_max:g}"
        parts.append(f"z {lo}〜{hi} Å")
    return L("、".join(parts) or "全原子", ", ".join(parts) or "all atoms")


def _choose(atoms: Atoms, step, what: str) -> np.ndarray:
    idx = np.flatnonzero(select(atoms, step.where))
    if step.count is not None:
        n = int(step.count)
    else:
        if not 0.0 <= step.fraction <= 1.0:
            raise RecipeError(L(f"割合は 0〜1 です (いま {step.fraction})", f"fraction must be 0..1 (now {step.fraction})"))
        n = int(round(step.fraction * len(idx)))
    if n <= 0:
        raise RecipeError(L(f"{what}原子の数が 0 です (条件 {describe(step.where)} に合う原子 {len(idx)} 個)。個数か割合を増やしてください",
                            f"{what}: the number of atoms is 0 ({len(idx)} atoms match {describe(step.where)}). Increase the count or fraction"))
    if n > len(idx):
        raise RecipeError(L(f"{what}: 条件 ({describe(step.where)}) に合う原子は {len(idx)} 個しかなく、{n} 個は選べません。条件を広げるか個数を減らしてください",
                            f"{what}: only {len(idx)} atoms match ({describe(step.where)}), cannot choose {n}. Widen the condition or reduce the count"))
    rng = np.random.default_rng(step.seed)
    return np.sort(rng.choice(idx, size=n, replace=False))


def _molecule(ref: MoleculeRef | Component) -> Atoms:
    from adit.mixture import _one_molecule

    return _one_molecule(Component(kind=ref.kind, ref=ref.ref))


def _mols(components: list[Component]) -> list[tuple[Atoms, int]]:
    if not components:
        raise RecipeError(L("成分がありません", "no components"))
    if any(c.count < 0 for c in components):
        raise RecipeError(L("個数に負の値があります", "a count is negative"))
    return [(_molecule(c), c.count) for c in components]


def _summary(mols: list[tuple[Atoms, int]]) -> str:
    return ", ".join(f"{m.get_chemical_formula()} × {n}" for m, n in mols if n > 0)


def _mass_g(mols: list[tuple[Atoms, int]]) -> float:
    return sum(float(m.get_masses().sum()) * n for m, n in mols) * AMU_TO_G


def in_plane_widths(cell) -> np.ndarray:
    cell = np.asarray(cell, dtype=float)
    area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    return np.array([area / np.linalg.norm(cell[1]), area / np.linalg.norm(cell[0])])


def fit_repeat(widths: np.ndarray, need: float, floor=(1, 1)) -> np.ndarray:
    return np.maximum(np.asarray(floor, dtype=int), np.ceil(need / widths - 1e-9).astype(int))


def cross_section_problem(cell, mols: list[tuple[Atoms, int]], min_distance: float) -> str | None:
    from adit.mixture import packing_need

    w = in_plane_widths(cell)
    need_min, thin_name, need_room, room_name = packing_need(mols, min_distance)
    if (w >= need_min - 1e-9).all():
        return None
    n = fit_repeat(w, need_room)
    return L(f"スラブの断面が狭すぎて、溶液の成分が周期の像と重ならずに入りません。周期の幅は a の方向 {w[0]:.2f} Å、b の方向 {w[1]:.2f} Å です。"
             f"{thin_name} がいちばん薄い向きでも「分子間の最短距離」{min_distance:g} Å を保つには、どちらも {need_min:.2f} Å 以上要ります。"
             f"この手順の前に「超格子」の手順を入れて、「繰り返し (a × b × c)」を {n[0]} × {n[1]} × 1 にしてください "
             f"(幅が {w[0] * n[0]:.2f} × {w[1] * n[1]:.2f} Å になり、いちばん大きい {room_name} も余裕をもって入ります)。"
             "「超格子」は「面で切る」より後ろに置いてください (ミラー指数はいまのセルの格子ベクトルに対する指数なので、"
             "先に超格子にすると同じ指数でも別の面になります)",
             f"the slab cross-section is too narrow for the solution to avoid its own periodic images. The periodic widths are "
             f"{w[0]:.2f} Å along a and {w[1]:.2f} Å along b, but {thin_name} needs at least {need_min:.2f} Å in both directions to keep "
             f"the minimum distance {min_distance:g} Å even in its thinnest orientation. Add a supercell step before this one with "
             f"repeat {n[0]} x {n[1]} x 1 (widths become {w[0] * n[0]:.2f} x {w[1] * n[1]:.2f} Å, enough room for the largest "
             f"component {room_name}). Put the supercell after the slab step, not before (Miller indices refer to the current cell "
             "vectors, so a supercell first turns the same indices into a different plane)")


def box_problem(cell, mols: list[tuple[Atoms, int]], min_distance: float, from_padding: bool) -> str | None:
    from adit.mixture import packing_need, perpendicular_widths

    w = perpendicular_widths(np.asarray(cell, dtype=float))
    need_min, thin_name, need_room, room_name = packing_need(mols, min_distance)
    if (w >= need_min - 1e-9).all():
        return None
    how = L("「余白」を大きくしてください", "increase the padding") if from_padding else \
        L("先に「超格子」の手順でセルを広げてください", "widen the cell with a supercell step first")
    return L(f"箱が狭すぎて、溶媒の成分が周期の像と重ならずに入りません。周期の幅は {w[0]:.2f} × {w[1]:.2f} × {w[2]:.2f} Å です。"
             f"{thin_name} がいちばん薄い向きでも「分子間の最短距離」{min_distance:g} Å を保つには、どの方向も {need_min:.2f} Å 以上要ります "
             f"(いちばん大きい {room_name} に余裕をみるなら {need_room:.2f} Å)。{how}",
             f"the box is too small for the solvent to avoid its own periodic images. The periodic widths are "
             f"{w[0]:.2f} x {w[1]:.2f} x {w[2]:.2f} Å, but {thin_name} needs at least {need_min:.2f} Å in every direction to keep the "
             f"minimum distance {min_distance:g} Å even in its thinnest orientation ({need_room:.2f} Å with room for the largest "
             f"component {room_name}). {how}")


def _fit_in_plane(a: Atoms, rep: tuple[int, int, int], st: Supercell) -> tuple[tuple[int, int, int], str]:
    from adit.mixture import packing_need

    _need_periodic(a, L("面内を溶液に合わせて広げる超格子", "a supercell fitted to the solution"))
    mols = _mols(list(st.fit_components))
    w = in_plane_widths(a.cell)
    _need_min, _thin, need, name = packing_need(mols, st.fit_min_distance)
    n = fit_repeat(w, need, floor=rep[:2])
    new = (int(n[0]), int(n[1]), rep[2])
    grown = L(f"面内を a {new[0]} 倍・b {new[1]} 倍にしました", f"repeated the surface {new[0]} times along a and {new[1]} times along b") \
        if (new[0], new[1]) != (rep[0], rep[1]) else \
        L(f"現在の繰り返し (a {new[0]} 倍・b {new[1]} 倍) で必要な幅があるため、断面は広げませんでした",
          f"kept the current repeats ({new[0]} along a and {new[1]} along b) because the cross-section was already wide enough")
    note = L(f"溶液の成分のうちいちばん大きい {name} が「分子間の最短距離」{st.fit_min_distance:g} Å を保って入るには周期の幅が {need:.2f} Å 要るので、"
             f"{grown} (幅 {w[0] * new[0]:.2f} × {w[1] * new[1]:.2f} Å)",
             f"the largest component {name} needs a periodic width of {need:.2f} Å to keep the minimum distance {st.fit_min_distance:g} Å, "
             f"so {grown} (widths {w[0] * new[0]:.2f} x {w[1] * new[1]:.2f} Å)")
    return new, note


def op_supercell(atoms: Atoms, st: Supercell) -> Atoms:
    a = prepared(atoms)
    note = ""
    if st.repeat is not None:
        rep = tuple(int(x) for x in st.repeat)
        if any(k < 1 for k in rep):
            raise RecipeError(L(f"繰り返しは 1 以上です: {rep}", f"repeats must be at least 1: {rep}"))
        if st.fit_components:
            rep, note = _fit_in_plane(a, rep, st)
        for k in range(3):
            if rep[k] > 1 and (not a.pbc[k] or np.linalg.norm(a.cell[k]) < _TOL):
                raise RecipeError(L(f"軸 {'abc'[k]} は周期ではないので繰り返せません (分子なら先に box)", f"axis {'abc'[k]} is not periodic, cannot repeat (use box first for a molecule)"))
        limit(len(a) * rep[0] * rep[1] * rep[2])
        out = renumber_images(a.repeat(rep), len(a), rep[0] * rep[1] * rep[2])
    else:
        _need_periodic(a, L("変換行列の超格子", "a supercell matrix"))
        P = np.array(st.matrix, dtype=int)
        det = int(round(np.linalg.det(P)))
        if det <= 0:
            raise RecipeError(L(f"変換行列の行列式が {det} です。1 以上 (右手系) になるようにしてください (負なら 1 行の符号を反転)",
                                f"the matrix determinant is {det}; it must be >= 1 (right-handed; flip the sign of one row if negative)"))
        limit(len(a) * det)
        from ase.build import make_supercell
        out = make_molecules_whole(renumber_images(make_supercell(a, P, order="cell-major"), len(a), det))
    if note:
        out.info = {**out.info, "recipe_note": note}
    return out


def renumber_images(out: Atoms, n_prim: int, n_images: int) -> Atoms:
    if "molecule" not in out.arrays or len(out) != n_prim * n_images:
        return out
    mol = out.get_array("molecule").copy()
    width = int(mol[:n_prim].max(initial=-1)) + 1
    if width > 0:
        img = np.repeat(np.arange(n_images), n_prim)
        out.set_array("molecule", np.where(mol >= 0, mol + img * width, mol))
    return out


def _fold_z(s: Atoms, tol: float = 0.2) -> Atoms:
    c = float(s.cell[2, 2])
    z = s.positions[:, 2]
    z[z > c - tol] -= c
    z -= z.min()
    return s


def list_terminations(atoms: Atoms, miller, tol: float = 0.2) -> list[str]:
    from ase.build import surface

    _need_periodic(atoms, L("面で切る手順", "slab"))
    one = _fold_z(surface(prepared(atoms), tuple(int(x) for x in miller), 1, vacuum=None, periodic=True), tol)
    groups = planes(one.positions[:, 2], tol)
    return [one[g].get_chemical_formula() for g in reversed(groups)]


def op_slab(atoms: Atoms, st: Slab) -> Atoms:
    from ase.build import surface

    _need_periodic(atoms, L("面で切る手順", "slab"))
    hkl = tuple(int(x) for x in st.miller)
    if hkl == (0, 0, 0):
        raise RecipeError(L("ミラー指数が (0 0 0) です", "Miller indices are (0 0 0)"))
    if st.layers < 1:
        raise RecipeError(L("層の数は 1 以上です", "layers must be at least 1"))
    if st.vacuum < 0:
        raise RecipeError(L("真空は 0 以上です", "vacuum must be >= 0"))
    a = prepared(atoms)
    try:
        one = _fold_z(surface(a, hkl, 1, vacuum=None, periodic=True))
    except Exception as ex:
        raise RecipeError(L(f"面 {hkl} で切れません: {ex}", f"cannot cut the plane {hkl}: {ex}")) from ex
    limit(len(one) * (st.layers + 1))
    groups1 = planes(one.positions[:, 2], 0.2)
    k = len(groups1)
    terms = [one[g].get_chemical_formula() for g in reversed(groups1)]
    if not 0 <= st.termination < k:
        listing = ", ".join(f"{i} = {f}" for i, f in enumerate(terms))
        raise RecipeError(L(f"終端の番号 {st.termination} はありません。この面の 1 周期の原子面は上から {listing} です",
                            f"no termination {st.termination}. The atomic planes of one period, from the top: {listing}"))
    period = np.asarray(one.cell[2]).copy()
    if st.termination == 0:
        s = _fold_z(surface(a, hkl, st.layers, vacuum=None, periodic=True))
    else:
        s = _fold_z(surface(a, hkl, st.layers + 1, vacuum=None, periodic=True))
        g = planes(s.positions[:, 2], 0.2)
        if len(g) != k * (st.layers + 1):
            raise RecipeError(L("原子面の数が周期と合わず、終端を選べません (面が波打っている構造)。termination = 0 で切り、原子を抜く手順で整えてください",
                                "the number of atomic planes does not match the period (corrugated planes); use termination 0 and remove atoms instead"))
        drop = np.concatenate(g[:k - st.termination] + g[len(g) - st.termination:])
        del s[drop]
        s.set_cell([s.cell[0], s.cell[1], period * st.layers], scale_atoms=False)
        s.positions[:, 2] -= s.positions[:, 2].min()
    if st.vacuum > 0:
        s.center(vacuum=st.vacuum, axis=2)
    s.pbc = (True, True, True)
    s.info.pop("adsorbate_info", None)
    return s


def op_vacuum(atoms: Atoms, st: Vacuum) -> Atoms:
    _need_periodic(atoms, L("真空層の手順", "vacuum"))
    if st.thickness < 0:
        raise RecipeError(L("真空の厚みは 0 以上です", "vacuum thickness must be >= 0"))
    a = unwrap_axis(prepared(atoms), st.axis)
    a.center(vacuum=st.thickness / 2.0, axis=st.axis)
    return a


def _orthogonal_matrix(cell: np.ndarray, max_multiple: int, nmax: int = 3) -> np.ndarray | None:
    uc = cell / np.linalg.norm(cell, axis=1)[:, None]
    if np.abs(uc @ uc.T - np.eye(3)).max() < 1e-5:
        return np.eye(3, dtype=int)
    rng = range(-nmax, nmax + 1)
    ints = np.array([(i, j, k) for i in rng for j in rng for k in rng if (i, j, k) != (0, 0, 0)])
    vec = ints @ cell
    ln = np.linalg.norm(vec, axis=1)
    unit = vec / ln[:, None]
    orth = np.abs(unit @ unit.T) < 1e-5
    best, best_key = None, None
    for u in np.argsort(ln, kind="stable"):
        if best_key is not None and best_key[0] == 1 and 3 * ln[u] > best_key[1]:
            break
        vs = np.flatnonzero(orth[u])
        for v in vs:
            ws = vs[orth[v, vs]]
            if not len(ws):
                continue
            dets = np.rint(np.cross(ints[v], ints[ws]) @ ints[u]).astype(int)  # det[u; v; w] = u · (v × w)
            ok = (dets != 0) & (np.abs(dets) <= max_multiple)
            for w, det in zip(ws[ok], dets[ok]):
                M = np.array([ints[u], ints[v], ints[w] * (1 if det > 0 else -1)])
                key = (abs(int(det)), float(ln[u] + ln[v] + ln[w]))
                if best_key is None or key < best_key:
                    best, best_key = M, key
    return best


def op_box(atoms: Atoms, st: Box) -> Atoms:
    a = prepared(atoms)
    if not any(a.pbc) or abs(np.linalg.det(np.asarray(a.cell))) < 1e-6 and not all(a.pbc):
        if st.padding < 0:
            raise RecipeError(L("padding は 0 以上です", "padding must be >= 0"))
        ext = np.ptp(a.positions, axis=0)
        if st.lengths is not None:
            size = np.asarray(st.lengths, dtype=float)
            if (size <= ext).any():
                raise RecipeError(L(f"箱 {tuple(size)} Å が構造の広がり {tuple(np.round(ext, 2))} Å より小さく、入りません",
                                    f"the box {tuple(size)} Å is smaller than the structure extent {tuple(np.round(ext, 2))} Å"))
        else:
            size = ext + 2 * st.padding
            if (size <= 0).any():
                raise RecipeError(L("padding を 0 より大きくしてください (1 原子の構造の箱の大きさが 0 になります)", "padding must be positive for a single atom"))
        a.set_cell(np.diag(size)); a.pbc = (True, True, True); a.center()
        return a
    if not all(a.pbc):
        raise RecipeError(L("一部の方向だけ周期の構造です。先に vacuum でその方向に真空を置いてください", "the structure is periodic in only some directions; add vacuum first"))
    cell = np.asarray(a.cell)
    P = _orthogonal_matrix(cell, st.max_multiple)
    if P is None:
        raise RecipeError(L(f"元のセルの {st.max_multiple} 倍までに、直交する格子の取り方が見つかりません。max_multiple を大きくしてください",
                            f"no orthogonal cell within {st.max_multiple} times the original; increase max_multiple"))
    det = int(round(np.linalg.det(P)))
    limit(len(a) * det)
    from ase.build import make_supercell
    out = renumber_images(make_supercell(a, P, order="cell-major"), len(a), det)
    rcell, _q = out.cell.standard_form()
    out.set_cell(rcell, scale_atoms=True)
    return make_molecules_whole(out)


def _rotate_to(v: np.ndarray, target: np.ndarray) -> np.ndarray:
    v = v / np.linalg.norm(v); t = target / np.linalg.norm(target)
    c = float(np.dot(v, t)); axis = np.cross(v, t); s = float(np.linalg.norm(axis))
    if s < 1e-12:
        return np.eye(3) if c > 0 else np.diag([1.0, -1.0, -1.0])
    k = axis / s
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + s * K + (1 - c) * K @ K


def op_adsorb(atoms: Atoms, st: Adsorb) -> Atoms:
    a = prepared(atoms)
    mol = _molecule(st.molecule)
    if not 0 <= st.down_atom < len(mol):
        raise RecipeError(L(f"下に向ける原子の番号 {st.down_atom} が分子 {mol.get_chemical_formula()} の範囲 (0〜{len(mol) - 1}) の外です",
                            f"down_atom {st.down_atom} is outside the molecule {mol.get_chemical_formula()} (0..{len(mol) - 1})"))
    p = mol.get_positions()
    v = p.mean(axis=0) - p[st.down_atom]
    p = p - p[st.down_atom]
    if np.linalg.norm(v) > 1e-8:
        p = p @ _rotate_to(v, np.array([0.0, 0.0, 1.0])).T
    if st.site is not None:
        info = a.info.get("adsorbate_info") or {}
        sites = info.get("sites") or {}
        if not sites:
            raise RecipeError(L("吸着サイトの名前は、ASE の面関数で作ったスラブ (土台が surface で、そのあと変換行列の超格子や面の切り直しをしていないもの) にだけあります。"
                                "above_atom (原子の番号) か xy で置き場所を指定してください",
                                "site names exist only for slabs made by ASE surface functions (surface base, not re-cut or matrix-supercelled). "
                                "Use above_atom (an atom index) or xy instead"))
        if st.site not in sites:
            raise RecipeError(L(f"この面に席 {st.site!r} はありません。使える席: {', '.join(sites)}", f"no site {st.site!r} on this surface. Available: {', '.join(sites)}"))
        xy = np.dot(np.asarray(sites[st.site], dtype=float), np.asarray(info.get("cell", np.asarray(a.cell)[:2, :2])))
        top = info.get("top layer atom index", int(np.argmax(a.positions[:, 2])))
        z = float(a.positions[top, 2])
    elif st.above_atom is not None:
        if not 0 <= st.above_atom < len(a):
            raise RecipeError(L(f"原子の番号 {st.above_atom} が範囲 (0〜{len(a) - 1}) の外です", f"atom index {st.above_atom} is out of range (0..{len(a) - 1})"))
        xy = a.positions[st.above_atom, :2]
        z = float(a.positions[st.above_atom, 2])
    else:
        xy = np.asarray(st.xy, dtype=float)
        z = float(a.positions[:, 2].max())
    target = np.array([xy[0], xy[1], z + st.height])
    ads = Atoms(symbols=mol.get_chemical_symbols(), positions=p + target, cell=a.cell, pbc=a.pbc)
    ads.set_array("molecule", np.zeros(len(ads), dtype=int))
    out = concat(a, ads)
    out.info = dict(a.info)
    return out


def op_remove(atoms: Atoms, st: Remove) -> Atoms:
    a = prepared(atoms)
    idx = _choose(a, st, L("原子を抜く手順", "remove"))
    del a[idx]
    return a


def op_substitute(atoms: Atoms, st: Substitute) -> Atoms:
    if st.to not in atomic_numbers:
        raise RecipeError(L(f"置き換える元素 {st.to!r} は元素記号ではありません", f"{st.to!r} is not a chemical symbol"))
    a = prepared(atoms)
    idx = _choose(a, st, L("元素の置換", "substitute"))
    numbers = a.get_atomic_numbers()
    numbers[idx] = atomic_numbers[st.to]
    a.set_atomic_numbers(numbers)
    return a


def op_solvent_layer(atoms: Atoms, st: SolventLayer) -> Atoms:
    why = slab_cell_problem(atoms)
    if why:
        raise RecipeError(L(f"溶液の層は、a, b が周期で c が面に垂直なスラブの上にだけ置けます ({why})。"
                            "slab 手順で切ったスラブか、surface の土台を使ってください。ファイルの構造なら、c を面に垂直に取り直してから読み込んでください",
                            f"a solvent layer needs a slab with periodic a, b and c normal to the surface ({why}). "
                            "Use a slab step or a surface base; for a file, re-define c to be normal to the surface"))
    for name in ("gap", "vacuum"):
        if getattr(st, name) < 0:
            raise RecipeError(L(f"{name} は 0 以上です", f"{name} must be >= 0"))
    if st.thickness is not None and st.thickness <= 0:
        raise RecipeError(L("厚みは 0 より大きくしてください", "thickness must be positive"))
    mols = _mols(st.components)
    if sum(n for _, n in mols) <= 0:
        raise RecipeError(L("個数が 0 です", "the counts are all zero"))
    limit(len(atoms) + sum(len(m) * n for m, n in mols))
    why = cross_section_problem(atoms.cell, mols, st.min_distance)
    if why:
        raise RecipeError(why)
    a = unwrap_axis(standardized(prepared(atoms)), 2)
    pos = a.get_positions(); pos[:, 2] -= pos[:, 2].min(); a.set_positions(pos)
    cell = np.asarray(a.cell).copy()
    area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    if st.thickness is not None:
        T = float(st.thickness)
    else:
        if st.density_g_cm3 <= 0:
            raise RecipeError(L("密度は 0 より大きくしてください", "density must be positive"))
        T = _mass_g(mols) / st.density_g_cm3 * 1e24 / area
    top = float(pos[:, 2].max())
    lo = top + st.gap
    hi = lo + T
    cell[2] = (0.0, 0.0, hi + st.gap + st.vacuum)
    a.set_cell(cell, scale_atoms=False)
    try:
        packed = pack_molecules(mols, cell, np.random.default_rng(st.seed), min_distance=st.min_distance,
                                steps=st.max_tries // 4, z_band=(lo, hi), obstacles=a.get_positions())
    except MixtureError as ex:
        raise RecipeError(str(ex)) from ex
    if not packed.ok:
        raise RecipeError(pack_failure_message(_summary(mols), packed, st.min_distance))
    out = concat(a, packed.to_atoms(cell))
    out.info = {"solvent_band": (lo, hi)}
    return out


def op_solvate(atoms: Atoms, st: Solvate) -> Atoms:
    a = prepared(atoms)
    from_padding = not any(a.pbc)
    if not any(a.pbc):
        if st.padding <= 0:
            raise RecipeError(L("padding は 0 より大きくしてください", "padding must be positive"))
        size = np.ptp(a.positions, axis=0) + 2 * st.padding
        a.set_cell(np.diag(size)); a.pbc = (True, True, True); a.center()
    elif not _periodic(a):
        raise RecipeError(L("一部の方向だけ周期の構造は溶媒和できません。先に vacuum でその方向に真空を置いてください",
                            "cannot solvate a structure periodic in only some directions; add vacuum first"))
    cell = np.asarray(a.cell).copy()
    comps = list(st.components)
    zero = [i for i, c in enumerate(comps) if c.count == 0]
    if len(zero) > 1:
        raise RecipeError(L("個数を密度から決める成分 (count = 0) は 1 つだけにしてください", "only one component may have count 0 (filled to the density)"))
    mols = _mols(comps)
    if zero:
        if st.density_g_cm3 <= 0:
            raise RecipeError(L("密度は 0 より大きくしてください", "density must be positive"))
        k = zero[0]
        vol_cm3 = abs(float(np.linalg.det(cell))) * 1e-24
        rest = float(a.get_masses().sum()) * AMU_TO_G + _mass_g([m for i, m in enumerate(mols) if i != k])
        m1 = float(mols[k][0].get_masses().sum()) * AMU_TO_G
        n = int((st.density_g_cm3 * vol_cm3 - rest) // m1)
        if n <= 0:
            raise RecipeError(L(f"箱 ({' × '.join(f'{x:.1f}' for x in np.linalg.norm(cell, axis=1))} Å) が溶質だけで密度 {st.density_g_cm3} g/cm³ に達し、"
                                f"{mols[k][0].get_chemical_formula()} を入れる余地がありません。padding を大きくしてください",
                                f"the box is already at {st.density_g_cm3} g/cm3 with the solute alone; no room for {mols[k][0].get_chemical_formula()}. Increase padding"))
        mols[k] = (mols[k][0], n)
    if sum(n for _, n in mols) <= 0:
        raise RecipeError(L("個数が 0 です", "the counts are all zero"))
    limit(len(a) + sum(len(m) * n for m, n in mols))
    why = box_problem(cell, mols, st.min_distance, from_padding)
    if why:
        raise RecipeError(why)
    try:
        packed = pack_molecules(mols, cell, np.random.default_rng(st.seed), min_distance=st.min_distance,
                                steps=st.max_tries // 4, obstacles=a.get_positions())
    except MixtureError as ex:
        raise RecipeError(str(ex)) from ex
    if not packed.ok:
        raise RecipeError(pack_failure_message(_summary(mols), packed, st.min_distance))
    return concat(a, packed.to_atoms(cell))


def op_fix(atoms: Atoms, st: Fix) -> Atoms:
    a = prepared(atoms)
    m = select(a, st.where)
    if st.bottom_layers < 0:
        raise RecipeError(L("層の数は 0 以上です", "bottom_layers must be >= 0"))
    if st.bottom_layers > 0:
        layers = planes(a.positions[:, 2], st.layer_tolerance)
        if st.bottom_layers > len(layers):
            raise RecipeError(L(f"層は {len(layers)} 枚しかありません (z の差 {st.layer_tolerance} Å で分けたとき)", f"there are only {len(layers)} layers"))
        mb = np.zeros(len(a), dtype=bool)
        mb[np.concatenate(layers[:st.bottom_layers])] = True
        m &= mb
    if not m.any():
        raise RecipeError(L(f"固定する条件 ({describe(st.where)}) に合う原子がありません", f"no atoms match the fix condition ({describe(st.where)})"))
    fixed = a.get_array("fixed").copy()
    fixed |= m
    a.set_array("fixed", fixed)
    return a


OPS = {"supercell": op_supercell, "slab": op_slab, "vacuum": op_vacuum, "box": op_box, "adsorb": op_adsorb,
       "remove": op_remove, "substitute": op_substitute, "solvent_layer": op_solvent_layer, "solvate": op_solvate, "fix": op_fix}
