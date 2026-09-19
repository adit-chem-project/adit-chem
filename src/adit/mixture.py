
from __future__ import annotations

from adit.progress import report

from adit.errors import AditError
import json
import math
from typing import Literal

import numpy as np
from ase import Atoms
from pydantic import BaseModel, Field

from adit.lang import L

AVOGADRO = 6.02214076e23
A3_TO_L = 1e-27  # 1 Å³ = 1e-27 L
AMU_TO_G = 1.66053907e-24


class MixtureError(AditError):
    pass


class Component(BaseModel):

    kind: Literal["preset", "smiles", "file"] = "preset"
    ref: str
    count: int = 1
    charge: int = 0
    label: str = ""

    def name(self) -> str:
        return self.label or self.ref


class MixtureSpec(BaseModel):
    components: list[Component] = Field(default_factory=list)
    box_a: float = 0.0
    density_g_cm3: float = 1.0
    min_distance: float = 2.0
    seed: int = 0
    max_tries: int = 2000
    cell: list[tuple[float, float, float]] | None = None
    z_band: tuple[float, float] | None = None

    def total_charge(self) -> int:
        return sum(c.count * c.charge for c in self.components)

    def to_ref(self) -> str:
        return self.model_dump_json(exclude_none=True)

    @classmethod
    def from_ref(cls, ref: str) -> "MixtureSpec":
        try:
            return cls.model_validate_json(ref)
        except ValueError as ex:
            raise MixtureError(L(f"混合物の指定を読めません: {ex}", f"cannot read the mixture specification: {ex}")) from ex

    def summary(self) -> str:
        return ", ".join(f"{c.name()} × {c.count}" for c in self.components)


def _one_molecule(c: Component) -> Atoms:
    from adit.structure import from_file, from_preset, from_smiles

    if c.kind == "preset":
        return from_preset(c.ref)
    if c.kind == "smiles":
        return from_smiles(c.ref)
    return from_file(c.ref)


def box_edge_for_density(mols: list[tuple[Atoms, int]], density_g_cm3: float) -> float:
    if density_g_cm3 <= 0:
        raise MixtureError(L("密度は 0 より大きい値にしてください", "density must be positive"))
    mass_g = sum(float(m.get_masses().sum()) * n for m, n in mols) * AMU_TO_G
    volume_cm3 = mass_g / density_g_cm3
    return float((volume_cm3 * 1e24) ** (1.0 / 3.0))  # cm^3 → Å³


def concentration_mol_per_l(count: int, box_a: float) -> float:
    return count / (AVOGADRO * box_a ** 3 * A3_TO_L) if box_a > 0 else float("nan")


def _random_rotation(rng: np.random.Generator) -> np.ndarray:
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a * a + b * b - c * c - d * d, 2 * (b * c - a * d), 2 * (b * d + a * c)],
                     [2 * (b * c + a * d), a * a - b * b + c * c - d * d, 2 * (c * d - a * b)],
                     [2 * (b * d - a * c), 2 * (c * d + a * b), a * a - b * b - c * c + d * d]])


def _min_image_dist(p: np.ndarray, others: np.ndarray, box: float) -> float:
    if len(others) == 0:
        return np.inf
    d = p[:, None, :] - others[None, :, :]
    d -= box * np.round(d / box)
    return float(np.sqrt((d ** 2).sum(axis=-1)).min())


def _candidate_pairs(pos: np.ndarray, mol_id: np.ndarray, box: float, cutoff: float):
    from ase.neighborlist import neighbor_list

    probe = Atoms(f"H{len(pos)}", positions=pos, cell=[box, box, box], pbc=True)
    i, j = neighbor_list("ij", probe, cutoff=cutoff)
    sel = (mol_id[i] != mol_id[j]) & (i < j)
    return i[sel], j[sel]


def _mic_vectors(pos: np.ndarray, i: np.ndarray, j: np.ndarray, box: float):
    dv = pos[j] - pos[i]
    dv -= box * np.round(dv / box)
    return dv, np.sqrt((dv * dv).sum(axis=1))


def min_inter_distance(pos: np.ndarray, mol_id: np.ndarray, box: float, cutoff: float) -> float:
    i, j = _candidate_pairs(pos, mol_id, box, cutoff)
    if len(i) == 0:
        return float("inf")
    return float(_mic_vectors(pos, i, j, box)[1].min())


def _relax(pos: np.ndarray, centers: np.ndarray, local: np.ndarray, mol_id: np.ndarray, box: float, r0: float, steps: int):
    damp, target, skin = 0.7, r0 * 1.03, 0.8
    idx_i = idx_j = None
    for it in range(max(50, steps)):
        if idx_i is None or it % 8 == 0:
            idx_i, idx_j = _candidate_pairs(pos, mol_id, box, target + skin)
        if len(idx_i) == 0:
            return pos, True
        dv, d = _mic_vectors(pos, idx_i, idx_j, box)
        close = (d < target) & (d > 1e-9)
        if not close.any():
            idx_i, idx_j = _candidate_pairs(pos, mol_id, box, target)
            if len(idx_i) == 0:
                return pos, True
            continue
        i, j, dv, d = idx_i[close], idx_j[close], dv[close], d[close]
        push = ((target - d) * damp / 2.0)[:, None] * (dv / d[:, None])
        shift = np.zeros_like(centers)
        np.add.at(shift, mol_id[i], -push)
        np.add.at(shift, mol_id[j], +push)
        big = np.linalg.norm(shift, axis=1)
        over = big > r0
        if over.any():
            shift[over] *= (r0 / big[over])[:, None]
        centers += shift
        centers -= box * np.floor(centers / box)
        pos = local + centers[mol_id]
    return pos, min_inter_distance(pos, mol_id, box, r0) >= r0 - 1e-6


def build_mixture(spec: MixtureSpec) -> Atoms:
    if not spec.components:
        raise MixtureError(L("成分がありません", "no components"))
    if any(c.count < 0 for c in spec.components):
        raise MixtureError(L("個数に負の値があります", "a count is negative"))
    mols = [(_one_molecule(c), c.count) for c in spec.components]
    total = sum(n for _, n in mols)
    if total <= 0:
        raise MixtureError(L("個数が 0 です", "the counts are all zero"))
    check_atom_limit(sum(len(m) * n for m, n in mols))
    if spec.cell is not None or spec.z_band is not None:
        if spec.cell is None:
            raise MixtureError(L("z の帯を使うときはセル (cell) も指定してください", "z_band needs a cell"))
        cell = np.asarray(spec.cell, dtype=float)
        packed = pack_molecules(mols, cell, np.random.default_rng(spec.seed), min_distance=spec.min_distance,
                                steps=spec.max_tries // 4, z_band=spec.z_band)
        if not packed.ok:
            raise MixtureError(pack_failure_message(spec.summary(), packed, spec.min_distance))
        return packed.to_atoms(cell)
    box = spec.box_a if spec.box_a > 0 else box_edge_for_density(mols, spec.density_g_cm3)
    need_min, thin_name, need_room, _room_name = packing_need(mols, spec.min_distance)
    if box < need_min:
        raise MixtureError(L(f"箱の一辺 {box:.2f} Å が狭すぎて、{thin_name} が周期の像と重ならずに入りません "
                             f"(いちばん薄い向きで測っても、分子間の最短距離 {spec.min_distance:g} Å を保つには一辺 {need_min:.2f} Å 以上 "
                             f"要ります。余裕をみるなら {need_room:.2f} Å)。一辺を大きくする、密度を下げる、個数を増やす、"
                             "最短距離を短くする、のどれかにしてください",
                             f"the box edge {box:.2f} Å is too small for {thin_name} to avoid its own periodic image "
                             f"(even in its thinnest direction it needs an edge of {need_min:.2f} Å to keep the minimum distance "
                             f"{spec.min_distance:g} Å; {need_room:.2f} Å with room to spare). Enlarge the edge, lower the density, "
                             "increase the counts, or lower the minimum distance"))
    rng = np.random.default_rng(spec.seed)

    k = max(1, int(math.ceil(total ** (1.0 / 3.0))))
    step = box / k
    sites = np.array([(a + 0.5, b + 0.5, c + 0.5) for a in range(k) for b in range(k) for c in range(k)]) * step
    rng.shuffle(sites)

    local_parts: list[np.ndarray] = []
    centers: list[np.ndarray] = []
    symbols: list[str] = []
    tags: list[int] = []
    mol_id: list[int] = []
    m_no = 0
    for ci in sorted(range(len(mols)), key=lambda x: -len(mols[x][0])):
        mol, n = mols[ci]
        base = mol.get_positions() - mol.get_positions().mean(axis=0)
        for _ in range(n):
            local_parts.append(base @ _random_rotation(rng).T)
            centers.append(sites[m_no % len(sites)] + rng.uniform(-0.3 * step, 0.3 * step, size=3))
            symbols += mol.get_chemical_symbols(); tags += [ci] * len(mol); mol_id += [m_no] * len(mol)
            m_no += 1
            report(m_no, total, L("分子を置いています", "placing molecules"))
    local = np.vstack(local_parts)
    cen = np.array(centers)
    ids = np.array(mol_id)
    pos = local + cen[ids]
    pos, ok = _relax(pos, cen, local, ids, box, spec.min_distance, spec.max_tries // 4)
    if not ok:
        near = min_inter_distance(pos, ids, box, spec.min_distance * 2)
        raise MixtureError(L(f"箱 {box:.1f} Å に {spec.summary()} を最短距離 {spec.min_distance} Å では詰められません (いちばん近い組で {near:.2f} Å)。"
                             "箱を大きくする、個数を減らす、最短距離を短くする、のどれかにしてください",
                             f"cannot pack {spec.summary()} into a {box:.1f} Å box with a minimum distance of {spec.min_distance} Å "
                             f"(closest pair {near:.2f} Å). Enlarge the box, reduce the counts, or lower the minimum distance"))
    atoms = Atoms(symbols=symbols, positions=pos, cell=[box, box, box], pbc=True, tags=tags)
    atoms.set_array("molecule", ids)
    return atoms



MAX_ATOMS = 50_000


def check_atom_limit(n: int, what: str = "") -> None:
    if n > MAX_ATOMS:
        raise MixtureError(L(f"{what}原子が {n:,} 個になり、上限 {MAX_ATOMS:,} 個を超えます (作る前に止めました)。個数や大きさを減らしてください",
                             f"{what}this would make {n:,} atoms, above the limit of {MAX_ATOMS:,} (stopped before building). Reduce the counts or the size"))


def salt_count(concentration_mol_l: float, volume_a3: float) -> tuple[int, float]:
    if concentration_mol_l < 0 or volume_a3 <= 0:
        raise MixtureError(L("濃度は 0 以上、体積は 0 より大きい値にしてください", "concentration must be >= 0 and volume > 0"))
    n = int(round(concentration_mol_l * AVOGADRO * volume_a3 * A3_TO_L))
    return n, n / (AVOGADRO * volume_a3 * A3_TO_L)


def perpendicular_widths(cell: np.ndarray) -> np.ndarray:
    vol = abs(float(np.linalg.det(cell)))
    return np.array([vol / np.linalg.norm(np.cross(cell[(k + 1) % 3], cell[(k + 2) % 3])) for k in range(3)])


def _directions(n: int = 1024) -> np.ndarray:
    k = np.arange(n) + 0.5
    z = k / n
    r = np.sqrt(np.clip(1.0 - z * z, 0.0, None))
    phi = np.pi * (1.0 + 5.0 ** 0.5) * k
    return np.stack([r * np.cos(phi), r * np.sin(phi), z], axis=1)


def molecule_widths(mol: Atoms) -> tuple[float, float]:
    p = np.asarray(mol.get_positions(), dtype=float)
    if len(p) < 2:
        return 0.0, 0.0
    thin = float(np.ptp(p @ _directions().T, axis=0).min())
    return thin, float(np.linalg.norm(p[:, None, :] - p[None, :, :], axis=-1).max())


def packing_need(mols: list[tuple[Atoms, int]], min_distance: float) -> tuple[float, str, float, str]:
    thin, long = (0.0, ""), (0.0, "")
    for m, n in mols:
        if n <= 0:
            continue
        t, l = molecule_widths(m)
        f = m.get_chemical_formula()
        thin = max(thin, (t, f)); long = max(long, (l, f))
    return thin[0] + min_distance, thin[1], long[0] + min_distance, long[1]


def vdw_volume(mols: list[tuple[Atoms, int]]) -> float:
    from ase.data import covalent_radii, vdw_radii

    v = 0.0
    for m, n in mols:
        if n <= 0:
            continue
        r = np.array([vdw_radii[z] if np.isfinite(vdw_radii[z]) else covalent_radii[z] + 0.8 for z in m.get_atomic_numbers()])
        v += n * float((4.0 / 3.0 * np.pi * r ** 3).sum())
    return v


def _grid_counts(widths: np.ndarray, total: int) -> np.ndarray:
    s = (float(np.prod(widths)) / total) ** (1.0 / 3.0)
    n = np.maximum(1, np.round(widths / s)).astype(int)
    while int(np.prod(n)) < total:
        k = int(np.argmax(widths / n)); n[k] += 1
    return n


class Packed:

    def __init__(self, symbols, positions, tags, mol_id, ok, near, band, cell=None, vdw=0.0):
        self.symbols, self.positions, self.tags, self.mol_id = symbols, positions, tags, mol_id
        self.ok, self.near, self.band = ok, near, band
        self.cell = None if cell is None else np.asarray(cell, dtype=float)
        self.vdw = float(vdw)

    def to_atoms(self, cell) -> Atoms:
        a = Atoms(symbols=self.symbols, positions=self.positions, cell=cell, pbc=True, tags=self.tags)
        a.set_array("molecule", np.asarray(self.mol_id, dtype=int))
        return a


def _room_note(packed: "Packed") -> str:
    if packed.cell is None:
        return ""
    cell = packed.cell
    area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
    if packed.band:
        thick = float(packed.band[1] - packed.band[0])
        vol = area * thick
        parts = [L(f"断面 {area:,.1f} Å²、帯の厚み {thick:,.1f} Å (帯の体積 {vol:,.0f} Å³)",
                   f"cross-section {area:,.1f} Å², band thickness {thick:,.1f} Å (band volume {vol:,.0f} Å³)")]
    else:
        vol = abs(float(np.linalg.det(cell)))
        parts = [L(f"断面 {area:,.1f} Å²、セルの体積 {vol:,.0f} Å³", f"cross-section {area:,.1f} Å², cell volume {vol:,.0f} Å³")]
    if packed.vdw > 0 and vol > 0:
        parts.append(L(f"成分のファンデルワールス体積の合計 {packed.vdw:,.0f} Å³ (体積の {100 * packed.vdw / vol:.0f} %)",
                       f"total van der Waals volume of the components {packed.vdw:,.0f} Å³ ({100 * packed.vdw / vol:.0f} % of the volume)"))
    return L("。いまの広さ: " + "、".join(parts), ". Current room: " + ", ".join(parts))


def pack_failure_message(summary: str, packed: "Packed", min_distance: float) -> str:
    where = L(f"z の帯 {packed.band[0]:.1f}〜{packed.band[1]:.1f} Å", f"the z band {packed.band[0]:.1f}-{packed.band[1]:.1f} Å") if packed.band \
        else L("セル", "the cell")
    band_note = L("スラブとの隙間 (gap) は帯の厚みを変えません。", "The slab gap does not change the band thickness. ") if packed.band else ""
    return L(f"{where}に {summary} を最短距離 {min_distance} Å では詰められません (いちばん近い組で {packed.near:.2f} Å){_room_note(packed)}。"
             f"{band_note}密度を下げる (帯なら厚くなる)、個数を減らす、最短距離を短くする、のどれかにしてください",
             f"cannot pack {summary} into {where} with a minimum distance of {min_distance} Å (closest pair {packed.near:.2f} Å){_room_note(packed)}. "
             f"{band_note}Lower the density (a thicker band), reduce the counts, or lower the minimum distance")


def pack_molecules(mols: list[tuple[Atoms, int]], cell, rng: np.random.Generator, *, min_distance: float = 2.0, steps: int = 500,
                   z_band: tuple[float, float] | None = None, obstacles: np.ndarray | None = None) -> Packed:
    from ase.neighborlist import neighbor_list

    cell = np.asarray(cell, dtype=float)
    if abs(np.linalg.det(cell)) < 1e-6:
        raise MixtureError(L("セルの体積が 0 です (格子ベクトルが同じ面に乗っています)", "the cell has zero volume"))
    inv = np.linalg.inv(cell)
    obs = np.zeros((0, 3)) if obstacles is None else np.asarray(obstacles, dtype=float).reshape(-1, 3)
    total = sum(n for _, n in mols)
    if z_band is not None:
        lo, hi = float(z_band[0]), float(z_band[1])
        if hi <= lo:
            raise MixtureError(L(f"z の帯の厚みが 0 以下です ({lo:.2f}〜{hi:.2f} Å)", f"the z band is empty ({lo:.2f}-{hi:.2f} Å)"))
        if abs(cell[0, 2]) > 1e-6 or abs(cell[1, 2]) > 1e-6 or abs(cell[2, 0]) > 1e-6 or abs(cell[2, 1]) > 1e-6:
            raise MixtureError(L("z の帯を使うには、a, b が xy 面内にあり c が z 方向 (面に垂直) のセルが要ります",
                                 "a z band needs a cell with a, b in the xy plane and c along z (normal to the plane)"))
        area = float(np.linalg.norm(np.cross(cell[0], cell[1])))
        widths = np.array([area / np.linalg.norm(cell[1]), area / np.linalg.norm(cell[0]), hi - lo])
    else:
        widths = perpendicular_widths(cell)
    n = _grid_counts(widths, total)
    sites = np.array([((a + 0.5) / n[0], (b + 0.5) / n[1], (c + 0.5) / n[2])
                      for a in range(n[0]) for b in range(n[1]) for c in range(n[2])])
    rng.shuffle(sites)

    local_parts, centers, symbols, tags, mol_id = [], [], [], [], []
    m_no = 0
    for ci in sorted(range(len(mols)), key=lambda x: -len(mols[x][0])):
        mol, cnt = mols[ci]
        base = mol.get_positions() - mol.get_positions().mean(axis=0)
        for _ in range(cnt):
            loc = base @ _random_rotation(rng).T
            f = sites[m_no % len(sites)] + rng.uniform(-0.3, 0.3, size=3) / n
            if z_band is not None:
                for _t in range(50):
                    ext = float(np.ptp(loc[:, 2]))
                    if ext <= hi - lo:
                        break
                    loc = base @ _random_rotation(rng).T
                else:
                    raise MixtureError(L(f"z の帯の厚み {hi - lo:.2f} Å が {mol.get_chemical_formula()} の大きさより薄く、入りません。帯を厚くしてください",
                                         f"the z band ({hi - lo:.2f} Å) is thinner than {mol.get_chemical_formula()}. Make the band thicker"))
                c = f[0] * cell[0] + f[1] * cell[1]
                c[2] = lo - loc[:, 2].min() + float(np.clip(f[2], 0.0, 1.0)) * (hi - lo - ext)
            else:
                c = f @ cell
            local_parts.append(loc); centers.append(c)
            symbols += mol.get_chemical_symbols(); tags += [ci] * len(mol); mol_id += [m_no] * len(mol)
            m_no += 1
            report(m_no, total, L("分子を置いています", "placing molecules"))
    local = np.vstack(local_parts); cen = np.array(centers); ids = np.array(mol_id)
    npos, nmol = len(local), m_no
    ids_all = np.concatenate([ids, np.full(len(obs), nmol)])

    def wrap(cen):
        fr = cen @ inv
        if z_band is not None:
            fr[:, :2] -= np.floor(fr[:, :2])
        else:
            fr -= np.floor(fr)
        return fr @ cell

    def pairs(pos, cutoff):
        probe = Atoms(numbers=np.ones(npos + len(obs), dtype=int), positions=np.vstack([pos, obs]), cell=cell, pbc=True)
        i, j, S = neighbor_list("ijS", probe, cutoff=cutoff)
        sel = (ids_all[i] != ids_all[j]) & (i < j)
        return i[sel], j[sel], S[sel] @ cell

    def clamp(cen):
        if z_band is None:
            return cen
        z = local[:, 2] + cen[ids, 2]
        zmin = np.full(nmol, np.inf); zmax = np.full(nmol, -np.inf)
        np.minimum.at(zmin, ids, z); np.maximum.at(zmax, ids, z)
        cen[:, 2] += np.where(zmin < lo, lo - zmin, 0.0) + np.where(zmax > hi, hi - zmax, 0.0)
        return cen

    r0 = float(min_distance)
    damp, target, skin = 0.7, r0 * 1.03, 0.8
    cen = clamp(wrap(cen))
    pos = local + cen[ids]
    ii = None
    for it in range(max(50, steps)):
        report(it + 1, max(50, steps), L("分子の重なりをほどいています", "relaxing overlaps"))
        if ii is None or it % 8 == 0:
            cen = wrap(cen); pos = local + cen[ids]
            ii, jj, sv = pairs(pos, target + skin)
        if len(ii) == 0:
            break
        allp = np.vstack([pos, obs])
        dv = allp[jj] - allp[ii] + sv
        d = np.sqrt((dv * dv).sum(axis=1))
        close = (d < target) & (d > 1e-9)
        if not close.any():
            cen = wrap(cen); pos = local + cen[ids]
            ii, jj, sv = pairs(pos, target)
            if len(ii) == 0:
                break
            continue
        i, j, dv, d = ii[close], jj[close], dv[close], d[close]
        push = ((target - d) * damp / 2.0)[:, None] * (dv / d[:, None])
        mob = j < npos
        shift = np.zeros_like(cen)
        np.add.at(shift, ids[i], -push * np.where(mob, 1.0, 2.0)[:, None])
        np.add.at(shift, ids[j[mob]], push[mob])
        big = np.linalg.norm(shift, axis=1)
        over = big > r0
        if over.any():
            shift[over] *= (r0 / big[over])[:, None]
        cen = clamp(cen + shift)
        pos = local + cen[ids]
    cen = wrap(cen); pos = local + cen[ids]
    ii, jj, sv = pairs(pos, 2 * r0)
    allp = np.vstack([pos, obs])
    d = np.sqrt(((allp[jj] - allp[ii] + sv) ** 2).sum(axis=1))
    near = float(d.min()) if len(d) else float("inf")
    return Packed(symbols, pos, tags, ids, near >= r0 - 1e-6, near, (lo, hi) if z_band is not None else None,
                  cell=cell, vdw=vdw_volume(mols))


def parse_mixture_text(text: str) -> list[Component]:
    out: list[Component] = []
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        head, _, rest = line.partition("*")
        ref = head.strip(); kind = "preset"
        for k in ("smiles", "file"):
            if ref.lower().startswith(k + ":"):
                kind, ref = k, ref[len(k) + 1:].strip()
        if not ref:
            raise MixtureError(L(f"成分の指定が空です: {raw!r}", f"empty component: {raw!r}"))
        toks = rest.split()
        try:
            count = int(toks[0]) if toks else 1
        except ValueError as ex:
            raise MixtureError(L(f"個数として読めません: {raw!r}", f"not a count: {raw!r}")) from ex
        charge, label = 0, ""
        for t in toks[1:]:
            k, _, v = t.partition("=")
            if k == "charge":
                try:
                    charge = int(v)
                except ValueError as ex:
                    raise MixtureError(L(f"電荷として読めません: {raw!r}", f"not a charge: {raw!r}")) from ex
            elif k == "label":
                label = v
        out.append(Component(kind=kind, ref=ref, count=count, charge=charge, label=label))
    return out


def mixture_text(components: list[Component]) -> str:
    lines = []
    for c in components:
        ref = c.ref if c.kind == "preset" else f"{c.kind}:{c.ref}"
        extra = (f" charge={c.charge:+d}" if c.charge else "") + (f" label={c.label}" if c.label else "")
        lines.append(f"{ref} * {c.count}{extra}")
    return "\n".join(lines)
