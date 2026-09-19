
from __future__ import annotations

import time
from dataclasses import dataclass

import numpy as np
from ase import Atoms

from adit.builder.bases import build_base, limit
from adit.builder.model import Recipe, RecipeError, op_label
from adit.builder.ops import OPS, make_molecules_whole, prepared, split_molecules
from adit.lang import L

_HINTS = {
    "base": ("土台の構造の時点で原子が重なっています。ファイルなら中身を、分子ならセルの大きさを確かめてください",
             "atoms already overlap in the base structure; check the file, or the cell size for a molecule"),
    "supercell": ("元のセルより大きい分子 (セルの外にはみ出した分子) を繰り返すと重なります。分子なら先に box で十分な箱に入れてください",
                  "repeating a molecule larger than its cell makes overlaps; put a molecule in a large enough box first"),
    "slab": ("終端 (termination) か層の数を変えてください", "change the termination or the number of layers"),
    "vacuum": ("真空を厚くしてください", "make the vacuum thicker"),
    "box": ("padding を大きくしてください", "increase the padding"),
    "adsorb": ("高さ (height) を上げるか、置き場所を変えてください。スラブの上に真空が無いと、周期の像のスラブの下面に当たります",
               "raise the height or move the site; without vacuum above the slab it hits the periodic image"),
    "remove": ("", ""), "substitute": ("", ""),
    "solvent_layer": ("隙間 (gap) か最短距離を大きくしてください", "increase the gap or the minimum distance"),
    "solvate": ("最短距離を大きくしてください", "increase the minimum distance"), "fix": ("", ""),
}


@dataclass
class StepLog:
    index: int
    op: str
    n_atoms: int
    min_distance: float
    seconds: float
    note: str = ""

    def line(self) -> str:
        d = "> 3 Å" if not np.isfinite(self.min_distance) else f"{self.min_distance:.2f} Å"
        return L(f"{self.index}. {op_label(self.op)}: {self.n_atoms} 原子、分子どうしの最短距離 {d}、{self.seconds:.2f} 秒{('、' + self.note) if self.note else ''}",
                 f"{self.index}. {op_label(self.op)}: {self.n_atoms} atoms, shortest distance between molecules {d}, {self.seconds:.2f} s{(', ' + self.note) if self.note else ''}")


def _pairs(atoms: Atoms, cutoff: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if any(atoms.pbc) and abs(np.linalg.det(np.asarray(atoms.cell))) > 1e-6:
        from ase.neighborlist import neighbor_list
        i, j, d = neighbor_list("ijd", atoms, cutoff)
        keep = i < j
        return i[keep], j[keep], d[keep]
    from scipy.spatial import cKDTree
    pos = atoms.get_positions()
    pr = cKDTree(pos).query_pairs(cutoff, output_type="ndarray").reshape(-1, 2)
    return pr[:, 0], pr[:, 1], np.linalg.norm(pos[pr[:, 0]] - pos[pr[:, 1]], axis=1)


def _nearest(atoms: Atoms, cutoff: float = 3.0) -> tuple[float, int, int, float]:
    i, j, d = _pairs(atoms, cutoff)
    if not len(d):
        return float("inf"), -1, -1, float("inf")
    k = int(np.argmin(d))
    mol = atoms.arrays.get("molecule")
    inter = d if mol is None else d[(mol[i] != mol[j]) | (mol[i] < 0)]
    return float(d[k]), int(i[k]), int(j[k]), float(inter.min()) if len(inter) else float("inf")


def check_step(atoms: Atoms, index: int, op: str) -> tuple[Atoms, StepLog]:
    from adit.validate import MIN_DISTANCE

    limit(len(atoms))
    if len(atoms) == 0:
        raise RecipeError(L(f"手順 {index} ({op_label(op)}) のあとで原子が 0 個になりました", f"no atoms left after step {index} ({op_label(op)})"))
    note = ""
    split = split_molecules(atoms)
    if split:
        atoms = make_molecules_whole(atoms)
        note = L(f"境界で切れた分子 {len(split)} 個をつなぎました", f"rejoined {len(split)} molecules split by the boundary")
    d, i, j, d_inter = _nearest(atoms)
    if d < MIN_DISTANCE:
        s = atoms.get_chemical_symbols()
        ja, en = _HINTS.get(op, ("", ""))
        raise RecipeError(L(f"手順 {index} ({op_label(op)}) のあとで、原子 {i + 1} ({s[i]}) と {j + 1} ({s[j]}) が {d:.2f} Å まで近づいています "
                            f"({MIN_DISTANCE} Å 未満は重なりとみなします)。{ja}",
                            f"after step {index} ({op_label(op)}), atoms {i + 1} ({s[i]}) and {j + 1} ({s[j]}) are {d:.2f} Å apart "
                            f"(below {MIN_DISTANCE} Å counts as overlapping). {en}"))
    return atoms, StepLog(index, op, len(atoms), d_inter, 0.0, note)


def order_notes(steps: list) -> dict[int, str]:
    notes: dict[int, list[str]] = {}
    slabs = [k for k, s in enumerate(steps, start=1) if s.op == "slab"]
    for k in slabs:
        before = [j for j, s in enumerate(steps[:k - 1], start=1) if s.op == "supercell"]
        if before:
            notes.setdefault(k, []).append(
                L(f"注意: 手順 {'、'.join(str(j) for j in before)} の超格子より後ろで切っています。ミラー指数はいまのセルの格子ベクトルに対する指数なので、"
                  "超格子のあとでは同じ指数でも別の面になります (面で切ってから超格子にすると、面はそのままで断面が広がります)",
                  f"note: this cut comes after the supercell in step {', '.join(str(j) for j in before)}. Miller indices refer to the "
                  "current cell vectors, so the same indices mean a different plane after a supercell (cutting first and then "
                  "repeating keeps the plane and widens the cross-section)"))
    for k in slabs[1:]:
        notes.setdefault(k, []).append(
            L(f"注意: 面で切る手順が {len(slabs)} 個あります (手順 {'、'.join(str(j) for j in slabs)})。このミラー指数は、手順 {slabs[0]} で切ったスラブのセルに対する指数です",
              f"note: there are {len(slabs)} slab steps (steps {', '.join(str(j) for j in slabs)}); these Miller indices refer to the "
              f"cell of the slab made in step {slabs[0]}"))
    return {k: L("、".join(v), ", ".join(v)) for k, v in notes.items()}


def build_recipe(recipe: Recipe) -> tuple[Atoms, list[StepLog]]:
    logs: list[StepLog] = []
    t0 = time.perf_counter()
    atoms = prepared(build_base(recipe.base))
    atoms, log = check_step(atoms, 0, "base")
    log.seconds = time.perf_counter() - t0
    logs.append(log)
    notes = order_notes(recipe.steps)
    for k, step in enumerate(recipe.steps, start=1):
        t0 = time.perf_counter()
        try:
            atoms = OPS[step.op](atoms, step)
        except RecipeError as ex:
            raise RecipeError(L(f"手順 {k} ({op_label(step.op)}): {ex}", f"step {k} ({op_label(step.op)}): {ex}")) from ex
        made = atoms.info.pop("recipe_note", "")
        atoms, log = check_step(atoms, k, step.op)
        log.note = L("、", ", ").join([x for x in (log.note, made, notes.get(k, "")) if x])
        log.seconds = time.perf_counter() - t0
        logs.append(log)
    return atoms, logs


def recipe_structure(recipe: Recipe | str, *, charge: int | None = None, multiplicity: int = 1):
    from adit.spec import AtomsData, Structure

    if isinstance(recipe, str):
        recipe = Recipe.from_ref(recipe)
    atoms, logs = build_recipe(recipe)
    fixed = [int(i) for i in np.flatnonzero(atoms.get_array("fixed"))] if "fixed" in atoms.arrays else []
    st = Structure(source="recipe", source_ref=recipe.to_ref(), atoms=AtomsData.from_ase(atoms),
                   charge=recipe.total_charge() if charge is None else charge, multiplicity=multiplicity, fixed_atoms=fixed)
    return st, logs
