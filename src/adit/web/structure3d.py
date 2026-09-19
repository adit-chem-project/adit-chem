
from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np
from ase import Atoms
from ase.data import covalent_radii
from ase.data.colors import jmol_colors

MAX_ATOMS = 2000
MAX_BOND_ATOMS = 600
BOND_FACTOR = 1.2

CELL_EDGES = ((0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7))


@dataclass
class Scene:
    positions: list[list[float]] = field(default_factory=list)
    radii: list[float] = field(default_factory=list)
    colors: list[str] = field(default_factory=list)              # "#rrggbb"
    symbols: list[str] = field(default_factory=list)
    bonds: list[list[int]] = field(default_factory=list)
    cell_lines: list[list[list[float]]] = field(default_factory=list)
    cell: list[list[float]] | None = None                          # rows = lattice vectors, for minimum-image measuring
    n_atoms: int = 0
    truncated_bonds: bool = False
    scale: float = 1.0

    def to_json(self) -> str:
        return json.dumps(self.__dict__, separators=(",", ":")).replace("<", "\\u003c")


def _hex(rgb) -> str:
    return "#" + "".join(f"{int(round(float(c) * 255)):02x}" for c in rgb[:3])


def scene_from_atoms(atoms: Atoms) -> Scene:
    scene = Scene()
    if atoms is None or len(atoms) == 0:
        return scene
    pos = np.asarray(atoms.get_positions(), dtype=float)
    numbers = list(atoms.get_atomic_numbers())
    center = pos.mean(axis=0)
    pos = pos - center
    scene.n_atoms = len(atoms)
    scene.positions = [[round(float(x), 4) for x in p] for p in pos]
    scene.radii = [round(float(covalent_radii[z]), 3) for z in numbers]
    scene.colors = [_hex(jmol_colors[z]) for z in numbers]
    scene.symbols = list(atoms.get_chemical_symbols())
    if any(atoms.pbc) and abs(np.linalg.det(np.asarray(atoms.cell, dtype=float))) > 0:
        cell = np.asarray(atoms.cell, dtype=float)
        corners = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)]) @ cell - center
        scene.cell_lines = [[[round(float(x), 4) for x in corners[a]],
                             [round(float(x), 4) for x in corners[b]]] for a, b in CELL_EDGES]
        scene.cell = [[round(float(x), 6) for x in row] for row in cell]
    if len(atoms) <= MAX_BOND_ATOMS:
        scene.bonds = _bonds(pos, numbers)
    else:
        scene.truncated_bonds = True
    far = float(np.linalg.norm(pos, axis=1).max()) if len(pos) else 1.0
    if scene.cell_lines:
        far = max(far, float(np.linalg.norm(np.array(scene.cell_lines).reshape(-1, 3), axis=1).max()))
    scene.scale = round(max(far, 1.0), 4)
    return scene


def _bonds(pos: np.ndarray, numbers) -> list[list[int]]:
    radii = np.array([covalent_radii[z] for z in numbers])
    out: list[list[int]] = []
    for i in range(len(pos) - 1):
        d = np.linalg.norm(pos[i + 1:] - pos[i], axis=1)
        limit = (radii[i] + radii[i + 1:]) * BOND_FACTOR
        for j in np.nonzero((d > 0.1) & (d <= limit))[0]:
            out.append([i, int(i + 1 + j)])
    return out


def summary_line(scene: Scene, periodic: bool) -> str:
    from adit.lang import L

    kind = L("周期系 (セルの辺も描いています)", "periodic (the cell edges are drawn)") if periodic else \
        L("非周期 (分子)", "non-periodic (a molecule)")
    bonds = (L("原子が多いので結合は描いていません", "too many atoms, so no bonds are drawn")
             if scene.truncated_bonds else L(f"結合 {len(scene.bonds)} 本", f"{len(scene.bonds)} bonds"))
    return L(f"原子 {scene.n_atoms} 個、{bonds}、{kind}。色は元素の色 (ASE の jmol)、大きさは共有結合半径です",
             f"{scene.n_atoms} atoms, {bonds}, {kind}. Colors are the ASE jmol element colors and the sizes are covalent radii")
