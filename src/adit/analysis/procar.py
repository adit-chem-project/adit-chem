"""Read VASP PROCAR files (projected band structure)."""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from adit.lang import L

_HEAD = re.compile(r"#\s*of\s*k-points:\s*(\d+)\s*#\s*of\s*bands:\s*(\d+)\s*#\s*of\s*ions:\s*(\d+)")
_KPT = re.compile(r"k-point\s+(\d+)\s*:(.*?)weight\s*=\s*([-\d.E+]+)")
_FLOAT = re.compile(r"[-+]?\d*\.\d+(?:[eE][-+]?\d+)?")
_BAND = re.compile(r"band\s+(\d+)\s*#\s*energy\s*([-\d.E+]+)\s*#\s*occ\.\s*([-\d.E+]+)")


@dataclass
class Procar:
    kpoints_frac: np.ndarray            # (nk, 3)
    weights: np.ndarray                 # (nk,)
    energies_ev: np.ndarray             # (nspin, nk, nbands)
    occupations: np.ndarray             # (nspin, nk, nbands)
    projections: np.ndarray             # (nspin, nk, nbands, nions, norbitals)
    orbitals: list[str] = field(default_factory=list)
    noncollinear: bool = False

    @property
    def n_spins(self) -> int:
        return self.energies_ev.shape[0]

    def by_ion(self, ion_index: int, spin: int = 0) -> np.ndarray:
        if not 1 <= ion_index <= self.projections.shape[3]:
            raise ValueError(L(f"原子 {ion_index} は範囲の外です (1〜{self.projections.shape[3]})",
                               f"atom {ion_index} is out of range (1-{self.projections.shape[3]})"))
        return self.projections[spin, :, :, ion_index - 1, :].sum(axis=-1)

    def by_orbital(self, name: str, spin: int = 0) -> np.ndarray:
        if name not in self.orbitals:
            raise ValueError(L(f"軌道 {name} がありません (あるのは {', '.join(self.orbitals)})",
                               f"no orbital {name} (available: {', '.join(self.orbitals)})"))
        return self.projections[spin, :, :, :, self.orbitals.index(name)].sum(axis=-1)

    def total(self, spin: int = 0) -> np.ndarray:
        return self.projections[spin].sum(axis=(-1, -2))


def _open(path: Path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return p.open(encoding="utf-8", errors="replace")


def read_procar(path: Path) -> Procar:
    """Read a PROCAR file, plain or gzipped."""
    with _open(path) as stream:
        lines = stream.read().splitlines()
    head = None
    for line in lines[:5]:
        m = _HEAD.search(line)
        if m:
            head = m
            break
    if head is None:
        raise ValueError(L(f"{Path(path).name} に見出し (# of k-points / bands / ions) がありません",
                           f"{Path(path).name} has no header (# of k-points / bands / ions)"))
    nk, nb, nion = (int(x) for x in head.groups())
    orbitals: list[str] = []
    kpoints, weights = [], []
    energies, occupations, projections = [], [], []
    spin_blocks: list[tuple[list, list, list]] = []
    cur_e, cur_o, cur_p = [], [], []
    noncollinear = False
    i = 0
    kpoint_seen = 0
    while i < len(lines):
        line = lines[i]
        mk = _KPT.search(line)
        if mk:
            if kpoint_seen == nk:
                spin_blocks.append((cur_e, cur_o, cur_p))
                cur_e, cur_o, cur_p = [], [], []
                kpoint_seen = 0
            kpoint_seen += 1
            if len(kpoints) < nk:
                coords = _FLOAT.findall(mk.group(2))
                if len(coords) != 3:
                    raise ValueError(L(f"{Path(path).name}: k 点の座標を 3 つ読めません: {mk.group(2)!r}",
                                       f"{Path(path).name}: cannot read 3 k-point coordinates: {mk.group(2)!r}"))
                kpoints.append([float(x) for x in coords])
                weights.append(float(mk.group(3)))
            k_e, k_o, k_p = [], [], []
            i += 1
            for _ in range(nb):
                while i < len(lines) and not _BAND.search(lines[i]):
                    i += 1
                mb = _BAND.search(lines[i])
                k_e.append(float(mb.group(2)))
                k_o.append(float(mb.group(3)))
                while i < len(lines) and not lines[i].lstrip().startswith("ion "):
                    i += 1
                if not orbitals:
                    orbitals = lines[i].split()[1:-1]
                nor = len(orbitals)
                rows = []
                i += 1
                for _ in range(nion):
                    rows.append([float(x) for x in lines[i].split()[1:1 + nor]])
                    i += 1
                k_p.append(rows)
                while i < len(lines) and not (_BAND.search(lines[i]) or _KPT.search(lines[i])):
                    if lines[i].lstrip().startswith("ion ") and "tot" not in lines[i]:
                        pass
                    i += 1
            cur_e.append(k_e)
            cur_o.append(k_o)
            cur_p.append(k_p)
            continue
        i += 1
    spin_blocks.append((cur_e, cur_o, cur_p))
    for e, o, p in spin_blocks:
        energies.append(e)
        occupations.append(o)
        projections.append(p)
    energies = np.array(energies, dtype=float)
    occupations = np.array(occupations, dtype=float)
    projections = np.array(projections, dtype=float)
    if energies.shape[1:] != (nk, nb):
        raise ValueError(L(f"{Path(path).name}: 読めた k 点とバンドの数 {energies.shape[1:]} が "
                           f"見出しの {(nk, nb)} と違います",
                           f"{Path(path).name}: the k-point and band counts read {energies.shape[1:]} "
                           f"differ from {(nk, nb)} in the header"))
    return Procar(np.array(kpoints), np.array(weights), energies, occupations, projections,
                  orbitals, noncollinear)


def summary_lines(pc: Procar) -> list[str]:
    nspin, nk, nb = pc.energies_ev.shape
    lines = [L(f"PROCAR: {nk} k 点 × {nb} バンド × 原子 {pc.projections.shape[3]} 個"
               f"{'、スピン 2 つ' if nspin == 2 else ''}",
               f"PROCAR: {nk} k-points x {nb} bands x {pc.projections.shape[3]} atoms"
               f"{', 2 spins' if nspin == 2 else ''}"),
             L(f"  軌道: {', '.join(pc.orbitals)}", f"  orbitals: {', '.join(pc.orbitals)}"),
             L(f"  全射影の和の平均 {pc.total().mean():.3f} (1 に足りない分は原子軌道で表しきれない成分です)",
               f"  mean of the sum of all projections {pc.total().mean():.3f} (the shortfall from 1 is the part atomic orbitals cannot represent)")]
    return lines
