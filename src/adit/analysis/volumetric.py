
from __future__ import annotations

from adit.errors import AditValueError
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ase import Atoms

from adit.lang import L

VASP_NAMES = ("LOCPOT", "CHGCAR", "PARCHG", "CHG")
CUBE_SUFFIXES = (".cube", ".cub")


class VolumetricError(AditValueError):
    pass


@dataclass
class Grid:
    values: np.ndarray
    atoms: Atoms
    kind: str            # "potential" / "density" / "unknown"
    unit: str
    source: Path
    origin: np.ndarray = field(default_factory=lambda: np.zeros(3))   # Å; grid point (0,0,0)
    periodic: bool = False                                              # True: the grid covers one cell without the end point

    @property
    def cell(self):
        return self.atoms.cell


def find_files(run_dir: Path | str) -> list[Path]:
    d = Path(run_dir).expanduser()
    if not d.is_dir():
        return []
    out = [p for p in sorted(d.iterdir()) if p.is_file()
           and (p.suffix.lower() in CUBE_SUFFIXES or p.name.startswith(VASP_NAMES))]
    return out


CUBE_UNITS = {"ev": ("eV", 1.0), "ry": ("eV", 13.605693122994), "hartree": ("eV", 27.211386245988),
              "ha": ("eV", 27.211386245988)}


def read_grid(path: Path | str, cube_unit: str = "") -> Grid:
    p = Path(path).expanduser()
    if not p.is_file():
        raise VolumetricError(L(f"ファイルがありません: {p}", f"file not found: {p}"))
    if p.suffix.lower() in CUBE_SUFFIXES:
        from ase.io import read

        try:
            with open(p, encoding="utf-8", errors="replace") as f:
                head = f.readline() + f.readline()
            dct = read(str(p), format="cube", read_data=True, full_output=True)
            values, atoms, origin = dct["data"], dct["atoms"], dct["origin"]
        except Exception as ex:
            raise VolumetricError(L(f"{p.name} を cube として読めません: {ex}",
                                    f"cannot read {p.name} as a cube file: {ex}")) from ex
        values = np.asarray(values, dtype=float)
        # pp.x (Quantum ESPRESSO, PP/src/cube.f90) writes the whole FFT grid, nr points at a/nr, without the end point
        periodic = _cube_is_periodic(head)
        extra = {"origin": np.asarray(origin, dtype=float), "periodic": periodic}
        if cube_unit:
            key = cube_unit.strip().lower()
            if key not in CUBE_UNITS:
                raise VolumetricError(L(f"cube の単位は {' / '.join(sorted(CUBE_UNITS))} から選んでください: {cube_unit!r}",
                                        f"the cube unit must be one of {' / '.join(sorted(CUBE_UNITS))}: {cube_unit!r}"))
            unit, factor = CUBE_UNITS[key]
            return Grid(values=values * factor, atoms=atoms, kind="potential", unit=unit, source=p, **extra)
        return Grid(values=values, atoms=atoms, kind="unknown", unit="", source=p, **extra)
    if p.name.startswith(VASP_NAMES):
        from ase.calculators.vasp import VaspChargeDensity

        try:
            vcd = VaspChargeDensity(str(p))
        except Exception as ex:
            raise VolumetricError(L(f"{p.name} を VASP の体積データとして読めません: {ex}",
                                    f"cannot read {p.name} as VASP volumetric data: {ex}")) from ex
        if not vcd.chg:
            raise VolumetricError(L(f"{p.name} に値がありません", f"{p.name} contains no values"))
        atoms = vcd.atoms[-1]
        values = np.asarray(vcd.chg[-1], dtype=float)
        if p.name.startswith("LOCPOT"):
            values = values * atoms.get_volume()
            return Grid(values=values, atoms=atoms, kind="potential", unit="eV", source=p, periodic=True)
        return Grid(values=values, atoms=atoms, kind="density", unit="e/Å³", source=p, periodic=True)
    raise VolumetricError(L(f"{p.name} は ADIT が読める体積データではありません (cube / LOCPOT / CHGCAR)",
                            f"{p.name} is not volumetric data ADIT can read (cube / LOCPOT / CHGCAR)"))


PERIODIC_CUBE_MARKS = ("pwscf", "quantum espresso")


def _cube_is_periodic(comment_lines: str) -> bool:
    text = comment_lines.lower()
    return any(mark in text for mark in PERIODIC_CUBE_MARKS)


def plane_average(grid: Grid, axis: int = 2) -> tuple[np.ndarray, np.ndarray]:
    if axis not in (0, 1, 2):
        raise VolumetricError(L("軸は 0 (a) / 1 (b) / 2 (c) のどれかです", "the axis must be 0 (a), 1 (b) or 2 (c)"))
    other = tuple(i for i in (0, 1, 2) if i != axis)
    mean = grid.values.mean(axis=other)
    length = float(np.linalg.norm(np.asarray(grid.cell)[axis]))
    positions = np.linspace(0.0, length, len(mean), endpoint=False)
    return positions, mean


def work_function(grid: Grid, fermi_ev: float | None, axis: int = 2) -> dict:
    positions, mean = plane_average(grid, axis)
    index = int(np.argmax(mean))
    vacuum = float(mean[index])
    out = {"axis": "abc"[axis], "vacuum_level_ev": vacuum, "vacuum_at_ang": float(positions[index]),
           "fermi_ev": fermi_ev, "work_function_ev": None, "source": grid.source.name,
           "flatness_ev": float(np.max(mean[max(0, index - 2):index + 3]) - np.min(mean[max(0, index - 2):index + 3])),
           "note": L("真空準位は、面平均が最大になる位置の値です。真空の領域が十分に広く平らかどうかは、図を見て利用者が確かめてください。",
                     "The vacuum level is the maximum of the plane-averaged profile; check on the figure that the vacuum region is wide and flat.")}
    if grid.unit != "eV":
        out["note"] = L(f"このファイルの値は eV ではありません ({grid.unit or '単位は不明'})。仕事関数は求めません。",
                        f"the values in this file are not in eV ({grid.unit or 'unknown unit'}); no work function is computed.")
        return out
    if fermi_ev is not None:
        out["work_function_ev"] = vacuum - float(fermi_ev)
    else:
        out["note"] += L(" フェルミ準位が読めないので、仕事関数は求めていません。",
                         " The Fermi level could not be read, so no work function is given.")
    return out


def write_csv(path: Path, positions: np.ndarray, mean: np.ndarray, unit: str) -> Path:
    import csv

    with open(path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["position_ang", f"mean_value{f' [{unit}]' if unit else ''}"])
        writer.writerows([[f"{x:.6f}", f"{y:.10g}"] for x, y in zip(positions, mean)])
    return path


def difference(grids: list, coefficients: list[float]) -> np.ndarray:
    if len(grids) != len(coefficients) or not grids:
        raise VolumetricError(L("格子データと係数の数が合いません", "the number of grids and coefficients differ"))
    first = grids[0]
    total = np.zeros_like(np.asarray(first.values, dtype=float))
    for grid, c in zip(grids, coefficients):
        values = np.asarray(grid.values, dtype=float)
        if values.shape != total.shape:
            raise VolumetricError(L(f"格子の分割が違います: {values.shape} と {total.shape} "
                                    "(同じ FFT 格子で計算し直してください)",
                                    f"the grids have different shapes: {values.shape} vs {total.shape} "
                                    "(recompute them on the same FFT grid)"))
        if not np.allclose(np.asarray(grid.atoms.cell), np.asarray(first.atoms.cell), atol=1e-6):
            raise VolumetricError(L("セルが違う格子データは引き算できません (同じセルで計算し直してください)",
                                    "grids with different cells cannot be combined (recompute them in the same cell)"))
        total += c * values
    return total


def macroscopic_average(positions: np.ndarray, values: np.ndarray, period_ang: float) -> np.ndarray:
    positions = np.asarray(positions, dtype=float)
    values = np.asarray(values, dtype=float)
    if positions.size != values.size or positions.size < 2:
        raise VolumetricError(L("位置と値の数が合いません", "positions and values must have the same length"))
    if period_ang <= 0:
        raise VolumetricError(L("周期 [Å] を正の値で指定してください (何を 1 周期と見るかは利用者が決めます)",
                                "give a positive period in Å (what counts as one period is your choice)"))
    step = float(positions[1] - positions[0])
    width = max(1, int(round(period_ang / step)))
    kernel = np.ones(width) / width
    extended = np.concatenate([values, values, values])
    smoothed = np.convolve(extended, kernel, mode="same")[len(values):2 * len(values)]
    return smoothed


class DensityGrid:

    def __init__(self, shape=(48, 48, 48)):
        self.shape = tuple(int(x) for x in shape)
        if min(self.shape) < 2:
            raise VolumetricError(L("分割数は各方向 2 以上です", "each dimension of the grid must be at least 2"))
        self.counts = np.zeros(self.shape, dtype=float)
        self.n_frames = 0
        self.cell = None
        self.reference = None

    def add(self, atoms, indices=None) -> None:
        cell = np.asarray(atoms.cell, dtype=float)
        if abs(np.linalg.det(cell)) <= 0:
            raise VolumetricError(L("3 次元の密度には周期セルが要ります", "a periodic cell is required for a 3D density"))
        if self.cell is None:
            self.cell, self.reference = cell, atoms.copy()
        pos = atoms.get_positions()
        if indices is not None:
            pos = pos[np.asarray(indices, dtype=int)]
        frac = np.linalg.solve(cell.T, pos.T).T % 1.0
        idx = np.minimum((frac * np.array(self.shape)).astype(int), np.array(self.shape) - 1)
        np.add.at(self.counts, (idx[:, 0], idx[:, 1], idx[:, 2]), 1.0)
        self.n_frames += 1

    def result(self) -> dict:
        if self.n_frames == 0 or self.cell is None:
            raise VolumetricError(L("フレームが 1 つも入っていません", "no frames were added"))
        voxel = abs(float(np.linalg.det(self.cell))) / float(np.prod(self.shape))
        values = self.counts / (self.n_frames * voxel)
        return {"values": values, "shape": self.shape, "n_frames": self.n_frames, "voxel_A3": voxel,
                "note": L("軌跡を通した数密度 [Å⁻³] です。分割数は利用者の指定です。cube に書き出すと、解析タブの「等値面」の枠や "
                          "VMD・VESTA・OVITO で等値面を見られます。",
                          "number density in Å⁻³ accumulated over the trajectory; the grid is yours. Written as a cube file, "
                          "the Isosurface box of the analysis tab, VMD, VESTA or OVITO can draw its isosurface.")}

    def write_cube(self, path, comment: str = "") -> "Path":
        from ase.io.cube import write_cube

        got = self.result()
        path = Path(path)
        with open(path, "w", encoding="utf-8") as f:
            write_cube(f, self.reference, data=got["values"], comment=comment or "adit 3D density [1/Ang^3]")
        return path
