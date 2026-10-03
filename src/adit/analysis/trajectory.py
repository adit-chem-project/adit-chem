
from __future__ import annotations

from adit.errors import AditValueError
import math
import re
from array import array
from collections.abc import Iterator, Sequence
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import chemical_symbols
from ase.units import create_units

from adit.lang import L

_QE_UNITS = create_units("2006")

MEMORY_BUDGET_MB = 1024
BYTES_PER_ATOM_FRAME = 24


class TrajectoryTooLarge(AditValueError):

    def __init__(self, message: str, suggested_stride: int):
        super().__init__(message)
        self.suggested_stride = suggested_stride


def _symbol(label: str) -> str:
    s = re.sub(r"[^A-Za-z]", " ", label).split()
    w = s[0] if s else label
    for n in (2, 1):
        cand = w[:n].capitalize()
        if cand in chemical_symbols:
            return cand
    raise ValueError(L(f"元素記号として読めません: {label!r}", f"not an element symbol: {label!r}"))


class Trajectory(Sequence):
    def __init__(self, path: Path | str, kind: str, *, cell=None, pbc=None, fmt: str | None = None,
                 comment_regex: re.Pattern | None = None, cells=None, _sel: tuple = (), _shared: dict | None = None):
        self.path = Path(path)
        self.kind = kind
        self.cell = cell
        self.cells = cells  # per-frame cells (xyz files without lattice, variable-cell MD)
        self.pbc = pbc
        self.fmt = fmt
        self.comment_regex = comment_regex
        self._sel = _sel
        self._shared = _shared if _shared is not None else {"total": None, "offsets": None, "comment_values": None}

    def __repr__(self) -> str:
        return f"Trajectory({self.path.name!r}, kind={self.kind!r}, sel={self._sel!r})"

    def total_frames(self) -> int:
        if self._shared["total"] is None:
            self._shared["total"] = self._count()
        return self._shared["total"]

    def _apply_sel(self, n: int) -> range:
        r = range(n)
        for s in self._sel:
            r = r[s]
        return r

    def _range(self) -> range:
        return self._apply_sel(self.total_frames())

    def __len__(self) -> int:
        return len(self._range())

    def __getitem__(self, i):
        if isinstance(i, slice):
            if i.step is not None and i.step <= 0:
                raise ValueError("Trajectory: step must be positive")
            return Trajectory(self.path, self.kind, cell=self.cell, pbc=self.pbc, fmt=self.fmt, comment_regex=self.comment_regex,
                              cells=self.cells, _sel=self._sel + (i,), _shared=self._shared)
        return self._read_at(self._range()[i])

    def _cell_at(self, k: int):
        if self.cells is not None and k < len(self.cells):
            return self.cells[k]
        return self.cell

    def __iter__(self) -> Iterator[Atoms]:
        r = self._range()
        if len(r) == 0:
            return
        start, step, last = r.start, r.step, r[-1]
        for k, atoms in self._frames(lambda k: k >= start and (k - start) % step == 0, last):
            if atoms is not None:
                yield atoms

    @property
    def comment_values(self) -> list[float]:
        if self._shared["comment_values"] is None:
            self._shared["total"] = self._count()
        return self._shared["comment_values"] or []

    def file_size(self) -> int:
        return self.path.stat().st_size

    def estimate_total_frames(self) -> int:
        if self._shared["total"] is not None:
            return self._shared["total"]
        per = self._first_frame_bytes()
        if not per:
            return self.total_frames()
        return max(1, math.ceil(self.file_size() / per))

    def estimate_len(self) -> int:
        return len(self._apply_sel(self.estimate_total_frames()))

    def first_natoms(self) -> int:
        for _, a in self._frames(lambda k: k == 0, 0):
            if a is not None:
                return len(a)
        return 0

    def _open(self):
        return open(self.path, "rb")

    def _count(self) -> int:
        if self.kind == "xyz":
            offsets = array("q")
            values: list[float] = []
            with self._open() as f:
                while True:
                    pos = f.tell()
                    head = f.readline()
                    if not head:
                        break
                    try:
                        n = int(head.split()[0])
                    except (ValueError, IndexError):
                        continue
                    comment = f.readline()
                    got = 0
                    for _ in range(n):
                        if not f.readline():
                            break
                        got += 1
                    if got < n:
                        break
                    offsets.append(pos)
                    if self.comment_regex is not None:
                        m = self.comment_regex.search(comment.decode("utf-8", "replace"))
                        if m:
                            values.append(float(m.group(1)))
            self._shared["offsets"] = offsets
            self._shared["comment_values"] = values
            return len(offsets)
        if self.kind == "xdatcar":
            n = 0
            with self._open() as f:
                for line in f:
                    if b"configuration=" in line:
                        n += 1
            return n
        return sum(1 for _ in self._frames(lambda k: False, None))

    def _first_frame_bytes(self) -> int | None:
        try:
            with self._open() as f:
                if self.kind == "xyz":
                    while True:
                        head = f.readline()
                        if not head:
                            return None
                        try:
                            n = int(head.split()[0]); break
                        except (ValueError, IndexError):
                            continue
                    start = f.tell() - len(head)
                    for _ in range(n + 1):
                        f.readline()
                    return max(1, f.tell() - start)
                if self.kind == "xdatcar":
                    first = None
                    for line in iter(f.readline, b""):
                        if b"configuration=" in line:
                            if first is not None:
                                return max(1, f.tell() - len(line) - first)
                            first = f.tell() - len(line)
                    return None
                marker = {"qe": b"ATOMIC_POSITIONS", "ase": b"<calculation>", "lammpsdump": b"ITEM: TIMESTEP"}.get(self.kind)
                if marker is None:
                    return None
                seen = []
                read = 0
                for line in iter(f.readline, b""):
                    read += len(line)
                    if marker in line:
                        seen.append(f.tell())
                        if len(seen) == 2:
                            return max(1, seen[1] - seen[0])
                    if read > 256 * 1024 * 1024:
                        return None
                return None
        except OSError:
            return None

    def _read_at(self, idx: int) -> Atoms:
        if self.kind == "xyz":
            self.total_frames()
            offsets = self._shared["offsets"]
            with self._open() as f:
                f.seek(offsets[idx])
                a = _parse_xyz_block(f, self.path, self._cell_at(idx), self.pbc)
            if a is None:
                raise IndexError(idx)
            return a
        for k, a in self._frames(lambda k: k == idx, idx):
            if a is not None:
                return a
        raise IndexError(idx)

    def _frames(self, want, last: int | None) -> Iterator[tuple[int, Atoms | None]]:
        gen = {"xyz": self._frames_xyz, "xdatcar": self._frames_xdatcar, "qe": self._frames_qe, "ase": self._frames_ase,
               "lammpsdump": self._frames_lammps}[self.kind]
        yield from gen(want, last)

    def _frames_lammps(self, want, last):
        types = self.fmt.split() if self.fmt else []
        k = -1
        with open(self.path, encoding="utf-8", errors="replace") as fd:
            while last is None or k < last:
                line = fd.readline()
                if not line:
                    return
                if not line.startswith("ITEM: TIMESTEP"):
                    continue
                fd.readline()
                fd.readline()  # ITEM: NUMBER OF ATOMS
                try:
                    n = int(fd.readline())
                except ValueError:
                    return
                head = fd.readline().split()  # ITEM: BOX BOUNDS [xy xz yz] pp pp pp
                b = [fd.readline().split() for _ in range(3)]
                cols = fd.readline().split()[2:]
                rows = [fd.readline() for _ in range(n)]
                if len(rows) < n or (n and len(rows[-1].split()) < len(cols)):
                    return
                k += 1
                if not want(k):
                    yield k, None
                    continue
                lo = np.array([float(x[0]) for x in b]); hi = np.array([float(x[1]) for x in b])
                xy = xz = yz = 0.0
                if "xy" in head:
                    xy, xz, yz = float(b[0][2]), float(b[1][2]), float(b[2][2])
                    lo[0] -= min(0.0, xy, xz, xy + xz); hi[0] -= max(0.0, xy, xz, xy + xz)
                    lo[1] -= min(0.0, yz); hi[1] -= max(0.0, yz)
                cell = np.array([[hi[0] - lo[0], 0, 0], [xy, hi[1] - lo[1], 0], [xz, yz, hi[2] - lo[2]]])
                flags = head[-3:]
                pbc = [f == "pp" for f in flags] if all(len(f) == 2 for f in flags) else True
                data = [r.split() for r in rows]
                if "id" in cols:
                    ic = cols.index("id")
                    data.sort(key=lambda w: int(w[ic]))
                if "element" in cols:
                    ie = cols.index("element")
                    syms = [w[ie] for w in data]
                elif "type" in cols and types:
                    it = cols.index("type")
                    syms = [types[int(w[it]) - 1] for w in data]
                else:
                    raise ValueError(L(f"{self.path.name} に元素 (element) の列が無く、型番号と元素の対応も分かりません",
                                       f"{self.path.name} has no element column and the type-to-element mapping is unknown"))
                for names, kind in ((("x", "y", "z"), "cart"), (("xu", "yu", "zu"), "cart"), (("xs", "ys", "zs"), "frac")):
                    if all(c in cols for c in names):
                        p = np.array([[float(w[cols.index(c)]) for c in names] for w in data])
                        break
                else:
                    raise ValueError(L(f"{self.path.name} に座標の列 (x y z / xu yu zu / xs ys zs) がありません",
                                       f"{self.path.name} has no coordinate columns (x y z / xu yu zu / xs ys zs)"))
                a = Atoms(symbols=syms, cell=cell, pbc=pbc)
                if kind == "frac":
                    a.set_scaled_positions(p)
                else:
                    a.set_positions(p - lo)
                yield k, a

    def _frames_xyz(self, want, last):
        k = -1
        with self._open() as f:
            while last is None or k < last:
                head = f.readline()
                if not head:
                    return
                try:
                    n = int(head.split()[0])
                except (ValueError, IndexError):
                    continue
                k += 1
                if want(k):
                    f.seek(f.tell() - len(head))
                    a = _parse_xyz_block(f, self.path, self._cell_at(k), self.pbc)
                    if a is None:
                        return
                    yield k, a
                else:
                    got = 0
                    for _ in range(n + 1):
                        if f.readline():
                            got += 1
                    if got < n + 1:
                        return
                    yield k, None

    def _frames_xdatcar(self, want, last):
        cell = np.eye(3)
        formula = ""
        k = -1
        with open(self.path, encoding="utf-8", errors="replace") as fd:
            while last is None or k < last:
                comment_line = fd.readline()
                if not comment_line:
                    return
                if "Direct configuration=" not in comment_line:
                    try:
                        scale = float(fd.readline())
                    except ValueError:
                        return
                    xx = [float(x) for x in fd.readline().split()]
                    yy = [float(y) for y in fd.readline().split()]
                    zz = [float(z) for z in fd.readline().split()]
                    cell = np.array([xx, yy, zz]) * scale
                    symbols = fd.readline().split()
                    numbers = [int(n) for n in fd.readline().split()]
                    total = sum(numbers)
                    formula = "".join(f"{s}{numbers[i]}" for i, s in enumerate(symbols))
                    fd.readline()
                k += 1
                rows = [fd.readline() for _ in range(total)]
                if len(rows[-1].split()) < 3 if rows else True:
                    return
                if want(k):
                    a = Atoms(formula, cell=cell, pbc=True)
                    a.set_scaled_positions(np.array([r.split()[:3] for r in rows], dtype=float))
                    yield k, a
                else:
                    yield k, None

    def _frames_qe(self, want, last):
        bohr = _QE_UNITS["Bohr"]
        alat = None
        nat = None
        cell = None
        pending = None  # (symbols, positions Å, cell)
        emitted = True
        k = -1
        with open(self.path, encoding="utf-8", errors="replace") as fd:
            for line in fd:
                if "lattice parameter (alat)" in line:
                    alat = float(line.split("=")[1].split()[0]) * bohr
                elif "number of atoms/cell" in line:
                    nat = int(line.split("=")[1].split()[0])
                elif "crystal axes: (cart. coord. in units of alat)" in line:
                    rows = [next(fd) for _ in range(3)]
                    cell = np.array([_floats_after_eq(r) for r in rows]) * alat
                elif "site n." in line and "positions (alat units)" in line:
                    syms, pos = [], []
                    for _ in range(nat):
                        r = next(fd)
                        syms.append(_symbol(r.split()[1]))
                        pos.append(_floats_after_eq(r))
                    pending = (syms, np.array(pos) * alat, cell.copy()); emitted = False
                elif line.startswith("CELL_PARAMETERS"):
                    rows = np.array([[float(x) for x in next(fd).split()[:3]] for _ in range(3)])
                    m = re.search(r"alat\s*=\s*([-\d.Ee+]+)", line)
                    if m:
                        cell = rows * float(m.group(1)) * bohr
                    elif "bohr" in line.lower():
                        cell = rows * bohr
                    else:
                        cell = rows
                elif line.startswith("ATOMIC_POSITIONS"):
                    unit = line.split("(")[1].split(")")[0].strip().lower() if "(" in line else "alat"
                    syms, pos = [], []
                    for _ in range(nat):
                        w = next(fd).split()
                        syms.append(_symbol(w[0])); pos.append([float(x) for x in w[1:4]])
                    p = np.array(pos)
                    if unit == "crystal":
                        p = p @ cell
                    elif unit == "bohr":
                        p = p * bohr
                    elif unit == "alat":
                        p = p * alat
                    pending = (syms, p, cell.copy()); emitted = False
                elif line.startswith("!") and "total energy" in line and pending is not None and not emitted:
                    emitted = True
                    k += 1
                    if want(k):
                        syms, p, c = pending
                        yield k, Atoms(symbols=syms, positions=p, cell=c, pbc=True)
                    else:
                        yield k, None
                    if last is not None and k >= last:
                        return

    def _frames_ase(self, want, last):
        from ase.io import iread
        for k, a in enumerate(iread(str(self.path), index=":", format=self.fmt)):
            if last is not None and k > last:
                return
            yield k, (a if want(k) else None)


def _floats_after_eq(line: str) -> list[float]:
    v = [float(x) for x in re.findall(r"-?\d+\.\d*(?:[EeDd][-+]?\d+)?", line.split("=", 1)[1].replace("D", "E").replace("d", "e"))]
    if len(v) < 3:
        raise ValueError(L(f"pw.x の出力の行を読めません: {line.strip()[:60]!r}", f"cannot read a line of the pw.x output: {line.strip()[:60]!r}"))
    return v[:3]


def _parse_xyz_block(f, path: Path, cell, pbc) -> Atoms | None:
    head = f.readline()
    while head and not head.split():
        head = f.readline()
    if not head:
        return None
    n = int(head.split()[0])
    f.readline()
    syms, pos = [], []
    for _ in range(n):
        raw = f.readline()
        if not raw:
            return None
        w = raw.split()
        try:
            syms.append(w[0].decode()); pos.append([float(w[1]), float(w[2]), float(w[3])])
        except (IndexError, ValueError) as ex:
            l = raw.decode("utf-8", "replace")
            raise ValueError(L(f"{path.name} のフレームを読めません (列が足りないか、数でない値): {l.strip()[:60]!r}",
                               f"cannot read a frame of {path.name} (missing columns or non-numbers): {l.strip()[:60]!r}")) from ex
    a = Atoms(symbols=syms, positions=pos)
    if cell is not None:
        a.set_cell(cell); a.pbc = pbc if pbc is not None else True
    return a


def check_budget(frames, n_atoms: int, budget_mb: float = MEMORY_BUDGET_MB, what: str = "MSD") -> dict:
    if isinstance(frames, Trajectory):
        est_frames = frames.estimate_len()
        size = frames.file_size()
    else:
        est_frames, size = len(frames), None
    need = est_frames * n_atoms * BYTES_PER_ATOM_FRAME
    limit = budget_mb * 1024 * 1024
    info = {"estimated_frames": int(est_frames), "n_atoms": int(n_atoms), "estimated_mb": need / 2**20, "budget_mb": float(budget_mb),
            "file_bytes": size}
    if need > limit:
        stride = math.ceil(need / limit)
        cur = frames._sel[-1].step if isinstance(frames, Trajectory) and frames._sel and frames._sel[-1].step else 1
        raise TrajectoryTooLarge(L(
            f"軌跡が大きすぎるので {what} を計算せずに止めました。見積もり: {est_frames:,} フレーム × {n_atoms:,} 原子 = 約 {need / 2**20:,.0f} MB "
            f"(上限 {budget_mb:,.0f} MB)。間引いてください: 例 --stride {stride * cur} (いまの {stride} 倍の間隔)。"
            f"上限そのものは --memory-mb で変えられます",
            f"the trajectory is too large, so {what} was not computed. Estimate: {est_frames:,} frames x {n_atoms:,} atoms = about {need / 2**20:,.0f} MB "
            f"(limit {budget_mb:,.0f} MB). Thin it out, e.g. --stride {stride * cur} ({stride} times the current spacing). "
            f"The limit itself can be changed with --memory-mb"), stride * cur)
    return info
