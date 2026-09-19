"""Vibrational normal modes (frequency + Cartesian displacement per atom) read from finished runs."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import chemical_symbols

from adit.lang import L

AU_TO_CM1 = 219474.6313705
AMU_TO_AU = 1822.888486209
SUPPORTED = ("dftbplus", "xtb", "gaussian", "orca")


@dataclass
class Modes:
    atoms: Atoms
    frequencies_cm1: list[float]
    vectors: np.ndarray                   # (n_modes, n_atoms, 3): Cartesian direction, scaled so the largest atom moves 1
    source: str
    notes: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.frequencies_cm1)

    def displaced(self, mode: int, amplitude_ang: float, phase: float) -> np.ndarray:
        """Positions with mode displaced by amplitude*sin(phase); amplitude is the largest atomic displacement in Å."""
        return self.atoms.get_positions() + amplitude_ang * np.sin(phase) * self.vectors[mode]

    def label(self, mode: int) -> str:
        return f"{mode + 1}: {self.frequencies_cm1[mode]:.1f} cm⁻¹"


def _normalize(vectors: np.ndarray) -> np.ndarray:
    out = np.array(vectors, dtype=float)
    for k in range(len(out)):
        big = float(np.linalg.norm(out[k], axis=1).max()) if out[k].size else 0.0
        if big > 1e-12:
            out[k] /= big
    return out


def modes_from_hessian(h_flat: np.ndarray, atoms: Atoms) -> tuple[list[float], np.ndarray]:
    """Frequencies [cm-1] and Cartesian mode vectors from a DFTB+ hessian.out (Hartree/Bohr^2, 3N x 3N)."""
    n = len(atoms) * 3
    h_flat = np.asarray(h_flat, dtype=float)
    if h_flat.size < n * n or not np.all(np.isfinite(h_flat[: n * n])):
        raise ValueError(L(f"hessian.out が途中で切れているか壊れています ({n}×{n} = {n * n} 個の数が要るところ、{h_flat.size} 個)",
                           f"hessian.out is truncated or broken ({n}x{n} = {n * n} numbers needed, {h_flat.size} found)"))
    h = h_flat[: n * n].reshape(n, n)
    h = 0.5 * (h + h.T)
    masses_au = np.repeat(atoms.get_masses() * AMU_TO_AU, 3)
    hw = h / np.sqrt(np.outer(masses_au, masses_au))
    lam, w = np.linalg.eigh(hw)
    freqs = [float(np.sign(x) * np.sqrt(abs(x)) * AU_TO_CM1) for x in lam]
    cart = (w / np.sqrt(masses_au)[:, None]).T.reshape(n, len(atoms), 3)
    return freqs, _normalize(cart)


_G98_FREQ = re.compile(r"^\s*Frequencies\s+--\s+(.*)$")
_G98_ROW = re.compile(r"^\s*(\d+)\s+(\d+)\s+((?:-?\d+\.\d+\s*)+)$")
_G98_ORIENT = re.compile(r"(Standard|Input) orientation:")


def _g98_orientation(lines: list[str], before: int) -> Atoms | None:
    start = None
    for k in range(before - 1, -1, -1):
        if _G98_ORIENT.search(lines[k]):
            start = k
            break
    if start is None:
        return None
    numbers, pos = [], []
    k = start + 5
    while k < len(lines):
        w = lines[k].split()
        if len(w) >= 6 and w[0].isdigit() and w[1].isdigit():
            numbers.append(int(w[1])); pos.append([float(x) for x in w[-3:]])
        elif pos:
            break
        k += 1
    return Atoms(numbers=numbers, positions=pos) if pos else None


def read_g98_modes(path: Path, atoms: Atoms | None = None) -> Modes | None:
    """Gaussian-style frequency output (Gaussian .log, xtb g98.out): 'Frequencies --' blocks with 3 modes per block."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    freqs: list[float] = []
    blocks: list[list[list[float]]] = []      # per block: rows of 3*k floats
    first_freq = None
    k = 0
    while k < len(lines):
        m = _G98_FREQ.match(lines[k])
        if not m:
            k += 1; continue
        if first_freq is None:
            first_freq = k
        vals = [float(x) for x in m.group(1).split()]
        freqs += vals
        rows: list[list[float]] = []
        k += 1
        while k < len(lines) and not re.match(r"^\s*Atom\s+AN\b", lines[k]):
            if _G98_FREQ.match(lines[k]) or not lines[k].strip():
                break
            k += 1
        if k < len(lines) and re.match(r"^\s*Atom\s+AN\b", lines[k]):
            k += 1
            while k < len(lines):
                r = _G98_ROW.match(lines[k])
                if not r:
                    break
                nums = [float(x) for x in r.group(3).split()]
                if len(nums) != 3 * len(vals):
                    break
                rows.append(nums); k += 1
        blocks.append(rows)
    if not freqs or first_freq is None:
        return None
    if atoms is None:
        atoms = _g98_orientation(lines, first_freq)
    if atoms is None:
        return None
    nat = len(atoms)
    vectors = []
    mode = 0
    for rows in blocks:
        per = len(rows[0]) // 3 if rows else 0
        for j in range(per):
            vec = np.zeros((nat, 3))
            for a, r in enumerate(rows[:nat]):
                vec[a] = r[3 * j: 3 * j + 3]
            vectors.append(vec); mode += 1
    if len(vectors) != len(freqs) or any(len(rows) != nat for rows in blocks):
        return None
    return Modes(atoms, freqs, _normalize(np.array(vectors)), Path(path).name)


_ORCA_FREQ = re.compile(r"^\s*(\d+):\s+(-?[\d.]+)\s+cm\*\*-1")


def _orca_coordinates(lines: list[str]) -> Atoms | None:
    start = None
    for k, l in enumerate(lines):
        if l.startswith("CARTESIAN COORDINATES (ANGSTROEM)"):
            start = k
    if start is None:
        return None
    symbols, pos = [], []
    for l in lines[start + 2:]:
        w = l.split()
        if len(w) == 4 and w[0] in chemical_symbols:
            symbols.append(w[0]); pos.append([float(x) for x in w[1:]])
        elif pos:
            break
    return Atoms(symbols=symbols, positions=pos) if pos else None


def read_orca_modes(path: Path, atoms: Atoms | None = None) -> Modes | None:
    """ORCA output: 'VIBRATIONAL FREQUENCIES' (3N values) and the 'NORMAL MODES' matrix (columns = modes)."""
    lines = Path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    freqs: list[float] = []
    matrix: dict[int, dict[int, float]] = {}
    in_freq = in_modes = False
    cols: list[int] = []
    for l in lines:
        if "VIBRATIONAL FREQUENCIES" in l:
            in_freq, freqs = True, []; continue
        if "NORMAL MODES" in l:
            in_freq, in_modes, matrix = False, True, {}; continue
        if in_freq:
            m = _ORCA_FREQ.match(l)
            if m:
                freqs.append(float(m.group(2)))
            continue
        if in_modes:
            w = l.split()
            if not w:
                continue
            if all(x.isdigit() for x in w) and len(w) <= 6 and not (len(w) > 1 and "." in l):
                cols = [int(x) for x in w]; continue
            if w[0].isdigit() and len(w) == len(cols) + 1:
                try:
                    vals = [float(x) for x in w[1:]]
                except ValueError:
                    continue
                row = int(w[0])
                for c, v in zip(cols, vals):
                    matrix.setdefault(c, {})[row] = v
                continue
            if matrix and l.strip() and not w[0].isdigit():
                in_modes = False
    if not freqs or not matrix:
        return None
    if atoms is None:
        atoms = _orca_coordinates(lines)
    if atoms is None:
        return None
    nat = len(atoms)
    n = 3 * nat
    if len(freqs) != n or len(matrix) != n:
        return None
    vectors, kept = [], []
    for j in range(n):
        col = np.array([matrix[j].get(i, 0.0) for i in range(n)]).reshape(nat, 3)
        if abs(freqs[j]) < 1e-6 and not np.any(col):
            continue                         # ORCA lists translations/rotations as zero modes
        vectors.append(col); kept.append(freqs[j])
    if not vectors:
        return None
    return Modes(atoms, kept, _normalize(np.array(vectors)), Path(path).name)


def unsupported_reason(code: str | None) -> str:
    names = "DFTB+ (hessian.out), xtb (g98.out), Gaussian, ORCA"
    return L(f"振動モードの表示は {names} に対応しています (このディレクトリのコード: {code or '不明'})",
             f"mode display supports {names} (this directory: {code or 'unknown'})")


def read_modes(run_dir: Path | str, code: str | None = None) -> Modes | None:
    """Modes of a run directory, or None when the code is not supported or nothing was found."""
    from adit.analysis.readers import detect_code, load_run

    d = Path(run_dir)
    code = code or detect_code(d)
    if code not in SUPPORTED:
        return None
    try:
        data = load_run(d, code)
        final = data.final
    except Exception:
        final = None
    if code == "dftbplus":
        hess = d / "hessian.out"
        if not hess.is_file() or final is None:
            return None
        freqs, vec = modes_from_hessian(np.array(hess.read_text(encoding="utf-8", errors="replace").split(), dtype=float), final)
        m = Modes(final.copy(), freqs, vec, "hessian.out")
        m.notes.append(L("hessian.out を質量重み付きで対角化した固有ベクトル (並進・回転は落としていない)",
                         "eigenvectors of the mass-weighted hessian.out (translations and rotations are not projected out)"))
        return m
    if code == "xtb":
        p = d / "g98.out"
        return read_g98_modes(p, final) if p.is_file() else None
    if code == "gaussian":
        from adit.analysis.readers_qc import find_output
        p = find_output(d)
        return read_g98_modes(p) if p is not None else None
    if code == "orca":
        p = d / "output.log"
        return read_orca_modes(p, final) if p.is_file() else None
    return None


__all__ = ["Modes", "modes_from_hessian", "read_g98_modes", "read_orca_modes", "read_modes", "unsupported_reason", "SUPPORTED"]
