"""Read Gaussian, GAMESS and Q-Chem outputs."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import chemical_symbols

HARTREE_EV = 27.211386245988
BOHR_ANG = 0.529177210903
GAMESS_IR_TO_KM_PER_MOL = 42.2561

OUTPUT_NAMES = ("output.log", "output.out", "out.log")


def find_output(run_dir: Path) -> Path | None:
    for name in OUTPUT_NAMES:
        p = Path(run_dir) / name
        if p.is_file():
            return p
    hits = sorted(Path(run_dir).glob("*.out")) + sorted(Path(run_dir).glob("*.log"))
    return hits[0] if hits else None


def _lines(path: Path):
    path = Path(path)
    if path.suffix == ".gz":
        import gzip

        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as stream:
            yield from stream
    else:
        with path.open(encoding="utf-8", errors="replace") as stream:
            yield from stream


def _floats(text: str) -> list[float]:
    return [float(x) for x in re.findall(r"-?\d+\.\d+", text)]


def _gamess_frequencies(text: str) -> list[float]:
    # GAMESS prints an imaginary frequency as its magnitude followed by "I"
    return [-abs(float(v)) if flag else float(v) for v, flag in re.findall(r"(-?\d+\.\d+)\s*(I\b)?", text)]


class QcOutput:

    def __init__(self, code: str, path: Path):
        self.code = code
        self.path = path
        self.energies_ev: list[float] = []
        self.frequencies_cm1: list[float] = []
        self.ir_intensities: list[float] = []
        self.raman_activities: list[float] = []
        self.frames: list[Atoms] = []
        self.charges: list[dict] = []


def _add_charges(out: QcOutput, values: list[float], definition: str, source: str) -> None:
    if values:
        out.charges.append({"definition": definition, "values": values, "source": source,
                            "unit": "e", "note": ""})


def read_gaussian(path: Path) -> QcOutput:
    """Read a Gaussian .log or .out file. Verified against Gaussian 16 output."""
    out = QcOutput("gaussian", Path(path))
    symbols: list[str] = []
    coords: list[list[float]] = []
    reading_geom = mulliken = False
    skip = 0
    charges: list[float] = []
    definition, source = "Mulliken", "Gaussian の Mulliken charges"
    for line in _lines(path):
        if line.startswith(" SCF Done:"):
            m = re.search(r"=\s*(-?\d+\.\d+)", line)
            if m:
                out.energies_ev.append(float(m.group(1)) * HARTREE_EV)
            continue
        if "Frequencies --" in line:
            out.frequencies_cm1 += _floats(line.split("--", 1)[1])
            continue
        if "IR Inten" in line:
            out.ir_intensities += _floats(line.split("--", 1)[1])
            continue
        if "Raman Activ" in line:
            out.raman_activities += _floats(line.split("--", 1)[1])
            continue
        if "Standard orientation:" in line or "Input orientation:" in line:
            reading_geom, skip, symbols, coords = True, 4, [], []
            continue
        if reading_geom:
            if skip > 0:
                skip -= 1
                continue
            if set(line.strip()) <= {"-"} and line.strip():
                reading_geom = False
                if symbols:
                    out.frames.append(Atoms(symbols=symbols, positions=np.array(coords)))
                continue
            parts = line.split()
            if len(parts) >= 6:
                symbols.append(chemical_symbols[int(parts[1])])
                coords.append([float(x) for x in parts[-3:]])
            continue
        if "Mulliken charges" in line and line.rstrip().endswith(":") and not line.lstrip().startswith("Sum"):
            if charges:
                _add_charges(out, charges, definition, source)
                charges = []
            summed = "hydrogens summed" in line
            definition = "Mulliken (水素を重原子にまとめた)" if summed else "Mulliken"
            source = "Gaussian の " + line.strip().rstrip(":")
            mulliken = True
            continue
        if mulliken:
            parts = line.split()
            if len(parts) == 3 and parts[0].isdigit():
                charges.append(float(parts[2]))
            elif charges:
                mulliken = False
                _add_charges(out, charges, definition, source)
                charges = []
            continue
    if charges:
        _add_charges(out, charges, definition, source)
    return _keep_last_per_definition(out)


def read_gamess(path: Path) -> QcOutput:
    out = QcOutput("gamess", Path(path))
    freqs: list[float] = []
    irs: list[float] = []
    ramans: list[float] = []
    symbols: list[str] = []
    coords: list[list[float]] = []
    reading_geom = in_bohr = False
    skip = 0
    charges: list[float] = []
    mulliken = False
    for line in _lines(path):
        if "FINAL" in line and "ENERGY IS" in line:
            m = re.search(r"IS\s+(-?\d+\.\d+)", line)
            if m:
                out.energies_ev.append(float(m.group(1)) * HARTREE_EV)
            continue
        if line.lstrip().startswith("FREQUENCY:"):
            freqs += _gamess_frequencies(line.split(":", 1)[1])
            continue
        if line.lstrip().startswith("IR INTENSITY:"):
            irs += _floats(line.split(":", 1)[1])
            continue
        if line.lstrip().startswith("RAMAN ACTIVITY:"):
            ramans += _floats(line.split(":", 1)[1])
            continue
        if "COORDINATES OF ALL ATOMS ARE (ANGS)" in line:
            reading_geom, skip, symbols, coords, in_bohr = True, 2, [], [], False
            continue
        if "COORDINATES (BOHR)" in line:
            reading_geom, skip, symbols, coords, in_bohr = True, 1, [], [], True
            continue
        if reading_geom:
            if skip > 0:
                skip -= 1
                continue
            parts = line.split()
            if len(parts) >= 5 and parts[0].isalpha():
                symbols.append(parts[0].capitalize())
                coords.append([float(x) for x in parts[-3:]])
                continue
            reading_geom = False
            if symbols:
                pos = np.array(coords) * (BOHR_ANG if in_bohr else 1.0)
                out.frames.append(Atoms(symbols=symbols, positions=pos))
            continue
        if "TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS" in line:
            mulliken, charges = True, []
            continue
        if mulliken:
            parts = line.split()
            if len(parts) == 6 and parts[0].isdigit() and parts[1].isalpha():
                charges.append(float(parts[3]))
            elif charges and not line.strip():
                mulliken = False
            continue
    n_zero = 6 if len(freqs) > 6 else 0
    out.frequencies_cm1 = freqs[n_zero:]
    out.ir_intensities = [v * GAMESS_IR_TO_KM_PER_MOL for v in irs[n_zero:]]
    out.raman_activities = ramans[n_zero:] if ramans else []
    if charges:
        _add_charges(out, charges, "Mulliken", "GAMESS の TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS")
    return _keep_last_per_definition(out)


def read_qchem(path: Path) -> QcOutput:
    out = QcOutput("qchem", Path(path))
    symbols: list[str] = []
    coords: list[list[float]] = []
    reading_geom = False
    skip = 0
    charges: list[float] = []
    mulliken = False
    for line in _lines(path):
        if "Total energy in the final basis set" in line:
            m = re.search(r"=\s*(-?\d+\.\d+)", line)
            if m:
                out.energies_ev.append(float(m.group(1)) * HARTREE_EV)
            continue
        if line.lstrip().startswith("Frequency:"):
            out.frequencies_cm1 += _floats(line.split(":", 1)[1])
            continue
        if line.lstrip().startswith("IR Intens:"):
            out.ir_intensities += _floats(line.split(":", 1)[1])
            continue
        if line.lstrip().startswith("Raman Intens:"):
            out.raman_activities += _floats(line.split(":", 1)[1])
            continue
        if "Standard Nuclear Orientation" in line:
            reading_geom, skip, symbols, coords = True, 2, [], []
            continue
        if reading_geom:
            if skip > 0:
                skip -= 1
                continue
            parts = line.split()
            if len(parts) == 5 and parts[0].isdigit():
                symbols.append(parts[1])
                coords.append([float(x) for x in parts[2:5]])
                continue
            reading_geom = False
            if symbols:
                out.frames.append(Atoms(symbols=symbols, positions=np.array(coords)))
            continue
        if "Ground-State Mulliken Net Atomic Charges" in line:
            mulliken, charges, skip = True, [], 3
            continue
        if mulliken:
            if skip > 0:
                skip -= 1
                continue
            parts = line.split()
            if len(parts) == 3 and parts[0].isdigit():
                charges.append(float(parts[2]))
                continue
            mulliken = False
            _add_charges(out, charges, "Mulliken", "Q-Chem の Ground-State Mulliken Net Atomic Charges")
            continue
    if mulliken and charges:
        _add_charges(out, charges, "Mulliken", "Q-Chem の Ground-State Mulliken Net Atomic Charges")
    return _keep_last_per_definition(out)


def _keep_last_per_definition(out: QcOutput) -> QcOutput:
    last: dict[str, dict] = {}
    for c in out.charges:
        last[c["definition"]] = c
    out.charges = list(last.values())
    return out


READERS = {"gaussian": read_gaussian, "gamess": read_gamess, "qchem": read_qchem}


def read_output(code: str, run_dir: Path) -> QcOutput | None:
    fn = READERS.get(code)
    if fn is None:
        return None
    path = find_output(run_dir)
    return fn(path) if path is not None else None
