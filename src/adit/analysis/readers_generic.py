
from __future__ import annotations

from pathlib import Path

import numpy as np

from adit.lang import L

CANDIDATES = (
    ("*.extxyz", "extxyz"), ("*.traj", "traj"), ("trajectory.xyz", "xyz"), ("*.xtc", None), ("*.trr", None), ("*.dcd", None),
    ("*.nc", "netcdftrajectory"), ("vasprun.xml", "vasprun"), ("*.log", "gaussian-out"), ("*.out", None),
    ("*.xyz", "xyz"), ("*.pdb", "proteindatabank"), ("*.gro", "gromacs"), ("*.cif", "cif"), ("*.gen", "gen"),
    ("POSCAR", "vasp"), ("CONTCAR", "vasp"),
)
MAX_CANDIDATE_BYTES = 512 * 2 ** 20


def ase_readable(path: Path, fmt: str | None = None):
    from ase.io import read
    from ase.io.formats import filetype

    try:
        kind = fmt or filetype(str(path), read=True)
    except Exception:
        return None
    try:
        frames = read(str(path), index=":", format=kind)
    except Exception:
        try:
            frames = read(str(path), format=kind)
        except Exception:
            return None
    frames = frames if isinstance(frames, list) else [frames]
    return (frames, kind) if frames and len(frames[0]) else None


def find_readable(run_dir: Path) -> tuple[list, str, Path] | None:
    from adit.analysis import trajectory_ext as ext

    run_dir = Path(run_dir)
    for pattern, fmt in CANDIDATES:
        for path in sorted(run_dir.glob(pattern)):
            if not path.is_file() or path.stat().st_size > MAX_CANDIDATE_BYTES:
                continue
            if ext.needs_external(path):
                try:
                    frames, lib = ext.read_external(path, ext.find_topology(run_dir, path))
                except ext.ExternalTrajectoryError:
                    continue
                return frames, lib, path
            got = ase_readable(path, fmt)
            if got is None:
                continue
            frames, kind = got
            return frames, kind, path
    return None


def energies_from(frames) -> list[float]:
    out = []
    for fr in frames:
        try:
            out.append(float(fr.get_potential_energy()))
        except Exception:
            return []
    return out


def read_generic(run_dir: Path, code: str):
    from adit.analysis.readers import RunData

    run_dir = Path(run_dir)
    found = find_readable(run_dir)
    if found is None:
        raise ValueError(L(
            f"{run_dir} には、ADIT が読める出力がありません ({code} には専用の読み取りがなく、"
            "ASE が読める軌跡・構造のファイルも見つかりませんでした)。"
            "extxyz や xyz など ASE が読める形で軌跡を書き出すと、構造から出せる解析はできます",
            f"{run_dir} has no output ADIT can read (there is no dedicated reader for {code}, and no "
            "trajectory or structure file that ASE can read). Write the trajectory in a format ASE reads "
            "(extxyz, xyz, ...) to enable the structure-based analyses"))
    frames, kind, path = found
    energies = energies_from(frames)
    notes = [L(f"{code} には専用の読み取りがないので、{path.name} を ASE の {kind} として読みました "
               f"({len(frames)} フレーム)。構造から出せる解析 (結合長・RDF・MSD・空間群・粉末回折) だけができます",
               f"there is no dedicated reader for {code}, so {path.name} was read with ASE as {kind} "
               f"({len(frames)} frames); only the structure-based analyses (bond lengths, RDF, MSD, space group, XRD) are available")]
    if not energies:
        notes.append(L("エネルギーの推移はこのファイルから読めません (温度・圧力も同じ)",
                       "no energies could be read from this file (the same applies to temperature and pressure)"))
    return RunData(code=code, run_dir=run_dir, frames=frames, energies_ev=energies,
                   frame_source=path.name, notes=notes)
