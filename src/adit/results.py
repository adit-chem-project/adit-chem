"""Summarise whether a run finished, and collect the output lines that explain a failure."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.data import covalent_radii, atomic_numbers

from adit.lang import L
from adit.structure import from_file

HARTREE_EV = 27.211386245988  # CODATA 2018


@dataclass
class RunSummary:
    converged: bool
    geometry_steps: int
    scc_iterations_last: int
    mermin_energy_hartree: float | None
    bonds: list[tuple[str, float]] = field(default_factory=list)  # ("O1-H2", 0.967)
    task_type: str | None = None
    finished: bool | None = None
    exit_code: int | None = None
    completion_assessed: bool = True
    convergence_assessed: bool = True

    def status_line(self) -> str:
        code_ok = self.exit_code == 0 if self.exit_code is not None else None
        bad_code = L(f" (終了コード {self.exit_code})", f" (exit code {self.exit_code})") if self.exit_code not in (None, 0) else ""
        if self.task_type == "geometry_optimization":
            if not self.convergence_assessed:
                return L("構造最適化: 最後まで走りました。収束したかは未判定です (このコードは収束の印を出しません)",
                         "geometry optimization: it ran to the end; convergence was not assessed (this code prints no convergence message)") + bad_code
            if self.converged:
                return L("構造最適化: 収束しました", "geometry optimization: converged") + bad_code
            return L("構造最適化: 収束しませんでした (出力に収束の印がありません)", "geometry optimization: not converged (no convergence message in the output)") + bad_code
        if code_ok is False:
            return L(f"計算: 終了コードが 0 ではありません ({self.exit_code})", f"calculation: exit code is not 0 ({self.exit_code})")
        if not self.completion_assessed:
            return L("計算: 正常終了は未判定です (出力と終了コードを確認してください)",
                     "calculation: normal termination was not assessed (check the output and exit code)")
        done = self.finished if self.finished is not None else self.converged
        if code_ok and done:
            return L("計算: 指定の計算が終わりました (終了コード 0)", "calculation: the requested calculation finished (exit code 0)")
        if done:
            return L("計算: 指定の計算が終わりました (出力の最後に正常終了の印があります)", "calculation: the requested calculation finished (normal-termination message in the output)")
        if code_ok:
            return L("計算: 終了コードは 0 ですが、出力に正常終了の印がありません", "calculation: exit code 0, but no normal-termination message in the output")
        return L("計算: 出力に正常終了の印がありません (途中で止まった可能性があります)", "calculation: no normal-termination message in the output (it may have stopped midway)")

    def text(self) -> str:
        steps_label = {"molecular_dynamics": L("MD のステップ数 (出力に記録された数)", "MD steps (as recorded in the output)")}.get(
            self.task_type or "", L("構造ステップ数", "geometry steps"))
        lines = [self.status_line()]
        if self.task_type not in ("vibrations", "single_point", "band_structure"):
            lines.append(f"{steps_label}: {self.geometry_steps}")
        if self.scc_iterations_last >= 0:
            lines.append(L(f"最終ステップの SCC 反復回数: {self.scc_iterations_last}", f"SCC iterations in the last step: {self.scc_iterations_last}"))
        if self.mermin_energy_hartree is not None:
            e = self.mermin_energy_hartree
            lines.append(L("全エネルギー (Mermin)", "total energy (Mermin)") + f": {e:.10f} Hartree = {e * HARTREE_EV:.4f} eV")
        lines.append(L("結合長 (共有結合半径の和 ×1.2 以内の原子対):", "bond lengths (atom pairs within 1.2 x the sum of covalent radii):"))
        lines += [f"  {n}: {d:.4f} Å" for n, d in self.bonds] or [L("  (該当なし)", "  (none)")]
        return "\n".join(lines)


_DONE_MARKS = {
    "dftbplus": ("DFTB+ running times",),
    "xtb": ("normal termination of xtb",),
    "espresso": ("JOB DONE",),
    "orca": ("ORCA TERMINATED NORMALLY",),
    "lammps": ("Total wall time",),
    "gromacs": ("Finished mdrun",),
    "cp2k": ("PROGRAM ENDED AT",),
    "mlip": ("adit-mlip: done",),
    "openmm": ("adit-openmm: done",),
    "abinit": ("Calculation completed.",),
    "psi4": ("*** Psi4 exiting successfully",),
}
_DONE_FILE = {"gromacs": "adit.log"}


_GENERATED_ALWAYS = {"spec.json", "README.txt", "submit.sh", "analyze.py", "code_version.txt", "stages.json",
                     "scan.json", "compare.json", "conversion.json"}


def not_run_yet(run_dir: Path | str) -> bool:
    d = Path(run_dir).expanduser()
    if not d.is_dir():
        return False
    from adit.provenance import read_provenance

    prov = read_provenance(d) or {}
    known = set(_GENERATED_ALWAYS)
    for key in ("inputs", "files", "generated_files"):
        known |= {item.get("name", "") for item in (prov.get(key) or [])}
    if not prov:
        # Without provenance only an ADIT-made directory (spec.json present) can be judged; other
        # directories are left to the readers so that foreign output names are not reported as "not run".
        outputs = [p for p in d.iterdir() if p.is_file() and p.name in
                   ("output.log", "log.lammps", "adit.log", "results.json", "detailed.out", "input.abo")]
        return (d / "spec.json").is_file() and not outputs
    rest = [p for p in d.rglob("*") if p.is_file()
            and p.relative_to(d).as_posix() not in known
            and not p.relative_to(d).as_posix().startswith("analysis/")]
    return not rest


# ABINIT 10.0.3 / Psi4 1.11
_ERROR_MARKS: dict[str, tuple[tuple[str, ...], tuple[str, ...]]] = {
    "dftbplus": (("output.log", "dftb.out"), ("ERROR!", "ERROR STOP")),
    "xtb": (("output.log", "xtb.log"), ("[ERROR]", "abnormal termination of xtb", "Fortran runtime error")),
    "espresso": (("output.log",), ("Error in routine", "%%%%%%%%%%%")),
    "cp2k": (("output.log", "adit.log"), ("[ABORT]",)),
    "lammps": (("log.lammps", "output.log"), ("ERROR:", "ERROR on proc")),
    "gromacs": (("output.log", "grompp.log", "md.log"), ("Fatal error:", "Error in user input:")),
    "abinit": (("output.log", "input.abo"), ("--- !ERROR",)),
    "psi4": (("output.dat", "output.log"), ("Psi4 encountered an error", "PsiException")),
}
MAX_ERROR_LINES = 6


def failure_lines(run_dir: Path | str, code: str = "") -> list[str]:
    d = Path(run_dir).expanduser()
    if not code:
        spec = d / "spec.json"
        if spec.is_file():
            import json

            try:
                code = json.loads(spec.read_text(encoding="utf-8")).get("method", {}).get("code", "")
            except ValueError:
                code = ""
    names, marks = _ERROR_MARKS.get(code, ((), ()))
    out: list[str] = []
    for name in names:
        path = d / name
        if not path.is_file():
            continue
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
        for number, line in enumerate(lines, 1):
            text = line.strip()
            if not any(mark in line for mark in marks) or not text.strip("%*| "):
                continue
            out.append(f"{name}:{number}  {text}")
            if text.endswith(":") or text == "--- !ERROR":
                for follow in lines[number:number + 6]:
                    body = follow.strip().strip("%*| ")
                    if body and not body.startswith(("src_file", "src_line", "mpi_rank", "message")):
                        out.append(f"{name}:{number + 1}  {body}")
                        break
            if len(out) >= MAX_ERROR_LINES:
                return out[:MAX_ERROR_LINES]
    return out


def failure_note(run_dir: Path | str, code: str = "") -> str:
    lines = failure_lines(run_dir, code)
    if not lines:
        return ""
    return L("計算が止まった印が出力にあります (コードが書いた行をそのまま写します): ",
             "the output contains a message that the calculation stopped (quoted as the code wrote it): ") + " / ".join(lines)


def _task_type(run_dir: Path) -> str | None:
    try:
        from adit.project import load_project
        return load_project(run_dir).task.type
    except Exception:
        return None


def _finished(run_dir: Path, code: str) -> bool | None:
    if code == "nwchem":
        return None
    if code == "vasp":
        x = run_dir / "vasprun.xml"
        return x.is_file() and x.read_text(encoding="utf-8", errors="replace").rstrip().endswith("</modeling>") if x.is_file() else None
    log = run_dir / _DONE_FILE.get(code, "output.log")
    if not log.is_file():
        return False
    with open(log, "rb") as f:
        f.seek(max(0, log.stat().st_size - 65536))
        text = f.read().decode("utf-8", "replace")
    return any(m in text for m in _DONE_MARKS[code])


def read_results_tag(path: Path) -> dict[str, np.ndarray]:
    out: dict[str, np.ndarray] = {}
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    i = 0
    while i < len(lines):
        m = re.match(r"^(\S+)\s*:(\w+):(\d+):(.*)$", lines[i])
        if not m:
            i += 1
            continue
        key, typ, rank, shape = m.group(1), m.group(2), int(m.group(3)), m.group(4).strip()
        n = int(np.prod([int(x) for x in shape.split(",")])) if rank > 0 else 1
        vals: list[float] = []
        i += 1
        while i < len(lines) and len(vals) < n and not re.match(r"^\S+\s*:\w+:", lines[i]):
            vals += [float(x) for x in lines[i].split()] if typ == "real" else [int(x) for x in lines[i].split()]
            i += 1
        out[key] = np.asarray(vals)
    return out


def bonds_of(atoms: Atoms, factor: float = 1.2) -> list[tuple[str, float]]:
    from adit.analysis.compute import bonds
    return bonds(atoms, factor)


_ = (covalent_radii, atomic_numbers)


def _summarize_extra(run_dir: Path, code: str) -> RunSummary:
    from adit.analysis.readers import load_run
    r = load_run(run_dir)
    log = {"lammps": "log.lammps", "gromacs": "adit.log", "cp2k": "output.log"}[code]
    text = (run_dir / log).read_text(encoding="utf-8", errors="replace") if (run_dir / log).is_file() else ""
    mark = {"lammps": r"Stopping criterion = (energy|force) tolerance", "gromacs": r"converged to Fmax", "cp2k": r"GEOMETRY OPTIMIZATION COMPLETED"}[code]
    converged = re.search(mark, text) is not None
    energy = r.energies_ev[-1] / HARTREE_EV if r.energies_ev else None
    final = r.final
    return RunSummary(converged, len(r.energies_ev), -1, energy, bonds_of(final) if final is not None else [])


def _summarize_mlip(run_dir: Path) -> RunSummary:
    import json

    from adit.analysis.readers import load_run

    r = load_run(run_dir)
    p = run_dir / "results.json"
    res = {}
    if p.is_file():
        try:
            obj = json.loads(p.read_text(encoding="utf-8", errors="replace"))
            res = obj if isinstance(obj, dict) else {}
        except ValueError:
            res = {}
    steps = int(res["steps"]) if isinstance(res.get("steps"), (int, float)) else len(r.energies_ev)
    converged = bool(res.get("converged")) if res.get("converged") is not None else bool(res)
    energy = r.energies_ev[-1] / HARTREE_EV if r.energies_ev else None
    return RunSummary(converged, steps, -1, energy, bonds_of(r.final) if r.final is not None else [])


def _summarize_openmm(run_dir: Path) -> RunSummary:
    from adit.analysis.readers import load_run

    r = load_run(run_dir)
    table = r.extra_tables.get("openmm", {})
    steps = int(table["steps"]) if isinstance(table.get("steps"), (int, float)) else len(r.energies_ev)
    energy = r.energies_ev[-1] / HARTREE_EV if r.energies_ev else None
    return RunSummary(False, steps, -1, energy, bonds_of(r.final) if r.final is not None else [], convergence_assessed=False)


def _summarize_abinit(run_dir: Path) -> RunSummary:
    from adit.analysis.readers import load_run

    r = load_run(run_dir)
    log = run_dir / "output.log"
    text = log.read_text(encoding="utf-8", errors="replace") if log.is_file() else ""
    converged = "gradients are converged" in text
    energy = r.energies_ev[-1] / HARTREE_EV if r.energies_ev else None
    return RunSummary(converged, len(r.energies_ev), -1, energy, bonds_of(r.final) if r.final is not None else [])


def _summarize_psi4(run_dir: Path) -> RunSummary:
    from adit.analysis.readers import load_run

    r = load_run(run_dir)
    table = r.extra_tables.get("psi4", {})
    converged = bool(table.get("converged"))
    steps = int(table.get("optimization_steps") or 0) or len(r.energies_ev)
    energy = r.energies_ev[-1] / HARTREE_EV if r.energies_ev else None
    return RunSummary(converged, steps, -1, energy, bonds_of(r.final) if r.final is not None else [])


def summarize_run(run_dir: Path | str, exit_code: int | None = None) -> RunSummary:
    """Summarise a run directory. Pass exit_code when the caller knows it."""
    run_dir = Path(run_dir)
    for marker, code, fn in (("INCAR", "vasp", _summarize_vasp), ("xtb.inp", "xtb", _summarize_xtb),
                             ("pw.in", "espresso", _summarize_espresso), ("orca.inp", "orca", _summarize_orca),
                             ("nwchem.nw", "nwchem", _summarize_nwchem),
                             ("in.lammps", "lammps", lambda d: _summarize_extra(d, "lammps")),
                             ("grompp.mdp", "gromacs", lambda d: _summarize_extra(d, "gromacs")),
                             ("cp2k.inp", "cp2k", lambda d: _summarize_extra(d, "cp2k")),
                             ("mlip_settings.json", "mlip", _summarize_mlip),
                             ("openmm_settings.json", "openmm", _summarize_openmm),
                             ("input.abi", "abinit", _summarize_abinit),
                             ("input.dat", "psi4", _summarize_psi4),
                             (None, "dftbplus", _summarize_dftb)):
        if marker is None or (run_dir / marker).is_file():
            s = fn(run_dir)
            s.task_type, s.finished, s.exit_code = _task_type(run_dir), _finished(run_dir, code), exit_code
            return s
    raise AssertionError("unreachable")


def _summarize_nwchem(run_dir: Path) -> RunSummary:
    # RunSummary labels its numeric field as Mermin energy. That is not the
    # NWChem SCF/DFT quantity; adit-analyze reports the latter accurately.
    return RunSummary(False, 0, -1, None, [], completion_assessed=False)


def _summarize_dftb(run_dir: Path) -> RunSummary:
    log = (run_dir / "output.log").read_text(encoding="utf-8", errors="replace") if (run_dir / "output.log").is_file() else ""
    steps = log.count("***  Geometry step:")
    last = log.rsplit("***  Geometry step:", 1)[-1] if steps else log
    scc_last = len(re.findall(r"^\s+\d+\s+-?\d\.\d+E[+-]\d+", last, re.M))
    converged = "Geometry converged" in log if steps else ("Total Energy" in log)
    energy = None
    tag = run_dir / "results.tag"
    if tag.is_file():
        data = read_results_tag(tag)
        if "mermin_energy" in data and data["mermin_energy"].size:
            energy = float(data["mermin_energy"][0])
    geom = run_dir / "geom.out.gen"
    atoms = from_file(geom if geom.is_file() else run_dir / "geometry.gen")
    if (run_dir / "geo_end.xyz").is_file():
        from adit.analysis.trajectory import Trajectory
        traj = Trajectory(run_dir / "geo_end.xyz", "xyz", cell=atoms.cell if any(atoms.pbc) else None, pbc=atoms.pbc if any(atoms.pbc) else None)
        if len(traj):
            atoms = traj[-1]
    return RunSummary(converged, steps, scc_last, energy, bonds_of(atoms))


def _summarize_vasp(run_dir: Path) -> RunSummary:
    log = (run_dir / "output.log").read_text(encoding="utf-8", errors="replace") if (run_dir / "output.log").is_file() else ""
    stdout_seen = bool(log)
    if not log and (run_dir / "OSZICAR").is_file():
        log = (run_dir / "OSZICAR").read_text(encoding="utf-8", errors="replace")
    f_lines = [l for l in log.splitlines() if " F= " in l]
    energy_ev = float(f_lines[-1].split("F=")[1].split()[0]) if f_lines else None
    last_block = log.rsplit(" F= ", 1)[0].rsplit("N       E", 1)[-1] if f_lines else log
    scc_last = len(re.findall(r"^\s*(?:DAV|RMM|CG)\s*:", last_block, re.M))
    incar = (run_dir / "INCAR").read_text(encoding="utf-8", errors="replace")
    converged = "reached required accuracy" in log or (bool(f_lines) and "NSW = 0" in incar)
    if not stdout_seen and not converged:
        m = re.search(r"^\s*NSW\s*=\s*(\d+)", incar, re.M)
        x = run_dir / "vasprun.xml"
        closed = x.is_file() and x.read_text(encoding="utf-8", errors="replace").rstrip().endswith("</modeling>")
        converged = bool(m and f_lines and closed and len(f_lines) < int(m.group(1)))
    geom = run_dir / "CONTCAR"
    atoms = from_file(geom if geom.is_file() and geom.stat().st_size > 0 else run_dir / "POSCAR")
    return RunSummary(converged, len(f_lines), scc_last, energy_ev / HARTREE_EV if energy_ev is not None else None, bonds_of(atoms))


def _summarize_xtb(run_dir: Path) -> RunSummary:
    log = (run_dir / "output.log").read_text(encoding="utf-8", errors="replace") if (run_dir / "output.log").is_file() else ""
    m = re.findall(r"TOTAL ENERGY\s+(-?\d+\.\d+) Eh", log)
    energy = float(m[-1]) if m else None
    mo = re.search(r"GEOMETRY OPTIMIZATION CONVERGED AFTER (\d+) ITERATIONS", log)
    steps = int(mo.group(1)) if mo else 0
    converged = bool(mo) or ("normal termination of xtb" in log and not (run_dir / "xtbopt.xyz").exists())
    scc = len(re.findall(r"^\s+\d+\s+-?\d+\.\d+\s+-?\d\.\d+E[+-]\d+", log.rsplit("SCC iter.", 1)[-1], re.M))
    geom = run_dir / "xtbopt.xyz"
    atoms = from_file(geom if geom.is_file() else run_dir / "struct.xyz")
    return RunSummary(converged, steps, scc, energy, bonds_of(atoms))


RY_EV = 13.605693122994


def _summarize_espresso(run_dir: Path) -> RunSummary:
    log = (run_dir / "output.log").read_text(encoding="utf-8", errors="replace") if (run_dir / "output.log").is_file() else ""
    en = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+) Ry", log, re.M)
    energy_ry = float(en[-1]) if en else None
    converged = "JOB DONE" in log and ("bfgs converged" in log or "calculation      = 'scf'" in (run_dir / "pw.in").read_text(encoding="utf-8", errors="replace"))
    scc = len(re.findall(r"^\s+iteration #\s*\d+", log.rsplit("!    total energy", 1)[0].rsplit("Self-consistent Calculation", 1)[-1], re.M))
    try:
        from ase.io import read
        atoms = read(run_dir / "output.log", format="espresso-out", index=-1)
    except Exception:
        atoms = from_file(run_dir / "pw.in")
    return RunSummary(converged, len(en), scc, energy_ry * RY_EV / HARTREE_EV if energy_ry is not None else None, bonds_of(atoms))


def _summarize_orca(run_dir: Path) -> RunSummary:
    log = (run_dir / "output.log").read_text(encoding="utf-8", errors="replace") if (run_dir / "output.log").is_file() else ""
    en = re.findall(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)", log)
    converged = "ORCA TERMINATED NORMALLY" in log and ("THE OPTIMIZATION HAS CONVERGED" in log or "Opt" not in (run_dir / "orca.inp").read_text(encoding="utf-8", errors="replace"))
    scc = len(re.findall(r"^\s*\d+\s+-\d+\.\d+\s+-?\d\.\d+e[+-]\d+", log.rsplit("SCF ITERATIONS", 1)[-1], re.M))
    geom = run_dir / "orca.xyz"
    atoms = from_file(geom if geom.is_file() else run_dir / "orca.inp") if geom.is_file() else None
    if atoms is None:
        from adit.spec import CalculationSpec
        atoms = CalculationSpec.load(run_dir / "spec.json").atoms
    return RunSummary(converged, len(en), scc, float(en[-1]) if en else None, bonds_of(atoms))
