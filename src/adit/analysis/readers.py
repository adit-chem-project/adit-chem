"""Read the output of each calculation code into the common RunData form. Missing values stay None."""

from __future__ import annotations

import json
import re
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import iread, read

from adit.analysis.trajectory import _QE_UNITS, Trajectory
from adit.lang import L

HARTREE_EV = 27.211386245988
RY_EV = 13.605693122994
BOHR_ANG = 0.529177210903
EA0_DEBYE = 2.541746473


@dataclass
class RunData:
    code: str
    run_dir: Path
    frames: Sequence[Atoms] = field(default_factory=list)
    energies_ev: list[float] = field(default_factory=list)
    temperatures_k: list[float] | None = None
    times_fs: list[float] | None = None
    eigenvalues_ev: np.ndarray | None = None
    eigen_weights: np.ndarray | None = None
    fermi_ev: float | None = None
    fermi_is_homo: bool = False
    force_max_ev_ang: float | None = None
    force_source: str = ""
    frequencies_cm1: list[float] | None = None
    ir_intensities: list[float] | None = None
    raman_activities: list[float] | None = None
    notes: list[str] = field(default_factory=list)
    frame_dt_fs: float | None = None
    frame_source: str = ""
    charges: list[dict] = field(default_factory=list)
    thermo: dict | None = None
    electronic: dict = field(default_factory=dict)
    series: dict[str, dict] = field(default_factory=dict)
    extra_tables: dict[str, dict] = field(default_factory=dict)

    @property
    def final(self) -> Atoms | None:
        return self.frames[-1] if self.frames else None


_INPUTS = (("dftb_in.hsd", "dftbplus"), ("INCAR", "vasp"), ("xtb.inp", "xtb"), ("pw.in", "espresso"), ("orca.inp", "orca"), ("nwchem.nw", "nwchem"),
           ("in.lammps", "lammps"), ("grompp.mdp", "gromacs"), ("cp2k.inp", "cp2k"), ("mlip_settings.json", "mlip"),
           ("openmm_settings.json", "openmm"), ("input.abi", "abinit"), ("input.dat", "psi4"),
           ("gaussian.gjf", "gaussian"), ("gamess.inp", "gamess"), ("qchem.in", "qchem"), ("openmx.dat", "openmx"),
           ("amber.in", "amber"), ("namd.conf", "namd"), ("grrm.com", "grrm"), ("dftb.inp", "dcdftbmd"))


def detect_code(run_dir: Path) -> str | None:
    for name, code in _INPUTS:
        if (run_dir / name).is_file():
            return code
    return None


def load_run(run_dir: Path | str, code: str | None = None) -> RunData:
    """Read a run directory. Pass code to read outputs when the input files are not present."""
    run_dir = Path(run_dir)
    known = {c for _, c in _INPUTS} | {"lammps", "gromacs", "cp2k", "openmm", "abinit", "psi4", "mlip"}
    if code:
        if code not in known:
            raise ValueError(L(f"知らないコードです: {code} (使えるもの: {', '.join(sorted(known))})",
                               f"unknown code: {code} (known: {', '.join(sorted(known))})"))
    else:
        code = detect_code(run_dir)
    if code is None:
        names = " / ".join(n for n, _ in _INPUTS)
        raise ValueError(L(f"{run_dir} に計算コードの入力がありません ({names})", f"{run_dir} has no input of a supported code ({names})"))
    if code in ("lammps", "gromacs", "cp2k", "openmm", "abinit", "psi4"):
        from adit.analysis import readers_extra
        return {"lammps": readers_extra.read_lammps, "gromacs": readers_extra.read_gromacs, "cp2k": readers_extra.read_cp2k,
                "openmm": readers_extra.read_openmm, "abinit": readers_extra.read_abinit,
                "psi4": readers_extra.read_psi4}[code](run_dir)
    special = {"dftbplus": _read_dftb, "vasp": _read_vasp, "xtb": _read_xtb, "espresso": _read_espresso,
               "orca": _read_orca, "nwchem": _read_nwchem, "mlip": _read_mlip}
    if code in special:
        return special[code](run_dir)
    if code in ("gaussian", "gamess", "qchem"):
        got = _read_qc(run_dir, code)
        if got is not None:
            return got
    if code == "openmx":
        got = _read_openmx(run_dir)
        if got is not None:
            return got
    from adit.analysis.readers_generic import read_generic

    return read_generic(run_dir, code)


def _read_qc(run_dir: Path, code: str) -> RunData | None:
    from adit.analysis.readers_qc import read_output

    got = read_output(code, run_dir)
    if got is None or not (got.energies_ev or got.frequencies_cm1 or got.frames):
        return None
    r = RunData(code, run_dir)
    r.energies_ev = got.energies_ev
    r.frames = got.frames
    r.frame_source = got.path.name
    r.frequencies_cm1 = got.frequencies_cm1 or None
    r.ir_intensities = got.ir_intensities or None
    r.raman_activities = got.raman_activities or None
    r.charges = got.charges
    if not got.frames:
        r.notes.append(L(f"{got.path.name} から構造を読めませんでした (エネルギーと振動数だけ読みました)",
                         f"could not read a structure from {got.path.name}; only energies and frequencies were read"))
    return r


def _read_openmx(run_dir: Path) -> RunData | None:
    from adit.analysis.readers_openmx import read_openmx
    from adit.analysis.readers_qc import find_output

    path = find_output(run_dir)
    if path is None:
        return None
    got = read_openmx(path)
    if not (got.energies_ev or got.frames or got.populations):
        return None
    r = RunData("openmx", run_dir)
    r.energies_ev = got.energies_ev
    r.frames = got.frames
    r.frame_source = path.name
    if got.eigenvalues_ev:
        r.eigenvalues_ev = np.array(got.eigenvalues_ev, dtype=float)
    if got.populations:
        r.charges = [{"definition": L("Mulliken の電子数", "Mulliken populations"),
                      "values": got.populations, "source": f"{path.name} の Mulliken populations",
                      "unit": "e", "note": L("電子の数です (価電子数を引いた電荷ではありません)",
                                             "these are electron counts, not charges relative to the valence")}]
    if got.homo_index is not None and got.eigenvalues_ev:
        r.notes.append(L(f"{path.name}: 最高被占準位は {got.homo_index} 番目です",
                         f"{path.name}: the highest occupied level is number {got.homo_index}"))
    if got.note_no_cell:
        r.notes.append(L(f"{path.name} にセルが書かれていないので、分率座標を Å に直せません (構造は読みませんでした)",
                         f"{path.name} has no cell, so the fractional coordinates cannot be converted to Å "
                         "(the structure was not read)"))
    return r


def _text(p: Path) -> str:
    return p.read_text(encoding="utf-8", errors="replace") if p.is_file() else ""


def _lines(p: Path) -> Iterator[str]:
    if not p.is_file():
        return
    with open(p, encoding="utf-8", errors="replace") as f:
        yield from f


def _grep(p: Path, pattern: str, group: int = 1) -> list[str]:
    rx = re.compile(pattern)
    out = []
    for line in _lines(p):
        m = rx.search(line)
        if m:
            out.append(m.group(group))
    return out


def _read_nwchem(d: Path) -> RunData:
    # Read only the total energy; NWChem geometry/frequency parsing is not validated.
    # Output label: https://nwchemgit.github.io/Sample.html
    result = RunData("nwchem", d)
    energy = re.compile(r"^\s*Total (?:SCF|DFT) energy\s*=\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+-]?\d+)?)\s*$")
    for line in _lines(d / "output.log"):
        match = energy.match(line)
        if match:
            result.energies_ev.append(float(match.group(1).replace("D", "E").replace("d", "e")) * HARTREE_EV)
    result.notes.append(L("NWChem の Total SCF/DFT energy だけを読みました。構造・振動数と正常終了は未判定です。",
                          "Only NWChem's Total SCF/DFT energy was read; structure, frequencies, and normal termination were not assessed."))
    return result


def _read_xyz_frames(path: Path, cell=None, pbc=None) -> list[Atoms]:
    return list(Trajectory(path, "xyz", cell=cell, pbc=pbc)) if Path(path).is_file() else []


def _last_xyz_frame(path: Path) -> list[Atoms]:
    if not path.is_file():
        return []
    t = Trajectory(path, "xyz")
    return [t[-1]] if len(t) else []


# ---------------- DFTB+ ----------------
def _read_dftb(d: Path) -> RunData:
    r = RunData("dftbplus", d)
    try:
        start = read(d / "geometry.gen")
    except Exception:
        start = None
    cell, pbc = (start.cell, start.pbc) if start is not None and any(start.pbc) else (None, None)
    md = d / "md.out"
    if md.is_file() and md.stat().st_size > 0:
        e_rx = re.compile(r"Total MD Energy:\s+-?\d+\.\d+ H\s+(-?\d+\.\d+) eV")
        t_rx = re.compile(r"MD Temperature:\s+[-\d.]+ au\s+(-?\d+\.\d+) K")
        energies, temps = [], []
        for line in _lines(md):
            m = e_rx.search(line)
            if m:
                energies.append(float(m.group(1))); continue
            m = t_rx.search(line)
            if m:
                temps.append(float(m.group(1)))
        r.energies_ev, r.temperatures_k = energies, temps
        hsd = _text(d / "dftb_pin.hsd")
        m = re.search(r"TimeStep\s*(\[(\w+)\])?\s*=\s*([-\d.E+e]+)", hsd)
        fs = None
        if m:
            v = float(m.group(3)); unit = (m.group(2) or "au").lower()
            fs = v if unit in ("fs", "femtosecond") else v * 1000 if unit in ("ps", "picosecond") else v / 41.341373
        m2 = re.search(r"MDRestartFrequency\s*=\s*(\d+)", hsd)
        every = int(m2.group(1)) if m2 else 1
        if fs is not None:
            r.times_fs = [i * every * fs for i in range(len(r.energies_ev))]
            r.frame_dt_fs = every * fs
        if (d / "geo_end.xyz").is_file():
            r.frames = Trajectory(d / "geo_end.xyz", "xyz", cell=cell, pbc=pbc)
            r.frame_source = "geo_end.xyz"
    else:
        r.energies_ev = [float(x) for x in _grep(d / "output.log", r"Total Energy:\s+-?\d+\.\d+ H\s+(-?\d+\.\d+) eV")]
        for name in ("geom.out.gen", "geometry.gen"):
            if (d / name).is_file():
                try:
                    r.frames = [read(d / name)]
                    break
                except Exception:
                    pass
    band = _text(d / "band.out")
    if band:
        eigs, ws = [], []
        w = 1.0
        for l in band.splitlines():
            if l.startswith("KPT"):
                m = re.search(r"KWEIGHT\s+([\d.Ee+-]+)", l); w = float(m.group(1)) if m else 1.0
                continue
            p = l.split()
            if len(p) == 3:
                try:
                    eigs.append(float(p[1])); ws.append(w)
                except ValueError:
                    pass
        if eigs:
            r.eigenvalues_ev, r.eigen_weights = np.array(eigs), np.array(ws)
    detailed = _text(d / "detailed.out")
    m = re.search(r"Fermi level:\s+[-\d.]+ H\s+(-?\d+\.\d+) eV", detailed)
    if m:
        r.fermi_ev = float(m.group(1))
    q = _dftb_charges(detailed, md_run=bool(r.temperatures_k))
    if q:
        r.charges.append(q)
    hess = d / "hessian.out"
    if hess.is_file() and r.final is not None:
        r.frequencies_cm1 = frequencies_from_hessian(np.array(_text(hess).split(), dtype=float), r.final)
        r.notes.append(L("振動数は hessian.out (Hartree/Bohr^2) を質量重み付きで対角化して求めた (並進・回転は落としていない)", "frequencies come from mass-weighted diagonalization of hessian.out (Hartree/Bohr^2); translations and rotations are not projected out"))
    return r


def _dftb_charges(detailed: str, md_run: bool = False) -> dict | None:
    lines = detailed.splitlines()
    for k, l in enumerate(lines):
        if l.strip().startswith("Atomic gross charges"):
            vals = []
            for l2 in lines[k + 2:]:
                w = l2.split()
                if len(w) != 2:
                    break
                try:
                    vals.append(float(w[1]))
                except ValueError:
                    break
            if not vals:
                return None
            note = L("DFTB+ の Mulliken 電荷 (価電子の数 − Mulliken の電子数。負なら電子が多い)",
                     "DFTB+ Mulliken charges (valence electrons minus Mulliken population; negative means excess electrons)")
            if md_run:
                note += L("。MD では detailed.out は最後に書かれたステップの値", "; for MD, detailed.out holds the last step written")
            return {"definition": "Mulliken", "values": vals, "source": f"detailed.out:{k + 1}", "unit": "e", "note": note}
    return None


def frequencies_from_hessian(h_flat: np.ndarray, atoms: Atoms) -> list[float]:
    n = len(atoms) * 3
    if h_flat.size < n * n or not np.all(np.isfinite(h_flat[: n * n])):
        raise ValueError(L(f"hessian.out が途中で切れているか壊れています ({n}×{n} = {n * n} 個の数が要るところ、{h_flat.size} 個)",
                           f"hessian.out is truncated or broken ({n}x{n} = {n * n} numbers needed, {h_flat.size} found)"))
    h = h_flat[: n * n].reshape(n, n)
    h = 0.5 * (h + h.T)
    masses_au = np.repeat(atoms.get_masses() * 1822.888486209, 3)
    hw = h / np.sqrt(np.outer(masses_au, masses_au))
    lam = np.linalg.eigvalsh(hw)
    au_to_cm1 = 219474.6313705
    return [float(np.sign(x) * np.sqrt(abs(x)) * au_to_cm1) for x in lam]


# ---------------- xtb ----------------
_XTB_THERMO = (
    ("total_free_energy_box", r"::\s*total free energy\s+(-?[\d.]+)\s*Eh", "total free energy"),
    ("zpe", r"::\s*zero point energy\s+(-?[\d.]+)\s*Eh", "zero point energy"),
    ("g_rrho_wo_zpve", r"::\s*G\(RRHO\) w/o ZPVE\s+(-?[\d.]+)\s*Eh", "G(RRHO) w/o ZPVE"),
    ("g_rrho_contrib", r"::\s*G\(RRHO\) contrib\.\s+(-?[\d.]+)\s*Eh", "G(RRHO) contrib."),
    ("total_energy", r"\|\s*TOTAL ENERGY\s+(-?[\d.]+)\s*Eh", "TOTAL ENERGY"),
    ("total_enthalpy", r"\|\s*TOTAL ENTHALPY\s+(-?[\d.]+)\s*Eh", "TOTAL ENTHALPY"),
    ("total_free_energy", r"\|\s*TOTAL FREE ENERGY\s+(-?[\d.]+)\s*Eh", "TOTAL FREE ENERGY"),
)


def _xtb_thermo(log: Path) -> dict | None:
    rx = [(k, re.compile(p), lab) for k, p, lab in _XTB_THERMO]
    found: dict[str, tuple[float, int, str]] = {}
    rows, in_table, table_line = [], False, None
    for no, line in enumerate(_lines(log), start=1):
        if "T/K" in line and "H(0)-H(T)+PV" in line:
            in_table, rows, table_line = True, [], no
            continue
        if in_table:
            w = line.split()
            if len(w) == 5:
                try:
                    rows.append([float(x) for x in w]); continue
                except ValueError:
                    pass
            if rows and line.strip().startswith("---"):
                in_table = False
            continue
        for k, p, lab in rx:
            m = p.search(line)
            if m:
                found[k] = (float(m.group(1)), no, lab)
                break
    if not any(k in found for k in ("zpe", "total_enthalpy", "total_free_energy", "g_rrho_contrib")):
        return None
    items = [{"key": k, "label": lab, "value_eh": v, "value_ev": v * HARTREE_EV, "line": f"output.log:{no}"}
             for k, (v, no, lab) in sorted(found.items(), key=lambda kv: kv[1][1])]
    if len(rows) >= 1:
        t, _hpv, h_t, ts, g_t = rows[-1] if len(rows) == 1 else rows[0]
        for key, lab, v in (("h_t", "H(T)", h_t), ("ts", "T*S", ts), ("g_t", "G(T)", g_t)):
            if len(rows) == 1:
                items.append({"key": key, "label": lab, "value_eh": v, "value_ev": v * HARTREE_EV, "line": f"output.log:{table_line}"})
    out = {"code": "xtb", "method": "RRHO (xtb)", "source": "output.log", "temperature_k": rows[0][0] if len(rows) == 1 else None,
           "items": items, "table_rows": [{"T_K": r_[0], "H0_minus_HT_plus_PV_eh": r_[1], "H_T_eh": r_[2], "TS_eh": r_[3], "G_T_eh": r_[4]} for r_ in rows]}
    if len(rows) > 1:
        out["note"] = L("温度の表に複数の温度があります。枠の中の値がどの温度のものかは出力からは読めません (表の各行は table_rows)",
                        "the temperature table has several temperatures; which one the boxed values belong to cannot be read from the output (rows in table_rows)")
    return out


def _xtb_json(d: Path, r: RunData) -> None:
    p = d / "xtbout.json"
    if not p.is_file():
        return
    text = _text(p)
    try:
        js = json.loads(text)
    except ValueError:
        r.notes.append(L("xtbout.json を JSON として読めません", "cannot parse xtbout.json as JSON"))
        return
    method = str(js.get("method", ""))

    def line_of(key: str) -> str:
        for no, l in enumerate(text.splitlines(), start=1):
            if f'"{key}"' in l:
                return f"xtbout.json:{no}"
        return "xtbout.json"

    q = js.get("partial charges")
    if isinstance(q, list) and q:
        if "GFN2" in method:
            definition = "Mulliken (GFN2-xTB)"
            note = L("GFN2-xTB が自己無撞着に解く原子の電荷 (Mulliken の分け方。GFN2-xTB の論文 Bannwarth ら J. Chem. Theory Comput. 15, 1652 (2019))",
                     "atomic charges solved self-consistently by GFN2-xTB (Mulliken partitioning; Bannwarth et al., J. Chem. Theory Comput. 15, 1652 (2019))")
        else:
            definition = f"{method or 'xtb'} partial charges"
            note = L("xtb が xtbout.json の partial charges に書く値 (この手法での電荷の種類は、xtb の文書で確かめていません)",
                     "the values xtb writes as 'partial charges' in xtbout.json (the kind of charge for this method is not confirmed in the xtb docs)")
        r.charges.append({"definition": definition, "values": [float(x) for x in q], "source": line_of("partial charges"), "unit": "e", "note": note})
    gap = js.get("HOMO-LUMO gap / eV")
    if gap is not None:
        r.electronic.update(homo_lumo_gap_ev=float(gap), gap_source=line_of("HOMO-LUMO gap / eV"))
    dip = js.get("dipole / a.u.")
    if isinstance(dip, list) and len(dip) == 3:
        v = np.array(dip, dtype=float) * EA0_DEBYE
        r.electronic.update(dipole_debye=[float(x) for x in v], dipole_norm_debye=float(np.linalg.norm(v)),
                            dipole_source=line_of("dipole / a.u.") + L(" (a.u. × 2.541746 で D に換算)", " (a.u. x 2.541746 to D)"))


def _read_xtb(d: Path) -> RunData:
    r = RunData("xtb", d)
    trj = d / "xtb.trj"
    log = d / "output.log"
    if trj.is_file():
        t = Trajectory(trj, "xyz", comment_regex=re.compile(r"energy:\s*(-?\d+\.\d+)"))
        r.frames, r.frame_source = t, "xtb.trj"
        r.energies_ev = [v * HARTREE_EV for v in t.comment_values]
        rx = re.compile(r"^\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(-?\d+\.\d+)\s*$")
        rows = [m.groups() for m in (rx.match(l) for l in _lines(log)) if m]
        if rows:
            r.times_fs = [float(a[0]) * 1000 for a in rows]
            r.temperatures_k = [float(a[4]) for a in rows]
            r.notes.append(L("温度は output.log の MD の表 (time, <Epot>, Ekin, <T>, T, Etot) の T", "temperature is the T column of the MD table in output.log (time, <Epot>, Ekin, <T>, T, Etot)"))
        m = re.search(r"^\s*dump\s*=\s*([\d.Ee+-]+)", _text(d / "xtb.inp"), re.M)
        if m:
            r.frame_dt_fs = float(m.group(1))
    else:
        cycles = _grep(log, r"^\s*\*\s*total energy\s*:\s*(-?\d+\.\d+)\s*Eh")
        found = cycles or _grep(log, r"TOTAL ENERGY\s+(-?\d+\.\d+) Eh")
        r.energies_ev = [float(x) * HARTREE_EV for x in found]
        if (d / "xtbopt.log").is_file():
            r.frames = Trajectory(d / "xtbopt.log", "xyz")
        else:
            for name in ("xtbopt.xyz", "struct.xyz"):
                if (d / name).is_file():
                    r.frames = _last_xyz_frame(d / name)
                    break
    vib = _text(d / "vibspectrum")
    if vib:
        fr, ir = [], []
        for l in vib.splitlines():
            p = l.split()
            if len(p) >= 4 and p[0].isdigit():
                try:
                    float(p[1])
                    fr.append(float(p[1])); ir.append(float(p[2]))
                except ValueError:
                    try:
                        fr.append(float(p[2])); ir.append(float(p[3]))
                    except ValueError:
                        pass
        r.frequencies_cm1, r.ir_intensities = fr, ir
        g = _grep(log, r"GRADIENT NORM\s+(-?[\d.]+)\s+Eh/")
        if g and r.final is not None:
            r.force_max_ev_ang = float(g[0]) / np.sqrt(len(r.final)) * HARTREE_EV / BOHR_ANG
            r.force_source = L(f"output.log の勾配のノルム {float(g[0]):.4g} Eh/Bohr を原子数の平方根で割った下限",
                               f"lower bound from the gradient norm {float(g[0]):.4g} Eh/Bohr in output.log divided by sqrt(number of atoms)")
    r.thermo = _xtb_thermo(log)
    _xtb_json(d, r)
    if "homo_lumo_gap_ev" not in r.electronic:
        gaps = []
        for no, line in enumerate(_lines(log), start=1):
            m = re.search(r"\|\s*HOMO-LUMO GAP\s+(-?[\d.]+)\s*eV", line)
            if m:
                gaps.append((float(m.group(1)), no))
        if gaps:
            r.electronic.update(homo_lumo_gap_ev=gaps[-1][0], gap_source=f"output.log:{gaps[-1][1]}")
    return r


# ---------------- VASP ----------------
def _read_vasp(d: Path) -> RunData:
    r = RunData("vasp", d)
    vr = d / "vasprun.xml"
    vr_ok = False
    try:
        energies = [float(a.get_potential_energy()) for a in iread(str(vr), index=":", format="vasp-xml")]
        r.energies_ev, vr_ok = energies, True
    except Exception as ex:
        r.notes.append(L(f"vasprun.xml を読めません: {ex}", f"cannot read vasprun.xml: {ex}"))
    incar = _text(d / "INCAR")
    xd = d / "XDATCAR"
    xtraj = Trajectory(xd, "xdatcar") if xd.is_file() else None
    if xtraj is not None and len(xtraj) > 0:
        r.frames, r.frame_source = xtraj, "XDATCAR"
    elif vr_ok:
        r.frames, r.frame_source = Trajectory(vr, "ase", fmt="vasp-xml"), "vasprun.xml"
    elif (d / "CONTCAR").is_file():
        try:
            r.frames = [read(d / "CONTCAR")]
        except Exception:
            pass
    temps = [float(t) for t in _grep(d / "OSZICAR", r"T=\s*([\d.]+)")]
    if temps:
        r.temperatures_k = temps
        m = re.search(r"POTIM\s*=\s*([\d.]+)", incar)
        if m:
            potim = float(m.group(1))
            r.times_fs = [i * potim for i in range(len(temps))]
            ibrion = re.search(r"^\s*IBRION\s*=\s*(-?\d+)", incar, re.M)
            if ibrion is None or int(ibrion.group(1)) == 0:
                nb = re.search(r"^\s*NBLOCK\s*=\s*(\d+)", incar, re.M)
                every = int(nb.group(1)) if nb and r.frame_source == "XDATCAR" else 1
                r.frame_dt_fs = potim * every
    dos = _text(d / "DOSCAR").splitlines()
    if len(dos) > 6:
        try:
            emax, emin, nedos, efermi = [float(x) for x in dos[5].split()[:4]]
            rows = [l.split() for l in dos[6: 6 + int(nedos)]]
            e = np.array([float(x[0]) for x in rows]); t = np.array([float(x[1]) for x in rows])
            r.eigenvalues_ev, r.eigen_weights, r.fermi_ev = e, t, efermi
            r.notes.append(L("DOS は DOSCAR の全 DOS をそのまま使う (重みは DOS の値)", "DOS is the total DOS of DOSCAR as is (weights are the DOS values)"))
        except (ValueError, IndexError):
            pass
    rx = re.compile(r"^\s*\d+\s+(f|f/i)\s*=\s*[\d.]+ THz\s+[\d.]+ 2PiTHz\s+([\d.]+) cm-1")
    fr = [m.groups() for m in (rx.match(l) for l in _lines(d / "OUTCAR")) if m]
    if fr:
        r.frequencies_cm1 = [float(v) if kind == "f" else -float(v) for kind, v in fr]
    return r


# ---------------- Quantum ESPRESSO ----------------
def _read_espresso(d: Path) -> RunData:
    r = RunData("espresso", d)
    log = d / "output.log"
    if not log.is_file():
        r.notes.append(L("output.log がありません", "output.log is missing"))
    else:
        r.frames, r.frame_source = Trajectory(log, "qe"), "output.log"
    energies, temps = [], []
    fermi = homo = None
    forces_block: list[str] | None = None
    forces_done = False
    new_config = False
    for line in _lines(log):
        if ("site n." in line and "positions (alat units)" in line) or line.startswith("ATOMIC_POSITIONS"):
            new_config = True
        elif line.startswith("!") and "total energy" in line:
            if new_config:
                energies.append(float(line.split("=")[1].split()[0]) * _QE_UNITS["Ry"])
                new_config = False
        elif "temperature" in line:
            m = re.match(r"^\s+temperature\s+=\s+([\d.]+) K", line)
            if m:
                temps.append(float(m.group(1)))
        if fermi is None and "the Fermi energy is" in line:
            m = re.search(r"the Fermi energy is\s+(-?[\d.]+) ev", line)
            if m:
                fermi = float(m.group(1))
        if homo is None and "highest occupied" in line:
            m = re.search(r"highest occupied(?:, lowest unoccupied)? level(?:s)? \(ev\):\s+(-?[\d.]+)", line)
            if m:
                homo = float(m.group(1))
        if not forces_done:
            if forces_block is None and "Forces acting on atoms" in line:
                forces_block = [line]
            elif forces_block is not None:
                forces_block.append(line)
                if "Total force" in line:
                    forces_done = True
    r.energies_ev = energies
    if temps:
        r.temperatures_k = temps
        m = re.search(r"dt\s*=\s*([\d.Ee+-]+)", _text(d / "pw.in"))
        if m:
            fs = float(m.group(1)) * 4.8378e-17 / 1e-15
            r.times_fs = [i * fs for i in range(len(temps))]
            r.frame_dt_fs = fs
    if fermi is not None:
        r.fermi_ev = fermi
    elif homo is not None:
        r.fermi_ev = homo; r.fermi_is_homo = True
    fm = _qe_first_forces_max("".join(forces_block)) if forces_block else None
    if fm is not None:
        r.force_max_ev_ang = fm
        r.force_source = L("output.log の最初の「Forces acting on atoms」", "the first 'Forces acting on atoms' block of output.log")
    dyn = d / "dynmat.out"
    if dyn.is_file():
        rx = re.compile(r"freq\s*\(\s*\d+\)\s*=\s*(-?[\d.]+)\s*\[THz\]\s*=\s*(-?[\d.]+)\s*\[cm-1\]")
        fr = [float(m.group(2)) for m in (rx.search(l) for l in _lines(dyn)) if m]
        if fr:
            r.frequencies_cm1 = fr
            r.notes.append(L("振動数は dynmat.out (dynmat.x が音響和則を当てたあとの Γ 点の値)", "frequencies are from dynmat.out (Gamma point, after dynmat.x applies the acoustic sum rule)"))
    return r


def _qe_first_forces_max(log: str) -> float | None:
    if "Forces acting on atoms" not in log:
        return None
    block = log.split("Forces acting on atoms", 1)[1].split("Total force", 1)[0]
    rows = re.findall(r"atom\s+\d+\s+type\s+\d+\s+force\s+=\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", block)
    if not rows:
        return None
    return float(max(np.linalg.norm([float(x) for x in row]) for row in rows)) * RY_EV / BOHR_ANG


# ---------------- ORCA ----------------
_ORCA_THERMO = (
    ("temperature", r"^\s*Temperature\s+\.\.\.\s+(-?[\d.]+)\s*K", "Temperature"),
    ("zpe", r"^\s*Zero point energy\s+\.\.\.\s+(-?[\d.]+)\s*Eh", "Zero point energy"),
    ("total_thermal_energy", r"^\s*Total thermal energy\s+(-?[\d.]+)\s*Eh", "Total thermal energy"),
    ("total_enthalpy", r"^\s*Total Enthalpy\s+\.\.\.\s+(-?[\d.]+)\s*Eh", "Total Enthalpy"),
    ("final_entropy_term", r"^\s*Final entropy term\s+\.\.\.\s+(-?[\d.]+)\s*Eh", "Final entropy term"),
    ("final_gibbs_free_energy", r"^\s*Final Gibbs free (?:energy|enthalpy)\s+\.\.\.\s+(-?[\d.]+)\s*Eh", "Final Gibbs free energy"),
    ("g_minus_eel", r"^\s*G-E\(el\)\s+\.\.\.\s+(-?[\d.]+)\s*Eh", "G-E(el)"),
)


def _read_orca(d: Path) -> RunData:
    r = RunData("orca", d)
    log = d / "output.log"
    energies, freqs = [], []
    in_vib = in_thermo = in_raman = False
    ramans: list[float] = []
    thermo: dict[str, tuple[float, int, str]] = {}
    rx_e = re.compile(r"FINAL SINGLE POINT ENERGY\s+(-?\d+\.\d+)")
    rx_f = re.compile(r"^\s*\d+:\s+(-?[\d.]+) cm\*\*-1")
    rx_raman = re.compile(r"^\s*\d+:\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$")
    rx_t = [(k, re.compile(p), lab) for k, p, lab in _ORCA_THERMO]
    for no, line in enumerate(_lines(log), start=1):
        m = rx_e.search(line)
        if m:
            energies.append(float(m.group(1)) * HARTREE_EV)
            continue
        if "RAMAN SPECTRUM" in line:
            in_raman, ramans = True, []
            continue
        if in_raman:
            m = rx_raman.match(line)
            if m:
                ramans.append(float(m.group(2)))
            elif ramans and line.strip() and not line.lstrip().startswith(("Mode", "-")):
                in_raman = False
            continue
        if "VIBRATIONAL FREQUENCIES" in line:
            in_vib, freqs = True, []
            continue
        if in_vib:
            if "NORMAL MODES" in line:
                in_vib = False
            else:
                m = rx_f.match(line)
                if m:
                    freqs.append(float(m.group(1)))
            continue
        if "THERMOCHEMISTRY" in line:
            in_thermo = True
            continue
        if in_thermo:
            for k, p, lab in rx_t:
                m = p.search(line)
                if m:
                    thermo[k] = (float(m.group(1)), no, lab)
                    break
    r.energies_ev = energies
    if freqs:
        r.frequencies_cm1 = freqs
    if ramans:
        r.raman_activities = ramans
    inp = _text(d / "orca.inp")
    # GOAT/DOCKER structures are ensembles of alternatives, not an MD trajectory.
    # File names and semantics: ORCA 6.1 official manuals:
    # https://www.faccts.de/docs/orca/6.1/tutorials/prop/goat.html
    # https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/docker.html
    if re.search(r"^!.*\bGOAT\b", inp, re.M | re.I):
        name = "orca.globalminimum.xyz"
        if (d / name).is_file():
            r.frames, r.frame_source = [read(d / name, index=0)], name
        r.energies_ev = []
        r.notes.append(L("GOAT の途中の最適化エネルギーは時系列として描いていません。配座のエネルギーと重みは output.log を確認してください",
                         "intermediate GOAT optimization energies are not plotted as a time series; inspect output.log for conformer energies and weights"))
    elif re.search(r"^%DOCKER\b", inp, re.M | re.I):
        name = "orca.docker.xyz"
        if (d / name).is_file():
            r.frames, r.frame_source = [read(d / name, index=0)], name
        r.energies_ev = []
        r.notes.append(L("DOCKER の配置候補は MD 軌跡ではありません。相互作用エネルギーと候補は output.log と orca.docker.xyz を確認してください",
                         "DOCKER poses are not an MD trajectory; inspect output.log and orca.docker.xyz for interaction energies and poses"))
    else:
        for name in ("trajectory.xyz", "orca.xyz"):
            if (d / name).is_file():
                r.frames, r.frame_source = Trajectory(d / name, "xyz"), name
                break
    ts = re.search(r"Timestep\s+([\d.Ee+-]+)_fs", inp)
    stride = re.search(r"Dump\s+Position\s+Stride\s+(\d+)", inp)
    if ts and r.frame_source == "trajectory.xyz":
        r.frame_dt_fs = float(ts.group(1)) * (int(stride.group(1)) if stride else 1)
    if any(k != "temperature" for k in thermo):
        t = thermo.pop("temperature", None)
        items = [{"key": k, "label": lab, "value_eh": v, "value_ev": v * HARTREE_EV, "line": f"output.log:{no}"}
                 for k, (v, no, lab) in sorted(thermo.items(), key=lambda kv: kv[1][1])]
        r.thermo = {"code": "orca", "method": "ORCA thermochemistry", "source": "output.log", "temperature_k": t[0] if t else None, "items": items}
    return r


EV_ANG3_GPA = 160.21766208  # 1 eV/Å³ [GPa]
_OPT_LOG = re.compile(r"^\w+:\s+(\d+)\s+\d+:\d+:\d+\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s*$")


def _mlip_md_log(p: Path) -> tuple[list[float], list[float], list[float], bool]:
    times, etot, temps, peratom = [], [], [], False
    for line in _lines(p):
        if "Etot" in line:
            peratom = "Etot/N" in line
            continue
        w = line.split()
        if len(w) >= 5:
            try:
                v = [float(x) for x in w[:5]]
            except ValueError:
                continue
            times.append(v[0] * 1000.0); etot.append(v[1]); temps.append(v[4])
    return times, etot, temps, peratom


def _read_mlip(d: Path) -> RunData:
    r = RunData("mlip", d)
    settings, results = {}, {}
    for name, target in (("mlip_settings.json", "settings"), ("results.json", "results")):
        p = d / name
        if p.is_file():
            try:
                obj = json.loads(_text(p))
            except ValueError:
                r.notes.append(L(f"{name} を JSON として読めません", f"cannot parse {name} as JSON"))
                continue
            if isinstance(obj, dict):
                if target == "settings":
                    settings = obj
                else:
                    results = obj
    if not results:
        r.notes.append(L("results.json がありません (run_mlip.py がまだ終わっていないかもしれません)",
                         "results.json is missing (run_mlip.py may not have finished)"))
    task = str(results.get("task") or settings.get("task") or "")
    for name, kind in (("trajectory.extxyz", "traj"), ("final.extxyz", "final"), ("structure.extxyz", "start")):
        p = d / name
        if p.is_file() and p.stat().st_size > 0:
            r.frames = Trajectory(p, "ase", fmt="extxyz") if kind == "traj" else [read(p, format="extxyz")]
            r.frame_source = name
            break
    md = d / "md.log"
    if task == "molecular_dynamics" and md.is_file():
        times, etot, temps, peratom = _mlip_md_log(md)
        if etot:
            r.energies_ev, r.temperatures_k, r.times_fs = etot, temps, times
            s = settings.get("md") or {}
            try:
                r.frame_dt_fs = float(s["timestep_fs"]) * int(s["dump_interval"])
            except (KeyError, TypeError, ValueError):
                r.frame_dt_fs = None
            r.notes.append(L("エネルギーと温度は md.log (ASE の MDLogger) の Etot と T。全エネルギー = ポテンシャル + 運動エネルギー",
                             "energies and temperatures are Etot and T of md.log (ASE MDLogger); the total energy is potential + kinetic"))
            if peratom:
                r.notes.append(L("md.log の値は原子 1 個あたりです (Etot/N)。そのまま並べています", "the md.log values are per atom (Etot/N) and are listed as they are"))
    opt = d / "opt.log"
    if not r.energies_ev and opt.is_file():
        rows = [m.groups() for m in (_OPT_LOG.match(l) for l in _lines(opt)) if m]
        if rows:
            r.energies_ev = [float(x[1]) for x in rows]
            r.notes.append(L("エネルギーの推移は opt.log (ASE の最適化の記録) の Energy の列 (格子も動かす設定では、その目的関数の値)",
                             "the energy trace is the Energy column of opt.log (ASE optimizer log; with cell relaxation it is the value of that objective)"))
    if not r.energies_ev and results.get("energy_ev") is not None:
        r.energies_ev = [float(results["energy_ev"])]
    if results.get("frequencies_cm1"):
        r.frequencies_cm1 = [float(x) for x in results["frequencies_cm1"]]
        r.notes.append(L("振動数は results.json (ASE の Vibrations。虚数は負の数で書かれています)", "frequencies are from results.json (ASE Vibrations; imaginary ones are written as negative)"))
        if results.get("fmax_ev_per_ang") is not None:
            r.force_max_ev_ang = float(results["fmax_ev_per_ang"])
            r.force_source = L("results.json の fmax_ev_per_ang (振動解析の基準の構造)", "fmax_ev_per_ang in results.json (the reference structure of the vibrational analysis)")
    t = {"task": task, "model_family": results.get("model_family") or settings.get("model_family"),
         "model": results.get("model") or settings.get("model"), "versions": results.get("versions"),
         "converged": results.get("converged"), "steps": results.get("steps"),
         "energy_ev": results.get("energy_ev"), "fmax_ev_per_ang": results.get("fmax_ev_per_ang"),
         "stress_voigt_ev_per_ang3": results.get("stress_voigt_ev_per_ang3"), "stress_note": results.get("stress_note"),
         "source": "results.json"}
    sv = results.get("stress_voigt_ev_per_ang3")
    if isinstance(sv, list) and len(sv) == 6:
        t["stress_voigt_gpa"] = [float(x) * EV_ANG3_GPA for x in sv]
        t["pressure_gpa"] = float(-np.mean([float(x) for x in sv[:3]]) * EV_ANG3_GPA)
        t["stress_convention"] = L("応力は ASE の約束 (圧縮で負)。圧力 = -(σxx + σyy + σzz) / 3",
                                   "stress follows the ASE sign convention (compression negative); pressure = -(sxx + syy + szz) / 3")
    r.extra_tables["mlip"] = t
    return r
