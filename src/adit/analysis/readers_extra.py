
from __future__ import annotations

import re
import shutil
from pathlib import Path

import numpy as np
from ase.io import read

from adit.analysis.readers import BOHR_ANG, HARTREE_EV, RunData, _grep, _lines, _text
from adit.analysis.trajectory import Trajectory
from adit.lang import L

KCAL_MOL_EV = 1.0 / 23.060547830619
KJ_MOL_EV = 1.0 / 96.485332123
ATM_BAR = 1.01325


def _spec(d: Path):
    try:
        from adit.project import load_project
        return load_project(d)
    except Exception:
        return None


# ---------------- LAMMPS ----------------
def _lammps_thermo(log: Path) -> tuple[list[dict], int | None]:
    rows: list[dict] = []
    cols = None
    first_line = None
    for no, line in enumerate(_lines(log), start=1):
        w = line.split()
        if w and w[0] == "Step":
            cols, first_line = w, first_line or no
            continue
        if cols is None:
            continue
        if line.startswith("Loop time"):
            cols = None
            continue
        if len(w) == len(cols):
            try:
                rows.append(dict(zip(cols, (float(x) for x in w))))
            except ValueError:
                pass
    return rows, first_line


def read_lammps(d: Path) -> RunData:
    r = RunData("lammps", d)
    inp = _text(d / "in.lammps")
    m = re.search(r"^\s*units\s+(\w+)", inp, re.M)
    units = m.group(1) if m else None
    t_fs = {"metal": 1000.0, "real": 1.0}.get(units)
    e_ev = {"metal": 1.0, "real": KCAL_MOL_EV}.get(units)
    p_bar = {"metal": 1.0, "real": ATM_BAR}.get(units)
    if units and t_fs is None:
        r.notes.append(L(f"単位系 {units} の換算を持っていないので、時間・エネルギー・圧力は読みません (metal と real だけ)",
                         f"no conversion for units {units}; time, energy and pressure are not read (metal and real only)"))
    is_md = re.search(r"^\s*run\s+[1-9]", inp, re.M) is not None and re.search(r"^\s*minimize\b", inp, re.M) is None
    rows, line = _lammps_thermo(d / "log.lammps")
    src = f"log.lammps:{line}" if line else "log.lammps"
    if rows and e_ev is not None:
        key = "TotEng" if "TotEng" in rows[0] and is_md else "PotEng"
        if key in rows[0]:
            r.energies_ev = [x[key] * e_ev for x in rows]
        ts = re.search(r"^\s*timestep\s+([\d.eE+-]+)", inp, re.M)
        if "Time" in rows[0]:
            r.times_fs = [x["Time"] * t_fs for x in rows]
        elif ts and is_md:
            r.times_fs = [x["Step"] * float(ts.group(1)) * t_fs for x in rows]
        if is_md and "Temp" in rows[0]:
            r.temperatures_k = [x["Temp"] for x in rows]
        if "Press" in rows[0]:
            r.series["pressure"] = {"values": [x["Press"] * p_bar for x in rows], "unit": "bar", "source": f"{src} (Press)"}
        if all(k in rows[0] for k in ("Pxy", "Pxz", "Pyz")):
            r.series["pressure_tensor"] = {"values": {k: [x[k] for x in rows] for k in ("Pxy", "Pxz", "Pyz")},
                                           "unit": "bar" if units == "metal" else "atm",
                                           "source": f"{src} (Pxy, Pxz, Pyz)"}
        r.notes.append(L(f"エネルギーは {src} の thermo の表の {key} ({units} 単位から eV に換算)",
                         f"energy is {key} of the thermo table in {src} (converted from {units} units to eV)"))
    dump = re.search(r"^\s*dump\s+\S+\s+\S+\s+custom\s+(\d+)\s+(\S+)", inp, re.M)
    els = re.search(r"^\s*dump_modify\s+\S+\s+element\s+(.+?)(?:\s+sort\s+\S+)?\s*$", inp, re.M)
    types = els.group(1).split() if els else []
    for name in ([dump.group(2)] if dump else []) + ["traj.lammpstrj", "forces.lammpstrj"]:
        if (d / name).is_file():
            r.frames, r.frame_source = Trajectory(d / name, "lammpsdump", fmt=" ".join(types)), name
            break
    else:
        if (d / "final.data").is_file():
            try:
                style = re.search(r"^\s*atom_style\s+(\w+)", inp, re.M)
                from ase.data import atomic_numbers
                zmap = {i + 1: atomic_numbers[e] for i, e in enumerate(types)} if types else None
                r.frames = [read(d / "final.data", format="lammps-data", atom_style=style.group(1) if style else "atomic",
                                 units=units or "metal", Z_of_type=zmap)]
                r.frame_source = "final.data"
            except Exception as ex:
                r.notes.append(L(f"final.data を読めません: {ex}", f"cannot read final.data: {ex}"))
    ts = re.search(r"^\s*timestep\s+([\d.eE+-]+)", inp, re.M)
    if is_md and dump and ts and t_fs is not None and r.frame_source == dump.group(2):
        r.frame_dt_fs = float(ts.group(1)) * t_fs * int(dump.group(1))
    attach_plumed(r, d)
    return r


# ---------------- GROMACS ----------------
def output_prefix(d: Path, default: str = "adit") -> str:
    text = _text(d / "submit.sh") + _text(d / "namd.conf")
    m = re.search(r"-deffnm\s+(\S+)", text) or re.search(r"^\s*outputName\s+(\S+)", text, re.M)
    if m:
        return m.group(1)
    from adit.config import OLD_APP_NAMES

    for name in (default, *OLD_APP_NAMES):
        if any(d.glob(f"{name}.*")):
            return name
    return default


def _gromacs_blocks(log: Path) -> list[dict]:
    blocks: list[dict] = []
    state = None
    names: list[str] | None = None
    cur: dict | None = None
    for no, line in enumerate(_lines(log), start=1):
        if "A V E R A G E S" in line:
            break
        s = line.split()
        if s == ["Step", "Time"]:
            state, cur = "step", None
            continue
        if state == "step" and len(s) == 2:
            try:
                cur = {"step": int(s[0]), "time_ps": float(s[1]), "values": {}, "line": no}
            except ValueError:
                cur = None
            state = "wait"
            continue
        if "Energies (kJ/mol)" in line:
            if cur is None:
                cur = {"step": None, "time_ps": None, "values": {}, "line": no}
            state, names = "names", None
            continue
        if state == "names":
            if not line.strip():
                if cur is not None and cur["values"]:
                    blocks.append(cur)
                state, cur = None, None
                continue
            if names is None:
                names = [line[i:i + 15].strip() for i in range(0, len(line.rstrip("\n")), 15)]
                names = [n for n in names if n]
            else:
                try:
                    vals = [float(x) for x in s]
                except ValueError:
                    names = None
                    continue
                cur["values"].update(zip(names, vals))
                names = None
    return blocks


def read_gromacs(d: Path) -> RunData:
    r = RunData("gromacs", d)
    mdp = _text(d / "grompp.mdp")
    m = re.search(r"^\s*integrator\s*=\s*(\S+)", mdp, re.M)
    integrator = m.group(1).lower() if m else ""
    is_md = integrator in ("md", "md-vv", "md-vv-avek", "sd", "bd")
    prefix = output_prefix(d)
    log = d / f"{prefix}.log"
    if not log.is_file():
        logs = [p for p in d.glob("*.log") if p.name not in ("grompp.log", "output.log")]
        log = logs[0] if logs else log
    blocks = _gromacs_blocks(log)
    if blocks:
        v0 = blocks[0]["values"]
        key = "Total Energy" if is_md and "Total Energy" in v0 else "Potential"
        r.energies_ev = [b["values"][key] * KJ_MOL_EV for b in blocks if key in b["values"]]
        r.notes.append(L(f"エネルギーは {log.name} の「{key}」(系全体の kJ/mol を 96.485 kJ/mol/eV で割って、系全体の eV に換算)",
                         f"energy is '{key}' in {log.name} (whole-system kJ/mol divided by 96.485 kJ/mol per eV to give eV for the whole system)"))
        if is_md and all(b["time_ps"] is not None for b in blocks):
            r.times_fs = [b["time_ps"] * 1000.0 for b in blocks]
        if is_md and "Temperature" in v0:
            r.temperatures_k = [b["values"].get("Temperature", np.nan) for b in blocks]
        if "Pressure (bar)" in v0:
            r.series["pressure"] = {"values": [b["values"].get("Pressure (bar)", np.nan) for b in blocks], "unit": "bar",
                                    "source": f"{log.name} (Pressure (bar))"}
        if "Density" in v0:
            r.series["density_log"] = {"values": [b["values"].get("Density", np.nan) / 1000.0 for b in blocks], "unit": "g/cm^3",
                                       "source": f"{log.name} (Density、kg/m^3 を g/cm^3 に)"}
        if is_md and "Conserved En." in v0:
            r.series["conserved"] = {"values": [b["values"]["Conserved En."] * KJ_MOL_EV for b in blocks], "unit": "eV",
                                     "source": f"{log.name} (Conserved En.)"}
    for name in (f"{prefix}.gro", "confout.gro"):
        if (d / name).is_file():
            try:
                r.frames, r.frame_source = [read(d / name, format="gromacs")], name
                break
            except Exception as ex:
                r.notes.append(L(f"{name} を読めません: {ex}", f"cannot read {name}: {ex}"))
    for name in (f"{prefix}.xtc", f"{prefix}.trr"):
        if not (d / name).is_file():
            continue
        from adit.analysis import trajectory_ext as ext

        try:
            frames, lib = ext.read_external(d / name, ext.find_topology(d, d / name))
        except ext.ExternalTrajectoryError as ex:
            r.notes.append(L(f"軌跡 {name} を読めません: {ex} "
                             f"RDF・MSD・書き出しは、いまは最後の構造 {r.frame_source or '-'} だけを使っています",
                             f"the trajectory {name} cannot be read: {ex} "
                             f"RDF, MSD and export currently use only the final structure {r.frame_source or '-'}"))
            break
        r.frames, r.frame_source = frames, name
        r.notes.append(L(f"軌跡 {name} を {lib} で読みました ({len(frames)} フレーム)。"
                         "座標は GROMACS の gmx trjconv の出力と 2.3e-06 Å 以内で一致することを確かめてあります",
                         f"the trajectory {name} was read with {lib} ({len(frames)} frames); the coordinates agree with "
                         "gmx trjconv output to within 2.3e-06 A"))
        if r.times_fs and len(r.times_fs) >= 2 and len(frames) >= 2:
            span = r.times_fs[-1] - r.times_fs[0]
            r.frame_dt_fs = span / (len(frames) - 1) if span > 0 else r.frame_dt_fs
        break
    if shutil.which("gmx"):
        r.notes.append(L("gmx energy -f adit.edr で、温度・圧力・密度などの全項目の時系列を取り出せます (エネルギーの記録 .edr は GROMACS のバイナリ形式)",
                         "gmx energy -f adit.edr extracts the time series of temperature, pressure, density and all other terms (.edr is GROMACS binary)"))
    attach_plumed(r, d)
    return r


# ---------------- CP2K ----------------
def _cp2k_charge_block(lines: list[str], k: int, total_prefix: str) -> list[float]:
    header = next((lines[j] for j in range(k + 1, min(k + 4, len(lines))) if "Net charge" in lines[j]), "")
    idx = -2 if "Spin moment" in header else -1
    vals = []
    started = False
    for l in lines[k + 1:]:
        s = l.split()
        if l.strip().startswith(total_prefix):
            break
        if not s or s[0].startswith("#") or "Net charge" in l:
            if started and not s:
                break
            continue
        try:
            int(s[0]); vals.append(float(s[idx])); started = True
        except (ValueError, IndexError):
            if started:
                break
    return vals


def read_cp2k(d: Path) -> RunData:
    r = RunData("cp2k", d)
    inp = _text(d / "cp2k.inp")
    m = re.search(r"^\s*PROJECT\s+(\S+)", inp, re.M | re.I)
    project = m.group(1) if m else output_prefix(d)
    m = re.search(r"^\s*RUN_TYPE\s+(\S+)", inp, re.M | re.I)
    run_type = m.group(1).upper() if m else "ENERGY_FORCE"
    cell_block = re.search(r"&CELL\b(.*?)&END CELL", inp, re.S | re.I)
    cell = pbc = None
    if cell_block:
        vecs = [re.search(rf"^\s*{c}\s+(\S+)\s+(\S+)\s+(\S+)", cell_block.group(1), re.M) for c in "ABC"]
        per = re.search(r"^\s*PERIODIC\s+(\S+)", cell_block.group(1), re.M | re.I)
        per = per.group(1).upper() if per else "XYZ"
        if all(vecs) and per != "NONE":
            cell = np.array([[float(x) for x in v.groups()] for v in vecs])
            pbc = [c in per for c in "XYZ"]
    log = d / "output.log"
    energies = []
    lines = [] if not log.is_file() else (list(_lines(log)) if log.stat().st_size < 64 * 2**20 else None)
    if lines is None:
        energies = [float(l.split()[-1]) * HARTREE_EV for l in _lines(log) if "ENERGY| Total FORCE_EVAL" in l]
        lines = []
    else:
        energies = [float(l.split()[-1]) * HARTREE_EV for l in lines if "ENERGY| Total FORCE_EVAL" in l]
    r.energies_ev = energies
    for title, prefix, definition in (("Mulliken Population Analysis", "# Total charge", "Mulliken (CP2K)"),
                                      ("Hirshfeld Charges", "Total Charge", "Hirshfeld (CP2K)")):
        ks = [k for k, l in enumerate(lines) if l.strip() == title]
        if ks:
            vals = _cp2k_charge_block(lines, ks[-1], prefix)
            if vals:
                r.charges.append({"definition": definition, "values": vals, "source": f"output.log:{ks[-1] + 1}", "unit": "e",
                                  "note": L("CP2K の output.log の最後の表の Net charge", "Net charge of the last table in CP2K output.log")})
    ener = d / f"{project}-1.ener"
    if run_type == "MD" and ener.is_file():
        rows = []
        for l in _lines(ener):
            s = l.split()
            if s and not s[0].startswith("#") and len(s) >= 6:
                try:
                    rows.append([float(x) for x in s[:6]])
                except ValueError:
                    pass
        if rows:
            a = np.array(rows)
            r.times_fs = a[:, 1].tolist()
            r.temperatures_k = a[:, 3].tolist()
            r.energies_ev = ((a[:, 2] + a[:, 4]) * HARTREE_EV).tolist()
            r.series["conserved"] = {"values": (a[:, 5] * HARTREE_EV).tolist(), "unit": "eV", "source": f"{ener.name} (Cons Qty)"}
            r.notes.append(L(f"MD のエネルギーは {ener.name} の運動エネルギー + ポテンシャルエネルギー (hartree を eV に)。保存量 (Cons Qty) は別の系列",
                             f"MD energy is kinetic + potential energy from {ener.name} (hartree to eV); the conserved quantity (Cons Qty) is a separate series"))
    pos = d / f"{project}-pos-1.xyz"
    if pos.is_file():
        r.frames, r.frame_source = Trajectory(pos, "xyz", cell=cell, pbc=pbc), pos.name
        if (d / f"{project}-1.cell").is_file():
            r.notes.append(L(f"セルが変わる計算です ({project}-1.cell)。軌跡には入力 (cp2k.inp) のセルを付けています",
                             f"the cell changes during this run ({project}-1.cell); the trajectory carries the input cell from cp2k.inp"))
    else:
        ks = [k for k, l in enumerate(lines) if "ATOMIC COORDINATES IN ANGSTROM" in l]
        if ks:
            syms, p = [], []
            for l in lines[ks[-1] + 3:]:
                s = l.split()
                if len(s) < 7:
                    break
                syms.append(s[2]); p.append([float(x) for x in s[4:7]])
            from ase import Atoms
            a = Atoms(syms, positions=p)
            if cell is not None:
                a.set_cell(cell); a.pbc = pbc
            r.frames, r.frame_source = [a], "output.log"
    ts = re.search(r"^\s*TIMESTEP\s+([\d.eE+-]+)", inp, re.M | re.I)
    each = re.search(r"&TRAJECTORY.*?&EACH\s+MD\s+(\d+)", inp, re.S | re.I)
    if run_type == "MD" and ts and r.frame_source == pos.name:
        r.frame_dt_fs = float(ts.group(1)) * (int(each.group(1)) if each else 1)
    mol = d / f"{project}-VIBRATIONS-1.mol"
    if mol.is_file():
        fr, it, sec = [], [], None
        for l in _lines(mol):
            t = l.strip()
            if t.startswith("["):
                sec = t.upper()
                continue
            if sec == "[FREQ]" and t:
                fr.append(float(t.split()[0]))
            elif sec == "[INT]" and t:
                it.append(float(t.split()[0]))
        if fr:
            r.frequencies_cm1 = fr
            r.ir_intensities = it if len(it) == len(fr) else None
            r.notes.append(L(f"振動数は {mol.name} (Molden 形式) の [FREQ]", f"frequencies are [FREQ] of {mol.name} (Molden format)"))
    return r


# ---------------- OpenMM ----------------
def _openmm_md_log(path: Path) -> tuple[list[str], list[list[float]]]:
    header: list[str] = []
    rows: list[list[float]] = []
    for line in _lines(path):
        text = line.strip()
        if not text:
            continue
        if not header:
            if text.startswith("#"):
                header = [c.strip().strip('"') for c in text.lstrip("#").split(",")]
            continue
        cells = text.split(",")
        if len(cells) != len(header):
            continue
        try:
            rows.append([float(c) for c in cells])
        except ValueError:
            continue
    return header, rows


def read_openmm(d: Path) -> RunData:
    import json

    r = RunData("openmm", d)
    settings, results = {}, {}
    for name, is_settings in (("openmm_settings.json", True), ("results.json", False)):
        p = d / name
        if not p.is_file():
            continue
        try:
            obj = json.loads(_text(p))
        except ValueError:
            r.notes.append(L(f"{name} を JSON として読めません", f"cannot parse {name} as JSON"))
            continue
        if isinstance(obj, dict):
            if is_settings:
                settings = obj
            else:
                results = obj
    if not results:
        r.notes.append(L("results.json がありません (run_openmm.py がまだ終わっていないかもしれません)",
                         "results.json is missing (run_openmm.py may not have finished)"))
    task = str(results.get("task") or settings.get("task") or "")
    md_log = d / "md.log"
    if task == "molecular_dynamics" and md_log.is_file():
        header, rows = _openmm_md_log(md_log)
        col = {name: i for i, name in enumerate(header)}

        def column(name: str) -> list[float] | None:
            i = col.get(name)
            return [row[i] for row in rows] if i is not None and rows else None

        total = column("Total Energy (kJ/mole)")
        if total:
            r.energies_ev = [x * KJ_MOL_EV for x in total]
            r.notes.append(L("エネルギーは md.log (OpenMM の StateDataReporter) の Total Energy (系全体の kJ/mol を 96.485 kJ/mol/eV で割って、系全体の eV に換算)",
                             "energy is Total Energy in md.log (OpenMM StateDataReporter); whole-system kJ/mol divided by 96.485 kJ/mol per eV to give eV for the whole system"))
            times_ps = column("Time (ps)")
            if times_ps:
                r.times_fs = [x * 1000.0 for x in times_ps]
            temps = column("Temperature (K)")
            if temps:
                r.temperatures_k = temps
            volume = column("Box Volume (nm^3)")
            if volume:
                r.series["volume"] = {"values": [x * 1000.0 for x in volume], "unit": "Å³",
                                      "source": "md.log (Box Volume (nm^3) を Å³ に)"}
            density = column("Density (g/mL)")
            if density:
                r.series["density_log"] = {"values": density, "unit": "g/cm^3", "source": "md.log (Density (g/mL))"}
            s = settings.get("md") or {}
            try:
                r.frame_dt_fs = float(s["timestep_fs"]) * int(s["dump_interval"])
            except (KeyError, TypeError, ValueError):
                r.frame_dt_fs = None
    if not r.energies_ev:
        for key in ("total_energy_ev", "potential_energy_ev"):
            if isinstance(results.get(key), (int, float)):
                r.energies_ev = [float(results[key])]
                r.notes.append(L(f"エネルギーは results.json の {key} の 1 点です", f"the energy is the single {key} value of results.json"))
                break
    traj = d / "trajectory.extxyz"
    final = d / "final.pdb"
    if traj.is_file() and traj.stat().st_size > 0:
        try:
            r.frames, r.frame_source = Trajectory(traj, "ase", fmt="extxyz"), traj.name
        except Exception as ex:
            r.notes.append(L(f"trajectory.extxyz を読めません: {ex}", f"cannot read trajectory.extxyz: {ex}"))
    if not r.frames and final.is_file():
        try:
            r.frames, r.frame_source = [read(final, format="proteindatabank")], final.name
        except Exception as ex:
            r.notes.append(L(f"final.pdb を読めません: {ex}", f"cannot read final.pdb: {ex}"))
    if (d / "trajectory.dcd").is_file() and r.frame_source != "trajectory.extxyz":
        from adit.analysis import trajectory_ext as ext

        try:
            frames, lib = ext.read_external(d / "trajectory.dcd", ext.find_topology(d, d / "trajectory.dcd"))
        except ext.ExternalTrajectoryError as ex:
            r.notes.append(L(f"軌跡 trajectory.dcd を読めません: {ex} "
                             "RDF・MSD・書き出しは最後の構造だけを使います "
                             "(openmm_settings.json の write_xyz_trajectory を true にして実行し直すと、extxyz でも書けます)",
                             f"the trajectory trajectory.dcd cannot be read: {ex} "
                             "RDF, MSD and export use only the final structure "
                             "(set write_xyz_trajectory to true in openmm_settings.json and rerun to get an extxyz as well)"))
        else:
            r.frames, r.frame_source = frames, "trajectory.dcd"
            r.notes.append(L(f"軌跡 trajectory.dcd を {lib} で読みました ({len(frames)} フレーム)",
                             f"the trajectory trajectory.dcd was read with {lib} ({len(frames)} frames)"))
    r.extra_tables["openmm"] = {"task": task, "input_format": results.get("input_format") or settings.get("input_format"),
                                "openmm_version": results.get("openmm_version"), "platform": results.get("platform"),
                                "ensemble": results.get("ensemble"), "steps": results.get("steps"),
                                "nonbonded_method": settings.get("nonbonded_method"),
                                "nonbonded_cutoff_nm": settings.get("nonbonded_cutoff_nm"),
                                "constraints": settings.get("constraints"),
                                "potential_energy_kj_per_mol": results.get("potential_energy_kj_per_mol"),
                                "total_energy_kj_per_mol": results.get("total_energy_kj_per_mol"),
                                "final_temperature_k": results.get("final_temperature_k"),
                                "degrees_of_freedom": results.get("degrees_of_freedom"),
                                "source": "results.json + openmm_settings.json"}
    r.notes.append(L("力場・原子型・電荷は外部のトポロジーが決めています。ADIT は照合していません。",
                     "the force field, atom types and charges come from the external topology; ADIT did not check them."))
    attach_plumed(r, d)
    return r


# ---------------- ABINIT ----------------
_ABINIT_STEP_ENERGY = re.compile(r"Total energy \(etotal\) \[Ha\]=\s*([-+0-9.EDed]+)")


def _abinit_outvars(path: Path) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    current: str | None = None
    started = False
    for line in _lines(path):
        if "-outvars: echo values of variables after computation" in line:
            started, out, current = True, {}, None
            continue
        if not started or not line.strip():
            continue
        words = line.split()
        if words[0] in ("-", "P"):
            words = words[1:]
        elif len(words[0]) > 1 and words[0][0] in "-P" and words[0][1].isalpha():
            words[0] = words[0][1:]
        if not words:
            continue
        if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", words[0]):
            current = words[0]
            out[current] = words[1:]
        elif current is not None:
            out[current] += words
    return out


def _abinit_floats(values: list[str]) -> list[float]:
    picked = []
    for v in values:
        try:
            picked.append(float(v.replace("D", "E").replace("d", "e")))
        except ValueError:
            break
    return picked


def read_abinit(d: Path) -> RunData:
    from ase import Atoms
    from ase.data import chemical_symbols

    r = RunData("abinit", d)
    log, abo = d / "output.log", d / "input.abo"
    energies = [float(x.replace("D", "E")) * HARTREE_EV for x in _grep(log, _ABINIT_STEP_ENERGY.pattern)]
    outvars = _abinit_outvars(abo) if abo.is_file() else {}
    if energies:
        r.energies_ev = energies
        r.notes.append(L("エネルギーは output.log の各段階の「Total energy (etotal) [Ha]」(eV に換算)",
                         "energies are 'Total energy (etotal) [Ha]' of each step in output.log (converted to eV)"))
    elif "etotal" in outvars:
        values = _abinit_floats(outvars["etotal"])
        if values:
            r.energies_ev = [values[0] * HARTREE_EV]
            r.notes.append(L("エネルギーは input.abo の最後の etotal (eV に換算)",
                             "the energy is the final etotal in input.abo (converted to eV)"))
    try:
        acell = _abinit_floats(outvars["acell"])
        rprim = _abinit_floats(outvars["rprim"])
        xred = _abinit_floats(outvars["xred"])
        znucl = _abinit_floats(outvars["znucl"])
        typat = [int(float(x)) for x in outvars["typat"]]
        cell = np.array(rprim[:9], dtype=float).reshape(3, 3) * np.array(acell[:3], dtype=float)[:, None] * BOHR_ANG
        scaled = np.array(xred[:3 * len(typat)], dtype=float).reshape(len(typat), 3)
        symbols = [chemical_symbols[int(round(znucl[i - 1]))] for i in typat]
        r.frames = [Atoms(symbols=symbols, scaled_positions=scaled, cell=cell, pbc=True)]
        r.frame_source = abo.name
        r.notes.append(L("最後の構造は input.abo の末尾の変数の表 (acell・rprim・xred) から組み立てました",
                         "the final structure was built from the variable table at the end of input.abo (acell, rprim, xred)"))
    except (KeyError, ValueError, IndexError):
        if abo.is_file():
            r.notes.append(L("input.abo から最後の構造を読めませんでした", "could not read the final structure from input.abo"))
    if _text(log).count("gradients are converged"):
        r.notes.append(L("構造最適化は ABINIT が「gradients are converged」と書いています",
                         "ABINIT reported 'gradients are converged' for the structural optimization"))
    return r


# ---------------- Psi4 ----------------
def _psi4_optimization_steps(log: Path) -> list[float]:
    energies: list[float] = []
    inside = False
    for line in _lines(log):
        if "==> Optimization Summary <==" in line:
            inside, energies = True, []
            continue
        if not inside:
            continue
        words = line.replace("~", "").split()
        if len(words) >= 3 and words[0].isdigit():
            try:
                energies.append(float(words[1]))
            except ValueError:
                continue
        elif energies and not line.strip():
            inside = False
    return energies


def read_psi4(d: Path) -> RunData:
    import json

    r = RunData("psi4", d)
    results = {}
    p = d / "results.json"
    if p.is_file():
        try:
            obj = json.loads(_text(p))
            results = obj if isinstance(obj, dict) else {}
        except ValueError:
            r.notes.append(L("results.json を JSON として読めません", "cannot parse results.json as JSON"))
    if not results:
        r.notes.append(L("results.json がありません (Psi4 がまだ終わっていないかもしれません)",
                         "results.json is missing (Psi4 may not have finished)"))
    steps = _psi4_optimization_steps(d / "output.log")
    if steps:
        r.energies_ev = [x * HARTREE_EV for x in steps]
        r.notes.append(L("エネルギーの推移は output.log の「Optimization Summary」の表 (各段階の Total Energy [Hartree] を eV に換算)",
                         "the energy trace is the 'Optimization Summary' table in output.log (Total Energy [hartree] of each step, converted to eV)"))
    elif isinstance(results.get("energy_hartree"), (int, float)):
        r.energies_ev = [float(results["energy_hartree"]) * HARTREE_EV]
        r.notes.append(L("エネルギーは results.json の energy_hartree (Psi4 の返り値。eV に換算)",
                         "the energy is energy_hartree in results.json (the Psi4 return value, converted to eV)"))
    if results.get("frequencies_cm1"):
        r.frequencies_cm1 = [float(x) for x in results["frequencies_cm1"]]
        r.notes.append(L("振動数は results.json (Psi4 の wfn.frequencies()。虚数は負の数)",
                         "frequencies are from results.json (Psi4's wfn.frequencies(); imaginary ones are negative)"))
    xyz = d / "final.xyz"
    if xyz.is_file():
        try:
            r.frames, r.frame_source = [read(xyz, format="xyz")], xyz.name
        except Exception as ex:
            r.notes.append(L(f"final.xyz を読めません: {ex}", f"cannot read final.xyz: {ex}"))
    r.extra_tables["psi4"] = {"task": results.get("task"), "method": results.get("method"), "basis": results.get("basis"),
                              "psi4_version": results.get("psi4_version"), "energy_hartree": results.get("energy_hartree"),
                              "geom_maxiter": results.get("geom_maxiter"), "converged": results.get("converged"),
                              "optimization_steps": len(steps) or None, "source": "results.json"}
    return r


def plumed_colvar_tables(d: Path) -> dict:
    tables: dict[str, dict] = {}
    for path in sorted(d.iterdir()):
        if not path.is_file() or path.name in ("plumed.dat", "plumed.log"):
            continue
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                head = f.readline()
        except OSError:
            continue
        if not head.startswith("#! FIELDS"):
            continue
        fields = head.split()[2:]
        sums = [0.0] * len(fields)
        first: list[float] | None = None
        last: list[float] | None = None
        rows = 0
        for line in _lines(path):
            if line.startswith("#") or not line.strip():
                continue
            words = line.split()
            if len(words) != len(fields):
                continue
            try:
                values = [float(x) for x in words]
            except ValueError:
                continue
            rows += 1
            sums = [a + b for a, b in zip(sums, values)]
            last = values
            if first is None:
                first = values
        if rows:
            tables[path.name] = {"fields": fields, "rows": rows,
                                 "first": dict(zip(fields, first or [])),
                                 "last": dict(zip(fields, last or [])),
                                 "mean": {k: v / rows for k, v in zip(fields, sums)}}
    return tables


def attach_plumed(r: RunData, d: Path) -> None:
    tables = plumed_colvar_tables(d)
    if not tables:
        return
    r.extra_tables["plumed"] = {"files": tables,
                                "source": L("PLUMED が書いた COLVAR 形式のファイル (#! FIELDS の列)",
                                            "files in PLUMED COLVAR format (columns of the #! FIELDS line)")}
    names = ", ".join(sorted(tables))
    r.notes.append(L(f"PLUMED の出力を読みました ({names})。列の意味は利用者の plumed.dat が決めるので、ADIT は値を並べるだけです。",
                     f"PLUMED output was read ({names}); the meaning of each column comes from your plumed.dat, so ADIT only lists the values."))
