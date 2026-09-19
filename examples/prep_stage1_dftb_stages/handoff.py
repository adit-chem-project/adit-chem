"""Carry the final structure, and velocities when continuing MD, from one run into the next. Standard library only, because this file is copied next to the run."""

from __future__ import annotations

import os
import re
import shutil
import sys

BOHR_ANG = 0.529177210903
CONVERTS = ("dftbplus", "xtb", "espresso", "orca")
VELOCITY_FILE = "velocities.dat"


PREFIXES = ("adit", "vista", "qcgui")


def existing_name(directory: str, pattern: str) -> str:
    for prefix in PREFIXES:
        name = pattern.format(prefix)
        if directory and os.path.isfile(os.path.join(directory, name)):
            return name
    return pattern.format(PREFIXES[0])


class HandoffError(Exception):
    pass


def copy_files(code: str, previous_task: str, velocities: bool, gromacs_conf: str = "conf.gro",
               previous_dir: str = "") -> dict[str, str]:
    md = previous_task == "molecular_dynamics"
    if code == "vasp" and previous_task in ("geometry_optimization", "molecular_dynamics"):
        return {"POSCAR": "CONTCAR"}
    if code == "cp2k" and previous_task in ("geometry_optimization", "molecular_dynamics"):
        return {"prev.restart": existing_name(previous_dir, "{}-1.restart")}
    if code == "lammps" and previous_task in ("geometry_optimization", "molecular_dynamics"):
        return {"data.lammps": "final.data"}
    if code == "gromacs" and previous_task in ("geometry_optimization", "molecular_dynamics", "single_point"):
        out = {gromacs_conf: existing_name(previous_dir, "{}.gro")}
        if velocities and md:
            out["prev.cpt"] = existing_name(previous_dir, "{}.cpt")
        return out
    if code == "xtb" and velocities and md:
        return {"mdrestart": "mdrestart"}
    return {}


def _floats(tokens) -> list[float]:
    return [float(t.replace("D", "E").replace("d", "e")) for t in tokens]


def _mat_vec(cell, frac):
    return [sum(frac[j] * cell[j][i] for j in range(3)) for i in range(3)]


def parse_gen(text: str) -> dict:
    lines = [l.split("#", 1)[0].strip() for l in text.splitlines()]
    lines = [l for l in lines if l]
    head = lines[0].split()
    n, kind = int(head[0]), head[1].upper()
    types = lines[1].split()
    symbols, pos = [], []
    for l in lines[2:2 + n]:
        p = l.split()
        symbols.append(types[int(p[1]) - 1])
        pos.append(_floats(p[2:5]))
    cell = None
    if kind in ("S", "F"):
        cell = [_floats(lines[2 + n + 1 + i].split()[:3]) for i in range(3)]
        if kind == "F":
            pos = [_mat_vec(cell, f) for f in pos]
    return {"symbols": symbols, "positions": pos, "cell": cell, "velocities": None}


def format_gen(symbols, positions, cell) -> str:
    types = []
    for s in symbols:
        if s not in types:
            types.append(s)
    lines = [f"{len(symbols)} {'S' if cell else 'C'}", " ".join(types)]
    for i, (s, p) in enumerate(zip(symbols, positions), 1):
        lines.append(f"{i} {types.index(s) + 1} {p[0]:.12f} {p[1]:.12f} {p[2]:.12f}")
    if cell:
        lines.append("0.0 0.0 0.0")
        lines += [f"{v[0]:.12f} {v[1]:.12f} {v[2]:.12f}" for v in cell]
    return "\n".join(lines) + "\n"


def last_xyz_frame(path: str) -> list[list[str]]:
    last = None
    with open(path, encoding="utf-8", errors="replace") as f:
        while True:
            line = f.readline()
            if not line:
                break
            s = line.strip()
            if not s:
                continue
            try:
                n = int(s.split()[0])
            except ValueError as ex:
                raise HandoffError(f"{path}: xyz の原子数の行が読めません: {s!r}") from ex
            f.readline()
            rows = [f.readline().split() for _ in range(n)]
            if len(rows) != n or any(len(r) < 4 for r in rows):
                break
            last = rows
    if last is None:
        raise HandoffError(f"{path}: フレームがありません")
    return last


def parse_poscar(text: str) -> dict:
    lines = text.splitlines()
    scale = float(lines[1].split()[0])
    cell = [[scale * x for x in _floats(lines[2 + i].split()[:3])] for i in range(3)]
    names = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    symbols = [n for n, c in zip(names, counts) for _ in range(c)]
    k = 7
    if lines[k].strip().lower().startswith("s"):
        k += 1  # Selective dynamics
    direct = lines[k].strip().lower()[:1] not in ("c", "k")
    k += 1
    n = len(symbols)
    pos = []
    for l in lines[k:k + n]:
        v = _floats(l.split()[:3])
        pos.append(_mat_vec(cell, v) if direct else [scale * x for x in v])
    k += n
    vel = None
    if k < len(lines):
        mode = lines[k].strip().lower()
        rows = lines[k + 1:k + 1 + n]
        if len(rows) == n and all(len(r.split()) >= 3 for r in rows):
            v = [_floats(r.split()[:3]) for r in rows]
            if mode.startswith("d"):
                v = None
            vel = v
    return {"symbols": symbols, "positions": pos, "cell": cell, "velocities": vel}


def _qe_block(path: str, key: str, n: int) -> tuple[str, list[list[str]]] | None:
    found = None
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith(key):
                rows = [next(f, "").split() for _ in range(n)]
                found = (line.strip(), rows)
    return found


def _qe_cell(header: str, rows, alat_bohr: float | None) -> list[list[float]]:
    v = [_floats(r[:3]) for r in rows]
    h = header.lower()
    if "angstrom" in h:
        return v
    if "bohr" in h:
        return [[x * BOHR_ANG for x in r] for r in v]
    m = re.search(r"alat\s*=\s*([-\d.Ee+]+)", h)
    a = float(m.group(1)) if m else alat_bohr
    if a is None:
        raise HandoffError("CELL_PARAMETERS (alat) の alat が分かりません")
    return [[x * a * BOHR_ANG for x in r] for r in v]


def read_pw_in_cell(path: str) -> list[list[float]] | None:
    lines = open(path, encoding="utf-8").read().splitlines()
    for i, l in enumerate(lines):
        if l.strip().upper().startswith("CELL_PARAMETERS"):
            return _qe_cell(l.strip(), [lines[i + 1 + j].split() for j in range(3)], None)
    return None


def read_pw_in_nat(path: str) -> int:
    m = re.search(r"\bnat\s*=\s*(\d+)", open(path, encoding="utf-8").read())
    if not m:
        raise HandoffError(f"{path}: nat がありません")
    return int(m.group(1))


def final_espresso(prev: str) -> dict:
    log = os.path.join(prev, "output.log")
    nat = read_pw_in_nat(os.path.join(prev, "pw.in"))
    if not os.path.isfile(log):
        raise HandoffError(f"{log} がありません")
    alat = None
    with open(log, encoding="utf-8", errors="replace") as f:
        for line in f:
            m = re.search(r"lattice parameter \(alat\)\s*=\s*([-\d.]+)\s*a\.u\.", line)
            if m:
                alat = float(m.group(1))
                break
    cellb = _qe_block(log, "CELL_PARAMETERS", 3)
    cell = _qe_cell(*cellb, alat) if cellb else read_pw_in_cell(os.path.join(prev, "pw.in"))
    posb = _qe_block(log, "ATOMIC_POSITIONS", nat)
    if posb is None:
        raise HandoffError(f"{log} に ATOMIC_POSITIONS がありません (最適化・MD の出力ではないか、途中で止まっています)")
    header, rows = posb
    h = header.lower()
    symbols = [r[0] for r in rows]
    raw = [_floats(r[1:4]) for r in rows]
    if "crystal" in h:
        pos = [_mat_vec(cell, v) for v in raw]
    elif "bohr" in h:
        pos = [[x * BOHR_ANG for x in v] for v in raw]
    elif "angstrom" in h:
        pos = raw
    else:  # alat
        pos = [[x * alat * BOHR_ANG for x in v] for v in raw]
    return {"symbols": [re.sub(r"\d+$", "", s) for s in symbols], "positions": pos, "cell": cell, "velocities": None,
            "source": "output.log"}


def final_cp2k(prev: str) -> dict:
    path = os.path.join(prev, existing_name(prev, "{}-1.restart"))
    if not os.path.isfile(path):
        raise HandoffError(f"{path} がありません")
    cell, symbols, pos, in_subsys, in_cell, in_coord, has_vel = [None, None, None], [], [], False, False, False, False
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            s = line.strip()
            u = s.upper()
            if u.startswith("&SUBSYS"):
                in_subsys = True
            elif u.startswith("&END SUBSYS"):
                break
            elif in_subsys and u.startswith("&CELL") and not u.startswith("&CELL_REF") and cell[0] is None:
                in_cell = True
            elif in_cell and u.startswith("&END"):
                in_cell = False
            elif in_cell and u.split()[:1] in (["A"], ["B"], ["C"]):
                cell["ABC".index(u.split()[0])] = _floats(s.split()[1:4])
            elif in_subsys and u.startswith("&COORD"):
                in_coord = True
            elif in_coord and u.startswith("&END"):
                in_coord = False
            elif in_coord:
                p = s.split()
                if p[0].upper() in ("UNIT", "SCALED"):
                    if p[0].upper() == "SCALED" or (len(p) > 1 and p[1].lower() != "angstrom"):
                        raise HandoffError(f"{path}: COORD が Å 以外で書かれています ({s})")
                    continue
                symbols.append(p[0])
                pos.append(_floats(p[1:4]))
            elif in_subsys and u.startswith("&VELOCITY"):
                has_vel = True
    if not pos:
        raise HandoffError(f"{path} に &COORD がありません")
    return {"symbols": [re.sub(r"[^A-Za-z].*$", "", s) for s in symbols], "positions": pos,
            "cell": cell if all(c is not None for c in cell) else None, "velocities": None, "has_velocities": has_vel,
            "source": os.path.basename(path)}


def final_structure(code: str, prev: str, previous_task: str) -> dict:
    j = lambda name: os.path.join(prev, name)  # noqa: E731
    md = previous_task == "molecular_dynamics"
    if code == "dftbplus":
        start = parse_gen(open(j("geometry.gen"), encoding="utf-8").read())
        if md:
            rows = last_xyz_frame(j("geo_end.xyz"))
            vel = [[x / 1000.0 for x in _floats(r[-3:])] for r in rows] if all(len(r) >= 7 for r in rows) else None  # Å/ps → Å/fs
            cell = parse_gen(open(j("geo_end.gen"), encoding="utf-8").read())["cell"] if os.path.isfile(j("geo_end.gen")) else start["cell"]
            return {"symbols": [r[0] for r in rows], "positions": [_floats(r[1:4]) for r in rows], "cell": cell, "velocities": vel,
                    "source": "geo_end.xyz"}
        if previous_task == "geometry_optimization":
            d = parse_gen(open(j("geom.out.gen"), encoding="utf-8").read())
            d["source"] = "geom.out.gen"
            return d
        start["source"] = "geometry.gen"
        return start
    if code == "vasp":
        name = "CONTCAR" if previous_task in ("geometry_optimization", "molecular_dynamics") else "POSCAR"
        d = parse_poscar(open(j(name), encoding="utf-8").read())
        if not md:
            d["velocities"] = None
        d["source"] = name
        return d
    if code in ("xtb", "orca"):
        name = {("xtb", True): "xtb.trj", ("xtb", False): "xtbopt.xyz", ("orca", True): "trajectory.xyz", ("orca", False): "orca.xyz"}[(code, md)]
        if previous_task not in ("geometry_optimization", "molecular_dynamics"):
            name = "struct.xyz" if code == "xtb" else None
        if name is None:
            return _orca_input_coords(j("orca.inp"))
        rows = last_xyz_frame(j(name))
        return {"symbols": [r[0] for r in rows], "positions": [_floats(r[1:4]) for r in rows], "cell": None, "velocities": None, "source": name}
    if code == "espresso":
        return final_espresso(prev)
    if code == "cp2k":
        return final_cp2k(prev)
    raise HandoffError(f"{code} の最終構造はこのファイルでは読みません")


def _orca_input_coords(path: str) -> dict:
    lines = open(path, encoding="utf-8").read().splitlines()
    k = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("* xyz"))
    rows = []
    for l in lines[k + 1:]:
        if l.strip() == "*":
            break
        rows.append(l.split())
    return {"symbols": [r[0] for r in rows], "positions": [_floats(r[1:4]) for r in rows], "cell": None, "velocities": None, "source": "orca.inp"}


def _check_same(symbols_here: list[str], symbols_prev: list[str], where: str) -> None:
    if [s.capitalize() for s in symbols_here] != [s.capitalize() for s in symbols_prev]:
        raise HandoffError(f"{where}: 前の段階の原子の並び ({len(symbols_prev)} 個) が、この段階の入力 ({len(symbols_here)} 個) と合いません")


def rewritten_at_run(code: str, velocities: bool = False) -> list[str]:
    if code == "dftbplus":
        return ["geometry.gen"] + ([VELOCITY_FILE] if velocities else [])
    if code == "xtb":
        return ["struct.xyz"]
    if code == "orca":
        return ["orca.inp"]
    if code == "espresso":
        return ["pw.in"]
    return []


def apply_stage(code: str, prev: str, here: str, previous_task: str, velocities: bool) -> list[str]:
    if code not in CONVERTS:
        return []
    d = final_structure(code, prev, previous_task)
    j = lambda name: os.path.join(here, name)  # noqa: E731
    if code == "dftbplus":
        cur = parse_gen(open(j("geometry.gen"), encoding="utf-8").read())
        _check_same(cur["symbols"], d["symbols"], "geometry.gen")
        cell = d["cell"] if cur["cell"] else None
        open(j("geometry.gen"), "w", encoding="utf-8").write(format_gen(d["symbols"], d["positions"], cell))
        done = ["geometry.gen"]
        if velocities:
            if not d["velocities"]:
                raise HandoffError(f"{prev} の {d.get('source')} に速度がありません")
            open(j(VELOCITY_FILE), "w", encoding="utf-8").write(
                "".join(f"{v[0] * 1000.0:.10f} {v[1] * 1000.0:.10f} {v[2] * 1000.0:.10f}\n" for v in d["velocities"]))  # Å/ps
            done.append(VELOCITY_FILE)
        return done
    if code == "xtb":
        cur = last_xyz_frame(j("struct.xyz"))
        _check_same([r[0] for r in cur], d["symbols"], "struct.xyz")
        body = "".join(f"{s} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n" for s, p in zip(d["symbols"], d["positions"]))
        open(j("struct.xyz"), "w", encoding="utf-8").write(f"{len(d['symbols'])}\nfrom {prev}/{d['source']}\n{body}")
        return ["struct.xyz"]
    if code == "orca":
        lines = open(j("orca.inp"), encoding="utf-8").read().splitlines()
        k = next(i for i, l in enumerate(lines) if l.strip().lower().startswith("* xyz"))
        e = next(i for i in range(k + 1, len(lines)) if lines[i].strip() == "*")
        _check_same([l.split()[0] for l in lines[k + 1:e]], d["symbols"], "orca.inp")
        lines[k + 1:e] = [f"  {s:2s} {p[0]:14.8f} {p[1]:14.8f} {p[2]:14.8f}" for s, p in zip(d["symbols"], d["positions"])]
        open(j("orca.inp"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
        return ["orca.inp"]
    if code == "espresso":
        lines = open(j("pw.in"), encoding="utf-8").read().splitlines()
        k = next(i for i, l in enumerate(lines) if l.strip().upper().startswith("ATOMIC_POSITIONS"))
        n = len(d["symbols"])
        rows = [lines[k + 1 + i].split() for i in range(n)]
        _check_same([re.sub(r"\d+$", "", r[0]) for r in rows], d["symbols"], "pw.in")
        lines[k] = "ATOMIC_POSITIONS angstrom"
        for i, (r, p) in enumerate(zip(rows, d["positions"])):
            lines[k + 1 + i] = " ".join([r[0], f"{p[0]:.10f}", f"{p[1]:.10f}", f"{p[2]:.10f}", *r[4:]])
        if d["cell"]:
            c = next((i for i, l in enumerate(lines) if l.strip().upper().startswith("CELL_PARAMETERS")), None)
            if c is not None:
                lines[c] = "CELL_PARAMETERS angstrom"
                for i in range(3):
                    lines[c + 1 + i] = " ".join(f"{x:.10f}" for x in d["cell"][i])
        open(j("pw.in"), "w", encoding="utf-8").write("\n".join(lines) + "\n")
        return ["pw.in"]
    return []


def run_stage_handoff(code: str, prev: str, previous_task: str, velocities: bool, files: dict[str, str], here: str = ".") -> None:
    for dest, src in files.items():
        s = os.path.join(prev, src)
        if not os.path.isfile(s):
            raise HandoffError(f"前の段階の出力 {s} がありません (前の段階が終わっていないか、失敗しています)")
        shutil.copyfile(s, os.path.join(here, dest))
    apply_stage(code, prev, here, previous_task, velocities)


def _main(argv: list[str]) -> int:
    import argparse
    import json

    ap = argparse.ArgumentParser(prog="handoff.py", description="前の段階の最終構造 (と速度) を、この段階の入力に入れる (ADIT が生成)")
    ap.add_argument("code")
    ap.add_argument("previous_dir")
    ap.add_argument("previous_task")
    ap.add_argument("--velocities", action="store_true")
    ap.add_argument("--files", default="{}", help="写すファイル (JSON: この段階での名前 → 前の段階での名前)")
    a = ap.parse_args(argv)
    try:
        run_stage_handoff(a.code, a.previous_dir, a.previous_task, a.velocities, json.loads(a.files))
    except (HandoffError, OSError, ValueError, StopIteration, IndexError) as ex:
        print(f"handoff.py: {ex}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))
