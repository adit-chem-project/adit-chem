"""Generate inputs that vary one setting, for convergence checks and lattice scans."""

from __future__ import annotations

from adit.errors import AditValueError
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from adit.lang import L
from adit.spec import CalculationSpec

SCAN_FILE = "scan.json"
TABLE_FILE = "scan_energies.csv"
SCALE = "scale"


class ScanError(AditValueError):
    pass


@dataclass(frozen=True)
class Scan:
    path: str            # "method.ecutwfc" / "kpoints.mesh" / "scale"
    values: list[str]

    @property
    def leaf(self) -> str:
        return self.path.split(".")[-1]

    def dir_name(self, value: str) -> str:
        safe = re.sub(r"[^0-9A-Za-z.+-]+", "_", value).strip("_") or "v"
        leaf = re.sub(r"[^0-9A-Za-z.+-]+", "_", self.leaf).strip("_") or "v"
        return f"{leaf}_{safe}"


def parse_scan(text: str) -> Scan:
    if "=" not in text:
        raise ScanError(L(f"変える項目と値は「項目=値,値,…」の形で書いてください (例 method.ecutwfc=30,40,50): {text!r}",
                          f"write it as item=value,value,... (e.g. method.ecutwfc=30,40,50): {text!r}"))
    path, _, rest = text.partition("=")
    path = path.strip()
    values = [v.strip() for v in rest.split(",") if v.strip()]
    if not (GEOM.match(path) or re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*(\.[A-Za-z_][A-Za-z0-9_]*)*", path)):
        raise ScanError(L(f"項目の名前の書き方が正しくありません: {path!r}", f"Invalid item name: {path!r}"))
    if len(values) < 2:
        raise ScanError(L("値はカンマで区切って 2 つ以上書いてください", "give at least two values, separated by commas"))
    if len(set(values)) != len(values):
        raise ScanError(L(f"同じ値が 2 回書かれています: {values}", f"The same value appears twice: {values}"))
    return Scan(path, values)


def _convert(value: str, current):
    if isinstance(current, (list, tuple)):
        parts = [p for p in re.split(r"[x×,\s]+", value) if p]
        return [type(current[0])(p) if current else p for p in parts]
    if isinstance(current, bool):
        return value.lower() in ("1", "true", "yes", "on")
    if isinstance(current, int):
        return int(value)
    if isinstance(current, float):
        return float(value)
    return value


GEOM = re.compile(r"geom\.(distance|angle|dihedral)\(([0-9,\s]+)\)\Z")


def _apply_geometry(spec: CalculationSpec, path: str, value: str) -> dict:
    m = GEOM.match(path)
    kind, raw = m.group(1), m.group(2)
    indices = [int(x) for x in raw.replace(" ", "").split(",") if x != ""]
    need = {"distance": 2, "angle": 3, "dihedral": 4}[kind]
    atoms = spec.atoms
    if len(indices) != need:
        raise ScanError(L(f"geom.{kind} には原子の番号を {need} 個書いてください (0 始まり)",
                          f"geom.{kind} needs {need} atom indices (0-based)"))
    if any(i < 0 or i >= len(atoms) for i in indices):
        raise ScanError(L(f"原子の番号が範囲の外です (0〜{len(atoms) - 1}): {indices}",
                          f"atom index out of range (0-{len(atoms) - 1}): {indices}"))
    try:
        x = float(value)
    except ValueError as ex:
        raise ScanError(L(f"{value!r} は数として読めません", f"{value!r} is not a number")) from ex
    if kind == "distance":
        if x <= 0:
            raise ScanError(L("距離は正の数にしてください", "the distance must be positive"))
        atoms.set_distance(indices[0], indices[1], x, fix=0)
    elif kind == "angle":
        atoms.set_angle(*indices, angle=x)
    else:
        atoms.set_dihedral(*indices, angle=x)
    from adit.spec import AtomsData

    data = spec.model_dump(mode="json")
    data["structure"]["atoms"] = AtomsData.from_ase(atoms).model_dump(mode="json")
    data["structure"]["source_ref"] = spec.structure.source_ref
    return data


def apply_value(spec: CalculationSpec, path: str, value: str) -> CalculationSpec:
    data = spec.model_dump(mode="json")
    if GEOM.match(path):
        data = _apply_geometry(spec, path, value)
    elif path == SCALE:
        s = float(value)
        if s <= 0:
            raise ScanError(L("倍率は正の数にしてください", "The scale factor must be positive"))
        if not spec.structure.periodic:
            raise ScanError(L("格子の大きさ (scale) を変えられるのは周期系だけです", "Scaling the cell (scale) works only for periodic systems"))
        atoms = data["structure"]["atoms"]
        atoms["cell"] = [[c * s for c in v] for v in atoms["cell"]]
        atoms["positions"] = [[c * s for c in v] for v in atoms["positions"]]
    else:
        keys = path.split(".")
        if keys[0] == "kpoints" and data.get("kpoints") is None:
            if not spec.structure.periodic:
                raise ScanError(L("分子 (非周期) では k 点を使いません", "Molecules (non-periodic) have no k-points"))
            data["kpoints"] = {"mode": "mesh"}
            from adit.spec import KPoints
            data["kpoints"] = KPoints(mode="mesh").model_dump(mode="json")
        node = data
        for k in keys[:-1]:
            if not isinstance(node, dict) or k not in node or node[k] is None:
                raise ScanError(L(f"計算設定に {path!r} という項目はありません ({k!r} の所で見つからなくなりました)",
                                  f"The calculation settings have no item {path!r} (nothing found at {k!r})"))
            node = node[k]
        leaf = keys[-1]
        if not isinstance(node, dict) or leaf not in node:
            have = ", ".join(sorted(node)) if isinstance(node, dict) else ""
            raise ScanError(L(f"計算設定に {path!r} という項目はありません。この階層にある項目: {have}",
                              f"The calculation settings have no item {path!r}. Items at this level: {have}"))
        try:
            node[leaf] = _convert(value, node[leaf])
        except (TypeError, ValueError) as ex:
            raise ScanError(L(f"{value!r} は {path} の値として読めません ({ex})", f"{value!r} is not a valid value for {path} ({ex})")) from ex
        if keys[:2] == ["kpoints", "mesh"] and isinstance(data.get("kpoints"), dict):
            data["kpoints"]["mode"] = "mesh"
        if keys[:2] == ["kpoints", "density"] and isinstance(data.get("kpoints"), dict):
            data["kpoints"]["mode"] = "density"
    data.setdefault("meta", {})
    data["meta"]["comment"] = (data["meta"].get("comment", "") + f" scan {path}={value}").strip()
    try:
        return CalculationSpec.model_validate(data)
    except Exception as ex:
        from pydantic import ValidationError as PydanticError
        if isinstance(ex, PydanticError):
            from adit.validate_types import friendly_pydantic
            raise ScanError(friendly_pydantic(ex)) from ex
        raise


def write_scan(spec: CalculationSpec, cfg, out_dir: Path | str, scan: Scan | list[Scan], *,
               overwrite: bool = False) -> list[Path]:
    """Write one directory per value. Nothing is written unless every value assembles. Several scans give the full grid of combinations."""
    import itertools

    from adit.project import build_project, write_project

    scans = [scan] if isinstance(scan, Scan) else list(scan)
    if not scans:
        raise ScanError(L("振る項目がありません", "no parameter to scan"))
    out = Path(out_dir).expanduser()
    combos = list(itertools.product(*[s.values for s in scans]))
    specs = []
    for combo in combos:
        one = spec
        for sc, v in zip(scans, combo):
            one = apply_value(one, sc.path, v)
        specs.append((combo, one))
    from adit.batch import check_output

    check_output(out, overwrite)
    out.mkdir(parents=True, exist_ok=True)
    names = ["__".join(sc.dir_name(v) for sc, v in zip(scans, combo)) for combo, _ in specs]
    for name, (_, s) in zip(names, specs):
        build_project(s, cfg, output_dir=out / name)
    dirs = []
    for name, (_, s) in zip(names, specs):
        d = out / name
        write_project(s, cfg, d, overwrite=overwrite)
        dirs.append(d)
    record = {"scans": [{"path": sc.path, "values": sc.values} for sc in scans],
              "dirs": [d.name for d in dirs],
              "combinations": [list(combo) for combo, _ in specs]}
    if len(scans) == 1:
        record["path"], record["values"] = scans[0].path, scans[0].values
    (out / SCAN_FILE).write_text(json.dumps(record, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return dirs


from adit.spec import (EV_PER_ANG3_IN_GPA as EV_ANG3_TO_GPA, HARTREE_PER_BOHR_IN_EV_PER_ANG as HARTREE_BOHR_TO_EV_ANG,
                        RY_PER_BOHR_IN_EV_PER_ANG as RY_BOHR_TO_EV_ANG)
MIN_EOS_POINTS = 5
COLUMNS = ["value", "dir", "natoms", "energy_ev", "diff_from_last_mev", "diff_per_atom_mev", "fmax_ev_ang", "pressure_gpa", "note"]


def final_force_and_pressure(d: Path, code: str) -> tuple[float | None, float | None]:
    import numpy as np

    def text(p: Path) -> str:
        return p.read_text(encoding="utf-8", errors="replace") if p.is_file() else ""

    fmax = pres = None
    if code == "dftbplus":
        t = text(d / "detailed.out")
        if "Total Forces" in t:
            block = t.split("Total Forces", 1)[1].split("\n\n", 1)[0]
            rows = re.findall(r"^\s*\d+\s+(-?[\d.]+(?:[EeDd][-+]?\d+)?)\s+(-?[\d.]+(?:[EeDd][-+]?\d+)?)\s+(-?[\d.]+(?:[EeDd][-+]?\d+)?)\s*$", block, re.M)
            if rows:
                f = np.array([[float(x.replace("D", "E").replace("d", "e")) for x in r] for r in rows])
                fmax = float(np.linalg.norm(f, axis=1).max()) * HARTREE_BOHR_TO_EV_ANG
        m = re.findall(r"^Pressure:.*?(-?[\d.]+E[-+]\d+)\s+Pa", t, re.M)
        if m:
            pres = float(m[-1]) / 1e9
    elif code == "espresso":
        t = text(d / "output.log")
        if "Forces acting on atoms" in t:
            block = t.rsplit("Forces acting on atoms", 1)[1].split("Total force", 1)[0]
            rows = re.findall(r"atom\s+\d+\s+type\s+\d+\s+force\s+=\s+(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", block)
            if rows:
                f = np.array([[float(x) for x in r] for r in rows])
                fmax = float(np.linalg.norm(f, axis=1).max()) * RY_BOHR_TO_EV_ANG
        m = re.findall(r"total\s+stress.*?P=\s*(-?[\d.]+)", t)
        if m:
            pres = float(m[-1]) / 10.0  # kbar → GPa
    elif code == "vasp" and (d / "vasprun.xml").is_file():
        try:
            from ase.io import read
            a = read(str(d / "vasprun.xml"), index=-1)
            fmax = float(np.linalg.norm(a.get_forces(), axis=1).max())
            try:
                s = a.get_stress(voigt=True)
                pres = float(-np.mean(s[:3]) * EV_ANG3_TO_GPA)
            except Exception:
                pass
        except Exception:
            pass
    return fmax, pres


def fit_eos(volumes: list[float], energies: list[float]) -> dict:
    from ase.eos import EquationOfState, birchmurnaghan

    v = np.asarray(volumes, dtype=float); e = np.asarray(energies, dtype=float)
    v0, e0, b = EquationOfState(list(v), list(e), eos="birchmurnaghan").fit(warn=False)
    params, bp = [float(e0), float(b), 4.0, float(v0)], None
    try:
        from scipy.optimize import curve_fit
        popt, _ = curve_fit(birchmurnaghan, v, e, p0=params, maxfev=20000)
        params, bp = [float(x) for x in popt], float(popt[2])
    except Exception:
        pass
    e0, b, v0 = params[0], params[1], params[3]
    return {"v0_A3": v0, "e0_ev": e0, "b0_gpa": b * EV_ANG3_TO_GPA, "bp": bp, "params": params}


@dataclass
class ScanResult:
    path: str
    rows: list[dict]
    eos: dict | None = None
    eos_note: str = ""
    figures: dict = None
    table: str = ""
    reference: str = "last"
    notes: list[str] | None = None

    def summary_text(self, include_table: bool = True) -> str:
        if self.notes:
            return "\n".join(self.notes)
        """include_table=False は値ごとの行を省く (GUI は表を別の部品に出すので、要約にはそれ以外だけを残す)。"""
        table = self.table if include_table else Path(self.table).name
        lines = [L(f"変えた項目: {self.path}   値の数: {len(self.rows)}   表: {table}",
                   f"Scanned item: {self.path}   Number of values: {len(self.rows)}   Table: {table}")]
        if include_table:
            ref_ja, ref_en = diff_label(self.reference)
            lines.append(L(f"値\tエネルギー [eV]\t{ref_ja} [meV/原子]\t力の最大値 [eV/Å]\t圧力 [GPa]\t注",
                           f"Value\tEnergy [eV]\t{ref_en} [meV/atom]\tMax force [eV/Å]\tPressure [GPa]\tNote"))
            for r in self.rows:
                lines.append("\t".join([r["value"], r["energy_ev"] or "-", r["diff_per_atom_mev"] or "-", r["fmax_ev_ang"] or "-",
                                        r["pressure_gpa"] or "-", r["note"]]))
        if self.eos:
            e = self.eos
            abc = " × ".join(f"{x:.4f}" for x in e["cell_lengths_A"])
            bp = f"、B0' {e['bp']:.2f}" if e.get("bp") is not None else ""
            bp_en = f", B0' {e['bp']:.2f}" if e.get("bp") is not None else ""
            lines.append(L(f"状態方程式 (Birch–Murnaghan) の当てはめ: 平衡の体積 {e['v0_A3']:.3f} Å³、体積弾性率 {e['b0_gpa']:.1f} GPa{bp}、"
                           f"エネルギーが最小になる倍率 {e['scale0']:.4f} (そのときの格子の長さ {abc} Å)",
                           f"Birch–Murnaghan fit: equilibrium volume {e['v0_A3']:.3f} Å³, bulk modulus {e['b0_gpa']:.1f} GPa{bp_en}, "
                           f"scale at the energy minimum {e['scale0']:.4f} (cell lengths {abc} Å)"))
        elif self.eos_note:
            lines.append(L("注: ", "Note: ") + self.eos_note)
        lines.append(L("どの値で十分とみなすかは、計算の目的に合わせて決めてください (ADIT は判断しません)",
                       "Decide which value is good enough for your purpose; ADIT does not judge this"))
        return "\n".join(lines)


def _reference_row(rows: list[dict], path: str) -> tuple[dict | None, str]:
    have = [r for r in rows if r["energy_ev"]]
    if not have:
        return None, "last"
    if path == SCALE:
        return min(have, key=lambda r: float(r["energy_ev"])), "min"
    return have[-1], "last"


def diff_label(reference: str) -> tuple[str, str]:
    return (("最小値との差", "Diff. from minimum") if reference == "min" else ("最後の値との差", "Diff. from last"))


def _is_number(s: str) -> bool:
    try:
        float(s); return True
    except ValueError:
        return False


def _analyze_grid(out: Path, meta: dict) -> "ScanResult":
    import csv as _csv

    from adit.analysis.readers import load_run

    paths = [sc["path"] for sc in meta["scans"]]
    columns = paths + ["dir", "energy_ev", "diff_from_min_mev", "note"]
    rows = []
    for combo, name in zip(meta["combinations"], meta["dirs"]):
        row = {p: v for p, v in zip(paths, combo)}
        row.update({"dir": name, "energy_ev": "", "diff_from_min_mev": "", "note": ""})
        try:
            run = load_run(out / name)
            if run.energies_ev:
                row["energy_ev"] = f"{run.energies_ev[-1]:.8f}"
            else:
                row["note"] = L("エネルギーが見つかりません (まだ計算していないかもしれません)",
                                "No energy found (perhaps not run yet)")
        except Exception as ex:
            row["note"] = L(f"読めません: {ex}", f"Cannot read: {ex}")
        rows.append(row)
    energies = [float(r["energy_ev"]) for r in rows if r["energy_ev"]]
    if energies:
        lowest = min(energies)
        for r in rows:
            if r["energy_ev"]:
                r["diff_from_min_mev"] = f"{(float(r['energy_ev']) - lowest) * 1000:.3f}"
    table = out / TABLE_FILE
    with open(table, "w", encoding="utf-8", newline="") as f:
        writer = _csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    lines = [L(f"振った項目: {', '.join(paths)} ({len(rows)} 通りの組み合わせ)",
               f"scanned: {', '.join(paths)} ({len(rows)} combinations)"),
             L(f"表: {table}", f"table: {table}"),
             L("エネルギーの差は、読めた中で最も低い値からの差 [meV] です。",
               "The energy difference is measured from the lowest value that could be read, in meV."),
             L("図は項目を 1 つだけ振ったときに描きます (2 つ以上の組み合わせでは表だけです)。",
               "Figures are drawn only when a single parameter is scanned; for a grid only the table is written.")]
    missing = [r["dir"] for r in rows if not r["energy_ev"]]
    if missing:
        lines.append(L(f"エネルギーを読めなかったディレクトリ: {', '.join(missing)}",
                       f"directories with no readable energy: {', '.join(missing)}"))
    return ScanResult(path=" × ".join(paths), rows=rows, table=str(table), figures={}, notes=lines)


def analyze_scan(out_dir: Path | str) -> ScanResult:
    from adit.analysis.readers import load_run
    from adit.project import load_project

    out = Path(out_dir).expanduser()
    meta = json.loads((out / SCAN_FILE).read_text(encoding="utf-8"))
    if "scans" in meta and len(meta["scans"]) > 1:
        return _analyze_grid(out, meta)
    rows, volumes = [], []
    for v, name in zip(meta["values"], meta["dirs"]):
        d = out / name
        row = {c: "" for c in COLUMNS}; row.update(value=v, dir=name)
        vol = None
        try:
            st = load_project(d).structure
            row["natoms"] = str(len(st.atoms.symbols))
            if st.periodic:
                vol = float(st.atoms.to_ase().get_volume())
        except Exception:
            pass
        try:
            run = load_run(d)
            if run.energies_ev:
                row["energy_ev"] = f"{run.energies_ev[-1]:.8f}"
                if not row["natoms"] and run.final is not None:
                    row["natoms"] = str(len(run.final))
            else:
                row["note"] = L("エネルギーが見つかりません (まだ計算していないかもしれません)", "No energy found (perhaps not run yet)")
            fmax, pres = final_force_and_pressure(d, run.code)
            row["fmax_ev_ang"] = f"{fmax:.6g}" if fmax is not None else ""
            row["pressure_gpa"] = f"{pres:.6g}" if pres is not None else ""
        except Exception as ex:
            row["note"] = L(f"読めません: {ex}", f"Cannot read: {ex}")
        rows.append(row); volumes.append(vol)
    base, reference = _reference_row(rows, meta["path"])
    for r in rows:
        if r["energy_ev"] and base is not None:
            diff = (float(r["energy_ev"]) - float(base["energy_ev"])) * 1000
            r["diff_from_last_mev"] = f"{diff:.3f}"
            if r["natoms"] and r["natoms"] == base["natoms"]:
                r["diff_per_atom_mev"] = f"{diff / int(r['natoms']):.3f}"
    with open(out / TABLE_FILE, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS); w.writeheader(); w.writerows(rows)
    res = ScanResult(path=meta["path"], rows=rows, figures={}, table=str(out / TABLE_FILE), reference=reference)

    if meta["path"] == SCALE:
        pts = [(float(r["value"]), vol, float(r["energy_ev"])) for r, vol in zip(rows, volumes) if r["energy_ev"] and vol]
        if len(pts) < MIN_EOS_POINTS:
            res.eos_note = L(f"状態方程式を当てはめるには、エネルギーのある点が {MIN_EOS_POINTS} つ以上要ります (いまは {len(pts)} つ)",
                             f"Fitting the equation of state needs at least {MIN_EOS_POINTS} points with an energy (now {len(pts)})")
        else:
            try:
                fit = fit_eos([p[1] for p in pts], [p[2] for p in pts])
                v_ref = float(np.median([vol / s ** 3 for s, vol, _ in pts]))
                scale0 = (fit["v0_A3"] / v_ref) ** (1 / 3)
                ref_dir = next(out / r["dir"] for r in rows if r["energy_ev"])
                s_ref = float(next(r["value"] for r in rows if r["energy_ev"]))
                lengths = load_project(ref_dir).structure.atoms.to_ase().cell.lengths() / s_ref
                fit.update(scale0=float(scale0), cell_lengths_A=[float(x * scale0) for x in lengths])
                if not (min(p[1] for p in pts) <= fit["v0_A3"] <= max(p[1] for p in pts)):
                    res.eos_note = L("当てはめた平衡の体積が、計算した範囲の外にあります。倍率の範囲を広げて計算し直してください",
                                     "The fitted equilibrium volume is outside the scanned range; widen the scale range and run again")
                res.eos = fit
                (out / "eos.json").write_text(json.dumps(fit, indent=2) + "\n", encoding="utf-8")
            except Exception as ex:
                res.eos_note = L(f"状態方程式を当てはめられませんでした: {ex}", f"Could not fit the equation of state: {ex}")
    _plot_scan(out, res, volumes)
    return res


def _plot_scan(out: Path, res: ScanResult, volumes: list) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    got = [(i, r) for i, r in enumerate(res.rows) if r["diff_per_atom_mev"] or r["diff_from_last_mev"]]
    if len(got) >= 2:
        per_atom = all(r["diff_per_atom_mev"] for _, r in got)
        numeric = all(_is_number(r["value"]) for _, r in got)
        x = [float(r["value"]) if numeric else i for i, r in got]
        y = [float(r["diff_per_atom_mev"] if per_atom else r["diff_from_last_mev"]) for _, r in got]
        fig, ax = plt.subplots(figsize=(6, 3.4))
        ax.plot(x, y, "o-", lw=1.2); ax.axhline(0, color="gray", lw=0.8, ls="--"); ax.grid(alpha=0.3)
        if not numeric:
            ax.set_xticks(x); ax.set_xticklabels([r["value"] for _, r in got])
        ref = "min" if res.reference == "min" else "last"
        ax.set_xlabel(res.path); ax.set_ylabel(f"E - E({ref}) [meV/atom]" if per_atom else f"E - E({ref}) [meV]")
        p = out / "scan_energy.png"; fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
        res.figures["scan_energy"] = str(p)
    if res.eos:
        import numpy as np
        from ase.eos import birchmurnaghan
        pts = [(vol, float(r["energy_ev"])) for r, vol in zip(res.rows, volumes) if r["energy_ev"] and vol]
        v = np.array([p[0] for p in pts]); e = np.array([p[1] for p in pts])
        fig, ax = plt.subplots(figsize=(6, 3.4))
        ax.plot(v, e, "o", label="calc.")
        vv = np.linspace(min(v.min(), res.eos["v0_A3"]), max(v.max(), res.eos["v0_A3"]), 200)
        ax.plot(vv, birchmurnaghan(vv, *res.eos["params"]), "-", label="Birch–Murnaghan")
        ax.axvline(res.eos["v0_A3"], color="gray", lw=0.8, ls="--")
        ax.set_xlabel("volume [Å$^3$]"); ax.set_ylabel("energy [eV]"); ax.legend(); ax.grid(alpha=0.3)
        ax.ticklabel_format(axis="y", useOffset=False, style="plain")
        p = out / "eos.png"; fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
        res.figures["eos"] = str(p)


def collect_scan(out_dir: Path | str) -> list[dict]:
    return analyze_scan(out_dir).rows


__all__ = ["Scan", "ScanError", "parse_scan", "apply_value", "write_scan", "collect_scan"]
