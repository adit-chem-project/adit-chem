
from __future__ import annotations

from adit.errors import AditValueError
import csv
import json
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from adit.lang import L

COMPARE_FILE = "compare.json"
EV_KJ_MOL = 96.485332123
EV_KCAL_MOL = 23.060547830619
BALANCE_TOL = 1e-9

_SKIP_PREFIX = ("meta.created", "meta.comment", "meta.provenance.generated_utc", "meta.provenance.generated_files",
                "meta.continued_from", "meta.template", "meta.stage", "handoff", "runtime", "structure")

_ENERGY_SOURCE = {
    "dftbplus": ("output.log の Total Energy (MD では md.out の Total MD Energy)", "Total Energy in output.log (Total MD Energy in md.out for MD)"),
    "xtb": ("output.log の total energy (MD では xtb.trj の energy)", "total energy in output.log (energy in xtb.trj for MD)"),
    "vasp": ("vasprun.xml の各ステップのエネルギー (ASE の読み。energy(sigma→0))", "energy of each step in vasprun.xml (ASE reader, energy(sigma->0))"),
    "espresso": ("output.log の「!    total energy」", "'!    total energy' in output.log"),
    "orca": ("output.log の FINAL SINGLE POINT ENERGY", "FINAL SINGLE POINT ENERGY in output.log"),
    "nwchem": ("output.log の Total SCF/DFT energy", "Total SCF/DFT energy in output.log"),
    "lammps": ("log.lammps の thermo の表 (MD は TotEng、それ以外は PotEng)", "thermo table in log.lammps (TotEng for MD, PotEng otherwise)"),
    "gromacs": ("GROMACS のエネルギーの表 (readers_extra の読み)", "GROMACS energy table (readers_extra)"),
    "cp2k": ("CP2K の output.log (readers_extra の読み)", "CP2K output.log (readers_extra)"),
    "openmm": ("md.log の Total Energy と results.json (run_openmm.py の生成したファイル)", "Total Energy in md.log and results.json (written by run_openmm.py)"),
    "abinit": ("output.log の Total energy (etotal) [Ha] (無ければ input.abo の最後の etotal)", "Total energy (etotal) [Ha] in output.log (or the final etotal in input.abo)"),
    "psi4": ("results.json の energy_hartree (入力の末尾の数行が書いた Psi4 の返り値)", "energy_hartree in results.json (the Psi4 return value written by the last lines of the input)"),
}
_CODE_G_KEYS = ("total_free_energy", "final_gibbs_free_energy")


class CompareError(AditValueError):
    pass


@dataclass
class Reaction:
    name: str
    terms: list[tuple[float, str]]  # (ν, dir)

    def text(self) -> str:
        return " ".join(f"{nu:+g}·{d}" for nu, d in self.terms)


def parse_compare(text: str) -> list[Reaction]:
    out = []
    for k, chunk in enumerate([c for c in text.split(";") if c.strip()], start=1):
        name, body = f"r{k}", chunk.strip()
        m = re.match(r"^\s*([A-Za-z_][\w.-]*)\s*=(.*)$", body)
        if m and ":" not in m.group(1):
            name, body = m.group(1), m.group(2)
        terms = []
        for t in [t for t in body.split(",") if t.strip()]:
            if ":" not in t:
                raise CompareError(L(f"項は「係数:ディレクトリ」の形で書いてください (例 -1:slab): {t.strip()!r}",
                                     f"write each term as coefficient:directory (e.g. -1:slab): {t.strip()!r}"))
            nu_s, d = t.split(":", 1)
            try:
                nu = float(nu_s)
            except ValueError as ex:
                raise CompareError(L(f"係数を数として読めません: {nu_s.strip()!r}", f"cannot read the coefficient as a number: {nu_s.strip()!r}")) from ex
            if nu == 0 or not np.isfinite(nu):
                raise CompareError(L(f"係数は 0 以外の有限の数にしてください: {t.strip()!r}", f"the coefficient must be a non-zero finite number: {t.strip()!r}"))
            terms.append((nu, d.strip()))
        if not terms:
            raise CompareError(L(f"反応 {name} に項がありません", f"reaction {name} has no terms"))
        out.append(Reaction(name, terms))
    if not out:
        raise CompareError(L("比べる組が書かれていません", "no reactions given"))
    names = [r.name for r in out]
    if len(set(names)) != len(names):
        raise CompareError(L(f"同じ名前の反応が 2 つあります: {names}", f"two reactions have the same name: {names}"))
    return out


def load_compare_file(path: Path | str) -> list[Reaction]:
    path = Path(path)
    try:
        js = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as ex:
        raise CompareError(L(f"{path.name} を JSON として読めません: {ex}", f"cannot read {path.name} as JSON: {ex}")) from ex
    if not isinstance(js, dict) or not isinstance(js.get("reactions"), list):
        raise CompareError(L(f"{path.name} に反応の一覧 (reactions) がありません", f"{path.name} has no list of reactions"))
    out = []
    for k, r in enumerate(js["reactions"], start=1):
        try:
            terms = [(float(t["nu"]), str(t["dir"])) for t in (r.get("terms") or [])]
        except (AttributeError, KeyError, TypeError, ValueError) as ex:
            raise CompareError(L(f"compare.json の {k} 番目の反応の項を読めません (各項は {{\"dir\": …, \"nu\": 数}}): {ex}",
                                 f"cannot read the terms of reaction {k} in compare.json (each term is {{\"dir\": ..., \"nu\": number}}): {ex}")) from ex
        if not terms or any(nu == 0 or not np.isfinite(nu) for nu, _ in terms):
            raise CompareError(L(f"compare.json の {k} 番目の反応に項が無いか、係数が 0 か数でない項があります",
                                 f"reaction {k} in compare.json has no terms, or a term with a zero or non-finite coefficient"))
        out.append(Reaction(str(r.get("name") or f"r{k}"), terms))
    if not out:
        raise CompareError(L("compare.json に reactions がありません", "compare.json has no reactions"))
    return out


def _flatten(x, prefix: str = "") -> dict[str, object]:
    out: dict[str, object] = {}
    if isinstance(x, dict):
        for k, v in x.items():
            out.update(_flatten(v, f"{prefix}.{k}" if prefix else str(k)))
    elif isinstance(x, list) and x and all(isinstance(v, dict) for v in x):
        for k, v in enumerate(x):
            name = v.get("name") if isinstance(v.get("name"), str) else str(k)
            out.update(_flatten({kk: vv for kk, vv in v.items() if kk not in ("name", "source")}, f"{prefix}.{name}"))
    else:
        out[prefix] = json.dumps(x, ensure_ascii=False) if isinstance(x, (list, tuple)) else x
    return out


def conditions_of(spec) -> dict[str, object]:
    data = spec.model_dump(mode="json")
    flat = _flatten({k: v for k, v in data.items() if k in ("method", "kpoints", "task", "meta")})
    flat = {k: v for k, v in flat.items() if not k.startswith(_SKIP_PREFIX) and not k.startswith("meta.app_version")}
    task = spec.task.type
    if task != "molecular_dynamics":
        flat = {k: v for k, v in flat.items() if not k.startswith("task.md.")}
    if task != "band_structure":
        flat = {k: v for k, v in flat.items() if not k.startswith("task.bands.")}
    flat["structure.charge"] = spec.structure.charge
    flat["structure.multiplicity"] = spec.structure.multiplicity
    flat["structure.periodic"] = bool(spec.structure.periodic)
    if spec.kpoints is not None and spec.structure.periodic and spec.kpoints.mode in ("mesh", "gamma"):
        cell = spec.structure.atoms.to_ase().cell
        b = 2 * np.pi * np.linalg.norm(np.asarray(cell.reciprocal()), axis=1)
        mesh = spec.kpoints.mesh if spec.kpoints.mode == "mesh" else (1, 1, 1)
        flat["kpoints.spacing_2pi_per_A"] = json.dumps([float(f"{x / n:.4g}") for x, n in zip(b, mesh)])
    return flat


@dataclass
class CompareResult:
    base: str
    runs: list[dict]
    reactions: list[dict]
    conditions: dict[str, dict[str, object]]
    differing: list[str]
    partial: list[str]
    files: dict[str, str] = field(default_factory=dict)
    figures: dict[str, str] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def summary_text(self) -> str:
        lines = [L(f"比べた計算: {len(self.runs)} 個、反応: {len(self.reactions)} 個   ディレクトリ: {self.base}",
                   f"runs compared: {len(self.runs)}, reactions: {len(self.reactions)}   directory: {self.base}")]
        for r in self.runs:
            e = f"{r['energy_ev']:.6f} eV" if r["energy_ev"] is not None else L("(エネルギーなし)", "(no energy)")
            lines.append(f"  {r['dir']}: {r['code'] or '?'} / {r['task'] or '?'} / {r['formula'] or '?'} ({r['natoms'] or '?'} atoms) / {e}")
        for x in self.reactions:
            if x["delta_e_ev"] is not None:
                differing = L(f" [条件が違う項目 {x['n_differing']} 個]", f" [{x['n_differing']} settings differ]")
                val = f"ΔE = {x['delta_e_ev']:+.6f} eV = {x['delta_e_kj_mol']:+.3f} kJ/mol = {x['delta_e_kcal_mol']:+.3f} kcal/mol{differing}"
            else:
                val = L("ΔE は出せません (エネルギーの無い計算がある)", "ΔE not available (a run has no energy)")
            lines.append(f"{x['name']}: {x['terms']}   {val}")
            if x.get("delta_g_code_ev") is not None:
                lines.append(L(f"  コードが出した G での差 ΣνG = {x['delta_g_code_ev']:+.6f} eV", f"  difference of the code-reported G: ΣνG = {x['delta_g_code_ev']:+.6f} eV"))
            if not x["balanced"]:
                lines.append(L(f"  組成が釣り合っていません (Σν·組成 = {x['imbalance']})", f"  composition is not balanced (Σν·composition = {x['imbalance']})"))
            lines.append(L(f"  条件が違う項目: {len(x['differing'])} 個 {x['differing'][:10]}{' …' if len(x['differing']) > 10 else ''}"
                           f"、一部の計算にだけある項目: {len(x['partial'])} 個",
                           f"  settings that differ: {len(x['differing'])} {x['differing'][:10]}{' ...' if len(x['differing']) > 10 else ''}"
                           f", settings present in only some runs: {len(x['partial'])}"))
        lines += [L("注: ", "note: ") + n for n in self.notes]
        lines.append(L("比べてよいか (条件の違いが結果に効くか) は、ADIT は判断しません", "ADIT does not judge whether these runs can be compared"))
        return "\n".join(lines)


def _formula(c: Counter) -> str:
    return "".join(f"{el}{'' if n == 1 else n}" for el, n in sorted(c.items()))


def _run_row(base: Path, d: str) -> tuple[dict, dict | None]:
    from adit import lang
    from adit.analysis.readers import load_run

    p = Path(d) if Path(d).is_absolute() else base / d
    row = {"dir": d, "path": str(p), "code": "", "task": "", "natoms": None, "formula": "", "composition": {}, "energy_ev": None,
           "energy_source": "", "code_g_ev": None, "code_g_label": "", "note": ""}
    cond = None
    if not p.is_dir():
        row["note"] = L("ディレクトリがありません", "directory not found")
        return row, None
    spec = None
    try:
        from adit.project import load_project
        spec = load_project(p)
        cond = conditions_of(spec)
        row["task"] = spec.task.type
    except Exception:
        row["note"] = L("spec.json が読めないので、条件は比べません", "spec.json cannot be read; settings are not compared")
    try:
        run = load_run(p)
    except Exception as ex:
        row["note"] = (row["note"] + " / " if row["note"] else "") + L(f"出力を読めません: {ex}", f"cannot read the output: {ex}")
        run = None
    atoms = None
    if run is not None:
        row["code"] = run.code
        if run.energies_ev:
            row["energy_ev"] = float(run.energies_ev[-1])
            ja, en = _ENERGY_SOURCE.get(run.code, (run.code, run.code))
            row["energy_source"] = en if lang.LANGUAGE == "en" else ja
        else:
            row["note"] = (row["note"] + " / " if row["note"] else "") + L("エネルギーが見つかりません (まだ計算していないかもしれません)", "no energy found (perhaps not run yet)")
        if run.thermo:
            for it in run.thermo.get("items", []):
                if it["key"] in _CODE_G_KEYS:
                    row["code_g_ev"], row["code_g_label"] = float(it["value_ev"]), f"{run.code}: {it['label']} ({it['line']})"
        try:
            atoms = run.final
        except Exception:
            atoms = None
    if atoms is None and spec is not None:
        atoms = spec.structure.atoms.to_ase()
    if atoms is not None:
        c = Counter(atoms.get_chemical_symbols())
        row.update(natoms=len(atoms), formula=_formula(c), composition=dict(c))
    if not row["code"] and spec is not None:
        row["code"] = spec.method.code
    return row, cond


def _diff_keys(conds: dict[str, dict | None], dirs: list[str]) -> tuple[list[str], list[str], dict[str, dict[str, object]]]:
    have = {d: conds[d] for d in dirs if conds.get(d) is not None}
    keys = sorted({k for c in have.values() for k in c})
    table = {k: {d: c.get(k, None) for d, c in have.items()} for k in keys}
    differing, partial = [], []
    for k in table:
        present = [json.dumps(have[d][k], sort_keys=True) for d in have if k in have[d]]
        if len(set(present)) >= 2:
            differing.append(k)
        if len(present) < len(have):
            partial.append(k)
    return differing, partial, table


def analyze_compare(base: Path | str, reactions: list[Reaction] | None = None) -> CompareResult:
    base = Path(base).expanduser()
    if reactions is None:
        if not (base / COMPARE_FILE).is_file():
            raise CompareError(L(f"{base} に {COMPARE_FILE} がありません (--compare \"名前=係数:ディレクトリ,…\" で直接書けます)",
                                 f"{base} has no {COMPARE_FILE} (you can also write --compare \"name=coef:dir,...\")"))
        reactions = load_compare_file(base / COMPARE_FILE)
    dirs = list(dict.fromkeys(d for r in reactions for _, d in r.terms))
    rows, conds = {}, {}
    for d in dirs:
        rows[d], conds[d] = _run_row(base, d)
    notes = []
    no_spec = [d for d in dirs if conds[d] is None]
    if no_spec:
        notes.append(L(f"spec.json が読めない計算は条件の比べから外しました: {', '.join(no_spec)}", f"runs without a readable spec.json are left out of the settings comparison: {', '.join(no_spec)}"))
    differing_all, partial_all, table_all = _diff_keys(conds, dirs)

    rx = []
    for r in reactions:
        es = [rows[d]["energy_ev"] for _, d in r.terms]
        de = float(sum(nu * e for (nu, _), e in zip(r.terms, es))) if all(e is not None for e in es) else None
        gs = [rows[d]["code_g_ev"] for _, d in r.terms]
        dg = float(sum(nu * g for (nu, _), g in zip(r.terms, gs))) if all(g is not None for g in gs) else None
        bal: Counter = Counter()
        comp_known = all(rows[d]["composition"] for _, d in r.terms)
        for nu, d in r.terms:
            for el, n in rows[d]["composition"].items():
                bal[el] += nu * n
        imb = {el: v for el, v in sorted(bal.items()) if abs(v) > BALANCE_TOL}
        differing, partial, _ = _diff_keys(conds, [d for _, d in r.terms])
        rx.append({"name": r.name, "terms": r.text(), "terms_list": [{"nu": nu, "dir": d} for nu, d in r.terms],
                   "delta_e_ev": de, "delta_e_kj_mol": de * EV_KJ_MOL if de is not None else None,
                   "delta_e_kcal_mol": de * EV_KCAL_MOL if de is not None else None, "delta_g_code_ev": dg,
                   "delta_g_code_kj_mol": dg * EV_KJ_MOL if dg is not None else None,
                   "balanced": comp_known and not imb, "composition_known": comp_known,
                   "imbalance": " ".join(f"{el}:{v:+g}" for el, v in imb.items()) if comp_known else L("組成が分からない計算がある", "composition unknown for a run"),
                   "differing": differing, "partial": partial, "n_differing": len(differing), "n_partial": len(partial)})

    res = CompareResult(base=str(base), runs=[rows[d] for d in dirs], reactions=rx, conditions=table_all,
                        differing=differing_all, partial=partial_all, notes=notes)
    _write(base, res)
    _plot(base, res)
    # the directory is copied to other machines: paths in the file are relative to it
    rel = lambda v: _relative_to(v, base)
    (base / "compare_summary.json").write_text(json.dumps({
        "runs": rel(res.runs), "reactions": res.reactions, "differing": res.differing, "partial": res.partial,
        "n_differing": len(res.differing), "n_partial": len(res.partial), "files": rel(res.files), "figures": rel(res.figures), "notes": res.notes,
        "units": {"energy": "eV", "kj_mol_per_ev": EV_KJ_MOL, "kcal_mol_per_ev": EV_KCAL_MOL},
        "rule": L("違う = その項目を持つ計算のあいだで値が 2 種類以上。一部だけ = その項目を持たない計算がある",
                  "differing = two or more distinct values among runs that have the setting; partial = some runs do not have the setting"),
    }, ensure_ascii=False, indent=2, default=str) + "\n", encoding="utf-8")
    res.files["summary"] = str(base / "compare_summary.json")
    return res


def _relative_to(obj, base: Path):
    if isinstance(obj, dict):
        return {k: _relative_to(v, base) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_relative_to(v, base) for v in obj]
    if isinstance(obj, str) and obj and Path(obj).is_absolute():
        try:
            return Path(obj).resolve().relative_to(Path(base).resolve()).as_posix()
        except ValueError:
            return obj
    return obj


def _cell(v) -> str:
    if v is None:
        return ""
    return json.dumps(v, ensure_ascii=False) if isinstance(v, (dict, list)) else str(v)


def _write(base: Path, res: CompareResult) -> None:
    cols = ["dir", "code", "task", "natoms", "formula", "energy_ev", "energy_source", "code_g_ev", "code_g_label", "note"]
    extra = res.differing + [k for k in res.partial if k not in res.differing]
    with open(base / "compare_runs.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(cols + extra)
        for r in res.runs:
            w.writerow([_cell(r[c]) for c in cols] + [_cell(res.conditions.get(k, {}).get(r["dir"])) for k in extra])
    with open(base / "compare_reactions.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["name", "terms", "delta_e_ev", "delta_e_kj_mol", "delta_e_kcal_mol", "delta_g_code_ev", "balanced", "imbalance",
                    "n_differing", "differing", "n_partial"])
        for x in res.reactions:
            w.writerow([x["name"], x["terms"], _cell(x["delta_e_ev"]), _cell(x["delta_e_kj_mol"]), _cell(x["delta_e_kcal_mol"]),
                        _cell(x["delta_g_code_ev"]), x["balanced"], x["imbalance"] if not x["balanced"] else "", x["n_differing"],
                        " ".join(x["differing"]), x["n_partial"]])
    dirs = [r["dir"] for r in res.runs]
    with open(base / "compare_conditions.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["setting", *dirs, "differs", "partial"])
        for k, vals in res.conditions.items():
            w.writerow([k, *[_cell(vals.get(d)) for d in dirs], k in res.differing, k in res.partial])
    res.files.update(runs=str(base / "compare_runs.csv"), reactions=str(base / "compare_reactions.csv"), conditions=str(base / "compare_conditions.csv"))


def _plot(base: Path, res: CompareResult) -> None:
    got = [x for x in res.reactions if x["delta_e_ev"] is not None]
    if not got:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(6, 0.5 + 0.45 * len(got) + 1.2))
    y = np.arange(len(got))
    v = [x["delta_e_ev"] for x in got]
    ax.barh(y, v, color=["tab:blue" if x["balanced"] else "tab:gray" for x in got])
    ax.set_yticks(y); ax.set_yticklabels([x["name"] for x in got]); ax.invert_yaxis()
    ax.axvline(0, color="k", lw=0.8)
    for yi, x in zip(y, got):
        ax.text(x["delta_e_ev"], yi, f" {x['delta_e_ev']:+.3f} eV ({x['delta_e_kj_mol']:+.1f} kJ/mol, {x['n_differing']} settings differ)",
                va="center", ha="left" if x["delta_e_ev"] >= 0 else "right", fontsize=8)
    ax.set_xlabel(r"$\Delta E = \sum \nu E$ [eV]  (gray: composition not balanced)"); ax.grid(alpha=0.3, axis="x")
    p = base / "compare_energy.png"; fig.savefig(p, dpi=110, bbox_inches="tight"); plt.close(fig)
    res.figures["compare_energy"] = str(p)


__all__ = ["Reaction", "CompareError", "CompareResult", "parse_compare", "load_compare_file", "analyze_compare", "conditions_of"]
