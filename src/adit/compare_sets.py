
from __future__ import annotations

from adit.errors import AditValueError
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np

from adit import batch
from adit.lang import L
from adit.spec import AtomsData, CalculationSpec, Structure

KINDS = ("adsorption", "reaction", "solvation")
COMPARE_FILE = "compare.json"
_NAME = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")
_MEMBER_KEYS = {"name", "structure", "spec", "charge", "multiplicity", "fixed_atoms", "nu"}


class CompareSetError(AditValueError):
    pass


def load_set(path: Path | str) -> tuple[dict, Path]:
    p = Path(path).expanduser()
    try:
        return json.loads(p.read_text(encoding="utf-8")), p.parent
    except FileNotFoundError as ex:
        raise CompareSetError(L(f"組の定義のファイルがありません: {p}", f"set file not found: {p}")) from ex
    except json.JSONDecodeError as ex:
        raise CompareSetError(L(f"{p} を JSON として読めません: {ex}", f"cannot read {p} as JSON: {ex}")) from ex


def _flat(d, prefix: str = "") -> dict:
    out = {}
    for k, v in (d or {}).items():
        key = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, dict) and v:
            out.update(_flat(v, key))
        else:
            out[key] = v
    return out


def _member_structure(base: CalculationSpec, m: dict, where: Path, name: str) -> Structure:
    from adit.structure import StructureError, from_file

    if not isinstance(m, dict):
        raise CompareSetError(L(f"{name}: {{…}} の形で書いてください", f"{name}: must be an object {{...}}"))
    bad = sorted(set(m) - _MEMBER_KEYS)
    if bad:
        raise CompareSetError(L(f"{name}: 使えない項目があります: {bad} (使えるのは {sorted(_MEMBER_KEYS)})", f"{name}: unknown items {bad} (allowed: {sorted(_MEMBER_KEYS)})"))
    if m.get("spec") and m.get("structure"):
        raise CompareSetError(L(f"{name}: structure と spec はどちらか一方にしてください", f"{name}: give either structure or spec, not both"))
    if m.get("spec"):
        p = Path(m["spec"]).expanduser()
        p = p if p.is_absolute() else where / p
        try:
            st = CalculationSpec.load(p).structure
        except Exception as ex:
            raise CompareSetError(L(f"{name}: {p} を読めません: {ex}", f"{name}: cannot read {p}: {ex}")) from ex
    elif m.get("structure") in (None, "", "base"):
        st = base.structure
    else:
        p = Path(m["structure"]).expanduser()
        p = p if p.is_absolute() else where / p
        try:
            atoms = from_file(p)
        except StructureError as ex:
            raise CompareSetError(f"{name}: {ex}") from ex
        st = Structure(source="file", source_ref=str(p), atoms=AtomsData.from_ase(atoms))
    upd = {k: m[k] for k in ("charge", "multiplicity", "fixed_atoms") if k in m}
    upd["velocities"] = None
    from pydantic import ValidationError as PydanticError

    from adit.validate_types import friendly_pydantic
    try:
        return Structure.model_validate({**st.model_dump(mode="json"), **upd})
    except PydanticError as ex:
        raise CompareSetError(f"{name}: " + friendly_pydantic(ex)) from ex
    except Exception as ex:
        raise CompareSetError(L(f"{name}: {ex}", f"{name}: {ex}")) from ex


def _member_spec(base: CalculationSpec, st: Structure, name: str, role: str, overrides: dict | None = None) -> tuple[CalculationSpec, list[str]]:
    from pydantic import ValidationError as PydanticError

    from adit.stages import _merge
    from adit.validate_types import friendly_pydantic

    data = base.model_dump(mode="json")
    data["structure"] = st.model_dump(mode="json")
    data["handoff"] = None
    changed = []
    if not st.periodic:
        if data.get("kpoints") is not None:
            data["kpoints"] = None
            changed.append(L("非周期なので k 点を使いません", "non-periodic, so no k-points"))
        if data["task"].get("relax_cell", "no") != "no":
            data["task"]["relax_cell"] = "no"
            changed.append(L("非周期なので格子を動かしません (task.relax_cell = no)", "non-periodic, so the cell is not relaxed (task.relax_cell = no)"))
    if overrides:
        data = _merge(data, overrides)
    data["meta"] = dict(data.get("meta") or {}, comment=(str((data.get("meta") or {}).get("comment", "")) + f" compare-set {role}:{name}").strip())
    try:
        return CalculationSpec.model_validate(data), changed
    except PydanticError as ex:
        raise CompareSetError(f"{name}: " + friendly_pydantic(ex)) from ex


def _check_overrides(ov, name: str) -> dict:
    if ov is None:
        return {}
    if not isinstance(ov, dict) or set(ov) - {"method"} or not isinstance(ov.get("method", {}), dict):
        raise CompareSetError(L(f"{name}: 上書きできるのは {{\"method\": {{…}}}} だけです", f"{name}: only {{\"method\": {{...}}}} can be overridden"))
    if "code" in ov.get("method", {}):
        raise CompareSetError(L(f"{name}: 組の中で計算コードは変えられません", f"{name}: the code cannot change within a set"))
    return ov


def plan_set(spec: CalculationSpec, data: dict, where: Path) -> dict:
    if not isinstance(data, dict) or data.get("kind") not in KINDS:
        raise CompareSetError(L(f"kind は {' / '.join(KINDS)} のどれかにしてください", f"kind must be one of {' / '.join(KINDS)}"))
    kind = data["kind"]
    members: list[tuple[str, str, CalculationSpec, list[str]]] = []
    intended: set[str] = set()
    if kind == "adsorption":
        roles = ("slab", "molecule", "adsorbed")
        missing = [r for r in roles if r not in data]
        if missing:
            raise CompareSetError(L(f"吸着の組には slab / molecule / adsorbed の 3 つが要ります (無いもの: {missing})",
                                    f"an adsorption set needs slab / molecule / adsorbed (missing: {missing})"))
        box = data.get("molecule_box", "as_is")
        if box not in ("as_is", "slab_cell"):
            raise CompareSetError(L("molecule_box は as_is か slab_cell です", "molecule_box must be as_is or slab_cell"))
        sts = {r: _member_structure(spec, data[r], where, r) for r in roles}
        atom_sets = [sts[r].atoms.model_dump_json() for r in roles]
        if len(set(atom_sets)) == 1:
            raise CompareSetError(L("スラブ・分子・吸着した構造がすべて同じです。3 つの構造をそれぞれ指定してください",
                                    "the slab, molecule and adsorbed structures are identical; provide each of the three structures"))
        if box == "slab_cell":
            slab = sts["slab"].atoms.to_ase()
            if not any(slab.pbc):
                raise CompareSetError(L("molecule_box = slab_cell には、周期のあるスラブが要ります", "molecule_box = slab_cell needs a periodic slab"))
            mol = sts["molecule"].atoms.to_ase()
            center = 0.5 * np.asarray(slab.cell).sum(axis=0)
            mol.positions = mol.positions - mol.positions.mean(axis=0) + center
            mol.cell, mol.pbc = slab.cell, slab.pbc
            sts["molecule"] = sts["molecule"].model_copy(update={"atoms": AtomsData.from_ase(mol)})
        for r in roles:
            s, ch = _member_spec(spec, sts[r], r, r)
            members.append((r, r, s, ch))
        reactions = [{"name": data.get("name") or "adsorption", "terms": [{"dir": "adsorbed", "nu": 1}, {"dir": "slab", "nu": -1}, {"dir": "molecule", "nu": -1}]}]
    elif kind == "reaction":
        items = data.get("members")
        if not isinstance(items, list) or len(items) < 2:
            raise CompareSetError(L("反応の組には members に 2 つ以上の計算が要ります", "a reaction set needs at least two members"))
        terms = []
        for i, it in enumerate(items, 1):
            name = str((it or {}).get("name") or "") if isinstance(it, dict) else ""
            if not _NAME.match(name):
                raise CompareSetError(L(f"{i} 番目の name は英数字と _ . - で書いてください (ディレクトリの名前になります): {name!r}",
                                        f"name of member {i} must be letters, digits, _ . - (it becomes a directory name): {name!r}"))
            if any(t["dir"] == name for t in terms):
                raise CompareSetError(L(f"{name}: 同じ name が 2 つあります (ディレクトリの名前になるので、別の名前にしてください)",
                                        f"{name}: this name is used twice (it becomes a directory name, so give different names)"))
            try:
                nu = float(it.get("nu"))
            except (TypeError, ValueError):
                nu = 0.0
            if nu == 0 or not np.isfinite(nu):
                raise CompareSetError(L(f"{name}: 係数 nu (生成物は正、反応物は負) を 0 以外の有限の数で書いてください", f"{name}: give a non-zero finite coefficient nu (products positive, reactants negative)"))
            s, ch = _member_spec(spec, _member_structure(spec, it, where, name), name, "member")
            members.append((name, "member", s, ch))
            terms.append({"dir": name, "nu": nu})
        reactions = [{"name": data.get("name") or "reaction", "terms": terms}]
    else:  # solvation
        if "solvent" not in data:
            raise CompareSetError(L("溶媒和の組には solvent (溶媒側の method の上書き) が要ります", "a solvation set needs solvent (method overrides for the solvated run)"))
        st = _member_structure(spec, {"structure": data.get("structure", "base")}, where, "structure")
        gas_ov, sol_ov = _check_overrides(data.get("gas"), "gas"), _check_overrides(data["solvent"], "solvent")
        g, gch = _member_spec(spec, st, "gas", "gas", gas_ov)
        s, sch = _member_spec(spec, st, "solvent", "solvent", sol_ov)
        if g.method == s.method:
            raise CompareSetError(L("溶媒側と気相側の method が同じです (solvent の上書きで何も変わっていません)", "the solvated and gas-phase methods are identical (the solvent overrides change nothing)"))
        members += [("gas", "gas", g, gch), ("solvent", "solvent", s, sch)]
        fg, fs = _flat({"method": g.model_dump(mode="json")["method"]}), _flat({"method": s.model_dump(mode="json")["method"]})
        intended = {k for k in set(fg) | set(fs) if fg.get(k) != fs.get(k)}
        reactions = [{"name": data.get("name") or "solvation", "terms": [{"dir": "solvent", "nu": 1}, {"dir": "gas", "nu": -1}]}]
    return {"kind": kind, "members": members, "reactions": reactions, "intended": intended}


def conditions(spec: CalculationSpec) -> dict:
    from adit.analysis.compare import conditions_of

    c = conditions_of(spec)
    if spec.kpoints is not None and spec.structure.periodic:
        c["kpoints.resolved_mesh"] = json.dumps(list(spec.kpoints.resolved_mesh(spec.structure.atoms.cell)))
    return c


def differences(members, intended: set[str]) -> tuple[list[dict], list[dict]]:
    conds = {d: conditions(s) for d, _, s, _ in members}
    keys = sorted({k for c in conds.values() for k in c})
    diff, partial = [], []
    for k in keys:
        vals = {d: c[k] for d, c in conds.items() if k in c}
        row = {"setting": k, "values": vals, "intended": k in intended}
        if len({json.dumps(v, sort_keys=True, default=str) for v in vals.values()}) >= 2:
            diff.append(row)
        if len(vals) < len(conds):
            partial.append({**row, "missing_in": [d for d in conds if d not in vals]})
    return diff, partial


def balance(members, reactions) -> list[dict]:
    specs = {d: s for d, _, s, _ in members}
    out = []
    for r in reactions:
        comp, charge = Counter(), 0.0
        for t in r["terms"]:
            s = specs[t["dir"]]
            for el, n in Counter(s.structure.atoms.symbols).items():
                comp[el] += t["nu"] * n
            charge += t["nu"] * s.structure.charge
        out.append({"name": r["name"], "imbalance": {el: v for el, v in sorted(comp.items()) if abs(v) > 1e-9}, "charge_imbalance": charge})
    return out


def _formula(s: CalculationSpec) -> str:
    return "".join(f"{el}{'' if n == 1 else n}" for el, n in sorted(Counter(s.structure.atoms.symbols).items()))


def _val(v) -> str:
    return v if isinstance(v, str) else json.dumps(v, default=str)


def write_compare_set(spec: CalculationSpec, cfg, out_dir: Path | str, data: dict, where: Path | str = ".", *, overwrite: bool = False) -> list[Path]:
    from adit import __version__

    out = Path(out_dir).expanduser()
    plan = plan_set(spec, data, Path(where))
    members, reactions = plan["members"], plan["reactions"]
    diff, partial = differences(members, plan["intended"])
    bal = balance(members, reactions)
    rx_text = " ".join(f"{t['nu']:+g}·{t['dir']}" for t in reactions[0]["terms"])
    items = []
    for d, role, s, changed in members:
        mine = [x for x in diff if d in x["values"]]
        lines = [L("== 比べる計算の組 ==", "== Set of runs to compare =="),
                 L(f"  この計算は {plan['kind']} の組の {d} です (組の説明と式は ../README.txt、ΔE = {rx_text})。",
                   f"  This run is {d} of a {plan['kind']} set (see ../README.txt; ΔE = {rx_text})."),
                 *[f"  - {c}" for c in changed]]
        if mine:
            lines.append(L(f"  組の中で値が違う条件: {len(mine)} 項目 (一覧は ../compare.json の differences)", f"  settings that differ within the set: {len(mine)} (see differences in ../compare.json)"))
        items.append(batch.Item(d, s, lines))
    dirs = batch.write_items(items, cfg, out, overwrite=overwrite)
    batch.write_json(out / COMPARE_FILE, {
        "kind": plan["kind"], "generated_by": f"adit {__version__}", "reactions": reactions,
        "members": [{"dir": d, "role": role, "formula": _formula(s), "natoms": len(s.structure.atoms.symbols), "periodic": s.structure.periodic,
                     "charge": s.structure.charge, "multiplicity": s.structure.multiplicity, "structure_source": s.structure.source_ref,
                     "kpoints_mesh": list(s.kpoints.resolved_mesh(s.structure.atoms.cell)) if s.kpoints is not None and s.structure.periodic else None,
                     "changed_for_this_run": changed} for d, role, s, changed in members],
        "differences": diff, "partial": partial, "balance": bal,
        "rule": L("違う = その項目を持つ計算のあいだで値が 2 種類以上。一部だけ = その項目を持たない計算がある。intended = 違えるのが目的の項目 (溶媒和の上書き)",
                  "differing = two or more distinct values among runs that have the setting; partial = some runs lack it; intended = overridden on purpose (solvation)"),
    })
    lines = [L(f"ADIT {__version__} が生成した、比べる計算の組です (種類: {plan['kind']}、計算コード: {spec.method.code})",
               f"Set of runs to compare, generated by ADIT {__version__} (kind: {plan['kind']}, code: {spec.method.code})"), "",
             L("  ディレクトリ    役割        組成      原子数  周期  電荷  多重度  k 点", "  directory       role        formula   atoms   pbc   chg   mult    k-mesh")]
    for d, role, s, _ in members:
        mesh = "x".join(str(k) for k in s.kpoints.resolved_mesh(s.structure.atoms.cell)) if s.kpoints is not None and s.structure.periodic else "-"
        lines.append(f"  {d:<15} {role:<11} {_formula(s):<9} {len(s.structure.atoms.symbols):<7} {'yes' if s.structure.periodic else 'no':<5} "
                     f"{s.structure.charge:<5} {s.structure.multiplicity:<7} {mesh}")
    lines += ["", L(f"比べる式: ΔE = {rx_text}   (係数 ν は生成物が正、反応物が負)", f"Difference: ΔE = {rx_text}   (ν: products positive, reactants negative)")]
    for b in bal:
        if b["imbalance"] or b["charge_imbalance"]:
            lines.append(L(f"  組成または電荷が釣り合っていません: Σν·組成 = {b['imbalance']}、Σν·電荷 = {b['charge_imbalance']:g} (この差は反応のエネルギーではありません)",
                           f"  composition or charge is not balanced: Σν·composition = {b['imbalance']}, Σν·charge = {b['charge_imbalance']:g} (this difference is not a reaction energy)"))
        else:
            lines.append(L("  組成と電荷は釣り合っています", "  composition and charge are balanced"))
    lines += ["", L(f"== 組の中で値が違う条件 ({len(diff)} 項目。ADIT は止めません。比べてよいかは判断しません) ==",
                    f"== Settings that differ within the set ({len(diff)}; ADIT does not stop, and does not judge whether they can be compared) ==")]
    for x in diff:
        tag = L(" (違えるのが目的の項目)", " (intended)") if x["intended"] else ""
        lines.append(f"  {x['setting']}{tag}: " + ", ".join(f"{d} = {_val(v)}" for d, v in x["values"].items()))
    if not diff:
        lines.append(L("  (なし)", "  (none)"))
    only = [x for x in partial if x["setting"] not in {y["setting"] for y in diff}]
    if only:
        lines.append(L(f"  一部の計算にだけある項目: {len(only)} 項目 (例: {', '.join(x['setting'] for x in only[:6])}。一覧は compare.json の partial)",
                       f"  settings present in only some runs: {len(only)} (e.g. {', '.join(x['setting'] for x in only[:6])}; see partial in compare.json)"))
    lines += ["", L("== 実行したあと ==", "== After running =="),
              L("  adit-analyze <このディレクトリ> --compare   (compare.json を読み、ΔE と条件の違いの表を書きます)",
                "  adit-analyze <this directory> --compare   (reads compare.json and writes ΔE and the table of differing settings)")]
    batch.write_top(out, cfg, spec, dirs, lines, "組の計算", "runs of the set")
    return dirs


__all__ = ["CompareSetError", "KINDS", "load_set", "plan_set", "write_compare_set", "differences", "balance"]
