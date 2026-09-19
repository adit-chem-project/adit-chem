
from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from pydantic import TypeAdapter

from adit.builder.model import (Adsorb, Base, Box, ClusterRef, Fix, MoleculeRef, PolymerRef, Recipe, Remove, Selection, Slab,
                               SolventLayer, Solvate, Step, Substitute, Supercell, TwoD, Vacuum, interface_steps as make_interface_steps,
                               op_label)
from adit.lang import L
from adit.mixture import Component, MixtureError, mixture_text, parse_mixture_text
from adit.structure import StructureError, pretty_formula

NEW_BASES = ("2d", "cluster", "polymer")
SLOW_BASES = ("polymer",)
TWOD_KINDS = [("graphene", "グラフェン型 (C₂、BN)", "Graphene-type (C₂, BN)"), ("mx2", "MX₂ 型 (MoS₂ など)", "MX₂ sheet (MoS₂ etc.)"),
              ("nanoribbon", "ナノリボン", "Nanoribbon"), ("nanotube", "ナノチューブ", "Nanotube")]
CLUSTER_KINDS = [("icosahedron", "正二十面体", "Icosahedron"), ("decahedron", "十面体", "Decahedron"), ("octahedron", "八面体", "Octahedron"),
                 ("wulff", "Wulff 形 (面のエネルギーから)", "Wulff shape (from surface energies)")]

ELECTROLYTE_EXAMPLE = [Component(kind="smiles", ref="O=C1OCCO1", count=32, label="C3H4O3"),
                       Component(kind="smiles", ref="[Li+]", count=1, charge=1, label="Li+"),
                       Component(kind="smiles", ref="F[P-](F)(F)(F)(F)F", count=1, charge=-1, label="PF6-")]
ADD_ORDER = ["slab", "solvent_layer", "fix", "supercell", "remove", "substitute", "adsorb", "vacuum", "box", "solvate"]
INTERFACE = "interface"

_STEP = TypeAdapter(Step)
_STEP_NO = re.compile(r"(?:手順|step) (\d+) \(")


class RecipeFormError(StructureError):
    pass


def default_step(op: str):
    return {"supercell": lambda: Supercell(repeat=(2, 2, 1)), "slab": lambda: Slab(miller=(1, 0, 0)), "vacuum": Vacuum, "box": Box,
            "adsorb": lambda: Adsorb(molecule=MoleculeRef(kind="preset", ref="H2O"), xy=(0.0, 0.0)),
            "remove": lambda: Remove(count=1), "substitute": lambda: Substitute(count=1, to="Al"),
            "solvent_layer": lambda: SolventLayer(components=[Component(ref="H2O", count=32, label="H2O")]),
            "solvate": lambda: Solvate(components=[Component(ref="H2O", count=0, label="H2O")]),
            "fix": lambda: Fix(bottom_layers=1)}[op]()


def interface_steps(base_is_slab: bool) -> list:
    return make_interface_steps(base_is_slab, ELECTROLYTE_EXAMPLE)


INTERFACE_NOTE = ("溶液の成分は例です。個数と密度 (または厚み) を決めてから「作る」を押してください。板の断面は、成分が周期境界で重ならない広さまで自動的に広がります。イオンの個数は「濃度から個数」でも決められます",
                  "the solution components are examples. Set their counts and the density (or thickness), then press Build. The slab cross-section is enlarged automatically until each component clears its periodic images. You can also set ion counts with \"Count from concentration\"")


def failed_step(message: str) -> int | None:
    m = _STEP_NO.search(message or "")
    k = int(m.group(1)) if m else 0
    return k if k >= 1 else None


_SUP = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")


def pretty_name(name: str) -> str:
    from ase.data import chemical_symbols

    m = re.fullmatch(r"([A-Z][A-Za-z0-9]*?)(\d*)([+-])", name)
    if not m:
        return pretty_formula(name)
    body, digits, sign = m.groups()
    if digits and body in chemical_symbols:
        return body + (digits + sign).translate(_SUP)
    return pretty_formula(body + digits) + sign.translate(_SUP)


def _fmt(x) -> str:
    return "" if x is None else f"{float(x):.15g}"


def describe_selection(sel: Selection) -> str:
    parts = []
    if sel.elements:
        parts.append(" ".join(sel.elements))
    if sel.z_min is not None or sel.z_max is not None:
        parts.append(f"z {_fmt(sel.z_min)}〜{_fmt(sel.z_max)} Å")
    if sel.indices:
        parts.append(L(f"番号 {len(sel.indices)} 個", f"{len(sel.indices)} indices"))
    return ", ".join(parts)


def _comps(s) -> str:
    return ", ".join(f"{pretty_name(c.name())} × {c.count if c.count else L('自動', 'auto')}" for c in s.components)


def summary(s) -> str:
    op = s.op
    if op == "supercell":
        repeat = " × ".join(str(x) for x in s.repeat) if s.repeat else L("変換行列", "matrix")
        return L(f"溶液に合わせて自動調整 (下限 {repeat})", f"auto-sized for solution (minimum {repeat})") if s.fit_components else repeat
    if op == "slab":
        return L(f"({' '.join(map(str, s.miller))})、{s.layers} 層", f"({' '.join(map(str, s.miller))}), {s.layers} layer{'s' if s.layers != 1 else ''}")
    if op == "vacuum":
        return f"{'abc'[s.axis]}, {s.thickness:g} Å"
    if op == "box":
        return L(f"余白 {s.padding:g} Å", f"padding {s.padding:g} Å")
    if op in ("remove", "substitute"):
        what = describe_selection(s.where) or L("全原子", "all atoms")
        n = f"{s.count}" if s.count is not None else f"{s.fraction:g}"
        amount = L(f"{what} から {n}", f"{n} of {what}")
        return amount if op == "remove" else f"{amount} → {s.to}"
    if op == "adsorb":
        where = s.site if s.site is not None else L(f"原子 {s.above_atom + 1} の上", f"above atom {s.above_atom + 1}") if s.above_atom is not None else f"xy ({s.xy[0]:g}, {s.xy[1]:g})"
        mol = pretty_formula(s.molecule.ref) if s.molecule.kind == "preset" else s.molecule.ref
        return f"{mol}, {where}, {s.height:g} Å"
    if op in ("solvent_layer", "solvate"):
        return _comps(s)
    if op == "fix":
        parts = [L(f"下から {s.bottom_layers} 層", f"bottom {s.bottom_layers} layer{'s' if s.bottom_layers != 1 else ''}")] if s.bottom_layers else []
        d = describe_selection(s.where)
        if d:
            parts.append(d)
        return ", ".join(parts) or L("全原子", "all atoms")
    return ""


def step_fields(k: int, s) -> dict[str, str]:
    p = f"st{k}_"
    f: dict[str, str] = {p + "op": s.op}

    def sel(w: Selection) -> None:
        f.update({p + "elements": " ".join(w.elements), p + "zmin": _fmt(w.z_min), p + "zmax": _fmt(w.z_max)})

    if s.op == "supercell":
        f[p + "mode"] = "matrix" if s.matrix is not None else "repeat"
        for i, v in enumerate(s.repeat or (1, 1, 1)):
            f[f"{p}r{i}"] = str(v)
        m = s.matrix or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        for i in range(3):
            for j in range(3):
                f[f"{p}m{i}{j}"] = str(m[i][j])
    elif s.op == "slab":
        for c, v in zip("hkl", s.miller):
            f[p + c] = str(v)
        f.update({p + "layers": str(s.layers), p + "vacuum": _fmt(s.vacuum), p + "term": str(s.termination)})
    elif s.op == "vacuum":
        f.update({p + "axis": str(s.axis), p + "thickness": _fmt(s.thickness)})
    elif s.op == "box":
        f.update({p + "padding": _fmt(s.padding), p + "maxm": str(s.max_multiple)})
    elif s.op in ("remove", "substitute"):
        sel(s.where)
        f.update({p + "pmode": "fraction" if s.fraction is not None else "count", p + "count": str(s.count if s.count is not None else 1),
                  p + "fraction": _fmt(s.fraction if s.fraction is not None else 0.1), p + "seed": str(s.seed)})
        if s.op == "substitute":
            f[p + "to"] = s.to
    elif s.op == "adsorb":
        place = "site" if s.site is not None else "above_atom" if s.above_atom is not None else "xy"
        xy = s.xy or (0.0, 0.0)
        f.update({p + "kind": s.molecule.kind, p + "ref": s.molecule.ref, p + "place": place, p + "site": s.site or "ontop",
                  p + "atom": str((s.above_atom or 0) + 1), p + "x": _fmt(xy[0]), p + "y": _fmt(xy[1]), p + "height": _fmt(s.height),
                  p + "down": str(s.down_atom + 1)})
    elif s.op in ("solvent_layer", "solvate"):
        f.update({p + "comps": mixture_text(s.components), p + "mind": _fmt(s.min_distance), p + "seed": str(s.seed),
                  p + "density": _fmt(s.density_g_cm3)})
        if s.op == "solvent_layer":
            f.update({p + "tmode": "density" if s.thickness is None else "thickness", p + "thickness": _fmt(s.thickness if s.thickness is not None else 20.0),
                      p + "gap": _fmt(s.gap), p + "vacuum": _fmt(s.vacuum), p + "conc": "1",
                      p + "conc_row": str(next((i for i, c in enumerate(s.components) if c.charge), 0))})
        else:
            f[p + "padding"] = _fmt(s.padding)
    elif s.op == "fix":
        f[p + "layers"] = str(s.bottom_layers)
        sel(s.where)
    return f


def _num(f: dict, name: str, what: str, *, integer: bool = False, optional: bool = False):
    v = (f.get(name) or "").strip()
    if not v:
        if optional:
            return None
        raise ValueError(L(f"{what} が空です", f"{what} is empty"))
    try:
        x = float(v)
    except ValueError as ex:
        raise ValueError(L(f"{what} を数値として読めません: {v!r}", f"{what} is not a number: {v!r}")) from ex
    if not np.isfinite(x):
        raise ValueError(L(f"{what} には有限の数を入れてください: {v!r}", f"{what} must be a finite number: {v!r}"))
    if integer:
        if x != int(x):
            raise ValueError(L(f"{what} には整数を入れてください: {v!r}", f"{what} must be an integer: {v!r}"))
        return int(x)
    return x


def _merge(k: int, s, f: dict):
    p = f"st{k}_"
    g = lambda n: (f.get(p + n) or "").strip()  # noqa: E731
    num = lambda n, what, **kw: _num(f, p + n, what, **kw)  # noqa: E731

    def where(old: Selection) -> dict:
        return {"elements": [e for e in re.split(r"[,\s、]+", g("elements")) if e], "indices": list(old.indices),
                "z_min": num("zmin", L("z の下限", "z min"), optional=True), "z_max": num("zmax", L("z の上限", "z max"), optional=True)}

    def comps() -> list:
        try:
            return [c.model_dump() for c in parse_mixture_text(f.get(p + "comps") or "")]
        except MixtureError as ex:
            raise ValueError(str(ex)) from ex

    op = s.op
    u: dict = {}
    if op == "supercell":
        if g("mode") == "matrix":
            u = {"repeat": None, "matrix": [[num(f"m{i}{j}", L("変換行列", "matrix"), integer=True) for j in range(3)] for i in range(3)]}
        else:
            u = {"repeat": tuple(num(f"r{i}", L("繰り返し", "repeat"), integer=True) for i in range(3)), "matrix": None}
    elif op == "slab":
        u = {"miller": tuple(num(c, L("ミラー指数", "Miller index"), integer=True) for c in "hkl"), "layers": num("layers", L("層の数", "layers"), integer=True),
             "vacuum": num("vacuum", L("真空", "vacuum")), "termination": num("term", L("終端", "termination"), integer=True)}
    elif op == "vacuum":
        u = {"axis": num("axis", L("軸", "axis"), integer=True), "thickness": num("thickness", L("隙間", "gap"))}
    elif op == "box":
        u = {"padding": num("padding", L("余白", "padding")), "max_multiple": num("maxm", L("倍数の上限", "max. multiple"), integer=True)}
    elif op in ("remove", "substitute"):
        frac = g("pmode") == "fraction"
        u = {"where": where(s.where), "count": None if frac else num("count", L("個数", "count"), integer=True),
             "fraction": num("fraction", L("割合", "fraction")) if frac else None, "seed": num("seed", L("乱数の種", "seed"), integer=True)}
        if op == "substitute":
            u["to"] = g("to")
    elif op == "adsorb":
        place = g("place") or "xy"
        u = {"molecule": {"kind": g("kind") or "preset", "ref": g("ref")},
             "site": g("site") if place == "site" else None,
             "above_atom": num("atom", L("原子の番号", "atom number"), integer=True) - 1 if place == "above_atom" else None,
             "xy": (num("x", "x"), num("y", "y")) if place == "xy" else None,
             "height": num("height", L("高さ", "height")), "down_atom": num("down", L("下に向ける原子", "atom facing down"), integer=True) - 1}
    elif op in ("solvent_layer", "solvate"):
        u = {"components": comps(), "min_distance": num("mind", L("分子間の最短距離", "min. distance")), "seed": num("seed", L("乱数の種", "seed"), integer=True),
             "density_g_cm3": num("density", L("密度", "density"))}
        if op == "solvent_layer":
            u.update(thickness=num("thickness", L("厚み", "thickness")) if g("tmode") == "thickness" else None,
                     gap=num("gap", L("隙間", "gap")), vacuum=num("vacuum", L("真空", "vacuum")))
        else:
            u["padding"] = num("padding", L("余白", "padding"))
    elif op == "fix":
        u = {"bottom_layers": num("layers", L("下から数えた層", "bottom layers"), integer=True), "where": where(s.where)}
    return type(s).model_validate({**s.model_dump(), **u})


def _first_msg(ex: Exception) -> str:
    errs = getattr(ex, "errors", None)
    if callable(errs):
        e = errs()[0] if errs() else {}
        return str(e.get("msg", ex)).removeprefix("Value error, ")
    return str(ex)


def raw_steps(f: dict) -> list:
    text = (f.get("recipe_steps") or "").strip() or "[]"
    try:
        data = json.loads(text)
        if not isinstance(data, list):
            raise ValueError("not a list")
    except ValueError as ex:
        raise RecipeFormError(L(f"組み立て手順の控え (recipe_steps) を読めません: {ex}", f"cannot read the stored build steps (recipe_steps): {ex}")) from ex
    out = []
    for k, d in enumerate(data, start=1):
        try:
            out.append(_STEP.validate_python(d))
        except ValueError as ex:
            raise RecipeFormError(L(f"手順 {k}: 読めません: {_first_msg(ex)}", f"step {k}: cannot read: {_first_msg(ex)}")) from ex
    return out


def step_count(f: dict) -> int:
    try:
        return len(json.loads((f.get("recipe_steps") or "").strip() or "[]"))
    except (ValueError, TypeError):
        return 0


def steps_from_form(f: dict) -> list:
    out = []
    for k, s in enumerate(raw_steps(f), start=1):
        if f.get(f"st{k}_op") == s.op:
            try:
                s = _merge(k, s, f)
            except ValueError as ex:
                raise RecipeFormError(L(f"手順 {k} ({op_label(s.op)}): {_first_msg(ex)}", f"step {k} ({op_label(s.op)}): {_first_msg(ex)}")) from ex
        out.append(s)
    return out


def write_steps(f: dict, steps: list) -> None:
    for key in [x for x in f if re.match(r"st\d+_", x)]:
        del f[key]
    f["recipe_steps"] = json.dumps([s.model_dump(mode="json", exclude_none=True) for s in steps], ensure_ascii=False)
    for k, s in enumerate(steps, start=1):
        f.update(step_fields(k, s))


def _compact(model, extra: dict) -> dict:
    d = model.model_dump(mode="json", exclude_defaults=True)
    full = model.model_dump(mode="json")
    for k in extra:
        if k not in d and k in full:
            d[k] = full[k]
    return d


def _extra(f: dict, src: str) -> dict:
    try:
        d = json.loads(f.get("base_extra") or "{}")
    except ValueError:
        return {}
    x = d.get(src) if isinstance(d, dict) else None
    return dict(x) if isinstance(x, dict) else {}


def _pos_int(f: dict, name: str, what: str, default: int) -> int:
    v = _num(f, name, what, integer=True, optional=True)
    return default if v is None else v


def base_ref(f: dict, src: str) -> dict:
    ex = _extra(f, src)
    d = dict(ex)
    s = lambda n, dflt="": (f.get(n) or dflt).strip()  # noqa: E731
    if src == "2d":
        k = s("td_kind", "graphene")
        d["kind"] = k
        v = _num(f, "td_vac", L("真空 (片側)", "vacuum (each side)"), optional=True)
        d["vacuum"] = 10.0 if v is None else v
        text = s("td_formula")
        if k in ("graphene", "mx2"):
            third = (d.get("size") or (1, 1, 1))[2]
            d["formula"] = text
            d["size"] = [_pos_int(f, "td_nx", L("繰り返し", "repeat"), 2), _pos_int(f, "td_ny", L("繰り返し", "repeat"), 2), int(third)]
        elif k == "nanotube":
            d.update(n=_pos_int(f, "td_n", "n", 6), m=_pos_int(f, "td_m", "m", 0), length=_pos_int(f, "td_len", L("長さ", "length"), 1), symbol=text or "C")
        else:
            d.update(n=_pos_int(f, "td_n", "n", 6), m=_pos_int(f, "td_m", "m", 0), ribbon_type=s("td_ribbon", "zigzag"),
                     saturated=s("td_sat", "yes") != "no")
        return _compact(TwoD.model_validate(d), ex)
    if src == "cluster":
        k = s("cl_kind", "icosahedron")
        d["kind"] = k; d["symbol"] = s("cl_el", "Cu")
        a = _num(f, "cl_a", L("格子定数 a", "lattice constant a"), optional=True)
        d["lattice_constant"] = a or None
        if k == "icosahedron":
            d["shells"] = _pos_int(f, "cl_shells", L("殻の数", "shells"), 3)
        elif k == "decahedron":
            d["p"], d["q"], d["r"] = _pos_int(f, "cl_p", "p", 2), _pos_int(f, "cl_q", "q", 2), _pos_int(f, "cl_r", "r", 0)
        elif k == "octahedron":
            d["length"], d["cutoff"] = _pos_int(f, "cl_len", L("稜の原子数", "atoms per edge"), 5), _pos_int(f, "cl_cut", L("頂点を切る層", "corner cut"), 0)
        else:
            d["size"], d["structure"] = _pos_int(f, "cl_size", L("原子数の目安", "target atoms"), 100), s("cl_struct", "fcc")
        return _compact(ClusterRef.model_validate(d), ex)
    unit = s("pl_unit")
    if not unit:
        raise ValueError(L("繰り返し単位の SMILES を入力してください (例 *CC*)", "enter the SMILES of the repeat unit (e.g. *CC*)"))
    d.update(unit=unit, n=_pos_int(f, "pl_n", L("重合度", "units"), 10), seed=_pos_int(f, "pl_seed", L("乱数の種", "seed"), 0))
    return _compact(PolymerRef.model_validate(d), ex)


def base_fields(src: str, ref: dict) -> dict[str, str]:
    f: dict[str, str] = {"base_extra": json.dumps({src: ref}, ensure_ascii=False)}
    if src == "2d":
        t = TwoD.model_validate(ref)
        f.update(td_kind=t.kind, td_formula=(t.symbol if t.symbol != "C" else "") if t.kind == "nanotube" else (t.formula if t.kind != "nanoribbon" else ""),
                 td_nx=str(t.size[0]), td_ny=str(t.size[1]), td_n=str(t.n), td_m=str(t.m), td_len=str(t.length), td_ribbon=t.ribbon_type,
                 td_sat="yes" if t.saturated else "no", td_vac=_fmt(t.vacuum))
    elif src == "cluster":
        c = ClusterRef.model_validate(ref)
        f.update(cl_kind=c.kind, cl_el=c.symbol, cl_shells=str(c.shells), cl_p=str(c.p), cl_q=str(c.q), cl_r=str(c.r), cl_len=str(c.length),
                 cl_cut=str(c.cutoff), cl_size=str(c.size), cl_struct=c.structure, cl_a=_fmt(c.lattice_constant) if c.lattice_constant else "0")
    elif src == "polymer":
        p = PolymerRef.model_validate(ref)
        f.update(pl_unit=p.unit, pl_n=str(p.n), pl_seed=str(p.seed))
    return f


BASE_DEFAULTS = {"td_kind": "graphene", "td_nx": "2", "td_ny": "2", "td_n": "6", "td_m": "0", "td_len": "1", "td_ribbon": "zigzag", "td_sat": "yes",
                 "td_vac": "10", "cl_kind": "icosahedron", "cl_el": "Cu", "cl_shells": "3", "cl_p": "2", "cl_q": "2", "cl_r": "0", "cl_len": "5",
                 "cl_cut": "0", "cl_size": "100", "cl_struct": "fcc", "cl_a": "0", "pl_n": "10", "pl_seed": "0", "recipe_steps": "[]",
                 "recipe_add": "slab"}

_SHA: dict[tuple, str] = {}


def base_from_form(f: dict) -> Base:
    from adit.web.forms import FormError, source_ref

    src = (f.get("source") or "preset").strip()
    if src == "fetch":
        from adit.web.forms import as_file_source
        try:
            f, _rec = as_file_source(f)
        except FormError as ex:
            raise RecipeFormError(_first_msg(ex)) from ex
        src = "file"
    try:
        if src in NEW_BASES:
            return Base(source=src, ref=base_ref(f, src))
        if src == "file":
            path = (f.get("file_path") or "").strip()
            if not path:
                raise ValueError(L("構造ファイルを指定してください", "choose a structure file"))
            if path == (f.get("file_sha_path") or "") and f.get("file_sha256"):
                return Base(source="file", ref=path, sha256=f["file_sha256"])
            p = Path(path).expanduser()
            if not p.is_file():
                raise ValueError(L(f"構造ファイル {p} がありません", f"structure file {p} was not found"))
            st = p.stat(); key = (str(p.resolve()), st.st_mtime_ns, st.st_size)
            if key not in _SHA:
                from adit.builder import file_sha256
                _SHA[key] = file_sha256(p)
            return Base(source="file", ref=path, sha256=_SHA[key])
        ref = source_ref(f, src)
        if not ref:
            raise ValueError(L(f"{src}: 何も指定されていません", f"{src}: nothing given"))
        return Base(source=src, ref=ref)
    except (ValueError, FormError) as ex:
        raise RecipeFormError(_first_msg(ex)) from ex


def recipe_from_form(f: dict) -> Recipe:
    return Recipe(base=base_from_form(f), steps=steps_from_form(f))


def recipe_mode(f: dict) -> bool:
    return step_count(f) > 0 or (f.get("source") or "") in NEW_BASES


def needs_build(f: dict) -> bool:
    return step_count(f) > 0 or (f.get("source") or "") in SLOW_BASES


def form_from_recipe(rec: Recipe) -> dict[str, str]:
    f: dict[str, str] = {"source": rec.base.source, "base_extra": "{}", "file_sha_path": "", "file_sha256": ""}
    if rec.base.source in NEW_BASES and isinstance(rec.base.ref, dict):
        try:
            f.update(base_fields(rec.base.source, rec.base.ref))
        except ValueError:
            pass
    if rec.base.source == "file" and rec.base.sha256:
        f.update(file_sha_path=str(rec.base.ref), file_sha256=rec.base.sha256)
    write_steps(f, list(rec.steps))
    f["recipe_q"] = str(rec.total_charge())
    return f


@dataclass
class Built:

    ref: str
    structure: object
    lines: list[str] = field(default_factory=list)


def structure_for(f: dict, built: Built | None, charge: int, multiplicity: int):
    from adit.builder import recipe_structure

    rec = recipe_from_form(f)
    ref = rec.to_ref()
    if built is not None and built.ref == ref:
        return built.structure, None, rec
    if needs_build(f):
        msg = (L("組み立て手順を変えました。「作る」を押すと作り直します", "the steps have changed; press Build to rebuild") if built is not None
               else L("「作る」を押すと、組み立て手順から構造を作ります", "press Build to make the structure from the steps"))
        raise RecipeFormError(msg)
    st, logs = recipe_structure(rec, charge=charge, multiplicity=multiplicity)
    return st, Built(ref, st, [x.line() for x in logs]), rec


_BEFORE: dict[str, object] = {}


def atoms_before(rec: Recipe, index: int):
    steps = rec.steps[:index]
    if rec.base.source in ("mixture", "polymer") or any(s.op in ("solvent_layer", "solvate") for s in steps):
        return None
    part = Recipe(base=rec.base, steps=steps)
    key = part.to_ref()
    if key not in _BEFORE:
        from adit.builder import build_recipe
        try:
            atoms = build_recipe(part)[0]
        except Exception:
            atoms = None
        if len(_BEFORE) > 32:
            _BEFORE.clear()
        _BEFORE[key] = atoms
    return _BEFORE[key]


def terminations(rec: Recipe, index: int) -> list[str] | None:
    from adit.builder import list_terminations
    s = rec.steps[index]
    if tuple(s.miller) == (0, 0, 0):
        return None
    key = f"terms:{Recipe(base=rec.base, steps=rec.steps[:index]).to_ref()}:{tuple(s.miller)}"
    if key not in _BEFORE:
        atoms = atoms_before(rec, index)
        try:
            _BEFORE[key] = None if atoms is None else list_terminations(atoms, s.miller)
        except Exception:
            _BEFORE[key] = None
    return _BEFORE[key]


def term_options(names: list[str] | None, cur: int) -> list[tuple[str, str]]:
    if names:
        opts = [(str(i), L(f"{i}: 最上面が {pretty_formula(n)}", f"{i}: {pretty_formula(n)} on top")) for i, n in enumerate(names)]
    else:
        opts = [(str(i), L(f"{i} (組成は「作る」のあとに分かります)", f"{i} (composition known after Build)") if i == 0 else str(i)) for i in range(max(6, cur + 1))]
    if str(cur) not in [v for v, _ in opts]:
        opts.append((str(cur), str(cur)))
    return opts


def area_before(rec: Recipe, index: int, built: Built | None) -> float | None:
    atoms = atoms_before(rec, index)
    if atoms is None and built is not None:
        atoms = built.structure.atoms.to_ase()
    if atoms is None or not all(atoms.pbc):
        return None
    return float(np.linalg.norm(np.cross(atoms.cell[0], atoms.cell[1])))


def base_is_slab(f: dict) -> bool:
    src = (f.get("source") or "").strip()
    if src in ("surface", "2d"):
        return True
    if src != "file":
        return False
    from adit.builder import slab_cell_problem
    from adit.structure import from_file
    try:
        atoms = from_file((f.get("file_path") or "").strip())
    except (StructureError, OSError, ValueError):
        return False
    if slab_cell_problem(atoms) is not None:
        return False
    z = np.sort(atoms.get_scaled_positions(wrap=True)[:, 2])
    gaps = np.diff(np.concatenate([z, [z[0] + 1.0]]))
    return float(gaps.max() * atoms.cell.lengths()[2]) > 4.0


def fill_count(step: SolventLayer, row: int, conc: float, area: float | None) -> tuple[SolventLayer, str, bool]:
    from adit.mixture import AMU_TO_G, _one_molecule, salt_count

    comps = list(step.components)
    if not 0 <= row < len(comps):
        return step, L("先に表で、個数を決める行 (イオンなど) を選んでください", "select the table row (such as an ion) whose count you want first"), False
    counts = [c.count for c in comps]

    def volume() -> float | None:
        if step.thickness is None:
            mass = sum(float(_one_molecule(c).get_masses().sum()) * n for c, n in zip(comps, counts)) * AMU_TO_G
            return mass / step.density_g_cm3 * 1e24
        return None if not area else area * step.thickness

    try:
        n = counts[row]
        v = None
        for _ in range(30):
            v = volume()
            if v is None:
                return step, L("断面積がまだ分かりません。先に「作る」を押すか、厚みを密度から決めてください",
                               "the cross-section is not known yet; press Build first, or choose the thickness from a density"), False
            new, _c = salt_count(conc, v)
            if new == n:
                break
            n = counts[row] = new
        if n < 1:
            return step, L(f"帯の体積 {v:.0f} Å³ では {conc:g} mol/L は 1 個未満です。成分を増やすか濃度を上げてください",
                           f"{conc:g} mol/L is less than one ion in a band of {v:.0f} Å³; add more components or raise the concentration"), False
        actual = salt_count(conc, v)[1]
    except Exception as ex:
        return step, str(ex), False
    comps[row] = comps[row].model_copy(update={"count": n})
    name = pretty_name(comps[row].name())
    msg = L(f"帯の体積 {v:.0f} Å³ に {name} を {n} 個 (実際の濃度 {actual:.2f} mol/L)", f"{n} × {name} in a band of {v:.0f} Å³ (actual {actual:.2f} mol/L)")
    if step.thickness is None:
        msg += L("。ほかの行の個数を変えたら、もう一度押してください", ". If you change another row's count, press again")
    return step.model_copy(update={"components": comps}), msg, True
