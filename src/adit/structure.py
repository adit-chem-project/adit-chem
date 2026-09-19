"""Obtain structures as ASE Atoms: presets, SMILES, files, bulk crystals, surfaces and mixtures."""

from __future__ import annotations

from adit.errors import AditError
import re

from pathlib import Path

from ase import Atoms

from adit.lang import L
from ase.build import molecule
from ase import build as ase_build
from ase.build import bulk
from ase.collections import g2
from ase.data import atomic_numbers, reference_states
from ase.io import read

from adit.spec import AtomsData, Structure


class StructureError(AditError):
    pass


def preset_names() -> list[str]:
    return sorted(g2.names)


_PRESET_ALIASES = {
    "H2O": ("water", "水"),
    "NH3": ("ammonia", "アンモニア"),
    "CH4": ("methane", "メタン"),
    "CH3OH": ("methanol", "メタノール"),
    "CH3CH2OH": ("ethanol", "エタノール"),
    "CH3OCH3": ("dimethyl ether", "ジメチルエーテル"),
    "CH3CHO": ("acetaldehyde", "アセトアルデヒド"),
    "CH3COCH3": ("acetone", "アセトン"),
    "CH3COOH": ("acetic acid", "酢酸"),
    "HCOOH": ("formic acid", "ギ酸"),
    "C6H6": ("benzene", "ベンゼン"),
}


def preset_search_text(name: str) -> str:
    return " ".join((name, pretty_formula(name), *_PRESET_ALIASES.get(name, ())))


def preset_matches(name: str, query: str) -> bool:
    q = query.strip().casefold()
    if not q:
        return True
    if q in name.casefold() or q in pretty_formula(name).casefold():
        return True
    return any(alias.casefold().startswith(q) for alias in _PRESET_ALIASES.get(name, ()))


_SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def pretty_formula(name: str) -> str:
    if not re.fullmatch(r"[A-Z][A-Za-z0-9]*", name):
        return name
    return re.sub(r"(?<=[A-Za-z])(\d+)", lambda m: m.group(1).translate(_SUB), name)


def from_preset(name: str) -> Atoms:
    if name not in g2.names:
        import difflib

        from adit.lang import L
        near = difflib.get_close_matches(name, list(g2.names), n=5, cutoff=0.5)
        near_s = L(f" 近い名前: {', '.join(near)}。", f" Close names: {', '.join(near)}.") if near else ""
        raise StructureError(L(f"プリセットに {name!r} はありません (プリセットは ASE の g2 セットの名前で、H2O、CH3CH2OH のような分子式が中心)。{near_s}"
                               f"一覧に無い分子は、種類を SMILES にして書くと作れます (例 エタノール = CCO)",
                               f"no preset named {name!r} (presets are the names of ASE's g2 set, mostly formulas such as H2O, CH3CH2OH).{near_s} "
                               f"For other molecules, switch the kind to SMILES (e.g. ethanol = CCO)"))
    return molecule(name)


# ---- 2. SMILES ----
def has_rdkit() -> bool:
    try:
        import rdkit  # noqa: F401
    except ImportError:
        return False
    return True


def from_smiles(smiles: str, *, seed: int = 0) -> Atoms:
    if not has_rdkit():
        raise StructureError(L("RDKit が無いので SMILES からは作れません", "RDKit is not installed, so SMILES cannot be used"))
    from rdkit import Chem
    from rdkit.Chem import AllChem

    if not smiles.strip():
        raise StructureError(L("SMILES が空です", "the SMILES string is empty"))
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise StructureError(L(f"SMILES を解釈できません: {smiles!r}", f"cannot parse SMILES: {smiles!r}"))
    if mol.GetNumAtoms() == 0:
        raise StructureError(L(f"SMILES に原子がありません: {smiles!r}", f"the SMILES has no atoms: {smiles!r}"))
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise StructureError(L(f"3 次元座標を作れませんでした: {smiles!r}", f"could not build 3D coordinates: {smiles!r}"))
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    positions = [tuple(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())]
    return Atoms(symbols=symbols, positions=positions)


def from_file(path: Path | str | None) -> Atoms:
    if path is None or not str(path).strip():
        raise StructureError(L("構造ファイルを指定してください", "choose a structure file"))
    path = Path(str(path).strip()).expanduser()
    if not path.is_file():
        what = L("ディレクトリで、ファイルではありません", "is a directory, not a file") if path.is_dir() else L("がありません", "was not found")
        raise StructureError(L(f"構造ファイル {path} {what}", f"structure file {path} {what}"))
    try:
        result = read(path)
    except Exception as ex:
        raise StructureError(L(f"構造ファイル {path} を読めません (形式が拡張子と合っているか確かめてください): {ex}",
                               f"cannot read the structure file {path} (check that the format matches the extension): {ex}")) from ex
    atoms = result[-1] if isinstance(result, list) else result
    if len(atoms) == 0:
        raise StructureError(L(f"構造ファイル {path} に原子がありません", f"the structure file {path} has no atoms"))
    return atoms


CRYSTAL_STRUCTURES = ["sc", "fcc", "bcc", "hcp", "diamond", "zincblende", "rocksalt", "cesiumchloride", "fluorite", "wurtzite"]


def default_bulk(symbol: str) -> tuple[str | None, float | None]:
    z = atomic_numbers.get(symbol)
    ref = reference_states[z] if z is not None and z < len(reference_states) else None
    if not ref:
        return None, None
    return ref.get("symmetry"), ref.get("a")


def from_bulk(ref: str) -> Atoms:
    tokens = ref.split()
    if not tokens:
        raise StructureError(L("バルク: 元素を指定してください (例 Si、NaCl rocksalt 5.64)", "bulk: give an element (e.g. Si, NaCl rocksalt 5.64)"))
    symbol, rest = tokens[0], tokens[1:]
    cubic = "cubic" in rest
    rest = [t for t in rest if t != "cubic"]
    structure = rest[0] if rest and rest[0] in CRYSTAL_STRUCTURES else None
    if structure:
        rest = rest[1:]
    a = None
    if rest:
        try:
            a = float(rest[0])
        except ValueError as ex:
            raise StructureError(L(f"バルク: 格子定数を数値として読めません: {rest[0]!r}", f"bulk: the lattice constant is not a number: {rest[0]!r}")) from ex
    try:
        return bulk(symbol, structure, a=a, cubic=cubic)
    except Exception as ex:
        raise StructureError(L(f"バルク: {ref!r} を作れません: {ex}", f"bulk: cannot build {ref!r}: {ex}")) from ex


def from_spacegroup(ref: str) -> Atoms:
    from ase.spacegroup import crystal

    tokens = ref.split()
    if len(tokens) < 3:
        raise StructureError(L("空間群: 「空間群 元素:x,y,z … cell=a[,b,c[,α,β,γ]]」の形で書いてください "
                               "(例: 225 Na:0,0,0 Cl:0.5,0,0 cell=5.64)",
                               "space group: write it as 'group element:x,y,z ... cell=a[,b,c[,alpha,beta,gamma]]' "
                               "(e.g. 225 Na:0,0,0 Cl:0.5,0,0 cell=5.64)"))
    group_text, rest = tokens[0], tokens[1:]
    cell_tokens = [t for t in rest if t.lower().startswith("cell=")]
    if not cell_tokens:
        raise StructureError(L("空間群: 格子定数を cell=a[,b,c[,α,β,γ]] の形で書いてください",
                               "space group: give the lattice with cell=a[,b,c[,alpha,beta,gamma]]"))
    try:
        cellpar = [float(x) for x in cell_tokens[-1].split("=", 1)[1].split(",") if x.strip()]
    except ValueError as ex:
        raise StructureError(L(f"空間群: 格子定数を数として読めません: {cell_tokens[-1]!r}",
                               f"space group: the lattice parameters are not numbers: {cell_tokens[-1]!r}")) from ex
    if len(cellpar) == 1:
        cellpar = cellpar * 3 + [90.0, 90.0, 90.0]
    elif len(cellpar) == 3:
        cellpar = cellpar + [90.0, 90.0, 90.0]
    elif len(cellpar) != 6:
        raise StructureError(L("空間群: 格子定数は 1 つ (立方)、3 つ (a,b,c)、6 つ (a,b,c,α,β,γ) のどれかで書いてください",
                               "space group: give 1 (cubic), 3 (a,b,c) or 6 (a,b,c,alpha,beta,gamma) lattice parameters"))
    symbols, basis = [], []
    for token in rest:
        if token.lower().startswith("cell="):
            continue
        if ":" not in token:
            raise StructureError(L(f"空間群: 原子は「元素:x,y,z」の形で書いてください: {token!r}",
                                   f"space group: write each atom as element:x,y,z: {token!r}"))
        element, _, coords = token.partition(":")
        try:
            xyz = [float(x) for x in coords.split(",")]
        except ValueError as ex:
            raise StructureError(L(f"空間群: 分率座標を数として読めません: {token!r}",
                                   f"space group: the fractional coordinates are not numbers: {token!r}")) from ex
        if len(xyz) != 3:
            raise StructureError(L(f"空間群: 分率座標は 3 つ書いてください: {token!r}",
                                   f"space group: give three fractional coordinates: {token!r}"))
        symbols.append(element)
        basis.append(xyz)
    if not symbols:
        raise StructureError(L("空間群: 原子を 1 つ以上書いてください", "space group: give at least one atom"))
    group: int | str = int(group_text) if group_text.isdigit() else group_text
    try:
        return crystal(symbols, basis=basis, spacegroup=group, cellpar=cellpar)
    except Exception as ex:
        raise StructureError(L(f"空間群: この組み合わせでは作れません ({ex})",
                               f"space group: cannot build this combination ({ex})")) from ex


SURFACE_FUNCTIONS = ["fcc100", "fcc110", "fcc111", "bcc100", "bcc110", "bcc111", "hcp0001", "hcp10m10", "diamond100", "diamond111"]


def from_surface(ref: str) -> Atoms:
    tokens = ref.split()
    if len(tokens) < 3:
        raise StructureError(L("スラブ: 「面 元素 nx×ny×層数」の 3 つが要ります (例: fcc111 Al 2x2x3)",
                               "slab: give facet, element and nx x ny x layers (e.g. fcc111 Al 2x2x3)"))
    facet, symbol, size_s = tokens[0], tokens[1], tokens[2]
    if facet not in SURFACE_FUNCTIONS:
        raise StructureError(L(f"スラブ: 面 {facet!r} には対応していません ({', '.join(SURFACE_FUNCTIONS)})",
                               f"slab: facet {facet!r} is not supported ({', '.join(SURFACE_FUNCTIONS)})"))
    try:
        size = tuple(int(x) for x in size_s.lower().split("x"))
        assert len(size) == 3 and all(n >= 1 for n in size)
    except (ValueError, AssertionError) as ex:
        raise StructureError(L(f"スラブ: 大きさは nx×ny×層数 の形で書いてください (例 2x2x3): {size_s!r}",
                               f"slab: write the size as nx x ny x layers (e.g. 2x2x3): {size_s!r}")) from ex
    kw: dict = {"vacuum": 10.0}
    for t in tokens[3:]:
        k, _, v = t.partition("=")
        if k not in ("vacuum", "a") or not v:
            raise StructureError(L(f"スラブ: 解釈できない指定: {t!r} (vacuum=<Å> か a=<Å>)", f"slab: cannot interpret {t!r} (use vacuum=<Å> or a=<Å>)"))
        try:
            kw[k] = float(v)
        except ValueError as ex:
            raise StructureError(L(f"スラブ: {k} を数値として読めません: {v!r}", f"slab: {k} is not a number: {v!r}")) from ex
    try:
        atoms = getattr(ase_build, facet)(symbol, size=size, **kw)
    except Exception as ex:
        raise StructureError(L(f"スラブ: {ref!r} を作れません: {ex}", f"slab: cannot build {ref!r}: {ex}")) from ex
    atoms.pbc = (True, True, True)
    return atoms


def build_structure(source: str, source_ref: str, *, charge: int = 0, multiplicity: int = 1) -> Structure:
    if source == "preset":
        atoms = from_preset(source_ref)
    elif source == "smiles":
        atoms = from_smiles(source_ref)
    elif source == "file":
        atoms = from_file(source_ref)
    elif source == "bulk":
        atoms = from_bulk(source_ref)
    elif source == "surface":
        atoms = from_surface(source_ref)
    elif source == "mixture":
        from adit.mixture import MixtureError, MixtureSpec, build_mixture

        try:
            atoms = build_mixture(MixtureSpec.from_ref(source_ref))
        except MixtureError as ex:
            raise StructureError(str(ex)) from ex
    elif source == "recipe":
        from adit.builder import recipe_structure

        st, _logs = recipe_structure(source_ref, charge=charge, multiplicity=multiplicity)
        return st
    else:
        raise StructureError(L(f"構造の作り方 {source!r} には対応していません", f"unknown structure source: {source!r}"))
    return Structure(source=source, source_ref=source_ref, atoms=AtomsData.from_ase(atoms),
                     charge=charge, multiplicity=multiplicity)
