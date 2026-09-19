
from __future__ import annotations

from dataclasses import dataclass

PLACE_LABELS: dict[str, tuple[str, str]] = {
    "structure": ("構造", "Structure"),
    "structure.atoms": ("構造", "Structure"),
    "structure.charge": ("全電荷", "Total charge"),
    "structure.multiplicity": ("スピン多重度", "Spin multiplicity"),
    "structure.fixed_atoms": ("固定原子", "Fixed atoms"),
    "structure.fixed_axes": ("固定原子", "Fixed atoms"),
    "structure.velocities": ("原子の速度 (前の MD から)", "Atom velocities (from the previous MD)"),
    "handoff": ("前の計算からの引き継ぎ", "Carried over from the previous run"),
    "method.code": ("計算コード", "Code"),
    "method.sk_set": ("Slater-Koster パラメータ", "Slater-Koster parameters"),
    "method.pseudo_set": ("擬ポテンシャルのセット", "Pseudopotential set"),
    "method.pseudo": ("元素ごとの UPF ファイル", "UPF file per element"),
    "method.potcar": ("POTCAR", "POTCAR"),
    "method.third_order": ("DFTB3", "DFTB3"),
    "kpoints": ("k 点", "k-points"),
    "task": ("計算の種類", "Calculation type"),
    "task.md": ("分子動力学 (MD)", "Molecular dynamics"),
    "task.bands": ("バンド計算", "Band structure"),
    "runtime": ("実行環境とリソース", "Runtime and resources"),
    "output_dir": ("出力ディレクトリ", "Output directory"),
    "method.xc": ("汎関数", "Functional"),
    "method.basis_file": ("基底関数のファイル", "Basis-set file"),
    "method.potential_file": ("擬ポテンシャルのファイル", "Pseudopotential file"),
    "method.basis": ("元素ごとの基底と擬ポテンシャル", "Basis and pseudopotential per element"),
    "method.potential": ("元素ごとの基底と擬ポテンシャル", "Basis and pseudopotential per element"),
    "method.dispersion": ("分散補正", "Dispersion correction"),
    "method.cutoff_ry": ("カットオフ [Ry]", "Cutoff [Ry]"),
    "method.rel_cutoff_ry": ("相対カットオフ [Ry]", "Relative cutoff [Ry]"),
    "method.eps_scf": ("SCF の収束の閾値 (EPS_SCF)", "SCF threshold (EPS_SCF)"),
    "method.max_scf": ("SCF の反復の上限 (MAX_SCF)", "Max SCF iterations (MAX_SCF)"),
    "method.uks": ("スピン分極", "Spin polarization"),
    "method.poisson_solver": ("ポアソン方程式の解き方", "Poisson solver"),
    "method.isolated_box_ang": ("分子の箱の一辺 [Å]", "Box edge for a molecule [Å]"),
    "method.extra_sections": ("追加の行 (節ごと)", "Extra lines (per section)"),
    # LAMMPS
    "method.units": ("単位系 (units)", "Units"),
    "method.atom_style": ("原子の形式 (atom_style)", "Atom style"),
    "method.data_file": ("data ファイル", "Data file"),
    "method.type_elements": ("型番号の元素", "Elements of atom types"),
    "method.pair_style": ("pair_style", "pair_style"),
    "method.pair_coeff": ("pair_coeff", "pair_coeff"),
    "method.potential_files": ("写すファイル", "Files to copy"),
    "method.style_commands": ("read_data の前の行", "Commands before read_data"),
    "method.extra_commands": ("pair_coeff の後の行", "Commands after pair_coeff"),
    "method.seed": ("乱数の種", "Random seed"),
    # GROMACS
    "method.topology_file": ("トポロジー (.top)", "Topology (.top)"),
    "method.structure_file": ("構造のファイル (.gro / .pdb)", "Structure file (.gro / .pdb)"),
    "method.coulombtype": ("静電相互作用 (coulombtype)", "Electrostatics (coulombtype)"),
    "method.rcoulomb_nm": ("カットオフ [nm]", "Cut-offs [nm]"),
    "method.rvdw_nm": ("カットオフ [nm]", "Cut-offs [nm]"),
    "method.constraints": ("拘束 (constraints)", "Constraints"),
    "method.pcoupl": ("圧力浴 (pcoupl)", "Barostat (pcoupl)"),
    "method.compressibility_per_bar": ("等温圧縮率 [1/bar]", "Compressibility [1/bar]"),
    "method.define": ("define", "define"),
    "method.checkpoint_file": ("前の段階の .cpt", "Previous stage .cpt"),
    "method.gen_seed": ("初速の乱数の種 (gen-seed)", "Velocity seed (gen-seed)"),
    "method.extra_mdp": ("追加の mdp", "Extra mdp options"),
    "method.ts_search": ("遷移状態の探索 (OptTS)", "Transition-state search (OptTS)"),
    "method.ts_calc_hess": ("最初にヘシアンを計算 (Calc_Hess)", "Compute the Hessian first (Calc_Hess)"),
    "method.ts_recalc_hess": ("ヘシアンを計算し直す間隔 (Recalc_Hess)", "Hessian recalculation interval (Recalc_Hess)"),
    "method.ts_freq": ("最後に振動数を計算 (Freq)", "Frequencies at the end (Freq)"),
    "method.irc": ("反応座標をたどる (IRC)", "Follow the reaction path (IRC)"),
    "method.irc_max_iter": ("IRC の反復の上限", "Max IRC iterations"),
    "method.irc_direction": ("IRC の向き", "IRC direction"),
    "method.model_family": ("機械学習ポテンシャルの種類", "Machine-learning potential"),
    "method.model": ("モデル", "Model"),
    "method.device": ("計算に使うデバイス (device)", "Device"),
    "method.dtype": ("数値の精度 (dtype)", "Precision (dtype)"),
    "compare": ("比べる計算の組", "Set of runs to compare"),
    "conformers": ("配座の候補", "Conformer candidates"),
    "conformers.rmsd": ("重複とみなす RMSD [Å]", "RMSD for duplicates [Å]"),
    "conformers.count": ("作る配座の数", "Number of conformers to embed"),
    "neb": ("反応経路 (NEB)", "Reaction path (NEB)"),
    "neb.images": ("中間の像の数", "Number of intermediate images"),
    "neb.end": ("終状態の構造", "Final-state structure"),
    "phonons": ("フォノン (有限変位)", "Phonons (finite displacements)"),
    "phonons.supercell": ("超格子の倍率", "Supercell"),
    "phonons.distance": ("変位の大きさ [Å]", "Displacement [Å]"),
    "elastic": ("弾性定数 (歪み)", "Elastic constants (strains)"),
    "elastic.strains": ("歪みの大きさ", "Strain magnitudes"),
}


def place_key(location: str) -> str:
    parts = location.split(".")
    for n in range(len(parts), 0, -1):
        key = ".".join(parts[:n])
        if key in PLACE_LABELS:
            ja, _ = PLACE_LABELS[key]
            if n < len(parts) and key not in ("structure", "task", "runtime") and parts[0] == "method":
                break
            return ja
    return parts[-1]


def place_label(location: str) -> str:
    from adit.lang import L

    parts = location.split(".")
    for n in range(len(parts), 0, -1):
        key = ".".join(parts[:n])
        if key in PLACE_LABELS:
            ja, en = PLACE_LABELS[key]
            if n < len(parts) and key not in ("structure", "task", "runtime") and parts[0] == "method":
                break
            return L(ja, en)
    return parts[-1]


def friendly_pydantic(ex) -> str:
    from adit.lang import L

    title = getattr(ex, "title", "") or ""
    prefix = next((p for suffix, p in (("Method", "method"), ("Task", "task"), ("MDSettings", "task.md"), ("Runtime", "runtime"),
                                        ("KPoints", "kpoints"), ("Structure", "structure")) if title.endswith(suffix)), "")
    lines, seen = [], set()
    for e in ex.errors():
        parts = [str(x) for x in e.get("loc", ()) if not str(x).startswith("function-")]
        loc = ".".join([prefix] * bool(prefix) + parts)
        if loc in seen:
            continue
        seen.add(loc)
        where = place_label(loc) if loc else L("入力", "input")
        lines.append(f"{where}: {_why(e)}")
    return "\n".join(lines) or str(ex)


def _why(err: dict) -> str:
    from adit.lang import L

    kind, given = err.get("type", ""), err.get("input")
    if kind == "literal_error":
        allowed = str(err.get("ctx", {}).get("expected", "")).replace("'", "")
        return L(f"{given!r} は使えません。使えるのは {allowed} です",
                 f"{given!r} is not allowed; allowed values are {allowed}")
    if kind == "missing":
        return L("この欄が必要です", "this field is required")
    if kind.endswith("_type") or kind.endswith("_parsing"):
        return L(f"値の型が違います ({err.get('msg', '')})。いまの値: {given!r}",
                 f"wrong type of value ({err.get('msg', '')}); the value was {given!r}")
    return f"{err.get('msg', '')}" + L(f" (いまの値: {given!r})", f" (the value was {given!r})")


@dataclass(frozen=True)
class ValidationError:
    location: str
    message: str

    def __str__(self) -> str:
        return f"{place_label(self.location)}: {self.message}"
