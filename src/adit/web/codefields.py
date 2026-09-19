
from __future__ import annotations

import re
from adit.lang import L

_SECTION = re.compile(r"^\[\s*([A-Za-z0-9_/]+)\s*\]$")


def parse_sections(text: str) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    cur: str | None = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        m = _SECTION.match(line)
        if m:
            cur = m.group(1).upper()
            out.setdefault(cur, [])
            continue
        if cur is None:
            raise ValueError(L(f"追加の行: 先に [セクションのパス] の見出しを書いてください (例 [FORCE_EVAL/DFT/SCF]): {raw.strip()!r}",
                               f"extra lines: write a [section path] header first (e.g. [FORCE_EVAL/DFT/SCF]): {raw.strip()!r}"))
        out[cur].append(line)
    return {k: "\n".join(v) for k, v in out.items() if v}


def sections_text(d: dict[str, str]) -> str:
    return "\n".join(f"[{k}]\n{v.strip()}" for k, v in d.items())


def parse_lines(text: str) -> list[str]:
    return [l.strip() for l in text.splitlines() if l.strip()]


def parse_words(text: str) -> list[str]:
    return text.replace(",", " ").split()


def _mdp_value(v: str):
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v


def parse_mdp(text: str) -> dict:
    out: dict = {}
    for raw in text.splitlines():
        line = re.split(r"[;#]", raw, maxsplit=1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise ValueError(L(f"追加の mdp: 「名前 = 値」の形ではありません: {raw.strip()!r}", f"extra mdp: not of the form name = value: {raw.strip()!r}"))
        k, _, v = line.partition("=")
        if not k.strip():
            raise ValueError(L(f"追加の mdp: 名前が空です: {raw.strip()!r}", f"extra mdp: empty name: {raw.strip()!r}"))
        out[k.strip()] = _mdp_value(v.strip())
    return out


def mdp_text(d: dict) -> str:
    return "\n".join(f"{k} = {('yes' if v else 'no') if isinstance(v, bool) else v}" for k, v in d.items())


def cp2k_candidates(cp2k_data: str, basis_file: str, potential_file: str, elements: list[str]) -> tuple[str, list[tuple[str, list[str], list[str]]]]:
    from adit.codes.cp2k_data import Cp2kData

    data = Cp2kData.discover(cp2k_data)
    if not data.found:
        return "", []
    return str(data.root), [(e, data.names_for(basis_file, e, potential=False), data.names_for(potential_file, e, potential=True)) for e in elements]


def cp2k_files(cp2k_data: str) -> tuple[list[str], list[str]]:
    from adit.codes.cp2k_data import Cp2kData

    data = Cp2kData.discover(cp2k_data)
    names = data.files()
    return [n for n in names if "BASIS" in n.upper()], [n for n in names if "POTENTIAL" in n.upper()]



def _num(text: str, what: str) -> float:
    import math

    try:
        x = float(text.strip())
    except ValueError as ex:
        raise ValueError(L(f"{what}: 数値として読めません: {text.strip()!r}", f"{what}: not a number: {text.strip()!r}")) from ex
    if not math.isfinite(x):
        raise ValueError(L(f"{what}: 有限の数を入れてください: {text.strip()!r}", f"{what}: enter a finite number: {text.strip()!r}"))
    return x


def hubbard_from_rows(rows) -> dict:
    from adit.spec import HubbardU

    out: dict = {}
    for el, orb, u, j in rows:
        el, orb, u, j = (str(x or "").strip() for x in (el, orb, u, j))
        if not el:
            if orb or u or j:
                raise ValueError(L("DFT+U: 元素を選んでいない行があります", "DFT+U: a row has no element"))
            continue
        if el in out:
            raise ValueError(L(f"DFT+U: {el} の行が 2 つあります", f"DFT+U: two rows for {el}"))
        if not orb:
            raise ValueError(L(f"DFT+U: {el} の殻を入れてください (例 3d、4f)", f"DFT+U: enter the shell for {el} (e.g. 3d, 4f)"))
        if not u:
            raise ValueError(L(f"DFT+U: {el} の U [eV] を入れてください (既定値はありません)", f"DFT+U: enter U [eV] for {el} (there is no default)"))
        out[el] = HubbardU(orbital=orb, u_ev=_num(u, f"DFT+U {el} U"), j_ev=_num(j, f"DFT+U {el} J") if j else 0.0)
    return out


def hubbard_rows(hubbard: dict) -> list[tuple[str, str, str, str]]:
    return [(e, h.orbital, f"{h.u_ev:g}", f"{h.j_ev:g}" if h.j_ev else "") for e, h in hubbard.items()]


def element_values_from_rows(rows, what: str) -> dict[str, float]:
    out: dict[str, float] = {}
    for el, text in rows:
        text = str(text or "").strip()
        if not text:
            continue
        v = _num(text, f"{what} ({el})")
        if v != 0:
            out[el] = v
    return out


STAGE_TYPES = ("geometry_optimization", "molecular_dynamics", "single_point", "vibrations", "band_structure")
STAGE_COLUMNS = ("name", "type", "ensemble", "thermostat", "temperature", "steps", "timestep", "pressure", "velocities", "extra")


def stage_column_titles() -> list[str]:
    return [L("名前", "Name"), L("計算の種類", "Task"), L("アンサンブル", "Ensemble"), L("熱浴", "Thermostat"), L("温度 [K]", "Temperature [K]"),
            L("ステップ数", "Steps"), L("時間刻み [fs]", "Time step [fs]"), L("圧力 [bar]", "Pressure [bar]"),
            L("速度", "Velocities"), L("ほかの項目 (JSON)", "Other items (JSON)")]


SHELLS = ["3d", "4d", "5d", "4f", "5f", "2p", "3p", "4p"]


def stage_choices(col: str) -> list[tuple[str, str]] | None:
    keep = ("", L("(画面のまま)", "(unchanged)"))
    if col == "type":
        names = {"geometry_optimization": L("構造最適化", "Geometry optimization"), "molecular_dynamics": L("分子動力学", "Molecular dynamics"),
                 "single_point": L("一点計算", "Single point"), "vibrations": L("振動解析", "Vibrations"), "band_structure": L("バンド計算", "Band structure")}
        return [keep] + [(k, names[k]) for k in STAGE_TYPES]
    if col == "ensemble":
        return [keep] + [(x, x) for x in ("NVT", "NVE", "NPT")]
    if col == "thermostat":
        return [keep] + [("berendsen", "Berendsen"), ("andersen", "Andersen"), ("nose_hoover", "Nosé-Hoover"), ("langevin", "Langevin"),
                         ("csvr", L("CSVR (速度再スケール)", "CSVR (velocity rescaling)"))]
    if col == "velocities":
        return stage_velocity_choices()
    return None


def carry_table() -> str:
    return L("MD の続き (前も今回も MD、同じ計算コード) で引き継ぐもの:\n"
             "  DFTB+・VASP  位置と速度\n"
             "  xtb  位置と mdrestart (座標と速度)\n"
             "  CP2K  位置・セルと restart (速度)。熱浴・圧力浴の状態は引き継ぎません\n"
             "  LAMMPS  final.data (位置と速度)\n"
             "  GROMACS  adit.gro と adit.cpt (速度と熱浴・圧力浴の状態)\n"
             "  Quantum ESPRESSO・ORCA  位置 (とセル) だけ。速度は出力に書かれないので引き継げません\n"
             "構造最適化の続きは最終構造から始めます。一点計算・振動解析・バンド計算は構造を動かさないので、前の入力の構造のままです。",
             "What an MD continuation (MD before and now, same code) carries over:\n"
             "  DFTB+, VASP  positions and velocities\n"
             "  xtb  positions and mdrestart (coordinates and velocities)\n"
             "  CP2K  positions, cell and the restart file (velocities); thermostat and barostat state are not carried over\n"
             "  LAMMPS  final.data (positions and velocities)\n"
             "  GROMACS  adit.gro and adit.cpt (velocities and thermostat/barostat state)\n"
             "  Quantum ESPRESSO, ORCA  positions (and cell) only; velocities are not written to their output\n"
             "An optimization continues from its final structure. Single points, vibrations and band structures do not move atoms, "
             "so the previous input structure is used.")


def stage_velocity_choices() -> list[tuple[str, str]]:
    return [("auto", L("自動 (MD → MD なら引き継ぐ)", "auto (carried over from MD to MD)")), ("yes", L("引き継ぐ", "carry over")),
            ("no", L("引き継がない", "do not carry over"))]


def stage_from_row(row: dict, index: int, base_type: str) -> dict | None:
    import json

    r = {k: str(row.get(k) or "").strip() for k in STAGE_COLUMNS}
    vel = r["velocities"] or "auto"
    if not any(r[k] for k in STAGE_COLUMNS if k != "velocities") and vel == "auto":
        return None
    where = L(f"{index} 段階目", f"stage {index}")
    task: dict = {}
    md: dict = {}
    kind = r["type"] or base_type
    if r["type"]:
        if r["type"] not in STAGE_TYPES:
            raise ValueError(L(f"{where}: 計算の種類 {r['type']!r} は使えません", f"{where}: unknown task {r['type']!r}"))
        task["type"] = r["type"]
    if r["ensemble"]:
        md["ensemble"] = r["ensemble"]
    if r["thermostat"]:
        md["thermostat"] = r["thermostat"]
    if r["temperature"]:
        md["temperature_k"] = _num(r["temperature"], L(f"{where}の温度 [K]", f"{where} temperature [K]"))
    if r["timestep"]:
        md["timestep_fs"] = _num(r["timestep"], L(f"{where}の時間刻み [fs]", f"{where} time step [fs]"))
    if r["pressure"]:
        md["pressure_bar"] = _num(r["pressure"], L(f"{where}の圧力 [bar]", f"{where} pressure [bar]"))
    if r["steps"]:
        n = _num(r["steps"], L(f"{where}のステップ数", f"{where} steps"))
        if n != int(n):
            raise ValueError(L(f"{where}のステップ数は整数にしてください", f"{where}: steps must be an integer"))
        if kind == "geometry_optimization":
            task["max_steps"] = int(n)
        elif kind == "molecular_dynamics":
            md["steps"] = int(n)
        else:
            raise ValueError(L(f"{where}: ステップ数を使うのは構造最適化と分子動力学だけです", f"{where}: steps apply only to optimization and MD"))
    if md:
        task["md"] = md
    out: dict = {"name": r["name"] or f"stage{index}"}
    if task:
        out["task"] = task
    if vel != "auto":
        out["velocities"] = vel == "yes"
    if r["extra"]:
        try:
            extra = json.loads(r["extra"])
        except json.JSONDecodeError as ex:
            raise ValueError(L(f"{where}のほかの項目: JSON として読めません ({ex.msg})", f"{where} other items: not valid JSON ({ex.msg})")) from ex
        if not isinstance(extra, dict):
            raise ValueError(L(f'{where}のほかの項目は {{"method": {{…}}}} の形にしてください', f'{where} other items must look like {{"method": {{...}}}}'))
        for k, v in extra.items():
            if k in ("task",) and isinstance(v, dict) and isinstance(out.get("task"), dict):
                out["task"] = _merge_dicts(v, out["task"])
            elif k not in out:
                out[k] = v
    return out


def _merge_dicts(a: dict, b: dict) -> dict:
    out = dict(a)
    for k, v in b.items():
        out[k] = _merge_dicts(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def row_from_stage(stage: dict) -> dict:
    import json

    task = dict(stage.get("task") or {})
    md = dict(task.pop("md", None) or {})
    kind = task.pop("type", "")
    row = {"name": str(stage.get("name", "")), "type": kind, "ensemble": md.pop("ensemble", ""), "thermostat": md.pop("thermostat", ""),
           "temperature": _g(md.pop("temperature_k", "")), "timestep": _g(md.pop("timestep_fs", "")), "pressure": _g(md.pop("pressure_bar", "")),
           "steps": "", "velocities": {True: "yes", False: "no"}.get(stage.get("velocities"), "auto")}
    if kind == "geometry_optimization" and "max_steps" in task:
        row["steps"] = _g(task.pop("max_steps"))
    elif "steps" in md:
        row["steps"] = _g(md.pop("steps"))
    rest = {k: v for k, v in stage.items() if k not in ("name", "task", "velocities", "dir")}
    if md:
        task["md"] = md
    if task:
        rest["task"] = task
    row["extra"] = json.dumps(rest, ensure_ascii=False) if rest else ""
    return row


def _g(v) -> str:
    return f"{v:g}" if isinstance(v, (int, float)) and not isinstance(v, bool) else str(v)


COMMON_LABELS = {
    "method.code": "計算コード", "task.type": "種類", "task.optimizer": "最適化アルゴリズム", "task.max_steps": "最大ステップ数 (MaxSteps)",
    "task.force_tolerance_ev_per_ang": "力の収束判定 [eV/Å]", "task.relax_cell": "セルの緩和 (周期系)", "task.md.ensemble": "アンサンブル",
    "task.md.thermostat": "熱浴", "task.md.temperature_k": "温度 [K]", "task.md.timestep_fs": "時間刻み [fs]", "task.md.steps": "MD ステップ数",
    "task.md.dump_interval": "軌跡の出力間隔 [ステップ]", "task.md.coupling_time_fs": "熱浴の緩和時間 [fs]", "task.md.pressure_bar": "圧力 [bar] (NPT)",
    "task.md.barostat_time_fs": "圧力浴の緩和時間 [fs] (NPT)", "task.bands.path": "k 点の経路 (バンド)", "task.bands.npoints": "経路上の k 点数",
    "task.bands.empty_bands": "空のバンド数 (pw.x)", "kpoints.mode": "サンプリング方法", "kpoints.mesh": "メッシュ (n1 n2 n3)", "kpoints.shift": "メッシュ (n1 n2 n3)",
    "kpoints.density": "k 点密度 [点/Å⁻¹]", "runtime.profile": "プロファイル", "runtime.nodes": "ノード数", "runtime.ncpus": "ノードあたりのコア数",
    "runtime.mpiprocs": "ノードあたりの MPI プロセス数", "runtime.omp_threads": "OpenMP スレッド数", "runtime.walltime": "制限時間 (HH:MM:SS)", "runtime.job_name": "ジョブ名",
}
METHOD_LABELS = {
    "dftbplus": {"sk_set": "Slater-Koster パラメータ", "scc": "SCC (自己無撞着電荷)", "scc_tolerance": "SCC の収束判定 (SccTolerance)", "max_scc_iterations": "SCC の反復上限 (MaxSccIterations)",
                 "third_order": "DFTB3 (ThirdOrderFull)", "dispersion": "分散力補正", "d3_params": "DFT-D3 係数 (BJ)", "filling_temperature": "電子温度 [K]",
                 "solvation_param_file": "溶媒のパラメータファイル (GBSA)"},
    "vasp": {"potcar_set": "POTCAR のセット", "potcar": "元素ごとの POTCAR 名", "binary": "実行ファイルの種類", "encut": "ENCUT [eV] (0 = 指定しない)",
             "ediff": "EDIFF [eV]", "nelm": "NELM", "ismear": "ISMEAR", "sigma": "SIGMA [eV]", "ispin": "ISPIN", "magmom": "MAGMOM", "ivdw": "IVDW",
             "ibrion": "IBRION (構造最適化の方法)", "algo": "ALGO", "prec": "PREC", "lreal": "LREAL", "magmom_by_element": "元素ごとの初期磁気モーメント [μB]",
             "hubbard": "DFT+U", "ldau_type": "LDAUTYPE", "extra_incar": "追加の INCAR 設定"},
    "xtb": {"gfn": "計算手法 (--gfn)", "accuracy": "精度 (--acc)", "etemp": "電子温度 [K] (--etemp)", "max_iterations": "SCC の最大反復回数",
            "opt_level": "最適化の収束レベル (--opt)", "solvation": "溶媒モデル (--alpb / --gbsa)", "solvent": "溶媒"},
    "espresso": {"pseudo_set": "擬ポテンシャルのセット", "pseudo": "元素ごとの UPF ファイル", "ecutwfc": "ecutwfc [Ry]", "ecutrho": "ecutrho [Ry] (0 = 指定しない)",
                 "conv_thr": "conv_thr [Ry]", "electron_maxstep": "electron_maxstep", "mixing_beta": "mixing_beta", "occupations": "occupations",
                 "smearing": "smearing", "degauss": "degauss [Ry]", "nspin": "nspin", "input_dft": "input_dft",
                 "starting_magnetization": "starting_magnetization (元素ごと)", "hubbard": "DFT+U", "hubbard_projector": "HUBBARD の射影",
                 "extra": "追加の変数 (名前空間.変数 = 値)"},
    "orca": {"method": "計算手法 (! 行)", "basis": "基底関数", "scf_convergence": "SCF の収束判定", "scf_maxiter": "SCF の最大反復回数",
             "maxcore_mb": "%maxcore [MB] (0 = 指定しない)", "extra_keywords": "追加のキーワード (! 行)", "extra_blocks": "追加の %ブロック",
             "solvation": "溶媒モデル (CPCM / SMD)", "solvent": "溶媒"},
    "cp2k": {"basis_file": "基底関数のファイル", "potential_file": "擬ポテンシャルのファイル", "basis": "元素ごとの基底と擬ポテンシャル",
             "potential": "元素ごとの基底と擬ポテンシャル", "xc": "汎関数", "dispersion": "分散補正", "cutoff_ry": "カットオフ [Ry]",
             "rel_cutoff_ry": "相対カットオフ [Ry]", "eps_scf": "SCF の収束の閾値 (EPS_SCF)", "max_scf": "SCF の反復の上限 (MAX_SCF)", "uks": "スピン分極",
             "poisson_solver": "ポアソン方程式の解き方", "isolated_box_ang": "分子の箱の一辺 [Å]", "magnetization_by_element": "元素ごとの MAGNETIZATION",
             "hubbard": "DFT+U", "plus_u_method": "PLUS_U_METHOD", "sccs_relative_permittivity": "溶媒の比誘電率 (SCCS)", "extra_sections": "追加の行 (節ごと)"},
    "lammps": {"units": "単位系 (units)", "atom_style": "原子の形式 (atom_style)", "data_file": "data ファイル", "type_elements": "型番号の元素",
               "pair_style": "pair_style", "pair_coeff": "pair_coeff", "potential_files": "写すファイル", "style_commands": "read_data の前の行",
               "extra_commands": "pair_coeff の後の行", "seed": "乱数の種"},
    "gromacs": {"topology_file": "トポロジー (.top)", "structure_file": "構造のファイル (.gro / .pdb)", "coulombtype": "静電相互作用 (coulombtype)",
                "rcoulomb_nm": "カットオフ [nm]", "rvdw_nm": "カットオフ [nm]", "constraints": "拘束 (constraints)", "pcoupl": "圧力浴 (pcoupl)",
                "compressibility_per_bar": "等温圧縮率 [1/bar]", "define": "define", "checkpoint_file": "前の段階の .cpt", "gen_seed": "初速の乱数の種 (gen-seed)",
                "extra_mdp": "追加の mdp"},
}


def label_for_path(code: str, path: str) -> str | None:
    parts = path.split(".")
    if parts[0] == "method" and len(parts) > 1:
        return COMMON_LABELS.get("method.code") if parts[1] == "code" else METHOD_LABELS.get(code, {}).get(parts[1])
    for n in range(len(parts), 0, -1):
        got = COMMON_LABELS.get(".".join(parts[:n]))
        if got:
            return got
    return None


def template_marks(spec) -> dict[str, str]:
    from adit.templates import flatten, template_origin

    t = spec.meta.template
    if not t:
        return {}
    code = spec.method.code
    values = t.get("values", {})
    origin = template_origin(spec)
    ok: dict[str, bool] = {}
    for path in values:
        lab = label_for_path(code, path)
        if lab:
            ok[lab] = ok.get(lab, True) and path in origin
    for path in flatten(spec.model_dump(mode="json")):
        lab = label_for_path(code, path)
        if lab in ok and path not in values:
            ok[lab] = False
    return {lab: t["name"] for lab, good in ok.items() if good}


# Spec fields the web form has no widget for. PrepOrigin keeps them from the loaded
# spec.json so that load -> preview -> generate does not silently reset them.
HIDDEN_METHOD_FIELDS = {
    "xtb": ("md_hmass", "md_shake", "md_sccacc"),
    "espresso": ("dipole_correction", "dipole_direction", "dipole_maxpos", "dipole_decrease", "dipole_amplitude"),
    "lammps": ("thermo_pressure_tensor",),
}
HIDDEN_META_FIELDS = ("comment", "stage")


def hidden_method_values(method) -> dict:
    fields = type(method).model_fields
    return {k: getattr(method, k) for k in HIDDEN_METHOD_FIELDS.get(method.code, ()) if getattr(method, k) != fields[k].default}


def hidden_meta_values(meta) -> dict:
    return {k: getattr(meta, k) for k in HIDDEN_META_FIELDS if getattr(meta, k)}


def uneven_shift(kpoints):
    # The form has one shift field for all three directions.
    if kpoints is None or len(set(kpoints.shift)) == 1:
        return None
    return tuple(kpoints.shift)


def _shift_text(shift) -> str:
    return "(" + ", ".join(f"{x:g}" for x in shift) + ")"


class PrepOrigin:

    def __init__(self, spec=None):
        self.continued_from = self.handoff = self.velocities = self.template = None
        self.plumed = None
        self.code, self.source_ref, self.symbols = "", "", []
        self.hidden_code, self.hidden_method, self.hidden_meta, self.kp_shift = "", {}, {}, None
        if spec is not None:
            self.take(spec)

    def take(self, spec, *, template_only: bool = False) -> None:
        self.template = spec.meta.template
        self.hidden_code, self.hidden_method = spec.method.code, hidden_method_values(spec.method)
        self.kp_shift = uneven_shift(spec.kpoints)
        if template_only:
            return
        self.hidden_meta = hidden_meta_values(spec.meta)
        self.plumed = spec.plumed
        self.continued_from, self.handoff = spec.meta.continued_from, spec.handoff
        self.velocities = spec.structure.velocities
        self.code, self.source_ref, self.symbols = spec.method.code, spec.structure.source_ref, list(spec.structure.atoms.symbols)

    def clear_continuation(self) -> None:
        self.continued_from = self.handoff = self.velocities = self.plumed = None
        self.code, self.source_ref, self.symbols = "", "", []

    @property
    def active(self) -> bool:
        return bool(self.continued_from or self.handoff or self.template or self.hidden_method or self.hidden_meta or self.kp_shift)

    def _shift_kept(self, spec) -> bool:
        kp = spec.kpoints
        if self.kp_shift is None or kp is None or len(set(kp.shift)) != 1:
            return False
        return f"{kp.shift[0]:g}" == f"{self.kp_shift[0]:g}"

    def same_structure(self, st) -> bool:
        return bool(self.source_ref) and st.source_ref == self.source_ref and list(st.atoms.symbols) == self.symbols

    def apply(self, spec):
        from adit.spec import Handoff

        meta_upd: dict = {}
        if self.template:
            meta_upd["template"] = self.template
        upd: dict = {}
        if (self.continued_from or self.handoff) and self.same_structure(spec.structure):
            meta_upd["continued_from"] = self.continued_from
            md = spec.task.type == "molecular_dynamics"
            same_code = spec.method.code == self.code
            h = self.handoff
            if h is not None:
                if not (md and same_code) and (h.velocities or h.files):
                    h = Handoff(previous_dir=h.previous_dir, previous_task=h.previous_task, previous_code=h.previous_code, at_run=h.at_run)
                upd["handoff"] = h
            if md and same_code:
                if self.velocities is not None and spec.structure.velocities is None:
                    upd["structure"] = spec.structure.model_copy(update={"velocities": self.velocities})
            elif spec.structure.velocities is not None:
                upd["structure"] = spec.structure.model_copy(update={"velocities": None})
        if self.plumed is not None and spec.plumed is None:
            from adit.codes.plumed import SUPPORTED_CODES

            if spec.method.code in SUPPORTED_CODES and spec.task.type == "molecular_dynamics":
                upd["plumed"] = self.plumed
        if self.hidden_method and spec.method.code == self.hidden_code:
            upd["method"] = spec.method.model_copy(update=dict(self.hidden_method))
        if self._shift_kept(spec):
            upd["kpoints"] = spec.kpoints.model_copy(update={"shift": self.kp_shift})
        for k, v in self.hidden_meta.items():
            if not getattr(spec.meta, k):
                meta_upd[k] = v
        if meta_upd:
            upd["meta"] = spec.meta.model_copy(update=meta_upd)
        return spec.model_copy(update=upd) if upd else spec

    def describe(self, spec=None) -> str:
        lines = []
        if self.continued_from:
            d = self.continued_from.get("dir", "")
            used = spec is not None and spec.meta.continued_from is not None
            vel = spec is not None and (spec.structure.velocities is not None or bool(spec.handoff and spec.handoff.velocities))
            if not used:
                lines.append(L(f"前の計算 {d} の続きでしたが、構造を作り直したので、続きの情報は使いません",
                               f"was a continuation of {d}, but the structure was rebuilt, so the restart data is not used"))
            else:
                lines.append(L(f"前の計算の続き: {d}", f"Continues the previous run: {d}")
                             + (L(" (速度も引き継ぎます)", " (velocities carried over)") if vel else L(" (速度は引き継ぎません)", " (velocities not carried over)")))
        if self.template:
            lines.append(L(f"研究室の雛形 {self.template.get('name')} から読み込みました。雛形の値のままの欄には、ラベルの下に小さく印を付けています",
                           f"Loaded from the group template {self.template.get('name')}; fields still at the template value are marked under their labels"))
        kept_shift = self.kp_shift is not None and (spec is None or (spec.kpoints is not None and tuple(spec.kpoints.shift) == self.kp_shift))
        if kept_shift:
            shown = _shift_text(self.kp_shift)
            lines.append(L(f"k 点シフトは {shown} のまま (シフトの欄には 1 つ目の成分だけ出ています)",
                           f"k-point shift kept at {shown} (the shift field shows the first component only)"))
        kept = [f"{k}={v}" for k, v in self.hidden_method.items()] if (spec is None or spec.method.code == self.hidden_code) else []
        kept += [f"meta.{k}" for k in self.hidden_meta]
        if kept:
            lines.append(L("画面に欄の無い設定は、読み込んだ値のまま: " + ", ".join(kept),
                           "Settings without a field on this page keep the loaded values: " + ", ".join(kept)))
        return "\n".join(lines)


def continuation_summary(spec) -> str:
    c = spec.meta.continued_from or {}
    out = [L(f"前の計算: {c.get('dir', '')} ({c.get('code', '')}、{c.get('task', '')})", f"Previous run: {c.get('dir', '')} ({c.get('code', '')}, {c.get('task', '')})"), ""]
    out.append(L("引き継いだもの", "Carried over"))
    out += [f"  - {x}" for x in c.get("carried", [])] or ["  -"]
    out.append(L("引き継がないもの", "Not carried over"))
    out += [f"  - {x}" for x in c.get("not_carried", [])] or [L("  (なし)", "  (none)")]
    return "\n".join(out)


def provenance_summary(prov: dict | None) -> tuple[str, list[str]]:
    if not prov:
        return "", []
    head = L(f"作成時の記録: ADIT {prov.get('adit_version')} / Python {prov.get('python')} / ASE {prov.get('ase')}",
             f"Provenance: adit {prov.get('adit_version')} / Python {prov.get('python')} / ASE {prov.get('ase')}")
    files = list(prov.get("files") or []) + list(prov.get("generated_files") or [])
    if files:
        head += L(f"、ファイルの SHA-256 {len(files)} 件", f", SHA-256 of {len(files)} files")
    lines = [f"{f['sha256']}  {f['name']}" for f in files]
    return head, lines
