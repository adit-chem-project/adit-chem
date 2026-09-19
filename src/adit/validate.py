"""Mechanical validation only: what makes a run impossible to start, never an opinion about quality."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
from ase.data import atomic_numbers

from adit.config import Config
from adit.spec import CalculationSpec
from adit.lang import L
from adit.validate_types import ValidationError

__all__ = ["ValidationError", "validate", "output_dir_errors", "MIN_DISTANCE"]

MIN_DISTANCE = 0.5


def output_dir_errors(output_dir: Path | str) -> list[ValidationError]:
    target = Path(output_dir).expanduser()
    parent = target.resolve().parent
    if target.exists() and not target.is_dir():
        return [ValidationError("output_dir", L(
            f"出力先がディレクトリではありません (同じ名前のファイルがあります): {target}",
            f"the output path is not a directory (a file with that name exists): {target}"))]
    if not parent.is_dir():
        return [ValidationError("output_dir", L(
            f"親ディレクトリがありません: {parent} (打ち間違いを防ぐため、ADIT は途中のディレクトリを勝手に作りません。"
            f"`mkdir -p {parent}` で作ってから、もう一度実行してください)",
            f"the parent directory does not exist: {parent} (ADIT does not create intermediate directories, to catch typos; "
            f"create it with `mkdir -p {parent}` and run again)"))]
    return []
_WALLTIME = re.compile(r"^\d{1,3}:\d{2}:\d{2}$")
_PLACEHOLDER = re.compile(r"<[^<>\n]*[^\x00-\x7f][^<>\n]*>|<[A-Za-z][A-Za-z0-9_-]*(?: [A-Za-z0-9_-]+)+>|/path/to\b")


def validate(spec: CalculationSpec, cfg: Config, *, output_dir: Path | str | None = None) -> list[ValidationError]:
    """Return the list of errors. An empty list means the settings can be generated."""
    errs: list[ValidationError] = []
    finite = _check_finite(spec)
    if finite:
        return finite
    errs += _check_structure(spec)
    errs += _check_numbers(spec)
    errs += _check_velocities_and_handoff(spec)
    if spec.plumed is not None:
        from adit.codes.plumed import validate_plumed
        errs += validate_plumed(spec)
    if output_dir is not None:
        errs += output_dir_errors(output_dir)
    if spec.runtime.profile.strip() and spec.runtime.profile not in cfg.profiles:
        errs.append(ValidationError("runtime.profile", L(f"プロファイル {spec.runtime.profile!r} が設定ファイルにありません", f"profile {spec.runtime.profile!r} is not in the settings file")))
    elif spec.runtime.profile in cfg.profiles:
        errs += _check_placeholders(spec.runtime.profile, cfg.profiles[spec.runtime.profile], spec.method.code)
        errs += _check_limits(spec.runtime.profile, cfg.profiles[spec.runtime.profile], spec.runtime)
        jn = spec.runtime.job_name.strip()
        if cfg.profiles[spec.runtime.profile].kind != "direct" and jn and not _JOB_NAME.match(jn):
            errs.append(ValidationError("runtime.job_name", L(
                f"ジョブ名 {jn!r} は、英字で始まる英数字と _ - . だけにしてください (空白や日本語はジョブスクリプトのヘッダで使えません)",
                f"job name {jn!r}: use only letters, digits, _ - . and start with a letter (spaces and non-ASCII cannot be used in the job header)")))
    from adit.codes import GENERATORS

    gen = GENERATORS.get(spec.method.code)
    if gen is None:
        errs.append(ValidationError("method.code", L(f"対応していないコードです: {spec.method.code!r}", f"unsupported code: {spec.method.code!r}")))
    elif not errs or all(e.location != "structure.atoms" for e in errs):
        errs += gen.validate(spec, cfg)
    return errs


def _check_finite(spec: CalculationSpec) -> list[ValidationError]:
    import math

    errs: list[ValidationError] = []

    def walk(node, path: str) -> None:
        if isinstance(node, bool):
            return
        if isinstance(node, float):
            if not math.isfinite(node):
                errs.append(ValidationError(path or "input", L(f"有限の数を入れてください (いまは {node})", f"enter a finite number (now {node})")))
        elif isinstance(node, dict):
            for k, v in node.items():
                walk(v, f"{path}.{k}" if path else str(k))
        elif isinstance(node, (list, tuple)):
            for v in node:
                walk(v, path)

    walk(spec.model_dump(exclude={"meta"}), "")
    seen, out = set(), []
    for e in errs:
        if e.location not in seen:
            seen.add(e.location); out.append(e)
    return out


_JOB_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")


def _check_placeholders(name: str, profile, code: str) -> list[ValidationError]:
    items: list[tuple[str, str]] = [("select_extra", profile.select_extra)]
    items += [(f"header_extra[{i}]", v) for i, v in enumerate(profile.header_extra)]
    items += [(f"modules[{i}]", v) for i, v in enumerate(profile.modules)]
    items += [(f"code_modules.{code}[{i}]", v) for i, v in enumerate(profile.code_modules.get(code, []))]
    items += [(f"env.{k}", v) for k, v in profile.env.items()]
    if code in profile.commands:
        items.append((f"commands.{code}", profile.commands[code]))
    errs = []
    for key, value in items:
        if value and _PLACEHOLDER.search(value):
            errs.append(ValidationError("runtime.profile", L(
                f"環境設定の仮の値のままです: [profiles.{name}] の {key} = {value!r}。<…> や /path/to の部分をお使いの環境の値に書き換えるか、"
                "要らない行なら行頭に # を付けて無効にしてください",
                f"a placeholder is still in the settings: {key} = {value!r} in [profiles.{name}]. Replace the <...> or /path/to part with your own value, "
                "or put # at the start of the line if it is not needed")))
    return errs


def walltime_seconds(text: str) -> int | None:
    if not _WALLTIME.match(text.strip()):
        return None
    h, m, s = (int(x) for x in text.strip().split(":"))
    return h * 3600 + m * 60 + s


def _check_limits(name: str, profile, r) -> list[ValidationError]:
    errs: list[ValidationError] = []
    over = L("環境設定の上限を超えています", "over the limit in the settings")
    if profile.nodes_max > 0 and r.nodes > profile.nodes_max:
        errs.append(ValidationError("runtime.nodes", L(
            f"{over}: [profiles.{name}] の nodes_max = {profile.nodes_max} に対して、ノード数 {r.nodes}",
            f"{over}: nodes_max = {profile.nodes_max} in [profiles.{name}], but nodes = {r.nodes}")))
    if profile.cores_max > 0 and r.nodes * r.ncpus > profile.cores_max:
        errs.append(ValidationError("runtime.ncpus", L(
            f"{over}: [profiles.{name}] の cores_max = {profile.cores_max} に対して、ノード数 {r.nodes} × ノードあたりのコア数 {r.ncpus} = {r.nodes * r.ncpus}",
            f"{over}: cores_max = {profile.cores_max} in [profiles.{name}], but nodes {r.nodes} x cores/node {r.ncpus} = {r.nodes * r.ncpus}")))
    if profile.walltime_max.strip():
        limit = walltime_seconds(profile.walltime_max)
        if limit is None:
            errs.append(ValidationError("runtime.profile", L(
                f"環境設定の [profiles.{name}] の walltime_max = {profile.walltime_max!r} が HH:MM:SS の形ではありません",
                f"walltime_max = {profile.walltime_max!r} in [profiles.{name}] is not of the form HH:MM:SS")))
        elif (walltime_seconds(r.walltime) or 0) > limit:
            errs.append(ValidationError("runtime.walltime", L(
                f"{over}: [profiles.{name}] の walltime_max = {profile.walltime_max} に対して、制限時間 {r.walltime}",
                f"{over}: walltime_max = {profile.walltime_max} in [profiles.{name}], but walltime = {r.walltime}")))
    return errs


def _check_structure(spec: CalculationSpec) -> list[ValidationError]:
    errs: list[ValidationError] = []
    st = spec.structure
    n = len(st.atoms.symbols)
    if n == 0:
        return [ValidationError("structure.atoms", L("原子が 1 つもありません", "no atoms"))]
    unknown = [s for s in st.atoms.symbols if s not in atomic_numbers]
    if unknown:
        return [ValidationError("structure.atoms", L(f"元素記号として解釈できません: {sorted(set(unknown))}", f"not chemical symbols: {sorted(set(unknown))}"))]

    if n >= 2:
        close = _close_pairs(st)
        for i, j, dij in close[:10]:
            errs.append(
                ValidationError(
                    "structure.atoms",
                    L(
                        f"原子 {i + 1} ({st.atoms.symbols[i]}) と {j + 1} ({st.atoms.symbols[j]}) の距離が "
                        f"{dij:.3f} Å で {MIN_DISTANCE} Å 未満",
                        f"atoms {i + 1} ({st.atoms.symbols[i]}) and {j + 1} ({st.atoms.symbols[j]}) are "
                        f"{dij:.3f} Å apart, less than {MIN_DISTANCE} Å",
                    ),
                )
            )
        if len(close) > 10:
            errs.append(ValidationError("structure.atoms", L(f"他に {len(close) - 10} 組が重なっています", f"{len(close) - 10} more overlapping pairs")))
    if st.charge is not None and sum(atomic_numbers[s] for s in st.atoms.symbols) - st.charge < 0:
        errs.append(ValidationError("structure.charge", L(f"電荷 {st.charge} では電子数が負になります", f"charge {st.charge} gives a negative number of electrons")))
    if st.multiplicity < 1:
        errs.append(ValidationError("structure.multiplicity", L("1 以上が必要です", "must be at least 1")))
    return errs


def _close_pairs(st) -> list[tuple[int, int, float]]:
    return close_pairs_atoms(st.atoms.to_ase(), MIN_DISTANCE)


def close_pairs_atoms(atoms, cutoff: float = MIN_DISTANCE) -> list[tuple[int, int, float]]:
    if any(atoms.pbc):
        cell = np.asarray(atoms.cell, dtype=float)
        if abs(np.linalg.det(cell)) < 1e-6:
            return []
        from ase.neighborlist import neighbor_list
        ii, jj, dd = neighbor_list("ijd", atoms, cutoff)
        best: dict[tuple[int, int], float] = {}
        for i, j, d in zip(ii, jj, dd):
            if i < j:
                best[(int(i), int(j))] = min(float(d), best.get((int(i), int(j)), float("inf")))
        return sorted((i, j, d) for (i, j), d in best.items())
    from scipy.spatial import cKDTree
    pos = np.asarray(atoms.positions, dtype=float)
    pairs = cKDTree(pos).query_pairs(cutoff, output_type="ndarray")
    return sorted((int(i), int(j), float(np.linalg.norm(pos[i] - pos[j]))) for i, j in pairs)


MAX_BAND_POINTS = 10_000
MAX_KPOINTS = 1_000_000


def _check_numbers(spec: CalculationSpec) -> list[ValidationError]:
    errs: list[ValidationError] = []
    r = spec.runtime
    for name in ("nodes", "ncpus", "mpiprocs", "omp_threads"):
        if getattr(r, name) < 1:
            errs.append(ValidationError(f"runtime.{name}", L("1 以上が必要です", "must be at least 1")))
    if not _WALLTIME.match(r.walltime):
        errs.append(ValidationError("runtime.walltime", L(f"HH:MM:SS の形ではありません: {r.walltime!r}", f"not of the form HH:MM:SS: {r.walltime!r}")))
    if not r.job_name.strip():
        errs.append(ValidationError("runtime.job_name", L("ジョブ名が空です", "job name is empty")))
    if not r.profile.strip():
        errs.append(ValidationError("runtime.profile", L("実行先プロファイルが選ばれていません", "no target profile selected")))
    t = spec.task
    if t.type == "geometry_optimization":
        if t.max_steps < 0:
            errs.append(ValidationError("task.max_steps", L("負のステップ数は指定できません", "negative number of steps is not allowed")))
        if t.force_tolerance_ev_per_ang <= 0:
            errs.append(ValidationError("task.force_tolerance_ev_per_ang", L("0 より大きい値が必要です", "must be greater than 0")))
    if t.type == "band_structure":
        if not spec.structure.periodic:
            errs.append(ValidationError("task.type", L("バンド計算は周期系だけです", "band structure is only for periodic systems")))
        if t.bands.npoints < 2:
            errs.append(ValidationError("task.bands.npoints", L("2 以上が必要です", "must be at least 2")))
        if t.bands.empty_bands < 0:
            errs.append(ValidationError("task.bands.empty_bands", L("0 以上が必要です", "must be at least 0")))
        if t.bands.path.strip() and spec.structure.periodic:
            errs += _check_band_path(spec)
        if t.bands.npoints > MAX_BAND_POINTS:
            errs.append(ValidationError("task.bands.npoints", L(
                f"{t.bands.npoints} 点は多すぎます ({MAX_BAND_POINTS} 点まで。生成に使うメモリが大きくなりすぎるため)。バンド図にはふつう数十〜数百点で足ります",
                f"{t.bands.npoints} points is too many (max {MAX_BAND_POINTS}; generation would use too much memory). A band plot usually needs tens to hundreds")))
    if t.type == "molecular_dynamics":
        md = t.md
        if md.steps < 1:
            errs.append(ValidationError("task.md.steps", L("1 以上が必要です", "must be at least 1")))
        if md.timestep_fs <= 0:
            errs.append(ValidationError("task.md.timestep_fs", L("0 より大きい値が必要です", "must be greater than 0")))
        if md.temperature_k < 0:
            errs.append(ValidationError("task.md.temperature_k", L("負の温度は指定できません", "negative temperature is not allowed")))
        if md.dump_interval < 1:
            errs.append(ValidationError("task.md.dump_interval", L("1 以上が必要です", "must be at least 1")))
        if md.coupling_time_fs <= 0:
            errs.append(ValidationError("task.md.coupling_time_fs", L("0 より大きい値が必要です", "must be greater than 0")))
        if md.ensemble == "NPT" and not spec.structure.periodic:
            errs.append(ValidationError("task.md.ensemble", L("NPT は周期系だけです", "NPT is only for periodic systems")))
        if md.ensemble == "NPT" and md.barostat_time_fs <= 0:
            errs.append(ValidationError("task.md.barostat_time_fs", L("0 より大きい値が必要です", "must be greater than 0")))
    errs += _check_periodic(spec)
    return errs


def _check_band_path(spec: CalculationSpec) -> list[ValidationError]:
    atoms = spec.structure.atoms.to_ase()
    try:
        if abs(np.linalg.det(np.asarray(atoms.cell, dtype=float))) < 1e-6:
            return []
        from adit.bandpath import band_path
        band_path(atoms, spec.task.bands.path.strip(), 2)
        return []
    except Exception:
        try:
            names = " ".join(sorted(atoms.cell.bandpath(npoints=2).special_points))
        except Exception:
            names = ""
        return [ValidationError("task.bands.path", L(
            f"経路 {spec.task.bands.path!r} を読めません。この格子で使える記号: {names} (記号をつなげて書き、途切れる所は , で区切ります。空欄なら標準の経路)",
            f"cannot read the path {spec.task.bands.path!r}. Symbols for this lattice: {names} (write them in a row, use , for a jump; leave empty for the standard path)"))]


def _check_periodic(spec: CalculationSpec) -> list[ValidationError]:
    errs: list[ValidationError] = []
    st, kp, t = spec.structure, spec.kpoints, spec.task
    n = len(st.atoms.symbols)
    bad = [i for i in st.fixed_atoms if i < 0 or i >= n]
    if bad:
        errs.append(ValidationError("structure.fixed_atoms", L(f"原子の番号が範囲外です (0〜{n - 1}): {bad}", f"atom index out of range (0..{n - 1}): {bad}")))
    bad_axes = [k for k in st.fixed_axes if not k.isdigit() or int(k) < 0 or int(k) >= n]
    if bad_axes:
        errs.append(ValidationError("structure.fixed_axes", L(f"原子の番号が範囲外です (0〜{n - 1}): {bad_axes}", f"atom index out of range (0..{n - 1}): {bad_axes}")))
    from adit.codes import GENERATORS
    gen = GENERATORS.get(spec.method.code)
    uses_k = gen is None or gen.uses_kpoints
    if not uses_k:
        kp = None
    if st.periodic:
        if kp is None and uses_k:
            errs.append(ValidationError("kpoints", L("周期系なので k 点の指定が必要です", "periodic system: k-points are required")))
        elif kp is not None:
            if kp.mode == "mesh" and any(k < 1 for k in kp.mesh):
                errs.append(ValidationError("kpoints.mesh", L(f"分割数は 1 以上にしてください: {kp.mesh}", f"mesh must be at least 1: {kp.mesh}")))
            if kp.mode == "density" and kp.density <= 0:
                errs.append(ValidationError("kpoints.density", L("0 より大きい値が必要です", "must be greater than 0")))
            if any(s not in (0.0, 0.5) for s in kp.shift):
                errs.append(ValidationError("kpoints.shift", L(f"シフトは 0 か 0.5 にしてください: {kp.shift}", f"shift must be 0 or 0.5: {kp.shift}")))
        cell = np.asarray(st.atoms.cell, dtype=float)
        if abs(np.linalg.det(cell)) < 1e-6:
            errs.append(ValidationError("structure.atoms", L("周期系なのに格子ベクトルが退化しています (体積 0)", "periodic system but the cell is degenerate (zero volume)")))
        elif kp is not None and (kp.mode != "density" or kp.density > 0) and (kp.mode != "mesh" or all(k >= 1 for k in kp.mesh)):
            mesh = kp.resolved_mesh(cell)
            total = int(np.prod([float(k) for k in mesh]))
            if total > MAX_KPOINTS:
                where = "kpoints.density" if kp.mode == "density" else "kpoints.mesh"
                errs.append(ValidationError(where, L(
                    f"k 点が {' × '.join(str(k) for k in mesh)} = {total:.3g} 点になり、多すぎます ({MAX_KPOINTS:.0e} 点まで。生成と計算に使うメモリが大きくなりすぎるため)",
                    f"{' x '.join(str(k) for k in mesh)} = {total:.3g} k-points is too many (max {MAX_KPOINTS:.0e}; generation and the calculation would use too much memory)")))
    else:
        if kp is not None:
            errs.append(ValidationError("kpoints", L("分子系 (非周期) に k 点は指定できません", "k-points cannot be given for a molecule (non-periodic)")))
        if t.relax_cell != "no":
            errs.append(ValidationError("task.relax_cell", L("分子系 (非周期) では格子を動かせません", "the cell cannot be relaxed for a molecule (non-periodic)")))
    return errs


_ORBITAL = re.compile(r"^([1-9])([spdf])$")


def orbital_l(orbital: str) -> int:
    return "spdf".index(orbital.strip().lower()[-1])


def check_hubbard(hubbard: dict, elements: list[str], location: str) -> list[ValidationError]:
    errs: list[ValidationError] = []
    absent = [e for e in hubbard if e not in elements]
    if absent:
        errs.append(ValidationError(location, L(f"U を付けた元素が構造にありません: {absent}", f"elements with U are not in the structure: {absent}")))
    bad = [f"{e}: {h.orbital!r}" for e, h in hubbard.items() if not _ORBITAL.match(h.orbital.strip().lower())]
    if bad:
        errs.append(ValidationError(location, L(f"殻は主量子数と s / p / d / f で書いてください (例 3d、4f): {bad}",
                                                f"write the shell as principal quantum number + s / p / d / f (e.g. 3d, 4f): {bad}")))
    return errs


def _check_velocities_and_handoff(spec: CalculationSpec) -> list[ValidationError]:
    errs: list[ValidationError] = []
    st, h, t = spec.structure, spec.handoff, spec.task
    from adit.codes import GENERATORS

    gen = GENERATORS.get(spec.method.code)
    if st.velocities is not None:
        if len(st.velocities) != len(st.atoms.symbols):
            errs.append(ValidationError("structure.velocities", L(f"速度の数 {len(st.velocities)} が原子数 {len(st.atoms.symbols)} と違います",
                                                                  f"{len(st.velocities)} velocities for {len(st.atoms.symbols)} atoms")))
        if t.type != "molecular_dynamics":
            errs.append(ValidationError("structure.velocities", L("速度は分子動力学 (MD) でしか使われません (計算の種類を MD にするか、速度を外してください)",
                                                                  "velocities are used only in molecular dynamics (choose MD or drop the velocities)")))
        elif gen is not None and not gen.writes_velocities:
            errs.append(ValidationError("structure.velocities", L(f"この計算コードの生成器 ({spec.method.code}) は速度を入力に書けません (書けるのは DFTB+ と VASP)",
                                                                  f"the {spec.method.code} generator cannot write velocities (DFTB+ and VASP can)")))
    if h is None:
        return errs
    if h.velocities and (t.type != "molecular_dynamics" or h.previous_task != "molecular_dynamics"):
        errs.append(ValidationError("handoff", L("前の計算の速度を引き継げるのは、前も今回も分子動力学 (MD) のときだけです",
                                                 "velocities can be carried over only from MD to MD")))
    if (h.velocities or h.files) and h.previous_code != spec.method.code:
        errs.append(ValidationError("handoff", L(f"前の計算 ({h.previous_code}) と計算コードが違うので、続きの情報は使えません (構造だけなら使えます)",
                                                 f"the previous run used {h.previous_code}, so its restart data cannot be used (the structure alone can)")))
    if not h.at_run:
        missing = [src for src in h.files.values() if not (Path(h.previous_dir).expanduser() / src).is_file()]
        if missing:
            errs.append(ValidationError("handoff", L(f"前の計算のディレクトリ {h.previous_dir} に、引き継ぐファイルがありません: {missing}",
                                                     f"files to carry over are missing from the previous directory {h.previous_dir}: {missing}")))
    return errs


def electron_parity_error(n_electrons: int, charge: int, multiplicity: int, *, location: str = "structure.multiplicity") -> ValidationError | None:
    unpaired = multiplicity - 1
    if n_electrons < 0:
        return None
    if unpaired > n_electrons:
        return ValidationError(location, L(f"不対電子 {unpaired} が電子数 {n_electrons} を超えています", f"{unpaired} unpaired electrons exceed the {n_electrons} electrons"))
    if (n_electrons - unpaired) % 2 != 0:
        hint = L("電子数が奇数なら多重度は 2, 4, … (不対電子が奇数個)、偶数なら 1, 3, …。周期系 (結晶) で基本セルの電子数が奇数なら、"
                 "多重度を 2 にするか、セルを 2 倍にして偶数にする方法があります",
                 "odd electron count → multiplicity 2, 4, …; even → 1, 3, …. For a crystal whose cell has an odd electron count, "
                 "set multiplicity 2 or double the cell")
        return ValidationError(location, L(f"電子数 {n_electrons} (電荷 {charge}) に対して多重度 {multiplicity} は成り立ちません (偶奇が合いません)。{hint}",
                                           f"multiplicity {multiplicity} is impossible for {n_electrons} electrons (charge {charge}): parity mismatch. {hint}"))
    return None
