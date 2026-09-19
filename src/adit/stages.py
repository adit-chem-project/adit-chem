
from __future__ import annotations

from adit.errors import AditValueError
import json
import re
import shutil
from dataclasses import dataclass, field
from pathlib import Path

from adit import handoff as hf
from adit.lang import L, pick
from adit.spec import CalculationSpec, Handoff

STAGES_FILE = "stages.json"
HANDOFF_SCRIPT = "handoff.py"
ALLOWED = ("task", "method", "runtime")
MD = "molecular_dynamics"
CARRY_VELOCITIES = ("dftbplus", "vasp", "xtb", "cp2k", "lammps", "gromacs")


class StageError(AditValueError):
    pass


@dataclass(frozen=True)
class Stage:
    name: str
    overrides: dict
    velocities: bool | None = None


@dataclass
class PlannedStage:
    dir: str
    spec: CalculationSpec
    pre_command: str = ""
    readme: list[str] = field(default_factory=list)
    uses_script: bool = False


def parse_stages(data) -> list[Stage]:
    items = data.get("stages") if isinstance(data, dict) else data
    if not isinstance(items, list) or not items:
        raise StageError(L('段階の一覧がありません。{"stages": [{"name": …, "task": {…}}, …]} の形で書いてください',
                           'no stages; write them as {"stages": [{"name": ..., "task": {...}}, ...]}'))
    out = []
    for i, it in enumerate(items, 1):
        if not isinstance(it, dict):
            raise StageError(L(f"{i} 段階目が {{…}} の形ではありません", f"stage {i} is not an object {{...}}"))
        # "dir" is what write_stages records, so its own stages.json reads back
        bad = [k for k in it if k not in ("name", "velocities", "dir", *ALLOWED)]
        if bad:
            raise StageError(L(f"{i} 段階目に使えない項目があります: {bad} (使えるのは name / task / method / runtime / velocities)",
                               f"stage {i} has unknown items: {bad} (allowed: name / task / method / runtime / velocities)"))
        vel = it.get("velocities")
        if vel is not None and not isinstance(vel, bool):
            raise StageError(L(f"{i} 段階目の velocities は true か false にしてください", f"velocities of stage {i} must be true or false"))
        for k in ALLOWED:
            if k in it and not isinstance(it[k], dict):
                raise StageError(L(f"{i} 段階目の {k} は {{…}} の形にしてください", f"{k} of stage {i} must be an object {{...}}"))
        out.append(Stage(str(it.get("name") or f"stage{i}"), {k: it[k] for k in ALLOWED if k in it}, vel))
    return out


def load_stages(path: Path | str) -> list[Stage]:
    try:
        return parse_stages(json.loads(Path(path).read_text(encoding="utf-8")))
    except json.JSONDecodeError as ex:
        raise StageError(L(f"{path} を JSON として読めません: {ex}", f"cannot read {path} as JSON: {ex}")) from ex


def _merge(a: dict, b: dict) -> dict:
    out = dict(a)
    for k, v in b.items():
        out[k] = _merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def stage_dir(i: int, name: str) -> str:
    safe = re.sub(r"[^0-9A-Za-z_-]+", "_", name).strip("_") or "stage"
    return f"stage_{i:02d}_{safe}"


def _describe(s: CalculationSpec) -> str:
    t = s.task
    if t.type == MD:
        md = t.md
        return L(f"MD ({md.ensemble}、{md.temperature_k:g} K、{md.steps} ステップ × {md.timestep_fs:g} fs、熱浴 {md.thermostat})",
                 f"MD ({md.ensemble}, {md.temperature_k:g} K, {md.steps} steps x {md.timestep_fs:g} fs, thermostat {md.thermostat})")
    return pick({"geometry_optimization": "構造最適化", "single_point": "一点計算", "vibrations": "振動解析", "band_structure": "バンド計算"},
                {"geometry_optimization": "geometry optimization", "single_point": "single point", "vibrations": "vibrations",
                 "band_structure": "band structure"}).get(t.type, t.type)


def plan_stages(spec: CalculationSpec, stages: list[Stage]) -> list[PlannedStage]:
    from pydantic import ValidationError as PydanticError

    from adit.validate_types import friendly_pydantic

    code = spec.method.code
    base = spec.model_dump(mode="json")
    keep = base.get("handoff") if (spec.handoff is not None and not spec.handoff.at_run) else None
    base["handoff"] = None
    n = len(stages)
    out: list[PlannedStage] = []
    for i, stg in enumerate(stages, 1):
        d = stage_dir(i, stg.name)
        data = _merge(base, stg.overrides)
        if i == 1 and keep is not None:
            data["handoff"] = keep
        if data["method"].get("code") != code:
            raise StageError(L(f"{i} 段階目 ({stg.name}) で計算コードを変えることはできません ({code} のまま)", f"stage {i} ({stg.name}) cannot change the code (stays {code})"))
        if data["runtime"].get("profile") != spec.runtime.profile:
            raise StageError(L(f"{i} 段階目 ({stg.name}) でプロファイルを変えることはできません (全段階で {spec.runtime.profile})",
                               f"stage {i} ({stg.name}) cannot change the profile (all stages use {spec.runtime.profile})"))
        data["meta"] = dict(data.get("meta") or {},
                            stage={"index": i, "count": n, "name": stg.name, "dir": d},
                            comment=(str((data.get("meta") or {}).get("comment", "")) + f" stage {i}/{n} {stg.name}").strip())
        try:
            s = CalculationSpec.model_validate(data)
        except PydanticError as ex:
            raise StageError(L(f"{i} 段階目 ({stg.name}): ", f"stage {i} ({stg.name}): ") + friendly_pydantic(ex)) from ex
        plan = PlannedStage(d, s)
        if i == 1:
            plan.readme = [L("== 段階 ==", "== Stage =="),
                           L(f"  {n} 段階のうち 1 段階目 ({stg.name}: {_describe(s)})。spec.json の構造から始めます。",
                             f"  Stage 1 of {n} ({stg.name}: {_describe(s)}); starts from the structure in spec.json.")]
        else:
            prev = out[-1]
            ptask, both_md = prev.spec.task.type, prev.spec.task.type == MD and s.task.type == MD
            vel = stg.velocities if stg.velocities is not None else (both_md and code in CARRY_VELOCITIES)
            if vel and not both_md:
                raise StageError(L(f"{i} 段階目 ({stg.name}): 速度を引き継げるのは、前の段階もこの段階も MD のときだけです",
                                   f"stage {i} ({stg.name}): velocities can be carried over only from MD to MD"))
            if vel and code not in CARRY_VELOCITIES:
                raise StageError(L(f"{i} 段階目 ({stg.name}): {code} は出力に速度を書かないので引き継げません (\"velocities\": false にしてください)",
                                   f"stage {i} ({stg.name}): {code} does not write velocities, so they cannot be carried over (set \"velocities\": false)"))
            files = hf.copy_files(code, ptask, vel)
            s = s.model_copy(update={"handoff": Handoff(previous_dir=f"../{prev.dir}", previous_task=ptask, previous_code=code,
                                                         at_run=True, velocities=vel, files=files)})
            plan.spec = s
            if code in hf.CONVERTS:
                args = [f"python3 ../{HANDOFF_SCRIPT}", code, f"../{prev.dir}", ptask]
                if vel:
                    args.append("--velocities")
                if files:
                    args.append(f"--files '{json.dumps(files)}'")
                plan.pre_command, plan.uses_script = " ".join(args), True
            elif files:
                plan.pre_command = " && ".join(f"cp ../{prev.dir}/{src} {dest}" for dest, src in files.items())
            what = (L(f"{', '.join(f'{src} → {dest}' for dest, src in files.items())} を写し", f"copies {', '.join(f'{src} -> {dest}' for dest, src in files.items())}")
                    if files else "")
            conv = L("構造を前の段階の最終構造に書き換え", "rewrites the structure with the final structure of the previous stage") if code in hf.CONVERTS else ""
            plan.readme = [
                L("== 段階 ==", "== Stage =="),
                L(f"  {n} 段階のうち {i} 段階目 ({stg.name}: {_describe(s)})。前の段階 ../{prev.dir} が終わってから実行します。",
                  f"  Stage {i} of {n} ({stg.name}: {_describe(s)}); run it after the previous stage ../{prev.dir} has finished."),
                L(f"  submit.sh は計算の前に、前の段階の出力から{what}{'、' if what and conv else ''}{conv}ます"
                  f"{' (python3 と ../handoff.py を使います)' if plan.uses_script else ''}。前の段階の出力が無ければ、そこで止まります。",
                  f"  Before the calculation, submit.sh {what}{'; ' if what and conv else ''}{conv}"
                  f"{' (uses python3 and ../handoff.py)' if plan.uses_script else ''}; it stops if the previous output is missing."),
                (L("  速度 (MD の続きの情報) も前の段階から引き継ぎます。", "  Velocities (MD restart data) are also carried over from the previous stage.") if vel else
                 L("  速度は引き継ぎません (MD なら新しく初速を作ります)。", "  Velocities are not carried over (MD generates new initial velocities).")),
                L("  spec.json の構造は 1 段階目の最初の構造のままです。実際の出発点は前の段階の出力です。",
                  "  The structure in spec.json is still the initial structure of stage 1; the actual starting point is the previous stage's output."),
            ]
        out.append(plan)
    return out


def _top_readme(plan: list[PlannedStage], profile) -> str:
    from adit import __version__

    s0 = plan[0].spec
    lines = [L(f"ADIT {__version__} が生成した、段階に分けた計算です ({s0.method.code}、{len(plan)} 段階)", f"Staged calculation generated by ADIT {__version__} ({s0.method.code}, {len(plan)} stages)"), ""]
    for p in plan:
        h = p.spec.handoff
        tail = (L("  ← 前の段階の最終構造" + ("と速度" if h.velocities else "") + "から", "  <- from the final structure" + (" and velocities" if h.velocities else "") + " of the previous stage")
                if h else "")
        lines.append(f"  {p.dir:<24} {_describe(p.spec)}{tail}")
    lines += ["", L("各段階のディレクトリは、通常の ADIT の生成したファイルと同じ形です (中の README.txt に、その段階の入力と出力の説明があります)。",
                    "Each stage directory has the usual ADIT layout (its README.txt explains its inputs and outputs)."),
              L("2 段階目からの submit.sh は、実行する直前に前の段階の出力から構造 (と速度) を写します。前の段階が終わっていなければ止まります。",
                "From stage 2 on, submit.sh copies the structure (and velocities) from the previous stage right before running; it stops if that stage has not finished.")]
    if any(p.uses_script for p in plan):
        lines.append(L("そのとき python3 で handoff.py (このディレクトリにあります。標準ライブラリだけで動きます) を使います。",
                       "It uses handoff.py (in this directory; standard library only) with python3."))
    lines += ["", L("== この PC で実行する ==", "== Run on this PC ==")]
    if profile.kind == "direct":
        lines += [L("  bash submit.sh   (段階を順に実行します。途中の段階が失敗したら、そこで止まります)", "  bash submit.sh   (runs the stages in order; stops at the first stage that fails)")]
    else:
        lines += [L("  この生成したファイルはクラスタのプロファイル用です。この PC で実行するなら、プロファイルを direct のものにして生成し直します。",
                    "  These files are for a cluster profile; to run on this PC, generate again with a direct profile.")]
    lines += ["", L("== クラスタで実行する (依存付きの投入の例。ADIT は投入しません) ==", "== Running on a cluster (example of dependent submission; ADIT does not submit) =="),
              L("  各段階の submit.sh をクラスタのプロファイルで生成したうえで、このディレクトリで次のように投入すると、",
                "  With the stage submit.sh files generated for a cluster profile, submitting as below from this directory"),
              L("  前の段階が正常に終わったときだけ次の段階が走り出します。", "  starts each stage only after the previous one has finished successfully."),
              "  PBS:"]
    for i, p in enumerate(plan, 1):
        dep = f" -W depend=afterok:$j{i - 1}" if i > 1 else ""
        lines.append(f"    j{i}=$(cd {p.dir} && qsub{dep} submit.sh)")
    lines.append("  Slurm:")
    for i, p in enumerate(plan, 1):
        dep = f" --dependency=afterok:$j{i - 1}" if i > 1 else ""
        lines.append(f"    j{i}=$(cd {p.dir} && sbatch --parsable{dep} submit.sh)")
    return "\n".join(lines) + "\n"


def _top_submit(plan: list[PlannedStage]) -> str:
    from adit import __version__

    dirs = " ".join(p.dir for p in plan)
    return (f"#!/bin/bash\n"
            + L(f"# ADIT {__version__} が生成。段階を順に実行する (前の段階が失敗したらそこで止まる)\n",
                f"# generated by ADIT {__version__}. Runs the stages in order (stops if a stage fails)\n")
            + 'cd "$(dirname "$0")"\nset -e\n'
            + f"for d in {dirs}; do\n  echo \"== $d\"\n  bash \"$d/submit.sh\"\ndone\n")


def write_stages(spec: CalculationSpec, cfg, out_dir: Path | str, stages: list[Stage], *, overwrite: bool = False) -> list[Path]:
    from adit.batch import check_output
    from adit.project import build_project, write_project

    out = Path(out_dir).expanduser()
    plan = plan_stages(spec, stages)
    check_output(out, overwrite)
    out.mkdir(parents=True, exist_ok=True)
    for p in plan:
        build_project(p.spec, cfg, output_dir=out / p.dir, pre_command=p.pre_command, extra_readme=p.readme)
    dirs = []
    for p in plan:
        write_project(p.spec, cfg, out / p.dir, overwrite=overwrite, pre_command=p.pre_command, extra_readme=p.readme)
        dirs.append(out / p.dir)
    if any(p.uses_script for p in plan):
        shutil.copyfile(Path(hf.__file__), out / HANDOFF_SCRIPT)
    (out / STAGES_FILE).write_text(json.dumps({"stages": [
        {"dir": p.dir, "name": p.spec.meta.stage["name"], "task": p.spec.task.model_dump(mode="json"),
         "velocities": bool(p.spec.handoff.velocities) if p.spec.handoff else False} for p in plan]}, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    profile = cfg.profile(spec.runtime.profile)
    if profile.kind == "direct":
        (out / "submit.sh").write_text(_top_submit(plan), encoding="utf-8", newline="\n")
        try:
            (out / "submit.sh").chmod(0o755)
        except OSError:
            pass
    (out / "README.txt").write_text(_top_readme(plan, profile), encoding="utf-8", newline="\n")
    return dirs


__all__ = ["Stage", "StageError", "parse_stages", "load_stages", "plan_stages", "write_stages", "stage_dir"]
