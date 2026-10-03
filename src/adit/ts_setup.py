
from __future__ import annotations

from adit.errors import AditValueError
import io
import json
from pathlib import Path

from adit.lang import L
from adit.spec import CalculationSpec

MODES = ("ts", "irc", "ts+irc")
SELLA_SCRIPT = "run_sella.py"
SELLA_SETTINGS = "sella_settings.json"
SELLA_STRUCTURE = "ts_guess.extxyz"
TBLITE_METHOD = {"1": "GFN1-xTB", "2": "GFN2-xTB"}
SELLA_URL = "https://github.com/zadorlab/sella"
TBLITE_DOC = "https://tblite.readthedocs.io/en/latest/users/ase.html"


class TsError(AditValueError):
    pass


def ts_spec(spec: CalculationSpec, mode: str) -> CalculationSpec:
    if spec.method.code != "orca":
        hint = L(" (xtb と機械学習ポテンシャルは --sella で ASE + Sella のスクリプトにできます)", " (for xtb and ML potentials, --sella writes an ASE + Sella script)")
        raise TsError(L(f"--ts は ORCA だけです (いまは {spec.method.code}){hint}", f"--ts supports ORCA only (now {spec.method.code}){hint}"))
    m = spec.method
    if mode == "ts":
        upd_m, upd_t = {"ts_search": True, "irc": False}, {"type": "geometry_optimization"}
    elif mode == "irc":
        upd_m, upd_t = {"ts_search": False, "irc": True}, {"type": "single_point"}
    else:
        raise TsError(L(f"--ts は {' / '.join(MODES)} のどれかです", f"--ts must be one of {' / '.join(MODES)}"))
    return spec.model_copy(update={"method": m.model_copy(update=upd_m), "task": spec.task.model_copy(update=upd_t), "handoff": None,
                                   "meta": spec.meta.model_copy(update={"comment": (spec.meta.comment + f" {mode}").strip()})})


def write_ts(spec: CalculationSpec, cfg, out_dir: Path | str, mode: str, *, overwrite: bool = False) -> list[Path]:
    from adit.project import write_project
    from adit.stages import Stage, write_stages

    out = Path(out_dir).expanduser()
    if mode in ("ts", "irc"):
        s = ts_spec(spec, mode)
        write_project(s, cfg, out, overwrite=overwrite)
        return [out]
    if mode != "ts+irc":
        raise TsError(L(f"--ts は {' / '.join(MODES)} のどれかです", f"--ts must be one of {' / '.join(MODES)}"))
    if spec.method.code != "orca":
        ts_spec(spec, "ts")
    stages = [Stage("ts", {"task": {"type": "geometry_optimization"}, "method": {"ts_search": True, "irc": False}}, velocities=False),
              Stage("irc", {"task": {"type": "single_point"}, "method": {"ts_search": False, "irc": True}}, velocities=False)]
    return write_stages(spec, cfg, out, stages, overwrite=overwrite)


SELLA_TEMPLATE = r'''#!/usr/bin/env python3
"""adit が生成: ASE + Sella で鞍点 (遷移状態) を探し、そこから IRC を前向き・後ろ向きにたどる。設定の値は sella_settings.json にある。
要るもの: ase、sella (pip install sella)、calculator のパッケージ (tblite か mace-torch / chgnet)。ADIT は要らない。
結果: ts.extxyz (探索の終わりの構造)、irc_forward.extxyz / irc_reverse.extxyz (IRC の終わり)、*.traj (各ステップ)、results.json"""
import json
from importlib import metadata

import numpy as np
from ase.io import read, write

with open("sella_settings.json", encoding="utf-8") as f:
    S = json.load(f)


def version(pkg):
    try:
        return metadata.version(pkg)
    except metadata.PackageNotFoundError:
        return "not installed"


VERSIONS = {p: version(p) for p in ("ase", "sella", "tblite", "mace-torch", "chgnet", "torch")}
print("adit-sella: " + ", ".join(f"{k} {v}" for k, v in VERSIONS.items()), flush=True)


def make_calculator():
    c = S["calculator"]
    if c["kind"] == "tblite":
        from tblite.ase import TBLite
        return TBLite(method=c["method"], charge=c["charge"], multiplicity=c["multiplicity"], accuracy=c["accuracy"])
    fam, model, kw = c["model_family"], c["model"], {}
    if fam in ("mace_mp", "mace_off"):
        from mace.calculators import mace_mp, mace_off
        if model:
            kw["model"] = model
        if c["device"]:
            kw["device"] = c["device"]
        if c["dtype"]:
            kw["default_dtype"] = c["dtype"]
        return mace_mp(dispersion=c["dispersion"], **kw) if fam == "mace_mp" else mace_off(**kw)
    if fam == "chgnet":
        from chgnet.model.dynamics import CHGNetCalculator
        return CHGNetCalculator(use_device=c["device"] or None)
    raise SystemExit(f"unknown calculator: {c}")


def fix(a):
    if S["fixed_atoms"]:
        from ase.constraints import FixAtoms
        a.set_constraint(FixAtoms(indices=S["fixed_atoms"]))
    return a


from sella import IRC, Sella

atoms = fix(read("ts_guess.extxyz"))
atoms.calc = make_calculator()
opt = Sella(atoms, trajectory="sella.traj")
converged = opt.run(fmax=S["fmax_ev_per_ang"], steps=S["steps"])
write("ts.extxyz", atoms)
results = {"versions": VERSIONS, "ts": {"converged": bool(converged), "energy_ev": float(atoms.get_potential_energy()),
                                         "fmax_ev_per_ang": float(np.linalg.norm(atoms.get_forces(), axis=1).max())}}
if S["irc"]:
    for direction in ("forward", "reverse"):
        a = fix(read("ts.extxyz"))
        a.calc = make_calculator()
        irc = IRC(a, trajectory=f"irc_{direction}.traj")
        ok = irc.run(fmax=S["fmax_ev_per_ang"], steps=S["steps"], direction=direction)
        write(f"irc_{direction}.extxyz", a)
        results[f"irc_{direction}"] = {"converged": bool(ok), "energy_ev": float(a.get_potential_energy())}
with open("results.json", "w", encoding="utf-8") as f:
    json.dump(results, f, indent=2)
print("adit-sella: done", json.dumps(results["ts"]), flush=True)
'''


def sella_settings(spec: CalculationSpec, *, irc: bool = True) -> dict:
    m, st, t = spec.method, spec.structure, spec.task
    if t.max_steps < 1:
        raise TsError(L("Sella の最大ステップ数 (task.max_steps) は 1 以上にしてください", "the Sella max steps (task.max_steps) must be at least 1"))
    if st.fixed_axes:
        raise TsError(L("Sella のスクリプトは軸ごとの固定を書きません (原子ごとの固定は書きます)", "the Sella script does not write per-axis constraints (per-atom ones are written)"))
    if m.code == "xtb":
        if m.gfn not in TBLITE_METHOD:
            raise TsError(L("tblite で使うのは GFN1 か GFN2 です (tblite の文書の method の例)", "tblite is used with GFN1 or GFN2 (the method examples in the tblite docs)"))
        if m.solvation != "none":
            raise TsError(L("tblite の溶媒の引数の書き方を確かめていないので、溶媒ありは書きません", "the tblite solvation argument was not verified, so solvation is not written"))
        calc = {"kind": "tblite", "method": TBLITE_METHOD[m.gfn], "charge": st.charge, "multiplicity": st.multiplicity, "accuracy": m.accuracy}
    elif m.code == "mlip":
        from adit.codes.mlip import MlipGenerator
        errs = MlipGenerator().validate(spec.model_copy(update={"task": t.model_copy(update={"type": "single_point"})}), None)
        if errs:
            raise TsError("\n".join(str(e) for e in errs))
        mf = MlipGenerator()._model_file(spec)
        calc = {"kind": "mlip", "model_family": m.model_family, "model": mf.name if mf else m.model.strip(), "device": m.device.strip(),
                "dtype": m.dtype, "dispersion": m.dispersion}
    else:
        raise TsError(L(f"Sella のスクリプトの calculator は xtb (tblite) と機械学習ポテンシャルだけです (いまは {m.code}。ORCA は --ts)",
                        f"the Sella script supports xtb (tblite) and ML potentials only (now {m.code}; for ORCA use --ts)"))
    return {"calculator": calc, "fmax_ev_per_ang": t.force_tolerance_ev_per_ang, "steps": t.max_steps, "irc": irc,
            "fixed_atoms": sorted(st.fixed_atoms)}


def write_sella(spec: CalculationSpec, cfg, out_dir: Path | str, *, irc: bool = True, overwrite: bool = False) -> list[Path]:
    import shutil

    from ase.io import write

    from adit import __version__, batch
    from adit.scripts.render import render_submit

    out = Path(out_dir).expanduser()
    settings = sella_settings(spec, irc=irc)
    profile = cfg.profile(spec.runtime.profile)
    batch.check_output(out, overwrite)
    buf = io.StringIO()
    write(buf, spec.atoms, format="extxyz")
    from adit.codes.base import command_values

    template = profile.command_for("mlip" if spec.method.code == "mlip" else "python", "python3")
    try:
        exe = template.format(**command_values(spec))
    except (KeyError, IndexError, ValueError) as ex:
        raise TsError(L(f"実行コマンド {template!r} に埋められない {{...}} があります ({ex})。使えるのは {{mpiprocs}} {{ntasks}} {{omp_threads}} {{binary}} です",
                        f"the run command {template!r} has a placeholder that cannot be filled ({ex}); the available ones are {{mpiprocs}} {{ntasks}} {{omp_threads}} {{binary}}")) from ex
    pkg = "tblite" if settings["calculator"]["kind"] == "tblite" else {"chgnet": "chgnet"}.get(settings["calculator"].get("model_family"), "mace-torch")
    readme = [
        L(f"ADIT {__version__} が生成した、ASE + Sella の遷移状態の探索{'と IRC' if irc else ''}のスクリプトです (ADIT は Sella に依存しません)",
          f"ASE + Sella transition-state search{' and IRC' if irc else ''} script generated by ADIT {__version__} (ADIT does not depend on Sella)"), "",
        L("  run_sella.py        スクリプト (Sella で鞍点を探し、IRC を前向き・後ろ向きにたどる)", "  run_sella.py        the script (Sella saddle-point search, then IRC forward and reverse)"),
        L("  sella_settings.json 値 (calculator、力の収束の基準 [eV/Å]、最大ステップ数、固定原子)", "  sella_settings.json values (calculator, force criterion [eV/Å], max steps, fixed atoms)"),
        L("  ts_guess.extxyz     出発の構造 (遷移状態の見当。spec.json の構造)", "  ts_guess.extxyz     starting structure (the transition-state guess; the structure in spec.json)"),
        "", L("== 実行する前に用意すること ==", "== Before running =="),
        f"  pip install ase sella {pkg}",
        L(f"  Sella: {SELLA_URL}" + (f"   tblite の ASE calculator: {TBLITE_DOC}" if pkg == "tblite" else ""), f"  Sella: {SELLA_URL}" + (f"   tblite ASE calculator: {TBLITE_DOC}" if pkg == "tblite" else "")),
        L("  Sella の探索と IRC の細かい設定 (IRC の刻み dx など) は Sella の既定のままです。", "  Sella search and IRC details (e.g. the IRC step dx) are Sella defaults."),
        "", L("== 実行する ==", "== Run =="), "  bash submit.sh",
        "", L("== 結果 ==", "== Results =="),
        L("  ts.extxyz / irc_forward.extxyz / irc_reverse.extxyz、results.json (エネルギー [eV]、収束したか)。鞍点かどうかの判定は ADIT はしません (振動解析で確かめます)",
          "  ts.extxyz / irc_forward.extxyz / irc_reverse.extxyz, results.json (energies [eV], converged or not); ADIT does not judge whether it is a saddle point (check with a vibrational analysis)"),
    ]
    out.mkdir(parents=True, exist_ok=True)
    files = {SELLA_SCRIPT: SELLA_TEMPLATE, SELLA_SETTINGS: json.dumps(settings, indent=2, ensure_ascii=False) + "\n", SELLA_STRUCTURE: buf.getvalue(),
             "README.txt": "\n".join(readme) + "\n", "submit.sh": render_submit(spec, profile, f"{exe} {SELLA_SCRIPT} > output.log 2>&1")}
    for name, text in files.items():
        (out / name).write_text(text, encoding="utf-8", newline="\n")
    (out / "submit.sh").chmod(0o755)
    spec.save(out / "spec.json")
    if spec.method.code == "mlip":
        from adit.codes.mlip import MlipGenerator
        mf = MlipGenerator()._model_file(spec)
        if mf is not None:
            shutil.copy2(mf, out / mf.name)
    return [out]


__all__ = ["TsError", "MODES", "ts_spec", "write_ts", "sella_settings", "write_sella"]
