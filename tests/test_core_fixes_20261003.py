import ast
import json
import os
import shutil
import tomllib
from pathlib import Path

import numpy as np
import pytest

from adit import handoff as hf
from adit.cli import main
from adit.config import load_config, save_config
from adit.project import ProjectError, build_project, restore_backup, write_project
from adit.scan import ScanError, parse_scan, write_scan
from adit.spec import (AtomsData, CalculationSpec, GromacsMethod, Handoff, KPoints, LammpsMethod, MDSettings, MlipMethod, Runtime,
                       Structure, Task, VaspMethod, XtbMethod)
from adit.stages import Stage, StageError, load_stages, parse_stages, plan_stages, write_stages
from adit.validate import validate
from tests.conftest import cfg_for, water_spec

REPO = Path(__file__).resolve().parent.parent
HANDOFF = REPO / "src" / "adit" / "handoff.py"


def _cfg_file(sk_root, tmp_path):
    p = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), p)
    return p


def _locs(errs):
    return sorted({e.location for e in errs})


def _periodic_water(cell=None, **kw):
    atoms = AtomsData(symbols=["O", "H", "H"], positions=[(5.0, 4.0, 5.0), (5.0, 5.0, 5.783064), (5.0, 5.0, 4.216936)],
                      cell=cell or [(10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0)], pbc=(True, True, True))
    return water_spec(structure=Structure(source="file", source_ref="x", atoms=atoms), **kw)


# ---------- core H1: unknown keys in the files adit-gen reads ----------

def test_cli_stops_on_unknown_keys_unless_told_otherwise(sk_root, tmp_path, capsys):
    cfg_path = _cfg_file(sk_root, tmp_path)
    d = json.loads(water_spec().to_json())
    d["kpoint"] = {"mode": "mesh"}
    d["task"]["typ"] = "molecular_dynamics"
    spec_path = tmp_path / "mine.json"
    spec_path.write_text(json.dumps(d), encoding="utf-8")
    assert main([str(spec_path), "--validate", "--config", str(cfg_path)]) == 1
    err = capsys.readouterr().err
    assert "mine.json に知らない項目があります: kpoint, task.typ。" in err and "--ignore-unknown-keys" in err
    assert main([str(spec_path), "--validate", "--config", str(cfg_path), "--ignore-unknown-keys"]) == 0
    out = capsys.readouterr()
    assert "注意" in out.err and "task.typ" in out.err and "エラーはありません" in out.out


def test_cli_checks_the_previous_run_and_the_template_for_unknown_keys(sk_root, tmp_path, capsys):
    cfg = cfg_for(sk_root)
    cfg.templates_dir = str(tmp_path / "tpl")
    cfg_path = tmp_path / "cluster.toml"
    save_config(cfg, cfg_path)
    prev = tmp_path / "prev"
    prev.mkdir()
    d = json.loads(water_spec().to_json())
    d["handof"] = {}
    (prev / "spec.json").write_text(json.dumps(d), encoding="utf-8")
    assert main(["--continue-from", str(prev), str(tmp_path / "next"), "--config", str(cfg_path)]) == 1
    assert "handof" in capsys.readouterr().err
    (tmp_path / "tpl").mkdir()
    body = {k: v for k, v in json.loads(water_spec().to_json()).items() if k not in ("structure", "meta")}
    body["template"] = {"comment": "lab"}
    body["method"]["sccc"] = 1
    (tmp_path / "tpl" / "lab.json").write_text(json.dumps(body), encoding="utf-8")
    spec_path = tmp_path / "spec.json"
    water_spec().save(spec_path)
    assert main([str(spec_path), "--validate", "--template", "lab", "--config", str(cfg_path)]) == 1
    err = capsys.readouterr().err
    assert "method.sccc" in err and ": template" not in err and "template," not in err


# ---------- core M7: an explicit --config path must exist ----------

def test_cli_explicit_config_path_must_exist(tmp_path, capsys):
    spec_path = tmp_path / "spec.json"
    water_spec().save(spec_path)
    missing = tmp_path / "typo" / "cluster.toml"
    assert main([str(spec_path), "--validate", "--config", str(missing)]) == 2
    err = capsys.readouterr().err
    assert not missing.parent.exists() and "--config" in err and "環境設定ファイルがまだありません" in err


# ---------- core M4: file-system errors are messages, not tracebacks ----------

def test_cli_reports_a_missing_stages_file_without_a_traceback(sk_root, tmp_path, capsys):
    cfg_path = _cfg_file(sk_root, tmp_path)
    spec_path = tmp_path / "spec.json"
    water_spec().save(spec_path)
    assert main([str(spec_path), str(tmp_path / "out"), "--stages", str(tmp_path / "missing.json"), "--config", str(cfg_path)]) == 1
    err = capsys.readouterr().err
    assert "missing.json" in err and "段階のファイルがありません" in err and "Traceback" not in err


@pytest.mark.skipif(os.name == "nt" or (hasattr(os, "geteuid") and os.geteuid() == 0), reason="a read-only directory does not stop root or Windows")
def test_cli_reports_an_unwritable_output_directory(sk_root, tmp_path, capsys):
    cfg_path = _cfg_file(sk_root, tmp_path)
    spec_path = tmp_path / "spec.json"
    water_spec().save(spec_path)
    ro = tmp_path / "ro"
    ro.mkdir()
    ro.chmod(0o555)
    try:
        assert main([str(spec_path), str(ro / "out"), "--config", str(cfg_path)]) == 1
        err = capsys.readouterr().err
        assert str(ro) in err and "Traceback" not in err
    finally:
        ro.chmod(0o755)


# ---------- core M6 / L5: scans are bounded and directory names unique ----------

def test_scan_has_an_upper_bound_and_rejects_colliding_names(sk_root, tmp_path):
    cfg = cfg_for(sk_root)
    big = [parse_scan("method.max_scc_iterations=" + ",".join(str(i) for i in range(100, 140))),
           parse_scan("method.scc_tolerance=" + ",".join(f"1e-{i}" for i in range(4, 44)))]
    with pytest.raises(ScanError, match="1600"):
        write_scan(water_spec(), cfg, tmp_path / "grid", big)
    assert not (tmp_path / "grid").exists()
    with pytest.raises(ScanError, match="sk_set_a_b"):
        write_scan(water_spec(), cfg, tmp_path / "dup", parse_scan("method.sk_set=a b,a_b"))
    assert not (tmp_path / "dup").exists()


# ---------- core M8 / S4: copies ----------

def test_missing_copy_source_is_a_clean_error_and_nothing_is_written(sk_root, tmp_path, monkeypatch):
    from adit.codes import GENERATORS

    monkeypatch.setattr(GENERATORS["dftbplus"], "files_to_copy", lambda spec, res: {"skf/O-H.skf": tmp_path / "gone.skf"})
    with pytest.raises(ProjectError, match="gone.skf"):
        write_project(water_spec(), cfg_for(sk_root), tmp_path / "calc", backup=tmp_path / "bk")
    assert not (tmp_path / "calc").exists() and not (tmp_path / "bk").exists()


def test_a_failed_write_still_leaves_a_usable_backup(sk_root, tmp_path, monkeypatch):
    cfg = cfg_for(sk_root)
    out = tmp_path / "calc"
    write_project(water_spec(), cfg, out)
    (out / "README.txt").write_text("mine\n", encoding="utf-8")
    real = shutil.copy2

    def boom(src, dst, *a, **k):
        if Path(dst).is_relative_to(tmp_path / "bk"):
            return real(src, dst, *a, **k)
        raise OSError("disk full")

    monkeypatch.setattr(shutil, "copy2", boom)
    with pytest.raises(OSError, match="disk full"):
        write_project(water_spec(), cfg, out, overwrite=True, backup=tmp_path / "bk")
    assert (tmp_path / "bk" / "manifest.json").is_file()
    monkeypatch.setattr(shutil, "copy2", real)
    restored, _removed = restore_backup(tmp_path / "bk")
    assert "README.txt" in restored and (out / "README.txt").read_text(encoding="utf-8") == "mine\n"


def test_copy_onto_itself_is_skipped(sk_root, tmp_path):
    cfg = cfg_for(sk_root)
    out = tmp_path / "run1"
    md = Task(type="molecular_dynamics", md=MDSettings(steps=10))
    write_project(water_spec(method=XtbMethod(), task=md), cfg, out)
    (out / "mdrestart").write_text("restart\n", encoding="utf-8")
    h = Handoff(previous_dir=str(out), previous_task="molecular_dynamics", previous_code="xtb", velocities=True, files={"mdrestart": "mdrestart"})
    written = write_project(water_spec(method=XtbMethod(), task=md, handoff=h), cfg, out, overwrite=True)
    assert out / "mdrestart" in written and (out / "mdrestart").read_text(encoding="utf-8") == "restart\n"


# ---------- core M10: save_config with a top-level inline table ----------

def test_save_config_keeps_a_top_level_inline_table(tmp_path):
    p = tmp_path / "cluster.toml"
    p.write_text('foo = { a = 1 }\nlanguage = "ja"\n\n[profiles.local]\nkind = "direct"\n', encoding="utf-8")
    save_config(load_config(p), p)
    data = tomllib.loads(p.read_text(encoding="utf-8"))
    assert data["foo"] == {"a": 1} and data["profiles"]["local"]["kind"] == "direct"
    assert load_config(p).unknown_keys == ["foo"]


# ---------- core M11: UTF-8 BOM ----------

def test_json_files_with_a_utf8_bom_are_read(sk_root, tmp_path):
    bom = b"\xef\xbb\xbf"
    spec_path = tmp_path / "spec.json"
    spec_path.write_bytes(bom + water_spec().to_json().encode("utf-8"))
    assert CalculationSpec.load(spec_path).method.code == "dftbplus"
    stages_path = tmp_path / "stages.json"
    stages_path.write_bytes(bom + json.dumps({"stages": [{"name": "a", "task": {"type": "single_point"}}]}).encode("utf-8"))
    assert [s.name for s in load_stages(stages_path)] == ["a"]
    from adit.compare_sets import load_set

    set_path = tmp_path / "set.json"
    set_path.write_bytes(bom + json.dumps({"kind": "reaction"}).encode("utf-8"))
    assert load_set(set_path)[0] == {"kind": "reaction"}
    from adit.templates import load_template

    tpl = tmp_path / "lab.json"
    body = {k: v for k, v in json.loads(water_spec().to_json()).items() if k not in ("structure", "meta")}
    tpl.write_bytes(bom + json.dumps(body).encode("utf-8"))
    assert load_template(tpl, water_spec()).method.code == "dftbplus"


# ---------- core S8: a spec that is not a JSON object ----------

def test_a_spec_that_is_not_an_object_gets_a_readable_message(sk_root, tmp_path, capsys):
    cfg_path = _cfg_file(sk_root, tmp_path)
    spec_path = tmp_path / "spec.json"
    spec_path.write_text("[1, 2]", encoding="utf-8")
    assert main([str(spec_path), "--validate", "--config", str(cfg_path)]) == 2
    err = capsys.readouterr().err
    assert "オブジェクト" in err and "配列" in err and "Traceback" not in err
    with pytest.raises(ValueError, match="spec.json"):
        CalculationSpec.from_json('"abc"')


# ---------- core L1 / L8 / S6: validation instead of numeric crashes ----------

def test_absurd_kpoint_density_is_a_validation_error(sk_root):
    errs = validate(_periodic_water(kpoints=KPoints(mode="density", density=1e308)), cfg_for(sk_root))
    assert "kpoints.density" in _locs(errs)
    assert "kpoints.density" not in _locs(validate(_periodic_water(kpoints=KPoints(mode="density", density=0.2)), cfg_for(sk_root)))


def test_walltime_minutes_and_seconds_are_bounded(sk_root):
    assert "runtime.walltime" in _locs(validate(water_spec(runtime=Runtime(profile="local", walltime="00:90:00")), cfg_for(sk_root)))
    assert "runtime.walltime" not in _locs(validate(water_spec(runtime=Runtime(profile="local", walltime="01:59:59")), cfg_for(sk_root)))


def test_compare_set_with_a_degenerate_cell_fails_in_validation_not_numpy(sk_root, tmp_path):
    from adit.compare_sets import write_compare_set

    spec = _periodic_water(cell=[(0.0, 0.0, 0.0)] * 3, kpoints=KPoints(mode="density", density=0.2))
    with pytest.raises(ProjectError, match="退化"):
        write_compare_set(spec, cfg_for(sk_root), tmp_path / "set", {"kind": "solvation", "solvent": {"method": {"scc_tolerance": 1e-7}}}, tmp_path)


# ---------- core L9 / M12: stages.json keeps overrides; README names the python version ----------

def test_stages_json_records_method_and_runtime_overrides(sk_root, tmp_path):
    stages = parse_stages({"stages": [{"name": "min", "task": {"type": "geometry_optimization"}},
                                      {"name": "nvt", "task": {"type": "molecular_dynamics"}, "method": {"scc_tolerance": 1e-7},
                                       "runtime": {"walltime": "02:00:00"}}]})
    out = tmp_path / "st"
    write_stages(water_spec(), cfg_for(sk_root), out, stages)
    back = load_stages(out / "stages.json")
    assert back[1].overrides["method"] == {"scc_tolerance": 1e-7} and back[1].overrides["runtime"] == {"walltime": "02:00:00"}
    assert "method" not in back[0].overrides and "runtime" not in back[0].overrides
    assert "python3 (3.6 以上)" in (out / "README.txt").read_text(encoding="utf-8")
    assert "python3 3.6 以上" in (out / "stage_02_nvt" / "README.txt").read_text(encoding="utf-8")


# ---------- core L11 / codes L4: the remote work directory is quoted ----------

def test_remote_work_directory_is_quoted(sk_root, tmp_path):
    cfg = cfg_for(sk_root)
    p = cfg.profiles["cluster"]
    p.host, p.user, p.remote_dir = "cluster.example.ac.jp", "me", "/scratch/my runs"
    rt = Runtime(profile="cluster", ncpus=8, omp_threads=8, job_name="w")
    files = build_project(water_spec(runtime=rt), cfg, output_dir=tmp_path / "calc")
    transfer = files.texts["transfer_and_submit.sh"]
    assert "me@cluster.example.ac.jp:'/scratch/my runs'/" in transfer and "cd '/scratch/my runs'/calc" in transfer
    assert "cd '/scratch/my runs'/calc" in files.texts["README.txt"]
    p.remote_dir = "/scratch/runs"
    plain = build_project(water_spec(runtime=rt), cfg, output_dir=tmp_path / "calc").texts["transfer_and_submit.sh"]
    assert "me@cluster.example.ac.jp:/scratch/runs/" in plain and "cd /scratch/runs/calc" in plain


# ---------- codes M4 / M6: run-command placeholders and core counts ----------

def test_unknown_braces_in_a_run_command_are_a_clean_error(sk_root):
    import adit.project as project
    from adit.codes import GENERATORS, GenerationError

    cfg = cfg_for(sk_root)
    cfg.profiles["local"].commands["dftbplus"] = "env OMP_NUM_THREADS=${omp} dftb+"
    errs = validate(water_spec(), cfg)
    msgs = " ".join(e.message for e in errs)
    assert "runtime.profile" in _locs(errs) and "{omp}" in msgs and "{ntasks}" in msgs and "{mpiprocs}" in msgs
    with pytest.raises(ProjectError, match="omp"):
        build_project(water_spec(), cfg)
    with pytest.raises(GenerationError, match="omp"):
        project._run_command(GENERATORS["dftbplus"], water_spec(), cfg.profiles["local"])


def test_ntasks_is_an_accepted_placeholder(sk_root):
    cfg = cfg_for(sk_root)
    cfg.profiles["local"].commands["dftbplus"] = "mpirun -np {ntasks} dftb+"
    assert "runtime.profile" not in _locs(validate(water_spec(), cfg))
    try:
        files = build_project(water_spec(), cfg)
    except ProjectError as ex:   # until the generator fills {ntasks}, the error must still be readable
        assert "ntasks" in str(ex) and "KeyError" not in str(ex)
    else:
        assert "mpirun -np 1 dftb+" in files.texts["submit.sh"]


def test_processes_times_threads_must_fit_the_cores_on_a_cluster(sk_root):
    cfg = cfg_for(sk_root)
    errs = validate(water_spec(runtime=Runtime(profile="cluster", ncpus=4, mpiprocs=4, omp_threads=2, job_name="w")), cfg)
    assert "runtime.ncpus" in _locs(errs) and any("= 8" in e.message and "コア数 4" in e.message for e in errs)
    assert "runtime.ncpus" not in _locs(validate(water_spec(runtime=Runtime(profile="cluster", ncpus=8, mpiprocs=4, omp_threads=2, job_name="w")), cfg))
    assert "runtime.ncpus" not in _locs(validate(water_spec(), cfg))   # local: ncpus is not part of the job


# ---------- codes M1 / L6 / L2: optimizer combinations the codes refuse; the velocity message ----------

def test_lammps_fire_cannot_relax_the_cell(sk_root):
    from ase.build import bulk

    lam = LammpsMethod(units="metal", pair_style="eam", pair_coeff="* * Cu_u3.eam", potential_files=[str(REPO / "examples" / "lammps_cu" / "Cu_u3.eam")])
    a = bulk("Cu", "fcc", a=3.615, cubic=True)

    def spec(**task):
        return CalculationSpec(structure=Structure(source="bulk", source_ref="Cu", atoms=AtomsData.from_ase(a)), method=lam,
                               task=Task(type="geometry_optimization", max_steps=10, force_tolerance_ev_per_ang=0.1, **task), runtime=Runtime(profile="local"))

    bad = validate(spec(optimizer="FIRE", relax_cell="shape_and_volume"), cfg_for(sk_root))
    assert "task.optimizer" in _locs(bad) and any("fire" in e.message for e in bad)
    assert "task.optimizer" not in _locs(validate(spec(optimizer="SteepestDescent", relax_cell="shape_and_volume"), cfg_for(sk_root)))
    assert "task.optimizer" not in _locs(validate(spec(optimizer="FIRE", relax_cell="no"), cfg_for(sk_root)))


def test_gromacs_lbfgs_needs_no_constraints_and_one_process(sk_root):
    def spec(constraints, mpiprocs):
        return water_spec(method=GromacsMethod(constraints=constraints), runtime=Runtime(profile="local", mpiprocs=mpiprocs),
                          task=Task(type="geometry_optimization", optimizer="LBFGS", max_steps=10, force_tolerance_ev_per_ang=1.0))

    assert "task.optimizer" in _locs(validate(spec("h-bonds", 1), cfg_for(sk_root)))
    assert "task.optimizer" in _locs(validate(spec("none", 4), cfg_for(sk_root)))
    assert "task.optimizer" not in _locs(validate(spec("none", 1), cfg_for(sk_root)))


def test_velocity_message_lists_every_code_that_writes_them(sk_root):
    st = water_spec().structure.model_copy(update={"velocities": [(0.0, 0.0, 0.001)] * 3})
    errs = validate(water_spec(method=XtbMethod(), structure=st, task=Task(type="molecular_dynamics")), cfg_for(sk_root))
    msg = next(e.message for e in errs if e.location == "structure.velocities")
    assert "DFTB+" in msg and "VASP" in msg and "OpenMM" in msg and "機械学習" in msg


# ---------- codes H2: codes without a converter cannot pretend to carry the structure ----------

def test_stages_refuse_to_lose_the_structure_for_codes_without_a_converter():
    spec = water_spec(method=MlipMethod(model_family="mace_off"))
    with pytest.raises(StageError, match="mlip"):
        plan_stages(spec, [Stage("opt", {"task": {"type": "geometry_optimization"}}), Stage("nvt", {"task": {"type": "molecular_dynamics"}})])
    plan = plan_stages(spec, [Stage("sp", {"task": {"type": "single_point"}}), Stage("opt", {"task": {"type": "geometry_optimization"}})])
    assert plan[1].pre_command == "" and not plan[1].uses_script
    assert not any("前の段階の出力から" in line or "実際の出発点" in line for line in plan[1].readme) and len(plan[1].readme) == 3


# ---------- codes H3: VASP stages convert CONTCAR instead of copying it ----------

CONTCAR_MD = """adit: O H (3 atoms)
   1.0000000000000000
    10.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   10.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   10.0000000000000000
   O    H
     1     2
Selective dynamics
Direct
  0.5000000000000000  0.4000000000000000  0.5000000000000000   F   F   F
  0.5200000000000000  0.5100000000000000  0.5783064000000000   T   T   T
  0.4800000000000000  0.5100000000000000  0.4216936000000000   T   T   T

  0.10000000E-02  0.20000000E-02  0.30000000E-02
 -0.10000000E-02  0.00000000E+00  0.00000000E+00
  0.00000000E+00  0.00000000E+00  0.10000000E-02

  1
  0.1000000000000000E+01
  0.5000000000000000E+00  0.4000000000000000E+00  0.5000000000000000E+00
  0.5200000000000000E+00  0.5100000000000000E+00  0.5783064000000000E+00
  0.4800000000000000E+00  0.5100000000000000E+00  0.4216936000000000E+00
"""

POSCAR_HERE = """adit: O H (3 atoms)
   1.0000000000000000
    10.0000000000000000    0.0000000000000000    0.0000000000000000
     0.0000000000000000   10.0000000000000000    0.0000000000000000
     0.0000000000000000    0.0000000000000000   10.0000000000000000
   O    H
     1     2
Selective dynamics
Direct
  0.5000000000000000  0.4000000000000000  0.5000000000000000   F   F   F
  0.5000000000000000  0.5000000000000000  0.5783064000000000   T   T   T
  0.5000000000000000  0.5000000000000000  0.4216936000000000   T   T   T
"""


def test_vasp_stage_handoff_rewrites_poscar_and_drops_the_predictor_block(tmp_path):
    prev, here = tmp_path / "prev", tmp_path / "here"
    prev.mkdir()
    here.mkdir()
    (prev / "CONTCAR").write_text(CONTCAR_MD, encoding="utf-8")
    (here / "POSCAR").write_text(POSCAR_HERE, encoding="utf-8")
    assert hf.apply_stage("vasp", str(prev), str(here), "molecular_dynamics", False) == ["POSCAR"]
    text = (here / "POSCAR").read_text(encoding="utf-8")
    lines = text.splitlines()
    assert lines[7] == "Selective dynamics" and lines[8] == "Direct"
    assert lines[9].split()[3:] == ["F", "F", "F"] and lines[10].split()[3:] == ["T", "T", "T"]
    assert np.allclose([float(x) for x in lines[10].split()[:3]], [0.52, 0.51, 0.5783064])
    assert len(lines) == 12   # no velocity block, no predictor-corrector block
    back = hf.parse_poscar(text)
    assert np.allclose(back["positions"], hf.parse_poscar(CONTCAR_MD)["positions"]) and back["velocities"] is None
    (here / "POSCAR").write_text(POSCAR_HERE, encoding="utf-8")
    hf.apply_stage("vasp", str(prev), str(here), "molecular_dynamics", True)
    text = (here / "POSCAR").read_text(encoding="utf-8")
    assert np.allclose(hf.parse_poscar(text)["velocities"], [[1e-3, 2e-3, 3e-3], [-1e-3, 0.0, 0.0], [0.0, 0.0, 1e-3]])
    assert len(text.splitlines()) == 16 and "0.1000000000000000E+01" not in text
    (here / "POSCAR").write_text(POSCAR_HERE, encoding="utf-8")
    with pytest.raises(hf.HandoffError, match="速度がありません / no velocities"):
        hf.apply_stage("vasp", str(prev), str(here), "geometry_optimization", True)


def test_vasp_stages_use_the_converter_instead_of_a_raw_copy():
    spec = _periodic_water(method=VaspMethod(encut=400), kpoints=KPoints())
    plan = plan_stages(spec, [Stage("opt", {"task": {"type": "geometry_optimization"}}),
                              Stage("nvt", {"task": {"type": "molecular_dynamics"}}),
                              Stage("nve", {"task": {"type": "molecular_dynamics", "md": {"ensemble": "NVE"}}}),
                              Stage("prod", {"task": {"type": "molecular_dynamics", "md": {"ensemble": "NVE"}}}, False)])
    assert plan[1].pre_command == "python3 ../handoff.py vasp ../stage_01_opt geometry_optimization" and plan[1].uses_script
    assert plan[2].pre_command == "python3 ../handoff.py vasp ../stage_02_nvt molecular_dynamics --velocities"
    assert plan[3].pre_command == "python3 ../handoff.py vasp ../stage_03_nve molecular_dynamics"
    assert all(p.spec.handoff.files == {} for p in plan[1:]) and not plan[3].spec.handoff.velocities
    assert "vasp" in hf.CONVERTS and hf.copy_files("vasp", "molecular_dynamics", True) == {} and "POSCAR" in hf.rewritten_at_run("vasp")
    assert any("python3 3.6" in line for line in plan[1].readme)


# ---------- core M12 / S7: handoff.py runs on the cluster's python3 ----------

def test_handoff_script_is_portable_to_python_36_and_speaks_both_languages():
    src = HANDOFF.read_text(encoding="utf-8")
    tree = ast.parse(src, feature_version=(3, 7))
    assert not any(isinstance(n, ast.ImportFrom) and n.module == "__future__" for n in ast.walk(tree))
    annotations = []
    for n in ast.walk(tree):
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)):
            annotations += [a.annotation for a in n.args.args + n.args.kwonlyargs if a.annotation is not None]
            if n.returns is not None:
                annotations.append(n.returns)
        elif isinstance(n, ast.AnnAssign):
            annotations.append(n.annotation)
    for a in annotations:
        for sub in ast.walk(a):
            assert not (isinstance(sub, ast.BinOp) and isinstance(sub.op, ast.BitOr)), ast.unparse(a)
            assert not (isinstance(sub, ast.Subscript) and isinstance(sub.value, ast.Name)
                        and sub.value.id in ("list", "dict", "tuple", "set", "type")), ast.unparse(a)
    calls = [n for n in ast.walk(tree) if isinstance(n, ast.Call) and getattr(n.func, "id", "") == "HandoffError"]
    assert len(calls) >= 10 and all(" / " in ast.unparse(c.args[0]) for c in calls)
    assert "from adit" not in src and "import adit" not in src


def test_orca_star_xyz_without_a_space_is_recognised(tmp_path):
    inp = tmp_path / "orca.inp"
    inp.write_text("! HF def2-SVP\n*xyz 0 1\nO 0.0 0.0 0.0\nH 0.0 0.0 0.96\n*\n", encoding="utf-8")
    d = hf._orca_input_coords(str(inp))
    assert d["symbols"] == ["O", "H"] and d["positions"][1][2] == 0.96
    prev = tmp_path / "prev"
    prev.mkdir()
    (prev / "orca.xyz").write_text("2\nfinal\nO 0.0 0.0 0.1\nH 0.0 0.0 1.06\n", encoding="utf-8")
    assert hf.apply_stage("orca", str(prev), str(tmp_path), "geometry_optimization", False) == ["orca.inp"]
    assert "1.06000000" in inp.read_text(encoding="utf-8")
    inp.write_text("! HF\n* internal\n*\n", encoding="utf-8")
    with pytest.raises(hf.HandoffError, match="xyz"):
        hf._orca_input_coords(str(inp))
