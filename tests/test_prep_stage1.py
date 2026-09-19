
from __future__ import annotations

import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import os

import pytest
from ase.build import bulk

from adit import handoff as hf
from adit.continuation import ContinuationError, continue_from
from adit.project import ProjectError, build_project, version_trap, write_project
from adit.spec import (AtomsData, CalculationSpec, Cp2kMethod, DftbMethod, EspressoMethod, Handoff, HubbardU, KPoints,
                        LammpsMethod, MDSettings, OrcaMethod, Runtime, Structure, Task, VaspMethod, XtbMethod)
from adit.stages import Stage, StageError, parse_stages, plan_stages, write_stages
from adit.templates import TemplateError, list_templates, load_template, save_template, template_origin
from tests.conftest import cfg_for, water_spec
from tests.test_cp2k import h2o as cp2k_h2o, make_data
from tests.test_espresso import make_fake_upf, si_spec

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"
REAL_SK = REPO / "slakos"
REAL_PSEUDO = REPO / "pseudo"
need_sk = pytest.mark.skipif(not (REAL_SK / "mio-1-1").is_dir(), reason="slakos/mio-1-1 が無い (リポジトリに入れていない)")


def errors(spec, cfg) -> list:
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg)
    return ex.value.errors


def md(**kw) -> Task:
    return Task(type="molecular_dynamics", md=MDSettings(**{"steps": 10, **kw}))


@pytest.fixture
def cfg(sk_root, tmp_path):
    c = cfg_for(sk_root)
    c.cp2k_data = str(make_data(tmp_path / "cp2k_data"))
    c.pseudo_root = str(make_fake_upf(tmp_path / "pseudo"))
    return c


def xtb_spec(**m):
    return water_spec(method=XtbMethod(**m), task=Task(type="single_point"))


def test_xtb_solvent_on_command_line(cfg):
    run = build_project(xtb_spec(solvation="alpb", solvent="Water"), cfg).texts["submit.sh"].rstrip().splitlines()[-1]
    assert "--alpb water" in run and run.endswith("> output.log 2>&1")
    run = build_project(xtb_spec(gfn="1", solvation="gbsa", solvent="chcl3"), cfg).texts["submit.sh"].rstrip().splitlines()[-1]
    assert "--gbsa chcl3" in run


@pytest.mark.parametrize("m, where", [
    (dict(solvation="alpb", solvent="hexandecane"), "method.solvent"),
    (dict(solvation="gbsa", solvent="dmf", gfn="1"), "method.solvent"),
    (dict(solvation="alpb", solvent="water", gfn="0"), "method.solvation"),
    (dict(solvation="alpb", solvent=""), "method.solvent"),
    (dict(solvation="none", solvent="water"), "method.solvent"),
])
def test_xtb_solvent_checks(cfg, m, where):
    assert any(e.location == where for e in errors(xtb_spec(**m), cfg))


def test_orca_solvent(cfg):
    from adit.codes.orca import solvent_names
    from adit.codes.orca_solvents import SOLVENTS

    assert len(SOLVENTS) == 184 and sum(r[1] for r in SOLVENTS) == 182 and sum(r[2] for r in SOLVENTS) == 179
    spec = lambda **m: water_spec(method=OrcaMethod(**m), task=Task(type="single_point"))  # noqa: E731
    inp = build_project(spec(solvation="smd", solvent="water"), cfg).texts["orca.inp"]
    assert "! HF def2-SVP SMD(water)" in inp
    assert "CPCM(dmso)" in build_project(spec(solvation="cpcm", solvent="dmso"), cfg).texts["orca.inp"]
    assert "diethylether" in solvent_names("smd") and "ammonia" in solvent_names("cpcm") and "ammonia" not in solvent_names("smd")
    for m, frag in ((dict(solvation="smd", solvent="ammonia"), "SMD"), (dict(solvation="cpcm", solvent="unobtainium"), "表にありません"),
                    (dict(solvation="cpcm", solvent="diethyl ether"), "diethylether")):
        errs = errors(spec(**m), cfg)
        assert any(e.location == "method.solvent" and frag in e.message for e in errs), errs


def test_dftb_solvation_param_file(cfg, tmp_path):
    p = tmp_path / "param_gbsa_h2o.txt"
    p.write_text("80.2\n", encoding="utf-8")
    files = build_project(water_spec(method=DftbMethod(sk_set="fake-1-0", solvation_param_file=str(p))), cfg)
    assert 'Solvation = GeneralisedBorn {\n    ParamFile = "param_gbsa_h2o.txt"\n  }' in files.texts["dftb_in.hsd"]
    assert files.copies["param_gbsa_h2o.txt"] == p
    assert any(e.location == "method.solvation_param_file" for e in errors(water_spec(method=DftbMethod(sk_set="fake-1-0", solvation_param_file=str(tmp_path / "no.txt"))), cfg))


def test_cp2k_sccs(cfg):
    inp = build_project(cp2k_h2o(m={"sccs_relative_permittivity": 78.36}), cfg).texts["cp2k.inp"]
    assert "    &SCCS\n      RELATIVE_PERMITTIVITY 78.36\n    &END SCCS" in inp
    assert "&SCCS" not in build_project(cp2k_h2o(), cfg).texts["cp2k.inp"]


def feo(method, **kw) -> CalculationSpec:
    a = bulk("FeO", "rocksalt", a=4.3)
    return CalculationSpec(structure=Structure(source="bulk", source_ref="FeO", atoms=AtomsData.from_ase(a)), method=method,
                           kpoints=KPoints(mode="mesh", mesh=(2, 2, 2)), task=kw.pop("task", Task(type="single_point")),
                           runtime=Runtime(profile="local"), **kw)


def test_vasp_magmom_and_ldau(cfg):
    m = VaspMethod(ispin=2, magmom_by_element={"Fe": 4.0}, hubbard={"Fe": HubbardU(orbital="3d", u_ev=5.3)})
    incar = build_project(feo(m), cfg).texts["INCAR"]
    for line in ("MAGMOM = 4 0", "LDAU = .TRUE.", "LDAUTYPE = 2", "LDAUL = 2 -1", "LDAUU = 5.3 0", "LDAUJ = 0 0"):
        assert line in incar, line
    assert "LMAXMIX" not in incar
    bands = build_project(feo(m, task=Task(type="band_structure")), cfg).texts
    assert "LMAXMIX = 4" in bands["INCAR"] and "LMAXMIX = 4" in bands["bands/INCAR"]


def test_vasp_per_atom_magmom_follows_poscar_order(cfg):
    a = bulk("FeO", "rocksalt", a=4.3).repeat((2, 1, 1))[[1, 0, 3, 2]]
    spec = CalculationSpec(structure=Structure(source="bulk", source_ref="x", atoms=AtomsData.from_ase(a)), method=VaspMethod(ispin=2, magmom=[0.1, 4.0, 0.2, -4.0]),
                           kpoints=KPoints(), runtime=Runtime(profile="local"))
    assert "MAGMOM = 0.1 0.2 4 -4" in build_project(spec, cfg).texts["INCAR"]


@pytest.mark.parametrize("m, where", [
    (VaspMethod(ispin=1, magmom_by_element={"Fe": 4.0}), "method.ispin"),
    (VaspMethod(ispin=2, magmom_by_element={"Ni": 1.0}), "method.magmom_by_element"),
    (VaspMethod(ispin=2, magmom=[1, 1], magmom_by_element={"Fe": 1.0}), "method.magmom_by_element"),
    (VaspMethod(hubbard={"Ni": HubbardU(orbital="3d", u_ev=6.2)}), "method.hubbard"),
    (VaspMethod(hubbard={"Fe": HubbardU(orbital="d", u_ev=5.3)}), "method.hubbard"),
])
def test_vasp_magnetism_checks(cfg, m, where):
    assert any(e.location == where for e in errors(feo(m), cfg))


def test_qe_starting_magnetization_and_hubbard(cfg):
    m = EspressoMethod(pseudo_set="fake-pbe", pseudo={"Si": "Si.pbe-fake.UPF"}, ecutwfc=40, nspin=2, starting_magnetization={"Si": 0.5},
                       hubbard={"Si": HubbardU(orbital="3p", u_ev=1.5)}, hubbard_projector="ortho-atomic")
    pw = build_project(si_spec(method=m), cfg).texts["pw.in"]
    assert "starting_magnetization(1) = 0.5" in pw
    assert pw.rstrip().endswith("HUBBARD {ortho-atomic}\nU Si-3p 1.5")
    for bad, where in ((m.model_copy(update={"nspin": 1, "hubbard": {}}), "method.nspin"),
                       (m.model_copy(update={"hubbard": {"Si": HubbardU(orbital="3p", u_ev=1, j_ev=0.5)}}), "method.hubbard"),
                       (m.model_copy(update={"starting_magnetization": {"O": 0.1}}), "method.starting_magnetization")):
        assert any(e.location == where for e in errors(si_spec(method=bad), cfg))


@pytest.mark.skipif(not (REAL_PSEUDO / "pslibrary").is_dir(), reason="pseudo/pslibrary が無い")
def test_qe_hubbard_shell_must_be_in_upf(sk_root):
    c = cfg_for(sk_root)
    c.pseudo_root = str(REAL_PSEUDO)
    m = EspressoMethod(pseudo_set="pslibrary", ecutwfc=30, hubbard={"Si": HubbardU(orbital="3d", u_ev=1.0)})
    errs = errors(si_spec(pseudo_set="pslibrary", method=m), c)
    assert any(e.location == "method.hubbard" and "PP_CHI" in e.message and "3P" in e.message for e in errs), errs
    build_project(si_spec(pseudo_set="pslibrary", method=m.model_copy(update={"hubbard": {"Si": HubbardU(orbital="3p", u_ev=1.0)}})), c)


def test_cp2k_magnetization_and_plus_u(cfg):
    m = {"uks": True, "magnetization_by_element": {"O": 2.0}, "hubbard": {"O": HubbardU(orbital="2p", u_ev=3.0, j_ev=0.5)}}
    inp = build_project(cp2k_h2o(m=m, structure=water_spec().structure.model_copy(update={"multiplicity": 3})), cfg).texts["cp2k.inp"]
    assert "    PLUS_U_METHOD MULLIKEN" in inp
    assert "    &KIND O\n      BASIS_SET DZVP-MOLOPT-SR-GTH\n      POTENTIAL GTH-PBE-q6\n      MAGNETIZATION 2\n      &DFT_PLUS_U\n        L 1\n        U_MINUS_J [eV] 2.5\n      &END DFT_PLUS_U\n    &END KIND" in inp
    assert any(e.location == "method.uks" for e in errors(cp2k_h2o(m={"magnetization_by_element": {"O": 2.0}}), cfg))


def test_continue_dftb_md_carries_positions_and_velocities(sk_root):
    d = EX / "dftb_md_water_generated"
    spec = continue_from(d)
    rows = hf.last_xyz_frame(str(d / "geo_end.xyz"))
    assert np.allclose(spec.structure.atoms.positions, [[float(x) for x in r[1:4]] for r in rows])
    assert np.allclose(spec.structure.velocities, [[float(x) / 1000 for x in r[-3:]] for r in rows])  # Å/ps → Å/fs
    assert spec.handoff.velocities and spec.meta.continued_from["dir"] == str(d) and spec.meta.continued_from["task"] == "molecular_dynamics"
    assert spec.meta.continued_from["spec_sha256"] == hashlib.sha256((d / "spec.json").read_bytes()).hexdigest()
    s2 = spec.model_copy(update={"method": DftbMethod(sk_set="fake-1-0")})
    files = build_project(s2, cfg_for(sk_root))
    hsd = files.texts["dftb_in.hsd"]
    assert "  Velocities [AA/ps] {\n    " + f"{float(rows[0][-3]):.10f}"[:8] in hsd
    assert "== 前の計算から引き継いだもの ==" in files.texts["README.txt"] and "geo_end.xyz" in files.texts["README.txt"]
    assert CalculationSpec.from_json(files.texts["spec.json"]) == s2
    nov = continue_from(d, velocities=False)
    assert nov.structure.velocities is None and any("引き継がない" in x for x in nov.meta.continued_from["not_carried"])


def test_continue_from_optimizations():
    s = continue_from(EX / "water_generated")
    gen = hf.parse_gen((EX / "water_generated" / "geom.out.gen").read_text(encoding="utf-8"))
    assert np.allclose(s.structure.atoms.positions, gen["positions"]) and s.structure.velocities is None
    x = continue_from(EX / "xtb_water_generated")
    assert np.allclose(x.structure.atoms.positions, [[float(v) for v in r[1:4]] for r in hf.last_xyz_frame(str(EX / "xtb_water_generated" / "xtbopt.xyz"))])
    v = continue_from(EX / "vasp_h2o_generated")
    c = hf.parse_poscar((EX / "vasp_h2o_generated" / "CONTCAR").read_text(encoding="utf-8"))
    assert v.structure.atoms.symbols == ["O", "H", "H"] and np.allclose(v.structure.atoms.positions[1], c["positions"][1])


def test_continue_with_new_conditions_and_across_codes():
    d = EX / "dftb_md_water_generated"
    prev = CalculationSpec.load(d / "spec.json")
    cond = prev.model_copy(update={"task": md(ensemble="NVE", temperature_k=250)})
    s = continue_from(d, cond)
    assert s.task.md.ensemble == "NVE" and s.task.md.temperature_k == 250 and s.structure.velocities is not None
    other = continue_from(d, water_spec(method=XtbMethod(), task=md(thermostat="berendsen")))
    assert other.structure.velocities is None and not other.handoff.velocities and any("xtb" in x for x in other.meta.continued_from["not_carried"])


def test_continue_qe_md_positions_only(tmp_path):
    s = continue_from(EX / "qe_md_si_generated")
    assert s.structure.velocities is None and any("速度" in x for x in s.meta.continued_from["not_carried"])
    assert not np.allclose(s.structure.atoms.positions, CalculationSpec.load(EX / "qe_md_si_generated" / "spec.json").structure.atoms.positions)


def test_continue_cp2k_lammps_gromacs(cfg, monkeypatch):
    monkeypatch.chdir(REPO)
    c = continue_from(EX / "cp2k_h2o_md_generated")
    from adit.handoff import existing_name

    restart = existing_name(str(EX / "cp2k_h2o_md_generated"), "{}-1.restart")
    assert restart.endswith("-1.restart")
    assert c.handoff.files == {"prev.restart": restart} and c.handoff.velocities
    inp = build_project(c.model_copy(update={"method": c.method.model_copy(update={"basis": {"O": "DZVP-MOLOPT-SR-GTH", "H": "DZVP-MOLOPT-SR-GTH"}})}), cfg)
    assert "&EXT_RESTART\n  RESTART_FILE_NAME prev.restart\n  RESTART_DEFAULT .FALSE.\n  RESTART_VEL .TRUE.\n&END EXT_RESTART" in inp.texts["cp2k.inp"]
    assert inp.copies["prev.restart"] == EX / "cp2k_h2o_md_generated" / restart
    lm = continue_from(EX / "lammps_cu_nvt_generated")
    assert lm.method.data_file.endswith("final.data") and lm.method.type_elements == ["Cu"] and lm.handoff.velocities
    text = build_project(lm, cfg).texts["in.lammps"]
    assert "velocity all create" not in text and "read_data data.lammps" in text
    gm = continue_from(EX / "gromacs_spce_nvt_generated")
    assert gm.method.checkpoint_file.endswith(".cpt") and gm.method.structure_file.endswith(".gro")
    f = build_project(gm, cfg)
    assert "continuation" in f.texts["grompp.mdp"] and "-t prev.cpt" in f.texts["submit.sh"]


def test_continue_needs_the_files(tmp_path):
    with pytest.raises(ContinuationError, match="spec.json"):
        continue_from(tmp_path)
    d = tmp_path / "half"
    d.mkdir()
    shutil.copy(EX / "dftb_md_water_generated" / "spec.json", d)
    shutil.copy(EX / "dftb_md_water_generated" / "geometry.gen", d)
    with pytest.raises(ContinuationError, match="geo_end.xyz"):
        continue_from(d)


def test_velocity_and_handoff_checks(cfg, tmp_path):
    st = water_spec().structure.model_copy(update={"velocities": [(0.0, 0.0, 0.0)] * 2})
    assert any(e.location == "structure.velocities" for e in errors(water_spec(structure=st, task=md(thermostat="berendsen")), cfg))
    st3 = water_spec().structure.model_copy(update={"velocities": [(0.0, 0.0, 0.0)] * 3})
    assert any(e.location == "structure.velocities" for e in errors(water_spec(structure=st3), cfg))
    assert any(e.location == "structure.velocities" for e in errors(water_spec(structure=st3, method=XtbMethod(), task=md()), cfg))
    h = Handoff(previous_dir=str(tmp_path), previous_task="molecular_dynamics", previous_code="xtb", velocities=True, files={"mdrestart": "mdrestart"})
    assert any(e.location == "handoff" and "mdrestart" in e.message for e in errors(water_spec(method=XtbMethod(), task=md(), handoff=h), cfg))


def test_vasp_poscar_velocities_roundtrip(cfg):
    a = bulk("FeO", "rocksalt", a=4.3).repeat((2, 1, 1))[[1, 0, 3, 2]]
    vel = [(0.01 * i, -0.02 * i, 0.003) for i in range(4)]
    spec = CalculationSpec(structure=Structure(source="bulk", source_ref="x", atoms=AtomsData.from_ase(a), velocities=vel), method=VaspMethod(),
                           kpoints=KPoints(), task=md(thermostat="andersen"), runtime=Runtime(profile="local"))
    back = hf.parse_poscar(build_project(spec, cfg).texts["POSCAR"])
    order = [0, 2, 1, 3]
    assert back["symbols"] == ["O", "O", "Fe", "Fe"] and np.allclose(back["velocities"], [vel[i] for i in order])


STAGES = {"stages": [{"name": "min", "task": {"type": "geometry_optimization", "max_steps": 20}},
                     {"name": "nvt", "task": {"type": "molecular_dynamics", "md": {"ensemble": "NVT", "thermostat": "berendsen", "steps": 20}}},
                     {"name": "nve", "task": {"type": "molecular_dynamics", "md": {"ensemble": "NVE", "steps": 20}}}]}


def test_stages_plan_and_files(cfg, tmp_path):
    out = tmp_path / "staged"
    dirs = write_stages(water_spec(), cfg, out, parse_stages(STAGES))
    assert [d.name for d in dirs] == ["stage_01_min", "stage_02_nvt", "stage_03_nve"]
    sub2 = (out / "stage_02_nvt" / "submit.sh").read_text(encoding="utf-8").rstrip().splitlines()[-1]
    assert sub2.startswith("python3 ../handoff.py dftbplus ../stage_01_min geometry_optimization && dftb+")
    sub3 = (out / "stage_03_nve" / "submit.sh").read_text(encoding="utf-8").rstrip().splitlines()[-1]
    assert sub3.startswith("python3 ../handoff.py dftbplus ../stage_02_nvt molecular_dynamics --velocities && ")
    hsd3 = (out / "stage_03_nve" / "dftb_in.hsd").read_text(encoding="utf-8")
    assert '<<< "velocities.dat"' in hsd3
    assert "Thermostat = None {}" in hsd3 and "InitialTemperature" not in hsd3
    assert "Velocities" not in (out / "stage_02_nvt" / "dftb_in.hsd").read_text(encoding="utf-8")
    top = (out / "submit.sh").read_text(encoding="utf-8")
    assert "for d in stage_01_min stage_02_nvt stage_03_nve; do" in top and "set -e" in top
    readme = (out / "README.txt").read_text(encoding="utf-8")
    assert "j2=$(cd stage_02_nvt && qsub -W depend=afterok:$j1 submit.sh)" in readme
    assert "j3=$(cd stage_03_nve && sbatch --parsable --dependency=afterok:$j2 submit.sh)" in readme
    assert (out / "handoff.py").read_text(encoding="utf-8") == Path(hf.__file__).read_text(encoding="utf-8")
    s3 = CalculationSpec.load(out / "stage_03_nve" / "spec.json")
    assert s3.handoff.at_run and s3.handoff.previous_dir == "../stage_02_nvt" and s3.meta.stage["index"] == 3
    assert json.loads((out / "stages.json").read_text(encoding="utf-8"))["stages"][2]["velocities"] is True


def test_stages_copy_codes(cfg, tmp_path, monkeypatch):
    lam = LammpsMethod(units="metal", pair_style="eam", pair_coeff="* * Cu_u3.eam", potential_files=[str(EX / "lammps_cu" / "Cu_u3.eam")])
    a = bulk("Cu", "fcc", a=3.615, cubic=True)
    spec = CalculationSpec(structure=Structure(source="bulk", source_ref="Cu", atoms=AtomsData.from_ase(a)), method=lam, runtime=Runtime(profile="local"))
    plan = plan_stages(spec, parse_stages({"stages": [{"name": "min", "task": {"type": "geometry_optimization", "max_steps": 10}},
                                                      {"name": "nvt", "task": {"type": "molecular_dynamics", "md": {"thermostat": "nose_hoover"}}},
                                                      {"name": "prod", "task": {"type": "molecular_dynamics", "md": {"thermostat": "nose_hoover"}}}]}))
    assert plan[1].pre_command == "cp ../stage_01_min/final.data data.lammps" and not plan[1].uses_script
    files = build_project(plan[2].spec, cfg, pre_command=plan[2].pre_command)
    assert "velocity all create" not in files.texts["in.lammps"]
    assert "velocity all create" in build_project(plan[1].spec, cfg).texts["in.lammps"]


def test_stage_errors(cfg):
    with pytest.raises(StageError, match="計算コード"):
        plan_stages(water_spec(), parse_stages({"stages": [{"task": {}}, {"method": {"code": "xtb"}}]}))
    with pytest.raises(StageError, match="MD"):
        plan_stages(water_spec(), [Stage("a", {"task": {"type": "geometry_optimization"}}), Stage("b", {"task": {"type": "molecular_dynamics"}}, True)])
    with pytest.raises(StageError, match="速度"):
        plan_stages(si_spec(), [Stage("a", {"task": {"type": "molecular_dynamics"}}), Stage("b", {"task": {"type": "molecular_dynamics"}}, True)])
    with pytest.raises(StageError, match="使えない項目"):
        parse_stages({"stages": [{"structure": {}}]})


@pytest.mark.skipif(os.name == "nt", reason="handoff.py を素の PATH で実行する試験 (走り先は Linux)")
def test_handoff_script_on_real_outputs(tmp_path):
    here = tmp_path / "here"
    here.mkdir()
    shutil.copy(EX / "dftb_md_water_generated" / "geometry.gen", here)
    script = Path(hf.__file__)
    r = subprocess.run(["python3", str(script), "dftbplus", str(EX / "dftb_md_water_generated"), "molecular_dynamics", "--velocities"],
                       cwd=here, capture_output=True, text=True, env={"PATH": "/usr/bin:/bin"})
    assert r.returncode == 0, r.stderr
    rows = hf.last_xyz_frame(str(EX / "dftb_md_water_generated" / "geo_end.xyz"))
    assert np.allclose(hf.parse_gen((here / "geometry.gen").read_text(encoding="utf-8"))["positions"], [[float(x) for x in r_[1:4]] for r_ in rows])
    assert np.allclose(np.loadtxt(here / "velocities.dat"), [[float(x) for x in r_[-3:]] for r_ in rows])
    qe = tmp_path / "qe"
    qe.mkdir()
    shutil.copy(EX / "qe_md_si_generated" / "pw.in", qe)
    hf.apply_stage("espresso", str(EX / "qe_md_si_generated"), str(qe), "molecular_dynamics", False)
    pw = (qe / "pw.in").read_text(encoding="utf-8")
    fin = hf.final_espresso(str(EX / "qe_md_si_generated"))
    assert "ATOMIC_POSITIONS angstrom" in pw and f"{fin['positions'][1][0]:.10f}" in pw
    r = subprocess.run(["python3", str(script), "dftbplus", str(tmp_path / "nothing"), "molecular_dynamics"], cwd=here, capture_output=True, text=True)
    assert r.returncode == 1 and "handoff.py:" in r.stderr


def test_templates_save_list_load(cfg, tmp_path, monkeypatch):
    monkeypatch.setenv("ADIT_CONFIG", str(tmp_path / "conf" / "cluster.toml"))
    cfg.templates_dir = str(tmp_path / "lab")
    base = water_spec(method=DftbMethod(sk_set="fake-1-0", scc_tolerance=1e-7), task=md(thermostat="nose_hoover"))
    p = save_template(base, "lab-dftb-md", cfg=cfg, comment="研究室の DFTB+ MD")
    data = json.loads(p.read_text(encoding="utf-8"))
    assert p.parent == tmp_path / "lab" and "structure" not in data and data["template"]["comment"] == "研究室の DFTB+ MD"
    with pytest.raises(TemplateError):
        save_template(base, "lab-dftb-md", cfg=cfg)
    (tmp_path / "conf" / "templates").mkdir(parents=True)
    (tmp_path / "conf" / "templates" / "broken.json").write_text("{}", encoding="utf-8")
    infos = {i.name: i for i in list_templates(cfg)}
    assert infos["lab-dftb-md"].code == "dftbplus" and infos["broken"].error
    other = Structure(source="smiles", source_ref="C", atoms=AtomsData(symbols=["C", "H", "H", "H", "H"],
                                                                         positions=[(0, 0, 0), (0.63, 0.63, 0.63), (-0.63, -0.63, 0.63), (-0.63, 0.63, -0.63), (0.63, -0.63, -0.63)]))
    s = load_template("lab-dftb-md", other, cfg)
    assert s.structure == other and s.method.scc_tolerance == 1e-7 and s.meta.template["name"] == "lab-dftb-md"
    origin = template_origin(s)
    assert origin["method.scc_tolerance"] == "lab-dftb-md" and origin["task.md.thermostat"] == "lab-dftb-md"
    changed = s.model_copy(update={"method": s.method.model_copy(update={"scc_tolerance": 1e-6})})
    assert "method.scc_tolerance" not in template_origin(changed) and "method.sk_set" in template_origin(changed)
    with pytest.raises(TemplateError):
        load_template("broken", other, cfg)


def test_provenance_in_spec_and_readme(cfg, tmp_path):
    out = tmp_path / "w"
    write_project(water_spec(), cfg, out)
    prov = json.loads((out / "spec.json").read_text(encoding="utf-8"))["provenance"]
    assert prov["ase"] == "3.29.0" and prov["python"] and prov["adit_version"] and prov["code_version_file"] == "code_version.txt"
    names = {f["name"]: f["sha256"] for f in prov["files"]}
    assert names["skf/H-O.skf"] == hashlib.sha256((out / "skf" / "H-O.skf").read_bytes()).hexdigest()
    readme = (out / "README.txt").read_text(encoding="utf-8")
    assert "== 作成時の記録 ==" in readme and names["skf/H-O.skf"] in readme
    sub = (out / "submit.sh").read_text(encoding="utf-8")
    assert "trap 'grep -m1 -E" in sub and sub.rstrip().splitlines()[-1] == "dftb+ > output.log 2>&1"
    cp = build_project(cp2k_h2o(), cfg)
    gen = json.loads(cp.texts["spec.json"])["provenance"]["generated_files"]
    assert {g["name"] for g in gen} == {"BASIS_adit", "POTENTIAL_adit"} and all(len(g["source_sha256"]) == 64 for g in gen)
    orca = json.loads(build_project(water_spec(method=OrcaMethod(), task=Task()), cfg).texts["spec.json"])["provenance"]
    assert orca["code_version_file"] is None


@pytest.mark.skipif(os.name == "nt", reason="sh で実行する試験 (走り先は Linux)")
def test_version_trap_keeps_exit_status(tmp_path):
    body = version_trap(("output.log", "xtb version")) + "\necho '   * xtb version 6.7.1 (edcfbbe)' > output.log\nfalse\n"
    (tmp_path / "s.sh").write_text("#!/bin/sh\n" + body, encoding="utf-8")
    r = subprocess.run(["sh", "s.sh"], cwd=tmp_path)
    assert r.returncode == 1 and (tmp_path / "code_version.txt").read_text(encoding="utf-8").strip() == "* xtb version 6.7.1 (edcfbbe)"


def test_espresso_version_probe_matches_pw_and_neb():
    from adit.codes.espresso import EspressoGenerator

    name, pattern = EspressoGenerator().version_probe(None)
    assert name == "output.log"
    assert re.search(pattern, "Program PWSCF v.7.5 starts on 1Jan2026")
    assert re.search(pattern, "Program NEB v.7.6 starts on 2Jan2026")


def test_cli_continue_stages_templates(sk_root, tmp_path, monkeypatch, capsys):
    from adit.cli import main
    from adit.config import save_config

    conf = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), conf)
    monkeypatch.setenv("ADIT_CONFIG", str(conf))
    spec = tmp_path / "spec.json"
    water_spec(task=md(thermostat="berendsen")).save(spec)
    prev = tmp_path / "prev"
    shutil.copytree(EX / "dftb_md_water_generated", prev)
    s = CalculationSpec.load(prev / "spec.json")
    s.model_copy(update={"method": DftbMethod(sk_set="fake-1-0")}).save(prev / "spec.json")
    assert main(["--continue-from", str(prev), str(tmp_path / "cont"), "--set", "task.md.steps=7"]) == 0
    c = CalculationSpec.load(tmp_path / "cont" / "spec.json")
    assert c.task.md.steps == 7 and c.structure.velocities is not None and "Velocities [AA/ps]" in (tmp_path / "cont" / "dftb_in.hsd").read_text(encoding="utf-8")
    st = tmp_path / "stages.json"
    st.write_text(json.dumps(STAGES), encoding="utf-8")
    assert main([str(spec), str(tmp_path / "staged"), "--stages", str(st)]) == 0
    assert (tmp_path / "staged" / "stage_03_nve" / "submit.sh").is_file()
    assert main([str(spec), "--save-template", "t1"]) == 0 and (tmp_path / "templates" / "t1.json").is_file()
    assert main(["--list-templates"]) == 0 and "t1" in capsys.readouterr().out
    assert main([str(spec), str(tmp_path / "fromtpl"), "--template", "t1"]) == 0
    assert CalculationSpec.load(tmp_path / "fromtpl" / "spec.json").meta.template["name"] == "t1"
    assert main(["--continue-from", str(tmp_path / "nothing"), str(tmp_path / "x")]) == 1


def test_old_spec_json_still_loads():
    for d in sorted(EX.glob("*_generated")):
        s = CalculationSpec.load(d / "spec.json")
        assert s.handoff is None and s.structure.velocities is None and s.meta.continued_from is None


def test_write_stages_does_not_create_missing_parents(cfg, tmp_path):
    from adit.project import ProjectError
    with pytest.raises(ProjectError) as ex:
        write_stages(water_spec(), cfg, tmp_path / "no" / "such" / "staged", parse_stages(STAGES))
    assert [e.location for e in ex.value.errors] == ["output_dir"] and not (tmp_path / "no").exists()
