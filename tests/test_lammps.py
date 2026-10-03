
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from ase.build import bulk

from adit.config import save_config
from adit.project import ProjectError, build_project, write_project
from adit.spec import AtomsData, CalculationSpec, KPoints, LammpsMethod, MDSettings, Runtime, Structure, Task
from tests.conftest import cfg_for, water_spec

REPO = Path(__file__).resolve().parent.parent
REAL_EAM = REPO / "examples" / "lammps_cu" / "Cu_u3.eam"


@pytest.fixture
def pot(tmp_path) -> Path:
    p = tmp_path / "pots" / "Cu_u3.eam"
    p.parent.mkdir()
    p.write_text("fake eam\n", encoding="utf-8")
    return p


def cu(pot, *, cell=None, **kw) -> CalculationSpec:
    atoms = bulk("Cu", "fcc", a=3.615, cubic=True)
    m = dict(units="metal", pair_style="eam", pair_coeff="* * Cu_u3.eam", potential_files=[str(pot)])
    m.update(kw.pop("m", {}))
    md = MDSettings(ensemble="NVT", thermostat="nose_hoover", temperature_k=300, timestep_fs=1.0, steps=100, dump_interval=10, coupling_time_fs=100)
    return CalculationSpec(structure=Structure(source="bulk", source_ref="Cu", atoms=AtomsData.from_ase(atoms)), method=LammpsMethod(**m),
                           task=kw.pop("task", Task(type="molecular_dynamics", md=md)),
                           runtime=Runtime(profile="local", mpiprocs=1, omp_threads=1, job_name="cu"), **kw)


def errors(spec, cfg) -> list:
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg)
    return ex.value.errors


def test_cu_nvt(sk_root, pot):
    files = build_project(cu(pot), cfg_for(sk_root))
    inp = files.texts["in.lammps"]
    for key in ["units metal", "atom_style atomic", "boundary p p p", "pair_style eam", "read_data data.lammps", "pair_coeff * * Cu_u3.eam",
                "velocity all create 300 12345 mom yes rot no dist gaussian", "timestep 0.001",
                "dump traj all custom 10 traj.lammpstrj id type element x y z", "dump_modify traj element Cu sort id",
                "fix integ all nvt temp 300 300 0.1", "run 100", "write_data final.data", "# 原子の型番号と元素: 1=Cu"]:
        assert key in inp, key
    assert inp.index("pair_style") < inp.index("read_data") < inp.index("pair_coeff")
    data = files.texts["data.lammps"]
    assert data.startswith("adit: atom types 1=Cu") and "4 atoms" in data and "1 atom types" in data and "Masses" in data
    assert files.copies == {"Cu_u3.eam": pot}
    assert files.texts["submit.sh"].rstrip().endswith("mpirun -np 1 lmp -in in.lammps -log log.lammps > output.log 2>&1")
    readme = files.texts["README.txt"]
    assert "LAMMPS" in readme and "1 = Cu" in readme and "conda install -c conda-forge lammps" in readme
    build_project(cu(pot, kpoints=KPoints(mode="mesh", mesh=(2, 2, 2))), cfg_for(sk_root))


@pytest.mark.parametrize("thermostat,lines", [
    ("berendsen", ["fix integ all nve", "fix tstat all temp/berendsen 300 300 0.1"]),
    ("csvr", ["fix integ all nve", "fix tstat all temp/csvr 300 300 0.1 12345"]),
    ("langevin", ["fix integ all nve", "fix tstat all langevin 300 300 0.1 12345"]),
])
def test_thermostats(sk_root, pot, thermostat, lines):
    md = MDSettings(ensemble="NVT", thermostat=thermostat, temperature_k=300, coupling_time_fs=100)
    inp = build_project(cu(pot, task=Task(type="molecular_dynamics", md=md)), cfg_for(sk_root)).texts["in.lammps"]
    for l in lines:
        assert l in inp, l


def test_npt_units_and_minimize(sk_root, pot):
    md = MDSettings(ensemble="NPT", thermostat="nose_hoover", temperature_k=300, timestep_fs=2.0, coupling_time_fs=100,
                    pressure_bar=1.01325, barostat_time_fs=1000)
    inp = build_project(cu(pot, m={"units": "real"}, task=Task(type="molecular_dynamics", md=md)), cfg_for(sk_root)).texts["in.lammps"]
    assert "timestep 2" in inp and "fix integ all npt temp 300 300 100 iso 1 1 1000" in inp
    md.thermostat = "berendsen"
    inp = build_project(cu(pot, task=Task(type="molecular_dynamics", md=md)), cfg_for(sk_root)).texts["in.lammps"]
    assert "fix pstat all press/berendsen iso 1.01325 1.01325 1" in inp
    opt = Task(type="geometry_optimization", optimizer="FIRE", max_steps=50, force_tolerance_ev_per_ang=0.1, relax_cell="shape_and_volume")
    with pytest.raises(ProjectError, match="min_style fire"):  # LAMMPS: fire does not support fix box/relax
        build_project(cu(pot, m={"units": "real"}, task=opt), cfg_for(sk_root))
    opt.relax_cell = "no"
    inp = build_project(cu(pot, m={"units": "real"}, task=opt), cfg_for(sk_root)).texts["in.lammps"]
    assert "box/relax" not in inp and "min_style fire" in inp and "minimize 0.0 2.306054783 50 500" in inp
    opt.relax_cell = "shape_and_volume"; opt.optimizer = "Rational"
    inp = build_project(cu(pot, m={"units": "real"}, task=opt), cfg_for(sk_root)).texts["in.lammps"]
    assert "fix relax all box/relax aniso 0.0" in inp and "min_style cg" in inp
    opt.relax_cell = "volume_only"; opt.optimizer = "Rational"
    inp = build_project(cu(pot, task=opt), cfg_for(sk_root)).texts["in.lammps"]
    assert "fix relax all box/relax iso 0.0" in inp and "min_style cg" in inp and "minimize 0.0 0.1 50 500" in inp
    sp = build_project(cu(pot, task=Task(type="single_point")), cfg_for(sk_root)).texts["in.lammps"]
    assert "run 0" in sp and "id type element x y z fx fy fz" in sp and "velocity" not in sp


def test_molecule_box_fixed_atoms_and_type_order(sk_root, pot):
    spec = water_spec(method=LammpsMethod(units="real", atom_style="charge", pair_style="lj/cut/coul/cut 10.0", pair_coeff="* * 0.1 3.0",
                                          type_elements=["H", "O"]),
                      task=Task(type="geometry_optimization", max_steps=10))
    spec = spec.model_copy(update={"structure": spec.structure.model_copy(update={"fixed_atoms": [0]})})
    files = build_project(spec, cfg_for(sk_root))
    inp, data = files.texts["in.lammps"], files.texts["data.lammps"]
    assert "boundary s s s" in inp and "group fixed id 1" in inp and "fix freeze fixed setforce 0.0 0.0 0.0" in inp
    assert "dump_modify traj element H O sort id" in inp and data.startswith("adit: atom types 1=H 2=O")
    assert "Atoms # charge" in data and "2 atom types" in data
    lo = [float(l.split()[0]) for l in data.splitlines() if l.endswith(("xlo xhi", "ylo yhi", "zlo zhi"))]
    assert lo == [0.0, 0.0, 0.0]


def test_external_data_file(sk_root, tmp_path):
    d = tmp_path / "water.data"
    d.write_text("LAMMPS data made elsewhere\n\n3 atoms\n2 bonds\n1 angles\n2 atom types\n1 bond types\n1 angle types\n\n"
                 "0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\nMasses\n\n1 15.999\n2 1.008\n\nAtoms # full\n\n"
                 "1 1 1 -0.8 5 5 5\n2 1 2 0.4 5.9 5 5\n3 1 2 0.4 4.7 5.9 5\n", encoding="utf-8")
    m = dict(units="real", atom_style="full", data_file=str(d), type_elements=["O", "H"], pair_style="lj/cut/coul/cut 10.0",
             pair_coeff="1 1 0.155 3.166\n2 2 0.0 0.0", style_commands="bond_style harmonic\nangle_style harmonic")
    spec = water_spec(method=LammpsMethod(**m), task=Task(type="single_point"))
    files = build_project(spec, cfg_for(sk_root))
    assert "data.lammps" not in files.texts and files.copies == {"data.lammps": d}
    inp = files.texts["in.lammps"]
    assert inp.index("bond_style harmonic") < inp.index("read_data") and "pair_coeff 2 2 0.0 0.0" in inp
    assert any(e.location == "method.type_elements" for e in errors(water_spec(method=LammpsMethod(**{**m, "type_elements": ["O"]})), cfg_for(sk_root)))
    assert any(e.location == "method.atom_style" for e in errors(water_spec(method=LammpsMethod(**{**m, "atom_style": "molecular"})), cfg_for(sk_root)))
    four = water_spec(structure=Structure(source="file", source_ref="x", atoms=AtomsData(symbols=["O", "H", "H", "H"],
                      positions=[(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)])), method=LammpsMethod(**m))
    assert any(e.location == "structure.atoms" for e in errors(four, cfg_for(sk_root)))


def test_things_lammps_cannot_do(sk_root, pot, tmp_path):
    cfg = cfg_for(sk_root)
    assert {e.location for e in errors(cu(pot, m={"units": "", "pair_style": ""}), cfg)} >= {"method.units", "method.pair_style"}
    assert any(e.location == "method.pair_coeff" for e in errors(cu(pot, m={"pair_coeff": ""}), cfg))
    assert any(e.location == "method.pair_coeff" and "パス" in e.message for e in errors(cu(pot, m={"pair_coeff": f"* * {pot}"}), cfg))
    assert any(e.location == "method.potential_files" for e in errors(cu(pot, m={"potential_files": [str(tmp_path / "nope.eam")]}), cfg))
    sp = cu(pot)
    assert any(e.location == "structure.charge" for e in errors(sp.model_copy(update={"structure": sp.structure.model_copy(update={"charge": 1})}), cfg))
    assert any(e.location == "task.md.thermostat" for e in errors(cu(pot, task=Task(type="molecular_dynamics", md=MDSettings(thermostat="andersen"))), cfg))
    assert any(e.location == "task.md.thermostat" for e in errors(cu(pot, task=Task(type="molecular_dynamics", md=MDSettings(ensemble="NPT", thermostat="csvr"))), cfg))
    assert any(e.location == "task.type" for e in errors(cu(pot, task=Task(type="vibrations")), cfg))
    assert any(e.location == "task.max_steps" for e in errors(cu(pot, task=Task(type="geometry_optimization", max_steps=0)), cfg))
    slab = sp.structure.model_copy(update={"atoms": sp.structure.atoms.model_copy(update={"pbc": (True, True, False)})})
    assert any(e.location == "structure.atoms" for e in errors(sp.model_copy(update={"structure": slab}), cfg))


def test_spec_roundtrip_and_cli(sk_root, pot, tmp_path, capsys):
    from adit.cli import main
    spec = cu(pot, m={"extra_commands": "neigh_modify every 1 delay 0", "seed": 7})
    assert CalculationSpec.from_json(spec.to_json()) == spec
    cfg_path = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), cfg_path)
    spec.save(tmp_path / "spec.json")
    assert main([str(tmp_path / "spec.json"), str(tmp_path / "calc"), "--config", str(cfg_path)]) == 0
    assert (tmp_path / "calc" / "Cu_u3.eam").is_file() and "neigh_modify every 1 delay 0" in (tmp_path / "calc" / "in.lammps").read_text(encoding="utf-8")


LMP = shutil.which("lmp") or ""


@pytest.mark.skipif(not (LMP and REAL_EAM.is_file()), reason="lmp か examples/lammps_cu/Cu_u3.eam が無い")
def test_lammps_real_run(sk_root, tmp_path):
    out = tmp_path / "cu"
    write_project(cu(REAL_EAM), cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=300,
                       env={**os.environ, "PATH": f"{Path(LMP).parent}:{os.environ.get('PATH', '')}"})
    assert r.returncode == 0, r.stderr
    assert "Total wall time" in (out / "output.log").read_text(encoding="utf-8")
    traj = (out / "traj.lammpstrj").read_text(encoding="utf-8")
    assert traj.count("ITEM: TIMESTEP") == 11 and "ITEM: ATOMS id type element x y z" in traj and (out / "final.data").is_file()
