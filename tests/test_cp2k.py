
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from adit.codes.cp2k_data import Cp2kData
from adit.config import save_config
from adit.project import ProjectError, build_project, write_project
from adit.spec import AtomsData, CalculationSpec, Cp2kMethod, KPoints, MDSettings, Runtime, Structure, Task
from adit.structure import from_bulk
from tests.conftest import cfg_for, water_spec

BASIS = """# fake basis file (書式は CP2K の BASIS_MOLOPT と同じ)
O  DZVP-MOLOPT-SR-GTH DZVP-MOLOPT-SR-GTH-q6
 1
 2 0 2 1 2 2 1
      1.0  0.1  0.1  0.1  0.1  0.1
O  SZV-MOLOPT-SR-GTH SZV-MOLOPT-SR-GTH-q6
 1
 2 0 1 1 1 1
      1.0  0.1  0.1
H  DZVP-MOLOPT-SR-GTH DZVP-MOLOPT-SR-GTH-q1
 1
 2 0 1 1 2 1
      1.0  0.1  0.1  0.1
H  BAD-q2 BAD-q2
 1
 2 0 0 1 1
      1.0  0.1
Si  DZVP-MOLOPT-SR-GTH DZVP-MOLOPT-SR-GTH-q4
 1
 2 0 2 1 2 2 1
      1.0  0.1  0.1  0.1  0.1  0.1
"""
POTENTIALS = """# fake potential file (書式は CP2K の GTH_POTENTIALS と同じ)
O GTH-PBE-q6 GTH-PBE
    2    4
     0.24455430    2   -16.66721480     2.48731132
    2
     0.22095592    1    18.33745811
     0.21133247    0
#
H GTH-PBE-q1 GTH-PBE
    1
     0.20000000    2    -4.17890044     0.72446331
    0
#
Si GTH-PBE-q4 GTH-PBE
    2    2
     0.44    1    -7.3
    2
     0.42    2     5.9    -1.2
                           3.2
     0.48    1     2.7
#
Si GTH-BLYP-q4 GTH-BLYP
    2    2
     0.44    1    -7.3
    0
"""


def make_data(root: Path) -> Path:
    root.mkdir(parents=True)
    (root / "BASIS_MOLOPT").write_text(BASIS, encoding="utf-8")
    (root / "GTH_POTENTIALS").write_text(POTENTIALS, encoding="utf-8")
    (root / "dftd3.dat").write_text("fake\n", encoding="utf-8")
    return root


@pytest.fixture
def cfg_cp2k(sk_root, tmp_path):
    cfg = cfg_for(sk_root)
    cfg.cp2k_data = str(make_data(tmp_path / "cp2k_data"))
    return cfg


def h2o(**kw) -> CalculationSpec:
    m = dict(xc="PBE", basis={"O": "DZVP-MOLOPT-SR-GTH", "H": "DZVP-MOLOPT-SR-GTH"}, cutoff_ry=280, rel_cutoff_ry=40, poisson_solver="MT",
             isolated_box_ang=8.0)
    m.update(kw.pop("m", {}))
    return water_spec(method=Cp2kMethod(**m), task=kw.pop("task", Task(type="single_point")),
                      runtime=Runtime(profile="local", mpiprocs=1, omp_threads=1, job_name="h2o"), **kw)


def si(**kw) -> CalculationSpec:
    m = dict(xc="PBE", potential={"Si": "GTH-PBE-q4"}, cutoff_ry=300, rel_cutoff_ry=50)
    m.update(kw.pop("m", {}))
    return CalculationSpec(structure=Structure(source="bulk", source_ref="Si", atoms=AtomsData.from_ase(from_bulk("Si"))),
                           method=Cp2kMethod(**m), kpoints=kw.pop("kpoints", KPoints(mode="mesh", mesh=(2, 2, 2))),
                           task=kw.pop("task", Task(type="single_point")), **kw)


def errors(spec, cfg) -> list:
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg)
    return ex.value.errors


def test_data_reader(tmp_path):
    d = Cp2kData(make_data(tmp_path / "d"))
    assert d.names_for("BASIS_MOLOPT", "O", potential=False) == ["DZVP-MOLOPT-SR-GTH", "SZV-MOLOPT-SR-GTH"]
    assert d.names_for("GTH_POTENTIALS", "Si", potential=True) == ["GTH-PBE-q4", "GTH-BLYP-q4"]
    si_pp = d.find("GTH_POTENTIALS", "Si", "GTH-PBE", potential=True)
    assert si_pp.names[0] == "GTH-PBE-q4" and si_pp.valence == 4 and "3.2" in si_pp.text
    assert d.find("GTH_POTENTIALS", "O", "GTH-PBE-q6", potential=True).valence == 6
    assert d.find("BASIS_MOLOPT", "H", "DZVP-MOLOPT-SR-GTH", potential=False).q == 1
    assert d.names_for("BASIS_MOLOPT", "Xe", potential=False) == []


def test_h2o_single_point(cfg_cp2k):
    files = build_project(h2o(), cfg_cp2k)
    inp = files.texts["cp2k.inp"]
    for key in ["RUN_TYPE ENERGY_FORCE", "BASIS_SET_FILE_NAME BASIS_adit", "POTENTIAL_FILE_NAME POTENTIAL_adit", "CUTOFF 280", "REL_CUTOFF 40",
                "EPS_SCF 1e-05", "MAX_SCF 50", "&XC_FUNCTIONAL PBE", "PERIODIC NONE", "POISSON_SOLVER MT", "&CENTER_COORDINATES",
                "A 8.0000000000 0.0000000000 0.0000000000", "&KIND O", "BASIS_SET DZVP-MOLOPT-SR-GTH", "POTENTIAL GTH-PBE-q6",
                "POTENTIAL GTH-PBE-q1", "MULTIPLICITY 1", "&FORCES ON"]:
        assert key in inp, key
    assert "&MOTION" not in inp and "&KPOINTS" not in inp and "STRESS_TENSOR" not in inp
    assert "DZVP-MOLOPT-SR-GTH-q6" in files.texts["BASIS_adit"] and "SZV" not in files.texts["BASIS_adit"] and "BAD" not in files.texts["BASIS_adit"]
    assert "GTH-PBE-q1" in files.texts["POTENTIAL_adit"] and "Si" not in files.texts["POTENTIAL_adit"]
    assert files.copies == {}
    assert files.texts["submit.sh"].rstrip().endswith("mpirun -np 1 cp2k.psmp -i cp2k.inp > output.log 2>&1")
    assert "CP2K" in files.texts["README.txt"] and "BASIS_adit" in files.texts["README.txt"]


def test_required_values_have_no_defaults(cfg_cp2k):
    errs = errors(h2o(m={"xc": "", "cutoff_ry": 0, "rel_cutoff_ry": 0, "poisson_solver": "", "isolated_box_ang": 0}), cfg_cp2k)
    locs = {e.location for e in errs}
    assert {"method.xc", "method.cutoff_ry", "method.rel_cutoff_ry", "method.poisson_solver", "method.isolated_box_ang"} <= locs


def test_names_are_picked_only_when_unique(cfg_cp2k):
    errs = errors(h2o(m={"basis": {}}), cfg_cp2k)
    msgs = [e.message for e in errs if e.location == "method.basis"]
    assert any("DZVP-MOLOPT-SR-GTH, SZV-MOLOPT-SR-GTH" in m for m in msgs)
    errs = errors(si(m={"potential": {}}), cfg_cp2k)
    assert any(e.location == "method.potential" and "GTH-BLYP-q4" in e.message for e in errs)
    errs = errors(h2o(m={"basis": {"O": "DZVP-MOLOPT-SR-GTH", "H": "NOPE"}}), cfg_cp2k)
    assert any(e.location == "method.basis" and "NOPE" in e.message for e in errs)


def test_basis_and_potential_valence_must_agree(cfg_cp2k):
    errs = errors(h2o(m={"basis": {"O": "DZVP-MOLOPT-SR-GTH", "H": "BAD-q2"}}), cfg_cp2k)
    assert any(e.location == "method.basis" and "-q2" in e.message for e in errs)


def test_spin_and_parity(cfg_cp2k):
    spec = h2o(m={"basis": {"O": "DZVP-MOLOPT-SR-GTH", "H": "DZVP-MOLOPT-SR-GTH"}})
    st = spec.structure
    errs = errors(spec.model_copy(update={"structure": st.model_copy(update={"charge": 1})}), cfg_cp2k)
    assert any(e.location == "structure.multiplicity" for e in errs)
    errs = errors(spec.model_copy(update={"structure": st.model_copy(update={"charge": 1, "multiplicity": 2})}), cfg_cp2k)
    assert [e.location for e in errs] == ["method.uks"]
    ok = spec.model_copy(update={"structure": st.model_copy(update={"charge": 1, "multiplicity": 2}),
                                 "method": spec.method.model_copy(update={"uks": True})})
    inp = build_project(ok, cfg_cp2k).texts["cp2k.inp"]
    assert "UKS .TRUE." in inp and "CHARGE 1" in inp and "MULTIPLICITY 2" in inp


def test_periodic_cell_opt_and_npt(cfg_cp2k):
    inp = build_project(si(task=Task(type="geometry_optimization", relax_cell="shape_and_volume", max_steps=30)), cfg_cp2k).texts["cp2k.inp"]
    for key in ["RUN_TYPE CELL_OPT", "STRESS_TENSOR ANALYTICAL", "&CELL_OPT", "MAX_ITER 30", "EXTERNAL_PRESSURE [bar] 0.0",
                "SCHEME MONKHORST-PACK 2 2 2", "PERIODIC XYZ", "MAX_FORCE 1.000000e-04"]:
        assert key in inp, key
    assert "&POISSON" not in inp and "CENTER_COORDINATES" not in inp
    md = MDSettings(ensemble="NPT", thermostat="csvr", temperature_k=500, timestep_fs=1.0, steps=10, dump_interval=2,
                    coupling_time_fs=50, pressure_bar=2.0, barostat_time_fs=500)
    inp = build_project(si(task=Task(type="molecular_dynamics", md=md)), cfg_cp2k).texts["cp2k.inp"]
    for key in ["RUN_TYPE MD", "ENSEMBLE NPT_I", "STEPS 10", "TIMESTEP 1", "TEMPERATURE 500", "TYPE CSVR", "TIMECON 50",
                "&BAROSTAT", "PRESSURE 2", "TIMECON 500", "MD 2", "&CELL ON", "STRESS_TENSOR ANALYTICAL"]:
        assert key in inp, key
    lang = build_project(si(task=Task(type="molecular_dynamics", md=MDSettings(ensemble="NVT", thermostat="langevin", coupling_time_fs=200))),
                         cfg_cp2k).texts["cp2k.inp"]
    assert "ENSEMBLE LANGEVIN" in lang and "GAMMA 0.005" in lang


def test_things_cp2k_cannot_do(cfg_cp2k):
    assert any(e.location == "task.relax_cell" for e in errors(si(task=Task(type="geometry_optimization", relax_cell="volume_only")), cfg_cp2k))
    assert any(e.location == "task.md.thermostat" for e in errors(si(task=Task(type="molecular_dynamics", md=MDSettings(thermostat="berendsen"))), cfg_cp2k))
    assert any(e.location == "task.md.thermostat" for e in errors(si(task=Task(type="molecular_dynamics", md=MDSettings(ensemble="NPT", thermostat="langevin"))), cfg_cp2k))
    assert any(e.location == "kpoints.shift" for e in errors(si(kpoints=KPoints(mode="mesh", mesh=(2, 2, 2), shift=(0.5, 0.5, 0.5))), cfg_cp2k))
    assert any(e.location == "task.type" for e in errors(si(task=Task(type="band_structure")), cfg_cp2k))
    assert any(e.location == "method.extra_sections" for e in errors(si(m={"extra_sections": {"FORCE_EVAL/DFT/NOPE": "X 1"}}), cfg_cp2k))


def test_fixed_atoms_and_extra_sections(cfg_cp2k):
    spec = h2o(m={"basis": {"O": "DZVP-MOLOPT-SR-GTH", "H": "DZVP-MOLOPT-SR-GTH"}, "extra_sections": {"force_eval/dft/scf": "SCF_GUESS RESTART"}},
               task=Task(type="geometry_optimization", optimizer="LBFGS", max_steps=0))
    spec = spec.model_copy(update={"structure": spec.structure.model_copy(update={"fixed_atoms": [0], "fixed_axes": {"1": (True, True, False)}})})
    inp = build_project(spec, cfg_cp2k).texts["cp2k.inp"]
    assert "OPTIMIZER LBFGS" in inp and "MAX_ITER" not in inp
    assert "LIST 1\n      COMPONENTS_TO_FIX XYZ" in inp and "LIST 2\n      COMPONENTS_TO_FIX Z" in inp
    assert "      SCF_GUESS RESTART\n    &END SCF" in inp


def test_without_data_directory(sk_root, tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path / "empty"))
    monkeypatch.delenv("CP2K_DATA_DIR", raising=False)
    cfg = cfg_for(sk_root)
    missing = errors(h2o(m={"basis": {}, "potential": {}}), cfg)
    assert any(e.location == "method.potential" for e in missing)
    assert len({e.message for e in missing}) == len(missing)
    assert {e.location for e in missing} >= {"method.basis", "method.potential"}
    files = build_project(h2o(m={"basis": {"O": "DZVP-MOLOPT-SR-GTH", "H": "DZVP-MOLOPT-SR-GTH"}, "potential": {"O": "GTH-PBE-q6", "H": "GTH-PBE-q1"}}), cfg)
    assert "BASIS_SET_FILE_NAME BASIS_MOLOPT" in files.texts["cp2k.inp"] and "BASIS_adit" not in files.texts
    assert "未確認" in files.texts["README.txt"]
    cfg.cp2k_data = str(tmp_path / "missing")
    assert any("cp2k_data" in e.message for e in errors(h2o(), cfg))


def test_spec_roundtrip_and_cli(cfg_cp2k, tmp_path, capsys):
    from adit.cli import main
    spec = h2o(m={"dispersion": "d3bj", "extra_sections": {"FORCE_EVAL/DFT/SCF": "SCF_GUESS ATOMIC"}})
    assert CalculationSpec.from_json(spec.to_json()) == spec
    cfg_path = tmp_path / "cluster.toml"
    save_config(cfg_cp2k, cfg_path)
    spec.save(tmp_path / "spec.json")
    assert main([str(tmp_path / "spec.json"), "--print", "cp2k.inp", "--config", str(cfg_path)]) == 0
    out = capsys.readouterr().out
    assert "TYPE DFTD3(BJ)" in out and "REFERENCE_FUNCTIONAL PBE" in out and "PARAMETER_FILE_NAME dftd3.dat" in out
    assert main([str(tmp_path / "spec.json"), str(tmp_path / "calc"), "--config", str(cfg_path)]) == 0
    assert (tmp_path / "calc" / "BASIS_adit").is_file()


def test_ot_and_surface_dipole_are_written_only_when_selected(cfg_cp2k):
    spec = h2o(m={"ot": True, "ot_minimizer": "DIIS", "ot_preconditioner": "FULL_SINGLE_INVERSE",
                  "ot_extra": "ENERGY_GAP 0.001", "surface_dipole_correction": True, "surf_dip_dir": "Z"})
    inp = build_project(spec, cfg_cp2k).texts["cp2k.inp"]
    for line in ("&OT", "MINIMIZER DIIS", "PRECONDITIONER FULL_SINGLE_INVERSE", "ENERGY_GAP 0.001",
                 "SURFACE_DIPOLE_CORRECTION .TRUE.", "SURF_DIP_DIR Z"):
        assert line in inp


CP2K = shutil.which("cp2k.psmp") or ""


@pytest.mark.skipif(not (CP2K and Cp2kData.discover().found), reason="cp2k.psmp or its data directory is missing")
def test_cp2k_real_run(sk_root, tmp_path):
    from ase.build import molecule
    cfg = cfg_for(sk_root)
    spec = CalculationSpec(structure=Structure(source="preset", source_ref="H2O", atoms=AtomsData.from_ase(molecule("H2O"))),
                           method=Cp2kMethod(xc="PBE", basis={"O": "DZVP-MOLOPT-SR-GTH", "H": "DZVP-MOLOPT-SR-GTH"},
                                             potential={"O": "GTH-PBE-q6", "H": "GTH-PBE-q1"}, cutoff_ry=280, rel_cutoff_ry=40,
                                             poisson_solver="MT", isolated_box_ang=6.0),
                           task=Task(type="single_point"), runtime=Runtime(profile="local", mpiprocs=1, omp_threads=1))
    out = tmp_path / "h2o"
    write_project(spec, cfg, out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=600,
                       env={**os.environ, "PATH": f"{Path(CP2K).parent}:{os.environ.get('PATH', '')}"})
    assert r.returncode == 0, r.stderr
    log = (out / "output.log").read_text(encoding="utf-8")
    line = next(l for l in log.splitlines() if "ENERGY| Total FORCE_EVAL" in l)
    assert abs(float(line.split()[-1]) - (-17.219399129357214)) < 1e-5
