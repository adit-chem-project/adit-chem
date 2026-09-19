
import os
import subprocess
from pathlib import Path

import pytest

from adit.codes.potcar import PotcarLibrary
from adit.codes.vasp import VaspGenerator, potcar_name
from adit.config import Profile
from adit.project import ProjectError, build_project, write_project
from adit.spec import AtomsData, CalculationSpec, KPoints, MDSettings, Runtime, Structure, Task, VaspMethod
from adit.structure import from_file
from tests.conftest import cfg_for

REPO = Path(__file__).resolve().parent.parent
gen = VaspGenerator()

FAKE_POTCARS = {
    "H": ("PAW_PBE H 15Jun2001", 1.0, 250.0),
    "O": ("PAW_PBE O 08Apr2002", 6.0, 400.0),
    "Si": ("PAW_PBE Si 05Jan2001", 4.0, 245.3),
    "Li_sv": ("PAW_PBE Li_sv 10Sep2004", 3.0, 499.0),
}


def make_fake_pp(root: Path, potcar_set="potpaw_PBE") -> Path:
    for name, (titel, zval, enmax) in FAKE_POTCARS.items():
        d = root / potcar_set / name
        d.mkdir(parents=True)
        (d / "POTCAR").write_text(f"  PAW_PBE {name.split('_')[0]} fake\n   {zval} ; ZVAL   = {zval}    mass and valenz\n"
                                  f"   TITEL  = {titel}\n   ENMAX  =  {enmax}; ENMIN  =  {enmax * 0.75} eV\n(data omitted)\n", encoding="utf-8")
    return root


def h2o_spec(**kw) -> CalculationSpec:
    s = 0.52918
    atoms = AtomsData(symbols=["O", "H", "H"], positions=[(0, 0, 0), (1.10 * s, -1.43 * s, 0), (1.10 * s, 1.43 * s, 0)],
                      cell=[(15 * s, 0, 0), (0, 15 * s, 0), (0, 0, 15 * s)], pbc=(True, True, True))
    st = Structure(source="file", source_ref="examples/vasp_h2o/POSCAR (数値を写した)", atoms=atoms, fixed_atoms=[0])
    d = dict(structure=st, method=VaspMethod(encut=400, prec="Normal", ismear=0, sigma=0.1), kpoints=KPoints(mode="gamma"),
             task=Task(type="geometry_optimization", max_steps=10, force_tolerance_ev_per_ang=0.02),
             runtime=Runtime(profile="cluster", ncpus=8, mpiprocs=8, omp_threads=1, job_name="h2o"))
    d.update(kw)
    return CalculationSpec(**d)


def si_spec(**kw) -> CalculationSpec:
    a = 5.5
    cell = [(0, 0.5 * a, 0.5 * a), (0.5 * a, 0, 0.5 * a), (0.5 * a, 0.5 * a, 0)]
    import numpy as np
    frac = np.array([(-0.125, -0.125, -0.125), (0.125, 0.125, 0.130)])
    pos = [tuple(float(x) for x in f @ np.array(cell)) for f in frac]
    atoms = AtomsData(symbols=["Si", "Si"], positions=pos, cell=cell, pbc=(True, True, True))
    st = Structure(source="file", source_ref="examples/vasp_cd_si/POSCAR (数値を写した)", atoms=atoms)
    d = dict(structure=st, method=VaspMethod(encut=240, ismear=0, sigma=0.1, extra_incar={"ISTART": 0, "ICHARG": 2}),
             kpoints=KPoints(mode="mesh", mesh=(11, 11, 11)),
             task=Task(type="geometry_optimization", max_steps=10, force_tolerance_ev_per_ang=0.0001),
             runtime=Runtime(profile="cluster", ncpus=8, mpiprocs=8, omp_threads=1, job_name="si"))
    d.update(kw)
    return CalculationSpec(**d)


@pytest.fixture
def cfg_pp(sk_root, tmp_path):
    cfg = cfg_for(sk_root)
    pp = make_fake_pp(tmp_path / "pp")
    cfg.profiles["cluster"].env["VASP_PP_PATH"] = str(pp)
    cfg.profiles["cluster"].commands["vasp"] = "mpirun -np {mpiprocs} /opt/vasp/bin/vasp_{binary}"
    return cfg


def test_potcar_library_headers(tmp_path):
    lib = PotcarLibrary(make_fake_pp(tmp_path), "potpaw_PBE")
    assert lib.names_for("Li") == ["Li_sv"] and lib.names_for("O") == ["O"] and lib.names_for("Xx") == []
    h = lib.header("Si")
    assert (h.titel, h.zval, h.enmax) == ("PAW_PBE Si 05Jan2001", 4.0, 245.3)
    assert PotcarLibrary.open_if_present(None, "potpaw_PBE") is None
    assert PotcarLibrary.open_if_present(str(tmp_path), "nope") is None


def test_potcar_name_precedence():
    m = VaspMethod(potcar={"O": "O_h"})
    assert potcar_name(m, "O") == "O_h" and potcar_name(m, "Li") == "Li_sv" and potcar_name(m, "Xx") == "Xx"


def test_h2o_files_match_wiki(cfg_pp):
    files = build_project(h2o_spec(), cfg_pp)
    inc, pos, kp = files.texts["INCAR"], files.texts["POSCAR"], files.texts["KPOINTS"]
    for key in ["ENCUT = 400", "PREC = Normal", "ISMEAR = 0", "SIGMA = 0.1", "NSW = 10", "IBRION = 2", "ISIF = 2", "EDIFFG = -0.02"]:
        assert key in inc, key
    lines = pos.splitlines()
    assert lines[5].split() == ["O", "H"] and lines[6].split() == ["1", "2"]
    assert lines[7].startswith("Selective") and lines[8].startswith("Direct")
    assert lines[9].split()[3:] == ["F", "F", "F"] and lines[10].split()[3:] == ["T", "T", "T"]
    assert kp.splitlines()[2] == "Gamma" and kp.splitlines()[3].split() == ["1", "1", "1"]
    spec_txt = files.texts["potcar.spec"]
    assert "O  O  'PAW_PBE O 08Apr2002'  6" in spec_txt and "H  H  'PAW_PBE H 15Jun2001'  1" in spec_txt
    assert files.copies == {}
    assert "POTCAR" not in files.texts
    assert files.texts["submit.sh"].rstrip().endswith("bash make_potcar.sh && mpirun -np 8 /opt/vasp/bin/vasp_std > output.log 2>&1")
    assert 'export VASP_PP_PATH="' in files.texts["submit.sh"]
    assert "ライセンス保持者" in files.texts["README.txt"]


def test_si_files_match_wiki(cfg_pp):
    files = build_project(si_spec(), cfg_pp)
    inc, kp = files.texts["INCAR"], files.texts["KPOINTS"]
    for key in ["ENCUT = 240", "ISTART = 0", "ICHARG = 2", "NSW = 10", "IBRION = 2", "ISIF = 2", "EDIFFG = -0.0001"]:
        assert key in inc, key
    assert kp.splitlines()[2].startswith("Monkhorst") and kp.splitlines()[3].split() == ["11", "11", "11"]
    assert files.texts["POSCAR"].splitlines()[5].split() == ["Si"]


@pytest.mark.parametrize("ensemble,thermostat,extra", [
    ("NVE", "berendsen", {}),
    ("NVT", "nose_hoover", {"SMASS": 2}),
    ("NPT", "langevin", {"PMASS": 800, "LANGEVIN_GAMMA_L": 4}),
])
def test_md_writes_xdatcar_interval(cfg_pp, ensemble, thermostat, extra):
    spec = si_spec(method=VaspMethod(extra_incar=extra), task=Task(type="molecular_dynamics",
        md=MDSettings(ensemble=ensemble, thermostat=thermostat, dump_interval=3)))
    files = build_project(spec, cfg_pp)
    assert "NBLOCK = 3" in files.texts["INCAR"]
    assert "NBLOCK" in files.texts["README.txt"]


@pytest.mark.parametrize("task_type", ["single_point", "geometry_optimization", "vibrations", "band_structure"])
def test_non_md_has_no_generated_nblock(cfg_pp, task_type):
    spec = si_spec(task=Task(type=task_type, md=MDSettings(dump_interval=3)))
    files = build_project(spec, cfg_pp)
    assert "NBLOCK" not in files.texts["INCAR"]


def test_gamma_centred_mesh_can_be_selected(cfg_pp):
    base = si_spec()
    spec = base.model_copy(update={"method": base.method.model_copy(update={"kpoints_centering": "gamma"}),
                                   "kpoints": KPoints(mode="mesh", mesh=(4, 4, 4))})
    assert build_project(spec, cfg_pp).texts["KPOINTS"].splitlines()[2] == "Gamma"


def test_without_library_generation_is_allowed_but_unverified(sk_root):
    cfg = cfg_for(sk_root)
    files = build_project(h2o_spec(), cfg)
    assert "(未確認)" in files.texts["potcar.spec"] and "未確認" in files.texts["README.txt"]
    charged = h2o_spec(structure=Structure(**{**h2o_spec().structure.model_dump(), "charge": 1}))
    with pytest.raises(ProjectError) as ex:
        build_project(charged, cfg)
    assert ex.value.errors[0].location == "structure.charge"


def test_library_checks_names_and_parity(cfg_pp):
    bad = h2o_spec(method=VaspMethod(potcar={"O": "O_nope"}))
    with pytest.raises(ProjectError) as ex:
        build_project(bad, cfg_pp)
    assert ex.value.errors[0].location == "method.potcar" and "O_nope" in ex.value.errors[0].message
    odd = h2o_spec(structure=Structure(**{**h2o_spec().structure.model_dump(), "multiplicity": 2}), method=VaspMethod(ispin=2))
    with pytest.raises(ProjectError) as ex:
        build_project(odd, cfg_pp)
    assert ex.value.errors[0].location == "structure.multiplicity"
    cation = h2o_spec(structure=Structure(**{**h2o_spec().structure.model_dump(), "charge": 1, "multiplicity": 2}), method=VaspMethod(ispin=2))
    inc = build_project(cation, cfg_pp).texts["INCAR"]
    assert "NELECT = 7" in inc and "NUPDOWN = 1" in inc and "ISPIN = 2" in inc


def test_molecule_without_box_is_error(sk_root):
    from tests.conftest import water_spec
    spec = water_spec(method=VaspMethod())
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_for(sk_root))
    assert any(e.location == "structure.atoms" and "周期セル" in e.message for e in ex.value.errors)


def test_extra_incar_and_relax_cell(cfg_pp):
    spec = si_spec(method=VaspMethod(encut=240, extra_incar={"lasph": True, "NCORE": 4, "GGA": "PS"}),
                   task=Task(type="geometry_optimization", relax_cell="shape_and_volume", max_steps=5))
    inc = build_project(spec, cfg_pp).texts["INCAR"]
    assert "LASPH = .TRUE." in inc and "NCORE = 4" in inc and "GGA = PS" in inc and "ISIF = 3" in inc


@pytest.mark.skipif(os.name == "nt", reason="生成した bash スクリプトを実行する試験 (走り先は Linux)")
def test_make_potcar_assembles_from_fake_library(cfg_pp, tmp_path):
    out = tmp_path / "h2o"
    write_project(h2o_spec(), cfg_pp, out)
    env = {**os.environ, "VASP_PP_PATH": cfg_pp.profiles["cluster"].env["VASP_PP_PATH"]}
    r = subprocess.run(["bash", "make_potcar.sh"], cwd=out, capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stderr
    potcar = (out / "POTCAR").read_text(encoding="utf-8")
    assert potcar.index("TITEL  = PAW_PBE O") < potcar.index("TITEL  = PAW_PBE H")
    assert (out / "potcar.used").read_text(encoding="utf-8").count("\n") == 2
    r2 = subprocess.run(["bash", "make_potcar.sh"], cwd=out, capture_output=True, text=True, env={**os.environ, "VASP_PP_PATH": ""})
    assert r2.returncode != 0 and "VASP_PP_PATH" in r2.stderr


def test_axis_constraints_and_ibrion(cfg_pp):
    st = h2o_spec().structure
    st2 = Structure(**{**st.model_dump(), "fixed_axes": {"1": (True, True, False), "2": (True, True, False)}})
    spec = h2o_spec(structure=st2, method=VaspMethod(encut=400, ibrion=1))
    files = build_project(spec, cfg_pp)
    pos = files.texts["POSCAR"].splitlines()
    assert pos[9].split()[3:] == ["F", "F", "F"] and pos[10].split()[3:] == ["T", "T", "F"] and pos[11].split()[3:] == ["T", "T", "F"]
    assert "IBRION = 1" in files.texts["INCAR"]


def test_axis_constraints_rejected_by_dftb(sk_root):
    from adit.codes.base import GenerationError
    from adit.codes.dftbplus import DftbPlusGenerator
    from adit.codes.sk_sets import SKSet
    from tests.conftest import water_spec
    st = water_spec().structure
    spec = water_spec(structure=Structure(**{**st.model_dump(), "fixed_axes": {"1": (True, True, False)}}))
    with pytest.raises(GenerationError):
        DftbPlusGenerator().dftb_in_hsd(spec, SKSet.from_dir(sk_root / "fake-1-0"))


@pytest.mark.parametrize("cell,movable,expected", [
    ([(0, 13, 0), (-12, 0, 0), (0, 0, 14)], (False, True, False), ["T", "F", "F"]),
    ([(12, 0, 0), (3, 13, 0), (0, 0, 14)], (True, True, False), ["T", "T", "F"]),
])
def test_axis_constraints_follow_cartesian_motion_after_element_sort(cfg_pp, cell, movable, expected):
    spec = h2o_spec()
    spec.structure.atoms.cell = cell
    spec.structure.atoms.symbols = ["H", "O", "H"]
    spec.structure.fixed_atoms = [0]
    spec.structure.fixed_axes = {"0": (True, False, True), "2": movable}
    before = spec.model_dump()
    pos = build_project(spec, cfg_pp).texts["POSCAR"].splitlines()
    assert pos[5].split() == ["H", "O"]
    assert [line.split()[3:] for line in pos[9:12]] == [["F"] * 3, expected, ["T"] * 3]
    assert spec.model_dump() == before


@pytest.mark.parametrize("cell", [
    [(12, 0, 0), (3, 13, 0), (0, 0, 14)],
    [(8, 8, 0), (-8, 8, 0), (0, 0, 14)],
])
def test_vasp_validation_rejects_axis_constraints_that_change_direction(cfg_pp, cell):
    spec = h2o_spec()
    spec.structure.atoms.cell = cell
    spec.structure.fixed_axes = {"1": (False, True, False)}
    with pytest.raises(ProjectError) as error:
        build_project(spec, cfg_pp)
    assert any(e.location == "structure.fixed_axes" and "Selective dynamics" in e.message for e in error.value.errors)
    # A fully fixed atom overrides any overlapping individual-axis entries.
    spec.structure.fixed_atoms = [0, 1]
    assert build_project(spec, cfg_pp).texts["POSCAR"].splitlines()[10].split()[3:] == ["F"] * 3


@pytest.mark.parametrize("task,extra", [
    (Task(type="geometry_optimization", relax_cell="shape_and_volume"), {}),
    (Task(type="molecular_dynamics", md=MDSettings(ensemble="NPT", thermostat="langevin")),
     {"LANGEVIN_GAMMA_L": 1, "PMASS": 100}),
    (Task(type="geometry_optimization"), {"ISIF": 4}),
    (Task(type="geometry_optimization"), {"ISIF": 4.0}),
])
def test_vasp_rejects_cartesian_axes_when_cell_directions_can_change(cfg_pp, task, extra):
    spec = h2o_spec(task=task, method=VaspMethod(extra_incar=extra))
    spec.structure.fixed_axes = {"1": (True, True, False)}
    with pytest.raises(ProjectError, match="Cartesian"):
        build_project(spec, cfg_pp)


def test_vasp_keeps_axis_constraints_when_only_cell_volume_changes(cfg_pp):
    spec = h2o_spec(task=Task(type="geometry_optimization", relax_cell="volume_only"))
    spec.structure.fixed_axes = {"1": (True, True, False)}
    assert build_project(spec, cfg_pp).texts["POSCAR"].splitlines()[10].split()[3:] == ["T", "T", "F"]


def test_parse_constraints():
    from adit.gui.panels.structure_panel import parse_constraints
    assert parse_constraints("1-2, 3:xy, 4:z", 5) == ([0, 1], {"2": (False, False, True), "3": (True, True, False)})
    with pytest.raises(ValueError):
        parse_constraints("3:q", 5)


@pytest.mark.parametrize("coupling_fs,timestep_fs,period", [(100.0, 2.0, 50), (100.0, 1.0, 100), (0.5, 1.0, 1)])
def test_csvr_period_is_a_number_of_md_steps(cfg_pp, coupling_fs, timestep_fs, period):
    spec = si_spec(task=Task(type="molecular_dynamics", md=MDSettings(
        ensemble="NVT", thermostat="csvr", coupling_time_fs=coupling_fs, timestep_fs=timestep_fs)))
    inc = build_project(spec, cfg_pp).texts["INCAR"]
    assert "MDALGO = 5" in inc and f"CSVR_PERIOD = {period}\n" in inc
    assert f"POTIM = {timestep_fs:g}" in inc
