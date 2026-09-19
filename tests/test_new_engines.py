"""Seven CLI-only engine subsets. Input syntax is traced in examples/new_engines/SOURCES.md.

No proprietary executable or parameter set is bundled; these are generation and
mechanical-validation tests, not claims of completed calculations.
"""

from pathlib import Path

import pytest

from tests.conftest import cfg_for, water_spec
from adit.cli import main as generate_main
from adit.config import save_config
from adit.convert import retarget_spec
from adit import lang
from adit.project import ProjectError, build_project, write_project
from adit.spec import (AmberMethod, AtomsData, CalculationSpec, GamessMethod,
                        GaussianMethod, GrrmMethod, KPoints, NamdMethod,
                        OpenmxMethod, OrcaMethod, QchemMethod, Runtime, Structure, Task, MDSettings)


def molecular(method, task="single_point"):
    return water_spec(method=method, task=Task(type=task), runtime=Runtime(profile="local"))


@pytest.mark.parametrize("method", [
    GaussianMethod(theory="HF", basis="6-31G(d)"),
    GamessMethod(gbasis="STO", ngauss=3),
    QchemMethod(theory="HF", basis="6-31G"),
    GrrmMethod(theory="HF", basis="6-31G"),
    OpenmxMethod(data_path="../DFT_DATA", pao={"H": "H5.0-s1p1"},
                 vps={"H": "H_CA19"}, valence={"H": 1}, xc="LDA", energycutoff_ry=220),
    AmberMethod(topology_file="system.prmtop", coordinates_file="system.rst7", cutoff_ang=9, igb=0),
    NamdMethod(structure_file="system.psf", coordinates_file="system.pdb", parameter_files=["force.prm"],
               exclude="scaled1-4", one_four_scaling=1, cutoff_ang=10, pairlistdist_ang=12, switching=False),
])
def test_new_method_json_round_trip(method):
    spec = molecular(method)
    assert CalculationSpec.from_json(spec.to_json()).method == method


@pytest.mark.parametrize(("method", "name", "fragment"), [
    (GaussianMethod(theory="HF", basis="6-31G(d)"), "gaussian.gjf", "# HF/6-31G(d) SP"),
    (GamessMethod(gbasis="STO", ngauss=3), "gamess.inp", "$BASIS GBASIS=STO NGAUSS=3 $END"),
    (QchemMethod(theory="HF", basis="6-31G"), "qchem.in", "METHOD HF\nBASIS 6-31G"),
    (GrrmMethod(theory="HF", basis="6-31G"), "grrm.com", "# MIN/HF/6-31G"),
])
def test_molecular_generators_round_trip_and_build(sk_root, tmp_path, method, name, fragment):
    task = "geometry_optimization" if method.code == "grrm" else "single_point"
    spec = molecular(method, task)
    spec.save(tmp_path / "spec.json")
    loaded = CalculationSpec.load(tmp_path / "spec.json")
    assert loaded.method == method
    files = build_project(loaded, cfg_for(sk_root))
    assert fragment in files.texts[name]
    assert "submit.sh" in files.texts and "README.txt" in files.texts
    assert "analyze.py" not in files.texts and "analyze.py" not in files.texts["README.txt"]
    assert "output.log" not in files.texts


@pytest.mark.parametrize("method", [GaussianMethod(), GamessMethod(), QchemMethod(), GrrmMethod()])
def test_molecular_missing_method_stops(sk_root, method):
    task = "geometry_optimization" if method.code == "grrm" else "single_point"
    with pytest.raises(ProjectError):
        build_project(molecular(method, task), cfg_for(sk_root))


def test_molecular_rejects_injection_and_lost_constraints(sk_root):
    with pytest.raises(ProjectError):
        build_project(molecular(GaussianMethod(theory="HF\n%chk=other", basis="6-31G")), cfg_for(sk_root))
    spec = molecular(QchemMethod(theory="HF", basis="6-31G"))
    spec.structure.fixed_atoms = [0]
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_for(sk_root))
    assert any(e.location == "structure.fixed_atoms" for e in ex.value.errors)


def test_openmx_molecular_and_periodic(sk_root):
    method = OpenmxMethod(data_path="../DFT_DATA", pao={"O": "O6.0-s2p2", "H": "H5.0-s2"},
                          vps={"O": "O_PBE", "H": "H_PBE"}, valence={"O": 6, "H": 1},
                          xc="GGA-PBE", energycutoff_ry=220)
    spec = molecular(method)
    inp = build_project(spec, cfg_for(sk_root)).texts["openmx.dat"]
    assert "scf.EigenvalueSolver Cluster" in inp
    assert "O 6.0-s2p2" not in inp
    assert "O O6.0-s2p2 O_PBE" in inp
    spec.structure.atoms = AtomsData(symbols=["O", "H", "H"],
        positions=[(0.0, -1.0, 0.0), (0.0, 0.0, 0.783064), (0.0, 0.0, -0.783064)],
        cell=[(12, 0, 0), (0, 12, 0), (0, 0, 12)], pbc=(True, True, True))
    spec.kpoints = KPoints(mode="mesh", mesh=(2, 2, 2))
    inp = build_project(spec, cfg_for(sk_root)).texts["openmx.dat"]
    assert "scf.EigenvalueSolver Band" in inp and "scf.Kgrid 2 2 2" in inp


def test_openmx_requires_user_pao_vps_valence(sk_root):
    with pytest.raises(ProjectError) as ex:
        build_project(molecular(OpenmxMethod(xc="GGA-PBE", energycutoff_ry=220)), cfg_for(sk_root))
    assert {"method.data_path", "method.pao", "method.vps", "method.valence"} <= {e.location for e in ex.value.errors}


def amber_files(tmp_path: Path) -> tuple[Path, Path]:
    top, crd = tmp_path / "system.prmtop", tmp_path / "system.rst7"
    top.write_text("%VERSION V0001\n%FLAG POINTERS\n%FORMAT(10I8)\n       3\n", encoding="ascii")
    crd.write_text("test\n     3\n 0.0 -1.0 0.0  0.0 0.0 0.783064\n 0.0 0.0 -0.783064\n", encoding="ascii")
    return top, crd


def test_amber_minimization_copies_user_files(sk_root, tmp_path):
    top, crd = amber_files(tmp_path)
    spec = molecular(AmberMethod(topology_file=str(top), coordinates_file=str(crd), cutoff_ang=9, igb=0),
                     "geometry_optimization")
    files = build_project(spec, cfg_for(sk_root))
    assert "imin=1" in files.texts["amber.in"] and "cut=9" in files.texts["amber.in"]
    assert files.copies["topology.prmtop"] == top
    assert files.copies["coordinates.rst7"] == crd
    assert "sander -O -i amber.in" in files.texts["submit.sh"]
    out = tmp_path / "amber_project"
    write_project(spec, cfg_for(sk_root), out)
    assert (out / "topology.prmtop").read_bytes() == top.read_bytes()
    assert (out / "coordinates.rst7").read_bytes() == crd.read_bytes()
    assert not (out / "analyze.py").exists()


def test_amber_atom_count_mismatch_stops(sk_root, tmp_path):
    top, crd = amber_files(tmp_path)
    crd.write_text("test\n     4\n", encoding="ascii")
    spec = molecular(AmberMethod(topology_file=str(top), coordinates_file=str(crd), cutoff_ang=9, igb=0),
                     "geometry_optimization")
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_for(sk_root))
    assert any(e.location == "structure.atoms" for e in ex.value.errors)


def test_amber_periodic_restart_needs_matching_box(sk_root, tmp_path):
    top, crd = amber_files(tmp_path)
    spec = molecular(AmberMethod(topology_file=str(top), coordinates_file=str(crd), cutoff_ang=9),
                     "geometry_optimization")
    spec.structure.atoms.cell = [(12, 0, 0), (0, 12, 0), (0, 0, 12)]
    spec.structure.atoms.pbc = (True, True, True)
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_for(sk_root))
    assert any(e.location == "method.coordinates_file" for e in ex.value.errors)
    crd.write_text(crd.read_text(encoding="ascii") + "12 12 12 90 90 90\n", encoding="ascii")
    assert "ntb=1" in build_project(spec, cfg_for(sk_root)).texts["amber.in"]


def namd_files(tmp_path: Path):
    psf, pdb, prm = tmp_path / "system.psf", tmp_path / "system.pdb", tmp_path / "force.prm"
    psf.write_text("PSF\n\n       3 !NATOM\n", encoding="ascii")
    pdb.write_text("ATOM      1  O   HOH A   1       0.000  -1.000   0.000\n"
                   "ATOM      2  H1  HOH A   1       0.000   0.000   0.783\n"
                   "ATOM      3  H2  HOH A   1       0.000   0.000  -0.783\n", encoding="ascii")
    prm.write_text("* test parameter file\n", encoding="ascii")
    return psf, pdb, prm


def test_namd_nve_copies_user_files(sk_root, tmp_path):
    psf, pdb, prm = namd_files(tmp_path)
    method = NamdMethod(structure_file=str(psf), coordinates_file=str(pdb), parameter_files=[str(prm)],
                        exclude="scaled1-4", one_four_scaling=1.0,
                        cutoff_ang=10, pairlistdist_ang=12, switching=True, switchdist_ang=8)
    spec = water_spec(method=method, task=Task(type="molecular_dynamics", md=MDSettings(
        ensemble="NVE", temperature_k=300, timestep_fs=1, steps=100, dump_interval=10)),
        runtime=Runtime(profile="local"))
    files = build_project(spec, cfg_for(sk_root))
    inp = files.texts["namd.conf"]
    assert "numsteps 100" in inp and "DCDfreq 10" in inp
    assert "switchdist 8" in inp and "parameters parameter_01.prm" in inp
    assert files.copies["topology.psf"] == psf and files.copies["coordinates.pdb"] == pdb
    assert files.copies["parameter_01.prm"] == prm
    out = tmp_path / "namd_project"
    write_project(spec, cfg_for(sk_root), out)
    assert (out / "topology.psf").read_bytes() == psf.read_bytes()
    assert (out / "parameter_01.prm").read_bytes() == prm.read_bytes()


def test_namd_requires_explicit_force_field_options(sk_root, tmp_path):
    psf, pdb, prm = namd_files(tmp_path)
    method = NamdMethod(structure_file=str(psf), coordinates_file=str(pdb), parameter_files=[str(prm)],
                        exclude="scaled1-4", cutoff_ang=10, pairlistdist_ang=12)
    spec = water_spec(method=method, task=Task(type="molecular_dynamics", md=MDSettings(ensemble="NVE")),
                      runtime=Runtime(profile="local"))
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_for(sk_root))
    assert {"method.one_four_scaling", "method.switching"} <= {e.location for e in ex.value.errors}


def test_cross_code_report_names_unapplied_settings(tmp_path):
    top, crd = amber_files(tmp_path)
    source = molecular(OrcaMethod(method="HF", basis="STO-3G"), "geometry_optimization")
    target = molecular(AmberMethod(topology_file=str(top), coordinates_file=str(crd),
                                   cutoff_ang=9, igb=0), "geometry_optimization")
    _, report = retarget_spec(source, target)
    assert "task.optimizer" in report["not_applied_by_target"]
    assert "task.force_tolerance_ev_per_ang" not in report["preserved"]
    assert "task.max_steps" in report["preserved"]
    assert "structure.charge" in report["not_applied_by_target"]
    assert "structure" not in report["preserved"]

    nve = water_spec(task=Task(type="molecular_dynamics", md=MDSettings(ensemble="NVE")),
                     runtime=Runtime(profile="local"))
    target_md = nve.model_copy(update={"method": NamdMethod()})
    _, report = retarget_spec(nve, target_md)
    assert "task.md.thermostat" in report["not_applied_by_target"]
    assert "task.md.temperature_k" in report["preserved"]


def test_new_engine_cli_does_not_run(sk_root, tmp_path):
    cfg_path, spec_path = tmp_path / "cluster.toml", tmp_path / "spec.json"
    save_config(cfg_for(sk_root), cfg_path)
    molecular(GaussianMethod(theory="HF", basis="6-31G(d)")).save(spec_path)
    out = tmp_path / "gaussian_project"
    assert generate_main([str(spec_path), str(out), "--config", str(cfg_path)]) == 0
    assert (out / "gaussian.gjf").is_file() and (out / "submit.sh").is_file()
    assert not (out / "output.log").exists()


def test_new_engine_generated_readme_is_clear_in_english(sk_root):
    old_language = lang.LANGUAGE
    try:
        lang.set_language("en")
        files = build_project(molecular(GaussianMethod(theory="HF", basis="6-31G(d)")), cfg_for(sk_root))
        readme = files.texts["README.txt"]
        assert "not editable in the GUI" in readme
        assert "does not yet parse it" in readme
        assert "analyze.py" not in readme
    finally:
        lang.set_language(old_language)


def test_grrm_readme_names_the_log_that_the_run_command_writes(sk_root):
    files = build_project(molecular(GrrmMethod(theory="HF", basis="6-31G"), "geometry_optimization"), cfg_for(sk_root))
    assert files.texts["submit.sh"].rstrip().endswith("grrm > output.log 2>&1")
    assert "output.log" in files.texts["README.txt"] and "grrm.log" not in files.texts["README.txt"]
