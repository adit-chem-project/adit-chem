
import pytest

from adit.config import Profile
from adit.project import ProjectError, build_project
from adit.spec import KPoints, OrcaMethod, Runtime, Structure, Task
from tests.conftest import cfg_for, water_spec


def test_orca_inp_matches_manual_forms(sk_root):
    spec = water_spec(method=OrcaMethod(method="B3LYP", basis="def2-TZVP", scf_convergence="TightSCF", scf_maxiter=200, maxcore_mb=2000,
                                        extra_keywords="D3BJ", extra_blocks="%basis\n  newGTO H \"def2-SVP\" end\nend"),
                      task=Task(type="geometry_optimization", max_steps=50, force_tolerance_ev_per_ang=0.0154266),
                      runtime=Runtime(profile="local", mpiprocs=4, omp_threads=1, job_name="w"))
    st = spec.structure
    spec = spec.model_copy(update={"structure": Structure(**{**st.model_dump(), "fixed_atoms": [0]})})
    files = build_project(spec, cfg_for(sk_root))
    inp = files.texts["orca.inp"]
    lines = inp.splitlines()
    assert lines[1] == "! B3LYP def2-TZVP Opt TightSCF D3BJ"
    for key in ["%pal nprocs 4 end", "%maxcore 2000", "%scf", "   MaxIter 200", "%geom", "   MaxIter 50", "   TolMaxG 3.000e-04",
                "   Constraints", "      { C 0 C }", "* xyz 0 1", "%basis"]:
        assert key in inp, key
    assert lines[-1] == "*" and lines[-4].split()[0] == "O"
    assert files.copies == {}
    assert files.texts["submit.sh"].rstrip().endswith("orca orca.inp > output.log 2>&1")


def test_orca_single_point_defaults(sk_root):
    inp = build_project(water_spec(method=OrcaMethod(), task=Task(type="single_point")), cfg_for(sk_root)).texts["orca.inp"]
    assert "! HF def2-SVP\n" in inp and "Opt" not in inp and "%geom" not in inp and "%pal" not in inp and "%maxcore" not in inp


def test_orca_goat_uses_documented_keyword_and_names_outputs(sk_root):
    spec = water_spec(method=OrcaMethod(method="XTB", basis="", goat=True), task=Task(type="single_point"))
    files = build_project(spec, cfg_for(sk_root))
    assert "! XTB GOAT\n" in files.texts["orca.inp"]
    assert "Opt" not in files.texts["orca.inp"]
    assert "orca.globalminimum.xyz" in files.texts["README.txt"]
    assert "orca.finalensemble.xyz" in files.texts["README.txt"]


def test_orca_goat_rejects_incompatible_task_and_fixed_atoms(sk_root):
    spec = water_spec(method=OrcaMethod(goat=True), task=Task(type="geometry_optimization"))
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_for(sk_root))
    assert any(e.location == "method.goat" for e in ex.value.errors)
    spec = water_spec(method=OrcaMethod(goat=True), task=Task(type="single_point"))
    spec.structure.fixed_atoms = [0]
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_for(sk_root))
    assert any(e.location == "structure.fixed_atoms" for e in ex.value.errors)


def test_orca_docker_copies_guest_and_uses_documented_block(sk_root, tmp_path):
    guest = tmp_path / "guest.xyz"
    guest.write_text("3\n0 1\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\nH 0.0 1.0 0.0\n", encoding="utf-8")
    spec = water_spec(method=OrcaMethod(method="XTB", basis="", docker_guest_file=str(guest)),
                      task=Task(type="single_point"))
    files = build_project(spec, cfg_for(sk_root))
    assert "! XTB\n" in files.texts["orca.inp"]
    assert '%DOCKER GUEST "guest.xyz" END' in files.texts["orca.inp"]
    assert files.copies["guest.xyz"] == guest
    assert "orca.docker.xyz" in files.texts["README.txt"]


def test_orca_docker_rejects_missing_guest_charge_and_basis(sk_root, tmp_path):
    guest = tmp_path / "guest.xyz"
    guest.write_text("1\ncomment\nH 0 0 0\n", encoding="utf-8")
    spec = water_spec(method=OrcaMethod(method="XTB", docker_guest_file=str(guest)), task=Task(type="single_point"))
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_for(sk_root))
    assert any(e.location == "method.docker_guest_file" for e in ex.value.errors)
    assert any(e.location == "method.method" for e in ex.value.errors)
    spec.method.basis = ""
    spec.method.docker_assume_neutral_singlet = True
    files = build_project(spec, cfg_for(sk_root))
    assert "ORCA の既定 (0, 1)" in files.texts["README.txt"]


def test_goat_and_docker_outputs_are_not_mistaken_for_md_trajectories(tmp_path):
    from ase.io import write
    from adit.analysis.readers import _read_orca

    atoms = water_spec().atoms
    (tmp_path / "orca.inp").write_text("! XTB GOAT\n* xyz 0 1\n*\n", encoding="utf-8")
    (tmp_path / "output.log").write_text("FINAL SINGLE POINT ENERGY -1.0\nFINAL SINGLE POINT ENERGY -2.0\n", encoding="utf-8")
    write(tmp_path / "orca.globalminimum.xyz", atoms, format="xyz")
    result = _read_orca(tmp_path)
    assert result.frame_source == "orca.globalminimum.xyz" and len(result.frames) == 1
    assert result.energies_ev == []
    (tmp_path / "orca.inp").write_text('%DOCKER GUEST "guest.xyz" END\n* xyz 0 1\n*\n', encoding="utf-8")
    write(tmp_path / "orca.docker.xyz", [atoms, atoms], format="xyz")
    result = _read_orca(tmp_path)
    assert result.frame_source == "orca.docker.xyz" and len(result.frames) == 1
    assert result.energies_ev == []


def test_orca_rejects_periodic_axes_and_empty(sk_root):
    from tests.test_periodic import tio2_spec
    with pytest.raises(ProjectError) as ex:
        build_project(tio2_spec(sk_set="mio-ext", method=OrcaMethod(), kpoints=KPoints()), cfg_for(sk_root))
    assert any("分子系だけ" in e.message for e in ex.value.errors)
    st = water_spec().structure
    with pytest.raises(ProjectError) as ex:
        build_project(water_spec(structure=Structure(**{**st.model_dump(), "fixed_axes": {"1": (True, True, False)}}), method=OrcaMethod()), cfg_for(sk_root))
    assert ex.value.errors[0].location == "structure.fixed_axes"
    with pytest.raises(ProjectError) as ex:
        build_project(water_spec(method=OrcaMethod(basis="  ")), cfg_for(sk_root))
    assert ex.value.errors[0].location == "method.basis"


def test_orca_command_from_profile(sk_root):
    cfg = cfg_for(sk_root)
    cfg.profiles["local"].commands["orca"] = "/opt/orca/orca"
    s = build_project(water_spec(method=OrcaMethod(), runtime=Runtime(profile="local", mpiprocs=8)), cfg).texts["submit.sh"]
    assert s.rstrip().endswith("/opt/orca/orca orca.inp > output.log 2>&1") and "mpirun" not in s


def test_orca_md_holds_fixed_atoms_with_a_cartesian_constraint(sk_root):
    from adit.spec import MDSettings
    spec = water_spec(method=OrcaMethod(), task=Task(type="molecular_dynamics", md=MDSettings(ensemble="NVT", thermostat="csvr", steps=20)))
    spec.structure.fixed_atoms = [2, 0]
    files = build_project(spec, cfg_for(sk_root))
    inp = files.texts["orca.inp"]
    md = inp[inp.index("%md"):]
    md = md[:md.index("\nend")]
    assert "   Constraint Add Cartesian 0\n   Constraint Add Cartesian 2\n" in md
    assert md.index("Constraint Add") < md.index("Run 20")
    assert "Constraint Add Cartesian" in files.texts["README.txt"]
    spec.structure.fixed_atoms = []
    assert "Constraint" not in build_project(spec, cfg_for(sk_root)).texts["orca.inp"]


@pytest.mark.parametrize("method,task", [
    (OrcaMethod(), Task(type="vibrations")),
    (OrcaMethod(irc=True), Task(type="single_point")),
])
def test_orca_rejects_fixed_atoms_where_they_would_be_ignored(sk_root, method, task):
    spec = water_spec(method=method, task=task)
    spec.structure.fixed_atoms = [0]
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg_for(sk_root))
    assert any(e.location == "structure.fixed_atoms" for e in ex.value.errors)
    spec.structure.fixed_atoms = []
    build_project(spec, cfg_for(sk_root))
