
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from ase.build import bulk, molecule
from ase.io import read

from adit.analysis.readers import load_run
from adit.codes.plumed import PLUMED_FILE, check_syntax, engine_has_plumed
from adit.project import ProjectError, build_project, write_project
from adit.spec import (AtomsData, CalculationSpec, DftbMethod, GromacsMethod, LammpsMethod, MDSettings,
                        OpenmmMethod, PlumedSettings, Runtime, Structure, Task)
from tests.conftest import cfg_for

REPO = Path(__file__).resolve().parent.parent
LINES = "d1: DISTANCE ATOMS=1,2\nPRINT ARG=d1 FILE=COLVAR STRIDE=10\n"


def md_task(steps=100, dump=50) -> Task:
    return Task(type="molecular_dynamics",
                md=MDSettings(ensemble="NVT", thermostat="nose_hoover", temperature_k=300.0,
                              timestep_fs=1.0, steps=steps, dump_interval=dump, coupling_time_fs=100.0))


def lammps_spec(tmp_path, plumed=None, task=None) -> CalculationSpec:
    atoms = bulk("Cu", "fcc", a=3.615, cubic=True) * (2, 2, 2)
    potential = REPO / "examples" / "lammps_cu" / "Cu_u3.eam"
    return CalculationSpec(
        structure=Structure(source="bulk", source_ref="Cu fcc", atoms=AtomsData.from_ase(atoms)),
        method=LammpsMethod(units="metal", atom_style="atomic", pair_style="eam",
                            pair_coeff="* * Cu_u3.eam", potential_files=[str(potential)]),
        task=task or md_task(),
        runtime=Runtime(profile="local", mpiprocs=1, omp_threads=1, job_name="cu"),
        plumed=plumed)


def gromacs_spec(plumed=None, task=None) -> CalculationSpec:
    src = REPO / "examples" / "gromacs_spce"
    atoms = read(src / "conf.gro")
    return CalculationSpec(
        structure=Structure(source="file", source_ref=str(src / "conf.gro"), atoms=AtomsData.from_ase(atoms)),
        method=GromacsMethod(topology_file=str(src / "topol.top"), structure_file=str(src / "conf.gro"),
                             coulombtype="PME", rcoulomb_nm=0.85, rvdw_nm=0.85, constraints="h-bonds"),
        task=task or md_task(),
        runtime=Runtime(profile="local", mpiprocs=1, omp_threads=1, job_name="spce"),
        plumed=plumed)


def openmm_spec(plumed=None, task=None) -> CalculationSpec:
    src = REPO / "examples" / "gromacs_spce"
    atoms = read(src / "conf.gro")
    from adit.codes.gromacs import gmx_top_dir
    top_dir = gmx_top_dir()
    return CalculationSpec(
        structure=Structure(source="file", source_ref=str(src / "conf.gro"), atoms=AtomsData.from_ase(atoms)),
        method=OpenmmMethod(input_format="gromacs", topology_file=str(src / "topol.top"),
                            coordinates_file=str(src / "conf.gro"), include_dir=str(top_dir) if top_dir else "",
                            nonbonded_method="PME", nonbonded_cutoff_nm=0.8, constraints="HBonds", platform="CPU"),
        task=task or md_task(),
        runtime=Runtime(profile="local", mpiprocs=1, omp_threads=1, job_name="spce"),
        plumed=plumed)


def errors(spec, cfg) -> list:
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg)
    return ex.value.errors


def places(errs) -> set[str]:
    return {e.location for e in errs}


def test_lammps_gets_fix_plumed(tmp_path, sk_root):
    out = tmp_path / "calc"
    write_project(lammps_spec(tmp_path, PlumedSettings(lines=LINES)), cfg_for(sk_root), out)
    assert (out / PLUMED_FILE).read_text(encoding="utf-8") == LINES
    assert "fix adit_plumed all plumed plumedfile plumed.dat outfile plumed.log" in (out / "in.lammps").read_text(encoding="utf-8")
    readme = (out / "README.txt").read_text(encoding="utf-8")
    assert "plumed.dat" in readme and "nm" in readme


def test_gromacs_mdrun_gets_plumed_option(tmp_path, sk_root):
    out = tmp_path / "calc"
    write_project(gromacs_spec(PlumedSettings(lines=LINES)), cfg_for(sk_root), out)
    assert (out / PLUMED_FILE).is_file()
    assert "-plumed plumed.dat" in (out / "submit.sh").read_text(encoding="utf-8")
    assert "PLUMED_KERNEL" in (out / "README.txt").read_text(encoding="utf-8")


@pytest.mark.skipif(not shutil.which("gmx"), reason="the OpenMM input needs gmx (for the GROMACS force field directory)")
def test_openmm_script_adds_plumed_force(tmp_path, sk_root):
    out = tmp_path / "calc"
    write_project(openmm_spec(PlumedSettings(lines=LINES)), cfg_for(sk_root), out)
    assert (out / PLUMED_FILE).is_file()
    assert '"plumed_file": "plumed.dat"' in (out / "openmm_settings.json").read_text(encoding="utf-8")
    assert "from openmmplumed import PlumedForce" in (out / "run_openmm.py").read_text(encoding="utf-8")


def test_input_file_is_copied(tmp_path, sk_root):
    user = tmp_path / "mine.dat"
    user.write_text(LINES, encoding="utf-8")
    out = tmp_path / "calc"
    write_project(lammps_spec(tmp_path, PlumedSettings(input_file=str(user))), cfg_for(sk_root), out)
    assert (out / PLUMED_FILE).read_text(encoding="utf-8") == LINES


def test_no_plumed_means_no_file(tmp_path, sk_root):
    out = tmp_path / "calc"
    write_project(lammps_spec(tmp_path), cfg_for(sk_root), out)
    assert not (out / PLUMED_FILE).exists()
    assert "plumed" not in (out / "in.lammps").read_text(encoding="utf-8")


def test_round_trip_spec(tmp_path):
    spec = lammps_spec(tmp_path, PlumedSettings(lines=LINES))
    assert CalculationSpec.from_json(spec.to_json()) == spec
    assert CalculationSpec.from_json(lammps_spec(tmp_path).to_json()).plumed is None


def test_requires_exactly_one_input(tmp_path, sk_root):
    cfg = cfg_for(sk_root)
    assert "plumed.input_file" in places(errors(lammps_spec(tmp_path, PlumedSettings()), cfg))
    both = PlumedSettings(input_file=str(tmp_path / "x.dat"), lines=LINES)
    assert "plumed.input_file" in places(errors(lammps_spec(tmp_path, both), cfg))
    missing = PlumedSettings(input_file=str(tmp_path / "missing.dat"))
    assert "plumed.input_file" in places(errors(lammps_spec(tmp_path, missing), cfg))


def test_only_md_and_supported_codes(tmp_path, sk_root):
    cfg = cfg_for(sk_root)
    task = Task(type="geometry_optimization", max_steps=10, optimizer="FIRE")
    assert "task.type" in places(errors(lammps_spec(tmp_path, PlumedSettings(lines=LINES), task), cfg))
    spec = lammps_spec(tmp_path, PlumedSettings(lines=LINES))
    spec = spec.model_copy(update={"method": DftbMethod(sk_set="fake"),
                                   "structure": spec.structure.model_copy(
                                       update={"atoms": AtomsData.from_ase(molecule("H2O"))})})
    assert "plumed" in places(errors(spec, cfg))


@pytest.mark.skipif(not shutil.which("plumed"), reason="plumed is not on PATH")
def test_broken_input_is_caught_by_plumed(tmp_path, sk_root):
    bad = PlumedSettings(lines="d1: NOSUCHACTION ATOMS=1,2\nPRINT ARG=d1 FILE=COLVAR\n")
    assert "plumed.lines" in places(errors(lammps_spec(tmp_path, bad), cfg_for(sk_root)))
    assert check_syntax(LINES, 32) is None


def test_engine_check_says_unknown_without_the_engine():
    assert engine_has_plumed("openmm") is None
    assert engine_has_plumed("lammps") in (True, False, None)


def _plumed_kernel() -> str:
    exe = shutil.which("plumed")
    if not exe:
        return ""
    lib = Path(exe).resolve().parent.parent / "lib" / "libplumedKernel.so"
    return str(lib) if lib.is_file() else ""


LAMMPS_OK = bool(shutil.which("lmp")) and engine_has_plumed("lammps") is True
GROMACS_OK = bool(shutil.which("gmx")) and engine_has_plumed("gromacs") is True and bool(_plumed_kernel())


def _openmm_plumed_python() -> str:
    exe = shutil.which("python3") or ""
    if not exe:
        return ""
    probe = subprocess.run([exe, "-c", "import openmm, openmmplumed"], capture_output=True)
    return exe if probe.returncode == 0 else ""


@pytest.mark.skipif(not LAMMPS_OK, reason="no lmp built with PLUMED")
def test_lammps_real_run_with_plumed(tmp_path, sk_root):
    out = tmp_path / "run"
    write_project(lammps_spec(tmp_path, PlumedSettings(lines=LINES), md_task(steps=50, dump=25)),
                  cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, (out / "output.log").read_text(encoding="utf-8")[-2000:]
    assert (out / "COLVAR").read_text(encoding="utf-8").startswith("#! FIELDS time d1")
    table = load_run(out).extra_tables["plumed"]["files"]["COLVAR"]
    assert table["fields"] == ["time", "d1"] and table["rows"] >= 2
    assert 0.2 < table["first"]["d1"] < 0.3


@pytest.mark.skipif(not GROMACS_OK, reason="no gmx with PLUMED support and libplumedKernel.so")
def test_gromacs_real_run_with_plumed(tmp_path, sk_root):
    out = tmp_path / "run"
    write_project(gromacs_spec(PlumedSettings(lines=LINES), md_task(steps=50, dump=25)), cfg_for(sk_root), out)
    env = {**os.environ, "PLUMED_KERNEL": _plumed_kernel(), "OMP_NUM_THREADS": "1"}
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=900, env=env)
    assert r.returncode == 0, (out / "output.log").read_text(encoding="utf-8")[-2000:]
    table = load_run(out).extra_tables["plumed"]["files"]["COLVAR"]
    assert table["fields"] == ["time", "d1"] and table["rows"] >= 2


@pytest.mark.skipif(not (_openmm_plumed_python() and shutil.which("gmx")),
                    reason="the python3 with openmm-plumed or gmx (for the force field directory) is missing")
def test_openmm_real_run_with_plumed(tmp_path, sk_root):
    out = tmp_path / "run"
    write_project(openmm_spec(PlumedSettings(lines=LINES), md_task(steps=50, dump=25)), cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=900,
                       env={**os.environ, "OPENMM_CPU_THREADS": "1"})
    assert r.returncode == 0, (out / "output.log").read_text(encoding="utf-8")[-2000:]
    table = load_run(out).extra_tables["plumed"]["files"]["COLVAR"]
    assert table["fields"] == ["time", "d1"] and table["rows"] >= 2


def test_prep_origin_keeps_plumed_for_supported_md(tmp_path):
    from adit.web.codefields import PrepOrigin

    spec = lammps_spec(tmp_path, PlumedSettings(lines=LINES))
    origin = PrepOrigin(spec)
    rebuilt = spec.model_copy(update={"plumed": None})
    assert origin.apply(rebuilt).plumed == spec.plumed


def test_prep_origin_drops_plumed_when_it_no_longer_applies(tmp_path):
    from adit.web.codefields import PrepOrigin

    spec = lammps_spec(tmp_path, PlumedSettings(lines=LINES))
    origin = PrepOrigin(spec)
    not_md = spec.model_copy(update={"plumed": None, "task": Task(type="geometry_optimization", max_steps=10)})
    assert origin.apply(not_md).plumed is None
    other_code = spec.model_copy(update={"plumed": None, "method": DftbMethod(sk_set="mio-1-1")})
    assert origin.apply(other_code).plumed is None


def test_conversion_does_not_carry_plumed(tmp_path):
    from adit.convert import retarget_spec

    source = lammps_spec(tmp_path, PlumedSettings(lines=LINES))
    target = gromacs_spec()
    converted, report = retarget_spec(source, target)
    assert converted.plumed is None
    assert "plumed" in report["not_applied_by_target"]
    assert "PLUMED" in report["rule"]
