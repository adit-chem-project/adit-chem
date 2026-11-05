
import shutil
import subprocess

import pytest
from ase.build import molecule

from adit.analysis.readers import load_run
from adit.project import ProjectError, build_project, write_project
from adit.results import summarize_run
from adit.spec import AtomsData, CalculationSpec, KPoints, Psi4Method, Runtime, Structure, Task
from tests.conftest import cfg_for


def spec_for(task=None, atoms=None, charge=0, multiplicity=1, kpoints=None, **m) -> CalculationSpec:
    atoms = atoms if atoms is not None else molecule("H2O")
    method = {"method": "scf", "basis": "cc-pVDZ", **m}
    return CalculationSpec(
        structure=Structure(source="preset", source_ref="H2O", atoms=AtomsData.from_ase(atoms),
                            charge=charge, multiplicity=multiplicity),
        method=Psi4Method(**method),
        task=task or Task(type="single_point"),
        kpoints=kpoints,
        runtime=Runtime(profile="local", mpiprocs=1, omp_threads=2, job_name="h2o"))


def errors(spec, cfg) -> list:
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg)
    return ex.value.errors


def places(errs) -> set[str]:
    return {e.location for e in errs}


def test_writes_molecule_basis_and_call(sk_root, tmp_path):
    out = tmp_path / "calc"
    write_project(spec_for(), cfg_for(sk_root), out)
    text = (out / "input.dat").read_text(encoding="utf-8")
    assert "molecule adit {" in text and "\n  0 1\n" in text
    assert "  units angstrom" in text and "  no_reorient" in text and "  no_com" in text and "  symmetry c1" in text
    assert "set basis cc-pVDZ" in text and "energy('scf', return_wfn=True)" in text
    assert "save_xyz_file('final.xyz', True)" in text
    assert "psi4 -i input.dat -o output.log -n 2 > stdout.log 2>&1" in (out / "submit.sh").read_text(encoding="utf-8")
    readme = (out / "README.txt").read_text(encoding="utf-8")
    assert "input.dat" in readme and "results.json" in readme


def test_task_chooses_the_psi4_call(sk_root, tmp_path):
    for kind, call in (("geometry_optimization", "optimize"), ("vibrations", "frequencies")):
        out = tmp_path / kind
        write_project(spec_for(Task(type=kind, max_steps=9)), cfg_for(sk_root), out)
        text = (out / "input.dat").read_text(encoding="utf-8")
        assert f"{call}('scf', return_wfn=True)" in text
        if kind == "geometry_optimization":
            assert "set geom_maxiter 9" in text and 'results["converged"] = True' in text
        else:
            assert "wfn.frequencies().to_array()" in text


def test_charge_multiplicity_reference_and_memory(sk_root, tmp_path):
    out = tmp_path / "calc"
    atoms = molecule("OH")
    write_project(spec_for(atoms=atoms, charge=0, multiplicity=2, reference="uhf", memory_mb=500),
                  cfg_for(sk_root), out)
    text = (out / "input.dat").read_text(encoding="utf-8")
    assert "\n  0 2\n" in text and "set reference uhf" in text and "memory 500 MB" in text


def test_extra_set_lines(sk_root, tmp_path):
    out = tmp_path / "calc"
    write_project(spec_for(extra_set={"scf_type": "df", "e_convergence": "10"}), cfg_for(sk_root), out)
    text = (out / "input.dat").read_text(encoding="utf-8")
    assert "set scf_type df" in text and "set e_convergence 10" in text


def test_round_trip_spec():
    spec = spec_for(reference="rks", memory_mb=1000, extra_set={"scf_type": "pk"})
    assert CalculationSpec.from_json(spec.to_json()) == spec


def test_requires_method_and_basis(sk_root):
    cfg = cfg_for(sk_root)
    assert {"method.method", "method.basis"} <= places(errors(spec_for(method="", basis=""), cfg))
    assert "method.basis" in places(errors(spec_for(basis="cc pVDZ"), cfg))


def test_open_shell_needs_an_open_shell_reference(sk_root):
    cfg = cfg_for(sk_root)
    assert "method.reference" in places(errors(spec_for(atoms=molecule("OH"), multiplicity=2), cfg))
    write = spec_for(atoms=molecule("OH"), multiplicity=2, reference="uhf")
    assert build_project(write, cfg) is not None


def test_electron_parity_is_checked(sk_root):
    assert "structure.multiplicity" in places(errors(spec_for(multiplicity=2, reference="uhf"), cfg_for(sk_root)))


def test_periodic_kpoints_and_unsupported_tasks(sk_root):
    cfg = cfg_for(sk_root)
    atoms = molecule("H2O")
    atoms.cell = [10, 10, 10]
    atoms.pbc = True
    assert "structure.atoms" in places(errors(spec_for(atoms=atoms), cfg))
    assert "kpoints" in places(errors(spec_for(kpoints=KPoints(mode="mesh", mesh=(2, 2, 2))), cfg))
    for kind in ("molecular_dynamics", "band_structure"):
        assert "task.type" in places(errors(spec_for(Task(type=kind)), cfg))


def test_fixed_atoms_are_stopped(sk_root):
    spec = spec_for()
    spec = spec.model_copy(update={"structure": spec.structure.model_copy(update={"fixed_atoms": [0]})})
    assert "structure.fixed_atoms" in places(errors(spec, cfg_for(sk_root)))


def test_extra_set_cannot_override_generated_options(sk_root):
    cfg = cfg_for(sk_root)
    assert "method.extra_set" in places(errors(spec_for(extra_set={"basis": "sto-3g"}), cfg))
    assert "method.extra_set" in places(errors(spec_for(extra_set={"1bad": "x"}), cfg))


def test_mpi_is_rejected(sk_root):
    spec = spec_for()
    spec = spec.model_copy(update={"runtime": spec.runtime.model_copy(update={"mpiprocs": 4})})
    assert "runtime.mpiprocs" in places(errors(spec, cfg_for(sk_root)))


PSI4 = shutil.which("psi4") or ""


def _run(spec, sk_root, out):
    write_project(spec, cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, (out / "stdout.log").read_text(encoding="utf-8")[-3000:]
    return load_run(out)


@pytest.mark.skipif(not PSI4, reason="psi4 is not on PATH")
def test_psi4_real_run_single_point(sk_root, tmp_path):
    run = _run(spec_for(), sk_root, tmp_path / "sp")
    assert run.energies_ev and run.energies_ev[-1] == pytest.approx(-76.026 * 27.211386245988, rel=1e-4)
    assert run.final is not None and len(run.final) == 3
    assert summarize_run(tmp_path / "sp").finished is True
    assert "adit-psi4: psi4" in (tmp_path / "sp" / "code_version.txt").read_text(encoding="utf-8")


@pytest.mark.skipif(not PSI4, reason="psi4 is not on PATH")
def test_psi4_real_run_optimization(sk_root, tmp_path):
    out = tmp_path / "opt"
    run = _run(spec_for(Task(type="geometry_optimization", max_steps=20)), sk_root, out)
    assert len(run.energies_ev) >= 2 and run.energies_ev[-1] < run.energies_ev[0]
    d = run.final.get_distance(0, 1)
    assert 0.93 < d < 0.96
    summary = summarize_run(out)
    assert summary.converged is True and summary.geometry_steps == len(run.energies_ev)


@pytest.mark.skipif(not PSI4, reason="psi4 is not on PATH")
def test_psi4_real_run_vibrations(sk_root, tmp_path):
    out = tmp_path / "vib"
    run = _run(spec_for(Task(type="vibrations")), sk_root, out)
    assert run.frequencies_cm1 and len(run.frequencies_cm1) >= 3
    assert max(run.frequencies_cm1) > 3000
