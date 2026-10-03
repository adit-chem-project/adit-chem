
from pathlib import Path

import pytest
from ase.io import read

from adit.native_import import import_native
from adit.native_verify import verify_generated_bundle
from adit.project import write_project
from adit.spec import AtomsData, CalculationSpec, GromacsMethod, MDSettings, Runtime, Structure, Task
from tests.conftest import cfg_for

REPO = Path(__file__).resolve().parent.parent
SRC = REPO / "examples" / "gromacs_spce"


def spec_for(task=None, **m) -> CalculationSpec:
    atoms = read(SRC / "conf.gro")
    method = {"topology_file": str(SRC / "topol.top"), "structure_file": str(SRC / "conf.gro"),
              "coulombtype": "PME", "rcoulomb_nm": 0.85, "rvdw_nm": 0.85, "constraints": "h-bonds", **m}
    return CalculationSpec(
        structure=Structure(source="file", source_ref=str(SRC / "conf.gro"), atoms=AtomsData.from_ase(atoms)),
        method=GromacsMethod(**method),
        task=task or Task(type="molecular_dynamics",
                          md=MDSettings(ensemble="NVT", thermostat="csvr", temperature_k=300.0, timestep_fs=2.0,
                                        steps=500, dump_interval=50, coupling_time_fs=100.0)),
        runtime=Runtime(profile="local", mpiprocs=1, omp_threads=1, job_name="spce"))


def test_generated_md_bundle_verifies(tmp_path, sk_root):
    out = tmp_path / "nvt"
    write_project(spec_for(), cfg_for(sk_root), out)
    result = verify_generated_bundle(out)
    fields = {x["field"] for x in result.preserved}
    assert not result.mismatched
    assert {"task.type", "task.md.steps", "task.md.timestep_fs", "task.md.temperature_k", "task.md.ensemble",
            "task.md.thermostat", "method.coulombtype", "method.rcoulomb_nm", "method.constraints"} <= fields
    assert {"structure.atoms.symbols", "structure.atoms.positions", "structure.atoms.cell"} <= fields
    assert any("topol.top" in str(x.get("file", "")) for x in result.unverifiable)


def test_units_are_converted_back(tmp_path, sk_root):
    out = tmp_path / "em"
    task = Task(type="geometry_optimization", optimizer="LBFGS", max_steps=200, force_tolerance_ev_per_ang=0.05)
    write_project(spec_for(task, constraints="none"), cfg_for(sk_root), out)   # GROMACS: L-BFGS cannot be combined with constraints
    native = import_native(out, code="gromacs")
    assert native.parsed["task.max_steps"] == 200
    assert native.parsed["task.force_tolerance_ev_per_ang"] == pytest.approx(0.05, rel=1e-6)
    assert native.spec.task.optimizer == "LBFGS"
    result = verify_generated_bundle(out)
    assert not result.mismatched and {"task.max_steps", "task.force_tolerance_ev_per_ang"} <= {x["field"] for x in result.preserved}


def test_npt_barostat_is_read_back(tmp_path, sk_root):
    out = tmp_path / "npt"
    task = Task(type="molecular_dynamics",
                md=MDSettings(ensemble="NPT", thermostat="csvr", temperature_k=300.0, timestep_fs=2.0, steps=100,
                              dump_interval=50, coupling_time_fs=100.0, pressure_bar=1.0, barostat_time_fs=2000.0))
    write_project(spec_for(task, pcoupl="C-rescale", compressibility_per_bar=4.5e-5), cfg_for(sk_root), out)
    native = import_native(out, code="gromacs")
    assert native.parsed["task.md.ensemble"] == "NPT" and native.parsed["method.pcoupl"] == "C-rescale"
    assert native.parsed["task.md.barostat_time_fs"] == pytest.approx(2000.0)
    assert native.parsed["method.compressibility_per_bar"] == pytest.approx(4.5e-5)
    assert not verify_generated_bundle(out).mismatched


def test_a_changed_input_is_caught(tmp_path, sk_root):
    out = tmp_path / "nvt"
    write_project(spec_for(), cfg_for(sk_root), out)
    mdp = out / "grompp.mdp"
    mdp.write_text(mdp.read_text(encoding="utf-8").replace("ref-t              = 300", "ref-t              = 350"), encoding="utf-8")
    result = verify_generated_bundle(out)
    assert [x["field"] for x in result.mismatched] == ["task.md.temperature_k"]
    assert result.mismatched[0]["expected"] == 300.0 and result.mismatched[0]["observed"] == 350.0


def test_unknown_mdp_option_is_reported(tmp_path, sk_root):
    out = tmp_path / "extra"
    write_project(spec_for(extra_mdp={"nstcalcenergy": 50}), cfg_for(sk_root), out)
    native = import_native(out, code="gromacs")
    assert native.spec is None
    assert any("nstcalcenergy" in x["text"] for x in native.unknown)


def test_detects_gromacs_without_being_told(tmp_path, sk_root):
    out = tmp_path / "nvt"
    write_project(spec_for(), cfg_for(sk_root), out)
    assert import_native(out).code == "gromacs"
