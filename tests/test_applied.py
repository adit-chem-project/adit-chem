from adit.applied import code_specific_unapplied, unapplied_settings
from adit.spec import MDSettings, OpenmmMethod, Psi4Method, Task, XtbMethod
from tests.conftest import water_spec


def test_psi4_optimizer_and_force_threshold_are_reported_as_unapplied():
    out = code_specific_unapplied(water_spec(method=Psi4Method(), task=Task(type="geometry_optimization")))
    assert {"task.optimizer", "task.force_tolerance_ev_per_ang"} <= set(out)
    assert "geom_maxiter" in out["task.optimizer"]["reason"]


def test_an_nve_run_lists_the_thermostat_and_pressure_fields():
    spec = water_spec(method=XtbMethod(), task=Task(type="molecular_dynamics", md=MDSettings(ensemble="NVE")))
    out = unapplied_settings(spec)
    for path in ("task.md.thermostat", "task.md.coupling_time_fs", "task.md.pressure_bar", "task.md.barostat_time_fs",
                 "task.max_steps", "task.relax_cell"):
        assert path in out, path
    assert out["task.max_steps"]["value"] == spec.task.max_steps


def test_a_force_field_code_cannot_take_charge_or_multiplicity():
    out = unapplied_settings(water_spec(method=OpenmmMethod(), task=Task(type="single_point")))
    assert out["structure.charge"]["value"] == 0 and "structure.multiplicity" in out
    assert "task.optimizer" in out                                   # a single point does not move atoms


def test_a_dftb_optimisation_of_a_molecule_reports_the_cell_and_the_shared_optimizer():
    out = unapplied_settings(water_spec(task=Task(type="geometry_optimization")))
    assert set(out) == {"task.relax_cell", "task.optimizer"}       # DFTB+ uses its own optimiser choice
    assert "dftbplus" in out["task.optimizer"]["reason"]
