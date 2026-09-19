
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from ase.build import bulk

from adit.analysis.readers import load_run
from adit.project import ProjectError, build_project, write_project
from adit.results import summarize_run
from adit.spec import AbinitMethod, AtomsData, CalculationSpec, KPoints, Runtime, Structure, Task
from tests.conftest import cfg_for


@pytest.fixture
def pseudo(tmp_path) -> Path:
    p = tmp_path / "fake" / "14si.pspnc"
    p.parent.mkdir(parents=True)
    p.write_text("fake pseudopotential for tests\n", encoding="utf-8")
    return p


def spec_for(pseudo, task=None, atoms=None, kpoints=None, **m) -> CalculationSpec:
    atoms = atoms if atoms is not None else bulk("Si", "diamond", a=5.43)
    method = {"pseudos": {"Si": str(pseudo)}, "ecut_ha": 8.0, "tolerance": "toldfe", "tolerance_value": 1e-8, **m}
    return CalculationSpec(
        structure=Structure(source="bulk", source_ref="Si diamond", atoms=AtomsData.from_ase(atoms)),
        method=AbinitMethod(**method),
        task=task or Task(type="single_point"),
        kpoints=kpoints if kpoints is not None else KPoints(mode="mesh", mesh=(4, 4, 4)),
        runtime=Runtime(profile="local", mpiprocs=1, omp_threads=1, job_name="si"))


def errors(spec, cfg) -> list:
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg)
    return ex.value.errors


def places(errs) -> set[str]:
    return {e.location for e in errs}


def variables(text: str) -> dict[str, list[str]]:
    out: dict[str, list[str]] = {}
    current = None
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        words = line.split()
        if words[0][0].isalpha():
            current = words[0]
            out[current] = words[1:]
        elif current is not None:
            out[current] += words
    return out


def test_writes_structure_pseudos_and_settings(pseudo, sk_root, tmp_path):
    out = tmp_path / "calc"
    write_project(spec_for(pseudo), cfg_for(sk_root), out)
    v = variables((out / "input.abi").read_text(encoding="utf-8"))
    assert v["natom"] == ["2"] and v["ntypat"] == ["1"] and v["znucl"] == ["14"] and v["typat"] == ["1", "1"]
    assert v["acell"] == ["1.0", "1.0", "1.0"]
    assert [float(x) for x in v["xred"]] == pytest.approx([0, 0, 0, 0.25, 0.25, 0.25])
    assert float(v["rprim"][1]) == pytest.approx(5.43 / 2 / 0.529177210903, rel=1e-9)
    assert v["ecut"] == ["8"] and v["toldfe"] == ["1e-08"] and v["nstep"] == ["30"]
    assert v["ngkpt"] == ["4", "4", "4"] and v["shiftk"] == ["0", "0", "0"]
    assert v["pseudos"] == ['"14si.pspnc"'] and v["pp_dirpath"] == ['"./"']
    assert (out / "14si.pspnc").is_file()
    assert "abinit input.abi > output.log 2>&1" in (out / "submit.sh").read_text(encoding="utf-8")
    readme = (out / "README.txt").read_text(encoding="utf-8")
    assert "input.abi" in readme and "14si.pspnc" in readme


def test_geometry_optimization_writes_ionmov(pseudo, sk_root, tmp_path):
    out = tmp_path / "calc"
    task = Task(type="geometry_optimization", max_steps=7, force_tolerance_ev_per_ang=0.05)
    write_project(spec_for(pseudo, task, tolerance="toldff", tolerance_value=5e-5), cfg_for(sk_root), out)
    v = variables((out / "input.abi").read_text(encoding="utf-8"))
    assert v["ionmov"] == ["2"] and v["ntime"] == ["7"]
    assert float(v["tolmxf"][0]) == pytest.approx(0.05 * 0.529177210903 / 27.211386245988)  # eV/Å → Ha/Bohr


def test_fixed_atoms_and_axes_become_iatfix(pseudo, sk_root, tmp_path):
    atoms = bulk("Si", "diamond", a=5.43, cubic=True)
    spec = spec_for(pseudo, Task(type="geometry_optimization", max_steps=5), atoms=atoms,
                    tolerance="toldff", tolerance_value=5e-5)
    spec = spec.model_copy(update={"structure": spec.structure.model_copy(
        update={"fixed_atoms": [0], "fixed_axes": {"2": (True, True, False)}})})
    out = tmp_path / "calc"
    write_project(spec, cfg_for(sk_root), out)
    v = variables((out / "input.abi").read_text(encoding="utf-8"))
    assert v["natfix"] == ["1"] and v["iatfix"] == ["1"]
    assert v["natfixz"] == ["1"] and v["iatfixz"] == ["3"]
    assert "natfixx" not in v and "natfixy" not in v


def test_cartesian_axes_are_converted_to_lattice_directions(pseudo, sk_root):
    # Cell with a2 rotated in the xy plane: fixing Cartesian z is still the third lattice direction,
    # fixing Cartesian x is not representable by iatfixx/y/z (which act on reduced coordinates).
    from ase import Atoms
    atoms = Atoms("Si2", positions=[(0, 0, 0), (1.3, 1.3, 1.3)],
                  cell=[(5.0, 0, 0), (2.5, 4.33, 0), (0, 0, 5.0)], pbc=True)
    task = Task(type="geometry_optimization", max_steps=5)
    spec = spec_for(pseudo, task, atoms=atoms, tolerance="toldff", tolerance_value=5e-5)
    ok = spec.model_copy(update={"structure": spec.structure.model_copy(update={"fixed_axes": {"1": (True, True, False)}})})
    v = variables(build_project(ok, cfg_for(sk_root)).texts["input.abi"])
    assert v["iatfixz"] == ["2"] and "natfixx" not in v and "natfixy" not in v
    bad = spec.model_copy(update={"structure": spec.structure.model_copy(update={"fixed_axes": {"1": (False, True, True)}})})
    assert "structure.fixed_axes" in places(errors(bad, cfg_for(sk_root)))


def test_extra_variables_are_written_and_checked(pseudo, sk_root, tmp_path):
    out = tmp_path / "calc"
    write_project(spec_for(pseudo, extra={"chksymtnons": "0"}), cfg_for(sk_root), out)
    assert variables((out / "input.abi").read_text(encoding="utf-8"))["chksymtnons"] == ["0"]
    assert "method.extra" in places(errors(spec_for(pseudo, extra={"ecut": "20"}), cfg_for(sk_root)))
    assert "method.extra" in places(errors(spec_for(pseudo, extra={"1bad": "0"}), cfg_for(sk_root)))


def test_round_trip_spec(pseudo):
    spec = spec_for(pseudo, occopt=3, tsmear_ha=0.01, ixc=11, nband=8)
    assert CalculationSpec.from_json(spec.to_json()) == spec


def test_requires_pseudos_ecut_and_tolerance(pseudo, sk_root):
    cfg = cfg_for(sk_root)
    spec = spec_for(pseudo, pseudos={}, ecut_ha=0.0, tolerance="")
    assert {"method.pseudos", "method.ecut_ha", "method.tolerance"} <= places(errors(spec, cfg))
    assert "method.pseudos" in places(errors(spec_for(pseudo, pseudos={"Si": "/no/such/file.psp8"}), cfg))


def test_molecular_and_charged_systems_are_stopped(pseudo, sk_root):
    cfg = cfg_for(sk_root)
    atoms = bulk("Si", "diamond", a=5.43)
    atoms.pbc = False
    assert "structure.atoms" in places(errors(spec_for(pseudo, atoms=atoms), cfg))
    spec = spec_for(pseudo)
    spec = spec.model_copy(update={"structure": spec.structure.model_copy(update={"charge": 1, "multiplicity": 2})})
    assert {"structure.charge", "structure.multiplicity"} <= places(errors(spec, cfg))


def test_unsupported_tasks_and_cell_relaxation(pseudo, sk_root):
    cfg = cfg_for(sk_root)
    for kind in ("molecular_dynamics", "vibrations", "band_structure"):
        assert "task.type" in places(errors(spec_for(pseudo, Task(type=kind)), cfg))
    task = Task(type="geometry_optimization", relax_cell="shape_and_volume", max_steps=5)
    assert "task.relax_cell" in places(errors(spec_for(pseudo, task, tolerance="toldff", tolerance_value=5e-5), cfg))


def test_toldfe_is_rejected_for_optimization(pseudo, sk_root):
    task = Task(type="geometry_optimization", max_steps=5)
    assert "method.tolerance" in places(errors(spec_for(pseudo, task), cfg_for(sk_root)))


ABINIT = shutil.which("abinit") or ""
PSEUDO_DIR = Path(os.environ.get("ADIT_ABINIT_PSEUDOS", ""))
SI_PSEUDO = PSEUDO_DIR / "14si.pspnc" if str(PSEUDO_DIR) else Path("/nonexistent")


@pytest.mark.skipif(not (ABINIT and SI_PSEUDO.is_file()),
                    reason="abinit か ADIT_ABINIT_PSEUDOS/14si.pspnc が無い")
def test_abinit_real_run_single_point(sk_root, tmp_path):
    out = tmp_path / "si"
    write_project(spec_for(SI_PSEUDO), cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, (out / "output.log").read_text(encoding="utf-8")[-3000:]
    run = load_run(out)
    assert len(run.energies_ev) == 1 and -260 < run.energies_ev[0] < -220
    assert run.final is not None and len(run.final) == 2 and run.final.pbc.all()
    assert summarize_run(out).finished is True
    assert "ABINIT" in (out / "code_version.txt").read_text(encoding="utf-8")


@pytest.mark.skipif(not (ABINIT and SI_PSEUDO.is_file()),
                    reason="abinit か ADIT_ABINIT_PSEUDOS/14si.pspnc が無い")
def test_abinit_real_run_relaxation(sk_root, tmp_path):
    atoms = bulk("Si", "diamond", a=5.43)
    atoms.positions[1] += [0.1, 0.0, 0.0]
    task = Task(type="geometry_optimization", max_steps=5, force_tolerance_ev_per_ang=0.05)
    spec = spec_for(SI_PSEUDO, task, atoms=atoms, tolerance="toldff", tolerance_value=5e-5,
                    extra={"chksymtnons": "0"})
    out = tmp_path / "si_opt"
    write_project(spec, cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=900)
    assert r.returncode == 0, (out / "output.log").read_text(encoding="utf-8")[-3000:]
    run = load_run(out)
    assert len(run.energies_ev) >= 2 and run.energies_ev[-1] < run.energies_ev[0]
    assert run.final is not None
    d = run.final.get_distance(0, 1, mic=True)
    assert 2.2 < d < 2.5
    summary = summarize_run(out)
    assert summary.converged is True and summary.geometry_steps == len(run.energies_ev)
