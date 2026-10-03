
import os
import shutil
import subprocess
from pathlib import Path

import pytest
from ase.io import read

from adit.analysis.readers import load_run
from adit.codes.gromacs import gmx_top_dir
from adit.project import ProjectError, build_project, write_project
from adit.results import summarize_run
from adit.spec import AtomsData, CalculationSpec, MDSettings, OpenmmMethod, Runtime, Structure, Task
from tests.conftest import cfg_for

GRO = """water
    3
    1SOL     OW    1   0.126   0.639   0.322
    1SOL    HW1    2   0.187   0.713   0.290
    1SOL    HW2    3   0.108   0.577   0.248
   3.00000   3.00000   3.00000
"""
PRMTOP = """%VERSION V0001
%FLAG POINTERS
%FORMAT(10I8)
       3       2       0       0       0       0       0       0       0       0
"""
RST7 = """water
    3
  1.2600000  6.3900000  3.2200000  1.8700000  7.1300000  2.9000000
  1.0800000  5.7700000  2.4800000
"""


@pytest.fixture
def src(tmp_path) -> Path:
    d = tmp_path / "src"
    (d / "myff.ff").mkdir(parents=True)
    (d / "myff.ff" / "forcefield.itp").write_text('[ defaults ]\n1 2 yes 0.5 0.5\n', encoding="utf-8")
    (d / "topol.top").write_text('#include "myff.ff/forcefield.itp"\n[ system ]\nx\n[ molecules ]\nSOL 1\n', encoding="utf-8")
    (d / "conf.gro").write_text(GRO, encoding="utf-8")
    (d / "system.prmtop").write_text(PRMTOP, encoding="utf-8")
    (d / "system.rst7").write_text(RST7, encoding="utf-8")
    return d


def spec_for(src, task=None, **m) -> CalculationSpec:
    atoms = read(src / "conf.gro")
    method = {"input_format": "gromacs", "topology_file": str(src / "topol.top"), "coordinates_file": str(src / "conf.gro"),
              "nonbonded_method": "PME", "nonbonded_cutoff_nm": 0.9, "constraints": "HBonds", **m}
    return CalculationSpec(
        structure=Structure(source="file", source_ref=str(src / "conf.gro"), atoms=AtomsData.from_ase(atoms)),
        method=OpenmmMethod(**method),
        task=task or Task(type="single_point"),
        runtime=Runtime(profile="local", mpiprocs=1, omp_threads=1, job_name="w"))


def md_task(ensemble="NVT", thermostat="langevin", **kw) -> Task:
    return Task(type="molecular_dynamics", md=MDSettings(ensemble=ensemble, thermostat=thermostat, temperature_k=300.0,
                                                         timestep_fs=2.0, steps=20, dump_interval=5, **kw))


def errors(spec, cfg) -> list:
    with pytest.raises(ProjectError) as ex:
        build_project(spec, cfg)
    return ex.value.errors


def places(errs) -> set[str]:
    return {e.location for e in errs}


def test_generates_script_settings_and_copies(src, sk_root, tmp_path):
    out = tmp_path / "calc"
    write_project(spec_for(src, md_task()), cfg_for(sk_root), out)
    assert (out / "run_openmm.py").is_file() and (out / "topol.top").is_file() and (out / "conf.gro").is_file()
    assert (out / "myff.ff" / "forcefield.itp").is_file()
    settings = (out / "openmm_settings.json").read_text(encoding="utf-8")
    assert '"nonbonded_method": "PME"' in settings and '"nonbonded_cutoff_nm": 0.9' in settings
    assert '"topology_file": "topol.top"' in settings
    assert '"task": "molecular_dynamics"' in settings and '"steps": 20' in settings
    assert "python3 run_openmm.py > output.log 2>&1" in (out / "submit.sh").read_text(encoding="utf-8")
    readme = (out / "README.txt").read_text(encoding="utf-8")
    assert "run_openmm.py" in readme and "conda install -c conda-forge openmm" in readme


def test_amber_copies_topology_and_restart(src, sk_root, tmp_path):
    out = tmp_path / "calc"
    spec = spec_for(src, input_format="amber", topology_file=str(src / "system.prmtop"),
                    coordinates_file=str(src / "system.rst7"), nonbonded_method="PME", nonbonded_cutoff_nm=0.9)
    write_project(spec, cfg_for(sk_root), out)
    assert (out / "system.prmtop").is_file() and (out / "system.rst7").is_file()
    assert '"input_format": "amber"' in (out / "openmm_settings.json").read_text(encoding="utf-8")


def test_round_trip_spec(src):
    spec = spec_for(src, md_task("NPT", "langevin", pressure_bar=1.0), barostat_interval_steps=50)
    assert CalculationSpec.from_json(spec.to_json()) == spec


def test_requires_format_files_and_nonbonded(src, sk_root):
    cfg = cfg_for(sk_root)
    spec = spec_for(src, input_format="", nonbonded_method="", topology_file=str(src / "missing.top"))
    assert {"method.input_format", "method.nonbonded_method", "method.topology_file"} <= places(errors(spec, cfg))


def test_cutoff_must_be_positive_unless_nocutoff(src, sk_root):
    cfg = cfg_for(sk_root)
    assert "method.nonbonded_cutoff_nm" in places(errors(spec_for(src, nonbonded_cutoff_nm=0.0), cfg))
    spec = spec_for(src, nonbonded_method="NoCutoff", nonbonded_cutoff_nm=0.9)
    assert "method.nonbonded_cutoff_nm" in places(errors(spec, cfg))


def test_periodic_method_needs_periodic_structure(src, sk_root, tmp_path):
    atoms = read(src / "conf.gro")
    atoms.pbc = False
    spec = spec_for(src)
    spec = spec.model_copy(update={"structure": spec.structure.model_copy(update={"atoms": AtomsData.from_ase(atoms)})})
    assert "method.nonbonded_method" in places(errors(spec, cfg_for(sk_root)))


def test_nonperiodic_method_rejected_for_periodic_structure(src, sk_root):
    assert "method.nonbonded_method" in places(errors(spec_for(src, nonbonded_method="NoCutoff"), cfg_for(sk_root)))


def test_atom_count_must_match_coordinates(src, sk_root):
    spec = spec_for(src)
    atoms = read(src / "conf.gro")[:2]
    spec = spec.model_copy(update={"structure": spec.structure.model_copy(update={"atoms": AtomsData.from_ase(atoms)})})
    assert "structure.atoms" in places(errors(spec, cfg_for(sk_root)))


def test_include_dir_required_when_topology_uses_gromacs_library(tmp_path, sk_root):
    d = tmp_path / "lib"
    d.mkdir()
    (d / "conf.gro").write_text(GRO, encoding="utf-8")
    (d / "topol.top").write_text('#include "oplsaa.ff/forcefield.itp"\n[ system ]\nx\n[ molecules ]\nSOL 1\n', encoding="utf-8")
    spec = spec_for(d)
    errs = errors(spec, cfg_for(sk_root))
    if gmx_top_dir() is None or not (gmx_top_dir() / "oplsaa.ff" / "forcefield.itp").is_file():
        pytest.skip("GROMACS の力場の置き場所が無いので、この検査は「未確認」になる")
    assert "method.include_dir" in places(errs)


def test_unsupported_tasks_and_thermostats(src, sk_root):
    cfg = cfg_for(sk_root)
    assert "task.type" in places(errors(spec_for(src, Task(type="vibrations")), cfg))
    assert "task.md.thermostat" in places(errors(spec_for(src, md_task(thermostat="berendsen")), cfg))
    assert "task.md.thermostat" in places(errors(spec_for(src, md_task(thermostat="csvr")), cfg))
    npt = spec_for(src, md_task("NPT", "langevin", pressure_bar=1.0), nonbonded_method="CutoffNonPeriodic")
    assert "method.nonbonded_method" in places(errors(npt, cfg))


def test_fixed_atoms_and_charge_rejected(src, sk_root):
    cfg = cfg_for(sk_root)
    spec = spec_for(src)
    spec = spec.model_copy(update={"structure": spec.structure.model_copy(update={"fixed_atoms": [0], "charge": -1})})
    assert {"structure.fixed_atoms", "structure.charge"} <= places(errors(spec, cfg))


def _openmm_python() -> str:
    exe = shutil.which("python3") or ""
    if not exe:
        return ""
    probe = subprocess.run([exe, "-c", "import openmm"], capture_output=True)
    return exe if probe.returncode == 0 else ""


OPENMM_PY = _openmm_python()
TOP_DIR = gmx_top_dir()
HAVE_SPCE = bool(TOP_DIR and (TOP_DIR / "spc216.gro").is_file() and (TOP_DIR / "oplsaa.ff" / "spce.itp").is_file())


@pytest.mark.skipif(not (OPENMM_PY and HAVE_SPCE), reason="openmm の python3 か GROMACS 同梱の SPC/E が無い")
def test_openmm_real_run_md(sk_root, tmp_path):
    s = tmp_path / "src"
    s.mkdir()
    shutil.copy(TOP_DIR / "spc216.gro", s / "conf.gro")
    (s / "topol.top").write_text('#include "oplsaa.ff/forcefield.itp"\n#include "oplsaa.ff/spce.itp"\n'
                                 "[ system ]\nSPC/E\n[ molecules ]\nSOL 216\n", encoding="utf-8")
    spec = spec_for(s, md_task(), include_dir=str(TOP_DIR), nonbonded_cutoff_nm=0.8, platform="CPU", write_xyz_trajectory=True)
    out = tmp_path / "md"
    write_project(spec, cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=900,
                       env={**os.environ, "OPENMM_CPU_THREADS": "1"})
    assert r.returncode == 0, (out / "output.log").read_text(encoding="utf-8")[-3000:]
    assert (out / "md.log").is_file() and (out / "trajectory.dcd").is_file() and (out / "final.pdb").is_file()
    run = load_run(out)
    assert len(run.energies_ev) == 4 and run.temperatures_k and len(run.temperatures_k) == 4
    assert run.times_fs == pytest.approx([10.0, 20.0, 30.0, 40.0])
    assert run.frame_dt_fs == 10.0
    assert run.frame_source == "trajectory.extxyz" and len(run.frames) == 4
    assert run.final is not None and len(run.final) == 648
    assert run.final.pbc.all() and run.final.cell.lengths()[0] == pytest.approx(18.6206, abs=1e-3)
    table = run.extra_tables["openmm"]
    assert table["ensemble"] == "NVT" and table["platform"] == "CPU" and table["steps"] == 20
    assert 100.0 < table["final_temperature_k"] < 600.0
    summary = summarize_run(out)
    assert summary.finished is True and summary.geometry_steps == 20


@pytest.mark.skipif(not (OPENMM_PY and HAVE_SPCE), reason="openmm の python3 か GROMACS 同梱の SPC/E が無い")
def test_openmm_real_run_md_without_text_trajectory(sk_root, tmp_path):
    s = tmp_path / "src"
    s.mkdir()
    shutil.copy(TOP_DIR / "spc216.gro", s / "conf.gro")
    (s / "topol.top").write_text('#include "oplsaa.ff/forcefield.itp"\n#include "oplsaa.ff/spce.itp"\n'
                                 "[ system ]\nSPC/E\n[ molecules ]\nSOL 216\n", encoding="utf-8")
    spec = spec_for(s, md_task(), include_dir=str(TOP_DIR), nonbonded_cutoff_nm=0.8, platform="CPU")
    out = tmp_path / "md_dcd"
    write_project(spec, cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=900,
                       env={**os.environ, "OPENMM_CPU_THREADS": "1"})
    assert r.returncode == 0, (out / "output.log").read_text(encoding="utf-8")[-3000:]
    assert not (out / "trajectory.extxyz").exists()
    run = load_run(out)
    assert run.frame_source == "final.pdb" and len(run.frames) == 1
    assert any("write_xyz_trajectory" in n for n in run.notes)


@pytest.mark.skipif(not (OPENMM_PY and HAVE_SPCE), reason="openmm の python3 か GROMACS 同梱の SPC/E が無い")
def test_openmm_real_run_minimization(sk_root, tmp_path):
    s = tmp_path / "src"
    s.mkdir()
    shutil.copy(TOP_DIR / "spc216.gro", s / "conf.gro")
    (s / "topol.top").write_text('#include "oplsaa.ff/forcefield.itp"\n#include "oplsaa.ff/spce.itp"\n'
                                 "[ system ]\nSPC/E\n[ molecules ]\nSOL 216\n", encoding="utf-8")
    task = Task(type="geometry_optimization", max_steps=20, force_tolerance_ev_per_ang=0.1)
    spec = spec_for(s, task, include_dir=str(TOP_DIR), nonbonded_cutoff_nm=0.8, platform="CPU")
    out = tmp_path / "em"
    write_project(spec, cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=900,
                       env={**os.environ, "OPENMM_CPU_THREADS": "1"})
    assert r.returncode == 0, (out / "output.log").read_text(encoding="utf-8")[-3000:]
    run = load_run(out)
    table = run.extra_tables["openmm"]
    assert table["potential_energy_kj_per_mol"] < 0
    assert run.energies_ev and run.energies_ev[-1] < 0
    summary = summarize_run(out)
    assert summary.convergence_assessed is False
    assert "未判定" in summary.status_line() or "not assessed" in summary.status_line()


def test_conversion_from_openmm_is_stopped(src):
    from adit.convert import ConversionError, retarget_spec
    from adit.spec import DftbMethod

    source = spec_for(src, md_task())
    target = source.model_copy(update={"method": DftbMethod(sk_set="mio-1-1")})
    with pytest.raises(ConversionError) as ex:
        retarget_spec(source, target)
    assert "OpenMM" in str(ex.value)


def test_conversion_to_openmm_lists_unapplied_settings(src):
    from adit.convert import retarget_spec
    from adit.spec import DftbMethod

    source = spec_for(src, md_task("NPT", "langevin", pressure_bar=1.0)).model_copy(update={"method": DftbMethod(sk_set="mio-1-1")})
    target = spec_for(src, md_task("NPT", "langevin", pressure_bar=1.0))
    _, report = retarget_spec(source, target)
    assert "task.md.barostat_time_fs" in report["not_applied_by_target"]
    assert {"structure.charge", "structure.multiplicity"} <= set(report["not_applied_by_target"])
    assert "task.md.pressure_bar" in report["preserved"]


def test_include_dir_is_written_as_an_absolute_path(src, tmp_path, monkeypatch):
    from adit.codes.openmm import OpenmmGenerator

    (tmp_path / "ff_top").mkdir()
    monkeypatch.chdir(tmp_path)
    assert OpenmmGenerator().settings(spec_for(src, include_dir="ff_top"))["include_dir"] == str((tmp_path / "ff_top").resolve())
    assert OpenmmGenerator().settings(spec_for(src))["include_dir"] == ""
