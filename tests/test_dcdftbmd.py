"""DCDFTBMD 2.0 official-input subset; no executable is bundled for a live run.

Source: https://www.chem.waseda.ac.jp/dcdftbmd/document/DCDFTBMD_2.0_en.pdf
"""

import pytest

from adit.project import ProjectError, build_project
from adit.spec import CalculationSpec, DcdftbmdMethod, HARTREE_PER_BOHR_IN_EV_PER_ANG, MDSettings, Runtime, Task
from tests.conftest import cfg_for, water_spec


def _method(tmp_path):
    files = {}
    for pair in ("O-O", "O-H", "H-O", "H-H"):
        path = tmp_path / f"{pair}.spl"
        path.write_text("test parameter\n", encoding="utf-8")
        files[pair] = str(path)
    return DcdftbmdMethod(scc=True, divide_and_conquer=False,
                          highest_angular_momentum={"O": 2, "H": 1}, sk_files=files)


def test_dcdftbmd_single_point_matches_manual_sections(sk_root, tmp_path):
    spec = water_spec(method=_method(tmp_path), task=Task(type="single_point"), runtime=Runtime(profile="local"))
    spec.save(tmp_path / "spec.json")
    spec = CalculationSpec.load(tmp_path / "spec.json")
    files = build_project(spec, cfg_for(sk_root))
    inp = files.texts["dftb.inp"]
    assert inp.startswith("SCC=TRUE DC=FALSE\n\nADIT calculation\n\n2\nO 2\nparams/O-O.spl params/O-H.spl\nH 1\nparams/H-O.spl params/H-H.spl\n\n3 0 1\n")
    assert "O0 " in inp and "H0 " in inp
    assert set(files.copies) == {f"params/{pair}.spl" for pair in ("O-O", "O-H", "H-O", "H-H")}
    assert "dftb_serial.00.x > output.log 2>&1" in files.texts["submit.sh"]


def test_dcdftbmd_opt_and_nvt_units(sk_root, tmp_path):
    m = _method(tmp_path)
    opt = build_project(water_spec(method=m, task=Task(type="geometry_optimization", optimizer="FIRE",
                                                         max_steps=12, force_tolerance_ev_per_ang=HARTREE_PER_BOHR_IN_EV_PER_ANG * 0.001),
                                   runtime=Runtime(profile="local")), cfg_for(sk_root)).texts["dftb.inp"]
    assert "OPT=(MAXITER=12 GRADCONV=0.001 OPTTYPE=5)" in opt
    md = build_project(water_spec(method=m, task=Task(type="molecular_dynamics", md=MDSettings(
        ensemble="NVT", thermostat="berendsen", steps=21, timestep_fs=0.5,
        dump_interval=3, temperature_k=310, coupling_time_fs=100)),
                                   runtime=Runtime(profile="local")), cfg_for(sk_root)).texts["dftb.inp"]
    assert "NSTEP=21" in md and "DELTAT=5e-16" in md and "PRINT=3" in md
    assert "NVT=TRUE NVTTYPE=4 BATHTEMP=310 TAUTEMP=1e-13" in md


def test_dcdftbmd_missing_parameters_and_unmapped_thermostat_stop_before_generation(sk_root, tmp_path):
    m = _method(tmp_path)
    m.sk_files.pop("H-O")
    with pytest.raises(ProjectError) as ex:
        build_project(water_spec(method=m), cfg_for(sk_root))
    assert any(e.location == "method.sk_files" for e in ex.value.errors)
    m = _method(tmp_path)
    with pytest.raises(ProjectError) as ex:
        build_project(water_spec(method=m, task=Task(type="molecular_dynamics", md=MDSettings(
            ensemble="NVT", thermostat="nose_hoover"))), cfg_for(sk_root))
    assert any(e.location == "task.md.thermostat" for e in ex.value.errors)


def test_dcdftbmd_writes_no_analysis_script(sk_root, tmp_path):
    spec = water_spec(method=_method(tmp_path), task=Task(type="single_point"), runtime=Runtime(profile="local"))
    files = build_project(spec, cfg_for(sk_root))
    assert "analyze.py" not in files.texts
    readme = files.texts["README.txt"]
    assert "analyze.py" not in readme and "自動解析は未対応" in readme
