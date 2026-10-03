
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from adit.codes.sk_sets import SKSet
from adit.codes.dftbplus import DftbPlusGenerator
from adit.project import ProjectError, build_project, write_project
from adit.spec import DftbMethod, MDSettings, OrcaMethod, Task, VaspMethod, XtbMethod
from tests.conftest import REAL_SK_ROOT, cfg_for, water_spec

DFTB = os.environ.get("ADIT_DFTB_EXE") or shutil.which("dftb+") or ""
XTB = os.environ.get("ADIT_XTB_EXE") or shutil.which("xtb") or ""


def md_task(**kw) -> Task:
    d = dict(ensemble="NVT", thermostat="berendsen", temperature_k=300, timestep_fs=0.5, steps=100, dump_interval=10, coupling_time_fs=100)
    d.update(kw)
    return Task(type="molecular_dynamics", md=MDSettings(**d))


def test_dftb_md_blocks(sk_root):
    skset = SKSet.from_dir(sk_root / "fake-1-0")
    gen = DftbPlusGenerator()
    hsd = gen.dftb_in_hsd(water_spec(task=md_task()), skset)
    for key in ["Driver = VelocityVerlet", "Steps = 100", "TimeStep [fs] = 0.5", "MDRestartFrequency = 10", "Thermostat = Berendsen",
                "Temperature [Kelvin] = 300", "Timescale [fs] = 100"]:
        assert key in hsd, key
    nve = gen.dftb_in_hsd(water_spec(task=md_task(ensemble="NVE")), skset)
    assert "Thermostat = None" in nve and "InitialTemperature [Kelvin] = 300" in nve
    an = gen.dftb_in_hsd(water_spec(task=md_task(thermostat="andersen", timestep_fs=1.0, coupling_time_fs=50)), skset)
    assert "Thermostat = Andersen" in an and "ReselectProbability = 0.02" in an
    nh = gen.dftb_in_hsd(water_spec(task=md_task(thermostat="nose_hoover", coupling_time_fs=100)), skset)
    assert "Thermostat = NoseHoover" in nh and "CouplingStrength [THz] = 10" in nh
    hess = gen.dftb_in_hsd(water_spec(task=Task(type="vibrations")), skset)
    assert "Driver = SecondDerivatives" in hess and "Atoms = 1:-1" in hess
    with pytest.raises(ProjectError) as ex:
        build_project(water_spec(task=md_task(thermostat="langevin")), cfg_for(sk_root))
    assert ex.value.errors[0].location == "task.md.thermostat"
    with pytest.raises(ProjectError) as ex:
        build_project(water_spec(task=md_task(ensemble="NPT")), cfg_for(sk_root))
    assert ex.value.errors[0].location == "task.md.ensemble"


def test_xtb_md_and_hess(sk_root):
    files = build_project(water_spec(method=XtbMethod(), task=md_task(steps=400, timestep_fs=0.5, dump_interval=10)), cfg_for(sk_root))
    ctl, cmd = files.texts["xtb.inp"], files.texts["submit.sh"].rstrip().splitlines()[-1]
    for key in ["$md", "temp=300", "time=0.2", "dump=5", "step=0.5", "nvt=true", "hmass=4", "shake=2", "sccacc=2"]:
        assert key in ctl, key
    assert "--md" in cmd and "--opt" not in cmd
    nve = build_project(water_spec(method=XtbMethod(), task=md_task(ensemble="NVE")), cfg_for(sk_root)).texts["xtb.inp"]
    assert "nvt=false" in nve
    hess = build_project(water_spec(method=XtbMethod(), task=Task(type="vibrations")), cfg_for(sk_root)).texts["submit.sh"]
    assert "--hess" in hess
    with pytest.raises(ProjectError):
        build_project(water_spec(method=XtbMethod(), task=md_task(thermostat="andersen")), cfg_for(sk_root))


def test_vasp_md_and_freq(sk_root, tmp_path):
    from tests.test_vasp import h2o_spec, make_fake_pp
    cfg = cfg_for(sk_root); cfg.profiles["cluster"].env["VASP_PP_PATH"] = str(make_fake_pp(tmp_path / "pp"))
    inc = build_project(h2o_spec(task=md_task(thermostat="andersen", steps=500, timestep_fs=1.0, coupling_time_fs=50)), cfg).texts["INCAR"]
    for key in ["IBRION = 0", "NSW = 500", "POTIM = 1", "TEBEG = 300", "TEEND = 300", "MDALGO = 1", "ANDERSEN_PROB = 0.02", "ISIF = 2"]:
        assert key in inc, key
    lang = build_project(h2o_spec(task=md_task(thermostat="langevin", coupling_time_fs=100)), cfg).texts["INCAR"]
    assert "MDALGO = 3" in lang and "LANGEVIN_GAMMA = 10 10" in lang
    nve = build_project(h2o_spec(task=md_task(ensemble="NVE")), cfg).texts["INCAR"]
    assert "SMASS = -3" in nve and "MDALGO" not in nve
    freq = build_project(h2o_spec(task=Task(type="vibrations")), cfg).texts["INCAR"]
    assert "IBRION = 5" in freq and "NFREE = 2" in freq and "POTIM = 0.015" in freq
    with pytest.raises(ProjectError) as ex:
        build_project(h2o_spec(task=md_task(thermostat="nose_hoover")), cfg)
    assert ex.value.errors[0].location == "method.extra_incar"
    ok = build_project(h2o_spec(task=md_task(thermostat="nose_hoover"), method=VaspMethod(encut=400, extra_incar={"SMASS": 2.0})), cfg).texts["INCAR"]
    assert "MDALGO = 2" in ok and "SMASS = 2" in ok


def test_espresso_md(sk_root, tmp_path):
    from tests.test_espresso import make_fake_upf, si_spec
    cfg = cfg_for(sk_root); cfg.pseudo_root = str(make_fake_upf(tmp_path / "pseudo"))
    inp = build_project(si_spec(task=md_task(thermostat="csvr", timestep_fs=2.0, steps=50, coupling_time_fs=40)), cfg).texts["pw.in"]
    for key in ["calculation      = 'md'", "nstep            = 50", "ion_dynamics     = 'verlet'", "ion_temperature  = 'svr'", "tempw            = 300", "nraise           = 20"]:
        assert key in inp, key
    assert "dt               = 41.3" in inp  # 2 fs = 41.34 Ry a.u.
    vib = build_project(si_spec(task=Task(type="vibrations")), cfg)
    assert "&inputph" in vib.texts["ph.in"] and "0.0 0.0 0.0" in vib.texts["ph.in"] and "asr = 'crystal'" in vib.texts["dynmat.in"]
    assert "ph.x -in ph.in" in vib.texts["submit.sh"] and "dynmat.x < dynmat.in" in vib.texts["submit.sh"]
    npt = build_project(si_spec(task=md_task(ensemble="NPT", thermostat="berendsen", pressure_bar=10)), cfg).texts["pw.in"]
    assert "calculation      = 'vc-md'" in npt and "cell_dynamics    = 'pr'" in npt and "press            = 0.01" in npt and "ion_dynamics     = 'beeman'" in npt


def test_orca_md_and_freq(sk_root):
    inp = build_project(water_spec(method=OrcaMethod(), task=md_task(thermostat="csvr", timestep_fs=0.5, steps=200, dump_interval=2, coupling_time_fs=10)), cfg_for(sk_root)).texts["orca.inp"]
    for key in ["! HF def2-SVP MD", "%md", "   Timestep 0.5_fs", "   Initvel 300_K", "   Thermostat CSVR 300_K Timecon 10_fs",
                '   Dump Position Stride 2 Filename "trajectory.xyz"', "   Run 200"]:
        assert key in inp, key
    freq = build_project(water_spec(method=OrcaMethod(), task=Task(type="vibrations")), cfg_for(sk_root)).texts["orca.inp"]
    assert "! HF def2-SVP Freq" in freq
    with pytest.raises(ProjectError):
        build_project(water_spec(method=OrcaMethod(), task=md_task(thermostat="langevin")), cfg_for(sk_root))


@pytest.mark.skipif(not (DFTB and (REAL_SK_ROOT / "mio-1-1").is_dir()), reason="dftb+ か mio-1-1 が無い")
def test_dftb_md_real_run(tmp_path):
    out = tmp_path / "md"
    write_project(water_spec(method=DftbMethod(sk_set="mio-1-1"), task=md_task(steps=50, dump_interval=5)), cfg_for(REAL_SK_ROOT), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=300, env={**os.environ, "PATH": f"{Path(DFTB).parent}:{os.environ.get('PATH', '')}"})
    assert r.returncode == 0
    assert (out / "md.out").read_text(encoding="utf-8").count("MD step:") == 11 and (out / "geo_end.xyz").is_file()


@pytest.mark.skipif(not XTB, reason="xtb が無い")
def test_xtb_hess_real_run(sk_root, tmp_path):
    out = tmp_path / "h"
    write_project(water_spec(method=XtbMethod(), task=Task(type="vibrations")), cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=300, env={**os.environ, "PATH": f"{Path(XTB).parent}:{os.environ.get('PATH', '')}"})
    assert r.returncode == 0 and (out / "vibspectrum").is_file()
