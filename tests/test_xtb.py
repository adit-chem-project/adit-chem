
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from adit.codes.xtb import XtbGenerator
from adit.project import ProjectError, build_project, write_project
from adit.spec import KPoints, Structure, Task, XtbMethod
from tests.conftest import cfg_for, water_spec, run_log_tail

XTB = os.environ.get("ADIT_XTB_EXE") or shutil.which("xtb") or ""
gen = XtbGenerator()


def test_xtb_files_and_command(sk_root):
    spec = water_spec(method=XtbMethod(gfn="1", accuracy=0.5, etemp=400, max_iterations=100, opt_level="tight"),
                      task=Task(type="geometry_optimization", max_steps=50))
    files = build_project(spec, cfg_for(sk_root))
    assert set(files.texts) >= {"struct.xyz", "xtb.inp", "submit.sh", "README.txt", "spec.json"} and files.copies == {}
    ctl = files.texts["xtb.inp"]
    assert "$scc" in ctl and "maxiterations=100" in ctl and "maxcycle=50" in ctl and ctl.rstrip().endswith("$end")
    cmd = files.texts["submit.sh"].rstrip().splitlines()[-1]
    for key in ["xtb struct.xyz", "--gfn 1", "--chrg 0", "--uhf 0", "--acc 0.5", "--etemp 400", "--input xtb.inp", "--opt tight", "--json"]:
        assert key in cmd, key
    assert files.texts["struct.xyz"].splitlines()[0].strip() == "3"


def test_xtb_gfnff_single_point_and_fix(sk_root):
    st = water_spec().structure
    spec = water_spec(structure=Structure(**{**st.model_dump(), "fixed_atoms": [0], "charge": 1, "multiplicity": 2}),
                      method=XtbMethod(gfn="ff"), task=Task(type="single_point"))
    files = build_project(spec, cfg_for(sk_root))
    cmd = files.texts["submit.sh"].rstrip().splitlines()[-1]
    assert "--gfnff" in cmd and "--chrg 1" in cmd and "--uhf 1" in cmd and "--opt" not in cmd
    assert "$fix" in files.texts["xtb.inp"] and "atoms: 1" in files.texts["xtb.inp"]


def test_xtb_rejects_periodic_and_axes(sk_root):
    from tests.test_periodic import tio2_spec
    with pytest.raises(ProjectError) as ex:
        build_project(tio2_spec(sk_set="mio-ext", method=XtbMethod(), kpoints=KPoints()), cfg_for(sk_root))
    assert any("分子系だけ" in e.message for e in ex.value.errors)
    st = water_spec().structure
    bad = water_spec(structure=Structure(**{**st.model_dump(), "fixed_axes": {"1": (True, True, False)}}), method=XtbMethod())
    with pytest.raises(ProjectError) as ex:
        build_project(bad, cfg_for(sk_root))
    assert ex.value.errors[0].location == "structure.fixed_axes"


@pytest.mark.skipif(not XTB, reason="xtb is missing")
def test_xtb_real_run(sk_root, tmp_path):
    from adit.results import summarize_run
    out = tmp_path / "w"
    write_project(water_spec(method=XtbMethod(), task=Task(type="geometry_optimization", max_steps=100)), cfg_for(sk_root), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=300,
                       env={**os.environ, "PATH": f"{Path(XTB).parent}:{os.environ.get('PATH', '')}"})
    assert r.returncode == 0, f"{r.stderr}\n{run_log_tail(out)}"
    log = (out / "output.log").read_text(encoding="utf-8")
    assert "GEOMETRY OPTIMIZATION CONVERGED" in log and (out / "xtbopt.xyz").is_file()
    s = summarize_run(out)
    assert s.converged and s.mermin_energy_hartree is not None and s.mermin_energy_hartree < -5.0
    assert any(abs(d - 0.96) < 0.03 for _, d in s.bonds)
