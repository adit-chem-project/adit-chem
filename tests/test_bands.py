
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from ase.build import bulk

from adit.bandpath import band_path, dftb_klines, qe_crystal_b, vasp_line_mode
from adit.project import ProjectError, build_project, write_project
from adit.spec import BandSettings, KPoints, OrcaMethod, Task, XtbMethod
from tests.conftest import REAL_SK_ROOT, cfg_for, water_spec

DFTB = os.environ.get("ADIT_DFTB_EXE") or shutil.which("dftb+") or ""
REPO = Path(__file__).resolve().parent.parent


def test_band_path_and_formats():
    kp = band_path(bulk("Si", "diamond", a=5.43), npoints=60)
    assert kp.path.startswith("GX") and len(kp.kpts) == 60 and kp.labels[0] == (0, "G")
    kl = dftb_klines(kp)
    assert kl[0].strip().startswith("KPointsAndWeights = Klines") and kl[1].strip().startswith("1  0.000000 0.000000 0.000000")
    assert sum(int(l.split()[0]) for l in kl[1:-1]) == 60
    vl = vasp_line_mode(kp, 10)
    assert "Line-mode" in vl and "fractional" in vl and vl.count("\n\n") == len(kp.segments())
    qb = qe_crystal_b(kp)
    assert qb.startswith("K_POINTS crystal_b\n") and int(qb.splitlines()[1]) == qb.count("\n") - 2


def test_dftb_band_files(sk_root):
    from tests.test_periodic import tio2_spec
    from tests.conftest import make_fake_skset
    make_fake_skset(sk_root, "ti-0-0", ["Ti", "O"])
    files = build_project(tio2_spec(sk_set="ti-0-0", task=Task(type="band_structure", bands=BandSettings(npoints=30))), cfg_for(sk_root))
    assert "bands/dftb_in.hsd" in files.texts and "bands/kpath.json" in files.texts
    b = files.texts["bands/dftb_in.hsd"]
    for key in ['<<< "../geometry.gen"', "ReadInitialCharges = Yes", "MaxSccIterations = 1", "ConvergentSccOnly = No", 'Prefix = "../skf/"', "Klines"]:
        assert key in b, key
    assert "SupercellFolding" in files.texts["dftb_in.hsd"] and "Driver {}" in files.texts["dftb_in.hsd"]
    assert files.texts["submit.sh"].rstrip().endswith("cp charges.bin bands/ && cd bands && dftb+ > output.log 2>&1")


def test_vasp_and_qe_band_files(sk_root, tmp_path):
    from tests.test_vasp import make_fake_pp, si_spec as vasp_si
    cfg = cfg_for(sk_root); cfg.profiles["cluster"].env["VASP_PP_PATH"] = str(make_fake_pp(tmp_path / "pp"))
    files = build_project(vasp_si(task=Task(type="band_structure", bands=BandSettings(npoints=40))), cfg)
    assert "ICHARG = 11" in files.texts["bands/INCAR"] and "Line-mode" in files.texts["bands/KPOINTS"] and "NSW = 0" in files.texts["INCAR"]
    assert "cp CHGCAR POTCAR bands/" in files.texts["submit.sh"]
    from tests.test_espresso import make_fake_upf, si_spec as qe_si
    cfg2 = cfg_for(sk_root); cfg2.pseudo_root = str(make_fake_upf(tmp_path / "pseudo"))
    files = build_project(qe_si(task=Task(type="band_structure", bands=BandSettings(npoints=40, empty_bands=6))), cfg2)
    b = files.texts["bands/pw.in"]
    assert "calculation      = 'bands'" in b and "K_POINTS crystal_b" in b and "nbnd             = 10" in b
    assert "outdir           = '../tmp'" in b and "verbosity        = 'high'" in b
    assert "cd bands && mpirun -np 2 pw.x -in pw.in" in files.texts["submit.sh"]


def test_band_structure_rejected_for_molecules_and_codes(sk_root):
    with pytest.raises(ProjectError) as ex:
        build_project(water_spec(task=Task(type="band_structure")), cfg_for(sk_root))
    assert any("周期系" in e.message for e in ex.value.errors)
    from tests.test_periodic import tio2_spec
    for m in (XtbMethod(), OrcaMethod()):
        with pytest.raises(ProjectError) as ex:
            build_project(tio2_spec(sk_set="mio-ext", method=m, kpoints=KPoints()), cfg_for(sk_root))
        assert any("バンド計算はありません" in e.message or "分子系だけ" in e.message for e in ex.value.errors)


def test_read_vasp_eigenval(tmp_path):
    from adit.analysis.bands import _vasp_eigenval
    p = tmp_path / "EIGENVAL"
    p.write_text("h\nh\nh\nh\nh\n  8  2  3\n\n 0.0 0.0 0.0 0.5\n 1 -5.0 1.0\n 2 6.0 0.0\n 3 8.0 0.0\n\n 0.5 0.0 0.5 0.5\n 1 -4.0 1.0\n 2 5.0 0.0\n 3 7.0 0.0\n", encoding="utf-8")
    e = _vasp_eigenval(p)
    assert np.allclose(e, [[-5, 6, 8], [-4, 5, 7]])


@pytest.mark.skipif(not (DFTB and (REAL_SK_ROOT / "mio-ext" / "README").is_file()), reason="dftb+ か mio-ext が無い")
def test_dftb_bands_real_run(tmp_path):
    from adit.analysis import AnalysisOptions, run_analysis
    from tests.test_periodic import tio2_spec
    out = tmp_path / "b"
    write_project(tio2_spec(task=Task(type="band_structure", bands=BandSettings(npoints=30))), cfg_for(REAL_SK_ROOT), out)
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=600, env={**os.environ, "PATH": f"{Path(DFTB).parent}:{os.environ.get('PATH', '')}"})
    assert r.returncode == 0, r.stderr
    res = run_analysis(out, AnalysisOptions())
    assert "bands" in res.figures and res.tables["bands"]["n_kpoints"] == 30 and res.tables["bands"]["gap_ev"] is not None
