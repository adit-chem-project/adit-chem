
import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk

from adit.analysis import AnalysisOptions, run_analysis
from adit.analysis import compute
from adit.analysis.readers import frequencies_from_hessian, load_run

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"


@pytest.mark.parametrize("language", ["ja", "en"])
@pytest.mark.parametrize("error", [None, 2e-6])
def test_diffusion_summary_does_not_attach_block_sem_to_full_run_d(language, error, monkeypatch):
    from adit import lang
    from adit.analysis.report import AnalysisResult

    monkeypatch.setattr(lang, "LANGUAGE", language)
    diffusion = {"D_cm2_s": 1e-5, "D_err_cm2_s": error}
    block_error = ({"d_err_cm2_s": error, "n_blocks": 5, "block_frames": 20,
                    "block_frame_counts": [20] * 5, "fit_range_fs": [10, 50]} if error is not None else None)
    result = AnalysisResult(code="xtb", run_dir="run", tables={"msd": {
        **diffusion, "species": None, "n_frames": 100, "last_A2": 0.5,
        "fit_range_fs": [10, 50], "fit_range_user": True,
        "by_element": {"H": diffusion, "O": diffusion},
        "D_error": block_error,
    }})
    text = result.summary_text()
    assert "+-" not in text
    assert "1.000e-05 cm^2/s" in text or "1.000e-05 cm²/s" in text
    if error is None:
        assert "±" not in text
    else:
        assert "1.000e-05 ± 2.0e-06 cm^2/s" not in text
        assert ("ブロック D 平均の標準誤差" if language == "ja" else "standard error of the mean block D") in text


def test_read_dftb_md():
    r = load_run(EX / "dftb_md_water_generated")
    assert r.code == "dftbplus" and len(r.energies_ev) == 21 and len(r.temperatures_k) == 21
    assert r.times_fs is not None and abs(r.times_fs[1] - 5.0) < 1e-6
    assert len(r.frames) == 21 and len(r.frames[0]) == 3
    assert r.eigenvalues_ev is not None and r.fermi_ev is not None


def test_read_dftb_opt_and_vib(tmp_path):
    r = load_run(EX / "water_generated")
    assert r.code == "dftbplus" and len(r.energies_ev) == 10 and r.final is not None
    h = (EX / "dftb_md_water_generated").parent / "xtb_vib_water_generated"  # dummy to keep path usage
    assert h.is_dir()


def test_frequencies_from_hessian_water_matches_modes():
    hess = (REPO / "tests" / "data" / "water_hessian.out")
    if not hess.is_file():
        pytest.skip("tests/data/water_hessian.out が無い")
    atoms = Atoms("OHH", positions=[(0, -1, 0), (0, 0, 0.783064), (0, 0, -0.783064)])
    f = frequencies_from_hessian(np.array(hess.read_text(encoding="utf-8").split(), dtype=float), atoms)
    big = sorted(x for x in f if x > 100)
    ref = [687.59, 830.43, 892.59, 1208.76, 1247.43, 1797.61]
    assert len(big) == 6 and all(abs(a - b) < 5 for a, b in zip(big, ref))


def test_read_xtb_vib():
    r = load_run(EX / "xtb_vib_water_generated")
    assert r.code == "xtb" and r.frequencies_cm1 and len(r.frequencies_cm1) == 9
    assert any(abs(f - 1656.71) < 0.1 for f in r.frequencies_cm1) and r.ir_intensities and max(r.ir_intensities) > 50


def test_read_xtb_optimization_energy_per_cycle():
    r = load_run(EX / "xtb_water_generated")
    hartree_ev = 27.211386245988
    assert len(r.energies_ev) == 9
    assert abs(r.energies_ev[-1] - (-5.0705444) * hartree_ev) < 1e-3
    assert abs(r.energies_ev[0] - (-4.9778649) * hartree_ev) < 1e-3
    assert r.final is not None


def test_read_espresso_md():
    r = load_run(EX / "qe_md_si_generated")
    assert r.code == "espresso" and len(r.energies_ev) == 20 and r.temperatures_k and len(r.temperatures_k) == 20
    assert r.times_fs is not None and abs(r.times_fs[1] - 2.0) < 1e-3


def test_run_analysis_writes_figures(tmp_path):
    import shutil
    d = tmp_path / "run"
    shutil.copytree(EX / "dftb_md_water_generated", d)
    res = run_analysis(d, AnalysisOptions(rdf=True, msd=True, dos=True, skip_frames=5))
    assert set(res.figures) >= {"energy", "temperature", "rdf", "msd", "dos"}
    assert all(Path(p).stat().st_size > 1000 for p in res.figures.values())
    assert (d / "analysis" / "summary.json").is_file() and "温度" in res.summary_text()
    assert res.tables["temperature"]["skipped"] == 5
    assert "unwrap_check" in res.tables["msd"]
    assert (d / "analysis" / "rdf.json").is_file()


def test_rdf_rmax_is_capped_at_half_cell_width(tmp_path):
    import shutil
    d = tmp_path / "run"
    shutil.copytree(EX / "qe_md_si_generated", d)
    frames = load_run(d).frames
    assert frames and any(frames[0].pbc), "周期系の軌跡が読めない (このテストが比べる相手が無い)"
    cell = np.asarray(frames[0].cell, dtype=float); vol = abs(np.linalg.det(cell))
    half = 0.5 * min(vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3))
    res = run_analysis(d, AnalysisOptions(rdf=True, rdf_rmax=8.0))
    t = res.tables["rdf"]
    assert t["rmax_requested"] == 8.0 and abs(t["rmax"] - min(8.0, half)) < 1e-9 and t["rmax"] < 8.0


def test_rdf_and_msd_synthetic():
    a = bulk("Si", "diamond", a=5.43, cubic=True)
    frames = [a.copy() for _ in range(5)]
    r, g = compute.rdf(frames, ("Si", "Si"), rmax=2.7, nbins=54)
    peak = r[np.argmax(g)]
    assert abs(peak - 2.35) < 0.1
    moving = []
    for k in range(11):
        b = Atoms("H", positions=[(0.1 * k, 0, 0)], cell=[10, 10, 10], pbc=True)
        moving.append(b)
    t, m, D = compute.msd(moving, None, [10.0 * k for k in range(11)])
    assert abs(m[-1] - 1.0) < 1e-9 and D is not None and D > 0
    wrap = [Atoms("H", positions=[(9.9, 0, 0)], cell=[10, 10, 10], pbc=True), Atoms("H", positions=[(0.1, 0, 0)], cell=[10, 10, 10], pbc=True)]
    _, m2, _ = compute.msd(wrap, None, None)
    assert abs(m2[-1] - 0.04) < 1e-9


def test_analyze_script_runs(tmp_path):
    import shutil, subprocess, sys
    d = tmp_path / "run"
    shutil.copytree(EX / "xtb_vib_water_generated", d)
    r = subprocess.run([sys.executable, "analyze.py"], cwd=d, capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr
    assert "振動数" in r.stdout and (d / "analysis" / "spectrum.png").is_file()


def test_read_vasp_real_outputs():
    r = load_run(EX / "vasp_h2o_generated")
    assert r.code == "vasp" and len(r.energies_ev) == 5 and abs(r.energies_ev[-1] - (-14.222575)) < 1e-4
    assert r.final is not None and r.eigenvalues_ev is not None and r.fermi_ev is not None
    assert len(r.frames) == 5


def test_summary_json_paths_are_relative_to_the_run_directory(tmp_path):
    import shutil
    d = tmp_path / "run"
    shutil.copytree(EX / "dftb_md_water_generated", d)
    res = run_analysis(d, AnalysisOptions(skip_frames=5))
    assert all(Path(p).is_absolute() for p in res.figures.values())
    js = json.loads((d / "analysis" / "summary.json").read_text(encoding="utf-8"))
    assert js["figures"] and all(not Path(p).is_absolute() and (d / p).is_file() for p in js["figures"].values())
    assert js["figures"]["energy"] == "analysis/energy.png"
    charge_file = js["tables"]["charges"][0]["file"]
    assert charge_file == "analysis/charges_Mulliken.csv" and (d / charge_file).is_file()
    assert str(d) not in json.dumps(js["figures"]) + json.dumps(js["tables"])
    elsewhere = tmp_path / "elsewhere"
    run_analysis(d, AnalysisOptions(skip_frames=5, out_dir=elsewhere))
    js = json.loads((elsewhere / "summary.json").read_text(encoding="utf-8"))
    assert js["figures"]["energy"] == "energy.png" and (elsewhere / "energy.png").is_file()


def test_shipped_example_summary_has_no_developer_paths():
    for p in EX.glob("*/analysis/summary.json"):
        assert "/home/" not in p.read_text(encoding="utf-8"), p
