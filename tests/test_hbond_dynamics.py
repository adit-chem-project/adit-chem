import csv
import json
import math
import shutil
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms

from adit.analysis import AnalysisOptions, run_analysis
from adit.analysis import hbond as HB

REPO = Path(__file__).resolve().parent.parent
PATTERN = [1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1]


def _dimer(present: bool, cell=None) -> Atoms:
    # water O-H ... O: donor O at origin, its H along +x, acceptor O 2.8 Å away (or 10 Å when the bond is absent)
    d = 2.8 if present else 10.0
    atoms = Atoms("OHHOHH", positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0], [d, 0, 0], [d + 0.3, 0.9, 0], [d + 0.3, -0.9, 0]])
    if cell is not None:
        atoms.set_cell(cell); atoms.pbc = True
    return atoms


def _expected(pattern, tau_max):
    h = np.array(pattern, dtype=float)
    n = len(h)
    ci, cc = [], []
    for tau in range(tau_max + 1):
        origins = [t0 for t0 in range(n - tau) if h[t0] > 0]
        ci.append(np.mean([h[t0 + tau] for t0 in origins]))
        cc.append(np.mean([float(all(h[t0:t0 + tau + 1] > 0)) for t0 in origins]))
    return np.array(ci), np.array(cc)


def test_count_frame_matches_the_geometry_and_the_periodic_image():
    got = HB.count_frame(_dimer(True), 3.5, 150.0)
    assert got["count"] == 1 and got["pairs"][0]["donor"] == 1 and got["pairs"][0]["hydrogen"] == 2 and got["pairs"][0]["acceptor"] == 4
    assert got["pairs"][0]["distance_A"] == pytest.approx(2.8) and got["pairs"][0]["angle_deg"] == pytest.approx(180.0)
    assert HB.count_frame(_dimer(False), 3.5, 150.0)["count"] == 0
    # 10 Å apart in a 12 Å box: the minimum image is 2 Å away on the other side, but the angle is then 0 deg
    far = HB.count_frame(_dimer(False, cell=np.diag([12.0, 12.0, 12.0])), 3.5, 150.0)
    assert far["count"] == 0
    tri = HB._geometry(_dimer(False, cell=np.diag([12.0, 12.0, 12.0])), HB.DEFAULT_DONORS, "H", 1.3, 3.5)
    assert len(tri[3]) >= 1 and tri[3].min() == pytest.approx(2.0)
    with pytest.raises(HB.HydrogenBondError):
        HB.count_frame(_dimer(True), 0.0, 150.0)


def test_lifetime_matches_the_brute_force_definition():
    frames = [_dimer(bool(p)) for p in PATTERN]
    t = HB.lifetime(frames, 3.5, 150.0, dt_fs=2.0)
    assert t["n_frames"] == len(PATTERN) and t["n_pairs"] == 1 and t["tau_max"] == len(PATTERN) // 2 and t["unit"] == "fs"
    ci, cc = _expected(PATTERN, t["tau_max"])
    assert t["intermittent"] == pytest.approx(ci, abs=1e-9)
    assert t["continuous"] == pytest.approx(cc, abs=1e-9)
    assert t["lag"] == pytest.approx([2.0 * k for k in range(t["tau_max"] + 1)])
    assert t["intermittent"][0] == 1.0 and t["continuous"][0] == 1.0
    assert all(c <= i + 1e-12 for c, i in zip(t["continuous"], t["intermittent"]))
    trap = getattr(np, "trapezoid", None) or np.trapz
    assert t["lifetime_continuous"]["integral"] == pytest.approx(float(trap(cc, np.arange(len(cc)) * 2.0)))
    assert t["lifetime_continuous"]["reached_one_over_e"] and t["lifetime_continuous"]["one_over_e"] > 0
    assert "MDAnalysis" in t["definition"] and t["source"].startswith("https://docs.mdanalysis.org")
    assert t["memory"]["actual_mb"] == pytest.approx(len(PATTERN) / 2 ** 20)
    t2 = HB.lifetime(frames, 3.5, 150.0, tau_max=3)
    assert t2["unit"] == "frame" and len(t2["lag"]) == 4 and t2["lag"] == [0.0, 1.0, 2.0, 3.0]


def test_lifetime_with_two_pairs_and_the_memory_guard():
    # two independent dimers: pair A present always, pair B follows PATTERN; the ratio averages per origin
    frames = []
    for p in PATTERN:
        a = _dimer(True)
        b = _dimer(bool(p)); b.translate([0, 20.0, 0])
        frames.append(a + b)
    t = HB.lifetime(frames, 3.5, 150.0)
    assert t["n_pairs"] == 2
    h = np.array(PATTERN, dtype=float)
    n = len(h)
    for tau in range(t["tau_max"] + 1):
        exp_i = np.mean([(1 + h[t0] * h[t0 + tau]) / (1 + h[t0]) for t0 in range(n - tau)])
        exp_c = np.mean([(1 + float(all(h[t0:t0 + tau + 1] > 0)) * h[t0]) / (1 + h[t0]) for t0 in range(n - tau)])
        assert t["intermittent"][tau] == pytest.approx(exp_i, abs=1e-9)
        assert t["continuous"][tau] == pytest.approx(exp_c, abs=1e-9)
    with pytest.raises(HB.HydrogenBondError) as ex:
        HB.lifetime(frames, 3.5, 150.0, budget_mb=1e-6)
    assert "--stride" in str(ex.value)
    with pytest.raises(HB.HydrogenBondError):
        HB.lifetime([_dimer(False), _dimer(False)], 3.5, 150.0)


def test_distance_angle_map_counts_every_triplet():
    frames = [_dimer(True)] * 3
    t = HB.distance_angle_map(frames, 4.0, bins=(40, 36))
    d = t["distribution"]
    # 4 hydrogens, each bonded to its nearest O with the other O within 4 Å: 4 triplets per frame, any angle
    assert t["n_frames"] == 3 and t["n_points"] == 12 and d.counts.sum() == 12
    area = np.outer(np.diff(np.linspace(0, 4.0, 41)), np.diff(np.linspace(0, 180, 37)))
    assert float((d.density * area).sum()) == pytest.approx(1.0)
    # the donor hydrogen (2.8 Å, 180 deg) lands in the last angle bin, once per frame
    assert d.counts[26:30, -1].sum() == 3 and d.counts[:, -1].sum() == 3
    assert t["peak"]["distance_A"] <= 4.0 and 0 <= t["peak"]["angle_deg"] <= 180
    assert HB.distance_angle_map([_dimer(False)], 4.0)["n_points"] == 0
    with pytest.raises(HB.HydrogenBondError):
        HB.distance_angle_map(frames, 0.0)


def test_run_analysis_writes_tables_figures_and_csv(tmp_path):
    src = REPO / "examples" / "openmm_spce_nvt_generated"
    d = tmp_path / "md"
    shutil.copytree(src, d, ignore=shutil.ignore_patterns("analysis"))
    res = run_analysis(d, AnalysisOptions(hbond="3.5,150", hbond_lifetime=True, hbond_cdf=4.0, stride=1))
    assert "hbond" in res.tables and "hbond_lifetime" in res.tables and "hbond_map" in res.tables, res.notes
    for k in ("hbond", "hbond_lifetime", "hbond_map"):
        assert Path(res.figures[k]).is_file()
    lt = res.tables["hbond_lifetime"]
    assert lt["intermittent"][0] == 1.0 and len(lt["lag"]) == lt["tau_max"] + 1
    with open(d / "analysis" / "hbond_lifetime.csv", newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == lt["tau_max"] + 1 and float(rows[0]["c_continuous"]) == 1.0
    mp = res.tables["hbond_map"]
    assert mp["rmax_A"] == 4.0 and len(mp["counts"]) == 60 and (d / "analysis" / "hbond_map.csv").is_file()
    js = json.loads((d / "analysis" / "summary.json").read_text(encoding="utf-8"))
    assert js["figures"]["hbond_lifetime"] == "analysis/hbond_lifetime.png" and js["figures"]["hbond_map"] == "analysis/hbond_map.png"
    assert any("寿命" in n or "lifetime" in n for n in res.notes)
    res2 = run_analysis(d, AnalysisOptions(hbond_lifetime=True))
    assert "hbond_lifetime" not in res2.tables and any("--hbond" in n for n in res2.notes)


def test_cli_options_and_gui_fields(tmp_path):
    from adit.analysis.cli import main as analyze
    from adit.gui import analysis_fields as AF
    src = REPO / "examples" / "openmm_spce_nvt_generated"
    d = tmp_path / "md"
    shutil.copytree(src, d, ignore=shutil.ignore_patterns("analysis"))
    assert analyze([str(d), "--hbond", "3.5,150", "--hbond-lifetime", "--hbond-cdf", "4"]) == 0
    assert (d / "analysis" / "hbond_lifetime.png").is_file() and (d / "analysis" / "hbond_map.png").is_file()
    opts = AnalysisOptions(hbond="3.5,150", hbond_lifetime=True, hbond_cdf=4.0)
    f = AF.fields_from_options(opts)
    assert f["hbond_lifetime"] == "on" and f["hbond_cdf"] == "4"
    back = AF.options_from_fields(f)
    assert back.hbond_lifetime and back.hbond_cdf == 4.0
    assert not AF.options_from_fields({"hbond_cdf": ""}).hbond_lifetime
    res = run_analysis(d, opts)
    secs = {s.key: s for s in AF.result_sections(res)}
    assert "hbond_lifetime" in secs and len(secs["hbond_lifetime"].rows) == 2 and math.isfinite(float(secs["hbond_lifetime"].rows[0][1]))
