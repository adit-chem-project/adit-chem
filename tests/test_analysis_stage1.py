
import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.io import read

from adit.analysis import AnalysisOptions, compute, run_analysis
from adit.analysis import trajectory as trj
from adit.analysis.readers import HARTREE_EV, load_run

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"


def _copy(name: str, tmp_path: Path) -> Path:
    d = tmp_path / name
    shutil.copytree(EX / name, d, ignore=shutil.ignore_patterns("analysis"))
    return d


def test_xtb_thermochemistry_charges_gap_dipole():
    r = load_run(EX / "xtb_vib_water_generated")
    th = {it["key"]: it for it in r.thermo["items"]}
    assert r.thermo["temperature_k"] == pytest.approx(298.15)
    assert th["zpe"]["value_eh"] == pytest.approx(0.008147568762)
    assert th["g_rrho_contrib"]["value_eh"] == pytest.approx(-0.010288936689)
    assert th["total_enthalpy"]["value_eh"] == pytest.approx(-4.965840055468)
    assert th["total_free_energy"]["value_eh"] == pytest.approx(-4.988153848206)
    assert th["h_t"]["value_eh"] == pytest.approx(0.120249e-01)
    assert th["total_enthalpy"]["value_ev"] == pytest.approx(-4.965840055468 * HARTREE_EV)
    assert th["zpe"]["line"].startswith("output.log:")
    lines = (EX / "xtb_vib_water_generated" / "output.log").read_text(encoding="utf-8").splitlines()
    assert "zero point energy" in lines[int(th["zpe"]["line"].split(":")[1]) - 1]
    q = r.charges[0]
    assert q["definition"] == "Mulliken (GFN2-xTB)" and q["values"] == pytest.approx([-0.50329938, 0.25164969, 0.25164969])
    assert q["source"].startswith("xtbout.json:")
    assert r.electronic["homo_lumo_gap_ev"] == pytest.approx(6.33829684)
    assert r.electronic["dipole_norm_debye"] == pytest.approx(1.13716233 * 2.541746473, rel=1e-6)


def test_xtb_optimization_has_no_thermo_but_has_gap():
    r = load_run(EX / "xtb_water_generated")
    assert r.thermo is None and r.electronic.get("homo_lumo_gap_ev") is not None


def test_dftb_gross_charges_with_line():
    r = load_run(EX / "water_generated")
    q = r.charges[0]
    assert q["definition"] == "Mulliken" and q["source"] == "detailed.out:13"
    assert q["values"] == pytest.approx([-0.59260702, 0.29630351, 0.29630351])
    md = load_run(EX / "dftb_md_water_generated").charges[0]
    assert "最後" in md["note"]


def test_qe_dynmat_frequencies():
    r = load_run(EX / "qe_si_phonon_generated")
    assert r.frequencies_cm1 == pytest.approx([0, 0, 0, 587.449126, 587.449126, 587.449126])


ORCA_THERMO = """\
FINAL SINGLE POINT ENERGY      -157.180901163000
--------------------------
THERMOCHEMISTRY AT 298.15K
--------------------------

Temperature         ...   298.15 K
Pressure            ...     1.00 atm
Zero point energy                ...      0.10746477 Eh      67.44 kcal/mol
Total thermal energy                   -157.07337829 Eh
Total Enthalpy                    ...   -157.07243408 Eh
Final entropy term                ...      0.03353675 Eh     21.04 kcal/mol
Final Gibbs free energy         ...   -157.10597083 Eh
G-E(el)                           ...      0.08017472 Eh     50.31 kcal/mol
MULLIKEN ATOMIC CHARGES
ORCA TERMINATED NORMALLY
"""


def test_orca_thermochemistry_lines_from_the_tutorial(tmp_path):
    (tmp_path / "orca.inp").write_text("! B3LYP def2-SVP Freq\n* xyz 0 1\nH 0 0 0\nH 0 0 0.74\n*\n", encoding="utf-8")
    (tmp_path / "output.log").write_text(ORCA_THERMO, encoding="utf-8")
    r = load_run(tmp_path)
    th = {it["key"]: it["value_eh"] for it in r.thermo["items"]}
    assert r.thermo["temperature_k"] == pytest.approx(298.15)
    assert th == pytest.approx({"zpe": 0.10746477, "total_thermal_energy": -157.07337829, "total_enthalpy": -157.07243408,
                                "final_entropy_term": 0.03353675, "final_gibbs_free_energy": -157.10597083, "g_minus_eel": 0.08017472})
    assert r.charges == []


def test_summary_json_has_new_tables(tmp_path):
    d = _copy("xtb_vib_water_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions())
    s = json.loads((d / "analysis" / "summary.json").read_text(encoding="utf-8"))
    assert {"charges", "thermochemistry", "electronic"} <= set(s["tables"])
    assert s["tables"]["charges"][0]["atoms"] == ["O1", "H2", "H3"] and abs(s["tables"]["charges"][0]["sum"]) < 1e-6
    txt = res.summary_text()
    assert "原子の電荷 (Mulliken (GFN2-xTB)" in txt and "TOTAL FREE ENERGY -4.988154 Eh" in txt and "HOMO-LUMO ギャップ 6.3383 eV" in txt


def test_qe_stream_matches_ase():
    for name in ("qe_md_si_generated", "qe_si_generated"):
        r = load_run(EX / name)
        ref = read(EX / name / "output.log", format="espresso-out", index=":")
        assert len(r.frames) == len(ref) and len(r.energies_ev) == len(ref)
        for a, b in zip(r.frames, ref):
            assert np.allclose(a.positions, b.positions, atol=1e-6) and np.allclose(a.cell, b.cell, atol=1e-6)
        assert np.allclose(r.energies_ev, [a.get_potential_energy() for a in ref])


def test_vasp_reads_xdatcar_trajectory():
    r = load_run(EX / "vasp_h2o_generated")
    ref = read(EX / "vasp_h2o_generated" / "XDATCAR", index=":")
    assert r.frame_source == "XDATCAR" and len(r.frames) == len(ref) == 5
    assert all(np.allclose(a.positions, b.positions) and np.allclose(a.cell, b.cell) for a, b in zip(r.frames, ref))


def test_trajectory_is_lazy_sequence():
    r = load_run(EX / "dftb_md_water_generated")
    t = r.frames
    assert isinstance(t, trj.Trajectory) and len(t) == 21
    sub = t[5::4]
    assert isinstance(sub, trj.Trajectory) and len(sub) == len(range(21)[5::4]) == 4
    all_frames = list(t)
    assert np.allclose(sub[1].positions, all_frames[9].positions) and np.allclose(t[-1].positions, all_frames[-1].positions)
    assert np.allclose(all_frames[0].positions, [[0, -1, 0], [0, 0, 0.783064], [0, 0, -0.783064]])
    assert [a.get_chemical_symbols() for a in sub] == [["O", "H", "H"]] * 4
    with pytest.raises(IndexError):
        t[21]
    assert r.frame_dt_fs == pytest.approx(5.0)


def test_estimate_before_reading(tmp_path):
    t = trj.Trajectory(EX / "dftb_md_water_generated" / "geo_end.xyz", "xyz")
    est = t.estimate_total_frames()
    assert t._shared["total"] is None, "見積もりのために全部を数えてはいけない"
    assert abs(est - 21) <= 2


def test_budget_stops_before_reading_with_fake_size(tmp_path, monkeypatch):
    d = _copy("dftb_md_water_generated", tmp_path)
    monkeypatch.setattr(trj.Trajectory, "file_size", lambda self: 50 * 2**30)
    read_calls = []
    real = trj.Trajectory._frames_xyz
    monkeypatch.setattr(trj.Trajectory, "_frames_xyz", lambda self, want, last: (read_calls.append(last), (yield from real(self, want, last)))[1])
    with pytest.raises(trj.TrajectoryTooLarge) as ei:
        run_analysis(d, AnalysisOptions(energy=False, temperature=False, bonds=False, msd=True, memory_budget_mb=1024))
    assert "--stride" in str(ei.value) and ei.value.suggested_stride > 1
    assert all(last == 0 for last in read_calls), "上限を超えたら、最初のフレーム以外を読んではいけない"
    est = trj.Trajectory(d / "geo_end.xyz", "xyz").estimate_total_frames()
    info = trj.check_budget(trj.Trajectory(d / "geo_end.xyz", "xyz")[:: ei.value.suggested_stride], 3, 1024)
    assert info["estimated_mb"] <= 1024 and est > 10**8


def test_cli_prints_the_stride_advice(tmp_path, monkeypatch, capsys):
    from adit.analysis.cli import main
    d = _copy("dftb_md_water_generated", tmp_path)
    monkeypatch.setattr(trj.Trajectory, "file_size", lambda self: 50 * 2**30)
    assert main([str(d), "--msd"]) == 1
    assert "--stride" in capsys.readouterr().err


def test_stride_on_small_file(tmp_path):
    d = _copy("dftb_md_water_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions(rdf=True, msd=True, stride=2))
    t = res.tables["trajectory"]
    assert t["n_frames_used"] == len(range(21)[::2]) == 11 and t["dt_used_fs"] == pytest.approx(10.0)
    assert res.tables["msd"]["dt_fs"] == pytest.approx(10.0) and res.tables["rdf"]["n_frames"] == 11
    with pytest.raises(ValueError):
        run_analysis(d, AnalysisOptions(stride=0))


def test_msd_fft_matches_direct_average_over_origins():
    rng = np.random.default_rng(1)
    pos = np.cumsum(rng.normal(size=(60, 5, 3)), axis=0)
    direct = np.array([np.mean([np.mean(np.sum((pos[k + m] - pos[k]) ** 2, axis=1)) for k in range(60 - m)]) for m in range(60)])
    assert np.allclose(compute.msd_fft(pos), direct)
    assert np.allclose(compute.msd_fft(pos, chunk_bytes=1), direct)


def test_msd_fit_range_default_and_user():
    t = np.arange(11) * 10.0
    m = 6 * 1e-4 * t  # D = 1e-4 Å^2/fs → 1e-5 cm^2/s
    assert compute.fit_range_fs(t, None) == (10.0, 50.0)
    assert compute.diffusion_fit(t, m, (20.0, 60.0)) == pytest.approx(1e-5)


def test_msd_by_element_and_fit_range_in_table(tmp_path):
    d = _copy("dftb_md_water_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions(msd=True, msd_fit_fs=(20.0, 60.0)))
    m = res.tables["msd"]
    assert m["fit_range_fs"] == [20.0, 60.0] and m["fit_range_user"] and set(m["by_element"]) == {"H", "O"}
    assert m["method"] == "fft_multiple_time_origins" and len(m["lag_fs"]) == len(m["msd_A2"]) == 21


def test_block_average_and_autocorrelation_time_ar1():
    rng = np.random.default_rng(2)
    phi, n = 0.8, 40000
    x = np.zeros(n)
    for k in range(1, n):
        x[k] = phi * x[k - 1] + rng.normal()
    ac = compute.autocorrelation_time(x)
    assert ac["tau_int"] == pytest.approx(4.5, rel=0.15) and ac["window_reached"]
    blocks = compute.block_average(x)
    assert blocks[0]["block_size"] == 1 and blocks[1]["block_size"] == 2 and blocks[0]["n_blocks"] == n
    plateau = [b["sem"] for b in blocks if 64 <= b["block_size"] <= 512]
    assert np.mean(plateau) == pytest.approx(ac["sem"], rel=0.25)
    assert blocks[0]["sem"] < 0.5 * np.mean(plateau)


def test_negative_autocorrelation_time_gives_no_error():
    x = np.array([1.0, -1.0] * 10)
    ac = compute.autocorrelation_time(x)
    assert ac["tau_int"] <= 0 and ac["sem"] is None


def test_md_stats_in_summary(tmp_path):
    d = _copy("qe_md_si_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions())
    st = res.tables["timeseries_stats"]
    assert set(st) == {"temperature", "energy"} and st["temperature"]["n"] == 20 and st["temperature"]["unit"] == "K"
    assert st["temperature"]["blocks"][0]["n_blocks"] == 20 and "blocking" in res.figures
    assert "温度の統計" in res.summary_text()


def _brute_rdf(frames, a, b, rmax, nbins):
    edges = np.linspace(0, rmax, nbins + 1)
    shell = 4 / 3 * np.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    g = np.zeros(nbins); n = np.zeros(nbins)
    for fr in frames:
        s = np.array(fr.get_chemical_symbols()); ia, ib = np.where(s == a)[0], np.where(s == b)[0]
        dd = fr.get_all_distances(mic=True)[np.ix_(ia, ib)]
        if a == b:
            dd = dd[~np.eye(len(ia), dtype=bool)]
        h = np.histogram(dd.ravel(), bins=edges)[0]
        npairs = len(ia) * len(ib) if a != b else len(ia) * (len(ia) - 1)
        g += h * fr.get_volume() / (npairs * shell); n += np.cumsum(h) / len(ia)
    return g / len(frames), n / len(frames)


def test_rdf_neighbor_list_matches_all_distances():
    rng = np.random.default_rng(3)
    L = 9.0
    frames = [Atoms(["O"] * 12 + ["H"] * 20, positions=rng.uniform(0, L, (32, 3)), cell=[L, L, L], pbc=True) for _ in range(4)]
    acc = compute.RDFAccumulator([("O", "O"), ("H", "O"), ("O", "H")], rmax=4.4, nbins=44)
    for fr in frames:
        acc.add(fr)
    for a, b in acc.pairs:
        g, n = _brute_rdf(frames, a, b, 4.4, 44)
        res = acc.result((a, b))
        assert np.allclose(res["g"], g) and np.allclose(res["n"], n)
    assert np.allclose(acc.result(("H", "O"))["n_reverse"], acc.result(("O", "H"))["n"])


def test_coordination_of_diamond_si():
    from ase.build import bulk
    a = bulk("Si", "diamond", a=5.43, cubic=True)
    acc = compute.RDFAccumulator([("Si", "Si")], rmax=2.6, nbins=26)
    acc.add(a)
    assert acc.result(("Si", "Si"))["n"][-1] == pytest.approx(4.0)


def test_unwrap_records_steps_near_half_cell():
    acc = compute.UnwrapAccumulator(3)
    for x in (1.0, 2.0, 6.2):
        acc.add(Atoms("H", positions=[(x, 0, 0)], cell=[10, 10, 10], pbc=True))
    assert acc.large_step_count == 1
    assert acc.max_step_fraction == pytest.approx(0.42)


def test_zdensity_slab():
    pos = [[1, 1, 2.05], [3, 3, 2.05], [1, 3, 8.05], [2, 2, 8.05]]
    a = Atoms(["Pt", "Pt", "O", "O"], positions=pos, cell=[4, 4, 10], pbc=True)
    z = compute.ZDensityAccumulator(0.5)
    z.add(a); z.add(a)
    r = z.result()
    assert r["c_perpendicular"] and r["axis_length_A"] == pytest.approx(10.0) and r["bin_A"] == pytest.approx(0.5)
    area = 16.0
    for el, zc in (("Pt", 2.25), ("O", 8.25)):
        v = r["number_density_A3"][el]
        assert np.sum(v) * area * r["bin_A"] == pytest.approx(2.0)
        assert r["z_A"][np.argmax(v)] == pytest.approx(zc)


def test_zdensity_and_rdf_json_from_analysis(tmp_path):
    d = _copy("qe_md_si_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions(rdf=True, zdens=True))
    assert "zdensity" in res.figures and not res.tables["zdensity"]["c_perpendicular"]
    rdf = json.loads((d / "analysis" / "rdf.json").read_text(encoding="utf-8"))
    assert set(rdf["Si-Si"]) >= {"r", "g", "n", "n_reverse", "r_upper"}
    zd = json.loads((d / "analysis" / "zdensity.json").read_text(encoding="utf-8"))
    assert len(zd["z_A"]) == len(zd["number_density_A3"]["Si"])


def _broken_water_frames(n=3):
    frames = []
    for k in range(n):
        o = np.array([0.1 + 0.01 * k, 5.0, 5.0])
        h1 = o + [-0.757, 0.586, 0.0]
        h2 = o + [0.757, 0.586, 0.0]
        a = Atoms("OHH", positions=[o, h1, h2], cell=[10, 10, 10], pbc=True)
        a.wrap()
        frames.append(a)
    return frames


def test_export_with_molecule_unwrapping(tmp_path):
    from adit.analysis.export import TrajectoryExporter
    ex = TrajectoryExporter(tmp_path / "export", unwrap_molecules=True)
    frames = _broken_water_frames()
    assert np.linalg.norm(frames[0].positions[1] - frames[0].positions[0]) > 5
    for fr in frames:
        ex.add(fr)
    info = ex.close(run_dir=tmp_path, code="dftbplus", source="geo_end.xyz", dt_frame_fs=2.5, stride=2, skip=0, n_total=6, rdf_cutoff=5.0)
    assert info["n_frames"] == 3 and info["n_molecules"] == 1 and info["periodic_networks"] == 0
    back = read(tmp_path / "export" / "trajectory.extxyz", index=":")
    assert len(back) == 3 and np.allclose(back[0].cell.lengths(), [10, 10, 10]) and all(back[0].pbc)
    for a in back:
        d = np.linalg.norm(a.positions[1:] - a.positions[0], axis=1)
        assert np.allclose(d, np.sqrt(0.757**2 + 0.586**2), atol=1e-6)
    assert len(read(tmp_path / "export" / "trajectory.pdb", index=":")) == 3
    assert len(read(tmp_path / "export" / "trajectory.xyz", index=":")) == 3
    vmd = (tmp_path / "export" / "view.vmd").read_text(encoding="utf-8")
    assert "mol new trajectory.xyz type xyz waitfor all" in vmd and "pbc set {10.000000 10.000000 10.000000 90.0000 90.0000 90.0000} -all" in vmd
    readme = (tmp_path / "export" / "export_README.txt").read_text(encoding="utf-8")
    assert "1000.0000 pm" in readme and "5 fs" in readme and "travis -p trajectory.xyz" in readme and "O 1, H 2" in readme
    compile((tmp_path / "export" / "ovito_pipeline.py").read_text(encoding="utf-8"), "ovito_pipeline.py", "exec")


@pytest.mark.parametrize("language", ["ja", "en"])
def test_nonperiodic_export_does_not_claim_unwrapping_or_cell(tmp_path, monkeypatch, language):
    from adit import lang
    from adit.analysis.export import TrajectoryExporter

    monkeypatch.setattr(lang, "LANGUAGE", language)
    ex = TrajectoryExporter(tmp_path / "export", unwrap_molecules=True)
    ex.add(Atoms("OHH", positions=[(0, 0, 0), (1.27, 0, 0), (-1.27, 0, 0)]))
    info = ex.close(run_dir=tmp_path, code="dftbplus", source="geo_end.xyz", dt_frame_fs=1.0,
                    stride=1, skip=0, n_total=1, rdf_cutoff=3.0)
    readme = (tmp_path / "export" / "export_README.txt").read_text(encoding="utf-8")
    pdb = (tmp_path / "export" / "trajectory.pdb").read_text(encoding="utf-8")
    assert info["unwrap_requested"] and not info["unwrap_molecules"]
    assert "CRYST1" not in pdb
    if language == "ja":
        assert "つなぎ直しは行っていません" in readme and "セル情報なし" in readme
        assert "伸びた結合は別グループ" in readme and "セルと周期境界はありません" in readme
        assert "拡張 xyz (セル・周期境界なし)" in readme
        assert "つなぎ直しました" not in readme
    else:
        assert "no periodic unwrapping was performed" in readme and "no cell information" in readme
        assert "Stretched bonds may be split" in readme and "no cell or periodic boundaries" in readme
        assert "largest group: 1 atom" in readme
        assert "were made whole" not in readme


def test_export_keeps_periodic_networks(tmp_path):
    from ase.build import bulk
    from adit.analysis.export import TrajectoryExporter
    a = bulk("Si", "diamond", a=5.43, cubic=True)
    ex = TrajectoryExporter(tmp_path / "e", unwrap_molecules=True)
    ex.add(a)
    info = ex.close(run_dir=tmp_path, code="vasp", source="XDATCAR", dt_frame_fs=None, stride=1, skip=0, n_total=1, rdf_cutoff=2.7)
    assert info["periodic_networks"] == 1
    assert np.allclose(read(tmp_path / "e" / "trajectory.extxyz").positions, a.positions)
    assert "不明" in (tmp_path / "e" / "export_README.txt").read_text(encoding="utf-8")


def test_cli_export_with_stride(tmp_path, capsys):
    from adit.analysis.cli import main
    d = _copy("qe_md_si_generated", tmp_path)
    assert main([str(d), "--export", "--stride", "2"]) == 0
    e = d / "analysis" / "export"
    assert len(read(e / "trajectory.extxyz", index=":")) == 10 and len(read(e / "trajectory.pdb", index=":")) == 10
    assert "4 fs" in (e / "export_README.txt").read_text(encoding="utf-8") and "書き出し:" in capsys.readouterr().out


def test_export_writes_vmd_tcl_ovito_templates_and_travis_answers(tmp_path):
    from adit.analysis.export import TrajectoryExporter
    ex = TrajectoryExporter(tmp_path / "e", unwrap_molecules=False)
    for fr in _broken_water_frames():
        fr.set_cell([10, 12, 10])
        ex.add(fr)
    info = ex.close(run_dir=tmp_path, code="dftbplus", source="geo_end.xyz", dt_frame_fs=2.5, stride=2, skip=0, n_total=3, rdf_cutoff=5.0,
                    select="element O and z < 10")
    names = {Path(f).name for f in info["files"]}
    assert {"vmd_load.tcl", "travis_rdf.in", "travis_cdf.in", "travis_msd.in", "travis_hbond.in", "travis_acf.in"} <= names
    assert all((tmp_path / "e" / n).is_file() for n in names)
    tcl = (tmp_path / "e" / "vmd_load.tcl").read_text(encoding="utf-8")
    assert "mol new trajectory.xyz type xyz waitfor all" in tcl and "pbc set {10.000000 12.000000 10.000000 90.0000 90.0000 90.0000} -all" in tcl
    assert "mol modselect 1 $m {(name O and z < 10)}" in tcl and "mol modstyle 1 $m VDW" in tcl and "color Display Background white" in tcl
    ovito = (tmp_path / "e" / "ovito_pipeline.py").read_text(encoding="utf-8")
    compile(ovito, "ovito_pipeline.py", "exec")
    for cls in ("CommonNeighborAnalysisModifier", "PolyhedralTemplateMatchingModifier", "AcklandJonesModifier", "WignerSeitzAnalysisModifier",
                "CalculateDisplacementsModifier", "AtomicStrainModifier", "ClusterAnalysisModifier", "SpatialBinningModifier",
                "TimeAveragingModifier", "ExpressionSelectionModifier"):
        assert f"om.{cls}(" in ovito
    assert "expression='(ParticleType == \"O\" && Position.Z < 10)'" in ovito
    assert "compute_com=True" in ovito and 'data.tables["binning[average]"]' in ovito
    rdf = (tmp_path / "e" / "travis_rdf.in").read_text(encoding="utf-8").splitlines()
    answers = [l for l in rdf if not l.startswith("!")]
    assert answers == ["", "no", "1000.0000", "1200.0000", "1000.0000", "", "", "rdf", ""]
    msd = (tmp_path / "e" / "travis_msd.in").read_text(encoding="utf-8")
    assert "\nmsd\n" in msd and "未確認" in msd
    readme = (tmp_path / "e" / "export_README.txt").read_text(encoding="utf-8")
    assert "travis -p trajectory.xyz -i travis_rdf.in" in readme and "vmd -e vmd_load.tcl" in readme
    assert "1 フレームあたり 5 fs" in readme and "acf (vacf ではありません)" in readme


def test_export_travis_answers_only_for_orthorhombic_cells(tmp_path):
    from adit.analysis.export import TrajectoryExporter
    ex = TrajectoryExporter(tmp_path / "cubic")
    ex.add(Atoms("OHH", positions=[(0, 0, 0), (0.96, 0, 0), (-0.24, 0.93, 0)], cell=[10, 10, 10], pbc=True))
    info = ex.close(run_dir=tmp_path, code="dftbplus", source="x", dt_frame_fs=None, stride=1, skip=0, n_total=1, rdf_cutoff=5.0)
    answers = [l for l in (tmp_path / "cubic" / "travis_rdf.in").read_text(encoding="utf-8").splitlines() if not l.startswith("!")]
    assert answers == ["", "", "1000.0000", "", "", "rdf", ""] and info["travis_answer_files"] == sorted(info["travis_answer_files"])
    ovito = (tmp_path / "cubic" / "ovito_pipeline.py").read_text(encoding="utf-8")
    assert "expression='ParticleType == \"<元素>\" && Position.Z < <z の上限 Å>'" in ovito
    ex = TrajectoryExporter(tmp_path / "hex")
    ex.add(Atoms("C2", positions=[(0, 0, 0), (1.42, 0, 0)], cell=[[2.46, 0, 0], [-1.23, 2.13, 0], [0, 0, 10]], pbc=True))
    info = ex.close(run_dir=tmp_path, code="vasp", source="x", dt_frame_fs=None, stride=1, skip=0, n_total=1, rdf_cutoff=5.0)
    assert info["travis_answer_files"] == [] and not (tmp_path / "hex" / "travis_rdf.in").exists()
    assert "直方体でないセル" in (tmp_path / "hex" / "export_README.txt").read_text(encoding="utf-8")
    ex = TrajectoryExporter(tmp_path / "mol")
    ex.add(Atoms("OHH", positions=[(0, 0, 0), (0.96, 0, 0), (-0.24, 0.93, 0)]))
    info = ex.close(run_dir=tmp_path, code="xtb", source="x", dt_frame_fs=None, stride=1, skip=0, n_total=1, rdf_cutoff=5.0)
    assert info["travis_answer_files"] == [] and "pbc set" not in (tmp_path / "mol" / "vmd_load.tcl").read_text(encoding="utf-8")


def test_cli_export_passes_the_selection(tmp_path):
    from adit.analysis.cli import main
    d = _copy("qe_md_si_generated", tmp_path)
    assert main([str(d), "--export", "--select", "index 1-4"]) == 0
    tcl = (d / "analysis" / "export" / "vmd_load.tcl").read_text(encoding="utf-8")
    assert "mol modselect 1 $m {index 0 to 3}" in tcl
