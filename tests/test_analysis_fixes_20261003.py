"""Fixes from the 2026-10-03 audit of the analysis module (H1, H2, M1-M14, L1, L2, L4, L6, S1, S3)."""

import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.build import bulk, molecule

from adit.analysis import AnalysisOptions, run_analysis
from adit.analysis import bands as B
from adit.analysis import compute
from adit.analysis import crest
from adit.analysis import heavy_setup as H
from adit.analysis import local_order as LO
from adit.analysis import msd_worker as W
from adit.analysis import readers, readers_extra, report
from adit.analysis.geometry_series import dihedral_series
from adit.analysis.hbond import criteria_note
from adit.analysis.readers_qc import read_gamess
from adit.analysis.thermo import ThermoOptions, compute_thermo, summary_lines

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"
A_CU = 3.615


def _copy(name: str, tmp_path: Path) -> Path:
    d = tmp_path / name
    shutil.copytree(EX / name, d, ignore=shutil.ignore_patterns("analysis"))
    return d


def _quiet():
    return dict(energy=False, temperature=False, bonds=False, spacegroup=False)


# ---------------- H1: every periodic image counts ----------------
def test_fcc_conventional_cell_counts_all_images():
    cu = bulk("Cu", "fcc", a=A_CU, cubic=True)  # 4 atoms, cell 3.615 Å < 2 x cutoff
    got = LO.coordination(cu, 3.0)
    assert got["coordination"] == [12, 12, 12, 12]
    st = LO.steinhardt(cu, 3.0)
    assert np.allclose(st["q4"], 0.1909, atol=1e-3) and np.allclose(st["q6"], 0.5745, atol=1e-3)
    assert LO.clusters(cu, 3.0)["n_clusters"] == 1
    adf = LO.angle_distribution([cu], cutoff=3.0)
    assert adf["n_angles"] == 4 * 12 * 11 // 2
    centers, counts = np.array(adf["angle_deg"]), np.array(adf["counts"])
    for angle in (60.0, 90.0, 120.0, 180.0):  # the angles between nearest neighbours of an FCC site
        assert counts[np.abs(centers - angle) <= 2.0].sum() > 0, angle


def test_primitive_one_atom_cell_sees_its_own_images():
    prim = bulk("Cu", "fcc", a=A_CU)
    assert LO.coordination(prim, 3.0)["coordination"] == [12]
    st = LO.steinhardt(prim, 3.0)
    assert st["q4"][0] == pytest.approx(0.1909, abs=1e-3) and st["q6"][0] == pytest.approx(0.5745, abs=1e-3)


def test_molecule_without_a_cell_still_works():
    w = molecule("H2O")
    assert LO.coordination(w, 1.2)["coordination"] == [2, 1, 1]


# ---------------- H2 / M1: DOS grid flag and ISPIN = 2 ----------------
def test_many_eigenvalues_are_broadened_and_a_grid_passes_through():
    rng = np.random.default_rng(0)
    eigs = np.sort(rng.normal(-5, 3, 250))
    x, y = compute.dos(eigs, np.ones(250), sigma=0.1)
    assert len(x) == 800 and np.trapezoid(y, x) == pytest.approx(250.0, rel=1e-3)
    assert y.max() > 1.0
    grid = np.linspace(-10, 10, 101)
    vals = np.exp(-grid ** 2)
    gx, gy = compute.dos(grid, vals, sigma=0.1, is_grid=True)
    assert np.array_equal(gx, grid) and np.array_equal(gy, vals)


def _write_doscar(d: Path, two_spins: bool, nedos: int = 101) -> None:
    d.mkdir(parents=True, exist_ok=True)
    (d / "INCAR").write_text("ISPIN = 2\n" if two_spins else "ISPIN = 1\n", encoding="utf-8")
    lines = ["2 2 0 0", "x", "x", "x", "x", f" 10.0 -10.0 {nedos} 0.5 1.0"]
    for e in np.linspace(-10, 10, nedos):
        up, down = 1.0 + 0.1 * e, 2.0
        lines.append(f" {e:.4f} {up:.4f} {down:.4f} 0.0 0.0" if two_spins else f" {e:.4f} {up:.4f} 0.0")
    (d / "DOSCAR").write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_doscar_sets_the_grid_flag_and_sums_both_spins(tmp_path):
    _write_doscar(tmp_path / "v2", two_spins=True)
    r = readers.load_run(tmp_path / "v2")
    assert r.dos_is_grid and r.fermi_ev == 0.5
    expected = 1.0 + 0.1 * np.linspace(-10, 10, 101) + 2.0
    assert np.allclose(r.eigen_weights, expected)
    assert any("ISPIN = 2" in n and ("和" in n or "summed" in n) for n in r.notes)
    _write_doscar(tmp_path / "v1", two_spins=False)
    r1 = readers.load_run(tmp_path / "v1")
    assert np.allclose(r1.eigen_weights, 1.0 + 0.1 * np.linspace(-10, 10, 101))


def test_doscar_grid_is_not_broadened_again_and_band_out_is(tmp_path):
    _write_doscar(tmp_path / "vasp", two_spins=False)
    res = run_analysis(tmp_path / "vasp", AnalysisOptions(dos=True, **_quiet()))
    assert res.tables["dos"] == {"n_eigen": 101, "fermi_ev": 0.5, "broadened": False, "sigma_ev": None}
    assert "そのまま" in res.summary_text() or "as is" in res.summary_text()
    d = _copy("dftb_md_water_generated", tmp_path)
    res2 = run_analysis(d, AnalysisOptions(dos=True, dos_sigma=0.2, **_quiet()))
    assert res2.tables["dos"]["broadened"] is True and res2.tables["dos"]["sigma_ev"] == 0.2
    assert "σ = 0.2 eV" in res2.summary_text()


# ---------------- M2: xtb MD table ----------------
def test_xtb_md_table_in_the_fortran_format(tmp_path):
    d = tmp_path / "xtb"; d.mkdir()
    (d / "xtb.inp").write_text("$md\n   temp=300\n   dump=50.0\n   step=1.0\n$end\n", encoding="utf-8")
    (d / "xtb.trj").write_text("".join(f"3\n energy: {-5.07 - 0.001 * k:.8f} gnorm: 0.1 xtb: 6.6.1\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n"
                                       for k in range(3)), encoding="utf-8")
    rows = ["%7i%8.2f%13.5f%9.4f%6.0f%6.0f%12.5f" % (200, 0.20, -5.07043, 0.01094, 301., 305., -5.05949),
            "%7i%8.2f%13.5f%9.4f%6.0f%6.0f%12.5f%14.6E" % (400, 0.40, -5.07050, 0.01100, 302., 310., -5.05950, 1.2e-3)]
    (d / "output.log").write_text("      time (ps)    <Epot>      Ekin   <T>   T     Etot\n" + "\n".join(rows) + "\n", encoding="utf-8")
    r = readers.load_run(d)
    assert r.temperatures_k == [305.0, 310.0]
    assert r.temperature_times_fs == pytest.approx([200.0, 400.0])
    assert r.frame_dt_fs == 50.0 and r.times_fs == [0.0, 50.0, 100.0]
    assert len(r.energies_ev) == 3


# ---------------- M3: dihedral sign ----------------
def test_dihedral_sign_matches_ase():
    rng = np.random.default_rng(1)
    signs = set()
    for _ in range(12):
        at = Atoms("CCCC", positions=rng.random((4, 3)) * 3)
        ours = dihedral_series([at], 0, 1, 2, 3)[0]
        ref = at.get_dihedral(0, 1, 2, 3)
        assert (ours - ref) % 360.0 == pytest.approx(0.0, abs=1e-6) or (ours - ref) % 360.0 == pytest.approx(360.0, abs=1e-6)
        signs.add(np.sign(ours))
    assert signs == {-1.0, 1.0}


# ---------------- M4 / M9 / M13: selection, absent element, dt unknown ----------------
def test_select_restricts_the_msd_and_the_summary_counts_atoms(tmp_path):
    d = _copy("dftb_md_water_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions(msd=True, select="element H", **_quiet()))
    m = res.tables["msd"]
    assert m["n_atoms_used"] == 2 and m["selection"] == "element H" and set(m["by_element"]) == set()
    assert m["D_cm2_s"] is not None
    assert "使った原子 2 個" in res.summary_text()
    res2 = run_analysis(d, AnalysisOptions(msd=True, msd_species="O", select="element H", **_quiet()))
    assert res2.tables["msd"]["D_cm2_s"] is None and res2.tables["msd"]["n_atoms_used"] == 0
    assert any("--select" in n and "O" in n for n in res2.notes)


def test_absent_element_gives_no_d_and_names_the_present_ones(tmp_path):
    d = _copy("dftb_md_water_generated", tmp_path)
    res = run_analysis(d, AnalysisOptions(msd=True, msd_species="Na", **_quiet()))
    m = res.tables["msd"]
    assert m["D_cm2_s"] is None and "Na" in m["reason"] and "H, O" in m["reason"]
    txt = res.summary_text()
    assert "MSD (平均二乗変位): 出せません" in txt and "注: MSD を出せません" in txt
    a = compute.msd_analysis(np.zeros((5, 2, 3)), "Xx", None, symbols=["Ar", "Ar"], dt_fs=1.0)
    assert a["D_cm2_s"] is None and a["n_atoms"] == 0 and "Ar" in a["reason"]


def _fake_run(tmp_path, monkeypatch, *, frame_dt_fs, n_frames=10, series_points=100):
    rng = np.random.default_rng(0)
    frames = []
    for k in range(n_frames):
        pos = rng.random((4, 3)) * 8.0
        frames.append(Atoms("Ar4", positions=pos, cell=[8.0, 8.0, 8.0], pbc=True))
    fake = readers.RunData("lammps", tmp_path, frames=frames, temperatures_k=[300.0 + rng.normal() for _ in range(series_points)],
                           times_fs=[float(k) for k in range(series_points)], frame_dt_fs=frame_dt_fs)
    monkeypatch.setattr(report, "load_run", lambda _d, _code=None: fake)
    return fake


def test_per_atom_d_is_skipped_when_dt_is_unknown(tmp_path, monkeypatch):
    _fake_run(tmp_path, monkeypatch, frame_dt_fs=None)
    res = run_analysis(tmp_path, AnalysisOptions(msd=True, msd_per_atom=True, energy=False, bonds=False, spacegroup=False))
    assert "msd_per_atom" not in res.tables
    assert any("原子ごとの D は出しません" in n or "no per-atom D" in n for n in res.notes)


# ---------------- L1: --skip in time ----------------
def test_skip_is_applied_in_time_to_a_series_with_its_own_stride(tmp_path, monkeypatch):
    _fake_run(tmp_path, monkeypatch, frame_dt_fs=10.0)
    res = run_analysis(tmp_path, AnalysisOptions(skip_frames=2, energy=False, bonds=False, spacegroup=False))
    t = res.tables["temperature"]
    assert t["skipped"] == 20 and t["n"] == 80  # 2 frames x 10 fs, not 2 points
    assert res.tables["timeseries_stats"]["temperature"]["n"] == 80


def test_skip_falls_back_to_points_without_a_time_axis():
    values, skipped = report._skip_series([1.0, 2.0, 3.0, 4.0], None, 2, 10.0)
    assert values.tolist() == [3.0, 4.0] and skipped == 2
    values, skipped = report._skip_series([1.0, 2.0, 3.0, 4.0], [0.0, 1.0, 2.0, 3.0], 1, None)
    assert values.tolist() == [2.0, 3.0, 4.0] and skipped == 1


# ---------------- M5: structure factor ----------------
def test_structure_factor_of_fcc_peaks_near_two_pi_over_d111():
    sc = bulk("Cu", "fcc", a=A_CU, cubic=True).repeat((4, 4, 4))
    rmax = 0.5 * sc.cell.lengths().min()
    acc = compute.RDFAccumulator([("Cu", "Cu")], rmax=rmax, nbins=400)
    acc.add(sc)
    rr = acc.result(("Cu", "Cu"))
    got = LO.structure_factor(rr["r"], rr["g"], len(sc) / sc.get_volume())
    q, s = np.array(got["q_1_A"]), np.array(got["s_q"])
    window = (q > 2.0) & (q < 4.0)
    peak = q[window][np.argmax(s[window])]
    assert peak == pytest.approx(2 * np.pi / (A_CU / np.sqrt(3)), abs=0.25)  # (111) at 3.01 1/Å; (200) at 3.48 merges in
    assert got["reliable_above_q"] == pytest.approx(2 * np.pi / rr["r"][-1])


def test_sq_is_written_with_provenance_and_the_q_min_caveat(tmp_path, monkeypatch):
    _fake_run(tmp_path, monkeypatch, frame_dt_fs=10.0)
    quiet = dict(energy=False, bonds=False, spacegroup=False)
    res = run_analysis(tmp_path, AnalysisOptions(rdf=True, structure_factor=True, **quiet))
    assert "sq" in res.figures and Path(res.figures["sq"]).is_file()
    sq = json.loads((tmp_path / "analysis" / "sq.json").read_text(encoding="utf-8"))
    assert set(sq) == {"Ar-Ar", "_provenance"}
    prov = sq["_provenance"]
    assert prov["n_atoms"] == 4 and prov["rmax_A"] == res.tables["rdf"]["rmax"] == pytest.approx(4.0) and "Faber-Ziman" in prov["convention"]
    assert prov["number_density_A3"] == pytest.approx(4 / 8.0 ** 3)
    t = res.tables["structure_factor"]
    assert t["pairs"] == ["Ar-Ar"] and t["reliable_above_q_1_A"] == pytest.approx(2 * np.pi / 4.0)
    txt = res.summary_text()
    assert f"q < 2π/rmax = {2 * np.pi / 4.0:.2f}" in txt
    res2 = run_analysis(tmp_path, AnalysisOptions(structure_factor=True, **quiet))
    assert "sq" not in res2.figures and any("--rdf" in n for n in res2.notes)
    # a molecule has no cell volume to normalize g(r) with, so no S(q)
    monkeypatch.undo()
    d = _copy("dftb_md_water_generated", tmp_path)
    res3 = run_analysis(d, AnalysisOptions(rdf=True, structure_factor=True, **_quiet()))
    assert "sq" not in res3.figures and any("周期系" in n or "periodic" in n for n in res3.notes)


# ---------------- M6 / M7: metals and spin-polarised bands ----------------
def _metal_bands():
    k = np.linspace(0, 1, 21)
    return np.stack([-3 + 0 * k, -1 + 2 * k, 4 + 0 * k], axis=1)  # the middle band crosses E_F = 0 at k = 0.5


def test_gap_details_reports_a_metal_without_vbm_or_cbm():
    g = B.gap_details(_metal_bands(), 0.0)
    assert g["gap_ev"] == 0.0 and g["metal"] is True and g["n_bands_crossing"] == 1
    assert "vbm_ev" not in g and "cbm_ev" not in g
    k = np.linspace(0, 1, 21)
    insulator = np.stack([-3 + 0 * k, -1 - 0.5 * k, 2 + 0.5 * k], axis=1)
    g2 = B.gap_details(insulator, 0.0)
    assert g2["metal"] is False and g2["gap_ev"] == pytest.approx(3.0)
    # a crossing in the second spin channel alone also makes it a metal
    g3 = B.gap_details(insulator, 0.0, energies_down_ev=_metal_bands())
    assert g3["metal"] is True


def _dftb_band_dir(tmp_path: Path, energies, energies_down=None) -> Path:
    d = tmp_path / "run"; (d / "bands").mkdir(parents=True)
    (d / "dftb_in.hsd").write_text("Geometry = {}\n", encoding="utf-8")
    (d / "detailed.out").write_text("Fermi level:                         0.0000000000 H            0.0000 eV\n", encoding="utf-8")
    nk = len(energies)
    kpts = [[0.0, 0.0, float(k) / (nk - 1) * 0.5] for k in range(nk)]
    (d / "bands" / "kpath.json").write_text(json.dumps({"kpts": kpts, "labels": [[0, "G"], [nk - 1, "X"]],
                                                          "cell": np.diag([5.0, 5.0, 5.0]).tolist()}), encoding="utf-8")
    lines = []
    for spin, block in ((1, energies), (2, energies_down)):
        if block is None:
            continue
        for ik, row in enumerate(block, start=1):
            lines.append(f" KPT {ik:12d}  SPIN {spin:12d}  KWEIGHT    1.0000000000000000")
            lines += [f"{b + 1:6d} {e:10.4f}  {2.0 if e <= 0 else 0.0:.5f}" for b, e in enumerate(row)]
            lines.append("")
    (d / "bands" / "band.out").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return d


def test_report_says_a_band_crosses_the_fermi_level(tmp_path):
    d = _dftb_band_dir(tmp_path, _metal_bands())
    res = run_analysis(d, AnalysisOptions(**_quiet()))
    b = res.tables["bands"]
    assert b["metal"] is True and b["gap_ev"] == 0.0 and b["n_spins"] == 1 and b["n_kpoints"] == 21
    assert "vbm_ev" not in b["gap_details"]
    txt = res.summary_text()
    assert "フェルミ準位を横切るバンド" in txt and "effective_mass" not in res.tables


def test_dftb_spin_polarised_band_out_reads_both_channels(tmp_path):
    k = np.linspace(0, 1, 11)
    up = np.stack([-2 - k, 1 + k], axis=1)
    down = np.stack([-1.5 - k, 2 + k], axis=1)
    d = _dftb_band_dir(tmp_path, up, down)
    bd = B.load_bands(d, "dftbplus", 0.0)
    assert bd.n_spins == 2 and np.allclose(bd.energies_ev, up) and np.allclose(bd.energies_down_ev, down)
    res = run_analysis(d, AnalysisOptions(**_quiet()))
    b = res.tables["bands"]
    assert b["n_spins"] == 2 and b["n_bands"] == 2 and b["gap_ev"] == pytest.approx(2.5)  # VBM -1.5 (down), CBM 1.0 (up)
    assert "スピン 2" in res.summary_text() and "bands" in res.figures
    # the same file without the SPIN field: spin is the outer loop, so 2 x nk blocks
    text = (d / "bands" / "band.out").read_text(encoding="utf-8")
    import re
    (d / "bands" / "band.out").write_text(re.sub(r"SPIN\s+\d+", "", text), encoding="utf-8")
    bd2 = B.load_bands(d, "dftbplus", 0.0)
    assert bd2.n_spins == 2 and np.allclose(bd2.energies_down_ev, down)


def test_eigenval_with_ispin_2_reads_both_columns(tmp_path):
    p = tmp_path / "EIGENVAL"
    p.write_text("h\nh\nh\nh\nh\n  8  2  2\n\n 0.0 0.0 0.0 0.5\n 1 -5.0 -4.5 1.0 1.0\n 2 6.0 1.0 0.0 0.0\n\n"
                 " 0.5 0.0 0.5 0.5\n 1 -4.0 -3.5 1.0 1.0\n 2 5.0 -0.5 0.0 0.0\n", encoding="utf-8")
    up, down = B._vasp_eigenval(p)
    assert np.allclose(up, [[-5, 6], [-4, 5]]) and np.allclose(down, [[-4.5, 1.0], [-3.5, -0.5]])
    g = B.gap_details(np.array(up), 0.0, energies_down_ev=np.array(down))
    assert g["metal"] is True  # the spin-down band 2 crosses E_F


def test_qe_spin_up_and_down_sections_are_split(tmp_path):
    def block(k, e1, e2):
        return f"\n          k = 0.0000 0.0000 {k:.4f} (   100 PWs)   bands (ev):\n\n   {e1:8.4f}  {e2:8.4f}\n"
    txt = ("\n ------ SPIN UP ------------\n" + block(0.0, -5.0, 1.0) + block(0.5, -4.0, 2.0)
           + "\n ------ SPIN DOWN ----------\n" + block(0.0, -4.5, 1.5) + block(0.5, -3.5, 2.5)
           + "\n     the Fermi energy is     0.0000 ev\n")
    p = tmp_path / "output.log"; p.write_text(txt, encoding="utf-8")
    up, down = B._qe_bands_out(p)
    assert up == [[-5.0, 1.0], [-4.0, 2.0]] and down == [[-4.5, 1.5], [-3.5, 2.5]]


# ---------------- M8 / L4: thermochemistry mode selection ----------------
WATER = Atoms("OHH", positions=[[0, 0, 0], [0.96, 0, 0], [-0.24, 0.93, 0]])
IDEAL = dict(model="ideal_gas", temperatures_k=(298.15,), pressure_pa=101325.0, symmetry_number=2, geometry="nonlinear", spin=0.0)


def test_imaginary_mode_stays_among_the_used_modes():
    freqs = [-500.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 1600.0, 3600.0]
    t = compute_thermo(freqs, WATER, -100.0, ThermoOptions(**IDEAL))
    assert not t["computed"] and t["n_imaginary_used"] == 1 and sorted(t["modes_used_cm1"]) == [-500.0, 1600.0, 3600.0]
    assert any("虚振動" in r or "imaginary" in r for r in t["reasons"])
    t2 = compute_thermo(freqs, WATER, -100.0, ThermoOptions(imaginary="ignore", **IDEAL))
    assert t2["computed"] and t2["n_imaginary_used"] == 1 and t2["warnings"] == []


def test_a_used_mode_below_50_cm1_is_listed():
    freqs = [3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 45.0, 1600.0, 3600.0]
    t = compute_thermo(freqs, WATER, -100.0, ThermoOptions(**IDEAL))
    assert t["computed"] and len(t["warnings"]) == 1 and "45.0" in t["warnings"][0]
    assert t["reasons"] == t["warnings"]
    assert any("45.0" in line for line in summary_lines(t))


def test_missing_ase_thermo_class_gives_a_reason(monkeypatch):
    import ase.thermochemistry as tc

    monkeypatch.delattr(tc, "QuasiHarmonicThermo", raising=False)
    o = ThermoOptions(model="quasi_harmonic", temperatures_k=(300.0,), exclude_lowest=6, qh_cutoff_cm1=100.0)
    t = compute_thermo([3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 1600.0, 3600.0, 3700.0], WATER, -100.0, o)
    assert not t["computed"] and "QuasiHarmonicThermo" in t["reasons"][0] and "3.28" in t["reasons"][0]


def test_steinhardt_without_spherical_harmonics_gives_a_reason(monkeypatch):
    import scipy.special

    monkeypatch.delattr(scipy.special, "sph_harm_y", raising=False)
    monkeypatch.delattr(scipy.special, "sph_harm", raising=False)
    with pytest.raises(LO.LocalOrderError, match="SciPy 1.15"):
        LO.steinhardt(bulk("Cu", "fcc", a=A_CU), 3.0)


# ---------------- M10 / M12 / L2 / L6: readers and notes ----------------
def test_cp2k_without_output_log_does_not_raise(tmp_path):
    (tmp_path / "cp2k.inp").write_text("&GLOBAL\n PROJECT adit\n RUN_TYPE ENERGY\n&END GLOBAL\n", encoding="utf-8")
    r = readers_extra.read_cp2k(tmp_path)
    assert r.code == "cp2k" and len(r.frames) == 0 and r.energies_ev == []


def test_crest_table_with_an_origin_column(tmp_path):
    p = tmp_path / "crest.log"
    p.write_text(" Erel/kcal        Etot weight/tot  conformer     set   degen     origin\n"
                 "       1   0.000   -14.44226    0.13456    0.26912       1       2     mtd1\n"
                 "       2   0.001   -14.44226    0.13456                                mtd3\n"
                 "       3   0.019   -14.44223    0.13036    0.13036       2       1     mtd2\n\n"
                 "T /K                                  :   298.15\n", encoding="utf-8")
    assert crest.read_log_degeneracy(p) == {1: 2, 2: 1}


def test_gamess_drops_the_modes_it_names(tmp_path):
    p = tmp_path / "gamess.out"
    p.write_text("     MODES    1 TO    5 ARE TAKEN AS ROTATIONS AND TRANSLATIONS.\n"
                 "       FREQUENCY:         1.20        0.80        0.00        0.00        0.00\n"
                 "    IR INTENSITY:      0.00000     0.00000     0.00000     0.00000     0.00000\n"
                 "       FREQUENCY:       650.00     2300.00     2400.00\n"
                 "    IR INTENSITY:      0.10000     0.20000     0.30000\n", encoding="utf-8")
    out = read_gamess(p)
    assert out.frequencies_cm1 == [650.0, 2300.0, 2400.0]
    assert len(out.ir_intensities) == 3 and out.ir_intensities[0] == pytest.approx(0.1 * 42.2561, rel=1e-3)


def test_hbond_note_warns_about_the_gromacs_angle_definition():
    note = criteria_note()
    assert "H–D–A" in note and ("同じ本数にはなりません" in note or "same count" in note)


# ---------------- M11: block error by default ----------------
def test_random_walk_gets_a_block_error_by_default():
    rng = np.random.default_rng(0)
    pos = np.cumsum(rng.normal(0, 0.1, (1000, 32, 3)), axis=0)
    a = compute.msd_analysis(pos, None, None, symbols=["Ar"] * 32, dt_fs=1.0, remove_drift=False)
    e = a["error"]
    assert e["fit_mode"] == "fraction_of_each_block" and e["n_blocks"] == 5
    assert e["d_err_cm2_s"] is not None and np.isfinite(e["d_err_cm2_s"]) and e["d_err_cm2_s"] > 0
    assert all(b["fit_range_fs"] == [pytest.approx(0.1 * b["duration_fs"]), pytest.approx(0.5 * b["duration_fs"])] for b in e["blocks"])
    assert a["D_cm2_s"] == pytest.approx(0.1 ** 2 / 2 * 0.1, rel=0.5)  # sigma^2 / (2 dt), Å²/fs -> cm²/s


# ---------------- M14 / testsci H1: the heavy job ----------------
def test_worker_source_comes_from_package_data():
    assert H.worker_source() == Path(W.__file__).read_text(encoding="utf-8")


def test_write_job_fills_dt_trajectory_and_cell(tmp_path):
    written = H.write_job(tmp_path, dt_fs=5.0, trajectory="geo_end.xyz", cell=np.diag([10.0, 11.0, 12.0]))
    assert {p.name for p in written} == {H.WORKER_FILE, H.JOB_FILE, H.CELL_FILE}
    job = (tmp_path / H.JOB_FILE).read_text(encoding="utf-8")
    assert "msd_worker.py geo_end.xyz --dt 5 " in job and f"--cell {H.CELL_FILE}" in job
    assert W.read_cell(tmp_path / H.CELL_FILE).tolist() == np.diag([10.0, 11.0, 12.0]).tolist()
    blank = H.write_job(tmp_path / "b", dt_fs=None)
    text = (tmp_path / "b" / H.JOB_FILE).read_text(encoding="utf-8")
    assert "<軌跡のファイル>" in text and "--dt <" in text and len(blank) == 2


def test_cli_write_msd_job_uses_the_run_data(tmp_path, capsys):
    from adit.analysis.cli import main

    d = _dftb_npt_dir(tmp_path, [10.0] * 12)  # periodic DFTB+ run: geo_end.xyz carries no lattice
    assert main([str(d), "--write-msd-job"]) == 0
    job = (d / H.JOB_FILE).read_text(encoding="utf-8")
    assert "msd_worker.py geo_end.xyz --dt 1 " in job and f"--cell {H.CELL_FILE}" in job
    assert (d / H.CELL_FILE).is_file() and (d / H.WORKER_FILE).is_file()
    assert str(d / H.JOB_FILE) in capsys.readouterr().out
    out = tmp_path / "msd_vanhove.json"
    W.main([str(d / "geo_end.xyz"), "--dt", "1", "--cell", str(d / H.CELL_FILE), "--taus", "3", "--out", str(out)])
    got = json.loads(out.read_text(encoding="utf-8"))
    assert got["warnings"] == [] and got["cell_A"] == np.diag([10.0, 10.0, 10.0]).tolist() and got["frames"] == 12
    W.main([str(d / "geo_end.xyz"), "--dt", "1", "--taus", "3", "--out", str(out)])
    got = json.loads(out.read_text(encoding="utf-8"))
    assert got["cell_A"] is None and got["warnings"] and "--cell" in got["warnings"][0]


def test_heavy_job_from_the_report_names_the_trajectory_and_cell(tmp_path):
    d = _dftb_npt_dir(tmp_path, [10.0] * 12)
    res = run_analysis(d, AnalysisOptions(msd=True, vanhove=5, heavy_limit_seconds=0.0, out_dir=tmp_path / "out", **_quiet()))
    assert H.CELL_FILE in {Path(p).name for p in res.tables["vanhove_job"]["files"]}
    job = (d / H.JOB_FILE).read_text(encoding="utf-8")
    assert "geo_end.xyz --dt 1 " in job and f"--cell {H.CELL_FILE}" in job


# ---------------- S1 / S3: barostat cells and variable-cell unwrapping ----------------
def _dftb_npt_dir(tmp_path: Path, lengths) -> Path:
    d = tmp_path / "npt"; d.mkdir()
    (d / "dftb_in.hsd").write_text("Geometry = GenFormat { <<< 'geometry.gen' }\n", encoding="utf-8")
    (d / "dftb_pin.hsd").write_text("Driver = VelocityVerlet {\n  TimeStep [fs] = 1.0\n  MDRestartFrequency = 1\n}\n", encoding="utf-8")
    (d / "geometry.gen").write_text("2  S\n Si\n 1 1 0.0 0.0 0.0\n 2 1 1.35 1.35 1.35\n 0.0 0.0 0.0\n"
                                    f" {lengths[0]} 0.0 0.0\n 0.0 {lengths[0]} 0.0\n 0.0 0.0 {lengths[0]}\n", encoding="utf-8")
    md, xyz = [], []
    for k, a in enumerate(lengths):
        md += [f"MD step: {k}", "Lattice vectors (A)",  # written per step only with a barostat (mainio.F90)
               f"  {a:.8E}  0.00000000E+00  0.00000000E+00", f"  0.00000000E+00  {a:.8E}  0.00000000E+00",
               f"  0.00000000E+00  0.00000000E+00  {a:.8E}",
               f"Volume:   {a ** 3 / 0.148:.6E} au^3   {a ** 3:.6E} A^3",
               "Potential Energy:  -3.9798793068 H  -108.2980 eV", "MD Kinetic Energy:  0.0028501338 H  0.0776 eV",
               "Total MD Energy:  -3.9770291730 H  -108.2205 eV", "MD Temperature:  0.0009500446 au  300.0000 K"]
        x = (0.95 + 0.1 * k) % 1.0 * a  # the first atom walks through the boundary while the box changes
        xyz.append(f"2\nMD iter: {k}\nSi {x:.8f} 0.0 0.0\nSi 1.35 1.35 1.35\n")
    (d / "md.out").write_text("\n".join(md) + "\n", encoding="utf-8")
    (d / "geo_end.xyz").write_text("".join(xyz), encoding="utf-8")
    return d


def test_dftb_barostat_cells_are_attached_per_frame(tmp_path):
    lengths = [10.0, 9.9, 9.8, 9.7]
    r = readers.load_run(_dftb_npt_dir(tmp_path, lengths))
    assert len(r.frames) == 4 and all(r.frames[k].cell[0, 0] == pytest.approx(lengths[k]) for k in range(4))
    assert all(r.frames[k].pbc.all() for k in range(4))
    assert any("Lattice vectors" in n for n in r.notes)


def test_variable_cell_is_unwrapped_in_fractional_coordinates(tmp_path):
    lengths = [10.0, 9.9, 9.8, 9.7]
    r = readers.load_run(_dftb_npt_dir(tmp_path, lengths))
    acc = compute.UnwrapAccumulator(len(r.frames))
    for fr in r.frames:
        acc.add(fr)
    pos, _ = acc.result()
    assert acc.variable_cell
    expected = [(0.95 + 0.1 * k) * lengths[k] for k in range(4)]  # unwrapped fraction x current cell length
    assert pos[:, 0, 0] == pytest.approx(expected)
    res = run_analysis(tmp_path / "npt", AnalysisOptions(msd=True, **_quiet()))
    assert res.tables["msd"]["unwrap_check"]["variable_cell"] is True
    assert any("分率座標" in n or "fractional" in n for n in res.notes)


def test_fixed_cell_unwrapping_is_unchanged():
    cell = np.diag([10.0, 10.0, 10.0])
    frames = [Atoms("Ar", positions=[[9.5 + 0.0, 0, 0]], cell=cell, pbc=True), Atoms("Ar", positions=[[0.4, 0, 0]], cell=cell, pbc=True),
              Atoms("Ar", positions=[[1.3, 0, 0]], cell=cell, pbc=True)]
    acc = compute.UnwrapAccumulator(3)
    for fr in frames:
        acc.add(fr)
    pos, _ = acc.result()
    assert pos[:, 0, 0] == pytest.approx([9.5, 10.4, 11.3]) and not acc.variable_cell
