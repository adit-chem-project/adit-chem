
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms, units
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write

from adit.analysis import AnalysisOptions, run_analysis
from adit.analysis import compare, neb, pdos, phonons, symmetry, thermo, uvvis
from adit.analysis.readers import HARTREE_EV, load_run

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"


def _copy(name: str, dest: Path) -> Path:
    d = dest / name
    shutil.copytree(EX / name, d, ignore=shutil.ignore_patterns("analysis"))
    return d


def test_parse_compare_forms():
    r = compare.parse_compare("ads=1:slab_mol,-1:slab,-1:mol; -0.5:C:\\runs\\o2, 1:h2o")
    assert [x.name for x in r] == ["ads", "r2"]
    assert r[0].terms == [(1.0, "slab_mol"), (-1.0, "slab"), (-1.0, "mol")]
    assert r[1].terms == [(-0.5, "C:\\runs\\o2"), (1.0, "h2o")]
    for bad in ("ads=slab", "x:slab", "0:slab", "", "a=1:x; a=1:y"):
        with pytest.raises(compare.CompareError):
            compare.parse_compare(bad)


def test_compare_energies_balance_and_conditions(tmp_path):
    for n in ("xtb_water_generated", "xtb_vib_water_generated", "water_generated", "qe_si_generated"):
        _copy(n, tmp_path)
    rx = compare.parse_compare("same=1:xtb_water_generated,-1:xtb_vib_water_generated;"
                               "cross=1:water_generated,-1:xtb_water_generated; unbal=1:water_generated")
    res = compare.analyze_compare(tmp_path, rx)
    e = {d: load_run(tmp_path / d).energies_ev[-1] for d in ("xtb_water_generated", "xtb_vib_water_generated", "water_generated")}
    same, cross, unbal = res.reactions
    assert same["delta_e_ev"] == pytest.approx(e["xtb_water_generated"] - e["xtb_vib_water_generated"])
    assert same["delta_e_kj_mol"] == pytest.approx(same["delta_e_ev"] * 96.485332123)
    assert same["delta_e_kcal_mol"] == pytest.approx(same["delta_e_ev"] * 23.060547830619)
    assert same["balanced"] and "task.type" in same["differing"] and "method.code" not in same["differing"]
    assert same["delta_g_code_ev"] is None
    assert cross["differing"] == ["method.code"] and cross["n_differing"] == 1
    assert "method.sk_set" in cross["partial"] and "method.gfn" in cross["partial"]
    assert not unbal["balanced"] and unbal["imbalance"] == "H:+2 O:+1"
    runs = {r["dir"]: r for r in res.runs}
    assert runs["water_generated"]["formula"] == "H2O" and runs["water_generated"]["natoms"] == 3
    assert runs["xtb_vib_water_generated"]["code_g_ev"] == pytest.approx(-4.988153848206 * HARTREE_EV)
    for k in ("runs", "reactions", "conditions", "summary"):
        assert Path(res.files[k]).is_file()
    assert Path(res.figures["compare_energy"]).is_file()
    js = json.loads((tmp_path / "compare_summary.json").read_text(encoding="utf-8"))
    assert js["reactions"][2]["imbalance"] == "H:+2 O:+1"
    assert js["files"]["runs"] == "compare_runs.csv" and js["figures"]["compare_energy"] == "compare_energy.png"
    assert all(not Path(r["path"]).is_absolute() for r in js["runs"])
    head = (tmp_path / "compare_runs.csv").read_text(encoding="utf-8").splitlines()[0]
    assert "method.code" in head and "task.type" in head
    from adit.project import load_project
    sp = load_project(tmp_path / "qe_si_generated")
    b = 2 * np.pi * np.linalg.norm(np.linalg.inv(np.array(sp.structure.atoms.cell)).T, axis=1)
    got = json.loads(compare.conditions_of(sp)["kpoints.spacing_2pi_per_A"])
    assert got == pytest.approx(b / np.array(sp.kpoints.mesh), rel=1e-3)


def test_compare_json_and_cli(tmp_path, capsys):
    from adit.analysis.cli import main as analyze
    for n in ("xtb_water_generated", "xtb_vib_water_generated"):
        _copy(n, tmp_path)
    (tmp_path / "compare.json").write_text(json.dumps({"reactions": [{"name": "same", "terms": [
        {"dir": "xtb_water_generated", "nu": 1}, {"dir": "xtb_vib_water_generated", "nu": -1}]}]}), encoding="utf-8")
    assert analyze([str(tmp_path), "--compare"]) == 0
    assert "same" in capsys.readouterr().out
    assert analyze([str(tmp_path), "--compare", "x=1:xtb_water_generated,-1:missing_dir"]) == 0
    assert "missing_dir" in (tmp_path / "compare_runs.csv").read_text(encoding="utf-8")
    empty = tmp_path / "empty"; empty.mkdir()
    assert analyze([str(empty), "--compare"]) == 1
    assert "compare.json" in capsys.readouterr().err


WATER = Atoms("OH2", positions=[[0, 0, 0.119], [0, 0.763, -0.477], [0, -0.763, -0.477]])


def _ig(**kw):
    base = dict(model="ideal_gas", temperatures_k=(298.15,), pressure_pa=100000.0, symmetry_number=2, geometry="nonlinear", spin=0.0)
    base.update(kw)
    return thermo.ThermoOptions(**base)


def test_thermo_ideal_gas_water_from_xtb():
    r = load_run(EX / "xtb_vib_water_generated")
    t = thermo.compute_thermo(r.frequencies_cm1, r.final, r.energies_ev[-1], _ig(), r.thermo)
    assert t["computed"] and t["n_modes_used"] == 3
    top3 = sorted(r.frequencies_cm1)[-3:]
    row = t["rows"][0]
    assert row["zpe_ev"] == pytest.approx(0.5 * sum(top3) * units.invcm)
    assert row["zpe_ev"] == pytest.approx(0.008147568762 * HARTREE_EV, rel=2e-3)
    assert row["h_corr_ev"] - row["u_corr_ev"] == pytest.approx(units.kB * 298.15)
    assert row["s_j_mol_k"] == pytest.approx(46.9633 * 4.184, rel=2e-3)
    assert row["g_total_ev"] == pytest.approx(r.energies_ev[-1] + row["g_corr_ev"])
    assert any(it["key"] == "total_free_energy" for it in t["code_values"])
    lower = thermo.compute_thermo(r.frequencies_cm1, r.final, None, _ig(pressure_pa=1e6))["rows"][0]
    assert row["s_ev_per_k"] - lower["s_ev_per_k"] == pytest.approx(units.kB * np.log(10))


def test_thermo_needs_inputs_and_no_defaults():
    f = [0.0] * 6 + [1600.0, 3600.0, 3700.0]
    t = thermo.compute_thermo(f, WATER, None, thermo.ThermoOptions(model="ideal_gas", temperatures_k=(298.15,)))
    assert not t["computed"] and "圧力" in t["reasons"][0] or "pressure" in t["reasons"][0]
    t = thermo.compute_thermo(f, WATER, None, thermo.ThermoOptions(model="harmonic", temperatures_k=(300.0,)))
    assert not t["computed"]
    t = thermo.compute_thermo([], WATER, None, _ig())
    assert not t["computed"]
    t = thermo.compute_thermo(f, WATER, None, thermo.ThermoOptions(model="harmonic", temperatures_k=(300.0,), exclude_lowest=5))
    assert not t["computed"] and ("0 cm^-1" in t["reasons"][0] or "0 cm⁻¹" in t["reasons"][0])
    t = thermo.compute_thermo(f, WATER, -10.0, thermo.ThermoOptions(model="harmonic", temperatures_k=(300.0, 600.0), exclude_lowest=6))
    assert t["computed"] and t["n_modes_used"] == 3 and len(t["rows"]) == 2
    r0 = t["rows"][0]
    assert r0["f_corr_ev"] == pytest.approx(r0["u_corr_ev"] - 300.0 * r0["s_ev_per_k"]) and r0["f_total_ev"] == pytest.approx(-10.0 + r0["f_corr_ev"])


def test_thermo_imaginary_modes_are_the_users_choice():
    f = [-300.0, 1600.0, 3600.0, 3700.0]
    o = dict(model="harmonic", temperatures_k=(300.0,), exclude_lowest=0)
    t = thermo.compute_thermo(f, WATER, None, thermo.ThermoOptions(**o))
    assert not t["computed"] and t["n_imaginary_used"] == 1 and "ignore" in t["reasons"][0]
    assert not thermo.compute_thermo(f, WATER, None, thermo.ThermoOptions(**o, imaginary="stop"))["computed"]
    t = thermo.compute_thermo(f, WATER, None, thermo.ThermoOptions(**o, imaginary="ignore"))
    assert t["computed"] and t["rows"][0]["zpe_ev"] == pytest.approx(0.5 * (1600 + 3600 + 3700) * units.invcm)
    t = thermo.compute_thermo([-50.0, 1.0, 2.0, 3.0, 4.0, 5.0, 1600.0, 3600.0, 3700.0], WATER, None, _ig())
    assert t["computed"] and t["n_imaginary_all"] == 1 and t["n_imaginary_used"] == 0


def test_thermo_quasi_harmonic_and_msrrho():
    f = [20.0, 1600.0, 3600.0]
    t = thermo.compute_thermo(f, WATER, None, thermo.ThermoOptions(model="quasi_harmonic", temperatures_k=(298.15,), exclude_lowest=0, qh_cutoff_cm1=100.0))
    assert t["computed"] and t["rows"][0]["zpe_ev"] == pytest.approx(0.5 * (100 + 1600 + 3600) * units.invcm)
    assert not thermo.compute_thermo(f, WATER, None, thermo.ThermoOptions(model="quasi_harmonic", temperatures_k=(298.15,), exclude_lowest=0))["computed"]
    t = thermo.compute_thermo(f, WATER, None, thermo.ThermoOptions(model="msrrho", temperatures_k=(298.15,), exclude_lowest=0, msrrho_tau_cm1=35.0))
    h = thermo.compute_thermo(f, WATER, None, thermo.ThermoOptions(model="harmonic", temperatures_k=(298.15,), exclude_lowest=0))
    assert t["computed"] and t["ase_class"] == "MSRRHOThermo"
    assert t["rows"][0]["s_ev_per_k"] < h["rows"][0]["s_ev_per_k"]


def test_thermo_in_run_analysis_and_cli(tmp_path):
    from adit.analysis.cli import main as analyze
    d = _copy("xtb_vib_water_generated", tmp_path)
    assert analyze([str(d), "--thermo", "ideal_gas", "--temperature", "298.15,400", "--pressure", "101325",
                    "--symmetry-number", "2", "--geometry", "nonlinear", "--spin", "0"]) == 0
    js = json.loads((d / "analysis" / "summary.json").read_text(encoding="utf-8"))
    t = js["tables"]["thermo_ase"]
    assert t["computed"] and [r["T_K"] for r in t["rows"]] == [298.15, 400.0] and (d / "analysis" / "thermo.csv").is_file()
    assert analyze([str(d), "--thermo", "ideal_gas"]) == 0
    js = json.loads((d / "analysis" / "summary.json").read_text(encoding="utf-8"))
    assert not js["tables"]["thermo_ase"]["computed"] and js["tables"]["thermo_ase"]["reasons"]


# ---------------- NEB ----------------
A_, B_, LX = 0.8, -0.3, 2.0


def _E(x):
    return A_ * np.sin(np.pi * x / LX) ** 2 + B_ * x / LX


def _dE(x):
    return A_ * np.pi / LX * np.sin(2 * np.pi * x / LX) + B_ / LX


def _images(n=7):
    out = []
    for x in np.linspace(0, LX, n):
        a = Atoms("H", positions=[[x, 0, 0]], cell=[10, 10, 10], pbc=True)
        a.calc = SinglePointCalculator(a, energy=float(_E(x)) - 5.0, forces=[[-float(_dE(x)), 0, 0]])
        out.append(a)
    return out


def test_neb_from_images_barrier():
    t = neb.neb_from_images(_images())
    xs = np.linspace(0, LX, 20001)
    assert t["reaction_energy_ev"] == pytest.approx(B_)
    assert t["forward_barrier_raw_ev"] == pytest.approx(_E(1.0))
    assert t["forward_barrier_ev"] == pytest.approx(_E(xs).max(), abs=0.02)
    assert t["backward_barrier_ev"] == pytest.approx(t["forward_barrier_ev"] - B_)
    assert t["path_A"] == pytest.approx(np.linspace(0, LX, 7))


def test_neb_vasp_dirs_oszicar_only(tmp_path):
    ims = _images()
    for i, a in enumerate(ims):
        d = tmp_path / f"{i:02d}"; d.mkdir()
        write(d / "POSCAR", a, format="vasp")
        e = a.get_potential_energy()
        (d / "OSZICAR").write_text(f"       N       E                     dE             d eps       ncg     rms          rms(c)\n"
                                   f"   1 F= {e:.8E} E0= {e:.8E}  d E =0.000000E+00\n", encoding="utf-8")
    shutil.copy(EX / "vasp_h2o_generated" / "spec.json", tmp_path / "spec.json")
    result = run_analysis(tmp_path)
    t = result.tables["neb"]
    assert t["kind"] == "vasp" and t["n_images"] == 7 and t["interpolated"] is False and t["forward_barrier_ev"] is None
    assert t["forward_barrier_raw_ev"] == pytest.approx(_E(1.0)) and t["reaction_energy_ev"] == pytest.approx(B_)
    assert t["path_A"] == pytest.approx(np.linspace(0, LX, 7)) and Path(t["figure"]).is_file()
    assert not any("vasprun.xml" in note for note in result.notes)
    (tmp_path / "06" / "OSZICAR").write_text("", encoding="utf-8")
    t = neb.analyze_neb(tmp_path, tmp_path)
    assert "energies_rel_ev" not in t and "06" in t["reasons"][0]


def _qe_neb(d: Path):
    x = np.linspace(0, 1, 7)
    e = _E(x * LX)
    (d / "adit.dat").write_text("".join(f"{a:18.10f}{b:18.10f}{0.05:18.10f}\n" for a, b in zip(x, e)), encoding="utf-8")
    xi = np.linspace(0, 1, 101)
    (d / "adit.int").write_text("".join(f"{a:18.10f}{b:18.10f}\n" for a, b in zip(xi, _E(xi * LX))), encoding="utf-8")
    fwd = _E(xi * LX).max()
    (d / "neb.out").write_text(f"\n     activation energy (->) = {fwd:10.6f} eV\n     activation energy (<-) = {fwd - B_:10.6f} eV\n\n", encoding="utf-8")
    return fwd


def test_neb_qe_files_and_run_analysis(tmp_path):
    d = tmp_path / "neb"; d.mkdir()
    fwd = _qe_neb(d)
    res = run_analysis(d)
    t = res.tables["neb"]
    assert res.code == "espresso"
    assert t["kind"] == "espresso" and t["interpolated"] and t["forward_barrier_ev"] == pytest.approx(fwd, abs=1e-9)
    assert t["forward_barrier_raw_ev"] == pytest.approx(_E(1.0)) and t["reaction_energy_ev"] == pytest.approx(B_)
    assert t["backward_barrier_ev"] == pytest.approx(fwd - B_)
    assert t["neb_x_activation"]["forward"]["value_ev"] == pytest.approx(round(fwd, 6)) and t["neb_x_activation"]["forward"]["line"] == "neb.out:2"
    assert Path(res.figures["neb"]).is_file() and "NEB" in res.summary_text()
    (d / "other.dat").write_text("0 0 0\n", encoding="utf-8"); (d / "other.int").write_text("0 0\n", encoding="utf-8")
    assert "reasons" in neb.analyze_neb(d, d)


def test_neb_qe_reports_maximum_steps(tmp_path):
    d = tmp_path / "qe_neb"
    d.mkdir()
    fwd = _qe_neb(d)
    with open(d / "neb.out", "a", encoding="utf-8") as fh:
        fh.write("\n     neb: reached the maximum number of steps\n\n     JOB DONE.\n")
    result = run_analysis(d)
    table = result.tables["neb"]
    assert table["forward_barrier_ev"] == pytest.approx(fwd, abs=1e-9)
    assert table["stopping_reason"] == "maximum_steps"
    assert any("最大反復数" in reason or "maximum number" in reason for reason in table["reasons"])


# ---------------- PDOS ----------------
def _doscar(path: Path, nions: int, e, ef: float, total_cols: list, atom_blocks: list):
    head6 = f"  {e[-1]:.4f}  {e[0]:.4f}  {len(e)}  {ef:.4f}  1.0"
    lines = [f"{nions:4d}{nions:4d}   {1 if atom_blocks else 0}   0", "  1.0  1e-10  1e-10  1e-10  1e-15", "  1.0E-004", "  CAR ", " synthetic", head6]
    lines += [" ".join(f"{v:.6f}" for v in [x, *row]) for x, row in zip(e, total_cols)]
    for blk in atom_blocks:
        lines.append(head6)
        lines += [" ".join(f"{v:.6f}" for v in [x, *row]) for x, row in zip(e, blk)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_pdos_doscar_lorbit11_spin(tmp_path):
    e = np.linspace(-2, 2, 5)
    atoms = Atoms("SiO", positions=[[0, 0, 0], [1.5, 0, 0]], cell=[4, 4, 4], pbc=True)
    write(tmp_path / "POSCAR", atoms, format="vasp")
    (tmp_path / "INCAR").write_text("ISPIN = 2\nLORBIT = 11\n", encoding="utf-8")
    blk = [[v for k in range(9) for v in (k + 1.0, 10.0 * (k + 1))] for _ in e]
    _doscar(tmp_path / "DOSCAR", 2, e, 0.5, [[1, 2, 0, 0]] * 5, [blk, blk])
    d = pdos.read_doscar(tmp_path)
    ch = d["channels"]
    assert d["spin"] == 2 and set(ch) == {("Si", "s"), ("Si", "p"), ("Si", "d"), ("O", "s"), ("O", "p"), ("O", "d")}
    assert ch[("Si", "s")]["up"] == pytest.approx([1.0] * 5) and ch[("Si", "p")]["up"] == pytest.approx([2 + 3 + 4.0] * 5)
    assert ch[("O", "d")]["down"] == pytest.approx([10.0 * (5 + 6 + 7 + 8 + 9)] * 5)
    res = run_analysis(tmp_path)
    t = res.tables["pdos"]
    assert t["shifted_to_fermi"] and t["fermi_ev"] == pytest.approx(0.5) and "Si_p_up" in t["columns"]
    csv0 = (tmp_path / "analysis" / "pdos.csv").read_text(encoding="utf-8").splitlines()
    assert csv0[0].startswith("energy_minus_ef_ev,total_up,total_down") and float(csv0[1].split(",")[0]) == pytest.approx(-2.5)


def test_pdos_doscar_lorbit10_and_missing(tmp_path):
    e = np.linspace(-1, 1, 3)
    write(tmp_path / "POSCAR", Atoms("Cu", cell=[3, 3, 3], pbc=True), format="vasp")
    _doscar(tmp_path / "DOSCAR", 1, e, 0.0, [[1, 0]] * 3, [[[1.0, 2.0, 3.0]] * 3])
    d = pdos.read_doscar(tmp_path)
    assert d["orbitals_per_atom"] == 3 and d["channels"][("Cu", "d")][""] == pytest.approx([3.0] * 3)
    d = pdos.read_doscar(EX / "vasp_h2o_generated")
    assert not d["channels"] and "LORBIT" in d["reasons"][0]
    h = _copy("vasp_h2o_generated", tmp_path)
    assert "pdos" not in run_analysis(h).tables
    assert "LORBIT" in " ".join(run_analysis(h, AnalysisOptions(pdos=True)).tables["pdos"]["reasons"])


def test_pdos_projwfc_qe(tmp_path):
    d = _copy("qe_si_generated", tmp_path)
    e = np.linspace(0, 10, 6)
    def put(name, head, cols):
        (d / name).write_text(head + "\n" + "".join(" ".join(f"{v:.3f}" for v in (x, *c)) + "\n" for x, c in zip(e, cols)), encoding="utf-8")
    put("adit.pdos_atm#1(Si)_wfc#1(s)", "# E (eV)  ldos(E)   pdos(E)", [[1.0, 1.0]] * 6)
    put("adit.pdos_atm#1(Si)_wfc#2(p)", "# E (eV)  ldos(E)   pz(E)     px(E)     py(E)", [[3.0, 1, 1, 1]] * 6)
    put("adit.pdos_atm#2(Si)_wfc#1(s)", "# E (eV)  ldos(E)   pdos(E)", [[1.5, 1.5]] * 6)
    put("adit.pdos_atm#2(Si)_wfc#2(p)", "# E (eV)  ldos(E)   pz(E)     px(E)     py(E)", [[3.0, 1, 1, 1]] * 6)
    put("adit.pdos_tot", "# E (eV)  dos(E)    pdos(E)", [[9.0, 8.5]] * 6)
    ef = load_run(d).fermi_ev
    res = run_analysis(d)
    t = res.tables["pdos"]
    assert t["channels"] == ["Si_p", "Si_s"] and t["fermi_ev"] == pytest.approx(ef) and t["total_source"] == "adit.pdos_tot"
    rows = (d / "analysis" / "pdos.csv").read_text(encoding="utf-8").splitlines()
    assert rows[0] == "energy_minus_ef_ev,total,Si_p,Si_s" and rows[1].split(",")[2:] == ["6", "2.5"]
    assert float(rows[1].split(",")[0]) == pytest.approx(0 - ef, abs=1e-4)
    assert any("highest occupied" in n or "最高被占準位" in n for n in res.notes)
    put("adit.pdos_atm#1(Si)_wfc#1(s)", "# E (eV)  ldosup(E)  ldosdw(E)  pdosup(E)  pdosdw(E)", [[1.0, 2.0, 1.0, 2.0]] * 6)
    assert "reasons" in pdos.read_projwfc(d, 0.0) and "channels" not in pdos.read_projwfc(d, 0.0)


# ---------------- UV-Vis ----------------
ORCA6_TABLE = """
----------------------------------------------------------------------------------------------------
                     ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
----------------------------------------------------------------------------------------------------
     Transition      Energy     Energy  Wavelength fosc(D2)      D2        DX        DY        DZ
                      (eV)      (cm-1)    (nm)                 (au**2)    (au)      (au)      (au)
----------------------------------------------------------------------------------------------------
  0-1A  ->  1-1A    2.544515   20522.9   487.3   0.034716965   0.55690   0.74464   0.00916  -0.04823
  0-1A  ->  2-1A    4.292647   34622.5   288.8   0.015206604   0.14459  -0.00926  -0.34082  -0.16837
  0-1A  ->  3-1A    4.506793   36349.7   275.1   0.144575201   1.30939   1.14340  -0.00032  -0.04505

"""
STDA_TABLE = """
-----------------------------------------------------------------------------
         ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-----------------------------------------------------------------------------
State   Energy  Wavelength   fosc         T2         TX        TY        TZ
        (cm-1)    (nm)                  (au**2)     (au)      (au)      (au)
-----------------------------------------------------------------------------
   1   36108.8    276.9   0.000000000   0.00000   0.00000   0.00000   0.00000
   2   40000.0    250.0   0.100000000   0.80000   0.89443   0.00000   0.00000

"""


def test_uvvis_read_formats(tmp_path):
    p = tmp_path / "output.log"
    p.write_text("dummy\n" + ORCA6_TABLE, encoding="utf-8")
    t = uvvis.read_orca_absorption(p)
    assert [r["energy_ev"] for r in t["transitions"]] == [2.544515, 4.292647, 4.506793]
    assert [r["fosc"] for r in t["transitions"]] == [0.034716965, 0.015206604, 0.144575201]
    assert t["transitions"][0]["label"] == "0-1A -> 1-1A" and t["line"] == "output.log:4" and t["format"].startswith("eV")
    p.write_text(ORCA6_TABLE + STDA_TABLE, encoding="utf-8")
    t = uvvis.read_orca_absorption(p)
    assert t["n_tables"] == 2 and t["transitions"][1]["energy_ev"] == pytest.approx(40000.0 / 8065.543937)
    assert t["transitions"][1]["wavelength_nm"] == 250.0 and t["transitions"][1]["fosc"] == 0.1
    p.write_text(ORCA6_TABLE.replace("ABSORPTION SPECTRUM", "SOC CORRECTED ABSORPTION SPECTRUM"), encoding="utf-8")
    assert uvvis.read_orca_absorption(p) is None


def test_uvvis_broadening_area():
    grid = np.linspace(-50, 60, 200001)
    for shape in ("gauss", "lorentz"):
        y = uvvis.broaden([2.0, 5.0], [0.1, 0.3], 0.4, shape, grid)
        assert np.trapezoid(y, grid) == pytest.approx(0.4, rel=5e-3)
    y = uvvis.broaden([2.0], [1.0], 0.4, "gauss", np.array([2.0, 2.2]))
    assert y[1] / y[0] == pytest.approx(0.5)


def test_uvvis_in_run_analysis_and_cli(tmp_path):
    from adit.analysis.cli import main as analyze
    d = tmp_path / "orca"; d.mkdir()
    (d / "orca.inp").write_text("! B3LYP def2-SVP\n%tddft nroots 3 end\n* xyz 0 1\nO 0 0 0\n*\n", encoding="utf-8")
    (d / "output.log").write_text(ORCA6_TABLE, encoding="utf-8")
    res = run_analysis(d)
    assert res.tables["uvvis"]["n_transitions"] == 3 and res.tables["uvvis"]["broadening"] is None and Path(res.figures["uvvis"]).is_file()
    assert analyze([str(d), "--uvvis", "gauss:0.3"]) == 0
    js = json.loads((d / "analysis" / "summary.json").read_text(encoding="utf-8"))
    assert js["tables"]["uvvis"]["broadening"] == {"shape": "gauss", "fwhm_ev": 0.3} and (d / "analysis" / "uvvis_spectrum.csv").is_file()
    assert analyze([str(d), "--uvvis", "box:0.3"]) == 1


def test_spacegroup_si_and_without_spglib(tmp_path, monkeypatch):
    pytest.importorskip("spglib")
    from ase.build import bulk
    t = symmetry.spacegroup(bulk("Si", "diamond", a=5.43))
    assert [r["number"] for r in t["results"]] == [227, 227, 227] and t["results"][0]["international"] == "Fd-3m"
    assert symmetry.spacegroup(WATER) is None
    d = _copy("qe_si_generated", tmp_path)
    assert run_analysis(d).tables["spacegroup"]["results"][0]["number"] == 227
    monkeypatch.setitem(sys.modules, "spglib", None)
    assert "spglib" in symmetry.spacegroup(bulk("Si", "diamond", a=5.43))["reason"]
    res = run_analysis(d)
    assert "spacegroup" not in res.tables and any("spglib" in n for n in res.notes)


# ---------------- phonopy ----------------
BAND_YAML = """nqpoint: 5
npath: 2
segment_nqpoint:
- 3
- 2
labels:
- [ 'G', 'X' ]
- [ 'X', 'L' ]
reciprocal_lattice:
- [     0.18416500,     0.00000000,     0.00000000 ] # a*
natom: 2
lattice:
- [     5.43000000,     0.00000000,     0.00000000 ] # a
points:
- symbol: Si # 1
  coordinates: [  0.000000000000000,  0.000000000000000,  0.000000000000000 ]
  mass: 28.085500

phonon:
- q-position: [    0.0000000,    0.0000000,    0.0000000 ]
  distance:    0.0000000
  band:
  - # 1
    frequency:   -0.0200000000
    eigenvector:
    - # atom 1
      - [  0.70710678118655,  0.00000000000000 ]
  - # 2
    frequency:   15.5000000000
- q-position: [    0.2500000,    0.0000000,    0.2500000 ]
  distance:    0.0500000
  band:
  - # 1
    frequency:    3.0000000000
  - # 2
    frequency:   15.0000000000
- q-position: [    0.5000000,    0.0000000,    0.5000000 ]
  distance:    0.1000000
  band:
  - # 1
    frequency:    4.0000000000
  - # 2
    frequency:   14.0000000000
- q-position: [    0.5000000,    0.0000000,    0.5000000 ]
  distance:    0.1000000
  band:
  - # 1
    frequency:    4.0000000000
  - # 2
    frequency:   14.0000000000
- q-position: [    0.5000000,    0.5000000,    0.5000000 ]
  distance:    0.1500000
  band:
  - # 1
    frequency:    3.5000000000
  - # 2
    frequency:   14.5000000000
"""


def test_phonopy_band_yaml_and_dos(tmp_path):
    (tmp_path / "band.yaml").write_text(BAND_YAML, encoding="utf-8")
    (tmp_path / "total_dos.dat").write_text("# Sigma = 0.100000\n" + "".join(f"{x:20.10f}{y:20.10f}\n" for x, y in [(0, 0), (5, 0.3), (15, 1.2)]), encoding="utf-8")
    b = phonons.read_band_yaml(tmp_path / "band.yaml")
    assert b["frequency"].shape == (5, 2) and b["segment_nqpoint"] == [3, 2] and b["labels"] == [("G", "X"), ("X", "L")]
    assert phonons._ticks(b) == [(0.0, "G"), (0.1, "X"), (0.15, "L")]
    t = phonons.analyze_phonopy(tmp_path, tmp_path)
    assert t["n_qpoints"] == 5 and t["min_frequency"] == pytest.approx(-0.02) and t["n_negative_values"] == 1 and t["dos_points"] == 3
    assert Path(t["figure_bands"]).is_file() and t["dos_comment"] == "Sigma = 0.100000"
    d = _copy("qe_si_generated", tmp_path)
    (d / "phonopy").mkdir(); (d / "phonopy" / "total_dos.dat").write_text("0 0\n5 1\n", encoding="utf-8")
    res = run_analysis(d)
    assert res.tables["phonopy"]["dos_points"] == 2 and Path(res.figures["phonon_dos"]).is_file()
