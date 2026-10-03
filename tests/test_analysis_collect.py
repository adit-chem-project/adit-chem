
import json
import shutil
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import write

from adit.analysis import AnalysisOptions, run_analysis
from adit.analysis import neb
from adit.analysis.readers import detect_code, load_run
from adit.results import summarize_run

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"


def _copy(name: str, dest: Path) -> Path:
    d = dest / name
    shutil.copytree(EX / name, d, ignore=shutil.ignore_patterns("analysis", "skf"))
    return d


def _mlip_dir(d: Path, task: str, results: dict, *, frames=None, md_log: str = "", opt_log: str = "") -> Path:
    d.mkdir(parents=True, exist_ok=True)
    atoms = Atoms("Ar8", positions=np.array([[x, y, z] for x in (0.0, 5.0) for y in (0.0, 5.0) for z in (0.0, 5.0)]),
                  cell=[10, 10, 10], pbc=True)
    (d / "mlip_settings.json").write_text(json.dumps({"model_family": "mace_mp", "model": "", "task": task,
                                                      "md": {"timestep_fs": 1.0, "dump_interval": 2, "steps": 8}}), encoding="utf-8")
    (d / "results.json").write_text(json.dumps({"task": task, "model_family": "mace_mp", "model": "(package default)",
                                                "versions": {"ase": "3.29.0", "mace-torch": "0.3.14"}, **results}), encoding="utf-8")
    write(d / "structure.extxyz", atoms, format="extxyz")
    write(d / "final.extxyz", atoms, format="extxyz")
    if frames:
        write(d / "trajectory.extxyz", frames, format="extxyz")
    if md_log:
        (d / "md.log").write_text(md_log, encoding="utf-8")
    if opt_log:
        (d / "opt.log").write_text(opt_log, encoding="utf-8")
    (d / "output.log").write_text("adit-mlip: ase 3.29.0, mace-torch 0.3.14\nadit-mlip: done {}\n", encoding="utf-8")
    return d


def _ar8(shift: float) -> Atoms:
    a = Atoms("Ar8", positions=np.array([[x, y, z] for x in (0.0, 5.0) for y in (0.0, 5.0) for z in (0.0, 5.0)]) + shift,
              cell=[10, 10, 10], pbc=True)
    a.calc = SinglePointCalculator(a, energy=-1.0 - shift, forces=np.zeros((8, 3)))
    return a


def test_mlip_single_point_is_detected_and_read(tmp_path):
    d = _mlip_dir(tmp_path / "sp", "single_point",
                  {"energy_ev": -12.5, "fmax_ev_per_ang": 0.031, "stress_voigt_ev_per_ang3": [-0.01, -0.01, -0.01, 0.0, 0.0, 0.0]})
    assert detect_code(d) == "mlip"
    r = load_run(d)
    assert r.code == "mlip" and r.energies_ev == [-12.5] and len(r.final) == 8
    t = r.extra_tables["mlip"]
    assert t["model_family"] == "mace_mp" and t["fmax_ev_per_ang"] == 0.031
    assert t["pressure_gpa"] == pytest.approx(0.01 * 160.21766208)
    assert t["stress_voigt_gpa"][0] == pytest.approx(-0.01 * 160.21766208)
    res = run_analysis(d, AnalysisOptions())
    assert res.code == "mlip" and res.tables["mlip"]["task"] == "single_point"
    js = json.loads((d / "analysis" / "summary.json").read_text(encoding="utf-8"))
    assert js["tables"]["mlip"]["versions"]["ase"] == "3.29.0"
    s = summarize_run(d)
    assert s.finished is True and s.mermin_energy_hartree == pytest.approx(-12.5 / 27.211386245988)


def test_mlip_md_trajectory_feeds_rdf_and_msd(tmp_path):
    frames = [_ar8(0.1 * k) for k in range(6)]
    log = ("Time[ps]      Etot[eV]     Epot[eV]     Ekin[eV]    T[K]\n"
           + "".join(f"{0.002 * k:<10.4f} {-1.0 - 0.1 * k:12.4f} {-1.2 - 0.1 * k:12.4f} {0.2:12.4f}   {300.0 + k:6.1f}\n" for k in range(6)))
    d = _mlip_dir(tmp_path / "md", "molecular_dynamics", {"energy_ev": -1.5, "fmax_ev_per_ang": 0.5}, frames=frames, md_log=log)
    r = load_run(d)
    assert r.frame_source == "trajectory.extxyz" and len(r.frames) == 6
    assert r.energies_ev[0] == pytest.approx(-1.0) and r.temperatures_k[-1] == pytest.approx(305.0)
    assert r.times_fs[1] == pytest.approx(2.0) and r.frame_dt_fs == pytest.approx(2.0)
    res = run_analysis(d, AnalysisOptions(rdf=True, msd=True))
    assert res.tables["rdf"]["n_frames"] == 6 and res.tables["msd"]["n_frames"] == 6
    assert res.tables["temperature"]["mean_k"] == pytest.approx(302.5)
    assert {"rdf", "msd", "energy", "temperature"} <= set(res.figures)


def test_mlip_optimization_and_vibrations(tmp_path):
    opt = ("      Step     Time          Energy          fmax\n"
           "LBFGS:    0 12:00:00      -10.100000        0.400000\n"
           "LBFGS:    1 12:00:01      -10.300000        0.050000\n")
    d = _mlip_dir(tmp_path / "opt", "geometry_optimization", {"energy_ev": -10.3, "fmax_ev_per_ang": 0.05, "converged": True, "steps": 1},
                  opt_log=opt)
    r = load_run(d)
    assert r.energies_ev == [-10.1, -10.3]
    assert summarize_run(d).converged is True
    v = _mlip_dir(tmp_path / "vib", "vibrations", {"energy_ev": -9.0, "fmax_ev_per_ang": 0.002, "frequencies_cm1": [-12.0, 1600.0, 3700.0]})
    rv = load_run(v)
    assert rv.frequencies_cm1 == [-12.0, 1600.0, 3700.0] and rv.force_max_ev_ang == pytest.approx(0.002)
    res = run_analysis(v, AnalysisOptions())
    assert res.tables["frequencies"][0] == -12.0 and "spectrum" in res.figures


def test_phonon_parent_dir_plots_and_asks_for_collect(tmp_path):
    d = _copy("prep_stage3_dftb_phonons_ase", tmp_path)
    res = run_analysis(d, AnalysisOptions())
    assert res.code == "phonons" and "phonon_bands" in res.figures
    t = res.tables["phonon_set"]
    assert t["backend"] == "ase" and t["dim"] == [2, 2, 2] and t["collected"] is True
    assert res.tables["phonopy"]["n_bands"] == 6
    (d / "band.yaml").unlink()
    (d / "total_dos.dat").unlink()
    res2 = run_analysis(d, AnalysisOptions())
    t2 = res2.tables["phonon_set"]
    assert t2["collected"] is False and any("phonon_collect.py" in r for r in t2["reasons"])
    assert "phonon_bands" not in res2.figures


def test_phonon_collect_option_rebuilds_band_yaml(tmp_path):
    d = _copy("prep_stage3_dftb_phonons_ase", tmp_path)
    (d / "band.yaml").unlink()
    res = run_analysis(d, AnalysisOptions(collect=True))
    assert (d / "band.yaml").is_file()
    assert res.tables["phonon_set"]["collected"] is True and "band.yaml" in res.tables["phonon_set"]["collect_summary"]
    assert "phonon_bands" in res.figures


def test_elastic_table_figure_and_collect(tmp_path):
    d = _copy("prep_stage3_dftb_elastic", tmp_path)
    res = run_analysis(d, AnalysisOptions())
    t = res.tables["elastic"]
    assert t["collected"] is True and t["c_gpa"][0][0] == pytest.approx(1133.3374, abs=1e-3)
    assert t["unstrained_pressure_gpa"] == pytest.approx(-4.08690079713397, abs=1e-6)
    assert t["max_abs_residual_gpa"] >= 0
    assert "elastic_stress_strain" in res.figures and Path(t["files"]["fits"]).is_file()
    rows = Path(t["files"]["fits"]).read_text(encoding="utf-8").splitlines()
    assert rows[0].startswith("i,j,c_gpa,intercept_gpa")
    rec = json.loads((d / "elastic.json").read_text(encoding="utf-8"))
    for k in ("c_gpa", "fits", "stress_voigt_gpa"):
        rec.pop(k)
    (d / "elastic.json").write_text(json.dumps(rec), encoding="utf-8")
    res2 = run_analysis(d, AnalysisOptions())
    assert res2.tables["elastic"]["collected"] is False and any("elastic_collect.py" in r for r in res2.tables["elastic"]["reasons"])
    res3 = run_analysis(d, AnalysisOptions(collect=True))
    assert res3.tables["elastic"]["c_gpa"][0][0] == pytest.approx(1133.3374, abs=1e-3)


def test_conformer_table(tmp_path):
    d = _copy("prep_stage2_xtb_conformers", tmp_path)
    res = run_analysis(d, AnalysisOptions())
    t = res.tables["conformers"]
    assert res.code == "conformers" and t["force_field"] == "MMFF94" and t["n_with_energy"] == 3
    rows = {r["dir"]: r for r in t["conformers"]}
    assert rows["conf_001"]["code"] == "xtb"
    assert rows["conf_002"]["relative_kcal_mol"] == pytest.approx(0.0, abs=1e-9)
    assert rows["conf_001"]["relative_ev"] * 1000 == pytest.approx(67.253, abs=0.01)
    assert rows["conf_001"]["ff_relative_kcal_mol"] == pytest.approx(0.0, abs=1e-9)
    csv_rows = Path(t["files"]["table"]).read_text(encoding="utf-8").splitlines()
    assert csv_rows[0].startswith("rank,dir,rdkit_id") and len(csv_rows) == 4
    assert "conformers" in res.figures


def test_neb_images_dirs(tmp_path):
    d = _copy("prep_stage3_dftb_neb_images", tmp_path)
    res = run_analysis(d, AnalysisOptions())
    t = res.tables["neb"]
    assert t["kind"] == "images" and t["n_images"] == 7 and "neb" in res.figures
    js = json.loads((d / "neb.json").read_text(encoding="utf-8"))
    assert t["interpolated"] is True and t["fit_energies_rel_ev"]
    assert any("NEB の最適化はしていません" in r or "no NEB optimization" in r for r in t["reasons"])
    assert t["climbing_image"] is False
    assert any("鞍点" in r or "saddle point" in r for r in t["reasons"])
    assert max(t["energies_rel_ev"]) == pytest.approx(t["forward_barrier_raw_ev"])
    assert len(js["path_length_ang"]) == 7
    for name in ("output.log", "results.tag", "detailed.out"):
        (d / "image_03" / name).unlink()
    res2 = run_analysis(d, AnalysisOptions())
    assert "energies_rel_ev" not in res2.tables["neb"]
    assert any("image_03" in r for r in res2.tables["neb"]["reasons"])
    shutil.rmtree(d / "image_03")
    res3 = run_analysis(d, AnalysisOptions())
    assert "neb" not in res3.tables and any("image_" in n for n in res3.notes)


CP2K_BAND = """ *******************************************************************************
 BAND TYPE                     =                                          CI-NEB
 BAND TYPE OPTIMIZATION        =                                              SD
 STEP NUMBER                   =                                               3
 NUMBER OF NEB REPLICA         =                                               5
 DISTANCES REP =        0.252266        0.229810        0.225933        0.227113
 ENERGIES [au] =      -24.815004      -24.813368      -24.808500      -24.803002
                      -24.810000
 BAND TOTAL ENERGY [au]        =                             -124.05087400000000
 *******************************************************************************
"""


def test_cp2k_band_output(tmp_path):
    d = tmp_path / "band"
    d.mkdir()
    (d / "cp2k.inp").write_text("&GLOBAL\n  PROJECT adit\n  RUN_TYPE BAND\n&END GLOBAL\n", encoding="utf-8")
    (d / "output.log").write_text(CP2K_BAND, encoding="utf-8")
    assert neb.find_neb(d) == {"kind": "cp2k"}
    t = neb.analyze_neb(d, d)
    assert t["kind"] == "cp2k" and t["n_images"] == 5 and t["n_replica_in_output"] == 5 and t["band_type"] == "CI-NEB"
    ha = 27.211386245988
    assert t["energies_rel_ev"][1] == pytest.approx((-24.813368 + 24.815004) * ha)
    assert t["reaction_energy_ev"] == pytest.approx((-24.810000 + 24.815004) * ha)
    assert t["forward_barrier_raw_ev"] == pytest.approx((-24.803002 + 24.815004) * ha)
    assert t["interpolated"] is False and len(t["replica_distances"]) == 4
    assert (d / "neb.png").is_file()


def test_compare_json_is_read_without_the_option(tmp_path):
    d = _copy("prep_stage2_dftb_adsorption", tmp_path)
    res = run_analysis(d, AnalysisOptions())
    t = res.tables["compare"]
    assert res.code == "compare" and len(t["runs"]) == 3 and t["reactions"][0]["name"] == "adsorption"
    assert t["reactions"][0]["delta_e_ev"] == pytest.approx(-488.6238 + 377.7017 + 110.9604, abs=1e-6)
    assert "compare_energy" in res.figures and Path(t["files"]["summary"]).is_file()
    res2 = run_analysis(d, AnalysisOptions(compare=False))
    assert "compare" not in res2.tables


def test_charges_csv(tmp_path):
    d = _copy("prep_stage2_dftb_adsorption", tmp_path) / "molecule"
    res = run_analysis(d, AnalysisOptions())
    q = res.tables["charges"][0]
    p = Path(q["file"])
    assert p.name == "charges_Mulliken.csv" and p.parent.name == "analysis"
    rows = p.read_text(encoding="utf-8").splitlines()
    assert rows[0] == "index,atom,charge_e" and len(rows) == len(q["values"]) + 1
    assert rows[1].split(",")[1] == q["atoms"][0]
    assert float(rows[1].split(",")[2]) == pytest.approx(q["values"][0], abs=1e-6)
    assert sum(float(r.split(",")[2]) for r in rows[1:]) == pytest.approx(q["sum"], abs=1e-5)
