import json
import math
from pathlib import Path

import numpy as np
import pytest

from adit.analysis import AnalysisOptions, run_analysis
from adit.analysis import crest
from adit.analysis.free_energy import KB_EV_K
from adit.analysis.readers import HARTREE_EV

E0 = -33.88484238
REL_KCAL = [0.0, 0.195, 0.580]
KCAL_PER_EH = HARTREE_EV * crest.EV_KCAL_MOL


def _xyz_block(energy_eh: float, shift: float) -> str:
    # C-O-H-H: rotate the hydroxyl hydrogen so the heavy-atom RMSD stays 0 while the geometry changes
    atoms = [("C", 0.0, 0.0, 0.0), ("O", 1.43, 0.0, 0.0), ("H", -0.5, 0.9, 0.0), ("H", 1.8, 0.9 * math.cos(shift), 0.9 * math.sin(shift))]
    return f"{len(atoms)}\n  {energy_eh:18.8f}\n" + "".join(f"{s} {x:.6f} {y:.6f} {z:.6f}\n" for s, x, y, z in atoms)


def write_crest(d: Path, *, energies_file: bool = True, log: bool = True) -> Path:
    d.mkdir(parents=True, exist_ok=True)
    es = [E0 + r / KCAL_PER_EH for r in REL_KCAL]
    (d / "crest_conformers.xyz").write_text("".join(_xyz_block(e, 0.5 * k) for k, e in enumerate(es)), encoding="utf-8")
    if energies_file:
        (d / "crest.energies").write_text("".join(f"  {k + 1}  {r:12.3f}\n" for k, r in enumerate(REL_KCAL)), encoding="utf-8")
    if log:
        (d / "crest.log").write_text(
            "       Erel/kcal        Etot weight/tot  conformer     set   degen     origin\n"
            "       1   0.000   -33.88484    0.06131    0.30454       1       5\n"
            "       2   0.004   -33.88483    0.06088\n"
            "       6   0.195   -33.88452    0.04409    0.13066       2       3\n"
            "       9   0.580   -33.88391    0.02307    0.06891       3       3\n"
            "\nT /K                                  :   298.15\n", encoding="utf-8")
    return d


def test_read_ensemble_energies_and_log(tmp_path):
    d = write_crest(tmp_path / "crest")
    ens = crest.read_ensemble(d / "crest_conformers.xyz")
    assert len(ens) == 3 and ens[0][0] == pytest.approx(E0) and len(ens[0][1]) == 4
    assert crest.read_energies(d / "crest.energies") == {1: 0.0, 2: 0.195, 3: 0.58}
    assert crest.read_log_degeneracy(d / "crest.log") == {1: 5, 2: 3, 3: 3}
    assert crest.read_log_degeneracy(d / "missing.log") == {}
    assert crest.find_crest_dir(tmp_path) == d and crest.find_crest_dir(d) == d
    assert crest.find_crest_dir(tmp_path / "nowhere") is None


def test_boltzmann_weights_formula():
    de = np.array([0.0, 0.01, 0.05])
    g = np.array([2, 1, 1])
    w = crest.boltzmann_weights(de, g, 300.0)
    x = g * np.exp(-de / (KB_EV_K * 300.0))
    assert w == pytest.approx(x / x.sum()) and w.sum() == pytest.approx(1.0)
    with pytest.raises(crest.CrestError):
        crest.boltzmann_weights(de, g, 0.0)


def test_analyze_without_temperature_has_no_weights(tmp_path):
    d = write_crest(tmp_path / "crest")
    out = tmp_path / "analysis"; out.mkdir()
    t = crest.analyze_crest(d, out)
    rows = t["conformers"]
    assert [r["relative_kcal_mol"] for r in rows] == pytest.approx(REL_KCAL, abs=1e-5)
    assert [r["degeneracy"] for r in rows] == [5, 3, 3]
    assert all("weight" not in r for r in rows) and "population_lowest" not in t
    assert all(r["rmsd_heavy_A"] == pytest.approx(0.0, abs=1e-9) for r in rows)
    assert rows[1]["relative_kj_mol"] == pytest.approx(rows[1]["relative_kcal_mol"] * crest.EV_KJ_MOL / crest.EV_KCAL_MOL)
    assert any("温度" in r or "temperature" in r for r in t["reasons"])
    assert Path(t["files"]["table"]).is_file() and Path(t["figure"]).is_file()
    head = Path(t["files"]["table"]).read_text(encoding="utf-8").splitlines()[0]
    assert "weight" in head and "rmsd_heavy_A" in head


def test_analyze_with_temperature(tmp_path):
    d = write_crest(tmp_path / "crest", energies_file=False, log=False)
    out = tmp_path / "analysis"; out.mkdir()
    t = crest.analyze_crest(d, out, 298.15)
    rows = t["conformers"]
    assert [r["degeneracy"] for r in rows] == [1, 1, 1]
    de = np.array([r["relative_ev"] for r in rows])
    x = np.exp(-de / (KB_EV_K * 298.15))
    assert [r["weight"] for r in rows] == pytest.approx(x / x.sum())
    assert t["population_lowest"] == pytest.approx(rows[0]["weight"])
    assert "k_B = 8.617333262e-5" in t["weight_formula"]
    lines = crest.summary_lines(t)
    assert "298.15" in lines[0] and len(lines) >= 5


def test_run_analysis_and_cli_on_a_crest_directory(tmp_path, capsys):
    d = write_crest(tmp_path / "crest")
    res = run_analysis(d, AnalysisOptions(conformer_temperature_k=300.0))
    assert res.tables["crest_conformers"]["temperature_k"] == 300.0
    assert Path(res.figures["crest_conformers"]).is_file()
    js = json.loads((d / "analysis" / "summary.json").read_text(encoding="utf-8"))
    assert len(js["tables"]["crest_conformers"]["conformers"]) == 3 and "weight" in js["tables"]["crest_conformers"]["conformers"][0]
    assert js["figures"]["crest_conformers"] == "analysis/crest_conformers.png"
    from adit.analysis.cli import main as analyze
    assert analyze([str(d), "--conformer-temperature", "298.15"]) == 0
    assert "298.15" in capsys.readouterr().out
    res = run_analysis(d)
    assert "weight" not in res.tables["crest_conformers"]["conformers"][0]


def test_parent_directory_with_crest_subdirectory(tmp_path):
    from adit.analysis.collections import detect_collection
    d = write_crest(tmp_path / "set" / "crest")
    assert detect_collection(d) == "crest"
    (tmp_path / "set" / "conformers.json").write_text(json.dumps({"conformers": []}), encoding="utf-8")
    res = run_analysis(tmp_path / "set")
    assert "conformers" in res.tables and "crest_conformers" in res.tables


def test_gui_section_and_fields_round_trip(tmp_path):
    from adit.gui import analysis_fields as AF
    d = write_crest(tmp_path / "crest")
    res = run_analysis(d, AnalysisOptions(conformer_temperature_k=250.0))
    secs = {s.key: s for s in AF.result_sections(res)}
    assert "crest_conformers" in secs and len(secs["crest_conformers"].rows) == 3 and len(secs["crest_conformers"].columns) == 6
    f = AF.fields_from_options(AnalysisOptions(conformer_temperature_k=250.0))
    assert f["conformer_temperature"] == "250"
    assert AF.options_from_fields(f).conformer_temperature_k == 250.0
    assert AF.options_from_fields({"conformer_temperature": ""}).conformer_temperature_k == 0.0
