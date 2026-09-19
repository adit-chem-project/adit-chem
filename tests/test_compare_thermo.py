import csv
import json
import shutil
from pathlib import Path

import pytest

from adit.analysis import AnalysisOptions, run_analysis
from adit.analysis import compare
from adit.analysis.thermo import ThermoOptions

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"


def _ig(**kw):
    base = dict(model="ideal_gas", temperatures_k=(298.15,), pressure_pa=100000.0, symmetry_number=2, geometry="nonlinear", spin=0.0)
    base.update(kw)
    return ThermoOptions(**base)


def _prepare(tmp_path: Path, name: str, thermo: ThermoOptions | None) -> Path:
    d = tmp_path / name
    shutil.copytree(EX / "xtb_vib_water_generated", d, ignore=shutil.ignore_patterns("analysis"))
    if thermo is not None:
        res = run_analysis(d, AnalysisOptions(thermo=thermo))
        assert res.tables["thermo_ase"]["computed"], res.tables["thermo_ase"].get("reasons")
    return d


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    tmp = tmp_path_factory.mktemp("compare_thermo")
    _prepare(tmp, "a", _ig())
    _prepare(tmp, "b", _ig())
    _prepare(tmp, "two_t", _ig(temperatures_k=(298.15, 400.0)))
    _prepare(tmp, "other_p", _ig(pressure_pa=101325.0))
    _prepare(tmp, "harm", ThermoOptions(model="harmonic", temperatures_k=(298.15,), exclude_lowest=6))
    _prepare(tmp, "none", None)
    return tmp


def _thermo_csv(d: Path) -> dict:
    with open(d / "analysis" / "thermo.csv", newline="", encoding="utf-8") as f:
        return next(csv.DictReader(f))


def test_delta_h_s_g_from_thermo_csv(runs):
    rx = compare.parse_compare("same=1:a,-1:b; scaled=2:a,-1:a")
    res = compare.analyze_compare(runs, rx)
    same, scaled = res.reactions
    assert same["thermo_note"] == "" and len(same["thermo"]) == 1
    t = same["thermo"][0]
    assert t["T_K"] == 298.15 and t["P_Pa"] == 100000.0 and t["g_label"] == "G"
    assert t["delta_h_ev"] == pytest.approx(0.0, abs=1e-9) and t["delta_s_ev_per_k"] == pytest.approx(0.0, abs=1e-12)
    assert t["delta_g_ev"] == pytest.approx(0.0, abs=1e-9)
    row = _thermo_csv(runs / "a")
    s = scaled["thermo"][0]
    assert s["delta_h_ev"] == pytest.approx(float(row["h_total_ev"]))
    assert s["delta_s_ev_per_k"] == pytest.approx(float(row["s_ev_per_k"]))
    assert s["delta_g_ev"] == pytest.approx(float(row["g_total_ev"]))
    assert s["delta_h_kj_mol"] == pytest.approx(s["delta_h_ev"] * compare.EV_KJ_MOL)
    assert s["delta_s_j_mol_k"] == pytest.approx(s["delta_s_ev_per_k"] * compare.EV_J_MOL_K)
    assert s["delta_g_ev"] == pytest.approx(s["delta_h_ev"] - 298.15 * s["delta_s_ev_per_k"], abs=1e-9)
    by_dir = {r["dir"]: r for r in res.runs}
    assert by_dir["a"]["thermo_source"] == "analysis/thermo.csv" and by_dir["a"]["thermo"][0]["T_K"] == 298.15
    head = (runs / "compare_reactions.csv").read_text(encoding="utf-8").splitlines()
    assert "delta_h_kj_mol" in head[0] and "delta_s_j_mol_k" in head[0] and "thermo_note" in head[0]
    with open(runs / "compare_reactions.csv", newline="", encoding="utf-8") as f:
        rows = {r["name"]: r for r in csv.DictReader(f)}
    assert rows["scaled"]["thermo_T_K"] == "298.15" and float(rows["scaled"]["delta_g_ev"]) == pytest.approx(s["delta_g_ev"])
    js = json.loads((runs / "compare_summary.json").read_text(encoding="utf-8"))
    assert js["reactions"][0]["thermo"][0]["delta_s_j_mol_k"] == pytest.approx(0.0, abs=1e-6)
    assert "thermo_source" in (runs / "compare_runs.csv").read_text(encoding="utf-8").splitlines()[0]
    text = res.summary_text()
    assert "ΔH" in text and "ΔS" in text and "298.15 K" in text


def test_mismatched_temperature_pressure_model_or_missing_table(runs):
    rx = compare.parse_compare("mis=1:a,-1:two_t; pres=1:a,-1:other_p; mix=1:a,-1:harm; missing=1:a,-1:none; f=2:harm,-1:harm")
    res = compare.analyze_compare(runs, rx)
    mis, pres, mix, missing, f = res.reactions
    assert mis["thermo"] == [] and ("温度が揃っていません (298.15 / 298.15, 400)" in mis["thermo_note"]
                                    or "temperatures do not match (298.15 / 298.15, 400)" in mis["thermo_note"])
    assert pres["thermo"] == [] and ("圧力" in pres["thermo_note"] or "pressures" in pres["thermo_note"])
    assert mix["thermo"] == [] and ("混ざって" in mix["thermo_note"] or "mixed" in mix["thermo_note"])
    assert missing["thermo"] == [] and "thermo.csv" in missing["thermo_note"] and "none" in missing["thermo_note"]
    assert f["thermo"][0]["g_label"] == "F" and f["thermo"][0]["delta_h_ev"] is None
    row = _thermo_csv(runs / "harm")
    assert f["thermo"][0]["delta_g_ev"] == pytest.approx(float(row["f_total_ev"]))
    assert "F" in f["thermo_note"]
    with open(runs / "compare_reactions.csv", newline="", encoding="utf-8") as fh:
        rows = {r["name"]: r for r in csv.DictReader(fh)}
    assert rows["mis"]["delta_h_ev"] == "" and rows["mis"]["thermo_note"] == mis["thermo_note"]
    assert rows["f"]["delta_g_label"] == "F"


def test_two_temperatures_when_all_runs_share_them(runs):
    res = compare.analyze_compare(runs, compare.parse_compare("t2=1:two_t,-1:two_t"))
    th = res.reactions[0]["thermo"]
    assert [t["T_K"] for t in th] == [298.15, 400.0]
    assert all(t["delta_g_ev"] == pytest.approx(0.0, abs=1e-9) for t in th)
    with open(runs / "compare_reactions.csv", newline="", encoding="utf-8") as fh:
        row = next(csv.DictReader(fh))
    assert row["thermo_T_K"] == "298.15;400.0"


def test_gui_compare_table_has_the_thermo_columns(runs):
    from adit.gui import analysis_fields as AF
    res = compare.analyze_compare(runs, compare.parse_compare("same=1:a,-1:b; mis=1:a,-1:two_t"))
    sec = AF.compare_sections(res)[0]
    assert len(sec.columns) == 12 and any("ΔH" in c for c in sec.columns) and any("ΔS" in c for c in sec.columns)
    same, mis = sec.rows
    assert same[6] == "298.15" and same[7] == "+0.000" and same[9] == "+0.000"
    assert "298.15" in mis[6] and mis[7] == "-" and mis[10] in ("釣り合う", "balanced")
