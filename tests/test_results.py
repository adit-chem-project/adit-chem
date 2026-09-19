from pathlib import Path

from adit.results import read_results_tag, summarize_run

REPO = Path(__file__).resolve().parent.parent


def test_read_results_tag():
    d = read_results_tag(REPO / "examples" / "water_generated" / "results.tag")
    assert abs(d["mermin_energy"][0] - (-4.0779379326)) < 1e-9
    assert d["number_of_electrons"][0] == 8.0
    assert d["forces"].size == 9


def test_summarize_water():
    s = summarize_run(REPO / "examples" / "water_generated")
    assert s.converged and s.geometry_steps == 10 and s.scc_iterations_last == 3
    assert abs(s.mermin_energy_hartree - (-4.0779379326)) < 1e-9
    names = dict(s.bonds)
    assert set(names) == {"O1-H2", "O1-H3"} and abs(names["O1-H2"] - 0.9672) < 1e-3
    assert "0.9672" in s.text()


def test_not_run_yet_without_provenance_needs_spec_json(tmp_path):
    from adit.results import not_run_yet
    d = tmp_path / "run"
    d.mkdir()
    (d / "OUTCAR").write_text("foreign output\n", encoding="utf-8")
    assert not_run_yet(d) is False
    (d / "spec.json").write_text("{}", encoding="utf-8")
    (d / "OUTCAR").unlink()
    assert not_run_yet(d) is True
    (d / "output.log").write_text("Total Energy: -4.0 H\n", encoding="utf-8")
    assert not_run_yet(d) is False
    assert not_run_yet(tmp_path / "missing") is False
