
import csv
import json
import shutil
from pathlib import Path

import pytest

from tests.conftest import cfg_for, water_spec
from adit.scan import ScanError, apply_value, collect_scan, parse_scan, write_scan
from adit.spec import DftbMethod

REPO = Path(__file__).resolve().parent.parent


def test_parse_scan():
    s = parse_scan("kpoints.mesh = 4x4x4, 6x6x6")
    assert s.path == "kpoints.mesh" and s.values == ["4x4x4", "6x6x6"] and s.dir_name("4x4x4") == "mesh_4x4x4"
    for bad in ("method.ecutwfc", "method.ecutwfc=30", "method.ecutwfc=30,30", "1bad=1,2"):
        with pytest.raises(ScanError):
            parse_scan(bad)


def test_apply_value_changes_only_that_field():
    spec = water_spec(method=DftbMethod(sk_set="fake-1-0"))
    new = apply_value(spec, "method.scc_tolerance", "1e-7")
    assert new.method.scc_tolerance == 1e-7 and new.structure == spec.structure and new.task == spec.task
    assert "scan method.scc_tolerance=1e-7" in new.meta.comment
    with pytest.raises(ScanError):
        apply_value(spec, "method.no_such", "1")
    with pytest.raises(ScanError):
        apply_value(spec, "method.scc_tolerance", "abc")


def test_scale_is_for_periodic_only_and_scales_cell_and_positions():
    from ase.build import bulk
    from adit.spec import AtomsData, Structure
    spec = water_spec(method=DftbMethod(sk_set="fake-1-0"))
    with pytest.raises(ScanError):
        apply_value(spec, "scale", "1.01")
    si = bulk("Si", "diamond", a=5.43)
    per = spec.model_copy(update={"structure": Structure(source="bulk", source_ref="Si", atoms=AtomsData.from_ase(si))})
    new = apply_value(per, "scale", "1.02")
    a0, a1 = per.structure.atoms.to_ase(), new.structure.atoms.to_ase()
    assert abs(a1.get_volume() / a0.get_volume() - 1.02 ** 3) < 1e-9
    assert abs(a1.positions[1][0] / a0.positions[1][0] - 1.02) < 1e-9
    k = apply_value(per, "kpoints.mesh", "4x4x4")
    assert k.kpoints.mode == "mesh" and tuple(k.kpoints.mesh) == (4, 4, 4)


def test_write_scan_generates_each_value(sk_root, tmp_path):
    spec = water_spec(method=DftbMethod(sk_set="fake-1-0"))
    dirs = write_scan(spec, cfg_for(sk_root), tmp_path / "scan", parse_scan("method.max_scc_iterations=50,100,200"))
    assert [d.name for d in dirs] == ["max_scc_iterations_50", "max_scc_iterations_100", "max_scc_iterations_200"]
    assert "MaxSccIterations = 200" in (dirs[2] / "dftb_in.hsd").read_text(encoding="utf-8")
    meta = json.loads((tmp_path / "scan" / "scan.json").read_text(encoding="utf-8"))
    assert meta["path"] == "method.max_scc_iterations" and len(meta["dirs"]) == 3


def test_write_scan_writes_nothing_if_one_value_is_bad(sk_root, tmp_path):
    spec = water_spec(method=DftbMethod(sk_set="fake-1-0"))
    with pytest.raises(Exception):
        write_scan(spec, cfg_for(sk_root), tmp_path / "scan", parse_scan("method.max_scc_iterations=50,0"))
    assert not (tmp_path / "scan").exists() or not any((tmp_path / "scan").iterdir())


def test_fit_eos_recovers_known_curve():
    from ase.eos import birchmurnaghan
    from adit.scan import EV_ANG3_TO_GPA, fit_eos
    v = [40 * s ** 3 for s in (0.97, 0.98, 0.99, 1.0, 1.01, 1.02, 1.03)]
    e = [birchmurnaghan(x, -10.0, 0.6, 4.2, 40.0) for x in v]
    fit = fit_eos(v, e)
    assert abs(fit["v0_A3"] - 40.0) < 1e-4 and abs(fit["b0_gpa"] - 0.6 * EV_ANG3_TO_GPA) < 0.1 and abs(fit["bp"] - 4.2) < 0.05


def test_analyze_scan_columns_from_dftb_outputs(tmp_path):
    from adit.scan import analyze_scan
    out = tmp_path / "scan"
    for name in ("x_1", "x_2"):
        shutil.copytree(REPO / "examples" / "dftb_tio2_generated", out / name)
    (out / "scan.json").write_text(json.dumps({"path": "method.x", "values": ["1", "2"], "dirs": ["x_1", "x_2"]}), encoding="utf-8")
    res = analyze_scan(out)
    r = res.rows[0]
    assert r["natoms"] and r["diff_per_atom_mev"] == "0.000" and r["fmax_ev_ang"] and r["pressure_gpa"]
    assert abs(float(r["pressure_gpa"]) - 0.627536e7 / 1e9) < 1e-6
    assert "scan_energy" in res.figures and "method.x" in res.summary_text()


def test_espresso_force_and_pressure_are_read():
    from adit.scan import final_force_and_pressure
    fmax, pres = final_force_and_pressure(REPO / "examples" / "qe_si_generated", "espresso")
    assert fmax is not None and fmax < 1e-6 and abs(pres - 2.011) < 1e-9


def test_scale_scan_with_too_few_points_explains_instead_of_fitting(tmp_path):
    from adit.scan import analyze_scan
    out = tmp_path / "scan"
    for name in ("scale_0.99", "scale_1.01"):
        shutil.copytree(REPO / "examples" / "dftb_tio2_generated", out / name)
    (out / "scan.json").write_text(json.dumps({"path": "scale", "values": ["0.99", "1.01"], "dirs": ["scale_0.99", "scale_1.01"]}), encoding="utf-8")
    res = analyze_scan(out)
    assert res.eos is None and "5" in res.eos_note


def test_difference_reference_is_minimum_for_scale_and_last_otherwise():
    from adit.scan import _reference_row
    rows = [{"value": v, "energy_ev": e} for v, e in (("0.98", "-10.0"), ("1.00", "-10.5"), ("1.02", "-10.2"), ("1.04", ""))]
    base, ref = _reference_row(rows, "scale")
    assert ref == "min" and base["value"] == "1.00"
    base, ref = _reference_row(rows, "method.ecutwfc")
    assert ref == "last" and base["value"] == "1.02"
    assert _reference_row([{"value": "1", "energy_ev": ""}], "scale") == (None, "last")


def test_collect_scan_reads_final_energies(tmp_path):
    out = tmp_path / "scan"
    for name in ("x_1", "x_2"):
        shutil.copytree(REPO / "examples" / "water_generated", out / name)
    (out / "x_3").mkdir()
    (out / "scan.json").write_text(json.dumps({"path": "method.x", "values": ["1", "2", "3"], "dirs": ["x_1", "x_2", "x_3"]}), encoding="utf-8")
    rows = collect_scan(out)
    assert rows[0]["energy_ev"] and rows[0]["energy_ev"] == rows[1]["energy_ev"] and rows[2]["energy_ev"] == "" and rows[2]["note"]
    assert rows[0]["diff_from_last_mev"] == "0.000"
    assert list(csv.DictReader(open(out / "scan_energies.csv", encoding="utf-8")))[2]["dir"] == "x_3"


def test_write_scan_checks_the_output_directory_before_creating_anything(sk_root, tmp_path):
    from adit.project import OutputNotEmpty, ProjectError
    spec = water_spec(method=DftbMethod(sk_set="fake-1-0"))
    scan = parse_scan("method.max_scc_iterations=50,100")
    with pytest.raises(ProjectError) as ex:
        write_scan(spec, cfg_for(sk_root), tmp_path / "no" / "such" / "scan", scan)
    assert [e.location for e in ex.value.errors] == ["output_dir"] and not (tmp_path / "no").exists()
    out = tmp_path / "scan"
    write_scan(spec, cfg_for(sk_root), out, scan)
    before = (out / "scan.json").read_text(encoding="utf-8")
    with pytest.raises(OutputNotEmpty):
        write_scan(spec, cfg_for(sk_root), out, parse_scan("method.max_scc_iterations=70,80"))
    assert (out / "scan.json").read_text(encoding="utf-8") == before
    write_scan(spec, cfg_for(sk_root), out, parse_scan("method.max_scc_iterations=70,80"), overwrite=True)
    assert json.loads((out / "scan.json").read_text(encoding="utf-8"))["values"] == ["70", "80"]


def test_batch_check_output_reports_a_missing_parent_and_a_file(tmp_path):
    from adit.batch import check_output
    from adit.project import ProjectError
    with pytest.raises(ProjectError) as ex:
        check_output(tmp_path / "no" / "such" / "dir", False)
    assert [e.location for e in ex.value.errors] == ["output_dir"]
    (tmp_path / "file").write_text("x", encoding="utf-8")
    with pytest.raises(ProjectError):
        check_output(tmp_path / "file", True)
    check_output(tmp_path / "new", False)
    assert not (tmp_path / "new").exists()
