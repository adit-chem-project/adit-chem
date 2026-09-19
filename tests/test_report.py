
import json
import zipfile
from pathlib import Path

import pytest

from adit import lang
from adit.project import write_project
from adit.report import (ReportError, bundle_files, check_lines, conditions_csv, condition_rows,
                          load_run_report, main, methods_markdown, tidy_version, write_bundle)
from adit.spec import KPoints, MDSettings, Task
from tests.conftest import cfg_for, water_spec


@pytest.fixture
def run_dir(tmp_path, sk_root) -> Path:
    out = tmp_path / "water_opt"
    write_project(water_spec(), cfg_for(sk_root), out)
    return out


@pytest.fixture(autouse=True)
def _japanese():
    before = lang.LANGUAGE
    lang.set_language("ja")
    yield
    lang.set_language(before)


def test_reads_spec_and_says_what_is_missing(run_dir):
    report = load_run_report(run_dir)
    assert report.code == "dftbplus"
    assert report.code_version == ""
    assert any("バージョンは未記録" in n for n in report.notes)
    assert report.provenance["adit_version"]


def test_missing_directory_and_spec_are_errors(tmp_path):
    with pytest.raises(ReportError):
        load_run_report(tmp_path / "nope")
    empty = tmp_path / "empty"
    empty.mkdir()
    with pytest.raises(ReportError):
        load_run_report(empty)


def test_tidy_version_strips_only_our_marker():
    assert tidy_version("adit-psi4: psi4 1.11") == "psi4 1.11"
    assert tidy_version(".Version 10.0.3 of ABINIT, released Apr 2024.") == "Version 10.0.3 of ABINIT, released Apr 2024"
    assert tidy_version("DFTB+ release 25.1") == "DFTB+ release 25.1"


def test_methods_text_has_structure_code_and_settings(run_dir):
    text = methods_markdown([load_run_report(run_dir)], "ja")
    assert "H2O (3 原子" in text and "非周期" in text
    assert "DFTB+ (バージョンは未記録)" in text
    assert "構造最適化を行いました" in text
    assert "| method.sk_set | fake-1-0 |" in text
    assert "task.force_tolerance_ev_per_ang" in text and "eV/Å" in text
    assert "spec.json` の記録であって" in text


def test_methods_text_after_a_run_uses_the_recorded_version(run_dir):
    (run_dir / "code_version.txt").write_text("DFTB+ release 25.1\n", encoding="utf-8")
    text = methods_markdown([load_run_report(run_dir)], "ja")
    assert "DFTB+ (DFTB+ release 25.1)" in text and "バージョンは計算の出力から写したもの" in text
    assert "バージョンは未記録" not in text


def test_both_languages(run_dir):
    text = methods_markdown([load_run_report(run_dir)], "both")
    assert "# 計算条件 (方法)" in text and "# Computational details (methods)" in text
    assert lang.LANGUAGE == "ja"


def test_unknown_language_is_rejected(run_dir):
    with pytest.raises(ReportError):
        methods_markdown([load_run_report(run_dir)], "fr")


def test_md_and_kpoints_rows(tmp_path, sk_root):
    from adit.spec import AtomsData

    spec = water_spec(task=Task(type="molecular_dynamics",
                                md=MDSettings(ensemble="NVT", thermostat="berendsen", temperature_k=300.0,
                                              timestep_fs=1.0, steps=100, dump_interval=10)),
                      kpoints=KPoints(mode="mesh", mesh=(2, 2, 2)))
    atoms = spec.atoms
    atoms.set_cell([8.0, 8.0, 8.0])
    atoms.pbc = True
    spec = spec.model_copy(update={"structure": spec.structure.model_copy(
        update={"atoms": AtomsData.from_ase(atoms)})})
    out = tmp_path / "md"
    write_project(spec, cfg_for(sk_root), out)
    rows = dict((name, value) for name, value, _ in condition_rows(load_run_report(out)))
    assert rows["task.md.temperature_k"] == "300" and rows["task.md.ensemble"] == "NVT"
    assert rows["kpoints.mesh"] == "2×2×2"
    assert "task.optimizer" not in rows


def test_zero_method_fields_are_not_shown_as_values(run_dir):
    rows = [name for name, _, _ in condition_rows(load_run_report(run_dir))]
    assert all(not name.endswith("_ev") or True for name in rows)
    values = {name: value for name, value, _ in condition_rows(load_run_report(run_dir))}
    assert all(value != "0" for name, value in values.items() if name.startswith("method."))


def test_conditions_csv_has_one_row_per_run(run_dir, tmp_path, sk_root):
    other = tmp_path / "water_sp"
    write_project(water_spec(task=Task(type="single_point")), cfg_for(sk_root), other)
    text = conditions_csv([load_run_report(run_dir), load_run_report(other)])
    lines = text.strip().splitlines()
    assert len(lines) == 3 and lines[0].startswith("run_dir,code,code_version,formula")
    assert "geometry_optimization" in lines[1] and "single_point" in lines[2]
    assert ",dftbplus,," in lines[1]


def test_bundle_zip_has_inputs_and_manifest(run_dir, tmp_path):
    dest = tmp_path / "pack.zip"
    write_bundle([load_run_report(run_dir)], dest)
    with zipfile.ZipFile(dest) as zf:
        names = zf.namelist()
        manifest = json.loads(zf.read("manifest.json"))
    assert "methods.md" in names and f"{run_dir.name}/spec.json" in names
    assert f"{run_dir.name}/dftb_in.hsd" in names and f"{run_dir.name}/submit.sh" in names
    assert not any(n.endswith("detailed.out") for n in names)
    files = {f["name"]: f for f in manifest["runs"][0]["files"]}
    assert files["dftb_in.hsd"]["state"] == "same" and len(files["dftb_in.hsd"]["sha256"]) == 64


def test_bundle_directory_and_no_overwrite(run_dir, tmp_path):
    dest = tmp_path / "pack"
    write_bundle([load_run_report(run_dir)], dest)
    assert (dest / "methods.md").is_file() and (dest / run_dir.name / "dftb_in.hsd").is_file()
    with pytest.raises(ReportError):
        write_bundle([load_run_report(run_dir)], dest)


def test_bundle_list_skips_missing_files(run_dir):
    names = bundle_files(load_run_report(run_dir))
    assert "spec.json" in names and "README.txt" in names
    assert "code_version.txt" not in names


def test_check_detects_changed_and_missing_inputs(run_dir):
    lines, verdict = check_lines(load_run_report(run_dir))
    assert verdict == "ok" and any("同じ" in l for l in lines)
    (run_dir / "dftb_in.hsd").write_text("# 人が書き換えた\n", encoding="utf-8")
    lines, verdict = check_lines(load_run_report(run_dir))
    assert verdict == "bad" and any("書き換え" in l and "dftb_in.hsd" in l for l in lines)
    assert any("書き換えられたファイル 1 件" in l for l in lines)
    (run_dir / "dftb_in.hsd").unlink()
    lines, verdict = check_lines(load_run_report(run_dir))
    assert verdict == "bad" and any("ありません" in l for l in lines)
    assert any("無くなったファイル 1 件" in l for l in lines)


def test_check_covers_copied_parameter_files(run_dir):
    copied = sorted((run_dir / "skf").glob("*.skf"))
    assert copied, "写したファイルが無いと試験にならない"
    copied[0].write_text("tampered\n", encoding="utf-8")
    lines, verdict = check_lines(load_run_report(run_dir))
    assert verdict == "bad" and any(copied[0].name in l and "書き換え" in l for l in lines)
    assert any("写したファイル" in l for l in lines)


def test_check_says_it_cannot_check_without_records(run_dir):
    data = json.loads((run_dir / "spec.json").read_text(encoding="utf-8"))
    for key in ("inputs", "files", "generated_files"):
        data["provenance"].pop(key, None)
    (run_dir / "spec.json").write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    lines, verdict = check_lines(load_run_report(run_dir))
    assert verdict == "unknown" and any("確かめられません" in l for l in lines)
    assert not any("2026-09-13 より前に生成したディレクトリです" in l for l in lines)
    assert main([str(run_dir), "--check"]) == 2


# ---- CLI ----
def test_cli_writes_files_and_reports_tampering(run_dir, tmp_path, capsys):
    out = tmp_path / "methods.md"
    csv_path = tmp_path / "conditions.csv"
    assert main([str(run_dir), "-o", str(out), "--csv", str(csv_path), "--lang", "both"]) == 0
    assert "# 計算条件 (方法)" in out.read_text(encoding="utf-8")
    assert csv_path.read_text(encoding="utf-8").startswith("run_dir,")
    assert main([str(run_dir), "--check"]) == 0
    (run_dir / "submit.sh").write_text("echo tampered\n", encoding="utf-8")
    assert main([str(run_dir), "--check"]) == 1
    assert "書き換え" in capsys.readouterr().out


def test_cli_reports_a_bad_directory(tmp_path, capsys):
    assert main([str(tmp_path / "missing")]) == 2
    assert "ディレクトリがありません" in capsys.readouterr().out


def test_analysis_summary_is_copied_when_it_exists(run_dir):
    (run_dir / "analysis").mkdir()
    (run_dir / "analysis" / "summary.txt").write_text("エネルギー: 1 点、最終値 -100.0 eV\n", encoding="utf-8")
    text = methods_markdown([load_run_report(run_dir)], "ja")
    assert "### 結果 (解析の要約)" in text and "最終値 -100.0 eV" in text
    assert "adit-analyze" in text


def test_units_come_from_the_code_table(tmp_path, sk_root):
    from adit.spec import EspressoMethod, KPoints
    from ase.build import bulk

    from adit.spec import AtomsData

    spec = water_spec(method=EspressoMethod(pseudo_set="fake", pseudo={"Si": "Si.upf"}, ecutwfc=40.0, conv_thr=1e-8),
                      kpoints=KPoints(mode="mesh", mesh=(2, 2, 2)), task=Task(type="single_point"))
    spec = spec.model_copy(update={"structure": spec.structure.model_copy(
        update={"source": "bulk", "source_ref": "Si", "atoms": AtomsData.from_ase(bulk("Si", "diamond", a=5.43))})})
    from adit.report import RunReport, condition_rows

    rows = {name: unit for name, _, unit in condition_rows(RunReport(run_dir=tmp_path, spec=spec))}
    assert rows["method.ecutwfc"] == "Ry" and rows["method.conv_thr"] == "Ry"


def test_fields_not_written_into_the_input_go_to_their_own_table(run_dir):
    text = methods_markdown([load_run_report(run_dir)], "ja")
    assert "### この計算コードの入力に書かれていない共通の欄" in text
    body = text.split("### この計算コードの入力に書かれていない共通の欄")[1]
    assert "task.relax_cell" in body and "非周期系にセルはありません" in body
    head = text.split("### この計算コードの入力に書かれていない共通の欄")[0]
    assert "| task.relax_cell |" not in head


def test_runtime_section_is_written(run_dir):
    text = methods_markdown([load_run_report(run_dir)], "ja")
    assert "### 実行のしかた (spec.json の記録)" in text
    assert "runtime.mpiprocs" in text and "runtime.omp_threads" in text
    assert "計算機の機種・OS・コンパイラ" in text


def test_csv_uses_machine_readable_booleans(run_dir, tmp_path, sk_root):
    from adit.spec import MlipMethod

    spec = water_spec(method=MlipMethod(model_family="mace_mp", dispersion=True), task=Task(type="single_point"))
    out = tmp_path / "mlip"
    write_project(spec, cfg_for(sk_root), out)
    text = conditions_csv([load_run_report(out)])
    assert "true" in text and "はい" not in text


def test_warns_when_the_inputs_cannot_be_identified(run_dir):
    data = json.loads((run_dir / "spec.json").read_text(encoding="utf-8"))
    data["provenance"].pop("inputs")
    data["method"]["sk_set"] = "no-such-set"
    (run_dir / "spec.json").write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    text = methods_markdown([load_run_report(run_dir)], "ja")
    assert "この計算の入力ファイルを特定できません" in text



def test_analysis_stops_when_the_calculation_has_not_run(run_dir, capsys):
    from adit.analysis.cli import main as analyze_main
    from adit.results import not_run_yet

    assert not_run_yet(run_dir) is True
    assert analyze_main([str(run_dir)]) == 2
    assert "まだ実行していないようです" in capsys.readouterr().err
    (run_dir / "output.log").write_text("Total Energy: -4.0 H\n", encoding="utf-8")
    assert not_run_yet(run_dir) is False


def test_methods_section_names_the_references_to_cite(run_dir):
    text = methods_markdown([load_run_report(run_dir)], "ja")
    assert "### 引用 (文献)" in text
    assert "引用: dftbplus_hourahine2025, dftbplus_hourahine2020, adit_" in text
    assert "未記録: Slater-Koster セット fake-1-0 の文献" in text
    english = methods_markdown([load_run_report(run_dir)], "en")
    assert "### References to cite" in english and "Cite: dftbplus_hourahine2025" in english


def test_bib_option_and_bundle_write_references(run_dir, tmp_path, capsys):
    bib = tmp_path / "references.bib"
    assert main([str(run_dir), "--bib", str(bib)]) == 0
    text = bib.read_text(encoding="utf-8")
    assert text.count("@article{dftbplus_hourahine2020,") == 1 and "@software{adit_" in text
    assert "% 出典: https://dftbplus.org/about/index.html" in text
    assert "文献 (BibTeX) を書きました" in capsys.readouterr().out
    dest = tmp_path / "pack.zip"
    write_bundle([load_run_report(run_dir)], dest)
    with zipfile.ZipFile(dest) as zf:
        assert "references.bib" in zf.namelist() and b"@article{dftbplus_hourahine2020," in zf.read("references.bib")
