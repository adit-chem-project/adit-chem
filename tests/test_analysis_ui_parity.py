import dataclasses
from pathlib import Path

import pytest

from adit.analysis import AnalysisOptions
from adit.gui import analysis_fields as AF

INTERNAL = {
    "energy", "temperature", "bonds", "vibrations", "bands",
    "out_dir", "collect", "vanhove_here",
    "msd_fit_fs", "skip_frames", "stride", "rdf_rmax", "msd_species", "dos_sigma",
    "zdens_bin_ang", "memory_budget_mb", "symprecs", "thermo", "uvvis_broadening",
    "rdf", "msd", "dos", "zdens", "stats", "pdos", "export", "export_unwrap",
    "msd_axes", "msd_remove_drift", "msd_error_blocks", "msd_per_atom",
    "vanhove", "vanhove_displacement",
    "effective_mass", "projected_bands", "optical",
}
TEMPLATE = Path(__file__).resolve().parents[1] / "src" / "adit" / "web" / "templates" / "analysis.html"


def _option_names() -> set[str]:
    return {f.name for f in dataclasses.fields(AnalysisOptions)} - INTERNAL


def test_every_option_is_reachable_from_the_fields():
    source = Path(AF.__file__).read_text(encoding="utf-8")
    unreachable = [name for name in _option_names() if f'"{name}"' not in source]
    assert not unreachable, f"画面から指定できない解析の条件: {unreachable}"


def test_desktop_panel_has_a_widget_for_every_field():
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from adit.gui.panels.analysis_panel import AnalysisPanel

    QApplication.instance() or QApplication([])
    panel = AnalysisPanel()
    fields = panel.fields()
    skip = {"msd_fit", "rmax_label", "sigma_label", "skip_label"}
    missing = [k for k in AF.LABELS
               if k not in fields and k not in skip and not k.startswith(("rx_", "compare", "open_", "more", "export"))]
    assert not missing, f"デスクトップ版に欄が無い: {missing}"


def test_web_form_has_every_field():
    html = TEMPLATE.read_text(encoding="utf-8")
    skip = {"more", "compare", "compare_base", "compare_run", "compare_rows", "rx_name", "rx_nu", "rx_dir",
            "open_export", "export", "msd_fit",
            "rmax_label", "sigma_label", "skip_label"}
    missing = [k for k in AF.LABELS if k not in skip and f'"{k}"' not in html]
    assert not missing, f"ウェブ版に欄が無い: {missing}"


def test_fields_round_trip_through_options():
    opts = AnalysisOptions(select="element O", coordination_cutoff=3.2, centrosymmetry_neighbors=12,
                           steinhardt_cutoff=3.0, cluster_cutoff=2.8, adf="O,H", adf_cutoff=1.3,
                           structure_factor=True, hbond="3.5,150", radius_of_gyration=True,
                           density_grid="32,32,32", voronoi=True, voronoi_face_threshold=0.01,
                           sasa="bondi,1.4", distances=["1,2"], angles=["2,1,3"], dihedrals=["1,2,3,4"],
                           rmsd_reference=2, rmsf=True, vacf=True, conductivity_charge=1.0,
                           conductivity_temperature_k=500.0, displacement_reference=1, strain_cutoff=3.0,
                           pca=3, cluster=2, fes_temperature_k=300.0, fes_bins=40, fes_unit="eV",
                           bands_window_ev=6.0, effective_mass_points=7, bader="ACF.dat", bader_valence="4,6",
                           xrd="CuKa", xrd_range=(10.0, 80.0), viscosity=True, plane_average="c",
                           cube_unit="ev", work_function=True, heavy_limit_seconds=120.0,
                           rdf_pairs=[("O", "H")], zdens_axis="a")
    back = AF.options_from_fields(AF.fields_from_options(opts))
    for name in _option_names():
        assert getattr(back, name) == getattr(opts, name), name


def test_new_labels_have_english_and_help():
    from adit.gui.help import help_for
    from adit.gui.i18n import _EN

    for key in ("select", "coordination", "voronoi", "sasa", "fes", "bader", "xrd", "viscosity", "pca"):
        ja = AF.LABELS[key][0]
        assert ja in _EN and _EN[ja]
        h = help_for(ja)
        assert h is not None and h.ja and h.en and not h.required


ONLY_DESKTOP = {
    "ASE GUI を外部プロセスで開く": "サーバー側の画面が開いてしまうので、ウェブ版では出さない",
    "書き出したフォルダを開く": "同上 (サーバー側のフォルダが開く)",
    "構造を作っている途中の中止": "ウェブ版はプレビューの中で作り終えるまで待つ (途中で止める画面を持たない)",
    "ファイルを選ぶダイアログ": "ブラウザからはサーバー側のパスを選べないので、ウェブ版は文字で書く",
}
ONLY_CLI = {
    "画面を持たない一括処理": "コマンドは複数のディレクトリをまとめて回せる (画面は 1 つずつ)",
}
SAME_BUT_DIFFERENT_CONTROL = {
    "SHA-256 の表示": "デスクトップは押して開くボタン、ウェブは <details>。出る中身は同じ",
    "GROMACS の構造を構造の欄にも使う": "デスクトップはボタン、ウェブはチェックボックス",
    "溶液の成分": "デスクトップは行を足す表、ウェブは 1 行 1 成分のテキスト欄",
}


def test_the_report_feature_is_in_all_three_entry_points():
    from adit.gui import report_fields as RF

    assert {"rep_dirs", "rep_lang", "rep_methods", "rep_conditions", "rep_results", "rep_bundle",
            "rep_check"} <= set(RF.LABELS)
    root = Path(__file__).resolve().parents[1] / "src" / "adit"
    assert (root / "gui" / "report_dialog.py").is_file()
    assert (root / "web" / "templates" / "report.html").is_file()
    server = (root / "web" / "server.py").read_text(encoding="utf-8")
    assert '"/report"' in server
    base = (root / "web" / "templates" / "base.html").read_text(encoding="utf-8")
    assert 'href="/report"' in base
    panel = (root / "gui" / "panels" / "analysis_panel.py").read_text(encoding="utf-8")
    assert "btn_report" in panel


def test_the_report_labels_have_english_and_help():
    from adit.gui import report_fields as RF
    from adit.gui.help import help_for
    from adit.gui.i18n import _EN

    for ja, en in RF.LABELS.values():
        assert ja in _EN and _EN[ja] == en, ja
        assert help_for(ja) is not None, ja


def test_both_screens_build_the_report_through_the_same_code():
    root = Path(__file__).resolve().parents[1] / "src" / "adit"
    dialog = (root / "gui" / "report_dialog.py").read_text(encoding="utf-8")
    server = (root / "web" / "server.py").read_text(encoding="utf-8")
    assert "R.build(R.request_from_fields" in dialog
    assert "RF.build(RF.request_from_fields" in server


def test_the_same_report_comes_out_of_the_two_screens(tmp_path):
    from adit.gui import report_fields as RF

    run = Path(__file__).resolve().parents[1] / "examples" / "dftb_md_water_generated"
    fields = {"rep_dirs": str(run), "rep_lang": "ja"}
    first = RF.build(RF.request_from_fields(fields))
    second = RF.build(RF.request_from_fields({**fields, "rep_results": str(tmp_path / "r.csv")}))
    assert first.methods_text == second.methods_text
    assert (tmp_path / "r.csv").is_file()


def test_draw_is_in_both_screens_now():
    root = Path(__file__).resolve().parents[1] / "src" / "adit"
    assert (root / "sketch.py").is_file()
    assert (root / "web" / "static" / "sketch.js").is_file()
    assert (root / "web" / "templates" / "draw.html").is_file()
    assert '"/draw"' in (root / "web" / "server.py").read_text(encoding="utf-8")
    assert 'href="/draw"' in (root / "web" / "templates" / "base.html").read_text(encoding="utf-8")
    assert "Draw" not in " ".join(ONLY_DESKTOP)
    lines = (root / "sketch.py").read_text(encoding="utf-8").splitlines()
    assert not [x for x in lines if x.startswith(("import ", "from ")) and "PySide6" in x]


def test_the_macos_bundle_is_prepared():
    repo = Path(__file__).resolve().parents[1]
    spec = (repo / "packaging" / "adit.spec").read_text(encoding="utf-8")
    assert "BUNDLE(" in spec and "darwin" in spec
    workflow = (repo / ".github" / "workflows" / "macos-app.yml").read_text(encoding="utf-8")
    assert "macos-14" in workflow and "macos-13" not in workflow  # no Intel build: the runner is gone
    doc = (repo / "docs" / "INSTALL.md").read_text(encoding="utf-8")
    assert "署名" in doc and "公証" in doc
    assert "実機" in doc          # the document says macOS was not checked on a real machine
