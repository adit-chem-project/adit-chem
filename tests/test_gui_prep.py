
import os
import shutil
from pathlib import Path

import pytest

pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PySide6.QtWidgets import QApplication, QLabel, QMessageBox  # noqa: E402

from tests.conftest import pbs_profile  # noqa: E402
from adit.codes.orca import solvent_names as orca_names  # noqa: E402
from adit.codes.xtb import solvent_names as xtb_names  # noqa: E402
from adit.config import default_config  # noqa: E402
from adit.gui.help import help_for  # noqa: E402
from adit.gui.i18n import _EN  # noqa: E402
from adit.gui.convert_dialog import ConvertDialog  # noqa: E402
from adit.gui.main_window import MainWindow  # noqa: E402
from adit.gui.panels.cp2k_panel import Cp2kMethodPanel  # noqa: E402
from adit.gui.panels.espresso_panel import EspressoMethodPanel  # noqa: E402
from adit.gui.panels.method_panel import CODES, DftbMethodPanel  # noqa: E402
from adit.gui.panels.orca_panel import OrcaMethodPanel  # noqa: E402
from adit.gui.panels.vasp_panel import VaspMethodPanel  # noqa: E402
from adit.gui.panels.xtb_panel import XtbMethodPanel  # noqa: E402
from adit.gui.prep_widgets import HubbardTable  # noqa: E402
from adit.spec import Cp2kMethod, DftbMethod, EspressoMethod, HubbardU, OrcaMethod, VaspMethod, XtbMethod  # noqa: E402

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"
NEW_LABELS = ["溶媒モデル (--alpb / --gbsa)", "溶媒モデル (CPCM / SMD)", "溶媒", "溶媒のパラメータファイル (GBSA)", "溶媒の比誘電率 (SCCS)",
              "元素ごとの初期磁気モーメント [μB]", "DFT+U", "LDAUTYPE", "starting_magnetization (元素ごと)", "HUBBARD の射影",
              "元素ごとの MAGNETIZATION", "PLUS_U_METHOD"]


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


@pytest.fixture
def boxes(monkeypatch):
    seen = []
    for name in ("information", "critical", "warning"):
        monkeypatch.setattr(QMessageBox, name, staticmethod(lambda *a, _n=name, **k: seen.append((_n, a[2] if len(a) > 2 else ""))))
    monkeypatch.setattr(QMessageBox, "question", staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))
    return seen


@pytest.fixture
def conf(tmp_path, monkeypatch):
    monkeypatch.setenv("ADIT_CONFIG", str(tmp_path / "conf" / "cluster.toml"))
    return tmp_path


def make_window(sk_root, tmp_path) -> MainWindow:
    cfg = default_config(sk_root=str(sk_root)); cfg.profiles["cluster"] = pbs_profile()
    cfg.templates_dir = str(tmp_path / "lab")
    win = MainWindow(cfg, tmp_path / "cluster.toml")
    win.method.sk_set.setCurrentText("fake-1-0")
    win.runtime.outdir.setText(str(tmp_path / "out"))
    win.refresh_preview()
    return win


def use_code(win: MainWindow, code: str) -> None:
    win.method.code.setCurrentIndex(list(CODES).index(code))


def test_convert_dialog_uses_shared_structure_converter(app, tmp_path):
    cfg = default_config(); cfg.templates_dir = str(tmp_path / "templates")
    dialog = ConvertDialog(lambda: None, cfg)
    source = tmp_path / "water.xyz"
    source.write_text("3\nH2O\nO 0 0 0\nH 0.8 0.6 0\nH -0.8 0.6 0\n", encoding="utf-8")
    output = tmp_path / "water.cif"
    dialog.struct_source.setText(str(source)); dialog.struct_output.setText(str(output)); dialog.cell.setText("15")
    dialog.convert()
    assert output.is_file() and "構造を書きました" in dialog.status.text()
    dialog.kind.setCurrentIndex(1)
    assert not dialog.run.isEnabled()


def test_md_fields_do_not_silently_clamp_invalid_input(app):
    from adit.gui.panels.task_panel import TaskPanel
    panel = TaskPanel()
    panel.temperature.setValue(-50); panel.timestep.setValue(999); panel.md_steps.setValue(0)
    md = panel.task().md
    assert md.temperature_k == -50 and md.timestep_fs == 999 and md.steps == 0


def test_vasp_advanced_fields_roundtrip(app, tmp_path):
    panel = VaspMethodPanel(default_config())
    method = VaspMethod(nelmin=5, lasph=True, lmaxmix=4, nbands=96, isym=0, idipol=3, ldipol=True,
                        dipol="0.5 0.5 0.4", kpoints_centering="gamma")
    panel.set_method(method)
    assert panel.method() == method


def test_xtb_solvent_roundtrip_and_candidates(app):
    p = XtbMethodPanel()
    p.solvation.setCurrentIndex(p.solvation.findData("gbsa")); p.gfn.setCurrentIndex(p.gfn.findData("1"))
    names = [p.solvent.itemData(i) for i in range(1, p.solvent.count())]
    assert names == xtb_names("gbsa", "1") and p.solvent.currentData() == ""
    p.gfn.setCurrentIndex(p.gfn.findData("2"))
    assert [p.solvent.itemData(i) for i in range(1, p.solvent.count())] == xtb_names("gbsa", "2")
    m = XtbMethod(gfn="2", solvation="alpb", solvent="water")
    p.set_method(m)
    assert p.method() == m
    p.set_method(XtbMethod())
    assert p.method().solvation == "none" and p.method().solvent == ""


def test_xtb_gfn0_has_no_candidates_but_keeps_value(app):
    p = XtbMethodPanel()
    p.set_method(XtbMethod(gfn="0", solvation="alpb", solvent="water"))
    assert xtb_names("alpb", "0") == [] and p.method().solvent == "water"


def test_orca_solvent_candidates_follow_model(app):
    p = OrcaMethodPanel()
    for model in ("cpcm", "smd"):
        p.solvation.setCurrentIndex(p.solvation.findData(model))
        assert [p.solvent.itemData(i) for i in range(1, p.solvent.count())] == orca_names(model)
    m = OrcaMethod(solvation="smd", solvent="water")
    p.set_method(m)
    assert p.method() == m


def test_dftb_solvation_file_and_periodic(app, tmp_path):
    p = DftbMethodPanel("")
    f = tmp_path / "param_gbsa_h2o.txt"; f.write_text("x\n", encoding="utf-8")
    p.solv_file.setText(str(f))
    assert p.method().solvation_param_file == str(f)
    p.set_method(DftbMethod(sk_set="", solvation_param_file=""))
    p.set_periodic(True)
    assert not p.solv_file.isEnabled()
    p.set_periodic(False)
    assert p.solv_file.isEnabled()


def test_hubbard_table_elements_and_no_default_u(app):
    t = HubbardTable()
    t.set_elements(["Fe", "O"])
    t.add_row()
    ec = t.rows[0][0]
    assert [ec.itemData(i) for i in range(ec.count())] == ["", "Fe", "O"]
    assert t.rows[0][2].text() == ""
    ec.setCurrentIndex(ec.findData("Fe")); t.rows[0][1].setCurrentText("3d")
    with pytest.raises(ValueError, match="U"):
        t.hubbard()
    t.rows[0][2].setText("5.3")
    assert t.hubbard() == {"Fe": HubbardU(orbital="3d", u_ev=5.3, j_ev=0.0)}
    t.set_hubbard({"Ni": HubbardU(orbital="3d", u_ev=6.0, j_ev=1.0)})
    assert "構造に無い" in t.rows[0][0].currentText() or "not in the structure" in t.rows[0][0].currentText()
    t.remove_row(0)
    assert t.hubbard() == {}


@pytest.mark.parametrize("panel_cls, method", [
    (VaspMethodPanel, VaspMethod(ispin=2, magmom_by_element={"Fe": 4.0}, hubbard={"Fe": HubbardU(orbital="3d", u_ev=5.3, j_ev=1.0)}, ldau_type=1)),
    (EspressoMethodPanel, EspressoMethod(nspin=2, starting_magnetization={"Fe": 0.5}, hubbard={"Fe": HubbardU(orbital="3d", u_ev=4.3)},
                                         hubbard_projector="ortho-atomic")),
    (Cp2kMethodPanel, Cp2kMethod(uks=True, magnetization_by_element={"Fe": 2.0}, hubbard={"Fe": HubbardU(orbital="3d", u_ev=4.0, j_ev=0.5)},
                                 plus_u_method="LOWDIN", sccs_relative_permittivity=78.4)),
])
def test_magnetism_and_dftu_roundtrip(app, panel_cls, method):
    cfg = default_config()
    p = panel_cls(cfg)
    p.set_context(["Fe", "O"], "local", cfg)
    p.set_method(method)
    got = p.method()
    for k in ("magmom_by_element", "starting_magnetization", "magnetization_by_element", "hubbard", "ldau_type", "hubbard_projector",
              "plus_u_method", "sccs_relative_permittivity"):
        if hasattr(method, k):
            assert getattr(got, k) == getattr(method, k), k


def test_qe_table_has_no_j_column(app):
    p = EspressoMethodPanel(default_config())
    p.set_context(["Fe"], "local")
    p.hubbard.add_row("Fe", "3d", "4", "")
    assert not p.hubbard.with_j and p.method().hubbard["Fe"].j_ev == 0.0


def test_cp2k_per_element_rows_are_full_width(app, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path); win.resize(1800, 1000); win.show()
    win.set_mode(win.MODE_SETTINGS); app.processEvents()
    use_code(win, "cp2k"); app.processEvents()
    lab = next(w for w in win.method.cp2k.findChildren(QLabel) if w.property("adit_key") == "元素ごとの基底と擬ポテンシャル")
    assert lab.width() >= lab.sizeHint().width()
    combos = list(win.method.cp2k.basis_widgets.values())
    assert combos and min(c.width() for c in combos) >= 150
    win.close()


@pytest.fixture
def prev_md(tmp_path) -> Path:
    d = tmp_path / "prev_md"
    shutil.copytree(EX / "dftb_md_water_generated", d)
    return d


def test_continue_dialog_carries_velocities(app, boxes, sk_root, tmp_path, prev_md):
    win = make_window(sk_root, tmp_path)
    dlg = win.continue_dialog()
    assert dlg.table.minimumHeight() >= dlg.table.fontMetrics().lineSpacing() * 10
    dlg.folder.setText(str(prev_md))
    assert "dftbplus" in dlg.found.text()
    dlg.make()
    assert dlg.spec is not None and "geo_end.xyz" in dlg.summary()
    win.apply_continuation(dlg.spec)
    win.method.sk_set.setCurrentText("fake-1-0"); win.refresh_preview()
    s = win.current_spec()
    assert s.structure.velocities is not None and s.handoff.velocities and s.meta.continued_from["dir"] == str(prev_md)
    assert "Velocities [AA/ps]" in win.build().texts["dftb_in.hsd"]
    assert str(prev_md) in win.preview.origin.text()
    win.task.type.setCurrentIndex(win.task.type.findData("geometry_optimization"))
    s2 = win.current_spec()
    assert s2.structure.velocities is None and not s2.handoff.velocities
    win.clear_origin()
    assert win.current_spec().handoff is None and not win.preview.origin_row.isVisibleTo(win.preview)


def test_continue_without_velocities_and_errors(app, boxes, sk_root, tmp_path, prev_md):
    win = make_window(sk_root, tmp_path)
    dlg = win.continue_dialog(); dlg.folder.setText(str(prev_md)); dlg.velocities.setChecked(False); dlg.make()
    assert dlg.spec.structure.velocities is None
    bad = win.continue_dialog(); bad.folder.setText(str(tmp_path)); bad.make()
    assert bad.spec is None and boxes and "spec.json" in boxes[-1][1]


def test_stages_dialog_generates(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    dlg = win.stages_dialog()
    dlg.set_rows([{"name": "min", "type": "geometry_optimization", "steps": "50"},
                  {"name": "nvt", "type": "molecular_dynamics", "ensemble": "NVT", "thermostat": "berendsen", "temperature": "300", "steps": "100"},
                  {"name": "nve", "type": "molecular_dynamics", "ensemble": "NVE", "velocities": "yes"}])
    data = dlg.stages_data("single_point")
    assert data["stages"][0] == {"name": "min", "task": {"type": "geometry_optimization", "max_steps": 50}}
    assert data["stages"][1]["task"]["md"] == {"ensemble": "NVT", "thermostat": "berendsen", "temperature_k": 300.0, "steps": 100}
    assert data["stages"][2]["velocities"] is True
    dlg.folder.setText(str(tmp_path / "staged"))
    assert str(tmp_path / "staged") in dlg.output_path.text()
    dlg.generate()
    assert [d.name for d in dlg.dirs] == ["stage_01_min", "stage_02_nvt", "stage_03_nve"]
    assert (tmp_path / "staged" / "submit.sh").is_file() and (tmp_path / "staged" / "stages.json").is_file()
    assert boxes[-1][0] == "information"
    back = win.stages_dialog(); back.load_file(tmp_path / "staged" / "stages.json")
    assert [r["name"] for r in back.rows()] == ["min", "nvt", "nve"] and back.rows()[1]["temperature"] == "300"


def test_stages_dialog_reports_bad_cells(app, boxes, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    dlg = win.stages_dialog()
    dlg.set_rows([{"name": "a", "temperature": "hot"}])
    dlg.folder.setText(str(tmp_path / "s")); dlg.generate()
    assert dlg.out_dir is None and boxes[-1][0] == "critical" and "hot" in boxes[-1][1]


def test_template_save_load_and_marks(app, boxes, sk_root, conf):
    win = make_window(sk_root, conf)
    win.task.type.setCurrentIndex(win.task.type.findData("molecular_dynamics")); win.task.temperature.setValue(350)
    dlg = win.template_save_dialog(); dlg.name.setText("lab-md"); dlg.comment.setText("研究室の MD"); dlg.save()
    assert dlg.path == conf / "lab" / "lab-md.json" and dlg.path.is_file()
    win.task.temperature.setValue(300); win.task.type.setCurrentIndex(0)
    ld = win.template_load_dialog()
    assert any(i.name == "lab-md" for i in ld.infos)
    assert win.load_template("lab-md")
    s = win.current_spec()
    assert s.task.md.temperature_k == 350 and s.meta.template["name"] == "lab-md" and s.structure.atoms.symbols == ["O", "H", "H"]
    marked = {w.property("adit_key") for w in win._marked}
    assert "温度 [K]" in marked and "SCC の収束判定 (SccTolerance)" in marked
    win.task.temperature.setValue(360); win.refresh_preview()
    assert "温度 [K]" not in {w.property("adit_key") for w in win._marked}
    temp_label = next(w for w in win.task.findChildren(QLabel) if w.property("adit_key") == "温度 [K]")
    assert "雛形" not in temp_label.text()


def test_settings_dialog_has_templates_dir(app, boxes, tmp_path):
    from adit.config import load_config
    from adit.gui.settings_dialog import SettingsDialog
    path = tmp_path / "cluster.toml"
    dlg = SettingsDialog(path)
    dlg.templates_dir.setText(str(tmp_path / "lab"))
    dlg._save()
    assert load_config(path).templates_dir == str(tmp_path / "lab")


def test_provenance_is_shown_small(app, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    assert win.preview.prov.text().startswith("作成時の記録: ADIT") and "SHA-256" in win.preview.prov.text()
    win.preview.btn_prov.setChecked(True)
    assert "skf/H-O.skf" in win.preview.prov_text.toPlainText()


def test_ribbon_has_prep_buttons(app, sk_root, tmp_path):
    win = make_window(sk_root, tmp_path)
    titles = win.ribbon.page_titles()
    file_actions = win.ribbon.actions_on_page(titles.index("ファイル"))
    run_actions = win.ribbon.actions_on_page(titles.index("実行"))
    assert win.act_continue in file_actions and win.act_template_load in file_actions and win.act_template_save in file_actions
    assert run_actions.index(win.act_stages) == run_actions.index(win.act_scan) + 1


def test_new_labels_have_help_and_english():
    for key in NEW_LABELS:
        h = help_for(key)
        assert h is not None and h.ja and h.en, key
        assert key in _EN, key


def test_english_labels(app, sk_root, tmp_path):
    from adit.gui import i18n
    from adit.lang import set_language
    i18n.set_language("en")
    try:
        win = make_window(sk_root, tmp_path)
        i18n.translate_widgets(win)
        use_code(win, "xtb")
        texts = {w.text() for w in win.method.xtb.findChildren(QLabel)}
        assert "Solvation model (--alpb / --gbsa)" in texts and "Solvent" in texts
        assert win.act_continue.iconText() == "Continue" and win.act_stages.iconText() == "Stages"
        dlg = win.stages_dialog()
        assert dlg.windowTitle() == "Staged calculation" and dlg.table.horizontalHeaderItem(4).text() == "Temperature [K]"
    finally:
        i18n.set_language("ja"); set_language("ja")
