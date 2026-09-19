
from __future__ import annotations

import html as h
import json
import re
import shutil
import threading
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

from adit.codes.orca import solvent_names as orca_names
from adit.codes.xtb import solvent_names as xtb_names
from adit.config import default_config
from adit.spec import Cp2kMethod, DftbMethod, EspressoMethod, HubbardU, OrcaMethod, VaspMethod, XtbMethod
from adit.web.forms import FormError, default_form, form_from_spec, spec_from_form
from adit.web.server import WebApp, serve
from tests.conftest import water_spec
from tests.test_web import _fields_without_label

REPO = Path(__file__).resolve().parent.parent
EX = REPO / "examples"
CHECKBOXES = ["scc", "third", "bulk_cubic", "box", "cp_uks", "gmx_use_conf", "cont_vel", "cont_keep", "stg_overwrite", "tmpl_overwrite"]


@pytest.fixture
def web(sk_root, tmp_path, monkeypatch):
    monkeypatch.setenv("ADIT_CONFIG", str(tmp_path / "conf" / "cluster.toml"))
    cfg = default_config(sk_root=str(sk_root)); cfg.templates_dir = str(tmp_path / "lab")
    app = WebApp(cfg, tmp_path / "cluster.toml")
    app.form["output_dir"] = str(tmp_path / "out")
    httpd = serve(app, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    yield app, f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown(); httpd.server_close()


def _post(url: str, fields) -> str:
    return urllib.request.urlopen(url, data=urllib.parse.urlencode(fields).encode(), timeout=60).read().decode()


def _error(page: str) -> str:
    m = re.search(r'<div class="error">(.*?)</div>', page, re.S)
    return h.unescape(m.group(1)) if m else ""


def _base(app, **kw) -> dict:
    f = {**default_form(output_dir=app.form["output_dir"]), "sk_set": "fake-1-0"}
    f.update(kw)
    return f


def test_every_checkbox_has_an_empty_hidden_twin(web):
    _, base = web
    page = urllib.request.urlopen(base + "/").read().decode()
    for name in CHECKBOXES:
        assert re.search(rf'<input type="hidden" name="{name}" value=""><input type="checkbox" name="{name}"', page), name


def test_unticking_scc_reaches_the_spec(web):
    app, base = web
    f = _base(app); f.pop("scc")
    _post(base + "/preview", list(f.items()) + [("scc", "")])
    assert app.spec is not None and app.spec.method.scc is False
    _post(base + "/preview", list(f.items()) + [("scc", ""), ("scc", "on")])
    assert app.spec.method.scc is True
    f2 = _base(app, source="bulk", bulk_el="Si", box=""); f2.pop("scc")
    _post(base + "/preview", list(f2.items()) + [("bulk_cubic", ""), ("bulk_cubic", "on"), ("scc", "")])
    assert app.spec.structure.source_ref.endswith("cubic") and app.spec.method.scc is False


def test_solvent_form_roundtrip():
    for m in (XtbMethod(gfn="1", solvation="gbsa", solvent="chcl3"), OrcaMethod(solvation="cpcm", solvent="water"),
              DftbMethod(sk_set="fake-1-0", solvation_param_file="/tmp/param_gbsa_h2o.txt")):
        s = water_spec(method=m)
        back = spec_from_form({**default_form(), **form_from_spec(s)})
        assert back.method == s.method, m.code


def test_solvent_reads_only_the_chosen_combination():
    f = {**default_form(), "code": "xtb", "gfn": "2", "xtb_solvation": "alpb", "xtb_solvent__alpb_2": "water", "xtb_solvent__gbsa_2": "dmso"}
    m = spec_from_form(f).method
    assert (m.solvation, m.solvent) == ("alpb", "water")
    f["xtb_solvation"] = "none"
    assert spec_from_form(f).method.solvent == ""


def test_solvent_candidates_in_page(web):
    app, base = web
    page = _post(base + "/preview", _base(app, code="xtb", xtb_solvation="gbsa", gfn="1"))
    sel = re.search(r'<select name="xtb_solvent__gbsa_1".*?</select>', page, re.S).group(0)
    assert sel.count("<option") == len(xtb_names("gbsa", "1")) + 1
    osel = re.search(r'<select name="orca_solvent__smd".*?</select>', page, re.S).group(0)
    assert osel.count("<option") == len(orca_names("smd")) + 1
    assert '#xtb_solvation option[value="gbsa"]:checked):has(#gfn option[value="1"]:checked) .xsolv[data-k="gbsa_1"]' in page


@pytest.mark.parametrize("method", [
    VaspMethod(ispin=2, magmom_by_element={"O": 1.0}, hubbard={"O": HubbardU(orbital="2p", u_ev=3.0, j_ev=0.5)}, ldau_type=4),
    EspressoMethod(pseudo_set="x", nspin=2, starting_magnetization={"O": 0.2}, hubbard={"O": HubbardU(orbital="2p", u_ev=2.0)}, hubbard_projector="wf"),
    Cp2kMethod(uks=True, magnetization_by_element={"O": 2.0}, hubbard={"O": HubbardU(orbital="2p", u_ev=4.0, j_ev=1.0)},
               plus_u_method="LOWDIN", sccs_relative_permittivity=78.4),
])
def test_magnetism_dftu_form_roundtrip(method):
    s = water_spec(method=method)
    back = spec_from_form({**default_form(), **form_from_spec(s)})
    for k in ("magmom_by_element", "starting_magnetization", "magnetization_by_element", "hubbard", "ldau_type", "hubbard_projector",
              "plus_u_method", "sccs_relative_permittivity"):
        if hasattr(method, k):
            assert getattr(back.method, k) == getattr(method, k), k


def test_dftu_needs_u_and_uses_structure_elements(web):
    app, base = web
    with pytest.raises(FormError, match="U"):
        spec_from_form({**default_form(), "code": "vasp", "vasp_hub1_el": "O", "vasp_hub1_orb": "2p", "vasp_hub1_u": ""})
    page = _post(base + "/preview", _base(app, code="vasp"))
    sel = re.search(r'<select name="vasp_hub1_el".*?</select>', page, re.S).group(0)
    assert re.findall(r'<option value="(\w*)"', sel) == ["", "O", "H"]
    assert 'name="vasp_mag__O"' in page and 'name="vasp_mag__Fe"' not in page
    assert _fields_without_label(page) == []


@pytest.fixture
def prev_md(tmp_path) -> Path:
    d = tmp_path / "prev_md"
    shutil.copytree(EX / "dftb_md_water_generated", d)
    return d


def test_continue_form_carries_velocities(web, prev_md):
    app, base = web
    page = _post(base + "/continue", list(_base(app, cont_dir=str(prev_md)).items()) + [("cont_vel", ""), ("cont_vel", "on")])
    assert "geo_end.xyz" in h.unescape(page) and app.spec is not None
    assert app.spec.structure.velocities is not None and app.spec.handoff.velocities
    assert app.form["file_path"].endswith("geo_end.xyz") and app.form["source"] == "file"
    page = _post(base + "/preview", {**app.form, "sk_set": "fake-1-0"})
    assert "Velocities [AA/ps]" in h.unescape(page) and str(prev_md) in h.unescape(page)
    _post(base + "/origin_clear", {**app.form, "sk_set": "fake-1-0"})
    assert app.spec.handoff is None


def test_continue_form_without_velocities_and_errors(web, prev_md, tmp_path):
    app, base = web
    _post(base + "/continue", list(_base(app, cont_dir=str(prev_md)).items()) + [("cont_vel", "")])
    assert app.spec.structure.velocities is None
    page = _post(base + "/continue", _base(app, cont_dir=str(tmp_path / "nothing")))
    assert "spec.json" in _error(page)


def test_stages_form_generates(web, tmp_path):
    app, base = web
    rows = {"stg_rows": "3", "stg1_name": "min", "stg1_type": "geometry_optimization", "stg1_steps": "50",
            "stg2_name": "nvt", "stg2_type": "molecular_dynamics", "stg2_ensemble": "NVT", "stg2_temperature": "300", "stg2_steps": "100",
            "stg3_name": "", "stg3_velocities": "auto", "stg_dir": str(tmp_path / "staged")}
    page = _post(base + "/stages", _base(app, **rows))
    assert not _error(page), _error(page)
    assert [p.name for p in sorted((tmp_path / "staged").glob("stage_*"))] == ["stage_01_min", "stage_02_nvt"]
    data = json.loads((tmp_path / "staged" / "stages.json").read_text(encoding="utf-8"))
    assert data["stages"][1]["task"]["md"]["temperature_k"] == 300.0
    page = _post(base + "/stages", _base(app, **rows, stages_action="add"))
    assert 'name="stg4_name"' in page and app.form["stg_rows"] == "4"
    page = _post(base + "/stages", _base(app, **{**rows, "stg1_temperature": "hot"}))
    assert "hot" in h.unescape(page)


def test_template_save_and_load(web, tmp_path):
    app, base = web
    page = _post(base + "/template_save", _base(app, task_type="molecular_dynamics", temperature="350", tmpl_save_name="lab-md", tmpl_comment="MD"))
    assert (tmp_path / "lab" / "lab-md.json").is_file(), _error(page)
    page = _post(base + "/template_save", _base(app, tmpl_save_name="lab-md"))
    assert "上書き" in h.unescape(page)
    page = _post(base + "/", {}) if False else urllib.request.urlopen(base + "/").read().decode()
    assert "lab-md — dftbplus" in h.unescape(page)
    page = _post(base + "/template_load", _base(app, preset="NH3", tmpl_name="lab-md"))
    assert app.spec.task.md.temperature_k == 350 and app.spec.meta.template["name"] == "lab-md"
    assert app.spec.structure.atoms.symbols.count("N") == 1
    assert 'class="tmark"' in page and "雛形 lab-md の値" in h.unescape(page)


def test_provenance_in_preview(web):
    app, base = web
    page = h.unescape(_post(base + "/preview", _base(app)))
    assert "作成時の記録: ADIT" in page and "skf/H-O.skf" in page


def test_english_page(web):
    from adit import lang
    app, base = web
    lang.set_language("en")
    try:
        page = h.unescape(_post(base + "/preview", _base(app, code="xtb")))
        for text in ("Solvation model (--alpb / --gbsa)", "Staged calculation", "Continue a previous calculation", "Group templates",
                     "Provenance: ADIT"):
            assert text in page, text
    finally:
        lang.set_language("ja")


def _hidden_cases():
    # Specs whose values the form cannot show; everything else is form-shaped already.
    from adit.spec import KPoints, LammpsMethod, Meta

    xtb = spec_from_form({**default_form(), "code": "xtb", "task_type": "molecular_dynamics"})
    qe = spec_from_form({**default_form(), "code": "espresso", "source": "bulk", "bulk_el": "Si", "kp_mode": "mesh", "k1": "4", "k2": "4", "k3": "1"})
    lmp = spec_from_form({**default_form(), "code": "lammps", "source": "bulk", "bulk_el": "Cu", "lmp_units": "metal", "lmp_pair_style": "eam"})
    return [
        xtb.model_copy(update={"method": XtbMethod(md_hmass=2.0, md_shake=0, md_sccacc=1.0),
                               "meta": Meta(comment="hello", stage={"index": 1, "name": "nvt"})}),
        qe.model_copy(update={"method": EspressoMethod(dipole_correction=True, dipole_direction=3, dipole_maxpos=0.9, dipole_decrease=0.1,
                                                       dipole_amplitude=0.01),
                              "kpoints": KPoints(mode="mesh", mesh=(4, 4, 1), shift=(0.5, 0.5, 0.0))}),
        lmp.model_copy(update={"method": LammpsMethod(units="metal", pair_style="eam", thermo_pressure_tensor=True)}),
    ]


def test_fields_without_a_form_field_survive_the_round_trip():
    from adit.web.codefields import PrepOrigin

    skip = {"meta": {"created", "app_version"}}
    for spec in _hidden_cases():
        plain = spec_from_form({**default_form(), **form_from_spec(spec)})
        assert plain.model_dump(exclude=skip) != spec.model_dump(exclude=skip), spec.method.code
        origin = PrepOrigin(spec)
        assert origin.active
        back = origin.apply(plain)
        assert back.model_dump(exclude=skip) == spec.model_dump(exclude=skip), spec.method.code
        note = origin.describe(back)
        assert "読み込んだ値のまま" in note or "k 点シフトは (0.5, 0.5, 0) のまま" in note
    qe = _hidden_cases()[1]
    origin = PrepOrigin(qe)
    assert "k 点シフトは (0.5, 0.5, 0) のまま" in origin.describe(origin.apply(spec_from_form({**default_form(), **form_from_spec(qe)})))
    changed = spec_from_form({**default_form(), **form_from_spec(qe), "kp_shift": "0"})
    assert origin.apply(changed).kpoints.shift == (0.0, 0.0, 0.0)


def test_loaded_spec_json_is_written_back_unchanged(web):
    from adit.spec import CalculationSpec

    app, base = web
    spec = _hidden_cases()[0]
    b = "----adittest"
    body = (f"--{b}\r\nContent-Disposition: form-data; name=\"spec_file\"; filename=\"spec.json\"\r\nContent-Type: application/json\r\n\r\n").encode() \
        + spec.model_dump_json().encode() + f"\r\n--{b}--\r\n".encode()
    req = urllib.request.Request(base + "/load", data=body, headers={"Content-Type": f"multipart/form-data; boundary={b}"})
    page = urllib.request.urlopen(req, timeout=60).read().decode()
    assert not _error(page), _error(page)
    assert app.files is not None
    skip = {"meta": {"created", "app_version"}}
    written = CalculationSpec.from_json(app.files.texts["spec.json"])
    assert written.model_dump(exclude=skip) == spec.model_dump(exclude=skip)
    assert "md_hmass=2.0" in h.unescape(page) and "meta.comment" in h.unescape(page)
    page = _post(base + "/preview", {**app.form, "xtb_etemp": "500"})
    assert app.spec.method.md_hmass == 2.0 and app.spec.method.etemp == 500.0 and app.spec.meta.stage == {"index": 1, "name": "nvt"}
    _post(base + "/origin_clear", dict(app.form))
    assert app.spec.method.md_hmass == 4.0 and app.spec.meta.comment == ""
