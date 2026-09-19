
from __future__ import annotations

import re
import threading
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

from adit.config import default_config
from adit.web.forms import default_form, form_from_spec, spec_from_form
from adit.web.server import WebApp, serve
from tests.test_cp2k import BASIS, POTENTIALS
from tests.test_web import _fields_without_label

REPO = Path(__file__).resolve().parent.parent
SPCE = REPO / "examples" / "gromacs_spce"


@pytest.fixture
def cp2k_data(tmp_path) -> Path:
    d = tmp_path / "cp2k_data"; d.mkdir()
    (d / "BASIS_MOLOPT").write_text(BASIS, encoding="utf-8"); (d / "GTH_POTENTIALS").write_text(POTENTIALS, encoding="utf-8")
    (d / "dftd3.dat").write_text("fake\n", encoding="utf-8")
    return d


@pytest.fixture
def web(sk_root, tmp_path, cp2k_data):
    cfg = default_config(sk_root=str(sk_root)); cfg.cp2k_data = str(cp2k_data)
    app = WebApp(cfg, tmp_path / "cluster.toml")
    app.form["output_dir"] = str(tmp_path / "out")
    httpd = serve(app, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    yield app, f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown(); httpd.server_close()


def _post(url: str, fields: dict) -> str:
    return urllib.request.urlopen(url, data=urllib.parse.urlencode(fields).encode()).read().decode()


def _error(html: str) -> str:
    m = re.search(r'<div class="error">(.*?)</div>', html, re.S)
    return m.group(1) if m else ""


def cp2k_form() -> dict:
    f = default_form()
    f.update(code="cp2k", cp_xc="PBE", cp_cutoff="400", cp_rel_cutoff="60", cp_poisson="MT", cp_box="10", cp_uks="",
             cp_basis_map="O = DZVP-MOLOPT-SR-GTH\nH = DZVP-MOLOPT-SR-GTH", cp_extra="[FORCE_EVAL/DFT/SCF]\nSCF_GUESS ATOMIC", cp_dispersion="d3")
    return f


def lammps_form(tmp_path) -> dict:
    eam = tmp_path / "Cu_u3.eam"; eam.write_text("fake eam\n", encoding="utf-8")
    f = default_form()
    f.update(source="bulk", bulk_el="Cu", bulk_cubic="on", code="lammps", lmp_units="metal", lmp_pair_style="eam", lmp_pair_coeff="* * Cu_u3.eam",
             lmp_files=str(eam), lmp_seed="777", task_type="molecular_dynamics", kp_mode="mesh", k1="4", k2="4", k3="4")
    return f


def gromacs_form() -> dict:
    f = default_form()
    f.update(code="gromacs", gmx_top=str(SPCE / "topol.top"), gmx_conf=str(SPCE / "conf.gro"), gmx_use_conf="on", gmx_coulomb="PME",
             gmx_rcoulomb="0.9", gmx_rvdw="0.9", gmx_extra="nstcalcenergy = 100\nlincs-order = 4")
    return f


@pytest.mark.parametrize("make", ["cp2k", "lammps", "gromacs"])
def test_form_roundtrip(make, tmp_path):
    f = {"cp2k": cp2k_form, "lammps": lambda: lammps_form(tmp_path), "gromacs": gromacs_form}[make]()
    s = spec_from_form(f)
    assert s.method.code == make
    s2 = spec_from_form({**f, **form_from_spec(s)})
    assert s2.model_dump(exclude={"meta"}) == s.model_dump(exclude={"meta"})


def test_values_reach_the_spec(tmp_path):
    s = spec_from_form(cp2k_form())
    assert s.method.extra_sections == {"FORCE_EVAL/DFT/SCF": "SCF_GUESS ATOMIC"} and s.method.basis["O"] == "DZVP-MOLOPT-SR-GTH"
    s = spec_from_form(lammps_form(tmp_path))
    assert s.structure.periodic and s.kpoints is None and s.method.seed == 777
    s = spec_from_form(gromacs_form())
    assert s.structure.source == "file" and len(s.structure.atoms.symbols) == 648 and s.method.extra_mdp == {"nstcalcenergy": 100, "lincs-order": 4}
    assert form_from_spec(s)["gmx_use_conf"] == "on"


def _load(base: str, spec_json: bytes) -> str:
    b = "----adittest"
    body = (f"--{b}\r\nContent-Disposition: form-data; name=\"spec_file\"; filename=\"spec.json\"\r\nContent-Type: application/json\r\n\r\n").encode() + spec_json + f"\r\n--{b}--\r\n".encode()
    req = urllib.request.Request(base + "/load", data=body, headers={"Content-Type": f"multipart/form-data; boundary={b}"})
    return urllib.request.urlopen(req).read().decode()


@pytest.mark.parametrize("make", ["cp2k", "lammps", "gromacs"])
def test_load_spec_json(web, tmp_path, make):
    app, base = web
    f = {"cp2k": cp2k_form, "lammps": lambda: lammps_form(tmp_path), "gromacs": gromacs_form}[make]()
    spec = spec_from_form(f)
    html = _load(base, spec.model_dump_json().encode())
    assert "spec.json を読み込みました" in html and app.form["code"] == make
    assert app.spec is not None and app.spec.model_dump(exclude={"meta"}) == spec.model_dump(exclude={"meta"})


def test_cp2k_candidates_and_preview(web):
    app, base = web
    html = _post(base + "/preview", {**app.form, **cp2k_form()})
    assert '<div class="status ok">生成できます</div>' in html, _error(html)
    assert "DZVP-MOLOPT-SR-GTH, SZV-MOLOPT-SR-GTH" in html
    assert "BASIS_SET DZVP-MOLOPT-SR-GTH" in html


def test_kpoints_group_is_hidden_for_classical_md(web):
    app, base = web
    html = urllib.request.urlopen(base + "/").read().decode()
    assert '<fieldset class="kp-group">' in html
    for c in ("lammps", "gromacs"):
        assert f'#main:has(#code option[value="{c}"]:checked) .kp-group {{ display:none; }}' in html
    for c in ("cp2k", "lammps", "gromacs"):
        assert f'class="only" data-code="{c}"' in html and f'.only[data-code~="{c}"]' in html


def test_required_field_empty_shows_reason(web, tmp_path):
    app, base = web
    f = lammps_form(tmp_path); f["lmp_pair_style"] = ""
    html = _post(base + "/preview", {**app.form, **f})
    assert "pair_style" in _error(html) and 'class="error-banner"' in html
    f = cp2k_form(); f["task_type"] = "band_structure"; f["cp_xc"] = ""
    err = _error(_post(base + "/preview", {**app.form, **f}))
    assert "バンド計算はありません" in err and "汎関数" in err


def test_labels_help_and_english(web):
    from adit import lang
    app, base = web
    html = urllib.request.urlopen(base + "/").read().decode()
    assert _fields_without_label(html) == []
    pill = '<span class="pill">必須</span></label>'
    assert re.search(r'<label for="cp_xc" class="required" title="[^"]+">汎関数' + re.escape(pill), html)
    assert re.search(r'<label for="gmx_top" class="required" title="[^"]+">トポロジー \(.top\)' + re.escape(pill), html)
    assert "青字は必須" not in html      # the pill replaced the legend
    lang.set_language("en")
    try:
        html = urllib.request.urlopen(base + "/").read().decode()
        for t in ("Functional", "Basis and pseudopotential per element", "Units", "Files to copy", "Topology (.top)", "Compressibility [1/bar]"):
            assert f">{t}</label>" in html or f'>{t}<span class="pill">required</span></label>' in html, t
    finally:
        lang.set_language("ja")


def test_narrow_width_rules(web):
    app, base = web
    html = urllib.request.urlopen(base + "/").read().decode()
    assert "fieldset, .cols > div { min-width:0; }" in html
    media = html[html.index("@media (max-width: 600px)"):]
    assert ".row { grid-template-columns: minmax(0,1fr); gap:2px; }" in media[:600]
    assert ".row { display:grid; grid-template-columns: 180px minmax(0,1fr);" in html
