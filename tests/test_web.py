
from __future__ import annotations

import json
import re
import threading
import urllib.parse
import urllib.request
from pathlib import Path

import pytest
from ase.io import read

from adit.config import default_config
from adit.spec import Task
from adit.web.forms import FormError, default_form, form_from_spec, spec_from_form
from adit.web.server import WebApp, serve
from tests.conftest import pbs_profile, water_spec

REPO = Path(__file__).resolve().parent.parent


def test_form_roundtrip_dftb():
    f = default_form(); f["sk_set"] = "fake-1-0"
    s = spec_from_form(f)
    assert s.method.code == "dftbplus" and s.kpoints is None and s.structure.atoms.symbols == ["O", "H", "H"]
    s2 = spec_from_form({**f, **form_from_spec(s)})
    assert s2.model_dump(exclude={"meta"}) == s.model_dump(exclude={"meta"})


def test_form_vasp_bulk_extra_incar():
    f = default_form(); f.update(source="bulk", bulk_el="Si", code="vasp", kp_mode="mesh", k1="4", k2="4", k3="4",
                                 extra_incar="NSW = 10\nLWAVE = .FALSE.", potcar_map="Si = Si_sv", task_type="molecular_dynamics", thermostat="langevin")
    s = spec_from_form(f)
    assert s.structure.periodic and s.kpoints.mesh == (4, 4, 4) and s.method.extra_incar == {"NSW": 10, "LWAVE": False}
    assert s.method.potcar == {"Si": "Si_sv"} and s.task.md.thermostat == "langevin"
    back = form_from_spec(s)
    assert back["bulk_el"] == "Si" and back["k1"] == "4" and "LWAVE = .FALSE." in back["extra_incar"]


def test_form_vasp_advanced_fields_roundtrip():
    f = default_form(); f.update(code="vasp", nelmin="5", lasph="on", lmaxmix="4", nbands="96", isym="0",
                                 idipol="3", ldipol="on", dipol="0.5 0.5 0.4", kpoints_centering="gamma")
    method = spec_from_form(f).method
    assert (method.nelmin, method.lasph, method.lmaxmix, method.nbands, method.isym) == (5, True, 4, 96, 0)
    assert method.idipol == 3 and method.ldipol and method.dipol == "0.5 0.5 0.4" and method.kpoints_centering == "gamma"
    back = form_from_spec(spec_from_form(f))
    assert back["isym"] == "0" and back["ldipol"] == "on" and back["kpoints_centering"] == "gamma"


def test_form_errors_are_readable():
    f = default_form(); f["fixed"] = "99"
    with pytest.raises(FormError):
        spec_from_form(f)
    f = default_form(); f["charge"] = "abc"
    with pytest.raises(FormError):
        spec_from_form(f)


def test_form_from_existing_spec_task_fields():
    s = water_spec(task=Task(type="molecular_dynamics"))
    f = form_from_spec(s)
    assert f["task_type"] == "molecular_dynamics" and f["md_steps"] == "1000" and f["source"] == "preset"


@pytest.fixture
def web(sk_root, tmp_path):
    cfg = default_config(sk_root=str(sk_root)); cfg.profiles["cluster"] = pbs_profile()
    app = WebApp(cfg, tmp_path / "cluster.toml")
    httpd = serve(app, "127.0.0.1", 0)
    t = threading.Thread(target=httpd.serve_forever, daemon=True); t.start()
    yield app, f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown(); httpd.server_close()


def _post(url: str, fields: dict) -> str:
    return urllib.request.urlopen(url, data=urllib.parse.urlencode(fields).encode()).read().decode()


def _error(html: str) -> str:
    m = re.search(r'<div class="error">(.*?)</div>', html, re.S)
    return m.group(1) if m else ""


def test_convert_page_and_structure_route(web, tmp_path):
    app, base = web
    html = urllib.request.urlopen(base + "/convert").read().decode()
    assert "構造ファイルの形式を変換" in html and "別の計算コード用に変換" in html
    source = tmp_path / "water.xyz"
    source.write_text("3\nC3H4O3ではなくH2O\nO 0 0 0\nH 0.8 0.6 0\nH -0.8 0.6 0\n", encoding="utf-8")
    output = tmp_path / "water.cif"
    html = _post(base + "/convert_structure", {"source": str(source), "struct_output": str(output), "cell": "15"})
    assert "構造を書きました" in html and output.is_file()

    html = _post(base + "/convert_structure", {"source": "", "struct_output": str(tmp_path / "missing.cif")})
    assert "変換元の構造ファイルを指定してください" in html

    app.spec = water_spec().model_copy(update={"structure": water_spec().structure.model_copy(update={"fixed_atoms": [0]})})
    current = tmp_path / "current.extxyz"
    html = _post(base + "/save_current_structure", {"current_output": str(current)})
    assert "現在の構造を書きました" in html and current.is_file()
    saved = read(current)
    assert len(saved) == 3 and saved.constraints
    assert "すでにあります" in _post(base + "/save_current_structure", {"current_output": str(current)})


def test_index_and_preview_and_generate(web, tmp_path):
    app, base = web
    html = urllib.request.urlopen(base + "/").read().decode()
    assert "fake-1-0" in html and 'name="task_type"' in html
    f = dict(app.form); f.update(sk_set="fake-1-0", task_type="geometry_optimization", output_dir=str(tmp_path / "out"))
    html = _post(base + "/preview", f)
    assert not _error(html) and "dftb_in.hsd" in html and "MaxSteps = 200" in html
    assert app.spec is not None and app.spec.task.type == "geometry_optimization"
    spec_json = json.loads(urllib.request.urlopen(base + "/spec.json").read())
    assert spec_json["method"]["sk_set"] == "fake-1-0"
    html = _post(base + "/generate", f)
    assert not _error(html) and (tmp_path / "out" / "dftb_in.hsd").is_file() and (tmp_path / "out" / "README.txt").is_file()
    assert _error(_post(base + "/generate", f))
    assert not _error(_post(base + "/generate", {**f, "overwrite": "on"}))


def test_preview_shows_validation_errors(web):
    app, base = web
    f = dict(app.form); f.update(sk_set="nodocs-0-1")
    html = _post(base + "/preview", f)
    assert _error(html) and "LICENSE" in _error(html)


def test_cluster_profile_is_not_run(web, tmp_path):
    app, base = web
    f = dict(app.form); f.update(sk_set="fake-1-0", profile="cluster", ncpus="8", omp="8", output_dir=str(tmp_path / "c"))
    html = _post(base + "/generate", f)
    assert not _error(html) and "#PBS" in (tmp_path / "c" / "submit.sh").read_text(encoding="utf-8")
    assert "transfer_and_submit.sh" in html        # the person runs it in a terminal


def test_analysis_page_and_figures(web, tmp_path):
    app, base = web
    import shutil
    run_dir = tmp_path / "md"
    shutil.copytree(REPO / "examples" / "dftb_md_water_generated", run_dir)
    html = _post(base + "/analysis", {"run_dir": str(run_dir), "rdf": "on"})
    assert not _error(html) and "dftbplus" in html
    figs = re.findall(r'src="/file\?path=([^"]+)"', html)
    assert figs
    data = urllib.request.urlopen(base + "/file?path=" + figs[0]).read()
    assert data[:4] == b"\x89PNG"
    with pytest.raises(urllib.error.HTTPError):
        urllib.request.urlopen(base + "/file?path=" + urllib.parse.quote(str(REPO / "README.md")))


def test_load_spec_json_multipart(web):
    app, base = web
    spec_bytes = (REPO / "examples" / "qe_si_generated" / "spec.json").read_bytes()
    b = "----adittest"
    body = (f"--{b}\r\nContent-Disposition: form-data; name=\"spec_file\"; filename=\"spec.json\"\r\nContent-Type: application/json\r\n\r\n").encode() + spec_bytes + f"\r\n--{b}--\r\n".encode()
    req = urllib.request.Request(base + "/load", data=body, headers={"Content-Type": f"multipart/form-data; boundary={b}"})
    html = urllib.request.urlopen(req).read().decode()
    assert app.form["code"] == "espresso" and app.form["source"] == "bulk" and app.form["bulk_el"] == "Si"
    assert 'id="ecutwfc"' in html


def test_load_cli_only_spec_does_not_replace_web_form(web):
    from adit.spec import GaussianMethod, Runtime, Task
    from tests.conftest import water_spec

    app, base = web
    previous_code = app.form["code"]
    spec_bytes = water_spec(method=GaussianMethod(theory="HF", basis="6-31G(d)"),
                            task=Task(type="single_point"), runtime=Runtime(profile="local")).to_json().encode()
    boundary = "----aditclitest"
    body = (f"--{boundary}\r\nContent-Disposition: form-data; name=\"spec_file\"; filename=\"spec.json\"\r\n"
            "Content-Type: application/json\r\n\r\n").encode() + spec_bytes + f"\r\n--{boundary}--\r\n".encode()
    req = urllib.request.Request(base + "/load", data=body,
                                 headers={"Content-Type": f"multipart/form-data; boundary={boundary}"})
    html = urllib.request.urlopen(req).read().decode()
    assert app.form["code"] == previous_code
    assert "adit-gen" in html


def _fields_without_label(html: str) -> list[str]:
    from html.parser import HTMLParser

    class P(HTMLParser):
        def __init__(self):
            super().__init__(); self.fields, self.fors = [], set()

        def handle_starttag(self, tag, attrs):
            a = dict(attrs)
            if tag == "label" and a.get("for"):
                self.fors.add(a["for"])
            elif tag in ("input", "select", "textarea") and a.get("type") != "hidden":
                self.fields.append(a)

    p = P(); p.feed(html)
    return [f.get("name") or "?" for f in p.fields
            if not f.get("aria-label") and f.get("id") not in p.fors and f.get("name") != "overwrite"]


def test_every_field_has_a_label(web):
    app, base = web
    html = urllib.request.urlopen(base + "/").read().decode()
    assert html.count("<input") > 50
    assert _fields_without_label(html) == []


def test_only_the_chosen_fields_are_shown_by_css(web):
    app, base = web
    html = urllib.request.urlopen(base + "/").read().decode()
    for attr, sel, v in [("source", "source", "bulk"), ("source", "source", "mixture"), ("code", "code", "vasp"),
                         ("task", "task_type", "molecular_dynamics"), ("kp", "kp_mode", "mesh")]:
        assert f'class="only" data-{attr}="{v}"' in html
        assert f'#main:has(#{sel} option[value="{v}"]:checked) .only[data-{attr}~="{v}"]' in html
    assert 'name="source" id="source"' in html and 'type="radio"' not in html
    assert "show('code'" not in html


def test_direct_profile_hides_cluster_rows(web):
    app, base = web
    html = urllib.request.urlopen(base + "/").read().decode()
    assert 'data-kind="direct"' in html and 'data-kind="pbs"' in html
    assert '#main:has(#profile option[data-kind="direct"]:checked) .cluster-only { display:none; }' in html
    for fid in ("nodes", "walltime", "job_name"):
        assert re.search(r'<div class="row cluster-only"><label for="' + fid + '"', html), fid
    for fid in ("ncpus", "omp", "mpiprocs"):
        assert not re.search(r'<div class="row cluster-only"><label for="' + fid + '"', html), fid


def test_default_cores_match_desktop():
    f = default_form()
    assert f["ncpus"] == "8" and f["omp"] == "8"
    f["sk_set"] = "fake-1-0"
    r = spec_from_form(f).runtime
    assert r.ncpus == 8 and r.omp_threads == 8
    r = spec_from_form({**f, "ncpus": "", "omp": ""}).runtime
    assert r.ncpus == 8 and r.omp_threads == 8


def test_sk_label_follows_language(web):
    from adit import lang
    app, base = web
    before = lang.LANGUAGE
    try:
        lang.set_language("ja")
        html = urllib.request.urlopen(base + "/").read().decode()
        assert "Slater-Koster パラメータ" in html and "Slater-Koster set" not in html and "制限時間 (HH:MM:SS)" in html
        lang.set_language("en")
        html = urllib.request.urlopen(base + "/").read().decode()
        assert "Slater-Koster set" in html and "Slater-Koster パラメータ" not in html and "Time limit (HH:MM:SS)" in html
    finally:
        lang.set_language(before)


def test_error_headline_like_desktop(web):
    from adit import lang
    app, base = web
    before = lang.LANGUAGE
    try:
        lang.set_language("ja")
        f = dict(app.form); f.update(sk_set="nodocs-0-1", multiplicity="2")
        html = _post(base + "/preview", f)
        m = re.search(r'<div class="error-banner" role="alert">(.*?)</div>', html, re.S)
        assert m and app.errors
        import html as h
        head = h.unescape(m.group(1))
        assert head.startswith("生成できません: " + str(app.errors[0]))
        assert ("ほか" in head) == (len(app.errors) > 1)
        assert html.index("error-banner") < html.index('<fieldset>')
    finally:
        lang.set_language(before)


def test_mixture_form(web, tmp_path):
    app, base = web
    f = dict(app.form); f.update(sk_set="fake-1-0", source="mixture", mixture="H2O * 6\nNH3 * 1", mix_box_mode="edge", mix_edge="9", kp_mode="gamma")
    html = _post(base + "/preview", f)
    assert not _error(html), _error(html)
    assert app.spec.structure.source == "mixture" and app.spec.structure.periodic and len(app.spec.structure.atoms.symbols) == 22
    back = form_from_spec(app.spec)
    assert "H2O * 6" in back["mixture"] and back["mix_box_mode"] == "edge"
