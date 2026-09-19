
from __future__ import annotations

import json
import re
import threading
import urllib.parse
import urllib.request
from html import unescape

import numpy as np
import pytest
from ase import Atoms
from ase.io import write

from tests.test_web import _fields_without_label, _post
from adit.builder import Recipe
from adit.config import default_config
from adit.web.server import WebApp, serve

SKEW_AB = [[11.0, 0.0, 0.0], [-4.4, 16.5, 0.0]]


@pytest.fixture
def web(sk_root, tmp_path):
    app = WebApp(default_config(sk_root=str(sk_root)), tmp_path / "cluster.toml")
    app.form.update(sk_set="fake-1-0", output_dir=str(tmp_path / "out"))
    httpd = serve(app, "127.0.0.1", 0)
    t = threading.Thread(target=httpd.serve_forever, daemon=True); t.start()
    yield app, f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown(); httpd.server_close()


def act(app, base, action: str, **fields) -> str:
    return _post(base + "/recipe", {**app.form, **fields, "recipe_action": action})


def ops(app) -> list[str]:
    return [d["op"] for d in json.loads(app.form["recipe_steps"])]


def list_rows(html: str) -> list[str]:
    return [unescape(re.sub(r"\s+", " ", x)).strip() for x in re.findall(r'<li[^>]*>(?:<input[^>]*>)?<a href="#step\d+">(.*?)</a>', html, re.S)]


def status(html: str) -> str:
    m = re.search(r'<div class="([^"]*) recipe-status" role="status">(.*?)</div>', html, re.S)
    return unescape(m.group(2)) if m else ""


def log_lines(html: str) -> list[str]:
    m = re.search(r'<pre class="recipe-log">(.*?)</pre>', html, re.S)
    return unescape(m.group(1)).splitlines() if m else []


def _error(html: str) -> str:
    m = re.search(r'<div class="error">(.*?)</div>', html, re.S)
    return unescape(m.group(1)) if m else ""


@pytest.fixture(autouse=True)
def _ja():
    from adit import lang
    before = lang.LANGUAGE
    lang.set_language("ja")
    yield
    lang.set_language(before)


def test_add_move_remove(web):
    app, base = web
    app.form.update(source="bulk", bulk_cubic="on")
    for op in ("supercell", "slab", "fix"):
        html = act(app, base, "add", recipe_add=op)
    assert ops(app) == ["supercell", "slab", "fix"]
    assert [r.split("  —  ")[0].split(" — ")[0] for r in list_rows(html)] == ["1. 超格子", "2. 面で切る", "3. 固定"]
    assert re.search(r'<details class="stepbox" id="step3" open>', html)
    html = act(app, base, "up:3")
    assert ops(app) == ["supercell", "fix", "slab"] and app.form["st2_op"] == "fix" and list_rows(html)[1].startswith("2. 固定")
    act(app, base, "down:1")
    assert ops(app) == ["fix", "supercell", "slab"]
    html = act(app, base, "del:2")
    assert ops(app) == ["fix", "slab"] and len(list_rows(html)) == 2 and "st3_op" not in app.form
    assert 'name="fixed" id="fixed" value="" disabled' in html
    html = _post(base + "/preview", app.form)
    assert app.spec is None and "作る" in unescape(html)
    act(app, base, "del:1"); html = act(app, base, "del:1")
    assert ops(app) == [] and "手順を足すと" in html
    _post(base + "/preview", app.form)
    assert app.spec is not None and app.spec.structure.source == "bulk"


def test_each_step_fields_to_json(web):
    app, base = web
    app.form.update(source="bulk", bulk_cubic="on")
    for op in ("supercell", "slab", "vacuum", "remove", "substitute", "adsorb", "solvent_layer", "solvate", "fix", "box"):
        act(app, base, "add", recipe_add=op)
    edits = {"st1_mode": "matrix", "st1_m01": "1", "st1_m10": "-1",
             "st2_h": "1", "st2_k": "1", "st2_l": "1", "st2_layers": "2", "st2_vacuum": "8",
             "st3_thickness": "20",
             "st4_elements": "Si", "st4_zmin": "1.5", "st4_count": "2", "st4_seed": "4",
             "st5_pmode": "fraction", "st5_fraction": "0.25", "st5_to": "Ge",
             "st6_kind": "smiles", "st6_ref": "[C-]#[O+]", "st6_place": "above_atom", "st6_atom": "3", "st6_height": "1.9", "st6_down": "1",
             "st7_comps": "H2O * 32 label=H2O\nsmiles:[Li+] * 2 charge=+1 label=Li+", "st7_tmode": "thickness", "st7_thickness": "15", "st7_gap": "2.5",
             "st9_layers": "2", "st9_elements": "Si"}
    html = act(app, base, "noop", **edits)
    steps = json.loads(app.form["recipe_steps"])
    rec = Recipe.model_validate({"base": {"source": "bulk", "ref": "Si cubic"}, "steps": steps})
    s = rec.steps
    assert s[0].matrix == [[1, 1, 0], [-1, 1, 0], [0, 0, 1]] and s[0].repeat is None
    assert tuple(s[1].miller) == (1, 1, 1) and s[1].layers == 2 and s[1].vacuum == 8.0
    assert s[2].thickness == 20.0 and s[2].axis == 2
    assert s[3].where.elements == ["Si"] and s[3].where.z_min == 1.5 and s[3].count == 2 and s[3].seed == 4
    assert s[4].fraction == 0.25 and s[4].count is None and s[4].to == "Ge"
    assert s[5].molecule.kind == "smiles" and s[5].above_atom == 2 and s[5].height == 1.9 and s[5].down_atom == 0 and s[5].xy is None
    assert [c.ref for c in s[6].components] == ["H2O", "[Li+]"] and s[6].thickness == 15.0 and s[6].gap == 2.5
    assert s[7].components[0].count == 0
    assert s[8].bottom_layers == 2 and s[8].where.elements == ["Si"]
    assert app.form["charge"] == "2"
    rows = list_rows(html)
    assert rows[0] == "1. 超格子 — 変換行列" and rows[6].startswith("7. 溶液の層 (界面) — H₂O × 32, Li⁺ × 2")
    html = act(app, base, "up:5", st4_zmin="abc")
    assert ops(app)[3:5] == ["remove", "substitute"]
    assert status(html).startswith("手順 4 (原子を抜く): z の下限 を数値として読めません")
    assert re.search(r'<details class="stepbox" id="step4" open>', html)


def _oblique_slab(tmp_path) -> str:
    cell = np.array(SKEW_AB + [[0.0, 0.0, 20.0]])
    sym, pos = [], []
    for layer, (el, z) in enumerate([("Ni", 0.0), ("O", 2.1), ("Ni", 4.2)]):
        for i in range(4):
            for j in range(6):
                sym.append(el); pos.append((i / 4 + (layer % 2) * 0.125) * cell[0] + (j / 6) * cell[1] + [0, 0, z])
    path = tmp_path / "oblique_slab.extxyz"
    write(path, Atoms(sym, positions=pos, cell=cell, pbc=True))
    return str(path)


def _load(base: str, spec_bytes: bytes) -> str:
    b = "----adittest"
    body = (f"--{b}\r\nContent-Disposition: form-data; name=\"spec_file\"; filename=\"spec.json\"\r\nContent-Type: application/json\r\n\r\n").encode() + spec_bytes + f"\r\n--{b}--\r\n".encode()
    req = urllib.request.Request(base + "/load", data=body, headers={"Content-Type": f"multipart/form-data; boundary={b}"})
    return urllib.request.urlopen(req).read().decode()


def test_skew_cell_interface_and_roundtrip(web, tmp_path):
    pytest.importorskip("rdkit")
    app, base = web
    app.form.update(source="file", file_path=_oblique_slab(tmp_path))
    html = act(app, base, "add", recipe_add="interface")
    assert ops(app) == ["supercell", "solvent_layer", "fix"]
    assert "板の断面は" in status(html)
    comps = app.form["st2_comps"].splitlines()
    assert [re.search(r"label=(\S+)", c).group(1) for c in comps] == ["C3H4O3", "Li+", "PF6-"]
    assert "Li⁺" in html and "PF₆⁻" in html
    assert app.form["st2_vacuum"] == "0" and app.form["st2_tmode"] == "density"
    got = []
    for row in ("1", "2", "1"):
        html = act(app, base, "conc:2", st2_density="1.2", st2_conc="1", st2_conc_row=row)
        assert "mol/L" in html
        got.append([int(re.search(r"\* (\d+)", c).group(1)) for c in app.form["st2_comps"].splitlines()])
    assert [g[1] if i != 1 else g[2] for i, g in enumerate(got)] == [2, 3, 3] and got[-1] == [32, 3, 3], got
    assert app.form["charge"] == "0"
    html = act(app, base, "build")
    assert status(html) == "できました", status(html)
    st = app.spec.structure
    a = st.atoms.to_ase()
    assert len(a) == 72 + 32 * 10 + 3 + 3 * 7
    s = a.get_chemical_symbols()
    assert (s.count("Li"), s.count("P"), s.count("F"), s.count("Ni")) == (3, 3, 18, 48)
    assert a.cell.lengths()[0] == pytest.approx(11.0000) and a.cell.angles()[2] == pytest.approx(104.9, abs=0.1)
    assert len(st.fixed_atoms) == 24 and set(np.array(s)[st.fixed_atoms]) == {"Ni"}
    lines = log_lines(html)
    assert len(lines) == 4 and lines[0].startswith("0. 土台") and "416" in lines[2]
    assert app.spec is not None
    rec = Recipe.from_ref(st.source_ref)
    assert rec.base.source == "file" and len(rec.base.sha256) == 64
    spec_bytes = urllib.request.urlopen(base + "/spec.json").read()
    app.form.update(source="preset", recipe_steps="[]"); app.built = None
    html = _load(base, spec_bytes)
    assert app.form["source"] == "file" and ops(app) == ["supercell", "solvent_layer", "fix"]
    assert "spec.json の原子座標を使っています" in status(html)
    assert app.spec is not None, _error(html)
    assert app.spec.structure == st
    assert Recipe.from_ref(app.spec.structure.source_ref) == rec
    html = _post(base + "/preview", {**app.form, "st2_gap": "2.5"})
    assert app.spec is None and "作り直します" in unescape(html)


def test_new_bases(web):
    app, base = web
    html = _post(base + "/preview", {**app.form, "source": "2d"})
    assert app.spec is not None and app.spec.structure.source == "recipe" and len(app.spec.structure.atoms.symbols) == 8
    assert log_lines(html) and log_lines(html)[0].startswith("0. 土台")
    _post(base + "/preview", {**app.form, "source": "cluster", "cl_shells": "2"})
    assert len(app.spec.structure.atoms.symbols) == 13 and not app.spec.structure.periodic
    act(app, base, "add", recipe_add="box")
    act(app, base, "build", st1_padding="4")
    assert app.spec is not None and app.spec.structure.periodic
    ref = app.spec.structure.source_ref
    _load(base, app.spec.model_dump_json().encode())
    assert app.form["source"] == "cluster" and app.form["cl_shells"] == "2" and app.spec.structure.source_ref == ref
    pytest.importorskip("rdkit")
    act(app, base, "del:1")
    html = _post(base + "/preview", {**app.form, "source": "polymer", "pl_unit": "*CC*", "pl_n": "3"})
    assert app.spec is None and "作る" in unescape(html)
    html = act(app, base, "build")
    assert status(html) == "できました" and len(app.spec.structure.atoms.symbols) == 20  # C6H14


def test_atom_limit_stops_before_building(web):
    app, base = web
    app.form.update(source="bulk", bulk_cubic="on")
    act(app, base, "add", recipe_add="supercell")
    html = act(app, base, "build", st1_r0="50", st1_r1="50", st1_r2="50")
    assert "上限" in status(html) and "作る前に止めました" in status(html) and app.spec is None


def test_server_answers_while_building(web, monkeypatch):
    import adit.builder
    app, base = web
    started, release = threading.Event(), threading.Event()
    real = adit.builder.recipe_structure

    def slow(*a, **k):
        started.set(); release.wait(30)
        return real(*a, **k)

    monkeypatch.setattr(adit.builder, "recipe_structure", slow)
    app.form.update(source="bulk", bulk_cubic="on")
    act(app, base, "add", recipe_add="supercell")
    out: list[str] = []
    t = threading.Thread(target=lambda: out.append(act(app, base, "build")), daemon=True); t.start()
    assert started.wait(10)
    page = urllib.request.urlopen(base + "/analysis", timeout=5).read().decode()
    assert "<html" in page
    release.set(); t.join(30)
    assert out and status(out[0]) == "できました"


def test_english_and_labels(web):
    from adit import lang
    app, base = web
    lang.set_language("en")
    app.form.update(source="bulk")
    act(app, base, "add", recipe_add="slab")
    html = act(app, base, "add", recipe_add="solvent_layer")
    assert "Build steps" in html and "Miller indices (h k l)" in html and "Vacuum per side [Å]" in html
    assert "Count from concentration" in html and "Electrode–electrolyte interface" in html and ">Build<" in html
    assert list_rows(html)[0].startswith("1. slab") and ">組み立て手順<" not in html
    assert "2D material / nanotube" in html and "Nanoparticle" in html
    for op in ("supercell", "vacuum", "box", "adsorb", "remove", "substitute", "solvate", "fix"):
        html = act(app, base, "add", recipe_add=op)
    assert _fields_without_label(html) == []


def test_enter_key_goes_to_preview(web):
    app, base = web
    html = urllib.request.urlopen(base + "/").read().decode()
    main = html[html.index('<form method="post" action="/preview"'):]
    first = re.search(r"<button[^>]*type=\"submit\"[^>]*>", main).group(0)
    assert 'formaction="/preview"' in first


def test_unticked_step_is_skipped_without_deleting_it(web):
    app, base = web
    app.form.update(source="bulk", bulk_cubic="on")
    act(app, base, "add", recipe_add="supercell")
    html = act(app, base, "add", recipe_add="fix")
    assert 'name="st2_enabled" value="1" checked' in html and 'name="fixed" id="fixed" value="" disabled' in html
    posted = {k: v for k, v in app.form.items() if k != "st2_enabled"}      # an unticked box is not posted
    html = _post(base + "/recipe", {**posted, "st1_r0": "2", "st1_r1": "1", "st1_r2": "1", "recipe_action": "build"})
    steps = json.loads(app.form["recipe_steps"])
    assert steps[1]["op"] == "fix" and steps[1]["enabled"] is False and app.form["st2_enabled"] == ""
    assert 'name="st2_enabled" value="1" aria-label' in html and "(無効: 飛ばします)" in list_rows(html)[1]
    assert 'name="fixed" id="fixed" value="" disabled' not in html
    assert any("2. 固定: 無効なので飛ばしました" in x for x in log_lines(html))
    _post(base + "/preview", app.form)
    assert app.spec is not None and len(app.spec.structure.atoms.symbols) == 16 and app.spec.structure.fixed_atoms == []
    act(app, base, "build", st2_enabled="1")
    assert json.loads(app.form["recipe_steps"])[1]["enabled"] is True
    _post(base + "/preview", app.form)
    assert app.spec is not None and app.spec.structure.fixed_atoms
