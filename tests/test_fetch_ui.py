"""Fetching structures through the CLI, the web form and the desktop panel (network answered from tests/data/fetch)."""

from __future__ import annotations

import io
import json
import os
import threading
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

from adit.cli import main
from adit.config import default_config, save_config
from adit.spec import CalculationSpec
from tests.conftest import cfg_for, pbs_profile, water_spec
from tests.test_fetch import DATA, ROUTES, _Resp

REAL_URLOPEN = urllib.request.urlopen


@pytest.fixture
def offline(monkeypatch):
    """Answer the database URLs from fixtures; local (test server) URLs go through; everything else is 404."""
    seen: list[str] = []

    def fake_urlopen(req, *args, **kwargs):
        url = req.full_url if isinstance(req, urllib.request.Request) else req
        if url.startswith("http://127.0.0.1"):
            return REAL_URLOPEN(req, *args, **kwargs)
        seen.append(url)
        if url == "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/ethanol/cids/JSON":
            return _Resp(json.dumps({"IdentifierList": {"CID": [962, 963]}}).encode())
        name = ROUTES.get(url)
        if name is None:
            raise urllib.error.HTTPError(url, 404, "Not Found", {}, io.BytesIO(b""))
        return _Resp((DATA / name).read_bytes())

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    return seen


# ---- CLI ----
def test_cli_fetch_replaces_the_structure(offline, sk_root, tmp_path, capsys, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    cfg_path = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), cfg_path)
    spec_path = tmp_path / "spec.json"
    water_spec().save(spec_path)
    out = tmp_path / "calc"
    assert main([str(spec_path), str(out), "--config", str(cfg_path), "--fetch", "pubchem:water"]) == 0
    err = capsys.readouterr().err
    assert "構造を取得しました: PubChem 962 Water H2O" in err and "pubchem_962.sdf" in err
    assert (tmp_path / "home" / "adit_runs" / "fetched" / "pubchem_962.sdf").is_file()
    assert (out / "pubchem_962.sdf").is_file()
    prov = json.loads((out / "spec.json").read_text(encoding="utf-8"))["provenance"]
    rec = prov["fetched_structure"]
    assert rec["database"] == "pubchem" and rec["id"] == "962" and rec["sha256"] and rec["license_url"]
    assert any(f["name"] == "pubchem_962.sdf" and f["sha256"] == rec["sha256"] for f in prov["files"])
    readme = (out / "README.txt").read_text(encoding="utf-8")
    assert "== 構造の出どころ (データベースから取得) ==" in readme and "public domain" in readme
    saved = CalculationSpec.load(out / "spec.json")
    assert saved.structure.source == "file" and saved.structure.atoms.symbols == ["O", "H", "H"] and saved.structure.fetched == rec


def test_cli_fetch_candidates_and_failures_stop(offline, sk_root, tmp_path, capsys):
    cfg_path = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), cfg_path)
    spec_path = tmp_path / "spec.json"
    water_spec().save(spec_path)
    assert main([str(spec_path), str(tmp_path / "a"), "--config", str(cfg_path), "--fetch", "pubchem:ethanol"]) == 1
    err = capsys.readouterr().err
    assert "候補が 1 件あります" in err and "pubchem:962" in err and not (tmp_path / "a").exists()
    assert main([str(spec_path), str(tmp_path / "b"), "--config", str(cfg_path), "--fetch", "cod:9999999"]) == 1
    assert "取得できません: COD に ID 9999999 の項目がありません" in capsys.readouterr().err
    assert main([str(spec_path), str(tmp_path / "c"), "--config", str(cfg_path), "--fetch", "mp:mp-149"]) == 1
    assert "mp_api_key" in capsys.readouterr().err
    assert main([str(spec_path), "--validate", "--config", str(cfg_path), "--fetch", "nowhere:1"]) == 1
    assert "データベース:名前または ID" in capsys.readouterr().err


def test_cli_fetch_uses_the_configured_mp_key(offline, sk_root, tmp_path, capsys, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    cfg = cfg_for(sk_root); cfg.mp_api_key = "KEY"
    cfg_path = tmp_path / "cluster.toml"; save_config(cfg, cfg_path)
    assert "mp_api_key" in cfg_path.read_text(encoding="utf-8")
    spec_path = tmp_path / "spec.json"
    from adit.spec import XtbMethod
    water_spec(method=XtbMethod()).save(spec_path)
    assert main([str(spec_path), "--validate", "--config", str(cfg_path), "--fetch", "mp:mp-149"]) in (0, 1)
    err = capsys.readouterr().err
    assert "構造を取得しました: Materials Project mp-149 Si Fd-3m" in err


# ---- web ----
def test_web_form_roundtrip_with_fetch(offline, tmp_path):
    from adit.fetch import fetch
    from adit.structure import from_fetched
    from adit.web.forms import FormError, default_form, form_from_spec, spec_from_form

    res = fetch("cod:1000041"); res.fetched.save(tmp_path)
    st = from_fetched(res.fetched)
    spec = water_spec().model_copy(update={"structure": st})
    f = form_from_spec(spec)
    assert f["source"] == "fetch" and f["fetch_db"] == "cod" and f["fetch_query"] == "1000041" and f["fetch_file"] == st.source_ref
    assert json.loads(f["fetch_record"])["id"] == "1000041"
    back = spec_from_form({**default_form(), **f, "sk_set": "fake-1-0"})
    assert back.structure.source == "file" and back.structure.source_ref == st.source_ref and back.structure.fetched == st.fetched
    assert back.structure.atoms.symbols[:2] == ["Na", "Na"] and back.structure.periodic
    with pytest.raises(FormError, match="「取得」を押して構造を取得してください"):
        spec_from_form({**default_form(), "source": "fetch", "fetch_query": "1000041"})


def test_web_recipe_keeps_the_record(offline, tmp_path):
    from adit.fetch import fetch
    from adit.web.forms import FormError, default_form, spec_from_form
    from adit.web.recipe_form import Built, recipe_from_form, write_steps
    from adit.builder import recipe_structure
    from adit.builder.model import Supercell

    res = fetch("cod:1000041"); res.fetched.save(tmp_path)
    f = {**default_form(), "source": "fetch", "fetch_db": "cod", "fetch_query": "1000041", "fetch_file": res.fetched.record["file"],
         "fetch_record": json.dumps(res.fetched.record), "sk_set": "fake-1-0"}
    write_steps(f, [Supercell(repeat=(2, 1, 1))])
    rec = recipe_from_form(f)
    assert rec.base.source == "file" and rec.base.ref == res.fetched.record["file"] and rec.base.sha256

    class State:
        built = None

    state = State()
    with pytest.raises(FormError, match="「作る」を押す"):
        spec_from_form(f, state)
    built_st, _logs = recipe_structure(rec)
    state.built = Built(rec.to_ref(), built_st, [])
    st = spec_from_form(f, state).structure
    assert st.source == "recipe" and len(st.atoms.symbols) == 16 and st.fetched["id"] == "1000041"


@pytest.fixture
def web(sk_root, tmp_path):
    from adit.web.server import WebApp, serve

    cfg = default_config(sk_root=str(sk_root)); cfg.profiles["cluster"] = pbs_profile()
    app = WebApp(cfg, tmp_path / "cluster.toml")
    app.fetch_dir = tmp_path / "fetched"
    httpd = serve(app, "127.0.0.1", 0)
    t = threading.Thread(target=httpd.serve_forever, daemon=True); t.start()
    yield app, f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown(); httpd.server_close()


def _post(url: str, fields: dict) -> str:
    return REAL_URLOPEN(url, data=urllib.parse.urlencode(fields).encode()).read().decode()


def test_web_fetch_button_then_candidates(offline, web, tmp_path):
    app, base = web
    html = REAL_URLOPEN(base + "/").read().decode()
    assert 'value="fetch"' in html and "データベースから取得" in html and 'formaction="/fetch_structure"' in html
    assert not offline, "opening the page must not contact any database"
    html = _post(base + "/fetch_structure", {"source": "fetch", "fetch_db": "cod", "fetch_query": "1000041", "sk_set": "fake-1-0"})
    assert "構造を取得しました" in html and "CC0 1.0" in html and app.spec is not None
    assert app.spec.structure.source == "file" and app.spec.structure.fetched["id"] == "1000041" and app.spec.structure.periodic
    assert (tmp_path / "fetched" / "cod_1000041.cif").is_file() and app.form["source"] == "fetch"
    assert 'name="fetch_record"' in html and "Crystallography Open Database (COD) 1000041" in html
    html = _post(base + "/preview", {**app.form, "charge": "0"})
    assert app.spec.structure.fetched["id"] == "1000041" and "生成できます" in html
    html = _post(base + "/fetch_structure", {"source": "fetch", "fetch_db": "pubchem", "fetch_query": "ethanol", "sk_set": "fake-1-0"})
    assert "候補が 1 件あります" in html and 'name="fetch_pick"' in html and 'value="pubchem:962"' in html
    html = _post(base + "/fetch_structure", {"source": "fetch", "fetch_db": "pubchem", "fetch_query": "ethanol", "fetch_pick": "pubchem:962",
                                             "sk_set": "fake-1-0"})
    assert "構造を取得しました: PubChem 962 Water H2O" in html and app.spec.structure.fetched["id"] == "962" and not app.fetch_candidates
    html = _post(base + "/fetch_structure", {"source": "fetch", "fetch_db": "cod", "fetch_query": "424242", "sk_set": "fake-1-0"})
    assert "取得できません: COD に ID 424242 の項目がありません" in html
    html = _post(base + "/fetch_structure", {"source": "fetch", "fetch_db": "cod", "fetch_query": "", "sk_set": "fake-1-0"})
    assert "名前か ID を入力してください" in html


# ---- desktop ----
pytest.importorskip("PySide6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def qapp():
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


def test_panel_fetch_flow(offline, qapp, tmp_path, monkeypatch):
    from adit.gui.panels.structure_panel import SOURCES, StructurePanel

    assert list(SOURCES)[-1] == "fetch" and list(SOURCES)[:3] == ["preset", "smiles", "file"]
    p = StructurePanel(); p.fetch_in_thread = False; p.fetch_dir = tmp_path
    assert [p.source.itemData(i) for i in range(p.source.count())][-1] == "fetch"
    p.set_source("fetch")
    assert p.structure() is None and "「取得」を押してください" in p.error() and not offline
    p.fetch_db.setCurrentIndex(p.fetch_db.findData("cod")); p.fetch_query.setText("1000041")
    assert "1000041" in p.fetch_query.placeholderText()
    p._fetch(); qapp.processEvents()
    s = p.structure()
    assert s is not None and s.source == "file" and s.source_ref == str(tmp_path / "cod_1000041.cif") and s.fetched["id"] == "1000041"
    assert "CC0 1.0" in p.fetch_note.text() and p.fetch_go.isEnabled()
    assert p.current_recipe().base.source == "file" and p.current_recipe().base.sha256
    p.charge.setValue(1)
    assert p.structure().fetched["id"] == "1000041" and p.structure().charge == 1
    # restoring a saved spec brings the fetch rows back
    q = StructurePanel(); q.set_structure(p.structure())
    assert q.current_source() == "fetch" and q.fetch_db.currentData() == "cod" and q.fetch_query.text() == "1000041"
    assert q.structure().fetched["id"] == "1000041" and "1000041" in q.fetch_note.text()
    q.set_source("preset")
    assert q.structure().fetched is None
    # several candidates: the user picks one in a dialog
    from PySide6.QtWidgets import QInputDialog
    asked: list = []

    def pick(parent, title, text, items, current=0, editable=True):
        asked.append(list(items)); return items[0], True

    monkeypatch.setattr(QInputDialog, "getItem", staticmethod(pick))
    p.fetch_db.setCurrentIndex(p.fetch_db.findData("pubchem")); p.fetch_query.setText("ethanol")
    p._fetch(); qapp.processEvents(); qapp.processEvents()
    assert asked == [["CID 962: Water H2O"]] and p.structure().fetched["id"] == "962" and p.structure().atoms.symbols == ["O", "H", "H"]
    # a failure leaves the previous structure and shows the reason
    p.fetch_db.setCurrentIndex(p.fetch_db.findData("cod")); p.fetch_query.setText("424242")
    p._fetch(); qapp.processEvents()
    assert p.fetch_note.text().startswith("取得できません: COD に ID 424242") and p.structure().fetched["id"] == "962"
