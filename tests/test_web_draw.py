import json
import threading
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

from adit.sketch import Sketch, ring_atoms, sketch_from_dict, sketch_to_dict

pytestmark = pytest.mark.skipif(__import__("importlib.util", fromlist=["util"]).find_spec("rdkit") is None,
                                reason="needs RDKit")

REPO = Path(__file__).resolve().parents[1]
JS = REPO / "src" / "adit" / "web" / "static" / "sketch.js"


@pytest.fixture
def web(tmp_path):
    from adit.config import default_config
    from adit.web.server import WebApp, serve

    app = WebApp(default_config(sk_root=str(REPO / "slakos")), tmp_path / "cluster.toml")
    srv = serve(app, port=0)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_port}"

    def get(path="/draw"):
        with urllib.request.urlopen(base + path, timeout=30) as r:
            return r.read().decode()

    def post(data):
        req = urllib.request.Request(base + "/draw", data=urllib.parse.urlencode(data).encode())
        with urllib.request.urlopen(req, timeout=60) as r:
            return r.geturl(), r.read().decode()

    try:
        yield app, get, post
    finally:
        srv.shutdown()


def _benzene_json() -> str:
    sk = Sketch()
    sk.atoms, sk.bonds = ring_atoms("benzene", 300.0, 200.0)
    return json.dumps(sketch_to_dict(sk))


def test_the_page_has_the_same_tools_as_the_desktop(web):
    from adit.sketch import ELEMENTS, TEMPLATES

    _, get, _ = web
    html = get()
    assert 'id="sketchpad"' in html
    for element in ELEMENTS:
        assert f'data-elem="{element}"' in html
    for ring in TEMPLATES:
        assert f'data-tool="ring:{ring}"' in html
    for tool in ("1", "2", "3", "eraser"):
        assert f'data-tool="{tool}"' in html
    assert 'data-charge="1"' in html and 'data-charge="-1"' in html
    assert "noscript" in html and "JavaScript" in html


def test_a_drawn_ring_becomes_the_same_smiles_as_in_the_desktop(web):
    _, _, post = web
    _, html = post({"action": "check", "sketch": _benzene_json()})
    assert "c1ccccc1" in html
    sk = sketch_from_dict(json.loads(_benzene_json()))
    assert sk.to_smiles() == "c1ccccc1"
    assert sk.formula() == "C6H6"


def test_using_the_molecule_fills_the_structure_form(web):
    app, _, post = web
    url, _ = post({"action": "use", "sketch": _benzene_json()})
    assert url.endswith("/")
    assert app.form["source"] == "smiles" and app.form["smiles"] == "c1ccccc1"


def test_loading_a_smiles_puts_atoms_on_the_canvas(web):
    _, _, post = web
    _, html = post({"action": "load", "from_smiles": "CCO"})
    assert "3 原子" in html or "3 atoms" in html
    data = json.loads(html.split('id="sketch_data" value="', 1)[1].split('"', 1)[0].replace("&#34;", '"'))
    assert len(data["atoms"]) == 3 and len(data["bonds"]) == 2


def test_unreadable_smiles_says_why(web):
    _, _, post = web
    _, html = post({"action": "load", "from_smiles": "これは SMILES ではない"})
    assert "解釈できません" in html or "cannot parse" in html


def test_a_drawing_rdkit_refuses_is_reported_not_silently_dropped(web):
    _, _, post = web
    bad = json.dumps({"atoms": [{"x": 0, "y": 0, "elem": "C"}, {"x": 40, "y": 0, "elem": "C"}],
                      "bonds": [{"a": 0, "b": 1, "order": 3}, {"a": 0, "b": 1, "order": 3}]})
    _, html = post({"action": "check", "sketch": bad})
    assert "受け付けません" in html or "does not accept" in html


def test_a_broken_payload_is_refused_with_a_reason(web):
    _, _, post = web
    _, html = post({"action": "check", "sketch": json.dumps({"atoms": [{"x": 0, "y": 0, "elem": "C"}],
                                                             "bonds": [{"a": 0, "b": 5, "order": 1}]})})
    assert "範囲の外" in html or "out of range" in html
    _, empty = post({"action": "check", "sketch": "{}"})
    assert "まだ何も描かれていません" in empty or "nothing has been drawn" in empty


def test_the_editor_loads_nothing_from_the_outside():
    js = JS.read_text(encoding="utf-8")
    assert "http://" not in js and "https://" not in js
    assert "import " not in js and "require(" not in js
    assert "sketchpad" in js and "ring:" in js


def test_the_desktop_and_the_web_share_the_same_chemistry():
    desktop = (REPO / "src" / "adit" / "gui" / "sketcher.py").read_text(encoding="utf-8")
    assert "from adit.sketch import" in desktop
    assert "Chem.MolToSmiles" not in desktop
    server = (REPO / "src" / "adit" / "web" / "server.py").read_text(encoding="utf-8")
    assert "from adit.sketch import" in server
