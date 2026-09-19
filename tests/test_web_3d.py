import json
import re
from pathlib import Path

import pytest
from ase.build import bulk, molecule

from adit.web.structure3d import BOND_FACTOR, MAX_BOND_ATOMS, scene_from_atoms, summary_line

JS = Path(__file__).resolve().parents[1] / "src" / "adit" / "web" / "static" / "viewer3d.js"


def test_water_scene_matches_the_desktop_rules():
    scene = scene_from_atoms(molecule("H2O"))
    assert scene.n_atoms == 3
    assert scene.bonds == [[0, 1], [0, 2]]
    assert scene.colors[0] == "#ff0d0d" and scene.colors[1] == "#ffffff"
    assert scene.radii[0] == pytest.approx(0.66, abs=0.01)
    assert scene.cell_lines == []
    for axis in range(3):
        assert abs(sum(p[axis] for p in scene.positions)) < 1e-3


def test_periodic_cell_gets_twelve_edges():
    scene = scene_from_atoms(bulk("Si", "diamond", a=5.43, cubic=True))
    assert len(scene.cell_lines) == 12
    assert scene.n_atoms == 8
    assert scene.scale > 4.0


def test_bond_rule_is_the_same_factor_as_the_desktop():
    from ase import Atoms
    from ase.data import covalent_radii

    limit = (covalent_radii[1] * 2) * BOND_FACTOR
    close = scene_from_atoms(Atoms("H2", positions=[[0, 0, 0], [limit * 0.9, 0, 0]]))
    far = scene_from_atoms(Atoms("H2", positions=[[0, 0, 0], [limit * 1.1, 0, 0]]))
    assert close.bonds == [[0, 1]] and far.bonds == []


def test_too_many_atoms_says_so_instead_of_freezing_the_browser():
    from ase import Atoms

    n = MAX_BOND_ATOMS + 1
    atoms = Atoms("H" * n, positions=[[i * 3.0, 0, 0] for i in range(n)])
    scene = scene_from_atoms(atoms)
    assert scene.truncated_bonds and scene.bonds == []
    assert "結合" in summary_line(scene, False) or "bonds" in summary_line(scene, False)


def test_json_is_safe_to_put_inside_a_script_tag():
    scene = scene_from_atoms(molecule("H2O"))
    text = scene.to_json()
    assert "<" not in text
    assert json.loads(text)["n_atoms"] == 3


def test_the_viewer_loads_nothing_from_the_outside():
    js = JS.read_text(encoding="utf-8")
    assert "http://" not in js and "https://" not in js
    assert "import" not in js and "require(" not in js
    assert "getContext" in js and "ADIT_SCENE" in js


def test_the_page_shows_the_viewer_and_says_what_happens_without_javascript(tmp_path):
    from adit.config import default_config
    from adit.web.server import WebApp

    repo = Path(__file__).resolve().parents[1]
    app = WebApp(default_config(sk_root=str(repo / "slakos")), tmp_path / "cluster.toml")
    status, err = app.preview({"source": "preset", "preset": "H2O", "code": "dftbplus",
                               "task_type": "single_point", "sk_set": "mio-1-1",
                               "output_dir": str(tmp_path / "out")})
    if err:
        pytest.skip(f"この環境では構造を組み立てられない: {err}")
    scene_json, note, periodic = app.scene()
    assert json.loads(scene_json)["n_atoms"] == 3
    assert "原子 3 個" in note or "3 atoms" in note
    assert periodic is False
    html = (repo / "src" / "adit" / "web" / "templates" / "index.html").read_text(encoding="utf-8")
    assert 'id="adit3d"' in html and "noscript" in html
    assert "JavaScript" in html
    assert "scene_json | safe" in html


def test_no_scene_before_a_structure_is_built(tmp_path):
    from adit.config import default_config
    from adit.web.server import WebApp

    app = WebApp(default_config(sk_root=""), tmp_path / "cluster.toml")
    assert app.scene() == ("", "", False)


# ---- measuring in the browser and the trajectory player (2026-09-19) ----
PLAYBACK_JS = JS.with_name("playback.js")


def test_scene_carries_the_cell_for_minimum_image_measuring():
    assert scene_from_atoms(molecule("H2O")).cell is None
    si = scene_from_atoms(bulk("Si"))
    assert si.cell is not None and len(si.cell) == 3 and len(si.cell[0]) == 3
    assert "cell" in json.loads(si.to_json())


def test_the_scripts_have_the_measuring_and_playing_hooks_and_load_nothing():
    js = JS.read_text(encoding="utf-8")
    for needle in ("ADIT_VIEWER", "shiftKey", "adit3d-send-fixed", "adit3d-measure", "function dihedral", "function mic"):
        assert needle in js
    pb = PLAYBACK_JS.read_text(encoding="utf-8")
    assert "http://" not in pb and "https://" not in pb and "import" not in pb and "require(" not in pb
    assert "ADIT_FRAMES_URL" in pb and "adit-play-slider" in pb
    html = (Path(__file__).resolve().parents[1] / "src" / "adit" / "web" / "templates" / "index.html").read_text(encoding="utf-8")
    assert 'id="adit3d-send-fixed"' in html and 'id="adit3d-clear"' in html and "Shift" in html


def test_frames_payload_thins_to_the_browser_limit(tmp_path):
    from adit.web.frames import frames_json, frames_payload

    repo = Path(__file__).resolve().parents[1]
    run = repo / "examples" / "dftb_md_water_generated"
    d = frames_payload(run)
    assert d["n_total"] == 21 and d["n_frames"] == 21 and d["stride"] == 1
    assert d["scene"]["n_atoms"] == 3 and len(d["positions"]) == 21 and len(d["positions"][0]) == 3
    assert d["times_fs"][5] == pytest.approx(25.0) and len(d["energies_ev"]) == 21 and len(d["temperatures_k"]) == 21
    assert d["is_md"] and d["source"] == "geo_end.xyz"
    small = frames_payload(run, max_frames=5)
    assert small["stride"] == 5 and small["n_frames"] == 5 and small["frame_index"] == [0, 5, 10, 15, 20]
    assert any("5 フレームに 1 回" in n for n in small["notes"])
    assert frames_payload(run, stride=2, skip=1)["frame_index"] == list(range(1, 21, 2))
    text = frames_json(run, max_frames=3)
    assert "<" not in text and json.loads(text)["n_frames"] == 3
    one = frames_payload(repo / "examples" / "water_generated")
    assert one["n_frames"] == 0 and any("1 つ" in n for n in one["notes"])


def test_analysis_page_serves_the_frames_only_for_the_analyzed_directory(tmp_path):
    import shutil
    import threading
    import urllib.error
    import urllib.parse
    import urllib.request

    from adit.config import default_config
    from adit.web.server import WebApp, serve

    repo = Path(__file__).resolve().parents[1]
    run_dir = tmp_path / "md"
    shutil.copytree(repo / "examples" / "dftb_md_water_generated", run_dir)
    app = WebApp(default_config(sk_root=""), tmp_path / "cluster.toml")
    httpd = serve(app, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_port}"
    try:
        with pytest.raises(urllib.error.HTTPError):                      # nothing analyzed yet
            urllib.request.urlopen(base + "/frames.json?dir=" + urllib.parse.quote(str(run_dir)))
        html = urllib.request.urlopen(base + "/analysis", data=urllib.parse.urlencode({"run_dir": str(run_dir), "stride": "2"}).encode()).read().decode()
        assert 'id="adit-play"' in html and "ADIT_FRAMES_URL" in html and "ADIT_VIEWER" in html and "adit-play-slider" in html
        m = re.search(r'window.ADIT_FRAMES_URL = "([^"]+)"', html)
        assert m
        data = json.loads(urllib.request.urlopen(base + m.group(1)).read().decode())
        assert data["n_frames"] == 11 and data["stride"] == 2                # frames 0, 2, ..., 20: the stride of the analysis is applied
        with pytest.raises(urllib.error.HTTPError):                      # any other directory is refused
            urllib.request.urlopen(base + "/frames.json?dir=" + urllib.parse.quote(str(repo / "examples" / "water_generated")))
    finally:
        httpd.shutdown(); httpd.server_close()
