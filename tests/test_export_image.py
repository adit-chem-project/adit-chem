import json
import shutil
import threading
import urllib.parse
import urllib.request
from urllib.error import HTTPError
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from ase.build import bulk, molecule

from adit.export_image import scene_to_svg, sketch_to_svg
from adit.sketch import Sketch, ring_atoms, sketch_to_dict
from adit.web.structure3d import scene_from_atoms

REPO = Path(__file__).resolve().parents[1]


def _parse(svg: str):
    root = ET.fromstring(svg)
    assert root.tag.endswith("svg")
    return root


def test_a_drawn_structure_becomes_a_valid_svg():
    sketch = Sketch()
    sketch.atoms, sketch.bonds = ring_atoms("benzene", 100.0, 100.0)
    svg = sketch_to_svg(sketch)
    root = _parse(svg)
    lines = [e for e in root.iter() if e.tag.endswith("line")]
    assert len(lines) == 9
    assert root.get("width") and root.get("viewBox")


def test_hetero_atoms_get_their_symbol_and_color():
    pytest.importorskip("rdkit")
    sketch = Sketch.from_smiles("CCO")
    svg = sketch_to_svg(sketch)
    root = _parse(svg)
    texts = [e.text for e in root.iter() if e.tag.endswith("text")]
    assert texts == ["O"]
    assert "#D02020" in svg or "#d02020" in svg.lower()


def test_charges_are_written_next_to_the_symbol():
    pytest.importorskip("rdkit")
    sketch = Sketch.from_smiles("[NH4+]")
    svg = sketch_to_svg(sketch)
    assert "N+" in svg


def test_an_empty_sketch_is_refused():
    with pytest.raises(ValueError, match="描かれていません"):
        sketch_to_svg(Sketch())


def test_the_3d_view_becomes_a_valid_svg():
    scene = scene_from_atoms(bulk("Si", "diamond", a=5.43, cubic=True))
    svg = scene_to_svg(scene)
    root = _parse(svg)
    circles = [e for e in root.iter() if e.tag.endswith("circle")]
    lines = [e for e in root.iter() if e.tag.endswith("line")]
    assert len(circles) == 8
    assert len(lines) >= 12
    assert svg.index("<circle") < svg.rindex("<circle")


def _circle_x_by_color(svg: str) -> dict[str, float]:
    import re

    return {m.group(2): float(m.group(1)) for m in re.finditer(r'<circle cx="([\d.]+)" cy="[\d.]+" r="[\d.]+" fill="(#[0-9a-fA-F]+)"', svg)}


def test_the_3d_svg_follows_the_given_rotation():
    scene = scene_from_atoms(molecule("CO"))                 # C and O differ only along one axis
    carbon, oxygen = scene.colors
    axis = next(k for k in range(3) if abs(scene.positions[0][k] - scene.positions[1][k]) > 0.5)
    identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    # a quarter turn that brings that axis to the screen x axis, and the same turn mirrored
    to_x = [[0.0] * 3 for _ in range(3)]; to_x[0][axis] = 1.0; to_x[1][(axis + 1) % 3] = 1.0; to_x[2][(axis + 2) % 3] = 1.0
    flipped = [[-v for v in to_x[0]], to_x[1], [-v for v in to_x[2]]]
    plain = _circle_x_by_color(scene_to_svg(scene, rotation=to_x))
    mirror = _circle_x_by_color(scene_to_svg(scene, rotation=flipped))
    assert (plain[carbon] < plain[oxygen]) == (scene.positions[0][axis] < scene.positions[1][axis])
    assert (mirror[carbon] < mirror[oxygen]) != (plain[carbon] < plain[oxygen])
    assert scene_to_svg(scene, rotation=identity) != scene_to_svg(scene)   # the default view is a tilted one


def test_the_molecule_svg_has_no_cell_lines():
    svg = scene_to_svg(scene_from_atoms(molecule("H2O")))
    root = _parse(svg)
    assert len([e for e in root.iter() if e.tag.endswith("circle")]) == 3
    assert len([e for e in root.iter() if e.tag.endswith("line")]) == 2 * 3


def test_the_mol_file_can_be_read_back():
    pytest.importorskip("rdkit")
    from rdkit import Chem

    from adit.export_image import sketch_to_mol

    mol = sketch_to_mol(Sketch.from_smiles("CC(=O)O"))
    back = Chem.MolFromMolBlock(mol)
    assert back is not None and Chem.MolToSmiles(back) == "CC(=O)O"


def test_the_analysis_figures_can_be_written_as_vectors(tmp_path):
    import subprocess
    import sys

    out = tmp_path / "o"
    r = subprocess.run([sys.executable, "-m", "adit.analysis.cli",
                        str(REPO / "examples" / "dftb_md_water_generated"), "-o", str(out),
                        "--figure-format", "svg,pdf"], capture_output=True, text=True, cwd=REPO, timeout=600)
    assert r.returncode == 0, r.stderr
    assert (out / "energy.png").is_file() and (out / "energy.svg").is_file() and (out / "energy.pdf").is_file()
    assert _parse((out / "energy.svg").read_text(encoding="utf-8")) is not None
    bad = subprocess.run([sys.executable, "-m", "adit.analysis.cli",
                          str(REPO / "examples" / "dftb_md_water_generated"), "-o", str(tmp_path / "x"),
                          "--figure-format", "docx"], capture_output=True, text=True, cwd=REPO, timeout=600)
    assert "図の形式は" in bad.stdout + bad.stderr


def test_the_browser_can_download_the_svg_and_the_mol(tmp_path):
    pytest.importorskip("rdkit")
    from adit.config import default_config
    from adit.web.server import WebApp, serve

    app = WebApp(default_config(sk_root=str(REPO / "slakos")), tmp_path / "cluster.toml")
    srv = serve(app, port=0)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_port}"
    try:
        sketch = Sketch()
        sketch.atoms, sketch.bonds = ring_atoms("benzene", 100.0, 100.0)
        payload = urllib.parse.urlencode({"action": "svg", "sketch": json.dumps(sketch_to_dict(sketch))}).encode()
        with urllib.request.urlopen(urllib.request.Request(base + "/draw", data=payload), timeout=60) as r:
            assert r.headers.get("Content-Type") == "image/svg+xml"
            assert "attachment" in r.headers.get("Content-Disposition", "")
            assert _parse(r.read().decode()) is not None
        html = urllib.request.urlopen(base + "/draw", timeout=30).read().decode()
        assert 'value="svg"' in html and 'value="mol"' in html
    finally:
        srv.shutdown()


def test_the_desktop_offers_copy_and_save():
    pytest.importorskip("PySide6")
    source = (REPO / "src" / "adit" / "gui").rglob("*.py")
    text = "\n".join(p.read_text(encoding="utf-8") for p in source)
    assert "copy_widget" in text and "save_dialog" in text
    panel = (REPO / "src" / "adit" / "gui" / "panels" / "structure_view_panel.py").read_text(encoding="utf-8")
    sketcher = (REPO / "src" / "adit" / "gui" / "sketcher.py").read_text(encoding="utf-8")
    assert "btn_copy" in panel and "btn_save" in panel
    assert "btn_copy" in sketcher and "btn_save" in sketcher


def test_the_desktop_can_copy_and_save_each_figure(tmp_path):
    pytest.importorskip("PySide6")
    from PySide6.QtWidgets import QApplication

    from adit.gui.panels.analysis_panel import AnalysisPanel

    QApplication.instance() or QApplication([])
    panel = AnalysisPanel()
    run_dir = tmp_path / "run"
    shutil.copytree(REPO / "examples" / "dftb_md_water_generated", run_dir)   # never analyse the example in place
    panel.run_dir.setText(str(run_dir))
    panel.figure_format.setText("svg")
    panel.options()
    out = tmp_path / "figs"
    panel.out_dir_override = str(out) if hasattr(panel, "out_dir_override") else None
    res = panel.run()
    assert res is not None and res.figures
    path = next(iter(res.figures.values()))
    others = panel.figure_siblings(path)
    assert "svg" in others, "--figure-format svg を指定したのに SVG が無い"
    panel._copy_figure(path)
    assert "コピー" in panel.summary.text() or "clipboard" in panel.summary.text()


def test_the_browser_lists_a_download_link_for_every_figure(tmp_path):
    from adit.config import default_config
    from adit.web.server import WebApp, serve

    app = WebApp(default_config(sk_root=str(REPO / "slakos")), tmp_path / "cluster.toml")
    srv = serve(app, port=0)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{srv.server_port}"
    try:
        run_dir = tmp_path / "run"
        shutil.copytree(REPO / "examples" / "dftb_md_water_generated", run_dir)   # never analyse the example in place
        form = {"run_dir": str(run_dir), "figure_format": "svg"}
        req = urllib.request.Request(base + "/analysis", data=urllib.parse.urlencode(form).encode())
        with urllib.request.urlopen(req, timeout=900) as r:
            html = r.read().decode()
        import re

        links = re.findall(r'href="(/file\?path=[^"]+download=1)"', html)
        assert links, "図を保存するリンクが無い"
        assert any(".svg" in link for link in links), "SVG のリンクが無い"
        target = next(link for link in links if ".svg" in link).replace("&amp;", "&")
        with urllib.request.urlopen(base + target, timeout=60) as r:
            assert r.headers.get("Content-Type") == "image/svg+xml"
            assert "attachment" in r.headers.get("Content-Disposition", "")
            assert _parse(r.read().decode()) is not None
        with pytest.raises(HTTPError):
            urllib.request.urlopen(base + "/file?path=" + urllib.parse.quote("/etc/passwd"), timeout=30)
    finally:
        srv.shutdown()
