"""Isosurfaces of volumetric grids: marching tetrahedra, the readers, the desktop box and the browser route."""

import json
import os
import re
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.io.cube import write_cube

from adit.analysis.isosurface import (CORNERS, MAX_TRIANGLES, TETRAHEDRA, IsosurfaceError, IsosurfaceTooLarge, check_budget,
                                      grid_stats, isosurface, isosurfaces, payload)
from adit.analysis.volumetric import Grid, read_grid

REPO = Path(__file__).resolve().parent.parent
BOX = 8.0


def _gaussian_grid(n=40, sigma=1.0, kind="s", shear=None):
    """A Gaussian (s) or x * Gaussian (p) sampled on n^3 points around the origin; ASE cube convention (cell = n * step)."""
    step = np.eye(3) * BOX / (n - 1) if shear is None else np.asarray(shear, dtype=float) * BOX / (n - 1)
    origin = np.array([-BOX / 2] * 3)
    idx = np.stack(np.meshgrid(np.arange(n), np.arange(n), np.arange(n), indexing="ij"), axis=-1).reshape(-1, 3).astype(float)
    xyz = origin + idx @ step
    r2 = (xyz ** 2).sum(axis=1)
    vals = np.exp(-r2 / (2 * sigma ** 2))
    if kind == "p":
        vals = xyz[:, 0] * vals
    atoms = Atoms("H", positions=[[0.0, 0.0, 0.0]], cell=step * n, pbc=True)
    return Grid(values=vals.reshape(n, n, n), atoms=atoms, kind="unknown", unit="", source=Path("test.cube"), origin=origin)


def test_the_six_tetrahedra_fill_the_cube_exactly_once():
    vols = []
    for tet in TETRAHEDRA:
        p = CORNERS[list(tet)].astype(float)
        vols.append(abs(np.linalg.det(p[1:] - p[0])) / 6)
    assert vols == pytest.approx([1 / 6] * 6)
    rng = np.random.default_rng(0)
    pts = rng.random((3000, 3))
    cover = np.zeros(len(pts))
    for tet in TETRAHEDRA:
        p = CORNERS[list(tet)].astype(float)
        lam = np.linalg.solve((p[1:] - p[0]).T, (pts - p[0]).T).T
        cover += (lam >= -1e-9).all(axis=1) & (lam.sum(axis=1) <= 1 + 1e-9)
    assert cover.min() >= 1 and cover.max() == 1


def test_sphere_of_a_gaussian_has_the_expected_radius_and_outward_normals():
    g = _gaussian_grid()
    level = 0.5
    m = isosurface(g, level)
    r = np.linalg.norm(m.vertices, axis=1)
    expected = np.sqrt(-2 * np.log(level))                  # exp(-r^2/2) = level
    assert m.n_triangles > 1000
    assert r.mean() == pytest.approx(expected, abs=0.01) and r.std() < 0.01
    assert (np.einsum("ij,ij->i", m.normals, m.vertices) > 0).all()             # away from the enclosed high values
    p = m.vertices[m.faces]
    fn = np.cross(p[:, 1] - p[:, 0], p[:, 2] - p[:, 0])
    assert (np.einsum("ij,ij->i", fn, m.normals[m.faces].mean(axis=1)) > 0).all()  # winding agrees with the normals
    assert np.linalg.norm(m.normals, axis=1) == pytest.approx(np.ones(len(m.normals)))
    assert m.stride == 1 and m.notes == []


def test_vertices_are_placed_right_on_a_sheared_grid():
    g = _gaussian_grid(shear=[[1.0, 0, 0], [0.5, 1.0, 0], [0, 0.3, 1.0]])
    m = isosurface(g, 0.5)
    r = np.linalg.norm(m.vertices, axis=1)
    assert r.mean() == pytest.approx(np.sqrt(-2 * np.log(0.5)), abs=0.01) and r.std() < 0.01


def test_positive_and_negative_lobes_of_a_p_like_field():
    g = _gaussian_grid(kind="p")
    both = isosurfaces(g, 0.2)
    assert [m.level for m in both] == [0.2, -0.2]
    assert both[0].vertices[:, 0].min() > 0 and both[1].vertices[:, 0].max() < 0
    assert both[0].n_triangles == both[1].n_triangles                         # the field is antisymmetric in x
    assert (np.einsum("ij,ij->i", both[1].normals, both[1].vertices - both[1].vertices.mean(axis=0)) > 0).mean() > 0.95
    only = isosurfaces(_gaussian_grid(), 0.5)                                  # no negative values: one surface
    assert len(only) == 1
    with pytest.raises(IsosurfaceError):
        isosurfaces(g, -0.2)
    with pytest.raises(IsosurfaceError):
        isosurface(g, float("nan"))


def test_periodic_grids_are_wrapped_so_the_surface_closes_at_the_cell_faces():
    n = 24
    vals = np.cos(2 * np.pi * np.arange(n) / n)[:, None, None] * np.ones((n, n, n))
    g = Grid(values=vals, atoms=Atoms("H", cell=np.eye(3) * 10, pbc=True), kind="density", unit="", source=Path("CHGCAR"), periodic=True)
    m = isosurface(g, 0.0)
    assert sorted(set(np.round(m.vertices[:, 0], 2))) == [2.5, 7.5]          # two sheets where cos crosses zero
    assert m.vertices[:, 1].max() == pytest.approx(10.0) and m.vertices[:, 2].max() == pytest.approx(10.0)   # up to the far face
    assert any("周期" in n_ or "periodic" in n_ for n_ in m.notes)
    open_ = isosurface(g, 0.0, periodic=False)
    assert open_.vertices[:, 1].max() == pytest.approx(10.0 * (n - 1) / n)


def test_the_budget_stops_before_extracting_and_suggests_a_stride():
    g = _gaussian_grid()
    info = check_budget(g.values, 0.5, 1, False, None)
    assert info["active_cubes"] > 0 and info["estimated_triangles"] > 0
    with pytest.raises(IsosurfaceTooLarge) as ex:
        isosurface(g, 0.5, max_triangles=500)
    assert ex.value.suggested_stride >= 2 and "500" in str(ex.value)
    m = isosurface(g, 0.5, stride=ex.value.suggested_stride, max_triangles=500)
    assert 0 < m.n_triangles <= 500 * 1.5                                     # the estimate is a ratio, so allow slack
    assert np.linalg.norm(m.vertices, axis=1).mean() == pytest.approx(np.sqrt(-2 * np.log(0.5)), abs=0.03)
    assert any("間引" in n_ or "every" in n_ for n_ in m.notes)
    assert MAX_TRIANGLES >= 10_000
    # the estimate is close to the real count (measured ratio of 6 triangles per active cube)
    real = isosurface(g, 0.5, max_triangles=None).n_triangles
    assert 0.7 < real / info["estimated_triangles"] < 1.4


def test_stats_show_the_range_and_the_reference_but_no_default():
    st = grid_stats(_gaussian_grid(kind="p"))
    assert st["has_negative"] and st["min"] == pytest.approx(-st["max"])
    assert st["tenth_of_abs_max"] == pytest.approx(0.1 * st["abs_max"])
    assert st["shape"] == (40, 40, 40) and st["n_points"] == 64000
    assert "default" not in st and "level" not in st


def test_payload_is_flat_and_rounded():
    m = isosurface(_gaussian_grid(n=20), 0.5)
    out = payload([m])
    assert out[0]["n_triangles"] == m.n_triangles and len(out[0]["tri"]) == 9 * m.n_triangles and len(out[0]["normal"]) == 3 * m.n_triangles
    assert out[0]["level"] == 0.5


# ---- readers ----
def _write_cube(path: Path, grid: Grid, comment: str = "adit test") -> Path:
    with open(path, "w", encoding="utf-8") as f:
        write_cube(f, grid.atoms, data=grid.values, origin=grid.origin, comment=comment)
    return path


def test_cube_reader_keeps_the_origin_and_marks_pp_x_cubes_periodic(tmp_path):
    g = _gaussian_grid(n=13)                                     # odd: a grid point sits on the origin, where the value is 1
    p = _write_cube(tmp_path / "orbital.cube", g)
    got = read_grid(p)
    assert got.origin == pytest.approx(g.origin, abs=1e-4) and got.periodic is False
    assert got.values.shape == (13, 13, 13) and got.values.max() == pytest.approx(1.0, abs=1e-3)
    m = isosurface(got, 0.5)
    assert np.linalg.norm(m.vertices, axis=1).mean() == pytest.approx(np.sqrt(-2 * np.log(0.5)), abs=0.05)
    q = _write_cube(tmp_path / "rho.cube", g, comment="Cubefile created from PWScf calculation")
    assert read_grid(q).periodic is True                         # pp.x writes nr points at a/nr: no end point
    chg = REPO / "examples" / "vasp_h2o"
    for name in ("CHGCAR", "PARCHG", "LOCPOT"):
        if (chg / name).is_file():
            assert read_grid(chg / name).periodic is True


def test_find_files_lists_cubes_and_vasp_grids(tmp_path):
    from adit.analysis.volumetric import find_files

    _write_cube(tmp_path / "wp-1-1-3-real.cube", _gaussian_grid(n=6))
    (tmp_path / "CHGCAR").write_text("")
    (tmp_path / "output.log").write_text("")
    assert [p.name for p in find_files(tmp_path)] == ["CHGCAR", "wp-1-1-3-real.cube"]


# ---- desktop ----
@pytest.fixture(scope="module")
def qapp():
    pytest.importorskip("PySide6")
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from PySide6.QtWidgets import QApplication

    return QApplication.instance() or QApplication([])


def test_viewer_draws_surfaces_sorted_with_the_atoms(qapp):
    from PySide6.QtGui import QColor

    from adit.gui.viewer3d import Viewer3D

    g = _gaussian_grid(n=24)
    m = isosurface(g, 0.5)
    atoms = g.atoms.copy(); atoms.pbc = False                          # centre on the atom, as the box does for non-periodic cubes
    v = Viewer3D(); v.resize(300, 300); v.set_atoms(atoms)
    plain = v.grab().toImage()
    v.set_surfaces([(m.triangles(), m.face_normals(), QColor(40, 110, 230, 140))])
    assert v.n_triangles() == m.n_triangles
    img = v.grab().toImage()
    c = img.pixelColor(150, 150)
    assert c.blue() > c.red() + 40                                      # the blue surface is over the centre
    assert plain.pixelColor(150, 150).blue() <= plain.pixelColor(150, 150).red() + 40
    v.set_selection([0]); assert v.selected == [0]                     # selection still works with a surface
    v.clear_surfaces(); assert v.n_triangles() == 0


def test_the_analysis_box_lists_files_shows_the_range_and_draws(qapp, tmp_path):
    from adit.gui.help import help_for
    from adit.gui.i18n import _EN
    from adit.gui.panels.analysis_panel import AnalysisPanel

    _write_cube(tmp_path / "wp-1-1-2-real.cube", _gaussian_grid(n=16, kind="p"))
    panel = AnalysisPanel()
    panel.set_run_dir(tmp_path, ["H"])
    assert not panel.iso_box.isHidden()
    iso = panel.iso
    assert iso.files() == [str(tmp_path / "wp-1-1-2-real.cube")]
    assert "10 %" in iso.range.text() and "最小" in iso.range.text()
    assert iso.level.text() == ""                                       # no default level
    assert iso.show_surface() is False and "等値を入れて" in iso.note.text()
    iso.level.setText("0.2")
    assert iso.show_surface() is True
    assert len(iso.meshes()) == 2 and iso.viewer.n_triangles() > 0
    assert iso.colour(0.2).name() != iso.colour(-0.2).name()
    iso.opacity.setValue(30); assert iso.colour(0.2).alpha() == round(255 * 0.3)
    iso.clear_surface(); assert iso.viewer.n_triangles() == 0
    panel.set_run_dir(tmp_path / "nothing", [])
    assert panel.iso_box.isHidden()
    for key in ("等値", "間引き (格子点)", "透明度", "等値面を表示"):
        h = help_for(key)
        assert h is not None and h.ja and h.en, key
        assert key in _EN
    assert help_for("等値").required


def test_the_box_suggests_a_stride_for_a_large_grid(qapp, tmp_path):
    from adit.gui.isosurface_panel import IsosurfacePanel

    _write_cube(tmp_path / "big.cube", _gaussian_grid(n=24))
    p = IsosurfacePanel(); p.max_triangles = 300
    p.set_run_dir(tmp_path); p.level.setText("0.5")
    assert p.show_surface() is False and not p.btn_use_stride.isHidden()
    p._use_suggested_stride()
    assert p.stride.value() >= 2 and p.viewer.n_triangles() > 0


# ---- browser ----
def test_web_payload_stats_then_surfaces_and_the_stride_suggestion(tmp_path):
    from adit.web.isosurface import MAX_WEB_TRIANGLES, isosurface_json, isosurface_payload

    _write_cube(tmp_path / "orb.cube", _gaussian_grid(n=20, kind="p"))
    d = isosurface_payload(tmp_path, "orb.cube")
    assert d["surfaces"] == [] and d["stats"]["has_negative"] and d["scene"]["n_atoms"] == 1
    d = isosurface_payload(tmp_path, "orb.cube", level="0.2")
    assert [s["level"] for s in d["surfaces"]] == [0.2, -0.2] and all(s["n_triangles"] > 0 for s in d["surfaces"])
    small = isosurface_payload(tmp_path, "orb.cube", level="0.2", max_triangles=100)
    assert small["surfaces"] == [] and small["suggested_stride"] >= 2 and small["notes"]
    assert MAX_WEB_TRIANGLES > MAX_TRIANGLES
    with pytest.raises(IsosurfaceError):
        isosurface_payload(tmp_path, "../orb.cube", level="0.2")
    with pytest.raises(IsosurfaceError):
        isosurface_payload(tmp_path, "orb.cube", level="abc")
    assert "<" not in isosurface_json(tmp_path, "orb.cube", level="0.2")


def test_the_scripts_have_the_surface_hooks_and_load_nothing():
    static = REPO / "src" / "adit" / "web" / "static"
    js = (static / "viewer3d.js").read_text(encoding="utf-8")
    assert "setSurfaces" in js and "projectSurfaces" in js
    iso = (static / "isosurface.js").read_text(encoding="utf-8")
    for bad in ("http://", "https://", "import ", "require("):
        assert bad not in iso
    assert "ADIT_ISO_URL" in iso and "adit-iso-level" in iso and "adit-iso-use-stride" in iso
    html = (REPO / "src" / "adit" / "web" / "templates" / "analysis.html").read_text(encoding="utf-8")
    assert 'id="adit-iso"' in html and 'id="adit3d-iso"' in html and "isosurface_js" in html


def test_analysis_page_serves_the_isosurface_only_for_the_analyzed_directory(tmp_path):
    import shutil
    import threading
    import urllib.error
    import urllib.parse
    import urllib.request

    from adit.config import default_config
    from adit.web.server import WebApp, serve

    run_dir = tmp_path / "run"
    shutil.copytree(REPO / "examples" / "water_generated", run_dir, ignore=shutil.ignore_patterns("analysis"))
    _write_cube(run_dir / "wp-1-1-1-real.cube", _gaussian_grid(n=17))
    app = WebApp(default_config(sk_root=""), tmp_path / "cluster.toml")
    httpd = serve(app, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_port}"
    try:
        with pytest.raises(urllib.error.HTTPError):
            urllib.request.urlopen(base + "/isosurface.json?dir=" + urllib.parse.quote(str(run_dir)) + "&file=wp-1-1-1-real.cube")
        html = urllib.request.urlopen(base + "/analysis", data=urllib.parse.urlencode({"run_dir": str(run_dir)}).encode()).read().decode()
        assert 'id="adit-iso"' in html and "ADIT_ISO_URL" in html and "ADIT_VIEWER" in html and "wp-1-1-1-real.cube" in html
        m = re.search(r'window.ADIT_ISO_URL = "([^"]+)"', html)
        assert m
        data = json.loads(urllib.request.urlopen(base + m.group(1) + "&file=wp-1-1-1-real.cube&level=0.5").read().decode())
        assert data["surfaces"][0]["n_triangles"] > 0 and data["stats"]["max"] == pytest.approx(1.0, abs=1e-3)
        with pytest.raises(urllib.error.HTTPError):
            urllib.request.urlopen(base + m.group(1) + "&file=nothing.cube&level=0.5")
    finally:
        httpd.shutdown(); httpd.server_close()
