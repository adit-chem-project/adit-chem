
from __future__ import annotations

import html as h
import re
import threading
import urllib.parse
import urllib.request

import pytest
from ase.build import molecule
from ase.io import write

from tests.test_web import _fields_without_label
from adit.config import default_config
from adit.spec import MlipMethod, OrcaMethod
from adit.web.forms import default_form, form_from_spec, spec_from_form
from adit.web.server import WebApp, serve
from tests.conftest import water_spec

UPF = """<UPF version="2.0.1">
  <PP_HEADER element="O" z_valence="6.0" wfc_cutoff="3.5000000000E+01" rho_cutoff="2.8000000000E+02"/>
  <PP_MESH>
</UPF>
"""


@pytest.fixture
def web(sk_root, tmp_path, monkeypatch):
    monkeypatch.setenv("ADIT_CONFIG", str(tmp_path / "conf" / "cluster.toml"))
    cfg = default_config(sk_root=str(sk_root))
    app = WebApp(cfg, tmp_path / "cluster.toml")
    app.form["output_dir"] = str(tmp_path / "out")
    httpd = serve(app, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    yield app, f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown(); httpd.server_close()


def _post(url: str, fields) -> str:
    return urllib.request.urlopen(url, data=urllib.parse.urlencode(fields).encode(), timeout=120).read().decode()


def _error(page: str) -> str:
    m = re.search(r'<div class="error">(.*?)</div>', page, re.S)
    return h.unescape(m.group(1)) if m else ""


def _base(app, **kw) -> dict:
    f = {**default_form(output_dir=app.form["output_dir"]), "sk_set": "fake-1-0"}
    f.update(kw)
    return f


def _carbon(app, **kw) -> dict:
    return _base(app, source="bulk", bulk_el="C", bulk_struct="diamond", **kw)


def test_batch_compare_reaction(web, tmp_path):
    app, base = web
    nh3 = tmp_path / "nh3.xyz"; write(nh3, molecule("NH3"))
    out = tmp_path / "set"
    page = _post(base + "/batch", _base(app, batch_kind="compare", cmp_kind="reaction", cmp_rx_rows="3",
                                        cmp_rx1_name="water", cmp_rx1_nu="-1", cmp_rx2_name="nh3", cmp_rx2_file=str(nh3),
                                        cmp_rx2_nu="1", batch_dir=str(out)))
    assert not _error(page), _error(page)
    assert (out / "compare.json").is_file() and (out / "water" / "dftb_in.hsd").is_file()
    assert "/compare?dir=" in page
    page = _post(base + "/batch", _base(app, batch_kind="compare", cmp_kind="reaction", cmp_rx_rows="3", batch_action="add_row"))
    assert 'name="cmp_rx4_name"' in page and app.form["cmp_rx_rows"] == "4"


def test_batch_conformers(web, tmp_path):
    pytest.importorskip("rdkit")
    app, base = web
    out = tmp_path / "conf"
    page = _post(base + "/batch", _base(app, source="smiles", smiles="CCO", batch_kind="conformers", conf_n="4", conf_rmsd="0.5",
                                        batch_dir=str(out)))
    assert not _error(page), _error(page)
    assert (out / "conformers.json").is_file() and (out / "scan.json").is_file()
    assert "/analysis?dir=" in page


def test_batch_neb_images_and_native_message(web, tmp_path):
    app, base = web
    end = tmp_path / "end.xyz"
    a = molecule("H2O"); a.positions[1] += (0.0, 0.0, 0.3); write(end, a)
    out = tmp_path / "neb"
    page = _post(base + "/batch", _base(app, batch_kind="neb", neb_end=str(end), neb_images="3", neb_mode="images", batch_dir=str(out)))
    assert not _error(page), _error(page)
    assert sorted(p.name for p in out.glob("image_*")) == ["image_00", "image_01", "image_02", "image_03", "image_04"]
    page = _post(base + "/batch", _base(app, batch_kind="neb", neb_end=str(end), neb_images="3", neb_mode="native", batch_dir=str(tmp_path / "neb2")))
    assert "VASP" in _error(page)


def test_batch_phonons_and_elastic(web, tmp_path):
    app, base = web
    out = tmp_path / "ph"
    page = _post(base + "/batch", _carbon(app, batch_kind="phonons", ph_dim="1x1x1", ph_backend="ase", batch_dir=str(out)))
    assert not _error(page), _error(page)
    assert len(list(out.glob("disp-*"))) == 12 and (out / "phonon_collect.py").is_file()
    assert "phonon_collect.py" in h.unescape(page)
    out2 = tmp_path / "el"
    off = {f"el_c{j}": "" for j in range(2, 7)}
    page = _post(base + "/batch", _carbon(app, batch_kind="elastic", el_strains="-0.01,0.01", el_c1="on", batch_dir=str(out2), **off))
    assert not _error(page), _error(page)
    assert sorted(p.name for p in out2.glob("e*") if p.is_dir()) == ["e0", "e1_+0.01", "e1_-0.01"] and (out2 / "elastic_collect.py").is_file()


def test_batch_ts_stages(web, tmp_path):
    app, base = web
    out = tmp_path / "ts"
    page = _post(base + "/batch", _base(app, code="orca", task_type="geometry_optimization", batch_kind="ts", ts_mode="ts+irc",
                                        batch_dir=str(out)))
    assert not _error(page), _error(page)
    assert (out / "stage_01_ts" / "orca.inp").is_file() and "OptTS" in (out / "stage_01_ts" / "orca.inp").read_text(encoding="utf-8")


def test_batch_errors_name_the_field(web, tmp_path):
    app, base = web
    page = _post(base + "/batch", _carbon(app, batch_kind="phonons", ph_dim="", batch_dir=str(tmp_path / "x")))
    assert "超格子の倍率" in _error(page)
    page = _post(base + "/batch", _base(app, batch_kind="conformers", conf_n="4", conf_rmsd="", batch_dir=str(tmp_path / "y")))
    assert "RMSD" in _error(page)


def test_mlip_form_roundtrip_and_page(web):
    app, base = web
    s = water_spec(method=MlipMethod(model_family="mace_mp", model="small", device="cpu", dtype="float64", dispersion=True, seed=7))
    back = spec_from_form({**default_form(), **form_from_spec(s)})
    assert back.method == s.method and back.kpoints is None
    page = _post(base + "/preview", _base(app, code="mlip", mlip_family="mace_mp", mlip_model="small", task_type="single_point"))
    assert not _error(page), _error(page)
    assert app.spec.method.model_family == "mace_mp" and "run_mlip.py" in app.files.texts
    page = h.unescape(page)
    assert "pip install ase mace-torch" in page
    assert '#main:has(#code option[value="mlip"]:checked) .kp-group { display:none; }' in page
    assert _fields_without_label(page) == []


def test_orca_ts_irc_roundtrip_and_rules(web):
    app, base = web
    for m in (OrcaMethod(ts_search=True, ts_calc_hess=False, ts_recalc_hess=5, ts_freq=False),
              OrcaMethod(irc=True, irc_max_iter=30, irc_direction="forward")):
        s = water_spec(method=m)
        back = spec_from_form({**default_form(), **form_from_spec(s)})
        assert back.method == m, m
    page = _post(base + "/preview", _base(app, code="orca", task_type="geometry_optimization", orca_ts="on", orca_recalc_hess="5"))
    assert app.spec.method.ts_search and app.spec.method.ts_recalc_hess == 5
    assert '#main:has(#orca_ts:checked) .only[data-orcats~="on"]' in page


def test_orca_cli_options_survive_web_roundtrip(web, tmp_path):
    app, base = web
    m = OrcaMethod(goat=True, docker_guest_file="guests.xyz", docker_assume_neutral_singlet=True)
    spec = water_spec(method=m)
    fields = {**default_form(), **form_from_spec(spec)}
    assert spec_from_form(fields).method == m
    guest = tmp_path / "guests.xyz"
    guest.write_text("1\n0 1\nH 0 0 0\n", encoding="utf-8")
    m = OrcaMethod(method="XTB", basis="", docker_guest_file=str(guest), docker_assume_neutral_singlet=True)
    fields = {**default_form(), **form_from_spec(water_spec(method=m))}
    page = _post(base + "/preview", fields)
    assert app.spec.method == m
    assert f'name="orca_docker_guest_file" value="{guest}"' in page


def test_documented_values_shown_and_used(web, tmp_path):
    app, base = web
    lib = tmp_path / "pseudo" / "SSSP_test"; lib.mkdir(parents=True)
    (lib / "O.upf").write_text(UPF, encoding="utf-8")
    (lib / "H.upf").write_text(UPF.replace('element="O"', 'element="H"'), encoding="utf-8")
    app.cfg.pseudo_root = str(tmp_path / "pseudo")
    page = h.unescape(_post(base + "/preview", _base(app, code="espresso", pseudo_set="SSSP_test", ecutwfc="0")))
    assert "O.upf:2 の記載: O wfc_cutoff = 35 Ry" in page, page[:0]
    assert 'value="ecutwfc=35"' in page
    page = _post(base + "/docvalue", {**app.form, "doc_put": "ecutwfc=35"})
    assert app.form["ecutwfc"] == "35" and app.spec.method.ecutwfc == 35.0
    page = _post(base + "/docvalue", {**app.form, "doc_put": "ecutwfc=../../etc/passwd"})
    assert app.form["ecutwfc"] == "35"


def test_english_page(web):
    from adit import lang
    app, base = web
    lang.set_language("en")
    try:
        page = h.unescape(_post(base + "/preview", _base(app, code="mlip", mlip_family="chgnet", task_type="single_point")))
        for text in ("Batch generation", "Machine-learning potential", "Transition-state search (OptTS)", "Supercell",
                     "Strain magnitudes", "What to generate"):
            assert text in page, text
    finally:
        lang.set_language("ja")
