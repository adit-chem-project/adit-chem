
from __future__ import annotations

import html as h
import http.client
import json
import os
import threading
import urllib.error
import urllib.parse
import urllib.request

import pytest

from adit.config import default_config
from adit.spec import CalculationSpec, Task
from adit.web.forms import FormError, _f, _i, default_form, form_from_spec, spec_from_form
from adit.web.server import WebApp, serve
from tests.conftest import make_fake_skset, water_spec


@pytest.fixture
def web(tmp_path):
    root = tmp_path / "slakos"
    make_fake_skset(root, "fake-1-0", ["H", "C", "N", "O", "Si", "Al"])
    app = WebApp(default_config(sk_root=str(root)), tmp_path / "cluster.toml")
    httpd = serve(app, "127.0.0.1", 0)
    t = threading.Thread(target=httpd.serve_forever, daemon=True); t.start()
    yield app, f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown(); httpd.server_close()


def _post(url: str, fields: dict) -> str:
    return urllib.request.urlopen(url, data=urllib.parse.urlencode(fields).encode(), timeout=60).read().decode()


def _error(page: str) -> str:
    import re
    m = re.search(r'<div class="error">(.*?)</div>', page, re.S)
    return h.unescape(m.group(1)) if m else ""


@pytest.mark.parametrize("v", ["nan", "inf", "-inf", "1e400", "NaN", "Infinity"])
def test_number_fields_reject_non_finite(v):
    with pytest.raises(FormError, match="有限|finite"):
        _f(v, 0.0)
    with pytest.raises(FormError, match="有限|finite"):
        _i(v, 0)


BAD = {
    "md_steps inf": dict(task_type="molecular_dynamics", md_steps="inf"),
    "md_steps 1e400": dict(task_type="molecular_dynamics", md_steps="1e400"),
    "surf_nx inf": dict(source="surface", surf_nx="inf"),
    "kp_density nan": dict(source="bulk", bulk_el="Si", kp_mode="density", kp_density="nan"),
    "band QQQ": dict(source="bulk", bulk_el="Si", kp_mode="mesh", task_type="band_structure", band_path="QQQ"),
    "surf_vac nan": dict(source="surface", surf_vac="nan"),
    "box nan": dict(box="on", box_size="nan"),
    "magmom nan": dict(code="vasp", magmom="1 nan"),
}


@pytest.mark.parametrize("name", list(BAD))
def test_preview_bad_values_show_error(web, name):
    app, base = web
    f = {**app.form, "sk_set": "fake-1-0", **BAD[name]}
    page = _post(base + "/preview", f)
    err = _error(page)
    assert err, name
    assert "Traceback" not in err and "想定外" not in err, err
    assert app.spec is None


def test_band_path_error_names_the_point(web):
    app, base = web
    err = _error(_post(base + "/preview", {**app.form, "sk_set": "fake-1-0", **BAD["band QQQ"]}))
    assert "Q" in err and ("バンドの経路" in err or "band path" in err)


@pytest.mark.parametrize("fields", [{"rmax": "abc"}, {"skip": "inf"}, {"sigma": "nan"}, {"rmax": "1e400"}])
def test_analysis_bad_options_show_error(web, tmp_path, fields):
    app, base = web
    page = _post(base + "/analysis", {"run_dir": str(tmp_path), "rdf": "on", **fields})
    assert _error(page)


def test_unexpected_exception_is_500_page_not_disconnect(web, monkeypatch):
    app, base = web

    def boom(*a, **k):
        raise RuntimeError("boom-for-test")
    monkeypatch.setattr(app, "choices", boom)
    with pytest.raises(urllib.error.HTTPError) as ei:
        urllib.request.urlopen(base + "/", timeout=60)
    assert ei.value.code == 500
    assert "boom-for-test" in ei.value.read().decode()


def test_pydantic_error_first_line_is_field_name(web):
    app, base = web
    err = _error(_post(base + "/preview", {**app.form, "sk_set": "fake-1-0", "code": "vasp", "ispin": "3"}))
    assert err and "validation error for" not in err.splitlines()[0], err


def _multipart(name: str, data: bytes) -> tuple[bytes, dict]:
    b = "----adittest"
    body = (f"--{b}\r\nContent-Disposition: form-data; name=\"{name}\"; filename=\"spec.json\"\r\nContent-Type: application/json\r\n\r\n").encode() \
        + data + f"\r\n--{b}--\r\n".encode()
    return body, {"Content-Type": f"multipart/form-data; boundary={b}"}


def test_load_old_spec_migrates_force_tolerance(web):
    app, base = web
    d = json.loads(water_spec().to_json())
    d.pop("version", None)
    d["task"].pop("force_tolerance_ev_per_ang")
    d["task"]["force_tolerance"] = 1e-3
    expected = CalculationSpec.from_json(json.dumps(d)).task.force_tolerance_ev_per_ang
    assert abs(expected - Task().force_tolerance_ev_per_ang) > 1e-6
    body, hdr = _multipart("spec_file", json.dumps(d).encode())
    urllib.request.urlopen(urllib.request.Request(base + "/load", data=body, headers=hdr), timeout=60).read()
    assert float(app.form["force_tol"]) == pytest.approx(expected, rel=1e-5)


@pytest.mark.parametrize("data", [b"not json", b"[1, 2]", b"\xff\xfe", b'{"version": "x"}'])
def test_load_broken_spec_shows_error(web, data):
    app, base = web
    body, hdr = _multipart("spec_file", data)
    page = urllib.request.urlopen(urllib.request.Request(base + "/load", data=body, headers=hdr), timeout=60).read().decode()
    assert "spec.json" in _error(page)


def test_box_roundtrip_keeps_pbc_and_kpoints():
    f = default_form(); f.update(sk_set="fake-1-0", box="on", box_size="12.5", kp_mode="mesh", k1="2", k2="2", k3="2")
    s = spec_from_form(f)
    assert s.structure.periodic and s.kpoints is not None and s.kpoints.mesh == (2, 2, 2)
    back = form_from_spec(s)
    assert back["box"] == "on" and float(back["box_size"]) == 12.5
    s2 = spec_from_form({**default_form(), **back})
    assert s2.model_dump(exclude={"meta"}) == s.model_dump(exclude={"meta"})


def test_loading_molecule_clears_previous_box():
    f = default_form(); f["sk_set"] = "fake-1-0"
    s = spec_from_form(f)
    back = form_from_spec(s)
    assert back["box"] == ""
    s2 = spec_from_form({**default_form(), "box": "on", **back})
    assert not s2.structure.periodic and s2.kpoints is None


def _raw_post(base: str, path: str, headers: dict, body: bytes) -> int:
    u = urllib.parse.urlparse(base)
    conn = http.client.HTTPConnection(u.hostname, u.port, timeout=60)
    try:
        conn.request("POST", path, body=body, headers=headers)
        r = conn.getresponse(); r.read()
        return r.status
    finally:
        conn.close()


def test_content_length_is_checked(web):
    _, base = web
    form = {"Content-Type": "application/x-www-form-urlencoded"}
    assert _raw_post(base, "/preview", {**form, "Content-Length": "abc"}, b"x=1") == 400
    assert _raw_post(base, "/preview", {**form, "Content-Length": "-1"}, b"x=1") == 400
    assert _raw_post(base, "/preview", {**form, "Content-Length": str(10 ** 12)}, b"x=1") == 413


def test_superscript_digits_do_not_crash_the_recipe_page(web):
    app, base = web
    app.form.update(sk_set="fake-1-0", source="bulk", bulk_el="Si")
    _post(base + "/recipe", {**app.form, "recipe_action": "add", "recipe_add": "slab"})
    assert app.form.get("st1_op") == "slab"
    page = _post(base + "/recipe", {**app.form, "recipe_action": "del:\u00b2"})
    assert "内部エラー" not in page and app.form.get("st1_op") == "slab"
    page = _post(base + "/recipe", {**app.form, "st1_term": "\u00b2", "recipe_action": "build"})
    assert "内部エラー" not in page


def test_unreadable_directory_is_reported_not_500(web, tmp_path):
    if os.geteuid() == 0:
        pytest.skip("root can read everything")
    _, base = web
    locked = tmp_path / "locked"; inner = locked / "inner"
    inner.mkdir(parents=True); locked.chmod(0)
    try:
        page = urllib.request.urlopen(base + "/analysis?dir=" + urllib.parse.quote(str(inner)), timeout=60).read().decode()
        assert "内部エラー" not in page and "解析できません" in h.unescape(page)
        page = _post(base + "/compare", {"base": str(inner), "action": "compare"})
        assert "ディレクトリがありません" in _error(page)
    finally:
        locked.chmod(0o700)


def test_status_endpoint_is_gone(web):
    _, base = web
    with pytest.raises(urllib.error.HTTPError) as ex:
        urllib.request.urlopen(base + "/status", timeout=60)
    assert ex.value.code == 404


def test_requests_that_change_state_run_one_at_a_time(web, monkeypatch):
    import time

    app, base = web
    active, peak, lock = [0], [0], threading.Lock()
    original = app.preview

    def slow_preview(form):
        with lock:
            active[0] += 1; peak[0] = max(peak[0], active[0])
        try:
            time.sleep(0.3)
            return original(form)
        finally:
            with lock:
                active[0] -= 1

    monkeypatch.setattr(app, "preview", slow_preview)
    fields = {**default_form(), "sk_set": "fake-1-0"}
    threads = [threading.Thread(target=_post, args=(base + "/preview", fields)) for _ in range(3)]
    for t in threads:
        t.start()
    for t in threads:
        t.join()
    assert peak[0] == 1
