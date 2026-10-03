"""Fixes from the 2026-10-03 audit of the web front end: request-origin checks, body handling, round trips, error texts."""

from __future__ import annotations

import html as h
import http.client
import json
import os
import re
import shutil
import socket
import subprocess
import threading
import time
import urllib.error
import urllib.parse
import urllib.request
import webbrowser
from html.parser import HTMLParser
from pathlib import Path

import pytest
from ase.build import bulk, fcc111, molecule

from adit.compat import upload_basename
from adit.config import default_config
from adit.spec import DftbMethod, Task
from adit.web import server as S
from adit.web.forms import FormError, _i, default_form, form_from_spec, spec_from_form
from adit.web.server import WebApp, serve
from adit.web.structure3d import scene_from_atoms
from tests.conftest import make_fake_skset, water_spec

REPO = Path(__file__).resolve().parent.parent
JS = REPO / "src" / "adit" / "web" / "static" / "viewer3d.js"


@pytest.fixture
def web(tmp_path):
    root = tmp_path / "slakos"
    make_fake_skset(root, "fake-1-0", ["H", "C", "N", "O", "Si", "Al"])
    app = WebApp(default_config(sk_root=str(root)), tmp_path / "cluster.toml")
    app.form["output_dir"] = str(tmp_path / "out")
    httpd = serve(app, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    yield app, f"http://127.0.0.1:{httpd.server_port}", httpd
    httpd.shutdown(); httpd.server_close()


def _request(base: str, method: str, path: str, body: bytes | None = None, headers: dict | None = None) -> tuple[int, dict, str]:
    u = urllib.parse.urlparse(base)
    conn = http.client.HTTPConnection(u.hostname, u.port, timeout=60)
    try:
        hdr = {"Content-Type": "application/x-www-form-urlencoded"} if body is not None else {}
        hdr.update(headers or {})
        conn.request(method, path, body=body, headers=hdr)
        r = conn.getresponse()
        return r.status, {k.lower(): v for k, v in r.getheaders()}, r.read().decode("utf-8", errors="replace")
    finally:
        conn.close()


def _post(url: str, fields) -> str:
    data = urllib.parse.urlencode(fields).encode()
    return urllib.request.urlopen(url, data=data, timeout=120).read().decode()


def _error(page: str) -> str:
    m = re.search(r'<div class="error">(.*?)</div>', page, re.S)
    return h.unescape(m.group(1)) if m else ""


def _multipart(name: str, data: bytes, filename: str = "spec.json") -> tuple[bytes, dict]:
    b = "----adit20261003"
    body = (f"--{b}\r\nContent-Disposition: form-data; name=\"{name}\"; filename=\"{filename}\"\r\nContent-Type: application/json\r\n\r\n").encode() \
        + data + f"\r\n--{b}--\r\n".encode()
    return body, {"Content-Type": f"multipart/form-data; boundary={b}"}


def _load(base: str, spec_json: bytes) -> str:
    body, hdr = _multipart("spec_file", spec_json)
    return urllib.request.urlopen(urllib.request.Request(base + "/load", data=body, headers=hdr), timeout=120).read().decode()


class _FormScraper(HTMLParser):
    """What a browser would submit for the main form: inputs, checked boxes, the selected (or first) option, textareas."""

    def __init__(self):
        super().__init__()
        self.items: list[tuple[str, str]] = []
        self._select = None
        self._textarea = None

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "input":
            kind, name = (a.get("type") or "text").lower(), a.get("name")
            if not name or kind in ("submit", "button", "file"):
                return
            if kind in ("checkbox", "radio"):
                if "checked" in a:
                    self.items.append((name, a.get("value") or "on"))
            else:
                self.items.append((name, a.get("value") or ""))
        elif tag == "select":
            self._select = [a.get("name"), None, None]
        elif tag == "option" and self._select is not None:
            v = a.get("value") or ""
            if self._select[1] is None:
                self._select[1] = v
            if "selected" in a:
                self._select[2] = v
        elif tag == "textarea":
            self._textarea = (a.get("name"), [])

    def handle_endtag(self, tag):
        if tag == "select" and self._select is not None:
            name, first, chosen = self._select
            if name:
                self.items.append((name, chosen if chosen is not None else (first or "")))
            self._select = None
        elif tag == "textarea" and self._textarea is not None:
            self.items.append((self._textarea[0], "".join(self._textarea[1])))
            self._textarea = None

    def handle_data(self, data):
        if self._textarea is not None:
            self._textarea[1].append(data)


def _browser_submission(page: str) -> list[tuple[str, str]]:
    p = _FormScraper()
    p.feed(page)
    return [(k, v) for k, v in p.items if k]


# ---- finding 1: requests from other origins and for other hosts are refused ----

def test_cross_origin_post_is_refused_and_writes_nothing(web, tmp_path):
    app, base, _ = web
    out = tmp_path / "evil_out"
    fields = {**default_form(), "sk_set": "fake-1-0", "output_dir": str(out)}
    body = urllib.parse.urlencode(fields).encode()
    status, _, text = _request(base, "POST", "/generate", body, {"Origin": "http://evil.example"})
    assert status == 403 and not out.exists()
    assert "CSRF" in text
    status, _, _ = _request(base, "POST", "/generate", body, {"Referer": "http://evil.example/attack.html"})
    assert status == 403 and not out.exists()
    status, _, _ = _request(base, "POST", "/generate", body, {"Origin": "null"})
    assert status == 403 and not out.exists()
    status, _, _ = _request(base, "POST", "/generate", body, {"Origin": f"http://evil.example:{urllib.parse.urlparse(base).port + 1}"})
    assert status == 403 and not out.exists()
    assert app.written is None


def test_same_origin_post_works(web, tmp_path):
    app, base, _ = web
    port = urllib.parse.urlparse(base).port
    out = tmp_path / "good_out"
    fields = {**default_form(), "sk_set": "fake-1-0", "output_dir": str(out)}
    body = urllib.parse.urlencode(fields).encode()
    status, _, page = _request(base, "POST", "/generate", body, {"Origin": f"http://127.0.0.1:{port}"})
    assert status == 200 and not _error(page) and (out / "dftb_in.hsd").is_file()
    status, _, page = _request(base, "POST", "/generate", body + b"&overwrite=on", {"Referer": f"http://localhost:{port}/"})
    assert status == 200 and not _error(page)
    status, _, page = _request(base, "POST", "/preview", body)        # scripts and curl send neither header
    assert status == 200 and not _error(page)


def test_host_header_must_name_this_server(web):
    _, base, _ = web
    port = urllib.parse.urlparse(base).port
    assert _request(base, "GET", "/", headers={"Host": "evil.example"})[0] == 403
    assert _request(base, "GET", "/", headers={"Host": f"evil.example:{port}"})[0] == 403
    assert _request(base, "GET", "/spec.json", headers={"Host": f"evil.example:{port}"})[0] == 403
    # from loopback a loopback name may carry another port: ssh -L 9000:127.0.0.1:8765 shows up as localhost:9000
    assert _request(base, "GET", "/", headers={"Host": f"127.0.0.1:{port + 1}"})[0] == 200
    assert _request(base, "GET", "/", headers={"Host": f"evil.example:{port + 1}"})[0] == 403
    assert _request(base, "GET", "/", headers={"Host": f"127.0.0.1:{port}"})[0] == 200
    assert _request(base, "GET", "/", headers={"Host": f"localhost:{port}"})[0] == 200
    assert _request(base, "GET", "/", headers={"Host": f"LOCALHOST:{port}"})[0] == 200
    # a link to adit-web from another page is an ordinary GET with a foreign Referer; it must still open
    assert _request(base, "GET", "/", headers={"Referer": "https://wiki.example/adit"})[0] == 200
    assert _request(base, "GET", "/", headers={"Origin": "https://wiki.example"})[0] == 403


def test_host_allowed_rules(monkeypatch):
    monkeypatch.setattr(S, "_OWN_NAMES", {"pc.example.org", "pc"})
    ok = dict(server_port=8765)
    assert S.host_allowed("127.0.0.1", 8765, **ok) and S.host_allowed("::1", 8765, **ok) and S.host_allowed("localhost", 8765, **ok)
    assert S.host_allowed("127.0.0.2", 8765, **ok)                                # any loopback address
    assert not S.host_allowed("127.0.0.1", 80, **ok) and not S.host_allowed("127.0.0.1", None, **ok)
    assert S.host_allowed("127.0.0.1", None, server_port=80)
    # ssh -L 9000:127.0.0.1:8765: the browser says localhost:9000 and the connection comes from loopback
    assert S.host_allowed("localhost", 9000, server_port=8765, peer_loopback=True)
    assert S.host_allowed("127.0.0.1", 9000, server_port=8765, peer_loopback=True)
    assert not S.host_allowed("localhost", 9000, server_port=8765)                 # not from loopback: the port must match
    assert not S.host_allowed("pc.example.org", 9000, server_port=8765, peer_loopback=True)   # only loopback names
    assert S.host_allowed("pc.example.org", 8765, **ok) and S.host_allowed("PC.example.org.", 8765, **ok) and S.host_allowed("pc", 8765, **ok)
    assert not S.host_allowed("evil.example", 8765, **ok) and not S.host_allowed("", 8765, **ok)
    assert S.host_allowed("192.168.1.5", 8765, server_port=8765, local_ip="192.168.1.5")          # the address the client connected to
    assert S.host_allowed("192.168.1.5", 8765, server_port=8765, local_ip="::ffff:192.168.1.5")   # same, through a dual-stack socket
    assert not S.host_allowed("10.0.0.9", 8765, server_port=8765, local_ip="192.168.1.5")
    assert S.host_allowed("lab-node", 8765, server_port=8765, bound_host="lab-node")
    assert S.split_authority("[::1]:8765") == ("::1", 8765) and S.split_authority("Evil.Example") == ("evil.example", None)
    assert S.split_authority("127.0.0.1:abc") is None and S.split_authority("") is None


def test_open_browser_url_still_opens(tmp_path, monkeypatch):
    seen = []
    monkeypatch.setattr(webbrowser, "open", lambda url: seen.append(url) or True)
    app = WebApp(default_config(sk_root=""), tmp_path / "cluster.toml")
    httpd = serve(app, "127.0.0.1", 0, open_browser=True)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    try:
        for _ in range(100):
            if seen:
                break
            time.sleep(0.1)
        assert seen and seen[0] == f"http://127.0.0.1:{httpd.server_port}/"
        assert urllib.request.urlopen(seen[0], timeout=60).status == 200
    finally:
        httpd.shutdown(); httpd.server_close()


# ---- finding 2: a select keeps a loaded value that is not among its options ----

def test_select_keeps_values_outside_the_options_through_a_browser_round_trip(web):
    app, base, _ = web
    f = {**default_form(), "source": "bulk", "bulk_el": "Si", "code": "vasp", "ivdw": "21", "ibrion": "5", "kp_mode": "mesh", "k1": "2", "k2": "2", "k3": "2"}
    spec = spec_from_form(f)
    assert spec.method.ivdw == 21 and spec.method.ibrion == 5
    page = _load(base, spec.to_json().encode())
    assert app.spec is not None and app.spec.method.ivdw == 21 and app.spec.method.ibrion == 5
    assert '<option value="21" selected>21</option>' in page and '<option value="5" selected>5</option>' in page
    page = _post(base + "/preview", _browser_submission(page))      # post the rendered form back unchanged
    assert not _error(page), _error(page)
    assert app.spec.method.ivdw == 21 and app.spec.method.ibrion == 5


def test_unknown_profile_stays_selected(web):
    app, base, _ = web
    app.form["profile"] = "rccs"
    page = urllib.request.urlopen(base + "/", timeout=60).read().decode()
    assert '<option value="rccs" selected>rccs (' in page
    assert dict(_browser_submission(page))["profile"] == "rccs"


# ---- finding 3: the DFTB+ seed survives load -> preview ----

def test_dftb_seed_survives_load_and_preview(web):
    app, base, _ = web
    spec = water_spec(method=DftbMethod(sk_set="fake-1-0", seed=12348))
    page = _load(base, spec.to_json().encode())
    assert app.spec is not None and app.spec.method.seed == 12348
    assert "seed=12348" in h.unescape(page)
    _post(base + "/preview", dict(app.form))
    assert app.spec.method.seed == 12348
    assert "RandomSeed = 12348" in app.files.texts["dftb_in.hsd"]


# ---- finding 4: the body is read before the lock, with a socket timeout ----

def test_a_stalled_upload_does_not_block_other_requests(web, monkeypatch):
    app, base, httpd = web
    port = urllib.parse.urlparse(base).port
    monkeypatch.setattr(httpd.RequestHandlerClass, "timeout", 6)
    s = socket.create_connection(("127.0.0.1", port), timeout=60)
    try:
        s.sendall((f"POST /preview HTTP/1.1\r\nHost: 127.0.0.1:{port}\r\nContent-Type: application/x-www-form-urlencoded\r\n"
                   "Content-Length: 100\r\n\r\nsource=preset").encode())
        time.sleep(0.3)
        t0 = time.monotonic()
        assert urllib.request.urlopen(base + "/", timeout=60).status == 200   # takes app.lock; must not wait for the stalled POST
        s.settimeout(0.05)
        try:
            early = s.recv(4096)
        except (TimeoutError, socket.timeout):
            early = b""
        assert early == b"", early                                            # the POST was still waiting for its body
        s.settimeout(60)
        answer = s.recv(4096)
        assert answer.split(b"\r\n", 1)[0].endswith(b"408 Request Timeout"), answer[:80]
        assert time.monotonic() - t0 < 40
    finally:
        s.close()


# ---- finding 5: body limits per route and a multipart parser without the extra copies ----

def test_oversize_body_is_refused_on_form_routes_but_uploads_keep_the_large_limit(web):
    _, base, _ = web
    status, _, text = _request(base, "POST", "/analysis", b"x=1", {"Content-Length": str(5 * 1024 * 1024)})
    assert status == 413 and str(S.MAX_FORM_BODY_BYTES) in text
    status, _, _ = _request(base, "POST", "/draw", b"x=1", {"Content-Length": str(5 * 1024 * 1024)})
    assert status == 413
    big = urllib.parse.urlencode({**default_form(), "sk_set": "fake-1-0", "pad": "x" * (5 * 1024 * 1024)}).encode()
    status, _, page = _request(base, "POST", "/preview", big)
    assert status == 200 and not _error(page)
    assert S.MAX_UPLOAD_BODY_BYTES == 256 * 1024 * 1024 and all(p in S.UPLOAD_PATHS for p in ("/load", "/preview", "/generate", "/scan"))


def test_multipart_parser_splits_fields_and_binary_files():
    b = "xYzBoundary"
    payload = b"\x00\x01--xYz\r\n\r\n--not-the-boundary\r\nend"
    body = (f"--{b}\r\nContent-Disposition: form-data; name=\"a\"\r\n\r\nhello\r\n"
            f"--{b}\r\nContent-Disposition: form-data; name=\"f\"; filename=\"x y.bin\"\r\nContent-Type: application/octet-stream\r\n\r\n").encode() \
        + payload + f"\r\n--{b}\r\nContent-Disposition: form-data; name=\"empty\"; filename=\"\"\r\n\r\n\r\n--{b}--\r\n".encode()
    parts = list(S._multipart_items(body, b))
    assert [bytes(p) for _, p in parts] == [b"hello", payload, b""]
    assert b'name="f"' in parts[1][0]
    assert S._boundary_of(f'multipart/form-data; boundary="{b}"') == b and S._boundary_of("multipart/form-data") is None
    lf_only = body.replace(b"\r\n", b"\n")
    assert [bytes(p) for _, p in S._multipart_items(lf_only, b)][0] == b"hello"


# ---- finding 6: other apps' cookies on the same host ----

@pytest.mark.parametrize("extra", ["foo(bar)=1", "a,b=1", "x=a b c"])
def test_other_apps_cookies_do_not_break_the_token(tmp_path, extra):
    app = WebApp(default_config(sk_root=""), tmp_path / "cluster.toml")
    httpd = serve(app, "127.0.0.1", 0, token="secret-token")
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_port}"
    name = S.cookie_name(httpd.server_port)
    try:
        assert _request(base, "GET", "/", headers={"Cookie": f"{name}=secret-token; {extra}"})[0] == 200
        assert _request(base, "GET", "/", headers={"Cookie": f"{extra}; {name}=secret-token"})[0] == 200
        assert _request(base, "GET", "/", headers={"Cookie": f"{extra}; {name}=wrong"})[0] == 403
    finally:
        httpd.shutdown(); httpd.server_close()


def test_cookie_value_parsing():
    assert S.cookie_value("a=1; adit_token_1=abc ; b=2", "adit_token_1") == "abc"
    assert S.cookie_value('adit_token_1="q=x"', "adit_token_1") == "q=x"
    assert S.cookie_value("foo(bar)=1; adit_token_1=t", "adit_token_1") == "t"
    assert S.cookie_value("", "adit_token_1") is None and S.cookie_value("other=1", "adit_token_1") is None


# ---- finding 7: /generate reports OS errors as text ----

@pytest.mark.skipif(os.name == "nt", reason="chmod does not lock a directory on Windows")
def test_generate_into_an_unwritable_parent_is_reported_not_500(web, tmp_path):
    if os.geteuid() == 0:
        pytest.skip("root can write everywhere")
    app, base, _ = web
    parent = tmp_path / "ro"; parent.mkdir(); parent.chmod(0o555)
    try:
        fields = {**default_form(), "sk_set": "fake-1-0", "output_dir": str(parent / "out")}
        status, _, page = _request(base, "POST", "/generate", urllib.parse.urlencode(fields).encode())
        assert status == 200 and "内部エラー" not in page
        assert _error(page) and ("Permission" in _error(page) or "許可" in _error(page) or "denied" in _error(page).lower())
        assert app.written is None
    finally:
        parent.chmod(0o755)


# ---- finding 8: frames and isosurfaces answer while the prepare page is busy ----

def test_frames_and_isosurface_do_not_wait_for_the_prepare_lock(web, tmp_path):
    app, base, _ = web
    assert "/frames.json" in S.ANALYSIS_PATHS and "/isosurface.json" in S.ANALYSIS_PATHS
    with app.lock:                      # a slow recipe build would hold this
        for path in ("/frames.json", "/isosurface.json"):
            status, headers, _ = _request(base, "GET", f"{path}?dir={urllib.parse.quote(str(tmp_path))}")
            assert status == 404 and headers["content-type"] == "text/plain; charset=utf-8"


# ---- finding 9: charset on the plain-text answers ----

def test_plain_text_errors_declare_utf8(web, tmp_path):
    app, base, _ = web
    for path in ("/structure.svg", "/frames.json?dir=x", "/isosurface.json?dir=x", "/file?path=x", "/no-such-page"):
        status, headers, _ = _request(base, "GET", path)
        assert status == 404 and headers["content-type"] == "text/plain; charset=utf-8", path
    app.analysis_dir = str(tmp_path)
    status, headers, text = _request(base, "GET", f"/isosurface.json?dir={urllib.parse.quote(str(tmp_path))}&stride=abc")
    assert status == 400 and headers["content-type"] == "text/plain; charset=utf-8"
    assert "invalid literal" not in text and ("間引き" in text or "stride" in text)
    status, _, _ = _request(base, "GET", f"/isosurface.json?dir={urllib.parse.quote(str(tmp_path))}&stride=0")
    assert status == 400


# ---- finding 10: odd input is an error page, not a 500 ----

def test_huge_recipe_index_and_odd_sketch_json_are_not_500(web):
    app, base, _ = web
    app.form.update(sk_set="fake-1-0", source="bulk", bulk_el="Si")
    _post(base + "/recipe", {**app.form, "recipe_action": "add", "recipe_add": "slab"})
    assert app.form.get("st1_op") == "slab"
    page = _post(base + "/recipe", {**app.form, "recipe_action": "del:" + "9" * 5000})
    assert "内部エラー" not in page and app.form.get("st1_op") == "slab"
    for sketch in ('{"atoms": "x"}', "[]", '{"atoms": [{"elem": 5}]}', '{"atoms": [[1, 2]]}'):
        status, _, page = _request(base, "POST", "/draw", urllib.parse.urlencode({"sketch": sketch, "action": "check"}).encode())
        assert status == 200 and "内部エラー" not in page, sketch


def test_upload_names_lose_nul_and_directories():
    assert upload_basename("a\x00b.xyz") == "ab.xyz"
    assert upload_basename("../../evil\x00.xyz") == "evil.xyz"
    assert upload_basename("C:\\Users\\x\\water.xyz") == "water.xyz"
    assert upload_basename("\x00") == "upload"


# ---- finding 11: integer fields reject fractions ----

def test_integer_fields_reject_fractions():
    with pytest.raises(FormError, match="整数|integer"):
        _i("2.9", 1)
    assert _i("2.0", 1) == 2 and _i("1e3", 0) == 1000 and _i("", 7) == 7 and _i("-3", 0) == -3
    with pytest.raises(FormError, match="整数|integer"):
        spec_from_form({**default_form(), "sk_set": "fake-1-0", "multiplicity": "2.9"})
    with pytest.raises(FormError, match="整数|integer"):
        spec_from_form({**default_form(), "sk_set": "fake-1-0", "task_type": "molecular_dynamics", "dump": "2.5"})


# ---- finding 12: the overwrite question counts files with a cutoff ----

def test_file_count_for_the_overwrite_question_has_a_cutoff(tmp_path):
    d = tmp_path / "many"; (d / "sub").mkdir(parents=True)
    for i in range(3):
        (d / f"f{i}").write_text("x", encoding="utf-8")
    (d / ".adit_backup").mkdir(); (d / ".adit_backup" / "old").write_text("x", encoding="utf-8")
    assert S._file_count_phrase(d) in ("3 ファイル", "3 files")
    for i in range(1005):
        (d / "sub" / f"g{i}").write_text("x", encoding="utf-8")
    assert S._file_count_phrase(d) in ("1000 以上のファイル", "more than 1000 files")
    assert S._file_count_phrase(d, limit=5) in ("5 以上のファイル", "more than 5 files")


# ---- finding 13: numbers survive the form round trip exactly ----

def test_floats_survive_the_form_round_trip_exactly():
    f = {**default_form(), "sk_set": "fake-1-0", "source": "bulk", "bulk_el": "Si", "bulk_a": "5.431020511", "force_tol": "0.005142206709048065"}
    s = spec_from_form(f)
    assert s.structure.source_ref == "Si 5.431020511" and s.task.force_tolerance_ev_per_ang == 0.005142206709048065
    back = form_from_spec(s)
    assert back["bulk_a"] == "5.431020511" and back["force_tol"] == "0.005142206709048065"
    s2 = spec_from_form({**f, **back})
    assert s2.structure.source_ref == s.structure.source_ref and s2.task.force_tolerance_ev_per_ang == s.task.force_tolerance_ev_per_ang
    assert float(default_form()["force_tol"]) == Task().force_tolerance_ev_per_ang
    assert spec_from_form({**default_form(), "sk_set": "fake-1-0"}).task.force_tolerance_ev_per_ang == Task().force_tolerance_ev_per_ang
    assert back["kp_shift"] == "0" and back["temperature"] == "300"


# ---- finding 14: missing fields take the page defaults ----

def test_missing_fields_fall_back_to_the_page_defaults():
    s = spec_from_form({"source": "preset", "preset": "H2O", "sk_set": "fake-1-0"})
    assert s.method.scc is True
    assert spec_from_form({"source": "preset", "preset": "H2O", "sk_set": "fake-1-0", "scc": ""}).method.scc is False
    o = spec_from_form({"source": "preset", "preset": "H2O", "code": "orca"}).method
    assert o.ts_calc_hess is True and o.ts_freq is True
    g = spec_from_form({"source": "preset", "preset": "H2O", "code": "gromacs"}).method
    assert g.gen_seed == 12345


# ---- finding 15: auto-named directories on Windows paths ----

def test_auto_named_directories_drop_a_trailing_backslash(tmp_path):
    app = WebApp(default_config(sk_root=""), tmp_path / "cluster.toml")
    app.form["output_dir"] = "C:\\runs\\x\\"
    assert app.scan_default_dir() == "C:\\runs\\x_scan" and app.stages_default_dir() == "C:\\runs\\x_stages"
    assert app.batch_default_dir() == "C:\\runs\\x_set"
    app.form["output_dir"] = "/runs/x/"
    assert app.scan_default_dir() == "/runs/x_scan"


# ---- finding 16: documented values must be finite ----

def test_doc_value_rejects_nan_and_inf(tmp_path):
    app = WebApp(default_config(sk_root=""), tmp_path / "cluster.toml")
    before = app.form.get("ecutwfc")
    assert app.put_doc_value("ecutwfc=nan") == "" and app.put_doc_value("ecutwfc=inf") == "" and app.form.get("ecutwfc") == before
    assert app.put_doc_value("ecutwfc=30") and app.form["ecutwfc"] == "30"


# ---- suspected items: /spec.json under the lock, continue_run on odd errors ----

def test_spec_json_waits_for_the_prepare_lock(web):
    app, base, _ = web
    done = []
    with app.lock:
        t = threading.Thread(target=lambda: done.append(_request(base, "GET", "/spec.json")[0]), daemon=True)
        t.start(); t.join(0.5)
        assert not done                 # blocked on app.lock while a preview would be rebuilding
    t.join(30)
    assert done == [404]


def test_continue_with_an_unexpected_error_keeps_the_typed_fields(web, tmp_path, monkeypatch):
    import adit.continuation as C

    app, base, _ = web
    monkeypatch.setattr(C, "continue_from", lambda *a, **k: (_ for _ in ()).throw(KeyError("velocities")))
    status, _, page = _request(base, "POST", "/continue", urllib.parse.urlencode({**app.form, "cont_dir": str(tmp_path), "job_name": "typed"}).encode())
    assert status == 200 and "内部エラー" not in page
    assert app.form["cont_dir"] == str(tmp_path) and app.form["job_name"] == "typed"
    assert "velocities" in h.unescape(page)


# ---- finding 18: unknown keys in a loaded spec.json are shown ----

def test_load_warns_about_unknown_keys_and_still_loads(web):
    app, base, _ = web
    d = json.loads(water_spec().to_json())
    d["task"]["force_tolerance_ev_per_angg"] = 0.01
    d["mystery"] = 1
    page = h.unescape(_load(base, json.dumps(d).encode()))
    assert app.spec is not None and app.form["code"] == "dftbplus"
    assert "task.force_tolerance_ev_per_angg" in page and "mystery" in page and ("無視" in page or "ignored" in page)
    page = h.unescape(_load(base, water_spec().to_json().encode()))
    assert "無視" not in page and "ignored" not in page


# ---- finding 19: measuring in the browser honours per-axis periodicity ----

def test_scene_json_carries_pbc():
    slab = fcc111("Al", (2, 2, 3), vacuum=10.0)
    assert list(slab.pbc) == [True, True, False]
    scene = scene_from_atoms(slab)
    assert scene.pbc == [True, True, False] and json.loads(scene.to_json())["pbc"] == [True, True, False]
    assert scene_from_atoms(bulk("Si")).pbc == [True, True, True]
    assert scene_from_atoms(molecule("H2O")).pbc is None and scene_from_atoms(molecule("H2O")).cell is None


def test_browser_measuring_does_not_wrap_across_the_vacuum(tmp_path):
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not installed")
    pos = [[0.0, 0.0, 1.0], [0.0, 0.0, 13.0], [0.5, 0.0, 1.0], [9.5, 0.0, 1.0]]
    cell = [[10.0, 0, 0], [0, 10.0, 0], [0, 0, 20.0]]
    script = tmp_path / "check.js"
    script.write_text(
        "var window = {}; var document = { getElementById: function () { return null; }, addEventListener: function () {} };\n"
        + JS.read_text(encoding="utf-8")
        + "\nvar G = window.ADIT_GEOMETRY, pos = %s, cell = %s, inv = G.inv3(cell), slab = [true, true, false];\n"
        "console.log(JSON.stringify([G.distance(pos, 0, 1, cell, inv, slab), G.distance(pos, 0, 1, cell, inv), G.distance(pos, 0, 1, cell, inv, null),"
        " G.distance(pos, 2, 3, cell, inv, slab), G.distance(pos, 0, 1, cell, inv, [false, false, false])]));\n" % (json.dumps(pos), json.dumps(cell)),
        encoding="utf-8")
    out = subprocess.run([node, str(script)], capture_output=True, text=True, timeout=30)
    assert out.returncode == 0, out.stderr
    got = json.loads(out.stdout.strip().splitlines()[-1])
    assert got[0] == pytest.approx(12.0)          # across the vacuum: no image along z
    assert got[1] == pytest.approx(8.0) and got[2] == pytest.approx(8.0)   # without pbc all axes wrap (old behaviour, still used by the players)
    assert got[3] == pytest.approx(1.0)           # the periodic axes still take the nearest image
    assert got[4] == pytest.approx(12.0)
