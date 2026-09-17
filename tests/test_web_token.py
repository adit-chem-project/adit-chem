from __future__ import annotations

import http.client
import threading
from urllib.parse import urlparse

import pytest

from adit.config import default_config
from adit.web.server import WebApp, main, serve, server_url
from tests.conftest import make_fake_skset


@pytest.fixture
def web(tmp_path):
    root = tmp_path / "slakos"
    make_fake_skset(root, "fake-1-0", ["H", "O"])
    app = WebApp(default_config(sk_root=str(root)), tmp_path / "cluster.toml")
    httpd = serve(app, "127.0.0.1", 0, token="secret-token")
    t = threading.Thread(target=httpd.serve_forever, daemon=True); t.start()
    yield f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown(); httpd.server_close()


def _status(url: str, body: bytes | None = None, headers: dict | None = None) -> tuple[int, dict]:
    u = urlparse(url)
    conn = http.client.HTTPConnection(u.hostname, u.port, timeout=60)
    try:
        conn.request("POST" if body is not None else "GET", u.path + ("?" + u.query if u.query else ""), body=body, headers=headers or {})
        r = conn.getresponse()
        r.read()
        return r.status, dict(r.headers)
    finally:
        conn.close()


def test_requests_without_token_are_refused(web):
    assert _status(web + "/")[0] == 403
    assert _status(web + "/preview", body=b"x=1")[0] == 403
    assert _status(web + "/?token=wrong")[0] == 403
    assert _status(web + "/", headers={"Cookie": "adit_token=wrong"})[0] == 403


def test_token_in_url_sets_cookie_and_strips_token(web):
    code, headers = _status(web + "/analysis?token=secret-token&dir=abc")
    assert code == 303
    assert headers["Location"] == "/analysis?dir=abc"
    assert "adit_token=secret-token" in headers["Set-Cookie"] and "HttpOnly" in headers["Set-Cookie"]
    assert _status(web + "/", headers={"Cookie": "adit_token=secret-token"})[0] == 200


def test_no_token_keeps_old_behaviour(tmp_path):
    app = WebApp(default_config(sk_root=str(tmp_path)), tmp_path / "cluster.toml")
    httpd = serve(app, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    try:
        assert _status(f"http://127.0.0.1:{httpd.server_port}/")[0] == 200
    finally:
        httpd.shutdown(); httpd.server_close()


def test_server_url():
    assert server_url("127.0.0.1", 8765) == "http://127.0.0.1:8765/"
    assert server_url("0.0.0.0", 80, "abc") == "http://0.0.0.0:80/?token=abc"
    assert server_url("::", 80) == "http://[::]:80/"


@pytest.mark.parametrize("argv, expect_token", [(["--host", "0.0.0.0"], True), ([], False),
                                                (["--host", "0.0.0.0", "--no-token"], False)])
def test_main_makes_a_token_off_localhost(monkeypatch, tmp_path, capsys, argv, expect_token):
    import adit.web.server as S

    seen = {}

    class Fake:
        server_port = 1

        def serve_forever(self):
            raise KeyboardInterrupt

        def server_close(self):
            pass

    def fake_serve(app, host, port, open_browser, token=None):
        seen["token"] = token
        return Fake()

    monkeypatch.setattr(S, "serve", fake_serve)
    assert main(argv + ["--config", str(tmp_path / "cluster.toml")]) == 0
    assert bool(seen["token"]) == expect_token
    assert ("token=" in capsys.readouterr().err) == expect_token
