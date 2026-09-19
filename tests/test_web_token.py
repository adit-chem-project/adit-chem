from __future__ import annotations

import http.client
import socket
import threading
from urllib.parse import quote, urlparse

import pytest

from adit.config import default_config
from adit.web.server import WebApp, cookie_name, main, redact_token, serve, server_lines, server_url
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


def _cookie(url: str) -> str:
    return cookie_name(urlparse(url).port)


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
    assert _status(web + "/", headers={"Cookie": f"{_cookie(web)}=wrong"})[0] == 403


def test_token_in_url_sets_cookie_and_strips_token(web):
    code, headers = _status(web + "/analysis?token=secret-token&dir=abc")
    assert code == 303
    assert headers["Location"] == "/analysis?dir=abc"
    assert f"{_cookie(web)}=secret-token" in headers["Set-Cookie"] and "HttpOnly" in headers["Set-Cookie"]
    assert "SameSite=Lax" in headers["Set-Cookie"]
    assert _status(web + "/", headers={"Cookie": f"{_cookie(web)}=secret-token"})[0] == 200


def test_post_with_token_in_url_is_handled_not_redirected(web):
    code, headers = _status(web + "/language?token=secret-token", body=b"")
    assert code == 303 and headers["Location"] == "/" and f"{_cookie(web)}=secret-token" in headers["Set-Cookie"]
    code, headers = _status(web + "/preview?token=secret-token", body=b"source=preset&preset=H2O&sk_set=fake-1-0",
                            headers={"Content-Type": "application/x-www-form-urlencoded"})
    assert code == 200 and f"{_cookie(web)}=secret-token" in headers["Set-Cookie"]


def test_non_ascii_token_is_refused_not_500(web):
    assert _status(web + "/?token=%E6%97%A5")[0] == 403
    assert _status(web + "/", headers={"Cookie": f"{_cookie(web)}=%E6%97%A5"})[0] == 403


def test_cookie_name_carries_the_port():
    assert cookie_name(8765) == "adit_token_8765"


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
    assert server_url("::", 80, "a b&c") == "http://[::]:80/?token=" + quote("a b&c", safe="")


def test_server_lines_for_wildcard_hosts(monkeypatch):
    monkeypatch.setattr(socket, "getfqdn", lambda: "pc.example.org")
    lines = server_lines("0.0.0.0", 8765, "abc")
    assert "http://127.0.0.1:8765/?token=abc" in lines[1] and "http://pc.example.org:8765/?token=abc" in lines[2]
    assert "IP" in lines[3] and len(lines) == 5
    one = server_lines("127.0.0.1", 8765)
    assert len(one) == 1 and one[0].endswith("http://127.0.0.1:8765/  (Ctrl-C で停止 / Ctrl-C to stop)")


def test_redact_token():
    assert redact_token("GET /?token=abc&dir=x HTTP/1.1") == "GET /?token=***&dir=x HTTP/1.1"
    assert redact_token("GET /analysis?dir=x HTTP/1.1") == "GET /analysis?dir=x HTTP/1.1"


def test_log_hides_the_token(web, monkeypatch, capsys):
    monkeypatch.setenv("ADIT_WEB_LOG", "1")
    assert _status(web + "/?token=secret-token")[0] == 303
    err = capsys.readouterr().err
    assert "token=***" in err and "secret-token" not in err


def test_ipv6_loopback_serves(tmp_path):
    app = WebApp(default_config(sk_root=str(tmp_path)), tmp_path / "cluster.toml")
    try:
        httpd = serve(app, "::1", 0)
    except OSError as ex:
        pytest.skip(f"no IPv6 loopback: {ex}")
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    try:
        conn = http.client.HTTPConnection("::1", httpd.server_port, timeout=60)
        conn.request("GET", "/"); r = conn.getresponse(); r.read()
        assert r.status == 200
        conn.close()
    finally:
        httpd.shutdown(); httpd.server_close()


def test_main_rejects_odd_token_and_busy_port(monkeypatch, tmp_path, capsys):
    import adit.web.server as S

    assert main(["--token", "日本語", "--config", str(tmp_path / "cluster.toml")]) == 2
    assert "--token" in capsys.readouterr().err

    def busy(app, host, port, open_browser, token=None):
        raise OSError(98, "Address already in use")

    monkeypatch.setattr(S, "serve", busy)
    assert main(["--port", "1", "--config", str(tmp_path / "cluster.toml")]) == 2
    assert "127.0.0.1:1" in capsys.readouterr().err


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
    monkeypatch.setattr(socket, "getfqdn", lambda: "pc.example.org")
    assert main(argv + ["--config", str(tmp_path / "cluster.toml")]) == 0
    assert bool(seen["token"]) == expect_token
    assert ("token=" in capsys.readouterr().err) == expect_token
