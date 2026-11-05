import threading
import urllib.request

import pytest

from adit.config import default_config, save_config
from adit.web.server import WebApp, serve


@pytest.mark.parametrize("theme,media,dark_rule", [("auto", True, True), ("light", False, False), ("dark", False, True)])
def test_the_web_pages_follow_the_theme_setting(tmp_path, monkeypatch, theme, media, dark_rule):
    monkeypatch.delenv("ADIT_THEME", raising=False)
    cfg = default_config(sk_root="")
    cfg.theme = theme
    path = tmp_path / "cluster.toml"
    save_config(cfg, path)
    httpd = serve(WebApp(cfg, path), "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    try:
        html = urllib.request.urlopen(f"http://127.0.0.1:{httpd.server_port}/convert").read().decode()
    finally:
        httpd.shutdown()
        httpd.server_close()
    style = html.split("</style>")[0]
    assert ("@media (prefers-color-scheme: dark)" in style) is media
    assert ("--bg-top:#1c1c1e" in style) is dark_rule
    assert f'var theme = "{theme}"' in html
