
from __future__ import annotations

import os
import re
import threading
import time
import urllib.parse
import urllib.request

import pytest

from adit import lang
from adit.config import Profile, default_config, load_config, save_config, set_top_level_value
from adit.project import _local_time, build_project
from adit.spec import OrcaMethod, VaspMethod, XtbMethod
from adit.validate import validate
from adit.web.server import WebApp, serve
from tests.conftest import cfg_for, make_fake_skset, water_spec


@pytest.fixture
def ja():
    before = lang.LANGUAGE
    lang.set_language("ja")
    yield
    lang.set_language(before)


def _start(app):
    httpd = serve(app, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd, f"http://127.0.0.1:{httpd.server_port}"


def _post(url, fields):
    return urllib.request.urlopen(url, data=urllib.parse.urlencode(fields).encode()).read().decode()


def test_web_reloads_config_and_language_keeps_other_lines(sk_root, tmp_path, ja):
    path = tmp_path / "cluster.toml"
    save_config(default_config(), path)
    app = WebApp(load_config(path), path)
    httpd, base = _start(app)
    try:
        assert "fake-1-0" not in urllib.request.urlopen(base + "/").read().decode()
        text = path.read_text(encoding="utf-8").replace('sk_root = ""', f'sk_root = "{sk_root.as_posix()}"  # 自分で書いた行') + "# 最後のコメント\n"
        time.sleep(0.02)
        path.write_text(text, encoding="utf-8")
        assert "fake-1-0" in urllib.request.urlopen(base + "/").read().decode()
        _post(base + "/language", {})
        after = path.read_text(encoding="utf-8")
        assert f'sk_root = "{sk_root.as_posix()}"  # 自分で書いた行' in after and "# 最後のコメント" in after
        assert 'language = "en"' in after and lang.LANGUAGE == "en"
        _post(base + "/language", {})
    finally:
        httpd.shutdown(); httpd.server_close()


def test_web_shows_broken_config_and_keeps_running(tmp_path, ja):
    path = tmp_path / "cluster.toml"
    save_config(default_config(), path)
    app = WebApp(load_config(path), path)
    httpd, base = _start(app)
    try:
        time.sleep(0.02)
        path.write_text('sk_root = "C:\\Users\\me\\slakos"\n', encoding="utf-8")
        html = urllib.request.urlopen(base + "/").read().decode()
        assert "環境設定ファイルを読めないため" in html and "/mnt/c/Users/" in html
    finally:
        httpd.shutdown(); httpd.server_close()


def test_set_top_level_value_inserts_before_tables(tmp_path):
    path = tmp_path / "c.toml"
    path.write_text('sk_root = "/x"\n\n[profiles.local]\nkind = "direct"\n', encoding="utf-8")
    set_top_level_value(path, "language", "en")
    text = path.read_text(encoding="utf-8")
    assert text.index('language = "en"') < text.index("[profiles.local]") and load_config(path).language == "en"


def test_the_generated_directory_is_cleared_on_code_change(sk_root, tmp_path, ja, monkeypatch):
    app = WebApp(cfg_for(sk_root), tmp_path / "none.toml")
    f = dict(app.form); f.update(sk_set="fake-1-0", output_dir=str(tmp_path / "out"))
    app.preview(f); msg, err = app.generate(False)
    assert not err and app.written is not None
    app.preview({**f, "code": "xtb"})
    assert app.written is None      # changing the code clears "written" until regenerated


def test_analysis_page_not_run_notice_and_md_defaults(sk_root, tmp_path, ja):
    from tests.conftest import Task
    from adit.project import write_project
    cfg = cfg_for(sk_root)
    out = tmp_path / "md"
    write_project(water_spec(task=Task(type="molecular_dynamics")), cfg, out)
    app = WebApp(cfg, tmp_path / "none.toml")
    httpd, base = _start(app)
    try:
        html = urllib.request.urlopen(base + "/analysis?dir=" + urllib.parse.quote(str(out))).read().decode()
        assert "解析はまだ実行していません" in html
        assert html.index("まだ実行していません") < html.index("<fieldset>")
        assert re.search(r'name="rdf" checked', html) and re.search(r'name="msd" checked', html)
        assert not (out / "analysis").exists()      # GET runs nothing
        html = urllib.request.urlopen(base + "/analysis", data=urllib.parse.urlencode({"run_dir": str(out), "action": "run"}).encode()).read().decode()
        assert "まだ実行していないようです" in html and (out / "analysis").is_dir()
    finally:
        httpd.shutdown(); httpd.server_close()


def test_placeholder_in_profile_is_an_error(sk_root, ja):
    cfg = cfg_for(sk_root)
    cfg.profiles["cluster"].env["VASP_PP_PATH"] = "<POTCAR ライブラリの親ディレクトリ>"
    from adit.spec import Runtime
    errs = validate(water_spec(runtime=Runtime(profile="cluster")), cfg)
    assert any("仮の値のままです" in str(e) and "VASP_PP_PATH" in str(e) for e in errs)
    cfg.profiles["cluster"].env["VASP_PP_PATH"] = "/opt/potcar"
    cfg.profiles["cluster"].commands["dftbplus"] = "dftb+ < in > out"
    assert not any("仮の値" in str(e) for e in validate(water_spec(runtime=Runtime(profile="cluster")), cfg))


def test_missing_elements_list_other_sets(sk_root, ja):
    make_fake_skset(sk_root, "withsi-1-0", ["H", "O", "Si"])
    spec = water_spec()
    spec.structure.atoms.symbols[0] = "Si"  # SiH2
    from adit.spec import DftbMethod
    spec.method = DftbMethod(sk_set="fake-1-0")
    msgs = [str(e) for e in validate(spec, cfg_for(sk_root))]
    assert any("無い元素があります" in m and "withsi-1-0" in m for m in msgs), msgs
    spec.structure.atoms.symbols[0] = "Al"
    msgs = [str(e) for e in validate(spec, cfg_for(sk_root))]
    assert any("xtb (GFN-xTB)" in m and "dftb.org" in m for m in msgs), msgs


def test_sk_root_unset_message_works_on_every_screen(ja):
    msgs = " ".join(str(e) for e in validate(water_spec(), cfg_for(None)))
    assert "sk_root" in msgs and "ウェブ版とコマンド行" in msgs and "cluster.toml" in msgs


# ---- 12, 13: README.txt ----
@pytest.mark.parametrize("method,word", [(VaspMethod(), "export PATH=<VASP"), (OrcaMethod(), "export PATH=<ORCA"), (XtbMethod(), "conda activate")])
def test_readme_path_step_per_code(sk_root, ja, method, word):
    from adit.spec import Runtime
    spec = water_spec(method=method, runtime=Runtime(profile="local"))
    if isinstance(method, VaspMethod):
        from tests.test_vasp import h2o_spec
        spec = h2o_spec(runtime=Runtime(profile="local"))
    r = build_project(spec, cfg_for(sk_root)).texts["README.txt"]
    local = r.split("== この PC で実行する ==", 1)[1].split("== 研究室", 1)[0]
    assert word in local
    if not isinstance(method, XtbMethod):
        assert "conda activate" not in local


def test_local_time():
    assert _local_time("2026-09-11T09:10:51+00:00") != "2026-09-11T09:10:51+00:00"
    assert "(" in _local_time("2026-09-11T09:10:51+00:00")
    assert _local_time("broken") == "broken"
