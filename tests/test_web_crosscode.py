"""Web trial of native import, read-back verification and cross-run audit."""

import json
import threading
import urllib.parse
import urllib.request

import pytest

from adit import lang
from adit.config import default_config
from adit.web.server import WebApp, serve
from tests.test_native_import_cli import _vasp_bundle
from tests.test_native_verify import make_bundle


@pytest.fixture
def web_crosscode(tmp_path, monkeypatch):
    monkeypatch.setattr(lang, "LANGUAGE", "ja")
    app = WebApp(default_config(), tmp_path / "cluster.toml")
    httpd = serve(app, "127.0.0.1", 0)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown()
    httpd.server_close()


def _get(base, route):
    return urllib.request.urlopen(base + route).read().decode("utf-8")


def _post(base, route, fields):
    body = urllib.parse.urlencode(fields).encode("utf-8")
    return urllib.request.urlopen(base + route, data=body).read().decode("utf-8")


def test_web_native_import_writes_review_without_touching_source(web_crosscode, tmp_path):
    source, output = tmp_path / "source", tmp_path / "review"
    _vasp_bundle(source)
    before = {p.name: p.read_bytes() for p in source.iterdir()}
    page = _get(web_crosscode, "/convert")
    assert "既存の計算入力を読み込む" in page and 'action="/native_import"' in page

    page = _post(web_crosscode, "/native_import", {"native_source": str(source), "native_output": str(output)})
    assert "確認用の下書き" in page
    assert (output / "draft_spec.json").is_file()
    assert json.loads((output / "import_report.json").read_text(encoding="utf-8"))["complete"]
    assert before == {p.name: p.read_bytes() for p in source.iterdir()}
    again = _post(web_crosscode, "/native_import", {"native_source": str(source), "native_output": str(output)})
    assert "新しいディレクトリ" in again


def test_web_native_import_shows_unknown_and_no_draft(web_crosscode, tmp_path):
    source, output = tmp_path / "source", tmp_path / "review"
    _vasp_bundle(source, extra="UNKNOWN_KEY = 17\n")
    page = _post(web_crosscode, "/native_import", {"native_source": str(source), "native_output": str(output)})
    assert "UNKNOWN_KEY = 17" in page and "下書きの計算設定 (draft_spec.json) は作っていません" in page
    assert not (output / "draft_spec.json").exists()


def test_web_verify_distinguishes_mapped_pass_from_mismatch(web_crosscode, tmp_path):
    bundle = tmp_path / "generated"
    bundle.mkdir()
    make_bundle(bundle)
    before = {p.name: p.read_bytes() for p in bundle.iterdir()}
    page = _post(web_crosscode, "/native_verify", {"verify_project": str(bundle)})
    assert "対応項目のみの点検結果: 問題なし" in page and "未確認" in page
    assert before == {p.name: p.read_bytes() for p in bundle.iterdir()}

    incar = bundle / "INCAR"
    incar.write_text(incar.read_text(encoding="utf-8").replace("EDIFF = 1e-08", "EDIFF = 2e-08"), encoding="utf-8")
    page = _post(web_crosscode, "/native_verify", {"verify_project": str(bundle)})
    assert "対応項目のみの点検結果: 問題あり" in page and "method.ediff" in page


def test_web_cross_run_audit_is_read_only(web_crosscode, tmp_path, monkeypatch):
    from types import SimpleNamespace

    paths_seen = []
    def fake(paths, *, require_md):
        paths_seen.append((paths, require_md))
        return SimpleNamespace(comparable=False, summary_text=lambda: "異なる組成を検出")
    monkeypatch.setattr("adit.analysis.audit_run_dirs", fake)
    page = _get(web_crosscode, "/analysis")
    assert 'action="/audit_runs"' in page
    page = _post(web_crosscode, "/audit_runs", {"audit_ref": str(tmp_path / "a"), "audit_other": str(tmp_path / "b"), "audit_kind": "energy"})
    assert "異なる組成を検出" in page
    assert paths_seen == [([str(tmp_path / "a"), str(tmp_path / "b")], False)]
    assert not (tmp_path / "a").exists() and not (tmp_path / "b").exists()


def test_web_crosscode_english_labels(web_crosscode, tmp_path, monkeypatch):
    monkeypatch.setattr(lang, "LANGUAGE", "en")
    page = _get(web_crosscode, "/convert")
    assert "Import existing calculation input" in page
    assert "Verify generated input" in page
    page = _get(web_crosscode, "/analysis")
    assert "Cross-run audit" in page
