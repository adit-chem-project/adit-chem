
from __future__ import annotations

import html as h
import json
import re
import shutil
import threading
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

import pytest

from adit import lang
from adit.config import default_config
from adit.web.server import WebApp, serve
from tests.conftest import make_fake_skset

REPO = Path(__file__).resolve().parent.parent


@pytest.fixture
def ja():
    before = lang.LANGUAGE
    lang.set_language("ja")
    yield
    lang.set_language(before)


@pytest.fixture
def web(tmp_path, ja):
    root = tmp_path / "slakos"
    make_fake_skset(root, "fake-1-0", ["H", "C", "N", "O", "Si", "Al"])
    app = WebApp(default_config(sk_root=str(root)), tmp_path / "cluster.toml")
    app.form["output_dir"] = str(tmp_path / "run")
    httpd = serve(app, "127.0.0.1", 0)
    t = threading.Thread(target=httpd.serve_forever, daemon=True); t.start()
    yield app, f"http://127.0.0.1:{httpd.server_port}"
    httpd.shutdown(); httpd.server_close()


def _get(url: str) -> str:
    return urllib.request.urlopen(url, timeout=60).read().decode()


def _post(url: str, fields: dict) -> str:
    return urllib.request.urlopen(url, data=urllib.parse.urlencode(fields).encode(), timeout=60).read().decode()


def _error(page: str) -> str:
    m = re.search(r'<div class="error">(.*?)</div>', page, re.S)
    return h.unescape(m.group(1)) if m else ""


def _banner(page: str, cls: str) -> str:
    m = re.search(rf'<div class="{cls}" role="\w+">(.*?)</div>', page, re.S)
    return h.unescape(m.group(1)) if m else ""


def _selected_item(page: str) -> str:
    sel = re.search(r'<select name="scan_item" id="scan_item".*?</select>', page, re.S).group(0)
    return re.search(r'<option value="([^"]+)"[^>]*selected', sel).group(1)


def _water(app, tmp_path, **kw) -> dict:
    f = dict(app.form); f.update(sk_set="fake-1-0", output_dir=str(tmp_path / "run"))
    f.update(kw)
    return f


def test_choices_match_desktop():
    pytest.importorskip("PySide6")
    from adit.gui.scan_dialog import choices_for as desktop
    from adit.web.scan_choices import choices_for as web_
    before = lang.LANGUAGE
    try:
        for language in ("ja", "en"):
            lang.set_language(language)
            for code in ("dftbplus", "vasp", "xtb", "espresso", "orca"):
                for periodic in (False, True):
                    a = [(c.path, c.label, c.example, c.purpose) for c in desktop(code, periodic)]
                    b = [(c.path, c.label, c.example, c.purpose) for c in web_(code, periodic)]
                    assert a == b, (language, code, periodic)
    finally:
        lang.set_language(before)


def test_all_choices_carry_the_conditions():
    from adit.web.scan_choices import OTHER, all_choices
    got = {c.path: (code, per) for c, code, per in all_choices()}
    assert got == {"method.ecutwfc": ("espresso", False), "method.encut": ("vasp", False), "kpoints.mesh": ("", True),
                   "kpoints.density": ("", True), "scale": ("", True), OTHER: ("", False)}


def test_scan_box_is_on_the_page_with_css_rules(web):
    app, base = web
    page = _get(base + "/")
    assert '<details class="scan"' in page and "1 つの条件を変えて一括生成" in page
    for label in ("変える項目", "項目の場所", "値", "保存先"):
        assert f">{label}</label>" in page, label
    assert 'option value="method.ecutwfc" data-code="espresso"' in page and 'option value="method.encut" data-code="vasp"' in page
    assert 'option value="scale" class="per"' in page
    assert '#main:has(#code option[value="espresso"]:checked) #scan_item option[data-code="espresso"] { display:block; }' in page
    assert '#main:has(#source option[value="bulk"]:checked) #scan_item option.per' in page
    assert '#main:has(#scan_item option[value="other"]:checked) .scan-path-row { display:grid; }' in page
    assert 'formaction="/scan"' in page
    assert _selected_item(page) == "other"
    assert f'name="scan_dir" id="scan_dir" value="{app.form["output_dir"]}_scan"' in page


def test_scan_fields_have_labels(web):
    from tests.test_web import _fields_without_label
    app, base = web
    assert _fields_without_label(_get(base + "/")) == []


def test_default_item_follows_code_and_periodicity(web, tmp_path):
    app, base = web
    page = _post(base + "/preview", _water(app, tmp_path, source="bulk", bulk_el="Si", kp_mode="mesh", k1="2", k2="2", k3="2"))
    assert not _error(page)
    assert _selected_item(page) == "kpoints.mesh"
    page = _post(base + "/preview", _water(app, tmp_path, source="preset"))
    assert _selected_item(page) == "other"


def test_scan_dir_follows_output_dir_until_edited(web, tmp_path):
    app, base = web
    _get(base + "/")
    auto = app.form["output_dir"] + "_scan"
    page = _post(base + "/preview", _water(app, tmp_path, output_dir=str(tmp_path / "other"), scan_dir=auto))
    assert f'value="{tmp_path / "other"}_scan"' in page
    page = _post(base + "/preview", _water(app, tmp_path, scan_dir=str(tmp_path / "mine")))
    assert f'name="scan_dir" id="scan_dir" value="{tmp_path / "mine"}"' in page


def test_scan_generates_directories(web, tmp_path):
    app, base = web
    out = tmp_path / "out_scan"
    page = _post(base + "/scan", _water(app, tmp_path, scan_item="other", scan_path="method.max_scc_iterations",
                                        scan_values="50, 100", scan_dir=str(out)))
    assert not _error(page), _error(page)
    msg = _banner(page, "message-banner")
    assert f"{out} の中に 2 個のディレクトリを作りました。" in msg and "max_scc_iterations_50" in msg and "次にすること" in msg
    assert f'href="/analysis?dir={urllib.parse.quote(str(out))}"' in page
    assert (out / "scan.json").is_file() and (out / "max_scc_iterations_100" / "dftb_in.hsd").is_file()
    assert "MaxSccIterations = 100" in (out / "max_scc_iterations_100" / "dftb_in.hsd").read_text(encoding="utf-8")
    page = _post(base + "/scan", _water(app, tmp_path, scan_item="other", scan_path="method.max_scc_iterations",
                                        scan_values="50,100", scan_dir=str(out)))
    assert "ファイルがあります" in _error(page) and "上書きを許可" in _error(page)
    page = _post(base + "/scan", _water(app, tmp_path, scan_item="other", scan_path="method.max_scc_iterations",
                                        scan_values="50,100", scan_dir=str(out), scan_overwrite="on"))
    assert not _error(page) and _banner(page, "message-banner")


def test_scan_kpoints_mesh_on_bulk(web, tmp_path):
    app, base = web
    out = tmp_path / "si_scan"
    page = _post(base + "/scan", _water(app, tmp_path, source="bulk", bulk_el="Si", scan_item="kpoints.mesh",
                                        scan_values="2x2x2,4x4x4", scan_dir=str(out)))
    assert not _error(page), _error(page)
    assert sorted(p.name for p in out.iterdir() if p.is_dir()) == ["mesh_2x2x2", "mesh_4x4x4"]


@pytest.mark.parametrize("kw, expect", [
    (dict(scan_item="other", scan_path="method.max_scc_iterations", scan_values="50"), "値はカンマで区切って 2 つ以上"),
    (dict(scan_item="other", scan_path="", scan_values="1,2"), "項目の場所を書いてください"),
    (dict(scan_item="other", scan_path="method.max_scc_iterations", scan_values="1,2", scan_dir=""), "保存先を指定してください"),
    (dict(scan_item="scale", scan_values="0.99,1.01"), "周期系だけ"),
    (dict(scan_item="other", scan_path="method.no_such", scan_values="1,2"), "という項目はありません"),
    (dict(scan_item="other", scan_path="method.max_scc_iterations", scan_values="a,b"), "として読めません"),
    (dict(scan_item="other", scan_path="method.max_scc_iterations", scan_values="50,0"), "生成できません"),
    (dict(scan_item="nonsense", scan_values="1,2"), "変える項目を選んでください"),
    (dict(scan_item="method.ecutwfc", scan_values="30,40"), "という項目はありません"),
])
def test_scan_bad_inputs_show_reason(web, tmp_path, kw, expect):
    app, base = web
    out = tmp_path / "bad_scan"
    fields = _water(app, tmp_path, scan_dir=str(out)); fields.update(kw)
    page = _post(base + "/scan", fields)
    err = _error(page)
    assert expect in err or expect in _banner(page, "error-banner"), (err, _banner(page, "error-banner"))
    assert _banner(page, "error-banner").startswith("生成できません: ")
    assert not _banner(page, "message-banner")
    box = re.search(r'<details class="scan".*?</details>', page, re.S).group(0)
    assert ' open>' in box.split("\n", 1)[0] and page.count('<div class="error">') == 1 and '<div class="error">' in box
    assert not out.exists() or not any(out.glob("*/dftb_in.hsd"))


def test_scan_with_broken_form_shows_preview_error(web, tmp_path):
    app, base = web
    page = _post(base + "/scan", _water(app, tmp_path, charge="abc", scan_item="other", scan_path="method.max_scc_iterations",
                                        scan_values="1,2", scan_dir=str(tmp_path / "x")))
    assert _error(page) and not (tmp_path / "x").exists()


def test_scan_to_a_file_path_is_refused(web, tmp_path):
    app, base = web
    f = tmp_path / "afile"; f.write_text("x", encoding="utf-8")
    page = _post(base + "/scan", _water(app, tmp_path, scan_item="other", scan_path="method.max_scc_iterations",
                                        scan_values="1,2", scan_dir=str(f)))
    assert "保存先にできません" in _error(page)


def _scan_copy(tmp_path: Path, path: str, values: list[str]) -> Path:
    out = tmp_path / "scan"
    dirs = [f"{path.split('.')[-1]}_{v}" for v in values]
    for d in dirs:
        shutil.copytree(REPO / "examples" / "dftb_tio2_generated", out / d)
    (out / "scan.json").write_text(json.dumps({"path": path, "values": values, "dirs": dirs}), encoding="utf-8")
    return out


def _table(page: str) -> tuple[list[str], list[list[str]]]:
    t = re.search(r'<table class="scan">(.*?)</table>', page, re.S).group(1)
    head = [h.unescape(x) for x in re.findall(r"<th(?:\s[^>]*)?>(.*?)</th>", t)]
    body = [[h.unescape(x) for x in re.findall(r"<td[^>]*>(.*?)</td>", r)] for r in re.findall(r"<tr>(.*?)</tr>", t.split("<tbody>")[1], re.S)]
    return head, body


def test_analysis_page_shows_scan_table_and_figure(web, tmp_path):
    from adit.analysis.report import figure_title
    app, base = web
    out = _scan_copy(tmp_path, "method.x", ["1", "2"])
    page = _post(base + "/analysis", {"run_dir": str(out), "rdf": "on", "rmax": "abc"})
    assert not _error(page), _error(page)
    head, body = _table(page)
    assert head == ["値", "エネルギー [eV]", "最後の値との差 [meV/原子]", "力の最大値 [eV/Å]", "圧力 [GPa]", "注"]
    assert len(body) == 2 and body[0][0] == "1" and body[0][2] == "0.000" and body[0][4] != "-"
    assert '<td class="num">' in page and '<th scope="col" class="num">' in page and '<th scope="col">値</th>' in page
    assert "1 つの条件だけを変えた一連の計算です" in page and 'name="rdf" disabled' in page.replace(" checked", "")
    summary = h.unescape(re.search(r"<legend>要約</legend><pre>(.*?)</pre>", page, re.S).group(1))
    assert "変えた項目: method.x" in summary and "scan_energies.csv" in summary and "\t" not in summary
    assert f"<legend>{figure_title('scan_energy')}</legend>" in page
    figs = re.findall(r'src="/file\?path=([^"]+)"', page)
    assert len(figs) == 1 and _get_bytes(base + "/file?path=" + figs[0])[:4] == b"\x89PNG"
    assert "計算結果が見つかりません" not in page


def _get_bytes(url: str) -> bytes:
    return urllib.request.urlopen(url, timeout=60).read()


def test_analysis_get_link_and_scale_reference(web, tmp_path):
    app, base = web
    out = _scan_copy(tmp_path, "scale", ["0.99", "1.01"])
    page = _get(base + "/analysis?dir=" + urllib.parse.quote(str(out)))   # GET only fills the form
    assert "解析はまだ実行していません" in page and '<table class="scan">' not in page
    assert f'name="run_dir" value="{out}"' in page and 'name="rdf" disabled' in page
    page = _post(base + "/analysis", {"run_dir": str(out)})
    head, body = _table(page)
    assert head[2] == "最小値との差 [meV/原子]" and len(body) == 2
    assert "5 つ以上要ります" in h.unescape(page)


def test_analysis_scan_table_in_english(web, tmp_path):
    app, base = web
    out = _scan_copy(tmp_path, "method.x", ["1", "2"])
    lang.set_language("en")
    page = _post(base + "/analysis", {"run_dir": str(out)})
    head, _ = _table(page)
    assert head == ["Value", "Energy [eV]", "Diff. from last [meV/atom]", "Max force [eV/Å]", "Pressure [GPa]", "Note"]
    assert "These runs differ in one setting only." in page and "Results per value" in page


def test_analysis_broken_scan_json_shows_reason(web, tmp_path):
    app, base = web
    out = tmp_path / "broken"; out.mkdir(); (out / "scan.json").write_text("{not json", encoding="utf-8")
    page = _post(base + "/analysis", {"run_dir": str(out)})
    assert _error(page).startswith("解析できません: ") and '<table class="scan">' not in page


def test_analysis_normal_dir_unchanged(web, tmp_path):
    app, base = web
    run_dir = tmp_path / "md"
    shutil.copytree(REPO / "examples" / "dftb_md_water_generated", run_dir)
    page = _post(base + "/analysis", {"run_dir": str(run_dir), "rdf": "on"})
    assert not _error(page) and '<table class="scan">' not in page and "1 つの条件だけを変えた一連の計算です" not in page
    assert 'name="rdf" checked>' in page
