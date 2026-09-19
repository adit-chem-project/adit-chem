"""Fetching structures from databases: every request is answered from tests/data/fetch (no network)."""

from __future__ import annotations

import io
import json
import urllib.error
import urllib.request
from pathlib import Path

import pytest

from adit import fetch as F
from adit.fetch import FetchError

DATA = Path(__file__).resolve().parent / "data" / "fetch"

# URL -> fixture file, in the shape each service documents (responses recorded on 2026-09-19, see tests/data/SOURCES.md).
ROUTES = {
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/water/cids/JSON": "pubchem_water_cids.json",
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/962/SDF?record_type=3d": "pubchem_962_3d.sdf",
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/962/property/MolecularFormula,Title/JSON": "pubchem_962_props.json",
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/962,963/property/MolecularFormula,IUPACName,Title/JSON": "pubchem_962_props.json",
    "https://www.crystallography.net/cod/1000041.cif": "cod_1000041.cif",
    "https://www.crystallography.net/cod/result?id=1000041&format=json": "cod_1000041_meta.json",
    "https://api.materialsproject.org/materials/summary/?material_ids=mp-149&_fields=material_id%2Cformula_pretty%2Cstructure%2Csymmetry%2Cdeprecated":
        "mp_summary_mp-149.json",
    "https://providers.optimade.org/providers.json": "providers.json",
    "https://providers.optimade.org/index-metadbs/oqmd/v1/links": "idx_oqmd_links.json",
    "https://oqmd.org/optimade/v1/structures?filter=chemical_formula_reduced%3D%22Si%22&page_limit=5"
    "&response_fields=chemical_formula_descriptive%2Cchemical_formula_reduced%2Cnsites": "oqmd_si.json",
    "https://oqmd.org/optimade/v1/structures/4061352": "oqmd_4061352.json",
}


class _Resp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.close()


@pytest.fixture
def routes(monkeypatch):
    seen: list[tuple[str, dict]] = []
    extra: dict[str, bytes] = {}

    def fake_urlopen(req, timeout=None):
        url = req.full_url if isinstance(req, urllib.request.Request) else req
        headers = dict(req.header_items()) if isinstance(req, urllib.request.Request) else {}
        seen.append((url, headers))
        if url in extra:
            return _Resp(extra[url])
        name = ROUTES.get(url)
        if name is None:
            raise urllib.error.HTTPError(url, 404, "Not Found", {}, io.BytesIO(b""))
        return _Resp((DATA / name).read_bytes())

    monkeypatch.setattr(urllib.request, "urlopen", fake_urlopen)
    return seen, extra


def test_pubchem_by_name(routes, tmp_path):
    seen, _ = routes
    res = F.fetch("pubchem:water")
    assert res.fetched is not None and not res.candidates
    a = res.fetched.atoms
    assert a.get_chemical_formula() == "H2O" and not any(a.pbc)
    rec = res.fetched.record
    assert rec["database"] == "pubchem" and rec["id"] == "962" and rec["query"] == "water" and rec["format"] == "sdf"
    assert rec["title"] == "Water H2O" and rec["url"].endswith("/compound/cid/962/SDF?record_type=3d")
    assert len(rec["sha256"]) == 64 and rec["bytes"] == (DATA / "pubchem_962_3d.sdf").stat().st_size
    assert "public domain" in rec["license"] and rec["license_url"].startswith("https://www.ncbi.nlm.nih.gov/")
    assert rec["fetched_utc"].endswith("+00:00")
    assert all(h.get("User-agent", "").startswith("ADIT/") for _, h in seen)
    path = res.fetched.save(tmp_path)
    assert path == tmp_path / "pubchem_962.sdf" and path.read_bytes() == (DATA / "pubchem_962_3d.sdf").read_bytes()
    assert rec["file"] == str(path) and rec["response_file"] == str(path)


def test_pubchem_by_cid_and_missing_3d(routes):
    seen, extra = routes
    res = F.fetch("pubchem:962")
    assert res.fetched.record["id"] == "962" and seen[0][0].endswith("/cid/962/SDF?record_type=3d")
    with pytest.raises(FetchError, match="取得できません: PubChem の CID 12345 には 3D 座標がありません"):
        F.fetch("pubchem:12345")
    with pytest.raises(FetchError, match="取得できません: PubChem に 'nosuchmolecule' という名前の化合物がありません"):
        F.fetch("pubchem:nosuchmolecule")


def test_pubchem_several_cids_lists_candidates(routes):
    _, extra = routes
    extra["https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/water/cids/JSON"] = json.dumps({"IdentifierList": {"CID": [962, 963]}}).encode()
    res = F.fetch("pubchem:water")
    assert res.fetched is None and [c.ref for c in res.candidates] == ["pubchem:962"]
    assert res.candidates[0].label == "CID 962: Water H2O"
    text = F.candidates_text(res)
    assert "候補が 1 件" in text and "pubchem:962" in text


def test_cod(routes, tmp_path):
    res = F.fetch("cod:1000041")
    a = res.fetched.atoms
    assert a.get_chemical_formula() == "Cl4Na4" and all(a.pbc) and abs(a.cell.lengths()[0] - 5.62) < 1e-6
    rec = res.fetched.record
    assert rec["database"] == "cod" and rec["id"] == "1000041" and rec["format"] == "cif"
    assert rec["license"].startswith("CC0 1.0") and rec["license_url"] == "https://creativecommons.org/publicdomain/zero/1.0/"
    assert rec["title"].startswith("Sodium chloride / Accuracy of an automatic diffractometer")
    assert rec["doi"] == "10.1107/S0365110X65002244" and rec["authors"].startswith("Abrahams")
    assert res.fetched.save(tmp_path).name == "cod_1000041.cif"
    with pytest.raises(FetchError, match="取得できません: COD に ID 7777777 の項目がありません"):
        F.fetch("cod:7777777")
    with pytest.raises(FetchError, match="COD の ID は数字です"):
        F.fetch("cod:NaCl")


def test_materials_project(routes, tmp_path):
    seen, _ = routes
    with pytest.raises(FetchError, match="取得できません: Materials Project の API キーが環境設定"):
        F.fetch("mp:mp-149")
    assert not seen
    res = F.fetch("mp:mp-149", mp_api_key="KEY123")
    assert seen[0][1]["X-api-key"] == "KEY123"
    a = res.fetched.atoms
    assert a.get_chemical_formula() == "Si2" and all(a.pbc) and abs(a.cell.volume - 40.33) < 0.05
    rec = res.fetched.record
    assert rec["database"] == "mp" and rec["id"] == "mp-149" and rec["title"] == "Si Fd-3m" and rec["deprecated"] is False
    assert rec["license"].startswith("CC BY 4.0") and rec["license_url"] == "https://creativecommons.org/licenses/by/4.0/"
    path = res.fetched.save(tmp_path)
    assert path == tmp_path / "mp-149.extxyz" and (tmp_path / "mp-149.json").is_file()
    from ase.io import read
    assert read(path).get_chemical_formula() == "Si2"
    assert F.fetch("mp:149", mp_api_key="k").fetched.record["id"] == "mp-149"


def test_materials_project_rejects_key_and_partial_occupancy(routes):
    _, extra = routes

    def forbidden(req, timeout=None):
        raise urllib.error.HTTPError(req.full_url, 403, "Forbidden", {}, io.BytesIO(b""))

    import urllib.request as ur
    real = ur.urlopen
    ur.urlopen = forbidden
    try:
        with pytest.raises(FetchError, match="API キーを受け付けませんでした"):
            F.fetch("mp:mp-149", mp_api_key="bad")
    finally:
        ur.urlopen = real
    doc = json.loads((DATA / "mp_summary_mp-149.json").read_text())
    doc["data"][0]["structure"]["sites"][0]["species"] = [{"element": "Si", "occu": 0.5}, {"element": "Ge", "occu": 0.5}]
    extra[next(u for u in ROUTES if "materialsproject" in u)] = json.dumps(doc).encode()
    with pytest.raises(FetchError, match="部分占有"):
        F.fetch("mp:mp-149", mp_api_key="k")


def test_optimade_search_then_fetch(routes, tmp_path):
    seen, _ = routes
    res = F.fetch("optimade:Si")
    assert res.fetched is None and len(res.candidates) == 2
    assert res.candidates[0].ref == "optimade:oqmd:4061352" and "oqmd:4061352" in res.candidates[0].label and "1 原子" in res.candidates[0].label
    urls = [u for u, _ in seen]
    assert urls[0] == "https://providers.optimade.org/providers.json"
    assert not any("exmpl" in u or "aiida" in u for u in urls)
    assert res.notes == []
    res2 = F.fetch(res.candidates[0].ref)
    a = res2.fetched.atoms
    assert a.get_chemical_formula() == "Si" and all(a.pbc) and abs(a.cell.volume - 14.43) < 0.05
    rec = res2.fetched.record
    assert rec["database"] == "optimade" and rec["provider"] == "oqmd" and rec["id"] == "4061352" and rec["query"] == "oqmd:4061352"
    assert "OQMD" in rec["license"] and rec["license_url"] == "https://oqmd.org"
    assert rec["title"] == "Si F m -3 m"
    assert res2.fetched.save(tmp_path).name == "optimade_oqmd_4061352.extxyz" and (tmp_path / "optimade_oqmd_4061352.json").is_file()


def test_optimade_formula_and_filter():
    assert F.reduced_formula("SiO2") == "O2Si" and F.reduced_formula("Na2Cl2") == "ClNa" and F.reduced_formula("H2O") == "H2O"
    assert F.optimade_filter("SiO2") == 'chemical_formula_reduced="O2Si"'
    assert F.optimade_filter('elements HAS ALL "Si","O"') == 'elements HAS ALL "Si","O"'
    with pytest.raises(FetchError, match="組成式として読めません"):
        F.reduced_formula("si o2")


def test_optimade_no_match_and_unreachable_provider(routes):
    _, extra = routes
    extra["https://oqmd.org/optimade/v1/structures?filter=chemical_formula_reduced%3D%22Xe%22&page_limit=5"
          "&response_fields=chemical_formula_descriptive%2Cchemical_formula_reduced%2Cnsites"] = b'{"data": []}'
    with pytest.raises(FetchError, match="どのデータベースにも chemical_formula_reduced=\"Xe\" に合う項目がありません"):
        F.fetch("optimade:Xe")
    with pytest.raises(FetchError, match="OPTIMADE の提供元一覧に 'nowhere' がありません"):
        F.fetch("optimade:nowhere:1")


def test_no_network_stops_with_reason(monkeypatch):
    def down(req, timeout=None):
        raise urllib.error.URLError("Name or service not known")

    monkeypatch.setattr(urllib.request, "urlopen", down)
    with pytest.raises(FetchError, match=r"取得できません: PubChem に接続できません \(Name or service not known\)"):
        F.fetch("pubchem:water")
    with pytest.raises(FetchError, match="取得できません: COD に接続できません"):
        F.fetch("cod:1000041")
    with pytest.raises(FetchError, match="取得できません: providers.optimade.org に接続できません"):
        F.fetch("optimade:Si")


def test_bad_refs():
    for ref in ("", "water", "icsd:1", "pubchem:"):
        with pytest.raises(FetchError, match="データベース:名前または ID"):
            F.fetch(ref)


def test_readme_and_summary_lines(routes, tmp_path):
    res = F.fetch("cod:1000041")
    res.fetched.save(tmp_path)
    lines = F.readme_lines(res.fetched.record)
    text = "\n".join(lines)
    assert "== 構造の出どころ (データベースから取得) ==" in text and "CC0 1.0" in text and "SHA-256" in text
    assert "cod_1000041.cif" in text and "DOI: 10.1107/S0365110X65002244" in text
    assert "Crystallography Open Database (COD) 1000041" in F.summary_line(res.fetched.record)


def test_structure_and_project_carry_the_record(routes, tmp_path, sk_root):
    from adit.project import build_project
    from adit.structure import from_fetched
    from tests.conftest import cfg_for, water_spec

    res = F.fetch("pubchem:water")
    res.fetched.save(tmp_path / "fetched")
    st = from_fetched(res.fetched, charge=0, multiplicity=1)
    assert st.source == "file" and st.source_ref.endswith("pubchem_962.sdf") and st.fetched["id"] == "962"
    spec = water_spec().model_copy(update={"structure": st})
    spec2 = spec.model_validate_json(spec.model_dump_json())
    assert spec2.structure.fetched == st.fetched
    files = build_project(spec2, cfg_for(sk_root), output_dir=tmp_path / "calc")
    assert "pubchem_962.sdf" in files.copies
    prov = json.loads(files.texts["spec.json"])["provenance"]
    assert prov["fetched_structure"]["sha256"] == st.fetched["sha256"]
    assert any(f["name"] == "pubchem_962.sdf" and f["sha256"] == st.fetched["sha256"] for f in prov["files"])
    readme = files.texts["README.txt"]
    assert "== 構造の出どころ (データベースから取得) ==" in readme and "public domain" in readme and "pubchem_962.sdf" in readme
