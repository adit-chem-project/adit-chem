"""Fetch structures from public databases (PubChem, COD, Materials Project, OPTIMADE) with urllib only."""

from __future__ import annotations

import hashlib
import io
import json
import math
import re
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from ase import Atoms

from adit import __version__
from adit.errors import AditError
from adit.lang import L

TIMEOUT = 30.0
OPTIMADE_TIMEOUT = 20.0
OPTIMADE_PAGE_LIMIT = 5
MAX_CANDIDATES = 40

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
COD_BASE = "https://www.crystallography.net/cod"
MP_BASE = "https://api.materialsproject.org"
OPTIMADE_PROVIDERS_URL = "https://providers.optimade.org/providers.json"

DATABASES: dict[str, tuple[str, str]] = {
    "pubchem": ("PubChem", "https://pubchem.ncbi.nlm.nih.gov/"),
    "cod": ("Crystallography Open Database (COD)", "https://www.crystallography.net/cod/"),
    "mp": ("Materials Project", "https://next-gen.materialsproject.org/"),
    "optimade": ("OPTIMADE", "https://www.optimade.org/"),
}

# License wording as published by each site (see docs/USAGE.md for the pages checked).
LICENSES: dict[str, tuple[str, str]] = {
    "pubchem": ("NCBI policy: information created by the US government is in the public domain; "
                "individual data sources may assert their own terms", "https://www.ncbi.nlm.nih.gov/home/about/policies/"),
    "cod": ("CC0 1.0 (dedicated to the public domain)", "https://creativecommons.org/publicdomain/zero/1.0/"),
    "mp": ("CC BY 4.0 (attribution to the Materials Project required)", "https://creativecommons.org/licenses/by/4.0/"),
    "optimade": ("depends on the provider; check the provider's terms of use", ""),
}


class FetchError(AditError):
    pass


@dataclass(frozen=True)
class Candidate:
    ref: str
    label: str


@dataclass
class Fetched:
    atoms: Atoms
    record: dict
    data: bytes
    filename: str
    structure_text: str | None = None
    structure_filename: str | None = None

    def save(self, directory: Path | str) -> Path:
        """Write the response (and a readable structure file when the response is not one) and return the structure path."""
        d = Path(directory).expanduser()
        d.mkdir(parents=True, exist_ok=True)
        raw = d / self.filename
        raw.write_bytes(self.data)
        path = raw
        if self.structure_text is not None and self.structure_filename:
            path = d / self.structure_filename
            path.write_text(self.structure_text, encoding="utf-8", newline="\n")
        self.record["file"] = str(path)
        self.record["response_file"] = str(raw)
        return path


@dataclass
class FetchResult:
    fetched: Fetched | None = None
    candidates: list[Candidate] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)


def default_fetch_dir() -> Path:
    return Path.home() / "adit_runs" / "fetched"


def database_label(db: str) -> str:
    return DATABASES.get(db, (db, ""))[0]


def parse_ref(ref: str) -> tuple[str, str]:
    db, sep, rest = (ref or "").strip().partition(":")
    db = db.strip().lower()
    if not sep or db not in DATABASES or not rest.strip():
        names = ", ".join(f"{k}:<...>" for k in DATABASES)
        raise FetchError(L(f"取得の指定は「データベース:名前または ID」の形で書いてください ({names})。例 pubchem:water、cod:1000041、mp:mp-149、optimade:Si",
                           f"write the fetch reference as database:name-or-id ({names}), e.g. pubchem:water, cod:1000041, mp:mp-149, optimade:Si"))
    return db, rest.strip()


def _cannot(reason: str) -> FetchError:
    return FetchError(L(f"取得できません: {reason}", f"cannot fetch: {reason}"))


def _get(url: str, *, what: str, headers: dict[str, str] | None = None, timeout: float = TIMEOUT) -> bytes:
    req = urllib.request.Request(url, headers={"User-Agent": f"ADIT/{__version__}", "Accept": "*/*", **(headers or {})})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.read()
    except urllib.error.HTTPError as ex:
        raise _cannot(L(f"{what} が HTTP {ex.code} を返しました ({url})", f"{what} returned HTTP {ex.code} ({url})")) from ex
    except urllib.error.URLError as ex:
        raise _cannot(L(f"{what} に接続できません ({ex.reason})", f"cannot connect to {what} ({ex.reason})")) from ex
    except (TimeoutError, OSError) as ex:
        raise _cannot(L(f"{what} との通信に失敗しました ({ex})", f"communication with {what} failed ({ex})")) from ex


def _json(data: bytes, what: str) -> dict:
    try:
        return json.loads(data.decode("utf-8"))
    except (ValueError, UnicodeDecodeError) as ex:
        raise _cannot(L(f"{what} の応答を JSON として読めません ({ex})", f"the response of {what} is not JSON ({ex})")) from ex


def _record(db: str, query: str, entry_id: str, url: str, data: bytes, fmt: str, *, title: str = "", provider: str = "",
            license_text: str | None = None, license_url: str | None = None, extra: dict | None = None) -> dict:
    lic, lic_url = LICENSES[db]
    rec = {
        "database": db, "database_name": database_label(db), "provider": provider, "query": query, "id": entry_id,
        "title": title, "url": url, "fetched_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data), "format": fmt,
        "license": license_text if license_text is not None else lic, "license_url": license_url if license_url is not None else lic_url,
    }
    if extra:
        rec.update(extra)
    return rec


def _read_atoms(text: str, fmt: str, what: str) -> Atoms:
    from ase.io import read

    try:
        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            result = read(io.StringIO(text), format=fmt)
    except Exception as ex:
        raise _cannot(L(f"{what} の応答を {fmt} として読めません ({ex})", f"the response of {what} could not be read as {fmt} ({ex})")) from ex
    atoms = result[-1] if isinstance(result, list) else result
    if len(atoms) == 0:
        raise _cannot(L(f"{what} の応答に原子がありません", f"the response of {what} has no atoms"))
    return atoms


def _extxyz(atoms: Atoms) -> str:
    from ase.io import write

    buf = io.StringIO()
    write(buf, atoms, format="extxyz")
    return buf.getvalue()


def _safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("_") or "entry"


# ---- PubChem (PUG REST) ----
def _pubchem_candidates(name: str, cids: list[int]) -> list[Candidate]:
    ids = ",".join(str(c) for c in cids[:MAX_CANDIDATES])
    url = f"{PUBCHEM_BASE}/compound/cid/{ids}/property/MolecularFormula,IUPACName,Title/JSON"
    props = _json(_get(url, what="PubChem"), "PubChem").get("PropertyTable", {}).get("Properties", [])
    out = []
    for p in props:
        cid = p.get("CID")
        label = " ".join(str(x) for x in (p.get("Title") or p.get("IUPACName") or "", p.get("MolecularFormula") or "") if x)
        out.append(Candidate(ref=f"pubchem:{cid}", label=f"CID {cid}: {label}"))
    return out


def fetch_pubchem(query: str) -> FetchResult:
    query = query.strip()
    if query.isdigit():
        cid = int(query)
    else:
        url = f"{PUBCHEM_BASE}/compound/name/{urllib.parse.quote(query, safe='')}/cids/JSON"
        try:
            data = _get(url, what="PubChem")
        except FetchError as ex:
            if "HTTP 404" in str(ex):
                raise _cannot(L(f"PubChem に {query!r} という名前の化合物がありません", f"PubChem has no compound named {query!r}")) from ex
            raise
        cids = [int(c) for c in _json(data, "PubChem").get("IdentifierList", {}).get("CID", [])]
        if not cids:
            raise _cannot(L(f"PubChem に {query!r} という名前の化合物がありません", f"PubChem has no compound named {query!r}"))
        if len(cids) > 1:
            return FetchResult(candidates=_pubchem_candidates(query, cids))
        cid = cids[0]
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/SDF?record_type=3d"
    try:
        data = _get(url, what="PubChem")
    except FetchError as ex:
        if "HTTP 404" in str(ex):
            raise _cannot(L(f"PubChem の CID {cid} には 3D 座標がありません (3D 配座は一部の化合物にだけ用意されています)",
                            f"PubChem CID {cid} has no 3D coordinates (3D conformers exist only for some compounds)")) from ex
        raise
    atoms = _read_atoms(data.decode("utf-8", errors="replace"), "sdf", "PubChem")
    title = ""
    try:
        props = _json(_get(f"{PUBCHEM_BASE}/compound/cid/{cid}/property/MolecularFormula,Title/JSON", what="PubChem"), "PubChem")
        p = (props.get("PropertyTable", {}).get("Properties") or [{}])[0]
        title = " ".join(str(x) for x in (p.get("Title", ""), p.get("MolecularFormula", "")) if x)
    except FetchError:
        pass
    rec = _record("pubchem", query, str(cid), url, data, "sdf", title=title)
    return FetchResult(fetched=Fetched(atoms=atoms, record=rec, data=data, filename=f"pubchem_{cid}.sdf"))


# ---- COD ----
def fetch_cod(query: str) -> FetchResult:
    cod_id = query.strip()
    if not cod_id.isdigit():
        raise _cannot(L(f"COD の ID は数字です (例 1000041): {query!r}", f"a COD ID is a number (e.g. 1000041): {query!r}"))
    url = f"{COD_BASE}/{cod_id}.cif"
    try:
        data = _get(url, what="COD")
    except FetchError as ex:
        if "HTTP 404" in str(ex):
            raise _cannot(L(f"COD に ID {cod_id} の項目がありません", f"COD has no entry {cod_id}")) from ex
        raise
    atoms = _read_atoms(data.decode("utf-8", errors="replace"), "cif", "COD")
    title, extra = "", {}
    try:
        meta = _json(_get(f"{COD_BASE}/result?id={cod_id}&format=json", what="COD"), "COD")
        if isinstance(meta, list) and meta:
            m = meta[0]
            title = " / ".join(str(x) for x in (m.get("chemname") or m.get("formula") or "", m.get("title") or "") if x)
            extra = {"authors": m.get("authors") or "", "journal": " ".join(str(x) for x in (m.get("journal") or "", m.get("year") or "") if x),
                     "doi": m.get("doi") or ""}
    except FetchError:
        pass
    rec = _record("cod", query, cod_id, url, data, "cif", title=title, extra=extra)
    return FetchResult(fetched=Fetched(atoms=atoms, record=rec, data=data, filename=f"cod_{cod_id}.cif"))


# ---- Materials Project (new API) ----
def _atoms_from_pymatgen(struct: dict, what: str) -> Atoms:
    try:
        matrix = struct["lattice"]["matrix"]
        sites = struct["sites"]
    except (KeyError, TypeError) as ex:
        raise _cannot(L(f"{what} の構造の形式を読めません ({ex})", f"cannot read the structure format of {what} ({ex})")) from ex
    symbols, positions = [], []
    for site in sites:
        species = site.get("species") or []
        if len(species) != 1 or abs(float(species[0].get("occu", 1.0)) - 1.0) > 1e-6:
            raise _cannot(L(f"{what} の構造に部分占有のサイトがあります。ADIT では扱えないので、別の項目を選んでください",
                            f"the {what} structure has a partially occupied site; ADIT cannot use it, choose another entry"))
        symbols.append(species[0]["element"])
        positions.append([float(x) for x in site["xyz"]])
    pbc = struct["lattice"].get("pbc") or (True, True, True)
    return Atoms(symbols=symbols, positions=positions, cell=matrix, pbc=tuple(bool(b) for b in pbc))


def fetch_mp(query: str, api_key: str = "") -> FetchResult:
    mp_id = query.strip()
    if mp_id.isdigit():
        mp_id = f"mp-{mp_id}"
    if not re.fullmatch(r"(mp|mvc)-\d+", mp_id):
        raise _cannot(L(f"Materials Project の ID は mp-149 のような形です: {query!r}", f"a Materials Project ID looks like mp-149: {query!r}"))
    if not api_key.strip():
        raise _cannot(L("Materials Project の API キーが環境設定 (cluster.toml の mp_api_key) にありません。"
                        "キーは https://next-gen.materialsproject.org/dashboard で発行できます",
                        "no Materials Project API key in the settings (mp_api_key in cluster.toml). "
                        "Get one at https://next-gen.materialsproject.org/dashboard"))
    q = urllib.parse.urlencode({"material_ids": mp_id, "_fields": "material_id,formula_pretty,structure,symmetry,deprecated"})
    url = f"{MP_BASE}/materials/summary/?{q}"
    try:
        data = _get(url, what="Materials Project", headers={"X-API-KEY": api_key.strip()})
    except FetchError as ex:
        if "HTTP 401" in str(ex) or "HTTP 403" in str(ex):
            raise _cannot(L("Materials Project が API キーを受け付けませんでした (環境設定の mp_api_key を確かめてください)",
                            "the Materials Project rejected the API key (check mp_api_key in the settings)")) from ex
        raise
    docs = _json(data, "Materials Project").get("data") or []
    if not docs:
        raise _cannot(L(f"Materials Project に {mp_id} がありません", f"the Materials Project has no {mp_id}"))
    doc = docs[0]
    atoms = _atoms_from_pymatgen(doc.get("structure") or {}, "Materials Project")
    sym = doc.get("symmetry") or {}
    title = " ".join(str(x) for x in (doc.get("formula_pretty") or "", sym.get("symbol") or "") if x)
    extra = {"deprecated": bool(doc.get("deprecated", False))}
    rec = _record("mp", query, mp_id, url, data, "json", title=title, extra=extra)
    return FetchResult(fetched=Fetched(atoms=atoms, record=rec, data=data, filename=f"{mp_id}.json",
                                       structure_text=_extxyz(atoms), structure_filename=f"{mp_id}.extxyz"))


# ---- OPTIMADE ----
_FORMULA_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def reduced_formula(text: str) -> str:
    """OPTIMADE chemical_formula_reduced: elements in alphabetical order, smallest integer proportions, 1 omitted."""
    s = text.replace(" ", "")
    if not s or "".join(m.group(0) for m in _FORMULA_TOKEN.finditer(s)) != s:
        raise _cannot(L(f"組成式として読めません: {text!r} (例 Si、SiO2、NaCl)", f"not a chemical formula: {text!r} (e.g. Si, SiO2, NaCl)"))
    counts: dict[str, int] = {}
    for m in _FORMULA_TOKEN.finditer(s):
        counts[m.group(1)] = counts.get(m.group(1), 0) + (int(m.group(2)) if m.group(2) else 1)
    g = 0
    for n in counts.values():
        g = math.gcd(g, n)
    return "".join(f"{el}{counts[el] // g if counts[el] // g != 1 else ''}" for el in sorted(counts))


def optimade_filter(query: str) -> str:
    q = query.strip()
    if any(op in q for op in ("=", "<", ">", " HAS ", " has ")):
        return q
    return f'chemical_formula_reduced="{reduced_formula(q)}"'


@dataclass(frozen=True)
class OptimadeDatabase:
    id: str
    name: str
    base_url: str
    provider: str
    homepage: str


def _providers() -> list[dict]:
    data = _json(_get(OPTIMADE_PROVIDERS_URL, what="providers.optimade.org"), "providers.optimade.org").get("data") or []
    return [p for p in data if (p.get("attributes") or {}).get("base_url") and p.get("id") != "exmpl"]


def _child_databases(provider: dict) -> list[OptimadeDatabase]:
    attrs = provider.get("attributes") or {}
    base = str(attrs["base_url"]).rstrip("/")
    links = _json(_get(f"{base}/v1/links", what=f"OPTIMADE ({provider.get('id')})", timeout=OPTIMADE_TIMEOUT),
                  f"OPTIMADE ({provider.get('id')})").get("data") or []
    out = []
    for link in links:
        a = link.get("attributes") or {}
        if a.get("link_type") == "child" and a.get("base_url"):
            out.append(OptimadeDatabase(id=str(link.get("id")), name=str(a.get("name") or link.get("id")), base_url=str(a["base_url"]).rstrip("/"),
                                        provider=str(provider.get("id")), homepage=str(attrs.get("homepage") or a.get("homepage") or "")))
    return out


def optimade_databases(provider_id: str | None = None) -> tuple[list[OptimadeDatabase], list[str]]:
    """All child databases of the registered providers (or of the provider whose id matches), with per-provider failure notes."""
    providers = _providers()
    if provider_id:
        exact = [p for p in providers if p.get("id") == provider_id]
        providers = exact or [p for p in providers if provider_id.startswith(str(p.get("id")) + "-") or provider_id.startswith(str(p.get("id")))]
        if not providers:
            raise _cannot(L(f"OPTIMADE の提供元一覧に {provider_id!r} がありません", f"{provider_id!r} is not in the OPTIMADE providers list"))
    dbs: list[OptimadeDatabase] = []
    notes: list[str] = []

    def one(p: dict):
        try:
            return _child_databases(p), ""
        except FetchError as ex:
            return [], f"{p.get('id')}: {ex}"

    with ThreadPoolExecutor(max_workers=8) as pool:
        for found, note in pool.map(one, providers):
            dbs.extend(found)
            if note:
                notes.append(note)
    if provider_id:
        dbs = [d for d in dbs if d.id == provider_id] or dbs
    return dbs, notes


def _atoms_from_optimade(attrs: dict, what: str) -> Atoms:
    species = {s.get("name"): s for s in attrs.get("species") or []}
    symbols = []
    for name in attrs.get("species_at_sites") or []:
        sp = species.get(name) or {}
        syms = sp.get("chemical_symbols") or []
        conc = sp.get("concentration") or [1.0]
        if len(syms) != 1 or abs(float(conc[0]) - 1.0) > 1e-6 or syms[0] in ("X", "vacancy"):
            raise _cannot(L(f"{what} の構造に部分占有か空孔のサイトがあります。ADIT では扱えないので、別の項目を選んでください",
                            f"the {what} structure has a partially occupied or vacancy site; ADIT cannot use it, choose another entry"))
        symbols.append(syms[0])
    positions = attrs.get("cartesian_site_positions") or []
    if not symbols or len(symbols) != len(positions):
        raise _cannot(L(f"{what} の構造にサイトの情報がありません", f"the {what} structure has no usable site information"))
    cell = attrs.get("lattice_vectors")
    dims = attrs.get("dimension_types") or ([1, 1, 1] if cell else [0, 0, 0])
    if cell is None or any(v is None for row in cell for v in row):
        return Atoms(symbols=symbols, positions=positions)
    return Atoms(symbols=symbols, positions=positions, cell=cell, pbc=tuple(bool(d) for d in dims))


def _optimade_label(db: OptimadeDatabase, entry: dict) -> str:
    a = entry.get("attributes") or {}
    parts = [f"{db.id}:{entry.get('id')}", str(a.get("chemical_formula_descriptive") or a.get("chemical_formula_reduced") or "")]
    if a.get("nsites") is not None:
        parts.append(L(f"{a['nsites']} 原子", f"{a['nsites']} atoms"))
    sg = a.get("space_group_symbol_hermann_mauguin")
    if sg:
        parts.append(str(sg))
    return "  ".join(p for p in parts if p)


def _optimade_search(db: OptimadeDatabase, flt: str) -> tuple[list[Candidate], str]:
    q = urllib.parse.urlencode({"filter": flt, "page_limit": OPTIMADE_PAGE_LIMIT,
                                "response_fields": "chemical_formula_descriptive,chemical_formula_reduced,nsites"})
    try:
        data = _json(_get(f"{db.base_url}/v1/structures?{q}", what=f"OPTIMADE ({db.id})", timeout=OPTIMADE_TIMEOUT), f"OPTIMADE ({db.id})")
    except FetchError as ex:
        return [], str(ex)
    entries = data.get("data") or []
    return [Candidate(ref=f"optimade:{db.id}:{e.get('id')}", label=_optimade_label(db, e)) for e in entries if e.get("id") is not None], ""


def fetch_optimade(query: str) -> FetchResult:
    query = query.strip()
    head, sep, tail = query.partition(":")
    if sep and tail.strip() and not any(op in query for op in ("=", "<", ">")):
        return _fetch_optimade_entry(head.strip(), tail.strip())
    flt = optimade_filter(query)
    dbs, notes = optimade_databases()
    if not dbs:
        raise _cannot(L("OPTIMADE のデータベースを 1 つも見つけられませんでした", "no OPTIMADE database could be reached") + ("\n" + "\n".join(notes) if notes else ""))
    candidates: list[Candidate] = []
    with ThreadPoolExecutor(max_workers=8) as pool:
        for found, note in pool.map(lambda d: _optimade_search(d, flt), dbs):
            candidates.extend(found)
            if note:
                notes.append(note)
    if not candidates:
        raise _cannot(L(f"OPTIMADE のどのデータベースにも {flt} に合う項目がありません ({len(dbs)} 件を検索)",
                        f"no OPTIMADE database has an entry matching {flt} ({len(dbs)} searched)") + ("\n" + "\n".join(notes) if notes else ""))
    return FetchResult(candidates=candidates[:MAX_CANDIDATES], notes=notes)


def _fetch_optimade_entry(db_id: str, entry_id: str) -> FetchResult:
    dbs, notes = optimade_databases(db_id)
    matches = [d for d in dbs if d.id == db_id]
    if not matches:
        raise _cannot(L(f"OPTIMADE のデータベース {db_id!r} が見つかりません", f"OPTIMADE database {db_id!r} was not found") + ("\n" + "\n".join(notes) if notes else ""))
    db = matches[0]
    url = f"{db.base_url}/v1/structures/{urllib.parse.quote(entry_id, safe='')}"
    data = _get(url, what=f"OPTIMADE ({db.id})", timeout=OPTIMADE_TIMEOUT)
    entry = _json(data, f"OPTIMADE ({db.id})").get("data")
    if isinstance(entry, list):
        entry = entry[0] if entry else None
    if not entry:
        raise _cannot(L(f"{db.id} に {entry_id} がありません", f"{db.id} has no entry {entry_id}"))
    attrs = entry.get("attributes") or {}
    atoms = _atoms_from_optimade(attrs, f"OPTIMADE ({db.id})")
    title = " ".join(str(x) for x in (attrs.get("chemical_formula_descriptive") or attrs.get("chemical_formula_reduced") or "",
                                      attrs.get("space_group_symbol_hermann_mauguin") or "") if x)
    lic = L(f"提供元 {db.name} の規約による ({db.homepage})", f"as set by the provider {db.name} ({db.homepage})") if db.homepage else None
    rec = _record("optimade", f"{db_id}:{entry_id}", str(entry.get("id")), url, data, "json", title=title, provider=db.id,
                  license_text=lic, license_url=db.homepage or None, extra={"provider_name": db.name, "base_url": db.base_url})
    name = _safe_name(f"optimade_{db.id}_{entry.get('id')}")
    return FetchResult(fetched=Fetched(atoms=atoms, record=rec, data=data, filename=f"{name}.json",
                                       structure_text=_extxyz(atoms), structure_filename=f"{name}.extxyz"), notes=notes)


# ---- entry point ----
def fetch(ref: str, *, mp_api_key: str = "") -> FetchResult:
    """Fetch by 'db:query'. Returns a structure, or candidates to choose from (the choice is the user's)."""
    db, rest = parse_ref(ref)
    if db == "pubchem":
        return fetch_pubchem(rest)
    if db == "cod":
        return fetch_cod(rest)
    if db == "mp":
        return fetch_mp(rest, mp_api_key)
    return fetch_optimade(rest)


def candidates_text(result: FetchResult) -> str:
    lines = [L(f"候補が {len(result.candidates)} 件あります。1 つ選んで、その指定で取得し直してください:",
               f"{len(result.candidates)} candidates; choose one and fetch again with its reference:")]
    lines += [f"  {c.ref:<40} {c.label}" for c in result.candidates]
    if result.notes:
        lines.append(L("応答が無かった提供元:", "providers that did not answer:"))
        lines += [f"  {n}" for n in result.notes]
    return "\n".join(lines)


def summary_line(rec: dict) -> str:
    name = rec.get("database_name") or rec.get("database", "")
    prov = f" ({rec['provider']})" if rec.get("provider") else ""
    title = f" {rec['title']}" if rec.get("title") else ""
    return L(f"{name}{prov} {rec.get('id', '')}{title} を {rec.get('fetched_utc', '')} に取得。ライセンス: {rec.get('license', '')}",
             f"{name}{prov} {rec.get('id', '')}{title}, fetched {rec.get('fetched_utc', '')}. License: {rec.get('license', '')}")


def readme_lines(rec: dict) -> list[str]:
    lines = [L("== 構造の出どころ (データベースから取得) ==", "== Structure origin (fetched from a database) =="),
             L(f"  データベース: {rec.get('database_name', '')}" + (f" / 提供元 {rec['provider']}" if rec.get("provider") else ""),
               f"  Database: {rec.get('database_name', '')}" + (f" / provider {rec['provider']}" if rec.get("provider") else "")),
             L(f"  指定: {rec.get('query', '')} / ID: {rec.get('id', '')}" + (f" / {rec['title']}" if rec.get("title") else ""),
               f"  Query: {rec.get('query', '')} / ID: {rec.get('id', '')}" + (f" / {rec['title']}" if rec.get("title") else "")),
             f"  URL: {rec.get('url', '')}",
             L(f"  取得日時: {rec.get('fetched_utc', '')} (UTC) / 応答の SHA-256: {rec.get('sha256', '')} ({rec.get('bytes', 0)} バイト)",
               f"  Fetched: {rec.get('fetched_utc', '')} (UTC) / SHA-256 of the response: {rec.get('sha256', '')} ({rec.get('bytes', 0)} bytes)"),
             L(f"  ライセンス: {rec.get('license', '')}" + (f" ({rec['license_url']})" if rec.get("license_url") else ""),
               f"  License: {rec.get('license', '')}" + (f" ({rec['license_url']})" if rec.get("license_url") else ""))]
    for key, ja, en in (("doi", "DOI", "DOI"), ("authors", "著者", "Authors"), ("journal", "掲載", "Journal")):
        if rec.get(key):
            lines.append(f"  {L(ja, en)}: {rec[key]}")
    if rec.get("deprecated"):
        lines.append(L("  注意: Materials Project でこの項目は deprecated (非推奨) の印が付いています", "  Note: this entry is marked deprecated by the Materials Project"))
    files = [Path(rec[k]).name for k in ("file", "response_file") if rec.get(k)]
    if files:
        lines.append(L(f"  写したファイル: {', '.join(dict.fromkeys(files))} (元の置き場所 {rec.get('file', '')})",
                       f"  Copied files: {', '.join(dict.fromkeys(files))} (original location {rec.get('file', '')})"))
    lines.append(L("  取得したデータをそのまま公開・再配布するときは、上のライセンスの条件 (出典の表示など) に従ってください",
                   "  When publishing or redistributing the fetched data, follow the license above (attribution etc.)"))
    return lines + [""]
