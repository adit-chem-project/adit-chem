
from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path

from adit.lang import L

MAX_SCAN_LINES = 3000
MAX_JSON_BYTES = 5 * 2 ** 20
_NUM = r"(-?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?)"


@dataclass(frozen=True)
class DocValue:
    element: str | None
    quantity: str
    value: float | str
    unit: str
    source: str
    line: int
    text: str

    def as_dict(self) -> dict:
        return asdict(self)

    def __str__(self) -> str:
        v = f"{self.value:g}" if isinstance(self.value, float) else str(self.value)
        return f"{self.element or '-'}  {self.quantity} = {v} {self.unit}".rstrip() + f"   ({self.source}:{self.line})"


def _f(s: str) -> float:
    return float(s.replace("D", "E").replace("d", "e"))


def _lines(path: Path, limit: int = MAX_SCAN_LINES):
    with open(path, encoding="utf-8", errors="replace") as f:
        for no, line in enumerate(f, 1):
            if no > limit:
                return
            yield no, line.rstrip("\n")


# ---- pw.x ----
def _upf_values(path: Path, element: str) -> list[DocValue]:
    out = []
    in_header = False
    for no, line in _lines(path):
        if "<PP_HEADER" in line:
            in_header = True
        if in_header:
            for key in ("wfc_cutoff", "rho_cutoff"):
                m = re.search(rf'\b{key}\s*=\s*"\s*{_NUM}\s*"', line)
                if m and _f(m.group(1)) > 0:
                    out.append(DocValue(element, key, _f(m.group(1)), "Ry", str(path), no, line.strip()))
        if in_header and ("</PP_HEADER>" in line or line.rstrip().endswith("/>")):
            in_header = False
        m = re.search(rf"Suggested minimum cutoff for wavefunctions:\s*{_NUM}\s*(Ry|Ha)", line, re.I)
        if m:
            out.append(DocValue(element, "suggested_wfc_cutoff", _f(m.group(1)), m.group(2), str(path), no, line.strip()))
        m = re.search(rf"Suggested minimum cutoff for charge density:\s*{_NUM}\s*(Ry|Ha)", line, re.I)
        if m:
            out.append(DocValue(element, "suggested_rho_cutoff", _f(m.group(1)), m.group(2), str(path), no, line.strip()))
        if "<PP_MESH" in line:
            break
    return out


def _sssp_values(dirs: list[Path], element: str, upf: str) -> list[DocValue]:
    out = []
    seen = set()
    for d in dirs:
        if not d.is_dir():
            continue
        for p in sorted(d.glob("*.json")):
            if p in seen or p.stat().st_size > MAX_JSON_BYTES:
                continue
            seen.add(p)
            try:
                text = p.read_text(encoding="utf-8", errors="replace")
                js = json.loads(text)
            except (OSError, json.JSONDecodeError):
                continue
            e = js.get(element) if isinstance(js, dict) else None
            if not isinstance(e, dict) or e.get("filename") != upf:
                continue
            m = re.search(rf'"{re.escape(element)}"\s*:', text)
            no = text.count("\n", 0, m.start()) + 1 if m else 1
            for key in ("cutoff_wfc", "cutoff_rho"):
                if isinstance(e.get(key), (int, float)):
                    out.append(DocValue(element, f"sssp_{key}", float(e[key]), "Ry", str(p), no,
                                        f'"{element}": {{"filename": "{upf}", "{key}": {e[key]}}}'))
    return out


def _espresso(spec, cfg) -> list[DocValue]:
    from adit.codes.espresso import upf_name
    from adit.codes.upf import UpfLibrary

    lib = UpfLibrary.open_if_present(cfg.pseudo_root or None, spec.method.pseudo_set)
    if lib is None:
        return []
    out = []
    for el in spec.elements:
        name = upf_name(spec.method, lib, el)
        if not name or not lib.has(name):
            continue
        out += _upf_values(lib.path(name), el)
        out += _sssp_values([lib.set_dir, lib.root], el, name)
    return out


# ---- VASP ----
def _vasp(spec, cfg) -> list[DocValue]:
    from adit.codes.potcar import PotcarLibrary
    from adit.codes.vasp import PP_ENV, potcar_name

    prof = cfg.profiles.get(spec.runtime.profile)
    lib = PotcarLibrary.open_if_present((prof.env.get(PP_ENV) if prof else None) or None, spec.method.potcar_set)
    if lib is None:
        return []
    out = []
    for el in spec.elements:
        name = potcar_name(spec.method, el)
        p = lib.set_dir / name / "POTCAR"
        if not p.is_file():
            continue
        for no, line in _lines(p, 80):
            for key in ("ENMAX", "ENMIN"):
                m = re.search(rf"\b{key}\s*=\s*{_NUM}", line)
                if m:
                    out.append(DocValue(el, key, _f(m.group(1)), "eV", str(p), no, line.strip()))
    return out


# ---- CP2K ----
def _cp2k(spec, cfg) -> list[DocValue]:
    from adit.codes.cp2k import _data, _pick

    data = _data(cfg)
    m, out = spec.method, []
    for el in spec.elements:
        for fname, potential in ((m.basis_file, False), (m.potential_file, True)):
            p = data.resolve(fname)
            name = _pick(m, data, el, potential=potential) if p else None
            if not p or not name:
                continue
            comments: list[tuple[int, str]] = []
            for no, line in _lines(p, 10 ** 7):
                s = line.strip()
                if s.startswith("#"):
                    comments.append((no, s))
                    continue
                code = s.split("#", 1)[0].split()
                if len(code) >= 2 and code[0] == el and name in code[1:]:
                    for cno, c in comments:
                        if c.strip("# ").strip():
                            out.append(DocValue(el, "comment", c.lstrip("#").strip(), "", str(p), cno, c))
                    if "#" in s:
                        out.append(DocValue(el, "comment", s.split("#", 1)[1].strip(), "", str(p), no, s))
                    break
                if s:
                    comments = []
    return out


# ---- DFTB+ ----
def _dftb(spec, cfg) -> list[DocValue]:
    if not cfg.sk_root:
        return []
    d = Path(cfg.sk_root).expanduser() / spec.method.sk_set
    out = []
    readme = d / "README"
    if readme.is_file():
        in_hub = False
        for no, line in _lines(readme, 10 ** 6):
            s = line.strip()
            if "hubbard derivative" in s.lower():
                in_hub = True
                continue
            m = re.match(rf"^([A-Z][a-z]?)\s*=\s*{_NUM}\s*$", s)
            if in_hub and m:
                if m.group(1) in spec.elements:
                    out.append(DocValue(m.group(1), "hubbard_derivative", _f(m.group(2)), "", str(readme), no, s))
                continue
            if in_hub and not s:
                in_hub = False
            m = re.search(rf"\bzeta\s*=\s*{_NUM}", s, re.I)
            if m:
                out.append(DocValue(None, "zeta", _f(m.group(1)), "", str(readme), no, s))
            for key in ("s6", "s8", "a1", "a2"):
                m = re.search(rf"(?<![A-Za-z0-9_]){key}\s*[=:]\s*{_NUM}", s)
                if m:
                    out.append(DocValue(None, key, _f(m.group(1)), "", str(readme), no, s))
    spinw = d / "spinw.txt"
    if spinw.is_file():
        cur, last = None, None
        for no, line in _lines(spinw, 10 ** 6):
            m = re.match(r"^\s*([A-Z][a-z]?)\s*:\s*$", line)
            if m:
                if cur in spec.elements and last:
                    out.append(DocValue(cur, "spin_constant", last[0], "", str(spinw), last[1], last[2]))
                cur, last = m.group(1), None
            elif line.strip() and cur:
                w = line.split()
                try:
                    last = (float(w[-1]), no, line.strip())
                except ValueError:
                    pass
        if cur in spec.elements and last:
            out.append(DocValue(cur, "spin_constant", last[0], "", str(spinw), last[1], last[2]))
    return out


_READERS = {"espresso": _espresso, "vasp": _vasp, "cp2k": _cp2k, "dftbplus": _dftb}


def documented_values(spec, cfg=None) -> list[DocValue]:
    if cfg is None:
        from adit.config import Config, ConfigError, load_config
        try:
            cfg = load_config()
        except (ConfigError, OSError):
            cfg = Config()
    reader = _READERS.get(spec.method.code)
    return reader(spec, cfg) if reader else []




@dataclass(frozen=True)
class DocSetting:
    code: str
    tasks: tuple[str, ...]
    group: str
    field: str
    key: str
    value: str
    source: str
    url: str
    retrieved: str
    note: str = ""
    note_en: str = ""
    source_en: str = ""


_VASP_PHONON_SRC = ("VASP tutorial: Phonons — Part 1 (Graphene) の INCAR の例",
                    "https://www.vasp.at/tutorials/latest/phonon/part1/", "2026-09-12")
_VASP_PHONON_SRC_EN = "the INCAR example in the VASP tutorial Phonons — Part 1 (Graphene)"
DOCUMENTED_SETTINGS: list[DocSetting] = [
    DocSetting("vasp", ("vibrations",), "vasp_incar", "prec", "PREC", "Accurate", *_VASP_PHONON_SRC,
               "2 階微分には精度の高い力が要る、と wiki の「Phonons from finite differences」にも書かれています",
               "the wiki page 'Phonons from finite differences' also states that accurate forces are needed for the second derivatives",
               source_en=_VASP_PHONON_SRC_EN),
    DocSetting("vasp", ("vibrations",), "vasp_incar", "ediff", "EDIFF", "1e-8", *_VASP_PHONON_SRC,
               "電子の収束が甘いと力が不正確になる (文書の言葉では spurious convergence を避けるため)",
               "if the electronic convergence is too loose the forces are imprecise (the text speaks of avoiding a spurious convergence)",
               source_en=_VASP_PHONON_SRC_EN),
    DocSetting("vasp", ("vibrations",), "vasp_incar", "nelmin", "NELMIN", "5", *_VASP_PHONON_SRC,
               "イオンの 1 ステップあたり最低 5 回は電子の反復を回す",
               "at least 5 electronic iterations per ionic step", source_en=_VASP_PHONON_SRC_EN),
    DocSetting("vasp", ("vibrations",), "vasp_incar", "lreal", "LREAL", ".FALSE.", *_VASP_PHONON_SRC,
               "非常に精度の高い計算では逆空間の投影のままにする、と文書が書いています",
               "the text says to keep the reciprocal-space projection scheme for very accurate calculations",
               source_en=_VASP_PHONON_SRC_EN),
]


def documented_settings(spec=None, code: str | None = None, task: str | None = None) -> list[DocSetting]:
    if spec is not None:
        code = code or spec.method.code
        task = task if task is not None else spec.task.type
    return [s for s in DOCUMENTED_SETTINGS if (code is None or s.code == code) and (not s.tasks or task is None or task in s.tasks)]


def format_settings(settings: list[DocSetting]) -> str:
    if not settings:
        return L("この計算コードについて、文書に書かれた設定の控えはありません",
                 "no documented settings are on file for this code")
    head = L("計算コードの公式の文書に、例として書かれている設定 (ADIT の推奨ではありません):",
             "settings written as examples in the official documentation of the code (not recommendations by ADIT):")
    lines = [L(f"  {s.key} = {s.value}   ({s.source}、{s.retrieved} 取得)\n      {s.url}",
               f"  {s.key} = {s.value}   ({s.source_en or s.source}, retrieved {s.retrieved})\n      {s.url}") for s in settings]
    return "\n".join([head, *lines])


def format_values(values: list[DocValue]) -> str:
    if not values:
        return L("パラメータのファイルに、読める値は書かれていませんでした (またはファイルがこの PC にありません)",
                 "no values found in the parameter files (or the files are not on this PC)")
    head = L("パラメータ集の文書に書かれた値 (ADIT は値を決めません。欄に入れるかは利用者が決めます):",
             "values written in the parameter documentation (ADIT does not choose values; you decide whether to use them):")
    return "\n".join([head, *[f"  {v}" for v in values]])


__all__ = ["DOCUMENTED_SETTINGS", "DocSetting", "DocValue", "documented_settings", "documented_values", "format_settings", "format_values"]
