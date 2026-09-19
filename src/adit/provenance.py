
from __future__ import annotations

import hashlib
import platform
from datetime import datetime, timezone
from pathlib import Path

from adit import __version__
from adit.lang import L

VERSION_FILE = "code_version.txt"
_CACHE: dict[tuple[str, int, int], str] = {}


def sha256_file(path: Path | str) -> str:
    p = Path(path)
    st = p.stat()
    key = (str(p.resolve()), st.st_size, st.st_mtime_ns)
    if key not in _CACHE:
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        _CACHE[key] = h.hexdigest()
    return _CACHE[key]


def sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def collect(copies: dict[str, Path], texts: dict[str, str], sources: dict[str, Path], probe: tuple[str, str] | None,
            replaced_at_run: set[str] | None = None) -> dict:
    import ase

    replaced_at_run = set(replaced_at_run or ())
    files = []
    for rel, src in sorted(copies.items()):
        p = Path(src)
        files.append({"name": rel, "sha256": sha256_file(p), "bytes": p.stat().st_size, "source": str(p)})
    generated = []
    for rel, src in sorted(sources.items()):
        if rel in texts and src is not None and Path(src).is_file():
            generated.append({"name": rel, "sha256": sha256_text(texts[rel]), "copied_from": str(src), "source_sha256": sha256_file(src)})
    inputs = [{"name": rel, "sha256": sha256_text(text)} for rel, text in sorted(texts.items())]
    for item in inputs:
        if item["name"] in replaced_at_run:
            item["replaced_at_run"] = True
    for item in files:
        if item["name"] in replaced_at_run:
            item["replaced_at_run"] = True
    return {
        "adit_version": __version__,
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "python": platform.python_version(),
        "ase": ase.__version__,
        "files": files,
        "generated_files": generated,
        "inputs": inputs,
        "code_version_file": VERSION_FILE if probe else None,
    }


def spec_json_with(spec, prov: dict) -> str:
    import json

    data = spec.model_dump(mode="json")
    data["provenance"] = prov
    return json.dumps(data, indent=2, ensure_ascii=False) + "\n"


def read_provenance(output_dir: Path | str) -> dict | None:
    import json

    p = Path(output_dir).expanduser() / "spec.json"
    try:
        return json.loads(p.read_text(encoding="utf-8")).get("provenance")
    except (OSError, ValueError, AttributeError):
        return None


def readme_lines(prov: dict) -> list[str]:
    lines = [L("== 作成時の記録 ==", "== Provenance =="),
             L(f"  ADIT {prov['adit_version']} / Python {prov['python']} / ASE {prov['ase']} / 生成 {prov['generated_utc']} (UTC)",
               f"  ADIT {prov['adit_version']} / Python {prov['python']} / ASE {prov['ase']} / generated {prov['generated_utc']} (UTC)")]
    if prov["files"] or prov["generated_files"]:
        lines.append(L("  パラメータなどのファイルの SHA-256 (中身から計算する照合用のハッシュ。同じ値なら同じファイル。spec.json の provenance にも同じもの):",
                       "  SHA-256 of parameter and other files (a fingerprint of the contents; equal values mean identical files; also under provenance in spec.json):"))
        for f in prov["files"]:
            lines.append(f"    {f['sha256']}  {f['name']}")
        for f in prov["generated_files"]:
            lines.append(f"    {f['sha256']}  {f['name']}  " + L(f"(元のファイル {f['copied_from']} の SHA-256: {f['source_sha256']})",
                                                                 f"(SHA-256 of the source {f['copied_from']}: {f['source_sha256']})"))
    if prov.get("code_version_file"):
        lines.append(L(f"  計算コードのバージョン: 走り終えると submit.sh が出力のバージョンの行を {prov['code_version_file']} に写します",
                       f"  Code version: when the run ends, submit.sh copies the version line of the output to {prov['code_version_file']}"))
    else:
        lines.append(L("  計算コードのバージョン: この計算コードの出力からバージョンの行を拾う方法を確かめていないので、記録しません",
                       "  Code version: not recorded (how to pick the version line from this code's output has not been verified)"))
    lines.append("")
    if prov.get("fetched_structure"):
        from adit.fetch import readme_lines as fetched_lines

        lines += fetched_lines(prov["fetched_structure"])
    return lines


def verify_inputs(run_dir: Path | str) -> list[dict]:
    d = Path(run_dir).expanduser()
    prov = read_provenance(d) or {}
    out: list[dict] = []
    for key, kind in (("inputs", "input"), ("files", "copied"), ("generated_files", "written")):
        for item in prov.get(key) or []:
            name, recorded = item.get("name", ""), item.get("sha256", "")
            if not name or not recorded:
                continue
            path = d / name
            if not path.is_file():
                out.append({"name": name, "kind": kind, "state": "missing", "recorded": recorded, "found": ""})
                continue
            found = sha256_file(path)
            if found == recorded:
                state = "same"
            elif item.get("replaced_at_run"):
                state = "replaced_at_run"
            else:
                state = "changed"
            out.append({"name": name, "kind": kind, "state": state, "recorded": recorded, "found": found})
    return out
