
import pathlib
import re

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
SRC = sorted((ROOT / "src/adit").rglob("*.py")) + sorted((ROOT / "src/adit/web/templates").glob("*.html"))
DOCS = [ROOT / "README.md"] + sorted((ROOT / "docs").glob("*.md"))


def _hits(files, pattern):
    out = []
    for p in files:
        for i, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1):
            if re.search(pattern, line):
                out.append(f"{p.relative_to(ROOT)}:{i}  {line.strip()[:90]}")
    return out


def test_no_fullwidth_parentheses_in_source():
    assert _hits(SRC, r"[（）]") == []


def test_no_fullwidth_parentheses_or_slash_in_docs():
    assert _hits(DOCS, r"[（）／]") == []


@pytest.mark.parametrize("word", ["centring", "centre", "licence", "polarisation", "normalisation"])
def test_american_spelling(word):
    files = SRC + [p for p in DOCS if p.name != "TRIAL_20260913.md"]
    assert _hits(files, rf"\b{word}\b") == []


@pytest.mark.parametrize("word,instead", [("庫", "ライブラリ"), ("圧浴", "圧力浴"), ("参照関数", "参照波動関数")])
def test_withdrawn_words(word, instead):
    assert _hits(SRC, word) == [], f"use {instead}, not {word}"


def test_tool_is_katakana():
    assert _hits(SRC, r'"道具"') == []


def test_superscript_units():
    assert _hits(SRC, r"Å\^[23]|Å\^-3") == []


def test_product_name_is_capitalised_in_english():
    pattern = r"(?<![A-Za-z0-9./_`'\"-])adit (does|is|has|settings|itself|writes|generated)"
    assert _hits(SRC + DOCS, pattern) == []


def test_translation_table_has_no_duplicate_keys():
    import ast
    from collections import Counter

    tree = ast.parse((ROOT / "src/adit/gui/i18n.py").read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Dict) and len(node.keys) > 20:
            keys = [k.value for k in node.keys if isinstance(k, ast.Constant)]
            assert [k for k, n in Counter(keys).items() if n > 1] == []


def test_plain_text_output_has_no_markdown_emphasis():
    import ast

    files = [ROOT / "src/adit/analysis/report.py", ROOT / "src/adit/analysis/transport.py",
             ROOT / "src/adit/analysis/xrd.py", ROOT / "src/adit/analysis/rate.py",
             ROOT / "src/adit/analysis/vanhove.py", ROOT / "src/adit/analysis/heavy_setup.py",
             ROOT / "src/adit/codes/espresso.py"]
    hits = []
    for path in files:
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not (isinstance(node, ast.Call) and getattr(node.func, "id", "") == "L"):
                continue
            for text in _string_parts(node):
                if "**" in text:
                    hits.append(f"{path.relative_to(ROOT)}:{node.lineno}  {text[:70]}")
    assert hits == []


def _string_parts(node):
    import ast

    out = []
    for arg in node.args:
        for sub in ast.walk(arg):
            if isinstance(sub, ast.Constant) and isinstance(sub.value, str):
                out.append(sub.value)
    return out
