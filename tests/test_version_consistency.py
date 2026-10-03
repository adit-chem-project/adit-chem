import re
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent


def _cff_version(path: Path) -> str:
    m = re.search(r"^version:\s*(\S+)\s*$", path.read_text(encoding="utf-8"), re.M)
    assert m, path
    return m.group(1).strip('"')


def test_the_version_is_the_same_everywhere():
    pyproject = re.search(r'^version = "([^"]+)"', (REPO / "pyproject.toml").read_text(encoding="utf-8"), re.M).group(1)
    fallback = re.search(r'__version__ = "([^"]+)"', (REPO / "src" / "adit" / "__init__.py").read_text(encoding="utf-8")).group(1)
    assert pyproject == fallback == _cff_version(REPO / "CITATION.cff") == _cff_version(REPO / "src" / "adit" / "CITATION.cff")
