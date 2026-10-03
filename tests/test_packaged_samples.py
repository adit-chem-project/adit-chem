from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PACKAGED = REPO / "src" / "adit" / "examples"
TOP = REPO / "examples"


def test_every_packaged_sample_is_a_copy_of_the_repository_example():
    # pip installs ship src/adit/examples/*/spec.json; they must not drift from examples/*/spec.json
    names = sorted(p.parent.name for p in PACKAGED.glob("*/spec.json"))
    assert len(names) >= 20
    stale = [n for n in names if (PACKAGED / n / "spec.json").read_bytes() != (TOP / n / "spec.json").read_bytes()]
    assert stale == [], f"copy examples/<name>/spec.json over src/adit/examples/<name>/spec.json: {stale}"
    for n in names:
        src = PACKAGED / n / "SOURCES.md"
        if src.is_file():
            assert src.read_bytes() == (TOP / n / "SOURCES.md").read_bytes(), n
