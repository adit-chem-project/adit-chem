import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src" / "adit"
ALLOWED = {"config.py", "handoff.py", "readers_extra.py", "continuation.py"}


def test_the_package_is_named_adit():
    assert SRC.is_dir() and not (ROOT / "src" / "vista").exists()
    text = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    assert 'name = "adit-chem"' in text
    for command in ("adit =", "adit-gen =", "adit-analyze =", "adit-web =", "adit-convert =", "adit-report ="):
        assert command in text


def test_no_old_name_is_left_in_the_source():
    hits = []
    for path in sorted(SRC.rglob("*.py")) + sorted(SRC.rglob("*.html")) + sorted(SRC.rglob("*.js")):
        if path.name in ALLOWED:
            continue
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if re.search(r"(?<![A-Za-z])vista|VISTA", line, re.I):
                hits.append(f"{path.relative_to(SRC)}:{n}: {line.strip()[:80]}")
    assert hits == [], "旧名が残っています"


def test_the_old_settings_file_is_still_read(tmp_path, monkeypatch):
    from adit.config import OLD_APP_NAMES, config_path

    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path))
    import types

    from adit import config as _config
    monkeypatch.setattr(_config, "os", types.SimpleNamespace(name="posix", environ=_config.os.environ))
    monkeypatch.delenv("ADIT_CONFIG", raising=False)
    monkeypatch.delenv("VISTA_CONFIG", raising=False)
    monkeypatch.delenv("QCGUI_CONFIG", raising=False)
    for old in OLD_APP_NAMES:
        (tmp_path / old).mkdir()
        (tmp_path / old / "cluster.toml").write_text("", encoding="utf-8")
        assert config_path() == tmp_path / old / "cluster.toml"
        (tmp_path / old / "cluster.toml").unlink()
    (tmp_path / "adit").mkdir()
    (tmp_path / "adit" / "cluster.toml").write_text("", encoding="utf-8")
    assert config_path() == tmp_path / "adit" / "cluster.toml"


def test_the_old_environment_variables_still_work(monkeypatch):
    from adit.config import env_var

    for name in ("ADIT_THEME", "VISTA_THEME", "QCGUI_THEME"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("QCGUI_THEME", "dark")
    assert env_var("THEME") == "dark"
    monkeypatch.setenv("VISTA_THEME", "light")
    assert env_var("THEME") == "light"
    monkeypatch.setenv("ADIT_THEME", "auto")
    assert env_var("THEME") == "auto"


def test_runs_made_before_the_rename_can_still_be_continued(tmp_path):
    from adit.handoff import existing_name

    (tmp_path / "vista-1.restart").write_text("", encoding="utf-8")
    assert existing_name(str(tmp_path), "{}-1.restart") == "vista-1.restart"
    (tmp_path / "adit-1.restart").write_text("", encoding="utf-8")
    assert existing_name(str(tmp_path), "{}-1.restart") == "adit-1.restart"
    assert existing_name("", "{}.gro") == "adit.gro"


def test_the_output_prefix_comes_from_the_job_script(tmp_path):
    from adit.analysis.readers_extra import output_prefix

    (tmp_path / "submit.sh").write_text("gmx mdrun -deffnm mine -ntmpi 1\n", encoding="utf-8")
    assert output_prefix(tmp_path) == "mine"
    (tmp_path / "submit.sh").write_text("gmx mdrun\n", encoding="utf-8")
    (tmp_path / "vista.log").write_text("", encoding="utf-8")
    assert output_prefix(tmp_path) == "vista"
    assert output_prefix(tmp_path / "empty") == "adit"


def test_the_packaging_uses_the_new_name():
    assert (ROOT / "packaging" / "adit.spec").is_file()
    assert not (ROOT / "packaging" / "vista.spec").exists()
    spec = (ROOT / "packaging" / "adit.spec").read_text(encoding="utf-8")
    assert 'name="ADIT"' in spec and "ADIT.app" in spec
    for workflow in ("windows-exe.yml", "macos-app.yml"):
        text = (ROOT / ".github" / "workflows" / workflow).read_text(encoding="utf-8")
        assert "adit.spec" in text and "vista" not in text.lower()


def test_the_font_licence_is_bundled_with_the_font():
    spec = (ROOT / "packaging" / "adit.spec").read_text(encoding="utf-8")
    assert "OFL-1.1-NotoSansCJK.txt" in spec
    assert (ROOT / "licenses" / "OFL-1.1-NotoSansCJK.txt").is_file()
    text = (ROOT / "licenses" / "OFL-1.1-NotoSansCJK.txt").read_text(encoding="utf-8")
    assert "SIL Open Font License" in text
    assert "the font cannot be bundled without the OFL text" in spec
    assert "ONEFILE = False" in spec


def test_the_licence_files_ride_along_in_the_executables():
    spec = (ROOT / "packaging" / "adit.spec").read_text(encoding="utf-8")
    assert "THIRD_PARTY_NOTICES.md" in spec and 'licenses_dir.glob("*.txt")' in spec
    notices = (ROOT / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
    for name in ("GPL-2.0.txt", "LGPL-2.1.txt", "GNU-FDL-1.2.txt", "CC-BY-SA-4.0.txt"):
        assert (ROOT / "licenses" / name).is_file(), name
        assert name in notices, name
