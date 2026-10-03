import os

import pytest

from adit.cli import main
from adit.config import save_config
from tests.conftest import cfg_for, water_spec


def test_cli_generate_validate_print(sk_root, tmp_path, capsys):
    cfg_path = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), cfg_path)
    spec_path = tmp_path / "spec.json"
    water_spec().save(spec_path)
    assert main([str(spec_path), "--validate", "--config", str(cfg_path)]) == 0
    assert "エラーはありません" in capsys.readouterr().out
    assert main([str(spec_path), "--print", "dftb_in.hsd", "--config", str(cfg_path)]) == 0
    assert 'O = "p"' in capsys.readouterr().out
    out = tmp_path / "calc"
    assert main([str(spec_path), str(out), "--config", str(cfg_path)]) == 0
    assert (out / "submit.sh").is_file() and (out / "skf" / "O-H.skf").is_file()
    assert main([str(spec_path), str(out), "--config", str(cfg_path)]) == 1
    assert main([str(spec_path), str(out), "--config", str(cfg_path), "--overwrite"]) == 0


def test_cli_reports_errors(sk_root, tmp_path, capsys):
    cfg_path = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), cfg_path)
    spec_path = tmp_path / "spec.json"
    from adit.spec import DftbMethod
    water_spec(method=DftbMethod(sk_set="nope-1-1")).save(spec_path)
    assert main([str(spec_path), "--validate", "--config", str(cfg_path)]) == 1
    assert "Slater-Koster パラメータ:" in capsys.readouterr().err
    assert main([str(tmp_path / "missing.json"), "--validate", "--config", str(cfg_path)]) == 2
    capsys.readouterr()
    assert main([str(spec_path), "--validate", "--config", str(tmp_path / "nocfg.toml")]) == 1
    err = capsys.readouterr().err
    assert (tmp_path / "nocfg.toml").is_file() and "環境設定ファイルを作りました" in err and "sk_root" in err
    assert "save_config" not in err and "default_config" not in err


def test_cli_first_run_writes_template(tmp_path, capsys, monkeypatch):
    from adit.spec import XtbMethod
    monkeypatch.setenv("HOME", str(tmp_path / "home")); monkeypatch.delenv("ADIT_CONFIG", raising=False)
    monkeypatch.delenv("QCGUI_CONFIG", raising=False); monkeypatch.delenv("XDG_CONFIG_HOME", raising=False)
    monkeypatch.setenv("APPDATA", str(tmp_path / "home" / "appdata"))
    spec_path = tmp_path / "spec.json"
    water_spec(method=XtbMethod()).save(spec_path)
    assert main([str(spec_path), str(tmp_path / "out")]) == 0
    cfg_file = (tmp_path / "home" / "appdata" / "adit" / "cluster.toml" if os.name == "nt"
                else tmp_path / "home" / ".config" / "adit" / "cluster.toml")
    err = capsys.readouterr().err
    assert cfg_file.is_file() and str(cfg_file) in err and "sk_root" in err


def test_cli_overwrite_message_names_the_option(sk_root, tmp_path, capsys):
    cfg_path = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), cfg_path)
    spec_path = tmp_path / "spec.json"
    water_spec().save(spec_path)
    out = tmp_path / "calc"
    assert main([str(spec_path), str(out), "--config", str(cfg_path)]) == 0
    capsys.readouterr()
    assert main([str(spec_path), str(out), "--config", str(cfg_path)]) == 1
    err = capsys.readouterr().err
    assert "--overwrite" in err and "overwrite=True" not in err


@pytest.mark.parametrize("line,hint", [('sk_root = "C:\\Users\\me\\slakos"', "/mnt/c/Users/"), ("sk_root = /home/x/slakos", '"..." で囲み')])
def test_cli_toml_mistake_is_explained(tmp_path, capsys, line, hint):
    cfg_path = tmp_path / "cluster.toml"
    cfg_path.write_text(line + "\n", encoding="utf-8")
    spec_path = tmp_path / "spec.json"
    water_spec().save(spec_path)
    assert main([str(spec_path), "--validate", "--config", str(cfg_path)]) == 2
    err = capsys.readouterr().err
    assert "1 行目" in err and hint in err and line in err
    assert cfg_path.read_text(encoding="utf-8") == line + "\n"
