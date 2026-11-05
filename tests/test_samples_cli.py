
import pytest

from adit.cli import main as gen_main
from adit.config import save_config
from adit.samples import copy_sample, examples_dir, list_samples
from tests.conftest import cfg_for


@pytest.mark.skipif(examples_dir() is None, reason="no examples/ (installed with pip)")
def test_samples_are_listed():
    samples = list_samples()
    assert samples and all(s.code and s.task and s.formula for s in samples)
    names = [s.name for s in samples]
    assert "water_generated" in names
    line = next(s for s in samples if s.name == "water_generated").line()
    assert "dftbplus" in line


@pytest.mark.skipif(examples_dir() is None, reason="no examples/")
def test_copy_sample_writes_a_spec(tmp_path):
    out = copy_sample("water_generated", tmp_path / "mine.json")
    assert out.is_file()
    from adit.spec import CalculationSpec

    assert CalculationSpec.load(out).method.code == "dftbplus"
    with pytest.raises(FileExistsError):
        copy_sample("water_generated", out)
    with pytest.raises(FileNotFoundError):
        copy_sample("no_such_sample", tmp_path / "x.json")


@pytest.mark.skipif(examples_dir() is None, reason="no examples/")
def test_cli_lists_and_copies(tmp_path, capsys):
    assert gen_main(["--list-samples"]) == 0
    assert "water_generated" in capsys.readouterr().out
    assert gen_main(["--sample", "water_generated", str(tmp_path / "mine.json")]) == 0
    assert (tmp_path / "mine.json").is_file()


def test_check_config_prints_the_settings(tmp_path, sk_root, capsys):
    path = tmp_path / "cluster.toml"
    save_config(cfg_for(sk_root), path)
    assert gen_main(["--check-config", "--config", str(path)]) == 0
    out = capsys.readouterr().out
    assert str(path) in out and "sk_root" in out and "プロファイル local" in out


def test_check_config_warns_about_unknown_keys(tmp_path, capsys):
    path = tmp_path / "cluster.toml"
    path.write_text('sk_root = "/tmp/sk"\n\n[profiles.local]\nkind = "direct"\n\ntemplates_dir = "/tmp/t"\n', encoding="utf-8")
    assert gen_main(["--check-config", "--config", str(path)]) == 0
    assert "profiles.local.templates_dir" in capsys.readouterr().err
