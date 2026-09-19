import os
import shutil
import subprocess
from pathlib import Path

import pytest

from adit.config import Config, ConfigError, Profile, config_path, default_config, load_config, save_config
from adit.project import ProjectError, build_project, load_project, write_project
from adit.spec import DftbMethod, Runtime
from tests.conftest import REAL_SK_ROOT, water_spec, pbs_profile

DFTB_EXE = os.environ.get("ADIT_DFTB_EXE") or shutil.which("dftb+") or ""
HAVE_REAL = (REAL_SK_ROOT / "mio-1-1").is_dir() and Path(DFTB_EXE).is_file()


@pytest.fixture
def cfg(sk_root) -> Config:
    c = default_config(sk_root=str(sk_root))
    c.profiles["cluster"] = pbs_profile()
    from tests.conftest import slurm_profile
    c.profiles["slurm"] = slurm_profile()
    return c


def test_config_roundtrip(tmp_path, monkeypatch, cfg):
    monkeypatch.setenv("ADIT_CONFIG", str(tmp_path / "c" / "cluster.toml"))
    assert config_path() == tmp_path / "c" / "cluster.toml"
    with pytest.raises(ConfigError):
        load_config()
    save_config(cfg)
    back = load_config()
    assert back == cfg
    assert back.source_path == cfg.source_path == config_path().absolute()
    assert "source_path" not in back.model_dump_json()
    assert "source_path" not in config_path().read_text(encoding="utf-8")
    assert back.profile("cluster").kind == "pbs" and back.profile("cluster").modules_for("dftbplus") == ["dftbplus/25.1"]
    assert back.profile("cluster").modules_for("vasp") == []
    assert list(default_config().profiles) == ["local"]
    with pytest.raises(ConfigError):
        back.profile("nope")


def test_save_config_keeps_comments_and_unknown_keys(tmp_path):
    path = tmp_path / "cluster.toml"
    path.write_text('# my settings\nsk_root = "/x"  # parameters\nlanguage = "ja"\nmy_note = 1\n\n'
                    '[extra_table]\nkeep = true\n\n[profiles.local]\nkind = "direct"\ntemplates_dir = "/t"\n\n'
                    '[profiles.local.env]\nOMP_STACKSIZE = "1G"\n', encoding="utf-8")
    cfg = load_config(path)
    assert cfg.unknown_keys == ["extra_table", "my_note", "profiles.local.templates_dir"]
    cfg.language = "en"
    cfg.profiles["local"].env = {"OMP_NUM_THREADS": "2"}
    cfg.profiles["cluster"] = pbs_profile()
    save_config(cfg, path)
    text = path.read_text(encoding="utf-8")
    assert text.startswith('# my settings\nsk_root = "/x"  # parameters\nlanguage = "en"\nmy_note = 1\n')
    assert "OMP_STACKSIZE" not in text and 'templates_dir = "/t"' in text
    back = load_config(path)
    assert back.language == "en" and back.profiles["local"].env == {"OMP_NUM_THREADS": "2"}
    assert back.profiles["cluster"] == pbs_profile() and back.unknown_keys == cfg.unknown_keys
    cfg.window_frame = "native"
    save_config(cfg, path)
    text = path.read_text(encoding="utf-8")
    assert text.index('window_frame = "native"') < text.index("[") and load_config(path).window_frame == "native"
    path.write_text("not = [valid\n", encoding="utf-8")
    save_config(cfg, path)
    assert load_config(path).language == "en"


def test_config_broken(tmp_path):
    p = tmp_path / "cluster.toml"; p.write_text("profiles = 3\n", encoding="utf-8")
    with pytest.raises(ConfigError):
        load_config(p)


def test_write_project_layout(tmp_path, cfg):
    out = tmp_path / "calc"
    spec = water_spec()
    written = write_project(spec, cfg, out)
    names = sorted(p.relative_to(out).as_posix() for p in written)
    assert names == sorted(["dftb_in.hsd", "geometry.gen", "submit.sh", "spec.json", "README.txt", "analyze.py",
                            "skf/H-H.skf", "skf/H-O.skf", "skf/O-H.skf", "skf/O-O.skf", "skf/LICENSE", "skf/README"])
    assert os.access(out / "submit.sh", os.X_OK)
    assert load_project(out) == spec
    readme = (out / "README.txt").read_text(encoding="utf-8")
    assert "bash submit.sh" in readme and "fake-1-0" in readme


def test_refuse_overwrite(tmp_path, cfg):
    out = tmp_path / "calc"
    write_project(water_spec(), cfg, out)
    with pytest.raises(ProjectError):
        write_project(water_spec(), cfg, out)
    write_project(water_spec(), cfg, out, overwrite=True)


def test_validation_blocks_writing(tmp_path, cfg):
    out = tmp_path / "calc"
    bad = water_spec(method=DftbMethod(sk_set="nope-9-9"))
    with pytest.raises(ProjectError) as ex:
        write_project(bad, cfg, out)
    assert ex.value.errors and ex.value.errors[0].location == "method.sk_set"
    assert not out.exists()


@pytest.mark.parametrize("name", ["../escape.txt", "a/../../escape.txt", "/tmp/escape.txt", "C:\\escape.txt", "..\\escape.txt"])
def test_file_names_cannot_leave_the_output_directory(tmp_path, cfg, name):
    with pytest.raises(ProjectError, match="escape"):
        build_project(water_spec(), cfg, output_dir=tmp_path / "calc", extra_texts={name: "x"})
    from adit.spec import Handoff
    (tmp_path / "prev").mkdir()
    (tmp_path / "prev" / "charges.bin").write_bytes(b"x")
    h = Handoff(previous_dir=str(tmp_path / "prev"), previous_task="single_point", previous_code="dftbplus",
                files={name: "charges.bin"})
    with pytest.raises(ProjectError, match="escape"):
        write_project(water_spec(handoff=h), cfg, tmp_path / "calc")
    assert not (tmp_path / "escape.txt").exists() and not (tmp_path / "calc").exists()


def test_output_path_that_is_a_file_is_reported(tmp_path, cfg):
    (tmp_path / "calc").write_text("not a directory\n", encoding="utf-8")
    with pytest.raises(ProjectError) as ex:
        write_project(water_spec(), cfg, tmp_path / "calc")
    assert [e.location for e in ex.value.errors] == ["output_dir"]
    assert (tmp_path / "calc").read_text(encoding="utf-8") == "not a directory\n"


def test_unknown_profile(tmp_path, cfg):
    with pytest.raises(ProjectError) as ex:
        build_project(water_spec(runtime=Runtime(profile="mars")), cfg)
    assert ex.value.errors[0].location == "runtime.profile"


def test_cluster_readme_has_manual_steps(cfg):
    spec = water_spec(runtime=Runtime(profile="cluster", ncpus=8, omp_threads=8, walltime="01:00:00", job_name="w"))
    files = build_project(spec, cfg)
    r = files.texts["README.txt"]
    assert "scp -r" in r and "qsub submit.sh" in r and "qstat" in r
    assert "#PBS -l select=1:ncpus=8:mpiprocs=1:ompthreads=8:jobtype=core" in files.texts["submit.sh"]
    slurm = build_project(water_spec(runtime=Runtime(profile="slurm", ncpus=8, mpiprocs=8, job_name="w")), cfg)
    assert "#SBATCH --job-name=w" in slurm.texts["submit.sh"] and "sbatch submit.sh" in slurm.texts["README.txt"]


@pytest.mark.skipif(not HAVE_REAL, reason="mio-1-1 か dftb+ が無い")
def test_local_submit_runs(tmp_path):
    cfg = default_config(sk_root=str(REAL_SK_ROOT))
    out = tmp_path / "water"
    write_project(water_spec(method=DftbMethod(sk_set="mio-1-1")), cfg, out)
    env = {**os.environ, "PATH": f"{Path(DFTB_EXE).parent}:{os.environ.get('PATH', '')}"}
    r = subprocess.run(["bash", "submit.sh"], cwd=out, capture_output=True, text=True, timeout=300, env=env)
    assert r.returncode == 0, r.stderr
    log = (out / "output.log").read_text(encoding="utf-8")
    assert "Geometry converged" in log and (out / "results.tag").is_file()


def test_the_transfer_and_submit_commands_are_written_for_a_cluster(tmp_path, cfg):
    from adit.spec import Runtime

    out = tmp_path / "calc"
    spec = water_spec(runtime=Runtime(profile="cluster", ncpus=8, walltime="01:00:00", job_name="w"))
    plain = build_project(spec, cfg, output_dir=out).texts["transfer_and_submit.sh"]
    assert "<ユーザー名>@<クラスタのホスト名>" in plain and "qsub submit.sh" in plain
    cfg.profiles["cluster"].host = "cluster.example.ac.jp"
    cfg.profiles["cluster"].user = "me"
    cfg.profiles["cluster"].remote_dir = "/work/me"
    filled = build_project(spec, cfg, output_dir=out).texts["transfer_and_submit.sh"]
    assert 'rsync -av "$HERE" me@cluster.example.ac.jp:/work/me/' in filled
    assert "cd /work/me/calc && qsub submit.sh" in filled and "<" not in filled
    assert "transfer_and_submit.sh" not in build_project(water_spec(), cfg, output_dir=out).texts   # local profile


def test_check_remote_script_is_written_for_a_cluster_and_never_submits(tmp_path, cfg):
    import subprocess

    from adit.spec import Runtime

    out = tmp_path / "calc"
    spec = water_spec(runtime=Runtime(profile="cluster", ncpus=8, walltime="01:00:00", job_name="w"))
    assert "check_remote.sh" not in build_project(water_spec(), cfg, output_dir=out).texts   # local profile
    plain = build_project(spec, cfg, output_dir=out).texts["check_remote.sh"]
    assert 'TARGET="<user>@<host>"' in plain and "*'<'*)" in plain
    for name, profile in (("cluster", "pbs"), ("slurm", "slurm")):
        p = cfg.profiles[name]
        p.host, p.user, p.remote_dir = "cluster.example.ac.jp", "me", "/work/me"
        p.commands = {"dftbplus": "mpirun -np {mpiprocs} dftb+"}
        files = build_project(water_spec(runtime=Runtime(profile=name, ncpus=8, job_name="w")), cfg, output_dir=out)
        s = files.texts["check_remote.sh"]
        assert "ssh -o BatchMode=yes -o ConnectTimeout=10" in s and 'TARGET="me@cluster.example.ac.jp"' in s
        assert 'REMOTE_DIR="/work/me"' in s and "mkdir -p '$REMOTE_DIR' && touch" in s
        assert "command -v dftb+" in s and "command -v mpirun" in s and "load dftbplus/25.1 || exit 1" in s
        if profile == "pbs":
            assert 'QUEUE="normal"' in s and "qstat -Q $QUEUE" in s and "module -s load" in s and "sbatch" not in s
        else:
            assert "sbatch --test-only" in s and "qsub" not in s and "module load dftbplus/25.1" in s
        assert "qsub submit.sh" not in s and "sbatch submit.sh" not in s
        (tmp_path / f"check_{profile}.sh").write_text(s, encoding="utf-8")
        if shutil.which("bash"):
            r = subprocess.run(["bash", "-n", str(tmp_path / f"check_{profile}.sh")], capture_output=True, text=True)
            assert r.returncode == 0, r.stderr
        assert "check_remote.sh" in files.texts["README.txt"]
