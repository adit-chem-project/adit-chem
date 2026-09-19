from adit.config import Profile
from adit.scripts.render import render_submit
from tests.conftest import water_spec
from adit.spec import Runtime


def test_direct_script():
    p = Profile(kind="direct")
    s = render_submit(water_spec(), p, "dftb+ > output.log 2>&1")
    assert s.startswith("#!/bin/bash\n")
    assert "export OMP_NUM_THREADS=4" in s and "#PBS" not in s
    assert "ulimit" not in s and "module" not in s
    assert s.rstrip().endswith("dftb+ > output.log 2>&1")


def test_pbs_script_is_generic():
    p = Profile(kind="pbs", code_modules={"dftbplus": ["dftbplus/25.1"]}, select_extra=":jobtype=core",
                header_extra=["#PBS -q normal"], submit_command="qsub")
    spec = water_spec(runtime=Runtime(profile="cluster", nodes=1, ncpus=8, mpiprocs=1, omp_threads=8,
                                      walltime="72:00:00", job_name="water"))
    s = render_submit(spec, p, "dftb+ > output.log 2>&1")
    lines = s.splitlines()
    assert lines[0] == "#!/bin/sh"
    assert "#PBS -N water" in lines
    assert "#PBS -l select=1:ncpus=8:mpiprocs=1:ompthreads=8:jobtype=core" in lines
    assert "#PBS -l walltime=72:00:00" in lines and "#PBS -j oe" in lines and "#PBS -q normal" in lines
    assert "module -s purge" in lines
    assert any(l.startswith("module -s load dftbplus/25.1 ||") for l in lines)
    assert "export OMP_NUM_THREADS=8" in lines
    assert ". /etc/profile" in s and "command -v module" in s
    assert s.index(". /etc/profile") < s.index("module -s purge")
    assert "qsub submit.sh" in s
    assert s.rstrip().endswith("dftb+ > output.log 2>&1")
    assert "sample.sh" not in s and "公式サンプル" not in s


def test_pbs_without_extras_or_modules():
    p = Profile(kind="pbs")
    s = render_submit(water_spec(), p, "dftb+")
    assert "jobtype" not in s
    assert "module -s load" not in s and "module -s purge" not in s and "qsub submit.sh" in s


def test_slurm_script():
    p = Profile(kind="slurm", code_modules={"dftbplus": ["dftbplus/25.1"]}, header_extra=["#SBATCH --partition=short"],
                env={"VASP_PP_PATH": "/pp"}, commands={"vasp": "srun vasp_{binary}"})
    spec = water_spec(runtime=Runtime(profile="slurm", nodes=2, ncpus=16, mpiprocs=16, omp_threads=1, walltime="01:00:00", job_name="w"))
    s = render_submit(spec, p, "srun vasp_std > output.log 2>&1")
    lines = s.splitlines()
    assert lines[0] == "#!/bin/bash"
    for key in ["#SBATCH --job-name=w", "#SBATCH --nodes=2", "#SBATCH --ntasks-per-node=16", "#SBATCH --cpus-per-task=1",
                "#SBATCH --time=01:00:00", "#SBATCH --partition=short"]:
        assert key in lines, key
    first_cmd = next(i for i, l in enumerate(lines) if l and not l.startswith("#"))
    assert all(not l.startswith("#SBATCH") for l in lines[first_cmd:])
    assert "sbatch submit.sh" in s and 'export VASP_PP_PATH="/pp"' in s and "module load dftbplus/25.1 ||" in s
    assert s.rstrip().endswith("srun vasp_std > output.log 2>&1")


def test_queue_is_read_from_header_extra():
    from adit.scripts.render import queue_of

    assert queue_of(Profile(kind="pbs", header_extra=["#PBS -q normal"])) == "normal"
    assert queue_of(Profile(kind="slurm", header_extra=["#SBATCH --partition=short"])) == "short"
    assert queue_of(Profile(kind="slurm", header_extra=["#SBATCH -p gpu", "#SBATCH --account=x"])) == "gpu"
    assert queue_of(Profile(kind="pbs")) == ""
