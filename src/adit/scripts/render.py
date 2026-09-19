
from __future__ import annotations

from jinja2 import Environment, PackageLoader, StrictUndefined

from adit import lang
from adit import __version__
from adit.config import Profile
from adit.spec import CalculationSpec

_env = Environment(
    loader=PackageLoader("adit.scripts", "templates"),
    undefined=StrictUndefined,
    keep_trailing_newline=True,
    trim_blocks=True,
    lstrip_blocks=True,
)

TEMPLATE_OF_KIND = {"direct": "direct.sh.j2", "pbs": "pbs.sh.j2", "slurm": "slurm.sh.j2"}
CHECK_TEMPLATE_OF_KIND = {"pbs": "check_pbs.sh.j2", "slurm": "check_slurm.sh.j2"}


def render_submit(spec: CalculationSpec, profile: Profile, run_command: str) -> str:
    tmpl = _env.get_template(TEMPLATE_OF_KIND[profile.kind])
    return tmpl.render(
        runtime=spec.runtime, profile=profile, run_command=run_command,
        modules=profile.modules_for(spec.method.code), env=profile.env,
        app_version=__version__, created=spec.meta.created, lang=lang.LANGUAGE,
    )


def queue_of(profile: Profile) -> str:
    """The queue / partition named in header_extra, or "" when none is written there."""
    import re

    for line in profile.header_extra:
        m = re.match(r"^\s*#(?:PBS\s+-q|SBATCH\s+(?:--partition=|-p\s+))\s*(\S+)", line)
        if m:
            return m.group(1)
    return ""


def render_check_remote(spec: CalculationSpec, profile: Profile, programs: list[str]) -> str:
    """The pre-submission connectivity check for a cluster profile. It never submits."""
    tmpl = _env.get_template(CHECK_TEMPLATE_OF_KIND[profile.kind])
    return tmpl.render(
        runtime=spec.runtime, profile=profile, modules=profile.modules_for(spec.method.code),
        programs=programs, queue=queue_of(profile),
        target=profile.target or "<user>@<host>", remote_dir=profile.remote_dir.strip() or "<remote_dir>",
        app_version=__version__, created=spec.meta.created, lang=lang.LANGUAGE,
    )
