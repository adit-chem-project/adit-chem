
from __future__ import annotations

import itertools
import sys
from pathlib import Path
from typing import TYPE_CHECKING

import pytest

from adit.spec import HARTREE_PER_BOHR_IN_EV_PER_ANG, AtomsData, CalculationSpec, DftbMethod, Runtime, Structure, Task

if TYPE_CHECKING:
    from adit.config import Config

REPO = Path(__file__).resolve().parent.parent
REAL_SK_ROOT = REPO / "slakos"

FAKE_SHELLS = {"H": "1s", "O": "2s 2p", "C": "2s 2p", "N": "2s 2p", "Ti": "3d 4s 4p", "Si": "3s 3p", "Al": "3s 3p"}


@pytest.fixture(autouse=True)
def _isolated_config(tmp_path_factory, monkeypatch):
    path = tmp_path_factory.mktemp("adit_config") / "cluster.toml"
    monkeypatch.setenv("ADIT_CONFIG", str(path))
    yield path


@pytest.fixture(autouse=True)
def _language_is_japanese_for_every_test():
    # The suite asserts Japanese text; a test that switches the language must not leak into the next one.
    from adit import lang

    i18n = sys.modules.get("adit.gui.i18n")
    saved_lang, saved_i18n = lang.LANGUAGE, getattr(i18n, "LANGUAGE", None)
    lang.set_language("ja")
    if i18n is not None:
        i18n.set_language("ja")
    yield
    i18n = sys.modules.get("adit.gui.i18n")
    if i18n is not None:
        i18n.set_language(saved_i18n or "ja")   # also sets adit.lang
    lang.LANGUAGE = saved_lang


@pytest.fixture(autouse=True)
def _delete_windows_after_test():
    """Release test widgets before their styles and signal references accumulate.

    Existing module-local fixtures override this fallback. Do not import Qt for
    non-GUI tests or create an application solely for cleanup.
    """
    yield
    widgets = sys.modules.get("PySide6.QtWidgets")
    if widgets is None:
        return
    app = widgets.QApplication.instance()
    if app is None:
        return
    from PySide6.QtCore import QCoreApplication, QEvent

    for window in app.topLevelWidgets():
        if window.parentWidget() is None:
            window.close()
            window.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


def make_fake_skset(root: Path, name: str, elements: list[str], *, docs: bool = True) -> Path:
    d = root / name
    d.mkdir(parents=True)
    for a, b in itertools.product(elements, repeat=2):
        body = "0.02 500\n"
        if a == b:
            body += f"<Documentation>\n  <Basis atom=\"1\">\n    <Shells>{FAKE_SHELLS[a]} </Shells>\n  </Basis>\n</Documentation>\n"
        (d / f"{a}-{b}.skf").write_text(body, encoding="utf-8")
    if docs:
        (d / "LICENSE").write_text("CC-BY-SA-4.0 (fake for tests)\n", encoding="utf-8")
        (d / "README").write_text("fake set for tests\n", encoding="utf-8")
        (d / "spinw.txt").write_text("".join(f"{e}:\n   -0.05\n\n" for e in elements), encoding="utf-8")
    return d


@pytest.fixture
def sk_root(tmp_path: Path) -> Path:
    root = tmp_path / "slakos"
    make_fake_skset(root, "fake-1-0", ["H", "C", "N", "O"])
    make_fake_skset(root, "nodocs-0-1", ["H", "O"], docs=False)
    return root


def cfg_for(sk_root) -> "Config":
    from adit.config import default_config
    c = default_config(sk_root=str(sk_root) if sk_root else "")
    c.profiles["cluster"] = pbs_profile()
    c.profiles["slurm"] = slurm_profile()
    return c


def pbs_profile():
    from adit.config import Profile
    return Profile(kind="pbs", code_modules={"dftbplus": ["dftbplus/25.1"]}, select_extra=":jobtype=core",
                   header_extra=["#PBS -q normal"], description="テスト用の PBS")


def slurm_profile():
    from adit.config import Profile
    return Profile(kind="slurm", code_modules={"dftbplus": ["dftbplus/25.1"]}, header_extra=["#SBATCH --partition=short"],
                   description="テスト用の Slurm")


def water_spec(**overrides) -> CalculationSpec:
    kw = dict(
        structure=Structure(
            source="preset",
            source_ref="H2O",
            atoms=AtomsData(
                symbols=["O", "H", "H"],
                positions=[(0.0, -1.0, 0.0), (0.0, 0.0, 0.783064), (0.0, 0.0, -0.783064)],
            ),
        ),
        method=DftbMethod(sk_set="fake-1-0"),
        task=Task(type="geometry_optimization", optimizer="Rational", max_steps=100, force_tolerance_ev_per_ang=1e-4 * HARTREE_PER_BOHR_IN_EV_PER_ANG),
        runtime=Runtime(profile="local", omp_threads=4, job_name="water"),
    )
    kw.update(overrides)
    return CalculationSpec(**kw)

def run_log_tail(directory, lines: int = 25) -> str:
    from pathlib import Path as _Path

    out = []
    for name in ("output.log", "run.log", "grompp.log"):
        p = _Path(directory) / name
        if p.is_file():
            text = p.read_text(encoding="utf-8", errors="replace").splitlines()
            out.append(f"--- {name} の末尾 {min(lines, len(text))} 行 ---")
            out += text[-lines:]
    return "\n".join(out) or "(output.log がありません)"
