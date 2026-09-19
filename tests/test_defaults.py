from __future__ import annotations

from adit.defaults import MISSING, changed_fields, default_at, format_value, value_at, with_default
from adit.spec import DftbMethod, KPoints, MDSettings, Task, VaspMethod
from tests.conftest import water_spec


def test_untouched_spec_has_no_changed_fields():
    spec = water_spec(task=Task())
    assert changed_fields(spec) == {}
    assert set(changed_fields(water_spec())) == {"種類", "最大ステップ数 (MaxSteps)"}     # the test helper's own choices; runtime is not marked
    assert default_at(spec, "task.md.steps") == (True, 1000) and default_at(spec, "method.sk_set")[0] is False
    assert value_at(spec, "kpoints.mode") is MISSING          # a molecule has no k-points block


def test_changed_fields_are_keyed_by_the_screen_label():
    spec = water_spec(task=Task(type="molecular_dynamics", md=MDSettings(steps=5000, temperature_k=300.0)),
                      method=DftbMethod(sk_set="fake-1-0", scc_tolerance=1e-6, third_order=True))
    ch = changed_fields(spec)
    assert set(ch) == {"種類", "MD ステップ数", "SCC の収束判定 (SccTolerance)", "DFTB3 (ThirdOrderFull)"}
    assert ch["MD ステップ数"].paths == ("task.md.steps",) and ch["MD ステップ数"].value == 5000 and ch["MD ステップ数"].default == 1000
    assert ch["種類"].default == "single_point"
    spec2 = spec.model_copy(update={"kpoints": KPoints(mode="mesh", mesh=(4, 4, 4))})
    ch2 = changed_fields(spec2)
    assert ch2["メッシュ (n1 n2 n3)"].paths == ("kpoints.mesh", "kpoints.shift") and ch2["サンプリング方法"].value == "mesh"
    vasp = water_spec(task=Task(), method=VaspMethod(encut=520.0, ediff=1e-4))
    assert set(changed_fields(vasp)) == {"ENCUT [eV] (0 = 指定しない)"}      # ediff equals its default


def test_with_default_puts_one_field_back():
    t = Task(type="molecular_dynamics", md=MDSettings(steps=5000, timestep_fs=0.5), max_steps=999)
    back = with_default(t, ["md", "steps"])
    assert back.md.steps == 1000 and back.md.timestep_fs == 0.5 and back.max_steps == 999 and back.type == "molecular_dynamics"
    assert with_default(t, ["max_steps"]).max_steps == 200
    m = with_default(DftbMethod(sk_set="x", dispersion="dftd3", d3_params={"s6": 1.0}), ["d3_params"])
    assert m.d3_params is None and m.dispersion == "dftd3"


def test_format_value_reads_naturally():
    assert format_value(True) == "オン" and format_value(None) == "なし" and format_value("") == "空欄"
    assert format_value(1e-5) == "1e-05" and format_value((1, 1, 1)) == "1 1 1" and format_value({"s6": 1.0}) == "s6 = 1"
