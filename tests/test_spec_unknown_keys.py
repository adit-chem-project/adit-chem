import json

from adit.spec import CalculationSpec, unknown_keys, unknown_keys_in_file
from tests.conftest import water_spec


def _data(**changes) -> dict:
    d = json.loads(water_spec().to_json())
    d.update(changes)
    return d


def test_a_clean_spec_has_no_unknown_keys():
    assert unknown_keys(_data()) == []


def test_typos_are_reported_with_their_path():
    d = _data(kpoint={"mode": "mesh"})
    d["task"]["typ"] = "molecular_dynamics"
    d["task"]["md"]["temprature_k"] = 300
    assert unknown_keys(d) == ["kpoint", "task.md.temprature_k", "task.typ"]


def test_method_fields_follow_the_code_of_the_union():
    d = _data()
    d["method"] = {"code": "vasp", "encut": 500, "ENCUT": 400}
    assert unknown_keys(d) == ["method.ENCUT"]


def test_free_form_dictionaries_and_provenance_are_not_unknown():
    d = _data(provenance={"adit": "0.1.0a1"})
    d["method"] = {"code": "vasp", "extra_incar": {"NCORE": 4}, "hubbard": {"Ti": {"orbital": "d", "u_ev": 3.0}}}
    assert unknown_keys(d) == []
    d["method"]["hubbard"]["Ti"]["bogus"] = 1
    assert unknown_keys(d) == ["method.hubbard.Ti.bogus"]


def test_version_1_keys_are_migrated_before_the_check():
    d = _data(version=1)
    d["task"]["force_tolerance"] = 1e-3
    assert unknown_keys(d) == []


def test_file_helper_accepts_a_bom_and_non_objects(tmp_path):
    p = tmp_path / "spec.json"
    p.write_bytes(b"\xef\xbb\xbf" + json.dumps(_data(kpoint={})).encode("utf-8"))
    assert unknown_keys_in_file(p) == ["kpoint"]
    p.write_text("[1, 2]", encoding="utf-8")
    assert unknown_keys_in_file(p) == []
    assert CalculationSpec.load  # the loader itself is unchanged
