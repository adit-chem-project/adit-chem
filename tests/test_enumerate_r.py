import json

import pytest

from adit.enumerate_r import EnumerateError, enumerate_substituents, parse_groups
from adit.structures_batch import STRUCTURES_FILE, write_enumerated
from tests.conftest import cfg_for, water_spec

pytest.importorskip("rdkit")


def test_parse_groups_reads_one_list_per_attachment_point():
    assert parse_groups(["1=[H],C,OC", "2=F"]) == {1: ["[H]", "C", "OC"], 2: ["F"]}


def test_every_combination_becomes_a_product_with_a_unique_name():
    products = enumerate_substituents("c1ccc(cc1[*:1])[*:2]", {1: ["[H]", "C"], 2: ["F", "Cl", "OC"]})
    assert len(products) == 6
    assert len({p.name for p in products}) == 6 and len({p.smiles for p in products}) == 6
    assert all(len(p.atoms) >= 11 for p in products)                 # benzene ring plus substituents
    assert products[0].substituents == {1: "[H]", 2: "F"}
    assert "F" in products[0].atoms.get_chemical_symbols()


@pytest.mark.parametrize("core, groups, text", [
    ("c1ccccc1", {1: ["C"]}, "[*:1]"),                                # no attachment point at all
    ("c1ccccc1[*:1]", {}, "1"),                                       # point without substituents
    ("c1ccccc1[*:1]", {1: ["C"], 2: ["C"]}, "2"),                     # substituents for a missing point
])
def test_mismatched_points_are_reported(core, groups, text):
    with pytest.raises(EnumerateError, match=text.replace("[", r"\[").replace("]", r"\]")):
        enumerate_substituents(core, groups)


def test_write_enumerated_makes_one_directory_per_combination(sk_root, tmp_path):
    spec = water_spec()
    dirs = write_enumerated(spec, cfg_for(sk_root), tmp_path / "enum", "C[*:1]", ["1=O,N,CC"])
    assert sorted(d.name for d in dirs) == ["r1-01", "r1-02", "r1-03"]
    record = json.loads((tmp_path / "enum" / STRUCTURES_FILE).read_text(encoding="utf-8"))
    assert record["core"] == "C[*:1]" and [s["smiles"] for s in record["structures"]] == ["CO", "CN", "C(CC)"]
    for d in dirs:
        spec_out = json.loads((d / "spec.json").read_text(encoding="utf-8"))
        assert spec_out["structure"]["source"] == "smiles" and (d / "submit.sh").is_file()
