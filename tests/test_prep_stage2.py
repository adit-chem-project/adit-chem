
import json
from pathlib import Path

import pytest
from ase.io import write

from tests.conftest import cfg_for, make_fake_skset, water_spec
from adit.compare_sets import CompareSetError, write_compare_set
from adit.spec import (AtomsData, CalculationSpec, Cp2kMethod, DftbMethod, EspressoMethod, KPoints, Structure, Task, VaspMethod,
                        XtbMethod)


def _slab_files(tmp: Path):
    from ase.build import graphene, molecule

    slab = graphene(size=(2, 2, 1), vacuum=7.5)
    slab.pbc = True
    mol = molecule("H2O")
    ads = slab.copy()
    m2 = mol.copy()
    m2.translate(slab.positions.mean(axis=0) + [0, 0, 3.0] - m2.positions.mean(axis=0))
    ads += m2
    write(tmp / "slab.extxyz", slab)
    write(tmp / "h2o.xyz", mol)
    write(tmp / "ads.extxyz", ads)
    return slab


def _slab_spec(slab) -> CalculationSpec:
    return water_spec(structure=Structure(source="file", source_ref="slab.extxyz", atoms=AtomsData.from_ase(slab)),
                      kpoints=KPoints(mode="mesh", mesh=(2, 2, 1)), task=Task(type="single_point"))


ADS = {"kind": "adsorption", "slab": {"structure": "slab.extxyz"}, "molecule": {"structure": "h2o.xyz"}, "adsorbed": {"structure": "ads.extxyz"}}


def test_adsorption_set_writes_compare_json_and_records_the_nonperiodic_molecule(sk_root, tmp_path):
    from adit.analysis.compare import load_compare_file

    slab = _slab_files(tmp_path)
    out = tmp_path / "out"
    dirs = write_compare_set(_slab_spec(slab), cfg_for(sk_root), out, ADS, tmp_path)
    assert [d.name for d in dirs] == ["slab", "molecule", "adsorbed"]
    js = json.loads((out / "compare.json").read_text(encoding="utf-8"))
    assert js["reactions"] == [{"name": "adsorption", "terms": [{"dir": "adsorbed", "nu": 1}, {"dir": "slab", "nu": -1}, {"dir": "molecule", "nu": -1}]}]
    rx = load_compare_file(out / "compare.json")
    assert rx[0].terms == [(1.0, "adsorbed"), (-1.0, "slab"), (-1.0, "molecule")]
    per = next(x for x in js["differences"] if x["setting"] == "structure.periodic")
    assert per["values"] == {"slab": True, "molecule": False, "adsorbed": True} and per["intended"] is False
    assert any(x["setting"].startswith("kpoints.") and x["missing_in"] == ["molecule"] for x in js["partial"])
    assert "KPointsAndWeights" not in (out / "molecule" / "dftb_in.hsd").read_text(encoding="utf-8")
    assert "KPointsAndWeights" in (out / "slab" / "dftb_in.hsd").read_text(encoding="utf-8")
    assert js["balance"] == [{"name": "adsorption", "imbalance": {}, "charge_imbalance": 0.0}]
    top = (out / "README.txt").read_text(encoding="utf-8")
    assert "structure.periodic" in top and "adit-analyze" in top and (out / "submit.sh").is_file()
    text = (out / "molecule" / "README.txt").read_text(encoding="utf-8")
    assert "比べる計算の組" in text or "Set of runs" in text


def test_adsorption_molecule_in_slab_cell(sk_root, tmp_path):
    slab = _slab_files(tmp_path)
    out = tmp_path / "out"
    write_compare_set(_slab_spec(slab), cfg_for(sk_root), out, dict(ADS, molecule_box="slab_cell"), tmp_path)
    js = json.loads((out / "compare.json").read_text(encoding="utf-8"))
    assert all(m["periodic"] for m in js["members"])
    assert not any(x["setting"] == "structure.periodic" for x in js["differences"])
    mol = CalculationSpec.load(out / "molecule" / "spec.json").structure.atoms.to_ase()
    assert abs(mol.cell.array - slab.cell.array).max() < 1e-9


def test_adsorption_rejects_three_identical_base_structures(sk_root, tmp_path):
    from adit.compare_sets import CompareSetError, plan_set
    data = {"kind": "adsorption", "slab": {}, "molecule": {}, "adsorbed": {}}
    with pytest.raises(CompareSetError, match="すべて同じ|identical"):
        plan_set(water_spec(), data, tmp_path)


def test_reaction_set_records_imbalance(sk_root, tmp_path):
    from ase.build import molecule

    for n in ("H2", "O2", "H2O"):
        write(tmp_path / f"{n}.xyz", molecule(n))
    base = water_spec(task=Task(type="single_point"))
    members = [{"name": "h2", "structure": "H2.xyz", "nu": -1}, {"name": "o2", "structure": "O2.xyz", "nu": -0.5}, {"name": "h2o", "structure": "H2O.xyz", "nu": 1}]
    write_compare_set(base, cfg_for(sk_root), tmp_path / "ok", {"kind": "reaction", "name": "water", "members": members}, tmp_path)
    js = json.loads((tmp_path / "ok" / "compare.json").read_text(encoding="utf-8"))
    assert js["balance"][0]["imbalance"] == {} and js["reactions"][0]["name"] == "water"
    members[1]["nu"] = -1
    write_compare_set(base, cfg_for(sk_root), tmp_path / "bad", {"kind": "reaction", "members": members}, tmp_path)
    js = json.loads((tmp_path / "bad" / "compare.json").read_text(encoding="utf-8"))
    assert js["balance"][0]["imbalance"] == {"O": -1.0}
    readme = (tmp_path / "bad" / "README.txt").read_text(encoding="utf-8")
    assert "釣り合っていません" in readme or "not balanced" in readme


def test_solvation_set_marks_the_intended_difference(sk_root, tmp_path):
    base = water_spec(method=XtbMethod(), task=Task(type="single_point"))
    out = tmp_path / "out"
    write_compare_set(base, cfg_for(sk_root), out, {"kind": "solvation", "solvent": {"method": {"solvation": "alpb", "solvent": "water"}}}, tmp_path)
    assert "--alpb water" in (out / "solvent" / "submit.sh").read_text(encoding="utf-8")
    assert "--alpb" not in (out / "gas" / "submit.sh").read_text(encoding="utf-8")
    js = json.loads((out / "compare.json").read_text(encoding="utf-8"))
    diff = {x["setting"]: x for x in js["differences"]}
    assert diff["method.solvation"]["intended"] and diff["method.solvent"]["intended"]
    assert js["reactions"][0]["terms"] == [{"dir": "solvent", "nu": 1}, {"dir": "gas", "nu": -1}]


@pytest.mark.parametrize("data,words", [
    ({"kind": "nope"}, "kind"),
    ({"kind": "adsorption", "slab": {}, "molecule": {}}, "adsorbed"),
    ({"kind": "solvation", "solvent": {"method": {}}}, "method"),
    ({"kind": "solvation", "solvent": {"method": {"code": "orca"}}}, "コード"),
    ({"kind": "reaction", "members": [{"name": "a", "nu": 0}, {"name": "b", "nu": 1}]}, "nu"),
    ({"kind": "reaction", "members": [{"name": "a b", "nu": 1}, {"name": "b", "nu": 1}]}, "name"),
])
def test_set_errors_are_mechanical(sk_root, tmp_path, data, words):
    base = water_spec(method=XtbMethod(), task=Task(type="single_point"))
    with pytest.raises(CompareSetError) as ei:
        write_compare_set(base, cfg_for(sk_root), tmp_path / "o", data, tmp_path)
    assert words in str(ei.value)
    assert not (tmp_path / "o").exists()


def _ethanol_spec(**kw):
    from adit.structure import from_smiles
    return water_spec(structure=Structure(source="smiles", source_ref="CCO", atoms=AtomsData.from_ase(from_smiles("CCO")), **kw),
                      method=XtbMethod(), task=Task(type="geometry_optimization", max_steps=100))


def test_conformers_ethanol(sk_root, tmp_path):
    pytest.importorskip("rdkit")
    from adit.conformers import write_conformers

    out = tmp_path / "conf"
    dirs = write_conformers(_ethanol_spec(), cfg_for(sk_root), out, n=10, rmsd=0.3, heavy_only=False)
    js = json.loads((out / "conformers.json").read_text(encoding="utf-8"))
    assert js["embedded"] == 10 and js["kept"] == len(dirs) >= 2 and js["rmsd_threshold_ang"] == 0.3 and js["seed"] == 12345
    assert len({round(c["energy_kcal_mol"], 4) for c in js["conformers"]}) >= 2
    heavy = write_conformers(_ethanol_spec(), cfg_for(sk_root), tmp_path / "heavy", n=10, rmsd=0.3)
    assert len(heavy) == 1
    kept = [c for c in js["conformers"] if c["kept"]]
    assert [c["dir"] for c in kept] == [d.name for d in dirs] and dirs[0].name == "conf_001"
    assert [c["energy_kcal_mol"] for c in kept] == sorted(c["energy_kcal_mol"] for c in kept)
    assert all(c["duplicate_of_rdkit_id"] is not None for c in js["conformers"] if not c["kept"])
    scan = json.loads((out / "scan.json").read_text(encoding="utf-8"))
    assert scan["path"] == "conformer" and scan["dirs"] == [d.name for d in dirs]
    assert "crest struct.xyz --gfn2 --chrg 0 --uhf 0 -T 4" in (out / "crest" / "run_crest.sh").read_text(encoding="utf-8")
    assert (out / "conformers.sdf").read_text(encoding="utf-8").count("$$$$") == len(dirs)
    one = write_conformers(_ethanol_spec(), cfg_for(sk_root), tmp_path / "one", n=10, rmsd=100.0)
    assert len(one) == 1


def test_conformer_checks(sk_root, tmp_path):
    pytest.importorskip("rdkit")
    from adit.conformers import ConformerError, write_conformers

    with pytest.raises(ConformerError, match="charge|電荷"):
        write_conformers(_ethanol_spec(charge=1), cfg_for(sk_root), tmp_path / "a", n=3, rmsd=0.5)
    file_spec = water_spec(method=XtbMethod())
    with pytest.raises(ConformerError):
        write_conformers(file_spec.model_copy(update={"structure": file_spec.structure.model_copy(update={"source": "file", "source_ref": "w.xyz"})}),
                         cfg_for(sk_root), tmp_path / "b", n=3, rmsd=0.5)
    with pytest.raises(ConformerError):
        write_conformers(_ethanol_spec(), cfg_for(sk_root), tmp_path / "c", n=0, rmsd=0.5)
    with pytest.raises(ConformerError, match="同じ座標|identical"):
        write_conformers(_ethanol_spec(), cfg_for(sk_root), tmp_path / "s0", n=5, rmsd=0.5, seed=0)
    assert not (tmp_path / "s0").exists()
    dftb = _ethanol_spec().model_copy(update={"method": DftbMethod(sk_set="fake-1-0")})
    write_conformers(dftb, cfg_for(sk_root), tmp_path / "d", n=3, rmsd=0.5)
    assert not (tmp_path / "d" / "crest").exists()


def test_documented_values_upf_and_sssp(sk_root, tmp_path):
    from adit.docvalues import documented_values

    d = tmp_path / "pseudo" / "p"
    d.mkdir(parents=True)
    (d / "Si.pbe-x.UPF").write_text('<UPF version="2.0.1">\n<PP_INFO>\n  Suggested minimum cutoff for wavefunctions:  40. Ry\n'
                                    '  Suggested minimum cutoff for charge density: 320. Ry\n</PP_INFO>\n<PP_HEADER\n   element="Si"\n   z_valence="4.0"\n'
                                    '   wfc_cutoff="4.4e1"\n   rho_cutoff="0.0"/>\n<PP_MESH>\n</UPF>\n', encoding="utf-8")
    (d / "SSSP_test.json").write_text(json.dumps({"Si": {"filename": "Si.pbe-x.UPF", "cutoff_wfc": 30.0, "cutoff_rho": 240.0}}, indent=1), encoding="utf-8")
    (tmp_path / "pseudo" / "other.json").write_text(json.dumps({"Si": {"filename": "Si.other.UPF", "cutoff_wfc": 99.0}}), encoding="utf-8")
    from ase.build import bulk
    spec = water_spec(structure=Structure(source="bulk", source_ref="Si", atoms=AtomsData.from_ase(bulk("Si"))),
                      method=EspressoMethod(pseudo_set="p", ecutwfc=30), kpoints=KPoints(mode="mesh", mesh=(2, 2, 2)), task=Task())
    cfg = cfg_for(sk_root)
    cfg.pseudo_root = str(tmp_path / "pseudo")
    v = {x.quantity: x for x in documented_values(spec, cfg)}
    assert v["suggested_wfc_cutoff"].value == 40.0 and v["suggested_wfc_cutoff"].line == 3
    assert v["suggested_rho_cutoff"].value == 320.0
    assert v["wfc_cutoff"].value == 44.0 and v["wfc_cutoff"].line == 9
    assert "rho_cutoff" not in v
    assert v["sssp_cutoff_wfc"].value == 30.0 and v["sssp_cutoff_wfc"].source.endswith("SSSP_test.json") and v["sssp_cutoff_wfc"].line == 2
    assert not any(x.value == 99.0 for x in documented_values(spec, cfg))


def test_documented_values_potcar_sk_cp2k_and_none(sk_root, tmp_path):
    from ase.build import bulk, molecule

    from adit.docvalues import documented_values

    pp = tmp_path / "pp" / "potpaw_PBE" / "Si"
    pp.mkdir(parents=True)
    (pp / "POTCAR").write_text("  PAW_PBE Si 05Jan2001\n   TITEL  = PAW_PBE Si 05Jan2001\n   POMASS =   28.085; ZVAL   =    4.000\n"
                               "   ENMAX  =  245.345; ENMIN  =  184.009 eV\n", encoding="utf-8")
    cfg = cfg_for(sk_root)
    cfg.profiles["local"].env["VASP_PP_PATH"] = str(tmp_path / "pp")
    vs = water_spec(structure=Structure(source="bulk", source_ref="Si", atoms=AtomsData.from_ase(bulk("Si"))), method=VaspMethod(),
                    kpoints=KPoints(), task=Task())
    v = {x.quantity: x for x in documented_values(vs, cfg)}
    assert v["ENMAX"].value == 245.345 and v["ENMIN"].value == 184.009 and v["ENMAX"].line == 4 and v["ENMAX"].unit == "eV"

    make_fake_skset(tmp_path / "sk", "doc-1-0", ["H", "O"])
    (tmp_path / "sk" / "doc-1-0" / "README").write_text("set\n\nList of all atomic Hubbard derivatives:\n H = -0.1857\n O = -0.1575\n\nzeta = 4.00 (HX)\n", encoding="utf-8")
    cfg2 = cfg_for(tmp_path / "sk")
    v = documented_values(water_spec(method=DftbMethod(sk_set="doc-1-0")), cfg2)
    got = {(x.element, x.quantity): (x.value, x.line) for x in v}
    assert got[("H", "hubbard_derivative")] == (-0.1857, 4) and got[("O", "hubbard_derivative")] == (-0.1575, 5)
    assert got[(None, "zeta")] == (4.0, 7) and got[("H", "spin_constant")][0] == -0.05

    data = tmp_path / "cp2k"
    data.mkdir()
    (data / "BASIS_T").write_text("# first comment\n# DZVP basis for H\nH DZVP-TEST-GTH  # inline note\n 1\n 1 0 0 1 1\n 1.0 1.0\n", encoding="utf-8")
    (data / "POT_T").write_text("H GTH-TEST-q1\n 1\n 0.2 2 -4.0 0.7\n 0\n", encoding="utf-8")
    cfg3 = cfg_for(sk_root)
    cfg3.cp2k_data = str(data)
    h2 = molecule("H2")
    cs = water_spec(structure=Structure(source="preset", source_ref="H2", atoms=AtomsData.from_ase(h2)),
                    method=Cp2kMethod(basis_file="BASIS_T", potential_file="POT_T", xc="PBE", cutoff_ry=100, rel_cutoff_ry=30, poisson_solver="MT", isolated_box_ang=8))
    texts = [(x.value, x.line) for x in documented_values(cs, cfg3)]
    assert ("first comment", 1) in texts and ("DZVP basis for H", 2) in texts and ("inline note", 3) in texts
    assert documented_values(water_spec(method=XtbMethod()), cfg) == []


def test_cli_documented_and_compare_set(sk_root, tmp_path, capsys):
    from adit.cli import main
    from adit.config import save_config

    cfgp = tmp_path / "c.toml"
    save_config(cfg_for(sk_root), cfgp)
    slab = _slab_files(tmp_path)
    sp = tmp_path / "spec.json"
    _slab_spec(slab).save(sp)
    (tmp_path / "set.json").write_text(json.dumps(ADS), encoding="utf-8")
    assert main([str(sp), "--config", str(cfgp), "--documented"]) == 0
    assert main([str(sp), str(tmp_path / "o"), "--config", str(cfgp), "--compare-set", str(tmp_path / "set.json")]) == 0
    assert (tmp_path / "o" / "compare.json").is_file()
    assert main([str(sp), str(tmp_path / "o"), "--config", str(cfgp), "--compare-set", str(tmp_path / "set.json")]) == 1
    with pytest.raises(SystemExit):
        main([str(sp), str(tmp_path / "p"), "--config", str(cfgp), "--compare-set", str(tmp_path / "set.json"), "--scan", "method.scc_tolerance=1e-5,1e-6"])
    with pytest.raises(SystemExit):
        main([str(sp), str(tmp_path / "p"), "--config", str(cfgp), "--conformers", "3"])
