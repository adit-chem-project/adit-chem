
import re
from pathlib import Path

import pytest

from adit import lang
from adit.citations import DISPERSION, FUNCTIONALS, adit_citation, bibtex_text, citations_for
from adit.codes.base import GENERATORS
from adit.spec import CalculationSpec
from tests.conftest import water_spec

ROOT = Path(__file__).resolve().parents[1]
EX = ROOT / "examples"
VERIFIED = ("dftbplus", "vasp", "espresso", "xtb", "orca", "cp2k", "lammps", "gromacs", "openmm", "psi4", "abinit", "nwchem")


@pytest.fixture(autouse=True)
def _japanese():
    before = lang.LANGUAGE
    lang.set_language("ja")
    yield
    lang.set_language(before)


def _all_entries():
    import importlib

    from adit.codes.espresso import PSEUDO_SET_CITATIONS
    from adit.codes.sk_sets import SK_SET_CITATIONS
    from adit.codes.xtb import GFN_CITATIONS
    from adit.codes.vasp import PAW_CITATIONS

    out = list(PAW_CITATIONS)
    for code in VERIFIED:
        out += getattr(importlib.import_module(type(GENERATORS[code]).__module__), "CITATIONS")
    for table in (SK_SET_CITATIONS, PSEUDO_SET_CITATIONS, GFN_CITATIONS, FUNCTIONALS, DISPERSION):
        for entries in table.values():
            out += entries
    return out


@pytest.mark.parametrize("code", VERIFIED)
def test_every_verified_code_has_citations_with_doi_and_source(code):
    import importlib

    entries = getattr(importlib.import_module(type(GENERATORS[code]).__module__), "CITATIONS")
    assert entries
    for c in entries:
        assert c.doi and c.source.startswith("http") and c.key.startswith(code.replace("espresso", "qe").split("_")[0])
        assert c.bibtex.startswith("@") and f"{{{c.key}," in c.bibtex and c.doi in c.bibtex


def test_bibtex_entries_are_well_formed_and_keys_unique():
    entries = _all_entries()
    keys = [c.key for c in entries]
    assert len(entries) > 40
    for c in entries:
        assert c.bibtex.count("{") == c.bibtex.count("}"), c.key
        assert re.match(r"@\w+\{" + re.escape(c.key) + ",", c.bibtex)
        assert "year" in c.bibtex and "title" in c.bibtex
    seen = {}
    for c in entries:
        assert seen.setdefault(c.key, c.bibtex) == c.bibtex, c.key


def test_dftb_water_cites_code_set_and_adit():
    spec = CalculationSpec.load(EX / "dftb_md_water_generated" / "spec.json")
    cs = citations_for(spec, EX / "dftb_md_water_generated")
    assert cs.keys[:2] == ["dftbplus_hourahine2025", "dftbplus_hourahine2020"] and "mio_elstner1998" in cs.keys
    assert cs.not_recorded == []
    spec = water_spec()
    cs = citations_for(spec)
    assert any("fake-1-0" in n for n in cs.not_recorded)


def test_dispersion_and_functional_follow_the_record():
    spec = CalculationSpec.load(EX / "cp2k_h2o_generated" / "spec.json")
    spec.method.dispersion = "d3bj"
    cs = citations_for(spec)
    assert {"cp2k_kuhne2020", "pbe_perdew1996", "d3_grimme2010", "d3bj_grimme2011"} <= set(cs.keys)
    spec.method.xc = "BLYP"
    cs = citations_for(spec)
    assert "pbe_perdew1996" not in cs.keys and any("method.xc = BLYP" in n for n in cs.not_recorded)
    spec = CalculationSpec.load(EX / "vasp_h2o_generated" / "spec.json")
    spec.method.ivdw = 13
    cs = citations_for(spec)
    assert {"vasp_kresse1999", "pbe_perdew1996", "d4_caldeweyher2019"} <= set(cs.keys)
    spec.method.extra_incar = {"GGA": "PS"}
    spec.method.ivdw = 20
    cs = citations_for(spec)
    assert "pbesol_perdew2008" in cs.keys and "pbe_perdew1996" not in cs.keys and any("IVDW = 20" in n for n in cs.not_recorded)


def test_espresso_reads_the_functional_from_the_upf_headers_or_says_so():
    spec = CalculationSpec.load(EX / "qe_si_generated" / "spec.json")
    cs = citations_for(spec, EX / "qe_si_generated")
    assert "pslibrary_dalcorso2014" in cs.keys and any("UPF" in n for n in cs.not_recorded)
    spec.method.input_dft = "pbesol"
    spec.method.pseudo_set = "SSSP_efficiency"
    cs = citations_for(spec)
    assert "sssp_prandini2018" in cs.keys and "pbesol_perdew2008" in cs.keys and cs.not_recorded == []


def test_xtb_and_orca_variants():
    spec = CalculationSpec.load(EX / "xtb_vib_water_generated" / "spec.json")
    assert citations_for(spec).keys == ["xtb_bannwarth2020", "xtb_bannwarth2019"]
    spec.method.gfn = "ff"
    assert "xtb_spicher2020" in citations_for(spec).keys
    spec.method.gfn = "0"
    assert any("GFN0" in n for n in citations_for(spec).not_recorded)
    from adit.spec import OrcaMethod

    spec.method = OrcaMethod(method="B3LYP", basis="def2-SVP", extra_keywords="D3BJ TightSCF")
    cs = citations_for(spec)
    assert {"orca_neese2012", "orca_neese2025", "b3lyp_becke1993", "b3lyp_stephens1994", "d3_grimme2010", "d3bj_grimme2011"} <= set(cs.keys)
    spec.method = OrcaMethod(method="PBE-D4", basis="def2-SVP")
    assert {"pbe_perdew1996", "d4_caldeweyher2019"} <= set(citations_for(spec).keys)


def test_codes_without_a_verified_page_are_reported_as_not_recorded():
    from adit.spec import GaussianMethod

    spec = water_spec()
    spec.method = GaussianMethod(theory="HF", basis="sto-3g")
    cs = citations_for(spec)
    assert cs.keys == [] and any("gaussian" in n for n in cs.not_recorded)


def test_adit_citation_comes_from_citation_cff(tmp_path):
    c = adit_citation()
    assert c is not None and c.bibtex.startswith("@software{adit_") and "ADIT" in c.bibtex and "url" in c.bibtex
    assert adit_citation(tmp_path / "missing.cff") is None
    text = bibtex_text([c, c], header="h")
    assert text.count("@software") == 1 and text.startswith("% h\n")
