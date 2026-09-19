"""BibTeX references for the code, parameter families and corrections a run records in spec.json."""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from adit.lang import L


@dataclass(frozen=True)
class Citation:
    """One BibTeX entry with the DOI and the official page that asks for it."""

    key: str
    bibtex: str
    doi: str = ""
    source: str = ""


@dataclass
class CitationSet:
    """References found for one run plus the items that are not on record."""

    entries: list[Citation] = field(default_factory=list)
    not_recorded: list[str] = field(default_factory=list)

    @property
    def keys(self) -> list[str]:
        return [c.key for c in self.entries]

    def add(self, *citations: Citation) -> None:
        for c in citations:
            if c.key not in self.keys:
                self.entries.append(c)


_D3_WORDS = {"D3": "d3", "D3ZERO": "d3", "D3BJ": "d3bj", "D3(BJ)": "d3bj", "D4": "d4"}
# IVDW and GGA values as documented on the VASP wiki (IVDW, GGA pages); without GGA the POTCAR's functional is used
_VASP_IVDW = {11: "d3", 12: "d3bj", 13: "d4"}
_VASP_GGA = {"PE": "PBE", "PS": "PBESOL", "B3": "B3LYP", "B5": "B3LYP"}


def _functional_by_name(name: str, cs: CitationSet, where: str) -> None:
    token = name.strip().upper()
    if not token:
        return
    base = token
    for word, disp in _D3_WORDS.items():
        if base.endswith("-" + word):
            base = base[: -len(word) - 1]
            _dispersion(disp, cs, where)
    if base in FUNCTIONALS:
        cs.add(*FUNCTIONALS[base])
    else:
        cs.not_recorded.append(L(f"{where} = {name.strip()} の文献", f"reference for {where} = {name.strip()}"))


def _dispersion(kind: str, cs: CitationSet, where: str) -> None:
    if kind in DISPERSION:
        if kind == "d3bj":  # the BJ damping paper builds on the D3 paper
            cs.add(*DISPERSION["d3"])
        cs.add(*DISPERSION[kind])
    else:
        cs.not_recorded.append(L(f"{where} = {kind} の文献", f"reference for {where} = {kind}"))


def _code_citations(code: str) -> tuple[Citation, ...]:
    from adit.codes.base import GENERATORS
    import importlib

    gen = GENERATORS.get(code)
    if gen is None:
        return ()
    module = importlib.import_module(type(gen).__module__)
    found = getattr(module, "CITATIONS", ())
    if isinstance(found, dict):
        return tuple(found.get(code, ()))
    return tuple(found)


def _upf_functionals(run_dir: Path | None) -> set[str]:
    if run_dir is None or not (run_dir / "pseudo").is_dir():
        return set()
    from adit.codes.upf import UpfError, read_upf_header

    out: set[str] = set()
    for p in sorted((run_dir / "pseudo").iterdir()):
        if p.suffix.upper() != ".UPF":
            continue
        try:
            out.add(read_upf_header(p).functional.strip().upper())
        except UpfError:
            continue
    return out


def _dftbplus(m, cs: CitationSet) -> None:
    from adit.codes.sk_sets import SK_SET_CITATIONS

    family = m.sk_set.split("-")[0].lower()
    if family in SK_SET_CITATIONS:
        cs.add(*SK_SET_CITATIONS[family])
    else:
        cs.not_recorded.append(L(f"Slater-Koster セット {m.sk_set} の文献", f"reference for the Slater-Koster set {m.sk_set}"))
    if m.dispersion == "dftd3":
        # dftb_in.hsd is written with Damping = BeckeJohnson (see codes/dftbplus.py)
        _dispersion("d3bj", cs, "method.dispersion")
    elif m.dispersion != "none":
        cs.not_recorded.append(L(f"method.dispersion = {m.dispersion} の文献", f"reference for method.dispersion = {m.dispersion}"))


def _espresso(m, cs: CitationSet, run_dir: Path | None) -> None:
    from adit.codes.espresso import PSEUDO_SET_CITATIONS

    name = m.pseudo_set.strip().lower()
    family = next((k for k in PSEUDO_SET_CITATIONS if k in name), None)
    if family:
        cs.add(*PSEUDO_SET_CITATIONS[family])
    else:
        cs.not_recorded.append(L(f"擬ポテンシャルのセット {m.pseudo_set or '(未指定)'} の文献",
                                 f"reference for the pseudopotential set {m.pseudo_set or '(not set)'}"))
    if m.input_dft.strip():
        _functional_by_name(m.input_dft, cs, "method.input_dft")
        return
    found = _upf_functionals(run_dir)
    if len(found) == 1 and next(iter(found)) in FUNCTIONALS:
        cs.add(*FUNCTIONALS[next(iter(found))])
    else:
        cs.not_recorded.append(L("汎関数の文献 (input_dft は未指定。汎関数は UPF のヘッダにだけ書かれています)",
                                 "reference for the functional (input_dft is not set; the functional is stated only in the UPF headers)"))


def _vasp(m, cs: CitationSet) -> None:
    from adit.codes.vasp import PAW_CITATIONS

    cs.add(*PAW_CITATIONS)
    extra = {k.strip().upper(): str(v).strip().upper() for k, v in m.extra_incar.items()}
    if "GGA" in extra:
        name = _VASP_GGA.get(extra["GGA"], "")
        if name:
            _functional_by_name(name, cs, "GGA")
        else:
            cs.not_recorded.append(L(f"GGA = {extra['GGA']} の文献", f"reference for GGA = {extra['GGA']}"))
    elif "PBE" in m.potcar_set.upper():
        # without a GGA tag VASP uses the functional the POTCAR set was made with
        _functional_by_name("PBE", cs, "method.potcar_set")
    else:
        cs.not_recorded.append(L(f"汎関数の文献 (POTCAR のセット {m.potcar_set} から判断できません)",
                                 f"reference for the functional (cannot be read from the POTCAR set {m.potcar_set})"))
    if m.ivdw is not None:
        kind = _VASP_IVDW.get(m.ivdw)
        if kind:
            _dispersion(kind, cs, f"IVDW = {m.ivdw}")
        else:
            cs.not_recorded.append(L(f"IVDW = {m.ivdw} の文献", f"reference for IVDW = {m.ivdw}"))


def _keywords(text: str, cs: CitationSet, where: str) -> None:
    for word in re.split(r"[\s,]+", text.upper()):
        if word in _D3_WORDS:
            _dispersion(_D3_WORDS[word], cs, where)


def citations_for(spec, run_dir: Path | str | None = None) -> CitationSet:
    """Collect the references a spec records: code, parameter family, functional and dispersion."""
    cs = CitationSet()
    d = Path(run_dir).expanduser() if run_dir is not None else None
    m = spec.method
    code = m.code
    entries = _code_citations(code)
    if entries:
        cs.add(*entries)
    else:
        cs.not_recorded.append(L(f"計算コード {code} の文献", f"reference for the code {code}"))
    if code == "xtb":
        from adit.codes.xtb import GFN_CITATIONS

        if m.gfn in GFN_CITATIONS:
            cs.add(*GFN_CITATIONS[m.gfn])
        else:
            cs.not_recorded.append(L(f"GFN{m.gfn} の文献", f"reference for GFN{m.gfn}"))
    elif code == "dftbplus":
        _dftbplus(m, cs)
    elif code == "espresso":
        _espresso(m, cs, d)
    elif code == "vasp":
        _vasp(m, cs)
    elif code == "cp2k":
        _functional_by_name(m.xc, cs, "method.xc")
        if m.dispersion in ("d3", "d3bj"):
            _dispersion(m.dispersion, cs, "method.dispersion")
    elif code == "orca":
        _functional_by_name(m.method, cs, "method.method")
        _keywords(m.extra_keywords, cs, "method.extra_keywords")
    elif code == "psi4":
        _functional_by_name(m.method, cs, "method.method")
    elif code == "nwchem":
        if m.xc.strip():
            _functional_by_name(m.xc, cs, "method.xc")
    elif code == "openmx":
        if m.xc == "GGA-PBE":
            _functional_by_name("PBE", cs, "method.xc")
        elif m.xc:
            cs.not_recorded.append(L(f"method.xc = {m.xc} の文献", f"reference for method.xc = {m.xc}"))
    return cs


_CFF_SEARCH = (Path(__file__).resolve().parents[2], Path(__file__).resolve().parent)


def _cff_scalar(text: str, key: str) -> str:
    m = re.search(rf"^{key}:\s*(.+?)\s*$", text, re.MULTILINE)
    return m.group(1).strip().strip('"') if m else ""


def adit_citation(cff_path: Path | str | None = None) -> Citation | None:
    """Build the @software entry for ADIT from CITATION.cff; None when the file is not found."""
    if cff_path is not None:
        candidates = [Path(cff_path)]
    else:
        candidates = [d / "CITATION.cff" for d in _CFF_SEARCH]
    path = next((p for p in candidates if p.is_file()), None)
    if path is None:
        return None
    text = path.read_text(encoding="utf-8", errors="replace")
    title = _cff_scalar(text, "title")
    version = _cff_scalar(text, "version")
    url = _cff_scalar(text, "repository-code")
    date = _cff_scalar(text, "date-released")
    authors = re.findall(r"^\s*-\s*name:\s*(.+?)\s*$", text, re.MULTILINE)
    author = " and ".join("{" + a.strip().strip('"') + "}" for a in authors) or "{The ADIT project}"
    year = date[:4] if date else ""
    key = f"adit_{version.replace('.', '_')}" if version else "adit"
    lines = [f"@software{{{key},", f"  author  = {{{author}}},", f"  title   = {{{title or 'ADIT'}}},"]
    if version:
        lines.append(f"  version = {{{version}}},")
    if year:
        lines.append(f"  year    = {{{year}}},")
    if date:
        lines.append(f"  date    = {{{date}}},")
    if url:
        lines.append(f"  url     = {{{url}}},")
    lines.append("}")
    return Citation(key=key, bibtex="\n".join(lines), source=str(path))


def bibtex_text(citations: list[Citation], *, header: str = "") -> str:
    """Join unique entries into one .bib file; each entry is preceded by its source page as a comment."""
    out = [f"% {header}"] if header else []
    seen: set[str] = set()
    for c in citations:
        if c.key in seen:
            continue
        seen.add(c.key)
        if c.source:
            out.append(f"% {L('出典', 'source')}: {c.source}")
        out.append(c.bibtex.strip())
        out.append("")
    return "\n".join(out) + ("\n" if out and out[-1] != "" else "")


# Functionals and dispersion corrections are looked up by the name recorded in spec.json; the
# table is a lookup, not a recommendation. Bibliographic data were checked against Crossref.
_D3_README = "https://github.com/dftd3/simple-dftd3#readme"
_D4_README = "https://github.com/dftd4/dftd4#readme"
FUNCTIONALS: dict[str, tuple[Citation, ...]] = {
    "PBE": (Citation("pbe_perdew1996", r"""@article{pbe_perdew1996,
  author  = {Perdew, John P. and Burke, Kieron and Ernzerhof, Matthias},
  title   = {Generalized Gradient Approximation Made Simple},
  journal = {Physical Review Letters},
  volume  = {77},
  number  = {18},
  pages   = {3865--3868},
  year    = {1996},
  doi     = {10.1103/PhysRevLett.77.3865}
}""", doi="10.1103/PhysRevLett.77.3865", source="https://doi.org/10.1103/PhysRevLett.77.3865"),
            Citation("pbe_perdew1997erratum", r"""@article{pbe_perdew1997erratum,
  author  = {Perdew, John P. and Burke, Kieron and Ernzerhof, Matthias},
  title   = {Generalized Gradient Approximation Made Simple [Phys. Rev. Lett. 77, 3865 (1996)]},
  journal = {Physical Review Letters},
  volume  = {78},
  number  = {7},
  pages   = {1396},
  year    = {1997},
  doi     = {10.1103/PhysRevLett.78.1396}
}""", doi="10.1103/PhysRevLett.78.1396", source="https://doi.org/10.1103/PhysRevLett.78.1396")),
    "PBESOL": (Citation("pbesol_perdew2008", r"""@article{pbesol_perdew2008,
  author  = {Perdew, John P. and Ruzsinszky, Adrienn and Csonka, G{\'a}bor I. and Vydrov, Oleg A. and Scuseria, Gustavo E. and Constantin, Lucian A. and Zhou, Xiaolan and Burke, Kieron},
  title   = {Restoring the Density-Gradient Expansion for Exchange in Solids and Surfaces},
  journal = {Physical Review Letters},
  volume  = {100},
  number  = {13},
  pages   = {136406},
  year    = {2008},
  doi     = {10.1103/PhysRevLett.100.136406}
}""", doi="10.1103/PhysRevLett.100.136406", source="https://doi.org/10.1103/PhysRevLett.100.136406"),),
    "B3LYP": (Citation("b3lyp_becke1993", r"""@article{b3lyp_becke1993,
  author  = {Becke, Axel D.},
  title   = {Density-functional thermochemistry. {III}. The role of exact exchange},
  journal = {The Journal of Chemical Physics},
  volume  = {98},
  number  = {7},
  pages   = {5648--5652},
  year    = {1993},
  doi     = {10.1063/1.464913}
}""", doi="10.1063/1.464913", source="https://doi.org/10.1063/1.464913"),
              Citation("b3lyp_stephens1994", r"""@article{b3lyp_stephens1994,
  author  = {Stephens, P. J. and Devlin, F. J. and Chabalowski, C. F. and Frisch, M. J.},
  title   = {Ab Initio Calculation of Vibrational Absorption and Circular Dichroism Spectra Using Density Functional Force Fields},
  journal = {The Journal of Physical Chemistry},
  volume  = {98},
  number  = {45},
  pages   = {11623--11627},
  year    = {1994},
  doi     = {10.1021/j100096a001}
}""", doi="10.1021/j100096a001", source="https://doi.org/10.1021/j100096a001")),
}
DISPERSION: dict[str, tuple[Citation, ...]] = {
    "d3": (Citation("d3_grimme2010", r"""@article{d3_grimme2010,
  author  = {Grimme, Stefan and Antony, Jens and Ehrlich, Stephan and Krieg, Helge},
  title   = {A consistent and accurate ab initio parametrization of density functional dispersion correction ({DFT-D}) for the 94 elements {H-Pu}},
  journal = {The Journal of Chemical Physics},
  volume  = {132},
  number  = {15},
  pages   = {154104},
  year    = {2010},
  doi     = {10.1063/1.3382344}
}""", doi="10.1063/1.3382344", source=_D3_README),),
    "d3bj": (Citation("d3bj_grimme2011", r"""@article{d3bj_grimme2011,
  author  = {Grimme, Stefan and Ehrlich, Stephan and Goerigk, Lars},
  title   = {Effect of the damping function in dispersion corrected density functional theory},
  journal = {Journal of Computational Chemistry},
  volume  = {32},
  number  = {7},
  pages   = {1456--1465},
  year    = {2011},
  doi     = {10.1002/jcc.21759}
}""", doi="10.1002/jcc.21759", source=_D3_README),),
    "d4": (Citation("d4_caldeweyher2019", r"""@article{d4_caldeweyher2019,
  author  = {Caldeweyher, Eike and Ehlert, Sebastian and Hansen, Andreas and Neugebauer, Hagen and Spicher, Sebastian and Bannwarth, Christoph and Grimme, Stefan},
  title   = {A generally applicable atomic-charge dependent {London} dispersion correction},
  journal = {The Journal of Chemical Physics},
  volume  = {150},
  number  = {15},
  pages   = {154122},
  year    = {2019},
  doi     = {10.1063/1.5090222}
}""", doi="10.1063/1.5090222", source=_D4_README),
           Citation("d4_caldeweyher2017", r"""@article{d4_caldeweyher2017,
  author  = {Caldeweyher, Eike and Bannwarth, Christoph and Grimme, Stefan},
  title   = {Extension of the {D3} dispersion coefficient model},
  journal = {The Journal of Chemical Physics},
  volume  = {147},
  number  = {3},
  pages   = {034112},
  year    = {2017},
  doi     = {10.1063/1.4993215}
}""", doi="10.1063/1.4993215", source=_D4_README),
           Citation("d4_caldeweyher2020", r"""@article{d4_caldeweyher2020,
  author  = {Caldeweyher, Eike and Mewes, Jan-Michael and Ehlert, Sebastian and Grimme, Stefan},
  title   = {Extension and evaluation of the {D4} {London}-dispersion model for periodic systems},
  journal = {Physical Chemistry Chemical Physics},
  volume  = {22},
  number  = {16},
  pages   = {8499--8512},
  year    = {2020},
  doi     = {10.1039/D0CP00502A}
}""", doi="10.1039/D0CP00502A", source=_D4_README)),
}
