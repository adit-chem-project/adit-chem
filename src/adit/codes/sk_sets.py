
from __future__ import annotations

from adit.citations import Citation
from adit.lang import L

import itertools
import re
from dataclasses import dataclass, field
from pathlib import Path

_SKF_NAME = re.compile(r"^([A-Z][a-z]?)-([A-Z][a-z]?)\.skf$")
_SHELLS = re.compile(r"<Shells>\s*([^<]*?)\s*</Shells>", re.IGNORECASE)
_ELEM_VALUE = re.compile(r"^\s*([A-Z][a-z]?)\s*=\s*(-?\d+(?:\.\d+)?)\s*$")
_ZETA = re.compile(r"zeta\s*=\s*(\d+(?:\.\d+)?)", re.IGNORECASE)
_ELEM_HEADER = re.compile(r"^\s*([A-Z][a-z]?)\s*:\s*$")
_L_OF = {"s": 0, "p": 1, "d": 2, "f": 3}
_NAME_OF_L = {v: k for k, v in _L_OF.items()}

REQUIRED_DOCS = ("LICENSE", "README")


class SKSetError(Exception):
    pass


@dataclass
class SKSet:
    name: str
    root: Path
    pairs: dict[tuple[str, str], Path] = field(default_factory=dict)

    @classmethod
    def from_dir(cls, root: Path | str, name: str | None = None) -> "SKSet":
        root = Path(root)
        if not root.is_dir():
            raise SKSetError(L(f"SK セットのディレクトリがありません: {root}", f"Slater-Koster set directory not found: {root}"))
        pairs: dict[tuple[str, str], Path] = {}
        for p in sorted(root.iterdir()):
            m = _SKF_NAME.match(p.name)
            if m:
                pairs[(m.group(1), m.group(2))] = p
        if not pairs:
            raise SKSetError(L(f"*.skf が 1 つもありません: {root}", f"no *.skf files in {root}"))
        return cls(name=name or root.name, root=root, pairs=pairs)

    @property
    def elements(self) -> list[str]:
        return sorted({a for (a, b) in self.pairs if a == b})

    def missing_pairs(self, elements: list[str]) -> list[tuple[str, str]]:
        return [(a, b) for a, b in itertools.product(elements, repeat=2) if (a, b) not in self.pairs]

    def required_files(self, elements: list[str]) -> dict[tuple[str, str], Path]:
        missing = self.missing_pairs(elements)
        if missing:
            raise SKSetError(L(f"セット {self.name} に無いペアがあります: {missing}", f"pairs missing from set {self.name}: {missing}"))
        return {(a, b): self.pairs[(a, b)] for a, b in itertools.product(elements, repeat=2)}

    def doc_files(self) -> dict[str, Path]:
        return {n: self.root / n for n in REQUIRED_DOCS if (self.root / n).is_file()}

    def max_angular_momentum(self, element: str) -> str:
        path = self.pairs.get((element, element))
        if path is None:
            raise SKSetError(L(f"{element}-{element}.skf がありません", f"{element}-{element}.skf not found"))
        text = path.read_text(encoding="utf-8", errors="replace")
        m = _SHELLS.search(text)
        if not m:
            raise SKSetError(L(f"{path.name} の文書に <Shells> が無く、角運動量を決められません", f"{path.name} has no <Shells> in its documentation; cannot decide the angular momentum"))
        shells = m.group(1).split()
        ls = [_L_OF[s[-1].lower()] for s in shells if s[-1].lower() in _L_OF]
        if not ls:
            raise SKSetError(L(f"{path.name} の <Shells> を解釈できません: {m.group(1)!r}", f"cannot parse <Shells> of {path.name}: {m.group(1)!r}"))
        return _NAME_OF_L[max(ls)]


    def hubbard_derivs(self) -> dict[str, float] | None:
        readme = self.root / "README"
        if not readme.is_file():
            return None
        lines = readme.read_text(encoding="utf-8", errors="replace").splitlines()
        start = next((i for i, l in enumerate(lines) if "hubbard derivative" in l.lower()), None)
        if start is None:
            return None
        out: dict[str, float] = {}
        for l in lines[start + 1:]:
            m = _ELEM_VALUE.match(l)
            if m:
                out[m.group(1)] = float(m.group(2))
            elif out and l.strip() == "":
                break
        return out or None

    def damping_exponent(self) -> float | None:
        readme = self.root / "README"
        if not readme.is_file():
            return None
        m = _ZETA.search(readme.read_text(encoding="utf-8", errors="replace"))
        return float(m.group(1)) if m else None

    def spin_constants(self) -> dict[str, float] | None:
        f = self.root / "spinw.txt"
        if not f.is_file():
            return None
        out: dict[str, float] = {}
        current: str | None = None
        rows: list[list[float]] = []

        def flush() -> None:
            if current and rows:
                n = len(rows)
                if any(len(r) != n for r in rows):
                    raise SKSetError(L(f"spinw.txt の {current} が正方行列ではありません", f"{current} in spinw.txt is not a square matrix"))
                out[current] = rows[-1][-1]

        for l in f.read_text(encoding="utf-8", errors="replace").splitlines():
            m = _ELEM_HEADER.match(l)
            if m:
                flush(); current, rows = m.group(1), []
            elif l.strip() and current:
                rows.append([float(x) for x in l.split()])
        flush()
        return out or None


def discover_sets(sk_root: Path | str) -> dict[str, SKSet]:
    sk_root = Path(sk_root)
    found: dict[str, SKSet] = {}
    if not sk_root.is_dir():
        return found
    for d in sorted(sk_root.iterdir()):
        if d.is_dir():
            try:
                found[d.name] = SKSet.from_dir(d)
            except SKSetError:
                continue
    return found


# "Required references" from the README of each set on dftb.org (github.com/dftbparams/<set>), keyed by the family name
# before the first "-" (mio-1-1 -> mio). Sets whose README names only theses or unpublished tables carry those entries too.
_SK_README = "https://dftb.org/parameters/download.html"
_MIO = (
    Citation("mio_elstner1998", r"""@article{mio_elstner1998,
  author  = {Elstner, M. and Porezag, D. and Jungnickel, G. and Elsner, J. and Haugk, M. and Frauenheim, Th. and Suhai, S. and Seifert, G.},
  title   = {Self-consistent-charge density-functional tight-binding method for simulations of complex materials properties},
  journal = {Physical Review B},
  volume  = {58},
  number  = {11},
  pages   = {7260--7268},
  year    = {1998},
  doi     = {10.1103/PhysRevB.58.7260}
}""", doi="10.1103/PhysRevB.58.7260", source=_SK_README),
    Citation("mio_niehaus2001", r"""@article{mio_niehaus2001,
  author  = {Niehaus, T. A. and Elstner, M. and Frauenheim, Th. and Suhai, S.},
  title   = {Application of an approximate density-functional method to sulfur containing compounds},
  journal = {Journal of Molecular Structure: THEOCHEM},
  volume  = {541},
  number  = {1-3},
  pages   = {185--194},
  year    = {2001},
  doi     = {10.1016/S0166-1280(00)00762-4}
}""", doi="10.1016/S0166-1280(00)00762-4", source=_SK_README),
    Citation("mio_gaus2011", r"""@article{mio_gaus2011,
  author  = {Gaus, Michael and Cui, Qiang and Elstner, Marcus},
  title   = {{DFTB3}: Extension of the Self-Consistent-Charge Density-Functional Tight-Binding Method ({SCC-DFTB})},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {7},
  number  = {4},
  pages   = {931--948},
  year    = {2011},
  doi     = {10.1021/ct100684s}
}""", doi="10.1021/ct100684s", source=_SK_README),
)
_THREEOB = (
    Citation("threeob_gaus2013", r"""@article{threeob_gaus2013,
  author  = {Gaus, Michael and Goez, Albrecht and Elstner, Marcus},
  title   = {Parametrization and Benchmark of {DFTB3} for Organic Molecules},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {9},
  number  = {1},
  pages   = {338--354},
  year    = {2013},
  doi     = {10.1021/ct300849w}
}""", doi="10.1021/ct300849w", source=_SK_README),
    Citation("threeob_gaus2014", r"""@article{threeob_gaus2014,
  author  = {Gaus, Michael and Lu, Xiya and Elstner, Marcus and Cui, Qiang},
  title   = {Parameterization of {DFTB3/3OB} for Sulfur and Phosphorus for Chemical and Biological Applications},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {10},
  number  = {4},
  pages   = {1518--1537},
  year    = {2014},
  doi     = {10.1021/ct401002w}
}""", doi="10.1021/ct401002w", source=_SK_README),
    Citation("threeob_lu2015", r"""@article{threeob_lu2015,
  author  = {Lu, Xiya and Gaus, Michael and Elstner, Marcus and Cui, Qiang},
  title   = {Parametrization of {DFTB3/3OB} for Magnesium and Zinc for Chemical and Biological Applications},
  journal = {The Journal of Physical Chemistry B},
  volume  = {119},
  number  = {3},
  pages   = {1062--1082},
  year    = {2015},
  doi     = {10.1021/jp506557r}
}""", doi="10.1021/jp506557r", source=_SK_README),
    Citation("threeob_kubillus2015", r"""@article{threeob_kubillus2015,
  author  = {Kubillus, Maximilian and Kuba{\v r}, Tom{\'a}{\v s} and Gaus, Michael and {\v R}ez{\'a}{\v c}, Jan and Elstner, Marcus},
  title   = {Parameterization of the {DFTB3} Method for {Br}, {Ca}, {Cl}, {F}, {I}, {K}, and {Na} in Organic and Biological Systems},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {11},
  number  = {1},
  pages   = {332--342},
  year    = {2015},
  doi     = {10.1021/ct5009137}
}""", doi="10.1021/ct5009137", source=_SK_README),
)
_PBC = (
    Citation("pbc_sieck2000", r"""@phdthesis{pbc_sieck2000,
  author = {Sieck, A.},
  title  = {Structure and physical properties of silicon clusters and of vacancy clusters in bulk silicon},
  school = {University of Paderborn},
  year   = {2000}
}""", source=_SK_README),
    Citation("pbc_koehler2001", r"""@article{pbc_koehler2001,
  author  = {K{\"o}hler, Christof and Hajnal, Zolt{\'a}n and De{\'a}k, P{\'e}ter and Frauenheim, Thomas and Suhai, S{\'a}ndor},
  title   = {Theoretical investigation of carbon defects and diffusion in $\alpha$-quartz},
  journal = {Physical Review B},
  volume  = {64},
  number  = {8},
  pages   = {085333},
  year    = {2001},
  doi     = {10.1103/PhysRevB.64.085333}
}""", doi="10.1103/PhysRevB.64.085333", source=_SK_README),
    Citation("pbc_koehler2006", r"""@article{pbc_koehler2006,
  author  = {K{\"o}hler, Christof and Frauenheim, Thomas},
  title   = {Molecular dynamics simulations of {CFx} (x=2,3) molecules at {Si3N4} and {SiO2} surfaces},
  journal = {Surface Science},
  volume  = {600},
  number  = {2},
  pages   = {453--460},
  year    = {2006},
  doi     = {10.1016/j.susc.2005.10.044}
}""", doi="10.1016/j.susc.2005.10.044", source=_SK_README),
    Citation("pbc_koehler2005", r"""@article{pbc_koehler2005,
  author  = {K{\"o}hler, Christof and Seifert, Gotthard and Frauenheim, Thomas},
  title   = {Density functional based calculations for {Fe$_n$} ($n \le 32$)},
  journal = {Chemical Physics},
  volume  = {309},
  number  = {1},
  pages   = {23--31},
  year    = {2005},
  doi     = {10.1016/j.chemphys.2004.03.034}
}""", doi="10.1016/j.chemphys.2004.03.034", source=_SK_README),
)
_MATSCI_TABLES = Citation("matsci_frenzel2009", r"""@misc{matsci_frenzel2009,
  author       = {Frenzel, J. and Oliveira, A. F. and Jardillier, N. and Heine, T. and Seifert, G.},
  title        = {Semi-relativistic, self-consistent charge Slater-Koster tables for density-functional based tight-binding ({DFTB}) for materials science simulations},
  howpublished = {TU Dresden},
  year         = {2004--2009}
}""", source=_SK_README)
_GUIMARAES = Citation("matsci_guimaraes2007", r"""@article{matsci_guimaraes2007,
  author  = {Guimar{\~a}es, Luciana and Enyashin, Andrey N. and Frenzel, Johannes and Heine, Thomas and Duarte, H{\'e}lio A. and Seifert, Gotthard},
  title   = {Imogolite Nanotubes: Stability, Electronic, and Mechanical Properties},
  journal = {ACS Nano},
  volume  = {1},
  number  = {4},
  pages   = {362--368},
  year    = {2007},
  doi     = {10.1021/nn700184k}
}""", doi="10.1021/nn700184k", source=_SK_README)
_MATSCI = (
    _MATSCI_TABLES,
    Citation("matsci_frenzel2005", r"""@article{matsci_frenzel2005,
  author  = {Frenzel, Johannes and Oliveira, Augusto F. and Duarte, Helio A. and Heine, Thomas and Seifert, Gotthard},
  title   = {Structural and Electronic Properties of Bulk Gibbsite and Gibbsite Surfaces},
  journal = {Zeitschrift f{\"u}r anorganische und allgemeine Chemie},
  volume  = {631},
  number  = {6-7},
  pages   = {1267--1271},
  year    = {2005},
  doi     = {10.1002/zaac.200500051}
}""", doi="10.1002/zaac.200500051", source=_SK_README),
    _GUIMARAES,
    Citation("matsci_luschtinetz2008", r"""@article{matsci_luschtinetz2008,
  author  = {Luschtinetz, Regina and Oliveira, Augusto F. and Frenzel, Johannes and Joswig, Jan-Ole and Seifert, Gotthard and Duarte, Helio A.},
  title   = {Adsorption of phosphonic and ethylphosphonic acid on aluminum oxide surfaces},
  journal = {Surface Science},
  volume  = {602},
  number  = {7},
  pages   = {1347--1359},
  year    = {2008},
  doi     = {10.1016/j.susc.2008.01.035}
}""", doi="10.1016/j.susc.2008.01.035", source=_SK_README),
    Citation("matsci_luschtinetz2009", r"""@article{matsci_luschtinetz2009,
  author  = {Luschtinetz, Regina and Frenzel, Johannes and Milek, Theodor and Seifert, Gotthard},
  title   = {Adsorption of Phosphonic Acid at the {TiO$_2$} Anatase (101) and Rutile (110) Surfaces},
  journal = {The Journal of Physical Chemistry C},
  volume  = {113},
  number  = {14},
  pages   = {5730--5740},
  year    = {2009},
  doi     = {10.1021/jp8110343}
}""", doi="10.1021/jp8110343", source=_SK_README),
)
_OB2 = (
    Citation("ob2_vuong2018", r"""@article{ob2_vuong2018,
  author  = {Vuong, Van Quan and Akkarapattiakal Kuriappan, Jissy and Kubillus, Maximilian and Kranz, Julian J. and Mast, Thilo and Niehaus, Thomas A. and Irle, Stephan and Elstner, Marcus},
  title   = {Parametrization and Benchmark of Long-Range Corrected {DFTB2} for Organic Molecules},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {14},
  number  = {1},
  pages   = {115--125},
  year    = {2018},
  doi     = {10.1021/acs.jctc.7b00947}
}""", doi="10.1021/acs.jctc.7b00947", source=_SK_README),
)
_TIORG = (
    Citation("tiorg_dolgonos2010", r"""@article{tiorg_dolgonos2010,
  author  = {Dolgonos, Grygoriy and Aradi, B{\'a}lint and Moreira, Ney H. and Frauenheim, Thomas},
  title   = {An Improved Self-Consistent-Charge Density-Functional Tight-Binding ({SCC-DFTB}) Set of Parameters for Simulation of Bulk and Molecular Systems Involving Titanium},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {6},
  number  = {1},
  pages   = {266--278},
  year    = {2010},
  doi     = {10.1021/ct900422c}
}""", doi="10.1021/ct900422c", source=_SK_README),
)
_TRANS3D = (
    Citation("trans3d_zheng2007", r"""@article{trans3d_zheng2007,
  author  = {Zheng, Guishan and Witek, Henryk A. and Bobadova-Parvanova, Petia and Irle, Stephan and Musaev, Djamaladdin G. and Prabhakar, Rajeev and Morokuma, Keiji and Lundberg, Marcus and Elstner, Marcus and K{\"o}hler, Christof and Frauenheim, Thomas},
  title   = {Parameter Calibration of Transition-Metal Elements for the Spin-Polarized Self-Consistent-Charge Density-Functional Tight-Binding ({DFTB}) Method: {Sc}, {Ti}, {Fe}, {Co}, and {Ni}},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {3},
  number  = {4},
  pages   = {1349--1367},
  year    = {2007},
  doi     = {10.1021/ct600312f}
}""", doi="10.1021/ct600312f", source=_SK_README),
)
_HALORG = _MIO + (
    Citation("halorg_kubar2013", r"""@article{halorg_kubar2013,
  author  = {Kuba{\v r}, Tom{\'a}{\v s} and Bodrog, Zolt{\'a}n and Gaus, Michael and K{\"o}hler, Christof and Aradi, B{\'a}lint and Frauenheim, Thomas and Elstner, Marcus},
  title   = {Parametrization of the {SCC-DFTB} Method for Halogens},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {9},
  number  = {7},
  pages   = {2939--2949},
  year    = {2013},
  doi     = {10.1021/ct4001922}
}""", doi="10.1021/ct4001922", source=_SK_README),
)
_MAGSIL = (
    _MATSCI_TABLES,
    Citation("magsil_lourenco2012", r"""@article{magsil_lourenco2012,
  author  = {Louren{\c c}o, Maicon P. and de Oliveira, Claudio and Oliveira, Augusto F. and Guimar{\~a}es, Luciana and Duarte, H{\'e}lio A.},
  title   = {Structural, Electronic, and Mechanical Properties of Single-Walled Chrysotile Nanotube Models},
  journal = {The Journal of Physical Chemistry C},
  volume  = {116},
  number  = {17},
  pages   = {9405--9411},
  year    = {2012},
  doi     = {10.1021/jp301048p}
}""", doi="10.1021/jp301048p", source=_SK_README),
    _GUIMARAES,
)
_ZNORG = (
    Citation("znorg_moreira2009", r"""@article{znorg_moreira2009,
  author  = {Moreira, Ney H. and Dolgonos, Grygoriy and Aradi, B{\'a}lint and da Rosa, Andreia L. and Frauenheim, Thomas},
  title   = {Toward an Accurate Density-Functional Tight-Binding Description of Zinc-Containing Compounds},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {5},
  number  = {3},
  pages   = {605--614},
  year    = {2009},
  doi     = {10.1021/ct800455a}
}""", doi="10.1021/ct800455a", source=_SK_README),
)
_AUORG = _MIO[:2] + (
    Citation("auorg_fihey2015", r"""@article{auorg_fihey2015,
  author  = {Fihey, Arnaud and Hettich, Christian and Touzeau, J{\'e}r{\'e}my and Maurel, Fran{\c c}ois and Perrier, Aur{\'e}lie and K{\"o}hler, Christof and Aradi, B{\'a}lint and Frauenheim, Thomas},
  title   = {{SCC-DFTB} parameters for simulating hybrid gold-thiolates compounds},
  journal = {Journal of Computational Chemistry},
  volume  = {36},
  number  = {27},
  pages   = {2075--2087},
  year    = {2015},
  doi     = {10.1002/jcc.24046}
}""", doi="10.1002/jcc.24046", source=_SK_README),
)
_SIBAND = (
    Citation("siband_markov2015a", r"""@article{siband_markov2015a,
  author  = {Markov, Stanislav and Aradi, Balint and Yam, Chi-Yung and Xie, Hang and Frauenheim, Thomas and Chen, Guanhua},
  title   = {Atomic Level Modeling of Extremely Thin Silicon-on-Insulator {MOSFETs} Including the Silicon Dioxide: Electronic Structure},
  journal = {IEEE Transactions on Electron Devices},
  volume  = {62},
  number  = {3},
  pages   = {696--704},
  year    = {2015},
  doi     = {10.1109/TED.2014.2387288}
}""", doi="10.1109/TED.2014.2387288", source=_SK_README),
    Citation("siband_markov2015b", r"""@article{siband_markov2015b,
  author  = {Markov, Stanislav and Penazzi, Gabriele and Kwok, YanHo and Pecchia, Alessandro and Aradi, Balint and Frauenheim, Thomas and Chen, GuanHua},
  title   = {Permittivity of Oxidized Ultra-Thin Silicon Films From Atomistic Simulations},
  journal = {IEEE Electron Device Letters},
  volume  = {36},
  number  = {10},
  pages   = {1076--1078},
  year    = {2015},
  doi     = {10.1109/LED.2015.2465850}
}""", doi="10.1109/LED.2015.2465850", source=_SK_README),
)
_BORG = (
    Citation("borg_grundkoetterstock2012", r"""@article{borg_grundkoetterstock2012,
  author  = {Grundk{\"o}tter-Stock, Bernhard and Bezugly, Viktor and Kunstmann, Jens and Cuniberti, Gianaurelio and Frauenheim, Thomas and Niehaus, Thomas A.},
  title   = {{SCC-DFTB} Parametrization for Boron and Boranes},
  journal = {Journal of Chemical Theory and Computation},
  volume  = {8},
  number  = {3},
  pages   = {1153--1163},
  year    = {2012},
  doi     = {10.1021/ct200722n}
}""", doi="10.1021/ct200722n", source=_SK_README),
)
SK_SET_CITATIONS: dict[str, tuple[Citation, ...]] = {
    "mio": _MIO, "3ob": _THREEOB, "pbc": _PBC, "matsci": _MATSCI, "ob2": _OB2, "tiorg": _TIORG, "trans3d": _TRANS3D,
    "halorg": _HALORG, "magsil": _MAGSIL, "znorg": _ZNORG, "auorg": _AUORG, "siband": _SIBAND, "borg": _BORG,
}
