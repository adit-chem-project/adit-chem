
from __future__ import annotations

import math
from pathlib import Path

import pydantic

from adit.lang import L
from adit.validate_types import friendly_pydantic
from adit.spec import (AtomsData, BandSettings, CalculationSpec, Cp2kMethod, DftbMethod, EspressoMethod, GromacsMethod, KPoints, LammpsMethod,
                        MDSettings, MlipMethod, OrcaMethod, Runtime, Structure, Task, VaspMethod, XtbMethod)
from adit.web import codefields as CF
from adit.structure import StructureError, build_structure
from adit.textparse import parse_constraints, parse_element_map, parse_extra_incar, parse_extra_namelist

CODES = {"dftbplus": "DFTB+", "vasp": "VASP", "xtb": "xtb (GFN-xTB)", "espresso": "Quantum ESPRESSO (pw.x)", "orca": "ORCA",
         "cp2k": "CP2K", "lammps": "LAMMPS", "gromacs": "GROMACS", "mlip": "MACE / CHGNet"}
NO_KPOINTS = ("lammps", "gromacs", "mlip")
TYPES = {"single_point": ("一点計算", "Single point"), "geometry_optimization": ("構造最適化", "Geometry optimization"),
         "molecular_dynamics": ("分子動力学", "Molecular dynamics"), "vibrations": ("振動解析", "Vibrations"),
         "band_structure": ("バンド計算", "Band structure")}
RELAX_CELL = {"no": ("動かさない", "fixed"), "shape_and_volume": ("形と体積", "shape and volume"), "volume_only": ("体積だけ", "volume only")}
THERMOSTATS = {"berendsen": "Berendsen", "andersen": "Andersen", "nose_hoover": "Nosé-Hoover", "langevin": "Langevin", "csvr": "CSVR"}
KP_MODES = {"gamma": ("Γ 点のみ", "Γ only"), "mesh": ("分割数を指定", "mesh"), "density": ("密度から決める", "from density")}
GFN = {"2": "GFN2-xTB", "1": "GFN1-xTB", "0": "GFN0-xTB", "ff": "GFN-FF"}
IBRION = {2: "2: CG", 1: "1: RMM-DIIS", 3: "3: damped MD"}


def code_label(code: str) -> str:
    if code == "mlip":
        return L("機械学習ポテンシャル (MACE・CHGNet)", "Machine-learning potential (MACE, CHGNet)")
    return CODES.get(code, code)


class FormError(Exception):
    pass


def _f(v: str | None, default: float) -> float:
    v = (v or "").strip()
    if not v:
        return default
    try:
        x = float(v)
    except ValueError as ex:
        raise FormError(L(f"数値として読めません: {v!r}", f"not a number: {v!r}")) from ex
    if not math.isfinite(x):
        raise FormError(L(f"有限の数を入れてください: {v!r}", f"enter a finite number: {v!r}"))
    return x


def _i(v: str | None, default: int) -> int:
    v = (v or "").strip()
    if not v:
        return default
    try:
        x = float(v)
    except ValueError as ex:
        raise FormError(L(f"整数として読めません: {v!r}", f"not an integer: {v!r}")) from ex
    if not math.isfinite(x):
        raise FormError(L(f"有限の整数を入れてください: {v!r}", f"enter a finite integer: {v!r}"))
    return int(x)


def _pydantic_text(ex: ValueError) -> str:
    if isinstance(ex, pydantic.ValidationError):
        return friendly_pydantic(ex)
    return str(ex)


def _check_band_path(atoms, path: str) -> None:
    from ase.dft.kpoints import parse_path_string

    if not path or not any(atoms.pbc):
        return
    try:
        names = [n for seg in parse_path_string(path) for n in seg]
        known = atoms.cell.bandpath(npoints=0).special_points
    except Exception:
        return
    bad = sorted({n for n in names if n not in known})
    if bad:
        raise FormError(L(f"バンドの経路: この格子に無い点の名前です: {', '.join(bad)} (使える名前: {', '.join(sorted(known))})",
                          f"band path: no such point for this lattice: {', '.join(bad)} (available: {', '.join(sorted(known))})"))


def _b(v: str | None) -> bool:
    return (v or "").strip().lower() in ("1", "on", "true", "yes")


def _s(v: str | None, default: str = "") -> str:
    return (v or default).strip()


def source_ref(f: dict[str, str], source: str) -> str:
    if source == "preset":
        ref = _s(f.get("preset"))
    elif source == "smiles":
        ref = _s(f.get("smiles"))
    elif source == "file":
        ref = _s(f.get("file_path"))
    elif source == "fetch":
        ref = _s(f.get("fetch_file"))
    elif source == "bulk":
        parts = [_s(f.get("bulk_el"), "Si")]
        if _s(f.get("bulk_struct")):
            parts.append(_s(f.get("bulk_struct")))
        a = _f(f.get("bulk_a"), 0.0)
        if a > 0:
            parts.append(f"{a:g}")
        if _b(f.get("bulk_cubic")):
            parts.append("cubic")
        ref = " ".join(parts)
    elif source == "surface":
        n = "x".join(str(_i(f.get(k), 2 if k != "surf_nz" else 3)) for k in ("surf_nx", "surf_ny", "surf_nz"))
        ref = f"{_s(f.get('surf_facet'), 'fcc111')} {_s(f.get('surf_el'), 'Al')} {n} vacuum={_f(f.get('surf_vac'), 10.0):g}"
    elif source == "mixture":
        from adit.mixture import MixtureError, MixtureSpec, parse_mixture_text

        try:
            comps = parse_mixture_text(f.get("mixture") or "")
        except MixtureError as ex:
            raise FormError(str(ex)) from ex
        edge = _f(f.get("mix_edge"), 0.0) if _s(f.get("mix_box_mode"), "density") == "edge" else 0.0
        ref = MixtureSpec(components=comps, box_a=edge, density_g_cm3=_f(f.get("mix_density"), 1.0), min_distance=_f(f.get("mix_min_dist"), 2.0),
                          seed=_i(f.get("mix_seed"), 0)).to_ref() if comps else ""
    else:
        raise FormError(L(f"未知の構造の指定: {source!r}", f"unknown structure source: {source!r}"))
    return ref


_SOURCE_LABELS = {
    "preset": ("プリセット", "Preset"), "smiles": ("SMILES", "SMILES"), "file": ("ファイル", "File"),
    "bulk": ("バルク", "Bulk"), "surface": ("スラブ", "Slab"), "mixture": ("溶液・混合物", "Solution / mixture"),
    "recipe": ("組み立て手順", "Recipe"), "2d": ("2 次元材料", "2D material"), "fetch": ("データベースから取得", "Fetch from a database"),
}


def fetch_record(f: dict[str, str]) -> dict | None:
    """The fetch record kept in the hidden fetch_record field (JSON), or None."""
    import json

    text = _s(f.get("fetch_record"))
    if not text:
        return None
    try:
        rec = json.loads(text)
    except ValueError:
        return None
    return rec if isinstance(rec, dict) and rec.get("file") else None


def as_file_source(f: dict[str, str]) -> tuple[dict[str, str], dict | None]:
    """A 'fetch' form is a 'file' form whose path is the fetched file; returns the form to use and the fetch record."""
    if _s(f.get("source")) != "fetch":
        return f, None
    rec = fetch_record(f)
    if rec is None:
        raise FormError(L("データベースから取得: 「取得」を押して構造を取得してください", "Fetch from a database: press Fetch to get the structure first"))
    return {**f, "source": "file", "file_path": _s(f.get("fetch_file")) or str(rec["file"])}, rec


def _source_label(source: str) -> str:
    ja, en = _SOURCE_LABELS.get(source, (source, source))
    return L(ja, en)


def structure_from_form(f: dict[str, str], state=None) -> Structure:
    from adit.web import recipe_form

    if recipe_form.recipe_mode(f):
        return _recipe_structure(f, state)
    f, fetched = as_file_source(f)
    source = _s(f.get("source"), "preset")
    if _s(f.get("code")) == "gromacs" and _b(f.get("gmx_use_conf")) and _s(f.get("gmx_conf")):
        source, f, fetched = "file", {**f, "file_path": _s(f.get("gmx_conf"))}, None
    ref = source_ref(f, source)
    if not ref:
        raise FormError(L(f"{_source_label(source)}: 何も指定されていません", f"{_source_label(source)}: nothing given"))
    charge = _i(f.get("charge"), 0)
    if source == "mixture":
        from adit.mixture import MixtureSpec
        charge = MixtureSpec.from_ref(ref).total_charge()
    try:
        ost = getattr(state, "origin_structure", None)
        if ost is not None and source == "file" and ref == ost.source_ref:
            st = ost.model_copy(update={"charge": charge, "multiplicity": _i(f.get("multiplicity"), 1), "fixed_atoms": [], "fixed_axes": {}})
        else:
            st = build_structure(source, ref, charge=charge, multiplicity=_i(f.get("multiplicity"), 1))
        atoms = st.atoms.to_ase()
        if _b(f.get("box")) and not any(atoms.pbc):
            size = _f(f.get("box_size"), 15.0)
            atoms.set_cell([size, size, size]); atoms.center(); atoms.pbc = True
            st = st.model_copy(update={"atoms": AtomsData.from_ase(atoms)})
        fixed_atoms, fixed_axes = parse_constraints(_s(f.get("fixed")), len(atoms))
        return st.model_copy(update={"fixed_atoms": fixed_atoms, "fixed_axes": fixed_axes, "fetched": fetched})
    except (StructureError, ValueError) as ex:
        raise FormError(_pydantic_text(ex)) from ex


def _recipe_structure(f: dict[str, str], state) -> Structure:
    from adit.web import recipe_form

    charge, mult = _i(f.get("charge"), 0), _i(f.get("multiplicity"), 1)
    try:
        f, fetched = as_file_source(f)
        st, new, rec = recipe_form.structure_for(f, getattr(state, "built", None), charge, mult)
        if new is not None and state is not None:
            state.built = new
        atoms = st.atoms.to_ase()
        upd: dict = {"charge": charge, "multiplicity": mult, "fetched": fetched}
        if _b(f.get("box")) and not any(atoms.pbc):
            size = _f(f.get("box_size"), 15.0)
            atoms.set_cell([size, size, size]); atoms.center(); atoms.pbc = True
            upd["atoms"] = AtomsData.from_ase(atoms)
        if not any(s.op == "fix" for s in rec.steps):
            upd["fixed_atoms"], upd["fixed_axes"] = parse_constraints(_s(f.get("fixed")), len(atoms))
        return st.model_copy(update=upd)
    except (StructureError, ValueError) as ex:
        raise FormError(_pydantic_text(ex)) from ex


HUB_ROWS_MAX = 20
MAG_PREFIXES = ("vasp_mag__", "qe_mag__", "cp_mag__")
HUB_PREFIXES = ("vasp_hub", "qe_hub", "cp_hub")
PREP_PREFIXES = MAG_PREFIXES + HUB_PREFIXES + ("xtb_solvent__", "orca_solvent__")


def drop_prep_fields(f: dict[str, str]) -> dict[str, str]:
    return {k: v for k, v in f.items() if not k.startswith(PREP_PREFIXES)}


def hub_rows(f: dict[str, str], prefix: str) -> list[tuple[str, str, str, str]]:
    return [(f.get(f"{prefix}{i}_el", ""), f.get(f"{prefix}{i}_orb", ""), f.get(f"{prefix}{i}_u", ""), f.get(f"{prefix}{i}_j", ""))
            for i in range(1, HUB_ROWS_MAX + 1) if f"{prefix}{i}_el" in f]


def _hubbard(f: dict[str, str], prefix: str) -> dict:
    try:
        return CF.hubbard_from_rows(hub_rows(f, prefix))
    except ValueError as ex:
        raise FormError(str(ex)) from ex


def _element_values(f: dict[str, str], prefix: str, what: str) -> dict[str, float]:
    try:
        return CF.element_values_from_rows([(k[len(prefix):], v) for k, v in f.items() if k.startswith(prefix)], what)
    except ValueError as ex:
        raise FormError(str(ex)) from ex


def _solvent(f: dict[str, str], field: str, model: str, key: str) -> str:
    return "" if model == "none" else _s(f.get(f"{field}__{key}"))


def _prep_fields_from_spec(m) -> dict[str, str]:
    f: dict[str, str] = {}

    def hub(prefix: str, table: dict) -> None:
        for i, (el, orb, u, j) in enumerate(CF.hubbard_rows(table), start=1):
            f.update({f"{prefix}{i}_el": el, f"{prefix}{i}_orb": orb, f"{prefix}{i}_u": u, f"{prefix}{i}_j": j})

    if isinstance(m, DftbMethod):
        f["dftb_solv_file"] = m.solvation_param_file
    elif isinstance(m, VaspMethod):
        f.update({f"vasp_mag__{e}": f"{v:g}" for e, v in m.magmom_by_element.items()}, ldau_type=str(m.ldau_type))
        hub("vasp_hub", m.hubbard)
    elif isinstance(m, XtbMethod):
        f["xtb_solvation"] = m.solvation
        if m.solvation != "none":
            f[f"xtb_solvent__{m.solvation}_{m.gfn}"] = m.solvent
    elif isinstance(m, EspressoMethod):
        f.update({f"qe_mag__{e}": f"{v:g}" for e, v in m.starting_magnetization.items()}, hubbard_projector=m.hubbard_projector)
        hub("qe_hub", m.hubbard)
    elif isinstance(m, OrcaMethod):
        f["orca_solvation"] = m.solvation
        if m.solvation != "none":
            f[f"orca_solvent__{m.solvation}"] = m.solvent
    elif isinstance(m, Cp2kMethod):
        f.update({f"cp_mag__{e}": f"{v:g}" for e, v in m.magnetization_by_element.items()}, cp_plus_u=m.plus_u_method,
                 cp_sccs=f"{m.sccs_relative_permittivity:g}")
        hub("cp_hub", m.hubbard)
    return f


def method_from_form(f: dict[str, str]):
    code = _s(f.get("code"), "dftbplus")
    if code == "dftbplus":
        disp = _s(f.get("dispersion"), "none")
        d3 = {k: _f(f.get(f"d3_{k}"), 0.0) for k in ("s6", "s8", "a1", "a2")} if disp == "dftd3" else None
        return DftbMethod(sk_set=_s(f.get("sk_set")), scc=_b(f.get("scc")), scc_tolerance=_f(f.get("scc_tol"), 1e-5),
                          max_scc_iterations=_i(f.get("max_scc"), 100), third_order=_b(f.get("third")), dispersion=disp, d3_params=d3,
                          filling_temperature=_f(f.get("etemp"), 0.0), solvation_param_file=_s(f.get("dftb_solv_file")))
    if code == "vasp":
        try:
            extra = parse_extra_incar(f.get("extra_incar") or "")
            potcar = parse_element_map(f.get("potcar_map") or "")
        except ValueError as ex:
            raise FormError(str(ex)) from ex
        prep = dict(magmom_by_element=_element_values(f, "vasp_mag__", "MAGMOM"), hubbard=_hubbard(f, "vasp_hub"), ldau_type=_i(f.get("ldau_type"), 2))
        magmom_s = _s(f.get("magmom"))
        magmom = [_f(x, 0.0) for x in magmom_s.split()] if magmom_s else None
        ivdw_s = _s(f.get("ivdw"), "none")
        return VaspMethod(potcar_set=_s(f.get("potcar_set"), "potpaw_PBE"), potcar=potcar, binary=_s(f.get("binary"), "std"),
                          encut=_f(f.get("encut"), 0.0), ediff=_f(f.get("ediff"), 1e-4), nelm=_i(f.get("nelm"), 60), ismear=_i(f.get("ismear"), 0),
                          sigma=_f(f.get("sigma"), 0.1), ispin=_i(f.get("ispin"), 1), magmom=magmom,
                          ivdw=None if ivdw_s == "none" else _i(ivdw_s, 0), ibrion=_i(f.get("ibrion"), 2), algo=_s(f.get("algo"), "Normal"),
                          prec=_s(f.get("prec"), "Normal"), lreal=_s(f.get("lreal"), "Auto"), extra_incar=extra,
                          nelmin=_i(f.get("nelmin"), 0), lasph=_b(f.get("lasph")), lmaxmix=_i(f.get("lmaxmix"), 0),
                          nbands=_i(f.get("nbands"), 0), isym=_i(f.get("isym"), 0) if _s(f.get("isym")) else None,
                          idipol=_i(f.get("idipol"), 0), ldipol=_b(f.get("ldipol")), dipol=_s(f.get("dipol")),
                          kpoints_centering=_s(f.get("kpoints_centering"), "monkhorst-pack"), **prep)
    if code == "xtb":
        model, gfn = _s(f.get("xtb_solvation"), "none"), _s(f.get("gfn"), "2")
        return XtbMethod(gfn=gfn, accuracy=_f(f.get("accuracy"), 1.0), etemp=_f(f.get("xtb_etemp"), 300.0),
                         max_iterations=_i(f.get("xtb_max_iter"), 250), opt_level=_s(f.get("opt_level"), "normal"),
                         solvation=model, solvent=_solvent(f, "xtb_solvent", model, f"{model}_{gfn}"))
    if code == "espresso":
        try:
            extra = parse_extra_namelist(f.get("qe_extra") or "")
            pseudo = parse_element_map(f.get("pseudo_map") or "")
        except ValueError as ex:
            raise FormError(str(ex)) from ex
        return EspressoMethod(pseudo_set=_s(f.get("pseudo_set")), pseudo=pseudo, ecutwfc=_f(f.get("ecutwfc"), 0.0), ecutrho=_f(f.get("ecutrho"), 0.0),
                              conv_thr=_f(f.get("conv_thr"), 1e-6), electron_maxstep=_i(f.get("qe_maxstep"), 100), mixing_beta=_f(f.get("mixing_beta"), 0.7),
                              occupations=_s(f.get("occupations"), "fixed"), smearing=_s(f.get("smearing"), "gaussian"), degauss=_f(f.get("degauss"), 0.0),
                              nspin=_i(f.get("nspin"), 1), input_dft=_s(f.get("input_dft")), extra=extra,
                              starting_magnetization=_element_values(f, "qe_mag__", "starting_magnetization"), hubbard=_hubbard(f, "qe_hub"),
                              hubbard_projector=_s(f.get("hubbard_projector"), "atomic"), assume_isolated=_s(f.get("assume_isolated")))
    if code == "orca":
        model = _s(f.get("orca_solvation"), "none")
        return OrcaMethod(method=_s(f.get("orca_method"), "HF"), basis=f["orca_basis"].strip() if "orca_basis" in f else "def2-SVP", scf_convergence=_s(f.get("scf_conv"), "NormalSCF"),
                          scf_maxiter=_i(f.get("scf_maxiter"), 125), maxcore_mb=_i(f.get("maxcore"), 0), extra_keywords=_s(f.get("extra_kw")),
                          extra_blocks=(f.get("extra_blocks") or "").strip(), solvation=model, solvent=_solvent(f, "orca_solvent", model, model),
                          ts_search=_b(f.get("orca_ts")), ts_calc_hess=_b(f.get("orca_calc_hess")), ts_recalc_hess=_i(f.get("orca_recalc_hess"), 0),
                          ts_freq=_b(f.get("orca_freq")), irc=_b(f.get("orca_irc")), irc_max_iter=_i(f.get("irc_max_iter"), 0),
                          irc_direction=_s(f.get("irc_direction"), "both"), goat=_b(f.get("orca_goat")),
                          docker_guest_file=_s(f.get("orca_docker_guest_file")),
                          docker_assume_neutral_singlet=_b(f.get("orca_docker_assume_neutral_singlet")))
    if code == "mlip":
        return MlipMethod(model_family=_s(f.get("mlip_family")), model=_s(f.get("mlip_model")), device=_s(f.get("mlip_device")),
                          dtype=_s(f.get("mlip_dtype")), dispersion=_b(f.get("mlip_d3")), seed=_i(f.get("mlip_seed"), 12345))
    try:
        if code == "cp2k":
            return Cp2kMethod(xc=_s(f.get("cp_xc")), basis_file=_s(f.get("cp_basis_file"), "BASIS_MOLOPT"), potential_file=_s(f.get("cp_potential_file"), "GTH_POTENTIALS"),
                              basis=parse_element_map(f.get("cp_basis_map") or ""), potential=parse_element_map(f.get("cp_potential_map") or ""),
                              dispersion=_s(f.get("cp_dispersion"), "none"), cutoff_ry=_f(f.get("cp_cutoff"), 0.0), rel_cutoff_ry=_f(f.get("cp_rel_cutoff"), 0.0),
                              eps_scf=_f(f.get("cp_eps_scf"), 1e-5), max_scf=_i(f.get("cp_max_scf"), 50), uks=_b(f.get("cp_uks")),
                              poisson_solver=_s(f.get("cp_poisson")), isolated_box_ang=_f(f.get("cp_box"), 0.0),
                              extra_sections=CF.parse_sections(f.get("cp_extra") or ""),
                              magnetization_by_element=_element_values(f, "cp_mag__", "MAGNETIZATION"), hubbard=_hubbard(f, "cp_hub"),
                              plus_u_method=_s(f.get("cp_plus_u"), "MULLIKEN"), sccs_relative_permittivity=_f(f.get("cp_sccs"), 0.0),
                              ot=_b(f.get("cp_ot")), ot_minimizer=_s(f.get("cp_ot_min")), ot_preconditioner=_s(f.get("cp_ot_pre")),
                              ot_extra=(f.get("cp_ot_extra") or "").strip(), surface_dipole_correction=_b(f.get("cp_surf_dip")),
                              surf_dip_dir=_s(f.get("cp_surf_dir")))
        if code == "lammps":
            return LammpsMethod(units=_s(f.get("lmp_units")), atom_style=_s(f.get("lmp_atom_style"), "atomic"), data_file=_s(f.get("lmp_data_file")),
                                type_elements=CF.parse_words(f.get("lmp_type_elements") or ""), pair_style=_s(f.get("lmp_pair_style")),
                                pair_coeff=(f.get("lmp_pair_coeff") or "").strip(), potential_files=CF.parse_lines(f.get("lmp_files") or ""),
                                style_commands=(f.get("lmp_style_cmds") or "").strip(), extra_commands=(f.get("lmp_extra_cmds") or "").strip(),
                                seed=_i(f.get("lmp_seed"), 12345))
        if code == "gromacs":
            return GromacsMethod(topology_file=_s(f.get("gmx_top")), structure_file=_s(f.get("gmx_conf")), coulombtype=_s(f.get("gmx_coulomb"), "Cut-off"),
                                 rcoulomb_nm=_f(f.get("gmx_rcoulomb"), 1.0), rvdw_nm=_f(f.get("gmx_rvdw"), 1.0), constraints=_s(f.get("gmx_constraints"), "none"),
                                 pcoupl=_s(f.get("gmx_pcoupl")), compressibility_per_bar=_f(f.get("gmx_compress"), 0.0), define=_s(f.get("gmx_define")),
                                 checkpoint_file=_s(f.get("gmx_cpt")), gen_seed=_i(f.get("gmx_gen_seed"), -1), extra_mdp=CF.parse_mdp(f.get("gmx_extra") or ""))
    except ValueError as ex:
        if isinstance(ex, pydantic.ValidationError):
            raise
        raise FormError(str(ex)) from ex
    raise FormError(L(f"未知のコード: {code!r}", f"unknown code: {code!r}"))


def task_from_form(f: dict[str, str]) -> Task:
    md = MDSettings(ensemble=_s(f.get("ensemble"), "NVT"), thermostat=_s(f.get("thermostat"), "berendsen"), temperature_k=_f(f.get("temperature"), 300.0),
                    timestep_fs=_f(f.get("timestep"), 1.0), steps=_i(f.get("md_steps"), 1000), dump_interval=_i(f.get("dump"), 10),
                    coupling_time_fs=_f(f.get("coupling"), 100.0), pressure_bar=_f(f.get("pressure"), 1.0), barostat_time_fs=_f(f.get("barostat_time"), 1000.0))
    bands = BandSettings(path=_s(f.get("band_path")), npoints=_i(f.get("band_npoints"), 60), empty_bands=_i(f.get("band_empty"), 4))
    return Task(type=_s(f.get("task_type"), "single_point"), optimizer=_s(f.get("optimizer"), "Rational"), max_steps=_i(f.get("max_steps"), 200),
                force_tolerance_ev_per_ang=_f(f.get("force_tol"), Task().force_tolerance_ev_per_ang), relax_cell=_s(f.get("relax_cell"), "no"),
                md=md, bands=bands)


def kpoints_from_form(f: dict[str, str]) -> KPoints:
    shift = _f(f.get("kp_shift"), 0.0)
    return KPoints(mode=_s(f.get("kp_mode"), "gamma"), mesh=tuple(_i(f.get(k), 1) for k in ("k1", "k2", "k3")),
                   shift=(shift, shift, shift), density=_f(f.get("kp_density"), 0.0))


def runtime_from_form(f: dict[str, str]) -> Runtime:
    return Runtime(profile=_s(f.get("profile"), "local"), nodes=_i(f.get("nodes"), 1), ncpus=_i(f.get("ncpus"), 8), mpiprocs=_i(f.get("mpiprocs"), 1),
                   omp_threads=_i(f.get("omp"), 8), walltime=_s(f.get("walltime"), "01:00:00"), job_name=_s(f.get("job_name"), "adit"))


def spec_from_form(f: dict[str, str], state=None) -> CalculationSpec:
    st = structure_from_form(f, state)
    try:
        uses_k = _s(f.get("code"), "dftbplus") not in NO_KPOINTS
        spec = CalculationSpec(structure=st, method=method_from_form(f), kpoints=kpoints_from_form(f) if st.periodic and uses_k else None,
                               task=task_from_form(f), runtime=runtime_from_form(f))
    except ValueError as ex:
        raise FormError(_pydantic_text(ex)) from ex
    if spec.task.type == "band_structure":
        _check_band_path(spec.atoms, spec.task.bands.path)
    return spec


def form_from_spec(spec: CalculationSpec) -> dict[str, str]:
    from adit.web import recipe_form

    st, m, t, r = spec.structure, spec.method, spec.task, spec.runtime
    f: dict[str, str] = {"source": st.source, "charge": str(st.charge), "multiplicity": str(st.multiplicity), "recipe_steps": "[]",
                         "base_extra": "{}", "file_sha_path": "", "file_sha256": "", "recipe_q": ""}
    src, sref, rec = st.source, st.source_ref, None
    if st.source == "recipe":
        from adit.builder import Recipe
        try:
            rec = Recipe.from_ref(st.source_ref)
        except StructureError:
            rec = None
        if rec is not None:
            f.update(recipe_form.form_from_recipe(rec))
            src, sref = rec.base.source, rec.base.ref if isinstance(rec.base.ref, str) else ""
    if src == "preset":
        f["preset"] = sref
    elif src == "smiles":
        f["smiles"] = sref
    elif src == "file" and st.fetched and str(st.fetched.get("file", "")) == sref:
        import json

        f.update(source="fetch", fetch_db=str(st.fetched.get("database", "pubchem")), fetch_query=str(st.fetched.get("query", "")),
                 fetch_file=sref, fetch_record=json.dumps(st.fetched, ensure_ascii=False))
    elif src == "file":
        f["file_path"] = sref
    elif src == "bulk":
        toks = sref.split()
        f["bulk_el"] = toks[0] if toks else ""
        f["bulk_cubic"] = "on" if "cubic" in toks else ""
        rest = [x for x in toks[1:] if x != "cubic"]
        for x in rest:
            try:
                f["bulk_a"] = f"{float(x):g}"
            except ValueError:
                f["bulk_struct"] = x
    elif src == "surface":
        toks = sref.split()
        if len(toks) >= 3:
            f["surf_facet"], f["surf_el"] = toks[0], toks[1]
            n = toks[2].split("x")
            if len(n) == 3:
                f["surf_nx"], f["surf_ny"], f["surf_nz"] = n
        for x in toks[3:]:
            if x.startswith("vacuum="):
                f["surf_vac"] = x[7:]
    elif src == "mixture":
        from adit.mixture import MixtureError, MixtureSpec, mixture_text

        try:
            mx = MixtureSpec.from_ref(sref)
            f.update(mixture=mixture_text(mx.components), mix_box_mode="edge" if mx.box_a > 0 else "density", mix_edge=f"{mx.box_a:g}",
                     mix_density=f"{mx.density_g_cm3:g}", mix_min_dist=f"{mx.min_distance:g}", mix_seed=str(mx.seed))
        except MixtureError:
            pass
    if st.source in ("preset", "smiles", "file") and st.periodic:
        f["box"], f["box_size"] = "on", f"{float(st.atoms.to_ase().cell.lengths()[0]):.15g}"
    else:
        f["box"] = ""
    parts =[str(i + 1) for i in st.fixed_atoms]
    parts += [f"{int(k) + 1}:{''.join(c for c, mv in zip('xyz', v) if not mv)}" for k, v in st.fixed_axes.items()]
    f["fixed"] = ",".join(parts)
    if rec is not None and any(s.op == "fix" for s in rec.steps):
        f["fixed"] = ""
    f["code"] = m.code
    if isinstance(m, DftbMethod):
        f.update(sk_set=m.sk_set, scc="on" if m.scc else "", scc_tol=f"{m.scc_tolerance:g}", max_scc=str(m.max_scc_iterations),
                 third="on" if m.third_order else "", dispersion=m.dispersion, etemp=f"{m.filling_temperature:g}")
        for k, v in (m.d3_params or {}).items():
            f[f"d3_{k}"] = f"{v:g}"
    elif isinstance(m, VaspMethod):
        f.update(potcar_set=m.potcar_set, binary=m.binary, encut=f"{m.encut:g}", ediff=f"{m.ediff:g}", nelm=str(m.nelm), ismear=str(m.ismear),
                 sigma=f"{m.sigma:g}", ispin=str(m.ispin), magmom=" ".join(f"{x:g}" for x in m.magmom) if m.magmom else "",
                 ivdw="none" if m.ivdw is None else str(m.ivdw), ibrion=str(m.ibrion), algo=m.algo, prec=m.prec, lreal=str(m.lreal),
                 nelmin=str(m.nelmin), lasph="on" if m.lasph else "", lmaxmix=str(m.lmaxmix), nbands=str(m.nbands),
                 isym="" if m.isym is None else str(m.isym), idipol=str(m.idipol), ldipol="on" if m.ldipol else "", dipol=m.dipol,
                 kpoints_centering=m.kpoints_centering,
                 extra_incar="\n".join(f"{k} = {'.TRUE.' if v is True else '.FALSE.' if v is False else v}" for k, v in m.extra_incar.items()),
                 potcar_map="\n".join(f"{e} = {n}" for e, n in m.potcar.items()))
    elif isinstance(m, XtbMethod):
        f.update(gfn=m.gfn, accuracy=f"{m.accuracy:g}", xtb_etemp=f"{m.etemp:g}", xtb_max_iter=str(m.max_iterations), opt_level=m.opt_level)
    elif isinstance(m, EspressoMethod):
        f.update(pseudo_set=m.pseudo_set, ecutwfc=f"{m.ecutwfc:g}", ecutrho=f"{m.ecutrho:g}", conv_thr=f"{m.conv_thr:g}", qe_maxstep=str(m.electron_maxstep),
                 mixing_beta=f"{m.mixing_beta:g}", occupations=m.occupations, smearing=m.smearing, degauss=f"{m.degauss:g}", nspin=str(m.nspin),
                 input_dft=m.input_dft, pseudo_map="\n".join(f"{e} = {n}" for e, n in m.pseudo.items()),
                 assume_isolated=m.assume_isolated,
                 qe_extra="\n".join(f"{ns}.{k} = {v}" for ns, d in m.extra.items() for k, v in d.items()))
    elif isinstance(m, OrcaMethod):
        f.update(orca_method=m.method, orca_basis=m.basis, scf_conv=m.scf_convergence, scf_maxiter=str(m.scf_maxiter), maxcore=str(m.maxcore_mb),
                 extra_kw=m.extra_keywords, extra_blocks=m.extra_blocks,
                 orca_ts="on" if m.ts_search else "", orca_calc_hess="on" if m.ts_calc_hess else "", orca_recalc_hess=str(m.ts_recalc_hess),
                 orca_freq="on" if m.ts_freq else "", orca_irc="on" if m.irc else "", irc_max_iter=str(m.irc_max_iter), irc_direction=m.irc_direction,
                 orca_goat="on" if m.goat else "", orca_docker_guest_file=m.docker_guest_file,
                 orca_docker_assume_neutral_singlet="on" if m.docker_assume_neutral_singlet else "")
    elif isinstance(m, MlipMethod):
        f.update(mlip_family=m.model_family, mlip_model=m.model, mlip_device=m.device, mlip_dtype=m.dtype,
                 mlip_d3="on" if m.dispersion else "", mlip_seed=str(m.seed))
    elif isinstance(m, Cp2kMethod):
        f.update(cp_xc=m.xc, cp_basis_file=m.basis_file, cp_potential_file=m.potential_file, cp_basis_map="\n".join(f"{e} = {n}" for e, n in m.basis.items()),
                 cp_potential_map="\n".join(f"{e} = {n}" for e, n in m.potential.items()), cp_dispersion=m.dispersion, cp_cutoff=f"{m.cutoff_ry:g}",
                 cp_rel_cutoff=f"{m.rel_cutoff_ry:g}", cp_eps_scf=f"{m.eps_scf:g}", cp_max_scf=str(m.max_scf), cp_uks="on" if m.uks else "",
                 cp_poisson=m.poisson_solver, cp_box=f"{m.isolated_box_ang:g}", cp_extra=CF.sections_text(m.extra_sections))
        f.update(cp_ot="on" if m.ot else "", cp_ot_min=m.ot_minimizer, cp_ot_pre=m.ot_preconditioner, cp_ot_extra=m.ot_extra,
                 cp_surf_dip="on" if m.surface_dipole_correction else "", cp_surf_dir=m.surf_dip_dir)
    elif isinstance(m, LammpsMethod):
        f.update(lmp_units=m.units, lmp_atom_style=m.atom_style, lmp_data_file=m.data_file, lmp_type_elements=" ".join(m.type_elements),
                 lmp_pair_style=m.pair_style, lmp_pair_coeff=m.pair_coeff, lmp_files="\n".join(m.potential_files), lmp_style_cmds=m.style_commands,
                 lmp_extra_cmds=m.extra_commands, lmp_seed=str(m.seed))
    elif isinstance(m, GromacsMethod):
        f.update(gmx_top=m.topology_file, gmx_conf=m.structure_file, gmx_coulomb=m.coulombtype, gmx_rcoulomb=f"{m.rcoulomb_nm:g}", gmx_rvdw=f"{m.rvdw_nm:g}",
                 gmx_constraints=m.constraints, gmx_pcoupl=m.pcoupl, gmx_compress=f"{m.compressibility_per_bar:g}", gmx_define=m.define,
                 gmx_cpt=m.checkpoint_file, gmx_gen_seed=str(m.gen_seed), gmx_extra=CF.mdp_text(m.extra_mdp),
                 gmx_use_conf="on" if st.source == "file" and st.source_ref == m.structure_file and m.structure_file else "")
    f.update(_prep_fields_from_spec(m))
    f.update(task_type=t.type, optimizer=t.optimizer, max_steps=str(t.max_steps), force_tol=f"{t.force_tolerance_ev_per_ang:g}", relax_cell=t.relax_cell,
             ensemble=t.md.ensemble, thermostat=t.md.thermostat, temperature=f"{t.md.temperature_k:g}", timestep=f"{t.md.timestep_fs:g}",
             md_steps=str(t.md.steps), dump=str(t.md.dump_interval), coupling=f"{t.md.coupling_time_fs:g}", pressure=f"{t.md.pressure_bar:g}",
             barostat_time=f"{t.md.barostat_time_fs:g}", band_path=t.bands.path, band_npoints=str(t.bands.npoints), band_empty=str(t.bands.empty_bands))
    kp = spec.kpoints or KPoints()
    f.update(kp_mode=kp.mode, k1=str(kp.mesh[0]), k2=str(kp.mesh[1]), k3=str(kp.mesh[2]), kp_shift=f"{kp.shift[0]:g}", kp_density=f"{kp.density:g}")
    f.update(profile=r.profile, nodes=str(r.nodes), ncpus=str(r.ncpus), mpiprocs=str(r.mpiprocs), omp=str(r.omp_threads), walltime=r.walltime, job_name=r.job_name)
    return f


def default_form(cfg_profile: str = "local", output_dir: str | Path = "") -> dict[str, str]:
    f = {"source": "preset", "preset": "H2O", "charge": "0", "multiplicity": "1", "box_size": "15", "bulk_el": "Si", "surf_facet": "fcc111",
         "surf_el": "Al", "surf_nx": "2", "surf_ny": "2", "surf_nz": "3", "surf_vac": "10",
         "mixture": "H2O * 32 label=water", "mix_box_mode": "density", "mix_density": "1.0", "mix_edge": "15", "mix_min_dist": "2.0", "mix_seed": "0",
         "code": "dftbplus", "scc": "on", "scc_tol": "1e-05", "max_scc": "100", "dispersion": "none", "etemp": "0",
         "d3_s6": "0", "d3_s8": "0", "d3_a1": "0", "d3_a2": "0",
         "potcar_set": "potpaw_PBE", "binary": "std", "encut": "0", "ediff": "0.0001", "nelm": "60", "ismear": "0", "sigma": "0.1", "ispin": "1",
         "ivdw": "none", "ibrion": "2", "algo": "Normal", "prec": "Normal", "lreal": "Auto", "nelmin": "0", "lmaxmix": "0",
         "nbands": "0", "idipol": "0", "kpoints_centering": "monkhorst-pack",
         "gfn": "2", "accuracy": "1", "xtb_etemp": "300", "xtb_max_iter": "250", "opt_level": "normal",
         "ecutwfc": "0", "ecutrho": "0", "conv_thr": "1e-06", "qe_maxstep": "100", "mixing_beta": "0.7", "occupations": "fixed", "smearing": "gaussian",
         "degauss": "0", "nspin": "1",
         "orca_method": "HF", "orca_basis": "def2-SVP", "scf_conv": "NormalSCF", "scf_maxiter": "125", "maxcore": "0",
         "orca_calc_hess": "on", "orca_freq": "on", "orca_recalc_hess": "0", "irc_max_iter": "0", "irc_direction": "both",
         "mlip_family": "", "mlip_model": "", "mlip_device": "", "mlip_dtype": "", "mlip_seed": "12345",
         "batch_kind": "compare", "cmp_kind": "adsorption", "cmp_box": "as_is", "cmp_rx_rows": "3",
         "conf_seed": "12345", "conf_ff": "MMFF94", "conf_iters": "200", "neb_mode": "native", "neb_interp": "idpp",
         "ph_backend": "auto", "el_c1": "on", "el_c2": "on", "el_c3": "on", "el_c4": "on", "el_c5": "on", "el_c6": "on",
         "ts_mode": "ts", "ts_irc": "on",
         "cp_basis_file": "BASIS_MOLOPT", "cp_potential_file": "GTH_POTENTIALS", "cp_dispersion": "none", "cp_cutoff": "0", "cp_rel_cutoff": "0",
         "cp_eps_scf": "1e-05", "cp_max_scf": "50", "cp_box": "0",
         "lmp_atom_style": "atomic", "lmp_seed": "12345",
         "gmx_coulomb": "Cut-off", "gmx_rcoulomb": "1", "gmx_rvdw": "1", "gmx_constraints": "none", "gmx_compress": "0", "gmx_gen_seed": "12345",
         "xtb_solvation": "none", "orca_solvation": "none", "ldau_type": "2", "hubbard_projector": "atomic", "cp_plus_u": "MULLIKEN", "cp_sccs": "0",
         "dftb_solv_file": "", "cont_dir": "", "cont_vel": "on", "stg_rows": "3",
         "task_type": "single_point", "optimizer": "Rational", "max_steps": "200", "force_tol": f"{Task().force_tolerance_ev_per_ang:g}", "relax_cell": "no",
         "ensemble": "NVT", "thermostat": "berendsen", "temperature": "300", "timestep": "1", "md_steps": "1000", "dump": "10", "coupling": "100",
         "pressure": "1", "barostat_time": "1000", "band_npoints": "60", "band_empty": "4",
         "kp_mode": "gamma", "k1": "1", "k2": "1", "k3": "1", "kp_shift": "0", "kp_density": "0",
         "profile": cfg_profile, "nodes": "1", "ncpus": "8", "mpiprocs": "1", "omp": "8", "walltime": "01:00:00", "job_name": "adit",
         "output_dir": str(output_dir)}
    from adit.web.recipe_form import BASE_DEFAULTS
    f.update(BASE_DEFAULTS)
    return f
