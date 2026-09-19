
from __future__ import annotations

import warnings
from dataclasses import asdict, dataclass

import numpy as np
from ase import Atoms, units

from adit.lang import L

MODELS = ("ideal_gas", "harmonic", "quasi_harmonic", "msrrho")
GEOMETRIES = ("linear", "nonlinear", "monatomic")
EV_J_MOL = 96485.33212
EV_KJ_MOL = EV_J_MOL / 1000


@dataclass
class ThermoOptions:
    model: str | None = None
    temperatures_k: tuple[float, ...] = ()
    pressure_pa: float | None = None
    symmetry_number: int | None = None
    geometry: str | None = None
    spin: float | None = None
    imaginary: str | None = None
    exclude_lowest: int | None = None
    qh_cutoff_cm1: float | None = None
    msrrho_tau_cm1: float | None = None


def missing_inputs(o: ThermoOptions) -> list[str]:
    r = []
    if o.model not in MODELS:
        return [L(f"モデルを選んでください ({' / '.join(MODELS)})", f"choose a model ({' / '.join(MODELS)})")]
    if not o.temperatures_k:
        r.append(L("温度 [K]", "temperature [K]"))
    elif any(t <= 0 for t in o.temperatures_k):
        r.append(L("温度は正の数 [K]", "temperatures must be positive [K]"))
    if o.model == "ideal_gas":
        if o.pressure_pa is None or o.pressure_pa <= 0:
            r.append(L("圧力 [Pa] (正の数)", "pressure [Pa] (positive)"))
        if o.symmetry_number is None or o.symmetry_number < 1:
            r.append(L("回転の対称数 (1 以上の整数)", "rotational symmetry number (integer, 1 or more)"))
        if o.geometry not in GEOMETRIES:
            r.append(L(f"分子の形 ({' / '.join(GEOMETRIES)})", f"geometry ({' / '.join(GEOMETRIES)})"))
        if o.spin is None or o.spin < 0:
            r.append(L("全スピン S (0 以上。一重項 0、二重項 0.5)", "total spin S (0 or more; singlet 0, doublet 0.5)"))
    else:
        if o.exclude_lowest is None or o.exclude_lowest < 0:
            r.append(L("除く低い振動の本数 (0 以上の整数。並進・回転にあたるモードなど)", "number of lowest modes to exclude (integer, 0 or more; e.g. translations and rotations)"))
    if o.model == "quasi_harmonic" and (o.qh_cutoff_cm1 is None or o.qh_cutoff_cm1 <= 0):
        r.append(L("準調和の下限の振動数 [cm⁻¹] (正の数)", "quasi-harmonic cutoff frequency [cm^-1] (positive)"))
    if o.model == "msrrho" and (o.msrrho_tau_cm1 is None or o.msrrho_tau_cm1 <= 0):
        r.append(L("準 RRHO の減衰の τ [cm⁻¹] (正の数)", "quasi-RRHO damping tau [cm^-1] (positive)"))
    if o.imaginary not in (None, "ignore", "stop"):
        r.append(L("虚振動の扱いは ignore か stop", "imaginary-mode handling must be ignore or stop"))
    return r


def _energies(freqs_cm1) -> list[complex]:
    return [complex(f * units.invcm) if f >= 0 else complex(0, -f * units.invcm) for f in freqs_cm1]


def _order(e: list[complex]) -> list[int]:
    return sorted(range(len(e)), key=lambda i: (e[i] ** 2).real)


def _to_cm1(e: complex) -> float:
    return float(e.real / units.invcm) if e.imag == 0 else -float(e.imag / units.invcm)


def compute_thermo(freqs_cm1: list[float], atoms: Atoms | None, energy_ev: float | None, o: ThermoOptions,
                   code_thermo: dict | None = None) -> dict:
    from ase import thermochemistry as tc

    out = {"model": o.model, "inputs": asdict(o), "computed": False, "reasons": [], "rows": [],
           "n_frequencies": len(freqs_cm1 or []), "n_imaginary_all": int(sum(1 for f in freqs_cm1 or [] if f < 0)),
           "imaginary_cm1_all": [float(f) for f in freqs_cm1 or [] if f < 0],
           "electronic_energy_ev": energy_ev, "code_values": (code_thermo or {}).get("items", []),
           "code_values_source": ({k: code_thermo.get(k) for k in ("code", "method", "source", "temperature_k")} if code_thermo else None),
           "source": f"ASE {__import__('ase').__version__} ase.thermochemistry"}
    miss = missing_inputs(o)
    if not freqs_cm1:
        miss.insert(0, L("振動数がありません (振動解析の出力が要ります)", "no frequencies (a vibrational analysis output is needed)"))
    if atoms is None and o.model in ("ideal_gas", "msrrho"):
        miss.append(L("構造がありません (慣性モーメントに要ります)", "no structure (needed for the moments of inertia)"))
    if miss:
        out["reasons"] = [L("入力が足りないので計算していません: ", "not computed; missing input: ") + ", ".join(miss)]
        return out
    e = _energies(freqs_cm1)
    idx = _order(e)
    n = len(atoms) if atoms is not None else None
    if o.model == "ideal_gas":
        nv = {"nonlinear": 3 * n - 6, "linear": 3 * n - 5, "monatomic": 0}[o.geometry]
        if nv < 0 or len(e) < nv:
            out["reasons"] = [L(f"振動の本数が足りません ({o.geometry} の {n} 原子なら {max(nv, 0)} 本、読めたのは {len(e)} 本)",
                                f"not enough modes ({max(nv, 0)} needed for {o.geometry} with {n} atoms, {len(e)} read)")]
            return out
        use = idx[len(idx) - nv:] if nv else []
        rule = L(f"振動数の 2 乗の大きい順に {nv} 本 ({o.geometry}、{n} 原子)", f"the {nv} modes with the largest squared frequency ({o.geometry}, {n} atoms)")
    else:
        if o.exclude_lowest > len(e):
            out["reasons"] = [L(f"除く本数 {o.exclude_lowest} が振動の本数 {len(e)} より多い", f"exclude_lowest {o.exclude_lowest} exceeds the number of modes {len(e)}")]
            return out
        use = idx[o.exclude_lowest:]
        rule = L(f"振動数の 2 乗の小さい方から {o.exclude_lowest} 本を除いた {len(use)} 本", f"{len(use)} modes after excluding the {o.exclude_lowest} with the smallest squared frequency")
    sel = [e[i] for i in sorted(use)]
    imag = [x for x in sel if x.imag != 0]
    zero = [x for x in sel if x.imag == 0 and abs(x.real) < 1e-12]
    out.update(mode_rule=rule, n_modes_used=len(sel), modes_used_cm1=[_to_cm1(x) for x in sel],
               n_imaginary_used=len(imag), imaginary_cm1_used=[_to_cm1(x) for x in imag],
               excluded_cm1=[_to_cm1(e[i]) for i in idx if i not in set(use)])
    if imag and o.imaginary is None:
        out["reasons"] = [L(f"使うモードに虚振動が {len(imag)} 本あります ({', '.join(f'{_to_cm1(x):.1f}' for x in imag)} cm⁻¹)。"
                            "除いて計算するか (ignore)、計算しないか (stop) を選んでください",
                            f"{len(imag)} imaginary mode(s) among the modes used ({', '.join(f'{_to_cm1(x):.1f}' for x in imag)} cm^-1); "
                            "choose ignore (drop them) or stop (do not compute)")]
        return out
    if imag and o.imaginary == "stop":
        out["reasons"] = [L(f"使うモードに虚振動が {len(imag)} 本あるので、指定どおり計算していません", f"{len(imag)} imaginary mode(s) among the modes used; not computed as requested")]
        return out
    if zero:
        out["reasons"] = [L(f"使うモードに 0 cm⁻¹ のものが {len(zero)} 本あります (エントロピーが発散します。除く本数を見直してください)",
                            f"{len(zero)} mode(s) at 0 cm^-1 among the modes used (the entropy diverges; revisit the number of excluded modes)")]
        return out
    ign = bool(imag) and o.imaginary == "ignore"
    mol = None
    if atoms is not None:
        mol = atoms.copy()
        mol.pbc = False
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        if o.model == "ideal_gas":
            th = tc.IdealGasThermo(sel, o.geometry, potentialenergy=0.0, atoms=mol, symmetrynumber=o.symmetry_number, spin=o.spin,
                                   vib_selection="all", ignore_imag_modes=ign)
            includes = L("並進・回転・振動・電子 (2S+1) の寄与。H = U + k_B T、G = H − T S。圧力は並進のエントロピーだけに効く",
                         "translational, rotational, vibrational and electronic (2S+1) contributions; H = U + k_B T, G = H - T S; pressure enters only the translational entropy")
        elif o.model == "harmonic":
            th = tc.HarmonicThermo(sel, potentialenergy=0.0, ignore_imag_modes=ign)
            includes = L("振動の寄与だけ (並進・回転を含まない)。F = U − T S", "vibrational contributions only (no translation or rotation); F = U - T S")
        elif o.model == "quasi_harmonic":
            th = tc.QuasiHarmonicThermo(sel, potentialenergy=0.0, ignore_imag_modes=ign, raise_to=o.qh_cutoff_cm1 * units.invcm)
            includes = L(f"振動の寄与だけ。{o.qh_cutoff_cm1:g} cm⁻¹ より低い振動数を {o.qh_cutoff_cm1:g} cm⁻¹ に上げて調和振動子で数える。F = U − T S",
                         f"vibrational contributions only; frequencies below {o.qh_cutoff_cm1:g} cm^-1 are raised to it before the harmonic treatment; F = U - T S")
        else:
            real = [x.real for x in sel if x.imag == 0 and x.real > 0]
            th = tc.MSRRHOThermo(real, atoms=mol, potentialenergy=0.0, tau=o.msrrho_tau_cm1, nu_scal=1.0, treat_int_energy=False)
            includes = L(f"振動の寄与だけ。低い振動のエントロピーを自由回転子と混ぜる (τ = {o.msrrho_tau_cm1:g} cm⁻¹、ASE の既定 treat_int_energy=False: 内部エネルギーは調和振動子のまま)。F = U − T S",
                         f"vibrational contributions only; low-mode entropies are interpolated with free rotors (tau = {o.msrrho_tau_cm1:g} cm^-1; ASE default treat_int_energy=False: harmonic internal energy); F = U - T S")
        zpe = float(th.get_ZPE_correction())
        rows = []
        for T in o.temperatures_k:
            U = float(th.get_internal_energy(T, verbose=False))
            if o.model == "ideal_gas":
                H = float(th.get_enthalpy(T, verbose=False))
                S = float(th.get_entropy(T, pressure=o.pressure_pa, verbose=False))
                row = {"T_K": T, "P_Pa": o.pressure_pa, "zpe_ev": zpe, "u_corr_ev": U, "h_corr_ev": H, "s_ev_per_k": S, "ts_ev": T * S, "g_corr_ev": H - T * S}
            else:
                S = float(th.get_entropy(T, verbose=False))
                row = {"T_K": T, "zpe_ev": zpe, "u_corr_ev": U, "s_ev_per_k": S, "ts_ev": T * S, "f_corr_ev": U - T * S}
            row["s_j_mol_k"] = S * EV_J_MOL
            for k in [k for k in row if k.endswith("_corr_ev")]:
                row[k.replace("_corr_ev", "_corr_kj_mol")] = row[k] * EV_KJ_MOL
                if energy_ev is not None:
                    row[k.replace("_corr_ev", "_total_ev")] = energy_ev + row[k]
            rows.append(row)
    out.update(computed=True, rows=rows, includes=includes, ase_class=type(th).__name__,
               note=L("*_corr_ev は電子エネルギーを含まない補正 (ZPE を含む)。*_total_ev = 電子エネルギー (出力の最終エネルギー) + 補正",
                      "*_corr_ev are corrections without the electronic energy (ZPE included); *_total_ev = electronic energy (final energy of the output) + correction"))
    return out


def rows_csv(table: dict, path) -> None:
    import csv
    rows = table.get("rows") or []
    if not rows:
        return
    cols = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)


def summary_lines(t: dict) -> list[str]:
    lab = {"ideal_gas": L("理想気体", "ideal gas"), "harmonic": L("調和振動子", "harmonic"), "quasi_harmonic": L("準調和", "quasi-harmonic"),
           "msrrho": L("準 RRHO (Grimme)", "quasi-RRHO (Grimme)")}.get(t.get("model"), str(t.get("model")))
    if not t.get("computed"):
        return [L(f"熱化学 (ASE、{lab}): ", f"thermochemistry (ASE, {lab}): ") + " ".join(t.get("reasons", []))]
    lines = [L(f"熱化学 (ASE の {t['ase_class']}): {t['includes']}。使ったモード {t['n_modes_used']} 本 ({t['mode_rule']})、"
               f"虚振動 {t['n_imaginary_all']} 本 (使うモードの中 {t['n_imaginary_used']} 本)",
               f"thermochemistry (ASE {t['ase_class']}): {t['includes']}. modes used {t['n_modes_used']} ({t['mode_rule']}), "
               f"imaginary {t['n_imaginary_all']} ({t['n_imaginary_used']} among the modes used)")]
    for r in t["rows"]:
        if "g_corr_ev" in r:
            v = f"ZPE {r['zpe_ev']:.4f} eV, H-E {r['h_corr_ev']:.4f} eV, S {r['s_j_mol_k']:.2f} J/(mol K), G-E {r['g_corr_ev']:.4f} eV"
            tot = f", G {r['g_total_ev']:.6f} eV" if "g_total_ev" in r else ""
            lines.append(f"  T = {r['T_K']:g} K, P = {r['P_Pa']:g} Pa: {v}{tot}")
        else:
            v = f"ZPE {r['zpe_ev']:.4f} eV, U-E {r['u_corr_ev']:.4f} eV, S {r['s_j_mol_k']:.2f} J/(mol K), F-E {r['f_corr_ev']:.4f} eV"
            tot = f", F {r['f_total_ev']:.6f} eV" if "f_total_ev" in r else ""
            lines.append(f"  T = {r['T_K']:g} K: {v}{tot}")
    return lines


__all__ = ["ThermoOptions", "MODELS", "GEOMETRIES", "compute_thermo", "missing_inputs", "summary_lines", "rows_csv"]
