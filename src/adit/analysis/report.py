"""Run an analysis: read the outputs, compute, and write figures and tables into <run_dir>/analysis/."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from adit import lang
from adit.lang import L
from ase.data import atomic_numbers, covalent_radii

from adit.analysis import compute
from adit.analysis.readers import RunData, load_run
from adit.analysis.trajectory import MEMORY_BUDGET_MB, Trajectory, check_budget
from adit.analysis.symmetry import DEFAULT_SYMPRECS
from adit.analysis.thermo import ThermoOptions
from adit.analysis import plotstyle

OUT_SUBDIR = "analysis"

FEW_FRAMES = 20
TEMP_DEVIATION = 0.30
RDF_CLIP_RATIO = 3.0
MIN_STATS_POINTS = 8


def figure_title(name: str) -> str:
    return {
        "energy": L("エネルギーの推移", "Energy"),
        "temperature": L("温度の推移", "Temperature"),
        "rdf": L("動径分布関数 (RDF)", "Radial distribution function (RDF)"),
        "msd": L("平均二乗変位 (MSD)", "Mean square displacement (MSD)"),
        "dos": L("状態密度 (DOS)", "Density of states (DOS)"),
        "spectrum": L("振動スペクトル", "Vibrational spectrum"),
        "bands": L("バンド図", "Band structure"),
        "scan_energy": L("値ごとのエネルギーの差", "Energy difference for each value"),
        "eos": L("エネルギーと体積 (状態方程式の当てはめ)", "Energy vs. volume (equation-of-state fit)"),
        "blocking": L("ブロック平均による平均値の誤差", "Error of the mean by block averaging"),
        "zdensity": L("z 方向の密度分布", "Density profile along z"),
        "coordination": L("配位数 n(r)", "Coordination number n(r)"),
        "pressure": L("圧力の推移", "Pressure"),
        "neb": L("NEB のエネルギーの曲線", "NEB energy profile"),
        "pdos": L("射影状態密度 (PDOS)", "Projected density of states (PDOS)"),
        "uvvis": L("UV-Vis スペクトル", "UV-Vis spectrum"),
        "phonon_bands": L("フォノン分散", "Phonon dispersion"),
        "phonon_dos": L("フォノンの状態密度", "Phonon density of states"),
        "compare_energy": L("組にした計算のエネルギー差", "Energy differences of the compared runs"),
        "elastic_stress_strain": L("応力と歪み (弾性定数の当てはめ)", "Stress vs. strain (elastic constant fits)"),
        "conformers": L("配座ごとの相対エネルギー", "Relative energy of each conformer"),
    }.get(name, name)


@dataclass
class AnalysisOptions:
    energy: bool = True
    temperature: bool = True
    bonds: bool = True
    rdf: bool = False
    rdf_pairs: list[tuple[str, str]] = field(default_factory=list)
    rdf_rmax: float = 8.0
    msd: bool = False
    msd_species: str | None = None
    msd_fit_fs: tuple[float, float] | None = None
    msd_remove_drift: bool = True
    msd_axes: str = "xyz"
    msd_error_blocks: int = 5
    dos: bool = False
    dos_sigma: float = 0.1
    vibrations: bool = True
    bands: bool = True
    bands_window_ev: float = 10.0
    skip_frames: int = 0
    stride: int = 1
    zdens: bool = False
    zdens_bin_ang: float = 0.2
    zdens_axis: str = "c"
    stats: bool = True
    export: bool = False
    export_unwrap: bool = False
    memory_budget_mb: float = MEMORY_BUDGET_MB
    thermo: ThermoOptions | None = None
    uvvis_broadening: tuple[str, float] | None = None
    pdos: bool = False
    spacegroup: bool = True
    symprecs: tuple[float, ...] = DEFAULT_SYMPRECS
    compare: bool = True
    collect: bool = False
    viscosity: bool = False
    plane_average: str | None = None
    cube_unit: str = ""
    work_function: bool = False
    out_dir: Path | str | None = None
    msd_per_atom: bool = False
    vanhove: int = 0
    vanhove_displacement: str = "mic"
    coordination_cutoff: float = 0.0
    centrosymmetry_neighbors: int = 0
    steinhardt_cutoff: float = 0.0
    cluster_cutoff: float = 0.0
    adf: str | None = None
    adf_cutoff: float = 0.0
    structure_factor: bool = False
    hbond: str = ""
    radius_of_gyration: bool = False
    density_grid: str = ""
    effective_mass: bool = True
    effective_mass_points: int = 5
    projected_bands: bool = True
    optical: bool = True
    bader: str = ""
    bader_valence: str = ""
    displacement_reference: int = 0
    strain_cutoff: float = 0.0
    voronoi: bool = False
    voronoi_face_threshold: float = 0.0
    sasa: str = ""
    pca: int = 0
    cluster: int = 0
    fes_temperature_k: float = 0.0
    fes_bins: int = 50
    fes_unit: str = "kJ/mol"
    plot_colors: str = ""
    plot_ticks: str = ""
    plot_grid: str = ""
    plot_spines: str = ""
    plot_line_width: float = 0.0
    plot_font_size: float = 0.0
    plot_dpi: int = 0
    figure_format: str = ""
    code: str = ""
    freq_scale: float = 0.0
    spectrum_measured: str = ""
    select: str = ""
    distances: list = field(default_factory=list)
    angles: list = field(default_factory=list)
    dihedrals: list = field(default_factory=list)
    rmsd_reference: int = 0
    rmsf: bool = False
    vacf: bool = False
    conductivity_charge: float = 0.0
    conductivity_temperature_k: float = 0.0
    vanhove_here: bool = False
    heavy_limit_seconds: float = 60.0
    xrd: str = ""
    xrd_range: tuple[float, float] = (5.0, 90.0)
    xrd_measured: Path | str | None = None


@dataclass
class AnalysisResult:
    code: str
    run_dir: str
    figures: dict[str, str] = field(default_factory=dict)
    tables: dict[str, object] = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def summary_text(self) -> str:
        lines = [L(f"コード: {self.code}   ディレクトリ: {self.run_dir}", f"code: {self.code}   directory: {self.run_dir}")]
        t = self.tables
        if "energy" in t:
            e = t["energy"]
            lines.append(L(f"エネルギー: {e['n']} 点、最終値 {e['last_ev']:.6f} eV (最小 {e['min_ev']:.6f}、最大 {e['max_ev']:.6f})",
                           f"energy: {e['n']} points; last {e['last_ev']:.6f} eV (min {e['min_ev']:.6f}, max {e['max_ev']:.6f})"))
        if "temperature" in t:
            x = t["temperature"]
            lines.append(L(f"温度: 平均 {x['mean_k']:.1f} K (標準偏差 {x['std_k']:.1f}、{x['n']} 点、先頭 {x['skipped']} 点を除く)",
                           f"temperature: mean {x['mean_k']:.1f} K (std {x['std_k']:.1f}, {x['n']} points, first {x['skipped']} skipped)"))
        for key, name_ja, name_en in (("temperature", "温度", "temperature"), ("energy", "エネルギー", "energy"), ("density", "密度", "density"),
                                      ("pressure", "圧力", "pressure"), ("conserved", "保存量", "conserved quantity"),
                                      ("density_log", "密度 (ログファイルの値)", "density (from the log file)")):
            s = t.get("timeseries_stats", {}).get(key)
            if s:
                tau = ""
                if s["tau_int_samples"] is not None:
                    tfs = f" = {s['tau_int_fs']:.3g} fs" if s.get("tau_int_fs") else ""
                    if s["sem_acf"] is not None:
                        tau = L(f"、積分自己相関時間 {s['tau_int_samples']:.3g} 点{tfs}、それからの平均値の標準誤差 {s['sem_acf']:.3g} {s['unit']}",
                                f", integrated autocorrelation time {s['tau_int_samples']:.3g} samples{tfs}, giving a standard error of the mean of {s['sem_acf']:.3g} {s['unit']}")
                    else:
                        tau = L(f"、積分自己相関時間 {s['tau_int_samples']:.3g} 点{tfs} (0 以下なので、これからの誤差は出しません)",
                                f", integrated autocorrelation time {s['tau_int_samples']:.3g} samples{tfs} (not positive, so no error is derived from it)")
                    if not s.get("acf_window_reached"):
                        tau += L(" (自動の窓の条件を満たす前に系列が尽きた)", " (the series ended before the automatic window condition was met)")
                b = ", ".join(f"{x['block_size']}:{x['sem']:.2g}" for x in s["blocks"][:8]) + (" …" if len(s["blocks"]) > 8 else "")
                lines.append(L(f"{name_ja}の統計: 平均 {s['mean']:.6g} {s['unit']} (標準偏差 {s['std']:.3g}、{s['n']} 点){tau}。"
                               f"ブロック平均の誤差 (束ねた点の数:誤差) {b}",
                               f"{name_en} statistics: mean {s['mean']:.6g} {s['unit']} (std {s['std']:.3g}, {s['n']} points){tau}. "
                               f"block-averaging error (block size:error) {b}"))
        if "trajectory" in t:
            x = t["trajectory"]
            dt = L(f"、1 フレーム {x['dt_used_fs']:g} fs", f", {x['dt_used_fs']:g} fs per frame") if x.get("dt_used_fs") else ""
            lines.append(L(f"軌跡 ({x['source']}): {x['n_frames_total']} フレームのうち {x['n_frames_used']} フレームを使用 "
                           f"(先頭 {x['skip']} フレームを除き、{x['stride']} フレームに 1 回){dt}",
                           f"trajectory ({x['source']}): {x['n_frames_used']} of {x['n_frames_total']} frames used "
                           f"(first {x['skip']} frames skipped, then every {x['stride']}th frame){dt}"))
        if "bonds" in t:
            if t["bonds"]:
                bonds = ", ".join(f"{n} {d:.3f} Å" for n, d in t["bonds"][:12]) + (" …" if len(t["bonds"]) > 12 else "")
            else:
                bonds = L("該当なし (共有結合半径の和の 1.2 倍以内にある原子の組がありません)", "none (no atom pair within 1.2 x the sum of covalent radii)")
            lines.append(L("結合長 (最終構造): ", "bond lengths (final structure): ") + bonds)
        for q in t.get("charges", []):
            names = q.get("atoms") or [str(i + 1) for i in range(len(q["values"]))]
            vals = ", ".join(f"{n} {v:+.3f}" for n, v in list(zip(names, q["values"]))[:12]) + (" …" if len(q["values"]) > 12 else "")
            f = L(f"。表は {q['file']}", f"; table in {q['file']}") if q.get("file") else ""
            lines.append(L(f"原子の電荷 ({q['definition']}、{q['source']}、単位 e): {vals} (合計 {q['sum']:+.4f})",
                           f"atomic charges ({q['definition']}, {q['source']}, in e): {vals} (sum {q['sum']:+.4f})") + f)
        if "thermochemistry" in t:
            th = t["thermochemistry"]
            temp = L(f"、{th['temperature_k']:g} K", f", {th['temperature_k']:g} K") if th.get("temperature_k") is not None else ""
            vals = ", ".join(f"{it['label']} {it['value_eh']:.6f} Eh" for it in th["items"] if it["key"] != "total_free_energy_box")
            lines.append(L(f"熱化学 ({th['code']}、{th['source']}{temp}): ", f"thermochemistry ({th['code']}, {th['source']}{temp}): ") + vals)
        if "electronic" in t:
            el = t["electronic"]
            parts = []
            if el.get("homo_lumo_gap_ev") is not None:
                parts.append(L(f"HOMO-LUMO ギャップ {el['homo_lumo_gap_ev']:.4f} eV ({el['gap_source']})", f"HOMO-LUMO gap {el['homo_lumo_gap_ev']:.4f} eV ({el['gap_source']})"))
            if el.get("dipole_norm_debye") is not None:
                x, y, z = el["dipole_debye"]
                parts.append(L(f"双極子モーメント {el['dipole_norm_debye']:.4f} D (x, y, z = {x:.4f}, {y:.4f}, {z:.4f})",
                               f"dipole moment {el['dipole_norm_debye']:.4f} D (x, y, z = {x:.4f}, {y:.4f}, {z:.4f})"))
            if parts:
                lines.append("、".join(parts) if lang.LANGUAGE != "en" else "; ".join(parts))
        if "rdf" in t:
            r = t["rdf"]
            lines.append(L(f"RDF (動径分布関数): {', '.join(r['pairs'])} ({r['n_frames']} フレーム、r ≤ {r['rmax']:.2f} Å)。配位数 n(r) は rdf.json の n と n_reverse",
                           f"RDF (radial distribution function): {', '.join(r['pairs'])} ({r['n_frames']} frames, r ≤ {r['rmax']:.2f} Å); coordination n(r) is in rdf.json as n and n_reverse"))
            nz = r.get("normalization") or {}
            if nz.get("kind"):
                v = f"{nz['mean_volume_A3']:.1f}" if nz.get("mean_volume_A3") is not None else "?"
                if nz["kind"] == "cell_volume":
                    lines.append(L(f"  g(r) の規格化: セルの体積 V = {v} Å³ (周期系)。rmax を変えても g(r) の値は変わりません",
                                   f"  g(r) normalization: the cell volume V = {v} Å³ (periodic); the value of g(r) does not change with rmax"))
                else:
                    lines.append(L(f"  g(r) の規格化: 半径 rmax = {nz['rmax_A']:.2f} Å の球 V = {v} Å³ (非周期)。rmax を変えると g(r) の絶対値が変わるので、"
                                   "別の rmax や別の計算とは比べられません。rmax に依らない量は rdf.json の density [Å⁻³] と n(r) です",
                                   f"  g(r) normalization: a sphere of rmax = {nz['rmax_A']:.2f} Å, V = {v} Å³ (non-periodic); the absolute value of g(r) "
                                   "changes with rmax, so it cannot be compared with another rmax or another run. The rmax-independent quantities are "
                                   "density [Å⁻³] and n(r) in rdf.json"))
        if "msd" in t:
            m = t["msd"]
            last_msd = f"{m['last_A2']:.3g}"
            diff = ""
            if m.get("D_cm2_s") is not None:
                diff = L(f"、拡散係数 {m['D_cm2_s']:.3e} cm²/s", f", diffusion coefficient {m['D_cm2_s']:.3e} cm^2/s")
            lines.append(L(f"MSD (平均二乗変位。{m['species'] or '全原子'}、{m['n_frames']} フレーム): 最終 {last_msd} Å²",
                           f"MSD (mean squared displacement; {m['species'] or 'all atoms'}, {m['n_frames']} frames): last {last_msd} Å²") + diff)
            per = [f"{el} {v['D_cm2_s']:.3e}"
                   for el, v in m.get("by_element", {}).items() if v.get("D_cm2_s") is not None]
            if len(m.get("by_element", {})) > 1 and per:
                label = (L("  参考・全元素の拡散係数 [cm²/s]: ", "  for reference, diffusion coefficients for all elements [cm^2/s]: ")
                         if m.get("species") else L("  元素ごとの拡散係数 [cm²/s]: ", "  diffusion coefficient by element [cm^2/s]: "))
                lines.append(label + ", ".join(per))
            if m.get("formula"):
                lines.append(L(f"  使った式: {m['formula']} ({m['dimension']} 次元、成分 {m['axes']}、複数の時間原点で平均した MSD)",
                               f"  formula used: {m['formula']} ({m['dimension']}D, components {m['axes']}, MSD averaged over time origins)"))
            d = m.get("drift") or {}
            if d.get("removed"):
                w = L("質量重心", "center of mass") if d.get("mass_weighted") else L("幾何中心 (元素が読めないため)", "geometric center (elements unreadable)")
                lines.append(L(f"  重心の移動: 除去しました ({w}、最初と最後で {d['displacement_A']:.3f} Å、最大 {d['max_displacement_A']:.3f} Å 動いた)",
                               f"  center-of-mass drift: removed ({w}; it moved {d['displacement_A']:.3f} Å between first and last frame, at most {d['max_displacement_A']:.3f} Å)"))
            elif d:
                lines.append(L("  重心の移動: 除去していません (--msd-keep-drift)", "  center-of-mass drift: not removed (--msd-keep-drift)"))
            if m.get("fit_range_fs"):
                a, b = m["fit_range_fs"]
                how = L("利用者の指定", "set by the user") if m.get("fit_range_user") else L(
                    f"既定: 最大の遅れ時間 {m['lag_fs'][-1]:g} fs の {m['fit_fraction'][0]:.0%}〜{m['fit_fraction'][1]:.0%}",
                    f"default: {m['fit_fraction'][0]:.0%}-{m['fit_fraction'][1]:.0%} of the maximum lag {m['lag_fs'][-1]:g} fs")
                sl = L(f"、この範囲の log MSD 対 log t の傾き {m['loglog_slope']:.2f} (拡散なら 1 に近い)",
                       f", slope of log MSD vs log t over this range {m['loglog_slope']:.2f} (close to 1 for diffusion)") if m.get("loglog_slope") is not None else ""
                lines.append(L(f"  当てはめ範囲: {a:g}〜{b:g} fs ({how}){sl}", f"  fit range: {a:g}-{b:g} fs ({how}){sl}"))
            e = m.get("D_error") or {}
            if e.get("d_err_cm2_s") is not None:
                counts = e.get("block_frame_counts") or [e["block_frames"]] * e["n_blocks"]
                sizes = ", ".join(str(x) for x in counts)
                a, b = e.get("fit_range_fs") or m["fit_range_fs"]
                lines.append(L(f"  参考・ブロック D 平均の標準誤差: {e['d_err_cm2_s']:.2e} cm²/s。"
                               f"全 {sum(counts)} フレームを {e['n_blocks']} ブロック ({sizes} フレーム) に分け、各ブロックの {a:g}〜{b:g} fs を当てはめ、"
                               f"ブロックごとの D の標本標準偏差 ÷ √{e['n_blocks']} としました。全軌跡から出した D 自体の厳密な誤差ではありません",
                               f"  for reference, standard error of the mean block D: {e['d_err_cm2_s']:.2e} cm^2/s. "
                               f"All {sum(counts)} frames were split into {e['n_blocks']} blocks ({sizes} frames); each block was fitted over {a:g}-{b:g} fs. "
                               f"This is the sample standard deviation of block D values / sqrt({e['n_blocks']}), not a rigorous error on D from the full trajectory."))
            elif e.get("reason"):
                lines.append(L(f"  ブロック誤差: 出せません。{e['reason']}", f"  block-based error: not available. {e['reason']}"))
        if "xtb_md" in t:
            x = t["xtb_md"]
            if x["found"]:
                got = ", ".join(f"{k} = {v}" for k, v in x["found"].items())
                lines.append(L(f"  xtb の MD の設定 (xtb.inp の $md に書かれた値): {got}",
                               f"  xtb MD settings (values written in the $md block of xtb.inp): {got}"))
            if x["missing"]:
                miss = "、".join(f"{k} (既定 {x['defaults'][k]})" for k in x["missing"])
                en = ", ".join(f"{k} (default {x['defaults'][k]})" for k in x["missing"])
                lines.append(L(f"  xtb.inp に書かれていない MD の設定: {miss}。xtb の既定で走った可能性があります "
                               "(hmass は水素の質量 [u]、shake は結合の拘束 0/1/2、sccacc は MD 中の SCC の精度)",
                               f"  MD settings not written in xtb.inp: {en}; the run may have used the xtb defaults "
                               "(hmass is the hydrogen mass [u], shake constrains bonds 0/1/2, sccacc is the SCC accuracy during MD)"))
        if "zdensity" in t:
            z = t["zdensity"]
            lines.append(L(f"z 方向の密度分布: 長さ {z['axis_length_A']:.3f} Å を {z['bin_A']:.3f} Å 刻み、{z['n_frames']} フレーム (zdensity.json)",
                           f"density profile along z: {z['axis_length_A']:.3f} Å in {z['bin_A']:.3f} Å bins, {z['n_frames']} frames (zdensity.json)"))
        if "frequencies" in t:
            f = [x for x in t["frequencies"] if x > 50]
            neg = [x for x in t["frequencies"] if x < -50]
            top = ", ".join(f"{x:.0f}" for x in f[-8:]) + " cm^-1"
            lines.append(L(f"振動数 (振動の速さ): {len(f)} 本 (50 cm⁻¹ 超)。高い順に 8 本まで: ",
                               f"frequencies: {len(f)} above 50 cm^-1; up to the eight highest: ") + top
                         + (L(f"。虚振動 {len(neg)} 本", f". {len(neg)} imaginary") if neg else ""))
        if "bands" in t:
            b = t["bands"]
            if b.get("gap_ev") is not None:
                gap = (L(f"、最高被占準位を基準にした最小の間隔 {b['gap_ev']:.3f} eV",
                         f", smallest gap relative to the highest occupied level: {b['gap_ev']:.3f} eV")
                       if b.get("reference_level") == "highest_occupied" else
                       L(f"、フェルミ準位をまたぐ最小の間隔 {b['gap_ev']:.3f} eV",
                         f", smallest gap across the Fermi level {b['gap_ev']:.3f} eV"))
            else:
                gap = ""
            lines.append(L(f"バンド: {b['n_kpoints']} k 点 × {b['n_bands']} 本 (経路 {' '.join(b['labels'])})",
                           f"bands: {b['n_kpoints']} k-points × {b['n_bands']} bands (path {' '.join(b['labels'])})") + gap)
        if "dos" in t:
            ef = L(f"、フェルミ準位 {t['dos']['fermi_ev']:.3f} eV", f", Fermi level {t['dos']['fermi_ev']:.3f} eV") if t['dos'].get('fermi_ev') is not None else ""
            lines.append(L(f"DOS (状態密度): {t['dos']['n_eigen']} 個の固有値", f"DOS (density of states): {t['dos']['n_eigen']} eigenvalues") + ef)
        from adit.analysis import neb as _neb, pdos as _pdos, phonons as _ph, symmetry as _sym, thermo as _th, uvvis as _uv
        from adit.analysis import collections as _coll
        for key, fn in (("thermo_ase", _th.summary_lines), ("pdos", _pdos.summary_lines), ("neb", _neb.summary_lines),
                        ("uvvis", _uv.summary_lines), ("spacegroup", _sym.summary_lines), ("phonopy", _ph.summary_lines),
                        ("phonon_set", _coll.phonon_set_lines), ("elastic", _coll.elastic_lines), ("conformers", _coll.conformer_lines)):
            if key in t:
                lines += fn(t[key])
        if "mlip" in t:
            m = t["mlip"]
            ver = ", ".join(f"{k} {v}" for k, v in (m.get("versions") or {}).items())
            st = L(f"、圧力 {m['pressure_gpa']:+.4f} GPa", f", pressure {m['pressure_gpa']:+.4f} GPa") if m.get("pressure_gpa") is not None else ""
            fm = L(f"、最後の構造の力の最大値 {m['fmax_ev_per_ang']:.4g} eV/Å", f", largest force on the final structure {m['fmax_ev_per_ang']:.4g} eV/Å") if m.get("fmax_ev_per_ang") is not None else ""
            lines.append(L(f"機械学習ポテンシャル ({m.get('model_family')}、モデル {m.get('model')}、{m.get('task')}): {ver}",
                           f"machine-learning potential ({m.get('model_family')}, model {m.get('model')}, {m.get('task')}): {ver}") + fm + st)
            if m.get("converged") is not None:
                lines.append(L(f"  ASE の最適化: {'収束の条件を満たしました' if m['converged'] else '収束の条件を満たしていません'} ({m.get('steps')} ステップ)",
                               f"  ASE optimization: {'the convergence criterion was met' if m['converged'] else 'the convergence criterion was not met'} ({m.get('steps')} steps)"))
            if m.get("stress_note"):
                lines.append("  " + str(m["stress_note"]))
        if "compare" in t:
            lines += str(t["compare"].get("summary", "")).splitlines()
        if "export" in t:
            x = t["export"]
            lines.append(L(f"書き出し: {x['dir']} ({x['n_frames']} フレーム。trajectory.extxyz / .xyz / .pdb、view.vmd、ovito_pipeline.py、export_README.txt)",
                           f"export: {x['dir']} ({x['n_frames']} frames; trajectory.extxyz / .xyz / .pdb, view.vmd, ovito_pipeline.py, export_README.txt)"))
        lines += [L("注: ", "note: ") + n for n in self.notes]
        return "\n".join(lines)


def _jsonable(x):
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        return x.tolist()
    if isinstance(x, np.generic):
        return x.item()
    return x


def _check_options(opts: AnalysisOptions) -> None:
    import math

    if opts.stride < 1:
        raise ValueError(L(f"間引きの間隔 (stride) は 1 以上です: {opts.stride}", f"stride must be 1 or more: {opts.stride}"))
    if opts.skip_frames < 0:
        raise ValueError(L(f"捨てる先頭のフレーム数 (skip) は 0 以上です: {opts.skip_frames}", f"the number of frames to skip must be 0 or more: {opts.skip_frames}"))
    for value, ja, en in ((opts.rdf_rmax, "動径分布関数の距離の上限 [Å]", "the RDF upper distance [Å]"),
                          (opts.dos_sigma, "状態密度を広げる幅 [eV]", "the DOS broadening width [eV]"),
                          (opts.zdens_bin_ang, "z 方向の密度の区間の幅 [Å]", "the z-density bin width [Å]"),
                          (opts.memory_budget_mb, "MSD のメモリの上限 [MB]", "the MSD memory limit [MB]"),
                          (opts.bands_window_ev, "バンド図の縦軸の幅 [eV]", "the band-plot energy window [eV]")):
        if not math.isfinite(value) or value <= 0:
            raise ValueError(L(f"{ja} は 0 より大きい数にしてください: {value}", f"{en} must be a number greater than 0: {value}"))


FIGURE_FORMATS = ("svg", "pdf", "eps")


def _extra_figure_formats(opts: "AnalysisOptions | None") -> list[str]:
    if opts is None or not opts.figure_format:
        return []
    out = []
    for item in str(opts.figure_format).replace(" ", "").split(","):
        if not item:
            continue
        if item.lower() not in FIGURE_FORMATS:
            raise ValueError(L(f"図の形式は {', '.join(FIGURE_FORMATS)} のどれかです: {item}",
                               f"the figure format must be one of {', '.join(FIGURE_FORMATS)}: {item}"))
        out.append(item.lower())
    return out


def _apply_plot_style(opts: "AnalysisOptions | None") -> None:
    if opts is None:
        return
    from adit.analysis.plotstyle import from_text

    from adit.analysis.plotstyle import set_grid_choice

    style = from_text(opts.plot_colors, opts.plot_ticks, opts.plot_grid, opts.plot_spines,
                      opts.plot_line_width, opts.plot_font_size, opts.plot_dpi)
    style.apply()
    set_grid_choice(style.grid)


def use_cjk_font() -> str:
    import logging
    import sys

    from matplotlib import font_manager, rcParams

    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

    roots = [Path(sys.prefix) / "fonts"]
    bundled = getattr(sys, "_MEIPASS", "")
    if bundled:
        roots.insert(0, Path(bundled) / "fonts")
    for path in [root / name for root in roots for name in ("NotoSansCJKjp-VF.ttf", "NotoSansCJKjp-Regular.otf")]:
        if path.is_file():
            try:
                font_manager.fontManager.addfont(str(path))
                family = font_manager.FontProperties(fname=str(path)).get_name()
            except Exception:
                return ""
            rcParams["font.family"] = [family, "DejaVu Sans"]
            return family
    for family in ("Noto Sans CJK JP", "IPAGothic", "TakaoGothic", "Yu Gothic", "Hiragino Sans"):
        if any(f.name == family for f in font_manager.fontManager.ttflist):
            rcParams["font.family"] = [family, "DejaVu Sans"]
            return family
    return ""


def run_analysis(run_dir: Path | str, opts: AnalysisOptions | None = None) -> AnalysisResult:
    import matplotlib
    matplotlib.use("Agg")
    _apply_plot_style(opts)
    import matplotlib.pyplot as plt

    font = use_cjk_font()

    opts = opts or AnalysisOptions()
    _check_options(opts)
    run_dir = Path(run_dir)
    try:
        data: RunData = load_run(run_dir, opts.code or None)
    except ValueError:
        from adit.analysis.collections import detect_collection
        from adit.analysis.neb import find_neb
        kind = detect_collection(run_dir)
        found_neb = find_neb(run_dir)
        if found_neb is None and kind is None:
            raise
        native_kind = found_neb.get("kind") if found_neb is not None else None
        data = RunData(native_kind if native_kind not in (None, "images") else (kind or "neb"), run_dir)
    if data.code == "vasp" and not (run_dir / "vasprun.xml").exists():
        from adit.analysis.neb import find_neb
        found = find_neb(run_dir)
        if found is not None and found.get("kind") == "vasp":
            prefixes = ("vasprun.xml を読めません:", "cannot read vasprun.xml:")
            data.notes = [note for note in data.notes if not note.startswith(prefixes)]
    out = Path(opts.out_dir).expanduser() if opts.out_dir else run_dir / OUT_SUBDIR
    out.mkdir(parents=True, exist_ok=True)
    res = AnalysisResult(code=data.code, run_dir=str(run_dir), notes=list(data.notes))
    msd_memory = None
    if opts.msd:
        if isinstance(data.frames, Trajectory):
            nat = data.frames.first_natoms()
        else:
            nat = len(data.frames[0]) if len(data.frames) else 0
        if nat:
            sel = data.frames[opts.skip_frames::opts.stride] if (opts.skip_frames or opts.stride > 1) else data.frames
            msd_memory = check_budget(sel, nat, opts.memory_budget_mb)
    n_all = len(data.frames)
    skip = opts.skip_frames if n_all > opts.skip_frames else 0
    frames = data.frames[skip::opts.stride] if (skip or opts.stride > 1) else data.frames
    spec = _load_spec(run_dir)
    task_type = spec.task.type if spec is not None else None
    is_md = task_type == "molecular_dynamics" or bool(data.temperatures_k)
    if data.code == "xtb" and is_md:
        xtb_md = _xtb_md_settings(run_dir)
        if xtb_md is not None:
            res.tables["xtb_md"] = xtb_md
    fermi_used = False
    dt_used = data.frame_dt_fs * opts.stride if data.frame_dt_fs else None

    def save(fig, name: str) -> None:
        p = out / f"{name}.png"
        fig.savefig(p, dpi=110, bbox_inches="tight")
        for suffix in _extra_figure_formats(opts):
            fig.savefig(out / f"{name}.{suffix}", bbox_inches="tight")
        plt.close(fig)
        res.figures[name] = str(p)

    if opts.energy and data.energies_ev:
        e = np.array(data.energies_ev)
        if len(e) >= 2:
            x = np.array(data.times_fs[: len(e)]) if data.times_fs and len(data.times_fs) >= len(e) else np.arange(len(e))
            fig, ax = plt.subplots(figsize=(6, 3.2))
            ax.plot(x, e, lw=1.2); ax.set_xlabel("time [fs]" if data.times_fs else "step"); ax.set_ylabel("energy [eV]"); plotstyle.grid(ax)
            ax.ticklabel_format(axis="y", useOffset=False, style="plain")
            save(fig, "energy")
        res.tables["energy"] = {"n": int(len(e)), "last_ev": float(e[-1]), "min_ev": float(e.min()), "max_ev": float(e.max())}
    if opts.temperature and data.temperatures_k:
        T = np.array(data.temperatures_k)
        if len(T) >= 2:
            x = np.array(data.times_fs[: len(T)]) if data.times_fs and len(data.times_fs) >= len(T) else np.arange(len(T))
            fig, ax = plt.subplots(figsize=(6, 3.2))
            ax.plot(x, T, lw=1.0); ax.set_xlabel("time [fs]" if data.times_fs else "step"); ax.set_ylabel("T [K]"); plotstyle.grid(ax)
            save(fig, "temperature")
        Ts = T[opts.skip_frames:] if len(T) > opts.skip_frames else T
        res.tables["temperature"] = {"n": int(len(Ts)), "mean_k": float(Ts.mean()), "std_k": float(Ts.std()), "skipped": int(len(T) - len(Ts))}
        target = _target_temperature(spec)
        if target and abs(Ts.mean() - target) > TEMP_DEVIATION * target:
            res.tables["temperature"]["target_k"] = target
            res.notes.append(L(f"MD の温度: 目標 {target:.0f} K、平均 {Ts.mean():.0f} K", f"MD temperature: target {target:.0f} K, mean {Ts.mean():.0f} K"))
    if "pressure" in data.series and len(data.series["pressure"]["values"]) >= 1:
        P = np.array(data.series["pressure"]["values"], dtype=float)
        if len(P) >= 2:
            x = np.array(data.times_fs[: len(P)]) if data.times_fs and len(data.times_fs) >= len(P) else np.arange(len(P))
            fig, ax = plt.subplots(figsize=(6, 3.2))
            ax.plot(x, P, lw=1.0); ax.set_xlabel("time [fs]" if data.times_fs else "step"); ax.set_ylabel("P [bar]"); plotstyle.grid(ax)
            save(fig, "pressure")
        Ps = P[opts.skip_frames:] if len(P) > opts.skip_frames else P
        res.tables["pressure"] = {"n": int(len(Ps)), "mean_bar": float(np.nanmean(Ps)), "std_bar": float(np.nanstd(Ps)),
                                  "skipped": int(len(P) - len(Ps)), "source": data.series["pressure"]["source"]}
    if opts.bonds and data.final is not None:
        res.tables["bonds"] = compute.bonds(data.final)
    _add_properties(res, data, out)
    res.tables.update(data.extra_tables)

    first = frames[0] if len(frames) else None
    periodic = first is not None and bool(any(first.pbc)) and first.cell.rank == 3
    rdf_acc = unwrap = zacc = exporter = None
    densities: list[float] = []
    rmax = opts.rdf_rmax
    if first is not None and (periodic or opts.rdf):
        cell = np.asarray(first.cell, dtype=float)
        if periodic:
            vol = abs(np.linalg.det(cell))
            widths = [vol / np.linalg.norm(np.cross(cell[(i + 1) % 3], cell[(i + 2) % 3])) for i in range(3)]
            rmax = min(rmax, 0.5 * min(widths))
    if opts.rdf and first is not None:
        elems = sorted(set(first.get_chemical_symbols()))
        pairs = opts.rdf_pairs or [(a, b) for i, a in enumerate(elems) for b in elems[i:]]
        rdf_acc = compute.RDFAccumulator(pairs, rmax)
    if opts.msd and first is not None:
        if msd_memory is not None:
            res.tables["msd_memory"] = msd_memory
        unwrap = compute.UnwrapAccumulator(len(frames))
    if opts.zdens and first is not None:
        if periodic:
            zacc = compute.ZDensityAccumulator(opts.zdens_bin_ang, axis={'a': 0, 'b': 1, 'c': 2}.get(opts.zdens_axis, 2))
        else:
            res.notes.append(L("z 方向の密度分布は、セルのある (周期系の) 軌跡だけで出します", "the density profile along z is only produced for trajectories with a cell (periodic systems)"))
    want_density = opts.stats and is_md and periodic and len(frames) >= 2
    if opts.plane_average or opts.work_function:
        _add_volumetric(res, data, run_dir, out, opts, save)
    if opts.export and first is not None:
        from adit.analysis.export import EXPORT_SUBDIR, TrajectoryExporter
        exporter = TrajectoryExporter(out / EXPORT_SUBDIR, unwrap_molecules=opts.export_unwrap)
    if any(x is not None for x in (rdf_acc, unwrap, zacc, exporter)) or want_density:
        for fr in frames:
            if rdf_acc is not None:
                rdf_acc.add(fr)
            if unwrap is not None:
                unwrap.add(fr)
            if zacc is not None:
                zacc.add(fr)
            if exporter is not None:
                exporter.add(fr)
            if want_density:
                densities.append(compute.density_g_cm3(fr))
    if len(data.frames) > 1:
        res.tables["trajectory"] = {"source": data.frame_source, "n_frames_total": int(n_all), "skip": int(skip), "stride": int(opts.stride),
                                    "n_frames_used": int(len(frames)), "dt_frame_fs": data.frame_dt_fs, "dt_used_fs": dt_used}

    if rdf_acc is not None:
        fig, ax = plt.subplots(figsize=(6, 3.4))
        table = {}
        g_all, g_outer = 0.0, 0.0
        for a, b in rdf_acc.pairs:
            rr = rdf_acc.result((a, b))
            r, g = rr["r"], rr["g"]
            ax.plot(r, g, lw=1.2, label=f"{a}-{b}")
            table[f"{a}-{b}"] = {"r": r.tolist(), "g": g.tolist(), "r_upper": rr["r_edges_upper"].tolist(), "n": rr["n"].tolist(),
                                 "n_reverse": rr["n_reverse"].tolist(), "density": rr["density"].tolist(),
                                 "density_reverse": rr["density_reverse"].tolist()}
            if len(g):
                g_all = max(g_all, float(np.max(g)))
                outer = g[r > 1.2 * (covalent_radii[atomic_numbers[a]] + covalent_radii[atomic_numbers[b]])]
                g_outer = max(g_outer, float(np.max(outer)) if outer.size else 0.0)
        ax.set_xlabel("r [Å]"); ax.set_ylabel("g(r)"); ax.legend(); plotstyle.grid(ax)
        clipped = None
        if g_outer > 0 and g_all > RDF_CLIP_RATIO * g_outer:
            clipped = 1.2 * g_outer
            ax.set_ylim(0, clipped)
            ax.text(0.99, 0.97, L(f"縦軸は {clipped:.1f} で切っています (結合距離の山の最大 {g_all:.0f}。値は rdf.json に全部)",
                                  f"y axis cut at {clipped:.1f} (bonded peak up to {g_all:.0f}; full data in rdf.json)"),
                    transform=ax.transAxes, ha="right", va="top", fontsize=8, color="0.3")
        save(fig, "rdf")
        fig, ax = plt.subplots(figsize=(6, 3.4))
        for a, b in rdf_acc.pairs:
            v = table[f"{a}-{b}"]
            ax.plot(v["r_upper"], v["n"], lw=1.2, label=f"{b} around {a}")
            if a != b:
                ax.plot(v["r_upper"], v["n_reverse"], lw=1.0, ls="--", label=f"{a} around {b}")
        ax.set_xlabel("r [Å]"); ax.set_ylabel("n(r)"); ax.legend(fontsize=8); plotstyle.grid(ax)
        save(fig, "coordination")
        norm = rdf_acc.normalization()
        periodic_norm = norm["kind"] == "cell_volume"
        norm["definition"] = L(
            "g(r) = (距離 r の殻にある組の数) × V / (組の総数 × 殻の体積)。周期系は V = セルの体積 (rmax を変えても g(r) は変わらない)。"
            "非周期 (分子・クラスター) は V = 半径 rmax の球の体積で、基準にする密度が rmax で決まるので rmax を変えると g(r) の絶対値が変わる。"
            "有限の系には遠方の一様な密度が無いため、この規格化を rmax に依らない形にはできない。"
            "rmax に依らない量は rdf.json の density [Å⁻³] (元素 a の原子 1 個から距離 r の位置にある元素 b の数密度。density_reverse は入れ替え) と n(r)。"
            "別の rmax・別の計算と比べるときは g(r) ではなく density か n(r) を使う",
            "g(r) = (pairs in the shell at r) x V / (total number of pairs x shell volume). For periodic systems V is the cell volume, so g(r) does "
            "not change with rmax. For non-periodic systems (molecules, clusters) V is the volume of a sphere of radius rmax, so the reference density "
            "is set by rmax and the absolute value of g(r) changes when rmax changes. A finite system has no uniform density far away, so this "
            "normalization cannot be made independent of rmax. The rmax-independent quantities are density [Å⁻³] in rdf.json (number density of "
            "element b at distance r from one atom of element a; density_reverse swaps them) and n(r). Compare density or n(r), not g(r), across "
            "different rmax or different runs")
        res.tables["rdf"] = {"pairs": list(table), "n_frames": int(rdf_acc.n_frames), "rmax": float(rmax), "rmax_requested": float(opts.rdf_rmax),
                             "ylim_clipped": clipped, "g_max": g_all, "file": str(out / "rdf.json"), "normalization": norm,
                             "coordination_definition": L("n: 1 つ目の元素の原子 1 個のまわりで r 以内にある 2 つ目の元素の原子の数 (フレーム平均)。n_reverse は元素を入れ替えたもの。r は各区間の上端 r_upper",
                                                          "n: number of atoms of the second element within r of one atom of the first element (frame average); n_reverse swaps the elements; r is the upper bin edge r_upper")}
        table["_normalization"] = norm
        if not periodic_norm:
            v = norm["mean_volume_A3"]
            res.notes.append(L(
                f"g(r) の規格化: 非周期系なので半径 rmax = {rmax:.2f} Å の球 (V = {v:.1f} Å³) を体積に使っています。rmax を変えると g(r) の絶対値が変わるので、"
                "別の rmax や別の計算の g(r) とは比べられません。rmax に依らない量は rdf.json の density [Å⁻³] と n(r) です (周期系はセルの体積で規格化するので rmax に依りません)",
                f"g(r) normalization: this system is not periodic, so the volume of a sphere of rmax = {rmax:.2f} Å (V = {v:.1f} Å³) is used. The absolute "
                "value of g(r) changes with rmax, so it cannot be compared with a g(r) from another rmax or another run; the rmax-independent quantities "
                "are density [Å⁻³] and n(r) in rdf.json (periodic systems are normalized by the cell volume and do not depend on rmax)"))
        if rmax < opts.rdf_rmax:
            res.notes.append(L(f"RDF はセルの幅の半分 ({rmax:.2f} Å) までしか数えられません (周期境界の最小像のため。指定は {opts.rdf_rmax:.1f} Å)",
                               f"the RDF can only be counted up to half the cell width ({rmax:.2f} Å) because of the minimum-image convention (requested {opts.rdf_rmax:.1f} Å)"))
        if clipped is not None:
            res.notes.append(L(f"RDF の図は縦軸を {clipped:.1f} で切っています。結合距離 (共有結合半径の和の 1.2 倍) より内側の山は最大 {g_all:.0f} です",
                               f"the RDF plot is cut at {clipped:.1f} on the y axis; the peaks inside the bonded distance (1.2 x the sum of covalent radii) reach {g_all:.0f}"))
        (out / "rdf.json").write_text(json.dumps(table), encoding="utf-8")
    if unwrap is not None and len(frames) < 2:
        res.notes.append(L(
            f"MSD を出せません: 使えるフレームが {len(frames)} 個しかありません"
            + (f" (読んだのは {data.frame_source})" if data.frame_source else "")
            + "。GROMACS の .xtc / .trr は読まないので、gmx trjconv で xyz などに書き出してから解析してください。",
            f"no MSD: only {len(frames)} usable frame(s)"
            + (f" (read from {data.frame_source})" if data.frame_source else "")
            + ". GROMACS .xtc / .trr are not read; convert them with gmx trjconv first."))
        unwrap = None
    if unwrap is not None:
        pos, syms = unwrap.result()
        drift = {"removed": False, "displacement_A": None, "max_displacement_A": None, "requested": bool(opts.msd_remove_drift)}
        if opts.msd_remove_drift:
            pos, info = compute.remove_com_drift(pos, syms)
            drift = {**info, "requested": True}
        nb = max(0, int(opts.msd_error_blocks))
        kw = dict(symbols=syms, dt_fs=dt_used, fit_fs=opts.msd_fit_fs, axes=opts.msd_axes, remove_drift=False, n_blocks=nb or 5)
        a = compute.msd_analysis(pos, opts.msd_species, None, error=bool(nb), **kw)
        t, m, D, rng = a["lag"], a["msd_A2"], a["D_cm2_s"], a["fit_range_fs"]
        by = {}
        fig, ax = plt.subplots(figsize=(6, 3.2))
        ax.plot(t, m, lw=1.4, color="k", label=opts.msd_species or "all atoms")
        elems = sorted(set(syms))
        if len(elems) > 1:
            for el in elems:
                ae = compute.msd_analysis(pos, el, None, error=bool(nb), **kw)
                by[el] = {"D_cm2_s": ae["D_cm2_s"], "last_A2": float(ae["msd_A2"][-1]), "n_atoms": int(syms.count(el)),
                          "msd_A2": ae["msd_A2"].tolist(), "loglog_slope": ae["loglog_slope"],
                          "D_err_cm2_s": (ae["error"] or {}).get("d_err_cm2_s"), "n_blocks": (ae["error"] or {}).get("n_blocks")}
                if opts.msd_species is None:
                    ax.plot(ae["lag"], ae["msd_A2"], lw=1.0, ls="--", label=el)
        if rng is not None:
            how = "user" if a["fit_range_user"] else f"{compute.DEFAULT_FIT_FRACTION[0]:.0%}-{compute.DEFAULT_FIT_FRACTION[1]:.0%} of max lag"
            ax.axvspan(rng[0], rng[1], color="tab:orange", alpha=0.15, label=f"fit range {rng[0]:g}-{rng[1]:g} fs ({how})")
        sub = f"{a['formula']}, axes {a['axes']}, " + (
            f"COM drift removed ({drift['displacement_A']:.3f} Å)" if drift.get("displacement_A") is not None else "COM drift NOT removed")
        ax.set_title(sub, fontsize=8, color="0.3")
        ax.set_xlabel("lag time [fs]" if dt_used else "lag [frames]"); ax.set_ylabel("MSD [Å²]"); plotstyle.grid(ax); ax.legend(fontsize=8)
        save(fig, "msd")
        if opts.msd_per_atom:
            idx = [i for i, sym in enumerate(syms) if sym == opts.msd_species] if opts.msd_species else list(range(len(syms)))
            comps = compute.parse_axes(opts.msd_axes)
            each = compute.msd_per_atom(pos[:, idx, :][:, :, comps], dt_used or 1.0, len(comps), opts.msd_fit_fs)
            res.tables["msd_per_atom"] = {
                "species": opts.msd_species, "atom_index": idx, "symbol": [syms[i] for i in idx],
                "D_cm2_s": each["d_cm2_s"], "fit_range_fs": each["fit_range_fs"], "spread": each["spread"],
                "note": L("原子 1 個ごとに、時間原点を全部使った MSD を直線に当てはめた D です。"
                          "速い・遅いの判定はしていません (1 原子の統計は全体より悪く、ばらつきます)。",
                          "D per atom, from a straight-line fit to its own multiple-time-origin MSD. "
                          "No fast/slow judgement is made; single-atom statistics are much noisier than the average.")}
            if each["d_cm2_s"]:
                good = [v for v in each["d_cm2_s"] if v is not None]
                if good:
                    res.notes.append(L(f"原子ごとの D: {len(good)} 原子、{min(good):.3g}〜{max(good):.3g} cm²/s "
                                       f"(平均 {sum(good) / len(good):.3g})。表は analysis/summary.json の msd_per_atom",
                                       f"D per atom: {len(good)} atoms, {min(good):.3g} to {max(good):.3g} cm^2/s "
                                       f"(mean {sum(good) / len(good):.3g}); the table is in msd_per_atom of analysis/summary.json"))
        if opts.vacf:
            _add_vacf(res, pos, syms, frames, opts, dt_used, save)
        if opts.conductivity_charge:
            _add_conductivity(res, spec, data, syms, frames, opts, D)
        if opts.vanhove:
            _add_vanhove(res, pos, syms, frames, opts, dt_used, save, run_dir)
        n_atoms = len(frames[0]) if frames else 0
        span_fs = (len(frames) - 1) * dt_used if dt_used else None
        periodic = bool(np.any(frames[0].pbc)) if frames else False
        scale = {"n_atoms": n_atoms, "n_frames": len(frames), "span_fs": span_fs, "periodic": periodic}
        res.tables["msd_scale"] = scale
        res.notes.append(L(
            "この拡散係数を出した計算の規模: 原子 {} 個、フレーム {} 枚{}、周期境界 {}。"
            "値の良し悪しは判定しません".format(
                n_atoms, len(frames),
                f"、全体で {span_fs:g} fs" if span_fs else "",
                "あり" if periodic else "なし (拡散は普通 周期境界のある系で測ります)"),
            "the size of the run behind this diffusion coefficient: {} atoms, {} frames{}, periodic boundaries {}. "
            "ADIT does not judge whether the value is good".format(
                n_atoms, len(frames),
                f", {span_fs:g} fs in total" if span_fs else "",
                "on" if periodic else "off (diffusion is normally measured with periodic boundaries)")))
        res.tables["msd"] = {"species": opts.msd_species, "last_A2": float(m[-1]), "D_cm2_s": D, "n_frames": len(frames),
                             "method": "fft_multiple_time_origins", "dt_fs": dt_used, "fit_range_fs": list(rng) if rng else None,
                             "fit_range_user": opts.msd_fit_fs is not None, "lag_fs": t.tolist(), "msd_A2": m.tolist(), "by_element": by,
                             "axes": a["axes"], "dimension": a["dimension"], "formula": a["formula"], "drift": drift,
                             "fit_fraction": a["fit_fraction"], "loglog_slope": a["loglog_slope"],
                             "unwrap_check": {"max_step_fraction_of_shortest_cell_width": unwrap.max_step_fraction,
                                              "steps_at_or_above_0_4": unwrap.large_step_count},
                             "D_err_cm2_s": (a["error"] or {}).get("d_err_cm2_s"), "D_error": a["error"],
                             "definition": L(
                                 "D は当てはめ範囲で MSD を直線に当てはめた傾きから。D_error は各ブロックから出した D の平均の標準誤差で、"
                                 "全軌跡の D 自体の厳密な誤差ではありません。当てはめ範囲が 1 ブロックに収まらないときは、収まる範囲まで上限を下げて"
                                 "求めます。loglog_slope は当てはめ範囲での log MSD 対 log t の傾き (拡散なら 1 に近い)",
                                 "D comes from a straight-line fit of the MSD over the fit range. D_error is the standard error of the mean "
                                 "block D, not a rigorous error on D from the full trajectory. loglog_slope is the slope of log MSD vs log t "
                                 "over the fit range (close to 1 for diffusion)")}
        if unwrap.large_step_count:
            res.notes.append(L(
                f"MSD の境界越え補正: {unwrap.large_step_count} 個のフレーム間で、最小像の移動が最短セル幅の 40 % 以上でした "
                f"(最大 {unwrap.max_step_fraction:.1%})。1 フレームの間に半セル以上動くと移動方向を一意に復元できないため、間引く前の軌跡でも確認してください",
                f"MSD boundary unwrapping: the minimum-image displacement was at least 40% of the shortest cell width between "
                f"{unwrap.large_step_count} pairs of frames (maximum {unwrap.max_step_fraction:.1%}). If an atom moves by half a cell or more "
                "between frames, its direction cannot be reconstructed uniquely; also inspect the trajectory before subsampling"))
        if drift.get("displacement_A") is not None and dt_used:
            res.notes.append(L(f"MSD の前に系全体の質量重心を各フレームから引きました (重心は最初と最後で {drift['displacement_A']:.3f} Å、"
                               f"最初のフレームからの最大 {drift['max_displacement_A']:.3f} Å 動いています)",
                               f"the center-of-mass of the whole system was subtracted from every frame before the MSD (it moved "
                               f"{drift['displacement_A']:.3f} Å between the first and last frame, at most {drift['max_displacement_A']:.3f} Å)"))
        elif not opts.msd_remove_drift:
            res.notes.append(L("MSD から重心の流れを落としていません (--msd-keep-drift)。重心が流れていると MSD に並進の分が残ります",
                               "the center-of-mass drift was NOT removed from the MSD (--msd-keep-drift); any drift of the center of mass stays in the MSD"))
        if D is not None and D < 0:
            res.notes.append(L("拡散係数が負の値です。軌跡が短く統計が足りない可能性があります (拡散係数は MSD が時間に比例して増える範囲で読みます)",
                               "the diffusion coefficient is negative; the trajectory may be too short for the statistics (read D where the MSD grows linearly in time)"))
        if dt_used is None and len(frames) > 1:
            res.notes.append(L("軌跡の 1 フレームあたりの時間が出力から読めないので、MSD の横軸はフレームの番号で、拡散係数は出しません",
                               "the time per trajectory frame cannot be read from the output, so the MSD is plotted against frame lag and no diffusion coefficient is given"))
    if (opts.rdf or opts.msd) and 1 < len(frames) < FEW_FRAMES:
        res.notes.append(L(f"軌跡は {len(frames)} フレームです (RDF と MSD は、フレームが少ないと 1 フレームの揺らぎで形が変わります)",
                           f"the trajectory has {len(frames)} frames (with few frames, RDF and MSD change shape with single-frame fluctuations)"))
    elif (opts.rdf or opts.msd) and len(frames) <= 1:
        res.notes.append(L(f"軌跡がありません ({len(frames)} フレーム)。RDF と MSD は MD の結果で意味を持ちます",
                           f"no trajectory ({len(frames)} frame); RDF and MSD are meaningful for MD results"))
    if zacc is not None and zacc.n_frames:
        z = zacc.result()
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 5), sharex=True)
        for el, v in z["number_density_A3"].items():
            ax1.plot(z["z_A"], v, lw=1.0, label=el)
        ax1.set_ylabel("number density [Å$^{-3}$]"); ax1.legend(fontsize=8); plotstyle.grid(ax1)
        ax2.plot(z["z_A"], z["mass_density_g_cm3"], lw=1.0, color="k")
        ax2.set_ylabel("density [g/cm$^3$]"); ax2.set_xlabel("z [Å]" if z["c_perpendicular"] else "height normal to the a-b plane [Å]")
        plotstyle.grid(ax2)
        save(fig, "zdensity")
        (out / "zdensity.json").write_text(json.dumps(_jsonable(z)), encoding="utf-8")
        res.tables["zdensity"] = {"axis_length_A": z["axis_length_A"], "bin_A": z["bin_A"], "n_frames": z["n_frames"], "elements": list(z["number_density_A3"]),
                                  "c_perpendicular": z["c_perpendicular"], "file": str(out / "zdensity.json")}
        if not z["c_perpendicular"]:
            res.notes.append(L("セルの c が a-b 面に垂直でないので、z 方向の密度分布の横軸は a-b 面に垂直な方向の高さ (分率座標の 3 番目 × 面の間隔) です",
                               "the cell vector c is not perpendicular to the a-b plane, so the density profile is along the normal of the a-b plane (third fractional coordinate x interplanar spacing)"))
    if opts.stats and is_md:
        series = {}
        if data.temperatures_k and len(data.temperatures_k) - opts.skip_frames >= MIN_STATS_POINTS:
            series["temperature"] = (data.temperatures_k[opts.skip_frames:], "K", _spacing(data.times_fs))
        if data.energies_ev and len(data.energies_ev) - opts.skip_frames >= MIN_STATS_POINTS:
            series["energy"] = (data.energies_ev[opts.skip_frames:], "eV", _spacing(data.times_fs))
        dens = [x for x in densities if x is not None]
        if len(dens) >= MIN_STATS_POINTS and np.ptp(dens) > 1e-9 * max(dens):
            series["density"] = (dens, "g/cm^3", dt_used)
        for key, s in data.series.items():
            if not isinstance(s["values"], (list, tuple)):
                continue
            v = [x for x in s["values"][opts.skip_frames:] if np.isfinite(x)]
            if len(v) >= MIN_STATS_POINTS and np.ptp(v) > 0:
                series[key] = (v, s["unit"], _spacing(data.times_fs))
        if series:
            st = {}
            fig, ax = plt.subplots(figsize=(6, 3.2))
            for key, (x, unit, dt) in series.items():
                s = compute.series_stats(x, dt)
                s["unit"] = unit
                st[key] = s
                bs = [b["block_size"] for b in s["blocks"]]
                sem = np.array([b["sem"] for b in s["blocks"]])
                err = np.array([b["sem_err"] for b in s["blocks"]])
                ref = sem[0] if sem[0] > 0 else 1.0
                ax.errorbar(bs, sem / ref, yerr=err / ref, marker="o", ms=3, lw=1.0, capsize=2, label=f"{key} [{unit}] (×{ref:.2g})")
            ax.set_xscale("log", base=2); ax.set_xlabel("block size (points per block)")
            ax.set_ylabel("error of the mean / unblocked value"); plotstyle.grid(ax); ax.legend(fontsize=8)
            save(fig, "blocking")
            res.tables["timeseries_stats"] = st
            res.tables["timeseries_stats_definition"] = L(
                "blocks: Flyvbjerg–Petersen のブロック平均 (J. Chem. Phys. 91, 461 (1989))。束ねるごとの平均値の標準誤差 sem とその誤差 sem_err。"
                "tau_int_samples: 積分自己相関時間 (1/2 + Σρ、Sokal の自動の窓 c = 5)、sem_acf = 標準偏差 × sqrt(2 τ / N)。どの値を採るかは判定しない",
                "blocks: Flyvbjerg-Petersen block averaging (J. Chem. Phys. 91, 461 (1989)); standard error of the mean sem and its error sem_err per blocking step. "
                "tau_int_samples: integrated autocorrelation time (1/2 + sum rho, Sokal automatic window c = 5); sem_acf = std x sqrt(2 tau / N). No value is chosen for you")
    if exporter is not None:
        extra = []
        if data.code == "gromacs":
            extra = [L("GROMACS の軌跡 (adit.xtc) は書き出していません。書き出したのは最後の構造 (adit.gro) の 1 フレームです。",
                       "The GROMACS trajectory (adit.xtc) is not exported; only the final structure (adit.gro) is written as one frame."),
                     L("軌跡は GROMACS で PDB にしてから、VMD・OVITO・TRAVIS で開いてください (出力のグループを聞かれたら 0 = System):",
                       "Convert the trajectory to PDB with GROMACS first, then open it in VMD, OVITO or TRAVIS (answer 0 = System when asked for the output group):"),
                     "  echo 0 | gmx trjconv -f adit.xtc -s adit.tpr -o traj.pdb", ""]
        info = exporter.close(run_dir=run_dir, code=data.code, source=data.frame_source or "?", dt_frame_fs=data.frame_dt_fs, stride=opts.stride,
                              skip=skip, n_total=n_all, rdf_cutoff=float(min(6.0, rmax)) if periodic else 6.0, extra_lines=extra)
        res.tables["export"] = info
    if opts.dos and data.eigenvalues_ev is not None:
        x, y = compute.dos(data.eigenvalues_ev, data.eigen_weights, sigma=opts.dos_sigma)
        fig, ax = plt.subplots(figsize=(6, 3.2))
        ax.plot(x, y, lw=1.2)
        if data.fermi_ev is not None:
            ax.axvline(data.fermi_ev, color="gray", ls="--", lw=0.8); fermi_used = True
        ax.set_xlabel("E [eV]"); ax.set_ylabel("DOS"); plotstyle.grid(ax)
        save(fig, "dos")
        res.tables["dos"] = {"n_eigen": int(len(data.eigenvalues_ev)), "fermi_ev": data.fermi_ev}
    if opts.vibrations and data.frequencies_cm1:
        if opts.freq_scale:
            data.frequencies_cm1 = [f * opts.freq_scale for f in data.frequencies_cm1]
            res.notes.append(L(
                f"振動数に補正係数 {opts.freq_scale:g} を掛けました (あなたが指定した値です。"
                "ADIT は既定値を持ちません。図・要約・frequencies.csv はすべて掛けたあとの値です)",
                f"the frequencies were scaled by {opts.freq_scale:g} (the value you gave; ADIT has no default). "
                "The figure, the summary and frequencies.csv all show the scaled values"))
        measured = _measured_spectrum(res, opts)
        x, y = compute.spectrum(data.frequencies_cm1, data.ir_intensities)
        fig, ax = plt.subplots(figsize=(6, 3.2))
        ax.plot(x, y, lw=1.0)
        for f in data.frequencies_cm1:
            if f > 1:
                ax.axvline(f, color="gray", lw=0.5, alpha=0.5)
        ylabel = ("IR intensity (unit as printed by the code), Gaussian FWHM 47 cm$^{-1}$" if data.ir_intensities
                  else "(equal weights), Gaussian FWHM 47 cm$^{-1}$")
        if measured is not None:
            ax.plot(measured[0], measured[1] * (y.max() if y.size else 1.0), lw=1.0, color="tab:gray",
                    ls="--", label="measured (scaled to max)")
            ax.legend(fontsize=8)
        ax.set_xlabel("wavenumber [cm$^{-1}$]"); ax.set_ylabel(ylabel); plotstyle.grid(ax)
        save(fig, "spectrum")
        res.tables["frequencies"] = [float(f) for f in data.frequencies_cm1]
        if data.ir_intensities:
            res.tables["ir_intensities"] = [float(v) for v in data.ir_intensities]
        rows = [("frequency_cm-1", "ir_intensity", "raman_activity")]
        n = len(data.frequencies_cm1)
        ir = list(data.ir_intensities or []) + [""] * n
        ram = list(data.raman_activities or []) + [""] * n
        off = n - len(data.raman_activities) if data.raman_activities and len(data.raman_activities) < n else 0
        for i, f in enumerate(data.frequencies_cm1):
            ram_value = ram[i - off] if off <= i < off + len(data.raman_activities or []) else ""
            rows.append((f"{f:.6f}", f"{ir[i]:.6f}" if isinstance(ir[i], float) else "",
                         f"{ram_value:.6f}" if isinstance(ram_value, float) else ""))
        path = out / "frequencies.csv"
        path.write_text("\n".join(",".join(str(x) for x in row) for row in rows) + "\n", encoding="utf-8")
        res.tables["frequencies_csv"] = str(path)
        if data.raman_activities:
            ram = list(data.raman_activities)
            freqs = list(data.frequencies_cm1)
            if len(ram) != len(freqs):
                freqs = freqs[len(freqs) - len(ram):]
            xr, yr = compute.spectrum(freqs, ram)
            fig, ax = plt.subplots(figsize=(6, 3.2))
            ax.plot(xr, yr, lw=1.0, color="tab:red", label="calculated")
            if measured is not None:
                ax.plot(measured[0], measured[1] * (yr.max() if yr.size else 1.0), lw=1.0, color="tab:gray",
                        ls="--", label="measured (scaled to max)")
                ax.legend(fontsize=8)
            ax.set_xlabel("wavenumber [cm$^{-1}$]")
            ax.set_ylabel("Raman activity (unit as printed by the code), Gaussian FWHM 47 cm$^{-1}$")
            plotstyle.grid(ax)
            save(fig, "raman")
            res.tables["raman_activities"] = [float(v) for v in ram]
            peak = int(np.argmax(ram))
            res.notes.append(L(
                f"ラマン活性を {len(ram)} 本読みました (最大は {freqs[peak]:.1f} cm⁻¹ の {ram[peak]:.3g})。"
                "赤外の強度とは別の量です。測定との一致は判定しません",
                f"read {len(ram)} Raman activities (the largest is {ram[peak]:.3g} at {freqs[peak]:.1f} cm⁻¹); "
                "this is a different quantity from the infrared intensity, and agreement with experiment is not judged"))
        res.notes.append(L(
            "振動数と強度の表を frequencies.csv に書きました (1 列目 cm⁻¹、2 列目 赤外の強度、3 列目 ラマン活性)。"
            "図は幅 47 cm⁻¹ (半値全幅) のガウス関数で広げたものです。強度の単位はコードが書いたままで、"
            "測定と重ねるときは縦軸を自分で合わせてください",
            "the frequencies and intensities are in frequencies.csv (cm^-1, IR intensity, Raman activity); "
            "the figure broadens them with a Gaussian of 47 cm^-1 FWHM. The intensity unit is whatever the code "
            "printed, so scale the vertical axis yourself when overlaying a measured spectrum"))
        if task_type == "vibrations" or (task_type is None and len(data.frames) <= 1):
            res.notes.append(_vibration_note(data, spec))
    if opts.bands:
        from adit.analysis.bands import load_bands, plot_bands
        bd = load_bands(run_dir, data.code, data.fermi_ev)
        if bd is not None:
            plot_bands(bd, out / "bands.png", opts.bands_window_ev)
            res.figures["bands"] = str(out / "bands.png")
            e = bd.energies_ev
            gap = None
            if bd.fermi_ev is not None:
                occ = e[e <= bd.fermi_ev]; emp = e[e > bd.fermi_ev]
                if occ.size and emp.size:
                    gap = float(emp.min() - occ.max())
            from adit.analysis.bands import gap_details

            details = gap_details(e, bd.fermi_ev, bd.kpts_frac, dict(bd.labels))
            res.tables["bands"] = {"n_kpoints": int(e.shape[0]), "n_bands": int(e.shape[1]), "labels": [l for _, l in bd.labels],
                                   "gap_ev": gap, "fermi_ev": bd.fermi_ev,
                                   "reference_level": "highest_occupied" if data.fermi_is_homo else "fermi",
                                   "gap_details": details}
            if details is not None:
                where = lambda key: (f"{details[key + '_label']} " if details.get(key + "_label") else "") + f"(k 点 {details[key + '_kpoint_index'] + 1})"
                res.notes.append(L(
                    f"バンドギャップ {details['gap_ev']:.3f} eV: 価電子帯の頂上は {where('vbm')}、伝導帯の底は {where('cbm')}。"
                    + ("同じ k 点にあります (直接遷移)" if details["direct"] else
                       f"別の k 点です (間接遷移)。同じ k 点で最小の間隔は {details['direct_gap_ev']:.3f} eV")
                    + "。与えた経路の中での話で、経路の外はこの計算では分かりません",
                    f"band gap {details['gap_ev']:.3f} eV: the valence-band maximum is at k-point "
                    f"{details['vbm_kpoint_index'] + 1}, the conduction-band minimum at k-point {details['cbm_kpoint_index'] + 1}"
                    + (" (same k-point: direct)" if details["direct"] else
                       f" (different k-points: indirect); the smallest gap at a single k-point is {details['direct_gap_ev']:.3f} eV")
                    + ". This is within the path you gave; nothing is known about k-points outside it"))
            fermi_used = fermi_used or bd.fermi_ev is not None
            if opts.effective_mass:
                _add_effective_mass(res, bd, opts)
    _add_electronic_extras(res, run_dir, opts)
    if (opts.coordination_cutoff or opts.centrosymmetry_neighbors or opts.steinhardt_cutoff or opts.cluster_cutoff
            or opts.adf is not None or opts.hbond or opts.radius_of_gyration or opts.density_grid):
        _add_local_order(res, data, frames, opts, dt_used, save, out)
    if opts.distances or opts.angles or opts.dihedrals or opts.rmsd_reference or opts.rmsf:
        _add_geometry_series(res, data, frames, opts, dt_used, save, out)
    if (opts.displacement_reference or opts.voronoi or opts.sasa or opts.pca):
        _add_structure_extras(res, data, frames, opts, save, out)
    fermi_used = _add_stage23(res, data, run_dir, out, opts) or fermi_used
    if fermi_used and data.fermi_is_homo:
        res.notes.append(L("フェルミ準位の代わりに最高被占準位 (highest occupied level) を使った", "the highest occupied level is used instead of the Fermi level"))
    res.tables = _jsonable(res.tables)
    (out / "summary.json").write_text(json.dumps({
        "code": res.code, "tables": _paths_relative_to(res.tables, run_dir, out),
        "figures": _paths_relative_to(res.figures, run_dir, out), "notes": res.notes},
        ensure_ascii=False, indent=2, default=float), encoding="utf-8")
    (out / "summary.txt").write_text(res.summary_text() + "\n", encoding="utf-8")
    if opts.viscosity:
        _add_viscosity(res, data, out, save)
    if opts.xrd:
        _add_xrd(res, data, out, opts, save)
    return res


def _paths_relative_to(obj, run_dir: Path, out: Path):
    # summary.json is copied with the directory to other machines: write paths relative to run_dir
    # (or to the analysis directory when -o put it elsewhere). AnalysisResult itself keeps absolute paths.
    if isinstance(obj, dict):
        return {k: _paths_relative_to(v, run_dir, out) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_paths_relative_to(v, run_dir, out) for v in obj]
    if isinstance(obj, str) and obj and Path(obj).is_absolute():
        for base in (run_dir, out):
            try:
                return Path(obj).resolve().relative_to(Path(base).resolve()).as_posix()
            except ValueError:
                continue
    return obj


def _spacing(times_fs: list[float] | None) -> float | None:
    return float(times_fs[1] - times_fs[0]) if times_fs and len(times_fs) >= 2 else None


def _charge_csv(out: Path, item: dict, used: set[str]) -> str:
    import csv
    import re as _re

    name = _re.sub(r"[^0-9A-Za-z._-]+", "_", str(item.get("definition") or "charges")).strip("_") or "charges"
    base = name
    k = 2
    while name in used:
        name = f"{base}_{k}"; k += 1
    used.add(name)
    p = out / f"charges_{name}.csv"
    atoms = item.get("atoms") or [str(i + 1) for i in range(len(item["values"]))]
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["index", "atom", f"charge_{item.get('unit', 'e')}"])
        for i, (a, v) in enumerate(zip(atoms, item["values"]), start=1):
            w.writerow([i, a, f"{float(v):.6f}"])
    return str(p)


def _add_viscosity(res: AnalysisResult, data: RunData, out: Path, save) -> None:
    import matplotlib.pyplot as plt

    from adit.analysis.transport import TransportError, green_kubo_viscosity

    tensor = data.series.get("pressure_tensor")
    if not tensor:
        res.notes.append(L("粘度を出せません: 圧力テンソルの非対角成分 (pxy, pxz, pyz) が出力にありません。"
                           "いま対応しているのは LAMMPS だけで、生成のときに「thermo に圧力テンソルを足す」欄 "
                           "(method.thermo_pressure_tensor) を有効にすると出ます。GROMACS・OpenMM からは読めません",
                           "no viscosity: the output has no off-diagonal pressure components (pxy, pxz, pyz). "
                           "Only LAMMPS is supported: enable method.thermo_pressure_tensor when generating. "
                           "GROMACS and OpenMM outputs cannot be read for this"))
        return
    if not data.times_fs or data.final is None:
        res.notes.append(L("粘度を出せません: 時刻または最終構造 (体積) を読めません",
                           "no viscosity: the times or the final structure (for the volume) could not be read"))
        return
    temperatures = data.temperatures_k or []
    if not temperatures:
        res.notes.append(L("粘度を出せません: 温度の時系列がありません", "no viscosity: there is no temperature series"))
        return
    volume = float(data.final.get_volume()) if data.final.pbc.all() else 0.0
    values = tensor["values"]
    try:
        result = green_kubo_viscosity(data.times_fs, [values.get("Pxy"), values.get("Pxz"), values.get("Pyz")],
                                      volume, float(np.mean(temperatures)), pressure_unit=tensor.get("unit", "bar"))
    except TransportError as ex:
        res.notes.append(str(ex))
        return
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.plot(result.times_fs, result.running_pa_s * 1e3)
    ax.set_xlabel("lag time [fs]")
    ax.set_ylabel("running integral of eta [mPa s]")
    save(fig, "viscosity")
    table = result.as_dict()
    msd = res.tables.get("msd") or {}
    if msd.get("loglog_slope") is not None:
        table["msd_loglog_slope"] = msd["loglog_slope"]
        table["note_with_msd"] = L(
            f"同じ軌跡の MSD の log-log の傾きは {msd['loglog_slope']:.2f} です (拡散なら 1 に近い)。"
            "傾きが 0 に近い系 (固体) では、この粘度の値を液体の粘度として読まないでください。",
            f"the log-log slope of the MSD for the same trajectory is {msd['loglog_slope']:.2f} (close to 1 for diffusion); "
            "if it is near 0 (a solid), do not read this number as a liquid viscosity.")
    res.tables["viscosity"] = table
    if table.get("note_with_msd"):
        res.notes.append(table["note_with_msd"])
    res.notes.append(L(f"粘度 (Green-Kubo、積分の最後の値): {result.eta_pa_s * 1e3:.4g} mPa·s "
                       f"(温度 {result.temperature_k:.1f} K、体積 {result.volume_ang3:.1f} Å³、{result.points} 点)。{result.note}",
                       f"viscosity (Green-Kubo, value at the end of the integral): {result.eta_pa_s * 1e3:.4g} mPa s "
                       f"(temperature {result.temperature_k:.1f} K, volume {result.volume_ang3:.1f} A^3, {result.points} points). {result.note}"))


def _selected_atoms(res: AnalysisResult, frames, syms, opts: AnalysisOptions) -> list[int]:
    if opts.select:
        from adit.analysis.select import SelectionError, describe, select

        try:
            idx = select(frames[0], opts.select).tolist()
        except SelectionError as ex:
            res.notes.append(str(ex))
            return []
        res.notes.append(L(f"選び方 {opts.select!r}: ", f"selection {opts.select!r}: ") + describe(frames[0], opts.select))
        return idx
    if opts.msd_species:
        return [i for i, sym in enumerate(syms) if sym == opts.msd_species]
    return list(range(len(syms)))


def _add_effective_mass(res: AnalysisResult, bd, opts: AnalysisOptions) -> None:
    from adit.analysis.effective_mass import at_band_edges

    try:
        masses = at_band_edges(bd.x_axis(), bd.energies_ev, bd.fermi_ev, opts.effective_mass_points)
    except ValueError as ex:
        res.notes.append(L(f"有効質量を出せません: {ex}", f"cannot compute the effective mass: {ex}"))
        return
    if not masses:
        return
    res.tables["effective_mass"] = [
        {"kind": m.kind, "m_over_me": m.m_over_me, "band_index": m.band_index,
         "kpoint_index": m.kpoint_index, "n_points": m.n_points, "r_squared": m.r_squared,
         "k_window_inv_ang": m.k_window_inv_ang} for m in masses]
    for m in masses:
        who = L("正孔 (価電子帯の頂上)", "holes (valence-band maximum)") if m.kind == "hole" else \
            L("電子 (伝導帯の底)", "electrons (conduction-band minimum)")
        res.notes.append(L(
            f"有効質量 {who}: m*/m_e = {m.m_over_me:+.3f} ("
            f"{m.n_points} 点、幅 {m.k_window_inv_ang:.4f} 1/Å、決定係数 {m.r_squared:.4f})。"
            "経路に沿った向きの値で、テンソルではありません",
            f"effective mass for {who}: m*/m_e = {m.m_over_me:+.3f} ("
            f"{m.n_points} points, width {m.k_window_inv_ang:.4f} 1/Å, R^2 {m.r_squared:.4f}); "
            "this is along the path direction, not a tensor"))


def _add_electronic_extras(res: AnalysisResult, run_dir: Path, opts: AnalysisOptions) -> None:
    if opts.projected_bands:
        from adit.analysis.projected_bands import find_filproj, read_filproj

        found = find_filproj(run_dir)
        if found is not None:
            try:
                pb = read_filproj(found)
            except ValueError as ex:
                res.notes.append(L(f"射影バンドを読めません: {ex}", f"cannot read the projected bands: {ex}"))
            else:
                elements = sorted({p.element for p in pb.projections})
                res.tables["projected_bands"] = {
                    "n_kpoints": pb.n_kpoints, "n_bands": pb.n_bands,
                    "total_mean": float(pb.total().mean()),
                    "by_element": {el: float(pb.by_element(el).mean()) for el in elements},
                    "orbitals": [{"atom": p.atom_index, "element": p.element, "label": p.label,
                                  "l": p.l, "mean_weight": float(p.weights.mean())}
                                 for p in pb.projections]}
                per_el = ", ".join(f"{el} {pb.by_element(el).mean():.3f}" for el in elements)
                res.notes.append(L(
                    f"射影バンド ({found.name}): {pb.n_kpoints} k 点 × {pb.n_bands} バンド、"
                    f"原子軌道 {len(pb.projections)} 個。元素ごとの重みの平均 {per_el}、"
                    f"全射影の和の平均 {pb.total().mean():.3f} (1 に足りない分は原子軌道で表しきれない成分です)",
                    f"projected bands ({found.name}): {pb.n_kpoints} k-points x {pb.n_bands} bands, "
                    f"{len(pb.projections)} atomic orbitals; mean weight per element {per_el}, "
                    f"mean of all projections {pb.total().mean():.3f} (the rest is not representable by atomic orbitals)"))
    if opts.projected_bands:
        from adit.analysis.procar import read_procar

        for name in ("PROCAR", "PROCAR.gz"):
            path = run_dir / name
            if not path.is_file():
                continue
            try:
                pc = read_procar(path)
            except (ValueError, IndexError) as ex:
                res.notes.append(L(f"PROCAR を読めません: {ex}", f"cannot read the PROCAR: {ex}"))
                break
            nspin, nk, nb = pc.energies_ev.shape
            res.tables["procar"] = {"n_kpoints": nk, "n_bands": nb, "n_ions": int(pc.projections.shape[3]),
                                    "n_spins": nspin, "orbitals": pc.orbitals,
                                    "total_mean": float(pc.total().mean())}
            res.notes.append(L(
                f"PROCAR ({name}): {nk} k 点 × {nb} バンド × 原子 {pc.projections.shape[3]} 個"
                + ("、スピン 2 つ" if nspin == 2 else "")
                + f"、軌道 {', '.join(pc.orbitals)}。全射影の和の平均 {pc.total().mean():.3f}",
                f"PROCAR ({name}): {nk} k-points x {nb} bands x {pc.projections.shape[3]} ions"
                + (", two spins" if nspin == 2 else "")
                + f", orbitals {', '.join(pc.orbitals)}; mean of all projections {pc.total().mean():.3f}"))
            break
    if opts.optical:
        from adit.analysis.optical import read_epsilon

        try:
            op = read_epsilon(run_dir)
        except (FileNotFoundError, ValueError):
            op = None
        if op is not None:
            eps0 = op.static_dielectric()
            alpha = op.absorption_cm()
            res.tables["optical"] = {
                "n_energies": int(op.energy_ev.size), "energy_max_ev": float(op.energy_ev.max()),
                "static_eps1": [float(x) for x in eps0],
                "absorption_max_cm": [float(x) for x in alpha.max(axis=0)]}
            res.notes.append(L(
                f"光学 (epsilon.x、{op.energy_ev.size} 点、0〜{op.energy_ev.max():.2f} eV): "
                f"静的誘電率 eps1 は x {eps0[0]:.3f} / y {eps0[1]:.3f} / z {eps0[2]:.3f}、"
                f"吸収係数の最大は {alpha.max():.3e} 1/cm。k 点の数と広がりの幅で値が変わります",
                f"optics (epsilon.x, {op.energy_ev.size} points, 0-{op.energy_ev.max():.2f} eV): "
                f"static eps1 x {eps0[0]:.3f} / y {eps0[1]:.3f} / z {eps0[2]:.3f}, "
                f"largest absorption {alpha.max():.3e} 1/cm; the values depend on the k-point mesh and the broadening"))
    if opts.bader:
        from adit.analysis.bader import read_acf

        path = Path(opts.bader)
        if not path.is_file():
            res.notes.append(L(f"ACF.dat がありません: {path}", f"ACF.dat not found: {path}"))
            return
        valence = None
        if opts.bader_valence.strip():
            try:
                valence = [float(x) for x in opts.bader_valence.replace(" ", "").split(",") if x]
            except ValueError:
                res.notes.append(L("価電子数を数として読めません (例 --bader-valence 4,6,1)",
                                   "cannot read the valence electrons (e.g. --bader-valence 4,6,1)"))
                return
        try:
            bc = read_acf(path)
            net = bc.net_charge(valence) if valence is not None else None
        except ValueError as ex:
            res.notes.append(L(f"Bader 電荷を読めません: {ex}", f"cannot read the Bader charges: {ex}"))
            return
        res.tables["bader"] = {"electrons": [float(x) for x in bc.electrons],
                               "total_electrons": bc.total_electrons,
                               "net_charge": ([float(x) for x in net] if net is not None else None)}
        head = L(f"Bader 電荷 ({path.name}、原子 {bc.electrons.size} 個): 電子数 ",
                 f"Bader charges ({path.name}, {bc.electrons.size} atoms): electrons ")
        body = ", ".join(f"{e:.3f}" for e in bc.electrons[:12]) + (" ..." if bc.electrons.size > 12 else "")
        res.notes.append(head + body)
        if net is not None:
            res.notes.append(L("価電子との差 (正なら電子が減っています): ",
                               "difference from the valence electrons (positive means fewer electrons): ")
                             + ", ".join(f"{x:+.3f}" for x in net[:12]) + (" ..." if net.size > 12 else ""))
        else:
            res.notes.append(L("価電子数を --bader-valence で渡すと、価電子との差も出します (ADIT は推測しません)",
                               "pass the valence electrons with --bader-valence to also get the difference (ADIT does not guess them)"))


def _add_local_order(res: AnalysisResult, data: RunData, frames, opts: AnalysisOptions, dt_fs, save, out: Path) -> None:
    import matplotlib.pyplot as plt

    from adit.analysis import local_order as LO
    from adit.analysis.geometry_series import radius_of_gyration
    from adit.analysis.hbond import HydrogenBondError, count_series, criteria_note
    from adit.analysis.volumetric import DensityGrid, VolumetricError

    last = frames[-1] if len(frames) else None
    if last is None:
        return
    try:
        if opts.coordination_cutoff:
            got = LO.coordination(last, opts.coordination_cutoff)
            res.tables["coordination"] = got
            res.notes.append(L(f"配位数 (カットオフ {got['cutoff_A']:g} Å、最終構造): 平均 {got['mean']:.2f}、"
                               f"範囲 {min(got['coordination'])}〜{max(got['coordination'])}",
                               f"coordination number (cutoff {got['cutoff_A']:g} Å, final structure): mean {got['mean']:.2f}, "
                               f"range {min(got['coordination'])}-{max(got['coordination'])}"))
        if opts.centrosymmetry_neighbors:
            got = LO.centrosymmetry(last, opts.centrosymmetry_neighbors)
            res.tables["centrosymmetry"] = got
            values = np.array(got["centrosymmetry_A2"], dtype=float)
            res.notes.append(L(f"中心対称性パラメータ (相手 {got['n_neighbors']} 個): 最小 {np.nanmin(values):.4g}、"
                               f"最大 {np.nanmax(values):.4g} Å² (どこからを欠陥と呼ぶかは判定していません)",
                               f"centrosymmetry parameter ({got['n_neighbors']} neighbors): min {np.nanmin(values):.4g}, "
                               f"max {np.nanmax(values):.4g} Å² (no threshold for calling an atom a defect)"))
        if opts.steinhardt_cutoff:
            got = LO.steinhardt(last, opts.steinhardt_cutoff)
            res.tables["steinhardt"] = got
            res.notes.append(L(f"Steinhardt q4 = {np.nanmean(got['q4']):.4f}、q6 = {np.nanmean(got['q6']):.4f} "
                               f"(平均、カットオフ {got['cutoff_A']:g} Å)。どの値がどの構造かは言いません",
                               f"Steinhardt q4 = {np.nanmean(got['q4']):.4f}, q6 = {np.nanmean(got['q6']):.4f} "
                               f"(means, cutoff {got['cutoff_A']:g} Å); ADIT does not name the structures"))
        if opts.cluster_cutoff:
            got = LO.clusters(last, opts.cluster_cutoff)
            res.tables["clusters"] = got
            res.notes.append(L(f"かたまり (カットオフ {got['cutoff_A']:g} Å): {got['n_clusters']} 個、"
                               f"いちばん大きいもので {got['largest']} 原子",
                               f"clusters (cutoff {got['cutoff_A']:g} Å): {got['n_clusters']}, largest {got['largest']} atoms"))
        if opts.adf is not None:
            center, _, outer = (opts.adf or "").partition(",")
            got = LO.angle_distribution(frames, center=center or None, cutoff=opts.adf_cutoff, outer=outer or None)
            res.tables["angle_distribution"] = got
            fig, ax = plt.subplots(figsize=(6, 3.0))
            ax.plot(got["angle_deg"], got["counts"], lw=1.2)
            ax.set_xlabel("angle [deg]"); ax.set_ylabel("counts"); plotstyle.grid(ax)
            save(fig, "angle_distribution")
            top = got["angle_deg"][int(np.argmax(got["counts"]))] if got["n_angles"] else float("nan")
            res.notes.append(L(f"結合角の分布 ({center or '全元素'}、カットオフ {opts.adf_cutoff:g} Å): {got['n_angles']} 個の角度、"
                               f"いちばん多いのは {top:.0f} 度のあたり (帰属はしていません)",
                               f"angle distribution ({center or 'all elements'}, cutoff {opts.adf_cutoff:g} Å): "
                               f"{got['n_angles']} angles, most frequent near {top:.0f} degrees (no assignment is made)"))
    except LO.LocalOrderError as ex:
        res.notes.append(str(ex))
    if opts.hbond:
        try:
            distance, _, angle = opts.hbond.partition(",")
            got = count_series(frames, float(distance), float(angle))
        except (ValueError, HydrogenBondError) as ex:
            res.notes.append(L(f"水素結合を数えられません: {ex} 書き方は --hbond 3.5,150 (距離 [Å], 角度 [度]) です。",
                               f"cannot count hydrogen bonds: {ex} write it as --hbond 3.5,150 (distance in Å, angle in degrees).")
                             + criteria_note())
        else:
            res.tables["hbond"] = got
            fig, ax = plt.subplots(figsize=(6, 2.8))
            x = np.arange(len(got["counts"])) * (dt_fs if dt_fs else 1.0)
            ax.plot(x, got["counts"], lw=1.2)
            ax.set_xlabel("time [fs]" if dt_fs else "frame"); ax.set_ylabel("hydrogen bonds"); plotstyle.grid(ax)
            save(fig, "hbond")
            c = got["criteria"]
            res.notes.append(L(f"水素結合 (D–A ≤ {c['donor_acceptor_A']:g} Å、D–H···A ≥ {c['angle_deg']:g} 度): "
                               f"平均 {got['mean']:.2f} 本" + (f" (標準偏差 {got['std']:.2f})" if got["std"] else "")
                               + "。条件は指定された値です。" + criteria_note(),
                               f"hydrogen bonds (D-A <= {c['donor_acceptor_A']:g} Å, D-H...A >= {c['angle_deg']:g} deg): "
                               f"mean {got['mean']:.2f}" + (f" (sd {got['std']:.2f})" if got["std"] else "")
                               + ". The criteria are the ones you gave. " + criteria_note()))
    if opts.radius_of_gyration:
        got = radius_of_gyration(frames)
        res.tables["radius_of_gyration"] = got
        fig, ax = plt.subplots(figsize=(6, 2.8))
        x = np.arange(len(got["rg_A"])) * (dt_fs if dt_fs else 1.0)
        ax.plot(x, got["rg_A"], lw=1.2)
        ax.set_xlabel("time [fs]" if dt_fs else "frame"); ax.set_ylabel("Rg [Å]"); plotstyle.grid(ax)
        save(fig, "radius_of_gyration")
        res.notes.append(L(f"慣性半径 Rg: 平均 {np.mean(got['rg_A']):.3f} Å、最後 {got['rg_A'][-1]:.3f} Å (質量で重み付け)",
                           f"radius of gyration: mean {np.mean(got['rg_A']):.3f} Å, last {got['rg_A'][-1]:.3f} Å (mass weighted)"))
    if opts.density_grid:
        try:
            shape = tuple(int(v) for v in opts.density_grid.split(","))
            grid = DensityGrid(shape)
            for fr in frames:
                grid.add(fr)
            path = grid.write_cube(out / "density.cube")
        except (ValueError, VolumetricError) as ex:
            res.notes.append(L(f"3 次元の密度を出せません: {ex}", f"cannot accumulate the 3D density: {ex}"))
        else:
            res.tables["density_grid"] = {**{k: v for k, v in grid.result().items() if k != "values"}, "cube": str(path)}
            res.notes.append(L(f"3 次元の数密度を {path.name} に書きました ({'×'.join(map(str, shape))} の格子、"
                               f"{len(frames)} フレーム)。VMD・VESTA・OVITO で開いてください",
                               f"the 3D number density was written to {path.name} (grid {'x'.join(map(str, shape))}, "
                               f"{len(frames)} frames); open it in VMD, VESTA or OVITO"))


def _add_fes(res: AnalysisResult, table: dict, opts: AnalysisOptions, save) -> None:
    import matplotlib.pyplot as plt

    from adit.analysis.free_energy import free_energy_surface

    names = list(table)[:2]
    values = np.array([table[n] for n in names]).T
    try:
        fes = free_energy_surface(values if len(names) == 2 else values[:, 0],
                                  opts.fes_temperature_k, bins=opts.fes_bins, unit=opts.fes_unit)
    except ValueError as ex:
        res.notes.append(L(f"自由エネルギー面を出せません: {ex}", f"cannot compute the free-energy surface: {ex}"))
        return
    finite = fes.free_energy[np.isfinite(fes.free_energy)]
    res.tables["free_energy_surface"] = {
        "variables": names, "temperature_k": fes.temperature_k, "unit": fes.unit,
        "max": float(finite.max()), "empty_bins": int((fes.counts == 0).sum()),
        "centers": [[float(x) for x in c] for c in fes.centers],
        "free_energy": [[float(v) for v in row] for row in np.atleast_2d(fes.free_energy)]}
    filled = int((fes.counts > 0).sum())
    per_bin = float(fes.counts.sum()) / filled if filled else 0.0
    res.tables["free_energy_surface"]["points_per_filled_bin"] = per_bin
    res.notes.append(L(
        f"自由エネルギー面 A = -kT ln P ({', '.join(names)}、{fes.temperature_k:g} K、{fes.unit}): "
        f"最小を 0 にしたときの最大 {finite.max():.3f} {fes.unit}、点の入っていない区間 "
        f"{int((fes.counts == 0).sum())} 個、点の入っている区間 1 つあたり {per_bin:.2f} 点。"
        "収束したかどうかは判定しません",
        f"free-energy surface A = -kT ln P ({', '.join(names)}, {fes.temperature_k:g} K, {fes.unit}): "
        f"largest value {finite.max():.3f} {fes.unit} with the minimum set to zero, "
        f"{int((fes.counts == 0).sum())} empty bins, {per_bin:.2f} points per filled bin; "
        "ADIT does not judge convergence"))
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    if len(names) == 1:
        ax.plot(fes.centers[0], fes.free_energy, lw=1.4)
        ax.set_xlabel(names[0]); ax.set_ylabel(f"A [{fes.unit}]")
    else:
        masked = np.ma.masked_invalid(fes.free_energy.T)
        im = ax.pcolormesh(fes.centers[0], fes.centers[1], masked, shading="auto")
        fig.colorbar(im, ax=ax, label=f"A [{fes.unit}]")
        ax.set_xlabel(names[0]); ax.set_ylabel(names[1])
    plotstyle.grid(ax)
    save(fig, "free_energy")


def _measured_spectrum(res: AnalysisResult, opts: AnalysisOptions):
    if not opts.spectrum_measured:
        return None
    path = Path(opts.spectrum_measured).expanduser()
    if not path.is_file():
        res.notes.append(L(f"測定データのファイルがありません: {path}", f"measured spectrum not found: {path}"))
        return None
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        text = line.strip().replace(",", " ")
        if not text or text.startswith(("#", "%", "//")):
            continue
        parts = text.split()
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except (ValueError, IndexError):
            continue
    if len(rows) < 2:
        res.notes.append(L(f"{path.name} から測定の点を 2 つ以上 読めません (1 列目 波数 [cm⁻¹]、2 列目 強度)",
                           f"cannot read two or more points from {path.name} (column 1 wavenumber in cm^-1, column 2 intensity)"))
        return None
    arr = np.array(sorted(rows))
    y = arr[:, 1]
    top = float(np.abs(y).max())
    if top <= 0:
        res.notes.append(L(f"{path.name} の強度がすべて 0 です", f"all intensities in {path.name} are zero"))
        return None
    res.notes.append(L(
        f"測定したスペクトル {path.name} を重ねました ({arr.shape[0]} 点、"
        f"{arr[0, 0]:.1f}〜{arr[-1, 0]:.1f} cm⁻¹)。縦軸は互いの最大値で合わせただけで、一致の良し悪しは判定しません",
        f"overlaid the measured spectrum {path.name} ({arr.shape[0]} points, "
        f"{arr[0, 0]:.1f}-{arr[-1, 0]:.1f} cm^-1); the vertical axes are matched at their maxima only, "
        "and ADIT does not judge the agreement"))
    return arr[:, 0], y / top


def _add_structure_extras(res: AnalysisResult, data: RunData, frames, opts: AnalysisOptions, save, out: Path) -> None:
    import matplotlib.pyplot as plt

    if not frames:
        res.notes.append(L("構造がありません", "there are no structures"))
        return
    last = frames[-1]
    cell = np.asarray(last.cell, dtype=float) if bool(np.any(last.pbc)) else None
    if opts.displacement_reference:
        from adit.analysis.displacement import displacements, local_strain

        k = opts.displacement_reference - 1
        if not 0 <= k < len(frames):
            res.notes.append(L(f"変位の基準フレーム {opts.displacement_reference} がありません (全 {len(frames)} フレーム)",
                               f"reference frame {opts.displacement_reference} for the displacement does not exist "
                               f"({len(frames)} frames)"))
        else:
            ref = frames[k]
            u = displacements(ref.get_positions(), last.get_positions(), cell)
            mag = np.linalg.norm(u, axis=1)
            res.tables["displacement"] = {"reference_frame": opts.displacement_reference,
                                          "mean_ang": float(mag.mean()), "max_ang": float(mag.max()),
                                          "per_atom_ang": [float(x) for x in mag]}
            res.notes.append(L(
                f"変位 (フレーム {opts.displacement_reference} から最後まで、原子 {mag.size} 個): "
                f"平均 {mag.mean():.4f} Å、最大 {mag.max():.4f} Å。原子の並びが同じ順であることが前提です",
                f"displacement (frame {opts.displacement_reference} to the last, {mag.size} atoms): "
                f"mean {mag.mean():.4f} Å, largest {mag.max():.4f} Å; this assumes the atom order is unchanged"))
            if opts.strain_cutoff:
                ls = local_strain(ref.get_positions(), last.get_positions(), opts.strain_cutoff, cell, cell)
                ok = ls.neighbors >= 3
                if not ok.any():
                    res.notes.append(L(f"カットオフ {opts.strain_cutoff:g} Å の中に近傍が 3 個以上ある原子がありません",
                                       f"no atom has three or more neighbors within {opts.strain_cutoff:g} Å"))
                else:
                    res.tables["local_strain"] = {
                        "cutoff_ang": opts.strain_cutoff,
                        "volumetric_mean": float(ls.volumetric[ok].mean()),
                        "shear_mean": float(ls.shear[ok].mean()), "shear_max": float(ls.shear[ok].max()),
                        "residual_mean_ang": float(ls.residual[ok].mean())}
                    res.notes.append(L(
                        f"局所ひずみ (カットオフ {opts.strain_cutoff:g} Å、近傍 3 個以上の {int(ok.sum())} 原子): "
                        f"体積ひずみ 平均 {ls.volumetric[ok].mean():+.5f}、せん断 平均 {ls.shear[ok].mean():.5f}・"
                        f"最大 {ls.shear[ok].max():.5f}、当てはめの残差 平均 {ls.residual[ok].mean():.4f} Å",
                        f"local strain (cutoff {opts.strain_cutoff:g} Å, {int(ok.sum())} atoms with three or more "
                        f"neighbors): volumetric mean {ls.volumetric[ok].mean():+.5f}, shear mean "
                        f"{ls.shear[ok].mean():.5f} and max {ls.shear[ok].max():.5f}, fit residual "
                        f"{ls.residual[ok].mean():.4f} Å"))
    if opts.voronoi:
        try:
            from adit.analysis.voronoi import voronoi as _voronoi
        except ImportError:
            res.notes.append(L("Voronoi には scipy が要ります", "the Voronoi analysis needs scipy"))
        else:
            try:
                vr = _voronoi(last.get_positions(), cell, opts.voronoi_face_threshold)
            except (ValueError, Exception) as ex:
                res.notes.append(L(f"Voronoi を計算できません: {ex}", f"cannot compute the Voronoi analysis: {ex}"))
            else:
                ok = np.isfinite(vr.volume)
                res.tables["voronoi"] = {"closed": int(ok.sum()), "n_atoms": int(vr.volume.size),
                                         "volume_mean_ang3": float(vr.volume[ok].mean()),
                                         "faces_mean": float(vr.faces[ok].mean()),
                                         "face_threshold_ang2": vr.threshold_ang2}
                extra = ""
                if cell is not None and ok.all():
                    extra = L(f"、体積の合計 {vr.volume.sum():.3f} Å³ (セルは {abs(np.linalg.det(cell)):.3f} Å³)",
                              f", total volume {vr.volume.sum():.3f} Å³ (the cell is {abs(np.linalg.det(cell)):.3f} Å³)")
                res.notes.append(L(
                    f"Voronoi (閉じた多面体 {int(ok.sum())} / {vr.volume.size} 原子): 体積 平均 "
                    f"{vr.volume[ok].mean():.4f} Å³、面の数 平均 {vr.faces[ok].mean():.2f} "
                    f"(面積 {vr.threshold_ang2:g} Å² より大きい面だけ)",
                    f"Voronoi (closed cells for {int(ok.sum())} of {vr.volume.size} atoms): mean volume "
                    f"{vr.volume[ok].mean():.4f} Å³, mean face count {vr.faces[ok].mean():.2f} "
                    f"(faces larger than {vr.threshold_ang2:g} Å² only)") + extra)
    if opts.sasa:
        from adit.analysis.sasa import RADII_SETS, radii_for, sasa as _sasa

        name, _, probe_text = opts.sasa.partition(",")
        name = name.strip() or "bondi"
        try:
            probe = float(probe_text) if probe_text.strip() else 1.4
            radii = radii_for(list(last.get_chemical_symbols()), name)
            source = RADII_SETS[name][0] if name in RADII_SETS else name
            sr = _sasa(last.get_positions(), radii, probe, radii_source=source)
        except ValueError as ex:
            res.notes.append(L(f"溶媒接触表面積を計算できません: {ex} 書き方は --sasa bondi,1.4 です。",
                               f"cannot compute the solvent-accessible surface area: {ex} write it as --sasa bondi,1.4."))
        else:
            res.tables["sasa"] = {"total_ang2": sr.total_ang2, "probe_ang": sr.probe_ang,
                                  "radii_set": name, "n_points": sr.n_points,
                                  "per_atom_ang2": [float(x) for x in sr.per_atom_ang2]}
            res.notes.append(L(
                f"溶媒接触表面積 (Shrake–Rupley、プローブ {sr.probe_ang:g} Å、点 {sr.n_points} 個、"
                f"半径は {source}): 合計 {sr.total_ang2:.2f} Å²、表に出ていない原子 "
                f"{int((sr.per_atom_ang2 == 0).sum())} 個。半径とプローブ半径で値が変わります",
                f"solvent-accessible surface area (Shrake-Rupley, probe {sr.probe_ang:g} Å, {sr.n_points} points, "
                f"radii from {source}): total {sr.total_ang2:.2f} Å², {int((sr.per_atom_ang2 == 0).sum())} buried atoms; "
                "the value depends on the radii and the probe"))
    if opts.pca:
        from adit.analysis.frames_pca import kmeans, pca as _pca

        if len(frames) < 3:
            res.notes.append(L("主成分分析には 3 フレーム以上が要ります", "the PCA needs three or more frames"))
            return
        coords = np.array([f.get_positions() for f in frames])
        try:
            pr = _pca(coords, opts.pca)
        except ValueError as ex:
            res.notes.append(L(f"主成分分析を計算できません: {ex}", f"cannot compute the PCA: {ex}"))
            return
        table = {"explained": [float(x) for x in pr.explained],
                 "projections": [[float(v) for v in row] for row in pr.projections]}
        res.notes.append(L(
            "主成分分析 (重ね合わせたあと、成分 " + str(pr.explained.size) + " 個): 寄与率 "
            + ", ".join(f"{v * 100:.1f} %" for v in pr.explained)
            + f" (合わせて {pr.explained.sum() * 100:.1f} %)",
            f"PCA after superposition ({pr.explained.size} components): explained variance "
            + ", ".join(f"{v * 100:.1f}%" for v in pr.explained)
            + f" ({pr.explained.sum() * 100:.1f}% together)"))
        cl = None
        if opts.cluster:
            try:
                cl = kmeans(pr.projections, opts.cluster)
            except ValueError as ex:
                res.notes.append(L(f"分類できません: {ex}", f"cannot cluster the frames: {ex}"))
            else:
                table["cluster_sizes"] = [int(x) for x in cl.sizes]
                table["cluster_representatives"] = [int(x) + 1 for x in cl.representative]
                res.notes.append(L(
                    f"分類 (k-means、群 {int(cl.sizes.size)} 個): 大きさ "
                    + ", ".join(str(int(x)) for x in cl.sizes)
                    + "、中心に近いフレーム " + ", ".join(str(int(i) + 1) for i in cl.representative)
                    + "。群の数は利用者が決めた値で、いくつが正しいかは判定しません",
                    f"clustering (k-means, {int(cl.sizes.size)} groups): sizes "
                    + ", ".join(str(int(x)) for x in cl.sizes)
                    + ", frames closest to each center " + ", ".join(str(int(i) + 1) for i in cl.representative)
                    + "; the number of groups is the one you gave, and ADIT does not judge it"))
        res.tables["pca"] = table
        if pr.projections.shape[1] >= 2:
            fig, ax = plt.subplots(figsize=(4.6, 4.0))
            colors = None if cl is None else cl.labels
            ax.scatter(pr.projections[:, 0], pr.projections[:, 1], c=colors, s=14, cmap="tab10")
            ax.set_xlabel(f"PC1 ({pr.explained[0] * 100:.1f} %)")
            ax.set_ylabel(f"PC2 ({pr.explained[1] * 100:.1f} %)")
            plotstyle.grid(ax)
            save(fig, "pca")


def _add_geometry_series(res: AnalysisResult, data: RunData, frames, opts: AnalysisOptions, dt_fs, save, out: Path) -> None:
    import csv

    import matplotlib.pyplot as plt

    from adit.analysis import geometry_series as G
    from adit.analysis.select import SelectionError, select

    if len(frames) < 2:
        res.notes.append(L("フレームが 1 つしかないので、時系列は出せません",
                           "there is only one frame, so no time series can be computed"))
        return
    x = np.arange(len(frames), dtype=float) * (dt_fs if dt_fs else 1.0)
    xlabel = "time [fs]" if dt_fs else "frame"
    table: dict[str, list] = {}
    fig, ax = plt.subplots(figsize=(6, 3.2))
    drawn = False
    try:
        for text in opts.distances:
            i, j = G.parse_atom_list(text, 2)
            values = G.distance_series(frames, i, j)
            table[f"distance_{text}_A"] = values.tolist()
            ax.plot(x, values, lw=1.2, label=f"d({text}) [Å]"); drawn = True
        for text in opts.angles:
            i, j, k = G.parse_atom_list(text, 3)
            values = G.angle_series(frames, i, j, k)
            table[f"angle_{text}_deg"] = values.tolist()
            ax.plot(x, values, lw=1.2, ls="--", label=f"angle({text}) [deg]"); drawn = True
        for text in opts.dihedrals:
            i, j, k, l = G.parse_atom_list(text, 4)
            values = G.dihedral_series(frames, i, j, k, l)
            table[f"dihedral_{text}_deg"] = values.tolist()
            ax.plot(x, values, lw=1.2, ls=":", label=f"dihedral({text}) [deg]"); drawn = True
    except G.GeometryError as ex:
        plt.close(fig)
        res.notes.append(str(ex))
        return
    if drawn:
        ax.set_xlabel(xlabel); ax.set_ylabel("value"); plotstyle.grid(ax); ax.legend(fontsize=8)
        save(fig, "geometry")
    else:
        plt.close(fig)
    if opts.fes_temperature_k and table:
        _add_fes(res, table, opts, save)
    try:
        idx = select(frames[0], opts.select) if opts.select else None
    except SelectionError as ex:
        res.notes.append(str(ex))
        idx = None
    if opts.rmsd_reference:
        ref = max(0, int(opts.rmsd_reference) - 1)
        values = G.rmsd_series(frames, idx, reference=ref)
        table["rmsd_A"] = values.tolist()
        fig2, ax2 = plt.subplots(figsize=(6, 2.8))
        ax2.plot(x, values, lw=1.2); ax2.set_xlabel(xlabel); ax2.set_ylabel("RMSD [Å]"); plotstyle.grid(ax2)
        save(fig2, "rmsd")
        res.notes.append(L(f"RMSD (基準は {ref + 1} 番目のフレーム、Kabsch で重ね合わせ): 最後の値 {values[-1]:.3f} Å、"
                           f"最大 {values.max():.3f} Å (落ち着いたかどうかは判定していません)",
                           f"RMSD (reference: frame {ref + 1}, after a Kabsch superposition): last {values[-1]:.3f} Å, "
                           f"largest {values.max():.3f} Å. No judgement about equilibration is made"))
    if opts.rmsf:
        got = G.rmsf(frames, idx)
        res.tables["rmsf"] = got
        fig3, ax3 = plt.subplots(figsize=(6, 2.8))
        ax3.bar(range(len(got["rmsf_A"])), got["rmsf_A"], color="tab:green")
        ax3.set_xlabel("atom"); ax3.set_ylabel("RMSF [Å]"); ax3.grid(alpha=0.3, axis="y")
        save(fig3, "rmsf")
        top = int(np.argmax(got["rmsf_A"]))
        res.notes.append(L(f"RMSF (原子ごとの揺らぎ): 最大は {got['symbol'][top]} (番号 {got['index'][top]}) の "
                           f"{got['rmsf_A'][top]:.3f} Å、平均 {float(np.mean(got['rmsf_A'])):.3f} Å",
                           f"RMSF: the largest is {got['rmsf_A'][top]:.3f} Å on {got['symbol'][top]} (index {got['index'][top]}), "
                           f"mean {float(np.mean(got['rmsf_A'])):.3f} Å"))
    if table:
        path = out / "geometry_series.csv"
        with open(path, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow([xlabel] + list(table))
            for row in range(len(x)):
                writer.writerow([f"{x[row]:.6g}"] + [f"{table[k][row]:.6g}" for k in table])
        res.tables["geometry_series"] = {"columns": list(table), "csv": str(path), "n_frames": len(frames),
                                         "note": L("指定した原子の距離 [Å]・角度 [度]・二面角 [度] の時系列です。"
                                                   "距離は周期系なら最小像で測っています。",
                                                   "time series of the distances (Å), angles (deg) and dihedrals (deg) you asked for; "
                                                   "distances use the minimum image in periodic systems")}
        res.notes.append(L(f"幾何の時系列を {path.name} に書きました ({len(table)} 列)",
                           f"the geometry time series was written to {path.name} ({len(table)} columns)"))


def _add_vacf(res: AnalysisResult, pos, syms, frames, opts: AnalysisOptions, dt_fs, save) -> None:
    import matplotlib.pyplot as plt

    from adit.analysis import vacf as V

    if dt_fs is None:
        res.notes.append(L("フレームの間隔が分からないので、速度自己相関は求められません",
                           "the time between frames is unknown, so no velocity autocorrelation can be computed"))
        return
    idx = _selected_atoms(res, frames, syms, opts)
    if not idx:
        return
    stored = []
    for fr in frames:
        v = fr.get_velocities()
        if v is None or not np.any(v):
            stored = []
            break
        stored.append(v[idx])
    try:
        if stored and len(stored) == len(frames):
            velocities, source = np.array(stored, dtype=float), "velocities"
        else:
            velocities, source = V.velocities_from_positions(pos[:, idx, :], dt_fs), "finite_difference"
        result = V.vacf(velocities, dt_fs, source=source)
    except V.VacfError as ex:
        res.notes.append(str(ex))
        return
    fig, ax = plt.subplots(figsize=(6, 3.0))
    ax.plot(result.times_fs, result.vacf, lw=1.2)
    ax.axhline(0.0, color="0.6", lw=0.8)
    ax.set_xlabel("lag time [fs]"); ax.set_ylabel("normalized VACF"); plotstyle.grid(ax)
    save(fig, "vacf")
    fig2, ax2 = plt.subplots(figsize=(6, 3.0))
    ax2.plot(result.freq_cm1, result.spectrum, lw=1.0)
    ax2.set_xlim(0, min(result.freq_cm1[-1], 4000)); ax2.set_xlabel("wavenumber [cm^-1]")
    ax2.set_ylabel("intensity [arb.]"); plotstyle.grid(ax2)
    save(fig2, "vacf_spectrum")
    res.tables["vacf"] = {**result.as_dict(), "species": opts.msd_species, "n_atoms": len(idx)}
    res.notes += V.notes(result)


def _add_conductivity(res: AnalysisResult, spec, data, syms, frames, opts: AnalysisOptions, d_cm2_s) -> None:
    from adit.analysis.conductivity import ConductivityError, nernst_einstein, summary_line

    if d_cm2_s is None:
        res.notes.append(L("拡散係数が出ていないので、イオン伝導度は求められません",
                           "no diffusion coefficient is available, so no ionic conductivity is computed"))
        return
    last = frames[-1] if len(frames) else None
    if last is None or not any(last.pbc):
        res.notes.append(L("非周期の系なのでセルの体積が決まらず、イオン伝導度は求められません "
                           "(伝導度は単位体積あたりの量です)",
                           "the system is not periodic, so there is no cell volume and no ionic conductivity can be computed "
                           "(conductivity is a quantity per unit volume)"))
        return
    volume = float(last.get_volume())
    n_ions = sum(1 for sym in syms if sym == opts.msd_species) if opts.msd_species else len(syms)
    temperature = opts.conductivity_temperature_k
    how = L("指定された温度", "the temperature you gave")
    if not temperature and spec is not None and spec.task.type == "molecular_dynamics":
        temperature, how = float(spec.task.md.temperature_k), L("入力の設定温度", "the target temperature in the input")
    if not temperature and data.temperatures_k:
        temperature, how = float(np.mean(data.temperatures_k)), L("出力の平均温度", "the mean temperature in the output")
    try:
        c = nernst_einstein(d_cm2_s, n_ions, opts.conductivity_charge, volume, temperature)
    except ConductivityError as ex:
        res.notes.append(str(ex))
        return
    res.tables["conductivity"] = {**c.as_dict(), "species": opts.msd_species, "temperature_source": how}
    res.notes.append(summary_line(c) + L(f" (温度は{how})", f" (temperature from {how})"))


def _plot_vanhove_file(res: AnalysisResult, path: Path, save) -> None:
    import matplotlib.pyplot as plt

    try:
        got = json.loads(path.read_text(encoding="utf-8"))
    except ValueError as ex:
        res.notes.append(L(f"{path.name} を読めません: {ex}", f"cannot read {path.name}: {ex}"))
        return
    vh = got.get("vanhove") or {}
    times = np.array(vh.get("times_fs") or [])
    if len(times):
        fig, ax = plt.subplots(figsize=(6, 3.2))
        ax.plot(times, np.array(vh["d_fit_A2_fs"]) * 1e-1, "o-", ms=3, lw=1.0, label="fit of P(r, tau)")
        ax.plot(times, np.array(vh["d_direct_A2_fs"]) * 1e-1, "s--", ms=3, lw=1.0, label="<r^2> / (6 tau)")
        ax.set_xlabel("lag time [fs]"); ax.set_ylabel("D [cm^2/s]"); plotstyle.grid(ax); ax.legend(fontsize=8)
        ax.set_title(f"displacement distribution ({vh.get('displacement', '?')})", fontsize=8, color="0.3")
        save(fig, "vanhove")
        fig2, ax2 = plt.subplots(figsize=(6, 2.6))
        ax2.plot(times, vh.get("alpha2") or [], lw=1.2, color="tab:red"); ax2.axhline(0.0, color="0.6", lw=0.8)
        ax2.set_xlabel("lag time [fs]"); ax2.set_ylabel("alpha2"); plotstyle.grid(ax2)
        ax2.set_title("non-Gaussian parameter (0 = Gaussian)", fontsize=8, color="0.3")
        save(fig2, "vanhove_alpha2")
    per_atom = [v for v in (got.get("d_per_atom_cm2_s") or []) if v is not None]
    if per_atom:
        fig3, ax3 = plt.subplots(figsize=(6, 2.6))
        ax3.bar(range(len(per_atom)), per_atom, color="tab:blue")
        ax3.set_xlabel("atom"); ax3.set_ylabel("D [cm^2/s]"); ax3.grid(alpha=0.3, axis="y")
        save(fig3, "d_per_atom")
    res.tables["vanhove"] = {**vh, "source": str(path), "computed_here": False}
    res.tables["msd_per_atom"] = {"D_cm2_s": got.get("d_per_atom_cm2_s"), "symbol": got.get("symbols"),
                                  "atom_index": got.get("atom_index"), "fit_range_fs": got.get("fit_range_fs"),
                                  "source": str(path), "computed_here": False}
    res.notes.append(L(f"{path.name} を読みました (計算はこの場でしていません)。"
                       + (f"原子ごとの D: {len(per_atom)} 原子、{min(per_atom):.3g}〜{max(per_atom):.3g} cm²/s。" if per_atom else "")
                       + (f"変位の分布からの D: 当てはめ {np.nanmean(np.array(vh['d_fit_A2_fs']) * 1e-1):.3g} cm²/s、"
                          f"⟨r²⟩/(6τ) {np.mean(np.array(vh['d_direct_A2_fs']) * 1e-1):.3g} cm²/s" if len(times) else ""),
                       f"read {path.name} (nothing was computed here)."
                       + (f" D per atom: {len(per_atom)} atoms, {min(per_atom):.3g} to {max(per_atom):.3g} cm^2/s." if per_atom else "")
                       + (f" D from the displacement distribution: fit {np.nanmean(np.array(vh['d_fit_A2_fs']) * 1e-1):.3g} cm^2/s, "
                          f"<r^2>/(6 tau) {np.mean(np.array(vh['d_direct_A2_fs']) * 1e-1):.3g} cm^2/s" if len(times) else "")))
    if vh.get("truncated_from_fs"):
        res.notes.append(L(f"注意: 遅れ時間 {vh['truncated_from_fs']:g} fs から先は、変位の {vh.get('near_cap_fraction', 0):.1%} が"
                           f"最小像で測れる上限 {vh.get('half_width_A', float('nan')):.2f} Å の 9 割を超えています (D は低めに出ます)",
                           f"note: from lag {vh['truncated_from_fs']:g} fs, {vh.get('near_cap_fraction', 0):.1%} of displacements exceed "
                           f"90 % of {vh.get('half_width_A', float('nan')):.2f} Å, the minimum-image limit (D is biased low)"))


def _add_vanhove(res: AnalysisResult, pos, syms, frames, opts: AnalysisOptions, dt_fs, save, run_dir: Path) -> None:
    import matplotlib.pyplot as plt

    from adit.analysis import vanhove as vh

    from adit.analysis import heavy_setup

    done = run_dir / heavy_setup.RESULT_FILE
    if done.is_file():
        _plot_vanhove_file(res, done, save)
        return
    if dt_fs is None:
        res.notes.append(L("フレームの間隔が分からないので、変位の分布からは D を出せません",
                           "the time between frames is unknown, so no D can be obtained from the displacement distribution"))
        return
    last = frames[-1] if len(frames) else None
    cell = np.asarray(last.cell, dtype=float) if last is not None and any(last.pbc) else None
    if cell is None and opts.vanhove_displacement == "mic":
        res.notes.append(L("非周期の系なので、変位の分布は巻き戻した座標で求めます",
                           "the system is not periodic, so the displacement distribution uses unwrapped coordinates"))
    idx = _selected_atoms(res, frames, syms, opts)
    if not idx:
        return
    taus = vh.default_taus(pos.shape[0], count=int(opts.vanhove))
    kind = opts.vanhove_displacement if cell is not None else "unwrapped"
    from adit.analysis import heavy_setup

    guess = heavy_setup.estimate(pos.shape[0], len(idx), len(taus), kind)
    if not opts.vanhove_here and guess["seconds"] > opts.heavy_limit_seconds:
        written = heavy_setup.write_job(run_dir, dt_fs=dt_fs, species=opts.msd_species, taus=len(taus),
                                        displacement=kind)
        res.notes.append(L(
            f"変位の分布は、この場では計算しませんでした (見積もり {guess['seconds']:.0f} 秒。上限 {opts.heavy_limit_seconds:.0f} 秒)。"
            f"実行用のファイルを置きました: {', '.join(p.name for p in written)}。"
            f"人が実行して {heavy_setup.RESULT_FILE} ができたら、もう一度 adit-analyze すると図になります "
            "(その場で計算させるなら --vanhove-here)",
            f"the displacement distribution was not computed here (estimated {guess['seconds']:.0f} s, limit "
            f"{opts.heavy_limit_seconds:.0f} s). Files to run it were written: {', '.join(p.name for p in written)}. "
            f"Run them, and once {heavy_setup.RESULT_FILE} exists, adit-analyze turns it into figures "
            "(use --vanhove-here to compute it in place)"))
        res.tables["vanhove_job"] = {"estimate": guess, "files": [str(p) for p in written],
                                     "result_file": heavy_setup.RESULT_FILE}
        return
    try:
        result = vh.van_hove_self(pos[:, idx, :], cell, taus, dt_fs, displacement=kind)
    except vh.VanHoveError as ex:
        res.notes.append(str(ex))
        return
    fig, ax = plt.subplots(figsize=(6, 3.2))
    ax.plot(result.times, result.d_fit * 1e-1, "o-", ms=3, lw=1.0, label="fit of P(r, tau)")
    ax.plot(result.times, result.d_direct * 1e-1, "s--", ms=3, lw=1.0, label="<r^2> / (6 tau)")
    ax.set_xlabel("lag time [fs]"); ax.set_ylabel("D [cm^2/s]"); plotstyle.grid(ax); ax.legend(fontsize=8)
    ax.set_title(f"displacement distribution ({kind}), {opts.msd_species or 'all atoms'}", fontsize=8, color="0.3")
    save(fig, "vanhove")
    fig2, ax2 = plt.subplots(figsize=(6, 2.6))
    ax2.plot(result.times, result.alpha2, lw=1.2, color="tab:red")
    ax2.axhline(0.0, color="0.6", lw=0.8)
    ax2.set_xlabel("lag time [fs]"); ax2.set_ylabel("alpha2"); plotstyle.grid(ax2)
    ax2.set_title("non-Gaussian parameter (0 = Gaussian)", fontsize=8, color="0.3")
    save(fig2, "vanhove_alpha2")
    res.tables["vanhove"] = {**result.as_dict(), "species": opts.msd_species, "n_atoms": len(idx),
                             "definition": L(
                                 "P(r, τ) は遅れ時間 τ での変位の分布。fit はガウスの形を当てはめた D、direct は ⟨r²⟩/(6τ)。"
                                 "ガウス (フィックの拡散) なら 2 つは一致します。α₂ = 3⟨r⁴⟩/(5⟨r²⟩²) − 1 は 0 からのずれで"
                                 "非ガウス性を測る量です (Rahman 1964)。ADIT はどれが妥当かを判定しません。",
                                 "P(r, tau) is the distribution of displacements at lag tau. 'fit' is D from fitting the Gaussian form, "
                                 "'direct' is <r^2>/(6 tau); they agree for Gaussian (Fickian) motion. alpha2 = 3<r^4>/(5<r^2>^2) - 1 "
                                 "measures the deviation from Gaussian (Rahman 1964). ADIT does not judge which one to use.")}
    res.notes += vh.notes(result)


def _add_xrd(res: AnalysisResult, data: RunData, out: Path, opts: AnalysisOptions, save) -> None:
    import matplotlib.pyplot as plt

    from adit.analysis.xrd import XrdError, powder_pattern, read_measured, scale_to_100, write_csv

    if data.final is None:
        res.notes.append(L("最終構造を読めないので粉末回折は計算しません", "no final structure was read, so no powder pattern is computed"))
        return
    try:
        pattern = powder_pattern(data.final, opts.xrd, tuple(opts.xrd_range))
    except XrdError as ex:
        res.notes.append(str(ex))
        return
    fig, ax = plt.subplots(figsize=(7, 3.4))
    ax.vlines(pattern.two_theta, 0.0, pattern.intensity, color="#0a84ff",
              label="calculated (final structure)")
    measured_note = ""
    if opts.xrd_measured:
        try:
            mx, my = read_measured(opts.xrd_measured)
        except XrdError as ex:
            res.notes.append(str(ex))
        else:
            ax.plot(mx, scale_to_100(my), color="#888", lw=1.0, label="measured (scaled to a maximum of 100)")
            measured_note = L(f" 測定 {Path(opts.xrd_measured).name} ({len(mx)} 点) を重ねました。"
                              "どちらも最大を 100 にそろえてあります。一致・不一致は判定していません。",
                              f" The measurement {Path(opts.xrd_measured).name} ({len(mx)} points) is overlaid; "
                              "both are scaled to a maximum of 100 and no agreement is assessed.")
            ax.legend(fontsize=8)
    ax.set_xlim(*opts.xrd_range)
    ax.set_xlabel("2-theta [degrees]")
    ax.set_ylabel("intensity (max = 100)")
    save(fig, "xrd")
    table = write_csv(out / "xrd.csv", pattern)
    entry = pattern.as_dict()
    entry["csv"] = str(table)
    if opts.xrd_measured:
        entry["measured"] = str(opts.xrd_measured)
    res.tables["xrd"] = entry
    strongest = pattern.two_theta[int(np.argmax(pattern.intensity))]
    res.notes.append(L(f"粉末 X 線回折: ピーク {len(pattern.two_theta)} 本 (最も強いのは 2θ = {strongest:.2f} 度)。{pattern.note}{measured_note}",
                       f"powder XRD: {len(pattern.two_theta)} peaks (the strongest at 2-theta = {strongest:.2f} degrees). {pattern.note}{measured_note}"))


def _add_volumetric(res: AnalysisResult, data: RunData, run_dir: Path, out: Path, opts: AnalysisOptions, save) -> None:
    import matplotlib.pyplot as plt

    from adit.analysis.volumetric import VolumetricError, find_files, plane_average, read_grid, work_function, write_csv

    files = find_files(run_dir)
    if not files:
        res.notes.append(L("体積データ (LOCPOT・CHGCAR・*.cube) がこのディレクトリにありません",
                           "no volumetric data (LOCPOT, CHGCAR, *.cube) in this directory"))
        return
    axis = {"a": 0, "b": 1, "c": 2}.get((opts.plane_average or "c").lower(), 2)
    for path in files:
        try:
            grid = read_grid(path, cube_unit=opts.cube_unit)
        except VolumetricError as ex:
            res.notes.append(str(ex))
            continue
        positions, mean = plane_average(grid, axis)
        name = f"plane_average_{path.name}"
        fig, ax = plt.subplots(figsize=(6, 3.2))
        ax.plot(positions, mean)
        ax.set_xlabel(f"position along {'abc'[axis]} [Ang]")
        ax.set_ylabel(f"{path.name}" + (f" [{grid.unit}]" if grid.unit else ""))
        save(fig, name)
        table = write_csv(out / f"{name}.csv", positions, mean, grid.unit)
        entry = {"file": path.name, "axis": "abc"[axis], "unit": grid.unit, "kind": grid.kind,
                 "points": int(len(mean)), "min": float(mean.min()), "max": float(mean.max()), "csv": str(table)}
        if opts.work_function:
            entry["work_function"] = work_function(grid, data.fermi_ev, axis)
        res.tables.setdefault("volumetric", []).append(entry)
        res.notes.append(L(f"{path.name}: {'abc'[axis]} 軸に垂直な面で平均しました ({len(mean)} 点。表 {table.name})",
                           f"{path.name}: averaged over planes normal to {'abc'[axis]} ({len(mean)} points; table {table.name})"))


def _add_properties(res: AnalysisResult, data: RunData, out: Path) -> None:
    syms = data.final.get_chemical_symbols() if data.final is not None else None
    charges = []
    used: set[str] = set()
    for q in data.charges:
        item = dict(q)
        item["sum"] = float(np.sum(q["values"]))
        if syms is not None and len(syms) == len(q["values"]):
            item["atoms"] = [f"{s}{i + 1}" for i, s in enumerate(syms)]
        item["file"] = _charge_csv(out, item, used)
        charges.append(item)
    if charges:
        res.tables["charges"] = charges
    if data.thermo:
        res.tables["thermochemistry"] = data.thermo
        if data.thermo.get("note"):
            res.notes.append(data.thermo["note"])
    if data.electronic:
        res.tables["electronic"] = dict(data.electronic)


def _add_stage23(res: AnalysisResult, data: RunData, run_dir: Path, out: Path, opts: AnalysisOptions) -> bool:
    from adit.analysis import neb, pdos, phonons, symmetry, thermo, uvvis

    fermi_used = False
    if opts.thermo is not None:
        t = thermo.compute_thermo(data.frequencies_cm1 or [], data.final, data.energies_ev[-1] if data.energies_ev else None, opts.thermo, data.thermo)
        if t["computed"]:
            thermo.rows_csv(t, out / "thermo.csv")
            t["file"] = str(out / "thermo.csv")
        res.tables["thermo_ase"] = t
    nb = neb.analyze_neb(run_dir, out)
    if nb is not None:
        res.tables["neb"] = nb
        if nb.get("figure"):
            res.figures["neb"] = nb["figure"]
    pd = pdos.analyze_pdos(run_dir, data.code, data.fermi_ev, out)
    if pd is not None and ("channels" in pd or opts.pdos):
        res.tables["pdos"] = pd
        if pd.get("figure"):
            res.figures["pdos"] = pd["figure"]
        fermi_used = bool(pd.get("shifted_to_fermi")) and data.code == "espresso"
    elif opts.pdos and pd is None:
        res.notes.append(L("PDOS のファイル (VASP の DOSCAR、QE の projwfc.x の *.pdos_atm#…) がありません",
                           "no PDOS files (VASP DOSCAR, QE projwfc.x *.pdos_atm#...)"))
    if data.code == "orca":
        uv = uvvis.analyze_uvvis(run_dir, out, opts.uvvis_broadening)
        if uv is not None:
            res.tables["uvvis"] = uv
            res.figures["uvvis"] = uv["figure"]
    if opts.spacegroup:
        sg = symmetry.spacegroup(data.final, opts.symprecs)
        if sg is not None:
            if "results" in sg:
                res.tables["spacegroup"] = sg
            else:
                res.notes.append(sg["reason"])
    _add_collections(res, run_dir, out, opts)
    ph = phonons.analyze_phonopy(run_dir, out)
    if ph is not None:
        res.tables["phonopy"] = ph
        for k in ("bands", "dos"):
            if ph.get(f"figure_{k}"):
                res.figures[f"phonon_{k}"] = ph[f"figure_{k}"]
    return fermi_used


def _add_collections(res: AnalysisResult, run_dir: Path, out: Path, opts: AnalysisOptions) -> None:
    from adit.analysis import collections as coll

    kind = coll.detect_collection(run_dir)
    if kind == "phonons":
        res.tables["phonon_set"] = coll.analyze_phonon_set(run_dir, opts.collect)
    elif kind == "elastic":
        t = coll.analyze_elastic(run_dir, out, opts.collect)
        res.tables["elastic"] = t
        if t.get("figure_stress_strain"):
            res.figures["elastic_stress_strain"] = t["figure_stress_strain"]
        res.notes += t.get("moduli_summary") or []
        if t.get("moduli_error"):
            res.notes.append(t["moduli_error"])
    elif kind == "neb_images" and "neb" not in res.tables:
        res.notes.append(L("neb.json はありますが、像のディレクトリ (image_00 から連番) が揃っていないので、反応経路の曲線は出していません",
                           "there is a neb.json, but the image directories (image_00 upwards, consecutive) are not complete, so no reaction-path curve is produced"))
    elif kind == "conformers":
        t = coll.analyze_conformers(run_dir, out)
        res.tables["conformers"] = t
        if t.get("figure"):
            res.figures["conformers"] = t["figure"]
    if opts.compare and (run_dir / coll.COMPARE_FILE).is_file():
        from adit.analysis.compare import CompareError, analyze_compare
        try:
            c = analyze_compare(run_dir)
        except (CompareError, ValueError, OSError) as ex:
            res.notes.append(L(f"{coll.COMPARE_FILE} を読めないので、組にして比べる表は出していません: {ex}",
                               f"no comparison table; {coll.COMPARE_FILE} cannot be read: {ex}"))
            return
        res.tables["compare"] = {"base": c.base, "runs": c.runs, "reactions": c.reactions, "differing": c.differing,
                                 "partial": c.partial, "files": c.files, "figures": c.figures, "notes": c.notes,
                                 "summary": c.summary_text()}
        res.figures.update(c.figures)


XTB_MD_DEFAULTS = {"hmass": "4", "shake": "2", "sccacc": "2.0"}


def _xtb_md_settings(run_dir: Path) -> dict | None:
    p = run_dir / "xtb.inp"
    if not p.is_file():
        return None
    try:
        text = p.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    found: dict[str, str] = {}
    block = ""
    for raw in text.splitlines():
        line = raw.split("#")[0].strip()
        if not line:
            continue
        if line.startswith("$"):
            block = line[1:].strip().lower()
            continue
        if block == "md" and "=" in line:
            k, _, v = line.partition("=")
            k = k.strip().lower()
            if k in XTB_MD_DEFAULTS:
                found[k] = v.strip()
    return {"source": str(p), "found": found, "missing": [k for k in XTB_MD_DEFAULTS if k not in found], "defaults": dict(XTB_MD_DEFAULTS)}


def _load_spec(run_dir: Path):
    try:
        from adit.project import load_project
        return load_project(run_dir)
    except Exception:
        return None


def _target_temperature(spec) -> float | None:
    if spec is None or spec.task.type != "molecular_dynamics":
        return None
    md = spec.task.md
    return float(md.temperature_k) if md.ensemble in ("NVT", "NPT") and md.temperature_k > 0 else None


def _vibration_note(data: RunData, spec) -> str:
    usual = L("振動解析は、最適化済みの構造 (力がほぼ 0) で行うのが普通です", "vibrational analysis is normally done on an optimized structure (forces close to zero)")
    if data.force_max_ev_ang is None:
        return usual + L(" (この出力からは元の構造の力を読めません)", " (the forces on the starting structure cannot be read from this output)")
    tol = spec.task.force_tolerance_ev_per_ang if spec is not None else None
    s = L(f"振動解析の元の構造で、原子にかかる力の最大値は {data.force_max_ev_ang:.3g} eV/Å 以上 ({data.force_source})",
          f"on the starting structure of the vibrational analysis, the largest force on an atom is at least {data.force_max_ev_ang:.3g} eV/Å ({data.force_source})")
    if tol:
        s += L(f"。この計算条件の構造最適化の収束の目安は {tol:.3g} eV/Å", f"; the convergence threshold for geometry optimization in these settings is {tol:.3g} eV/Å")
    return s + L("。", ". ") + usual


_ANALYZE_JA = '''#!/usr/bin/env python
"""adit が生成した解析スクリプト。このディレクトリの計算結果を読み、analysis/ に図 (PNG) と summary.txt を書く。
使い方:  python analyze.py [--rdf] [--msd [元素]] [--dos] [--zdens] [--export] [--skip N] [--stride N] [--rmax R] [--sigma S]
ADIT が入った Python 環境で実行する (pip install adit)。"""
import argparse
from pathlib import Path

from adit.analysis import AnalysisOptions, run_analysis

ap = argparse.ArgumentParser()
ap.add_argument("--rdf", action="store_true", help="動径分布関数と配位数 (元素の全組み合わせ)")
ap.add_argument("--msd", nargs="?", const="", default=None, metavar="元素", help="平均二乗変位と拡散係数 (元素を省くと全原子)")
ap.add_argument("--dos", action="store_true", help="状態密度")
ap.add_argument("--zdens", action="store_true", help="z 方向の密度分布 (周期系)")
ap.add_argument("--export", action="store_true", help="軌跡を analysis/export/ に書き出す (TRAVIS・OVITO・VMD 用)")
ap.add_argument("--skip", type=int, default=0, help="軌跡の先頭を捨てる数 (平衡化)")
ap.add_argument("--stride", type=int, default=1, help="軌跡を N フレームおきに使う (大きな軌跡の間引き)")
ap.add_argument("--rmax", type=float, default=8.0)
ap.add_argument("--sigma", type=float, default=0.1, help="DOS のガウス幅 [eV]")
a = ap.parse_args()
opts = AnalysisOptions(rdf=a.rdf, msd=a.msd is not None, msd_species=(a.msd or None), dos=a.dos, zdens=a.zdens, export=a.export,
                       skip_frames=a.skip, stride=a.stride, rdf_rmax=a.rmax, dos_sigma=a.sigma)
res = run_analysis(Path(__file__).resolve().parent, opts)
print(res.summary_text())
print("図:", ", ".join(res.figures.values()) or "(無し)")
'''

_ANALYZE_EN = '''#!/usr/bin/env python
"""Analysis script generated by ADIT. Reads the results in this directory and writes figures (PNG) and summary.txt to analysis/.
Usage:  python analyze.py [--rdf] [--msd [ELEMENT]] [--dos] [--zdens] [--export] [--skip N] [--stride N] [--rmax R] [--sigma S]
Run it in a Python environment that has ADIT installed (pip install adit)."""
import argparse
from pathlib import Path

from adit.analysis import AnalysisOptions, run_analysis

ap = argparse.ArgumentParser()
ap.add_argument("--rdf", action="store_true", help="radial distribution functions and coordination numbers (all element pairs)")
ap.add_argument("--msd", nargs="?", const="", default=None, metavar="ELEMENT", help="mean square displacement and diffusion coefficient (all atoms if no element)")
ap.add_argument("--dos", action="store_true", help="density of states")
ap.add_argument("--zdens", action="store_true", help="density profile along z (periodic systems)")
ap.add_argument("--export", action="store_true", help="write the trajectory to analysis/export/ (for TRAVIS, OVITO, VMD)")
ap.add_argument("--skip", type=int, default=0, help="number of leading frames to discard (equilibration)")
ap.add_argument("--stride", type=int, default=1, help="use every N-th frame (to thin out large trajectories)")
ap.add_argument("--rmax", type=float, default=8.0)
ap.add_argument("--sigma", type=float, default=0.1, help="Gaussian width of the DOS [eV]")
a = ap.parse_args()
opts = AnalysisOptions(rdf=a.rdf, msd=a.msd is not None, msd_species=(a.msd or None), dos=a.dos, zdens=a.zdens, export=a.export,
                       skip_frames=a.skip, stride=a.stride, rdf_rmax=a.rmax, dos_sigma=a.sigma)
res = run_analysis(Path(__file__).resolve().parent, opts)
print(res.summary_text())
print("figures:", ", ".join(res.figures.values()) or "(none)")
'''


def analyze_script_text() -> str:
    return _ANALYZE_EN if lang.LANGUAGE == "en" else _ANALYZE_JA
