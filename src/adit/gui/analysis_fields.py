
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

from adit.lang import L

ROW_LIMIT = 20

LABELS: dict[str, tuple[str, str]] = {
    "more": ("詳しい条件", "More options"),
    "stride": ("間引き (N フレームおき)", "Stride (every N-th frame)"),
    "msd_fit": ("MSD の当てはめ範囲 [fs]", "MSD fit range [fs]"),
    "msd_axes": ("MSD を取る成分", "Components used for the MSD"),
    "msd_keep_drift": ("重心の流れを落とさない", "Keep the center-of-mass drift"),
    "msd_blocks": ("誤差を出すブロックの数", "Number of blocks for the error"),
    "msd_per_atom": ("原子 1 個ごとの拡散係数も出す", "Also report a diffusion coefficient per atom"),
    "vanhove": ("変位の分布 (van Hove) も見る", "Displacement distribution (van Hove)"),
    "vanhove_taus": ("見る遅れ時間の点数", "Number of lag times"),
    "vanhove_displacement": ("変位の測り方", "How displacements are measured"),
    "zdens": ("z 方向の密度分布", "Density profile along z"),
    "zdens_bin": ("z 密度の区間の幅 [Å]", "z-density bin width [Å]"),
    "stats": ("時系列の統計 (ブロック平均と自己相関時間)", "Time-series statistics"),
    "memory_mb": ("MSD のメモリの上限 [MB]", "Memory limit for the MSD [MB]"),
    "pdos": ("PDOS のファイルが無いときに理由を書く", "Report when there are no PDOS files"),
    "symprec": ("空間群の許容誤差 [Å]", "Space-group tolerances [Å]"),
    "th_model": ("熱化学のモデル", "Thermochemistry model"),
    "th_temps": ("熱化学の温度 [K]", "Thermochemistry temperatures [K]"),
    "th_pressure": ("熱化学の圧力 [Pa]", "Thermochemistry pressure [Pa]"),
    "th_sigma": ("回転の対称数", "Rotational symmetry number"),
    "th_geometry": ("分子の形", "Molecular geometry"),
    "th_spin": ("全スピン S", "Total spin S"),
    "th_imag": ("虚振動の扱い", "Imaginary modes"),
    "th_exclude": ("除く低い振動の本数", "Lowest modes to exclude"),
    "th_qh": ("準調和の下限 [cm⁻¹]", "Quasi-harmonic cutoff [cm⁻¹]"),
    "th_tau": ("準 RRHO の τ [cm⁻¹]", "Quasi-RRHO τ [cm⁻¹]"),
    "uv_shape": ("UV-Vis の広げ方", "UV-Vis broadening"),
    "uv_fwhm": ("UV-Vis の半値全幅 [eV]", "UV-Vis FWHM [eV]"),
    "export": ("TRAVIS・VMD・OVITO 用に書き出す", "Export for TRAVIS, VMD and OVITO"),
    "export_unwrap": ("分子を周期境界でつなぎ直す", "Make molecules whole"),
    "open_export": ("書き出したフォルダを開く", "Open the export folder"),
    "compare": ("組にして比べる…", "Compare runs…"),
    "compare_base": ("基準のディレクトリ", "Base directory"),
    "compare_run": ("比べる", "Compare"),
    "compare_rows": ("比べる組", "Runs to combine"),
    "rx_name": ("反応の名前", "Reaction name"),
    "rx_nu": ("係数 ν", "Coefficient ν"),
    "rx_dir": ("計算結果のディレクトリ", "Run directory"),
    "select": ("原子の選び方", "Atom selection"),
    "rdf_pairs": ("RDF の元素の組", "RDF element pairs"),
    "zdens_axis": ("密度分布を取る軸", "Axis for the density profile"),
    "coordination": ("配位数のカットオフ [Å]", "Coordination cutoff [Å]"),
    "centrosymmetry": ("中心対称性の相手の数", "Centrosymmetry neighbors"),
    "steinhardt": ("Steinhardt q4・q6 のカットオフ [Å]", "Steinhardt q4 and q6 cutoff [Å]"),
    "clusters": ("かたまりのカットオフ [Å]", "Cluster cutoff [Å]"),
    "adf": ("結合角の分布の元素", "Elements for the angle distribution"),
    "adf_cutoff": ("結合角の分布のカットオフ [Å]", "Angle-distribution cutoff [Å]"),
    "sq": ("構造因子 S(q) を出す", "Structure factor S(q)"),
    "hbond": ("水素結合 (距離 [Å], 角度 [度])", "Hydrogen bonds (distance [Å], angle [deg])"),
    "hbond_lifetime": ("水素結合の寿命 (存在の自己相関)", "Hydrogen-bond lifetime (presence autocorrelation)"),
    "hbond_cdf": ("水素結合の距離×角度の分布 (距離の上限 [Å])", "Hydrogen-bond distance-angle map (upper distance [Å])"),
    "rg": ("慣性半径 Rg の時系列", "Radius of gyration over time"),
    "density_grid": ("3 次元の数密度の格子", "3D number-density grid"),
    "voronoi": ("Voronoi の体積と面の数", "Voronoi volumes and face counts"),
    "voronoi_face": ("数える面の面積の下限 [Å²]", "Smallest face area to count [Å²]"),
    "sasa": ("溶媒接触表面積 (半径の表, プローブ [Å])", "Solvent-accessible surface (radii set, probe [Å])"),
    "distances": ("距離の時系列 (原子の組)", "Distance series (atom pairs)"),
    "angles": ("角度の時系列 (原子の組)", "Angle series (atom triples)"),
    "dihedrals": ("二面角の時系列 (原子の組)", "Dihedral series (atom quadruples)"),
    "rmsd_reference": ("RMSD の基準フレーム", "Reference frame for the RMSD"),
    "rmsf": ("原子ごとの揺らぎ (RMSF)", "Fluctuation per atom (RMSF)"),
    "vacf": ("速度自己相関と振動スペクトル", "Velocity autocorrelation (VACF) and spectrum"),
    "conductivity_charge": ("イオンの電荷 (伝導度)", "Ion charge (conductivity)"),
    "conductivity_temperature": ("伝導度の温度 [K]", "Temperature for the conductivity [K]"),
    "displacement": ("変位の基準フレーム", "Reference frame for the displacement"),
    "strain": ("局所ひずみのカットオフ [Å]", "Local-strain cutoff [Å]"),
    "pca": ("主成分分析の成分の数", "Number of principal components"),
    "cluster": ("分類する群の数 (k-means)", "Number of groups (k-means)"),
    "fes": ("自由エネルギー面の温度 [K]", "Temperature for the free-energy surface [K]"),
    "fes_bins": ("自由エネルギー面の区間の数", "Bins for the free-energy surface"),
    "fes_unit": ("自由エネルギーの単位", "Unit for the free energy"),
    "conformer_temperature": ("CREST の配座の重みの温度 [K]", "Temperature for the CREST conformer weights [K]"),
    "bands_window": ("バンド図の縦軸の幅 [eV]", "Band-plot energy window [eV]"),
    "effective_mass_points": ("有効質量に使う k 点の数", "k-points used for the effective mass"),
    "bader": ("Bader の ACF.dat", "Bader ACF.dat"),
    "bader_valence": ("Bader の価電子数", "Valence electrons for Bader"),
    "xrd": ("粉末回折の線源", "Radiation for the powder pattern"),
    "xrd_range": ("粉末回折の 2θ の範囲 [度]", "2-theta range for the powder pattern [deg]"),
    "xrd_measured": ("重ねる測定データのファイル", "Measured pattern to overlay"),
    "viscosity": ("粘度 (Green-Kubo)", "Viscosity (Green-Kubo)"),
    "plane_average": ("面平均を取る軸", "Axis for the plane average"),
    "cube_unit": ("cube の値の単位", "Unit of the cube values"),
    "work_function": ("仕事関数を出す", "Work function"),
    "heavy_limit": ("その場で計算する上限 [秒]", "Limit for computing here [s]"),
    "code": ("計算コードの指定", "Code (when there is no input file)"),
    "freq_scale": ("振動数の補正係数", "Frequency scaling factor"),
    "spectrum_measured": ("重ねる測定したスペクトル", "Measured spectrum to overlay"),
    "plot_colors": ("図の線の色", "Line colors"),
    "plot_ticks": ("目盛りの向き", "Tick direction"),
    "plot_grid": ("目盛り線", "Grid lines"),
    "plot_spines": ("グラフを囲む枠", "Axis frame"),
    "plot_line_width": ("線の太さ [pt]", "Line width [pt]"),
    "plot_font_size": ("図の文字の大きさ [pt]", "Font size in figures [pt]"),
    "plot_dpi": ("図の解像度 [dpi]", "Figure resolution [dpi]"),
    "figure_format": ("図の追加の形式", "Extra figure formats"),
}

THERMO_MODELS = [("", "(計算しない)", "(not computed)"), ("ideal_gas", "理想気体 (IdealGasThermo)", "Ideal gas (IdealGasThermo)"),
                 ("harmonic", "調和振動子 (HarmonicThermo)", "Harmonic (HarmonicThermo)"),
                 ("quasi_harmonic", "準調和 (QuasiHarmonicThermo)", "Quasi-harmonic (QuasiHarmonicThermo)"),
                 ("msrrho", "準 RRHO (MSRRHOThermo)", "Quasi-RRHO (MSRRHOThermo)")]
GEOMETRIES = [("", "(未指定)", "(not set)"), ("linear", "直線 (linear)", "Linear"), ("nonlinear", "非直線 (nonlinear)", "Nonlinear"),
              ("monatomic", "単原子 (monatomic)", "Monatomic")]
IMAGINARY = [("", "(未指定)", "(not set)"), ("ignore", "除いて計算する (ignore)", "Drop them and compute (ignore)"),
             ("stop", "計算しない (stop)", "Do not compute (stop)")]
UV_SHAPES = [("", "(遷移の棒だけ)", "(transition sticks only)"), ("gauss", "ガウス関数", "Gaussian"), ("lorentz", "ローレンツ関数", "Lorentzian")]
AXES = [("c", "c 軸", "c axis"), ("a", "a 軸", "a axis"), ("b", "b 軸", "b axis")]
PLANE_AXES = [("", "(面平均しない)", "(no plane average)"), ("c", "c 軸", "c axis"), ("a", "a 軸", "a axis"), ("b", "b 軸", "b axis")]
FES_UNITS = [("kJ/mol", "kJ/mol", "kJ/mol"), ("kcal/mol", "kcal/mol", "kcal/mol"), ("eV", "eV", "eV")]
CUBE_UNITS = [("", "(書かない)", "(not stated)"), ("ev", "eV", "eV"), ("ry", "Ry", "Ry"), ("hartree", "Hartree", "Hartree")]
TICK_DIRECTIONS = [("", "(既定)", "(default)"), ("in", "内向き", "Inwards"), ("out", "外向き", "Outwards"),
                   ("inout", "内外", "Both")]
GRID_CHOICES = [("", "(既定)", "(default)"), ("both", "縦横", "Both axes"), ("x", "縦線だけ", "Vertical only"),
                ("y", "横線だけ", "Horizontal only"), ("none", "引かない", "None")]
SPINE_CHOICES = [("", "(既定)", "(default)"), ("all", "四方", "All four"), ("left-bottom", "左と下だけ", "Left and bottom"),
                 ("none", "枠なし", "None")]

PLACEHOLDERS: dict[str, tuple[str, str]] = {
    "msd_fit": ("空欄なら最大ずれ時間の 10〜50 %", "Empty = 10–50% of the maximum lag time"),
    "msd_blocks": ("既定 5。0 なら誤差を出しません", "Default 5; 0 disables the error"),
    "vanhove_taus": ("既定 100。重ければ実行用のファイルを置きます", "Default 100; if it is heavy, files to run it are written"),
    "memory_mb": ("空欄なら 1024", "Empty = 1024"),
    "symprec": ("空欄なら 1e-5, 1e-3, 1e-1 を並べる", "Empty = list 1e-5, 1e-3, 1e-1"),
    "th_temps": ("カンマで複数 (例 298.15, 400)", "Comma separated (e.g. 298.15, 400)"),
    "select": ("例: element O and z < 10", "e.g. element O and z < 10"),
    "rdf_pairs": ("例: O-H; O-O。空欄なら全部の組", "e.g. O-H; O-O. Empty = every pair"),
    "coordination": ("空欄なら出しません", "Empty = not computed"),
    "centrosymmetry": ("FCC なら 12、BCC なら 8", "12 for FCC, 8 for BCC"),
    "steinhardt": ("空欄なら出しません", "Empty = not computed"),
    "clusters": ("空欄なら出しません", "Empty = not computed"),
    "adf": ("例: O または O,H", "e.g. O or O,H"),
    "hbond": ("例: 3.5,150。既定値はありません", "e.g. 3.5,150; there is no default"),
    "hbond_cdf": ("例: 4.0。空欄なら出しません", "e.g. 4.0; empty = not computed"),
    "density_grid": ("例: 48,48,48", "e.g. 48,48,48"),
    "voronoi_face": ("空欄なら 0 (全部数える)", "Empty = 0 (count every face)"),
    "sasa": ("例: bondi,1.4", "e.g. bondi,1.4"),
    "distances": ("例: 1,2; 3,4", "e.g. 1,2; 3,4"),
    "angles": ("例: 2,1,3", "e.g. 2,1,3"),
    "dihedrals": ("例: 1,2,3,4", "e.g. 1,2,3,4"),
    "rmsd_reference": ("空欄なら出しません", "Empty = not computed"),
    "conductivity_charge": ("空欄なら出しません", "Empty = not computed"),
    "conductivity_temperature": ("空欄なら MD の温度", "Empty = the MD temperature"),
    "displacement": ("空欄なら出しません", "Empty = not computed"),
    "strain": ("変位と一緒に使います", "Use together with the displacement"),
    "pca": ("空欄なら出しません", "Empty = not computed"),
    "cluster": ("主成分分析が要ります", "The PCA is required"),
    "fes": ("空欄なら出しません", "Empty = not computed"),
    "fes_bins": ("空欄なら 50", "Empty = 50"),
    "conformer_temperature": ("空欄なら重みを出しません", "Empty = no weights"),
    "bands_window": ("空欄なら 10", "Empty = 10"),
    "effective_mass_points": ("空欄なら 5", "Empty = 5"),
    "bader": ("外部の bader が書いたファイル", "written by the external bader program"),
    "bader_valence": ("例: 4,6,1", "e.g. 4,6,1"),
    "xrd": ("例: CuKa。空欄なら出しません", "e.g. CuKa; empty = not computed"),
    "xrd_range": ("空欄なら 5,90", "Empty = 5,90"),
    "xrd_measured": ("1 列目 2θ、2 列目 強度", "first column 2-theta, second column intensity"),
    "heavy_limit": ("空欄なら 60", "Empty = 60"),
    "code": ("空欄なら入力ファイルから判定", "Empty = detected from the input file"),
    "freq_scale": ("例 0.96。空欄なら掛けません", "e.g. 0.96; empty = not applied"),
    "spectrum_measured": ("1 列目 波数 [cm⁻¹]、2 列目 強度", "column 1 wavenumber in cm^-1, column 2 intensity"),
    "plot_colors": ("例: black,#1f77b4,tab:red", "e.g. black,#1f77b4,tab:red"),
    "plot_line_width": ("空欄なら既定", "Empty = default"),
    "plot_font_size": ("空欄なら既定", "Empty = default"),
    "plot_dpi": ("空欄なら既定", "Empty = default"),
    "figure_format": ("例: svg,pdf。空欄なら PNG だけ", "e.g. svg,pdf; empty = PNG only"),
}

CHECKS = ("rdf", "msd", "dos", "zdens", "stats", "pdos", "export_unwrap",
          "msd_keep_drift", "msd_per_atom", "vanhove", "sq", "rg", "voronoi", "rmsf", "vacf",
          "viscosity", "work_function")


def lab(key: str) -> str:
    ja, en = LABELS[key]
    return L(ja, en)


def ph(key: str) -> str:
    ja, en = PLACEHOLDERS[key]
    return L(ja, en)


def choices(items: list[tuple[str, str, str]]) -> list[tuple[str, str]]:
    return [(v, L(ja, en)) for v, ja, en in items]


class FieldError(ValueError):
    pass


def _text(f: dict, key: str) -> str:
    return str(f.get(key) or "").strip()


def _num(f: dict, key: str, kind=float, *, label_key: str | None = None, positive: bool = False, nonneg: bool = False):
    t = _text(f, key)
    if not t:
        return None
    name = lab(label_key or key)
    try:
        v = kind(t)
    except ValueError:
        what = L("整数", "an integer") if kind is int else L("数", "a number")
        raise FieldError(L(f"{name}: {what}として読めません ({t!r})", f"{name}: cannot read {t!r} as {what}")) from None
    if not math.isfinite(v):
        raise FieldError(L(f"{name}: 有限の数にしてください ({t})", f"{name}: must be a finite number ({t})"))
    if positive and not v > 0:
        raise FieldError(L(f"{name}: 正の値にしてください ({t})", f"{name}: must be positive ({t})"))
    if nonneg and v < 0:
        raise FieldError(L(f"{name}: 0 以上にしてください ({t})", f"{name}: must be 0 or more ({t})"))
    return v


def _floats(f: dict, key: str, *, positive: bool = True) -> tuple[float, ...]:
    t = _text(f, key)
    if not t:
        return ()
    out = []
    for part in [p.strip() for p in t.replace("、", ",").split(",") if p.strip()]:
        try:
            v = float(part)
        except ValueError:
            raise FieldError(L(f"{lab(key)}: 数をカンマで区切って書いてください ({t!r})", f"{lab(key)}: write numbers separated by commas ({t!r})")) from None
        if not math.isfinite(v) or (positive and not v > 0):
            raise FieldError(L(f"{lab(key)}: 正の有限の値にしてください ({part})", f"{lab(key)}: values must be positive and finite ({part})"))
        out.append(v)
    return tuple(out)


def _on(f: dict, key: str) -> bool:
    v = f.get(key)
    return v is True or (isinstance(v, str) and v.lower() in ("on", "1", "true", "yes"))


def options_from_fields(f: dict):
    from adit.analysis import AnalysisOptions
    from adit.analysis.symmetry import DEFAULT_SYMPRECS
    from adit.analysis.thermo import ThermoOptions
    from adit.analysis.trajectory import MEMORY_BUDGET_MB

    rmax = _num(f, "rmax", label_key="rmax_label", positive=True) if _text(f, "rmax") else None
    sigma = _num(f, "sigma", label_key="sigma_label", positive=True) if _text(f, "sigma") else None
    skip = _num(f, "skip", int, label_key="skip_label", nonneg=True)
    stride = _num(f, "stride", int, positive=True)
    lo, hi = _num(f, "msd_fit_from", label_key="msd_fit", nonneg=True), _num(f, "msd_fit_to", label_key="msd_fit", positive=True)
    if (lo is None) != (hi is None):
        raise FieldError(L(f"{lab('msd_fit')}: 始めと終わりの両方を入れてください (両方空欄なら最大ずれ時間の 10〜50 %)",
                           f"{lab('msd_fit')}: give both the start and the end (leave both empty for 10–50% of the maximum lag time)"))
    if lo is not None and not lo < hi:
        raise FieldError(L(f"{lab('msd_fit')}: 始め ({lo:g}) を終わり ({hi:g}) より小さくしてください",
                           f"{lab('msd_fit')}: the start ({lo:g}) must be smaller than the end ({hi:g})"))
    zbin = _num(f, "zdens_bin", positive=True)
    mem = _num(f, "memory_mb", positive=True)
    symprecs = _floats(f, "symprec") or DEFAULT_SYMPRECS
    thermo = None
    model = _text(f, "th_model")
    if model:
        thermo = ThermoOptions(model=model, temperatures_k=_floats(f, "th_temps"), pressure_pa=_num(f, "th_pressure"),
                               symmetry_number=_num(f, "th_sigma", int), geometry=_text(f, "th_geometry") or None,
                               spin=_num(f, "th_spin"), imaginary=_text(f, "th_imag") or None,
                               exclude_lowest=_num(f, "th_exclude", int), qh_cutoff_cm1=_num(f, "th_qh"), msrrho_tau_cm1=_num(f, "th_tau"))
    uv = None
    shape, fwhm = _text(f, "uv_shape"), _num(f, "uv_fwhm", positive=True)
    if shape and fwhm is None:
        raise FieldError(L(f"{lab('uv_fwhm')}: 広げ方を選んだときは幅も入れてください", f"{lab('uv_fwhm')}: give a width when a broadening is chosen"))
    if fwhm is not None and not shape:
        raise FieldError(L(f"{lab('uv_shape')}: 幅を入れたときは広げ方も選んでください", f"{lab('uv_shape')}: choose a shape when a width is given"))
    if shape:
        uv = (shape, fwhm)
    return AnalysisOptions(energy=True, temperature=True, bonds=True, rdf=_on(f, "rdf"), rdf_rmax=rmax if rmax is not None else 8.0,
                           msd=_on(f, "msd"), msd_species=_text(f, "msd_species") or None, dos=_on(f, "dos"),
                           dos_sigma=sigma if sigma is not None else 0.1, vibrations=True, skip_frames=skip or 0, stride=stride or 1,
                           msd_fit_fs=(lo, hi) if lo is not None else None, zdens=_on(f, "zdens"), zdens_bin_ang=zbin or 0.2,
                           stats=_on(f, "stats"), export=_on(f, "export"), export_unwrap=_on(f, "export_unwrap"),
                           memory_budget_mb=mem or MEMORY_BUDGET_MB, thermo=thermo, uvvis_broadening=uv, pdos=_on(f, "pdos"), symprecs=symprecs,
                           msd_axes=_text(f, "msd_axes") or "xyz", msd_remove_drift=not _on(f, "msd_keep_drift"),
                           msd_error_blocks=(_num(f, "msd_blocks", int, label_key="msd_blocks", nonneg=True)
                                             if _text(f, "msd_blocks") else 5),
                           msd_per_atom=_on(f, "msd_per_atom"),
                           vanhove=(_num(f, "vanhove_taus", int, label_key="vanhove_taus", positive=True) or 100)
                                   if _on(f, "vanhove") else 0,
                           vanhove_displacement=_text(f, "vanhove_displacement") or "mic",
                           **_extra_options(f))


def _pairs(f: dict, key: str) -> list:
    return [part.strip() for part in _text(f, key).replace("\n", ";").split(";") if part.strip()]


def _element_pairs(f: dict, key: str) -> list:
    # "O-H; O-O" → [("O", "H"), ("O", "O")]
    out = []
    for part in _pairs(f, key):
        a, sep, b = part.partition("-")
        if not sep or not a.strip() or not b.strip():
            raise FieldError(L(f"{lab(key)}: 元素の組は O-H のように書いてください ({part})",
                               f"{lab(key)}: write element pairs like O-H ({part})"))
        out.append((a.strip(), b.strip()))
    return out


def _extra_options(f: dict) -> dict:
    xrd_range = _floats(f, "xrd_range")
    if xrd_range and len(xrd_range) != 2:
        raise FieldError(L(f"{lab('xrd_range')}: 2 つの数をカンマで区切って入れてください",
                           f"{lab('xrd_range')}: give two numbers separated by a comma"))
    fes_bins = _num(f, "fes_bins", int, positive=True)
    return {
        "select": _text(f, "select"),
        "rdf_pairs": _element_pairs(f, "rdf_pairs"),
        "zdens_axis": _text(f, "zdens_axis") or "c",
        "coordination_cutoff": _num(f, "coordination", positive=True) or 0.0,
        "centrosymmetry_neighbors": _num(f, "centrosymmetry", int, positive=True) or 0,
        "steinhardt_cutoff": _num(f, "steinhardt", positive=True) or 0.0,
        "cluster_cutoff": _num(f, "clusters", positive=True) or 0.0,
        "adf": _text(f, "adf") or None,
        "adf_cutoff": _num(f, "adf_cutoff", positive=True) or 0.0,
        "structure_factor": _on(f, "sq"),
        "hbond": _text(f, "hbond"),
        "hbond_lifetime": _on(f, "hbond_lifetime"),
        "hbond_cdf": _num(f, "hbond_cdf", positive=True) or 0.0,
        "radius_of_gyration": _on(f, "rg"),
        "density_grid": _text(f, "density_grid"),
        "voronoi": _on(f, "voronoi"),
        "voronoi_face_threshold": _num(f, "voronoi_face", nonneg=True) or 0.0,
        "sasa": _text(f, "sasa"),
        "distances": _pairs(f, "distances"),
        "angles": _pairs(f, "angles"),
        "dihedrals": _pairs(f, "dihedrals"),
        "rmsd_reference": _num(f, "rmsd_reference", int, positive=True) or 0,
        "rmsf": _on(f, "rmsf"),
        "vacf": _on(f, "vacf"),
        "conductivity_charge": _num(f, "conductivity_charge") or 0.0,
        "conductivity_temperature_k": _num(f, "conductivity_temperature", positive=True) or 0.0,
        "displacement_reference": _num(f, "displacement", int, positive=True) or 0,
        "strain_cutoff": _num(f, "strain", positive=True) or 0.0,
        "pca": _num(f, "pca", int, positive=True) or 0,
        "cluster": _num(f, "cluster", int, positive=True) or 0,
        "fes_temperature_k": _num(f, "fes", positive=True) or 0.0,
        "fes_bins": fes_bins or 50,
        "fes_unit": _text(f, "fes_unit") or "kJ/mol",
        "conformer_temperature_k": _num(f, "conformer_temperature", positive=True) or 0.0,
        "bands_window_ev": _num(f, "bands_window", positive=True) or 10.0,
        "effective_mass_points": _num(f, "effective_mass_points", int, positive=True) or 5,
        "bader": _text(f, "bader"),
        "bader_valence": _text(f, "bader_valence"),
        "xrd": _text(f, "xrd"),
        "xrd_range": tuple(xrd_range) if xrd_range else (5.0, 90.0),
        "xrd_measured": _text(f, "xrd_measured") or None,
        "viscosity": _on(f, "viscosity"),
        "plane_average": _text(f, "plane_average") or None,
        "cube_unit": _text(f, "cube_unit"),
        "work_function": _on(f, "work_function"),
        "heavy_limit_seconds": _num(f, "heavy_limit", positive=True) or 60.0,
        "code": _text(f, "code"),
        "freq_scale": _num(f, "freq_scale", positive=True) or 0.0,
        "spectrum_measured": _text(f, "spectrum_measured"),
        "plot_colors": _text(f, "plot_colors"),
        "plot_ticks": _text(f, "plot_ticks"),
        "plot_grid": _text(f, "plot_grid"),
        "plot_spines": _text(f, "plot_spines"),
        "plot_line_width": _num(f, "plot_line_width", nonneg=True) or 0.0,
        "plot_font_size": _num(f, "plot_font_size", nonneg=True) or 0.0,
        "plot_dpi": _num(f, "plot_dpi", int, nonneg=True) or 0,
        "figure_format": _text(f, "figure_format"),
    }


def fields_from_options(o) -> dict[str, str]:
    from adit.analysis.symmetry import DEFAULT_SYMPRECS
    from adit.analysis.trajectory import MEMORY_BUDGET_MB

    def g(x) -> str:
        return "" if x is None else f"{x:g}" if isinstance(x, float) else str(x)

    f = {"rmax": g(o.rdf_rmax), "msd_species": o.msd_species or "", "sigma": g(o.dos_sigma), "skip": str(o.skip_frames),
         "stride": str(o.stride), "msd_fit_from": g(o.msd_fit_fs[0]) if o.msd_fit_fs else "", "msd_fit_to": g(o.msd_fit_fs[1]) if o.msd_fit_fs else "",
         "zdens_bin": g(o.zdens_bin_ang), "memory_mb": "" if o.memory_budget_mb == MEMORY_BUDGET_MB else g(o.memory_budget_mb),
         "symprec": "" if tuple(o.symprecs) == tuple(DEFAULT_SYMPRECS) else ", ".join(f"{x:g}" for x in o.symprecs),
         "uv_shape": o.uvvis_broadening[0] if o.uvvis_broadening else "", "uv_fwhm": g(o.uvvis_broadening[1]) if o.uvvis_broadening else "",
         "msd_axes": o.msd_axes, "msd_blocks": str(o.msd_error_blocks),
         "vanhove_taus": str(o.vanhove or 100), "vanhove_displacement": o.vanhove_displacement}
    for key, on in (("msd_keep_drift", not o.msd_remove_drift), ("msd_per_atom", o.msd_per_atom), ("vanhove", bool(o.vanhove))):
        f[key] = "on" if on else ""
    t = o.thermo
    f.update(th_model=(t.model or "") if t else "", th_temps=", ".join(f"{x:g}" for x in t.temperatures_k) if t else "",
             th_pressure=g(t.pressure_pa) if t else "", th_sigma=g(t.symmetry_number) if t else "", th_geometry=(t.geometry or "") if t else "",
             th_spin=g(t.spin) if t else "", th_imag=(t.imaginary or "") if t else "", th_exclude=g(t.exclude_lowest) if t else "",
             th_qh=g(t.qh_cutoff_cm1) if t else "", th_tau=g(t.msrrho_tau_cm1) if t else "")
    for k in CHECKS:
        f[k] = "on" if getattr(o, k, False) else ""
    f.update(_extra_fields(o))
    return f


def _extra_fields(o) -> dict[str, str]:
    def g(x) -> str:
        return "" if x in (None, 0, 0.0) else f"{x:g}" if isinstance(x, float) else str(x)

    f = {"select": o.select, "rdf_pairs": "; ".join(f"{a}-{b}" for a, b in o.rdf_pairs),
         "zdens_axis": o.zdens_axis, "coordination": g(o.coordination_cutoff),
         "centrosymmetry": g(o.centrosymmetry_neighbors), "steinhardt": g(o.steinhardt_cutoff),
         "clusters": g(o.cluster_cutoff), "adf": o.adf or "", "adf_cutoff": g(o.adf_cutoff),
         "hbond": o.hbond, "hbond_cdf": g(o.hbond_cdf), "density_grid": o.density_grid, "voronoi_face": g(o.voronoi_face_threshold),
         "sasa": o.sasa, "distances": "; ".join(o.distances), "angles": "; ".join(o.angles),
         "dihedrals": "; ".join(o.dihedrals), "rmsd_reference": g(o.rmsd_reference),
         "conductivity_charge": g(o.conductivity_charge),
         "conductivity_temperature": g(o.conductivity_temperature_k),
         "displacement": g(o.displacement_reference), "strain": g(o.strain_cutoff),
         "pca": g(o.pca), "cluster": g(o.cluster), "fes": g(o.fes_temperature_k),
         "fes_bins": "" if o.fes_bins == 50 else str(o.fes_bins), "fes_unit": o.fes_unit,
         "conformer_temperature": g(o.conformer_temperature_k),
         "bands_window": "" if o.bands_window_ev == 10.0 else g(o.bands_window_ev),
         "effective_mass_points": "" if o.effective_mass_points == 5 else str(o.effective_mass_points),
         "bader": o.bader, "bader_valence": o.bader_valence, "xrd": o.xrd,
         "xrd_range": "" if tuple(o.xrd_range) == (5.0, 90.0) else ",".join(f"{x:g}" for x in o.xrd_range),
         "xrd_measured": str(o.xrd_measured) if o.xrd_measured else "",
         "plane_average": o.plane_average or "", "cube_unit": o.cube_unit,
         "heavy_limit": "" if o.heavy_limit_seconds == 60.0 else g(o.heavy_limit_seconds),
         "code": o.code, "freq_scale": g(o.freq_scale), "spectrum_measured": o.spectrum_measured,
         "plot_colors": o.plot_colors, "plot_ticks": o.plot_ticks, "plot_grid": o.plot_grid,
         "plot_spines": o.plot_spines, "plot_line_width": g(o.plot_line_width),
         "plot_font_size": g(o.plot_font_size), "plot_dpi": g(o.plot_dpi),
         "figure_format": o.figure_format}
    for key, on in (("sq", o.structure_factor), ("rg", o.radius_of_gyration), ("voronoi", o.voronoi), ("hbond_lifetime", o.hbond_lifetime),
                    ("rmsf", o.rmsf), ("vacf", o.vacf), ("viscosity", o.viscosity),
                    ("work_function", o.work_function)):
        f[key] = "on" if on else ""
    return f


_DETAIL_DEFAULTS = {"stride": "1", "msd_fit_from": "", "msd_fit_to": "", "zdens": "", "zdens_bin": "0.2", "stats": "on", "memory_mb": "",
                    "pdos": "", "symprec": "", "th_model": "", "uv_shape": "", "uv_fwhm": ""}


def details_changed(f: dict) -> bool:
    return any(str(f.get(k, v) or "").strip() != v for k, v in _DETAIL_DEFAULTS.items())


MSD_AXES = [("xyz", "3 次元 (xyz、MSD = 6Dt)", "3D (xyz, MSD = 6Dt)"),
            ("xy", "面内 2 次元 (xy、MSD = 4Dt)", "in-plane 2D (xy, MSD = 4Dt)"),
            ("z", "1 次元 (z、MSD = 2Dt)", "1D (z, MSD = 2Dt)")]
VANHOVE_DISPLACEMENTS = [("mic", "最小像 (セルの幅の半分で頭打ち)", "minimum image (capped at half the cell width)"),
                         ("unwrapped", "巻き戻した座標の差", "difference of unwrapped coordinates")]

LABELS.update({"rmax_label": ("r の最大値 [Å]", "r max [Å]"), "sigma_label": ("ガウス幅 [eV]", "Gaussian width [eV]"),
               "skip_label": ("平衡化として捨てるフレーム数", "Frames to skip as equilibration")})


def reactions_from_rows(rows: list[tuple[str, str, str]]):
    from adit.analysis.compare import CompareError, Reaction

    order: list[str] = []
    terms: dict[str, list[tuple[float, str]]] = {}
    name = ""
    for i, (n, nu_s, d) in enumerate(rows, start=1):
        n, nu_s, d = (n or "").strip(), (nu_s or "").strip(), (d or "").strip()
        if not (n or nu_s or d):
            continue
        name = n or name or "r1"
        if not d:
            raise FieldError(L(f"{lab('compare_rows')} の {i} 行目: ディレクトリが空です", f"{lab('compare_rows')} row {i}: the directory is empty"))
        try:
            nu = float(nu_s)
        except ValueError:
            raise FieldError(L(f"{lab('compare_rows')} の {i} 行目: 係数を数として読めません ({nu_s!r})", f"{lab('compare_rows')} row {i}: cannot read the coefficient {nu_s!r} as a number")) from None
        if nu == 0:
            raise FieldError(L(f"{lab('compare_rows')} の {i} 行目: 係数が 0 です", f"{lab('compare_rows')} row {i}: the coefficient is 0"))
        if name not in terms:
            order.append(name); terms[name] = []
        terms[name].append((nu, d))
    if not order:
        raise CompareError(L("比べる組が入っていません (係数とディレクトリを 1 行以上)", "no runs given (at least one row with a coefficient and a directory)"))
    return [Reaction(n, terms[n]) for n in order]


def rows_from_reactions(reactions) -> list[tuple[str, str, str]]:
    out = []
    for r in reactions:
        for k, (nu, d) in enumerate(r.terms):
            out.append((r.name if k == 0 else "", f"{nu:g}", d))
    return out


@dataclass
class Section:
    key: str
    title: str
    columns: list[str]
    rows: list[list[str]]
    numeric: list[bool] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    total_rows: int = 0
    full_file: str = ""

    def more_text(self) -> str:
        if self.total_rows <= len(self.rows):
            return ""
        return L(f"先頭 {len(self.rows)} 行を表示 (全 {self.total_rows} 行)。全体は {self.full_file}",
                 f"first {len(self.rows)} of {self.total_rows} rows shown; all rows are in {self.full_file}")


def _cut(sec: Section, full_file: str) -> Section:
    sec.total_rows = len(sec.rows)
    if len(sec.rows) > ROW_LIMIT:
        sec.rows = sec.rows[:ROW_LIMIT]
        sec.full_file = full_file
    if not sec.numeric:
        sec.numeric = [False] * len(sec.columns)
    return sec


def _g(x, fmt: str = ".6g") -> str:
    if x is None:
        return "-"
    try:
        return format(float(x), fmt)
    except (TypeError, ValueError):
        return str(x)


def result_sections(res) -> list[Section]:
    t = res.tables
    summary_json = str(Path(res.run_dir) / "analysis" / "summary.json")
    out: list[Section] = []
    for k, q in enumerate(t.get("charges", [])):
        names = q.get("atoms") or [str(i + 1) for i in range(len(q["values"]))]
        rows = [[n, f"{v:+.4f}"] for n, v in zip(names, q["values"])]
        sec = _cut(Section(f"charges{k}", L(f"原子の電荷 ({q['definition']})", f"Atomic charges ({q['definition']})"),
                           [L("原子", "Atom"), L("電荷 [e]", "Charge [e]")], rows, [False, True],
                           [L(f"定義: {q['definition']}。{q.get('note') or ''}", f"definition: {q['definition']}. {q.get('note') or ''}").rstrip(". 。"),
                            L(f"出典: {q['source']} (ファイル:行)", f"source: {q['source']} (file:line)"),
                            L(f"合計 {q['sum']:+.4f} e", f"sum {q['sum']:+.4f} e")]), summary_json + f" (tables.charges[{k}])")
        out.append(sec)
    if "thermochemistry" in t:
        th = t["thermochemistry"]
        rows = [[it["label"], f"{it['value_eh']:.6f}", f"{it['value_ev']:.4f}", it.get("line", "")] for it in th["items"]
                if it["key"] != "total_free_energy_box"]
        temp = L(f"、{th['temperature_k']:g} K", f", {th['temperature_k']:g} K") if th.get("temperature_k") is not None else ""
        out.append(_cut(Section("thermochemistry", L(f"熱化学 ({th['code']} が出した値{temp})", f"Thermochemistry (values from {th['code']}{temp})"),
                                [L("項目", "Item"), L("値 [Eh]", "Value [Eh]"), L("値 [eV]", "Value [eV]"), L("出典 (ファイル:行)", "Source (file:line)")],
                                rows, [False, True, True, False], [L(f"方法: {th.get('method', '')}", f"method: {th.get('method', '')}")]), summary_json))
    if "thermo_ase" in t:
        out.append(_thermo_ase(t["thermo_ase"]))
    if "electronic" in t:
        el = t["electronic"]
        rows = []
        if el.get("homo_lumo_gap_ev") is not None:
            rows.append([L("HOMO-LUMO ギャップ", "HOMO-LUMO gap"), f"{el['homo_lumo_gap_ev']:.4f} eV", el.get("gap_source", "")])
        if el.get("dipole_norm_debye") is not None:
            x, y, z = el["dipole_debye"]
            rows.append([L("双極子モーメント", "Dipole moment"), f"{el['dipole_norm_debye']:.4f} D (x, y, z = {x:.4f}, {y:.4f}, {z:.4f})",
                         el.get("dipole_source", "")])
        if rows:
            out.append(_cut(Section("electronic", L("電子状態", "Electronic structure"), [L("量", "Quantity"), L("値", "Value"), L("出典 (ファイル:行)", "Source (file:line)")],
                                    rows), summary_json))
    if "timeseries_stats" in t:
        rows = []
        for key, s in t["timeseries_stats"].items():
            tau = "-" if s.get("tau_int_samples") is None else (f"{s['tau_int_samples']:.3g}" + (f" ({s['tau_int_fs']:.3g} fs)" if s.get("tau_int_fs") else ""))
            blocks = ", ".join(f"{b['block_size']}:{b['sem']:.2g}" for b in s["blocks"])
            rows.append([_series_name(key), f"{s['mean']:.6g} {s['unit']}", f"{s['std']:.3g}", str(s["n"]), tau, _g(s.get("sem_acf"), ".3g"), blocks])
        out.append(_cut(Section("timeseries_stats", L("時系列の統計", "Time-series statistics"),
                                [L("量", "Quantity"), L("平均", "Mean"), L("標準偏差", "Std"), L("点の数", "Points"),
                                 L("積分自己相関時間 [点]", "Integrated autocorrelation time [samples]"),
                                 L("平均値の標準誤差 (自己相関時間から)", "Std. error of the mean (from the autocorrelation time)"),
                                 L("ブロック平均の誤差 (束ねた点の数:誤差)", "Block-averaging error (block size:error)")],
                                rows, [False, True, True, True, True, True, False], [t.get("timeseries_stats_definition", "")]), summary_json))
    if "trajectory" in t:
        x = t["trajectory"]
        rows = [[L("読んだファイル", "File read"), str(x["source"])], [L("フレーム数 (ファイル全体)", "Frames in the file"), str(x["n_frames_total"])],
                [L("先頭で捨てたフレーム", "Leading frames skipped"), str(x["skip"])], [L("間引き (N フレームおき)", "Stride (every N-th frame)"), str(x["stride"])],
                [L("使ったフレーム数", "Frames used"), str(x["n_frames_used"])],
                [L("1 フレームの時間 (元 → 間引き後)", "Time per frame (original → after stride)"),
                 f"{_g(x.get('dt_frame_fs'), 'g')} fs → {_g(x.get('dt_used_fs'), 'g')} fs" if x.get("dt_frame_fs") else L("不明", "unknown")]]
        m = t.get("msd_memory")
        if m:
            rows.append([L("MSD の座標のメモリ (見積もり / 上限)", "Memory for MSD coordinates (estimate / limit)"), f"{m['estimated_mb']:.3g} MB / {m['budget_mb']:g} MB"])
        out.append(_cut(Section("trajectory", L("軌跡", "Trajectory"), [L("項目", "Item"), L("値", "Value")], rows), summary_json))
    if "msd" in t:
        m = t["msd"]
        rows = [[m["species"] or L("全原子", "all atoms"), "-", _g(m.get("D_cm2_s"), ".3e"), f"{m['last_A2']:.4g}"]]
        for el, v in m.get("by_element", {}).items():
            rows.append([el, str(v["n_atoms"]), _g(v.get("D_cm2_s"), ".3e"), f"{v['last_A2']:.4g}"])
        notes = []
        if m.get("fit_range_fs"):
            a, b = m["fit_range_fs"]
            formula = m.get("formula") or f"MSD = {2 * m.get('dimension', 3)} D t + c"
            setting = (L("利用者の指定", "set by the user") if m.get("fit_range_user") else
                       L("既定: 最大ずれ時間の 10〜50 %", "default: 10–50% of the maximum lag time"))
            notes.append(L(f"当てはめ範囲 {a:g}〜{b:g} fs ({setting})。使った式: {formula}。複数の時間原点で平均した MSD",
                           f"fit range {a:g}-{b:g} fs ({setting}); formula used: {formula}; MSD averaged over time origins"))
        else:
            notes.append(L("1 フレームの時間が分からないので拡散係数は出していません", "no diffusion coefficient (time per frame unknown)"))
        out.append(_cut(Section("msd", L("拡散係数 (MSD から)", "Diffusion coefficients (from the MSD)"),
                                [L("対象", "Atoms"), L("原子数", "Count"), "D [cm²/s]", L("最後の MSD [Å²]", "Last MSD [Å²]")], rows,
                                [False, True, True, True], notes), summary_json))
    if "spacegroup" in t:
        sg = t["spacegroup"]
        rows = [[f"{r['symprec_A']:g}", r["international"] or "?", "-" if r["number"] is None else str(r["number"]), r["pointgroup"] or "?"] for r in sg["results"]]
        out.append(_cut(Section("spacegroup", L("空間群", "Space group"), [L("許容誤差 [Å]", "Tolerance [Å]"), L("記号", "Symbol"), L("番号", "Number"), L("点群", "Point group")],
                                rows, [True, False, True, False], [L(f"{sg['structure']}、spglib {sg['spglib']}。どれを採るかは判断しません",
                                                                     f"{sg['structure']}, spglib {sg['spglib']}; no value is chosen for you")]), summary_json))
    if "uvvis" in t:
        uv = t["uvvis"]
        rows = [[r["label"], f"{r['energy_ev']:.4f}", f"{r['wavelength_nm']:.1f}", f"{r['fosc']:.5f}"] for r in uv["transitions"]]
        out.append(_cut(Section("uvvis", L("UV-Vis の遷移", "UV-Vis transitions"), [L("遷移", "Transition"), "E [eV]", "λ [nm]", L("振動子強度 f", "Oscillator strength f")],
                                rows, [False, True, True, True], [L(f"出典: {uv['source']}", f"source: {uv['source']}")] + list(uv.get("reasons", []))),
                        uv.get("file_transitions", summary_json)))
    if "hbond_lifetime" in t:
        h = t["hbond_lifetime"]
        u = h["unit"]

        def v(x) -> str:
            return "-" if x is None else f"{x:.4g}"

        rows = [[L("intermittent (切れて戻っても数える)", "intermittent (re-formed bonds count)"), v(h["lifetime_intermittent"]["integral"]),
                 v(h["lifetime_intermittent"]["one_over_e"]), v(h["lifetime_intermittent"]["last_value"])],
                [L("continuous (切れたら終わり)", "continuous (ends at the first break)"), v(h["lifetime_continuous"]["integral"]),
                 v(h["lifetime_continuous"]["one_over_e"]), v(h["lifetime_continuous"]["last_value"])]]
        out.append(_cut(Section("hbond_lifetime", L("水素結合の寿命", "Hydrogen-bond lifetime"),
                                [L("定義", "Definition"), L(f"C(τ) の積分 [{u}]", f"Integral of C(τ) [{u}]"), L(f"1/e の時間 [{u}]", f"1/e time [{u}]"),
                                 L("τ_max での C", "C at τ_max")], rows, [False, True, True, True],
                                [L(f"τ_max = {h['tau_max']} {u}、{h['n_frames']} フレーム、{h['n_pairs']} 組。", f"tau_max = {h['tau_max']} {u}, {h['n_frames']} frames, {h['n_pairs']} pairs. ")
                                 + h["definition"], L(f"出典: {h['source']}", f"source: {h['source']}")]), h.get("file", summary_json)))
    if "crest_conformers" in t:
        c = t["crest_conformers"]
        weighted = c.get("temperature_k") is not None and all("weight" in r for r in c.get("conformers", []))
        rows = [[str(r["index"]), _g(r.get("relative_kcal_mol"), ".3f"), _g(r.get("relative_kj_mol"), ".2f"), str(r["degeneracy"]),
                 _g(r.get("rmsd_heavy_A"), ".3f")] + ([_g(r.get("weight"), ".4f")] if weighted else []) for r in c.get("conformers", [])]
        cols = [L("番号", "Index"), "E−E(lowest) [kcal/mol]", "E−E(lowest) [kJ/mol]", L("縮退度", "Degeneracy"), L("重原子 RMSD [Å]", "Heavy-atom RMSD [Å]")]
        if weighted:
            cols.append(L(f"重み ({c['temperature_k']:g} K)", f"Weight ({c['temperature_k']:g} K)"))
        notes = [c.get("rmsd_note", "")] + list(c.get("reasons", []))
        if weighted:
            notes.append(c.get("weight_formula", ""))
        out.append(_cut(Section("crest_conformers", L("CREST の配座", "CREST conformers"), cols, rows, [True] * len(cols),
                                [n for n in notes if n]), (c.get("files") or {}).get("table", summary_json)))
    if "export" in t:
        x = t["export"]
        rows = [[L("書き出し先", "Folder"), x["dir"]], [L("フレーム数", "Frames"), str(x["n_frames"])],
                [L("分子をつなぎ直した", "Molecules made whole"), L("はい", "yes") if x.get("unwrap_molecules") else L("いいえ", "no")]]
        if "n_molecules" in x:
            rows.append([L("分子の数 (最大の原子数)", "Molecules (largest, atoms)"), f"{x['n_molecules']} ({x['largest_molecule_atoms']})"])
        out.append(_cut(Section("export", L("書き出し", "Export"), [L("項目", "Item"), L("値", "Value")], rows), summary_json))
    return out


def _series_name(key: str) -> str:
    return {"temperature": L("温度", "temperature"), "energy": L("エネルギー", "energy"), "density": L("密度", "density"),
            "pressure": L("圧力", "pressure"), "conserved": L("保存量", "conserved quantity"), "density_log": L("密度 (ログの値)", "density (log value)")}.get(key, key)


def _thermo_ase(th: dict) -> Section:
    title = L("熱化学 (ASE で計算)", "Thermochemistry (computed with ASE)")
    if not th.get("computed"):
        return _cut(Section("thermo_ase", title, [], [], [], list(th.get("reasons", []))), "")
    gas = any("g_corr_ev" in r for r in th["rows"])
    cols = ["T [K]"] + (["P [Pa]"] if gas else []) + ["ZPE [eV]", L("H−E [eV]", "H−E [eV]") if gas else L("U−E [eV]", "U−E [eV]"),
                                                     "S [J/(mol K)]", "G−E [eV]" if gas else "F−E [eV]", "G [eV]" if gas else "F [eV]"]
    rows = []
    for r in th["rows"]:
        tot = r.get("g_total_ev" if gas else "f_total_ev")
        rows.append([f"{r['T_K']:g}"] + ([f"{r['P_Pa']:g}"] if gas else []) +
                    [f"{r['zpe_ev']:.4f}", f"{r['h_corr_ev' if gas else 'u_corr_ev']:.4f}", f"{r['s_j_mol_k']:.2f}",
                     f"{r['g_corr_ev' if gas else 'f_corr_ev']:.4f}", _g(tot, ".6f")])
    notes = [f"{th['ase_class']} ({th['source']}): {th['includes']}",
             L(f"使ったモード {th['n_modes_used']} 本 ({th['mode_rule']})。虚振動 {th['n_imaginary_all']} 本 (使うモードの中 {th['n_imaginary_used']} 本)",
               f"modes used {th['n_modes_used']} ({th['mode_rule']}); imaginary {th['n_imaginary_all']} ({th['n_imaginary_used']} among the modes used)"),
             th.get("note", "")]
    if th.get("file"):
        notes.append(L(f"表の全体: {th['file']}", f"full table: {th['file']}"))
    return _cut(Section("thermo_ase", title, cols, rows, [True] * len(cols), [n for n in notes if n]), th.get("file", ""))


def compare_sections(cres) -> list[Section]:
    out = []
    rows = []
    for x in cres.reactions:
        de = x["delta_e_ev"]
        bal = L("釣り合う", "balanced") if x["balanced"] else (L("釣り合わない: ", "not balanced: ") + str(x["imbalance"]))
        th = x.get("thermo") or []

        def tcol(k: str, fmt: str, first: bool = False) -> str:
            if th:
                return "; ".join(_g(t.get(k), fmt) for t in th)
            return (x.get("thermo_note") or "-") if first else "-"

        rows.append([x["name"], x["terms"], _g(de, "+.6f"), _g(x.get("delta_e_kj_mol"), "+.3f"), _g(x.get("delta_e_kcal_mol"), "+.3f"),
                     _g(x.get("delta_g_code_ev"), "+.6f"), tcol("T_K", "g", True), tcol("delta_h_kj_mol", "+.3f"), tcol("delta_s_j_mol_k", "+.3f"),
                     tcol("delta_g_kj_mol", "+.3f"), bal, f"{x['n_differing']}: " + " ".join(x["differing"]) if x["differing"] else "0"])
    out.append(_cut(Section("compare_reactions", L("反応ごと (ΔE = ΣνE)", "Per reaction (ΔE = ΣνE)"),
                            [L("名前", "Name"), L("組 (ν·ディレクトリ)", "Terms (ν·directory)"), "ΔE [eV]", "ΔE [kJ/mol]", "ΔE [kcal/mol]",
                             L("コードの G の差 ΣνG [eV]", "Difference of code G, ΣνG [eV]"), L("熱化学の T [K]", "Thermo T [K]"),
                             "ΔH [kJ/mol]", "ΔS [J/(mol K)]", L("ΔG (F) [kJ/mol]", "ΔG (F) [kJ/mol]"), L("組成の釣り合い", "Composition balance"),
                             L("条件が違う項目", "Settings that differ")],
                            rows, [False, False, True, True, True, True, True, True, True, True, False, False],
                            [L("ν は生成物が正、反応物が負。組成の釣り合いは Σν·(元素ごとの原子数)", "ν is positive for products and negative for reactants; balance is Σν·(atoms per element)"),
                             L("ΔH・ΔS・ΔG は各計算の analysis/thermo.csv (ASE の熱化学) の H・S・G の ΣνX。温度 (と圧力) が全部の計算で揃っているときだけ出します。"
                               "振動だけのモデルでは H は無く、G の列は F = U − TS の差です",
                               "ΔH, ΔS and ΔG are ΣνX of H, S and G from analysis/thermo.csv (ASE thermochemistry) of each run, given only when every run has the same temperatures (and pressures). "
                               "Vibration-only models have no H, and the G column is then the difference of F = U − TS")]),
                    cres.files.get("reactions", "")))
    rows = [[r["dir"], r["code"] or "-", r["task"] or "-", r["formula"] or "-", "-" if r["natoms"] is None else str(r["natoms"]),
             _g(r["energy_ev"], ".6f"), r["energy_source"] or "-", r["note"] or ""] for r in cres.runs]
    out.append(_cut(Section("compare_runs", L("計算ごと", "Per run"),
                            [L("ディレクトリ", "Directory"), L("コード", "Code"), L("計算の種類", "Task"), L("組成", "Formula"), L("原子数", "Atoms"),
                             L("最終エネルギー [eV]", "Final energy [eV]"), L("エネルギーの出典", "Energy source"), L("注", "Note")],
                            rows, [False, False, False, False, True, True, False, False]), cres.files.get("runs", "")))
    dirs = [r["dir"] for r in cres.runs]
    rows = [[k] + [_cell(cres.conditions.get(k, {}), d) for d in dirs] for k in cres.differing]
    notes = [L(f"条件が違う項目 {len(cres.differing)} 個、一部の計算にだけある項目 {len(cres.partial)} 個 (全項目は {cres.files.get('conditions', '')})",
               f"{len(cres.differing)} settings differ, {len(cres.partial)} are present in only some runs (all settings in {cres.files.get('conditions', '')})")]
    notes += [L("注: ", "note: ") + n for n in cres.notes]
    notes.append(L("比べてよいか (条件の違いが結果に効くか) は、ADIT は判断しません", "ADIT does not judge whether these runs can be compared"))
    out.append(_cut(Section("compare_conditions", L("条件が違う項目", "Settings that differ"), [L("項目", "Setting")] + dirs, rows, [], notes),
                    cres.files.get("conditions", "")))
    return out


def _cell(vals: dict, d: str) -> str:
    if d not in vals:
        return L("(無し)", "(absent)")
    v = vals[d]
    return "-" if v is None else str(v)


def read_export_readme(res) -> tuple[str, str]:
    x = res.tables.get("export") if res is not None else None
    if not x:
        return "", ""
    p = Path(x["dir"]) / "export_README.txt"
    try:
        return x["dir"], p.read_text(encoding="utf-8")
    except OSError:
        return x["dir"], ""


__all__ = ["LABELS", "ROW_LIMIT", "FieldError", "Section", "options_from_fields", "fields_from_options", "reactions_from_rows",
           "rows_from_reactions", "result_sections", "compare_sections", "read_export_readme", "lab", "ph", "choices"]
