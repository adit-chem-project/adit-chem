
from __future__ import annotations

from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import write

from adit.analysis.compute import MoleculeUnwrapper
from adit.lang import L

EXPORT_SUBDIR = "export"
FILES = ("trajectory.extxyz", "trajectory.xyz", "trajectory.pdb", "view.vmd", "vmd_load.tcl", "ovito_pipeline.py", "export_README.txt")
# TRAVIS function keywords (its ">>> List of functions <<<" menu); one answer file is written per keyword
TRAVIS_FUNCTIONS = ("rdf", "cdf", "msd", "hbond", "acf")
_PDB_ATOM = "ATOM  %5d %4s %4s %4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n"


class TrajectoryExporter:

    def __init__(self, out_dir: Path, *, unwrap_molecules: bool = False):
        self.dir = Path(out_dir)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.unwrap = unwrap_molecules
        self._fe = open(self.dir / "trajectory.extxyz", "w", encoding="utf-8", newline="\n")
        self._fx = open(self.dir / "trajectory.xyz", "w", encoding="utf-8", newline="\n")
        self._fp = open(self.dir / "trajectory.pdb", "w", encoding="utf-8", newline="\n")
        self.n = 0
        self.first: Atoms | None = None
        self.last_cell = None
        self.cell_varies = False
        self.mol: MoleculeUnwrapper | None = None

    def add(self, fr: Atoms) -> None:
        if self.first is None:
            self.first = fr.copy()
            self.mol = MoleculeUnwrapper(fr)
        a = self.mol.apply(fr) if self.unwrap and self.mol is not None else fr
        periodic = bool(any(a.pbc)) and a.cell.rank == 3
        if periodic:
            c = np.asarray(a.cell, dtype=float)
            if self.last_cell is not None and not np.allclose(c, self.last_cell, atol=1e-6):
                self.cell_varies = True
            self.last_cell = c
        b = a.copy()
        b.calc = None
        write(self._fe, b, format="extxyz")
        self._fx.write(f"{len(a)}\n{L('フレーム', 'frame')} {self.n} (ADIT)\n")
        for s, p in zip(a.get_chemical_symbols(), a.get_positions()):
            self._fx.write(f"{s:<2s} {p[0]:16.8f} {p[1]:16.8f} {p[2]:16.8f}\n")
        self._write_pdb_model(a, periodic)
        self.n += 1

    def _write_pdb_model(self, a: Atoms, periodic: bool) -> None:
        p = a.get_positions()
        if periodic:
            cp = a.cell.cellpar()
            _, rot = a.cell.standard_form()
            p = p @ rot.T
            self._fp.write("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1\n" % tuple(cp))
        self._fp.write(f"MODEL     {self.n + 1}\n")
        labels = self.mol.labels if self.mol is not None else np.zeros(len(a), int)
        for i, (s, (x, y, z)) in enumerate(zip(a.get_chemical_symbols(), p)):
            self._fp.write(_PDB_ATOM % ((i + 1) % 100000, s, "MOL ", (labels[i] + 1) % 10000, x, y, z, 1.0, 0.0, s.upper()))
        self._fp.write("ENDMDL\n")

    def close(self, *, run_dir: Path, code: str, source: str, dt_frame_fs: float | None, stride: int, skip: int, n_total: int | None,
              rdf_cutoff: float, extra_lines: list[str] | None = None, select: str = "") -> dict:
        self.extra_lines = list(extra_lines or [])
        for f in (self._fe, self._fx, self._fp):
            f.close()
        info = {"dir": str(self.dir), "files": [str(self.dir / f) for f in FILES], "n_frames": self.n, "stride": stride, "skip": skip,
                "dt_frame_fs": dt_frame_fs, "cell_varies": self.cell_varies, "unwrap_requested": self.unwrap,
                "unwrap_molecules": bool(self.unwrap and self.mol is not None and self.mol.active)}
        a = self.first
        if a is None:
            return info
        periodic = bool(any(a.pbc)) and a.cell.rank == 3
        if self.mol is not None:
            info.update(n_molecules=self.mol.n_molecules, largest_molecule_atoms=int(self.mol.sizes.max()) if len(self.mol.sizes) else 0,
                        periodic_networks=int(self.mol.network.sum()))
        if periodic:
            info["cell_A"] = [float(x) for x in a.cell.cellpar()]
        sel = _selection_texts(select)
        info["select"] = select
        (self.dir / "view.vmd").write_text(_vmd_text(a, periodic, self.cell_varies), encoding="utf-8", newline="\n")
        (self.dir / "vmd_load.tcl").write_text(_vmd_tcl_text(a, periodic, self.cell_varies, sel), encoding="utf-8", newline="\n")
        (self.dir / "ovito_pipeline.py").write_text(_ovito_text(rdf_cutoff, sel), encoding="utf-8", newline="\n")
        travis = _travis_answer_files(a, periodic, self.cell_varies)
        for name, text in travis.items():
            (self.dir / name).write_text(text, encoding="utf-8", newline="\n")
        info["files"] += [str(self.dir / name) for name in travis]
        info["travis_answer_files"] = sorted(travis)
        text = _readme_text(a, info, run_dir=run_dir, code=code, source=source, n_total=n_total)
        if self.extra_lines:
            text += "\n".join(self.extra_lines) + "\n"
        (self.dir / "export_README.txt").write_text(text, encoding="utf-8", newline="\n")
        return info


def _selection_texts(select: str) -> dict[str, str]:
    """The ADIT --select expression rewritten for VMD and OVITO; a note replaces what cannot be rewritten."""
    out = {"adit": select.strip(), "vmd": "", "ovito": "", "ovito_note": ""}
    if not out["adit"]:
        return out
    from adit.analysis.select import SelectionError, to_ovito, to_vmd

    try:
        out["vmd"] = to_vmd(out["adit"])
    except SelectionError as ex:
        out["vmd_note"] = str(ex)
    try:
        out["ovito"] = to_ovito(out["adit"])
    except SelectionError as ex:
        out["ovito_note"] = str(ex)
    return out


def _vmd_text(a: Atoms, periodic: bool, cell_varies: bool) -> str:
    lines = [L("# ADIT が書いた VMD の読み込みスクリプト。使い方: このディレクトリで  vmd -e view.vmd",
               "# VMD script written by ADIT. Usage: in this directory,  vmd -e view.vmd"),
             "mol new trajectory.xyz type xyz waitfor all"]
    if periodic:
        cp = a.cell.cellpar()
        if cell_varies:
            lines.append(L("# セルはフレームごとに変わります。ここでは最初のフレームのセルを全フレームに付けます (各フレームのセルは trajectory.extxyz / trajectory.pdb)",
                           "# The cell changes from frame to frame. The first frame's cell is applied to all frames here (per-frame cells are in trajectory.extxyz / trajectory.pdb)"))
        lines += [L("# セル: a b c [Å] と alpha beta gamma [度] (PBCTools。VMD 1.8.6 から同梱)", "# cell: a b c [Å] and alpha beta gamma [deg] (PBCTools, bundled since VMD 1.8.6)"),
                  "pbc set {%.6f %.6f %.6f %.4f %.4f %.4f} -all" % tuple(cp), "pbc box"]
    lines += ["mol modstyle 0 [molinfo top] CPK", ""]
    return "\n".join(lines)


def _vmd_tcl_text(a: Atoms, periodic: bool, cell_varies: bool, sel: dict[str, str]) -> str:
    # Command names and arguments follow the VMD User's Guide (mol, color, display, axes, render) and the PBCTools plugin page
    elements = list(dict.fromkeys(a.get_chemical_symbols()))
    lines = [L("# ADIT が書いた VMD の Tcl スクリプト (代表的な表示まで)。使い方: このディレクトリで  vmd -e vmd_load.tcl",
               "# VMD Tcl script written by ADIT (loads the trajectory and sets up a typical display). Usage: in this directory,  vmd -e vmd_load.tcl"),
             L("# ADIT では実行を確かめていません。コマンドは VMD の User's Guide (mol / color / display / axes / render) と PBCTools の説明で確かめたものです",
               "# Not run by ADIT. Commands were checked against the VMD User's Guide (mol / color / display / axes / render) and the PBCTools page"),
             "", "mol new trajectory.xyz type xyz waitfor all", "set m [molinfo top]"]
    if periodic:
        cp = a.cell.cellpar()
        if cell_varies:
            lines.append(L("# セルはフレームごとに変わります。最初のフレームのセルを全フレームに付けます (各フレームのセルは trajectory.pdb を読めば付きます)",
                           "# The cell changes between frames; the first frame's cell is applied to all frames (load trajectory.pdb to get per-frame cells)"))
        lines += ["pbc set {%.6f %.6f %.6f %.4f %.4f %.4f} -all" % tuple(cp), "pbc box -style lines",
                  L("# 分子をセルの中へ戻すなら (分子単位で。fragment は結合でつながった単位):",
                    "# to wrap molecules back into the cell (per fragment, i.e. per bonded unit):"),
                  "# pbc wrap -compound fragment -all"]
    lines += ["", L("# 表示 0: 全原子を CPK、元素名で色分け", "# representation 0: all atoms as CPK, colored by name"),
              "mol modstyle 0 $m CPK", "mol modcolor 0 $m Name", "mol modselect 0 $m all", ""]
    rep = 1
    if sel["vmd"]:
        lines += [L(f"# 表示 {rep}: ADIT の --select \"{sel['adit']}\" と同じ原子 (VMD の選択式に直したもの。番号は 0 から) を VDW で",
                    f"# representation {rep}: the atoms of ADIT's --select \"{sel['adit']}\" (rewritten in VMD's selection language; indices from 0) as VDW"),
                  "mol addrep $m", f"mol modselect {rep} $m {{{sel['vmd']}}}", f"mol modstyle {rep} $m VDW", f"mol modcolor {rep} $m Name", ""]
        rep += 1
    first = elements[0]
    lines += [L(f"# 表示 {rep} (例。使うなら先頭の # を外す): 元素 {first} だけを VDW で。選択式の書き方は VMD の User's Guide の Selection の章",
                f"# representation {rep} (example; remove the leading # to use): element {first} only, as VDW. Selection syntax: VMD User's Guide, Selections chapter"),
              "# mol addrep $m", f"# mol modselect {rep} $m {{name {first}}}", f"# mol modstyle {rep} $m VDW", "",
              L(f"# 表示 {rep + 1} (例): フレームごとに結合を引き直す DynamicBonds。引数は 距離のカットオフ [Å]、結合の太さ、分割数",
                f"# representation {rep + 1} (example): DynamicBonds recomputed every frame; arguments are the distance cutoff [Å], bond radius and resolution"),
              "# mol addrep $m", f"# mol modselect {rep + 1} $m {{all}}", f"# mol modstyle {rep + 1} $m DynamicBonds <{L('距離のカットオフ Å', 'distance cutoff in Å')}> 0.3 12", "",
              L("# 画面: 白い背景、平行投影、座標軸を消す", "# display: white background, orthographic projection, axes off"),
              "color Display Background white", "display projection Orthographic", "axes location Off", "",
              L("# 画像に描くなら (例。最後のフレームへ移って Tachyon で描く):", "# to render an image (example: go to the last frame and render with Tachyon):"),
              "# animate goto end", "# render TachyonInternal frame_last.tga", ""]
    return "\n".join(lines)


_OVITO_TEMPLATE = '''#!/usr/bin/env python
"""{doc}"""
import sys
from pathlib import Path

from ovito.io import export_file, import_file
import ovito.modifiers as om

here = Path(__file__).resolve().parent
cutoff = float(sys.argv[1]) if len(sys.argv) > 1 else {cutoff:.4f}  # Å
pipeline = import_file(str(here / "trajectory.extxyz"))  # {cell_comment}
if hasattr(om, "RadialDistributionFunctionModifier"):  # {new_comment}
    pipeline.modifiers.append(om.RadialDistributionFunctionModifier(cutoff=cutoff, number_of_bins=200, partial=True))
else:  # {old_comment}
    pipeline.modifiers.append(om.CoordinationAnalysisModifier(cutoff=cutoff, number_of_bins=200, partial=True))
data = pipeline.compute()
key = next((k for k in ("rdf", "coordination-rdf") if k in data.tables), None)
if key is not None:
    export_file(pipeline, str(here / "ovito_rdf.*.txt"), "txt/table", key=key, multiple_frames=True)
    print("{wrote_rdf}", here / "ovito_rdf.*.txt")
if "Coordination" in data.particles.keys():
    export_file(pipeline, str(here / "ovito_coordination.xyz"), "xyz", multiple_frames=True,
                columns=["Particle Type", "Position.X", "Position.Y", "Position.Z", "Coordination"])
    print("{wrote_cn}", here / "ovito_coordination.xyz")
else:
    print("{no_cn}")

# ----------------------------------------------------------------------------------------------------
# {opt_title}
# {opt_note}
# {opt_ref}
# ----------------------------------------------------------------------------------------------------

# {s_title}
# {s_note}
# pipeline.modifiers.append(om.CommonNeighborAnalysisModifier(mode=om.CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff))
# pipeline.modifiers.append(om.PolyhedralTemplateMatchingModifier(rmsd_cutoff=<{ptm_rmsd}>, output_orientation=False))
# pipeline.modifiers.append(om.AcklandJonesModifier())
# export_file(pipeline, str(here / "ovito_structure_counts.txt"), "txt/attr", multiple_frames=True,
#             columns=["SourceFrame", "CommonNeighborAnalysis.counts.FCC", "CommonNeighborAnalysis.counts.HCP",
#                      "CommonNeighborAnalysis.counts.BCC", "CommonNeighborAnalysis.counts.ICO", "CommonNeighborAnalysis.counts.OTHER"])

# {ws_title}
# {ws_note}
# from ovito.pipeline import FileSource
# ws = om.WignerSeitzAnalysisModifier(per_type_occupancies=False, output_displaced=False)
# ws.reference = FileSource()
# ws.reference.load(str(here / "<{ws_ref}>"))
# pipeline.modifiers.append(ws)
# export_file(pipeline, str(here / "ovito_defects.txt"), "txt/attr", multiple_frames=True,
#             columns=["SourceFrame", "WignerSeitz.vacancy_count", "WignerSeitz.interstitial_count"])

# {d_title}
# {d_note}
# pipeline.modifiers.append(om.CalculateDisplacementsModifier(reference_frame=0, minimum_image_convention=True))
# pipeline.modifiers.append(om.AtomicStrainModifier(cutoff=<{strain_cutoff}>, reference_frame=0, output_strain_tensors=False))
# export_file(pipeline, str(here / "ovito_displacement_strain.xyz"), "xyz", multiple_frames=True,
#             columns=["Particle Type", "Position.X", "Position.Y", "Position.Z", "Displacement Magnitude", "Shear Strain", "Volumetric Strain"])

# {c_title}
# {c_note}
# pipeline.modifiers.append(om.ClusterAnalysisModifier(cutoff=<{cluster_cutoff}>, sort_by_size=True, compute_com=True, unwrap_particles=True))
# export_file(pipeline, str(here / "ovito_clusters.*.txt"), "txt/table", key="clusters", multiple_frames=True)

# {b_title}
# {b_note}
# pipeline.modifiers.append(om.SpatialBinningModifier(property="<{bin_property}>", direction=om.SpatialBinningModifier.Direction.Z,
#                                                     bin_count=<{bin_count}>, reduction_operation=om.SpatialBinningModifier.Operation.Mean))
# pipeline.modifiers.append(om.TimeAveragingModifier(operate_on="table:binning"))
# data = pipeline.compute()
# export_file(data.tables["binning[average]"], str(here / "ovito_profile_average.txt"), "txt/table")

# {e_title}
# {e_note}
# pipeline.modifiers.append(om.ExpressionSelectionModifier(expression='{expression}'))
# export_file(pipeline, str(here / "ovito_selection_count.txt"), "txt/attr", multiple_frames=True, columns=["SourceFrame", "ExpressionSelection.count"])
'''


def _ovito_text(cutoff: float, sel: dict[str, str] | None = None) -> str:
    sel = sel or {"adit": "", "ovito": "", "ovito_note": ""}
    if sel["ovito"]:
        expression = sel["ovito"]
        e_note = L(f"ADIT の --select \"{sel['adit']}\" を OVITO の式に直したもの (元素は ParticleType の名前、番号は ParticleIndex で 0 から)",
                   f"ADIT's --select \"{sel['adit']}\" rewritten as an OVITO expression (elements by ParticleType name, indices as ParticleIndex from 0)")
    else:
        expression = 'ParticleType == "<' + L("元素", "element") + '>" && Position.Z < <' + L("z の上限 Å", "upper z in Å") + ">"
        e_note = (sel["ovito_note"] if sel.get("ovito_note") else
                  L("式の書き方 (ParticleType、Position.X/Y/Z、ParticleIndex、&& || == != < <= > >=) は OVITO の Expression selection の説明を見てください",
                    "see OVITO's Expression selection page for the syntax (ParticleType, Position.X/Y/Z, ParticleIndex, && || == != < <= > >=)"))
    return _OVITO_TEMPLATE.format(
        opt_title=L("以下は任意の解析の雛形です。使う行の先頭の \"# \" を外し、<...> を埋めてください。数値は入れていません (系によるため)",
                    "Optional analyses (templates). Remove the leading \"# \" of the lines you want and fill in the <...> holes. No numerical values are given (they depend on the system)"),
        opt_note=L("どれも ADIT では実行していません。クラス名と引数名は OVITO 3.16 の Python リファレンスで確かめました",
                   "None of these has been run by ADIT. Class and argument names were checked against the OVITO 3.16 Python reference"),
        opt_ref="https://www.ovito.org/docs/current/python/modules/ovito_modifiers.html",
        s_title=L("構造同定 (結晶の局所構造: FCC / HCP / BCC / ICO ...)。3 つのうち 1 つを使う。結果は原子ごとの Structure Type と個数の属性",
                  "Structure identification (local crystal structure: FCC / HCP / BCC / ICO ...). Use one of the three; output is the per-atom Structure Type and the count attributes"),
        s_note=L("CNA (adaptive cutoff) / PTM (RMSD の閾値。リファレンスの既定値は 0.1) / Ackland-Jones (引数なし)",
                 "CNA (adaptive cutoff) / PTM (RMSD cutoff; the reference's default is 0.1) / Ackland-Jones (no arguments)"),
        ptm_rmsd=L("RMSD の閾値", "RMSD cutoff"),
        ws_title=L("Wigner-Seitz の欠陥解析 (空孔と格子間原子の数)。参照構造 (完全結晶) は人が用意して指定する",
                   "Wigner-Seitz defect analysis (vacancy and interstitial counts). The reference (perfect crystal) is supplied by you"),
        ws_note=L("参照構造のファイルは OVITO が読める形式 (extxyz、POSCAR、LAMMPS data など)。原子数と並びが参照のサイトに対応している必要はない",
                  "The reference file can be any format OVITO reads (extxyz, POSCAR, LAMMPS data ...); atoms need not be ordered like the reference sites"),
        ws_ref=L("参照構造のファイル名", "reference structure file"),
        d_title=L("変位ベクトルと原子ひずみ (最初のフレームを参照にする。reference_frame を変えれば別のフレーム)",
                  "Displacement vectors and atomic strain (frame 0 is the reference; change reference_frame for another frame)"),
        d_note=L("AtomicStrain の cutoff は近傍を数える半径 [Å]。出力は Displacement Magnitude、Shear Strain、Volumetric Strain",
                 "cutoff of AtomicStrain is the neighbor radius [Å]; outputs are Displacement Magnitude, Shear Strain and Volumetric Strain"),
        strain_cutoff=L("近傍の半径 Å", "neighbor radius in Å"),
        c_title=L("クラスタ解析 (カットオフ以内でつながった原子の集まり)。サイズ順に並べ、重心も出す",
                  "Cluster analysis (atoms connected within the cutoff), sorted by size, with centers of mass"),
        c_note=L("結果はフレームごとの表 clusters (Cluster Identifier、Cluster Size、Center of Mass) と属性 ClusterAnalysis.cluster_count",
                 "Output: the per-frame table clusters (Cluster Identifier, Cluster Size, Center of Mass) and the attribute ClusterAnalysis.cluster_count"),
        cluster_cutoff=L("カットオフ Å", "cutoff in Å"),
        b_title=L("空間ビニング + 時間平均 (z 方向に区切って原子の量を平均した分布。例: 電荷や速度の成分)",
                  "Spatial binning + time averaging (a profile along z of a per-atom quantity, e.g. a charge or velocity component)"),
        b_note=L("property は原子ごとの量の名前 (Position.X、Velocity.Z、Charge など、軌跡に入っている列)。時間平均した表は binning[average]",
                 "property is the name of a per-atom quantity present in the trajectory (Position.X, Velocity.Z, Charge ...); the averaged table is binning[average]"),
        bin_property=L("原子ごとの量", "per-atom property"), bin_count=L("区間の数", "number of bins"),
        e_title=L("式で原子を選ぶ (Expression selection)。選んだ数は属性 ExpressionSelection.count。続けて om.DeleteSelectedModifier() などを足せる",
                  "Select atoms by an expression (Expression selection). The count is the attribute ExpressionSelection.count; e.g. om.DeleteSelectedModifier() can follow"),
        e_note=e_note, expression=expression,
        doc=L("ADIT が書いた OVITO の Python スクリプト。ovito の Python パッケージ (MIT、pip install ovito) で動く。ADIT では実行を確かめていない。\n"
              "trajectory.extxyz を読み、元素の組ごとの動径分布関数と、原子ごとの配位数 (カットオフ以内の原子の数) を求めて書き出す。\n"
              "使い方: python ovito_pipeline.py [カットオフ Å]\n"
              "末尾に、構造同定・欠陥解析・変位とひずみ・クラスタ・空間ビニング・式による選択の雛形 (コメントアウト) がある。",
              "OVITO Python script written by ADIT. Runs with the ovito Python package (MIT, pip install ovito). Not run by ADIT.\n"
              "Reads trajectory.extxyz and writes partial radial distribution functions and per-atom coordination numbers (atoms within the cutoff).\n"
              "Usage: python ovito_pipeline.py [cutoff in Å]\n"
              "Commented-out templates for structure identification, defect analysis, displacements and strain, clusters, spatial binning and expression selection follow at the end."),
        cutoff=cutoff,
        cell_comment=L("拡張 xyz の Lattice= と pbc= をセルとして読む", "the Lattice= and pbc= keys of extended XYZ become the cell"),
        new_comment=L("新しいバージョンの名前", "name in recent versions"),
        old_comment=L("古いバージョンの名前 (原子ごとの配位数 Coordination も出す)", "name in older versions (also outputs per-atom Coordination)"),
        wrote_rdf=L("動径分布関数:", "RDF:"), wrote_cn=L("配位数:", "coordination:"),
        no_cn=L("このバージョンの OVITO は原子ごとの配位数を出しません。ADIT の analysis/rdf.json の n と n_reverse (配位数の曲線) を使ってください",
                "this OVITO version does not output per-atom coordination; use n and n_reverse (coordination curves) in ADIT's analysis/rdf.json"))


_TRAVIS_FUNCTION_WORDS = {
    "rdf": ("動径分布関数", "radial distribution function"),
    "cdf": ("組み合わせ分布関数 (距離と角度など 2 つの量の 2 次元分布)", "combined distribution function (2D distribution of two quantities, e.g. distance and angle)"),
    "msd": ("平均二乗変位", "mean square displacement"),
    "hbond": ("水素結合の動力学", "hydrogen bond dynamics"),
    "acf": ("自己相関関数 (速度自己相関はこの中)", "autocorrelation functions (the velocity autocorrelation is here)"),
}


def _travis_answer_files(a: Atoms, periodic: bool, cell_varies: bool) -> dict[str, str]:
    """Answer files for `travis -i`: the questions up to the function menu, for an orthorhombic fixed cell."""
    if not periodic:
        return {}
    cp = a.cell.cellpar()
    if not np.allclose(cp[3:], 90.0, atol=1e-3):
        return {}
    same = bool(np.allclose(cp[:3], cp[0], atol=1e-6))
    out: dict[str, str] = {}
    for fn in TRAVIS_FUNCTIONS:
        what = L(*_TRAVIS_FUNCTION_WORDS[fn])
        lines = [
            L(f"! ADIT が書いた TRAVIS の答えファイル ({fn}: {what})。使い方: このディレクトリで  travis -p trajectory.xyz -i travis_{fn}.in",
              f"! TRAVIS answer file written by ADIT ({fn}: {what}). Usage: in this directory,  travis -p trajectory.xyz -i travis_{fn}.in"),
            L("! 書式: 1 行が 1 つの答え、空行は既定値、! で始まる行はコメント (TRAVIS の Quick Start Guide とソースで確認)",
              "! Format: one answer per line, an empty line takes the default, lines starting with ! are comments (confirmed in the TRAVIS Quick Start Guide and source)"),
            L("! この並びは 2026-07-23 版のソースと 2018 年の公式チュートリアルの例に合わせたもので、ADIT では実行して確かめていません (未確認)。",
              "! The order follows the source of the 2026-07-23 version and the example in the 2018 official tutorial; ADIT has not run it (unverified)."),
            L("! 質問が版で変わることがあります。ファイルが尽きると TRAVIS はキーボード入力に切り替わるので、残りは対話で答え、",
              "! Questions change between versions. When the file runs out TRAVIS switches to keyboard input; answer the rest interactively,"),
            L("! 一度 -i なしで走らせて TRAVIS が書く input.txt を次回の -i に使うのが公式の勧めです。",
              "! and reuse the input.txt that TRAVIS writes in a run without -i, as the official tutorial recommends."),
            "! Use the advanced mode until the analysis selection menu (y/n)? [no]", "",
        ]
        if cell_varies:
            lines.append(L("! セルはフレームごとに変わりますが、ここでは最初のフレームのセルを固定して答えます (NPT なら advanced mode で時間依存のセルを選んでください)",
                           "! The cell changes between frames; the first frame's cell is given as fixed here (for NPT choose a time-dependent cell in the advanced mode)"))
        if same:
            lines += ["! Are the 3 cell vectors of the same size (yes/no)? [yes]", "",
                      "! Enter length of cell vector in pm:", f"{cp[0] * 100:.4f}"]
        else:
            lines += ["! Are the 3 cell vectors of the same size (yes/no)? [yes]", "no"]
            for axis, length in zip("XYZ", cp[:3]):
                lines += [f"! Enter length of {axis} cell vector in pm:", f"{length * 100:.4f}"]
        lines += [L("! (金属や希ガスの元素があると、ここで元素ごとに「結合の認識から外すか」を聞かれます。その分の行は入れていません)",
                    "! (with metal or noble-gas elements TRAVIS asks per element whether to exclude it from bond recognition here; those lines are not included)"),
                  "! Create images of the structural formulas (y/n)? [no]", "",
                  "! Accept these molecules (y) or change something (n)? [yes]", "",
                  "! Which functions to compute (comma separated)?", fn,
                  "! Use the advanced mode for the main part (y/n)? [no]", "",
                  L(f"! ここから先 ({fn} の設定、時間の刻み、分子と原子の選択、フレームの範囲) は対話で答えてください。目安は export_README.txt の TRAVIS の節",
                    f"! From here on (settings of {fn}, time step, molecules and atoms, frame range) answer interactively; see the TRAVIS section of export_README.txt"), ""]
        out[f"travis_{fn}.in"] = "\n".join(lines)
    return out


def _travis_readme_lines(info: dict, dt: float | None, stride: int, periodic: bool) -> list[str]:
    files = info.get("travis_answer_files") or []
    out = [L("  TRAVIS の答えファイル (-i):", "  TRAVIS answer files (-i):")]
    if files:
        out += [L(f"           {', '.join(files)}  (関数ごとに 1 つ。travis -p trajectory.xyz -i travis_rdf.in のように使う)",
                  f"           {', '.join(files)}  (one per function; use as travis -p trajectory.xyz -i travis_rdf.in)"),
                L("           書式は 1 行 1 答、空行 = 既定値、! で始まる行はコメント。中身はセルの大きさ [pm] と関数の選択までで、",
                  "           format: one answer per line, empty line = default, lines starting with ! are comments. They cover the cell size [pm] and the function choice;"),
                L("           ファイルが尽きると TRAVIS はキーボード入力に切り替わるので、残りは対話で答えます (下の目安)",
                  "           when the file runs out TRAVIS switches to keyboard input, so answer the rest interactively (see below)")]
    elif periodic:
        out.append(L("           書いていません (直方体でないセルは advanced mode でベクトルと角度を答える必要があり、その手順は未確認です)",
                     "           not written (a non-orthorhombic cell needs the advanced mode with vectors and angles, and that sequence is unverified)"))
    else:
        out.append(L("           書いていません (周期境界の無い軌跡。TRAVIS はそれでもセルの大きさを聞くので、分子が入る大きさを pm で答えてください)",
                     "           not written (non-periodic trajectory; TRAVIS still asks for a cell size, answer a size in pm that holds the molecule)"))
    out += [L("           質問の並び (2026-07-23 版のソースを読んだもので、実行では未確認。版によって変わります):",
              "           question sequence (read from the source of the 2026-07-23 version, not verified by running; changes between versions):"),
            L("             1. advanced mode を使うか [no]  2. 3 つのセルベクトルが同じ長さか [yes] → 長さを pm で",
              "             1. use the advanced mode [no]  2. are the 3 cell vectors the same size [yes] -> length in pm"),
            L("             3. 構造式の画像を作るか [no]  4. 認識した分子を受け入れるか [yes]  5. 計算する関数 (rdf, cdf, msd, hbond, acf, sdf, vdf, ...)",
              "             3. draw structural formulas [no]  4. accept the recognized molecules [yes]  5. functions to compute (rdf, cdf, msd, hbond, acf, sdf, vdf, ...)"),
            L("             6. 本体で advanced mode を使うか [no]  7. (msd・acf など時間の解析なら) フレーム間の時間 [fs]",
              "             6. use the advanced mode for the main part [no]  7. (for time-dependent analyses such as msd, acf) time between frames [fs]")]
    if dt:
        out.append(L(f"                → このデータでは 1 フレームあたり {dt * stride:g} fs", f"                -> {dt * stride:g} fs per frame for this data"))
    out += [L("             8. 参照分子と観測分子の番号 (分子の種類が 2 つ以上のとき)、原子の指定 (\"*\" = 全部、#2 = 質量中心)、",
              "             8. reference and observed molecule numbers (with 2 or more molecule types), atoms (\"*\" = all, #2 = center of mass),"),
            L("                rdf の最大半径 [pm] と区間の数、msd の相関の深さ [フレーム] など  9. 開始フレーム [1]・読むフレーム数 [all]・間引き [1]",
              "                max radius [pm] and bin count of the rdf, correlation depth [frames] of the msd, etc.  9. first frame [1], frames to read [all], stride [1]"),
            L("           速度自己相関は関数名 acf (vacf ではありません)。-i を付けた実行では input.txt は書かれず、答えが足りないと終了コード 0 のまま止まります",
              "           the velocity autocorrelation is the function acf (not vacf). A run with -i does not write input.txt, and stops with exit code 0 when answers run out"),
            L("           TRAVIS の引用: Brehm et al., J. Chem. Phys. 152, 164105 (2020) と Brehm, Kirchner, J. Chem. Inf. Model. 51, 2007 (2011) (公式サイトの Cite の指示)",
              "           citing TRAVIS: Brehm et al., J. Chem. Phys. 152, 164105 (2020) and Brehm, Kirchner, J. Chem. Inf. Model. 51, 2007 (2011) (per the Cite page of the official site)")]
    return out


def _readme_text(a: Atoms, info: dict, *, run_dir: Path, code: str, source: str, n_total: int | None) -> str:
    syms = a.get_chemical_symbols()
    counts = {s: syms.count(s) for s in dict.fromkeys(syms)}
    elems = ", ".join(f"{s} {n}" for s, n in counts.items())
    periodic = bool(any(a.pbc)) and a.cell.rank == 3
    dt = info.get("dt_frame_fs")
    stride = info["stride"]
    out = [L("ADIT が書き出した軌跡 (TRAVIS・OVITO・VMD 用)", "Trajectory exported by ADIT (for TRAVIS, OVITO and VMD)"), "",
           L(f"元の計算: {run_dir} ({code}、{source})", f"source calculation: {run_dir} ({code}, {source})"),
           L(f"フレーム数: {info['n_frames']} (元の軌跡 {n_total if n_total is not None else '?'} フレームのうち、先頭 {info['skip']} フレームを除き、{stride} フレームに 1 回)",
             f"frames: {info['n_frames']} (of {n_total if n_total is not None else '?'} in the original trajectory; first {info['skip']} skipped, every {stride})")]
    if dt:
        out.append(L(f"1 フレームあたりの時間: {dt * stride:g} fs (元の書き出しの間隔 {dt:g} fs × 間引き {stride})",
                     f"time per frame: {dt * stride:g} fs (original output interval {dt:g} fs x stride {stride})"))
    else:
        out.append(L("1 フレームあたりの時間: 不明 (出力から読めません。MD でない計算かもしれません)",
                     "time per frame: unknown (not readable from the output; perhaps not an MD run)"))
    out.append(L(f"元素: {elems} (合計 {len(a)} 原子)", f"elements: {elems} ({len(a)} atoms in total)"))
    if periodic:
        c = a.cell.cellpar()
        out += [L("セル (最初のフレーム):", "cell (first frame):"),
                f"  a = {c[0]:.6f} Å = {c[0] * 100:.4f} pm", f"  b = {c[1]:.6f} Å = {c[1] * 100:.4f} pm", f"  c = {c[2]:.6f} Å = {c[2] * 100:.4f} pm",
                f"  alpha = {c[3]:.4f}°, beta = {c[4]:.4f}°, gamma = {c[5]:.4f}°",
                L("  ベクトル [Å]:", "  vectors [Å]:")] + [f"    {v[0]:14.8f} {v[1]:14.8f} {v[2]:14.8f}" for v in np.asarray(a.cell)]
        if info.get("cell_varies"):
            out.append(L("  セルはフレームごとに変わります。trajectory.extxyz と trajectory.pdb には各フレームのセルがあり、trajectory.xyz には入りません",
                         "  the cell changes between frames; trajectory.extxyz and trajectory.pdb hold each frame's cell, trajectory.xyz has none"))
    else:
        out.append(L("セル: なし (周期境界の無い分子の計算)", "cell: none (non-periodic molecular calculation)"))
    if "n_molecules" in info:
        largest = info["largest_molecule_atoms"]
        s = L(f"距離による原子グループの推定 (最初のフレームで共有結合半径の和 × 1.2 以内): {info['n_molecules']} 個、最大 {largest} 原子。伸びた結合は別グループになる場合があります",
              f"distance-based atom groups (within 1.2 x the sum of covalent radii in the first frame): {info['n_molecules']}; "
              f"largest group: {largest} {'atom' if largest == 1 else 'atoms'}. Stretched bonds may be split into separate groups")
        out.append(s)
        if info["unwrap_molecules"]:
            out.append(L("  書き出す前に、分子を周期境界でつなぎ直しました (分子の中心がセルの中に来る位置へ)。区切りは最初のフレームのものを最後まで使います",
                         "  molecules were made whole across periodic boundaries before writing (each placed with its center inside the cell); the partition of the first frame is kept throughout"))
            if info.get("periodic_networks"):
                out.append(L(f"  周期的につながった集まり {info['periodic_networks']} 個 (結晶など) はつなぎ直さず、元の位置のままです",
                             f"  {info['periodic_networks']} periodically connected groups (e.g. crystals) were left as they were"))
        elif not periodic and info["unwrap_requested"]:
            out.append(L("  非周期軌跡のため、周期境界でのつなぎ直しは行っていません",
                         "  no periodic unwrapping was performed because this trajectory is non-periodic"))
        elif periodic:
            out.append(L("  つなぎ直しはしていません (座標は計算コードが書いたまま)。つなぎ直すには adit-analyze --export --unwrap-molecules",
                         "  molecules were not made whole (coordinates as written by the code); use adit-analyze --export --unwrap-molecules"))
        out.append(L("  trajectory.pdb の残基番号は、この推定グループの番号 (1 から)",
                     "  residue numbers in trajectory.pdb are these estimated group numbers (from 1)"))
    out += ["", L("ファイル:", "files:"),
            (L("  trajectory.extxyz  拡張 xyz (コメント行に Lattice= と pbc=)。OVITO・ASE で開く", "  trajectory.extxyz  extended XYZ (Lattice= and pbc= in the comment line); open with OVITO or ASE")
             if periodic else L("  trajectory.extxyz  拡張 xyz (セル・周期境界なし)。OVITO・ASE で開く",
                                "  trajectory.extxyz  extended XYZ (no cell or periodic boundaries); open with OVITO or ASE")),
            L("  trajectory.xyz     素の xyz (セルなし)。TRAVIS・VMD で開く", "  trajectory.xyz     plain XYZ (no cell); open with TRAVIS or VMD"),
            (L("  trajectory.pdb     複数モデルの PDB (各モデルに CRYST1)", "  trajectory.pdb     multi-model PDB (CRYST1 in each model)")
             if periodic else L("  trajectory.pdb     複数モデルの PDB (セル情報なし)", "  trajectory.pdb     multi-model PDB (no cell information)")),
            (L("  view.vmd           VMD の読み込みとセル", "  view.vmd           VMD loading and cell")
             if periodic else L("  view.vmd           VMD の読み込み", "  view.vmd           VMD loading")),
            L("  vmd_load.tcl       VMD の読み込みと代表的な表示 (CPK、選んだ原子の VDW、白背景、平行投影。例の行はコメントアウト)",
              "  vmd_load.tcl       VMD loading plus a typical display (CPK, selected atoms as VDW, white background, orthographic; example lines are commented out)"),
            L("  ovito_pipeline.py  OVITO の Python スクリプト (ADIT では実行を確かめていない)。動径分布関数と配位数のあとに、構造同定 (CNA / PTM / Ackland-Jones)、",
              "  ovito_pipeline.py  OVITO Python script (not run by ADIT). After the RDF and coordination, commented-out templates for structure identification (CNA / PTM / Ackland-Jones),"),
            L("                     Wigner-Seitz 欠陥解析、変位と原子ひずみ、クラスタ解析、空間ビニング + 時間平均、式による選択の雛形 (コメントアウト) がある",
              "                     Wigner-Seitz defect analysis, displacements and atomic strain, cluster analysis, spatial binning + time averaging and expression selection follow")]
    if info.get("select"):
        out.append(L(f"                     ADIT の --select \"{info['select']}\" は、vmd_load.tcl では VMD の選択式に、ovito_pipeline.py では OVITO の式に直して入れてあります",
                     f"                     ADIT's --select \"{info['select']}\" is rewritten as a VMD selection in vmd_load.tcl and as an OVITO expression in ovito_pipeline.py"))
    for name in info.get("travis_answer_files") or []:
        fn = name[len("travis_"):-len(".in")]
        out.append(L(f"  {name:<18} TRAVIS の答えファイル ({fn}: {L(*_TRAVIS_FUNCTION_WORDS[fn])})。travis -p trajectory.xyz -i {name}",
                     f"  {name:<18} TRAVIS answer file ({fn}: {L(*_TRAVIS_FUNCTION_WORDS[fn])}); travis -p trajectory.xyz -i {name}"))
    out += ["", L("開き方:", "how to open:"),
            L("  TRAVIS:  travis -p trajectory.xyz", "  TRAVIS:  travis -p trajectory.xyz")]
    if periodic:
        c = a.cell.cellpar()
        out.append(L(f"           対話でセルの大きさを聞かれたら pm で答える: {c[0] * 100:.2f} {c[1] * 100:.2f} {c[2] * 100:.2f}"
                     + (" (直方体でないセルの答え方は TRAVIS のバージョンの説明を見てください)" if not np.allclose(c[3:], 90.0) else ""),
                     f"           when asked for the cell size, answer in pm: {c[0] * 100:.2f} {c[1] * 100:.2f} {c[2] * 100:.2f}"
                     + (" (for non-orthogonal cells, see the documentation of your TRAVIS version)" if not np.allclose(c[3:], 90.0) else "")))
    if dt:
        out.append(L(f"           時間の刻みを聞かれたら 1 フレームあたり {dt * stride:g} fs", f"           when asked for the time step, {dt * stride:g} fs per frame"))
    out += [L("           対話の答えを 1 行ずつ書いたファイルを -i で渡すと、同じ解析を繰り返せます: travis -p trajectory.xyz -i <答えのファイル> (TRAVIS の man page の -i)",
              "           a file with the interactive answers, one per line, can be passed with -i to repeat the analysis: travis -p trajectory.xyz -i <answers file> (-i in the TRAVIS man page)")]
    out += _travis_readme_lines(info, dt, stride, periodic)
    out += [L("  VMD:     vmd -e view.vmd   (このディレクトリで。表示まで整えるなら vmd -e vmd_load.tcl)",
              "  VMD:     vmd -e view.vmd   (in this directory; vmd -e vmd_load.tcl also sets up the display)"),
            (L("  OVITO:   GUI で trajectory.extxyz を開く (セルと周期境界も読まれる)。または pip install ovito のあと python ovito_pipeline.py",
               "  OVITO:   open trajectory.extxyz in the GUI (cell and periodicity are read), or pip install ovito and run python ovito_pipeline.py")
             if periodic else L("  OVITO:   GUI で trajectory.extxyz を開く (セルと周期境界はありません)。または pip install ovito のあと python ovito_pipeline.py",
                                "  OVITO:   open trajectory.extxyz in the GUI (no cell or periodic boundaries), or pip install ovito and run python ovito_pipeline.py")),
            ""]
    return "\n".join(out)
