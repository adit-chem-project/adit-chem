#!/usr/bin/env python
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
