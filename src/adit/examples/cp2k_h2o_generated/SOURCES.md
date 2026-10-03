# examples/cp2k_h2o_generated —— CP2K の生成器の出力 (一点計算) と実際に実行した結果

ADIT の生成器 (`codes/cp2k.py`) が作った水分子の一点計算と、CP2K 2026.2 (conda-forge `cp2k`、`cp2k.psmp`) で実行した結果。

| 項目 | 内容 |
|---|---|
| 構造 | ASE `molecule("H2O")`。非周期。箱は一辺 6 Å の立方体 (画面の「分子の箱」)、`&TOPOLOGY &CENTER_COORDINATES` で中心へ |
| 条件 | PBE、DZVP-MOLOPT-SR-GTH (O, H)、GTH-PBE-q6 / GTH-PBE-q1、CUTOFF 280 Ry、REL_CUTOFF 40 Ry (どちらも CP2K の既定と同じ値を例として入れた)、POISSON_SOLVER MT、PERIODIC NONE |
| 基底・擬ポテンシャル | CP2K の data ディレクトリの `BASIS_MOLOPT` と `GTH_POTENTIALS` から、使う 4 項目だけを `BASIS_adit` / `POTENTIAL_adit` に写したもの (CP2K の配布物、GPL-2.0-or-later) |
| キーワードの出典 | CP2K 2026.2 の `cp2k.psmp --xml` が書く入力の説明 (cp2k_input.xml。https://manual.cp2k.org/ と同じ内容) |
| 実行 | `mpirun -np 1 cp2k.psmp -i cp2k.inp` (OpenMP 1 スレッド)。2026-09-12 |
| 結果 | 正常終了。所要時間 4.0 s。SCF は 10 回で収束。全エネルギー -17.219399129357214 Hartree (`ENERGY| Total FORCE_EVAL`)。力は output.log の `FORCES|` の行 |

同じ入力を `tests/test_cp2k.py::test_cp2k_real_run` が実行し、全エネルギーがこの値と 1e-5 Hartree 以内で一致することを確かめている。

## 名前について (2026-09-14)

この計算は**改名前 (VISTA) に実行した**ものです。改名 (VISTA → ADIT) に合わせて、
ファイル名と中の文字列を ADIT に直しました (いまの ADIT が生成するものと同じ形にするため)。
**計算そのもの (数値・収束・使ったソフトのバージョン) は実行したときのまま**です。
