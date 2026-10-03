# examples/xtb_water_generated —— 生成器が作った xtb の入力と、実際に実行した結果

adit の生成器 (`codes/xtb.py`) が作った水分子の構造最適化 (GFN2-xTB) と、xtb 6.7.1 (conda-forge) で実行した結果。
xtb にはパラメータファイルが無いので、生成したファイルはこれで完結する。

| 項目 | 内容 |
|---|---|
| コマンド行の引数の出典 | https://xtb-docs.readthedocs.io/en/latest/commandline.html と `xtb --help` (6.7.1)。`--gfn`, `--chrg`, `--uhf`, `--acc`, `--etemp`, `--opt [LEVEL]`, `--input`, `--parallel`, `--json` |
| xcontrol の出典 | https://xtb-docs.readthedocs.io/en/latest/xcontrol.html。`$scc maxiterations`、`$opt maxcycle`、`$fix atoms:`(1 始まり)、`$end` |
| xtb のライセンス | LGPL-3.0 (https://github.com/grimme-lab/xtb)。本体は conda-forge の `xtb` から |
| 結果 | 9 反復で収束。全エネルギー -5.070544 Eh (output.log)。最適化後の構造は xtbopt.xyz |

引用要件は xtb の文書に従う (GFN2-xTB: Bannwarth, Ehlert, Grimme, J. Chem. Theory Comput. 2019, 15, 1652)。
