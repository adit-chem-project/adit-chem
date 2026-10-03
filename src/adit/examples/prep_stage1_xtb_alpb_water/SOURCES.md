# examples/prep_stage1_xtb_alpb_water —— xtb の溶媒 (ALPB、水) の生成したファイルと実際に実行した結果

XtbMethod(solvation="alpb", solvent="water") の一点計算。submit.sh の xtb の行に `--alpb water` が入る。xtb 6.7.1 (conda-forge)、2026-09-12。

| 項目 | 値 |
|---|---|
| 溶媒の模型 (output.log) | ALPB、Solvent water |
| 溶媒和の自由エネルギー Gsolv | -0.029383426097 Eh |
| 全エネルギー | -4.997357957784 Eh (同じ構造の気相は -4.977864911515 Eh) |

溶媒の名前の一覧 (codes/xtb.py の SOLVENTS) は、xtb --help と公式文書の表にある 35 通りの名前を、手法 (GFN1 / GFN2 / GFN-FF) × 模型 (ALPB / GBSA) ごとに
1 回ずつ実行して、xtb が止まらなかったものだけにした。
