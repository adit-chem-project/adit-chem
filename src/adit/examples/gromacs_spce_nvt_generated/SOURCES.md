# examples/gromacs_spce_nvt_generated —— GROMACS の生成器の出力 (NVT) と実際に実行した結果

`examples/gromacs_spce_em_generated` の最後の構造 (`adit.gro`) を構造のファイルに指定して作った NVT。GROMACS 2026.3 (conda-forge)。

| 項目 | 内容 |
|---|---|
| 入力 | トポロジーは `examples/gromacs_spce/topol.top`、構造は `examples/gromacs_spce_em_generated/adit.gro` |
| mdp | integrator md、dt 0.002 ps、nsteps 250 (0.5 ps)、tcoupl V-rescale (画面の熱浴 csvr)、tau-t 0.1 ps、ref-t 300 K、gen-vel yes (gen-seed -1)、PME、0.85 nm |
| 実行 | grompp → `gmx mdrun -deffnm adit -ntmpi 1 -ntomp 1`。2026-09-12 |
| 結果 | 正常終了。所要時間 0.41 s。250 ステップの平均: 温度 307.94 K、ポテンシャルエネルギー -9832.91 kJ/mol、圧力 140.9 bar (adit.log の AVERAGES) |
| 出力 | `adit.xtc` (軌跡、25 ステップごと)、`adit.edr` (エネルギー)、`adit.gro` (最後の構造)、`adit.cpt` (次の段へ渡す続きの情報) |

gen-seed が -1 (GROMACS の既定。毎回ちがう乱数) なので、実行し直すと値は一致しない。0.5 ps は平衡化には短く、値は「走る」ことの確認用。

## 名前について (2026-09-14)

この計算は**改名前 (VISTA) に実行した**ものです。改名 (VISTA → ADIT) に合わせて、
ファイル名と中の文字列を ADIT に直しました (いまの ADIT が生成するものと同じ形にするため)。
**計算そのもの (数値・収束・使ったソフトのバージョン) は実行したときのまま**です。
