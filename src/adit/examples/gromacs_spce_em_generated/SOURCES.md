# examples/gromacs_spce_em_generated —— GROMACS の生成器の出力 (エネルギー最小化) と実際に実行した結果

ADIT の生成器 (`codes/gromacs.py`) が作った、SPC/E 水 216 分子の箱のエネルギー最小化と、GROMACS 2026.3 (conda-forge `gromacs`、
混合精度、thread-MPI) で実行した結果。続きの NVT は `examples/gromacs_spce_nvt_generated`。

| 項目 | 内容 |
|---|---|
| 入力 | `examples/gromacs_spce/topol.top` と `conf.gro` (GROMACS 同梱の spc216.gro。出典は同ディレクトリの SOURCES.md) |
| mdp | integrator steep、nsteps 50、emtol 1000 kJ mol⁻¹ nm⁻¹ (画面の力の収束基準 1.0364 eV/Å を換算)、coulombtype PME、rcoulomb = rvdw = 0.85 nm |
| mdp の出典 | https://manual.gromacs.org/current/user-guide/mdp-options.html |
| 実行 | `gmx grompp -f grompp.mdp -c conf.gro -p topol.top -o adit.tpr` → `gmx mdrun -deffnm adit -ntmpi 1 -ntomp 1`。2026-09-12 |
| 結果 | 正常終了。所要時間 0.30 s。1 ステップで最大の力が 1000 kJ mol⁻¹ nm⁻¹ を下回った (最初から 859.84、原子 103)。ポテンシャルエネルギー -9953.6836 kJ/mol |

カットオフを 0.9 nm にすると grompp が「カットオフが箱の短い辺の半分より長い」で止まった。最小化では近傍リストの半径が
カットオフの 1.05 倍になり (マニュアルの verlet-buffer-tolerance)、0.9 × 1.05 = 0.945 nm が箱 1.862 nm の半分 0.931 nm を超えるため。
生成器の検査は、この事実をそのまま判定に使う (`tests/test_gromacs.py::test_things_grompp_would_reject`)。

## 名前について (2026-09-14)

この計算は**改名前 (VISTA) に実行した**ものです。改名 (VISTA → ADIT) に合わせて、
ファイル名と中の文字列を ADIT に直しました (いまの ADIT が生成するものと同じ形にするため)。
**計算そのもの (数値・収束・使ったソフトのバージョン) は実行したときのまま**です。
