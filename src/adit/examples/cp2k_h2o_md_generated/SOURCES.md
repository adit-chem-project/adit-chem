# examples/cp2k_h2o_md_generated —— CP2K の生成器の出力 (数ステップの MD) と実際に実行した結果

`examples/cp2k_h2o_generated` と同じ電子状態の条件で、NVT の MD を 5 ステップだけ実行したもの。CP2K 2026.2 (conda-forge)。

| 項目 | 内容 |
|---|---|
| 条件 | ENSEMBLE NVT、THERMOSTAT CSVR (TIMECON 100 fs)、TEMPERATURE 300 K、TIMESTEP 0.5 fs、STEPS 5、軌跡は毎ステップ (XYZ) |
| 実行 | `mpirun -np 1 cp2k.psmp -i cp2k.inp` (OpenMP 1 スレッド)。2026-09-12 |
| 結果 | 正常終了。所要時間 10.2 s。`adit-1.ener` の保存量 (Cons Qty) は -17.217974062 → -17.217957779 Hartree (5 ステップで 1.6e-5 Hartree の変化)。温度は 300 K の初速から 104.9 K へ |
| 出力 | `adit-pos-1.xyz` (軌跡、6 フレーム)、`adit-1.ener` (時刻・運動エネルギー・温度・ポテンシャルエネルギー・保存量)、`adit-1.restart` |

3 原子の分子で 2.5 fs しか実行していないので、温度は初速の分配がまだ振動へ移る途中の値である。値は「走る」ことの確認用。

## 名前について (2026-09-14)

この計算は**改名前 (VISTA) に実行した**ものです。改名 (VISTA → ADIT) に合わせて、
ファイル名と中の文字列を ADIT に直しました (いまの ADIT が生成するものと同じ形にするため)。
**計算そのもの (数値・収束・使ったソフトのバージョン) は実行したときのまま**です。
