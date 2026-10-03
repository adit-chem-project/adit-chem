# examples/lammps_cu_nvt_generated —— LAMMPS の生成器の出力と実際に実行した結果

ADIT の生成器 (`codes/lammps.py`) が作った Cu の NVT と、LAMMPS (22 Jul 2025 update 5、conda-forge `lammps`) で実行した結果。

| 項目 | 内容 |
|---|---|
| 構造 | ASE `bulk("Cu", "fcc", a=3.615, cubic=True).repeat((3, 3, 3))` (108 原子、周期、立方体 10.845 Å) |
| 相互作用 | `pair_style eam` / `pair_coeff * * Cu_u3.eam`。ポテンシャルは `examples/lammps_cu/Cu_u3.eam` (出典は同ディレクトリの SOURCES.md) |
| 条件 | units metal、NVT (Nosé-Hoover、`fix nvt`)、300 K、1 fs × 300 ステップ、熱浴の時定数 100 fs、軌跡は 50 ステップごと |
| コマンドの出典 | https://docs.lammps.org/ の units、read_data、pair_style eam、velocity、fix nvt、dump custom、dump_modify element / sort、write_data の各ページ |
| 実行 | `mpirun -np 1 lmp -in in.lammps -log log.lammps` (1 プロセス、OpenMP 1 スレッド、nice 19、メモリ上限 1.5 GB)。2026-09-12 |
| 結果 | 正常終了。所要時間 1.7 s (submit.sh 全体)。300 ステップ目: 温度 252.67 K、ポテンシャルエネルギー -379.02162 eV (108 原子の合計)、全エネルギー -375.52698 eV |
| 軌跡 | `traj.lammpstrj` (dump custom、列は id type element x y z、7 フレーム)。`final.data` は最後の構造 |

温度の値は、300 ステップ (0.3 ps) では熱浴の時定数 (0.1 ps) の 3 倍しか経っていないので、揺らいでいる途中の値である。
この例は「生成したファイルがエラーなく走る」ことを確かめるためのもので、物性の値として使うものではない。
