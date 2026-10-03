# examples/openmm_spce_nvt_generated —— OpenMM で SPC/E 水 216 分子の NVT MD (2026-09-13 に実行)

ADIT が生成した OpenMM の計算ディレクトリと、それを実際に実行した結果。**ADIT は力場・原子型・電荷を作らない。**
トポロジーと座標は `examples/gromacs_spce/`(GROMACS 同梱の SPC/E の箱)をそのまま写したもので、出典とライセンスは
そのディレクトリの `SOURCES.md` にある。

| 項目 | 値 |
|---|---|
| 計算 | NVT MD、200 ステップ、刻み 2 fs、300 K、Langevin (時定数 100 fs → 衝突頻度 1/100 fs⁻¹)、50 ステップごとに記録 |
| 非結合 | PME、カットオフ 0.8 nm、`constraints = HBonds`、`rigidWater = True` (OpenMM の既定) |
| 実行 | openmm 8.4.0.dev-4768436 (conda-forge)、Platform = CPU、`OPENMM_CPU_THREADS=1`。乱数の種 12345 を Langevin に渡すので、同じ環境で実行し直すと同じ値になる (生成し直して実行しても `results.json` が一致することを確かめた)。入力の書式は 8.6.1 でも確かめている (`openmm-plumed` を入れたときに、この env の openmm が 8.6.1 から 8.4.0.dev へ下がった) |
| 力場 | GROMACS 同梱の `oplsaa.ff/forcefield.itp` と `oplsaa.ff/spce.itp` を `include_dir` から読む (同梱しない) |
| 結果 | ポテンシャルエネルギー -9980 kJ/mol (水 1 分子あたり -46.2 kJ/mol)、最後の温度 295.3 K、自由度 1293 |

自由度 1293 = 648 原子 × 3 − 剛体水の拘束 648 − 重心運動の除去 3。1 分子あたりのポテンシャルエネルギーは
SPC/E の文献値 (Berendsen, Grigera, Straatsma, J. Phys. Chem. 91, 6269 (1987)) の水の配置エネルギーと同じ桁で、
200 ステップの値なので平衡値としては扱わない。

## ファイル

| ファイル | 中身 |
|---|---|
| `run_openmm.py` / `openmm_settings.json` | ADIT が書いたスクリプトと設定の値 (計算の種類によらず同じスクリプト) |
| `topol.top` / `conf.gro` | 写したトポロジーと座標 (`examples/gromacs_spce/`) |
| `md.log` | OpenMM の `StateDataReporter` の CSV (ステップ・時刻・エネルギー・温度・体積・密度) |
| `results.json` / `final.pdb` | 最後のエネルギーと温度、最後の構造 |
| `trajectory.dcd` | OpenMM の `DCDReporter` の軌跡 (VMD・OVITO・TRAVIS 用。ADIT は読まない) |
| `trajectory.extxyz` | `write_xyz_trajectory` を有効にしたときに書くテキストの軌跡 (ADIT の RDF・MSD はこちらを読む) |
| `analysis/` | `adit-analyze . --rdf --msd O` の出力 (4 フレームしかないので、値は書式の確認用) |

## 出典

- OpenMM の使い方と API: https://docs.openmm.org/latest/userguide/application/02_running_sims.html 、
  https://docs.openmm.org/latest/api-python/app.html (引数は openmm 8.6.1 の定義で確かめ、8.4.0.dev でも実際に実行した)
- OpenMM の文献: Eastman ほか, PLoS Comput. Biol. 13, e1005659 (2017), doi:10.1371/journal.pcbi.1005659
