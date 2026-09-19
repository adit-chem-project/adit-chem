# examples/prep_stage1_dftb_stages —— 段階に分けた MD (最小化 → NVT → NVE) の生成したファイルと実際に実行した結果

`adit.stages.write_stages` (コマンド行なら `adit-gen spec.json out/ --stages stages.json`) で作ったもの。構造はテスト用の歪んだ水 (tests/conftest.py の water_spec)、
SK セットは mio-1-1、DFTB+ 25.1 (conda-forge、開発バージョンの表示)。`bash submit.sh` を 1 回実行した (OMP 1 スレッド、2026-09-19 に現在の生成器で作り直して実行し直した)。

| 段 | 条件 | 結果 |
|---|---|---|
| stage_01_min | 構造最適化 (最大 100 ステップ) | 正常終了。geom.out.gen |
| stage_02_nvt | NVT、Berendsen 300 K (時定数 20 fs)、0.5 fs × 40 ステップ | 開始の構造は stage_01 の geom.out.gen と 10 桁一致 (handoff.py が geometry.gen を書き換え) |
| stage_03_nve | NVE、0.5 fs × 20 ステップ、速度は stage_02 から | 最初の温度 314.6884 K = stage_02 の最後の温度 314.6884 K。全エネルギーの変化 1.5e-6 Hartree (20 ステップ) |

最初の生成では stage_03 が止まった: 速度 (Velocities) を与えたときに `Thermostat = None { InitialTemperature … }` を書くと、DFTB+ 25.1 は
「読まなかった項目」として止まる。生成器を直し (速度があるときは `Thermostat = None {}`)、生成し直して上の結果を得た。
