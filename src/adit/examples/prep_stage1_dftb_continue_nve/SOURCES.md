# examples/prep_stage1_dftb_continue_nve —— 前の計算の続き (DFTB+ の MD → NVE) の生成したファイルと実際に実行した結果

`adit.continuation.continue_from("examples/dftb_md_water_generated", 条件)` で作ったもの (コマンド行なら `adit-gen --continue-from …`)。
条件は元の spec.json の MD を NVE・10 ステップ・毎ステップ出力にしたもの。位置と速度は元の geo_end.xyz の最後のフレーム
(速度は Å/ps → Å/fs に直して spec.json の structure.velocities、dftb_in.hsd には Velocities [AA/ps] で書く)。DFTB+ 25.1、2026-09-12。

| 項目 | 値 |
|---|---|
| 元の MD の最後の温度 (examples/dftb_md_water_generated/md.out) | 5885.4585 K |
| 続きの最初の温度 (md.out) | 5885.4585 K |

2 つが一致したので、速度は単位も含めて引き継がれている。温度が高いのは元の例 (歪んだ水の MD) の値で、続きの確認とは関係しない。
