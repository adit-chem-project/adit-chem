# examples/prep_stage1_cp2k_dftu_h2o —— CP2K の DFT+U (KIND/DFT_PLUS_U) の生成したファイルと実際に実行した結果

tests/conftest.py の歪んだ水、PBE、DZVP-MOLOPT-SR-GTH / GTH-PBE、CUTOFF 280 Ry、REL_CUTOFF 40 Ry、箱 8 Å (MT)。
O の 2p に U = 2.0 eV (KIND O の &DFT_PLUS_U で L 1、U_MINUS_J [eV] 2.0。DFT 節に PLUS_U_METHOD MULLIKEN)。仕組みの確認用の値。CP2K 2026.2、2026-09-12。

| 項目 | 値 (output.log) |
|---|---|
| 読んだ U | "DFT+U| Method MULLIKEN"、"U(eff) = (U - J) value in [eV]: 2.000" |
| DFT+U のエネルギー | 0.04085351432522 Hartree |
| 全エネルギー | -17.084127263102282 Hartree |

## 名前について (2026-09-14)

この計算は**改名前 (VISTA) に実行した**ものです。改名 (VISTA → ADIT) に合わせて、
ファイル名と中の文字列を ADIT に直しました (いまの ADIT が生成するものと同じ形にするため)。
**計算そのもの (数値・収束・使ったソフトのバージョン) は実行したときのまま**です。
