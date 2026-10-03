# examples/prep_stage1_qe_dftu_si —— pw.x の DFT+U (HUBBARD カード) と starting_magnetization の生成したファイルと実際に実行した結果

Si (ダイヤモンド構造、pslibrary の Si.pbe-n-rrkjus_psl.1.0.0.UPF)、ecutwfc 20 Ry、k 点 2×2×2、nspin = 2、smearing (gaussian 0.01 Ry)、
starting_magnetization(Si) = 0.2、HUBBARD {atomic} で U Si-3p 1.0 eV。仕組みの確認用の値で、物理的な意味は持たせていない。QE 7.5 (conda-forge)、2026-09-12。

| 項目 | 値 (output.log) |
|---|---|
| 読んだ U | "Hubbard projectors: atomic"、"U(Si-3p) = 1.0000" (Dudarev) |
| 磁化 | 0.72 Bohr mag/cell から始まり 0.00 に (Si は非磁性) |
| 全エネルギー | -22.57798452 Ry (Hubbard energy 0.07069816 Ry)、6 回で収束 |

tmp/ は pw.x の途中のファイル (波動関数と電荷密度)。
