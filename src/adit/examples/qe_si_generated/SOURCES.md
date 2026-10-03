# examples/qe_si_generated —— 生成器が作った pw.x の入力と、実際に実行した結果

adit の生成器 (`codes/espresso.py`) が作ったダイヤモンド Si (2 原子、ASE の参照状態 a = 5.43 Å) の構造緩和と、
Quantum ESPRESSO 7.5 (conda-forge) の pw.x で実行した結果。`pseudo/` (UPF) はリポジトリに入れていない。

| 項目 | 内容 |
|---|---|
| 入力変数の出典 | https://www.quantum-espresso.org/Doc/INPUT_PW.html (calculation, nstep, forc_conv_thr, ecutwfc, ecutrho, occupations, smearing, degauss, nspin, tot_charge, tot_magnetization, conv_thr, electron_maxstep, cell_dofree, if_pos) |
| 書き出し | ASE 3.29 の `ase.io.espresso.write_espresso_in` |
| 公式の例 | QE の `PW/examples/example01` (Si の scf: ibrav 2, celldm 10.20, ecutwfc 18, conv_thr 1e-8。擬ポテンシャル Si.pz-vbc)。ここでは PBE の pslibrary を使い、ecutwfc / ecutrho は UPF に書かれた推奨値 (44 / 175 Ry) にした |
| 擬ポテンシャル | pslibrary 1.0.0 の `Si.pbe-n-rrkjus_psl.1.0.0.UPF` (https://pseudopotentials.quantum-espresso.org/upf_files/ から。pslibrary は GNU GPL: https://dalcorso.github.io/pslibrary/) |
| 結果 | JOB DONE。全エネルギー -22.83956699 Ry。bfgs は 1 SCF、0 ステップで収束 (対称な初期構造なので力が 0) |

QE 自体は GPL-2.0 (https://github.com/QEF/q-e)。引用要件は https://www.quantum-espresso.org/ を参照。
