# examples/abinit_si_generated —— ABINIT でダイヤモンド構造 Si の SCF 一点計算 (2026-09-13 に実行)

ADIT が生成した ABINIT の入力と、それを実際に実行した結果 (テキストのファイルだけ)。

| 項目 | 値 |
|---|---|
| 構造 | ダイヤモンド構造の Si、格子定数 5.43 Å の基本格子 (2 原子)。`ase.build.bulk` で作った |
| 条件 | `ecut 8 Hartree`、`ngkpt 4 4 4`(シフトなし)、`toldfe 1e-8`、`nstep 30`。汎関数は擬ポテンシャルのもの (`ixc` を書いていない) |
| 実行 | ABINIT 10.0.3 (conda-forge)、1 プロセス・1 スレッド |
| 結果 | `etotal = -8.8652901202 Hartree` (= -241.2368 eV、Si 2 原子) |

**`ecut 8 Hartree` と `4×4×4` は書式を見せるための値で、収束させた条件ではありません。**
実際の計算では、`adit-gen --scan method.ecut_ha=...` などで収束を確かめてください。

## 擬ポテンシャルは入っていません

この計算には ABINIT のテスト用擬ポテンシャル `14si.pspnc`(Troullier–Martins、ABINIT のリポジトリの
`tests/Pspdir/PseudosTM_pwteter/14si.pspnc`)を使いました。**ADIT は擬ポテンシャルを同梱しないので、
このディレクトリにも置いていません。**同じ計算を再現するには、そのファイルを入手して
`input.abi` と同じディレクトリに置いてください (入力の `pp_dirpath "./"` が同じ場所を指しています)。

    https://raw.githubusercontent.com/abinit/abinit/master/tests/Pspdir/PseudosTM_pwteter/14si.pspnc

`spec.json` の `method.pseudos` には、生成したときに指定したファイルの場所 (`/home/<user>/pseudos/14si.pspnc`) が
残っています。別の場所に置いたときは、その値を書き換えてから生成し直してください。

## ファイル

| ファイル | 中身 |
|---|---|
| `input.abi` | ADIT が書いた ABINIT の入力 (構造は Bohr 単位の `rprim` と `xred`) |
| `input.abo` | ABINIT の主出力。末尾の変数の表から、ADIT は最後の構造 (`acell`・`rprim`・`xred`) を読みます |
| `output.log` | 標準出力。ADIT はここから段ごとの `Total energy (etotal) [Ha]` を読みます |
| `code_version.txt` | `submit.sh` が残した ABINIT のバージョン |

ABINIT が書く大きなファイル (`*_GSR.nc`、`*_WFK`、`*_DEN` など) は、この例には含めていません (実行すると作られます)。

## 出典

- 入力変数の説明: https://docs.abinit.org/variables/ (ABINIT 10.0.3 で実際に実行して確かめました)
- ABINIT の文献: Gonze ほか, Comput. Phys. Commun. 248, 107042 (2020), doi:10.1016/j.cpc.2019.107042
- 擬ポテンシャル: 上の URL の ABINIT 配布物に含まれるテスト用データ
