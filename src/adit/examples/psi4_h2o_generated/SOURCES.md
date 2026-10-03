# examples/psi4_h2o_generated —— Psi4 で水分子の構造最適化 (2026-09-13 に実行)

ADIT が生成した Psi4 の入力と、それを実際に実行した結果。

| 項目 | 値 |
|---|---|
| 構造 | ASE の `molecule("H2O")` の座標 (非周期、電荷 0、多重度 1) |
| 条件 | 手法 `scf`、基底 `cc-pVDZ`、`geom_maxiter 20`。参照関数は書いていないので Psi4 の既定 (RHF) |
| 実行 | Psi4 1.11 (conda-forge)、1 スレッド (`psi4 -i input.dat -o output.log -n 1`) |
| 結果 | 4 段で収束、`energy_hartree = -76.0270327837`、最適化後の O-H = 0.946 Å、H-O-H = 104.6° |

**手法と基底は書式を見せるための値で、ADIT の推奨ではありません。**
O-H 0.946 Å は RHF/cc-pVDZ の水でよく知られた値で、測定 (0.958 Å) より短くなります。

## ADIT がしていること・していないこと

- 入力 (`input.dat`) は psithon (Psi4 の入力言語)。分子の並びは ADIT の構造そのままで、
  `no_reorient`・`no_com`・`symmetry c1` を付けて Psi4 が座標を回転・並進しないようにしています
- 入力の末尾の数行が `results.json` と `final.xyz` を書きます。解析はその 2 つと、
  `output.log` の「Optimization Summary」の表 (段ごとのエネルギー) を読みます
- 手法・基底・参照関数・収束の条件は**利用者が決めます**。ADIT は選びません

## ファイル

| ファイル | 中身 |
|---|---|
| `input.dat` | Psi4 の入力 |
| `output.log` | Psi4 の出力 (`-o` で指定)。段ごとのエネルギーの表と、`adit-psi4:` の 2 行が入ります |
| `stdout.log` | 画面に出る分 (optking の途中経過)。Psi4 が止まったときの理由もここに出ます |
| `results.json` / `final.xyz` | 入力の末尾が書いた結果と、最適化後の構造 |

## 出典

- Psi4 の入門と入力の書き方: https://psicode.org/psi4manual/master/index.html 、
  https://psicode.org/psi4manual/master/psithoninput.html (Psi4 1.11 で実際に実行して確かめました)
- Psi4 の文献: Smith ほか, J. Chem. Phys. 152, 184108 (2020), doi:10.1063/5.0006002
