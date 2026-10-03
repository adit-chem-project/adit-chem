# examples/plumed_lammps_cu_generated —— LAMMPS の MD に PLUMED をつないだ例 (2026-09-13 に実行)

`examples/lammps_cu_nvt_generated` と同じ Cu の NVT MD に、**利用者が書いた PLUMED の入力**を足したもの。
**ADIT は集合変数もバイアスも作りません。**`plumed.dat` の中身はこの例のために人が書いたものです。

| 項目 | 値 |
|---|---|
| 構造・力場 | 面心立方 Cu (`examples/lammps_cu` の EAM。出典はそのディレクトリの `SOURCES.md`) |
| 計算 | NVT MD、300 ステップ、刻み 1 fs、300 K、Nosé–Hoover (時定数 100 fs)、50 ステップごとに記録 |
| PLUMED | `d1: DISTANCE ATOMS=1,2` と `cn: COORDINATION GROUPA=1-32 GROUPB=1-32 R_0=3.0`、`PRINT ... FILE=COLVAR STRIDE=50` |
| つなぎ方 | `in.lammps` の `run` の直前に `fix adit_plumed all plumed plumedfile plumed.dat outfile plumed.log` |
| 実行 | LAMMPS 22 Jul 2025 (conda-forge、PLUMED パッケージ入り)、PLUMED 2.10、1 プロセス・1 スレッド |

## 単位に注意 (実際に実行して確かめた)

LAMMPS の `metal` 単位では距離は Å ですが、**PLUMED の既定の単位は nm** です。
`COLVAR` の `d1` は 0.2556 前後で、これは 2.556 Å (Cu の最近接距離) を nm で書いた値です。
PLUMED の入力に `UNITS` を書けば変えられます。この違いは ADIT が換算せず、README.txt にも注意として入れています。

## ファイル

| ファイル | 中身 |
|---|---|
| `plumed.dat` | 利用者が書いた PLUMED の入力 (ADIT がそのまま写したもの) |
| `COLVAR` | PLUMED が書いた集合変数の記録 (先頭が `#! FIELDS`)。ADIT の解析はこの列の名前と値を要約に入れます |
| `plumed.log` | PLUMED 自身の記録 (`fix plumed` の `outfile`) |
| `in.lammps` / `data.lammps` / `Cu_u3.eam` | LAMMPS の入力・構造・力場 |
| `log.lammps` / `output.log` / `traj.lammpstrj` / `final.data` | LAMMPS の出力 |

## 出典

- PLUMED: https://www.plumed.org/doc-v2.10/user-doc/html/index.html 、
  LAMMPS 側のつなぎ方: https://docs.lammps.org/fix_plumed.html
- PLUMED の文献は、使った機能ごとに `plumed.log` の最後に並びます (そこに出たものを引用してください)
- LAMMPS: https://docs.lammps.org/Intro_citing.html
