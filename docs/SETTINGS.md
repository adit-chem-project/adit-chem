# 環境設定

計算コードの置き場所と、どこで実行するか (プロファイル) を決めるファイルの説明です。
解析の条件は [解析の詳細](ANALYSIS.md) にあります。

## 環境設定 (cluster.toml)

```toml
sk_root = "/home/<ユーザー名>/slakos"   # Slater-Koster パラメータの置き場所 (直下に mio-1-1/ など)
pseudo_root = "/home/<ユーザー名>/pseudo" # Quantum ESPRESSO の擬ポテンシャル (直下に <セット名>/*.UPF)
cp2k_data = ""                       # CP2K の data ディレクトリ (BASIS_MOLOPT など)。空なら CP2K_DATA_DIR、cp2k の隣の share/cp2k/data の順に探します
templates_dir = ""                   # 研究室の雛形 (構造を除いた spec.json) の置き場所。空なら、このファイルと同じ場所の templates/ だけを見ます
default_profile = "local"
enable_run = true                    # false にすると、この PC では計算を実行しない設定になります (生成と解析はできます)
language = "ja"                      # "ja" または "en" (環境変数 ADIT_LANG が優先)
theme = "auto"                       # "auto" / "light" / "dark" (環境変数 ADIT_THEME が優先)
window_frame = "auto"                # ウィンドウの枠。"auto" (Linux では ADIT が描く) / "custom" (ADIT が描く) / "native" (OS に任せる)。環境変数 ADIT_FRAME が優先

[profiles.local]
kind = "direct"
description = "この PC で bash submit.sh を実行"
```

### 環境変数

環境設定ファイルより優先されます。一時的に切り替えたいときに使います。

| 環境変数 | 意味 |
|---|---|
| `ADIT_CONFIG` | 読む環境設定ファイルのパス (既定は `~/.config/adit/cluster.toml`。Windows は `%APPDATA%\adit\cluster.toml`) |
| `ADIT_LANG` | 表示言語 `ja` / `en` |
| `ADIT_THEME` | 画面のテーマ `auto` / `light` / `dark` |
| `ADIT_FRAME` | ウィンドウの枠 `auto` / `custom` / `native` |
| `ADIT_SHELL` | ワークスペースのターミナルで使うシェル (例 `ADIT_SHELL="bash --norc"`) |

### クラスタで実行する場合

環境設定ファイルに `kind = "pbs"` または `kind = "slurm"` のプロファイルを追加すると、`submit.sh` にそのスケジューラのヘッダが付きます。
テンプレートが書くのは汎用の項目 (ジョブ名、ノード数、コア数、MPI プロセス数、OpenMP スレッド数、制限時間) だけで、
サイト固有の値はすべて環境設定から埋めます (既定は空)。
下の例の `<...>` はそのままでは動かないので、クラスタの管理者か研究室の先輩に確かめて置き換えてください。置き換え忘れがあると、生成の前の検証で「仮の値のままです」と表示されます。
使わない行は、行頭に `#` を付けて無効にしたままにします。

`host` と `remote_dir` を書いておくと、生成した `transfer_and_submit.sh` に、そのまま貼れる `rsync` と `ssh`、投入コマンドが書き出されます (**ADIT は実行しません**)。
`cores_max` / `nodes_max` / `walltime_max` を書いておくと、画面や `spec.json` の runtime がそれを超えたとき、「仮の値のままです」と同じ形で生成の前に止まります (既定は無し = 照合しない)。

```toml
[profiles.remote]
kind = "pbs"                        # または "slurm"
description = "所属のクラスタ"
select_extra = ""                   # PBS の select 行の末尾に付ける文字列 (サイトが要求する場合。例 ":jobtype=core")
header_extra = ["#PBS -q <キュー名>"]  # ヘッダに追加する行 (queue / partition / account など。Slurm なら "#SBATCH --partition=<パーティション名>")
submit_command = "qsub"             # 投入コマンド名 (README.txt に書くだけで、ADIT は実行しません)
status_command = "qstat -u $USER"   # 状態確認のコマンド (同上)
host = "<クラスタのホスト名>"          # 転送と投入のコマンドに使います (空でも生成できます。その場合は <...> のまま出ます)
user = "<ログイン名>"                 # 空なら手元のログイン名を使う想定で、host だけを書きます
remote_dir = "<クラスタでの作業ディレクトリ>"  # 転送先
cores_max = 0                       # そのキューで使えるコア数の上限 (ノード数 × ノードあたりのコア数と照合)。0 なら照合しません
nodes_max = 0                       # ノード数の上限。0 なら照合しません
walltime_max = ""                   # 制限時間の上限 (HH:MM:SS)。空なら照合しません。超えていると生成の前の検証で止まります
modules = []                        # どの計算コードでも module load するもの (例 ["intel-mpi"])。下の code_modules の前に読み込みます

[profiles.remote.code_modules]   # 計算コードごとに module load するもの (クラスタで module avail と打つと一覧が出ます)
dftbplus = ["<DFTB+ の module 名>"]

# vasp = ["<MPI の module 名>", "<数値計算ライブラリの module 名>"]

# [profiles.remote.env]          # submit.sh の先頭で export する環境変数 (VASP を使うときだけ。# を外して使います)

# VASP_PP_PATH = "<POTCAR ライブラリの親ディレクトリ>"

[profiles.remote.commands]       # 計算コードごとの実行コマンド。{mpiprocs} {omp_threads} {binary} が埋められます
dftbplus = "dftb+"

# vasp = "mpirun -np {mpiprocs} <VASP の bin>/vasp_{binary}"

# orca = "<ORCA を展開したディレクトリ>/orca"   # ORCA は絶対パスで呼び、mpirun は付けません
```

生成したあとの手順は次の 3 つです (`transfer_and_submit.sh` には、この 3 つが `host` と `remote_dir` を埋めた形で書かれています)。

```bash
rsync -av <生成したディレクトリ> <クラスタ>:<作業ディレクトリ>/       # 1. 転送 (scp -r でも同じです)
ssh <クラスタ> 'cd <作業ディレクトリ>/<名前> && qsub submit.sh'     # 2. 投入 (Slurm なら sbatch)
ssh <クラスタ> 'qstat -u $USER'                                      # 3. 状態確認 (Slurm なら squeue)
```

`#!/bin/sh` で `module` コマンドが定義されていない環境があるため、`submit.sh` は `module` を使う前に `/etc/profile` を読み、それでも見つからなければ理由を表示して止まります。

クラスタ向けのプロファイルでは `check_remote.sh` も生成されます。手元の PC で `bash check_remote.sh` と打つと、`ssh -o BatchMode=yes` でのログイン、`module load`、実行ファイルの有無 (`command -v`)、
`remote_dir` に書けるか、Slurm なら `sbatch --test-only submit.sh`、PBS なら `qstat -Q <キュー名>` を順に試し、1 行ずつ OK / NG を出します。**投入はしません** (`host` と `remote_dir` が要ります)。

## 各計算コードで用意するもの

| 計算コード | 本体 | パラメータ |
|---|---|---|
| DFTB+ | conda-forge の `dftbplus` など。`dftb+` を PATH に | Slater-Koster パラメータを https://dftb.org/parameters/download.html から取得し、`sk_root` の下に置きます (CC BY-SA 4.0)。必要な元素ペアの skf ファイルと LICENSE / README が出力ディレクトリの `skf/` にコピーされます |
| VASP | ライセンスを持つ利用者が用意します。実行パスを `commands.vasp` に | POTCAR ライブラリを `potpaw_PBE/<名前>/POTCAR` の階層で置き、親ディレクトリを `env.VASP_PP_PATH` に。**POTCAR 自体は出力に含めません** (`potcar.spec` と `make_potcar.sh` を書き、実行時に連結します) |
| xtb | conda-forge の `xtb` を PATH に (LGPL-3.0) | 不要 (計算手法に内蔵) |
| Quantum ESPRESSO | conda-forge の `qe` などで `pw.x` を PATH に | UPF ファイルを `pseudo_root/<セット名>/` に置きます (pslibrary、SSSP https://www.materialscloud.org/discover/sssp など。ライセンスは配布元で確認してください)。必要な元素の UPF と付属文書が `pseudo/` にコピーされます |
| ORCA | 登録して入手し、実行パスを `commands.orca` に (再配布不可) | 不要 |
| CP2K | conda-forge の `cp2k` などで `cp2k.psmp` を PATH に (GPL-2.0-or-later) | CP2K に付いている data ディレクトリ (BASIS_MOLOPT、GTH_POTENTIALS など)。場所は環境設定の `cp2k_data` に書きます (空なら環境変数 `CP2K_DATA_DIR`、PATH にある cp2k の隣の `share/cp2k/data` の順に探します)。使う項目だけが `BASIS_adit` / `POTENTIAL_adit` に写されます |
| LAMMPS | conda-forge の `lammps` などで `lmp` を PATH に (GPL-2.0) | 力場のファイル (EAM、ReaxFF の ffield、機械学習ポテンシャルのモデル) か、外部で作った data ファイル。ADIT は力場の係数を決めません |
| GROMACS | conda-forge の `gromacs` などで `gmx` を PATH に (LGPL-2.1) | 外部で作ったトポロジー (.top / .itp) と構造 (.gro / .pdb)。CHARMM-GUI、acpype、pdb2gmx などで作ります |
| 機械学習ポテンシャル (MACE・CHGNet) | 生成した `run_mlip.py` を実行する環境に `pip install ase mace-torch` (または `chgnet`)。PyTorch が入るので数 GB になります。**ADIT 自身は使いません** | 学習済みモデル。空欄ならパッケージの既定のモデル (初回に自動でダウンロードされます)。モデルのファイルを指定すると出力ディレクトリに写します。ライセンスはモデルごとに違うので配布元で確認してください |

### CP2K・LAMMPS・GROMACS の画面

計算コードのプルダウンで選ぶと、そのコードの欄が出ます。「必須」の印が付いた青字のラベルは必須の欄です (ラベルにカーソルを合わせると説明が出ます)。

- **CP2K**: 汎関数、カットオフと相対カットオフ [Ry] は既定を持たないので、自分で決めます (CP2K のマニュアルの「CUTOFF と REL_CUTOFF の収束」の手順)。元素ごとの基底と擬ポテンシャルは、data ディレクトリのファイルから読んだ候補がプルダウンに並びます。候補が 1 つだけの元素は空欄のままでそれを使います。分子や表面のように周期でない方向があるときは「ポアソン方程式の解き方」を選び、箱の無い分子なら「分子の箱の一辺」を入れます。画面に無い設定は「追加の行 (節ごと)」に `[FORCE_EVAL/DFT/SCF]` のような見出しを書き、その下に行を書きます。
- **LAMMPS**: 単位系 (metal / real) と pair_style、pair_coeff を書き、力場のファイルを「写すファイル」に入れます (入力の中ではファイル名だけで書きます)。data ファイルを空にすると構造から `data.lammps` を書きます (atomic か charge 形式)。外部の data ファイルを使うときは、型番号の元素を順に入れ、構造の群にも同じ系を読み込みます。k 点の群は出ません。
- **GROMACS**: トポロジー (.top) と構造のファイル (.gro / .pdb) を指定します。**構造の正本は構造のファイルです。**「構造の群にも読み込む」を押すと、構造の群の作り方が「ファイル」になって同じファイルを読み込みます (原子数の確認と 3D 表示に使います。ウェブ版では「構造の欄にもこのファイルを使う」に印を付けます)。NPT では「圧力浴 (pcoupl)」と等温圧縮率が必須です (GROMACS のマニュアルの水の例は 4.5e-5 /bar)。NVT → NPT → 本計算のように段階を重ねるときは、前の段階の `adit.cpt` を「前の段階の .cpt」に指定します。k 点の群は出ません。

計算の種類のうち、そのコードの生成器に無いもの (LAMMPS と GROMACS の振動解析とバンド計算、CP2K のバンド計算) は、選ぶと「生成できません」の欄に理由が出ます。ほかのコードと同じ扱いです。

ORCA は本体が再配布できないため、入力ファイルの内容だけをテストしています。VASP はクラスタで、DFTB+、xtb、pw.x、CP2K、LAMMPS、GROMACS、OpenMM、Psi4、ABINIT、PLUMED は手元の PC で実行して確認しています。

## 出力ディレクトリ

```
<出力ディレクトリ>/
  submit.sh      実行スクリプト (ローカル: bash submit.sh / クラスタ: qsub または sbatch submit.sh)
  spec.json      計算設定一式。ADIT で開くと復元できます
  README.txt     実行手順と出力ファイルの見方
  analyze.py     解析スクリプト (analysis/ に図と要約を書きます)。解析に対応した計算コードのときだけ
  transfer_and_submit.sh   転送と投入のコマンド。プロファイルが pbs / slurm のときだけ
  check_remote.sh          投入前の疎通確認 (ssh、module、実行ファイル、作業ディレクトリ、キュー)。同上。投入はしません
  dftb_in.hsd, geometry.gen, skf/       DFTB+ の場合
  INCAR, POSCAR, KPOINTS, potcar.spec, make_potcar.sh   VASP の場合
  struct.xyz, xtb.inp                   xtb の場合
  pw.in, pseudo/                        Quantum ESPRESSO の場合
  orca.inp                              ORCA の場合
```

