# 各 OS のインストール方法

ADIT を入れて、最初の計算を実行するまでの手順です。

## いちばん簡単な方法: 実行ファイルを使う (Windows / macOS)

**Python も WSL も要りません。**

1. [Releases](https://github.com/adit-chem-project/adit-chem/releases) を開きます
2. Windows なら `ADIT-windows-x64.zip`、macOS (Apple Silicon) なら `ADIT-arm64.zip` を
   ダウンロードします (Intel の Mac 用は配っていません。下の pip で入れてください)
3. 展開し、Windows は `ADIT.exe`、macOS は `ADIT.app` をダブルクリックします

署名を付けていないため、初回だけ警告が出ます。Windows は「詳細情報」→「実行」、
macOS は Finder で右クリック →「開く」を選んでください。

コマンドライン (`adit-gen` など) を使いたいときは、同じフォルダの `adit-cli.exe` (macOS は `adit-cli`) を呼びます。

同じ PC で計算まで実行したい場合は、次の節へ進んでください。

## 同じ PC で計算まで実行したい人へ (Windows、この順に進めてください)

Linux を使ったことがなくても、この節の手順を上から順に行えば、最初の入力ファイルを作って計算を実行するところまで進めます。
コマンドは **1 行ずつ**コピーしてターミナル (文字でコンピュータに指示を出す黒い画面) に貼り付け、Enter を押します。
`#` から行末までは説明なので、貼り付けなくてもかまいません。`<...>` は自分の値に置き換えます。

### 1. WSL (Windows の中で Linux を動かす仕組み) を入れる

1. スタートメニューで「PowerShell」を右クリックし、「管理者として実行」を選びます
2. 次の 1 行を入力して Enter を押し、終わったら PC を再起動します

   ```powershell
   wsl --install
   ```

3. 再起動後に「Ubuntu」の画面が開くので、Linux 用のユーザー名とパスワードを決めます (Windows のものとは別です。パスワードは入力しても画面に表示されません)

以後は、スタートメニューの「Ubuntu」で開く画面 (ターミナル) で作業します。

### 2. miniforge (conda) を入れる

conda は、Python や計算コードを「環境」という箱に分けて入れる仕組みです。miniforge はその配布版で、conda-forge (DFTB+、xtb、Quantum ESPRESSO を配布している場所) を既定で使います。

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh      # 質問には Enter と yes で答えます
```

終わったら、ターミナルを一度閉じて開き直します。行頭に `(base)` と出ていれば入っています。

> **venv ではなく conda を勧める理由**: ADIT は Python だけで動きますが、この PC で計算を実行するには DFTB+ などの計算コードが要ります。
> conda ならそれらを同じ手順で入れられます。WSL の Ubuntu には最初 `python` コマンドが無く `python3` だけですが、conda の環境を有効にすると `python` が使えます。
> 計算はクラスタでしか実行せず、手元では入力の生成だけをする場合は、`python3 -m venv adit-env` と `source adit-env/bin/activate` の venv でもかまいません。

### 3. ADIT と DFTB+ を入れる環境を作る

```bash
conda create -n adit -c conda-forge "python>=3.11" dftbplus rdkit font-ttf-noto-cjk git
conda activate adit                         # 行頭が (adit) に変わります。ターミナルを開くたびに打ちます
```

`font-ttf-noto-cjk` はデスクトップ版の画面で日本語を表示するためのフォントです (WSL には最初、日本語のフォントがありません)。
xtb と Quantum ESPRESSO は、DFTB+ と同じ環境に入れると DFTB+ のバージョンが古いものに下がることがあるため、別の環境に入れます (使うものだけで十分です)。

```bash
conda create -n xtb -c conda-forge xtb
conda create -n qe -c conda-forge qe
conda activate --stack xtb                   # adit の環境を有効にしたまま、xtb も使えるようにする例
```

### 4. ADIT を入れる

```bash
git clone https://github.com/adit-chem-project/adit-chem.git
cd adit-chem && pip install -e ".[gui,smiles,analysis]"
```

`[gui,smiles,analysis]` は、デスクトップ版の画面、SMILES と Draw、空間群の解析に要るパッケージをまとめて入れる指定です。

### 5. パラメータを入手して置く

計算コードによっては、元素ごとのパラメータのファイルが要ります。ADIT には入っていないので、配布元から入手します。
Windows のブラウザでダウンロードしたファイルは、WSL の中からは `/mnt/c/Users/<Windows のユーザー名>/Downloads/` に見えます。

**DFTB+ の Slater-Koster パラメータ** (元素の組ごとの `.skf` ファイル。ライセンスは CC BY-SA 4.0 です)

入手先: https://dftb.org/parameters/download.html (水や有機分子なら mio が入り口です。mio が扱う元素は H, C, N, O, S, P)

```bash
mkdir -p ~/slakos                                                   # 置き場所を作る (~ は自分のホームフォルダ)
tar -xf /mnt/c/Users/<Windows のユーザー名>/Downloads/mio-1-1.tar.xz -C ~/slakos
ls ~/slakos/mio-1-1                                                 # H-H.skf などが並べば置けています
```

展開後は次の形になります。セットのフォルダ (`mio-1-1/`) の中に `.skf` と `LICENSE`、`README` が直接入っている形です。

```
~/slakos/
  mio-1-1/
    H-H.skf  H-O.skf  O-H.skf  O-O.skf  C-C.skf  …  LICENSE  README
```

**Quantum ESPRESSO の擬ポテンシャル (UPF ファイル)**

入手先の例: SSSP https://www.materialscloud.org/discover/sssp (ライセンスは配布元で確認してください)。
展開して、`.UPF` ファイルが直接入ったフォルダを、たとえば `~/pseudo/SSSP_efficiency/` として置きます。フォルダの名前が ADIT の「擬ポテンシャルのセット」になります。

```
~/pseudo/
  SSSP_efficiency/
    Si.pbe-n-rrkjus_psl.1.0.0.UPF  O.pbe-n-kjpaw_psl.0.1.UPF  …
```

xtb は計算手法にパラメータが入っているので不要です。VASP の POTCAR と ORCA は、[環境設定](SETTINGS.md#各計算コードで用意するもの) の「各計算コードで用意するもの」を見てください。

### 6. 環境設定ファイルに置き場所を書く

はじめて `adit-web` か `adit-gen` を実行すると、環境設定ファイル `~/.config/adit/cluster.toml` が作られ、その場所が表示されます。
次のコマンドで開き (nano は端末で使う簡単な文字エディタです。Ctrl+O → Enter で保存、Ctrl+X で終了)、置き場所を書きます。

```bash
echo $HOME                                   # 自分のホームフォルダ (例 /home/taro) が出ます
nano ~/.config/adit/cluster.toml
```

```toml
sk_root = "/home/<ユーザー名>/slakos"        # 5 で作った置き場所。セットのフォルダ (mio-1-1) の 1 つ上
pseudo_root = "/home/<ユーザー名>/pseudo"    # Quantum ESPRESSO を使うとき
```

パスは必ず `"..."` で囲みます。Windows の書き方 (`C:\Users\...`) はそのままでは読めません。WSL の中では `C:\Users\taro` は `/mnt/c/Users/taro` です。

### 7. 最初の入力を作って実行する

```bash
adit-web --open                             # ブラウザが開かなければ、Windows のブラウザで http://127.0.0.1:8765/ を開きます
```

ブラウザの画面で、構造「プリセット」の H2O、計算コード DFTB+、Slater-Koster パラメータ mio-1-1、プロファイル local のまま「プレビュー」→「生成」と進み、ターミナルで `bash submit.sh` を実行したあと「解析へ進む」を押し、解析のページで「解析を実行」を押します。
ターミナルだけで行う場合は、画面で保存した計算設定 (spec.json) から次のように作れます。

```bash
adit-gen spec.json ~/adit_runs/water       # 入力ファイル一式を ~/adit_runs/water に書きます
cd ~/adit_runs/water && bash submit.sh      # 計算を実行します (記録は output.log)
python analyze.py                            # 図と要約を analysis/ に書きます
```

各ディレクトリの `README.txt` に、そのディレクトリでの手順が書かれています。

## 実行ファイル (.exe / .app) を自分で作る

配布している実行ファイルは、この手順で作っています。`v` で始まるタグ (例 `v0.1.0`) をプッシュすると、
GitHub Actions が Windows の `.exe` と macOS の `.app` を作り、Releases に添付します。
手元で作るときは、下の手順です。この設定は Linux で実際に作り、CLI が動くところまで確かめてあります (2026-09-13)。

### Windows

#### 1 Windows の PC で作る

PowerShell で、リポジトリを置いた場所に移動してから:

```powershell
py -3.12 -m venv build-env
build-env\Scripts\activate
pip install -e ".[gui,smiles,analysis]" pyinstaller
pyinstaller packaging\adit.spec --noconfirm
```

`dist\adit\` に**実行ファイルが 2 つ**出来ます。フォルダごと配ってください (中のファイルが要ります)。

| ファイル | 使いみち |
|---|---|
| `ADIT.exe` | 画面 (デスクトップ版)。黒い画面が出ない代わりに、**標準出力を持たないので CLI には使えません** |
| `adit-cli.exe` | CLI。`adit-cli.exe gen ...` / `analyze` / `report` / `convert` / `web` |

2 つに分けるのは Windows の決まりのためです。画面用の実行ファイル (GUI サブシステム) には標準出力が無く、
**1 つにすると CLI の表示が何も出ません** (実際に作って確かめました)。
1 つの .exe にまとめたいときは `packaging/adit.spec` の `ONEFILE = False` を `True` に変えます
(起動のたびに一時フォルダへ展開するので、開くまで数秒かかります)。

#### 2 出来た実行ファイルの確認 (作った人が必ず通す)

```powershell
dist\adit\ADIT.exe                             # 画面が開く。日本語が □ になっていないか
dist\adit\adit-cli.exe gen --list-samples      # サンプルの一覧が出るか (同梱の examples を読めているか)
dist\adit\adit-cli.exe gen --sample water_generated mine.json
dist\adit\adit-cli.exe gen mine.json out\run1  # 入力・submit.sh・README.txt が書けるか
dist\adit\adit-cli.exe analyze out\run1        # 実行していないディレクトリでは理由を言って止まるか (終了コード 2)
dist\adit\adit-cli.exe web --open              # ブラウザで http://127.0.0.1:8765 が開くか (--open が無いとブラウザは開きません)
```

**画面の日本語が □ (豆腐) になるとき**は、フォントが同梱されていません。conda の環境で
`conda install -c conda-forge font-ttf-noto-cjk` を入れてから作り直すか、`ADIT_FONT_DIR` に
`NotoSansCJKjp-VF.ttf` のあるフォルダを指定して作り直します。

#### 3 同梱するもの・しないもの

| 同梱する | 理由 |
|---|---|
| `adit/scripts/templates/*.j2`、`adit/web/templates/*.html` | submit.sh と画面の雛形 |
| `examples/*/spec.json` | `gen --list-samples` / `--sample` が読む |
| `fonts/NotoSansCJKjp-VF.ttf` (あれば) | 画面と図の日本語 |

| 同梱しない | 理由 |
|---|---|
| DFTB+・xtb・Quantum ESPRESSO などの計算コード | 配布条件が別。**Windows では計算を実行しない** (生成した入力を Linux のサーバーへ転送して使う) ので、そもそも要らない |
| Slater-Koster セット、擬ポテンシャル、POTCAR | ライセンス上、ADIT が配ってはいけない |

#### 4 分かっている制限

- **Windows では計算を実行できません。**ワークスペースのターミナル (PowerShell) から `ssh` でクラスタに入り、
  `transfer_and_submit.sh` のコマンドで転送して投入します。exe で出来るのは入力の生成と、手元にある出力の解析です
- **大きさは 300〜500 MB** (PySide6 と matplotlib と SciPy を同梱するため)。`ONEFILE = True` でも縮みません
- **署名していません。**SmartScreen が「発行元不明」と警告します。配るなら署名するか、受け取る人に
  「詳細情報 → 実行」を案内してください
- **RDKit を入れた環境で作ると SMILES から構造を作れます**。入れずに作ると、その欄だけが使えません
- Linux で確かめたのは `gen --list-samples` / `--sample` / 入力の生成 / `analyze` (未実行のディレクトリで終了コード 2) /
  `report` / `convert verify` の 6 つと、**画面が起動して窓を作るところまで** (`QT_QPA_PLATFORM=offscreen`)。
  **Windows で生成したファイルは未確認**です (Windows の実行ファイルは Windows 上でしか作れないため)。
  上の「2 出来た実行ファイルの確認」を必ず通してください
- 作ってみて分かったこと 4 つ (どれも spec に入れてあります):
  1. **ASE は形式ごとのモジュールを名前で動的に読み込む**ので、`hiddenimports` に `ase.io` を入れないと
     `.gen` の書き出しが `UnknownFileTypeError` で落ちる
  2. Linux では **conda の `libOpenGL.so.0` を同梱しないと画面が起動しない** (Windows の PySide6 は
     OpenGL の DLL を自分で持っているので、この処理は Linux でだけ働く)
  3. **画面用の実行ファイルは標準出力を持たない**ので、CLI 用に console 版 (`adit-cli.exe`) を別に作る
  4. 窓の大きさの既定 (1800x1000) が **1536x864 の画面からはみ出した**ので、画面の 95 % に収めるようにした


### macOS

#### 作る

GitHub の Actions から手で動かします。

`v` で始まるタグ (例 `v0.1.0a1`) をプッシュすると自動で作られ、Releases に `ADIT-arm64.zip` が添付されます。
手で動かすこともできます。

1. リポジトリの **Actions** → **macos-app** → **Run workflow**
2. Apple Silicon (arm64) の 1 つが走ります
3. 終わると成果物 (artifact) に `adit-macos-arm64` が出ます。中身は
   - `ADIT-arm64.zip` … 二重クリックで開く `ADIT.app` (`ditto` でまとめたもの。Finder の権限が保たれます)
   - `dist/adit/` … 画面を使わない人向けの一式 (`adit-cli` が入っています)

手元の Mac で作るなら:

```bash
python -m pip install -e ".[gui,smiles,analysis]" pyinstaller
pyinstaller packaging/adit.spec --noconfirm
open dist/ADIT.app          # 画面
./dist/adit/adit-cli gen --list-samples   # コマンド
```

#### 署名と公証をしていません

Apple の開発者登録 (有料) が要るため、**署名 (codesign) も公証 (notarization) もしていません。**
そのため、受け取った人が初めて開くときは次のどちらかが要ります。

- Finder で `ADIT.app` を**右クリック → 開く** → 出てくる確認で「開く」
- または `xattr -dr com.apple.quarantine /Applications/ADIT.app`

**この手順を配布物の案内に必ず書いてください。**書かないと「壊れているから開けません」と表示され、
利用者は原因が分かりません。

署名するときは、Apple Developer Program に登録したうえで

```bash
codesign --deep --force --options runtime --sign "Developer ID Application: <名前> (<チーム ID>)" dist/ADIT.app
xcrun notarytool submit dist/ADIT-arm64.zip --apple-id <id> --team-id <team> --password <app 用パスワード> --wait
xcrun stapler staple dist/ADIT.app
```

#### 確認 (実機で最初に通すこと)

Windows の 6 点 (上の「2 出来た実行ファイルの確認」) に、別の Mac で確かめる 1 点を足した 7 点です。**通ったものに印を付けて、この文書を更新してください。**

1. `ADIT.app` を二重クリックして画面が出る (日本語が豆腐 □ にならない)
2. 構造を作り、入力を生成できる (`adit-cli gen --sample water_generated` → 生成 → `README.txt` がある)
3. 解析が図まで出る (`adit-cli analyze <ディレクトリ>`)
4. Draw が開き、描いた分子が SMILES になる (RDKit が同梱されているか)
5. ウェブ版が開く (`adit-cli web` → ブラウザで 127.0.0.1:8765)
6. 3D 表示が回る (画面の「構造」タブと、ウェブ版の「構造 (3D)」)
7. 別の Mac (作った機械ではない Mac) に写して、上の 1〜6 が通る

**7 が最も大事です。**作った機械では動いて、他の機械では足りないものがあって動かない、が起きます。

#### Intel の Mac

**配っているのは Apple Silicon (arm64) 版だけです。**GitHub が Intel の実行環境 (`macos-13`) の提供を
終えたため、CI では作れません。Intel の Mac では pip で入れてください。
1 つにまとめた universal2 は、PySide6 と RDKit が両対応の wheel を配っていないと作れないので、していません。

## インストールと起動 (Linux / macOS / 慣れている人向け)

Python 3.11 以上が必要です。仮想環境 (conda か venv) を作り、pip でインストールします。Windows で初めての人は上の節から進めてください。
依存パッケージは ASE、pydantic、Jinja2、matplotlib、SciPy、tomli-w と、デスクトップ版には PySide6、pyqtdarktheme-fork、
ワークスペースのターミナル用の pyte と ptyprocess (Windows は pywinpty) です。SMILES と Draw (分子を描く機能) を使うには RDKit も必要です。

```bash
python3 -m venv adit-env                # conda を使うなら上の節の 3 のとおり
source adit-env/bin/activate            # Windows (PowerShell) では  adit-env\Scripts\Activate.ps1
git clone https://github.com/adit-chem-project/adit-chem.git
cd adit-chem && pip install -e ".[gui,smiles,analysis]"
```

配布名は `adit-chem` です (PyPI の `adit` は別のパッケージです)。インポート名とコマンド名は `adit` です。

| 使い方 | コマンド | 備考 |
|---|---|---|
| デスクトップ版 | `adit` | PySide6 が必要です |
| ウェブ版 (ブラウザ) | `adit-web --open` | PySide6 は不要です。`http://127.0.0.1:8765/` が開きます |
| コマンドライン (生成) | `adit-gen spec.json <出力ディレクトリ>` | GUI で保存した計算設定から生成します |
| コマンドライン (変換) | `adit-convert structure ...` / `adit-convert calculation ...` | 構造形式を変換するか、共通条件を保って別の計算コードの入力を生成します |
| コマンドライン (解析) | `adit-analyze <計算結果のディレクトリ> --rdf --msd --dos` | 図と要約を `analysis/` に書きます |

### 構造形式と計算コードの変換

構造ファイルは、ASE が対応する形式の間で変換できます。入力と出力の形式は通常、ファイル名から判定されます。
デスクトップ版では「ファイル」タブの「変換」、ウェブ版では上部の「変換」から同じ機能を使えます。判定できないファイル名では、ASE の形式名を入力形式・出力形式の欄に指定できます。
どちらの画面でも、準備画面でプレビューした現在の構造を extended XYZ として保存できます。セル、周期境界、固定原子・軸固定、設定済みの初速度も保持されます。CLI で `spec.json` や生成ディレクトリから構造を変換するときも同じ情報を読みます。

```bash
adit-convert structure data.lammps POSCAR
adit-convert structure POSCAR structure.data
adit-convert structure POSCAR structure.xyz
adit-convert structure trajectory.extxyz final.cif --frame -1
adit-convert combine molecule_a.xyz molecule_b.xyz combined.xyz
adit-convert openbabel source.sdf target.mol2 --input-format sdf --output-format mol2
adit-convert dock6 prepared/dock.in ready/
```

GOAT、ORCA DOCKER、DCDFTBMD 2.0、DOCK6 の入力整理、Open Babel の変換入口は、[対応ソフトウェア](SOFTWARE.md) に手順と制限を書いています。これらはコマンドラインだけで、デスクトップ版とウェブ版の欄はまだありません。
