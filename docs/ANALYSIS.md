# 解析の詳細

解析タブとコマンド (`adit-analyze`) で選べる項目の中身です。
環境設定は [環境設定](SETTINGS.md) にあります。

## 解析

各出力ディレクトリの `analyze.py` を計算後に実行すると (`python analyze.py --rdf --msd --dos`)、`analysis/` に図 (PNG) と要約 (summary.txt / summary.json) が書かれます。
GUI の「解析」タブ、ウェブ版の「解析」ページ、`adit-analyze` も同じ処理です。数値の良し悪しの判断はしません。

| 解析 | 読むファイル | 出力 |
|---|---|---|
| エネルギー・温度の推移 | DFTB+ md.out / output.log、xtb output.log と xtb.trj、VASP vasprun.xml と OSZICAR、pw.x output.log、ORCA output.log | 折れ線グラフ |
| 動径分布関数 (RDF) | 軌跡 (geo_end.xyz、xtb.trj、vasprun.xml / XDATCAR、pw.x output.log、ORCA trajectory.xyz) | 元素の組ごとの g(r) |
| 平均二乗変位 (MSD) と拡散係数 | 同上 (周期系は境界をまたぐ移動を補正) | MSD の図と、既定では最大ずれ時間の 10〜50 % を直線フィットして求めた D [cm²/s] |
| 状態密度 (DOS) | DFTB+ band.out と detailed.out のフェルミ準位、VASP DOSCAR | ガウス関数で広げた DOS |
| 振動数とスペクトル | DFTB+ hessian.out (質量重み付きヘシアンの対角化)、xtb vibspectrum、VASP OUTCAR、ORCA output.log | 振動数の一覧とスペクトル (IR 強度があれば重み付き) |
| 結合長 | 最終構造 | 共有結合半径の和の 1.2 倍以内にある原子対 |
| バンド図 | `bands/kpath.json` と DFTB+ band.out、pw.x output.log、VASP EIGENVAL | 高対称点のラベル付きバンド図と、フェルミ準位をまたぐ最小の間隔 |
| 原子の電荷 | DFTB+ detailed.out、xtb xtbout.json、CP2K output.log (Mulliken・Hirshfeld) | 原子ごとの表。電荷の定義と出典 (ファイル:行) を添える |
| 熱化学 (コードの値) | xtb・ORCA の output.log | ZPE・H・G などをそのまま (行番号つき) |
| 熱化学 (ASE で計算) | 振動数と最終構造 | 理想気体・調和・準調和・準 RRHO。温度などは利用者が入れる (既定値は無く、足りなければ計算しない理由を出す)。thermo.csv |
| HOMO-LUMO ギャップ・双極子 | xtb xtbout.json など | 値と出典 |
| 時系列の統計 | MD の温度・エネルギー・密度・圧力 | ブロック平均による平均値の誤差 (blocking.png)、積分自己相関時間 |
| 配位数 n(r) | RDF と同じ軌跡 | coordination.png、rdf.json の n と n_reverse |
| z 方向の密度分布 | 周期系の軌跡 | 元素ごとの数密度と質量密度 (zdensity.png / zdensity.json) |
| 圧力 | LAMMPS log.lammps、GROMACS のエネルギーの表 | pressure.png |
| 元素ごとの拡散係数 | MSD と同じ軌跡 | 元素ごとの D [cm²/s] (当てはめ範囲は指定できる) |
| NEB | VASP の像のディレクトリ (00, 01, …)、QE neb.x の .dat / .int、CP2K BAND の output.log | エネルギーの曲線 (neb.png) と障壁。CP2K は像の値のみ |
| PDOS | VASP DOSCAR、QE projwfc.x の *.pdos_atm#… | pdos.png / pdos.csv |
| UV-Vis | ORCA の TD-DFT の吸収の表 | 遷移の表 (uvvis_transitions.csv) と、広げたスペクトル (uvvis.png) |
| 空間群 | 最終構造 (spglib が入っているとき) | 許容誤差ごとの空間群 |
| フォノン分散・DOS | phonopy の band.yaml、total_dos.dat | phonon_bands.png、phonon_dos.png |
| 組にして比べる表 | 複数の計算のディレクトリ | ΣνE、組成の釣り合い、条件が違う項目 (compare_*.csv、compare_energy.png) |

`adit-analyze` の主なオプション (`adit-analyze --help-all` に全部):

| オプション | 意味 |
|---|---|
| `--stride N` | 軌跡を N フレームおきに使う (大きな軌跡の間引き。RDF・MSD・z 密度・書き出しに効く) |
| `--msd-fit T0 T1` | 拡散係数を当てはめる時間の範囲 [fs] (既定は最大の遅れ時間の 10〜50 %) |
| `--zdens [Å]` | z 方向の密度分布 (数を続けると区間の幅。既定 0.2 Å。周期系だけ) |
| `--export` | 軌跡を `analysis/export/` に書き出す (extxyz・xyz・pdb、VMD の view.vmd と vmd_load.tcl、OVITO の ovito_pipeline.py、TRAVIS の答えファイル travis_*.in、export_README.txt。下の「書き出したファイル」) |
| `--unwrap-molecules` | 書き出す前に分子を周期境界でつなぎ直す (`--export` を含む) |
| `--memory-mb MB` | MSD で座標を持つメモリの上限 (既定 1024)。超える軌跡は読まずに止まり、間引きの間隔を示す |
| `--compare [組]` | 組にして比べる表。`"ads=1:slab_mol,-1:slab,-1:mol"` (係数:ディレクトリ。反応は `;` で区切る)。省くと compare.json を読む |
| `--thermo MODEL` ほか | 熱化学 (ASE)。`--temperature` `--pressure` `--symmetry-number` `--geometry` `--spin` `--imaginary` `--exclude-lowest` `--qh-cutoff` `--msrrho-tau` をモデルに合わせて入れる |
| `--uvvis SHAPE:FWHM` | ORCA の UV-Vis を広げる形と半値全幅 [eV] (例 `gauss:0.3`) |
| `--pdos` | PDOS のファイルが無いときも理由を書く (あれば指定しなくても描く) |
| `--symprec Å[,Å…]` | 空間群の許容誤差 (既定は 1e-5、1e-3、1e-1 を並べる) |

```bash
adit-analyze out/md --msd --stride 10 --msd-fit 1000 5000 --zdens 0.5
adit-analyze out/md --export --unwrap-molecules
adit-analyze out/vib --thermo ideal_gas --temperature 298.15 --pressure 100000 --symmetry-number 2 --geometry nonlinear --spin 0
adit-analyze runs/ --compare "ads=1:slab_mol,-1:slab,-1:mol"
```

画面では、解析タブ (ウェブ版は解析のページ) の「詳しい条件」を開くと、上のオプションと同じ欄があります (間引き、MSD の当てはめ範囲、
z 方向の密度分布、時系列の統計、熱化学の各欄、UV-Vis の広げ方、空間群の許容誤差)。熱化学の欄には既定値を入れていません。
結果は要約の下に表 (原子の電荷、熱化学、電子状態、時系列の統計、軌跡、拡散係数、空間群、UV-Vis の遷移) と図で出ます。表は見出しを押すと畳めます。
長い表は先頭 20 行だけを出し、全体のファイルの場所を添えます。

- 軌跡が大きすぎて止まったときは、理由の文と「間引きを N にする」ボタンが出ます。押すと間引きの欄に N が入るので、もう一度「解析を実行」を押します
- MSD では周期境界を越えた移動を最小像から復元します。間引き後の隣接フレーム間で移動が最短セル幅の 40 % 以上になった場合は、移動方向を一意に復元できなくなる半セル幅へ近いため、要約と `summary.json` に注意を残します。
- MSD の D は指定した時間範囲の直線フィットから求めます。ブロックごとの D を使う参考誤差は、すべてのブロックで同じ時間範囲を使える場合だけ表示します。短い軌跡では D が出てもブロック誤差は出ないことがあります。参考誤差はブロック D の平均の標準誤差であり、全軌跡から求めた D の厳密な誤差ではありません。
- 「TRAVIS・VMD・OVITO 用に書き出す」で `analysis/export/` に書き出し、export_README.txt の中身を表示します。
  デスクトップ版は「書き出したフォルダを開く」でファイルマネージャを開きます (ウェブ版は場所を表示するだけ)
- 「組にして比べる…」で、比べる計算のディレクトリと係数を行で入れます (名前が空の行は上の行と同じ反応)。
  デスクトップ版は小さな画面、ウェブ版は別のページ (`/compare`)。基準のディレクトリに compare.json があれば読み込めます

### 書き出したファイル (TRAVIS・OVITO・VMD 用) / Exported files for TRAVIS, OVITO and VMD

`--export` は `analysis/export/` に次を書きます。**どれも ADIT では実行していません** (TRAVIS・OVITO・VMD は手元に無い前提)。
雛形の中に数値は入れず、`<...>` の穴か、出典を添えた「例」にしています。

| ファイル | 中身 | 確かめ方 |
|---|---|---|
| `trajectory.extxyz` / `.xyz` / `.pdb` | 軌跡 (extxyz にはセル、pdb には各フレームの CRYST1) | ASE で読み戻す検査 |
| `view.vmd` | VMD の読み込みとセル (最小) | 検査で文字列を照合 |
| `vmd_load.tcl` | VMD の読み込みに加えて代表的な表示: 全原子 CPK、`--select` の原子を VDW、白背景、平行投影。DynamicBonds と画像の描画は例 (コメントアウト) | コマンド名と引数は VMD User's Guide (mol / color / display / axes / render) と PBCTools の説明で照合 |
| `ovito_pipeline.py` | 動径分布関数と配位数 (実行部) のあとに、**コメントアウトした雛形**: 構造同定 (CNA / PTM / Ackland-Jones)、Wigner-Seitz 欠陥解析 (参照構造は `<...>`)、変位ベクトルと原子ひずみ、クラスタ解析、空間ビニング + 時間平均、式による選択 | クラス名と引数名は OVITO 3.16 の Python リファレンス (`ovito.modifiers`) で照合 |
| `travis_rdf.in` / `travis_cdf.in` / `travis_msd.in` / `travis_hbond.in` / `travis_acf.in` | TRAVIS の `-i` に渡す答えファイル。セルの大きさ [pm] と関数の選択まで | 書式 (1 行 1 答、空行 = 既定値、`!` = コメント) は TRAVIS の Quick Start Guide とソースで確認。**質問の並びはソースを読んだもので、実行では未確認** |
| `export_README.txt` | 上の説明と、TRAVIS の残りの質問の目安 (時間の刻み、分子と原子の選択、フレームの範囲) | — |

TRAVIS の答えファイルは、**直方体で固定のセル**のときだけ書きます (直方体でないセルは advanced mode の手順が未確認、非周期は TRAVIS がセルを聞くので手で答える)。
ファイルが尽きると TRAVIS はキーボード入力に切り替わるので、残りは対話で答えられます。公式の勧めどおり、一度 `-i` なしで走らせて TRAVIS が書く `input.txt` を次回に使うのが確実です。
速度自己相関は TRAVIS では関数名 `acf` です (`vacf` ではありません)。

`--select` の式 (`element O and z < 10` など) は、`vmd_load.tcl` では VMD の選択式 (`name O and z < 10`。番号は 0 から)、
`ovito_pipeline.py` では OVITO の式 (`ParticleType == "O" && Position.Z < 10`) に直して入れます。`within` は OVITO の式には直せないので、その旨を書きます。

`--export` writes the files above into `analysis/export/`. None of them has been run by ADIT. Templates contain no numerical values: they use
`<...>` holes or examples with their source. The VMD commands were checked against the VMD User's Guide and the PBCTools page, the OVITO class
and argument names against the OVITO 3.16 Python reference, and the TRAVIS answer-file format (one answer per line, empty line = default,
`!` = comment) against the TRAVIS Quick Start Guide and source. The order of the TRAVIS questions was read from the source and is not verified
by running; answer files are written only for orthorhombic fixed cells, and TRAVIS switches to keyboard input when the file runs out.
A `--select` expression is rewritten as a VMD selection in `vmd_load.tcl` and as an OVITO expression in `ovito_pipeline.py` (`within` cannot be rewritten for OVITO).

## 収束の確認と格子定数 (1 つの条件を変えて一括生成)

平面波 DFT (VASP、Quantum ESPRESSO) の結果は、化学的な条件とは別に、数値計算の細かさの設定で変わります。
代表は **カットオフエネルギー** (電子の波をどこまで細かく表すか。QE の `ecutwfc`、VASP の `ENCUT`) と **k 点** (結晶の中の電子の状態を何点で代表させるか) です。
小さすぎると答えがずれ、大きすぎると計算時間が膨らむので、値を上げていき、エネルギーがほとんど変わらなくなる所を探します。これが収束の確認です。
格子の大きさを少しずつ変えてエネルギーを並べると、平衡の格子定数と体積弾性率 (どれだけ潰れにくいか) も求まります。

ADIT は、1 つの条件だけを変えた入力を値ごとのディレクトリにまとめて作り、実行したあとの結果を表と図にします。

**画面から**

1. いつもどおり構造と計算の条件を決めます。
2. 「実行」メニューの「1 つの条件を変えて一括生成…」を開き、変える項目 (カットオフ / k 点の分割数 / k 点の密度 / 格子の大きさ / その他) と値 (カンマ区切り) を入れて「生成」を押します。
3. できたディレクトリをそれぞれ実行します (この PC なら各ディレクトリで `bash submit.sh`、クラスタならジョブとして。大きな系はクラスタで)。
4. すべて終わったら、解析タブで「解析を実行」を押します (ディレクトリの欄には保存先が入っています)。

**コマンド行から**

```bash
adit-gen spec.json out/ --scan method.ecutwfc=30,40,50,60     # QE のカットオフ [Ry]
adit-gen spec.json out/ --scan kpoints.mesh=4x4x4,6x6x6,8x8x8  # k 点の分割数
adit-gen spec.json out/ --scan scale=0.97,0.98,0.99,1.00,1.01,1.02,1.03   # 格子の大きさ (周期系)

# それぞれを実行したあと
adit-analyze out/ --scan
```

**表と図の読み方** (`out/scan_energies.csv`、`out/scan_energy.png`、格子の大きさなら `out/eos.png` と `out/eos.json`)

| 列 | 意味 |
|---|---|
| 最後の値との差 [meV/原子] | いちばん最後に書いた値 (ふつう最も細かい条件) のエネルギーとの差を、原子の数で割ったもの。格子の大きさを変えたときだけは「最小値との差」で、エネルギーがいちばん低い倍率を 0 にします |
| 力の最大値 [eV/Å]、圧力 [GPa] | 最終構造の値 (DFTB+、QE、VASP の出力から読めるときだけ) |
| 状態方程式の当てはめ | 格子の大きさを 5 点以上変えたときだけ。Birch–Murnaghan の式で、平衡の体積・体積弾性率・エネルギーが最小になる倍率とそのときの格子の長さを出します |

どの値で十分とみなすかは、計算の目的によって違うので、ADIT は判断しません。
よく使われる目安は 1 原子あたり 1 meV 程度ですが、反応エネルギーのように差を見る計算ではもっと緩くてよいことが多く、力や応力 (格子定数、振動) を使う計算ではもっと厳しくする必要があります。

**注意すること**

- ふつうはカットオフを先に決め、次に k 点を決めます。
- QE の `ecutrho` (電荷密度のカットオフ) を 0 (指定しない) にすると、pw.x は `ecutwfc` の 4 倍を使います。ウルトラソフトや PAW の擬ポテンシャルでは、もっと大きな値を求められることがあります。
- VASP で体積を変える計算 (格子の大きさの一括生成や格子の緩和) では、カットオフが低いと見かけの応力 (Pulay 応力) が出て、格子定数が小さめに出ます。POTCAR の既定より 3 割ほど高い `ENCUT` がよく勧められます。
- 金属は k 点の収束が遅いので、占有の広げ方 (smearing) とその幅も合わせて決めます。
- 値の数だけ計算が増えます。

## ウェブ版

```bash
adit-web --open        # http://127.0.0.1:8765/
```

デスクトップ版と同じ生成と解析を、ブラウザのフォームから行います。「プレビュー」で生成ファイルの内容を確認し、「生成」でサーバー (adit-web を動かしている PC) 上の出力ディレクトリに書き出します。
「解析」ページでは計算結果のディレクトリを解析し、図と要約をそのページに表示します。spec.json の読み込みとダウンロードもできます。

依存は Python の標準ライブラリと Jinja2 だけです。既定では 127.0.0.1 でのみ待ち受けます。
画面の状態 (入力中のフォーム、プレビュー、解析の図) は adit-web 1 つにつき 1 組なので、同時に使うのは 1 つのブラウザからにしてください。

`--host 0.0.0.0` のように他の PC から届くアドレスで起動すると、起動のたびに合言葉 (トークン) が作られ、端末に `http://…:8765/?token=…` の形で表示されます。この URL を一度開くと、そのブラウザでは以後合言葉なしで使えます。合言葉は `--token` で自分で決めることもでき (使える文字は英数字と `_` `-`)、`--no-token` で外せます。通信は暗号化されないので、学内などの信頼できるネットワークの中だけで使ってください。

- 合言葉付きの URL はブラウザの履歴に残ります。共有の PC では、使い終わったら履歴を消してください。
- adit-web を起動し直すと合言葉が変わり、前の URL では開けなくなります (`--token` で決めた場合は変わりません)。
- 通信は暗号化されません (http)。合言葉も平文で流れます。
- `--host 0.0.0.0` (すべてのアドレス) で起動すると、自分用の URL と他の PC 用の URL の 2 行が表示されます。`--host ::` なら IPv6 でも待ち受けます。

開発に加わる人向けの手順 (試験の走らせ方、計算コードの追加) は [CONTRIBUTING.md](../CONTRIBUTING.md) にあります。
