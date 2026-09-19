ADIT 0.0.1 が生成した GROMACS の計算ディレクトリです (2026-09-12 02:13 (UTC+09:00))

計算の種類: 分子動力学 (MD) / 元素: O H / 原子数: 648
計算コード: GROMACS
実行先: プロファイル local (この PC で直接実行する)
  (プロファイル = ADIT の環境設定に書いた「どこで・どう実行するか」の組。ADIT の画面の「プロファイル」で選びます)

== このディレクトリのファイル ==
  submit.sh     ジョブスクリプト (計算を実行する手順を書いたシェルスクリプト。下の「実行する」で使います)
  spec.json     ADIT で決めた設定一式。ADIT の「ファイル」→「計算設定 (spec.json) を開く…」で読み込むと、同じ設定を画面に戻せます
  analyze.py    解析スクリプト (下の「結果の見方」を参照)
  grompp.mdp    計算の設定 (mdp。gmx grompp が読みます)
  topol.top     トポロジー (topol.top を写したもの。分子の種類と数、力場)
  conf.gro     構造 (adit.gro を写したもの)
                GROMACS に付いている力場から読むもの: oplsaa.ff/forcefield.itp, oplsaa.ff/spce.itp
                submit.sh が gmx grompp (設定とトポロジーをまとめて adit.tpr にする) → gmx mdrun (計算) の順に実行します

== 実行する前に用意すること ==
  この入力は時間刻み dt = 2 fs、constraints = none を使います。時間刻みと拘束の組み合わせは、使用する力場の文書に照らして利用者が決めます。
  初速を作る乱数種 gen-seed = -1 は grompp.mdp と spec.json の両方に記録されます。

== この PC で実行する ==
  以下はターミナル (コマンドを 1 行ずつ打ち込んで PC を操作する画面。Windows なら WSL の Ubuntu) で行います。
  1. このディレクトリへ移動します (cd = 作業する場所を変えるコマンド)
     <...> を、いまこのディレクトリを置いている場所のパスに置き換えてください。別の PC へ写した場合は、写し先のパスを使います。
       cd '<この計算ディレクトリのパス>'
  2. gmx が PATH (コマンドを探す場所の一覧) にあるか確かめます
       command -v gmx
     場所 (例: /home/.../bin/gmx) が 1 行出れば準備できています。
     何も出なければ、計算ソフトを入れた conda 環境 (ソフトごとに分けたインストール先) を有効にしてから、
     もう一度確かめます。例: conda activate <環境名>  (入れていなければ: conda install -c conda-forge gromacs)
  3. 実行します
       bash submit.sh
     画面に何も出なくても動いています。記録は output.log に書かれ、終わると次のコマンドを打てる状態に戻ります。
     途中経過は別のターミナルで  tail -f output.log  (Ctrl+C で表示だけ止まり、計算は続きます)

== 研究室のクラスタで実行したいとき ==
  クラスタ = 研究室や計算センターが共同で使う計算機の集まり。計算はジョブスケジューラ (PBS や Slurm。計算の順番待ちを
  管理するソフト) に預けて実行します。預けた 1 件の計算を「ジョブ」と呼びます。
  いまの submit.sh はこの PC 用です。クラスタで使うには、先に次の 2 つを行います。
  a. ADIT の環境設定ファイルに、kind = "pbs" か "slurm" のプロファイルを足します。環境設定ファイルの場所:
       /home/<user>/.config/adit/cluster.toml
     いちばん短い書き方は次のとおりです。ファイルの最後に書き足し、<...> を自分の値に置き換えます (# から行末までは説明で、
     消してもかまいません)。キュー名や module 名 (クラスタで計算ソフトを使えるようにするための名前) はクラスタごとに違うので、
     管理者か研究室の先輩に確かめてください。
       [profiles.remote]
       kind = "pbs"                                # Slurm のクラスタなら "slurm"
       header_extra = ["#PBS -q <キュー名>"]        # Slurm なら ["#SBATCH --partition=<パーティション名>"]
       [profiles.remote.code_modules]
       gromacs = ["<module 名>"]                   # 計算ソフトを使えるようにする module の名前 (クラスタで module avail と打つと一覧が出ます)
     ほかの項目 (投入コマンド、環境変数、実行コマンドなど) は ADIT の README.md の「クラスタで実行する場合」にあります。
  b. ADIT でプロファイルをそれに切り替えて生成し直します。submit.sh が qsub (PBS) / sbatch (Slurm) で預ける
     ジョブスクリプトになり、この README.txt にも、そのクラスタでの投入と確認のコマンドが入ります。
  そのあとの流れは次のとおりです (<...> は自分の値に置き換えます)。
  1. このディレクトリごとクラスタへ写します (手元の PC のターミナルで。scp / rsync = ネットワーク越しにファイルを写すコマンド)
     <この計算ディレクトリのパス> は、手元でいま置いている場所に置き換えてください。
       scp -r '<この計算ディレクトリのパス>' <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/
     または
       rsync -av '<この計算ディレクトリのパス>' <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/
  2. クラスタにログインして、写したディレクトリへ移動します (ssh = 別の計算機にログインするコマンド)
       ssh <ユーザー名>@<クラスタのホスト名>
       cd <クラスタでの作業ディレクトリ>/gromacs_spce_nvt_generated
  3. ジョブを預けます (投入)。ジョブ番号が表示されれば受け付けられています
       qsub submit.sh     (PBS のとき)   /   sbatch submit.sh   (Slurm のとき)
  4. 状態を確かめます ($USER は自分のユーザー名に置き換わります)。一覧から消えたら終わっています
       qstat -u $USER     (PBS のとき)   /   squeue -u $USER    (Slurm のとき)
  5. 終わったら結果を手元に取り戻します (手元の PC のターミナルで)
       scp -r <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/gromacs_spce_nvt_generated <手元の置き場所>/
     または
       rsync -av <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/gromacs_spce_nvt_generated <手元の置き場所>/

== 結果の見方 ==
  grompp.log    gmx grompp の出力 (トポロジーの誤りや注意はここに出ます)
  output.log / adit.log   gmx mdrun の出力 (エネルギーの表)
  adit.edr     エネルギーの記録。gmx energy -f adit.edr で温度・圧力・密度などを取り出せます
  adit.gro     最後の構造
  adit.xtc     MD の軌跡 (dump の間隔ごと。VMD や MDAnalysis で開けます)
  adit.cpt     チェックポイント (続きから実行するための状態)。次の段階 (NVT → NPT → 本計算) を ADIT で作るとき「前の段階の .cpt」に指定します
  analyze.py    ADIT の解析タブと同じ処理で、図と要約を analysis/ に書きます。ADIT が入った Python で実行します
       python analyze.py --rdf --msd   (動径分布関数と拡散係数も求めます)

== 作成時の記録 ==
  ADIT 0.1.0a1 / Python 3.14.7 / ASE 3.29.0 / 生成 2026-09-19T01:29:36+00:00 (UTC)
  パラメータなどのファイルの SHA-256 (中身から計算する照合用のハッシュ。同じ値なら同じファイル。spec.json の provenance にも同じもの):
    5f636d8165bd99daf3f3bb0dbeddb36ede461c8547f4a3ecac921c04f0166c6f  conf.gro
    f86096c40f9cb9c75da0b5ce22f7920532a832884b8b77e4da58f32517ee3f77  topol.top
  計算コードのバージョン: 走り終えると submit.sh が出力のバージョンの行を code_version.txt に写します
