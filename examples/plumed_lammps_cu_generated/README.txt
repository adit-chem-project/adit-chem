ADIT 0.0.1 が生成した LAMMPS の計算ディレクトリです (2026-09-12 02:13 (UTC+09:00))

計算の種類: 分子動力学 (MD) / 元素: Cu / 原子数: 108
計算コード: LAMMPS
実行先: プロファイル local (この PC で直接実行する)
  (プロファイル = ADIT の環境設定に書いた「どこで・どう実行するか」の組。ADIT の画面の「プロファイル」で選びます)

== このディレクトリのファイル ==
  submit.sh     ジョブスクリプト (計算を実行する手順を書いたシェルスクリプト。下の「実行する」で使います)
  spec.json     ADIT で決めた設定一式。ADIT の「ファイル」→「計算設定 (spec.json) を開く…」で読み込むと、同じ設定を画面に戻せます
  analyze.py    解析スクリプト (下の「結果の見方」を参照)
  in.lammps     LAMMPS の入力 (単位系 metal (時間 ps、エネルギー eV、長さ Å、圧力 bar))
  data.lammps   構造 (atomic 形式。ASE で書いたもの。斜めのセルは LAMMPS の向きに回して書かれます)
                原子の型番号と元素: 1 = Cu
  Cu_u3.eam   力場・モデルのファイル (pair_coeff などから名前で読みます)
  plumed.dat   PLUMED の入力 (利用者が書いたもの)

== 実行する前に用意すること ==
  PLUMED の入力 plumed.dat は利用者が書いたものです。ADIT は集合変数もバイアスも作らず、内容も確かめていません(読めるかどうかだけ、plumed があるときに確かめています)。
  LAMMPS は PLUMED パッケージを入れてビルドしたものが要ります (fix plumed)。
  PLUMED の既定の単位は nm・kJ/mol・ps で、計算コードの単位系とは別です (入力の UNITS で変えられます)。原子の番号は 1 始まりです。
  書き方と機能の説明は https://www.plumed.org/doc-v2.10/user-doc/html/index.html。

== この PC で実行する ==
  以下はターミナル (コマンドを 1 行ずつ打ち込んで PC を操作する画面。Windows なら WSL の Ubuntu) で行います。
  1. このディレクトリへ移動します (cd = 作業する場所を変えるコマンド)
     <...> を、いまこのディレクトリを置いている場所のパスに置き換えてください。別の PC へ写した場合は、写し先のパスを使います。
       cd '<この計算ディレクトリのパス>'
  2. lmp が PATH (コマンドを探す場所の一覧) にあるか確かめます
       command -v lmp
     場所 (例: /home/.../bin/lmp) が 1 行出れば準備できています。
     何も出なければ、計算ソフトを入れた conda 環境 (ソフトごとに分けたインストール先) を有効にしてから、
     もう一度確かめます。例: conda activate <環境名>  (入れていなければ: conda install -c conda-forge lammps)
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
       lammps = ["<module 名>"]                   # 計算ソフトを使えるようにする module の名前 (クラスタで module avail と打つと一覧が出ます)
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
       cd <クラスタでの作業ディレクトリ>/plumed_lammps_cu_generated
  3. ジョブを預けます (投入)。ジョブ番号が表示されれば受け付けられています
       qsub submit.sh     (PBS のとき)   /   sbatch submit.sh   (Slurm のとき)
  4. 状態を確かめます ($USER は自分のユーザー名に置き換わります)。一覧から消えたら終わっています
       qstat -u $USER     (PBS のとき)   /   squeue -u $USER    (Slurm のとき)
  5. 終わったら結果を手元に取り戻します (手元の PC のターミナルで)
       scp -r <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/plumed_lammps_cu_generated <手元の置き場所>/
     または
       rsync -av <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/plumed_lammps_cu_generated <手元の置き場所>/

== 結果の見方 ==
  output.log / log.lammps   LAMMPS の出力 (thermo の表: ステップ、時刻、温度、ポテンシャルエネルギー pe、運動エネルギー、全エネルギー、圧力、体積)
  traj.lammpstrj  MD の軌跡 (dump custom 形式、id type element x y z。OVITO や VMD でそのまま開けます)。final.data は最後の構造
  plumed.dat が書き出すファイル (COLVAR、HILLS など) の名前と中身は、その入力しだいです。
  plumed.log   PLUMED 自身のログ (fix plumed の outfile)
  ADIT の解析は、PLUMED の COLVAR 形式のファイル (先頭が #! FIELDS) があれば、その列を読んで要約に入れます。
  analyze.py    ADIT の解析タブと同じ処理で、図と要約を analysis/ に書きます。ADIT が入った Python で実行します
       python analyze.py --rdf --msd   (動径分布関数と拡散係数も求めます)

== 作成時の記録 ==
  ADIT 0.1.0a1 / Python 3.14.7 / ASE 3.29.0 / 生成 2026-09-19T01:29:40+00:00 (UTC)
  パラメータなどのファイルの SHA-256 (中身から計算する照合用のハッシュ。同じ値なら同じファイル。spec.json の provenance にも同じもの):
    3436c491a4c75ea8b7141adbc6ee382a118f5fdb47f609c2a660fc1eb772599f  Cu_u3.eam
  計算コードのバージョン: 走り終えると submit.sh が出力のバージョンの行を code_version.txt に写します
