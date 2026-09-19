ADIT 0.0.1 が生成した DFTB+ の計算ディレクトリです (2026-09-10 17:59 (UTC+09:00))

計算の種類: 分子動力学 (MD) / 元素: O H / 原子数: 3
計算コード: DFTB+
実行先: プロファイル local (この PC で直接実行する)
  (プロファイル = ADIT の環境設定に書いた「どこで・どう実行するか」の組。ADIT の画面の「プロファイル」で選びます)

== このディレクトリのファイル ==
  submit.sh     ジョブスクリプト (計算を実行する手順を書いたシェルスクリプト。下の「実行する」で使います)
  spec.json     ADIT で決めた設定一式。ADIT の「ファイル」→「計算設定 (spec.json) を開く…」で読み込むと、同じ設定を画面に戻せます
  analyze.py    解析スクリプト (下の「結果の見方」を参照)
  dftb_in.hsd   DFTB+ の入力 (計算の設定。DFTB+ はこのディレクトリで起動すると自動でこれを読みます)
  geometry.gen  構造 (原子の種類と座標。gen 形式、長さの単位は Å)
  skf/          Slater-Koster ファイル (セット mio-1-1。DFTB+ が使う、元素の組ごとのパラメータ) と LICENSE、README: H-H.skf, H-O.skf, LICENSE, O-H.skf, O-O.skf, README

== この PC で実行する ==
  以下はターミナル (コマンドを 1 行ずつ打ち込んで PC を操作する画面。Windows なら WSL の Ubuntu) で行います。
  1. このディレクトリへ移動します (cd = 作業する場所を変えるコマンド)
     <...> を、いまこのディレクトリを置いている場所のパスに置き換えてください。別の PC へ写した場合は、写し先のパスを使います。
       cd '<この計算ディレクトリのパス>'
  2. dftb+ が PATH (コマンドを探す場所の一覧) にあるか確かめます
       command -v dftb+
     場所 (例: /home/.../bin/dftb+) が 1 行出れば準備できています。
     何も出なければ、計算コードを入れた conda 環境 (コードごとに分けたインストール先) を有効にしてから、
     もう一度確かめます。例: conda activate <環境名>  (入れていなければ: conda install -c conda-forge dftbplus)
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
     消してもかまいません)。キュー名や module 名 (クラスタで計算コードを使えるようにするための名前) はクラスタごとに違うので、
     管理者か研究室の先輩に確かめてください。
       [profiles.remote]
       kind = "pbs"                                # Slurm のクラスタなら "slurm"
       header_extra = ["#PBS -q <キュー名>"]        # Slurm なら ["#SBATCH --partition=<パーティション名>"]
       [profiles.remote.code_modules]
       dftbplus = ["<module 名>"]                   # 計算コードを使えるようにする module の名前 (クラスタで module avail と打つと一覧が出ます)
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
       cd <クラスタでの作業ディレクトリ>/dftb_md_water_generated
  3. ジョブを預けます (投入)。ジョブ番号が表示されれば受け付けられています
       qsub submit.sh     (PBS のとき)   /   sbatch submit.sh   (Slurm のとき)
  4. 状態を確かめます ($USER は自分のユーザー名に置き換わります)。一覧から消えたら終わっています
       qstat -u $USER     (PBS のとき)   /   squeue -u $USER    (Slurm のとき)
  5. 終わったら結果を手元に取り戻します (手元の PC のターミナルで)
       scp -r <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/dftb_md_water_generated <手元の置き場所>/
     または
       rsync -av <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/dftb_md_water_generated <手元の置き場所>/

== 結果の見方 ==
  output.log    実行ログ (SCC = 電荷を自己無撞着に決める反復。その様子と、各ステップのエネルギー。単位は Hartree)
  detailed.out  最後のステップのエネルギーの内訳、Mulliken 電荷 (原子ごとの電荷の目安)、原子にかかる力
  results.tag   全エネルギーなどを、プログラムで読みやすい形で書いたもの
  geo_end.xyz   MD の軌跡 (MDRestartFrequency ステップごとの構造と速度)
  md.out        各ステップのエネルギーと温度
  dftb_pin.hsd  省略した設定を既定値で埋めた入力。実際に使われた設定を確かめられます
  analyze.py    ADIT の解析タブと同じ処理で、図と要約を analysis/ に書きます。ADIT が入った Python で実行します
       python analyze.py --rdf --msd   (動径分布関数と拡散係数も求めます)

== 失敗したときに文献が挙げる対処 ==
  ADIT が出力を走査して見つけた印 (解析の要約の「失敗の原因の候補」) ごとに、その計算コードの文書が挙げる対処を写したものです。
  ADIT の推奨ではありません。数値や手順は出典のとおりで、系に合うかはご自身で判断してください。
  - 「SCC is NOT converged, maximal SCC iterations exceeded」: DFTB+ recipes は、収束の判定が SccTolerance (既定 1e-5 電子)、反復の上限が MaxSccIterations (既定 100) で、上限に達すると「そこまでの電荷で全エネルギーを計算して警告を出して止まる」と説明しています。出典: https://dftbplus-recipes.readthedocs.io/en/latest/basics/firstcalc.html
  - 同上: DFTB+ manual は、静的な構造でない計算 (Driver が空でない) では ConvergentSccOnly で「収束しなくても先へ進む」ようにできると書いています (MaxSCCIterations の項)。出典: https://github.com/dftbplus/dftbplus/blob/main/doc/dftb+/manual/dftbp.tex
  - 「!!! Geometry did NOT converge!」: manual の MaxSteps は「収束する前にこの回数で最適化を止める」上限で、-1 で事実上無制限です。出典: https://github.com/dftbplus/dftbplus/blob/main/doc/dftb+/manual/dftbp.tex
  - 「WARNING! -> Current stacksize not set to unlimited」: DFTB+ 自身の警告で、ulimit -s unlimited を勧めています (クラスタ用の submit.sh には入れてあります)。出典: DFTB+ の出力そのもの
  - 「DUE TO TIME LIMIT」(Slurm) / 「=>> PBS: job killed: walltime … exceeded limit …」(PBS): 制限時間で止められた印。出典: Slurm slurmstepd https://github.com/SchedMD/slurm/blob/master/src/slurmd/slurmstepd/req.c、OpenPBS pbs_mom https://github.com/openpbs/openpbs/blob/master/src/resmom/mom_main.c
  - 「oom_kill event … OOM Killed」(Slurm) / 「=>> PBS: job killed: mem …」(PBS) / 「Out of memory: Killed process」(Linux) / 行頭近くの「Killed」: メモリ不足で強制終了された印。出典: Slurm https://github.com/SchedMD/slurm/blob/master/src/plugins/task/cgroup/task_cgroup_memory.c、Linux https://github.com/torvalds/linux/blob/master/mm/oom_kill.c
  - 「MPI_ABORT was invoked on rank …」(Open MPI) / 「application called MPI_Abort」「BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES」(MPICH): MPI のプロセスが異常終了した印。原因はその前の計算コードの行にあります。出典: Open MPI help-mpi-api.txt (配布物に同梱)、MPICH https://github.com/pmodels/mpich/blob/main/src/pm/hydra/mpiexec/pmiserv_cb.c

== 作成時の記録 ==
  ADIT 0.1.0a1 / Python 3.14.7 / ASE 3.29.0 / 生成 2026-09-19T05:32:34+00:00 (UTC)
  パラメータなどのファイルの SHA-256 (中身から計算する照合用のハッシュ。同じ値なら同じファイル。spec.json の provenance にも同じもの):
    d43946902e4d120105025683b54c8696856cbd2b459cfd16df166164774359c8  skf/H-H.skf
    e357795fc4da4e6317fa29a1ffe57c13d17b96d8929fc840648cb2430acd1a6d  skf/H-O.skf
    9348ddfd44da5a127c59141981954746a860ec8e03e0412cf3af7134af0f97e2  skf/LICENSE
    0bc4d6d3c013069ac5964cb58d52e50ea05eb7c5e2097fa76dd8255b535727bd  skf/O-H.skf
    423de839173a667cc1fd75bd111e6dfa46435b874544a5ce4441b84a31b80b43  skf/O-O.skf
    92f2a22a702ed3e8548e376a28c47bbba9fe556df8b735d42d4cbb8f67b63e15  skf/README
  計算コードのバージョン: 走り終えると submit.sh が出力のバージョンの行を code_version.txt に写します
