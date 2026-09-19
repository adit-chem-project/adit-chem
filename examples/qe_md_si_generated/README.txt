ADIT 0.0.1 が生成した Quantum ESPRESSO (pw.x) の計算ディレクトリです (2026-09-10 18:01 (UTC+09:00))

計算の種類: 分子動力学 (MD) / 元素: Si / 原子数: 2
計算コード: Quantum ESPRESSO (pw.x)
実行先: プロファイル local (この PC で直接実行する)
  (プロファイル = ADIT の環境設定に書いた「どこで・どう実行するか」の組。ADIT の画面の「プロファイル」で選びます)

== このディレクトリのファイル ==
  submit.sh     ジョブスクリプト (計算を実行する手順を書いたシェルスクリプト。下の「実行する」で使います)
  spec.json     ADIT で決めた設定一式。ADIT の「ファイル」→「計算設定 (spec.json) を開く…」で読み込むと、同じ設定を画面に戻せます
  analyze.py    解析スクリプト (下の「結果の見方」を参照)
  pw.in         pw.x の入力 (&CONTROL &SYSTEM などの設定と、元素・座標・k 点 = 周期系で電子の状態を計算する波数空間の点)
  pseudo/       この計算に必要な擬ポテンシャル (内殻電子の効果をまとめた、元素ごとのファイル。UPF 形式): Si.pbe-n-rrkjus_psl.1.0.0.UPF

== この PC で実行する ==
  以下はターミナル (コマンドを 1 行ずつ打ち込んで PC を操作する画面。Windows なら WSL の Ubuntu) で行います。
  1. このディレクトリへ移動します (cd = 作業する場所を変えるコマンド)
     <...> を、いまこのディレクトリを置いている場所のパスに置き換えてください。別の PC へ写した場合は、写し先のパスを使います。
       cd '<この計算ディレクトリのパス>'
  2. pw.x が PATH (コマンドを探す場所の一覧) にあるか確かめます
       command -v pw.x
     場所 (例: /home/.../bin/pw.x) が 1 行出れば準備できています。
     何も出なければ、計算コードを入れた conda 環境 (コードごとに分けたインストール先) を有効にしてから、
     もう一度確かめます。例: conda activate <環境名>  (入れていなければ: conda install -c conda-forge qe)
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
       espresso = ["<module 名>"]                   # 計算コードを使えるようにする module の名前 (クラスタで module avail と打つと一覧が出ます)
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
       cd <クラスタでの作業ディレクトリ>/qe_md_si_generated
  3. ジョブを預けます (投入)。ジョブ番号が表示されれば受け付けられています
       qsub submit.sh     (PBS のとき)   /   sbatch submit.sh   (Slurm のとき)
  4. 状態を確かめます ($USER は自分のユーザー名に置き換わります)。一覧から消えたら終わっています
       qstat -u $USER     (PBS のとき)   /   squeue -u $USER    (Slurm のとき)
  5. 終わったら結果を手元に取り戻します (手元の PC のターミナルで)
       scp -r <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/qe_md_si_generated <手元の置き場所>/
     または
       rsync -av <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/qe_md_si_generated <手元の置き場所>/

== 結果の見方 ==
  output.log    pw.x が画面に出す文字を保存したもの ('!    total energy' の行が全エネルギー。単位は Ry = リュードベリ)
                各ステップの座標 (ATOMIC_POSITIONS) と温度も output.log に出ます
  tmp/          波動関数と電荷密度 (続きの計算に使う途中のファイル。大きくなることがあります)
  analyze.py    ADIT の解析タブと同じ処理で、図と要約を analysis/ に書きます。ADIT が入った Python で実行します
       python analyze.py --rdf --msd   (動径分布関数と拡散係数も求めます)

== 失敗したときに文献が挙げる対処 ==
  ADIT が出力を走査して見つけた印 (解析の要約の「失敗の原因の候補」) ごとに、その計算コードの文書が挙げる対処を写したものです。
  ADIT の推奨ではありません。数値や手順は出典のとおりで、系に合うかはご自身で判断してください。
  - 「convergence NOT achieved after N iterations」: aiida-quantumespresso の PwBaseWorkChain は mixing_beta を 0.8 倍にして、前の計算から再開します (数値は aiida-quantumespresso の実装値: delta_factor_mixing_beta = 0.8)。relax の途中で起きた場合は、出力の構造に更新して最初から計算し直します。出典: https://github.com/aiidateam/aiida-quantumespresso/blob/main/src/aiida_quantumespresso/workflows/pw/base.py
  - 「Maximum CPU time exceeded」: 出力の構造に更新し、CONTROL の restart_mode = 'restart' で続きを計算します。同ワークチェーンは max_seconds を制限時間の 0.95 倍に設定しています (数値は aiida-quantumespresso の実装値: delta_factor_max_seconds = 0.95)。出典: https://github.com/aiidateam/aiida-quantumespresso/blob/main/src/aiida_quantumespresso/workflows/pw/base.py
  - 対角化の失敗 (「S matrix not positive definite」「problems computing cholesky」「zhegvd failed」「too many bands are not converged」「eigenvectors failed to converge」「Error in routine broyden」): ELECTRONS の diagonalization を david / paro / cg の別のものに順に切り替え、全部試したら止めます。出典: https://github.com/aiidateam/aiida-quantumespresso/blob/main/src/aiida_quantumespresso/workflows/pw/base.py
  - 「history already reset at previous step」(BFGS の履歴の失敗) / 「bfgs failed」: relax では ion_dynamics = 'damp' にする (対称性の低い系では、まず構造を標準偏差 0.01 の乱数で揺らす: 数値は aiida-quantumespresso の実装値 rattle_stdev = 0.01)。vc-relax では trust_radius_min を 0.1 倍にし (数値は aiida-quantumespresso の実装値 delta_factor_trust_radius_min = 0.1)、それでも駄目なら ion_dynamics = 'damp' と cell_dynamics = 'damp-w'。出典: https://github.com/aiidateam/aiida-quantumespresso/blob/main/src/aiida_quantumespresso/workflows/pw/base.py
  - 構造最適化が nstep 内に収束しない: 出力の構造に更新し、電荷密度を読んで再開します。出典: https://github.com/aiidateam/aiida-quantumespresso/blob/main/src/aiida_quantumespresso/workflows/pw/base.py
  - 「Not enough space allocated for radial FFT」(vc-relax で体積が大きく縮んだ): CELL の cell_factor を 2 倍にして (既定 2 → 4) 最初から計算し直します。出典: https://github.com/aiidateam/aiida-quantumespresso/blob/main/src/aiida_quantumespresso/workflows/pw/base.py
  - relax の途中で SCF の精度が上がらなくなった: IONS の upscale を 0.3 倍にして (最小 1) 再開します (数値は aiida-quantumespresso の実装値: delta_factor_upscale = 0.3)。出典: https://github.com/aiidateam/aiida-quantumespresso/blob/main/src/aiida_quantumespresso/workflows/pw/base.py
  - 最高のバンドが占有されている (バンド数の不足): nbnd を「5 % か 4 本の大きいほう」だけ増やします (数値は aiida-quantumespresso の実装値: delta_factor_nbnd = 0.05、delta_minimum_nbnd = 4)。出典: https://github.com/aiidateam/aiida-quantumespresso/blob/main/src/aiida_quantumespresso/workflows/pw/base.py
  - 「DUE TO TIME LIMIT」(Slurm) / 「=>> PBS: job killed: walltime … exceeded limit …」(PBS): 制限時間で止められた印。出典: Slurm slurmstepd https://github.com/SchedMD/slurm/blob/master/src/slurmd/slurmstepd/req.c、OpenPBS pbs_mom https://github.com/openpbs/openpbs/blob/master/src/resmom/mom_main.c
  - 「oom_kill event … OOM Killed」(Slurm) / 「=>> PBS: job killed: mem …」(PBS) / 「Out of memory: Killed process」(Linux) / 行頭近くの「Killed」: メモリ不足で強制終了された印。出典: Slurm https://github.com/SchedMD/slurm/blob/master/src/plugins/task/cgroup/task_cgroup_memory.c、Linux https://github.com/torvalds/linux/blob/master/mm/oom_kill.c
  - 「MPI_ABORT was invoked on rank …」(Open MPI) / 「application called MPI_Abort」「BAD TERMINATION OF ONE OF YOUR APPLICATION PROCESSES」(MPICH): MPI のプロセスが異常終了した印。原因はその前の計算コードの行にあります。出典: Open MPI help-mpi-api.txt (配布物に同梱)、MPICH https://github.com/pmodels/mpich/blob/main/src/pm/hydra/mpiexec/pmiserv_cb.c

== 作成時の記録 ==
  ADIT 0.1.0a1 / Python 3.14.7 / ASE 3.29.0 / 生成 2026-09-19T05:32:36+00:00 (UTC)
  パラメータなどのファイルの SHA-256 (中身から計算する照合用のハッシュ。同じ値なら同じファイル。spec.json の provenance にも同じもの):
    669fb75395a9d26973b0ea1ce8223bbcb30d3396c5d48bf5e794d1243c52375a  pseudo/Si.pbe-n-rrkjus_psl.1.0.0.UPF
  計算コードのバージョン: 走り終えると submit.sh が出力のバージョンの行を code_version.txt に写します
