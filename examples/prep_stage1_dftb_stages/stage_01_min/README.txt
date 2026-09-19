ADIT 0.0.1 が生成した DFTB+ の計算ディレクトリです (2026-09-12 02:59 (JST))

計算の種類: 構造最適化 / 元素: O H / 原子数: 3
計算コード: DFTB+
実行先: プロファイル local (この PC で直接走らせる)
  (プロファイル = ADIT の環境設定に書いた「どこで・どう走らせるか」の組。ADIT の画面の「プロファイル」で選びます)

== このディレクトリのファイル ==
  submit.sh     ジョブスクリプト (計算を走らせる手順を書いたシェルスクリプト。下の「走らせる」で使います)
  spec.json     ADIT で決めた設定一式。ADIT の「ファイル」→「計算設定 (spec.json) を開く…」で読み込むと、同じ設定を画面に戻せます
  analyze.py    解析スクリプト (下の「結果の見方」を参照)
  dftb_in.hsd   DFTB+ の入力 (計算の設定。DFTB+ はこのディレクトリで起動すると自動でこれを読みます)
  geometry.gen  構造 (原子の種類と座標。gen 形式、長さの単位は Å)
  skf/          Slater-Koster ファイル (DFTB+ が使う、元素の組ごとのパラメータ) と LICENSE、README: H-H.skf, H-O.skf, LICENSE, O-H.skf, O-O.skf, README

== 段 ==
  3 段のうち 1 段目 (min: 構造最適化)。spec.json の構造から始めます。

== この PC で走らせる ==
  以下はターミナル (コマンドを 1 行ずつ打ち込んで PC を操作する画面。Windows なら WSL の Ubuntu) で行います。
  1. このディレクトリへ移動します (cd = 作業する場所を変えるコマンド)
       cd /home/<user>/adit/examples/prep_stage1_dftb_stages/stage_01_min
  2. dftb+ が PATH (コマンドを探す場所の一覧) にあるか確かめます
       command -v dftb+
     場所 (例: /home/.../bin/dftb+) が 1 行出れば準備できています。
     何も出なければ、計算ソフトを入れた conda 環境 (ソフトごとに分けたインストール先) を有効にしてから、
     もう一度確かめます。例: conda activate <環境名>  (入れていなければ: conda install -c conda-forge dftbplus)
  3. 走らせます
       bash submit.sh
     画面に何も出なくても動いています。記録は output.log に書かれ、終わると次のコマンドを打てる状態に戻ります。
     途中経過は別のターミナルで  tail -f output.log  (Ctrl+C で表示だけ止まり、計算は続きます)

== 研究室のクラスタで走らせたいとき ==
  クラスタ = 研究室や計算センターが共同で使う計算機の集まり。計算はジョブスケジューラ (PBS や Slurm。計算の順番待ちを
  管理するソフト) に預けて走らせます。預けた 1 件の計算を「ジョブ」と呼びます。
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
       dftbplus = ["<module 名>"]                   # 計算ソフトを使えるようにする module の名前 (クラスタで module avail と打つと一覧が出ます)
     ほかの項目 (投入コマンド、環境変数、実行コマンドなど) は ADIT の README.md の「クラスタで実行する場合」にあります。
  b. ADIT でプロファイルをそれに切り替えて生成し直します。submit.sh が qsub (PBS) / sbatch (Slurm) で預ける
     ジョブスクリプトになり、この README.txt にも、そのクラスタでの投入と確認のコマンドが入ります。
  そのあとの流れは次のとおりです (<...> は自分の値に置き換えます)。
  1. このディレクトリごとクラスタへ写します (手元の PC のターミナルで。scp / rsync = ネットワーク越しにファイルを写すコマンド)
       scp -r /home/<user>/adit/examples/prep_stage1_dftb_stages/stage_01_min <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/
     または
       rsync -av /home/<user>/adit/examples/prep_stage1_dftb_stages/stage_01_min <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/
  2. クラスタにログインして、写したディレクトリへ移動します (ssh = 別の計算機にログインするコマンド)
       ssh <ユーザー名>@<クラスタのホスト名>
       cd <クラスタでの作業ディレクトリ>/stage_01_min
  3. ジョブを預けます (投入)。ジョブ番号が表示されれば受け付けられています
       qsub submit.sh     (PBS のとき)   /   sbatch submit.sh   (Slurm のとき)
  4. 状態を確かめます ($USER は自分のユーザー名に置き換わります)。一覧から消えたら終わっています
       qstat -u $USER     (PBS のとき)   /   squeue -u $USER    (Slurm のとき)
  5. 終わったら結果を手元に取り戻します (手元の PC のターミナルで)
       scp -r <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/stage_01_min <手元の置き場所>/
     または
       rsync -av <ユーザー名>@<クラスタのホスト名>:<クラスタでの作業ディレクトリ>/stage_01_min <手元の置き場所>/

== 結果の見方 ==
  output.log    実行ログ (SCC = 電荷を自己無撞着に決める反復 の様子と、各ステップのエネルギー)
  detailed.out  最後のステップのエネルギーの内訳、Mulliken 電荷 (原子ごとの電荷の目安)、原子にかかる力
  results.tag   全エネルギーなどを、プログラムで読みやすい形で書いたもの
  geom.out.gen / geom.out.xyz   最適化後の構造 (.xyz は分子ビューアで開けます)
  dftb_pin.hsd  省略した設定を既定値で埋めた入力。実際に使われた設定を確かめられます
  analyze.py    ADIT の解析タブと同じ処理で、図と要約を analysis/ に書きます。ADIT が入った Python で走らせます
       python analyze.py

== 引用 ==
  Slater-Koster セット mio-1-1 は CC BY-SA 4.0 で配布されています。論文などでは、セットの README
  (skf/README) に書かれた文献の引用が求められます。DFTB+ 自体の引用先は output.log の先頭に表示されます。

== 来歴 ==
  ADIT 0.0.1 / Python 3.14.7 / ASE 3.29.0 / 生成 2026-09-11T17:59:19+00:00 (UTC)
  パラメータなどのファイルの SHA-256 (中身から計算する指紋。同じ値なら同じファイル。spec.json の provenance にも同じもの):
    d43946902e4d120105025683b54c8696856cbd2b459cfd16df166164774359c8  skf/H-H.skf
    e357795fc4da4e6317fa29a1ffe57c13d17b96d8929fc840648cb2430acd1a6d  skf/H-O.skf
    9348ddfd44da5a127c59141981954746a860ec8e03e0412cf3af7134af0f97e2  skf/LICENSE
    0bc4d6d3c013069ac5964cb58d52e50ea05eb7c5e2097fa76dd8255b535727bd  skf/O-H.skf
    423de839173a667cc1fd75bd111e6dfa46435b874544a5ce4441b84a31b80b43  skf/O-O.skf
    92f2a22a702ed3e8548e376a28c47bbba9fe556df8b735d42d4cbb8f67b63e15  skf/README
  計算コードの版: 走り終えると submit.sh が出力の版の行を code_version.txt に写します
