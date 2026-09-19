ADIT 0.1.0a1 が生成した、フォノン (有限変位) の計算です (dftbplus、phonopy)

  超格子 2x2x2 (16 原子)、変位 1 個、変位の大きさ 0.01 Å
  k 点は分割数の指定 (4x4x4) なので、同じ分割数を超格子にも使っています (k 点の間隔は超格子の倍率だけ細かくなります)。
  phonopy 4.5.0 で作りました (phonopy_disp.yaml)。基本セルの選び方と変位の大きさは、指定のないものは phonopy の既定です。
  状態密度: q 点のメッシュ 10x10x10 (phonopy の既定のテトラヘドロン法)

== 実行したあと ==
  python phonon_collect.py   (ADIT が入った Python で。力を集めて band.yaml を書きます。phonopy で作ったなら phonopy も要ります)
  そのあと adit-analyze <このディレクトリ> でフォノン分散の図 (band.yaml と total_dos.dat を読みます)
  振動数の単位は THz。負の値は虚振動数 (phonopy の慣習)。安定かどうかは ADIT は判断しません。

== この PC で実行する ==
  bash submit.sh   (変位の計算を順に実行します。1 つずつ実行するなら、各ディレクトリで bash submit.sh)

== クラスタで実行する (ADIT は投入しません) ==
  各ディレクトリの submit.sh をクラスタのプロファイルで生成したうえで、このディレクトリで次のように投入します (互いに独立)。
  PBS:   for d in disp-001; do (cd $d && qsub submit.sh); done
  Slurm: for d in disp-001; do (cd $d && sbatch submit.sh); done
