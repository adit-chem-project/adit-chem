ADIT 0.1.0a1 が生成した、弾性定数のための歪みの計算です (dftbplus、原子を動かさない)

  成分 (Voigt の番号 1=xx 2=yy 3=zz 4=yz 5=xz 6=xy): [1, 4]、歪みの大きさ: [-0.005, 0.005] (4〜6 は工学歪み)、計算の数 5 (歪みなしの e0 を含む)

== 実行したあと ==
  python elastic_collect.py   (応力を集め、成分ごとに δ = 0 を含む点で直線を当てた傾き C_ij [GPa] を elastic_constants.csv に書きます)
  値の良し悪し (歪みの大きさが線形の範囲か) は ADIT は判断しません。当てはめの残差は elastic.json に書きます。

== この PC で実行する ==
  bash submit.sh   (歪みの計算を順に実行します。1 つずつ実行するなら、各ディレクトリで bash submit.sh)

== クラスタで実行する (ADIT は投入しません) ==
  各ディレクトリの submit.sh をクラスタのプロファイルで生成したうえで、このディレクトリで次のように投入します (互いに独立)。
  PBS:   for d in e0 e1_-0.005 e1_+0.005 e4_-0.005 e4_+0.005; do (cd $d && qsub submit.sh); done
  Slurm: for d in e0 e1_-0.005 e1_+0.005 e4_-0.005 e4_+0.005; do (cd $d && sbatch submit.sh); done
