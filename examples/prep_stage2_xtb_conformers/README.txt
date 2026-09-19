ADIT 0.1.0a1 が生成した、配座の候補の計算です (xtb)

  元の分子: CCO
  RDKit の ETKDG で 10 個を頼み 10 個を作り (乱数の種 12345)、MMFF94 で最適化 (反復の上限 200)。
  RMSD 0.3 Å 以下 (全原子) を重複として 7 個を外し、3 個を計算にしました。
  力場のエネルギーは並べるだけです (どの配座を使うかは ADIT は判断しません)。

  ディレクトリ  力場のエネルギー [kcal/mol]  最低との差
  conf_001           -1.5171              0.0000
  conf_002           -1.3369              0.1802
  conf_003           -1.3369              0.1802

  conformers.json  全部の配座 (重複として外したものと、その相手と RMSD を含む)
  conformers.sdf   残した配座 (分子ビューアで開けます)
  crest/           CREST の入力 (ADIT は実行しません。crest/README.txt)

== 実行したあと ==
  adit-analyze <このディレクトリ> --scan   (配座ごとの最終エネルギーの表 scan_energies.csv)

== この PC で実行する ==
  bash submit.sh   (配座の計算を順に実行します。1 つずつ実行するなら、各ディレクトリで bash submit.sh)

== クラスタで実行する (ADIT は投入しません) ==
  各ディレクトリの submit.sh をクラスタのプロファイルで生成したうえで、このディレクトリで次のように投入します (互いに独立)。
  PBS:   for d in conf_001 conf_002 conf_003; do (cd $d && qsub submit.sh); done
  Slurm: for d in conf_001 conf_002 conf_003; do (cd $d && sbatch submit.sh); done
