ADIT 0.1.0a1 が生成した、比べる計算の組です (種類: adsorption、計算コード: dftbplus)

  ディレクトリ    役割        組成      原子数  周期  電荷  多重度  k 点
  slab            slab        C8        8       yes   0     1       4x4x1
  molecule        molecule    H2O       3       no    0     1       -
  adsorbed        adsorbed    C8H2O     11      yes   0     1       4x4x1

比べる式: ΔE = +1·adsorbed -1·slab -1·molecule   (係数 ν は生成したファイルが正、反応物が負)
  組成と電荷は釣り合っています

== 組の中で値が違う条件 (1 項目。ADIT は止めません。比べてよいかは判断しません) ==
  structure.periodic: slab = true, molecule = false, adsorbed = true
  一部の計算にだけある項目: 7 項目 (例: kpoints, kpoints.density, kpoints.mesh, kpoints.mode, kpoints.resolved_mesh, kpoints.shift。一覧は compare.json の partial)

== 実行したあと ==
  adit-analyze <このディレクトリ> --compare   (compare.json を読み、ΔE と条件の違いの表を書きます)

== この PC で実行する ==
  bash submit.sh   (組の計算を順に実行します。1 つずつ実行するなら、各ディレクトリで bash submit.sh)

== クラスタで実行する (ADIT は投入しません) ==
  各ディレクトリの submit.sh をクラスタのプロファイルで生成したうえで、このディレクトリで次のように投入します (互いに独立)。
  PBS:   for d in slab molecule adsorbed; do (cd $d && qsub submit.sh); done
  Slurm: for d in slab molecule adsorbed; do (cd $d && sbatch submit.sh); done
