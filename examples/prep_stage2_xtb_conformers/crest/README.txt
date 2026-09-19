CREST の入力 (ADIT は実行しません)

  struct.xyz   出発の構造 (力場のエネルギーが最も低い配座)
  run_crest.sh crest struct.xyz --gfn2 --chrg 0 --uhf 0 -T 1
               電荷・不対電子の数・溶媒・並列数は spec.json の xtb の設定から写したものです

インストール方法: conda install -c conda-forge crest   (xtb と同じ conda-forge)
結果: crest_conformers.xyz (重複を除いた配座)、crest_best.xyz (最もエネルギーの低いもの)、crest.energies (相対エネルギー)。文書: https://crest-lab.github.io/crest-docs/
