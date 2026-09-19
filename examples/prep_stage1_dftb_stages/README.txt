ADIT 0.1.0a1 が生成した、段階に分けた計算です (dftbplus、3 段階)

  stage_01_min             構造最適化
  stage_02_nvt             MD (NVT、300 K、40 ステップ × 0.5 fs、熱浴 berendsen)  ← 前の段階の最終構造から
  stage_03_nve             MD (NVE、300 K、20 ステップ × 0.5 fs、熱浴 berendsen)  ← 前の段階の最終構造と速度から

各段階のディレクトリは、通常の ADIT の生成したファイルと同じ形です (中の README.txt に、その段階の入力と出力の説明があります)。
2 段階目からの submit.sh は、実行する直前に前の段階の出力から構造 (と速度) を写します。前の段階が終わっていなければ止まります。
そのとき python3 で handoff.py (このディレクトリにあります。標準ライブラリだけで動きます) を使います。

== この PC で実行する ==
  bash submit.sh   (段階を順に実行します。途中の段階が失敗したら、そこで止まります)

== クラスタで実行する (依存付きの投入の例。ADIT は投入しません) ==
  各段階の submit.sh をクラスタのプロファイルで生成したうえで、このディレクトリで次のように投入すると、
  前の段階が正常に終わったときだけ次の段階が走り出します。
  PBS:
    j1=$(cd stage_01_min && qsub submit.sh)
    j2=$(cd stage_02_nvt && qsub -W depend=afterok:$j1 submit.sh)
    j3=$(cd stage_03_nve && qsub -W depend=afterok:$j2 submit.sh)
  Slurm:
    j1=$(cd stage_01_min && sbatch --parsable submit.sh)
    j2=$(cd stage_02_nvt && sbatch --parsable --dependency=afterok:$j1 submit.sh)
    j3=$(cd stage_03_nve && sbatch --parsable --dependency=afterok:$j2 submit.sh)
