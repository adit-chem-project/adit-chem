# ADIT を Python から使う

コマンド (`adit-gen` / `adit-analyze` / `adit-report` / `adit-convert`) ではなく、**自分のスクリプトから
部品として呼ぶ**ときの地図です。

**約束**: ここに書いた名前は `adit.<モジュール>` から直接 import できます。内部用の関数 (`_` で始まるもの) は
いつ変わるか分からないので使わないでください。

## 1 最短の流れ

```python
from pathlib import Path

from adit.config import load_config
from adit.spec import CalculationSpec
from adit.project import write_project
from adit.analysis.report import AnalysisOptions, run_analysis
from adit.report import load_run_report, methods_markdown

spec = CalculationSpec.load("spec.json")        # 条件を読む
cfg = load_config()                             # ~/.config/adit/cluster.toml (ADIT_CONFIG で差し替え)
write_project(spec, cfg, Path("out/run1"))      # 入力・submit.sh・README.txt を書く (投入はしない)

# --- ここで人が計算を実行する (ADIT はジョブを投入しない) ---

res = run_analysis("out/run1", AnalysisOptions(msd=True, out_dir=Path("scratch/analysis1")))
print(res.summary_text())                       # 図と表は out_dir に書かれる
print(methods_markdown([load_run_report("out/run1")], "ja"))   # 論文の「方法」の節
```

## 2 主な関数の一覧

| したいこと | 呼ぶもの | 戻り値 |
|---|---|---|
| 条件を読む・書く | `adit.spec.CalculationSpec.load / save / from_json` | `CalculationSpec` (pydantic) |
| 条件の間違いを調べる | `adit.validate.validate(spec, cfg)` | `list[ValidationError]` (空なら問題なし) |
| 入力を組み立てる (書かない) | `adit.project.build_project(spec, cfg)` | `ProjectFiles` (`.texts` / `.copies`) |
| 入力を書く | `adit.project.write_project(spec, cfg, out_dir)` | 書いたファイルの `list[Path]` |
| 1 つの値を振る | `adit.scan.write_scan(spec, cfg, out, scans)` | 作ったディレクトリ |
| 構造を差し替えて一括 | `adit.structures_batch.write_structures / write_enumerated` | 作ったディレクトリ |
| 出力を読む | `adit.analysis.readers.load_run(run_dir)` | `RunData` (エネルギー・温度・軌跡・最終構造) |
| 解析する (図と表) | `adit.analysis.report.run_analysis(run_dir, opts)` | `AnalysisResult` |
| 走り終えたかを見る | `adit.results.summarize_run(run_dir)` / `not_run_yet` / `failure_lines` | `RunSummary` / `bool` / `list[str]` |
| 報告・再現パッケージ | `adit.report.load_run_report / methods_markdown / conditions_csv / results_csv / write_bundle / check_lines` | 文字列・パス |
| 既存の入力を読み戻す | `adit.native_import.import_native(source, code=None)` | `NativeImportResult` (`.spec` は読めないとき `None`、`.report()` に読み取れなかった欄) |
| 構造を作る・変える | `adit.structure.from_smiles / from_file / from_spacegroup`、`adit.convert` | `ase.Atoms` |
| 原子を選ぶ | `adit.analysis.select.select(atoms, "element O and z < 10")` | `numpy` の添字 (0 始まり) |
| 距離・角・二面角・RMSD・Rg の時系列 | `adit.analysis.geometry_series.*` | `numpy` の配列 |
| 配位数・中心対称性・Steinhardt・かたまり・S(q) | `adit.analysis.local_order.*` | `numpy` の配列 / `dict` |
| 水素結合を数える | `adit.analysis.hbond.count_series(frames, distance, angle)` | 1 フレームごとの本数 (しきい値は必須) |
| CREST の配座を読む・重みを付ける | `adit.analysis.crest.analyze_crest(crest_dir, out_dir, temperature_k=None)` / `read_ensemble` / `read_energies` / `boltzmann_weights(relative_ev, degeneracy, temperature_k)` | 配座の表 (`dict`) / `list[(E [Eh], Atoms)]` / `dict` / 重みの配列 (温度を渡したときだけ) |
| 速度自己相関と振動スペクトル | `adit.analysis.vacf.vacf(velocities, dt_fs, *, max_lag_fraction=0.5, source="velocities", window="hann")` (`velocities` は (フレーム, 原子, 3) の配列、Å/fs) | `VacfResult` |
| 有効質量 | `adit.analysis.effective_mass.at_band_edges(kdist, energies, fermi)` | `list[EffectiveMass]` |
| 射影バンド (QE の projwfc.x) | `adit.analysis.projected_bands.read_filproj(path)` | `ProjectedBands` |
| 光学 (QE の epsilon.x) | `adit.analysis.optical.read_epsilon(run_dir)` | `Optical` (n, k, 吸収係数, 反射率, EELS) |
| Bader 電荷 (外部の bader の ACF.dat) | `adit.analysis.bader.read_acf(path)` | `BaderCharges` (分割は計算しない) |
| 重い解析を計算機へ出す | `adit.analysis.heavy_setup.write_job(run_dir, ...)` | 実行するスクリプトのパス |
| VASP の PROCAR | `adit.analysis.procar.read_procar(path)` | `Procar` (射影・固有値・k 点) |
| Gaussian・GAMESS・Q-Chem の出力 | `adit.analysis.readers_qc.read_gaussian / read_gamess / read_qchem` | `QcOutput` (エネルギー・構造・振動数・赤外・ラマン・電荷) |
| OpenMX の出力 | `adit.analysis.readers_openmx.read_openmx(path)` | `OpenmxOutput` |
| VASP の WAVEDER | `adit.analysis.waveder.read_waveder(path)` / `dielectric_imag(...)` | `Waveder` / eps2 の配列 |
| 変位・局所ひずみ | `adit.analysis.displacement.displacements / local_strain` | `numpy` の配列 / `LocalStrain` |
| 自由エネルギー面・同時分布 | `adit.analysis.free_energy.free_energy_surface / distribution_2d / read_colvar` | `FreeEnergySurface` / `Distribution2D` |
| Voronoi の体積 | `adit.analysis.voronoi.voronoi(positions, cell)` | `VoronoiResult` |
| 溶媒接触表面積 | `adit.analysis.sasa.sasa(positions, radii, probe)` | `SasaResult` |
| フレームの主成分分析・分類 | `adit.analysis.frames_pca.pca / kmeans` | `PcaResult` / `ClusterResult` |

## 3 例外

すべての例外は `adit.errors.AditError` を継承します。**外部要因 (ファイルが無い、メモリ不足) と ADIT 起因の
失敗を 1 行で分けられる**ようにするためです。

```python
from adit.errors import AditError, AditValueError

try:
    write_project(spec, cfg, out)
except AditError as ex:          # ADIT が「これでは動かない」と判断して止めたもの
    print(f"ADIT: {ex}")
```

`AditValueError` は `AditError` と `ValueError` の両方を継承するので、`except ValueError` で書いた
既存のコードもそのまま動きます。

| 例外 | いつ | モジュール |
|---|---|---|
| `GenerationError` | 入力を組み立てられない | `adit.codes.base` |
| `OutputNotEmpty` | 出力先に既にファイルがある | `adit.project` |
| `StructureError` | 構造を作れない・読めない | `adit.structure` |
| `ScanError` | 振る値の指定が読めない | `adit.scan` |
| `StructuresError` / `EnumerateError` | 一括生成・置換基の列挙 | `adit.structures_batch` / `adit.enumerate_r` |
| `TemplateError` | 雛形の保存・読み出し | `adit.templates` |
| `ConfigMissing` / `ConfigError` | 環境設定ファイル | `adit.config` |
| `ReportError` | 報告を作れない (spec.json が無いなど) | `adit.report` |
| `TrajectoryTooLarge` | 軌跡がメモリの上限を超える | `adit.analysis.trajectory` |
| `TransportError` / `VolumetricError` / `XrdError` / `RateError` | 粘度・体積データ・粉末回折・速度定数 | `adit.analysis.*` |
| `HandoffError` | 段階の引き継ぎ (**このモジュールだけ `AditError` を継承しない**) | `adit.handoff` |

`adit/handoff.py` は計算機の上に写して単体で実行するので、ADIT を import しません。例外もその場で定義しています。

## 4 覚えておくこと

- **ジョブは投入しません。** `write_project` は `submit.sh` を書くだけで、実行するのは人です
- **解析は既定で計算のディレクトリに書きます。** 汚したくないときは `AnalysisOptions(out_dir=...)` を指定します
- **環境設定ファイルは 1 つです。** 試験などで差し替えるときは環境変数 `ADIT_CONFIG` を使います
 (`load_config(path)` に直接渡しても構いません)
- **言語は `adit.lang.set_language("ja" / "en")`** で切り替わります。エラー文も表の見出しも同じ設定を見ます
- **GUI を import しないでください。** `adit.gui` は PySide6 を要求します。ウェブ版・CLI・解析は import しません
