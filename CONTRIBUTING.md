# 開発に加わる人へ

## 動かす

```bash
python -m pip install -e ".[dev,gui,smiles,analysis]"
pytest -q                     # 全件 (GUI は画面なしで走ります)
```

計算ソフト (DFTB+、xtb、Quantum ESPRESSO、CP2K、LAMMPS、GROMACS…) が PATH にあると、
実際に実行する試験も動きます。無ければその試験だけ skip されます。

## ADIT が守っている 3 つの原則

1. **画面にロジックを置きません。**画面は条件 (`CalculationSpec`) を組み立てるだけで、
   入力ファイルを作るのは `codes/` と `project.py` です
2. **ジョブは投入しません。**`qsub` / `sbatch` は人が打ちます。手元で直接実行する方法だけがあります
3. **利用者の化学的判断に口を出しません。**判定するのは「このままでは動かない」機械的な事実だけです。
   収束したか、値が妥当かは言いません

## 書くときの決まり

- **キーワードと書式は、公式文書か実物の出力で確かめてから書きます。**記憶で書きません
  (確かめたものは `tests/data/SOURCES.md` に出典を残します)
- **3 つの入口 (コマンド・デスクトップ版・ウェブ版) に同じ機能を出します。**
  解析の条件を足したら `gui/analysis_fields.py` + `gui/panels/analysis_panel.py` + `gui/help.py` +
  `gui/i18n.py` + `web/templates/analysis.html` を同時に直します (見張りは `tests/test_analysis_ui_parity.py`)
- **画面の文言は日本語と英語の両方で書きます。**`gui/help.py` の説明と `gui/i18n.py` の対訳を必ず足します
- 数値を出すときは**単位と出典**を添えます

## 計算コードを追加する

`codes/<name>.py` に `InputGenerator` を実装して `register()` し、`spec.py` に Method を追加し、
デスクトップ版のパネル (`gui/panels/`) とウェブ版のフォームを追加します。

## 試験

新しい機能には、**答えの分かる小さな例**での試験を付けてください
(例: FCC で配位数 12、自由電子で m*/m_e = 1、孤立原子の表面積 4πr²)。
