# 重い解析: 大きな軌跡を計算機で処理する

ADIT は計算を投入しません。**重い解析も同じ扱い**にしました。

```
ADIT が用意する        人が実行する              ADIT が図にする
msd_worker.py   ──▶   bash msd_run.sh    ──▶   adit-analyze <ディレクトリ>
msd_run.sh           (→ msd_vanhove.json)     (図・要約・報告)
```

フォノン (`phonon_collect.py`) や弾性定数と同じ流儀です。画面の中で数十分かかる計算を回さないので、
**GUI が固まらず、ノート PC のメモリも使い切りません。**

## 何が「重い」のか

`adit-analyze --vanhove` は、実行する前に**算術で見積もります** (実際には計算しません)。

```
変位の数 ≒ τ の本数 × フレーム数 × 原子数 × 0.75
時間 ≒ 変位の数 × 100 ns + τ の本数 × 5 ms     (2026-09-13 に手元で測った係数)
座標のメモリ = フレーム数 × 原子数 × 24 バイト
```

見積もりが **60 秒** (`--heavy-limit` で変更可) を超えると、その場では計算せず、実行用のファイルを置きます。

| 系 | 変位の数 | 見積もり | どうなるか |
|---|---|---|---|
| 18 原子 × 1,000 フレーム × τ 100 本 | 135 万 | 0.6 秒 | その場で計算 |
| 18 原子 × 11,001 フレーム × τ 100 本 | 1,485 万 | 2 秒 | その場で計算 |
| 1,000 原子 × 10,000 フレーム × τ 100 本 | 7.5 億 | 76 秒 | ファイルを置く |
| 18 原子 × 11,001 フレーム × τ 3,000 本 | 4.5 億 | 45 秒 + 当てはめ 15 秒 | ファイルを置く |

## 使い方

```bash
adit-analyze run/ --msd --vanhove              # 見積もって、重ければ msd_worker.py と msd_run.sh を置く
bash run/msd_run.sh                             # 人が実行する (クラスタならジョブとして投入)
adit-analyze run/ --msd --vanhove              # 2 回目。msd_vanhove.json を読んで図にする
adit-analyze run/ --msd --vanhove --vanhove-here   # 小さい系なら、その場で計算してもよい
```

`msd_worker.py` は **numpy と scipy だけで動きます** (ADIT を import しません。`handoff.py` と同じ決まりで、
試験 `tests/test_heavy_setup.py` が見張っています)。クラスタに ADIT が入っていなくても走ります。

```bash
# 1 行 = 1 原子 (ラベル x y z) の独自形式をそのまま読む場合
python msd_worker.py Na.dat --natoms 18 --cell dftb.inp --dt 10 --taus 100 --out msd_vanhove.json

# ASE が読める形式 (extxyz、vasprun.xml など) はそのまま
python msd_worker.py trajectory.extxyz --dt 10 --species Na
```

## 出てくるもの

`msd_vanhove.json` に入るのは**数値だけ** (判定はしません)。

| 鍵 | 中身 |
|---|---|
| `d_per_atom_cm2_s` | 原子 1 個ごとの D。時間原点を全部使った MSD の傾きから |
| `d_from_mean_msd_cm2_s` | 全原子の平均の MSD から出した D |
| `vanhove.d_fit_A2_fs` | 変位の分布にガウスの形を当てはめた D(τ) |
| `vanhove.d_direct_A2_fs` | 同じ変位から ⟨r²⟩/(6τ) で出した D(τ)。ガウスなら上と一致します |
| `vanhove.alpha2` | 非ガウス因子 α₂ = 3⟨r⁴⟩/(5⟨r²⟩²) − 1 (Rahman 1964)。0 ならガウス |
| `vanhove.truncated_from_fs` | 最小像の上限で変位が頭打ちになり始めた遅れ時間 (なければ null) |

**最小像の頭打ちに注意。**`--vanhove-displacement mic` (既定) は折り返した座標の最小像を使うので、
**セルの最小の幅の半分**より大きい変位を測れません。変位の 1 % がその 9 割を超えたら、その遅れ時間を
`truncated_from_fs` に記録し、要約でも注意します。`--vanhove-displacement unwrapped` なら頭打ちしません。

## 判定しないこと

- どの D を採るか (MSD の傾き / 分布の当てはめ / ⟨r²⟩/(6τ))
- どこからを「非ガウス」と呼ぶか
- どの原子が「速い」か
- 当てはめ範囲が妥当か (既定は遅れ時間の 10〜50 %。`--msd-fit` で変えられます)
