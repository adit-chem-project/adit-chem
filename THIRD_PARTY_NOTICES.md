# 同梱している第三者のファイルと、そのライセンス

**ADIT 自身のコード (`src/adit/`、`tests/*.py`、文書) は MIT です** (`LICENSE`)。
一方、`examples/` と `tests/data/` には、**他のプロジェクトから取ってきたファイル**が入っています。
それらは**元のライセンスのまま**で、MIT にはなりません。この文書がその一覧です。
各ディレクトリの `SOURCES.md` に、取得日・URL・sha256 まで書いてあります。

ライセンスの全文は `licenses/` にあります。

## PyPI に出す配布物 (adit-chem) には入っていません

`pip install adit-chem` で入るもの (wheel / sdist) は **ADIT 自身のコードだけ**です。
下の第三者のファイルは、GitHub のリポジトリにだけ入っています。

## 計算コードの配布物から取ったファイル

| ファイル | 出典 | ライセンス | 全文 |
|---|---|---|---|
| `examples/lammps_cu/Cu_u3.eam` と、その写し `examples/lammps_cu_nvt_generated/Cu_u3.eam`、`examples/plumed_lammps_cu_generated/Cu_u3.eam` | LAMMPS の `potentials/` | GPL-2.0 | `licenses/GPL-2.0.txt` |
| `examples/cp2k_h2o_generated/BASIS_adit`、`POTENTIAL_adit` ほか CP2K の例 | CP2K の `data/` (BASIS_MOLOPT、GTH_POTENTIALS の一部) | GPL-2.0-or-later | `licenses/GPL-2.0.txt` |
| `examples/gromacs_spce/conf.gro` と、その写し `examples/gromacs_spce_em_generated/`、`examples/gromacs_spce_nvt_generated/`、`examples/openmm_spce_nvt_generated/` の `conf.gro` | GROMACS の `share/gromacs/top/spc216.gro` | LGPL-2.1 | `licenses/LGPL-2.1.txt` |
| `examples/prep_stage1_dftb_continue_nve/skf/`、`examples/prep_stage1_dftb_stages/*/skf/`、`examples/prep_stage2_dftb_adsorption/*/skf/`、`examples/prep_stage3_dftb_elastic/*/skf/`、`examples/prep_stage3_dftb_neb_images/*/skf/`、`examples/prep_stage3_dftb_phonons_ase/*/skf/`、`examples/prep_stage3_dftb_phonons_phonopy/*/skf/` の `.skf`、`LICENSE`、`README` (合計 76 個の `.skf`) | dftb.org の Slater-Koster セット mio-1-1 (https://dftb.org/parameters/download.html) | CC BY-SA 4.0 | `licenses/CC-BY-SA-4.0.txt` (各 `skf/LICENSE` にも同じ全文) |
| `examples/xtb_water_generated` の一部 | xtb の配布物 | LGPL-3.0 | (xtb の配布物を参照) |
| `examples/qe_si_generated` の一部 | Quantum ESPRESSO の配布物 | GPL-2.0 | `licenses/GPL-2.0.txt` |

同じディレクトリの `topol.top` (`examples/gromacs_spce/` とその写し) は、GROMACS 同梱の力場を `#include` するだけの短いファイルで、ADIT 側で書いたものです。

**ADIT はこれらのファイルを作りません。**利用者が用意したものを写すだけで、例は「動く形」を示すために置いています。

## 文書・レシピから取ったもの

| ファイル | 出典 | ライセンス | 全文 |
|---|---|---|---|
| `examples/vasp_h2o`、`examples/vasp_cd_si` | VASP 公式 wiki の例 (原文のまま転載) | GNU Free Documentation License 1.2 | `licenses/GNU-FDL-1.2.txt` |
| `examples/dftb_tio2`、`examples/water`、`examples/prep_stage2_dftb_adsorption` | DFTB+ の公式レシピ集 | CC BY-SA 4.0 | `licenses/CC-BY-SA-4.0.txt` |

## 書式を確かめるために置いた、他のプロジェクトの試験ファイル

| ファイル | 出典 | ライセンス |
|---|---|---|
| `tests/data/gaussian_dvb_raman.out.gz`、`gamess_dvb_ir.out.gz`、`qchem_dvb_raman.out.gz`、`orca_dvb_raman.out.gz` | [cclib](https://github.com/cclib/cclib) の `data/` | BSD-3-Clause |
| `tests/data/PROCAR.simple`、`PROCAR.new_format_5.4.4.gz`、`WAVEDER.gz` | [pymatgen](https://github.com/materialsproject/pymatgen) の `test-files/` | MIT |

**これらは「読み取りが本当に合っているか」を、実物の書式で確かめるために置いています。**
出典と確かめた内容は `tests/data/SOURCES.md` にあります。

## データベースから取得した構造 (利用者が取得するもの)

「データベースから取得」(`adit-gen --fetch`) で手元に来るファイルは ADIT の配布物ではなく、**利用者がそのサイトから取得したもの**です。
ADIT は取得元・ID・取得日時・応答の SHA-256 とライセンスの表記を spec.json (`provenance.fetched_structure`) と README.txt に記録し、
取得したファイルを出力ディレクトリに写します。**公開・再配布するときは、その表記の条件に従ってください。**
各サイトの表記は 2026-09-19 に確認しました。

| 取得元 | 表記 | 確認したページ |
|---|---|---|
| PubChem | NCBI の方針: 米国政府が作った情報はパブリックドメイン。提供元 (submitter) が特許・著作権などを主張する場合があり、NCBI は個別の許諾を出せないとしている。PubChem の案内も「提供元のライセンス表記を確かめること」としている | https://www.ncbi.nlm.nih.gov/home/about/policies/ 、https://pubchem.ncbi.nlm.nih.gov/docs/downloads |
| COD (Crystallography Open Database) | 「All data in the COD and the database itself are dedicated to the public domain and licensed under the CC0 License」 | https://www.crystallography.net/cod/ |
| Materials Project | 「By downloading Content from Materials Project, User agrees to accept the Creative Commons Attribution 4.0 License … provided proper attribution is given to the Materials Project」 | https://next-gen.materialsproject.org/about/terms (ブラウザ以外からは 403 になるので、同じ文面のある https://legacy.materialsproject.org/terms で確認) |
| OPTIMADE (横断検索) | 提供元ごとに異なる。ADIT は提供元の名前とホームページを記録するので、そこで条件を確かめてください | https://providers.optimade.org/ |

PubChem の利用方針は「1 秒に 5 回まで」の要求で、ADIT は 1 回の取得で 3 回 (名前 → CID、SDF、化合物名) までしか要求しません。

## 同梱していないもの (配れないもの)

- **VASP の POTCAR**: ライセンス保持者以外に配布できません。ADIT は**書きも写しもしません**
- **Slater-Koster パラメータ (slakos/)**: 利用者が dftb.org から入れます。セットごとの置き場所 `slakos/` はリポジトリには入れていません (`.gitignore`)。上の表にある `examples/prep_stage*/**/skf/` の mio-1-1 の写しだけが例外です
- **擬ポテンシャル (UPF、psp8 など)**: 同上
- **計算コード本体** (DFTB+、VASP、Gaussian、ORCA ほか): 入れていません

## 実行ファイル (.exe / .app) に入るもの

PyInstaller で作る配布物には、**日本語のフォント (Noto Sans CJK JP、OFL-1.1)** が入ります。
**その全文 `licenses/OFL-1.1-NotoSansCJK.txt` も一緒に入ります** (`packaging/adit.spec` がそのように組み立てます。
全文が無いときは、フォントを入れずに組み立てが止まります)。

同じ配布物に、`LICENSE`、この `THIRD_PARTY_NOTICES.md`、`licenses/*.txt` も入ります
(受け取った人が、手元で権利の表示を読めるように)。

| 入るもの | ライセンス | 全文 |
|---|---|---|
| Noto Sans CJK JP (フォント) | SIL Open Font License 1.1 | `fonts/OFL-1.1-NotoSansCJK.txt` |
| Python、ASE、pydantic、Jinja2、matplotlib、SciPy、NumPy | それぞれの配布物のライセンス (PSF / LGPL-2.1+ / MIT / BSD-3 ほか) | 各パッケージの配布物に含まれます |
| tomli-w | MIT | 同上 |
| PySide6 (画面のとき) | **LGPL-3.0** | Qt / PySide6 の配布物を参照 |
| PyQtDarkTheme-fork (画面のテーマ) | MIT | 各パッケージの配布物に含まれます |
| pyte (ワークスペースのターミナルの画面) | **LGPL-3.0** | 同上 |
| ptyprocess (Linux / macOS のターミナル) | ISC | 同上 |
| pywinpty (Windows のターミナル) | MIT | 同上 |
| RDKit (SMILES と Draw。入れて作ったときだけ) | BSD-3-Clause | 同上 |
| spglib (空間群。入れて作ったときだけ) | BSD-3-Clause | 同上 |

**PySide6 (Qt) と pyte は LGPL-3.0 です。**画面つきの実行ファイルを配るときは、LGPL の求めに応じて
「利用者が Qt を差し替えられること」を満たす必要があります (PyInstaller の 1 ディレクトリ形式なら、
Qt の共有ライブラリが別ファイルとして入るので差し替えられます。**1 つにまとめた形式 (onefile) は避けてください**)。
`packaging/adit.spec` の `ONEFILE` は False のままにしてあります。

## 商標について

VASP、Gaussian、ORCA、GAMESS、Q-Chem、GRRM、Materials Studio などは、それぞれの権利者の商標または登録商標です。
ADIT はこれらのソフトの**入力ファイルを書き、出力を読むだけ**で、これらの提供元とは無関係であり、
提供元からの承認も受けていません。
