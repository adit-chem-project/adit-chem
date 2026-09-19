# 試験に使う実物のファイルの出典

すべて、書式を**実物で**確かめるために置いてあります。

| ファイル | 出典 | ライセンス | いつ・何を確かめたか |
|---|---|---|---|
| `si.proj.projwfc_up` | 手元で Quantum ESPRESSO 7.x の projwfc.x を Si 2 原子で実行 | — (自分で作った) | 2026-09-13。filproj の並び |
| `epsr_si.dat` / `epsi_si.dat` / `eels_si.dat` | 手元で QE の epsilon.x を Si (ノルム保存の擬ポテンシャル) で実行 | — (自分で作った) | 2026-09-13。列の並びと、-Im(1/eps) が QE 自身の値と一致すること |
| `bader_ACF.dat` | 手元で bader 1.0.5 を、電子数 6 と 2 のガウス関数を置いた cube で実行 | — (自分で作った) | 2026-09-13。ACF.dat の列と、分けられた電子数 (5.9999 / 1.9987) |
| `PROCAR.simple` / `PROCAR.new_format_5.4.4.gz` | [pymatgen](https://github.com/materialsproject/pymatgen) の `test-files/io/vasp/outputs/` | MIT | 2026-09-13。素の形式と lm 分解 + 位相 (VASP 5.4.4)。読んだ値は pymatgen の `Procar` と一致 |
| `WAVEDER.gz` | [pymatgen](https://github.com/materialsproject/pymatgen) の `test-files/io/vasp/outputs/WAVEDER` を gzip したもの | MIT | 2026-09-13。Fortran の書式なしレコードの並び。読んだ配列は pymatgen の `Waveder` と完全一致 |
| `gaussian_dvb_raman.out.gz` / `gamess_dvb_ir.out.gz` / `qchem_dvb_raman.out.gz` / `orca_dvb_raman.out.gz` | [cclib](https://github.com/cclib/cclib) の `data/` (Gaussian 16、GAMESS-US 2018、Q-Chem 5.4、ORCA 6.0 の実際の出力) | BSD-3 | 2026-09-13。エネルギー・構造・振動数・赤外の強度・**ラマン活性**・Mulliken 電荷。読んだ値は cclib 自身の読み取りと一致 (GAMESS の赤外だけ単位の換算が要る) |
| `water_hessian.out` | DFTB+ の実行 (2026-09-10) | — (自分で作った) | hessian.out の並び |

| `fetch/pubchem_water_cids.json` / `fetch/pubchem_962_3d.sdf` / `fetch/pubchem_962_props.json` | PubChem PUG REST の実際の応答 (`compound/name/water/cids/JSON`、`compound/cid/962/SDF?record_type=3d`、`compound/cid/962/property/MolecularFormula,IUPACName,Title/JSON`) | NCBI の方針 (米国政府の著作物はパブリックドメイン) | 2026-09-19。名前 → CID、3D SDF、化合物名の応答の形 |
| `fetch/cod_1000041.cif` / `fetch/cod_1000041_meta.json` | COD の実際の応答 (`https://www.crystallography.net/cod/1000041.cif`、`result?id=1000041&format=json`) | CC0 1.0 | 2026-09-19。CIF と JSON の項目名 (formula, title, authors, doi …) |
| `fetch/providers.json` / `fetch/idx_oqmd_links.json` / `fetch/oqmd_si.json` / `fetch/oqmd_4061352.json` | OPTIMADE の実際の応答 (`https://providers.optimade.org/providers.json` を 3 提供元に間引いたもの、`index-metadbs/oqmd/v1/links`、`oqmd.org/optimade/v1/structures?filter=chemical_formula_reduced="Si"` を 2 件に間引いたもの、その 1 件目を単体の応答の形にしたもの) | providers.json は OPTIMADE 連合 (MIT)、OQMD の項目は OQMD の規約 | 2026-09-19。提供元一覧 → links → structures の 3 段の形と、lattice_vectors / species / cartesian_site_positions |
| `fetch/mp_summary_mp-149.json` | **実物ではない。**Materials Project の OpenAPI (`https://api.materialsproject.org/openapi.json` の `TypedStructureDict` / `TypedLatticeDict` / `TypedSiteDict`) に合わせて自分で書いた応答 (API キーが無く実物を取れなかった) | — (自分で作った) | 2026-09-19。`data[0].structure.lattice.matrix` と `sites[].species[].element` / `xyz` の形だけ |

**開発機に VASP は無いので、VASP の出力は公開されている本物のファイルで確かめています。**

## リポジトリに置いていないが、実物で確かめたもの

ライセンスが ADIT (MIT) と合わないので**置いていません**。書式は 2026-09-13 に実物で確かめました。

| 何 | 出典 | ライセンス | 試験の実行方法 |
|---|---|---|---|
| OpenMX の `.out` | OpenMX 付属の `work/input_example/Methane.out`・`H2O.out` ([FermiQ/openmx-square](https://github.com/FermiQ/openmx-square) で公開) | GPL-3 | `ADIT_OPENMX_OUT=<file> pytest tests/test_readers_openmx.py` |
| VASP の `vasprun.xml` + `WAVEDER` (誘電関数の突き合わせ) | pymatgen の `fixtures/reproduce_eps` | MIT だが 2.7 MB と大きい | 手元でだけ実施 (結果は `analysis/waveder.py` の説明に記録) |
