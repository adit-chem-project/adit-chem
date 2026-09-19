# 対応ソフトウェア

**計算コードごとに、ADIT がどこまで対応しているかの一覧です。**
「対応している」とは、その計算の入力ファイルを ADIT が作れる、という意味です。
計算コード本体・力場・基底関数・擬ポテンシャルは同梱していません (利用者が用意します)。
計算の実行と投入も行いません。

対応している計算の種類は、README の対応表が一覧です。このページは、**表に収まらない細かい範囲と制限**を書いています。

| コード / Code | 対応する計算 / Supported task | `method` に必要な値 / Required `method` fields | 主な入力 / Main input |
|---|---|---|---|
| Gaussian | 一重項分子の一点計算 / Singlet molecular single point | `theory`, `basis` | `gaussian.gjf` |
| US GAMESS | 一重項分子の RHF 一点計算 / Singlet molecular RHF single point | `gbasis` (`STO`, `N21`, `N31`, `N311`), `ngauss` | `gamess.inp` |
| Q-Chem | 一重項分子の一点計算 / Singlet molecular single point | `theory` (`METHOD`), `basis` (`BASIS`) | `qchem.in` |
| GRRM17 | 一重項分子の MIN または FREQ / Singlet molecular MIN or FREQ | `theory`, `basis`; `task.type` は `geometry_optimization` または `vibrations` / set `task.type` accordingly | `grrm.com` |
| OpenMX | 中性一重項の非周期・3 次元周期の一点計算 / Neutral singlet, nonperiodic or fully periodic single point | `data_path`, 元素ごとの `pao`, `vps`, `valence`, `xc`, `energycutoff_ry` / per-element PAO, VPS and valence | `openmx.dat` |
| Amber (`sander`) | 非周期・3 次元周期のエネルギー最小化 / Nonperiodic or fully periodic energy minimization | `topology_file` (テキスト `prmtop`), `coordinates_file` (テキスト `rst7`/`inpcrd`), `cutoff_ang`; 非周期系では `igb` を明示 / explicitly set `igb` for nonperiodic systems | `amber.in` と外部ファイル / plus supplied files |
| NAMD 3 | 非周期・3 次元周期の NVE MD / Nonperiodic or fully periodic NVE MD | `structure_file` (`.psf`), `coordinates_file` (`.pdb`), `parameter_files`, `exclude`, `one_four_scaling`, `cutoff_ang`, `pairlistdist_ang`; switching 使用時は `switchdist_ang` / `switchdist_ang` when switching is enabled | `namd.conf` と外部ファイル / plus supplied files |

`spec.json` の既存の `method` を次の形に置き換えます。下の文字列と数値は**書式を示す例であり、推奨条件ではありません**。利用する基底・擬ポテンシャル・力場が実際の実行環境にあるかは利用者が確認してください。

Replace the existing `method` object in `spec.json` with one of the shapes below. The names and numbers **illustrate syntax, not recommended settings**. Check that the selected basis, pseudopotentials, and force field are available on the execution host.

```json
{"code":"gaussian","theory":"HF","basis":"6-31G(d)"}
{"code":"gamess","gbasis":"STO","ngauss":3}
{"code":"qchem","theory":"HF","basis":"6-31G"}
{"code":"grrm","theory":"HF","basis":"6-31G"}
{"code":"openmx","data_path":"../DFT_DATA","pao":{"H":"H5.0-s1p1"},"vps":{"H":"H_CA19"},"valence":{"H":1},"xc":"LDA","energycutoff_ry":220}
{"code":"amber","topology_file":"system.prmtop","coordinates_file":"system.rst7","cutoff_ang":9,"igb":0}
{"code":"namd","structure_file":"system.psf","coordinates_file":"system.pdb","parameter_files":["force.prm"],"exclude":"scaled1-4","one_four_scaling":1,"cutoff_ang":10,"pairlistdist_ang":12,"switching":true,"switchdist_ang":8}
```

各行は独立した `method` オブジェクトの例です。OpenMX の例は H を含む系の書式のみを示します。実際には**構造にある全元素**の `pao`・`vps`・`valence` が必要です。周期系なら Spec の `kpoints` も指定し、shift はゼロにします。`DATA.PATH` は**計算を実行する場所**から読める必要があり、ADIT は PAO/VPS を同梱しません。OpenMX の初期スピンは元素ごとに指定された価電子数を半分ずつ両スピンに置く中性・非スピン分極の形です。磁性・帯電系は未対応です。

Each line is a separate example `method` object. The OpenMX example only shows the syntax for a system containing H. In a real calculation, supply `pao`, `vps`, and `valence` for **every element in the structure**. Periodic systems also need `kpoints` in the spec, with zero shift. `DATA.PATH` must be readable **where the calculation runs**; ADIT does not package PAO/VPS data. This subset splits the user-specified valence count equally between initial up/down spins for a neutral, non-spin-polarized system. Magnetic and charged systems are not supported.

Amber と NAMD は、利用者が完成させたトポロジーと座標をコピーします。ADIT は原子型・電荷・力場を作りません。Spec と座標ファイルの原子数および座標の順序が一致するか、1 mÅ の許容差で検査します。PSF/prmtop の原子型と力場の内容、Amber の周期セル、NAMD のパラメータの完全性は利用者が確かめてください。Amber は `task.max_steps` を `maxcyc` に書きますが、共通 Spec の `optimizer` と力の閾値を Amber の条件に変換しません。NAMD は共通 MD の `steps`・`temperature_k`・`timestep_fs`・`dump_interval` を書き、温度から初速度を作る NVE のみです。

Amber and NAMD copy topology and coordinates prepared by the user; ADIT does not assign atom types, charges, or force fields. It checks atom count and coordinate order against the spec, with a 1 mÅ tolerance. You must still check the PSF/prmtop atom types and force field, Amber periodic cell, and completeness of NAMD parameters. Amber writes `task.max_steps` as `maxcyc` but does not translate the common `optimizer` or force threshold into Amber settings. NAMD writes the common MD `steps`, `temperature_k`, `timestep_fs`, and `dump_interval`; this subset only runs NVE with initial velocities drawn from the specified temperature.

GRRM17 のこの入口は Gaussian 連携の MIN/FREQ だけです。GRRM とライセンス済みの Gaussian を用意し、GRRM の公式手順に従って連携コマンドを設定してください。共通 Spec の最適化器・最大反復数・力の閾値は GRRM のオプションに移しません。異なるコードへ同じ手法名・基底名を入力しても、実装や既定値を含む計算条件全体が一致するわけではありません。

This GRRM17 entry point covers only MIN/FREQ with Gaussian integration. Provide GRRM and licensed Gaussian and configure their linkage as described in the GRRM manual. The common optimizer, maximum iteration count, and force threshold are not translated into GRRM options. Matching method and basis names across codes does not make all implementation details and defaults equivalent.

コード間変換では、適用しない共通条件を `conversion.json` の `not_applied_by_target` に記録します。GRRM17 と Amber の最小化条件、NAMD の NVE で使わない熱浴・圧力条件、Amber/NAMD で外部トポロジーに委ねる電荷とスピン多重度が対象です。変換後の入力と元の条件を照合してください。

Cross-code conversion records common settings that are not applied in `conversion.json` under `not_applied_by_target`. This covers GRRM17 and Amber minimization settings, thermostat/pressure fields unused by NAMD NVE, and charge/spin fields left to external Amber/NAMD topologies. Compare the generated input against the original conditions.

計算コード本体やライセンスのない環境では、生成ファイルの構文と機械的な検査までを確認しています。7 コードとも**実行と出力解析は未確認・未実装**であり、生成したファイルに `analyze.py` は入りません。各入力キーワードと実行方法の出典は [SOURCES.md](../examples/new_engines/SOURCES.md) にあります。Materials Studio・Schrödinger Suite・Winmostar には対応していません。

Without the executables and required licenses, verification is limited to input syntax and mechanical checks. **Live runs and output analysis are not yet verified or implemented** for these seven codes, so their projects do not contain `analyze.py`. Input-keyword and launch references are in [SOURCES.md](../examples/new_engines/SOURCES.md). Materials Studio, Schrödinger Suite, and Winmostar are not supported.

## NWChem の最初の対応範囲

`spec.json` の `method` に `{"code":"nwchem","theory":"dft","basis":"cc-pvdz","xc":"b3lyp"}` のように指定し、`task.type` は `single_point` にします。SCF の場合は `theory` を `scf` にし、`xc` は空にします。この値は**記法の例**であり、ADIT の推奨ではありません。`adit-gen spec.json out/` は `nwchem.nw`、投入用の `submit.sh`、`spec.json`、`README.txt`、`analyze.py` を作りますが、投入・実行しません。並列実行では環境設定の `commands.nwchem` に起動コマンドを明示してください。実行後の `adit-analyze out/` は `output.log` の `Total SCF energy` / `Total DFT energy` だけを Hartree から eV に換算して読みます。正常終了や構造・振動数は判定しません。NWChem の basis library に指定した基底関数があるかは、実行環境で確認してください。NWChem と ORCA に同じ手法名・基底名を書いても、その他の既定条件まで同一にはなりません。

In `spec.json`, set `method` to, for example, `{"code":"nwchem","theory":"dft","basis":"cc-pvdz","xc":"b3lyp"}` and `task.type` to `single_point`. For SCF, set `theory` to `scf` and leave `xc` empty. These values illustrate syntax; they are **not ADIT recommendations**. `adit-gen spec.json out/` writes `nwchem.nw`, a `submit.sh` for the user to run, `spec.json`, `README.txt` and `analyze.py`; it never submits or runs the job. For parallel runs, explicitly set the launch command in `commands.nwchem`. After execution, `adit-analyze out/` reads only `Total SCF energy` or `Total DFT energy` from `output.log` and converts hartree to eV. It does not assess termination, geometry, or frequencies. Check that the selected basis exists in your NWChem installation. Using the same method and basis names in NWChem and ORCA does not make their other default settings identical.

入力の根拠 / Input references: [geometry](https://nwchemgit.github.io/Geometry.html), [basis](https://nwchemgit.github.io/Basis.html), [charge](https://nwchemgit.github.io/Charge.html), [SCF](https://nwchemgit.github.io/Hartree-Fock-Theory-for-Molecules.html), [DFT](https://nwchemgit.github.io/Density-Functional-Theory-for-Molecules.html), [task](https://nwchemgit.github.io/TASK.html). 出力の見出し / Output label: [NWChem sample output](https://nwchemgit.github.io/Sample.html).

## OpenMM の対応範囲

`spec.json` の `method` を `{"code":"openmm","input_format":"gromacs","topology_file":"topol.top","coordinates_file":"conf.gro","include_dir":"<GROMACS の share/gromacs/top>","nonbonded_method":"PME","nonbonded_cutoff_nm":0.9,"constraints":"HBonds"}` のように書きます。`input_format` は `gromacs` (`.top` と `.gro`) か `amber` (`prmtop` と `rst7`/`inpcrd`) です。`task.type` は `single_point`、`geometry_optimization`(`LocalEnergyMinimizer`)、`molecular_dynamics` に対応し、振動解析とバンドは止めます。MD の統計集団は NVE・NVT・NPT で、熱浴は OpenMM にある langevin (`LangevinMiddleIntegrator`)・nose_hoover (`NoseHooverIntegrator`)・andersen (`AndersenThermostat`) だけです。berendsen と csvr は OpenMM に無いので止めます。**上の値は書式の例であり、ADIT の推奨ではありません。**

`adit-gen spec.json out/` は `run_openmm.py`、`openmm_settings.json`、写したトポロジーと座標を書きます。投入も実行もしません。実行する環境に openmm が要ります (確かめたのは 8.6.1 と 8.4.0.dev。ADIT 自身は openmm を import しません)。非結合の扱い・カットオフ・拘束は既定を持たず、利用者が明示します (OpenMM の既定 `NoCutoff` を黙って使いません)。熱浴の時定数は衝突頻度 (1 / 時定数) として渡します。NPT の圧力浴は `MonteCarloBarostat` で、体積を試す間隔は `barostat_interval_steps` (OpenMM の既定 25 ステップ) です。共通 Spec の圧力浴の時定数はこの間隔に機械的に換算できないため、入力には使わず README に「使っていない」と書きます。固定原子・固定軸、電荷・多重度、続きの計算 (チェックポイント) は対応しません。

実行後の `adit-analyze out/` は `md.log` (OpenMM の `StateDataReporter` の CSV) からエネルギー・温度・体積・密度を、`results.json` から最後のエネルギーと温度を読みます。軌跡の `trajectory.dcd` はバイナリなので読みません (VMD・OVITO・TRAVIS でそのまま開けます)。ADIT の RDF・MSD に軌跡が要るときは `write_xyz_trajectory` を有効にすると、テキストの `trajectory.extxyz` も書きます。実際に実行して確かめた例は `examples/openmm_spce_nvt_generated` (openmm 8.4.0.dev、SPC/E 水 216 分子の NVT。同じ入力を 8.6.1 でも実行しています)。

In `spec.json`, write `method` as shown above. `input_format` is `gromacs` (`.top` plus `.gro`) or `amber` (`prmtop` plus `rst7`/`inpcrd`). Supported `task.type` values are `single_point`, `geometry_optimization` (`LocalEnergyMinimizer`), and `molecular_dynamics`; vibrations and band structure are rejected. MD ensembles are NVE, NVT, and NPT, and the thermostats are limited to those OpenMM provides: langevin (`LangevinMiddleIntegrator`), nose_hoover (`NoseHooverIntegrator`), and andersen (`AndersenThermostat`); berendsen and csvr are rejected because OpenMM has neither. **These values illustrate syntax; they are not ADIT recommendations.**

`adit-gen spec.json out/` writes `run_openmm.py`, `openmm_settings.json`, and copies of the topology and coordinates; it never submits or runs the job. The execution host needs openmm (verified with 8.6.1 and 8.4.0.dev; ADIT itself never imports openmm). The nonbonded method, cut-off, and constraints have no defaults and must be stated by the user, so the OpenMM default (`NoCutoff`) is never applied silently. The common coupling time is passed as a collision frequency (1 / coupling time). NPT uses `MonteCarloBarostat` with a volume-move interval of `barostat_interval_steps` (the OpenMM default is 25 steps); the common barostat time has no mechanical conversion to that interval, so it is not used, and the README says so. Fixed atoms and axes, total charge and multiplicity, and restarts from checkpoints are not supported.

After the run, `adit-analyze out/` reads energies, temperature, volume, and density from `md.log` (the OpenMM `StateDataReporter` CSV) and the final energy and temperature from `results.json`. The binary `trajectory.dcd` is not parsed (open it directly in VMD, OVITO, or TRAVIS). Enable `write_xyz_trajectory` to also write a text `trajectory.extxyz`, which adit's RDF and MSD read. A verified run is in `examples/openmm_spce_nvt_generated` (openmm 8.4.0.dev, 216 SPC/E waters, NVT; the same input was also run with 8.6.1).

入力と出力の根拠 / References: [running simulations](https://docs.openmm.org/latest/userguide/application/02_running_sims.html), [app layer API](https://docs.openmm.org/latest/api-python/app.html), [library API](https://docs.openmm.org/latest/api-python/library.html). 引数は openmm 8.6.1 の定義でも確かめました / Signatures were also verified against openmm 8.6.1.

## ABINIT の対応範囲

`method` を `{"code":"abinit","pseudos":{"Si":"/…/14si.pspnc"},"ecut_ha":8,"tolerance":"toldfe","tolerance_value":1e-8}` のように書き、`kpoints` も指定します。対応するのは **3 次元周期系の SCF 一点計算と構造最適化 (セルは動かさない)** だけで、分子 (非周期)、MD、振動 (anaddb)、バンド、帯電系、スピン分極は理由を示して止めます。構造は `acell 1 1 1` と Bohr 単位の `rprim`・`xred` で書き、固定原子は `iatfix`、軸ごとの固定は `iatfixx/y/z` になります (番号は ABINIT の約束どおり 1 始まり)。ABINIT は `iatfixx/y/z` を格子ベクトル方向 (還元座標) で解釈するので、Cartesian の軸固定は VASP の Selective dynamics と同じ変換で格子方向へ直し、直せないセル (固定する方向が格子ベクトルと合わないとき) は止めます。**擬ポテンシャルは利用者のファイルを写すだけで、ADIT は同梱も検証もしません。**汎関数も選びません (`ixc` を書かなければ、擬ポテンシャルに書かれた汎関数を ABINIT が使います)。`ecut`・収束の種類と閾値・占有の仕方は既定を持ちません。収束の種類は `toldfe`・`toldff`・`tolrff`・`tolvrs`・`tolwfr` の 5 つから 1 つを選びます (ABINIT は 1 つしか受け取りません)。構造最適化では ABINIT が `toldfe` を受け付けないので、`toldff`・`tolrff`・`tolvrs` (または `tolwfr`) から選ばせます (実際に実行して確かめた制限。ABINIT の文書も構造最適化では `toldfe` を避けるよう書いています)。画面に無い変数は `extra` にそのまま書けます (この生成器が書く変数と重なると止めます)。

解析は `output.log` の段階ごとの `Total energy (etotal) [Ha]`(無ければ `input.abo` の最後の `etotal`)と、`input.abo` 末尾の変数の表から組み立てた最後の構造を読みます。**ASE 3.29.0 の `abinit-out` の読みは、この出力の座標と合わなかったので使っていません。**実際に実行した例は `examples/abinit_si_generated` (ABINIT 10.0.3、ダイヤモンド構造 Si)。

Write `method` as shown above and also give `kpoints`. The supported range is **single-point SCF and geometry optimization (fixed cell) for fully periodic systems**; molecules, MD, phonons (anaddb), band structures, charged cells and spin polarization are rejected with a reason. The structure is written as `acell 1 1 1` with `rprim` and `xred` in bohr; fixed atoms become `iatfix` and per-axis constraints `iatfixx/y/z` (1-based, as ABINIT expects). ABINIT interprets `iatfixx/y/z` along the lattice vectors (reduced coordinates), so Cartesian axis constraints are converted to lattice directions the same way as VASP's Selective dynamics, and a cell where the fixed directions do not match lattice vectors is rejected. **Pseudopotentials are copied from the files you name; ADIT neither ships nor validates them,** and it selects no functional (without `ixc`, ABINIT uses the functional recorded in the pseudopotential). `ecut`, the convergence criterion and its threshold, and the occupation scheme have no defaults. The criterion is one of `toldfe`, `toldff`, `tolrff`, `tolvrs` and `tolwfr` (ABINIT accepts exactly one). For a geometry optimization ABINIT rejects `toldfe`, so ADIT requires `toldff`, `tolrff` or `tolvrs` (or `tolwfr`) (a limit found by running it; the ABINIT documentation also says to avoid `toldfe` when the geometry is optimized). Variables without a field can be passed through `extra` (clashes with what the generator writes are rejected).

Analysis reads the per-step `Total energy (etotal) [Ha]` from `output.log` (or the final `etotal` in `input.abo`) and rebuilds the final structure from the variable table at the end of `input.abo`. **ASE 3.29.0's `abinit-out` reader disagreed with those coordinates, so it is not used.** A verified run is in `examples/abinit_si_generated` (ABINIT 10.0.3, diamond Si).

## Psi4 の対応範囲

`method` を `{"code":"psi4","method":"scf","basis":"cc-pVDZ"}` のように書きます。`task.type` は `single_point`(`energy`)、`geometry_optimization`(`optimize`)、`vibrations`(`frequencies`)に対応し、周期系・MD・バンド・固定原子は止めます。入力は psithon (Psi4 の入力言語) で、分子の並びは ADIT の構造そのまま (`no_reorient`・`no_com`・`symmetry c1`)。多重度が 1 でないときは開殻の参照関数 (`uhf` / `rohf` / `uks`) を明示させます (Psi4 の既定 RHF は閉殻だけ)。電子数と多重度が合わない組み合わせは、ほかのコードと同じ共通の検査で止まります。入力の末尾の数行が `results.json` と `final.xyz` を書き、解析はそれと `output.log` の「Optimization Summary」の表 (段階ごとのエネルギー) を読みます。構造最適化の収束は、Psi4 の `optimize()` が収束しないと例外で止まることを使って判定します (`results.json` に `converged` があれば収束)。実際に実行した例は `examples/psi4_h2o_generated` (Psi4 1.11、水の RHF/cc-pVDZ 最適化)。

Write `method` as shown above. Supported `task.type` values are `single_point` (`energy`), `geometry_optimization` (`optimize`) and `vibrations` (`frequencies`); periodic systems, MD, band structures and fixed atoms are rejected. The input is psithon with the atoms in ADIT's order (`no_reorient`, `no_com`, `symmetry c1`). A multiplicity other than 1 requires an explicit open-shell reference (`uhf`, `rohf`, `uks`), because the Psi4 default RHF is closed-shell only; an electron-count/multiplicity mismatch is caught by the same shared check as the other codes. The last lines of the input write `results.json` and `final.xyz`; analysis reads those plus the 'Optimization Summary' table in `output.log`. Convergence of an optimization is decided from the fact that Psi4's `optimize()` raises instead of returning when it fails (`converged` in `results.json`). A verified run is in `examples/psi4_h2o_generated` (Psi4 1.11, RHF/cc-pVDZ water).

## PLUMED の対応範囲

PLUMED は単体の計算コードではないので、`method` ではなく Spec の `plumed` に書きます: `{"plumed":{"lines":"d1: DISTANCE ATOMS=1,2\nPRINT ARG=d1 FILE=COLVAR STRIDE=10"}}`、または `{"plumed":{"input_file":"/…/plumed.dat"}}`(どちらか一方)。**ADIT は集合変数もバイアスも作りません。**利用者の入力を `plumed.dat` として生成先へ写し、MD の入力に読み込みの 1 行を足すだけです。つなげるのは LAMMPS・GROMACS・OpenMM の MD だけで、ほかのコードや MD 以外の計算では止めます。

つなぎ方はどれも実際に実行して確かめました: LAMMPS は `fix adit_plumed all plumed plumedfile plumed.dat outfile plumed.log`、GROMACS は `gmx mdrun … -plumed plumed.dat`、OpenMM は `openmm-plumed` の `PlumedForce`。**GROMACS 2026.3 では、`-plumed` があっても実行時に環境変数 `PLUMED_KERNEL` が `libplumedKernel.so` を指していないと止まりました。**この事実は README.txt に入れています。`plumed` が PATH にあるときは、生成の前に `plumed driver --parse-only` で入力を読ませ、読めない入力を理由付きで止めます (中身の物理には口を出しません)。`lmp -h` / `gmx mdrun -h` でその実行ファイルが PLUMED を使えるかも見ますが、確かめられないときは何も表示せず、止めません。PLUMED の既定の単位は nm・kJ/mol・ps で、計算コードの単位系とは別です (実際に実行した例で Cu の 2.556 Å が `0.2556` と書かれることを確かめました)。解析は `#! FIELDS` で始まる COLVAR 形式のファイルを探し、列の名前・行数・最初と最後と平均を要約に入れます (列の意味は解釈しません)。実際に実行した例は `examples/plumed_lammps_cu_generated`。

PLUMED is not a standalone engine, so it lives in the spec's `plumed` field rather than in `method`: either `lines` (the PLUMED input as text) or `input_file` (a file to copy), never both. **ADIT invents no collective variables or bias terms**; it copies your input to `plumed.dat` and adds one line to the MD input. Only LAMMPS, GROMACS and OpenMM MD runs are supported; anything else is rejected.

Every integration was verified by running it: `fix adit_plumed all plumed plumedfile plumed.dat outfile plumed.log` for LAMMPS, `gmx mdrun … -plumed plumed.dat` for GROMACS, and `PlumedForce` from `openmm-plumed` for OpenMM. **With GROMACS 2026.3, `-plumed` alone was not enough: the run stopped unless `PLUMED_KERNEL` pointed at `libplumedKernel.so`** — the README says so. When `plumed` is on PATH, ADIT has it parse the input (`plumed driver --parse-only`) before generating and stops with the reason if it cannot; it never comments on the physics. It also checks `lmp -h` / `gmx mdrun -h` for PLUMED support, and stays silent when that cannot be determined. PLUMED's default units are nm, kJ/mol and ps regardless of the engine (the example shows Cu's 2.556 Å written as `0.2556`). Analysis finds files in COLVAR format (`#! FIELDS`) and reports the column names, the row count, and the first/last/mean values without interpreting them. A verified run is in `examples/plumed_lammps_cu_generated`.


## 追加の入口 (ORCA GOAT・DOCKER ほか)

### ORCA GOAT

`method.code="orca"`, `task.type="single_point"`, `method.goat=true` を `spec.json` に指定し、`adit-gen spec.json out/` で生成します。`method.method="XTB"` と `method.basis=""` は [ORCA 6.1 GOAT チュートリアル](https://www.faccts.de/docs/orca/6.1/tutorials/prop/goat.html)と同じ入力の例です。計算後は `orca.globalminimum.xyz` と `orca.finalensemble.xyz`、エネルギーと重みは `output.log` を確認してください。ORCA 本体は同梱しておらず、ORCA での実行は確認していません。

Set `method.code="orca"`, `task.type="single_point"`, and `method.goat=true` in `spec.json`, then run `adit-gen spec.json out/`. For the documented ORCA 6.1 tutorial form, set `method.method="XTB"` and `method.basis=""`. After the run, inspect `orca.globalminimum.xyz`, `orca.finalensemble.xyz`, and the energies and weights in `output.log`. ORCA is not bundled, and running this input with ORCA has not been verified.

### ORCA DOCKER / Host–guest docking

ホストは `spec.json` の構造、ゲストは `method.docker_guest_file` の XYZ で指定します。`task.type="single_point"`、`method.method="XTB"`、`method.basis=""` にし、XYZ の各フレームの 2 行目に電荷と多重度を整数で書きます (例 `0 1`)。書かれていないフレームを ORCA の既定 (0, 1) として扱う場合だけ、`method.docker_assume_neutral_singlet=true` を明示します。ADIT はゲストを `guest.xyz` として写し、`%DOCKER GUEST "guest.xyz" END` を生成します。[ORCA 6.1 マニュアル](https://www.faccts.de/docs/orca/6.1/manual/contents/structurereactivity/docker.html)に従い、DOCKER は XTB/GFN-xTB/GFN-FF に限ります。結果は `orca.docker.xyz` と `output.log` を確認してください。

The host is the structure in `spec.json`; set `method.docker_guest_file` to the guest XYZ. Use `task.type="single_point"`, `method.method="XTB"`, and `method.basis=""`. The second line of each XYZ frame must contain integer charge and multiplicity, for example `0 1`. Only if you intend ORCA's default (0, 1) for frames without these values, explicitly set `method.docker_assume_neutral_singlet=true`. ADIT copies the guest as `guest.xyz` and writes `%DOCKER GUEST "guest.xyz" END`. ORCA 6.1 restricts DOCKER to XTB/GFN-xTB/GFN-FF. Inspect `orca.docker.xyz` and `output.log` after running.

### DCDFTBMD 2.0

`method.code="dcdftbmd"` に加え、`method.scc`、`method.divide_and_conquer`、構造中の各元素の `method.highest_angular_momentum` (s=1、p=2、d=3、f=4)、**順序つき**の各元素対の `method.sk_files` (`.spl`)を明示します。例: O と H なら `O-O`、`O-H`、`H-O`、`H-H` が必要です。`adit-gen spec.json out/` は `dftb.inp`、`.spl` の複製、`submit.sh`、`spec.json`、`README.txt` を生成します。DFTB+ の `.skf` は変換しません。対応するのは一点計算、最急降下法/FIRE の構造最適化、NVE と Berendsen NVT です。NVT の温度と時定数は共通 MD 条件から単位だけ換算します。NPT、開殻、軸固定、その他の熱浴は生成前に止めます。実行ファイルと `.spl` は利用者が入手してください。[公式 2.0 マニュアル](https://www.chem.waseda.ac.jp/dcdftbmd/document/DCDFTBMD_2.0_en.pdf)。

For `method.code="dcdftbmd"`, explicitly set `method.scc`, `method.divide_and_conquer`, `method.highest_angular_momentum` for each element (s=1, p=2, d=3, f=4), and `method.sk_files` for every **ordered** element pair. For O and H, specify `O-O`, `O-H`, `H-O`, and `H-H` `.spl` files. `adit-gen spec.json out/` creates `dftb.inp`, copies the `.spl` files, and creates `submit.sh`, `spec.json` and `README.txt`; it does not convert DFTB+ `.skf` files. Supported tasks are single point, SteepestDescent/FIRE optimization, NVE, and Berendsen NVT. Temperature and coupling time are converted only in units. NPT, open-shell systems, axis constraints, and other thermostats stop before generation. Obtain the executable and `.spl` files yourself. See the [official 2.0 manual](https://www.chem.waseda.ac.jp/dcdftbmd/document/DCDFTBMD_2.0_en.pdf).

### DOCK6

`adit-convert dock6 prepared/dock.in ready/` は、**利用者が完成させた** `dock.in` と、そこから参照される配位子、球、格子、既知の定義ファイルをまとめます。必要に応じて `--asset relative/path` を繰り返して追加してください。パスは `dock.in` と同じ作業ディレクトリからの相対パスにし、出力先は入力ディレクトリの外の空の場所にします。元の `dock.in` は変更せず、複製の SHA-256 と `run.sh` を書きます。受容体・配位子の準備、電荷、原子型、スコア条件を作る機能ではありません。[DOCK 6.13 マニュアル](https://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm)の入力を完成させてから使ってください。

`adit-convert dock6 prepared/dock.in ready/` packages a **completed, user-authored** `dock.in`, its ligand, spheres, grid files, and recognized definition files. Repeat `--asset relative/path` for other dependencies. Paths must be relative to the `dock.in` working directory, and the empty output directory must be outside the input directory. The input deck is unchanged; ADIT writes SHA-256 hashes and `run.sh`. This does not prepare the receptor or ligand, assign charges or atom types, or choose scoring settings. Complete the input using the [DOCK 6.13 manual](https://dock.compbio.ucsf.edu/DOCK_6/dock6_manual.htm) first.

### Open Babel

`adit-convert openbabel source.sdf target.mol2 --input-format sdf --output-format mol2` は、PATH 上の `obabel` を明示的に使用します。入力と出力は同じファイルにできず、既存の出力は `--overwrite` を付けない限り変更しません。変換結果は一時ファイルに書いてから移動します。失敗時は以前の出力を保ち、途中ファイルの場所を表示します。結合次数・電荷・セル・座標が目的どおりか確認してください。[Open Babel の CLI 文書](https://openbabel.org/docs/Command-line_tools/babel.html)。

`adit-convert openbabel source.sdf target.mol2 --input-format sdf --output-format mol2` explicitly invokes `obabel` from PATH. Input and output cannot be the same file; an existing output is untouched unless `--overwrite` is given. ADIT stages the result before moving it into place. On failure it preserves the previous output and reports the partial file. Check bond orders, charges, cell, and coordinates for your use case. See the [Open Babel CLI documentation](https://openbabel.org/docs/Command-line_tools/babel.html).

## 計算コード間で条件を突き合わせる

`adit-convert calculation` は、変換の記録 `conversion.json` に `field_audit` を追加します。
これは値の**出典と適用状況**の記録であり、異なる計算コードが同じ物理モデル・
数値解を与えるという認証ではありません。

`adit-convert calculation` adds `field_audit` to `conversion.json`. It records
the **origin and input applicability** of settings, not a certification that
different engines implement the same physical model or produce the same result.

| status | 日本語 | English |
| --- | --- | --- |
| `preserved` | 元の Spec の値を保持 | Value preserved from the source spec |
| `changed_by_request` | 明示指定で雛形の値に変更 | Changed to a target-template value at the user's request |
| `omitted_by_request` | 明示指定で除外 | Omitted at the user's request |
| `retained_in_spec_not_input` | Spec に残るが生成入力には適用されない | Retained in the spec but not applied to generated input |
| `not_used_by_task` | この計算種類では使わない | Not used by this task type |
| `from_target_template` | 変換せず雛形から取得 | Taken from the target template without translation |
| `not_translated` | コード固有の方法を対応付けない | Engine-specific method settings are not mapped |

座標・原子記号・セルなどの大きな配列の値は `field_audit` に複写しません。
正確な値は生成先の `spec.json` と構造ファイルで確認できます。
入力に書かれていない暗黙の既定値、力場と擬ポテンシャルの物理的な同等性、
熱浴アルゴリズムの同等性は、この記録では保証しません。

Large arrays such as coordinates, symbols, and cell vectors are not copied into
`field_audit`; inspect the generated `spec.json` and structure files for their
exact values. The audit does not guarantee implicit engine defaults, physical
equivalence of force fields or pseudopotentials, or equivalence of thermostat
algorithms.

## 計算結果の機械的点検

```text
adit-analyze run_lammps --audit-with run_vasp --audit-kind msd
adit-analyze run_lammps --audit-with run_vasp --audit-kind energy
```

このコマンドは点検結果を標準出力へ示すだけで、MSD やエネルギーの図は作りません。`--msd` など通常の解析オプションと同時には指定できません。解析図が必要なら各計算を別々に解析してください。

This command reports the audit on standard output only; it does not produce MSD or energy plots. It cannot be combined with ordinary analysis options such as `--msd`. Analyze each run separately when you need figures.

MSD の点検では、元素記号の列・周期境界・セル・軌跡・フレーム間隔を確認します。同じ元素どうしの原子 ID は記号だけでは区別できないため、初期座標も周期境界を考慮して照合します。異なる場合は「同じ初期構造からの比較」の点検を通しません。異なる初期構造を用いる統計比較が妥当かどうかは判定しません。
エネルギーの点検では、同じ組成・周期境界・セルで、全系の eV として読めるかを確認します。軌跡の有無やフレーム数・原子順はエネルギー用の停止条件ではありません。出力は読み取り専用で、
失敗した点検は終了状態 1 を返します。異なる計算法のエネルギーや MSD を
**科学的に直接比較してよいか**は判定しません。
異なるコードのエネルギーでは、運動・ポテンシャル・電子の寄与など、出力が含む内訳を照合していないため、直接比較の点検は通しません。

For MSD, the command checks the element-symbol sequence, periodic boundaries,
cell geometry, trajectories, and frame spacing. Symbols cannot distinguish
same-element atom identities, so initial positions are matched with periodic
boundaries accounted for. A mismatch fails the check for runs starting from the
same structure; this does not judge statistical comparisons from different
initial structures. For energy, it checks matching composition, periodic boundaries, and cells, and whether energies are available in eV for the whole simulated system; missing
trajectories, frame-count differences, and atom order do not block this mode.
It does not write files and exits with status 1 if a mechanical check fails.
It does **not** decide whether the two methods' energies or MSDs are
scientifically comparable.
Direct energy comparison across different engines fails the check because the
contributions included in their reported energies have not been matched.

## 出力ログの警告と失敗の印の走査

解析の要約 (`adit-analyze`、解析タブ、ウェブ版の解析ページ、`analysis/summary.json` の `tables.diagnostics`) には、出力ログにある警告の件数と、失敗の原因の候補が載ります。
探す文字列は各コードが実際に書くものだけで、対応するコードは DFTB+、VASP、Quantum ESPRESSO (pw.x)、xtb、ORCA、CP2K、LAMMPS、GROMACS です。
ほかのコードでは、ジョブスケジューラ (PBS / Slurm) と MPI と OS の印 (制限時間、メモリ不足、MPI の異常終了) だけを探します。
見つけた行を写して分類するだけで、対処は判断しません。各コードの文書が挙げる対処は、生成した `README.txt` の末尾に出典付きで載せています。詳しくは [解析の詳細](ANALYSIS.md)。
