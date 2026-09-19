
from __future__ import annotations

from adit.config import env_var
import os

from PySide6.QtWidgets import (QAbstractButton, QComboBox, QGroupBox, QLabel, QLineEdit, QMenu, QMenuBar, QPlainTextEdit, QTabBar, QTabWidget, QToolBar,
                               QToolButton, QWidget)

LANGUAGE = "ja"

_EN: dict[str, str] = {
    "構造": "Structure", "プリセット": "Preset", "ファイル": "File", "バルク": "Bulk", "スラブ": "Slab",
    "3D 構造を生成": "Build 3D structure", "参照…": "Browse…", "全電荷": "Total charge", "スピン多重度": "Spin multiplicity",
    "固定原子": "Fixed atoms", "立方晶セル": "Cubic cell", "層": "layers", "真空層 [Å]": "Vacuum [Å]",
    "周期セルに入れる": "Put in a periodic cell", "分子を立方体の周期セルに置きます。VASP と pw.x では周期セルが必須です": "Place the molecule in a cubic periodic cell. VASP and pw.x require a periodic cell", "一辺 [Å]": "edge [Å]",
    "ase gui で構造を確認": "View structure with ase gui", "例 CCO": "e.g. CCO",
    "1 始まりの番号 (例 1-4,7)。軸を指定するなら 7:xy (x, y を固定)。空欄なら全原子を動かします":
        "Fixed atoms (1-based, e.g. 1-4,7). Axes only: 7:xy (freeze x,y). Empty = all free",
    "SMILES (RDKit が未導入のため使えません)": "SMILES (RDKit not installed)",
    "0 なら元素の既定値 (ASE の参照状態)": "0 = element default (ASE reference state)",
    "計算手法": "Method", "計算コード": "Code", "Slater-Koster パラメータ": "Slater-Koster set", "SCC (自己無撞着電荷)": "SCC (self-consistent charges)",
    "分散力補正": "Dispersion", "DFT-D3 係数 (BJ)": "DFT-D3 parameters (BJ)", "電子温度 [K]": "Electronic temperature [K]",
    "POTCAR のセット": "POTCAR set", "元素ごとの POTCAR 名": "POTCAR name per element", "実行ファイルの種類": "Binary", "IBRION (構造最適化の方法)": "IBRION (optimization method)",
    "ENCUT [eV] (0 = 指定しない)": "ENCUT [eV] (0 = omit)", "追加の INCAR 設定": "Additional INCAR settings",
    "画面にない INCAR のキーを 1 行に 1 つ、KEY = value の形で (例 NCORE = 4)": "INCAR keys not in the GUI, one per line: KEY = value (e.g. NCORE = 4)",
    "空欄なら指定しない。原子ごとに空白区切り (例 1 1 -1)": "Empty = omit. One value per atom (e.g. 1 1 -1)",
    "計算手法 (--gfn)": "Method (--gfn)", "精度 (--acc)": "Accuracy (--acc)", "電子温度 [K] (--etemp)": "Electronic temperature [K] (--etemp)",
    "SCC の最大反復回数": "Max SCC iterations", "最適化の収束レベル (--opt)": "Optimization level (--opt)", "GFN-FF (力場)": "GFN-FF (force field)",
    "パラメータファイルは不要です (計算手法に内蔵)。分子系のみ。電荷とスピン多重度は構造パネルの値をコマンドライン引数で渡します":
        "No parameter files are needed (built into the method). Molecules only. Charge and multiplicity are passed from the structure panel as command-line arguments",
    "擬ポテンシャルのセット": "Pseudopotential set", "元素ごとの UPF ファイル": "UPF file per element", "ecutwfc [Ry]": "ecutwfc [Ry]",
    "ecutrho [Ry] (0 = 指定しない)": "ecutrho [Ry] (0 = omit)", "追加の変数 (名前空間.変数 = 値)": "Additional variables (namelist.variable = value)",
    "画面にない変数を 1 行に 1 つ、名前空間.変数 = 値 の形で (例 system.nbnd = 20)": "Variables not in the GUI, one per line: namelist.var = value (e.g. system.nbnd = 20)",
    "空欄なら UPF の汎関数を使用 (例 PBE, PBEsol)": "Empty = functional from the UPF (e.g. PBE, PBEsol)",
    "k 点 (周期系)": "k-points (periodic)", "サンプリング方法": "Sampling", "Γ 点のみ": "Γ only", "メッシュを指定": "Mesh", "密度から自動": "From density",
    "メッシュ (n1 n2 n3)": "Mesh (n1 n2 n3)", "シフト": "shift", "k 点密度 [点/Å⁻¹]": "k-point density [points per Å⁻¹]", "分子系 (非周期) では使わない": "Not used for molecules",
    "計算の種類": "Task", "種類": "Type", "一点計算": "Single point", "構造最適化": "Geometry optimization", "最適化アルゴリズム": "Optimizer",
    "分子動力学": "Molecular dynamics", "振動解析": "Vibrations", "アンサンブル": "Ensemble", "熱浴": "Thermostat", "温度 [K]": "Temperature [K]",
    "時間刻み [fs]": "Time step [fs]", "MD ステップ数": "MD steps", "軌跡の出力間隔 [ステップ]": "Trajectory interval [steps]",
    "熱浴の緩和時間 [fs]": "Thermostat relaxation time [fs]", "圧力 [bar] (NPT)": "Pressure [bar] (NPT)", "圧力浴の緩和時間 [fs] (NPT)": "Barostat relaxation time [fs] (NPT)",
    "CSVR (速度再スケール)": "CSVR (velocity rescaling)", "バンド計算": "Band structure", "k 点の経路 (バンド)": "k-path (bands)",
    "経路上の k 点数": "k-points along the path", "空のバンド数 (pw.x)": "Empty bands (pw.x)",
    "空欄なら格子の標準経路 (例 GXWKGLUWLK,UX)": "Empty = standard path for the lattice (e.g. GXWKGLUWLK,UX)",
    "力の収束判定 [eV/Å]": "Force convergence [eV/Å]", "セルの緩和 (周期系)": "Cell relaxation (periodic)", "固定": "fixed", "形と体積": "shape and volume", "体積だけ": "volume only",
    "実行環境とリソース": "Execution target and resources", "プロファイル": "Profile", "ノード数": "Nodes", "ノードあたりのコア数": "Cores per node",
    "ノードあたりの MPI プロセス数": "MPI processes per node", "OpenMP スレッド数": "OpenMP threads", "制限時間 (HH:MM:SS)": "Time limit (HH:MM:SS)",
    "ジョブ名": "Job name", "出力ディレクトリ": "Output directory", "プロファイルが無い (cluster.toml を確認)": "No such profile (check cluster.toml)",
    "計算設定 (spec.json) を開く…": "Open calculation settings (spec.json)…", "環境設定を再読み込み": "Reload preferences", "実行ボタンを表示": "Show the Run button",
    "生成": "Generate", "操作": "Actions", "生成できます": "Ready to generate",
    "先に「生成」を押してください": "Press “Generate” first", "実行中": "Running",
    "クラスタ用 (PBS / Slurm) の入力は、このアプリからは実行しません": "Cluster (PBS / Slurm) input is not run from this app",
    "Windows では実行できません (生成した入力を Linux のサーバーに転送して使います)": "Cannot run on Windows (transfer the generated input to a Linux server)",
    "環境設定ファイルを作りました": "Settings file created", "表示言語を切り替えます (環境設定の language に保存。再起動後に反映)": "Switch the display language (saved as language in the preferences; effective after restart)", "挿入": "Insert", "表示": "View", "設定": "Preferences",
    "ヘルプ": "Help", "計算設定 (spec.json) を保存…": "Save calculation settings (spec.json)…", "終了": "Quit", "構造ファイルを読み込む…": "Load a structure file…",
    "溶液に成分を追加": "Add a component to the solution", "テーマ": "Theme", "システムに従う": "Follow the system", "ライト": "Light", "ダーク": "Dark",
    "ADIT について": "About ADIT", "README を開く": "Open the README", "言語 (再起動後に反映)": "Language (after restart)", "テーマ (再起動後に反映)": "Theme (after restart)",
    "ウィンドウの枠 (再起動後に反映)": "Window frame (after restart)", "自動 (Linux では ADIT が描く)": "Automatic (drawn by ADIT on Linux)",
    "ADIT が描く (ボタンにカーソルで色が点く)": "Drawn by ADIT (buttons light up on hover)", "OS に任せる": "Use the system frame",
    "フォルダを選ぶ…": "Choose folder…", "指定": "Reference",
    "Slater-Koster パラメータのセット (mio-1-1 など) を入れたフォルダを選びます。選んだ場所は環境設定 (sk_root) に保存されます":
        "Choose the folder holding Slater-Koster parameter sets (such as mio-1-1). Saved as sk_root in the preferences",
    "UPF のセット (SSSP など) のフォルダを並べた場所を選びます。選んだ場所は環境設定 (pseudo_root) に保存されます":
        "Choose the folder holding UPF sets (such as SSSP). Saved as pseudo_root in the preferences",
    "閉じる": "Close", "最小化": "Minimize", "最大化 / 元の大きさ": "Maximize / restore",
    "右側のタブ": "Right tabs", "溶液・混合物": "Solution / mixture", "成分を追加": "Add component", "削除": "Remove",
    "成分の種類": "Kind", "指定 (名前 / SMILES / ファイル)": "Reference (name / SMILES / path)", "個数": "Count", "電荷": "Charge", "表示名": "Label",
    "密度から自動 [g/cm³]": "From density [g/cm³]", "一辺を指定 [Å]": "Edge [Å]", "セル": "Box", "分子間の最短距離 [Å]": "Min. distance [Å]",
    "Draw: 分子を描いて SMILES にします (RDKit が要ります)": "Draw: sketch the molecule and turn it into SMILES (needs RDKit)", "構造の作り方": "Source", "共通設定": "Common settings", "環境設定…": "Preferences…",
    "元に戻す (1 つ前の設定)": "Undo (previous settings)", "1 つ前の設定に戻す": "Back to the previous settings",
    "探す": "Find", "やり直す": "Redo", "環境設定ファイル (cluster.toml) を編集": "Open and edit the settings file (cluster.toml)",
    "青字は必須項目です。ラベルにカーソルを合わせると説明が表示されます": "Blue labels are required. Hover a label for an explanation",
    "左ドラッグで回転、ホイールで拡大縮小、右ドラッグで移動、ダブルクリックでリセット": "Left-drag to rotate, wheel to zoom, right-drag to pan, double-click to reset",
    "分子系 (非周期) では使いません": "Not used for molecules", "生成ファイル": "Generated files", "実行": "Run",
    "この PC で実行": "Run on this PC", "視点": "Viewpoint", "繰り返し数": "Repeat", "ASE GUI で開く": "Open in ASE GUI", "斜め": "oblique", "z 軸から": "along z",
    "y 軸から": "along y", "x 軸から": "along x", "まだ実行していません": "Not run yet", "(元素の既定値)": "(element default)",
    "上書きの確認": "Confirm overwrite", "生成できません": "Cannot generate", "生成しました": "Generated", "読めません": "Cannot read", "設定を読めません": "Cannot read the settings",
    "生成ファイルのプレビュー": "Preview of generated files", "解析": "Analysis",
    "計算条件": "Calculation settings", "メインの面": "Main pane", "計算結果のディレクトリ": "Run directory", "エネルギーの推移": "Energy over time", "温度の推移 (MD)": "Temperature over time (MD)",
    "結合長 (最終構造)": "Bond lengths (final structure)", "動径分布関数 (RDF)": "Radial distribution function (RDF)", "r の最大値 [Å]": "Maximum r [Å]",
    "平均二乗変位 (MSD) と拡散係数": "MSD and diffusion coefficient", "元素": "element", "状態密度 (DOS)": "Density of states (DOS)", "ガウス幅 [eV]": "Gaussian width [eV]",
    "振動数とスペクトル": "Frequencies and spectrum", "平衡化として捨てるフレーム数": "Frames to skip (equilibration)", "解析を実行": "Run analysis",
    "計算結果のディレクトリ (生成すると自動で入ります)": "Run directory (filled in automatically after generation)", "(全原子)": "(all atoms)",
    "計算手法 (! 行)": "Method (! line)", "基底関数": "Basis set", "SCF の収束判定": "SCF convergence", "SCF の最大反復回数": "Max SCF iterations",
    "%maxcore [MB] (0 = 指定しない)": "%maxcore [MB] (0 = omit)", "追加のキーワード (! 行)": "Additional keywords (! line)", "追加の %ブロック": "Additional %blocks",
    "! 行に追加するキーワード (例 D3BJ RIJCOSX)": "Keywords added to the ! line (e.g. D3BJ RIJCOSX)",
    "%ブロックをそのまま追加 (例 %basis newGTO ... end)": "% blocks appended verbatim (e.g. %basis newGTO ... end)",
    "分子系のみ。電荷とスピン多重度は構造パネルの値を使います。並列数は実行環境の MPI プロセス数 (%pal nprocs)。ORCA 本体は登録して入手し、実行パスは環境設定に書きます":
        "Molecules only. Charge and multiplicity from the structure panel. Parallelism = MPI processes of the target (%pal nprocs). Obtain ORCA by registration; set its path in the settings",
    "conv_thr [Ry]": "conv_thr [Ry]", "mixing_beta": "mixing_beta", "保存しました。表示言語は次回の起動から: ": "Saved. The language applies from the next start: ",
    "開く": "Open", "保存": "Save", "再読み込み": "Reload", "元に戻す": "Undo", "構造ファイル": "Structure file",
    "計算設定": "Calculation settings", "アプリ": "Application", "ボタン": "Buttons", "入力": "Input", "この PC": "This PC",
    "環境設定": "Preferences", "自動": "Automatic", "ADIT が描く": "Drawn by ADIT", "OS の枠": "System frame",
    "ASE に収録された分子 (G2 集など) から選びます": "Choose a molecule from ASE's built-in set (G2 and others)",
    "例 [Na+]、CCO": "e.g. [Na+], CCO", "構造ファイルのパス": "Path to a structure file",
    "SMILES を入力して Enter。「Draw」で描くこともできます": "Type a SMILES string and press Enter, or draw the molecule with “Draw”",
    "この行の分子を描いて SMILES にします (RDKit が要ります)": "Draw this row's molecule and turn it into SMILES (requires RDKit)",
    "プリセットは一覧から選び、SMILES は入力するか「Draw」で描き、ファイルは「参照…」で選びます":
        "Pick a preset from the list, type a SMILES string or draw it with “Draw”, or choose a file with “Browse…”",
    "選んだ行 (入力欄を触った行) を消します": "Remove the selected row (the last row you clicked or typed in)",
    "分子を描いて、選んだ行 (無ければ新しい行) に SMILES で入れます (RDKit が要ります)":
        "Draw a molecule and put its SMILES in the selected row, or in a new row if none is selected (requires RDKit)",
    "2 次元材料・ナノチューブ": "2D material / nanotube", "ナノ粒子": "Nanoparticle", "ポリマー": "Polymer",
    "グラフェン型 (C₂、BN)": "Graphene-type (C₂, BN)", "MX₂ 型 (MoS₂ など)": "MX₂ sheet (MoS₂ etc.)", "ナノリボン": "Nanoribbon", "ナノチューブ": "Nanotube",
    "正二十面体": "Icosahedron", "十面体": "Decahedron", "八面体": "Octahedron", "Wulff 形 (面のエネルギーから)": "Wulff shape (from surface energies)",
    "組み立て手順": "Build steps", "超格子の取り方": "Supercell type", "繰り返し (a × b × c)": "Repeats (a × b × c)", "変換行列": "Matrix",
    "ミラー指数 (h k l)": "Miller indices (h k l)", "層の数 (周期)": "Layers (repeat units)", "真空 (片側) [Å]": "Vacuum per side [Å]",
    "終端 (最上面)": "Termination (top plane)", "軸": "Axis", "周期の像との隙間 [Å]": "Gap to periodic image [Å]", "余白 [Å]": "Padding [Å]",
    "倍数の上限": "Max. cell multiple", "元素 (条件)": "Elements", "z の範囲 [Å]": "z range [Å]", "選ぶ数": "Number to pick",
    "置き換え先の元素": "Replace with", "分子": "Molecule", "置き場所": "Position", "吸着サイトの名前": "Site", "原子の番号 (1 始まり)": "Atom number (from 1)",
    "xy 座標 [Å]": "x, y [Å]", "高さ [Å]": "Height [Å]", "下に向ける原子 (1 始まり)": "Atom facing down (from 1)", "厚みの決め方": "Thickness from",
    "隙間 [Å]": "Gap [Å]", "真空 [Å] (0 なら両側が溶液)": "Vacuum [Å] (0 = solution on both sides)", "濃度から個数": "Count from concentration",
    "溶質の周りの余白 [Å]": "Padding around solute [Å]", "密度 [g/cm³]": "Density [g/cm³]", "下から数えた層": "Bottom layers",
    "汎関数": "Functional", "基底関数のファイル": "Basis-set file", "擬ポテンシャルのファイル": "Pseudopotential file",
    "元素ごとの基底と擬ポテンシャル": "Basis and pseudopotential per element", "分散補正": "Dispersion correction",
    "カットオフ [Ry]": "Cutoff [Ry]", "相対カットオフ [Ry]": "Relative cutoff [Ry]", "SCF の収束の閾値 (EPS_SCF)": "SCF threshold (EPS_SCF)",
    "SCF の反復の上限 (MAX_SCF)": "Max SCF iterations (MAX_SCF)", "スピン分極": "Spin polarization", "ポアソン方程式の解き方": "Poisson solver",
    "分子の箱の一辺 [Å]": "Box edge for a molecule [Å]", "追加の行 (節ごと)": "Extra lines (per section)",
    "(3 方向とも周期なら不要)": "(not needed if periodic in all three directions)",
    "例 PBE、BLYP、PBE0 (CP2K の XC_FUNCTIONAL の名前)": "e.g. PBE, BLYP, PBE0 (a CP2K XC_FUNCTIONAL name)",
    "[セクションのパス] の行のあとに、その節の末尾に足す行 (例\n[FORCE_EVAL/DFT/SCF]\nSCF_GUESS ATOMIC)":
        "Lines after a [section path] line are appended to that section (e.g.\n[FORCE_EVAL/DFT/SCF]\nSCF_GUESS ATOMIC)",
    "基底関数と擬ポテンシャルは CP2K の data ディレクトリのファイルから読み、使う項目だけを生成したファイルに写します。data ディレクトリの場所は環境設定 (cp2k_data) に書きます":
        "Basis sets and pseudopotentials are read from the files in the CP2K data directory, and only the entries used are copied into the output. Set the data directory as cp2k_data in the settings",
    "単位系 (units)": "Units", "原子の形式 (atom_style)": "Atom style", "data ファイル": "Data file", "型番号の元素": "Elements of atom types",
    "写すファイル": "Files to copy", "read_data の前の行": "Commands before read_data", "pair_coeff の後の行": "Commands after pair_coeff", "乱数の種": "Random seed",
    "(選んでください)": "(choose one)",
    "空欄なら構造から data.lammps を書きます (atomic か charge のとき)": "Empty = data.lammps is written from the structure (atomic or charge style)",
    "型番号 1, 2, … の元素を順に空白区切りで (例 O H)。空欄なら構造の元素の順": "Elements of types 1, 2, … in order, space separated (e.g. O H). Empty = order in the structure",
    "例 eam、eam/alloy、reaxff NULL、mace no_domain_decomposition": "e.g. eam, eam/alloy, reaxff NULL, mace no_domain_decomposition",
    "pair_coeff の行 (行頭の pair_coeff は省けます。例 * * Cu_u3.eam)": "pair_coeff lines (the leading pair_coeff may be omitted; e.g. * * Cu_u3.eam)",
    "生成先へ写すファイル (力場・モデル) を 1 行に 1 つ。入力の中ではファイル名だけで書きます": "Files (force field, model) copied into the output, one per line; refer to them by file name only in the input",
    "例 bond_style harmonic / special_bonds lj/coul 0 0 0.5": "e.g. bond_style harmonic / special_bonds lj/coul 0 0 0.5",
    "例 kspace_style pppm 1e-4 / neigh_modify every 1": "e.g. kspace_style pppm 1e-4 / neigh_modify every 1",
    "相互作用は外部のもの (力場のファイル、data ファイル、機械学習ポテンシャルのモデル) をそのまま使います。全電荷とスピン多重度は使いません (0 と 1 のまま)。k 点はありません":
        "Interactions come from outside as is (force-field files, data files, machine-learning models). Total charge and multiplicity are not used (leave 0 and 1). There are no k-points",
    "pair_style": "pair_style", "pair_coeff": "pair_coeff", "define": "define",
    "トポロジー (.top)": "Topology (.top)", "構造のファイル (.gro / .pdb)": "Structure file (.gro / .pdb)", "静電相互作用 (coulombtype)": "Electrostatics (coulombtype)",
    "カットオフ [nm]": "Cut-offs [nm]", "拘束 (constraints)": "Constraints", "圧力浴 (pcoupl)": "Barostat (pcoupl)", "等温圧縮率 [1/bar]": "Compressibility [1/bar]",
    "前の段階の .cpt": "Previous stage .cpt", "初速の乱数の種 (gen-seed)": "Velocity seed (gen-seed)", "追加の mdp": "Extra mdp options",
    "(NPT のときに選びます)": "(choose for NPT)", "構造の群にも読み込む": "Also load into the structure group",
    "構造の群の作り方を「ファイル」にして、このファイルを読み込みます (原子数の確認と 3D 表示に使います)":
        "Sets the structure source to File and loads this file (used to check the atom count and for the 3D view)",
    "CHARMM-GUI、acpype、pdb2gmx などで作った .top": "A .top made with CHARMM-GUI, acpype, pdb2gmx, etc.",
    ".gro か .pdb (トポロジーと同じ原子の並び)": ".gro or .pdb (same atom order as the topology)",
    "例 -DPOSRES (位置の拘束)": "e.g. -DPOSRES (position restraints)",
    "前の段階 (NVT → NPT → 本計算) の .cpt。空欄なら初速を作ります": ".cpt of the previous stage (NVT → NPT → production). Empty = new velocities",
    "画面にない mdp の項目を 1 行に 1 つ、名前 = 値 の形で (例 nstcalcenergy = 100)": "mdp options not in the form, one per line: name = value (e.g. nstcalcenergy = 100)",
    "構造の正本は上の構造のファイル (.gro / .pdb) です。構造の群にも同じファイルを読み込んでください (原子数の確認と表示に使います)。全電荷とスピン多重度はトポロジーが決めるので使いません (0 と 1 のまま)。k 点はありません":
        "The structure file above (.gro / .pdb) is authoritative. Load the same file in the structure group too (used to check the atom count and for display). Total charge and multiplicity come from the topology and are not used (leave 0 and 1). There are no k-points",
}
_EN.update({
    "溶媒モデル (--alpb / --gbsa)": "Solvation model (--alpb / --gbsa)", "溶媒モデル (CPCM / SMD)": "Solvation model (CPCM / SMD)",
    "溶媒": "Solvent", "なし": "None", "溶媒のパラメータファイル (GBSA)": "Solvation parameter file (GBSA)",
    "空欄なら溶媒なし (param_gbsa_<溶媒>.txt の形のファイル)": "Empty = no solvent (a file like param_gbsa_<solvent>.txt)",
    "溶媒の比誘電率 (SCCS)": "Solvent permittivity (SCCS)",
    "元素ごとの初期磁気モーメント [μB]": "Initial magnetic moment per element [μB]", "DFT+U": "DFT+U", "LDAUTYPE": "LDAUTYPE",
    "starting_magnetization (元素ごと)": "starting_magnetization (per element)", "HUBBARD の射影": "HUBBARD projector",
    "元素ごとの MAGNETIZATION": "MAGNETIZATION per element", "PLUS_U_METHOD": "PLUS_U_METHOD",
})
_EN.update({
    "機械学習ポテンシャル (MACE・CHGNet)": "Machine-learning potential (MACE, CHGNet)",
    "機械学習ポテンシャルの種類": "Machine-learning potential", "モデル": "Model", "計算に使うデバイス (device)": "Device", "数値の精度 (dtype)": "Precision (dtype)",
    "D3 を足す (mace_mp の dispersion)": "Add D3 (the dispersion option of mace_mp)",
    "空欄ならパッケージの既定のモデル (例 small、medium-mpa-0、またはモデルのファイルのパス)":
        "Empty = the default model of the package (e.g. small, medium-mpa-0, or the path of a model file)",
    "空欄ならパッケージの既定 (例 cpu、cuda)": "Empty = the package default (e.g. cpu, cuda)",
    "(パッケージの既定)": "(package default)",
    "分子でも周期系でも使えます。全電荷とスピン多重度は使いません (0 と 1 のまま)。k 点はありません。"
    "計算は生成した run_mlip.py (ASE) が行い、ADIT 自身はモデルのパッケージを使いません":
        "Works for molecules and periodic systems. Total charge and multiplicity are not used (leave 0 and 1). There are no k-points. "
        "The generated run_mlip.py (ASE) does the calculation; ADIT itself does not use the model packages",
    "遷移状態の探索 (OptTS)": "Transition-state search (OptTS)", "最初にヘシアンを計算 (Calc_Hess)": "Compute the Hessian first (Calc_Hess)",
    "ヘシアンを計算し直す間隔 (Recalc_Hess)": "Hessian recalculation interval (Recalc_Hess)", "最後に振動数を計算 (Freq)": "Frequencies at the end (Freq)",
    "反応座標をたどる (IRC)": "Follow the reaction path (IRC)", "IRC の反復の上限": "Max IRC iterations", "IRC の向き": "IRC direction",
    "! OptTS にする (計算の種類は構造最適化)": "Use ! OptTS (with a geometry optimization)",
    "! Freq を足す": "Add ! Freq", "! Freq IRC にする (計算の種類は一点計算)": "Use ! Freq IRC (with a single point)",
    "遷移状態の探索は計算の種類が構造最適化のとき、IRC は一点計算のときに書かれます。"
    "両方を続けて実行するなら「実行」→「遷移状態と IRC」で段階に分けて生成します":
        "The transition-state search is written when the task is a geometry optimization, the IRC when it is a single point. "
        "To run both in turn, generate them as stages with Run → Transition state and IRC",
    "まとめて作る": "Batch generation", "組にして比べる": "Compare a set", "配座の候補": "Conformers", "反応経路 (NEB)": "Reaction path",
    "フォノン": "Phonons", "弾性定数": "Elastic constants", "遷移状態と IRC": "TS and IRC",
    "保存先": "Output directory", "組の種類": "Set kind", "分子の箱": "Molecule box", "スラブの構造": "Slab structure",
    "分子の構造": "Molecule structure", "吸着した構造": "Adsorbed structure", "画面の構造": "Screen structure",
    "スラブの固定原子": "Fixed atoms of the slab",
    "吸着した構造の固定原子": "Fixed atoms of the adsorbed structure", "反応に出てくる計算": "Runs of the reaction",
    "溶媒側で変える手法の項目": "Method items changed for the solvated run", "気相側で変える手法の項目": "Method items changed for the gas-phase run",
    "組の定義ファイル (set.json)": "Set file (set.json)",
    "作る配座の数": "Number of conformers to embed", "重複とみなす RMSD [Å]": "RMSD for duplicates [Å]", "力場": "Force field",
    "反復の上限": "Max iterations", "全原子で測る": "All atoms", "残す数": "Keep", "配座を作る SMILES": "SMILES for the conformers",
    "始状態の構造": "Initial-state structure", "終状態の構造": "Final-state structure", "中間の像の数": "Number of intermediate images",
    "NEB の作り方": "NEB mode", "補間": "Interpolation", "climbing image": "Climbing image",
    "超格子の倍率": "Supercell", "変位の大きさ [Å]": "Displacement [Å]", "変位の作り方": "How displacements are made",
    "DOS の q 点のメッシュ": "DOS q-point mesh", "DOS の幅 [THz]": "DOS width [THz]",
    "歪みの大きさ": "Strain magnitudes", "歪みを与える成分": "Strain components", "作るもの": "What to generate",
    "IRC もたどる (Sella)": "Also follow the IRC (Sella)",
    "行を追加": "Add row", "入れる": "Use", "生成": "Generate",
    "軌道変換法 (OT)": "Orbital transformation (OT)", "OT MINIMIZER": "OT MINIMIZER", "OT PRECONDITIONER": "OT PRECONDITIONER",
    "OT の追加行": "Additional OT lines", "SURFACE_DIPOLE_CORRECTION": "SURFACE_DIPOLE_CORRECTION", "SURF_DIP_DIR": "SURF_DIP_DIR",
    "SPRING (VASP)": "SPRING (VASP)", "K_SPRING (CP2K)": "K_SPRING (CP2K)", "opt_scheme (QE)": "opt_scheme (QE)",
})
_EN["指定"] = "Molecule"
_EN["設定"] = "Settings"

from adit.gui import analysis_fields as _AF  # noqa: E402  

from adit.gui import report_fields as _RF  # noqa: E402
from adit.web import prep23 as _P23  # noqa: E402

for _ja, _en in [*_AF.LABELS.values(), *_AF.PLACEHOLDERS.values(), *_RF.LABELS.values(), *_P23.LABELS.values(),
                 *((ja, en) for items in (_AF.THERMO_MODELS, _AF.GEOMETRIES, _AF.IMAGINARY, _AF.UV_SHAPES) for _, ja, en in items)]:
    _EN.setdefault(_ja, _en)
_EN.setdefault("最大ステップ数 (MaxSteps)", "Maximum steps (MaxSteps)")
_EN.setdefault("SCC の収束判定 (SccTolerance)", "SCC tolerance (SccTolerance)")
_EN.setdefault("SCC の反復上限 (MaxSccIterations)", "Maximum SCC iterations (MaxSccIterations)")


def set_language(lang: str) -> None:
    global LANGUAGE
    LANGUAGE = "en" if str(lang).lower().startswith("en") else "ja"
    from adit import lang as _lang
    _lang.set_language(LANGUAGE)


def language_from_env(default: str = "ja") -> str:
    return env_var("LANG", default)


def tr(text: str) -> str:
    if LANGUAGE == "en":
        return _EN.get(text, text)
    return text


def _tr_tip(action) -> str:
    # Keep the " (Ctrl+Z)" suffix out of the lookup so the tooltip before it still translates.
    from adit.gui.palette import split_shortcut

    base, suffix = split_shortcut(action.toolTip(), action)
    return tr(base) + suffix


def translate_widgets(root: QWidget) -> None:
    if LANGUAGE != "en":
        return
    for w in [root, *root.findChildren(QWidget)]:
        if isinstance(w, QGroupBox):
            w.setTitle(tr(w.title()))
        if isinstance(w, (QLabel, QAbstractButton)):
            w.setText(tr(w.text()))
        if isinstance(w, QComboBox):
            for i in range(w.count()):
                w.setItemText(i, tr(w.itemText(i)))
        if isinstance(w, (QLineEdit, QPlainTextEdit)):
            w.setPlaceholderText(tr(w.placeholderText()))
        if isinstance(w, QTabWidget):
            for i in range(w.count()):
                w.setTabText(i, tr(w.tabText(i)))
        if isinstance(w, QTabBar) and not isinstance(w.parentWidget(), QTabWidget):
            for i in range(w.count()):
                w.setTabText(i, tr(w.tabText(i)))
        if isinstance(w, QToolButton) and w.defaultAction() is not None:
            a = w.defaultAction()
            a.setText(tr(a.text())); a.setIconText(tr(a.iconText())); a.setToolTip(_tr_tip(a))
        if isinstance(w, (QToolBar, QMenu, QMenuBar)):
            for a in w.actions():
                a.setText(tr(a.text())); a.setToolTip(_tr_tip(a))
                if a.menu() is not None:
                    a.menu().setTitle(tr(a.menu().title()))
        if w.toolTip():
            w.setToolTip(tr(w.toolTip()))
    if hasattr(root, "windowTitle"):
        root.setWindowTitle(tr(root.windowTitle()))

_EN.update({
    "原子の選び方": "Atom selection",
    "RDF の元素の組": "RDF element pairs",
    "密度分布を取る軸": "Axis for the density profile",
    "配位数のカットオフ [Å]": "Coordination cutoff [Å]",
    "中心対称性の相手の数": "Centrosymmetry neighbors",
    "Steinhardt q4・q6 のカットオフ [Å]": "Steinhardt q4 and q6 cutoff [Å]",
    "かたまりのカットオフ [Å]": "Cluster cutoff [Å]",
    "結合角の分布の元素": "Elements for the angle distribution",
    "結合角の分布のカットオフ [Å]": "Angle-distribution cutoff [Å]",
    "構造因子 S(q) を出す": "Structure factor S(q)",
    "水素結合 (距離 [Å], 角度 [度])": "Hydrogen bonds (distance [Å], angle [deg])",
    "慣性半径 Rg の時系列": "Radius of gyration over time",
    "3 次元の数密度の格子": "3D number-density grid",
    "Voronoi の体積と面の数": "Voronoi volumes and face counts",
    "数える面の面積の下限 [Å²]": "Smallest face area to count [Å²]",
    "溶媒接触表面積 (半径の表, プローブ [Å])": "Solvent-accessible surface (radii set, probe [Å])",
    "距離の時系列 (原子の組)": "Distance series (atom pairs)",
    "角度の時系列 (原子の組)": "Angle series (atom triples)",
    "二面角の時系列 (原子の組)": "Dihedral series (atom quadruples)",
    "RMSD の基準フレーム": "Reference frame for the RMSD",
    "原子ごとの揺らぎ (RMSF)": "Fluctuation per atom (RMSF)",
    "速度自己相関と振動スペクトル": "Velocity autocorrelation (VACF) and spectrum",
    "イオンの電荷 (伝導度)": "Ion charge (conductivity)",
    "伝導度の温度 [K]": "Temperature for the conductivity [K]",
    "変位の基準フレーム": "Reference frame for the displacement",
    "局所ひずみのカットオフ [Å]": "Local-strain cutoff [Å]",
    "主成分分析の成分の数": "Number of principal components",
    "分類する群の数 (k-means)": "Number of groups (k-means)",
    "自由エネルギー面の温度 [K]": "Temperature for the free-energy surface [K]",
    "自由エネルギー面の区間の数": "Bins for the free-energy surface",
    "自由エネルギーの単位": "Unit for the free energy",
    "バンド図の縦軸の幅 [eV]": "Band-plot energy window [eV]",
    "有効質量に使う k 点の数": "k-points used for the effective mass",
    "Bader の ACF.dat": "Bader ACF.dat",
    "Bader の価電子数": "Valence electrons for Bader",
    "粉末回折の線源": "Radiation for the powder pattern",
    "粉末回折の 2θ の範囲 [度]": "2-theta range for the powder pattern [deg]",
    "重ねる測定データのファイル": "Measured pattern to overlay",
    "粘度 (Green-Kubo)": "Viscosity (Green-Kubo)",
    "面平均を取る軸": "Axis for the plane average",
    "cube の値の単位": "Unit of the cube values",
    "仕事関数を出す": "Work function",
    "その場で計算する上限 [秒]": "Limit for computing here [s]",
})

_EN.update({
    "計算コードの指定": "Code (when there is no input file)",
    "振動数の補正係数": "Frequency scaling factor",
    "重ねる測定したスペクトル": "Measured spectrum to overlay",
})

_EN.update({
    "報告をまとめる": "Build a report",
    "計算ディレクトリ (1 行に 1 つ)": "Run directories (one per line)",
    "報告の言語": "Report language",
    "方法の節 (Markdown)": "Methods section (Markdown)",
    "条件の表 (CSV)": "Conditions table (CSV)",
    "結果の表 (CSV)": "Results table (CSV)",
    "再現パッケージ (.zip)": "Reproducibility bundle (.zip)",
    "入力のハッシュを照合する": "Verify the input fingerprints",
    "報告を作る": "Build",
})

_EN.update({
    "図の線の色": "Line colors",
    "目盛りの向き": "Tick direction",
    "目盛り線": "Grid lines",
    "グラフを囲む枠": "Axis frame",
    "線の太さ [pt]": "Line width [pt]",
    "図の文字の大きさ [pt]": "Font size in figures [pt]",
    "図の解像度 [dpi]": "Figure resolution [dpi]",
})

_EN["図の追加の形式"] = "Extra figure formats"
