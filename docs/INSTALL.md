# Installation for each OS

These are the steps from installing ADIT to running your first calculation.

## The easiest way: Using the executables (Windows / macOS)

**You don't need Python or WSL.**

1. Open [Releases](https://github.com/adit-chem-project/adit-chem/releases)
2. Download `ADIT-windows-x64.zip` for Windows, or `ADIT-arm64.zip` for macOS (Apple Silicon)
   (there is no build for Intel Macs; install with pip as described below)
3. Extract it, then double-click `ADIT.exe` on Windows or `ADIT.app` on macOS

The executables are not signed, so a warning appears the first time only. On Windows, choose "More info" → "Run anyway".
On macOS, right-click the app in Finder and choose "Open".

To use the command line (`adit-gen` and so on), call `adit-cli.exe` in the same folder on Windows.
The macOS zip contains only `ADIT.app`, so call `ADIT.app/Contents/MacOS/adit-cli` in a terminal
(for example: `./ADIT.app/Contents/MacOS/adit-cli gen --list-samples`).

If you also want to run calculations on the same PC, go on to the next section.

## Running calculations on the same PC (Windows; follow these steps in order)

Even if you have never used Linux, following this section from top to bottom gets you as far as creating your first input files and running a calculation.
Copy the commands **one line at a time**, paste them into the terminal (the black window where you give the computer instructions as text), and press Enter.
The part from `#` to the end of a line is only an explanation, so you don't have to paste it. Replace `<...>` with your own values.

### 1. Install WSL (which runs Linux inside Windows)

1. In the Start menu, right-click "PowerShell" and choose "Run as administrator"
2. Type the following line and press Enter. When it finishes, restart the PC

   ```powershell
   wsl --install
   ```

3. After the restart, an "Ubuntu" window opens. Choose a user name and password for Linux (these are separate from your Windows ones, and the password does not appear on screen as you type)

From here on, work in the window (the terminal) that "Ubuntu" in the Start menu opens.

### 2. Install miniforge (conda)

conda is a system that installs Python and calculation codes into separate boxes called "environments". miniforge is a distribution of it that uses conda-forge (where DFTB+, xtb and Quantum ESPRESSO are distributed) by default.

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh      # answer the questions with Enter and yes
```

When it finishes, close the terminal and open it again. If `(base)` appears at the start of the line, it is installed.

> **Why conda rather than venv**: ADIT runs on Python alone, but running calculations on this PC needs a calculation code such as DFTB+.
> With conda, you can install those in the same way. WSL's Ubuntu starts with only `python3` and no `python` command, but activating a conda environment gives you `python`.
> If you only run calculations on a cluster and only generate input locally, a venv made with `python3 -m venv adit-env` and `source adit-env/bin/activate` is also fine.

### 3. Create an environment for ADIT and DFTB+

```bash
conda create -n adit -c conda-forge "python>=3.11" dftbplus rdkit font-ttf-noto-cjk git
conda activate adit                         # the start of the line changes to (adit); type this each time you open a terminal
```

`font-ttf-noto-cjk` is the font for showing Japanese in the desktop app (WSL has no Japanese font at first).
Install xtb and Quantum ESPRESSO in separate environments, because installing them alongside DFTB+ can downgrade DFTB+ to an older version (you only need the ones you will use).

```bash
conda create -n xtb -c conda-forge xtb
conda create -n qe -c conda-forge qe
conda activate --stack xtb                   # example: make xtb available too, with the adit environment still active
```

### 4. Install ADIT

```bash
git clone https://github.com/adit-chem-project/adit-chem.git
cd adit-chem && pip install -e ".[gui,smiles,analysis]"
```

`[gui,smiles,analysis]` installs the packages needed for the desktop app's GUI, SMILES and Draw, and space-group analysis.

### 5. Get the parameters and put them in place

Some calculation codes need parameter files for each element. ADIT does not include them, so get them from their distributors.
Files downloaded with a browser on Windows appear inside WSL at `/mnt/c/Users/<your Windows user name>/Downloads/`.

**Slater-Koster parameters for DFTB+** (`.skf` files for each element pair, licensed under CC BY-SA 4.0)

Where to get them: https://dftb.org/parameters/download.html (for water and organic molecules, start with mio; the elements mio covers are H, C, N, O, S, P)

```bash
mkdir -p ~/slakos                                                   # create the location (~ is your home folder)
tar -xf /mnt/c/Users/<your Windows user name>/Downloads/mio-1-1.tar.xz -C ~/slakos
ls ~/slakos/mio-1-1                                                 # if you see H-H.skf and so on, they are in place
```

After extraction it looks like this: the `.skf` files, `LICENSE` and `README` sit directly inside the folder of the set (`mio-1-1/`).

```
~/slakos/
  mio-1-1/
    H-H.skf  H-O.skf  O-H.skf  O-O.skf  C-C.skf  …  LICENSE  README
```

**Pseudopotentials for Quantum ESPRESSO (UPF files)**

One source: SSSP https://www.materialscloud.org/discover/sssp (check the license with the distributor).
Extract it and place the folder that directly contains the `.UPF` files, for example as `~/pseudo/SSSP_efficiency/`. The folder name becomes the "Pseudopotential set" in ADIT.

```
~/pseudo/
  SSSP_efficiency/
    Si.pbe-n-rrkjus_psl.1.0.0.UPF  O.pbe-n-kjpaw_psl.0.1.UPF  …
```

xtb needs no parameter files, because the parameters are built into the method. For VASP POTCAR files and ORCA, see "What to prepare for each code" in [Settings](SETTINGS.md#what-to-prepare-for-each-code).

### 6. Write the locations in the settings file

The first time you run `adit-web` or `adit-gen`, the settings file `~/.config/adit/cluster.toml` is created and its location is shown.
Open it with the commands below and write the locations (nano is a simple text editor for the terminal. Ctrl+O and then Enter saves; Ctrl+X exits).

```bash
echo $HOME                                   # shows your home folder (for example /home/taro)
nano ~/.config/adit/cluster.toml
```

```toml
sk_root = "/home/<user name>/slakos"        # the location created in step 5; one level above the folder of the set (mio-1-1)
pseudo_root = "/home/<user name>/pseudo"    # when you use Quantum ESPRESSO
```

Always enclose paths in `"..."`. Windows-style paths (`C:\Users\...`) cannot be read as they are. Inside WSL, `C:\Users\taro` is `/mnt/c/Users/taro`.

### 7. Create your first input and run it

```bash
adit-web --open                             # if no browser opens, open http://127.0.0.1:8765/ in a browser on Windows
```

In the browser, choose H2O under the "Preset" structure, the code DFTB+, the "Slater-Koster set" mio-1-1, and leave "Profile" as local. Go through "Preview" → "Generate", run `bash submit.sh` in the terminal, then press "Go to analysis", and on the analysis page press "Run analysis".
To work only in the terminal, you can create the run from calculation settings (spec.json) saved on screen, like this:

```bash
adit-gen spec.json ~/adit_runs/water       # writes the full set of input files to ~/adit_runs/water
cd ~/adit_runs/water && bash submit.sh      # runs the calculation (the log is output.log)
python analyze.py                            # writes figures and a summary to analysis/
```

The `README.txt` in each directory describes the steps for that directory.

## Building the executables (.exe / .app) yourself

The distributed executables are built with this procedure. When you push a tag starting with `v` (for example `v0.1.0`),
GitHub Actions builds the Windows `.exe` and the macOS `.app` and attaches them to Releases.
To build locally, follow the steps below. This configuration has actually been built on Linux and checked as far as the CLI running (2026-09-13).

### Windows

#### 1 Build on a Windows PC

In PowerShell, go to where you put the repository, then run:

```powershell
py -3.12 -m venv build-env
build-env\Scripts\activate
pip install -e ".[gui,smiles,analysis]" pyinstaller
pyinstaller packaging\adit.spec --noconfirm
```

`dist\adit\` will contain **two executables**. Distribute the whole folder (the files inside it are needed).

| File | Use |
|---|---|
| `ADIT.exe` | The GUI (desktop app). No black window appears, but in exchange it **has no standard output, so it cannot be used for the CLI** |
| `adit-cli.exe` | The CLI. `adit-cli.exe gen ...` / `analyze` / `report` / `convert` / `web` |

There are two because of how Windows works. An executable for the GUI (GUI subsystem) has no standard output,
and **with a single executable the CLI shows nothing** (we confirmed this by actually building it).
To combine them into one .exe, change `ONEFILE = False` to `True` in `packaging/adit.spec`
(it extracts to a temporary folder at every launch, so opening takes several seconds).

#### 2 Checking the built executables (whoever builds them must go through this)

```powershell
dist\adit\ADIT.exe                             # does the GUI open? Is any Japanese text shown as □?
dist\adit\adit-cli.exe gen --list-samples      # is the list of samples shown? (can it read the bundled examples?)
dist\adit\adit-cli.exe gen --sample water_generated mine.json
dist\adit\adit-cli.exe gen mine.json out\run1  # are the input, submit.sh and README.txt written?
dist\adit\adit-cli.exe analyze out\run1        # in a directory not yet run, does it give the reason and stop? (exit code 2)
dist\adit\adit-cli.exe web --open              # does http://127.0.0.1:8765 open in a browser? (without --open, no browser opens)
```

**If Japanese on screen appears as □ (tofu)**, the font is not bundled. In a conda environment, run
`conda install -c conda-forge font-ttf-noto-cjk` and rebuild, or set `ADIT_FONT_DIR` to
a folder containing `NotoSansCJKjp-VF.ttf` and rebuild.

#### 3 What is bundled and what is not

| Bundled | Why |
|---|---|
| `adit/scripts/templates/*.j2`, `adit/web/templates/*.html` | Templates for submit.sh and the screens |
| `examples/*/spec.json` | Read by `gen --list-samples` / `--sample` |
| `fonts/NotoSansCJKjp-VF.ttf` (if present) | Japanese on screen and in figures |

| Not bundled | Why |
|---|---|
| Calculation codes such as DFTB+, xtb and Quantum ESPRESSO | They have separate distribution terms. **Calculations are not run on Windows** (the generated input is transferred to a Linux server), so they are not needed in the first place |
| Slater-Koster sets, pseudopotentials, POTCAR | Under their licenses, ADIT must not distribute them |

#### 4 Known limitations

- **Calculations cannot be run on Windows.** From the Workspace terminal (PowerShell), log in to the cluster with `ssh`,
  then transfer and submit with the commands in `transfer_and_submit.sh`. What the exe can do is generate input and analyze output that is on the local machine
- **The size is 300 to 500 MB** (because PySide6, matplotlib and SciPy are bundled). `ONEFILE = True` does not make it smaller
- **It is not signed.** SmartScreen warns of an "Unknown publisher". If you distribute it, sign it or tell recipients
  to choose "More info" → "Run anyway"
- **Built in an environment with RDKit, it can create structures from SMILES.** Built without RDKit, only that field is unavailable
- What was checked on Linux is six things, `gen --list-samples` / `--sample` / generating input / `analyze` (exit code 2 on a directory not yet run) /
  `report` / `convert verify`, plus **the GUI starting and creating a window** (`QT_QPA_PLATFORM=offscreen`).
  **Files produced on Windows have not been checked** (Windows executables can only be built on Windows).
  Be sure to go through "2 Checking the built executables" above
- Four things we learned by building it (all are handled in the spec):
  1. **ASE loads its per-format modules dynamically by name**, so unless `ase.io` is in `hiddenimports`,
     writing `.gen` fails with `UnknownFileTypeError`
  2. On Linux, **the GUI does not start unless conda's `libOpenGL.so.0` is bundled** (PySide6 on Windows
     has its own OpenGL DLL, so this step only takes effect on Linux)
  3. **An executable for the GUI has no standard output**, so a separate console version (`adit-cli.exe`) is built for the CLI
  4. The default window size (1800x1000) **did not fit on a 1536x864 screen**, so it now fits within 95 % of the screen


### macOS

#### Building

Run it by hand from GitHub Actions.

When you push a tag starting with `v` (for example `v0.1.0a1`), it is built automatically and `ADIT-arm64.zip` is attached to Releases.
You can also run it by hand.

1. In the repository, go to **Actions** → **macos-app** → **Run workflow**
2. One Apple Silicon (arm64) job runs
3. When it finishes, `adit-macos-arm64` appears under the artifacts. It contains
   - `ADIT-arm64.zip` … `ADIT.app`, which opens with a double-click (packed with `ditto`, which preserves Finder permissions)
   - `dist/adit/` … the set for people who don't use the GUI (it contains `adit-cli`)

To build on your own Mac:

```bash
python -m pip install -e ".[gui,smiles,analysis]" pyinstaller
pyinstaller packaging/adit.spec --noconfirm
open dist/ADIT.app          # GUI
./dist/adit/adit-cli gen --list-samples   # command line
```

#### No code signing or notarization

This requires a paid Apple developer registration, so **the app has neither code signing (codesign) nor notarization.**
For that reason, recipients need one of the following the first time they open it.

- In Finder, **right-click `ADIT.app` → Open**, then choose "Open" in the confirmation that appears
- Or run `xattr -dr com.apple.quarantine /Applications/ADIT.app`

**Always include these steps in the instructions that come with the distribution.** Otherwise macOS says the app is damaged and can't be opened,
and users have no way of knowing why.

To sign it, enroll in the Apple Developer Program, then run:

```bash
codesign --deep --force --options runtime --sign "Developer ID Application: <name> (<team ID>)" dist/ADIT.app
xcrun notarytool submit dist/ADIT-arm64.zip --apple-id <id> --team-id <team> --password <app-specific password> --wait
xcrun stapler staple dist/ADIT.app
```

#### Checks (run these first on a real machine)

There are 7 checks: the 6 for Windows (in "2 Checking the built executables" above) plus 1 on another Mac. **Mark the ones that pass and update this document.**

1. Double-clicking `ADIT.app` shows the GUI (Japanese does not turn into tofu □)
2. A structure can be built and input generated (`adit-cli gen --sample water_generated` → generate → `README.txt` exists)
3. Analysis gets as far as figures (`adit-cli analyze <directory>`)
4. Draw opens and a drawn molecule becomes SMILES (is RDKit bundled?)
5. The web app opens (`adit-cli web` → 127.0.0.1:8765 in a browser)
6. The 3D view rotates (the "Structure" tab in the GUI, and "Structure (3D)" in the web app)
7. After copying it to another Mac (not the machine it was built on), checks 1 to 6 pass

**Check 7 is the most important.** It can work on the machine it was built on and still fail on other machines because something is missing.

#### Intel Macs

**Only the Apple Silicon (arm64) build is distributed.** GitHub has stopped providing its Intel runner (`macos-13`),
so CI cannot build it. On an Intel Mac, install with pip.
We do not make a combined universal2 build, because that is only possible if PySide6 and RDKit ship wheels for both architectures.

## Installing and starting (Linux / macOS / experienced users)

Python 3.11 or later is required. Create a virtual environment (conda or venv) and install with pip. If you are on Windows and new to this, start from the sections above.
The dependencies are ASE, pydantic, Jinja2, matplotlib, SciPy and tomli-w, plus, for the desktop app, PySide6, pyqtdarktheme-fork,
and pyte and ptyprocess for the Workspace terminal (pywinpty on Windows). To use SMILES and Draw (the molecule-drawing feature), you also need RDKit.

```bash
python3 -m venv adit-env                # if you use conda, follow step 3 in the section above
source adit-env/bin/activate            # on Windows (PowerShell):  adit-env\Scripts\Activate.ps1
git clone https://github.com/adit-chem-project/adit-chem.git
cd adit-chem && pip install -e ".[gui,smiles,analysis]"
```

The distribution name is `adit-chem` (`adit` on PyPI is a different package). The import name and the command name are `adit`.

| Use | Command | Notes |
|---|---|---|
| Desktop app | `adit` | Requires PySide6 |
| Web app (browser) | `adit-web --open` | Does not need PySide6. Opens `http://127.0.0.1:8765/` |
| Command line (generate) | `adit-gen spec.json <output directory>` | Generates from calculation settings saved in the GUI |
| Command line (convert) | `adit-convert structure ...` / `adit-convert calculation ...` | Converts a structure format, or generates input for another calculation code while keeping the shared settings |
| Command line (analyze) | `adit-analyze <run directory> --rdf --msd --dos` | Writes figures and a summary to `analysis/` |

### Converting structure formats and calculation codes

Structure files can be converted between the formats ASE supports. The input and output formats are normally inferred from the file names.
The same feature is under "Convert" on the "File" tab in the desktop app, and under "Convert" at the top in the web app. For file names that cannot be inferred, you can give ASE format names in the "Input format (optional)" and "Output format (optional)" fields.
In both, you can save the current structure previewed on the preparation screen ("Prepare" in the web app) as extended XYZ. The cell, periodic boundaries, fixed atoms and axis constraints, and any velocities already set are also kept. The CLI reads the same information when it converts a structure from `spec.json` or a generated directory.

```bash
adit-convert structure data.lammps POSCAR
adit-convert structure POSCAR structure.data
adit-convert structure POSCAR structure.xyz
adit-convert structure trajectory.extxyz final.cif --frame -1
adit-convert crystal "225 Na:0,0,0 Cl:0.5,0,0 cell=5.64" nacl.cif        # builds a crystal from the space group and the atoms' fractional coordinates
adit-convert surface "fcc111 Al 2x2x3 vacuum=10" slab.xyz                # builds a surface slab from the plane, element and size
adit-convert combine molecule_a.xyz molecule_b.xyz combined.xyz
adit-convert openbabel source.sdf target.mol2 --input-format sdf --output-format mol2
adit-convert dock6 prepared/dock.in ready/
```

GOAT, ORCA DOCKER, DCDFTBMD 2.0, DOCK6 input preparation and the Open Babel conversion entry point are covered, with steps and limits, in [Supported software](SOFTWARE.md). They are command-line only; the desktop and web apps have no fields for them yet.
