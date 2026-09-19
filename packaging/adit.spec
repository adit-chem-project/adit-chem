# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller build configuration for the Windows .exe, the macOS .app and the Linux executable.

    pyinstaller packaging/adit.spec      # writes dist/adit/

Bundled: Jinja templates, example spec.json files, an optional CJK font with the full text of its
licence, LICENSE, THIRD_PARTY_NOTICES.md and licenses/*.txt.
Not bundled: the simulation codes themselves, Slater-Koster sets, pseudopotentials and POTCAR files.
"""

import os
from pathlib import Path

from PyInstaller.utils.hooks import collect_data_files, collect_submodules

ONEFILE = False          # Keep this False: PySide6 (Qt) is LGPL-3.0 and users must be able to replace
                         # the Qt libraries, which the one-directory layout allows.

# Two executables are produced:
#   ADIT.exe      the window application (console=False, so it has no standard output on Windows)
#   adit-cli.exe  the command line (console=True): gen / analyze / report / convert / web

ROOT = Path(os.getcwd())
datas = collect_data_files("adit")                     # Jinja templates
datas += collect_data_files("ase", include_py_files=False)

examples = ROOT / "examples"
if examples.is_dir():                                    # sample settings only (spec.json)
    datas += [(str(p), str(Path("examples") / p.relative_to(examples).parent))
              for p in examples.rglob("spec.json")]

# CJK font. When it is bundled, the full text of the SIL Open Font License 1.1 must go with it.
for name in ("NotoSansCJKjp-VF.ttf", "NotoSansCJKjp-Regular.otf"):
    font = Path(os.environ.get("ADIT_FONT_DIR", Path(os.sys.prefix) / "fonts")) / name
    if font.is_file():
        datas += [(str(font), "fonts")]
        ofl = ROOT / "licenses" / "OFL-1.1-NotoSansCJK.txt"
        if not ofl.is_file():   # without the licence text the font may not be redistributed
            raise SystemExit(f"{ofl} is missing; the font cannot be bundled without the OFL text")
        datas += [(str(ofl), "fonts")]
        break

# Ship the third-party licence texts with the build
for extra in (ROOT / "LICENSE", ROOT / "THIRD_PARTY_NOTICES.md"):
    if extra.is_file():
        datas += [(str(extra), ".")]
licenses_dir = ROOT / "licenses"
if licenses_dir.is_dir():
    datas += [(str(p), "licenses") for p in licenses_dir.glob("*.txt")]

# ASE imports its per-format modules by name, so PyInstaller cannot find them on its own.
hiddenimports = (collect_submodules("adit") + collect_submodules("ase.io") + collect_submodules("ase.calculators")
                 + ["scipy.spatial.transform._rotation_groups"])

# Linux only: bundle the conda OpenGL/EGL libraries when they are present.
binaries = []
if os.sys.platform.startswith("linux"):
    for lib in ("libOpenGL.so.0", "libEGL.so.1", "libGLdispatch.so.0", "libGLX.so.0"):
        found = Path(os.sys.prefix) / "lib" / lib
        if found.is_file():
            binaries.append((str(found), "."))

a = Analysis(
    [str(ROOT / "packaging" / "launch.py")],
    pathex=[str(ROOT / "src")],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    excludes=["tkinter", "pytest", "IPython", "notebook"],
    noarchive=False,
)
pyz = PYZ(a.pure)

if ONEFILE:
    exe = EXE(pyz, a.scripts, a.binaries, a.datas, [], name="ADIT", console=False,
              icon=None, upx=False, disable_windowed_traceback=False)
    cli = EXE(pyz, a.scripts, a.binaries, a.datas, [], name="adit-cli", console=True,
              icon=None, upx=False, disable_windowed_traceback=False)
else:
    exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name="ADIT", console=False,
              icon=None, upx=False, disable_windowed_traceback=False)
    cli = EXE(pyz, a.scripts, [], exclude_binaries=True, name="adit-cli", console=True,
              icon=None, upx=False, disable_windowed_traceback=False)
    coll = COLLECT(exe, cli, a.binaries, a.datas, strip=False, upx=False, name="adit")

# macOS only: also produce a double-clickable .app. It is neither signed nor notarised,
# so the first launch needs right-click -> Open in Finder (docs/INSTALL.md).
if os.sys.platform == "darwin" and not ONEFILE:
    app = BUNDLE(coll, name="ADIT.app", icon=None, bundle_identifier="org.adit.app",
                 info_plist={
                     "CFBundleName": "ADIT",
                     "CFBundleDisplayName": "ADIT",
                     "CFBundleShortVersionString": os.environ.get("ADIT_VERSION", "0.1.0a1"),
                     "NSHighResolutionCapable": True,
                     "LSApplicationCategoryType": "public.app-category.education",
                 })
