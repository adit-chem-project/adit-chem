"""Convert structure formats, and convert a CalculationSpec from one calculation code to another."""

from __future__ import annotations

from adit.errors import AditValueError
import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import warnings
from pathlib import Path

from adit import lang
from adit.config import ConfigError, config_path, ensure_config, env_var
from adit.lang import L
from adit.project import OutputNotEmpty, ProjectError, write_project
from adit.provenance import sha256_file
from adit.spec import CalculationSpec
from adit.templates import TemplateError, load_template

REPORT_FILE = "conversion.json"


class ConversionError(AditValueError):
    pass


class StructureDataLossWarning(UserWarning):
    """A structure writer cannot retain data present in the source."""


def convert_with_openbabel(source: Path | str, output: Path | str, *,
                           input_format: str | None = None, output_format: str | None = None,
                           overwrite: bool = False, executable: str = "obabel") -> Path:
    """Use the explicitly selected Open Babel backend; do not infer chemistry options."""
    # Open Babel CLI: https://openbabel.org/docs/Command-line_tools/babel.html
    # Bond typing and format-specific metadata are Open Babel's responsibility;
    # this operation does not promise a lossless structural round trip.
    src, dst = Path(source).expanduser(), Path(output).expanduser()
    if not src.is_file():
        raise ConversionError(L(f"入力ファイルがありません: {src}", f"input file does not exist: {src}"))
    if src.resolve() == dst.resolve():
        raise ConversionError(L("入力ファイルと出力ファイルを同じファイルにはできません", "the input and output must not be the same file"))
    if dst.exists() and not overwrite:
        raise ConversionError(L(f"出力先 {dst} はすでにあります (上書きしません)", f"output {dst} already exists (not overwritten)"))
    if not dst.suffix and not output_format:
        raise ConversionError(L("出力形式を --output-format で指定してください", "specify the output format with --output-format"))
    if not src.suffix and not input_format:
        raise ConversionError(L("入力形式を --input-format で指定してください", "specify the input format with --input-format"))
    exe = shutil.which(executable)
    if exe is None:
        raise ConversionError(L(f"Open Babel の実行ファイル {executable!r} が見つかりません", f"Open Babel executable {executable!r} was not found"))
    dst.parent.mkdir(parents=True, exist_ok=True)
    fd, staged_name = tempfile.mkstemp(prefix=".adit-obabel-", suffix=dst.suffix, dir=dst.parent)
    os.close(fd)
    staged = Path(staged_name)
    command = [exe]
    if input_format:
        command += ["-i", input_format]
    command.append(str(src))
    if output_format:
        command += ["-o", output_format]
    command += ["-O", str(staged)]
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=False)
    except OSError as ex:
        raise ConversionError(L(f"Open Babel を実行できません。元の出力は変更していません。途中ファイル: {staged}。{ex}",
                                f"cannot run Open Babel. The previous output was not changed. Staged file: {staged}. {ex}")) from ex
    if result.returncode != 0 or not staged.is_file() or staged.stat().st_size == 0:
        detail = (result.stderr or result.stdout).strip()
        raise ConversionError(L(f"Open Babel の変換に失敗しました。元の出力は変更していません。途中ファイル: {staged}。{detail}",
                                f"Open Babel conversion failed. The previous output was not changed. Staged file: {staged}. {detail}"))
    try:
        if dst.exists() and not overwrite:
            raise ConversionError(L(f"出力先 {dst} が変換中に作成されました (上書きしません)。途中ファイル: {staged}",
                                    f"output {dst} was created during conversion (not overwritten). Staged file: {staged}"))
        os.replace(staged, dst)
    except OSError as ex:
        raise ConversionError(L(f"変換結果を {dst} に移せません。途中ファイル: {staged}。{ex}",
                                f"cannot move the converted file to {dst}. Staged file: {staged}. {ex}")) from ex
    return dst


def combine_xyz(sources: list[Path | str], output: Path | str, *, overwrite: bool = False) -> Path:
    """Concatenate single-frame, non-periodic XYZ files without changing coordinates."""
    # Equivalent in scope to CMMDE's combinexyz preparation command:
    # https://git.mki.or.id/CoreDev/CMMDE/src/commit/5e0e0a58482e371279e87c373536390ec9d11309/bin/cmmdepre.py
    import numpy as np
    from ase import Atoms
    from ase.io import read, write
    from adit.mixture import MAX_ATOMS

    if len(sources) < 2:
        raise ConversionError(L("結合する XYZ ファイルを 2 つ以上指定してください",
                                "provide at least two XYZ files to combine"))
    dst = Path(output).expanduser()
    if dst.suffix.lower() != ".xyz":
        raise ConversionError(L("結合後のファイル名は .xyz にしてください",
                                "the combined output filename must end in .xyz"))
    paths = [Path(item).expanduser() for item in sources]
    if any(path.suffix.lower() != ".xyz" for path in paths):
        raise ConversionError(L("結合できる入力は単一フレームの .xyz ファイルだけです",
                                "only single-frame .xyz files can be combined"))
    if any(path.resolve() == dst.resolve() for path in paths):
        raise ConversionError(L("出力先を、入力ファイルのどれかと同じファイルにはできません",
                                "the output must not be one of the input files"))
    if dst.exists() and not overwrite:
        raise ConversionError(L(f"出力先 {dst} はすでにあります (上書きしません)",
                                f"output {dst} already exists (not overwritten)"))
    combined = Atoms()
    for path in paths:
        try:
            with path.open(encoding="utf-8") as stream:
                try:
                    n = int(stream.readline().strip())
                except ValueError as ex:
                    raise ConversionError(L(f"{path} の先頭行には原子数を整数で書いてください",
                                            f"the first line of {path} must contain an integer atom count")) from ex
                if n < 1:
                    raise ConversionError(L(f"{path} に原子がありません", f"{path} contains no atoms"))
                if len(combined) + n > MAX_ATOMS:
                    raise ConversionError(L(f"結合後は {len(combined) + n:,} 原子になり、上限 {MAX_ATOMS:,} 原子を超えます",
                                            f"the combined structure would have {len(combined) + n:,} atoms, above the {MAX_ATOMS:,}-atom limit"))
                comment = stream.readline()
                if any(key in comment for key in ("Lattice=", "pbc=", "Properties=")):
                    raise ConversionError(L(f"{path} にはセルなどの拡張 XYZ 情報があります。combine はそれを保持できないため、通常の XYZ だけを受け付けます",
                                            f"{path} contains extended XYZ metadata such as a cell; combine cannot preserve it and accepts plain XYZ only"))
                if any(not stream.readline() for _ in range(n)):
                    raise ConversionError(L(f"{path} の原子行が足りません", f"{path} has fewer atom lines than declared"))
                if any(line.strip() for line in stream):
                    raise ConversionError(L(f"{path} には複数フレームか余分な行があります。1 フレームずつ分けてください",
                                            f"{path} contains multiple frames or extra lines; split it into single-frame files"))
            atoms = read(path, format="xyz")
            if len(atoms) != n or not np.isfinite(atoms.positions).all():
                raise ConversionError(L(f"{path} の原子数か座標が正しくありません",
                                        f"{path} has an invalid atom count or coordinates"))
            combined += atoms
        except ConversionError:
            raise
        except (OSError, UnicodeError) as ex:
            raise ConversionError(L(f"XYZ ファイル {path} を開くか読むことができません",
                                    f"cannot open or read the XYZ file {path}")) from ex
        except (ValueError, TypeError, IndexError, StopIteration) as ex:
            raise ConversionError(L(f"XYZ を {path} から読めません。原子行の元素と座標を確認してください",
                                    f"cannot read XYZ from {path}; check the elements and coordinates in the atom lines")) from ex
    dst.parent.mkdir(parents=True, exist_ok=True)
    try:
        write(dst, combined, format="xyz")
    except (OSError, ValueError) as ex:
        raise ConversionError(L(f"XYZ を {dst} に書けません: {ex}",
                                f"cannot write XYZ to {dst}: {ex}")) from ex
    return dst


def _spec_export_atoms(spec: CalculationSpec):
    # Restore Spec-only constraints and Å/fs velocities for ASE export.
    import numpy as np
    from ase.constraints import FixAtoms, FixCartesian
    from ase.units import fs

    atoms = spec.atoms.copy()
    # Apply velocities before constraints: this is a lossless export, not an
    # MD update that projects momenta onto the unconstrained directions.
    if spec.structure.velocities is not None:
        velocities = np.asarray(spec.structure.velocities, dtype=float)
        if velocities.shape != (len(atoms), 3) or not np.isfinite(velocities).all():
            raise ConversionError(L("初速度は原子ごとに有限の x, y, z 成分 [Å/fs] が必要です",
                                    "initial velocities require finite x, y, z components in Å/fs for each atom"))
        atoms.set_velocities(velocities / fs)
    constraints = []
    if spec.structure.fixed_atoms:
        constraints.append(FixAtoms(indices=spec.structure.fixed_atoms))
    for key, movable in sorted(spec.structure.fixed_axes.items(), key=lambda item: int(item[0])):
        constraints.append(FixCartesian(int(key), mask=[not value for value in movable]))
    atoms.set_constraint(constraints)
    return atoms


def _structure_losses(atoms, fmt: str) -> list[str]:
    # Known omissions in common ASE structure writers; no chemical judgment.
    from ase.constraints import FixAtoms, FixCartesian

    constrained = any(
        bool(len(c.index)) if isinstance(c, FixAtoms) else
        bool(c.mask.any()) if isinstance(c, FixCartesian) else True
        for c in atoms.constraints
    )
    losses = []
    if constrained and fmt in {"cif", "lammps-data", "xyz", "gen", "proteindatabank"}:
        losses.append(L("固定原子・軸固定", "fixed-atom / fixed-axis constraints"))
    if atoms.has("momenta") and fmt in {"cif", "xyz", "gen", "proteindatabank"}:
        losses.append(L("初速度", "initial velocities"))
    return losses


def _structure_format(path: Path, explicit: str | None, *, reading: bool) -> str | None:
    if explicit:
        return explicit
    from ase.io.formats import (
        PEEK_BYTES, UnknownFileTypeError, get_compression, match_magic, open_with_compression,
    )

    filename, _compression = get_compression(str(path))
    name = Path(filename)
    if name.name.upper() in ("POSCAR", "CONTCAR"):
        return "vasp"
    if name.suffix.lower() in (".data", ".lammps"):
        if reading:
            # A RuNNer input.data or LAMMPS dump can have the same suffix.
            # Retain ASE's content-based detection before applying our alias.
            with open_with_compression(str(path), "rb") as stream:
                header = stream.read(PEEK_BYTES)
            for candidate in (header, header.replace(b"\r\n", b"\n")):
                try:
                    return match_magic(candidate).name
                except UnknownFileTypeError:
                    continue
        return "lammps-data"
    return None


def write_spec_structure(spec: CalculationSpec, output: Path | str, *, overwrite: bool = False) -> Path:
    if not str(output).strip():
        raise ConversionError(L("構造の保存先を指定してください", "Choose an output file for the structure"))
    dst = Path(output).expanduser()
    if dst.exists() and not overwrite:
        raise ConversionError(L(f"出力先 {dst} はすでにあります (上書きしません)",
                                f"output {dst} already exists (not overwritten)"))
    try:
        from ase.io import write
        atoms = _spec_export_atoms(spec)
        dst.parent.mkdir(parents=True, exist_ok=True)
        write(dst, atoms, format="extxyz")
    except Exception as ex:
        raise ConversionError(L(f"構造を {dst} に書けません: {ex}", f"cannot write the structure to {dst}: {ex}")) from ex
    return dst


def retarget_spec(source: CalculationSpec, target_conditions: CalculationSpec, *,
                  use_target_thermostat: bool = False, no_velocities: bool = False) -> tuple[CalculationSpec, dict]:
    if source.method.code == target_conditions.method.code:
        raise ConversionError(L(f"変換元と変換先がどちらも {source.method.code} です",
                                f"source and target are both {source.method.code}"))
    # The structure in a Spec is not necessarily the structure actually used
    # by a force-field run backed by an external topology/data file. Copying it
    # into a DFT input would silently discard per-atom charges, bonds and even
    # potentially the actual coordinate ordering.
    if source.method.code == "lammps" and source.method.data_file.strip():
        raise ConversionError(L(
            "変換元の LAMMPS は外部 data ファイルを使います。Spec の構造だけでは、実際に使う座標・原子電荷・結合情報との一致を確認できません。このコード間変換は停止します。data ファイルの内容を点検し、構造と条件を明示した Spec を用意してください。",
            "The source LAMMPS run uses an external data file. The spec alone cannot verify its agreement with the actual coordinates, atomic charges, or topology. Cross-code conversion is stopped; inspect the data file and prepare an explicit spec for the structure and settings."))
    if source.method.code in {"gromacs", "amber", "namd", "openmm"}:
        code_name = {"gromacs": "GROMACS", "amber": "Amber", "namd": "NAMD", "openmm": "OpenMM"}[source.method.code]
        raise ConversionError(L(
            f"変換元の {code_name} は外部の構造とトポロジーを使います。Spec だけでは実際の原子順・電荷・結合との一致を確かめられないため、このコード間変換は停止します。",
            f"The source {code_name} run uses external structure and topology files. The spec alone cannot verify actual atom ordering, charges, or bonds, so this cross-code conversion is stopped."))
    task = source.task
    structure = source.structure
    changes = {}
    omitted = []
    notes = []
    if use_target_thermostat:
        if task.type != "molecular_dynamics" or task.md.ensemble == "NVE":
            raise ConversionError(L(
                "--use-target-thermostat は変換元が NVT または NPT の MD のときだけ使えます",
                "--use-target-thermostat requires an NVT or NPT MD source calculation"))
        if "thermostat" not in target_conditions.task.md.model_fields_set:
            raise ConversionError(L(
                "変換先の雛形に task.md.thermostat を明示してください。熱浴の種類を ADIT が選ぶことはしません",
                "Set task.md.thermostat explicitly in the target template; ADIT does not choose a thermostat"))
        thermostat = target_conditions.task.md.thermostat
        task = task.model_copy(update={"md": task.md.model_copy(update={"thermostat": thermostat})})
        changes["task.md.thermostat"] = {"source": source.task.md.thermostat, "target": thermostat}
        notes.append(L(
            f"明示指定により、熱浴を変換先の雛形から取りました: {source.task.md.thermostat} → {thermostat}。温度・時間刻み・ステップ数・出力間隔・時定数・圧力は変換元の値を保持しています。熱浴の変更は、同じ運動方程式の再現を意味しません。",
            f"As explicitly requested, the thermostat came from the target template: {source.task.md.thermostat} → {thermostat}. Temperature, timestep, step count, output interval, coupling times, and pressure retain their source values. Changing the thermostat does not reproduce the same equations of motion."))
    if no_velocities:
        structure = structure.model_copy(update={"velocities": None})
        if source.structure.velocities is not None:
            omitted.append("structure.velocities")
            notes.append(L(
                "明示指定により初速度を引き継いでいません。同じ初速度から始まる MD の比較にはなりません。変換先コードの初速度の設定を確認してください。",
                "Initial velocities were omitted as explicitly requested. The MD runs will not start with the same velocities; check the target code's velocity initialization settings."))
    comment = source.meta.comment
    target_note = f"converted {source.method.code} -> {target_conditions.method.code}"
    meta = target_conditions.meta.model_copy(update={"comment": (comment + " " + target_note).strip()})
    out = target_conditions.model_copy(update={
        "structure": structure,
        "task": task,
        "runtime": source.runtime,
        "handoff": None,
        "meta": meta,
    })
    from adit.applied import code_specific_unapplied

    not_applied = dict(code_specific_unapplied(out))
    if out.method.code == "vasp":
        from adit.codes.vasp import md_unapplied_settings
        for field, detail in md_unapplied_settings(out).items():
            notes.append(L(
                f"{field} = {detail['value']:g} fs は spec.json に保持していますが、VASP 入力には適用していません。{detail['reason']} コード間変換では、これらの INCAR パラメータは変換先の雛形から取っています。",
                f"{field} = {detail['value']:g} fs is retained in spec.json but is not applied to VASP input. {detail['reason']} For this conversion, these INCAR parameters came from the target template."))
    if source.plumed is not None:
        not_applied["plumed"] = {"value": "plumed.dat",
                                 "reason": L("PLUMED の入力は原子の番号で書くため、原子の並びが変わりうるコード間変換では引き継ぎません。変換先の並びで書き直してください。",
                                             "A PLUMED input refers to atoms by index, so it is not carried across a conversion that may reorder atoms; rewrite it for the target ordering.")}
        notes.append(L(
            "変換元には PLUMED の設定がありましたが、変換先の Spec には入れていません (原子の番号が変わりうるため)。",
            "The source had PLUMED settings; they were not copied to the target spec because atom indices may change."))
    report = {
        "source_code": source.method.code,
        "target_code": out.method.code,
        "preserved": (
            ([f"structure.{key}" for key in type(structure).model_fields
              if (key != "velocities" or structure.velocities is not None) and f"structure.{key}" not in not_applied]
             if omitted or any(k.startswith("structure.") for k in not_applied) else ["structure"])
            + ([f"task.{key}" for key in type(task).model_fields if key != "md" and f"task.{key}" not in not_applied]
               + [f"task.md.{key}" for key in type(task.md).model_fields
                  if not (use_target_thermostat and key == "thermostat") and f"task.md.{key}" not in not_applied]
               if use_target_thermostat or not_applied else ["task"]) + ["runtime"]),
        "taken_from_target_template": ["method", "kpoints"] + (["task.md.thermostat"] if use_target_thermostat else []),
        "changes": changes,
        "omitted": omitted,
        "not_applied_by_target": not_applied,
        "preserved_in_spec_only": list(not_applied),
        "options": {"use_target_thermostat": use_target_thermostat, "no_velocities": no_velocities},
        "not_translated": [f"method ({source.method.code} -> {out.method.code})"],
        "rule": L(
            "構造とコード共通の計算条件は、以下に記した明示変更と変換先で適用されない項目を除き、保持しました。力場、電子状態計算の手法、カットオフ、擬ポテンシャルなど、コード固有の条件は相互変換せず、変換先の雛形の値を使いました。",
            "The structure and code-independent calculation settings were preserved except for the explicit changes and unapplied settings listed below. Code-specific settings such as force fields, electronic-structure methods, cutoffs, and pseudopotentials were not translated; values came from the target template.",
        ) + (" " + " ".join(notes) if notes else ""),
    }
    from adit.conversion_audit import field_audit
    report["field_audit"] = field_audit(source, target_conditions, out, report)
    return out, report


def _read_without_crlf(src: Path, frame: int, fmt: str | None):
    import tempfile

    from ase.io import read

    raw = src.read_bytes()
    if b"\r\n" not in raw:
        raise
    with tempfile.TemporaryDirectory() as tmp:
        copy = Path(tmp) / src.name
        copy.write_bytes(raw.replace(b"\r\n", b"\n"))
        return read(copy, index=frame, format=fmt)


def convert_structure(source: Path | str, output: Path | str, *, input_format: str | None = None,
                      output_format: str | None = None, frame: int = -1, cell: tuple[float, float, float] | None = None,
                      overwrite: bool = False) -> Path:
    if not str(source).strip():
        raise ConversionError(L("変換元の構造ファイルを指定してください", "Choose a source structure file"))
    if not str(output).strip():
        raise ConversionError(L("変換後の構造ファイルを指定してください", "Choose an output structure file"))
    src, dst = Path(source).expanduser(), Path(output).expanduser()
    if dst.exists() and not overwrite:
        raise ConversionError(L(f"出力先 {dst} はすでにあります (上書きしません)",
                                f"output {dst} already exists (not overwritten)"))
    spec_file = src / "spec.json" if src.is_dir() else src
    try:
        if spec_file.name == "spec.json":
            atoms = _spec_export_atoms(CalculationSpec.load(spec_file))
        else:
            from ase.io import read
            fmt = _structure_format(src, input_format, reading=True)
            try:
                atoms = read(src, index=frame, format=fmt)
            except Exception:
                atoms = _read_without_crlf(src, frame, fmt)
    except Exception as ex:
        raise ConversionError(L(f"構造を {src} から読めません: {ex}", f"cannot read a structure from {src}: {ex}")) from ex
    if cell is not None:
        if any(not math.isfinite(x) or x <= 0 for x in cell):
            raise ConversionError(L("セルの長さには有限で 0 より大きい値を指定してください",
                                    "cell lengths must be finite and greater than zero"))
        atoms.set_cell(cell)
        atoms.center()
        atoms.pbc = True
    try:
        from ase.io.formats import filetype
        fmt = _structure_format(dst, output_format, reading=False) or filetype(dst, read=False)
    except Exception as ex:
        raise ConversionError(L(f"出力形式を判定できません: {ex}", f"cannot determine the output format: {ex}")) from ex
    if fmt == "espresso-in":
        raise ConversionError(L(
            "Quantum ESPRESSO の入力には元素ごとの擬ポテンシャルと計算条件が必要です。構造形式の変換では指定できません。QE の設定を含む spec.json を adit-gen へ渡すか、変換先の雛形を使って adit-convert calculation を実行してください。",
            "Quantum ESPRESSO input requires pseudopotentials for each element and calculation settings, which structure conversion cannot supply. Use adit-gen with a QE spec.json, or adit-convert calculation with a target template."))
    if fmt == "vasp" and atoms.cell.rank < 3:
        raise ConversionError(L("POSCAR には 3 方向のセル (格子ベクトル) が必要ですが、入力構造にセルがありません。--cell A または --cell A,B,C でセルの長さ [Å] を指定してください",
                                "POSCAR needs a three-dimensional cell, but the input has none. Give cell lengths in Å with --cell A or --cell A,B,C"))
    try:
        from ase.io import write
        if fmt == "vasp":
            from adit.vasp_constraints import prepare_vasp_constraints
            atoms = prepare_vasp_constraints(atoms)
        dst.parent.mkdir(parents=True, exist_ok=True)
        # Without Masses, ASE reads type 1/2 as H/He rather than the source
        # elements. This is a structure export, not a molecular topology.
        options = {"masses": True, "atom_style": "atomic", "velocities": atoms.has("momenta")} if fmt == "lammps-data" else {}
        write(dst, atoms, format=fmt, **options)
    except Exception as ex:
        raise ConversionError(L(f"構造を {dst} に書けません: {ex}", f"cannot write the structure to {dst}: {ex}")) from ex
    losses = _structure_losses(atoms, fmt)
    if losses:
        omitted = ", ".join(losses)
        warnings.warn(L(
            f"{fmt} 出力には {omitted} を保存していません。これらを保持するには extended XYZ (.extxyz) を使ってください。",
            f"The {fmt} output does not retain {omitted}. Use extended XYZ (.extxyz) to preserve these data."),
            StructureDataLossWarning, stacklevel=2)
    return dst


def _parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="adit-convert", description=L(
        "構造形式を変換するか、共通条件を保って別の計算コードの入力一式を生成します",
        "convert a structure format, or generate input for another code while preserving shared settings"))
    sub = ap.add_subparsers(dest="command", required=True)
    st = sub.add_parser("structure", help=L("構造ファイルの形式を変換します", "convert a structure file format"))
    st.add_argument("source"); st.add_argument("output")
    st.add_argument("--input-format", help=L("入力形式を明示します (通常は拡張子から判定)", "input format (normally inferred from the filename)"))
    st.add_argument("--output-format", help=L("出力形式を明示します (通常は拡張子から判定)", "output format (normally inferred from the filename)"))
    st.add_argument("--frame", type=int, default=-1, help=L("軌跡から読むフレーム番号 (既定: 最後)", "trajectory frame index (default: last)"))
    st.add_argument("--cell", help=L("セルのない構造に与えるセルの長さ [Å]: A または A,B,C", "cell lengths in Å for a structure without a cell: A or A,B,C"))
    st.add_argument("--overwrite", action="store_true", help=L("既存の出力を上書きします", "overwrite an existing output"))
    crystal = sub.add_parser("crystal", help=L("空間群と分率座標から結晶を作り、構造ファイルに書きます",
                                               "build a crystal from a space group and fractional coordinates, and write a structure file"))
    crystal.add_argument("spec", help=L("「空間群 元素:x,y,z … cell=a[,b,c[,α,β,γ]]」 (例: 225 Na:0,0,0 Cl:0.5,0,0 cell=5.64)",
                                        "'group element:x,y,z ... cell=a[,b,c[,alpha,beta,gamma]]' (e.g. 225 Na:0,0,0 Cl:0.5,0,0 cell=5.64)"))
    crystal.add_argument("output", help=L("書き出す構造ファイル (拡張子で形式が決まります: .cif / .xyz / POSCAR など)",
                                          "structure file to write (the format follows the extension: .cif / .xyz / POSCAR ...)"))
    crystal.add_argument("--overwrite", action="store_true", help=L("すでにあるファイルに上書きします", "overwrite an existing file"))
    surface = sub.add_parser("surface", help=L("表面スラブを作り、構造ファイルに書きます",
                                               "build a surface slab and write a structure file"))
    surface.add_argument("spec", help=L("「面 元素 nx×ny×層数 [vacuum=Å] [a=Å]」 (例: fcc111 Al 2x2x3 vacuum=10)",
                                        "'facet element nx x ny x layers [vacuum=A] [a=A]' (e.g. fcc111 Al 2x2x3 vacuum=10)"))
    surface.add_argument("output", help=L("書き出す構造ファイル (拡張子で形式が決まります)",
                                          "structure file to write (the format follows the extension)"))
    surface.add_argument("--overwrite", action="store_true", help=L("すでにあるファイルに上書きします", "overwrite an existing file"))
    merge = sub.add_parser("combine", help=L("単一フレームの XYZ を座標を動かさずに結合します",
                                            "combine single-frame XYZ files without moving coordinates"))
    merge.add_argument("paths", nargs="+", help=L("入力 XYZ を 2 つ以上、最後に出力 XYZ",
                                                  "two or more input XYZ files, followed by the output XYZ"))
    merge.add_argument("--overwrite", action="store_true", help=L("既存の出力を上書きします", "overwrite an existing output"))
    babel = sub.add_parser("openbabel", help=L("Open Babel で構造形式を変換します (結合次数などの解釈は Open Babel に従います)",
                                                "convert a structure with Open Babel (bond typing and interpretation follow Open Babel)"))
    babel.add_argument("source"); babel.add_argument("output")
    babel.add_argument("--input-format", help=L("Open Babel の入力形式名", "Open Babel input format name"))
    babel.add_argument("--output-format", help=L("Open Babel の出力形式名", "Open Babel output format name"))
    babel.add_argument("--overwrite", action="store_true", help=L("既存の出力を上書きします", "overwrite an existing output"))
    dock = sub.add_parser("dock6", help=L("利用者が作成した DOCK6 の dock.in と入力ファイルをまとめます",
                                             "package a user-authored DOCK6 dock.in and its input files"))
    dock.add_argument("source", help=L("既存の dock.in", "existing dock.in"))
    dock.add_argument("output", help=L("空の出力ディレクトリ", "empty output directory"))
    dock.add_argument("--asset", action="append", default=[], help=L("追加で複製する、dock.in と同じ作業ディレクトリ内の相対パス",
                                                                "additional relative path inside the dock.in working directory to copy"))
    native = sub.add_parser("import", help=L("既存の計算入力を読み取り、出典付きの点検記録と下書きの計算設定 (draft_spec.json) を作ります",
                                              "inspect existing native input and write a sourced report and draft spec"))
    native.add_argument("source", help=L("既存の入力ディレクトリまたは主入力ファイル",
                                          "existing input directory or primary input file"))
    native.add_argument("output", help=L("新しい点検記録の出力先ディレクトリ", "new audit-report directory"))
    native.add_argument("--code", choices=("vasp", "espresso", "qe", "lammps", "gromacs"),
                        help=L("自動判定できないときの入力コード", "input code when it cannot be detected"))
    verify = sub.add_parser("verify", help=L("生成済み入力のうち対応付けた条件を spec.json と読み戻し照合します",
                                              "re-import mapped native settings and compare them with spec.json"))
    verify.add_argument("project", help=L("spec.json を含む生成ディレクトリ", "generated directory containing spec.json"))
    verify.add_argument("--report", help=L("新しい JSON 点検記録ファイル (省略時は書き込みなし)",
                                            "new JSON audit-report file (nothing is written if omitted)"))
    calc = sub.add_parser("calculation", help=L("別の計算コードの入力一式を生成します", "generate input for another calculation code"))
    calc.add_argument("source", help=L("変換元の spec.json または生成ディレクトリ", "source spec.json or generated directory"))
    calc.add_argument("target_template", help=L("変換先コードの研究室の雛形名または雛形ファイル", "target-code template name or file"))
    calc.add_argument("output")
    calc.add_argument("--config", help=L("cluster.toml の場所", "path to cluster.toml"))
    calc.add_argument("--use-target-thermostat", action="store_true", help=L(
        "NVT/NPT の熱浴だけを変換先の雛形から取ります (温度や時定数は変換元を保持)",
        "use the target template's thermostat for NVT/NPT (retain source temperature and coupling times)"))
    calc.add_argument("--no-velocities", action="store_true", help=L(
        "初速度を明示的に除外し、その事実を変換記録に残します",
        "omit initial velocities and record the omission in the conversion report"))
    calc.add_argument("--accept-import-defaults", action="store_true", help=L(
        "読み込んだ下書きの計算設定 (draft_spec.json) の元入力にない既定値を確認済みとして、変換を続けます",
        "continue after reviewing defaults absent from an imported draft spec"))
    calc.add_argument("--overwrite", action="store_true", help=L("空でない出力ディレクトリを上書きします", "overwrite a non-empty output directory"))
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _parser().parse_args(argv)
    if args.command in ("crystal", "surface"):
        from ase.io import write as ase_write

        from adit.structure import StructureError, from_spacegroup, from_surface

        out = Path(args.output).expanduser()
        if out.exists() and not args.overwrite:
            print(L(f"すでにあります: {out} (上書きしてよければ --overwrite を付けてください)",
                    f"already exists: {out} (add --overwrite to replace it)"), file=sys.stderr)
            return 1
        try:
            atoms = from_spacegroup(args.spec) if args.command == "crystal" else from_surface(args.spec)
        except StructureError as ex:
            print(str(ex), file=sys.stderr)
            return 1
        try:
            ase_write(out, atoms)
        except Exception as ex:
            print(L(f"{out} に書けません: {ex}", f"cannot write {out}: {ex}"), file=sys.stderr)
            return 1
        print(L(f"{atoms.get_chemical_formula()} ({len(atoms)} 原子) を {out} に書きました",
                f"wrote {atoms.get_chemical_formula()} ({len(atoms)} atoms) to {out}"))
        return 0
    if args.command == "verify":
        from adit.native_verify import verify_generated_bundle
        from adit.native_workflow import verification_passed
        result = verify_generated_bundle(args.project)
        checked = verification_passed(result)
        if args.report:
            dst = Path(args.report).expanduser().resolve()
            project = Path(args.project).expanduser().resolve()
            if dst.exists() or dst.is_relative_to(project):
                print(L(f"点検記録は既存ファイルや生成ディレクトリの中に書けません: {dst}",
                        f"the audit report cannot overwrite an existing file or be written inside the generated directory: {dst}"), file=sys.stderr)
                return 1
            try:
                dst.parent.mkdir(parents=True, exist_ok=True)
                dst.write_text(json.dumps(result.report(), ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            except OSError as ex:
                print(str(ex), file=sys.stderr)
                return 1
        print(L(
            f"読み戻して照合した項目: {len(result.preserved)} 件一致、{len(result.mismatched)} 件不一致、{len(result.unverifiable)} 件未確認。計算全体の同等性は判定していません。",
            f"Re-imported mapped fields: {len(result.preserved)} matched, {len(result.mismatched)} mismatched, {len(result.unverifiable)} unverified. Whole-calculation equivalence was not assessed."))
        print(L(
            f"対応項目のみの点検結果: {'通過' if checked else '不通過'}。未確認項目は成功終了でも検証済みを意味しません。",
            f"Mapped-field check only: {'passed' if checked else 'failed'}. Unverified fields remain unverified even when this command exits successfully."))
        for item in result.mismatched[:5]:
            print(L(f"  不一致: {item['field']}", f"  mismatch: {item['field']}"), file=sys.stderr)
        if not result.mismatched:
            for item in result.unverifiable[:5]:
                reason = item.get("reason") or item.get("detail") or ""
                print(L(f"  未確認: {item.get('field', '?')}{': ' + str(reason) if reason else ''}",
                        f"  not verified: {item.get('field', '?')}{': ' + str(reason) if reason else ''}"), file=sys.stderr)
            codes = "VASP / Quantum ESPRESSO / LAMMPS / GROMACS"
            print(L(f"  (読み戻して照合できるのは {codes} の入力だけです。ほかのコードでは必ず「未確認」になり、"
                    "それは計算に問題があることを意味しません)",
                    f"  (only {codes} inputs can be read back; for other codes the result is always 'not verified', "
                    "which does not mean anything is wrong with the run)"), file=sys.stderr)
        if not checked:
            for item in result.unverifiable[:3]:
                print(L(f"  未確認: {item['field']}", f"  not verified: {item['field']}"), file=sys.stderr)
        if args.report:
            print(L(f"点検記録: {dst}", f"audit report: {dst}"))
        return 0 if checked else 1
    if args.command == "import":
        from adit.native_workflow import NativeWorkflowError, write_import_review
        src = Path(args.source).expanduser().resolve()
        try:
            result, dst = write_import_review(src, args.output, code=args.code)
        except (NativeWorkflowError, OSError, ValueError) as ex:
            print(str(ex), file=sys.stderr)
            return 1
        if result.spec is None:
            print(L(f"読み取れない条件があります。理由を {dst / 'import_report.json'} に記録しました。下書きの計算設定 (draft_spec.json) は作っていません。",
                    f"Some settings could not be imported. Reasons are in {dst / 'import_report.json'}; no draft spec was written."), file=sys.stderr)
            issues = [*result.unknown, *result.unsupported]
            for issue in issues[:3]:
                line_ref = f"{issue.get('file', src)}:{issue.get('line', '?')}"
                source_text = str(issue.get("text", "")).strip().replace("\n", " ")[:120]
                detail = source_text or L("入力を確認してください", "check this input")
                print(f"  {line_ref}: {detail}", file=sys.stderr)
            if len(issues) > 3:
                print(L(f"  ほか {len(issues) - 3} 件は import_report.json を参照してください。",
                        f"  See import_report.json for {len(issues) - 3} more issue(s)."), file=sys.stderr)
            return 1
        print(L(f"点検記録と下書きの計算設定 (draft_spec.json) を {dst} に作りました。既定値と外部パラメータを確認してから使ってください。",
                f"Wrote an audit report and draft spec in {dst}. Review defaults and external parameters before using it."))
        return 0
    if args.command == "dock6":
        from adit.dock6 import Dock6Error, package_dock6
        try:
            p = package_dock6(args.source, args.output, assets=args.asset)
        except (Dock6Error, OSError) as ex:
            print(str(ex), file=sys.stderr)
            return 1
        print(L(f"DOCK6 の入力を {p.resolve()} にまとめました。計算条件と入力ファイルを確認してください。",
                f"packaged DOCK6 input in {p.resolve()}. Review the calculation settings and input files."))
        return 0
    if args.command == "openbabel":
        try:
            p = convert_with_openbabel(args.source, args.output, input_format=args.input_format,
                                       output_format=args.output_format, overwrite=args.overwrite)
        except ConversionError as ex:
            print(str(ex), file=sys.stderr)
            return 1
        print(L(f"Open Babel で構造を書きました: {p.resolve()}。結合次数・電荷・周期セルなどが保持されたか確認してください。",
                f"wrote the structure with Open Babel: {p.resolve()}. Check whether bond orders, charges, and periodic cells were preserved."))
        return 0
    if args.command == "combine":
        try:
            if len(args.paths) < 3:
                raise ConversionError(L("入力 XYZ を 2 つ以上、そのあとに出力 XYZ を指定してください",
                                        "give at least two input XYZ files followed by the output XYZ"))
            p = combine_xyz(args.paths[:-1], args.paths[-1], overwrite=args.overwrite)
        except ConversionError as ex:
            print(str(ex), file=sys.stderr)
            return 1
        print(L(f"XYZ を書きました: {p.resolve()}", f"wrote XYZ: {p.resolve()}"))
        return 0
    if args.command == "structure":
        try:
            cell = None
            if args.cell:
                try:
                    values = [float(x.strip()) for x in args.cell.split(",")]
                except ValueError as ex:
                    raise ConversionError(L("--cell には数値を指定してください", "--cell requires numeric values")) from ex
                if len(values) == 1:
                    values *= 3
                if len(values) != 3:
                    raise ConversionError(L("--cell は A または A,B,C の形で指定してください", "--cell must be A or A,B,C"))
                cell = tuple(values)
            with warnings.catch_warnings(record=True) as notices:
                warnings.simplefilter("always", StructureDataLossWarning)
                p = convert_structure(args.source, args.output, input_format=args.input_format,
                                      output_format=args.output_format, frame=args.frame, cell=cell, overwrite=args.overwrite)
            for notice in notices:
                print(str(notice.message), file=sys.stderr)
        except (ConversionError, ValueError) as ex:
            print(str(ex), file=sys.stderr); return 1
        print(L(f"構造を書きました: {p.resolve()}", f"wrote the structure: {p.resolve()}")); return 0
    try:
        cfg, _path, _created = ensure_config(Path(args.config).expanduser() if args.config else config_path())
        lang.set_language(env_var("LANG", cfg.language))
        src_path = Path(args.source).expanduser()
        source_file = src_path / "spec.json" if src_path.is_dir() else src_path
        if source_file.name == "draft_spec.json" and not args.accept_import_defaults:
            raise ConversionError(L(
                "読み込んだ draft_spec.json には元入力に無い既定値があります。隣の import_report.json と外部パラメータを確認し、続ける場合だけ --accept-import-defaults を指定してください。",
                "An imported draft_spec.json contains defaults absent from the native input. Review the adjacent import_report.json and external parameters; pass --accept-import-defaults only if you want to continue."))
        source = CalculationSpec.load(source_file)
        target = load_template(args.target_template, source, cfg)
        converted, report = retarget_spec(source, target, use_target_thermostat=args.use_target_thermostat,
                                           no_velocities=args.no_velocities)
        report["input_sources"] = {
            "source_spec": {"path": str(source_file.resolve()),
                            "sha256": sha256_file(source_file)},
            "target_template_reference": args.target_template,
            "accepted_import_defaults": bool(args.accept_import_defaults),
        }
        import_report = source_file.with_name("import_report.json") if source_file.name == "draft_spec.json" else None
        if import_report is not None and import_report.is_file():
            report["input_sources"]["import_report"] = {
                "path": str(import_report.resolve()),
                "sha256": sha256_file(import_report),
            }
        report_text = json.dumps(report, ensure_ascii=False, indent=2) + "\n"
        extra = [L("== 計算コード間の変換 ==", "== Conversion between calculation codes =="),
                 f"  {source.method.code} -> {converted.method.code}", f"  {report['rule']}",
                 L(f"  詳細: {REPORT_FILE}", f"  Details: {REPORT_FILE}"), ""]
        written = write_project(converted, cfg, args.output, overwrite=args.overwrite,
                                extra_readme=extra, extra_texts={REPORT_FILE: report_text})
    except (ConversionError, ConfigError, TemplateError, ProjectError, OSError, ValueError) as ex:
        print(str(ex), file=sys.stderr)
        if isinstance(ex, ProjectError):
            if any(e.location == "structure.velocities" for e in ex.errors):
                print(L("初速度を引き継がずに生成する場合だけ --no-velocities を指定してください。",
                        "Use --no-velocities only if you intend to generate input without the source velocities."), file=sys.stderr)
            if any(e.location == "task.md.thermostat" for e in ex.errors) and not args.use_target_thermostat:
                print(L("熱浴を変更する場合は、変換先の雛形に task.md.thermostat を指定し、--use-target-thermostat を付けてください。",
                        "To change the thermostat, set task.md.thermostat in the target template and pass --use-target-thermostat."), file=sys.stderr)
        return 1
    print(L(f"{converted.method.code} の入力を {Path(args.output).resolve()} に生成しました ({len(written)} ファイル)",
            f"generated {converted.method.code} input in {Path(args.output).resolve()} ({len(written)} files)"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
