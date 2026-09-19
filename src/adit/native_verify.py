"""Read-only comparison of mapped native values against a generated spec."""
# This checks represented fields, never whole-calculation equivalence.
# VASP element grouping follows codes/vasp.py and https://vasp.at/wiki/POSCAR.
# LAMMPS time units: https://docs.lammps.org/units.html (converted by native_import).
# QE input units: https://www.quantum-espresso.org/Doc/INPUT_PW.html.
from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
import math
from pathlib import Path

from adit.lang import L
from adit.native_import import MAX_ATOMS, MAX_BYTES, NativeImportError, import_native
from adit.spec import CalculationSpec


@dataclass
class NativeVerificationResult:
    preserved: list[dict] = field(default_factory=list)
    mismatched: list[dict] = field(default_factory=list)
    unverifiable: list[dict] = field(default_factory=list)
    source: dict = field(default_factory=dict)
    native_report: dict = field(default_factory=dict)

    def report(self) -> dict:
        return {"preserved": self.preserved, "mismatched": self.mismatched,
                "unverifiable": self.unverifiable, "source": self.source,
                "native_import": self.native_report,
                "scope": L("元の spec.json と入力に明示された対応項目だけの比較です。計算全体の同等性や実行可能性は判定しません。",
                           "Only explicitly mapped native fields are compared with the source spec. Whole-calculation equivalence and executability are not assessed."),
                "tolerances": {"structure_angstrom_absolute": 1e-8,
                               "numeric_relative": 1e-8, "numeric_absolute": 1e-14}}


def _same(a, b, *, structure=False):
    if isinstance(a, bool) or isinstance(b, bool):
        return type(a) is type(b) and a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not math.isfinite(a) or not math.isfinite(b):
            return False
        if isinstance(a, int) and isinstance(b, int):
            return a == b
        return math.isclose(a, b, rel_tol=0 if structure else 1e-8, abs_tol=1e-8 if structure else 1e-14)
    if isinstance(a, (tuple, list)) and isinstance(b, (tuple, list)):
        return len(a) == len(b) and all(_same(x, y, structure=structure) for x, y in zip(a, b))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(_same(a[k], b[k], structure=structure) for k in a)
    return type(a) is type(b) and a == b


def _get(data, key):
    for part in key.split("."):
        data = data[int(part)] if isinstance(data, list) else data[part]
    return data


def verify_generated_bundle(project_dir: Path | str) -> NativeVerificationResult:
    """Compare safe mapped fields, without modifying or executing the bundle."""
    result = NativeVerificationResult()
    root = Path(project_dir).expanduser().resolve()
    path = root / "spec.json"
    try:
        if not path.is_file() or path.stat().st_size > MAX_BYTES:
            raise NativeImportError(L("spec.json がないか、読み込み上限を超えています",
                                      "spec.json is missing or exceeds the import size limit"))
        with open(path, "rb") as stream:
            raw = stream.read(MAX_BYTES + 1)
        if len(raw) > MAX_BYTES:
            raise NativeImportError(L("spec.json が読み込み上限を超えています",
                                      "spec.json exceeds the import size limit"))
        source = CalculationSpec.from_json(raw.decode("utf-8-sig"))
        if not 0 < len(source.structure.atoms.symbols) <= MAX_ATOMS:
            raise NativeImportError(L("元の spec.json が原子数の上限を超えています",
                                      "the source spec exceeds the atom-count limit"))
        result.source = {"file": str(path), "sha256": hashlib.sha256(raw).hexdigest(), "code": source.method.code}
        native = import_native(root, code=source.method.code)
    except (OSError, UnicodeError, ValueError, KeyError, TypeError) as exc:
        result.unverifiable.append({"field": "bundle", "reason": L(
            "比較用の spec.json または計算入力を読めません。ファイルと形式を確認してください。",
            "Cannot read spec.json or the calculation input for verification. Check the files and their formats."),
            "detail": str(exc)})
        return result
    result.native_report = native.report()
    expected = source.model_dump(mode="json")
    compared = set()

    def compare(key, a, b, provenance, *, structure=False):
        item = {"field": key, "expected": a, "observed": b, "provenance": provenance}
        (result.preserved if _same(a, b, structure=structure) else result.mismatched).append(item)
        compared.add(key)

    # Method and task fields are compared only when the parser recorded a
    # direct mapping; inferred defaults are deliberately excluded.
    for key, provenance in native.provenance.items():
        if not key.startswith(("method.", "task.")):
            continue
        if key == "method.pseudo" and not source.method.pseudo:
            result.unverifiable.append({"field": key, "reason": L(
                "元の spec.json に擬ポテンシャルの指定がなく、生成時の設定から補われています。ファイル名を元の spec.json と照合できません。",
                "The source spec has no pseudopotential mapping; generation supplied it from configuration, so its filenames cannot be checked against the source spec."),
                "provenance": provenance})
            continue
        try:
            a = _get(expected, key)
            if key == "method.kpoints_centering" and source.kpoints is not None and source.kpoints.mode == "gamma":
                a = "gamma"
        except (KeyError, IndexError, TypeError):
            result.unverifiable.append({"field": key, "reason": L("元の spec.json に対応する明示項目がありません", "no matching field is present in source spec"), "provenance": provenance})
            continue
        compare(key, a, native.parsed[key], provenance)

    if native.spec is not None:
        actual = native.spec
        if "task.type" not in compared:
            task_provenance = (native.provenance.get("native.control.calculation") or
                               native.provenance.get("native.NSW") or native.provenance.get("native.run"))
            if task_provenance:
                if source.task.type in {"single_point", "geometry_optimization", "molecular_dynamics"}:
                    compare("task.type", source.task.type, actual.task.type, task_provenance)
                else:
                    result.unverifiable.append({"field": "task.type", "reason": L(
                        "この入力は複数段階の計算の一部です。前段階の計算種別と元の spec.json の全体タスクは直接照合しません。",
                        "This input is one stage of a multi-stage calculation. Its stage task is not compared directly with the source spec's overall task."),
                        "provenance": task_provenance})
        kp_source = source.kpoints
        if kp_source is not None and actual.kpoints is not None:
            # Compare actual mesh and fractional shift, not gamma/mesh/density
            # spelling. Density resolution uses the recorded source cell.
            try:
                mesh = list(kp_source.resolved_mesh(source.structure.atoms.cell))
                compare("kpoints.mesh", mesh, list(actual.kpoints.resolved_mesh(actual.structure.atoms.cell)),
                        native.provenance.get("kpoints.mesh", native.provenance.get("kpoints", {})))
                shift = list(kp_source.shift) if kp_source.mode != "gamma" else [0., 0., 0.]
                compare("kpoints.shift", shift, list(actual.kpoints.shift),
                        native.provenance.get("kpoints", native.provenance.get("kpoints.mesh", {})))
            except (ValueError, TypeError, ArithmeticError) as exc:
                result.unverifiable.append({"field": "kpoints", "reason": L(
                    "k 点のメッシュまたはシフトを照合できません。",
                    "The k-point mesh or shift could not be verified."), "detail": str(exc)})
        # Do not claim coordinates stored in an external data file are verified
        # by possibly stale coordinates embedded in spec.json.
        external = source.method.code == "lammps" and bool(source.method.data_file)
        if external:
            result.unverifiable.append({"field": "structure.atoms", "reason": L("元の spec.json は外部の data ファイルを参照します。埋め込み座標との一致を保証できません。", "source spec references external data; embedded coordinates cannot verify that file")})
        else:
            atoms = source.structure.atoms
            order = list(range(len(atoms.symbols)))
            if source.method.code == "vasp":
                order = [i for element in source.elements for i, symbol in enumerate(atoms.symbols) if symbol == element]
            cell = atoms.cell
            rotated = source.method.code == "lammps" and any(abs(cell[i][j]) > 1e-12 for i in range(3) for j in range(3) if i != j)
            if rotated:
                result.unverifiable.append({"field": "structure.atoms", "reason": L("LAMMPS のセル回転を伴う原子対応はこの検査の対象外です。", "LAMMPS cell rotation / atom correspondence is outside this check")})
            else:
                for name, a in (("symbols", [atoms.symbols[i] for i in order]), ("positions", [list(atoms.positions[i]) for i in order]), ("cell", [list(v) for v in cell]), ("pbc", list(atoms.pbc))):
                    key = "structure.atoms." + name
                    provenance = native.provenance.get(key)
                    if provenance is None:
                        card = "ATOMIC_POSITIONS" if name in {"symbols", "positions"} else "CELL_PARAMETERS"
                        provenance = native.provenance.get("native." + card)
                    if provenance is None and source.method.code == "lammps":
                        provenance = {"file": actual.structure.source_ref, "basis": "atomic data with validated explicit IDs and type labels"}
                    b = actual.structure.atoms.model_dump(mode="json")[name]
                    compare(key, a, b, provenance or {}, structure=True)
    else:
        result.unverifiable.append({"field": "structure / task", "reason": L("未対応の指定があるため入力全体の読み込みは成立していません。", "unsupported input prevents complete parsing")})

    for issue in native.unknown + native.unsupported:
        result.unverifiable.append({"field": "native input", **issue})
    for item in native.not_applied:
        result.unverifiable.append(item)
    for dependency in native.unresolved_dependencies:
        result.unverifiable.append({"field": "dependency", **dependency})
    # These fields cannot be established from the native subset, regardless
    # of whether the corresponding source values happen to equal defaults.
    result.unverifiable.append({"field": "unmapped fields", "fields": sorted(set(native.inferred_defaults) - compared),
                               "reason": L("元入力に対応のない既定値は検証していません。", "defaults without explicit native mappings are not verified")})
    return result
