"""Read-only, deliberately limited native-input inspection."""
# Keyword sources: https://vasp.at/wiki/INCAR, https://vasp.at/wiki/KPOINTS,
# https://vasp.at/wiki/POSCAR, https://www.quantum-espresso.org/Doc/INPUT_PW.html,
# https://docs.lammps.org/read_data.html, https://docs.lammps.org/run.html.
# Benign output labels: https://vasp.at/wiki/SYSTEM; QE CONTROL variables are
# defined in INPUT_PW above. Output-directory references are reported, not opened.
# No input text is executed and no referenced potential or script is opened.
from __future__ import annotations

from adit.errors import AditValueError
from dataclasses import dataclass, field
import hashlib
from io import StringIO
import math
from pathlib import Path
import re
import shlex

from adit.lang import L
from adit.spec import AtomsData, CalculationSpec, Structure

MAX_BYTES = 8 * 1024 * 1024
MAX_ATOMS = 10000


class NativeImportError(AditValueError):
    pass


_ISSUE_REASONS_JA = {
    "unsupported INCAR syntax": "この INCAR 行の書式には対応していません",
    "only single-point NSW=0 / IBRION=-1 is supported": "NSW=0・IBRION=-1 の一点計算のみ読み込めます",
    "unmapped INCAR keyword": "この INCAR キーワードには対応していません",
    "Selective dynamics import is unsupported": "Selective dynamics の読み込みには対応していません",
    "velocities / trailing sections are unsupported": "速度または POSCAR 末尾の追加セクションには対応していません",
    "only automatic Gamma / Monkhorst-Pack meshes are supported": "自動生成の Gamma または Monkhorst-Pack メッシュのみ読み込めます",
    "only positive meshes without shifts or trailing sections are supported": "正のメッシュで、シフトと末尾の追加セクションがない場合のみ読み込めます",
    "unsupported namelist syntax": "この namelist 行の書式には対応していません",
    "duplicate assignment": "同じ項目が重複して指定されています",
    "output flags require logical values": "出力フラグには論理値が必要です",
    "unmapped namelist variable": "この namelist 変数には対応していません",
    "only scf with explicit ibrav=0 is supported": "ibrav=0 を明示した SCF 計算のみ読み込めます",
    "duplicate card": "同じカードが重複しています",
    "explicit angstrom cell and angstrom / crystal positions required": "セルは angstrom、座標は angstrom または crystal で単位を明示してください",
    "CELL_PARAMETERS must use angstrom": "CELL_PARAMETERS は angstrom 単位で指定してください",
    "unique plain element labels are required": "重複のない元素記号を指定してください",
    "only gamma / automatic k points supported": "gamma または automatic の k 点のみ読み込めます",
    "unmapped QE card": "この QE カードには対応していません",
    "explicit k points required": "K_POINTS を明示してください",
    "MD requires explicit generated velocity create; restart / supplied velocities are unsupported": "MD では生成器と同じ velocity create の明示指定が必要です。restart や外部から与えた速度には対応していません",
    "temperature ramps / unequal initial and target temperatures are unsupported": "温度ランプや初期温度と目標温度が異なる設定には対応していません",
    "only generated NVE / Nose-Hoover NVT without extra fix options supported": "生成器と同じ NVE または追加オプションのない Nose-Hoover NVT のみ読み込めます",
    "variables / continuation are unsupported": "変数展開と行継続には対応していません",
    "duplicate command": "同じコマンドが重複しています",
    "unmapped LAMMPS command": "この LAMMPS コマンドには対応していません",
    "explicit supported task, atom_style atomic and units metal / real required": "対応する計算タスク、atom_style atomic、units metal または real を明示してください",
    "only analytic lj/cut supported; potential dependencies are not inferred": "解析的な lj/cut のみ読み込めます。ポテンシャルへの外部依存は推測しません",
    "explicit pair_coeff required": "pair_coeff を明示してください",
    "only numeric lj/cut coefficients supported": "数値で指定した lj/cut の係数のみ読み込めます",
    "units, atom_style and boundary must precede read_data": "units・atom_style・boundary は read_data より前に指定してください",
    "explicit p / f boundary flags required": "境界条件には p または f を各方向に明示してください",
    "one bundle-relative data filename required": "入力一式の内部を指す data ファイル名を 1 つ指定してください",
    "data symlink leaves input bundle": "data ファイルのシンボリックリンクが入力一式の外を指しています",
    "duplicate atom type in generated title": "生成時の見出しで原子タイプが重複しています",
    "unsupported LAMMPS data section": "この LAMMPS data セクションには対応していません",
    "only five-column atomic data rows supported": "原子行は 5 列の atomic 形式だけ読み込めます",
    "invalid atom type": "原子タイプが不正です",
    "mass must be positive": "原子質量は正の値にしてください",
    "duplicate mass / atom type": "原子タイプの質量指定が重複しています",
    "Masses label disagrees with generated type title": "Masses の元素ラベルが生成時の見出しと一致しません",
    "Masses rows require explicit '# Element' labels": "Masses の各行に「# 元素記号」のラベルを明示してください",
    "atom IDs must cover 1..N exactly": "原子 ID は 1 から N までを重複なく指定してください",
    "all atom types require explicit element labels": "全原子タイプに元素ラベルを明示してください",
    "commands differ from the supported generated sequence; options / ordering / initialization cannot be preserved": "コマンドの内容・順序・初期化が対応する生成入力と異なります。そのまま保持できないため読み込みません",
}


@dataclass
class NativeImportResult:
    code: str
    spec: CalculationSpec | None = None
    provenance: dict = field(default_factory=dict)
    unknown: list[dict] = field(default_factory=list)
    unsupported: list[dict] = field(default_factory=list)
    files: dict = field(default_factory=dict)
    parsed: dict = field(default_factory=dict)
    inferred_defaults: dict = field(default_factory=dict)
    not_applied: list[dict] = field(default_factory=list)
    unresolved_dependencies: list[dict] = field(default_factory=list)

    def report(self) -> dict:
        return {"code": self.code, "complete": self.spec is not None,
                "ready_to_regenerate": False,
                "provenance": self.provenance, "unknown": self.unknown,
                "unsupported": self.unsupported, "files": self.files,
                "parsed": self.parsed, "inferred_defaults": self.inferred_defaults,
                "not_applied": self.not_applied,
                "unresolved_dependencies": self.unresolved_dependencies,
                "note": L("読み込みは再生成の同一性や実行可能性を保証しません。既定値とパラメータを確認してください。",
                          "Import does not guarantee identical regeneration or executability. Review defaults and parameter files.")}

    def issue(self, path, line, text, reason, *, unknown=False):
        if reason in _ISSUE_REASONS_JA:
            reason = L(_ISSUE_REASONS_JA[reason], reason)
        (self.unknown if unknown else self.unsupported).append(
            {"file": str(path), "line": line, "text": text, "reason": reason})

    def record(self, key, value, path, line, text, *, end_line=None):
        if key in self.parsed:
            self.issue(path, line, text, L("同じ項目が重複して指定されています", "duplicate assignment"))
        self.parsed[key] = value
        self.provenance[key] = {"file": str(path), "line": line, "text": text}
        if end_line is not None:
            self.provenance[key]["end_line"] = end_line


def _read(path, result):
    if not path.is_file() or path.stat().st_size > MAX_BYTES:
        raise NativeImportError(L(f"入力がないか読み込み上限を超えています: {path}",
                                  f"input missing or exceeds import size limit: {path}"))
    with open(path, "rb") as stream:
        raw = stream.read(MAX_BYTES + 1)
    if len(raw) > MAX_BYTES:
        raise NativeImportError(L("入力が読み込み上限を超えています", "input exceeds the import size limit"))
    result.files[str(path)] = {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
    return raw.decode("utf-8-sig")


def _number(text):
    value = float(text.replace("D", "e").replace("d", "e"))
    if not math.isfinite(value):
        raise ValueError("non-finite value")
    return int(value) if value.is_integer() else value


def _scalar(text):
    text = text.strip()
    if text.lower() in (".true.", "true", "t"):
        return True
    if text.lower() in (".false.", "false", "f"):
        return False
    if text[:1] in ("'", '"') and text[-1:] == text[:1]:
        return text[1:-1]
    try:
        return _number(text)
    except ValueError:
        return text


def _count(value):
    if not isinstance(value, int) or isinstance(value, bool) or not 0 < value <= MAX_ATOMS:
        raise NativeImportError(L("原子数が不正か読み込み上限を超えています", "invalid atom count or import atom limit exceeded"))


def _structure(atoms, path, result):
    import numpy as np
    _count(len(atoms))
    if not np.isfinite(atoms.positions).all() or not np.isfinite(atoms.cell).all():
        raise NativeImportError(L("構造に有限でない座標またはセル成分があります", "the structure contains non-finite coordinates or cell components"))
    if atoms.constraints or atoms.has("momenta"):
        result.issue(path, 1, "structure", L("固定・速度の読み込みは未対応", "constraint / velocity import is unsupported"))
    return Structure(source="file", source_ref=str(path), atoms=AtomsData.from_ase(atoms))


def _finish(result, structure, method, kpoints=None, task=None):
    if result.unknown or result.unsupported:
        return
    spec = CalculationSpec(structure=structure, method=method, kpoints=kpoints, **({"task": task} if task else {}))
    def leaves(data, prefix=""):
        for key, value in data.items():
            name = f"{prefix}.{key}" if prefix else key
            if isinstance(value, dict) and value:
                yield from leaves(value, name)
            else:
                yield name, value
    # Every non-source default is visible, including scheduler values; none is
    # represented as originating in the imported calculation.
    result.inferred_defaults = {k: v for k, v in leaves(spec.model_dump(mode="json"))
                               if k not in result.provenance and not k.startswith(("structure.atoms.", "meta."))}
    result.spec = spec


def _vasp(path, result):
    from ase.io import read
    root = path if path.is_dir() else path.parent
    incar, poscar, kp = (root / name for name in ("INCAR", "POSCAR", "KPOINTS"))
    result.unresolved_dependencies.append({"kind": "POTCAR library", "file": "POTCAR", "source_file": str(poscar),
        "line": 6, "reason": L("POTCAR の内容・種類・ライブラリは読み込まず、確認もしていません。", "POTCAR contents, variants and library are neither read nor verified.")})
    method = {"code": "vasp"}
    result.parsed["method.code"] = "vasp"
    allowed = {"ENCUT", "EDIFF", "NELM", "ISMEAR", "SIGMA", "ISPIN", "ALGO", "PREC", "LREAL", "NELMIN", "LASPH", "LMAXMIX", "NBANDS", "ISYM"}
    for no, raw in enumerate(_read(incar, result).splitlines(), 1):
        for chunk in re.split("[;]", re.split("[#!]", raw, maxsplit=1)[0]):
            if not chunk.strip():
                continue
            match = re.fullmatch(r"\s*([A-Za-z][A-Za-z0-9_]*)\s*=\s*(.+?)\s*", chunk)
            if not match:
                result.issue(incar, no, raw, "unsupported INCAR syntax")
                continue
            key, value = match[1].upper(), _scalar(match[2])
            if key in allowed:
                method[key.lower()] = value
                result.record("method." + key.lower(), value, incar, no, raw)
            elif key == "SYSTEM":
                # SYSTEM is a title only, not a chemical structure declaration.
                method.setdefault("extra_incar", {})[key] = str(value)
                result.record("method.extra_incar.SYSTEM", str(value), incar, no, raw)
            elif key in {"NSW", "IBRION"}:
                result.record("native." + key, value, incar, no, raw)
                if (key == "NSW" and value != 0) or (key == "IBRION" and value != -1):
                    result.issue(incar, no, raw, "only single-point NSW=0 / IBRION=-1 is supported")
            else:
                result.issue(incar, no, raw, "unmapped INCAR keyword", unknown=True)
    text = _read(poscar, result)
    lines = text.splitlines()
    if len(lines) < 8 or all(x.isdigit() for x in lines[5].split()):
        raise NativeImportError(L("POSCAR に元素記号を明示してください (VASP 5 形式)", "POSCAR requires explicit element symbols (VASP 5 format)"))
    from ase.data import atomic_numbers
    elements = lines[5].split()
    counts = [int(x) for x in lines[6].split()]
    if len(counts) != len(elements) or any(x not in atomic_numbers for x in elements) or any(n <= 0 for n in counts):
        raise NativeImportError(L("POSCAR の元素記号と正の原子数を対応させてください", "POSCAR element symbols must match positive atom counts"))
    scale = [_number(x) for x in lines[1].split()]
    if len(scale) != 1 or scale[0] == 0:
        raise NativeImportError(L("POSCAR のスケール係数は 0 以外の値を 1 つだけ指定してください", "only a single nonzero POSCAR scale factor is supported"))
    _count(sum(counts))
    if lines[7].strip().lower().startswith("s"):
        result.issue(poscar, 8, lines[7], "Selective dynamics import is unsupported")
        return
    if len(lines) < 8 + sum(counts):
        raise NativeImportError(L("POSCAR の原子座標行が足りません", "POSCAR has fewer atomic-position rows than declared"))
    if any(x.strip() for x in lines[8 + sum(counts):]):
        result.issue(poscar, 9 + sum(counts), "POSCAR trailing data", "velocities / trailing sections are unsupported")
        return
    atoms = read(StringIO(text), format="vasp")
    structure = _structure(atoms, poscar, result)
    result.record("native.POSCAR.scale", scale[0], poscar, 2, lines[1])
    result.record("structure.atoms.cell", structure.atoms.model_dump(mode="json")["cell"], poscar, 3,
                  "\n".join(lines[2:5]), end_line=5)
    result.record("structure.atoms.symbols", structure.atoms.symbols, poscar, 6,
                  "\n".join(lines[5:7]), end_line=7)
    result.record("structure.atoms.positions", structure.atoms.model_dump(mode="json")["positions"], poscar, 9,
                  "\n".join(lines[8:8 + sum(counts)]), end_line=8 + sum(counts))
    result.provenance["structure.atoms.positions"]["coordinate_mode"] = lines[7]
    result.provenance["structure.atoms.positions"]["depends_on"] = ["native.POSCAR.scale", "structure.atoms.cell"]
    klines = _read(kp, result).splitlines()
    if len(klines) < 4 or klines[1].strip() != "0" or klines[2].strip().lower()[:1] not in {"g", "m"}:
        result.issue(kp, 2, "KPOINTS", "only automatic Gamma / Monkhorst-Pack meshes are supported")
        return
    mesh = tuple(int(x) for x in klines[3].split())
    shift = tuple(float(x) for x in klines[4].split()) if len(klines) > 4 and klines[4].strip() else (0., 0., 0.)
    if len(mesh) != 3 or min(mesh) < 1 or len(shift) != 3 or any(x != 0 for x in shift) or any(x.strip() for x in klines[5:]):
        result.issue(kp, 4, klines[3], "only positive meshes without shifts or trailing sections are supported")
        return
    method["kpoints_centering"] = "gamma" if klines[2].strip().lower().startswith("g") else "monkhorst-pack"
    result.record("kpoints.mesh", list(mesh), kp, 4, klines[3])
    result.record("method.kpoints_centering", method["kpoints_centering"], kp, 3, klines[2])
    _finish(result, structure, method, {"mode": "mesh", "mesh": mesh, "shift": shift})


def _qe(path, result):
    from ase.io import read
    from ase.data import atomic_numbers
    path = path / "pw.in" if path.is_dir() else path
    text = _read(path, result)
    method = {"code": "espresso"}
    values, section, cards = {}, None, []
    allowed = {"system": {"ecutwfc", "ecutrho", "occupations", "smearing", "degauss", "nspin", "input_dft"},
               "electrons": {"conv_thr", "electron_maxstep", "mixing_beta"}}
    for no, raw in enumerate(text.splitlines(), 1):
        # Quotes protect comments and commas. Single-line namelists are accepted.
        chunks = re.findall(r"'[^']*'|\"[^\"]*\"|[^'\"]+", raw)
        clean = ""
        for chunk in chunks:
            if chunk.startswith(("'", '"')):
                clean += chunk
            else:
                clean += chunk.split("!", 1)[0]
                if "!" in chunk:
                    break
        clean = clean.strip()
        if not clean:
            continue
        if clean.startswith("&"):
            match = re.match(r"&(\w+)", clean)
            if section or not match:
                raise NativeImportError(L("namelist の開始位置または名前が不正です", "invalid namelist start or name"))
            section = match[1].lower()
            clean = clean[match.end():].strip()
        if section:
            tokens = re.split(r",(?=(?:[^']*'[^']*')*[^']*$)", clean)
            for token in tokens:
                token = token.strip()
                end = token.endswith("/")
                if end:
                    token = token[:-1].strip()
                if token:
                    match = re.fullmatch(r"([\w()]+)\s*=\s*(.+)", token)
                    if not match:
                        result.issue(path, no, raw, "unsupported namelist syntax")
                    else:
                        key, value = match[1].lower(), _scalar(match[2])
                        full = section + "." + key
                        if full in values:
                            result.issue(path, no, raw, "duplicate assignment")
                        values[full] = value
                        if key in allowed.get(section, set()):
                            method[key] = value
                            result.record("method." + key, value, path, no, raw)
                        elif full in {"system.ibrav", "system.nat", "system.ntyp", "control.calculation"}:
                            result.record("native." + full, value, path, no, raw)
                        elif full in {"control.tstress", "control.tprnfor"}:
                            if not isinstance(value, bool):
                                result.issue(path, no, raw, "output flags require logical values")
                            method.setdefault("extra", {}).setdefault("control", {})[key] = value
                            result.record("method.extra.control." + key, value, path, no, raw)
                        elif full in {"control.prefix", "control.outdir", "control.pseudo_dir"}:
                            result.record("native." + full, value, path, no, raw)
                            result.not_applied.append({"field": full, "value": value, "file": str(path), "line": no,
                                "reason": L("元入力のファイル配置。再生成の配置やパラメータのライブラリには適用しません。",
                                            "Source file layout; not applied to regenerated paths or parameter libraries.")})
                        else:
                            result.issue(path, no, raw, "unmapped namelist variable", unknown=True)
                if end:
                    section = None
            continue
        cards.append((no, clean))
    if section:
        raise NativeImportError(L("namelist の終端「/」がありません", "namelist is missing its terminating '/'"))
    _count(values.get("system.nat"))
    _count(values.get("system.ntyp"))
    if values.get("system.ibrav") != 0 or values.get("control.calculation", "scf") != "scf":
        result.issue(path, 1, "calculation / ibrav", "only scf with explicit ibrav=0 is supported")
        return
    kpoints = None
    seen = set()
    i = 0
    while i < len(cards):
        no, header = cards[i]
        name = header.split()[0].upper()
        if name in seen:
            result.issue(path, no, header, "duplicate card")
        seen.add(name)
        i += 1
        count = {"CELL_PARAMETERS": 3, "ATOMIC_POSITIONS": values["system.nat"], "ATOMIC_SPECIES": values["system.ntyp"]}.get(name)
        if count is not None:
            rows = cards[i:i + count]
            if len(rows) != count:
                raise NativeImportError(L("QE カードのデータ行が足りません", "QE card has fewer data rows than declared"))
            if name in {"CELL_PARAMETERS", "ATOMIC_POSITIONS"} and not re.search(r"\b(angstrom|crystal)\b", header, re.I):
                result.issue(path, no, header, "explicit angstrom cell and angstrom / crystal positions required")
            if name == "CELL_PARAMETERS" and "angstrom" not in header.lower():
                result.issue(path, no, header, "CELL_PARAMETERS must use angstrom")
            if name == "ATOMIC_SPECIES":
                pseudo = {}
                for rno, row in rows:
                    parts = row.split()
                    if len(parts) != 3 or parts[0] not in atomic_numbers or parts[0] in pseudo:
                        result.issue(path, rno, row, "unique plain element labels are required")
                    else:
                        pseudo[parts[0]] = parts[2]
                        result.unresolved_dependencies.append({"kind": "UPF", "file": parts[2], "element": parts[0],
                            "source_file": str(path), "line": rno,
                            "reason": L("ファイル名だけを読み取りました。UPF の存在・内容は未確認で、複製していません。",
                                        "Only the filename was imported. UPF existence and contents are unverified; the file was not copied.")})
                method["pseudo"] = pseudo
                result.record("method.pseudo", pseudo, path, no, header, end_line=rows[-1][0])
            result.record("native." + name, [row for _, row in rows], path, no, header, end_line=rows[-1][0])
            i += count
        elif name == "K_POINTS":
            if "gamma" in header.lower():
                kpoints = {"mode": "gamma"}
            elif "automatic" in header.lower() and i < len(cards):
                parts = [int(x) for x in cards[i][1].split()]
                if len(parts) != 6 or min(parts[:3]) < 1 or any(x not in (0, 1) for x in parts[3:]):
                    raise NativeImportError(L("automatic k 点のメッシュまたはシフトが不正です", "invalid automatic k-point mesh or shift"))
                kpoints = {"mode": "mesh", "mesh": parts[:3], "shift": [x / 2 for x in parts[3:]]}
                i += 1
            else:
                result.issue(path, no, header, "only gamma / automatic k points supported")
            result.record("kpoints", kpoints, path, no, header)
        else:
            result.issue(path, no, header, "unmapped QE card", unknown=True)
    if kpoints is None:
        result.issue(path, 1, "K_POINTS", "explicit k points required")
    if result.unknown or result.unsupported:
        return
    structure = _structure(read(StringIO(text), format="espresso-in"), path, result)
    _finish(result, structure, method, kpoints)


def _lammps_md(commands, result, path, method):
    # Inverse of the generated NVE / Nose-Hoover subset; no velocity guessing.
    # https://docs.lammps.org/velocity.html
    # https://docs.lammps.org/units.html
    # https://docs.lammps.org/fix_nve.html
    # https://docs.lammps.org/fix_nh.html
    run = commands.get("run", [])
    if len(run) != 1 or not run[0].isdigit() or int(run[0]) <= 0:
        return None
    velocity = commands.get("velocity", [])
    if len(velocity) != 10 or velocity[:2] != ["all", "create"] or velocity[4:] != ["mom", "yes", "rot", "no", "dist", "gaussian"]:
        result.issue(path, 1, "velocity", "MD requires explicit generated velocity create; restart / supplied velocities are unsupported")
        return None
    temperature, seed = _number(velocity[2]), _number(velocity[3])
    if temperature <= 0 or not isinstance(seed, int) or not 0 < seed < 900000000:
        raise NativeImportError(L("MD の温度または乱数種が不正です", "invalid MD temperature or random seed"))
    dt = commands.get("timestep", [])
    every = commands.get("thermo", [])
    if len(dt) != 1 or len(every) != 1 or not every[0].isdigit() or int(every[0]) <= 0:
        raise NativeImportError(L("時間刻みと正の出力間隔を明示してください", "an explicit timestep and positive output interval are required"))
    factor = 1000 if method.get("units") == "metal" else 1
    timestep = _number(dt[0]) * factor
    if timestep <= 0:
        raise NativeImportError(L("MD の時間刻みは正の値にしてください", "MD timestep must be positive"))
    md = {"steps": int(run[0]), "temperature_k": temperature, "timestep_fs": timestep, "dump_interval": int(every[0])}
    fix = commands.get("fix", [])
    if fix == ["integ", "all", "nve"]:
        md["ensemble"] = "NVE"
    elif len(fix) == 7 and fix[:4] == ["integ", "all", "nvt", "temp"]:
        start, stop, damping = (_number(x) for x in fix[4:])
        if start != stop or start != temperature or damping <= 0:
            result.issue(path, 1, "fix", "temperature ramps / unequal initial and target temperatures are unsupported")
            return None
        md.update(ensemble="NVT", thermostat="nose_hoover", coupling_time_fs=damping * factor)
    else:
        result.issue(path, 1, "fix", "only generated NVE / Nose-Hoover NVT without extra fix options supported")
        return None
    method["seed"] = seed
    mapping = {"steps": "run", "temperature_k": "velocity", "timestep_fs": "timestep", "dump_interval": "thermo",
               "ensemble": "fix", "thermostat": "fix", "coupling_time_fs": "fix"}
    for name, value in md.items():
        src = result.provenance["native." + mapping[name]]
        result.record("task.md." + name, value, path, src["line"], src["text"])
    src = result.provenance["native.velocity"]
    result.record("method.seed", seed, path, src["line"], src["text"])
    result.record("task.type", "molecular_dynamics", path, result.provenance["native.run"]["line"], "run " + run[0])
    return {"type": "molecular_dynamics", "md": md}


def _lammps(path, result):
    from ase.io import read
    from ase.data import atomic_numbers
    path = path / "in.lammps" if path.is_dir() else path
    method = {"code": "lammps"}
    commands = {}
    command_rows = []
    for no, raw in enumerate(_read(path, result).splitlines(), 1):
        tokens = shlex.split(raw, comments=True)
        if not tokens:
            continue
        key, args = tokens[0], tokens[1:]
        command_rows.append(tokens)
        if any(c in raw.split("#", 1)[0] for c in "$&"):
            result.issue(path, no, raw, "variables / continuation are unsupported")
        if key in commands:
            result.issue(path, no, raw, "duplicate command")
        commands[key] = args
        if key in {"units", "atom_style", "pair_style", "pair_coeff"}:
            value = " ".join(args)
            method[key] = value
            result.record("method." + key, value, path, no, raw)
        elif key in {"read_data", "boundary", "run", "thermo_style", "thermo", "velocity", "timestep", "dump", "dump_modify", "fix", "write_data"}:
            result.record("native." + key, args, path, no, raw)
        else:
            result.issue(path, no, raw, "unmapped LAMMPS command", unknown=True)
    task = _lammps_md(commands, result, path, method)
    if (commands.get("run") != ["0"] and task is None) or commands.get("atom_style") != ["atomic"] or commands.get("units") not in (["metal"], ["real"]):
        result.issue(path, 1, "run / atom_style / units", "explicit supported task, atom_style atomic and units metal / real required")
    # Analytic lj/cut has no external potential references to guess or copy.
    if not commands.get("pair_style") or commands["pair_style"][0] != "lj/cut":
        result.issue(path, 1, "pair_style", "only analytic lj/cut supported; potential dependencies are not inferred")
    if not commands.get("pair_coeff"):
        result.issue(path, 1, "pair_coeff", "explicit pair_coeff required")
    else:
        coeff = commands["pair_coeff"]
        if len(coeff) not in {4, 5} or any(not re.fullmatch(r"\d*\*?\d*", x) or not x for x in coeff[:2]):
            result.issue(path, 1, "pair_coeff", "only numeric lj/cut coefficients supported")
        else:
            for value in coeff[2:]:
                _number(value)
    order = list(commands)
    if "read_data" in order and any(key in order and order.index(key) > order.index("read_data") for key in ("units", "atom_style", "boundary")):
        result.issue(path, 1, "read_data", "units, atom_style and boundary must precede read_data")
    boundary = commands.get("boundary")
    if boundary is None or len(boundary) != 3 or any(x not in {"p", "f"} for x in boundary):
        result.issue(path, 1, "boundary", "explicit p / f boundary flags required")
    ref = commands.get("read_data", [])
    if len(ref) != 1 or Path(ref[0]).is_absolute() or ".." in Path(ref[0]).parts:
        result.issue(path, 1, "read_data", "one bundle-relative data filename required")
        return
    data = path.parent / ref[0]
    if not data.resolve().is_relative_to(path.parent.resolve()):
        result.issue(path, 1, ref[0], "data symlink leaves input bundle")
        return
    text = _read(data, result)
    # The generated data title is an explicit type mapping, not a mass guess.
    header_map = {}
    title = text.splitlines()[0] if text.splitlines() else ""
    title_match = re.fullmatch(r"adit: atom types ((?:\d+=[A-Z][a-z]?\s*)+) \(written by ASE\)", title)
    if title_match:
        pairs = re.findall(r"(\d+)=([A-Z][a-z]?)", title_match[1])
        if len({int(k) for k, _ in pairs}) != len(pairs):
            result.issue(data, 1, title, "duplicate atom type in generated title")
        header_map = {int(k): v for k, v in pairs}
        result.record("native.atom_type_title", header_map, data, 1, title)
    match = re.search(r"^\s*(\d+)\s+atoms\s*$", text, re.M)
    _count(int(match[1]) if match else None)
    labels, in_mass, in_atoms = {}, False, False
    atom_ids = []
    for no, raw in enumerate(text.splitlines(), 1):
        plain = raw.split("#", 1)[0].strip()
        if no > 1 and plain and plain[0].isalpha() and plain not in {"Masses", "Atoms"}:
            result.issue(data, no, raw, "unsupported LAMMPS data section")
        if raw.strip().split("#", 1)[0].strip() == "Masses":
            in_mass = True
            in_atoms = False
            continue
        if plain == "Atoms":
            in_atoms, in_mass = True, False
            continue
        if in_atoms and plain:
            parts = plain.split()
            if len(parts) != 5:
                result.issue(data, no, raw, "only five-column atomic data rows supported")
            else:
                atom_ids.append(int(parts[0]))
                if not 0 < int(parts[1]) <= MAX_ATOMS:
                    result.issue(data, no, raw, "invalid atom type")
                for value in parts[2:]:
                    _number(value)
        if in_mass and raw.strip():
            if raw.lstrip()[0].isalpha():
                in_mass = False
                continue
            match = re.fullmatch(r"\s*(\d+)\s+\S+\s+#\s*([A-Z][a-z]?)\s*", raw)
            mass_tokens = plain.split()
            if match is None and len(mass_tokens) == 2 and mass_tokens[0].isdigit() and int(mass_tokens[0]) in header_map:
                match = re.fullmatch(r"(\d+) ([A-Z][a-z]?)", mass_tokens[0] + " " + header_map[int(mass_tokens[0])])
            if len(mass_tokens) == 2 and _number(mass_tokens[1]) <= 0:
                result.issue(data, no, raw, "mass must be positive")
            if match and match[2] in atomic_numbers:
                if int(match[1]) in labels:
                    result.issue(data, no, raw, "duplicate mass / atom type")
                labels[int(match[1])] = match[2]
                if int(match[1]) in header_map and header_map[int(match[1])] != match[2]:
                    result.issue(data, no, raw, "Masses label disagrees with generated type title")
                result.record(f"method.type_elements.{int(match[1]) - 1}", match[2], data, no, raw)
                if int(match[1]) in header_map:
                    result.provenance[f"method.type_elements.{int(match[1]) - 1}"]["depends_on"] = ["native.atom_type_title"]
            else:
                result.issue(data, no, raw, "Masses rows require explicit '# Element' labels")
    types = re.search(r"^\s*(\d+)\s+atom types\s*$", text, re.M)
    declared_atoms = int(re.search(r"^\s*(\d+)\s+atoms\s*$", text, re.M)[1])
    if sorted(atom_ids) != list(range(1, declared_atoms + 1)):
        result.issue(data, 1, "Atoms", "atom IDs must cover 1..N exactly")
    if not types or len(labels) != int(types[1]) or sorted(labels) != list(range(1, len(labels) + 1)):
        result.issue(data, 1, "Masses", "all atom types require explicit element labels")
    if result.unknown or result.unsupported:
        return
    atoms = read(StringIO(text), format="lammps-data", atom_style="atomic", units=method["units"],
                 Z_of_type={key: atomic_numbers[value] for key, value in labels.items()})
    atoms.pbc = [x == "p" for x in boundary]
    method["data_file"] = str(data.resolve())
    method["type_elements"] = [labels[key] for key in sorted(labels)]
    structure = _structure(atoms, data, result)
    # MD and generated output commands are accepted only if the generator can
    # reproduce the entire command sequence, including initialization and IDs.
    if task or any(k in commands for k in ("thermo_style", "thermo", "dump", "dump_modify", "fix", "velocity", "timestep", "write_data")):
        from adit.codes.lammps import LammpsGenerator
        candidate = CalculationSpec(structure=structure, method=method, **({"task": task} if task else {}))
        generated = [shlex.split(row, comments=True) for row in LammpsGenerator().in_lammps(candidate).splitlines()]
        generated = [row for row in generated if row]
        normalize = lambda rows: [[_scalar(token) for token in row] for row in rows]
        if normalize(generated) != normalize(command_rows):
            result.issue(path, 1, "command sequence", "commands differ from the supported generated sequence; options / ordering / initialization cannot be preserved")
    _finish(result, structure, method, task=task)


_GMX_THERMOSTAT = {"berendsen": "berendsen", "nose-hoover": "nose_hoover", "andersen": "andersen", "v-rescale": "csvr"}
_GMX_METHOD_KEYS = {"coulombtype": "coulombtype", "rcoulomb": "rcoulomb_nm", "rvdw": "rvdw_nm",
                    "constraints": "constraints", "pcoupl": "pcoupl", "compressibility": "compressibility_per_bar",
                    "define": "define", "gen-seed": "gen_seed"}
_GMX_IGNORED = {"cutoff-scheme", "tc-grps", "pcoupltype", "nstlog", "nstenergy", "gen-temp", "gen-vel", "continuation"}


def _gromacs(path, result):
    from ase.io import read

    root = path if path.is_dir() else path.parent
    mdp_path = path if path.is_file() and path.suffix == ".mdp" else root / "grompp.mdp"
    text = _read(mdp_path, result)
    values: dict[str, tuple] = {}
    for no, raw in enumerate(text.splitlines(), 1):
        line = re.split(r"[;#]", raw, maxsplit=1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            result.issue(mdp_path, no, raw, L("「名前 = 値」の形ではありません", "not of the form name = value"))
            continue
        key, _, value = line.partition("=")
        key, value = key.strip().lower(), value.strip()
        if key in values:
            result.issue(mdp_path, no, raw, L("同じ項目が重複して指定されています", "duplicate assignment"))
        values[key] = (value, no, raw)
    method: dict = {"code": "gromacs"}
    task: dict = {}
    md: dict = {}

    def take(key, target, convert=lambda x: x, *, into=method):
        if key not in values:
            return None
        value, no, raw = values[key]
        try:
            converted = convert(value)
        except (TypeError, ValueError):
            result.issue(mdp_path, no, raw, L("値を数として読めません", "the value is not a number"))
            return None
        into[target.split(".")[-1]] = converted
        prefix = {id(method): "method", id(task): "task", id(md): "task.md"}[id(into)]
        result.record(f"{prefix}.{target.split('.')[-1]}", converted, mdp_path, no, raw)
        return converted

    integrator = values.get("integrator", ("", 1, ""))[0].lower()
    if integrator in ("steep", "l-bfgs"):
        task["type"] = "geometry_optimization"
        task["optimizer"] = "LBFGS" if integrator == "l-bfgs" else "SteepestDescent"
        result.record("task.type", "geometry_optimization", mdp_path, *values["integrator"][1:][::-1][::-1])
        take("nsteps", "max_steps", int, into=task)
        if "emtol" in values:  # kJ/mol/nm → eV/Å
            value, no, raw = values["emtol"]
            task["force_tolerance_ev_per_ang"] = float(value) / 964.8533212
            result.record("task.force_tolerance_ev_per_ang", task["force_tolerance_ev_per_ang"], mdp_path, no, raw)
    elif integrator in ("md", "sd"):
        steps = int(values.get("nsteps", ("0", 1, ""))[0])
        task["type"] = "single_point" if steps == 0 and "dt" not in values else "molecular_dynamics"
        result.record("task.type", task["type"], mdp_path, values.get("integrator", ("", 1, ""))[1],
                      values.get("integrator", ("", 1, ""))[2])
        if task["type"] == "molecular_dynamics":
            take("nsteps", "steps", int, into=md)
            take("dt", "timestep_fs", lambda v: float(v) * 1000.0, into=md)
            take("nstxout-compressed", "dump_interval", int, into=md)
            take("ref-t", "temperature_k", float, into=md)
            take("tau-t", "coupling_time_fs", lambda v: float(v) * 1000.0, into=md)
            take("ref-p", "pressure_bar", float, into=md)
            take("tau-p", "barostat_time_fs", lambda v: float(v) * 1000.0, into=md)
            pcoupl = values.get("pcoupl", ("no", 1, ""))[0].lower()
            tcoupl = values.get("tcoupl", ("no", 1, ""))[0].lower()
            md["ensemble"] = "NPT" if pcoupl not in ("no", "") else ("NVE" if tcoupl in ("no", "") and integrator != "sd" else "NVT")
            result.record("task.md.ensemble", md["ensemble"], mdp_path, values.get("pcoupl", values.get("tcoupl", ("", 1, "")))[1],
                          values.get("pcoupl", values.get("tcoupl", ("", 1, "")))[2])
            if integrator == "sd":
                md["thermostat"] = "langevin"
                result.record("task.md.thermostat", "langevin", mdp_path, *values["integrator"][1:])
            elif tcoupl not in ("no", ""):
                mapped = _GMX_THERMOSTAT.get(tcoupl)
                if mapped is None:
                    result.issue(mdp_path, values["tcoupl"][1], values["tcoupl"][2],
                                 L("共通の欄に対応のない熱浴です", "thermostat has no equivalent in the shared fields"))
                else:
                    md["thermostat"] = mapped
                    result.record("task.md.thermostat", mapped, mdp_path, *values["tcoupl"][1:])
    else:
        result.issue(mdp_path, values.get("integrator", ("", 1, ""))[1], values.get("integrator", ("", 1, ""))[2],
                     L("対応していない integrator です (steep / l-bfgs / md / sd)",
                       "unsupported integrator (steep / l-bfgs / md / sd)"))
    for key, field_name in _GMX_METHOD_KEYS.items():
        if key not in values:
            continue
        if key in ("pcoupl", "tcoupl") and values[key][0].strip().lower() in ("no", ""):
            continue
        if field_name in ("rcoulomb_nm", "rvdw_nm", "compressibility_per_bar"):
            take(key, field_name, float)
        elif field_name == "gen_seed":
            take(key, field_name, int)
        else:
            take(key, field_name)
    for key in values:
        if key in _GMX_METHOD_KEYS or key in _GMX_IGNORED:
            continue
        if key in ("integrator", "nsteps", "dt", "emtol", "nstxout-compressed", "ref-t", "tau-t", "ref-p", "tau-p", "tcoupl"):
            continue
        result.issue(mdp_path, values[key][1], values[key][2], L("共通の欄に対応のない mdp の項目です",
                                                                 "mdp option without an equivalent in the shared fields"), unknown=True)
    structure = None
    for name in ("conf.gro", "conf.pdb"):
        candidate = root / name
        if candidate.is_file():
            try:
                atoms = read(candidate)
            except Exception as exc:
                result.issue(candidate, 1, name, L("構造を読めません", "cannot read the structure"))
                result.unsupported[-1]["detail"] = str(exc)
                break
            _read(candidate, result)
            structure = _structure(atoms, candidate, result)
            result.provenance["structure.atoms.symbols"] = {"file": str(candidate), "line": 1, "text": name}
            result.provenance["structure.atoms.positions"] = {"file": str(candidate), "line": 1, "text": name}
            result.provenance["structure.atoms.cell"] = {"file": str(candidate), "line": 1, "text": name}
            result.provenance["structure.atoms.pbc"] = {"file": str(candidate), "line": 1, "text": name}
            method["structure_file"] = str(candidate)
            break
    if structure is None:
        result.issue(root, 1, "conf.gro", L("構造のファイル (conf.gro / conf.pdb) がありません",
                                            "no structure file (conf.gro / conf.pdb)"))
        return
    topology = root / "topol.top"
    if topology.is_file():
        _read(topology, result)
        method["topology_file"] = str(topology)
        result.unresolved_dependencies.append({
            "file": str(topology), "kind": "topology",
            "reason": L("トポロジー (原子型・電荷・力場) の中身は読み戻しません。ADIT はこのファイルを写すだけです。",
                        "The topology (atom types, charges, force field) is not parsed; ADIT only copies this file.")})
    else:
        result.issue(root, 1, "topol.top", L("トポロジー (topol.top) がありません", "no topology (topol.top)"))
        return
    from adit.spec import GromacsMethod, MDSettings, Task

    if md:
        task["md"] = MDSettings(**md)
    _finish(result, structure, GromacsMethod(**method), None, Task(**task) if task else None)


def import_native(source: Path | str, code: str | None = None) -> NativeImportResult:
    """Inspect a bundle. Unsupported data produces a report and no usable Spec."""
    path = Path(source).expanduser().resolve()
    if code is None:
        matches = [name for name, filename in (("vasp", "INCAR"), ("espresso", "pw.in"), ("lammps", "in.lammps"),
                                               ("gromacs", "grompp.mdp"))
                   if (path / filename).is_file()] if path.is_dir() else []
        if not path.is_dir():
            matches = [{"INCAR": "vasp", "POSCAR": "vasp", "KPOINTS": "vasp", "pw.in": "espresso", "in.lammps": "lammps",
                        "grompp.mdp": "gromacs"}.get(path.name)]
        if len(matches) != 1 or not matches[0]:
            raise NativeImportError(L("入力コードを明示してください", "specify the native input code explicitly"))
        code = matches[0]
    code = {"qe": "espresso"}.get(code, code)
    if code not in {"vasp", "espresso", "lammps", "gromacs"}:
        raise NativeImportError(L(f"未対応の入力コード: {code}", f"unsupported native input code: {code}"))
    result = NativeImportResult(code)
    try:
        {"vasp": _vasp, "espresso": _qe, "lammps": _lammps, "gromacs": _gromacs}[code](path, result)
    except NativeImportError as exc:
        result.issue(path, 1, "", str(exc))
        result.spec = None
    except (OSError, UnicodeError, ValueError, KeyError, IndexError, TypeError, StopIteration) as exc:
        result.issue(path, 1, "", L(
            "入力を読み込めません。ファイルの形式と数値を確認してください。",
            "Could not import the input. Check its format and numeric values."))
        result.unsupported[-1]["detail"] = str(exc)
        result.spec = None
    return result
