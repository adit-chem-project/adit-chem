
from __future__ import annotations

from adit.lang import L


def short_number(v) -> str:
    if isinstance(v, bool) or not isinstance(v, (int, float)):
        return str(v)
    if isinstance(v, int):
        return str(v)
    return min((repr(v), format(v, ".17g")), key=len)


def parse_indices(text: str, n: int) -> list[int]:
    return parse_constraints(text, n)[0]


def parse_constraints(text: str, n: int) -> tuple[list[int], dict[str, tuple[bool, bool, bool]]]:
    atoms: set[int] = set()
    axes: dict[str, tuple[bool, bool, bool]] = {}
    for part in text.replace(" ", "").split(","):
        if not part:
            continue
        rng, _, ax = part.partition(":")
        a, _, b = rng.partition("-")
        try:
            lo = int(a); hi = int(b) if b else lo
        except ValueError as ex:
            raise ValueError(L(f"原子の番号として読めません: {part!r}", f"not an atom index: {part!r}")) from ex
        if lo < 1 or hi > n or lo > hi:
            raise ValueError(L(f"原子の番号が範囲外です (1〜{n}): {part!r}", f"atom index out of range (1..{n}): {part!r}"))
        if ax:
            if not set(ax) <= set("xyz"):
                raise ValueError(L(f"固定する軸は x y z の組み合わせにしてください: {part!r}", f"axes to fix must be a combination of x y z: {part!r}"))
            move = tuple(c not in ax for c in "xyz")
            for i in range(lo - 1, hi):
                axes[str(i)] = move
        else:
            atoms.update(range(lo - 1, hi))
    return sorted(atoms), axes


def _value(v: str):
    up = v.upper()
    if up in (".TRUE.", "TRUE", "T"):
        return True
    if up in (".FALSE.", "FALSE", "F"):
        return False
    try:
        return int(v)
    except ValueError:
        try:
            return float(v)
        except ValueError:
            return v


def parse_extra_incar(text: str) -> dict:
    out: dict = {}
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" not in line:
            raise ValueError(L(f"「KEY = value」の形ではありません: {raw!r}", f"not of the form KEY = value: {raw!r}"))
        k, _, v = line.partition("=")
        k, v = k.strip().upper(), v.strip()
        if not k:
            raise ValueError(L(f"キーが空です: {raw!r}", f"empty key: {raw!r}"))
        out[k] = _value(v)
    return out


def parse_extra_namelist(text: str) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if "=" not in line or "." not in line.split("=", 1)[0]:
            raise ValueError(L(f"「ネームリスト.変数 = 値」の形ではありません (例 system.nbnd = 20): {raw!r}",
                               f"not of the form namelist.variable = value (e.g. system.nbnd = 20): {raw!r}"))
        k, _, v = line.partition("=")
        ns, _, key = k.strip().partition(".")
        val = _value(v.strip())
        if isinstance(val, str):
            val = val.strip("'\"")
        out.setdefault(ns.lower(), {})[key.strip().lower()] = val
    return out


def parse_element_map(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        sep = "=" if "=" in line else ":"
        if sep not in line:
            raise ValueError(L(f"「元素 = 名前」の形ではありません: {raw!r}", f"not of the form element = name: {raw!r}"))
        e, _, name = line.partition(sep)
        e, name = e.strip(), name.strip()
        if not e or not name:
            raise ValueError(L(f"元素か名前が空です: {raw!r}", f"empty element or name: {raw!r}"))
        out[e] = name
    return out
