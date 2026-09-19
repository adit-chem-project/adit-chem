"""Atom selection language: the part of the MDAnalysis syntax that contains no judgment."""

from __future__ import annotations

import re

import numpy as np
from ase import Atoms

from adit.errors import AditValueError
from adit.lang import L

_TOKEN = re.compile(r"\s*(\(|\)|<=|>=|<|>|,|[A-Za-z_][A-Za-z_0-9]*|-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)")


class SelectionError(AditValueError):
    pass


def parse_indices(text: str) -> list[int]:
    out: list[int] = []
    for part in text.replace(" ", "").split(","):
        if not part:
            continue
        if "-" in part[1:]:
            a, b = part.split("-", 1) if not part.startswith("-") else (part, part)
            try:
                lo, hi = int(a), int(b)
            except ValueError as ex:
                raise SelectionError(L(f"番号の範囲を読めません: {part!r}", f"cannot read the index range {part!r}")) from ex
            if lo < 1 or hi < lo:
                raise SelectionError(L(f"番号の範囲が逆か、1 より小さいです: {part!r}",
                                       f"the index range is reversed or below 1: {part!r}"))
            out += list(range(lo - 1, hi))
        else:
            try:
                value = int(part)
            except ValueError as ex:
                raise SelectionError(L(f"番号を読めません: {part!r}", f"cannot read the index {part!r}")) from ex
            if value < 1:
                raise SelectionError(L("番号は 1 から始まります", "indices start at 1"))
            out.append(value - 1)
    return sorted(set(out))


class _Parser:

    def __init__(self, text: str, atoms: Atoms):
        self.tokens = [m.group(1) for m in _TOKEN.finditer(text)]
        self.pos = 0
        self.atoms = atoms
        self.n = len(atoms)

    def peek(self) -> str | None:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def take(self) -> str:
        if self.pos >= len(self.tokens):
            raise SelectionError(L("選び方が途中で終わっています", "the selection ends unexpectedly"))
        self.pos += 1
        return self.tokens[self.pos - 1]

    def parse(self) -> np.ndarray:
        mask = self.parse_or()
        if self.peek() is not None:
            raise SelectionError(L(f"余分な文字があります: {' '.join(self.tokens[self.pos:])!r}",
                                   f"unexpected text: {' '.join(self.tokens[self.pos:])!r}"))
        return mask

    def parse_or(self) -> np.ndarray:
        mask = self.parse_and()
        while (self.peek() or "").lower() == "or":
            self.take()
            mask = mask | self.parse_and()
        return mask

    def parse_and(self) -> np.ndarray:
        mask = self.parse_unit()
        while (self.peek() or "").lower() == "and":
            self.take()
            mask = mask & self.parse_unit()
        return mask

    def parse_unit(self) -> np.ndarray:
        token = self.take()
        low = token.lower()
        if low == "not":
            return ~self.parse_unit()
        if token == "(":
            mask = self.parse_or()
            if self.peek() != ")":
                raise SelectionError(L("括弧が閉じていません", "a parenthesis is not closed"))
            self.take()
            return mask
        if low == "element":
            return self._elements()
        if low == "index":
            return self._indices()
        if low in ("x", "y", "z"):
            return self._coordinate(low)
        if low == "within":
            return self._within()
        if low == "all":
            return np.ones(self.n, dtype=bool)
        from ase.data import chemical_symbols

        if token in chemical_symbols:
            self.pos -= 1
            return self._elements()
        raise SelectionError(L(f"分からない言葉です: {token!r} (使えるのは element / index / x y z / within / and / or / not / all)",
                               f"unknown word {token!r} (use element, index, x/y/z, within, and, or, not, all)"))

    def _elements(self) -> np.ndarray:
        from ase.data import chemical_symbols

        wanted = []
        while (t := self.peek()) is not None and t in chemical_symbols:
            wanted.append(self.take())
        if not wanted:
            raise SelectionError(L("element のあとに元素記号が要ります", "element must be followed by chemical symbols"))
        syms = np.array(self.atoms.get_chemical_symbols())
        return np.isin(syms, wanted)

    def _indices(self) -> np.ndarray:
        parts = []
        while (t := self.peek()) is not None and (t == "," or t == "-" or re.fullmatch(r"-?\d+", t)):
            parts.append(self.take())
        idx = parse_indices("".join(parts))
        if idx and max(idx) >= self.n:
            raise SelectionError(L(f"原子は {self.n} 個しかありません (番号 {max(idx) + 1} は範囲外)",
                                   f"there are only {self.n} atoms (index {max(idx) + 1} is out of range)"))
        mask = np.zeros(self.n, dtype=bool)
        mask[idx] = True
        return mask

    def _coordinate(self, axis: str) -> np.ndarray:
        op = self.take()
        if op not in ("<", "<=", ">", ">="):
            raise SelectionError(L(f"{axis} のあとは < <= > >= のどれかです", f"{axis} must be followed by <, <=, > or >="))
        try:
            value = float(self.take())
        except ValueError as ex:
            raise SelectionError(L(f"{axis} {op} のあとに数が要ります", f"a number must follow {axis} {op}")) from ex
        col = self.atoms.get_positions()[:, "xyz".index(axis)]
        return {"<": col < value, "<=": col <= value, ">": col > value, ">=": col >= value}[op]

    def _within(self) -> np.ndarray:
        try:
            radius = float(self.take())
        except ValueError as ex:
            raise SelectionError(L("within のあとに距離 [Å] が要ります", "within must be followed by a distance in Å")) from ex
        if (self.take() or "").lower() != "of":
            raise SelectionError(L("書き方は「within 5 of element O」です", "write it as: within 5 of element O"))
        other = self.parse_unit()
        if not other.any():
            return np.zeros(self.n, dtype=bool)
        from adit.analysis.compute import _mic_step

        pos = self.atoms.get_positions()
        cell = np.asarray(self.atoms.cell, dtype=float)
        periodic = bool(np.any(self.atoms.pbc)) and abs(np.linalg.det(cell)) > 0
        mask = np.zeros(self.n, dtype=bool)
        targets = pos[other]
        for start in range(0, self.n, 512):
            block = pos[start:start + 512]
            diff = block[:, None, :] - targets[None, :, :]
            flat = diff.reshape(-1, 3)
            if periodic:
                flat = _mic_step(flat, cell)
            d = np.linalg.norm(flat, axis=1).reshape(len(block), len(targets))
            mask[start:start + 512] = (d <= radius).any(axis=1)
        return mask


def select(atoms: Atoms, expression: str | None) -> np.ndarray:
    """Return the 0-based indices of the selected atoms. An empty selection means all atoms."""
    if expression is None or not expression.strip():
        return np.arange(len(atoms))
    mask = _Parser(expression, atoms).parse()
    return np.flatnonzero(mask)


def describe(atoms: Atoms, expression: str | None) -> str:
    idx = select(atoms, expression)
    syms = np.array(atoms.get_chemical_symbols())[idx]
    kinds = ", ".join(f"{s} {int((syms == s).sum())}" for s in sorted(set(syms.tolist())))
    return L(f"選んだ原子: {len(idx)} 個 ({kinds or '無し'}) / 全 {len(atoms)} 個",
             f"selected atoms: {len(idx)} ({kinds or 'none'}) out of {len(atoms)}")


class _Translator:
    """Rewrite an ADIT selection as a VMD selection or an OVITO expression (same grammar, no atoms needed)."""

    def __init__(self, text: str, target: str):
        self.tokens = [m.group(1) for m in _TOKEN.finditer(text)]
        self.pos = 0
        self.target = target

    def peek(self) -> str | None:
        return self.tokens[self.pos] if self.pos < len(self.tokens) else None

    def take(self) -> str:
        if self.pos >= len(self.tokens):
            raise SelectionError(L("選び方が途中で終わっています", "the selection ends unexpectedly"))
        self.pos += 1
        return self.tokens[self.pos - 1]

    def parse(self) -> str:
        text = self.parse_or()
        if self.peek() is not None:
            raise SelectionError(L(f"余分な文字があります: {' '.join(self.tokens[self.pos:])!r}",
                                   f"unexpected text: {' '.join(self.tokens[self.pos:])!r}"))
        return text

    def parse_or(self) -> str:
        parts = [self.parse_and()]
        while (self.peek() or "").lower() == "or":
            self.take()
            parts.append(self.parse_and())
        joiner = " or " if self.target == "vmd" else " || "
        return joiner.join(parts) if len(parts) == 1 else "(" + joiner.join(parts) + ")"

    def parse_and(self) -> str:
        parts = [self.parse_unit()]
        while (self.peek() or "").lower() == "and":
            self.take()
            parts.append(self.parse_unit())
        joiner = " and " if self.target == "vmd" else " && "
        return parts[0] if len(parts) == 1 else "(" + joiner.join(parts) + ")"

    def parse_unit(self) -> str:
        from ase.data import chemical_symbols

        token = self.take()
        low = token.lower()
        if low == "not":
            inner = self.parse_unit()
            # muparser has no documented unary not; a comparison with 0 uses only documented operators
            return f"not {inner}" if self.target == "vmd" else f"({inner}) == 0"
        if token == "(":
            inner = self.parse_or()
            if self.peek() != ")":
                raise SelectionError(L("括弧が閉じていません", "a parenthesis is not closed"))
            self.take()
            return inner if inner.startswith("(") else f"({inner})"
        if low == "element" or token in chemical_symbols:
            if token in chemical_symbols:
                self.pos -= 1
            return self._elements()
        if low == "index":
            return self._indices()
        if low in ("x", "y", "z"):
            return self._coordinate(low)
        if low == "within":
            return self._within()
        if low == "all":
            return "all" if self.target == "vmd" else "1"
        raise SelectionError(L(f"分からない言葉です: {token!r} (使えるのは element / index / x y z / within / and / or / not / all)",
                               f"unknown word {token!r} (use element, index, x/y/z, within, and, or, not, all)"))

    def _elements(self) -> str:
        from ase.data import chemical_symbols

        wanted = []
        while (t := self.peek()) is not None and t in chemical_symbols:
            wanted.append(self.take())
        if not wanted:
            raise SelectionError(L("element のあとに元素記号が要ります", "element must be followed by chemical symbols"))
        if self.target == "vmd":
            return "name " + " ".join(wanted)
        parts = [f'ParticleType == "{s}"' for s in wanted]
        return parts[0] if len(parts) == 1 else "(" + " || ".join(parts) + ")"

    def _indices(self) -> str:
        parts = []
        while (t := self.peek()) is not None and (t == "," or t == "-" or re.fullmatch(r"-?\d+", t)):
            parts.append(self.take())
        idx = parse_indices("".join(parts))
        if not idx:
            raise SelectionError(L("index のあとに番号が要ります", "index must be followed by numbers"))
        runs: list[tuple[int, int]] = []
        for i in idx:
            if runs and runs[-1][1] == i - 1:
                runs[-1] = (runs[-1][0], i)
            else:
                runs.append((i, i))
        if self.target == "vmd":
            return "index " + " ".join(f"{a}" if a == b else f"{a} to {b}" for a, b in runs)
        terms = [f"ParticleIndex == {a}" if a == b else f"(ParticleIndex >= {a} && ParticleIndex <= {b})" for a, b in runs]
        return terms[0] if len(terms) == 1 else "(" + " || ".join(terms) + ")"

    def _coordinate(self, axis: str) -> str:
        op = self.take()
        if op not in ("<", "<=", ">", ">="):
            raise SelectionError(L(f"{axis} のあとは < <= > >= のどれかです", f"{axis} must be followed by <, <=, > or >="))
        try:
            value = float(self.take())
        except ValueError as ex:
            raise SelectionError(L(f"{axis} {op} のあとに数が要ります", f"a number must follow {axis} {op}")) from ex
        if self.target == "vmd":
            return f"{axis} {op} {value:g}"
        return f"Position.{axis.upper()} {op} {value:g}"

    def _within(self) -> str:
        try:
            radius = float(self.take())
        except ValueError as ex:
            raise SelectionError(L("within のあとに距離 [Å] が要ります", "within must be followed by a distance in Å")) from ex
        if (self.take() or "").lower() != "of":
            raise SelectionError(L("書き方は `within 5 of element O` です", "write it as `within 5 of element O`"))
        other = self.parse_unit()
        if self.target != "vmd":
            raise SelectionError(L("within は OVITO の式に直せません (距離による選択は OVITO の Expand selection などで作ってください)",
                                   "within cannot be written as an OVITO expression (build distance-based selections with e.g. Expand selection in OVITO)"))
        return f"(within {radius:g} of {other})"


def to_vmd(expression: str) -> str:
    """Rewrite an ADIT selection in VMD's selection language (indices become 0-based, elements become names)."""
    return _Translator(expression, "vmd").parse()


def to_ovito(expression: str) -> str:
    """Rewrite an ADIT selection as an OVITO expression-selection string (raises SelectionError for `within`)."""
    return _Translator(expression, "ovito").parse()
