# Progress reports from long builds and analyses, kept per thread so callers need no callback arguments.

from __future__ import annotations

import threading
from typing import Callable

Reporter = Callable[[int, int, str], None]

_local = threading.local()


def set_reporter(cb: Reporter | None) -> None:
    _local.cb = cb


def report(done: int, total: int, what: str = "") -> None:
    cb = getattr(_local, "cb", None)
    if cb is not None:
        cb(int(done), int(total), what)


__all__ = ["set_reporter", "report"]
