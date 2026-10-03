
from __future__ import annotations

import os
import shlex
import sys
from pathlib import PurePosixPath, PureWindowsPath


def readme_shell_path(path: str, *, windows: bool | None = None) -> str:
    windows = os.name == "nt" if windows is None else windows
    if not windows:
        return shlex.quote(path)
    w = PureWindowsPath(path)
    if len(w.drive) == 2 and w.drive[1] == ":":
        return shlex.quote(str(PurePosixPath("/mnt", w.drive[0].lower(), *w.parts[1:])))
    return shlex.quote(w.as_posix())


def upload_basename(name: str) -> str:
    base = PureWindowsPath(name.replace("\x00", "").replace("/", "\\")).name
    return base if base not in ("", ".", "..") else "upload"


def ensure_printable_stdio() -> None:
    for stream in (sys.stdout, sys.stderr):
        enc = (getattr(stream, "encoding", None) or "").lower().replace("-", "").replace("_", "")
        if enc in ("utf8", "utf8sig") or not hasattr(stream, "reconfigure"):
            continue
        try:
            stream.reconfigure(errors="backslashreplace")
        except (ValueError, OSError):
            pass
