"""Entry point of the frozen executable.

With no arguments it opens the desktop window; with arguments it behaves as the command line:

    ADIT.exe                       desktop window
    ADIT.exe web                   local web server
    ADIT.exe gen spec.json out/    same as adit-gen
    ADIT.exe analyze out/run1      same as adit-analyze
    ADIT.exe report out/run1       same as adit-report
    ADIT.exe convert structure ... same as adit-convert
"""

import multiprocessing
import sys


def _console_executable() -> str | None:
    # The windowed ADIT.exe has no stdout (sys.stdout is None); the console build adit-cli next to it does
    from pathlib import Path

    exe = Path(sys.executable)
    cli = exe.with_name("adit-cli" + exe.suffix)
    return str(cli) if cli.is_file() and cli.resolve() != exe.resolve() else None


def main() -> int:
    multiprocessing.freeze_support()
    commands = {"gen": "adit.cli", "analyze": "adit.analysis.cli", "web": "adit.web.server",
                "report": "adit.report", "convert": "adit.convert"}
    if len(sys.argv) > 1 and sys.argv[1] in commands:
        if sys.stdout is None:
            cli = _console_executable()
            if cli:
                import subprocess

                return subprocess.call([cli, *sys.argv[1:]])
        import importlib

        module = importlib.import_module(commands[sys.argv[1]])
        sys.argv = [f"adit-{sys.argv[1]}"] + sys.argv[2:]
        return int(module.main() or 0)
    from adit.gui.app import main as gui_main

    return int(gui_main() or 0)


if __name__ == "__main__":
    sys.exit(main())
