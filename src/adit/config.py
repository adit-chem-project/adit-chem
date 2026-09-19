"""Read and write cluster.toml. One profile is one run target. ADIT_CONFIG overrides the path."""

from __future__ import annotations

from adit.errors import AditError
import os
import re
import tomllib
from pathlib import Path
from typing import Literal

import tomli_w
from pydantic import BaseModel, Field, PrivateAttr, ValidationError

APP_NAME = "adit"


SUBMIT_DEFAULT = {"direct": "bash", "pbs": "qsub", "slurm": "sbatch"}
STATUS_DEFAULT = {"direct": "", "pbs": "qstat -u $USER", "slurm": "squeue -u $USER"}


class Profile(BaseModel):

    kind: Literal["direct", "pbs", "slurm"]
    modules: list[str] = Field(default_factory=list)
    code_modules: dict[str, list[str]] = Field(default_factory=dict)
    env: dict[str, str] = Field(default_factory=dict)
    commands: dict[str, str] = Field(default_factory=dict)
    select_extra: str = ""
    header_extra: list[str] = Field(default_factory=list)
    submit_command: str = ""
    status_command: str = ""
    description: str = ""
    host: str = ""          # cluster host name, for the transfer and submit commands
    user: str = ""          # login name on that host (empty: the local one)
    remote_dir: str = ""    # work directory on that host

    def modules_for(self, code: str) -> list[str]:
        return list(self.modules) + list(self.code_modules.get(code, []))

    def command_for(self, code: str, default: str) -> str:
        return self.commands.get(code) or default

    @property
    def submit(self) -> str:
        return self.submit_command or SUBMIT_DEFAULT[self.kind]

    @property
    def target(self) -> str:
        """`user@host`, `host`, or "" when no host is configured."""
        if not self.host.strip():
            return ""
        user = self.user.strip()
        return f"{user}@{self.host.strip()}" if user else self.host.strip()

    @property
    def status(self) -> str:
        return self.status_command or STATUS_DEFAULT[self.kind]

class Config(BaseModel):
    _source_path: Path | None = PrivateAttr(default=None)
    _unknown_keys: list[str] = PrivateAttr(default_factory=list)

    sk_root: str = ""
    pseudo_root: str = ""
    cp2k_data: str = ""
    templates_dir: str = ""
    default_profile: str = "local"
    enable_run: bool = True
    language: str = "ja"
    theme: str = "auto"
    window_frame: str = "auto"
    profiles: dict[str, Profile] = Field(default_factory=dict)

    @property
    def unknown_keys(self) -> list[str]:
        return list(self._unknown_keys)

    @property
    def source_path(self) -> Path | None:
        return self._source_path

    def profile(self, name: str) -> Profile:
        if name not in self.profiles:
            from adit.lang import L

            raise ConfigError(L(f"プロファイル {name!r} が環境設定ファイルにありません (あるもの: {sorted(self.profiles)})",
                                f"profile {name!r} is not in the settings file (available: {sorted(self.profiles)})"))
        return self.profiles[name]


class ConfigError(AditError):
    pass


OLD_APP_NAMES = ("vista", "qcgui")
OLD_APP_NAME = OLD_APP_NAMES[0]
OLD_ENV_PREFIXES = ("VISTA_", "QCGUI_")


def env_var(name: str, default: str = "") -> str:
    value = os.environ.get(f"ADIT_{name}")
    if value:
        return value
    for prefix in OLD_ENV_PREFIXES:
        value = os.environ.get(f"{prefix}{name}")
        if value:
            return value
    return default


def config_path() -> Path:
    env = env_var("CONFIG")
    if env:
        return Path(env).expanduser()
    if os.name == "nt":
        base = Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming"))
    else:
        base = Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config"))
    new = base / APP_NAME / "cluster.toml"
    if new.exists():
        return new
    for old_name in OLD_APP_NAMES:
        old = base / old_name / "cluster.toml"
        if old.exists():
            return old
    return new


def default_config(sk_root: str = "") -> Config:
    from adit.lang import L

    return Config(
        sk_root=sk_root,
        default_profile="local",
        profiles={
            "local": Profile(kind="direct", description=L("この PC で bash submit.sh を実行 (計算ソフトがインストール済みで PATH にあること)",
                                                          "run bash submit.sh on this PC (the program must be installed and on PATH)")),
        },
    )


class ConfigMissing(ConfigError):
    pass


_TOML_HINTS = {
    "Invalid hex value": ("\\ のあとの文字が特別な意味 (\\u や \\x で始まる文字の番号) に読まれました。Windows のパスの \\ は / か \\\\ にします",
                          "the character after \\ was read as an escape (\\u or \\x); write the \\ of a Windows path as / or \\\\"),
    "Unescaped '\\' in a string": ("文字列の中の \\ は \\\\ と 2 つ重ねるか、/ にします", "inside a string write \\ as \\\\ or use /"),
    "Illegal character in string": ("文字列の中に使えない文字があります。\\ は \\\\ と 2 つ重ねるか、/ にします", "illegal character in a string; write \\ as \\\\ or use /"),
    "Invalid value": ("値の書き方が正しくありません。パスなどの文字列は \"…\" で囲みます", "the value is not valid; strings such as paths must be enclosed in \"…\""),
    "Expected '=' after a key": ("項目の名前のあとに = がありません", "missing = after the key"),
    "Invalid statement": ("この行は「名前 = 値」の形になっていません", "this line is not of the form name = value"),
    "Unclosed string": ("文字列の \" が閉じていません", "the string's closing \" is missing"),
    "Cannot overwrite a value": ("同じ項目が 2 回書かれています", "the same key is written twice"),
    "Cannot declare": ("同じ見出し ([…]) が 2 回書かれています", "the same [table] header is written twice"),
}


def _english(text: str = "") -> bool:
    import re

    env = env_var("LANG")
    if env:
        return env.lower().startswith("en")
    m = re.search(r'^\s*language\s*=\s*["\']?(en|ja)', text, re.M)
    return bool(m and m.group(1) == "en")


def describe_config_error(ex: Exception, text: str, path: Path) -> str:
    en = _english(text)
    if isinstance(ex, tomllib.TOMLDecodeError):
        msg = getattr(ex, "msg", str(ex))
        line, col = getattr(ex, "lineno", None), getattr(ex, "colno", None)
        if line is None:
            m = re.search(r"at line (\d+), column (\d+)", str(ex))
            if m:
                line, col = int(m.group(1)), int(m.group(2))
        src = text.splitlines()[line - 1] if line and line <= len(text.splitlines()) else ""
        ja_hint, en_hint = next((v for k, v in _TOML_HINTS.items() if msg.startswith(k)), ("", ""))
        if en:
            out = [f"Cannot read the settings file {path}" + (f" (line {line}, column {col})" if line else "") + ".",
                   *([f"  line {line}:  {src.strip()}"] if src else []),
                   f"  {en_hint + ' ' if en_hint else ''}({msg})",
                   '  Enclose paths in "...", and write \\ as / or \\\\. For example:',
                   '    sk_root = "/home/<user>/slakos"',
                   '    sk_root = "/mnt/c/Users/<user>/slakos"     (inside WSL, C:\\Users\\<user> is /mnt/c/Users/<user>)',
                   '    sk_root = "C:/Users/<user>/slakos"         (Windows without WSL; "C:\\\\Users\\\\<user>\\\\slakos" also works)',
                   "  Fix the file and run again."]
        else:
            out = [f"環境設定ファイル {path} を読めません" + (f" ({line} 行目 {col} 文字目)" if line else "") + "。",
                   *([f"  {line} 行目:  {src.strip()}"] if src else []),
                   f"  {ja_hint + ' ' if ja_hint else ''}({msg})",
                   '  パスは "..." で囲み、\\ は / か \\\\ にします。例:',
                   '    sk_root = "/home/<ユーザー名>/slakos"',
                   '    sk_root = "/mnt/c/Users/<ユーザー名>/slakos"     (WSL の中では、Windows の C:\\Users\\<ユーザー名> は /mnt/c/Users/<ユーザー名> です)',
                   '    sk_root = "C:/Users/<ユーザー名>/slakos"         (WSL を使わない Windows。"C:\\\\Users\\\\<ユーザー名>\\\\slakos" でも可)',
                   "  直してから、もう一度実行してください。"]
        return "\n".join(out)
    if isinstance(ex, ValidationError):
        items = []
        for e in ex.errors():
            loc = ".".join(str(x) for x in e.get("loc", ()))
            items.append(f"  {loc or ('(whole file)' if en else '(全体)')}: {e.get('msg', '')}")
        head = (f"The settings file {path} is valid TOML, but these values cannot be used:" if en
                else f"環境設定ファイル {path} の書き方は正しいのですが、次の項目の値が使えません:")
        tail = ('  (true / false are written without quotes; kind is "direct", "pbs" or "slurm")' if en
                else '  (true / false は引用符なしで書きます。kind は "direct" / "pbs" / "slurm" のどれかです)')
        return "\n".join([head, *items, tail])
    return f"{path}: {ex}"


def read_config_text(path: Path) -> str:
    raw = path.read_bytes()
    try:
        return raw.decode("utf-8-sig")
    except UnicodeDecodeError as ex:
        en = _english()
        raise ConfigError(
            f"The settings file {path} is not saved as UTF-8 (byte {ex.start}). Save it again with the encoding UTF-8 "
            "(in Windows Notepad: File > Save as > Encoding: UTF-8)." if en else
            f"環境設定ファイル {path} が UTF-8 で保存されていません ({ex.start} バイト目)。文字コードを UTF-8 にして保存し直してください "
            "(Windows のメモ帳なら「名前を付けて保存」の「エンコード」で UTF-8 を選びます)。") from ex


def find_unknown_keys(data: dict) -> list[str]:
    known_top = set(Config.model_fields)
    known_profile = set(Profile.model_fields)
    out = [key for key in data if key not in known_top]
    profiles = data.get("profiles")
    if isinstance(profiles, dict):
        for name, body in profiles.items():
            if isinstance(body, dict):
                out += [f"profiles.{name}.{key}" for key in body if key not in known_profile]
    return sorted(out)


def load_config(path: Path | None = None) -> Config:
    path = path or config_path()
    if not path.is_file():
        en = _english()
        raise ConfigMissing(f"the settings file does not exist yet: {path}" if en else f"環境設定ファイルがまだありません: {path}")
    text = read_config_text(path)
    try:
        data = tomllib.loads(text)
        cfg = Config.model_validate(data)
        cfg._source_path = path.absolute()
        cfg._unknown_keys = find_unknown_keys(data)
        return cfg
    except (tomllib.TOMLDecodeError, ValidationError) as ex:
        raise ConfigError(describe_config_error(ex, text, path)) from ex


def ensure_config(path: Path | None = None) -> tuple[Config, Path, bool]:
    path = path or config_path()
    created = False
    if not path.exists():
        save_config(default_config(), path)
        created = True
    return load_config(path), path, created


def first_run_message(path: Path) -> str:
    from adit.lang import L

    return L(f"環境設定ファイルを作りました: {path}\n"
             "  DFTB+ を使うなら、このファイルの sk_root = \"\" の \"\" の中に、Slater-Koster パラメータを置いたフォルダを書いてください\n"
             "  (例: sk_root = \"/home/<ユーザー名>/slakos\"。その下に mio-1-1/ などのセットのフォルダがある場所)。\n"
             "  Quantum ESPRESSO を使うなら、同じように pseudo_root に UPF の置き場所を書きます。xtb だけなら書かなくてかまいません。",
             f"Created the settings file: {path}\n"
             "  For DFTB+, write the folder holding the Slater-Koster parameters between the quotes of sk_root = \"\"\n"
             "  (e.g. sk_root = \"/home/<user>/slakos\", the folder that contains set folders such as mio-1-1/).\n"
             "  For Quantum ESPRESSO, write the UPF folder in pseudo_root the same way. For xtb alone, nothing is needed.")


def set_top_level_value(path: Path, key: str, value) -> None:
    import re

    text = read_config_text(path)
    try:
        tomllib.loads(text)
    except tomllib.TOMLDecodeError as ex:
        raise ConfigError(describe_config_error(ex, text, path)) from ex
    lines = text.splitlines()
    _replace_top_level(lines, key, value)
    out = "\n".join(lines) + "\n"
    tomllib.loads(out)
    path.write_text(out, encoding="utf-8")


def _first_table(lines: list[str]) -> int:
    return next((i for i, s in enumerate(lines) if s.lstrip().startswith("[")), len(lines))


def _replace_top_level(lines: list[str], key: str, value) -> None:
    import re

    end = _first_table(lines)
    new_line = tomli_w.dumps({key: value}).strip()
    literal = new_line.split("=", 1)[1].strip()
    pat = re.compile(rf"""^(\s*{re.escape(key)}\s*=\s*)("[^"]*"|'[^']*'|[^\s#]+)(.*)$""")
    idx = next((i for i in range(end) if pat.match(lines[i])), None)
    if idx is not None:
        m = pat.match(lines[idx])
        lines[idx] = m.group(1) + literal + m.group(3)
    else:
        at = end
        while at > 0 and not lines[at - 1].strip():
            at -= 1
        lines.insert(at, new_line)


def _merged_config_text(data: dict, path: Path) -> str | None:
    # Keep the user's comments and unknown keys above the first table, and unknown keys inside
    # [profiles.*]; the tables themselves are rewritten from the model (comments in them are lost).
    if not path.is_file():
        return None
    try:
        text = read_config_text(path)
        old = tomllib.loads(text)
    except (OSError, tomllib.TOMLDecodeError):
        return None
    lines = text.splitlines()
    top = lines[:_first_table(lines)]
    tables = {k: v for k, v in data.items() if isinstance(v, dict)}
    if any(k in tomllib.loads("\n".join(top)) for k in tables):
        return None
    for key, value in data.items():
        if key not in tables:
            _replace_top_level(top, key, value)
    merged: dict = {name: body for name, body in old.items() if isinstance(body, dict) and name not in tables}
    known_profile = set(Profile.model_fields)
    for name, table in tables.items():
        if name == "profiles":
            old_profiles = old.get("profiles") if isinstance(old.get("profiles"), dict) else {}
            table = {pname: {**{k: v for k, v in (old_profiles.get(pname) or {}).items() if k not in known_profile}, **body}
                     for pname, body in table.items()}
        merged[name] = table
    out = "\n".join(top).rstrip("\n") + "\n\n" + tomli_w.dumps(merged)
    tomllib.loads(out)
    return out


def save_config(cfg: Config, path: Path | None = None) -> Path:
    path = path or config_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    data = cfg.model_dump(mode="json")
    text = _merged_config_text(data, path)
    path.write_text(text if text is not None else tomli_w.dumps(data), encoding="utf-8")
    cfg._source_path = path.absolute()
    return path


def unknown_keys_message(cfg: Config) -> str:
    keys = cfg.unknown_keys
    if not keys:
        return ""
    from adit.lang import L

    return L(f"注意: 設定ファイル {cfg.source_path} の次の項目は ADIT が知らないので使っていません: {', '.join(keys)}。"
             "[profiles.…] の見出しより下に書いた項目は、その節の中に入ります (TOML の決まり)。見出しより上へ移してください。",
             f"Note: these entries of the settings file {cfg.source_path} are unknown to ADIT and were not used: {', '.join(keys)}. "
             "In TOML, anything written below a [profiles.…] header belongs to that section; move them above the header.")
