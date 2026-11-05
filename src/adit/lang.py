"""Message language. L(ja, en) returns one of the two according to the current setting."""

from __future__ import annotations


from adit.config import env_var

LANGUAGE = "ja" if env_var("LANG").lower().startswith("ja") else "en"


def set_language(lang: str) -> None:
    global LANGUAGE
    LANGUAGE = "ja" if str(lang).lower().startswith("ja") else "en"


def L(ja: str, en: str) -> str:
    return en if LANGUAGE == "en" else ja


def pick(ja, en):
    return en if LANGUAGE == "en" else ja
