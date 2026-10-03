"""Message language. L(ja, en) returns one of the two according to the current setting."""

from __future__ import annotations


from adit.config import env_var

LANGUAGE = "en" if env_var("LANG").lower().startswith("en") else "ja"


def set_language(lang: str) -> None:
    global LANGUAGE
    LANGUAGE = "en" if str(lang).lower().startswith("en") else "ja"


def L(ja: str, en: str) -> str:
    return en if LANGUAGE == "en" else ja


def pick(ja, en):
    return en if LANGUAGE == "en" else ja
