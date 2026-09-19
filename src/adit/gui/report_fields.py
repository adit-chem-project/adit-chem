
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from adit.lang import L

LABELS: dict[str, tuple[str, str]] = {
    "report": ("報告をまとめる", "Build a report"),
    "rep_dirs": ("計算ディレクトリ (1 行に 1 つ)", "Run directories (one per line)"),
    "rep_lang": ("報告の言語", "Report language"),
    "rep_methods": ("方法の節 (Markdown)", "Methods section (Markdown)"),
    "rep_conditions": ("条件の表 (CSV)", "Conditions table (CSV)"),
    "rep_results": ("結果の表 (CSV)", "Results table (CSV)"),
    "rep_bundle": ("再現パッケージ (.zip)", "Reproducibility bundle (.zip)"),
    "rep_check": ("入力の照合用のハッシュを照合する", "Verify the input fingerprints"),
    "rep_run": ("報告を作る", "Build"),
}

LANGUAGES = [("", "(画面と同じ)", "(same as the screen)"), ("ja", "日本語", "Japanese"),
             ("en", "英語", "English"), ("both", "日本語と英語", "Japanese and English")]


class ReportFieldError(ValueError):
    pass


@dataclass
class ReportRequest:
    run_dirs: list[Path] = field(default_factory=list)
    language: str = ""
    methods_path: Path | None = None
    conditions_csv: Path | None = None
    results_csv: Path | None = None
    bundle: Path | None = None
    check: bool = False


@dataclass
class ReportOutcome:
    lines: list[str] = field(default_factory=list)
    written: list[Path] = field(default_factory=list)
    methods_text: str = ""
    check_verdict: str = ""                           # "ok" / "bad" / "unknown" / ""


def lab(key: str) -> str:
    ja, en = LABELS[key]
    return L(ja, en)


def _paths(text: str) -> list[Path]:
    out = []
    for line in str(text or "").replace(",", "\n").splitlines():
        s = line.strip()
        if s:
            out.append(Path(s).expanduser())
    return out


def request_from_fields(f: dict) -> ReportRequest:
    dirs = _paths(f.get("rep_dirs", ""))
    if not dirs:
        raise ReportFieldError(L(f"{lab('rep_dirs')}: 計算ディレクトリを 1 つ以上 入れてください",
                                 f"{lab('rep_dirs')}: give at least one run directory"))
    missing = [str(d) for d in dirs if not d.is_dir()]
    if missing:
        raise ReportFieldError(L(f"{lab('rep_dirs')}: ディレクトリがありません: {', '.join(missing)}",
                                 f"{lab('rep_dirs')}: no such directory: {', '.join(missing)}"))
    language = str(f.get("rep_lang", "") or "")
    if language not in ("", "ja", "en", "both"):
        raise ReportFieldError(L(f"{lab('rep_lang')}: ja / en / both のどれかにしてください",
                                 f"{lab('rep_lang')}: choose ja, en or both"))

    def one(key: str) -> Path | None:
        text = str(f.get(key, "") or "").strip()
        return Path(text).expanduser() if text else None

    return ReportRequest(run_dirs=dirs, language=language, methods_path=one("rep_methods"),
                         conditions_csv=one("rep_conditions"), results_csv=one("rep_results"),
                         bundle=one("rep_bundle"),
                         check=bool(f.get("rep_check") in (True, "on", "1", "true", "yes")))


def build(req: ReportRequest) -> ReportOutcome:
    from adit.report import (ReportError, check_lines, conditions_csv, load_run_report,
                              methods_markdown, results_csv, write_bundle)

    out = ReportOutcome()
    for path in (req.methods_path, req.conditions_csv, req.results_csv):
        if path is not None and path.exists():
            raise ReportError(L(f"すでにあります: {path} (別の名前を指定してください。上書きはしません)",
                                f"already exists: {path} (choose another name; nothing is overwritten)"))
    reports = [load_run_report(d) for d in req.run_dirs]
    out.methods_text = methods_markdown(reports, req.language)
    if req.methods_path is not None:
        req.methods_path.parent.mkdir(parents=True, exist_ok=True)
        req.methods_path.write_text(out.methods_text, encoding="utf-8")
        out.written.append(req.methods_path)
    for path, text in ((req.conditions_csv, conditions_csv), (req.results_csv, results_csv)):
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(text(reports), encoding="utf-8")
            out.written.append(path)
    if req.bundle is not None:
        try:
            out.written.append(write_bundle(reports, req.bundle, language=req.language))
        except ReportError as ex:
            out.lines.append(str(ex))
    if req.check:
        verdicts = []
        for report in reports:
            lines, verdict = check_lines(report)
            out.lines += lines
            verdicts.append(verdict)
        out.check_verdict = ("bad" if "bad" in verdicts else "unknown" if "unknown" in verdicts else "ok")
        out.lines.append({"ok": L("照合用のハッシュ: すべて一致しました", "fingerprints: everything matches"),
                          "bad": L("照合用のハッシュ: 生成したときと違うファイルがあります", "fingerprints: some files differ from when they were generated"),
                          "unknown": L("照合用のハッシュ: 記録が無くて確かめられないものがあります",
                                       "fingerprints: some files have no record and cannot be verified")}[out.check_verdict])
    if out.written:
        out.lines.append(L("書き出したファイル: " + ", ".join(str(p) for p in out.written),
                           "files written: " + ", ".join(str(p) for p in out.written)))
    return out
