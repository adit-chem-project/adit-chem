
from __future__ import annotations

import argparse
import email
import email.policy
import hmac
import html
import io
import json
import os
import re
import secrets
import shutil
import socket

from adit.compat import ensure_printable_stdio, upload_basename
import sys
import threading
import webbrowser
from datetime import datetime
from http import HTTPStatus
from http.cookies import SimpleCookie
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, quote, urlencode, urlparse

from jinja2 import Environment, PackageLoader, StrictUndefined

from adit import __version__, lang
from adit.analysis import AnalysisOptions, run_analysis
from adit.codes.sk_sets import discover_sets
from adit.config import (Config, ConfigError, config_path, ensure_config, env_var, first_run_message, load_config, save_config,
                          set_top_level_value)
from adit.gui.help import help_for
from adit.lang import L
from adit.project import OutputNotEmpty, ProjectError, ProjectFiles, build_project, load_project, write_project
from adit.spec import CalculationSpec
from adit.structure import CRYSTAL_STRUCTURES, FETCH_DATABASES, SURFACE_FUNCTIONS, has_rdkit, preset_names, preset_search_text
from adit.web import forms
from adit.web.forms import FormError, default_form, form_from_spec, spec_from_form

from ase.data import chemical_symbols  # noqa: E402

ELEMENTS = list(chemical_symbols[1:104])
_env = Environment(loader=PackageLoader("adit.web", "templates"), undefined=StrictUndefined, autoescape=True)
from adit.analysis.report import figure_title  # noqa: E402  
_env.globals["figure_title"] = figure_title
from adit.gui import analysis_fields as AF  # noqa: E402  
def _label_table(key: str) -> tuple[str, str]:
    from adit.gui import report_fields as _RF

    return AF.LABELS[key] if key in AF.LABELS else _RF.LABELS[key]


_env.globals["lab"] = lambda key: L(*_label_table(key))
_env.globals["hlp"] = lambda key: help_for(_label_table(key)[0])

RESULT_FILES = ("output.log", "md.out", "detailed.out", "geo_end.xyz", "hessian.out", "band.out", "vasprun.xml", "OSZICAR", "OUTCAR",
                "xtbopt.xyz", "xtb.trj", "vibspectrum", "trajectory.xyz", "dynmat.out", "bands/output.log")


def _mtime(path: Path) -> tuple[int, int] | None:
    try:
        st = path.stat()
        return (st.st_mtime_ns, st.st_size)
    except OSError:
        return None


class WebApp:

    def __init__(self, cfg: Config, cfg_path: Path):
        self.cfg, self.cfg_path = cfg, cfg_path
        # Requests run one at a time per page group, because the form, spec, figures and recipe caches are shared.
        # Two locks so the analysis pages still answer while a slow recipe build holds the prepare page.
        self.lock = threading.RLock()
        self.analysis_lock = threading.RLock()
        self._cfg_mtime = _mtime(cfg_path)
        self.config_error = ""
        self.form: dict[str, str] = default_form(cfg.default_profile, Path.home() / "adit_runs" / f"run_{datetime.now():%Y%m%d_%H%M%S}")
        self.spec: CalculationSpec | None = None
        self.files: ProjectFiles | None = None
        self.errors: list = []
        self.written: Path | None = None
        self.written_kind: str = ""
        self.written_exe: str = ""
        self.written_code: str = ""
        self.figures: dict[str, str] = {}
        self.upload_dir = Path.home() / "adit_runs" / "uploads"
        self.fetch_dir = Path.home() / "adit_runs" / "fetched"
        self.fetch_candidates: list = []
        self.scan_out: Path | None = None
        self._scan_auto = ""
        self.built = None
        self.recipe_status: tuple[str, str] | None = None
        self.recipe_open: int | None = None
        self.recipe_note: tuple[int, str, bool] | None = None
        from adit.web.codefields import PrepOrigin
        self.origin = PrepOrigin()
        self.origin_structure = None
        self.stages_out: Path | None = None
        self._stages_auto = ""
        self.batch_out = None
        self._batch_auto = ""

    def refresh_config(self) -> None:
        m = _mtime(self.cfg_path)
        if m is None or m == self._cfg_mtime:
            return
        self._cfg_mtime = m
        try:
            cfg = load_config(self.cfg_path)
        except (ConfigError, OSError) as ex:
            self.config_error = str(ex)
            return
        self.cfg, self.config_error = cfg, ""
        lang.set_language(env_var("LANG", cfg.language))

    def set_language(self, new: str) -> None:
        lang.set_language(new)
        self.cfg.language = new
        try:
            if self.cfg_path.is_file():
                set_top_level_value(self.cfg_path, "language", new)
            else:
                save_config(self.cfg, self.cfg_path)
        except (OSError, ConfigError):
            return
        self._cfg_mtime = _mtime(self.cfg_path)

    def choices(self) -> dict:
        sets = discover_sets(Path(self.cfg.sk_root).expanduser()) if self.cfg.sk_root else {}
        return {
            "codes": forms.CODES, "types": forms.TYPES, "relax_cell": forms.RELAX_CELL, "thermostats": forms.THERMOSTATS,
            "kp_modes": forms.KP_MODES, "gfn": forms.GFN, "ibrion": forms.IBRION, "presets": preset_names(), "has_rdkit": has_rdkit(),
            "crystal_structures": CRYSTAL_STRUCTURES, "surfaces": SURFACE_FUNCTIONS, "sk_sets": sorted(sets), "profiles": self.cfg.profiles,
            "optimizers": ["Rational", "LBFGS", "FIRE", "SteepestDescent"], "ensembles": ["NVT", "NVE", "NPT"],
            "opt_levels": ["crude", "sloppy", "loose", "normal", "tight", "verytight", "extreme"],
            "orca_methods": ["HF", "PBE", "B3LYP", "PBE0", "TPSS", "wB97X-D3", "MP2", "DLPNO-CCSD(T)"],
            "orca_basis": ["def2-SVP", "def2-TZVP", "def2-TZVPP", "def2-QZVPP", "cc-pVDZ", "cc-pVTZ", "ma-def2-SVP"],
        }

    def render(self, name: str, **ctx) -> str:
        tmpl = _env.get_template(name)
        return tmpl.render(T=L, HELP=help_for, lang=lang.LANGUAGE, version=__version__, cfg=self.cfg, cfg_path=str(self.cfg_path), app=self, **ctx)

    def scene(self) -> tuple[str, str, bool]:
        from adit.web.structure3d import scene_from_atoms, summary_line

        if self.spec is None:
            return "", "", False
        try:
            atoms = self.spec.structure.atoms.to_ase()   # AtomsData → ase.Atoms
            scene = scene_from_atoms(atoms)
            periodic = bool(getattr(self.spec.structure, "periodic", False))
        except Exception:
            return "", "", False
        if scene.n_atoms == 0:
            return "", "", False
        return scene.to_json(), summary_line(scene, periodic), periodic

    @staticmethod
    def sketch_js() -> str:
        from pathlib import Path as _Path

        path = _Path(__file__).with_name("static") / "sketch.js"
        try:
            return path.read_text(encoding="utf-8")
        except OSError:
            return ""

    @staticmethod
    def viewer_js() -> str:
        from pathlib import Path as _Path

        path = _Path(__file__).with_name("static") / "viewer3d.js"
        try:
            return path.read_text(encoding="utf-8")
        except OSError:
            return ""

    def preview(self, form: dict[str, str]) -> tuple[str, str]:
        self.form = {**default_form(self.cfg.default_profile, self.form.get("output_dir", "")), **form}
        self.form.pop("recipe_action", None)
        self._fill_step_fields(self.form)
        self._sync_recipe_charge(self.form)
        self.spec, self.files, self.errors = None, None, []
        try:
            spec = spec_from_form(self.form, self)
        except FormError as ex:
            return "", str(ex)
        except Exception as ex:
            return "", _unexpected(ex)
        spec = self.origin.apply(spec)
        self.spec = spec
        if self.written is not None and spec.method.code != self.written_code:
            self.written, self.written_kind, self.written_exe, self.written_code = None, "", "", ""
        try:
            self.files = build_project(spec, self.cfg, output_dir=self.form.get("output_dir") or None)
        except ProjectError as ex:
            self.errors = list(ex.errors)
            return "", str(ex)
        except Exception as ex:
            self.spec = None
            return "", _unexpected(ex)
        return L("生成できます", "Ready to generate"), ""

    def fetch_structure(self, form: dict[str, str]) -> tuple[str, str, str]:
        """Fetch on the user's request only; the chosen entry then feeds the normal file route. Returns (status, error, message)."""
        from adit.fetch import FetchError, candidates_text, fetch, summary_line

        self.form = {**default_form(self.cfg.default_profile, self.form.get("output_dir", "")), **form}
        self.form["source"] = "fetch"
        pick = (form.get("fetch_pick") or "").strip()
        db, query = (form.get("fetch_db") or "pubchem").strip(), (form.get("fetch_query") or "").strip()
        if pick and self.fetch_candidates and any(c.ref == pick for c in self.fetch_candidates):
            ref = pick
        elif query:
            ref = f"{db}:{query}"
        else:
            return "", L("データベースから取得: 名前か ID を入力してください", "Fetch from a database: enter a name or ID"), ""
        try:
            result = fetch(ref, mp_api_key=self.cfg.mp_api_key)
        except FetchError as ex:
            return "", str(ex), ""
        if result.fetched is None:
            self.fetch_candidates = list(result.candidates)
            return "", candidates_text(result), ""
        self.fetch_candidates = []
        result.fetched.save(self.fetch_dir)
        self.form["fetch_file"] = str(result.fetched.record["file"])
        self.form["fetch_record"] = json.dumps(result.fetched.record, ensure_ascii=False)
        self.form["file_sha_path"] = self.form["file_sha256"] = ""
        status_text, err = self.preview(self.form)
        note = L(f"構造を取得しました: {summary_line(result.fetched.record)}", f"Fetched the structure: {summary_line(result.fetched.record)}")
        if result.notes:
            note += "\n" + "\n".join(result.notes)
        return status_text, err, note

    def fetch_info(self) -> str:
        from adit.fetch import summary_line
        from adit.web.forms import fetch_record

        rec = fetch_record(self.form)
        return summary_line(rec) if rec else ""

    def generate(self, overwrite: bool) -> tuple[str, str]:
        if self.spec is None:
            return "", L("先にプレビューしてください", "Preview first")
        out = self.form.get("output_dir", "").strip()
        if not out:
            return "", L("出力ディレクトリを指定してください。", "specify an output directory.")
        try:
            written = write_project(self.spec, self.cfg, out, overwrite=overwrite)
        except OutputNotEmpty as ex:
            return "", str(ex) + "\n" + L("中のファイルを上書きしてよければ、「上書きを許可」に印を付けてから、もう一度「生成」を押してください。",
                                          "To overwrite the files inside, tick \"Allow overwrite\" and press \"Generate\" again.")
        except ProjectError as ex:
            return "", str(ex)
        self.written = Path(out).expanduser().resolve()
        profile = self.cfg.profile(self.spec.runtime.profile)
        self.written_kind = profile.kind
        self.written_code = self.spec.method.code
        from adit.codes import GENERATORS
        from adit.runner import executable_of

        self.written_exe = executable_of(GENERATORS[self.spec.method.code].run_command(self.spec, profile))
        return L(f"{len(written)} ファイルを {self.written} に書きました。実行の手順は README.txt", f"wrote {len(written)} files to {self.written}. See README.txt"), ""

    def scan_periodic_guess(self) -> bool:
        if self.spec is not None:
            return self.spec.structure.periodic
        return self.form.get("source") in ("bulk", "surface", "mixture", "2d") or forms._b(self.form.get("box"))

    def scan_default_dir(self) -> str:
        base = (self.form.get("output_dir") or "").strip().rstrip("/")
        return base + "_scan" if base else ""

    def scan_dir_value(self) -> str:
        v = (self.form.get("scan_dir") or "").strip()
        auto = self.scan_default_dir()
        if not v or v == self._scan_auto:
            v = auto
        self._scan_auto = auto
        return v

    def scan_item_value(self) -> str:
        from adit.web.scan_choices import choices_for
        paths = [c.path for c in choices_for(self.form.get("code", "dftbplus"), self.scan_periodic_guess())]
        cur = self.form.get("scan_item", "")
        return cur if cur in paths else paths[0]

    def scan(self, overwrite: bool) -> tuple[str, str]:
        from pydantic import ValidationError as PydanticError
        from adit.scan import ScanError, parse_scan, write_scan
        from adit.web.scan_choices import OTHER, all_choices

        if self.spec is None:
            return "", L("先にプレビューしてください", "Preview first")
        f = self.form
        item = (f.get("scan_item") or "").strip()
        if item == OTHER:
            path = (f.get("scan_path") or "").strip()
            if not path:
                return "", L("項目の場所を書いてください (例 method.scc_tolerance)。", "Type the setting path (e.g. method.scc_tolerance).")
        elif item in [c.path for c, _, _ in all_choices()]:
            path = item
        else:
            return "", L("変える項目を選んでください。", "Choose the parameter to change.")
        folder = (f.get("scan_dir") or "").strip()
        if not folder:
            return "", L("保存先を指定してください。", "Choose an output directory.")
        out = Path(folder).expanduser()
        try:
            scan = parse_scan(f"{path}={f.get('scan_values') or ''}")
            if out.exists() and not out.is_dir():
                return "", L(f"{out} は保存先にできません (同じ名前のファイルがあります)。", f"{out} cannot be used as the output directory (a file with that name exists).")
            if out.is_dir() and any(out.iterdir()) and not overwrite:
                return "", L(f"{out} は空ではありません。中のファイルを上書きしてよければ、一括生成の欄の「上書きを許可」に印を付けてから、もう一度「生成」を押してください。",
                             f"{out} is not empty. To overwrite the files inside, tick \"Allow overwrite\" in the parameter scan box and press \"Generate\" again.")
            dirs = write_scan(self.spec, self.cfg, out, scan, overwrite=overwrite)
        except (ScanError, ProjectError, ConfigError, ValueError, OSError) as ex:
            if isinstance(ex, PydanticError):
                from adit.validate_types import friendly_pydantic
                return "", friendly_pydantic(ex)
            return "", str(ex)
        except Exception as ex:
            return "", _unexpected(ex)
        self.scan_out = out.resolve()
        names = "\n".join(f"  {d.name}" for d in dirs)
        return L(
            f"{out} の中に {len(dirs)} 個のディレクトリを作りました。\n{names}\n\n"
            "次にすること\n"
            "1. それぞれの計算を実行します。この PC なら各ディレクトリで bash submit.sh、"
            "クラスタならジョブとして投入します (手順は各ディレクトリの README.txt にあります)。\n"
            "2. すべて終わったら、解析ページで「解析を実行」を押します。下の「解析へ進む」から開くと、計算結果のディレクトリ欄にこの保存先が入ります。",
            f"Created {len(dirs)} directories in {out}:\n{names}\n\n"
            "Next steps\n"
            "1. Run each calculation: on this PC, run bash submit.sh in each directory; on a cluster, submit each one as a job "
            "(see README.txt in each directory).\n"
            "2. When they have all finished, click \"Run analysis\" on the Analysis page. Open it with \"Go to analysis\" below "
            "and the Run directory field is already filled in."), ""

    def adopt_origin(self, spec: CalculationSpec) -> None:
        from adit.web.codefields import PrepOrigin
        self.origin = PrepOrigin(spec)
        self.origin_structure = spec.structure if (spec.meta.continued_from or spec.handoff) else None

    def clear_origin(self) -> None:
        from adit.web.codefields import PrepOrigin
        self.origin, self.origin_structure = PrepOrigin(), None

    def continue_run(self, form: dict[str, str]) -> tuple[str, str, str]:
        from adit.continuation import continue_from
        from adit.web.codefields import continuation_summary
        from adit.web.recipe_form import write_steps

        self.form = {**default_form(self.cfg.default_profile, self.form.get("output_dir", "")), **form}
        d = (self.form.get("cont_dir") or "").strip()
        if not d:
            return "", L("前の計算のディレクトリを指定してください。", "Choose the previous run directory."), ""
        cond = None
        if forms._b(self.form.get("cont_keep")):
            _, err = self.preview(dict(self.form))
            if self.spec is None:
                return "", L("いまの画面の条件を組み立てられません: ", "cannot build the current settings: ") + err, ""
            cond = self.spec
        try:
            spec = continue_from(d, cond, velocities=forms._b(self.form.get("cont_vel")))
        except (ValueError, OSError) as ex:
            return "", forms._pydantic_text(ex) if isinstance(ex, ValueError) else str(ex), ""
        keep = {k: self.form.get(k, "") for k in ("cont_dir", "cont_vel", "cont_keep")}
        self.adopt_origin(spec)
        self.form = {**forms.drop_prep_fields(self.form), **form_from_spec(spec), **keep}
        write_steps(self.form, [])
        prev = Path(spec.meta.continued_from["dir"])
        out = (self.form.get("output_dir") or "").strip()
        if out and Path(out).expanduser().resolve() == prev.resolve():
            self.form["output_dir"] = str(prev) + "_cont"
        status, err = self.preview(dict(self.form))
        return status, err, continuation_summary(spec) + "\n\n" + L(
            "出力ディレクトリは、前の計算とは別の場所にしてください (同じ場所に書くと前の出力を上書きします)。",
            "Use an output directory other than the previous run's (writing there would overwrite its output).")

    def batch_default_dir(self) -> str:
        base = (self.form.get("output_dir") or "").strip().rstrip("/")
        suffix = {"compare": "_set", "conformers": "_conf", "neb": "_neb", "phonons": "_phonon", "elastic": "_elastic",
                  "ts": "_ts"}.get(self.form.get("batch_kind") or "compare", "_batch")
        return base + suffix if base else ""

    def batch_dir_value(self) -> str:
        v = (self.form.get("batch_dir") or "").strip()
        auto = self.batch_default_dir()
        if not v or v == self._batch_auto:
            v = auto
        self._batch_auto = auto
        return v

    def batch_rx_rows(self) -> list[dict[str, str]]:
        from adit.web import prep23 as P
        return P.reaction_rows(self.form)

    def batch(self, overwrite: bool) -> tuple[str, str]:
        from adit.web import prep23 as P

        if self.spec is None:
            return "", L("先にプレビューしてください", "Preview first")
        kind = (self.form.get("batch_kind") or "compare").strip()
        folder = (self.form.get("batch_dir") or "").strip() or self.batch_default_dir()
        try:
            out = Path(folder).expanduser() if folder else Path("")
            if not folder:
                return "", L("保存先を指定してください。", "Choose an output directory.")
            if out.exists() and not out.is_dir():
                return "", L(f"{out} は保存先にできません (同じ名前のファイルがあります)。", f"{out} cannot be used as the output directory (a file with that name exists).")
            if out.is_dir() and any(out.iterdir()) and not overwrite:
                return "", L(f"{out} は空ではありません。中のファイルを上書きしてよければ、「まとめて作る」の欄の「上書きを許可」に印を付けてから、もう一度「生成」を押してください。",
                             f"{out} is not empty. To overwrite the files inside, tick \"Allow overwrite\" in the batch generation box and press \"Generate\" again.")
            res = P.run_batch(kind, self.spec, self.cfg, out, self.form, overwrite=overwrite)
        except ImportError as ex:
            return "", L(f"必要なパッケージがありません: {ex}", f"a required package is missing: {ex}")
        except (ProjectError, ConfigError, ValueError, OSError) as ex:
            return "", forms._pydantic_text(ex) if isinstance(ex, ValueError) else str(ex)
        except Exception as ex:
            return "", _unexpected(ex)
        self.batch_out = res
        return P.result_message(res), ""

    def put_doc_value(self, text: str) -> str:
        from adit.web.prep23 import DOC_FIELDS

        field, _, value = (text or "").partition("=")
        field, value = field.strip(), value.strip()
        fields = {f for _g, f, _u in DOC_FIELDS.values() if f}
        try:
            float(value)
        except ValueError:
            return ""
        if field not in fields:
            return ""
        self.form[field] = value
        return L(f"{field} に {value} を入れました (文書に書かれていた値です)。",
                 f"Put {value} in {field} (the value written in the documentation).")

    def stages_default_dir(self) -> str:
        base = (self.form.get("output_dir") or "").strip().rstrip("/")
        return base + "_stages" if base else ""

    def stages_dir_value(self) -> str:
        v = (self.form.get("stg_dir") or "").strip()
        auto = self.stages_default_dir()
        if not v or v == self._stages_auto:
            v = auto
        self._stages_auto = auto
        return v

    def stage_rows(self) -> list[dict]:
        from adit.web.codefields import STAGE_COLUMNS
        try:
            n = max(1, min(forms._i(self.form.get("stg_rows"), 3), 30))
        except FormError:
            n = 3
        return [{c: self.form.get(f"stg{i}_{c}", "auto" if c == "velocities" else "") or ("auto" if c == "velocities" else "")
                 for c in STAGE_COLUMNS} for i in range(1, n + 1)]

    def stages(self, overwrite: bool) -> tuple[str, str]:
        from adit.stages import StageError, parse_stages, write_stages
        from adit.web.codefields import stage_from_row

        if self.spec is None:
            return "", L("先にプレビューしてください", "Preview first")
        folder = (self.form.get("stg_dir") or "").strip() or self.stages_default_dir()
        if not folder:
            return "", L("保存先を指定してください。", "Choose an output directory.")
        out = Path(folder).expanduser()
        try:
            stages = []
            for row in self.stage_rows():
                s = stage_from_row(row, len(stages) + 1, self.spec.task.type)
                if s is not None:
                    stages.append(s)
            parsed = parse_stages({"stages": stages})
            if out.exists() and not out.is_dir():
                return "", L(f"{out} は保存先にできません (同じ名前のファイルがあります)。", f"{out} cannot be used as the output directory (a file with that name exists).")
            if out.is_dir() and any(out.iterdir()) and not overwrite:
                return "", L(f"{out} は空ではありません。中のファイルを上書きしてよければ、段階の欄の「上書きを許可」に印を付けてから、もう一度「段階に分けて生成」を押してください。",
                             f"{out} is not empty. To overwrite the files inside, tick \"Allow overwrite\" in the staged calculation box and press \"Generate stages\" again.")
            dirs = write_stages(self.spec, self.cfg, out, parsed, overwrite=overwrite)
        except (StageError, ProjectError, ConfigError, ValueError, OSError) as ex:
            return "", forms._pydantic_text(ex) if isinstance(ex, ValueError) else str(ex)
        except Exception as ex:
            return "", _unexpected(ex)
        self.stages_out = out.resolve()
        names = "\n".join(f"  {d.name}" for d in dirs)
        return L(f"{out} の中に {len(dirs)} 段階のディレクトリを作りました。\n{names}\n\n"
                 "次にすること\n"
                 "1. この PC なら、保存先で bash submit.sh を実行すると段階を順に実行します。クラスタなら、保存先の README.txt にある "
                 "依存付きの投入の例 (前の段階が終わってから次の段階が始まる) を使います。\n"
                 "2. 各段階が終わったら、解析ページでその段階のディレクトリを解析します。",
                 f"Created {len(dirs)} stage directories in {out}:\n{names}\n\n"
                 "Next steps\n"
                 "1. On this PC, run bash submit.sh in the output directory to run the stages in order. On a cluster, use the dependent "
                 "submission example in README.txt there (each stage starts after the previous one finishes).\n"
                 "2. When a stage has finished, analyze its directory on the Analysis page."), ""

    def template_load(self, form: dict[str, str]) -> tuple[str, str, str]:
        from adit.templates import TemplateError, load_template

        self.form = {**default_form(self.cfg.default_profile, self.form.get("output_dir", "")), **form}
        name = (self.form.get("tmpl_name") or "").strip()
        if not name:
            return "", L("雛形を選んでください。", "Choose a template."), ""
        try:
            st = spec_from_form(self.form, self).structure
        except FormError as ex:
            return "", L("先に構造を作ってください (雛形は構造を持たないため): ", "make a structure first (templates have no structure): ") + str(ex), ""
        try:
            t = load_template(name, st, self.cfg)
        except TemplateError as ex:
            return "", str(ex), ""
        keep = {k: v for k, v in self.form.items() if k.startswith(("tmpl_", "cont_", "stg", "scan_")) or k == "output_dir"}
        self.form = {**forms.drop_prep_fields(self.form), **form_from_spec(t), **keep}
        self.origin.take(t, template_only=True)
        status, err = self.preview(dict(self.form))
        return status, err, L(f"雛形 {name} を読み込みました (構造はそのまま)。雛形の値のままの欄には、ラベルの下に小さく印を付けています。",
                              f"Loaded the template {name} (structure kept). Fields still at the template value are marked under their labels.")

    def template_save(self, overwrite: bool) -> tuple[str, str]:
        from adit.templates import TemplateError, save_template

        if self.spec is None:
            return "", L("先にプレビューしてください", "Preview first")
        try:
            p = save_template(self.spec, (self.form.get("tmpl_save_name") or "").strip(), cfg=self.cfg,
                              comment=(self.form.get("tmpl_comment") or "").strip(), overwrite=overwrite)
        except (TemplateError, OSError) as ex:
            msg = str(ex)
            if not overwrite and ("すでにあります" in msg or "already exists" in msg):
                msg += "\n" + L("上書きしてよければ、雛形の欄の「上書きを許可」に印を付けてから、もう一度保存してください。",
                                "To overwrite it, tick \"Allow overwrite\" in the template box and save again.")
            return "", msg
        return L(f"雛形を保存しました: {p}", f"saved the template: {p}"), ""

    @staticmethod
    def _fill_step_fields(f: dict) -> None:
        from adit.web import recipe_form as R
        try:
            steps = R.raw_steps(f)
        except R.RecipeFormError:
            return
        for k, st in enumerate(steps, start=1):
            if f.get(f"st{k}_op") != st.op:
                f.update(R.step_fields(k, st))

    @staticmethod
    def _sync_recipe_charge(f: dict) -> None:
        from adit.web import recipe_form as R
        if not R.recipe_mode(f):
            f["recipe_q"] = ""
            return
        try:
            q = str(R.recipe_from_form(f).total_charge())
        except Exception:
            return
        if q != (f.get("recipe_q") or ""):
            f["recipe_q"] = f["charge"] = q

    def recipe_action(self, form: dict[str, str]) -> tuple[str, str]:
        from adit.builder import Recipe, recipe_structure
        from adit.structure import StructureError
        from adit.web import recipe_form as R

        self.form = {**default_form(self.cfg.default_profile, self.form.get("output_dir", "")), **form}
        f = self.form
        action = (f.pop("recipe_action", "") or "").strip()
        self._fill_step_fields(f)
        self.recipe_status, self.recipe_note = None, None
        verb, _, arg = action.partition(":")
        k = int(arg) if arg.isdecimal() else 0
        try:
            steps = R.steps_from_form(f)
        except StructureError as ex:
            self.recipe_status, self.recipe_open = ("ng", str(ex)), R.failed_step(str(ex))
            return "", ""
        keep: dict[str, str] = {}
        if verb == "add":
            op = (f.get("recipe_add") or "").strip()
            if op == R.INTERFACE:
                new = R.interface_steps(R.base_is_slab(f))
                steps += new
                self.recipe_open = len(steps) - 1
                self.recipe_status = ("hint", L(*R.INTERFACE_NOTE))
            elif op in R.ADD_ORDER:
                steps.append(R.default_step(op)); self.recipe_open = len(steps)
        elif verb in ("up", "down") and 1 <= k <= len(steps):
            t = k - 1 + (-1 if verb == "up" else 1)
            if 0 <= t < len(steps):
                steps[k - 1], steps[t] = steps[t], steps[k - 1]
                self.recipe_open = t + 1
        elif verb == "del" and 1 <= k <= len(steps):
            steps.pop(k - 1)
            self.recipe_open = min(k, len(steps)) or None
        elif verb == "conc" and 1 <= k <= len(steps) and steps[k - 1].op == "solvent_layer":
            keep = {n: f.get(n, "") for n in (f"st{k}_conc", f"st{k}_conc_row")}
            self.recipe_open = k
            try:
                conc = forms._f(f.get(f"st{k}_conc"), 1.0); row = forms._i(f.get(f"st{k}_conc_row"), 0)
            except FormError as ex:
                self.recipe_note = (k, str(ex), False)
            else:
                area = None
                if steps[k - 1].thickness is not None:
                    try:
                        area = R.area_before(Recipe(base=R.base_from_form(f), steps=steps), k - 1, self.built)
                    except StructureError:
                        area = None
                steps[k - 1], msg, ok = R.fill_count(steps[k - 1], row, conc, area)
                self.recipe_note = (k, msg, ok)
        elif verb == "build":
            pass
        R.write_steps(f, steps)
        f.update({n: v for n, v in keep.items() if v})
        self._sync_recipe_charge(f)
        self.spec, self.files, self.errors = None, None, []
        if verb != "build":
            return "", ""
        try:
            rec = Recipe(base=R.base_from_form(f), steps=steps)
            st, logs = recipe_structure(rec, charge=forms._i(f.get("charge"), 0), multiplicity=forms._i(f.get("multiplicity"), 1))
        except (StructureError, FormError) as ex:
            self.built, self.recipe_status = None, ("ng", str(ex))
            self.recipe_open = R.failed_step(str(ex)) or self.recipe_open
            return "", ""
        except Exception as ex:
            self.built = None
            self.recipe_status = ("ng", L(f"作れませんでした ({type(ex).__name__}: {ex})", f"could not build ({type(ex).__name__}: {ex})"))
            return "", ""
        self.built = R.Built(rec.to_ref(), st, [x.line() for x in logs])
        self.recipe_status = ("ok", L("できました", "Done"))
        return self.preview(dict(f))

    def adopt_loaded(self, spec: CalculationSpec) -> None:
        from adit.web import recipe_form as R
        self.built, self.recipe_status, self.recipe_note, self.recipe_open = None, None, None, None
        if spec.structure.source != "recipe":
            return
        try:
            ref = R.recipe_from_form(self.form).to_ref()
        except Exception:
            return
        self.built = R.Built(ref, spec.structure, [])
        self.recipe_status = ("hint", L("spec.json の原子座標を使っています。手順から作り直すなら「作る」を押してください",
                                        "using the coordinates saved in spec.json; press Build to rebuild from the steps"))
        self.recipe_open = 1 if R.step_count(self.form) else None

    def recipe_view(self) -> dict:
        from adit.builder import op_label
        from adit.structure import StructureError
        from adit.web import recipe_form as R

        f = self.form
        err = ""
        try:
            steps = R.steps_from_form(f)
        except StructureError as ex:
            err = str(ex)
            try:
                steps = R.raw_steps(f)
            except StructureError:
                steps = []
        rec, current = None, ""
        try:
            from adit.builder import Recipe
            rec = Recipe(base=R.base_from_form(f), steps=steps)
            current = rec.to_ref()
        except (StructureError, ValueError):
            rec = None
        items = []
        for k, s in enumerate(steps, start=1):
            it = {"k": k, "op": s.op, "label": op_label(s.op), "summary": R.summary(s), "p": f"st{k}_", "note": "", "note_ok": True,
                  "auto_fit": s.op == "supercell" and bool(s.fit_components),
                  "terms": [], "rows": []}
            if s.op == "slab":
                names = R.terminations(rec, k - 1) if rec is not None else None
                raw_term = (f.get(f"st{k}_term") or "").strip()
                try:
                    term = forms._i(raw_term, 0) if raw_term.lstrip("-").isdecimal() else s.termination
                except FormError:
                    term = s.termination
                it["terms"] = R.term_options(names, term)
            if s.op == "solvent_layer":
                it["rows"] = [(str(i), f"{i + 1}: {R.pretty_name(c.name())}") for i, c in enumerate(s.components)]
            if self.recipe_note and self.recipe_note[0] == k:
                it["note"], it["note_ok"] = self.recipe_note[1], self.recipe_note[2]
            items.append(it)
        built_ok = self.built is not None and current and self.built.ref == current
        status = self.recipe_status
        if err and not (status and status[0] == "ng"):
            status = ("ng", err)
        if status is None or (status[0] == "ok" and not built_ok):
            if built_ok:
                status = ("ok", L("できました", "Done"))
            elif R.needs_build(f):
                status = ("hint", L("組み立て手順を変えました。「作る」を押すと作り直します", "the steps have changed; press Build to rebuild") if self.built is not None
                          else L("「作る」を押すと、組み立て手順から構造を作ります", "press Build to make the structure from the steps"))
        lines = self.built.lines if (self.built is not None and R.recipe_mode(f)) else []
        open_k = self.recipe_open if self.recipe_open and 1 <= self.recipe_open <= len(items) else (1 if items else None)
        if status and status[0] == "ng" and R.failed_step(status[1]):
            open_k = R.failed_step(status[1])
        return {"items": items, "status": status, "lines": lines, "open": open_k, "has_fix": any(s.op == "fix" for s in steps),
                "add_opts": [(op, op_label(op)) for op in R.ADD_ORDER] + [(R.INTERFACE, L("電極と電解質の界面 (面 + 断面の自動調整 + 溶液 + 固定)",
                                                                                       "Electrode–electrolyte interface (surface cut + auto-sized cross-section + solution + fixed layer)"))]}



    def analyze(self, run_dir: str, opts: AnalysisOptions):
        res = run_analysis(run_dir, opts)
        self.figures = dict(res.figures)
        return res


def _is_md(run_dir: str) -> bool:
    try:
        return load_project(run_dir).task.type == "molecular_dynamics"
    except Exception:
        return False


def _is_scan_dir(d: str) -> bool:
    from adit.scan import SCAN_FILE
    try:
        return bool((d or "").strip()) and (Path(d.strip()).expanduser() / SCAN_FILE).is_file()
    except (OSError, ValueError):
        return False


COMPARE_MIN_ROWS = 6
COMPARE_MAX_ROWS = 200


def RF_choices(RF):
    return [(value, L(ja, en)) for value, ja, en in RF.LANGUAGES]


def _figure_files(res) -> dict:
    out: dict[str, list[tuple[str, str]]] = {}
    for name, path in (getattr(res, "figures", None) or {}).items():
        items = [("png", path)]
        for suffix in ("svg", "pdf", "eps"):
            other = Path(path).with_suffix("." + suffix)
            if other.is_file():
                items.append((suffix, str(other)))
        out[name] = items
    return out


def analysis_form_values(opts: AnalysisOptions | None, form: dict | None) -> dict[str, str]:
    base = AF.fields_from_options(opts or AnalysisOptions())
    if form is None:
        return base
    return {**base, **{k: "" for k in AF.CHECKS}, **{k: v for k, v in form.items() if isinstance(v, str)}}


def compare_base_of(run_dir: str) -> str:
    try:
        p = Path(run_dir.strip()).expanduser()
        return str(p.parent) if p.is_dir() else run_dir.strip()
    except (OSError, ValueError):
        return ""


def compare_rows_from_form(fields: dict[str, str]) -> list[tuple[str, str, str]]:
    try:
        n = min(max(int(fields.get("n_rows") or 0), 0), COMPARE_MAX_ROWS)
    except ValueError:
        n = 0
    rows = [(fields.get(f"rx_name_{i}", ""), fields.get(f"rx_nu_{i}", ""), fields.get(f"rx_dir_{i}", "")) for i in range(n)]
    while rows and not any(x.strip() for x in rows[-1]):
        rows.pop()
    return rows


def compare_rows_from_json(base: str) -> list[tuple[str, str, str]]:
    from adit.analysis.compare import COMPARE_FILE, load_compare_file
    try:
        p = Path(base.strip()).expanduser() / COMPARE_FILE
        return AF.rows_from_reactions(load_compare_file(p)) if base.strip() and p.is_file() else []
    except Exception:
        return []


def _prep_context(app: WebApp, prep: dict | None) -> dict:
    from adit.codes.orca import solvent_names as orca_names
    from adit.codes.xtb import SOLVENTS, solvent_names as xtb_names
    from adit.templates import list_templates, template_dirs
    from adit.web import codefields as CF

    f = app.form
    elements = list(app.spec.elements) if app.spec is not None else []
    for k, v in f.items():
        if k.startswith(forms.MAG_PREFIXES) and (v or "").strip():
            e = k.split("__", 1)[1]
            if e not in elements:
                elements.append(e)
    hub = {}
    for p in forms.HUB_PREFIXES:
        last = max([i for i in range(1, forms.HUB_ROWS_MAX + 1) if any((f.get(f"{p}{i}_{c}") or "").strip() for c in ("el", "orb", "u", "j"))] or [0])
        hub[p] = min(forms.HUB_ROWS_MAX, max(3, last + 2))
    prov_head, prov_lines = "", []
    if app.files is not None and "spec.json" in app.files.texts:
        try:
            prov = json.loads(app.files.texts["spec.json"]).get("provenance")
        except (ValueError, AttributeError):
            prov = None
        prov_head, prov_lines = CF.provenance_summary(prov)
    try:
        templates = list_templates(app.cfg)
    except OSError:
        templates = []
    return {"prep": prep or {}, "tmarks": CF.template_marks(app.spec) if app.spec is not None else {},
            "origin_note": app.origin.describe(app.spec) if app.origin.active else "", "prov_head": prov_head, "prov_lines": prov_lines,
            "prep_elements": elements, "hub_rows": hub, "shells": CF.SHELLS,
            "xtb_solvents": {f"{m}_{g}": xtb_names(m, g) for m, g in SOLVENTS}, "orca_solvents": {m: orca_names(m) for m in ("cpcm", "smd")},
            "stage_rows": app.stage_rows(), "stage_cols": list(zip(CF.STAGE_COLUMNS, CF.stage_column_titles())),
            "stage_choices": {c: CF.stage_choices(c) for c in CF.STAGE_COLUMNS if CF.stage_choices(c) is not None},
            "stages_dir": app.stages_dir_value(), "templates": templates, "template_folders": ", ".join(str(d) for d in template_dirs(app.cfg)),
            "carry_table": CF.carry_table(), **_prep23_context(app)}


def _prep23_context(app: WebApp) -> dict:
    from adit.web import prep23 as P

    pip = P.mlip_pip_line((app.form.get("mlip_family") or "").strip())
    try:
        docs = P.doc_lines(app.spec, app.cfg) if app.spec is not None else {}
    except (OSError, ValueError):
        docs = {}
    return {"batch_kinds": P.kind_names(), "batch_notes": {k: P.kind_note(k) for k in P.BATCH_KINDS},
            "compare_kinds": P.compare_kind_names(), "box_kinds": P.box_names(), "ff_kinds": P.ff_names(),
            "neb_modes": P.neb_mode_names(), "interp_kinds": P.interp_names(), "backend_kinds": P.backend_names(),
            "qe_opts": P.qe_opt_names(), "ts_modes": P.ts_mode_names(), "rx_rows": app.batch_rx_rows(),
            "batch_dir": app.batch_dir_value(), "mlip_family_opts": [("", L("(選んでください)", "(choose one)")), ("mace_mp", "MACE-MP (mace_mp)"),
                                                                     ("mace_off", "MACE-OFF (mace_off)"), ("chgnet", "CHGNet (chgnet)")],
            "mlip_pip": pip, "doc_values": docs, "doc_max": P.DOC_MAX_LINES, "doc_more": P.doc_more}


def scan_table(rows: list[dict], reference: str) -> dict:
    from adit.scan import diff_label
    have = [r for r in rows if r["energy_ev"]]
    per_atom = bool(have) and all(r["diff_per_atom_mev"] for r in have)
    diff_key = "diff_per_atom_mev" if per_atom else "diff_from_last_mev"
    ref_ja, ref_en = diff_label(reference)
    cols = [("value", L("値", "Value")), ("energy_ev", L("エネルギー [eV]", "Energy [eV]")),
            (diff_key, L(f"{ref_ja} [meV/原子]", f"{ref_en} [meV/atom]") if per_atom else L(f"{ref_ja} [meV]", f"{ref_en} [meV]")),
            ("fmax_ev_ang", L("力の最大値 [eV/Å]", "Max force [eV/Å]")), ("pressure_gpa", L("圧力 [GPa]", "Pressure [GPa]")),
            ("note", L("注", "Note"))]
    num = [k not in ("note", "value") for k, _ in cols]
    body = [[(r.get(k) or ("" if k == "note" else "-"), n) for (k, _), n in zip(cols, num)] for r in rows]
    return {"head": [(t, n) for (_, t), n in zip(cols, num)], "body": body}


def error_headline(message: str, errors: list, prefix: str) -> str:
    if not message:
        return ""
    first = str(errors[0]) if errors else (message.strip().splitlines() or [""])[0]
    more = L(f" (ほか {len(errors) - 1} 件)", f" (+{len(errors) - 1} more)") if len(errors) > 1 else ""
    return prefix + first + more


def _unexpected(ex: Exception) -> str:
    return L(f"入力を組み立てられません (想定外のエラー。値を見直してください): {type(ex).__name__}: {ex}",
             f"cannot build the input (unexpected error; check the values): {type(ex).__name__}: {ex}")


def _error_page(ex: Exception) -> str:
    why = html.escape(f"{type(ex).__name__}: {ex}")
    title = html.escape(L("内部エラー", "Internal error"))
    body = html.escape(L("処理中に想定外のエラーが起きました。入力の値を見直して、もう一度試してください。",
                         "An unexpected error occurred. Check the values and try again."))
    back = html.escape(L("最初の画面へ戻る", "Back to the start page"))
    return (f'<!doctype html><html><head><meta charset="utf-8"><title>{title}</title></head><body>'
            f'<h1>{title}</h1><p>{body}</p><div class="error">{why}</div><p><a href="/">{back}</a></p></body></html>')


def _executable_of(run_command: str) -> str:
    tail = run_command.split("&&")[-1].strip()
    words = tail.split()
    if not words:
        return ""
    if words[0] in ("mpirun", "mpiexec", "srun"):
        rest = [w for w in words[1:] if not w.startswith("-") and not w.isdigit()]
        return rest[0] if rest else ""
    return words[0]


# ---- HTTP ----
MAX_BODY_BYTES = 256 * 1024 * 1024


class BadRequest(Exception):

    def __init__(self, status: HTTPStatus, message: str):
        super().__init__(message)
        self.status = status


def _body_length(handler: BaseHTTPRequestHandler) -> int:
    raw = (handler.headers.get("Content-Length") or "0").strip()
    try:
        length = int(raw)
    except ValueError:
        length = -1
    if length < 0:
        raise BadRequest(HTTPStatus.BAD_REQUEST, L(f"Content-Length が数値ではありません: {raw!r}", f"Content-Length is not a number: {raw!r}"))
    if length > MAX_BODY_BYTES:
        raise BadRequest(HTTPStatus(413), L(f"送られた本文が大きすぎます ({length} バイト、上限 {MAX_BODY_BYTES})",
                                            f"request body too large ({length} bytes, limit {MAX_BODY_BYTES})"))
    return length


def _parse_body(handler: BaseHTTPRequestHandler) -> tuple[dict[str, str], dict[str, tuple[str, bytes]]]:
    length = _body_length(handler)
    body = handler.rfile.read(length) if length else b""
    ctype = handler.headers.get("Content-Type", "")
    fields: dict[str, str] = {}
    files: dict[str, tuple[str, bytes]] = {}
    if ctype.startswith("multipart/form-data"):
        msg = email.message_from_bytes(b"Content-Type: " + ctype.encode() + b"\r\nMIME-Version: 1.0\r\n\r\n" + body, policy=email.policy.HTTP)
        for part in msg.iter_parts():
            name = part.get_param("name", header="content-disposition")
            if not name:
                continue
            data = part.get_payload(decode=True) or b""
            filename = part.get_filename()
            if filename:
                files[name] = (filename, data)
            else:
                fields[name] = data.decode("utf-8", errors="replace")
    else:
        for k, v in parse_qs(body.decode("utf-8", errors="replace"), keep_blank_values=True).items():
            fields[k] = v[-1]
    return fields, files


LOCAL_HOSTS = ("127.0.0.1", "localhost", "::1")
READ_ONLY_PATHS = ("/file", "/spec.json")
ANALYSIS_PATHS = ("/analysis", "/compare", "/report", "/audit_runs")
WILDCARD_HOSTS = ("", "0.0.0.0", "::")
TOKEN_CHARS = re.compile(r"[A-Za-z0-9_-]+")


def cookie_name(port: int) -> str:
    # One cookie per port, so two servers on the same PC do not overwrite each other's login.
    return f"adit_token_{port}"


def new_token() -> str:
    return secrets.token_urlsafe(24)


def redact_token(text: str) -> str:
    return re.sub(r"token=[^&\s]*", "token=***", text)


def make_handler(app: WebApp, token: str | None = None):
    class Handler(BaseHTTPRequestHandler):
        server_version = f"adit-web/{__version__}"
        _pending_cookie: str | None = None

        def log_message(self, fmt, *args):
            if os.environ.get("ADIT_WEB_LOG"):
                super().log_message(fmt, *[redact_token(str(a)) for a in args])

        def end_headers(self) -> None:
            if self._pending_cookie is not None:
                self.send_header("Set-Cookie", self._pending_cookie)
                self._pending_cookie = None
            super().end_headers()

        def _send(self, text: str, status: HTTPStatus = HTTPStatus.OK, ctype: str = "text/html; charset=utf-8") -> None:
            data = text.encode("utf-8")
            self.send_response(status); self.send_header("Content-Type", ctype); self.send_header("Content-Length", str(len(data)))
            self.end_headers(); self.wfile.write(data)

        def _send_bytes(self, data: bytes, ctype: str, filename: str | None = None) -> None:
            self.send_response(HTTPStatus.OK); self.send_header("Content-Type", ctype); self.send_header("Content-Length", str(len(data)))
            if filename:
                self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
            self.end_headers(); self.wfile.write(data)

        def _redirect(self, where: str) -> None:
            self.send_response(HTTPStatus.SEE_OTHER); self.send_header("Location", where); self.send_header("Content-Length", "0"); self.end_headers()

        def _page(self, status_text: str = "", error: str = "", message: str = "", prefix: str | None = None, message_head: str = "",
                  scan_link: str = "", scan_error: str = "", prep: dict | None = None) -> None:
            ch = app.choices()
            ja = lang.LANGUAGE != "en"
            if prefix is None:
                prefix = L("生成できません: ", "cannot generate: ")
            cfg_problem = app.config_error
            from adit.web import prep23 as P23
            from adit.web import recipe_form as R
            only_rules = [("source", "source", ["preset", "smiles", "file", "bulk", "surface", "mixture", "2d", "cluster", "polymer", "fetch"]),
                          ("tdkind", "td_kind", [k for k, _, _ in R.TWOD_KINDS]), ("clkind", "cl_kind", [k for k, _, _ in R.CLUSTER_KINDS]),
                          ("task", "task_type", list(ch["types"])),
                          ("code", "code", list(ch["codes"])), ("kp", "kp_mode", list(ch["kp_modes"])), ("disp", "dispersion", ["dftd3"]),
                          ("mixbox", "mix_box_mode", ["density", "edge"]),
                          ("batch", "batch_kind", list(P23.BATCH_KINDS)), ("cmpkind", "cmp_kind", list(P23.COMPARE_KINDS))]
            opts = {
                "only_rules": only_rules,
                "error_head": (L("環境設定ファイルを読めないため、前の設定のまま動いています: ", "cannot read the settings file; still using the previous settings: ")
                               + cfg_problem.splitlines()[0]) if cfg_problem else (error_headline(error, app.errors if error else [], prefix)
                                                                                    or error_headline(scan_error, [], prefix)
                                                                                    or error_headline(next((v for k, v in (prep or {}).items()
                                                                                                            if k.endswith("_error") and v), ""), [], "")),
                "presets_opts": [(x, preset_search_text(x)) for x in ch["presets"]],
                "bulk_struct_opts": [("", L("(元素の既定)", "(element default)"))] + [(x, x) for x in ch["crystal_structures"]],
                "surf_opts": [(x, x) for x in ch["surfaces"]],
                "type_opts": [(k, v[0] if ja else v[1]) for k, v in ch["types"].items()],
                "relax_opts": [(k, v[0] if ja else v[1]) for k, v in ch["relax_cell"].items()],
                "ens_opts": [(x, x) for x in ch["ensembles"]],
                "thermo_opts": list(ch["thermostats"].items()),
                "opt_opts": [(x, x) for x in ch["optimizers"]],
                "code_opts": [(k, forms.code_label(k)) for k in ch["codes"]],
                "sk_opts": [(x, x) for x in ch["sk_sets"]],
                "ibrion_opts": [(str(k), v) for k, v in ch["ibrion"].items()],
                "gfn_opts": list(ch["gfn"].items()),
                "opt_level_opts": [(x, x) for x in ch["opt_levels"]],
                "kp_opts": [(k, v[0] if ja else v[1]) for k, v in ch["kp_modes"].items()],
                "profile_opts": [(k, f"{k} ({p.kind})") for k, p in ch["profiles"].items()],
                "td_kind_opts": [(k, ja_t if ja else en) for k, ja_t, en in R.TWOD_KINDS],
                "cl_kind_opts": [(k, ja_t if ja else en) for k, ja_t, en in R.CLUSTER_KINDS],
                "elements_opts": [(x, x) for x in ELEMENTS],
                "rc": app.recipe_view(),
                "fetch_db_opts": [(k, name) for k, (name, _url) in FETCH_DATABASES.items()],
                "fetch_candidates": app.fetch_candidates,
                "fetch_info": app.fetch_info(),
            }
            from adit.web import codefields as CF
            f = app.form
            opts["cp2k_basis_files"], opts["cp2k_potential_files"] = CF.cp2k_files(app.cfg.cp2k_data)
            opts["cp2k_root"], opts["cp2k_cands"] = CF.cp2k_candidates(app.cfg.cp2k_data, f.get("cp_basis_file") or "BASIS_MOLOPT",
                                                                       f.get("cp_potential_file") or "GTH_POTENTIALS",
                                                                       app.spec.elements if app.spec is not None else [])
            if cfg_problem:
                error = cfg_problem + ("\n\n" + error if error else "")
            from adit.web.scan_choices import OTHER, all_choices
            file_periodic = bool(app.spec is not None and app.spec.structure.source == "file" and app.spec.structure.periodic)
            scene_json, scene_note, scene_periodic = app.scene()
            self._send(app.render("index.html", form=app.form, choices=ch, files=app.files, spec=app.spec, status_text=status_text,
                                  scene_json=scene_json, scene_note=scene_note, scene_periodic=scene_periodic,
                                  viewer_js=app.viewer_js(),
                                  error=error, message=message, written=app.written,
                                  message_head=message_head, scan_link=scan_link, scan_error=scan_error, scan_choices=all_choices(), scan_other=OTHER,
                                  scan_item=app.scan_item_value(), scan_dir=app.scan_dir_value(), file_periodic=file_periodic,
                                  no_kpoints=forms.NO_KPOINTS, **opts, **_prep_context(app, prep)))

        def _convert_page(self, message: str = "", error: str = "", fields: dict | None = None,
                          native_result=None, native_review: Path | None = None,
                          verification=None, verification_ok: bool = False) -> None:
            from adit.templates import list_templates
            templates = [t for t in list_templates(app.cfg) if not t.error]
            self._send(app.render("convert.html", templates=templates, message=message, error=error, fields=fields or {},
                                  native_result=native_result, native_review=native_review,
                                  verification=verification, verification_ok=verification_ok))

        def do_GET(self) -> None:
            self._guarded(self._do_GET)

        def do_POST(self) -> None:
            self._guarded(self._do_POST)

        def _authorized(self) -> bool:
            if not token:
                return True
            want = token.encode("utf-8")
            name = cookie_name(self.server.server_port)
            u = urlparse(self.path)
            given = parse_qs(u.query).get("token", [""])[-1]
            if given and hmac.compare_digest(given.encode("utf-8"), want):
                cookie = f"{name}={token}; Path=/; HttpOnly; SameSite=Lax"
                if self.command == "POST":
                    # A redirect would drop the request body; answer directly and set the cookie on that response.
                    self._pending_cookie = cookie
                    return True
                rest = urlencode([(k, val) for k, vals in parse_qs(u.query, keep_blank_values=True).items() for val in vals if k != "token"])
                self.send_response(HTTPStatus.SEE_OTHER)
                self.send_header("Set-Cookie", cookie)
                self.send_header("Location", u.path + ("?" + rest if rest else "")); self.send_header("Content-Length", "0"); self.end_headers()
                return False
            cookie = SimpleCookie(self.headers.get("Cookie", ""))
            if name in cookie and hmac.compare_digest(cookie[name].value.encode("utf-8"), want):
                return True
            self._send(L("この画面を開くには、adit-web を起動した端末に表示された URL (token=… 付き) を開いてください。",
                         "Open the URL (with token=...) printed in the terminal where adit-web was started."),
                       HTTPStatus.FORBIDDEN, "text/plain; charset=utf-8")
            return False

        def _guarded(self, fn) -> None:
            try:
                if not self._authorized():
                    return
                path = urlparse(self.path).path
                if path in READ_ONLY_PATHS:
                    fn()
                else:
                    with (app.analysis_lock if path in ANALYSIS_PATHS else app.lock):
                        fn()
            except (BrokenPipeError, ConnectionResetError):
                raise
            except Exception as ex:
                import traceback
                traceback.print_exc(file=sys.stderr)
                try:
                    self._send(_error_page(ex), HTTPStatus.INTERNAL_SERVER_ERROR)
                except OSError:
                    pass

        # -- GET
        def _do_GET(self) -> None:
            app.refresh_config()
            u = urlparse(self.path); q = {k: v[-1] for k, v in parse_qs(u.query).items()}
            if u.path == "/":
                self._page()
            elif u.path == "/convert":
                self._convert_page()
            elif u.path == "/spec.json":
                if app.spec is None:
                    self._send(L("先にプレビューしてください", "Preview first"), HTTPStatus.NOT_FOUND, "text/plain; charset=utf-8"); return
                self._send_bytes(app.spec.model_dump_json(indent=2).encode(), "application/json", "spec.json")
            elif u.path == "/analysis":
                d = q.get("dir") or str(app.written or "")
                md = _is_md(d) if d else False
                self._analysis_page(d, AnalysisOptions(rdf=md, msd=md) if d else None, "")
            elif u.path == "/structure.svg":
                from adit.export_image import scene_to_svg
                from adit.web.structure3d import scene_from_atoms

                if app.spec is None:
                    self._send(L("先にプレビューしてください", "Preview first"), HTTPStatus.NOT_FOUND, "text/plain")
                    return
                try:
                    svg = scene_to_svg(scene_from_atoms(app.spec.structure.atoms.to_ase()))
                except Exception as ex:
                    self._send(str(ex), HTTPStatus.BAD_REQUEST, "text/plain")
                    return
                self._send_download(svg.encode("utf-8"), "image/svg+xml", "structure.svg")
            elif u.path == "/draw":
                self._draw_page()
            elif u.path == "/report":
                self._report_page({"rep_dirs": q.get("dir", "")}, None, "")
            elif u.path == "/compare":
                base = q.get("dir") or ""
                self._compare_page(base, compare_rows_from_json(base), None, "")
            elif u.path == "/file":
                p = q.get("path", "")
                allowed = set(app.figures.values())
                for shown in list(allowed):
                    allowed |= {str(Path(shown).with_suffix("." + x)) for x in ("svg", "pdf", "eps")}
                if p not in allowed or not Path(p).is_file():
                    self._send("not found", HTTPStatus.NOT_FOUND, "text/plain"); return
                kinds = {".png": "image/png", ".svg": "image/svg+xml", ".pdf": "application/pdf",
                         ".eps": "application/postscript"}
                kind = kinds.get(Path(p).suffix.lower(), "application/octet-stream")
                if q.get("download"):
                    self._send_download(Path(p).read_bytes(), kind, Path(p).name)
                else:
                    self._send_bytes(Path(p).read_bytes(), kind)
            else:
                self._send("not found", HTTPStatus.NOT_FOUND, "text/plain")

        # -- POST
        def _do_POST(self) -> None:
            app.refresh_config()
            u = urlparse(self.path)
            try:
                fields, files = _parse_body(self)
            except BadRequest as ex:
                self._send(str(ex), ex.status, "text/plain; charset=utf-8"); return
            if u.path == "/preview":
                self._save_upload(fields, files)
                status_text, err = app.preview(fields)
                self._page(status_text=status_text, error=err)
            elif u.path == "/convert_structure":
                from adit.convert import ConversionError, convert_structure
                try:
                    cell = None
                    if (fields.get("cell") or "").strip():
                        values = [float(x.strip()) for x in fields["cell"].split(",")]
                        if len(values) == 1: values *= 3
                        if len(values) != 3: raise ConversionError(L("セルは A または A,B,C の形で入力してください", "enter the cell as A or A,B,C"))
                        cell = tuple(values)
                    p = convert_structure(fields.get("source", ""), fields.get("struct_output", ""), cell=cell,
                                          input_format=(fields.get("input_format") or "").strip() or None,
                                          output_format=(fields.get("output_format") or "").strip() or None,
                                          overwrite=forms._b(fields.get("struct_overwrite")))
                    self._convert_page(L(f"構造を書きました: {p}", f"Wrote the structure: {p}"), fields=fields)
                except (ConversionError, ValueError, OSError) as ex:
                    self._convert_page(error=str(ex), fields=fields)
            elif u.path == "/save_current_structure":
                try:
                    if app.spec is None:
                        raise ValueError(L("先に準備画面で「プレビュー」を押してください", "Press Preview on the Prepare page first"))
                    from adit.convert import write_spec_structure
                    p = write_spec_structure(app.spec, fields.get("current_output", ""),
                                             overwrite=forms._b(fields.get("current_overwrite")))
                    self._convert_page(L(f"現在の構造を書きました: {p}", f"Saved the current structure: {p}"), fields=fields)
                except (ValueError, OSError) as ex:
                    self._convert_page(error=str(ex), fields=fields)
            elif u.path == "/convert_calculation":
                try:
                    if app.spec is None:
                        raise ValueError(L("先に準備画面で「プレビュー」を押してください。現在の構造と共通条件をそこから取ります",
                                           "Press Preview on the Prepare page first; the current structure and shared settings are taken from it"))
                    from adit.convert import REPORT_FILE, retarget_spec
                    from adit.project import write_project
                    from adit.templates import load_template
                    target = load_template(fields.get("target_template", ""), app.spec, app.cfg)
                    converted, report = retarget_spec(app.spec, target)
                    written = write_project(converted, app.cfg, fields.get("calc_output", ""), overwrite=forms._b(fields.get("calc_overwrite")),
                                            extra_readme=[L("== 計算コード間の変換 ==", "== Conversion between calculation codes =="),
                                                          f"  {app.spec.method.code} -> {converted.method.code}", f"  {report['rule']}",
                                                          L(f"  詳細: {REPORT_FILE}", f"  Details: {REPORT_FILE}"), ""],
                                            extra_texts={REPORT_FILE: json.dumps(report, ensure_ascii=False, indent=2) + "\n"})
                    self._convert_page(L(f"{converted.method.code} の入力を生成しました ({len(written)} ファイル): {fields.get('calc_output','')}",
                                         f"Generated {converted.method.code} input ({len(written)} files): {fields.get('calc_output','')}"), fields=fields)
                except Exception as ex:
                    self._convert_page(error=str(ex), fields=fields)
            elif u.path == "/native_import":
                from adit.native_workflow import NativeWorkflowError, write_import_review
                try:
                    result, review = write_import_review(fields.get("native_source", ""), fields.get("native_output", ""),
                                                         code=(fields.get("native_code") or "").strip() or None)
                except (NativeWorkflowError, OSError, ValueError) as ex:
                    self._convert_page(error=str(ex), fields=fields)
                    return
                self._convert_page(fields=fields, native_result=result, native_review=review)
            elif u.path == "/native_verify":
                from adit.native_verify import verify_generated_bundle
                from adit.native_workflow import verification_passed
                project = (fields.get("verify_project") or "").strip()
                if not project:
                    self._convert_page(error=L("生成入力のディレクトリを指定してください", "Specify a generated-input directory."), fields=fields)
                    return
                result = verify_generated_bundle(project)
                self._convert_page(fields=fields, verification=result, verification_ok=verification_passed(result))
            elif u.path == "/generate":
                self._save_upload(fields, files)
                status_text, err = app.preview(fields)
                if err:
                    self._page(status_text=status_text, error=err); return
                msg, err = app.generate(overwrite=forms._b(fields.get("overwrite")))
                self._page(status_text=status_text, error=err, message=msg)
            elif u.path == "/fetch_structure":
                status_text, err, msg = app.fetch_structure(fields)
                self._page(status_text=status_text, error=err, message=msg, prefix="" if not msg else None)
            elif u.path == "/scan":
                self._save_upload(fields, files)
                status_text, err = app.preview(fields)
                if err:
                    self._page(status_text=status_text, error=err); return
                msg, err = app.scan(overwrite=forms._b(fields.get("scan_overwrite")))
                link = f"/analysis?dir={quote(str(app.scan_out))}" if msg and app.scan_out else ""
                self._page(status_text=status_text, message_head=msg, scan_link=link, scan_error=err)
            elif u.path == "/recipe":
                self._save_upload(fields, files)
                status_text, err = app.recipe_action(fields)
                self._page(status_text=status_text, error=err)
            elif u.path == "/continue":
                self._save_upload(fields, files)
                status_text, err, summary = app.continue_run(fields)
                if not summary:
                    self._page(prep={"cont_error": err, "cont_open": True}); return
                self._page(status_text=status_text, error=err, message_head=summary)
            elif u.path == "/stages":
                self._save_upload(fields, files)
                if (fields.pop("stages_action", "") or "").strip() == "add":
                    app.form = {**default_form(app.cfg.default_profile, app.form.get("output_dir", "")), **fields}
                    app.form["stg_rows"] = str(len(app.stage_rows()) + 1)
                    self._page(prep={"stages_open": True}); return
                status_text, err = app.preview(fields)
                if err:
                    self._page(status_text=status_text, error=err, prep={"stages_open": True}); return
                msg, err = app.stages(overwrite=forms._b(fields.get("stg_overwrite")))
                self._page(status_text=status_text, message_head=msg, prep={"stages_error": err, "stages_open": True})
            elif u.path == "/batch":
                self._save_upload(fields, files)
                if (fields.pop("batch_action", "") or "").strip() == "add_row":
                    app.form = {**default_form(app.cfg.default_profile, app.form.get("output_dir", "")), **fields}
                    app.form["cmp_rx_rows"] = str(len(app.batch_rx_rows()) + 1)
                    self._page(prep={"batch_open": True}); return
                status_text, err = app.preview(fields)
                if err:
                    self._page(status_text=status_text, error=err, prep={"batch_open": True}); return
                msg, err = app.batch(overwrite=forms._b(fields.get("batch_overwrite")))
                link = ""
                if msg and app.batch_out is not None:
                    res = app.batch_out
                    link = (f"/compare?dir={quote(res.compare_base)}" if res.compare_base
                            else (f"/analysis?dir={quote(res.analysis_dir)}" if res.analysis_dir else ""))
                self._page(status_text=status_text, message_head=msg, scan_link=link, prep={"batch_error": err, "batch_open": True})
            elif u.path == "/docvalue":
                self._save_upload(fields, files)
                note = app.put_doc_value(fields.get("doc_put", ""))
                status_text, err = app.preview({**fields, **{k: v for k, v in app.form.items() if k in ("ecutwfc", "ecutrho", "encut",
                                                                                                       "d3_s6", "d3_s8", "d3_a1", "d3_a2")}})
                self._page(status_text=status_text, error=err, message=note)
            elif u.path == "/template_load":
                self._save_upload(fields, files)
                status_text, err, msg = app.template_load(fields)
                if not msg:
                    self._page(prep={"tmpl_error": err, "tmpl_open": True}); return
                self._page(status_text=status_text, error=err, prep={"tmpl_msg": msg, "tmpl_open": True})
            elif u.path == "/template_save":
                self._save_upload(fields, files)
                status_text, err = app.preview(fields)
                msg, serr = app.template_save(overwrite=forms._b(fields.get("tmpl_overwrite")))
                self._page(status_text=status_text, error=err, prep={"tmpl_msg": msg, "tmpl_error": serr, "tmpl_open": True})
            elif u.path == "/origin_clear":
                app.clear_origin()
                status_text, err = app.preview(fields)
                self._page(status_text=status_text, error=err)
            elif u.path == "/load":
                name_data = files.get("spec_file")
                if not name_data:
                    self._page(error=L("spec.json が選ばれていません", "no spec.json chosen"), prefix=""); return
                try:
                    spec = CalculationSpec.from_json(name_data[1].decode("utf-8-sig"))
                except Exception as ex:
                    why = forms._pydantic_text(ex) if isinstance(ex, ValueError) else f"{type(ex).__name__}: {ex}"
                    self._page(error=L(f"spec.json を読めません: {why}", f"cannot read the spec: {why}"), prefix=""); return
                if spec.method.code not in forms.CODES:
                    self._page(error=L(
                        f"{spec.method.code} の設定はウェブ版では編集できません。spec.json と adit-gen を使ってください。",
                        f"{spec.method.code} settings cannot be edited in the web app. Use spec.json with adit-gen instead."), prefix="")
                    return
                app.form = {**forms.drop_prep_fields(app.form), **form_from_spec(spec)}
                from adit.web.recipe_form import write_steps
                if spec.structure.source != "recipe":
                    write_steps(app.form, [])
                app.adopt_loaded(spec)
                app.adopt_origin(spec)
                status_text, err = app.preview(app.form)
                self._page(status_text=status_text, error=err, message=L("spec.json を読み込みました", "loaded spec.json"))
            elif u.path == "/analysis":
                self._analysis_post(fields)
            elif u.path == "/audit_runs":
                from adit.analysis import audit_run_dirs
                ref = (fields.get("audit_ref") or "").strip()
                other = (fields.get("audit_other") or "").strip()
                kind = (fields.get("audit_kind") or "msd").strip()
                audit_fields = {"audit_ref": ref, "audit_other": other, "audit_kind": kind}
                if not ref or not other:
                    error = L("比較する 2 つの計算ディレクトリを指定してください。",
                              "Specify both calculation directories to audit.")
                    self._analysis_page(ref, None, "", audit_error=error, audit_fields=audit_fields)
                    return
                if kind not in {"msd", "energy"}:
                    error = L("点検対象には MSD かエネルギーを選んでください",
                              "Choose MSD or energy as the audit kind.")
                    self._analysis_page(ref, None, "", audit_error=error, audit_fields=audit_fields)
                    return
                try:
                    audit = audit_run_dirs([ref, other], require_md=kind == "msd")
                except (OSError, ValueError) as ex:
                    self._analysis_page(ref, None, "", audit_error=L(f"横断点検ができません: {ex}",
                                                                   f"Cannot audit these runs: {ex}"), audit_fields=audit_fields)
                    return
                self._analysis_page(ref, None, "", audit_result=audit, audit_fields=audit_fields)
            elif u.path == "/draw":
                self._draw_post(fields)
            elif u.path == "/report":
                self._report_post(fields)
            elif u.path == "/compare":
                self._compare_post(fields)
            elif u.path == "/language":
                app.set_language("en" if lang.LANGUAGE != "en" else "ja")
                self._redirect("/")
            else:
                self._send("not found", HTTPStatus.NOT_FOUND, "text/plain")

        def _save_upload(self, fields: dict[str, str], files: dict[str, tuple[str, bytes]]) -> None:
            up = files.get("structure_file")
            if up and up[1]:
                app.upload_dir.mkdir(parents=True, exist_ok=True)
                dest = app.upload_dir / upload_basename(up[0])
                dest.write_bytes(up[1])
                fields["file_path"] = str(dest); fields["source"] = "file"
                fields["file_sha_path"] = fields["file_sha256"] = ""

        def _analysis_post(self, fields: dict[str, str]) -> None:
            run_dir = fields.get("run_dir", "")
            if _is_scan_dir(run_dir):
                self._analysis_page(run_dir, AnalysisOptions(), ""); return
            if fields.get("use_stride"):
                n = fields["use_stride"]
                self._analysis_page(run_dir, None, "", form={**fields, "stride": n},
                                    message=L(f"間引きの欄に {n} を入れました。「解析を実行」をもう一度押してください。",
                                              f"The stride is now {n}. Press \"Run analysis\" again.")); return
            try:
                opts = AF.options_from_fields({**fields, "export": "on" if fields.get("action") == "export" else ""})
            except AF.FieldError as ex:
                self._analysis_page(run_dir, None, L(f"解析の設定を読めません: {ex}", f"cannot read the analysis settings: {ex}"), form=fields); return
            self._analysis_page(run_dir, opts, "", form=fields)

        def _send_download(self, data: bytes, content_type: str, filename: str) -> None:
            self.send_response(HTTPStatus.OK)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Disposition", f'attachment; filename="{filename}"')
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)

        def _draw_page(self, sketch_json: str = "{}", smiles: str = "", formula: str = "",
                       smiles_error: str = "", error: str = "", message: str = "", from_smiles: str = "") -> None:
            import json as _json

            from adit.sketch import ELEMENT_COLORS, ELEMENTS, TEMPLATE_FORMULAS, TEMPLATES
            from adit.structure import pretty_formula

            rings = [(name, pretty_formula(TEMPLATE_FORMULAS[name])) for name in TEMPLATES]
            self._send(app.render("draw.html", elements=ELEMENTS, rings=rings,
                                  colors=_json.dumps(ELEMENT_COLORS), sketch_js=app.sketch_js(),
                                  sketch_json=sketch_json, smiles=smiles, formula=formula,
                                  smiles_error=smiles_error, error=error, message=message,
                                  from_smiles=from_smiles))

        def _draw_post(self, fields: dict) -> None:
            import json as _json

            from adit.sketch import Sketch, sketch_from_dict, sketch_to_dict

            action = (fields.get("action") or "check").strip()
            raw = fields.get("sketch") or "{}"
            if action == "load":
                text = (fields.get("from_smiles") or "").strip()
                if not text:
                    self._draw_page(raw, error=L("SMILES を入れてください", "enter a SMILES"))
                    return
                try:
                    sk = Sketch.from_smiles(text)
                except Exception as ex:
                    self._draw_page(raw, error=str(ex), from_smiles=text)
                    return
                data = sketch_to_dict(sk)
                for atom in data["atoms"]:
                    atom["x"] += 360.0
                    atom["y"] += 200.0
                self._draw_page(_json.dumps(data), from_smiles=text,
                                message=L(f"SMILES から {len(data['atoms'])} 原子を読み込みました",
                                          f"loaded {len(data['atoms'])} atoms from the SMILES"))
                return
            try:
                sk = sketch_from_dict(_json.loads(raw))
            except (ValueError, TypeError) as ex:
                self._draw_page("{}", error=str(ex))
                return
            if not sk.atoms:
                self._draw_page(raw, error=L("まだ何も描かれていません", "nothing has been drawn yet"))
                return
            try:
                smiles, formula = sk.to_smiles(), sk.formula()
            except Exception as ex:
                self._draw_page(raw, smiles_error=L(f"この描き方を RDKit が受け付けません: {ex}",
                                                    f"RDKit does not accept this drawing: {ex}"))
                return
            if action in ("svg", "mol"):
                from adit.export_image import sketch_to_mol, sketch_to_svg

                try:
                    text = sketch_to_svg(sk) if action == "svg" else sketch_to_mol(sk)
                except Exception as ex:
                    self._draw_page(raw, smiles_error=L(f"書き出せません: {ex}", f"cannot export: {ex}"))
                    return
                kind = "image/svg+xml" if action == "svg" else "chemical/x-mdl-molfile"
                name = "structure.svg" if action == "svg" else "structure.mol"
                self._send_download(text.encode("utf-8"), kind, name)
                return
            if action == "use":
                app.form = {**app.form, "source": "smiles", "smiles": smiles}
                self._redirect("/")
                return
            self._draw_page(raw, smiles=smiles, formula=formula)

        def _report_page(self, fields: dict, outcome, error: str) -> None:
            from adit.gui import report_fields as RF

            base = {k: "" for k in ("rep_dirs", "rep_lang", "rep_methods", "rep_conditions", "rep_results",
                                    "rep_bundle", "rep_check")}
            self._send(app.render("report.html", f={**base, **{k: v for k, v in fields.items() if isinstance(v, str)}},
                                  outcome=outcome, error=error, languages=RF_choices(RF)))

        def _report_post(self, fields: dict) -> None:
            from adit.gui import report_fields as RF
            from adit.report import ReportError

            try:
                outcome = RF.build(RF.request_from_fields(fields))
            except (RF.ReportFieldError, ReportError, OSError) as ex:
                self._report_page(fields, None, str(ex))
                return
            self._report_page(fields, outcome, "")

        def _compare_page(self, base: str, rows: list, result, error: str) -> None:
            from adit.analysis.compare import COMPARE_FILE
            try:
                has_json = bool(base.strip()) and (Path(base.strip()).expanduser() / COMPARE_FILE).is_file()
            except (OSError, ValueError):
                has_json = False
            rows = list(rows) + [("", "", "")] * max(0, COMPARE_MIN_ROWS - len(rows))
            self._send(app.render("compare.html", base=base, rows=rows, result=result, error=error, has_json=has_json,
                                  sections=AF.compare_sections(result) if result is not None else []))

        def _compare_post(self, fields: dict[str, str]) -> None:
            from adit.analysis.compare import COMPARE_FILE, CompareError, analyze_compare, load_compare_file
            base = (fields.get("base") or "").strip()
            rows = compare_rows_from_form(fields)
            action = fields.get("action") or "compare"
            if action == "add_rows":
                self._compare_page(base, rows + [("", "", "")] * 4, None, ""); return
            if action == "load_json":
                try:
                    rows = AF.rows_from_reactions(load_compare_file(Path(base).expanduser() / COMPARE_FILE))
                except (CompareError, ValueError, OSError, KeyError) as ex:
                    self._compare_page(base, rows, None, L(f"compare.json を読めません: {ex}", f"cannot read compare.json: {ex}")); return
                self._compare_page(base, rows, None, ""); return
            if not base or not os.path.isdir(os.path.expanduser(base)):
                self._compare_page(base, rows, None, L(f"{AF.lab('compare_base')}: ディレクトリがありません ({base!r})",
                                                       f"{AF.lab('compare_base')}: directory not found ({base!r})")); return
            try:
                cres = analyze_compare(Path(base).expanduser(), AF.reactions_from_rows(rows))
            except (AF.FieldError, CompareError, ValueError, OSError) as ex:
                self._compare_page(base, rows, None, L(f"比べられません: {ex}", f"cannot compare: {ex}")); return
            app.figures = dict(cres.figures or {})
            self._compare_page(base, rows, cres, "")

        def _analysis_page(self, run_dir: str, opts: AnalysisOptions | None, error: str, form: dict | None = None, message: str = "",
                           audit_result=None, audit_error: str = "", audit_fields: dict | None = None) -> None:
            from adit.analysis.trajectory import TrajectoryTooLarge
            res, notice, table, too_large, suggested = None, "", None, "", None
            scan = _is_scan_dir(run_dir)
            if scan and opts is not None:
                from adit.scan import analyze_scan
                try:
                    res = analyze_scan(Path(run_dir).expanduser())
                    app.figures = dict(res.figures or {})
                    table = scan_table(res.rows, res.reference)
                except Exception as ex:
                    error = L(f"解析できません: {ex}", f"cannot analyze: {ex}")
                    res = None
                opts = None
            elif opts is not None and run_dir:
                d = Path(run_dir).expanduser()
                if os.path.isdir(d) and not any(os.path.exists(d / f) for f in RESULT_FILES):
                    notice = L("計算結果が見つかりません (output.log などの出力ファイルがありません)。まだ実行していないようです。"
                               "下の要約は、出力が無いまま入力ファイルだけから求めたものです。",
                               "No calculation results found (no output files such as output.log). It seems the calculation has not been run yet; "
                               "the summary below is computed from the input files only.")
                try:
                    res = app.analyze(run_dir, opts)
                except TrajectoryTooLarge as ex:
                    error = L(f"解析できません: {ex}", f"cannot analyze: {ex}")
                    too_large, suggested = str(ex), ex.suggested_stride
                except Exception as ex:
                    error = L(f"解析できません: {ex}", f"cannot analyze: {ex}")
            f = analysis_form_values(opts, form)
            plain = res is not None and not scan
            export_dir, export_readme = AF.read_export_readme(res) if plain else ("", "")
            self._send(app.render("analysis.html", run_dir=run_dir, result=res, error=error, notice=notice, message=message, f=f,
                                  off=" disabled" if scan else "", more_open=bool(message or too_large) or AF.details_changed(f),
                                  sections=AF.result_sections(res) if plain else [], export_dir=export_dir, export_readme=export_readme,
                                  figure_files=_figure_files(res),
                                  too_large=too_large, suggested_stride=suggested, th_models=AF.choices(AF.THERMO_MODELS),
                                  geometries=AF.choices(AF.GEOMETRIES), imaginary=AF.choices(AF.IMAGINARY), uv_shapes=AF.choices(AF.UV_SHAPES),
                                  msd_axes=AF.choices(AF.MSD_AXES), vanhove_displacements=AF.choices(AF.VANHOVE_DISPLACEMENTS),
                                  axes=AF.choices(AF.AXES), plane_axes=AF.choices(AF.PLANE_AXES),
                                  fes_units=AF.choices(AF.FES_UNITS), cube_units=AF.choices(AF.CUBE_UNITS),
                                  tick_directions=AF.choices(AF.TICK_DIRECTIONS),
                                  grid_choices=AF.choices(AF.GRID_CHOICES),
                                  spine_choices=AF.choices(AF.SPINE_CHOICES),
                                  ph=AF.ph, compare_link="/compare?dir=" + quote(compare_base_of(run_dir)) if run_dir.strip() else "/compare",
                                  scan=scan, scan_table=table, scan_dir=str(Path(run_dir).expanduser()) if scan else "",
                                  audit_result=audit_result, audit_error=audit_error,
                                  audit_fields=audit_fields or {"audit_ref": run_dir, "audit_other": "", "audit_kind": "msd"}))

    return Handler


class _ThreadingHTTPServer6(ThreadingHTTPServer):
    address_family = socket.AF_INET6

    def server_bind(self) -> None:
        if self.server_address[0] in ("", "::"):
            try:
                self.socket.setsockopt(socket.IPPROTO_IPV6, socket.IPV6_V6ONLY, 0)
            except (OSError, AttributeError):
                pass
        super().server_bind()


def serve(app: WebApp, host: str = "127.0.0.1", port: int = 8765, open_browser: bool = False,
          token: str | None = None) -> ThreadingHTTPServer:
    cls = _ThreadingHTTPServer6 if ":" in host else ThreadingHTTPServer
    httpd = cls((host, port), make_handler(app, token))
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(server_url(browse_host(host), httpd.server_port, token))).start()
    return httpd


def browse_host(host: str) -> str:
    return "127.0.0.1" if host in WILDCARD_HOSTS else host


def server_url(host: str, port: int, token: str | None = None) -> str:
    shown = f"[{host}]" if ":" in host else host
    return f"http://{shown}:{port}/" + (f"?token={quote(token, safe='')}" if token else "")


def server_lines(host: str, port: int, token: str | None = None) -> list[str]:
    stop = "(Ctrl-C で停止 / Ctrl-C to stop)"
    if host not in WILDCARD_HOSTS:
        return [f"adit-web {__version__}: {server_url(host, port, token)}  {stop}"]
    return [f"adit-web {__version__}:",
            L("  この PC で開く:   ", "  on this PC:       ") + server_url("127.0.0.1", port, token),
            L("  他の PC から開く: ", "  from other PCs:   ") + server_url(socket.getfqdn(), port, token),
            L("  (繋がらなければ、ホスト名をこの PC の IP アドレスに置き換えてください)",
              "  (if it does not connect, replace the host name with this PC's IP address)"),
            f"  {stop}"]


def main(argv: list[str] | None = None) -> int:
    ensure_printable_stdio()
    ap = argparse.ArgumentParser(prog="adit-web", description=L("ブラウザで使う ADIT (準備 → 実行 → 解析)", "ADIT in the browser (prepare, run, analyze)"))
    ap.add_argument("--host", default="127.0.0.1", help=L("待ち受けるアドレス (既定 127.0.0.1)。これ以外では合言葉 (トークン) が必要になります",
                                                          "bind address (default 127.0.0.1); other addresses require a token"))
    ap.add_argument("--token", help=L("合言葉を自分で決める (省略すると、127.0.0.1 以外では毎回作ります)",
                                      "use this token (default: a new one each start, except on 127.0.0.1)"))
    ap.add_argument("--no-token", action="store_true", help=L("合言葉を使わない (信頼できるネットワークの中だけで)",
                                                              "do not require a token (trusted networks only)"))
    ap.add_argument("--port", type=int, default=8765)
    ap.add_argument("--config", help=L("cluster.toml のパス", "path to cluster.toml"))
    ap.add_argument("--open", action="store_true", help=L("起動後にブラウザを開きます", "open the browser after start"))
    args = ap.parse_args(argv)
    if args.token is not None and not TOKEN_CHARS.fullmatch(args.token):
        print(L("--token に使えるのは英数字と _ - だけです", "--token may contain only letters, digits, _ and -"), file=sys.stderr)
        return 2
    try:
        cfg, path, created = ensure_config(Path(args.config).expanduser() if args.config else config_path())
    except ConfigError as ex:
        print(str(ex), file=sys.stderr)
        return 2
    lang.set_language(env_var("LANG", cfg.language))
    if created:
        print(first_run_message(path) + "\n" + L("  書き換えたら、ブラウザの画面を読み込み直すだけで反映されます (adit-web を起動し直す必要はありません)。",
                                                "  After editing, just reload the page in the browser (no need to restart adit-web)."), file=sys.stderr)
    token = None if args.no_token else (args.token or (None if args.host in LOCAL_HOSTS else new_token()))
    if args.no_token and args.host not in LOCAL_HOSTS:
        print(L("警告: 合言葉なしのため、このアドレスに届く人は誰でも入力の生成と実行ができます", "warning: no token; anyone reaching this host can generate and run"), file=sys.stderr)
    app = WebApp(cfg, path)
    try:
        httpd = serve(app, args.host, args.port, args.open, token)
    except OSError as ex:
        print(L(f"そのアドレスでは待ち受けできません: {args.host}:{args.port} ({ex})", f"cannot listen on {args.host}:{args.port} ({ex})"), file=sys.stderr)
        return 2
    print("\n".join(server_lines(args.host, httpd.server_port, token)), file=sys.stderr)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        httpd.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
