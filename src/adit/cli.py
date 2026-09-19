
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from adit import lang
from adit.lang import L
from adit.config import ConfigError, config_path, ensure_config, env_var, first_run_message, unknown_keys_message
from adit.project import OutputNotEmpty, ProjectError, build_project, write_project
from adit.spec import CalculationSpec
from adit.validate import validate

_OVERWRITE_HINT = ("  上書きしてよければ、同じコマンドの最後に --overwrite を付けて実行してください。",
                   "  To overwrite, run the same command again with --overwrite at the end.")


def _apply_sets(spec: CalculationSpec, sets: list[str]) -> CalculationSpec:
    from adit.scan import apply_value

    for item in sets:
        if "=" not in item:
            raise ValueError(L(f"--set は「項目=値」の形で書いてください: {item!r}", f"--set must be item=value: {item!r}"))
        path, _, value = item.partition("=")
        comment = spec.meta.comment
        spec = apply_value(spec, path.strip(), value.strip())
        spec = spec.model_copy(update={"meta": spec.meta.model_copy(update={"comment": (comment + f" set {path.strip()}={value.strip()}").strip()})})
    return spec


def _run_group(args, spec: CalculationSpec, cfg, output_dir) -> int | None:
    if not (args.compare_set or args.conformers is not None or args.neb or args.phonons or args.elastic or args.ts or args.sella):
        return None
    try:
        if args.compare_set:
            from adit.compare_sets import load_set, write_compare_set
            data, where = load_set(args.compare_set)
            dirs = write_compare_set(spec, cfg, output_dir, data, where, overwrite=args.overwrite)
        elif args.conformers is not None:
            from adit.conformers import write_conformers
            dirs = write_conformers(spec, cfg, output_dir, n=args.conformers, rmsd=args.rmsd, seed=args.conf_seed, force_field=args.conf_ff,
                                    heavy_only=not args.conf_all_atoms, max_iters=args.conf_max_iters, keep=args.conf_keep, smiles=args.smiles,
                                    overwrite=args.overwrite)
        elif args.neb:
            from adit.neb_setup import NebOptions, write_neb
            opts = NebOptions(images=args.images, mode=args.neb_mode, interpolation=args.neb_interp, climb=args.climb,
                              vasp_spring=args.vasp_spring, cp2k_k_spring=args.cp2k_k_spring, qe_opt_scheme=args.qe_opt_scheme)
            dirs = write_neb(spec, cfg, output_dir, args.neb[0], args.neb[1], opts, where=Path.cwd(), overwrite=args.overwrite)
        elif args.phonons:
            from adit.phonon_setup import parse_dim, write_phonons
            dirs = write_phonons(spec, cfg, output_dir, parse_dim(args.phonons), distance=args.displacement, backend=args.phonon_backend,
                                 dos_mesh=parse_dim(args.phonon_dos_mesh) if args.phonon_dos_mesh else None,
                                 dos_width_thz=args.phonon_dos_width, overwrite=args.overwrite)
        elif args.elastic:
            from adit.elastic_setup import ElasticError, parse_values, write_elastic
            try:
                comps = [int(x) for x in args.elastic_components.split(",") if x.strip()]
            except ValueError as ex:
                raise ElasticError(L(f"成分は 1〜6 の整数をカンマで区切って書いてください: {args.elastic_components!r}",
                                     f"write components as integers 1..6 separated by commas: {args.elastic_components!r}")) from ex
            dirs = write_elastic(spec, cfg, output_dir, parse_values(args.elastic), comps, overwrite=args.overwrite)
        elif args.ts:
            from adit.ts_setup import write_ts
            dirs = write_ts(spec, cfg, output_dir, args.ts, overwrite=args.overwrite)
        else:
            from adit.ts_setup import write_sella
            dirs = write_sella(spec, cfg, output_dir, irc=not args.sella_no_irc, overwrite=args.overwrite)
    except OutputNotEmpty as ex:
        print(f"{ex}\n" + L(*_OVERWRITE_HINT), file=sys.stderr)
        return 1
    except ImportError as ex:
        print(L(f"必要なパッケージがありません: {ex}", f"a required package is missing: {ex}"), file=sys.stderr)
        return 1
    except (ProjectError, ConfigError, ValueError) as ex:
        print(str(ex), file=sys.stderr)
        return 1
    out = Path(output_dir).resolve()
    names = ", ".join(str(d.relative_to(out)) if out in d.resolve().parents else d.name for d in dirs)
    print(L(f"{len(dirs)} 個のディレクトリを {out} に書きました: {names}\n手順は {out / 'README.txt'}",
            f"wrote {len(dirs)} directories in {out}: {names}\nsee {out / 'README.txt'}"))
    return 0


def main(argv: list[str] | None = None) -> int:
    from adit.compat import ensure_printable_stdio

    ensure_printable_stdio()
    ap = argparse.ArgumentParser(prog="adit-gen", description=L("spec.json から計算ディレクトリ一式を生成します (GUI 不要)", "generate a calculation directory from spec.json (no GUI needed)"))
    ap.add_argument("spec", nargs="?", help=L("spec.json のパス (--continue-from で条件も前のままなら、ここに出力ディレクトリ)",
                                              "path to spec.json (with --continue-from and unchanged conditions, the output directory)"))
    ap.add_argument("output_dir", nargs="?", help=L("出力ディレクトリ (--validate のときは不要)", "output directory (not needed with --validate)"))
    ap.add_argument("--config", help=L("cluster.toml のパス (既定はユーザーの環境設定ファイル)", "path to cluster.toml (default: the user settings file)"))
    ap.add_argument("--run", action="store_true", help=L(
        "生成したあと、そのままこの PC で実行します (bash submit.sh)。画面と同じ判定で、"
        "実行できないときは理由を出して止まります。**クラスタへの投入はしません**",
        "after generating, run it here with bash submit.sh; the same checks as the screens apply and the reason is "
        "printed when it cannot run. ADIT never submits to a cluster"))
    ap.add_argument("--overwrite", action="store_true", help=L("出力ディレクトリが空でなくても上書きします", "overwrite even if the output directory is not empty"))
    ap.add_argument("--validate", action="store_true", help=L("検証だけ行い、生成しません", "validate only, do not generate"))
    ap.add_argument("--structures", metavar="FILE", action="append", default=[], help=L(
        "同じ条件で構造だけを差し替えた計算をまとめて作ります (ファイル名か * を使った書き方。複数回指定可)。"
        "1 つのファイルに複数フレームがあれば、フレームごとに 1 つの計算にします。例 --structures \"mols/*.xyz\"",
        "generate one calculation per structure, with all other settings unchanged (a file name or a pattern; repeatable). "
        "A file with several frames gives one calculation per frame. e.g. --structures \"mols/*.xyz\""))
    ap.add_argument("--core", metavar="SMILES", help=L(
        "骨格の SMILES ([*:1] のような印を置く)。--substituent と一緒に使い、置換基の組み合わせごとに計算を作ります",
        "SMILES of the core with attachment points such as [*:1]; use with --substituent to enumerate combinations"))
    ap.add_argument("--substituent", metavar="N=SMILES,SMILES,...", action="append", default=[], help=L(
        "印ごとの置換基 (例 --substituent \"1=[H],C,OC\")。複数回指定できます", "substituents per attachment point (repeatable), e.g. --substituent \"1=[H],C,OC\""))
    ap.add_argument("--scan", metavar="PATH=V1,V2,...", action="append", default=[], help=L(
        "1 つの値だけを振った入力を、値ごとのディレクトリにまとめて作ります (収束の確認・格子定数の探索)。"
        "例 method.ecutwfc=30,40,50 / kpoints.mesh=4x4x4,6x6x6 / scale=0.98,1.00,1.02 (周期系のセルを伸縮)。"
        "実行したあと adit-analyze <出力ディレクトリ> --scan で値とエネルギーの表を作ります。"
        "**2 回以上書くと、値のすべての組み合わせ (格子状) を作ります** (例 --scan method.ecutwfc=30,40 --scan kpoints.mesh=4x4x4,6x6x6)。"
        "構造の幾何も振れます: geom.distance(0,1)=0.9,1.0 / geom.angle(0,1,2)=100,105 / geom.dihedral(0,1,2,3)=0,30,60 (原子の番号は 0 始まり)",
        "generate one directory per value of a single parameter (convergence checks, lattice scans). "
        "e.g. method.ecutwfc=30,40,50 / kpoints.mesh=4x4x4,6x6x6 / scale=0.98,1.00,1.02 (scales a periodic cell). "
        "After running them, adit-analyze <output> --scan makes a value-energy table. "
        "**Give it more than once to scan a grid of all value combinations** (e.g. --scan method.ecutwfc=30,40 --scan kpoints.mesh=4x4x4,6x6x6). "
        "Geometry can be scanned too: geom.distance(0,1)=0.9,1.0 / geom.angle(0,1,2)=100,105 / geom.dihedral(0,1,2,3)=0,30,60 (0-based atom indices)"))
    ap.add_argument("--stages", metavar="STAGES.json", help=L(
        "段階に分けた計算 (例 最小化 → NVT → NPT → 本計算) を stage_01_… のディレクトリに並べて作ります。2 段階目からは前の段階の最終構造 (と速度) から始まります",
        "generate a staged calculation (e.g. minimization -> NVT -> NPT -> production) as stage_01_... directories; each stage starts from the previous one's final structure (and velocities)"))
    ap.add_argument("--continue-from", metavar="DIR", help=L("前の計算 (ADIT が生成して実行したディレクトリ) の最終構造 (MD なら速度も) から続きを作ります",
                                                            "continue from the final structure (and MD velocities) of a previous run generated by ADIT"))
    ap.add_argument("--no-velocities", action="store_true", help=L("--continue-from で、MD の続きでも速度を引き継がない (新しく初速を作る)",
                                                                 "with --continue-from, do not carry over MD velocities (generate new ones)"))
    ap.add_argument("--template", metavar="NAME", help=L("研究室の雛形 (構造抜きの spec.json) の条件を使い、構造だけを spec.json (か --continue-from の前の計算) から取ります",
                                                         "use the conditions of a group template (a spec.json without structure); the structure comes from spec.json (or --continue-from)"))
    ap.add_argument("--save-template", metavar="NAME", help=L("spec.json から構造を除いたものを雛形として保存します (置き場所は環境設定の templates_dir、無ければ ~/.config/adit/templates/)",
                                                              "save spec.json without its structure as a template (in templates_dir of the settings, or ~/.config/adit/templates/)"))
    ap.add_argument("--list-templates", action="store_true", help=L("雛形の一覧を表示します", "list the templates"))
    ap.add_argument("--list-samples", action="store_true", help=L(
        "同梱のサンプル (examples/) の一覧を表示します", "list the bundled samples (examples/)"))
    ap.add_argument("--sample", metavar="NAME", help=L(
        "サンプルの spec.json を写して、書き換えの出発点にします (写し先は出力ディレクトリの引数。--list-samples で名前が分かります)",
        "copy the spec.json of a sample as a starting point (the destination is the output argument; see --list-samples)"))
    ap.add_argument("--check-config", action="store_true", help=L(
        "環境設定ファイルの中身を確かめて表示します (プロファイル・置き場所・知らない項目)。生成はしません",
        "check and print the settings file (profiles, folders, unknown entries); nothing is generated"))
    ap.add_argument("--set", dest="sets", action="append", default=[], metavar="PATH=VALUE",
                    help=L("読み込んだ設定の 1 か所を変えます (繰り返し可。例 --set task.md.ensemble=NPT)", "change one item of the loaded settings (repeatable; e.g. --set task.md.ensemble=NPT)"))
    ap.add_argument("--print", dest="show", metavar="FILE", help=L("生成せずに、指定したファイル (例 dftb_in.hsd, INCAR, submit.sh) の中身を表示します", "do not generate; print the given file (e.g. dftb_in.hsd, INCAR, submit.sh)"))
    g = ap.add_argument_group(L("組・配座・経路・フォノン・弾性・遷移状態 (1 つだけ)", "sets, conformers, paths, phonons, elastic, transition states (one at a time)"))
    g.add_argument("--documented", action="store_true", help=L(
        "生成せずに、使うパラメータのファイル (UPF・SSSP の JSON・POTCAR・CP2K の基底・SK セットの README) に書かれた値を、ファイルと行つきで表示します",
        "do not generate; print values written in the parameter files used (UPF, SSSP JSON, POTCAR, CP2K basis, SK README) with file and line"))
    g.add_argument("--compare-set", metavar="SET.json", help=L("比べる計算の組 (吸着・反応・溶媒和) を同じ条件で生成し、compare.json を書きます (compare_sets.py に書式)",
                                                              "generate a set of runs to compare (adsorption / reaction / solvation) with shared settings and write compare.json"))
    g.add_argument("--conformers", type=int, metavar="N", help=L("RDKit の ETKDG で N 個の配座を作り、力場で最適化して重複を除き、配座ごとの計算にします (--rmsd が要ります)",
                                                                "embed N conformers with RDKit ETKDG, optimize with a force field, drop duplicates, one run per conformer (needs --rmsd)"))
    g.add_argument("--rmsd", type=float, metavar="Å", help=L("配座を重複とみなす RMSD [Å] (既定はありません)", "RMSD in Å below which conformers count as duplicates (no default)"))
    g.add_argument("--conf-seed", type=int, default=12345, metavar="S", help=L("ETKDG の乱数の種 (既定 12345。RDKit 2026.03 では 0 だと全部の配座が同じ座標になった)",
                                                                                "ETKDG random seed (default 12345; with 0, RDKit 2026.03 gave identical coordinates for all conformers)"))
    g.add_argument("--conf-ff", default="MMFF94", choices=("MMFF94", "MMFF94s", "UFF"), help=L("配座の最適化の力場 (既定 MMFF94)", "force field for the conformers (default MMFF94)"))
    g.add_argument("--conf-max-iters", type=int, default=200, metavar="N", help=L("力場の最適化の反復の上限 (既定 200 = RDKit の既定)", "max force-field iterations (default 200 = RDKit default)"))
    g.add_argument("--conf-all-atoms", action="store_true", help=L("RMSD を水素も含めた全原子で測ります (既定は水素を除く)", "measure the RMSD over all atoms (default: without hydrogens)"))
    g.add_argument("--conf-keep", type=int, metavar="N", help=L("残す配座の数の上限 (エネルギーの低い順)", "keep at most N conformers (lowest energy first)"))
    g.add_argument("--smiles", metavar="SMILES", help=L("配座を作る分子の SMILES (spec.json の構造が SMILES 由来でないとき)", "SMILES of the molecule for conformers (if the structure in spec.json is not from SMILES)"))
    g.add_argument("--neb", nargs=2, metavar=("START", "END"), help=L("始状態と終状態の構造のファイルから NEB の像を作ります (START に spec と書くと spec.json の構造)。--images が要ります",
                                                                     "make NEB images from the initial and final structure files (START = spec uses spec.json); needs --images"))
    g.add_argument("--images", type=int, metavar="N", help=L("NEB の中間の像の数", "number of intermediate NEB images"))
    g.add_argument("--climb", action="store_true", help=L("climbing image (QE: CI_scheme auto、CP2K: CI-NEB)。VASP 本体には無いので止めます", "climbing image (QE: CI_scheme auto, CP2K: CI-NEB); not in plain VASP, so it stops"))
    g.add_argument("--neb-mode", default="native", choices=("native", "images"), help=L("native: 計算コードの NEB (VASP / QE / CP2K)。images: 像ごとの一点計算 (どのコードでも)",
                                                                                        "native: the code's NEB (VASP / QE / CP2K); images: one single point per image (any code)"))
    g.add_argument("--neb-interp", default="idpp", choices=("idpp", "linear"), help=L("中間の像の作り方 (ASE の NEB.interpolate。既定 idpp)", "how to make the images (ASE NEB.interpolate; default idpp)"))
    g.add_argument("--vasp-spring", type=float, metavar="X", help=L("VASP の SPRING (書かなければ VASP の既定)", "VASP SPRING (VASP default if omitted)"))
    g.add_argument("--cp2k-k-spring", type=float, metavar="X", help=L("CP2K の K_SPRING (書かなければ CP2K の既定)", "CP2K K_SPRING (CP2K default if omitted)"))
    g.add_argument("--qe-opt-scheme", metavar="S", help=L("neb.x の opt_scheme (書かなければ neb.x の既定)", "neb.x opt_scheme (neb.x default if omitted)"))
    g.add_argument("--phonons", metavar="NxNxN", help=L("超格子 (例 2x2x2) と有限変位でフォノンの計算を並べます (DFTB+ / VASP / QE)", "phonons from supercell finite displacements, e.g. 2x2x2 (DFTB+ / VASP / QE)"))
    g.add_argument("--displacement", type=float, metavar="Å", help=L("変位の大きさ [Å] (省くと phonopy / ASE の既定)", "displacement in Å (phonopy / ASE default if omitted)"))
    g.add_argument("--phonon-backend", default="auto", choices=("auto", "phonopy", "ase"), help=L("変位の作り方 (auto: phonopy があれば phonopy)", "how to make displacements (auto: phonopy if installed)"))
    g.add_argument("--phonon-dos-mesh", metavar="NxNxN", help=L("状態密度の q 点のメッシュ (与えたときだけ total_dos.dat を書きます。既定はありません)", "q-point mesh for the phonon DOS (total_dos.dat only if given; no default)"))
    g.add_argument("--phonon-dos-width", type=float, metavar="THz", help=L("ASE で集めるときの DOS の幅 [THz] (ASE の DOS には必須)", "DOS width in THz for the ASE backend (required for the ASE DOS)"))
    g.add_argument("--elastic", metavar="S1,S2,...", help=L("弾性定数のための歪みの大きさ (例 --elastic=-0.01,0.01)。成分ごとに計算を並べます", "strain magnitudes for elastic constants (e.g. --elastic=-0.01,0.01); one run per component and strain"))
    g.add_argument("--elastic-components", default="1,2,3,4,5,6", metavar="1,2,...", help=L("歪みを与える Voigt の成分 (既定は 6 つ全部)", "Voigt components to strain (default all six)"))
    g.add_argument("--ts", choices=("ts", "irc", "ts+irc"), help=L("ORCA の遷移状態の探索 (OptTS)、IRC、または段階に分けた両方", "ORCA transition-state search (OptTS), IRC, or both as stages"))
    g.add_argument("--sella", action="store_true", help=L("ASE + Sella の遷移状態の探索と IRC のスクリプトを書きます (xtb は tblite、機械学習ポテンシャル)", "write an ASE + Sella TS search and IRC script (xtb via tblite, or an ML potential)"))
    g.add_argument("--sella-no-irc", action="store_true", help=L("Sella のスクリプトで IRC をしない", "no IRC in the Sella script"))
    args = ap.parse_args(argv)

    try:
        cfg, cfg_file, created = ensure_config(Path(args.config).expanduser() if args.config else config_path())
        warning = unknown_keys_message(cfg)
        if warning:
            print(warning, file=sys.stderr)
    except ConfigError as ex:
        print(str(ex), file=sys.stderr)
        return 2
    except OSError as ex:
        print(L(f"環境設定ファイルを読み書きできません: {ex}", f"cannot read or write the settings file: {ex}"), file=sys.stderr)
        return 2
    lang.set_language(env_var("LANG", cfg.language))
    if created:
        print(first_run_message(cfg_file), file=sys.stderr)

    from adit.templates import TemplateError, list_templates, load_template, save_template, template_dirs

    if args.list_samples:
        from adit.samples import list_samples, missing_message

        samples = list_samples()
        if not samples:
            print(missing_message())
            return 1
        print(L(f"サンプル {len(samples)} 件 (spec.json を写して書き換えると、自分の計算の出発点になります)",
                f"{len(samples)} samples (copy a spec.json and edit it as a starting point)"))
        print(L(f"  {'名前':<32} {'コード':<8} {'計算の種類':<20} {'組成':<12} 説明",
                f"  {'name':<34} {'code':<10} {'task':<22} {'formula':<12} note"))
        for sample in samples:
            print(sample.line())
        print(L("  写すには: adit-gen --sample <名前> <書き出す spec.json>",
                "  to copy one: adit-gen --sample <name> <spec.json to write>"))
        return 0
    if args.sample:
        from adit.samples import copy_sample

        destination = args.output_dir or args.spec
        if destination is None:
            ap.error(L("--sample には写し先のファイル名も指定してください", "--sample also needs a destination file name"))
        try:
            written = copy_sample(args.sample, destination, overwrite=args.overwrite)
        except (FileNotFoundError, FileExistsError) as ex:
            print(str(ex), file=sys.stderr)
            return 1
        print(L(f"サンプルの spec.json を写しました: {written}\n書き換えてから adit-gen {written} <出力ディレクトリ> で生成します",
                f"copied the sample spec.json: {written}\nedit it, then run adit-gen {written} <output>"))
        return 0
    if args.check_config:
        old_note = ""
        from adit.config import OLD_APP_NAMES

        text = str(cfg.source_path or "").replace("\\", "/")
        if cfg.source_path and any(f"/{name}/" in text for name in OLD_APP_NAMES):
            old_note = L(" (改名前の置き場所です。そのまま読めます。~/.config/adit/cluster.toml へ写しても構いません)",
                         " (this is the pre-rename location; it is still read. You may copy it to ~/.config/adit/cluster.toml)")
        print(L(f"環境設定ファイル: {cfg.source_path}{old_note}", f"settings file: {cfg.source_path}{old_note}"))
        for name, value in (("sk_root", cfg.sk_root), ("pseudo_root", cfg.pseudo_root),
                            ("cp2k_data", cfg.cp2k_data), ("templates_dir", cfg.templates_dir)):
            if not value:
                print(L(f"  {name}: (指定なし)", f"  {name}: (not set)"))
                continue
            exists = Path(value).expanduser().is_dir()
            print(L(f"  {name}: {value} {'(あります)' if exists else '(ありません)'}",
                    f"  {name}: {value} {'(exists)' if exists else '(missing)'}"))
        print(L(f"  既定のプロファイル: {cfg.default_profile}", f"  default profile: {cfg.default_profile}"))
        for name, profile in sorted(cfg.profiles.items()):
            codes = sorted(set(profile.commands) | set(profile.code_modules))
            print(L(f"  プロファイル {name}: 種類 {profile.kind}"
                    + (f"、コードごとの設定 {', '.join(codes)}" if codes else "、コードごとの設定なし")
                    + (f"、投入 {profile.submit}" if profile.kind != "direct" else ""),
                    f"  profile {name}: kind {profile.kind}"
                    + (f", per-code settings for {', '.join(codes)}" if codes else ", no per-code settings")
                    + (f", submit with {profile.submit}" if profile.kind != "direct" else "")))
        if cfg.default_profile not in cfg.profiles:
            print(L(f"  注意: 既定のプロファイル {cfg.default_profile!r} がありません",
                    f"  note: the default profile {cfg.default_profile!r} does not exist"), file=sys.stderr)
        message = unknown_keys_message(cfg)
        if message:
            print(message, file=sys.stderr)
        return 0
    if args.list_templates:
        infos = list_templates(cfg)
        print(L(f"雛形の置き場所: {', '.join(str(d) for d in template_dirs(cfg))}", f"template folders: {', '.join(str(d) for d in template_dirs(cfg))}"))
        for i in infos:
            print(f"  {i.name:<24} {i.code:<10} {i.error or i.comment}")
        if not infos:
            print(L("  (雛形はありません)", "  (no templates)"))
        return 0

    output_dir = args.output_dir
    try:
        if args.continue_from:
            from adit.continuation import ContinuationError, continue_from
            if args.spec is None:
                ap.error(L("--continue-from には出力ディレクトリが要ります", "--continue-from needs an output directory"))
            conditions = None
            if args.output_dir is None:
                output_dir = args.spec
            else:
                conditions = CalculationSpec.load(args.spec)
            if args.template:
                prev_st = CalculationSpec.load(Path(args.continue_from).expanduser() / "spec.json")
                conditions = load_template(args.template, prev_st, cfg)
            try:
                spec = continue_from(args.continue_from, conditions, velocities=not args.no_velocities)
            except ContinuationError as ex:
                print(str(ex), file=sys.stderr)
                return 1
        else:
            if args.spec is None:
                ap.error(L("spec.json を指定してください", "give a spec.json"))
            spec = CalculationSpec.load(args.spec)
            if args.template:
                spec = load_template(args.template, spec, cfg)
        spec = _apply_sets(spec, args.sets)
    except TemplateError as ex:
        print(str(ex), file=sys.stderr)
        return 1
    except Exception as ex:
        print(CalculationSpec.describe_error(ex), file=sys.stderr)
        return 2

    if args.save_template:
        try:
            p = save_template(spec, args.save_template, cfg=cfg, comment=spec.meta.comment, overwrite=args.overwrite)
        except TemplateError as ex:
            print(str(ex), file=sys.stderr)
            return 1
        print(L(f"雛形を保存しました: {p}", f"saved the template: {p}"))
        return 0

    if args.documented:
        from adit.docvalues import documented_values, format_values
        print(format_values(documented_values(spec, cfg)))
        return 0

    if args.validate or args.show:
        errs = validate(spec, cfg, output_dir=output_dir)
        if errs:
            print(L("エラー:", "Errors:"), file=sys.stderr)
            for e in errs:
                print(f"  - {e}", file=sys.stderr)
            return 1
        if args.show:
            try:
                files = build_project(spec, cfg, output_dir=output_dir)
            except ProjectError as ex:
                print(str(ex), file=sys.stderr)
                return 1
            if args.show not in files.texts:
                print(L(f"{args.show} は生成されません。生成されるのは: {', '.join(files.texts)}", f"{args.show} is not generated. Available: {', '.join(files.texts)}"), file=sys.stderr)
                return 1
            sys.stdout.write(files.texts[args.show])
            return 0
        print(L("エラーはありません (生成できます)", "no errors (can be generated)"))
        return 0

    if not output_dir:
        ap.error(L("出力ディレクトリを指定してください (検証だけなら --validate)", "an output directory is required (use --validate to only validate)"))
    if args.stages and args.scan:
        ap.error(L("--stages と --scan は同時に使えません", "--stages and --scan cannot be combined"))
    chosen = [n for n, v in (("--stages", args.stages), ("--scan", args.scan), ("--structures", bool(args.structures)),
                             ("--core", bool(args.core)), ("--compare-set", args.compare_set), ("--conformers", args.conformers is not None),
                             ("--neb", args.neb), ("--phonons", args.phonons), ("--elastic", args.elastic), ("--ts", args.ts), ("--sella", args.sella)) if v]
    if len(chosen) > 1:
        ap.error(L(f"{' と '.join(chosen)} は同時に使えません", f"{' and '.join(chosen)} cannot be combined"))
    if args.images is not None and not args.neb:
        ap.error(L("--images は --neb と一緒に指定してください", "--images must be given together with --neb"))
    if args.rmsd is not None and args.conformers is None:
        ap.error(L("--rmsd は --conformers と一緒に指定してください", "--rmsd must be given together with --conformers"))
    if args.conformers is not None and args.rmsd is None:
        ap.error(L("--conformers には --rmsd (重複とみなす RMSD [Å]。既定はありません) が要ります", "--conformers needs --rmsd (duplicate RMSD in Å; no default)"))
    if args.neb and args.images is None:
        ap.error(L("--neb には --images (中間の像の数) が要ります", "--neb needs --images (number of intermediate images)"))
    rc = _run_group(args, spec, cfg, output_dir)
    if rc is not None:
        return rc
    if args.stages:
        from adit.stages import StageError, load_stages, write_stages
        try:
            dirs = write_stages(spec, cfg, output_dir, load_stages(args.stages), overwrite=args.overwrite)
        except OutputNotEmpty as ex:
            print(f"{ex}\n" + L(*_OVERWRITE_HINT), file=sys.stderr)
            return 1
        except (StageError, ProjectError, ConfigError) as ex:
            print(str(ex), file=sys.stderr)
            return 1
        out = Path(output_dir).resolve()
        print(L(f"{len(dirs)} 段階のディレクトリを {out} に作りました: {', '.join(d.name for d in dirs)}\n手順は {out / 'README.txt'}",
                f"created {len(dirs)} stage directories in {out}: {', '.join(d.name for d in dirs)}\nsee {out / 'README.txt'}"))
        return 0
    if args.core or args.substituent:
        from adit.enumerate_r import EnumerateError
        from adit.structures_batch import write_enumerated

        if not (args.core and args.substituent):
            ap.error(L("--core と --substituent は一緒に指定してください", "--core and --substituent must be given together"))
        try:
            dirs = write_enumerated(spec, cfg, output_dir, args.core, args.substituent, overwrite=args.overwrite)
        except OutputNotEmpty as ex:
            print(f"{ex}\n" + L(*_OVERWRITE_HINT), file=sys.stderr)
            return 1
        except (EnumerateError, ProjectError, ConfigError, ValueError) as ex:
            print(str(ex), file=sys.stderr)
            return 1
        out = Path(output_dir).resolve()
        print(L(f"{len(dirs)} 個の置換体のディレクトリを {out} に作りました\n手順は {out / 'README.txt'}",
                f"created {len(dirs)} substituted-structure directories in {out}\nsee {out / 'README.txt'}"))
        return 0
    if args.structures:
        from adit.structures_batch import StructuresError, write_structures
        try:
            dirs = write_structures(spec, cfg, output_dir, args.structures, overwrite=args.overwrite)
        except OutputNotEmpty as ex:
            print(f"{ex}\n" + L(*_OVERWRITE_HINT), file=sys.stderr)
            return 1
        except (StructuresError, ProjectError, ConfigError, ValueError) as ex:
            print(str(ex), file=sys.stderr)
            return 1
        out = Path(output_dir).resolve()
        print(L(f"{len(dirs)} 個の構造のディレクトリを {out} に作りました: {', '.join(d.name for d in dirs[:5])}"
                + (" …" if len(dirs) > 5 else "") + f"\n手順は {out / 'README.txt'}",
                f"created {len(dirs)} structure directories in {out}: {', '.join(d.name for d in dirs[:5])}"
                + (" …" if len(dirs) > 5 else "") + f"\nsee {out / 'README.txt'}"))
        return 0
    if args.scan:
        from adit.scan import GEOM, ScanError, parse_scan, write_scan
        try:
            scans = [parse_scan(text) for text in args.scan]
            dirs = write_scan(spec, cfg, output_dir, scans, overwrite=args.overwrite)
        except OutputNotEmpty as ex:
            print(f"{ex}\n" + L(*_OVERWRITE_HINT), file=sys.stderr)
            return 1
        except (ScanError, ProjectError, ConfigError) as ex:
            print(str(ex), file=sys.stderr)
            return 1
        out = Path(output_dir).resolve()
        shown = ", ".join(d.name for d in dirs[:6]) + (" …" if len(dirs) > 6 else "")
        print(L(f"{len(dirs)} 個のディレクトリを {out} に作りました: {shown}\n"
                f"それぞれを実行したあと、adit-analyze {out} --scan で値とエネルギーの表 (scan_energies.csv) を作れます",
                f"created {len(dirs)} directories in {out}: {shown}\n"
                f"after running each, adit-analyze {out} --scan writes a value-energy table (scan_energies.csv)"))
        if any(GEOM.match(sc.path) for sc in scans) and spec.task.type in ("geometry_optimization", "molecular_dynamics"):
            print(L("注意: 振った幾何を固定していません。構造を動かす計算 (構造最適化・MD) では、その値のまま留まりません。"
                    "値を保ちたいときは、一点計算にするか、その原子を固定してください。",
                    "Note: the scanned geometry is not constrained. In a task that moves the atoms (optimization, MD) it will not stay "
                    "at the scanned value; use a single-point task or fix those atoms if you need it to."), file=sys.stderr)
        return 0
    try:
        written = write_project(spec, cfg, output_dir, overwrite=args.overwrite)
    except OutputNotEmpty as ex:
        print(f"{ex}\n" + L("  中のファイルを上書きしてよければ、同じコマンドの最後に --overwrite を付けて実行すると上書きします。",
                            "  To overwrite the files inside, run the same command again with --overwrite at the end."), file=sys.stderr)
        return 1
    except (ProjectError, ConfigError) as ex:
        print(str(ex), file=sys.stderr)
        return 1
    subdirs: dict[str, int] = {}
    for path in written:
        try:
            rel = Path(path).resolve().relative_to(Path(output_dir).resolve())
        except ValueError:
            continue
        if len(rel.parts) > 1:
            subdirs[rel.parts[0]] = subdirs.get(rel.parts[0], 0) + 1
    inside = (L("、うち " + "・".join(f"{name}/ の中に {n} 件" for name, n in sorted(subdirs.items())),
                ", of which " + ", ".join(f"{n} in {name}/" for name, n in sorted(subdirs.items())))
              if subdirs else "")
    print(L(f"{len(written)} ファイルを {Path(output_dir).resolve()} に書いた{inside}。実行の手順は README.txt",
            f"wrote {len(written)} files to {Path(output_dir).resolve()}{inside}. See README.txt for how to run"))
    if getattr(args, "run", False):
        return _run_here(Path(output_dir), cfg, args)
    return 0


def _run_here(output_dir: Path, cfg, args) -> int:
    from adit.runner import block_reason, run_and_wait, target_from_dir

    target = target_from_dir(output_dir, cfg)
    reason = block_reason(cfg, target, cfg_path=getattr(args, "config", None) or cfg.source_path)
    if reason:
        print(L(f"実行しませんでした: {reason}", f"did not run: {reason}"), file=sys.stderr)
        return 3
    print(L(f"実行します: bash submit.sh ({target.run_dir})   ※ 止めるときは Ctrl-C",
            f"running: bash submit.sh in {target.run_dir}   (Ctrl-C to stop)"))
    try:
        code = run_and_wait(target, on_line=lambda line: print(line, flush=True))
    except KeyboardInterrupt:
        print(L("中断しました (途中までの出力はディレクトリに残ります)",
                "interrupted (the partial output stays in the directory)"), file=sys.stderr)
        return 130
    except OSError as ex:
        print(L(f"実行できません: {ex}", f"cannot run: {ex}"), file=sys.stderr)
        return 3
    print(L(f"終了コード {code}。解析は adit-analyze {target.run_dir}",
            f"exit code {code}. To analyze: adit-analyze {target.run_dir}"))
    return 0 if code == 0 else 4


if __name__ == "__main__":
    sys.exit(main())
