
from __future__ import annotations

from adit.errors import AditValueError
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from adit import batch
from adit.lang import L
from adit.spec import AtomsData, CalculationSpec

FORCE_FIELDS = ("MMFF94", "MMFF94s", "UFF")
MAX_CONFORMERS = 5000
DEFAULT_SEED = 12345
CONF_FILE = "conformers.json"
SDF_FILE = "conformers.sdf"
CREST_DOC = "https://crest-lab.github.io/crest-docs/"


class ConformerError(AditValueError):
    pass


@dataclass
class Conformer:
    conf_id: int
    energy_kcal: float
    converged: bool
    kept: bool = False
    duplicate_of: int | None = None
    rmsd: float | None = None
    rank: int | None = None


def molecule_from(spec: CalculationSpec, smiles: str | None = None):
    from rdkit import Chem

    st = spec.structure
    if smiles is None and st.source == "file" and Path(st.source_ref).suffix.lower() in (".mol", ".sdf"):
        mol = Chem.MolFromMolFile(str(Path(st.source_ref).expanduser()), removeHs=False)
        if mol is None:
            raise ConformerError(L(f"{st.source_ref} を RDKit で読めません", f"RDKit cannot read {st.source_ref}"))
        return Chem.AddHs(mol, addCoords=True), f"file:{st.source_ref}"
    text = smiles if smiles is not None else (st.source_ref if st.source == "smiles" else None)
    if text is None:
        raise ConformerError(L("配座を作るには結合の情報が要ります。SMILES から作った構造にするか (--smiles で渡してもよい)、mol / sdf のファイルを使ってください",
                               "conformers need bonding information; use a structure made from SMILES (or pass --smiles), or a mol / sdf file"))
    mol = Chem.MolFromSmiles(text)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ConformerError(L(f"SMILES を解釈できません: {text!r}", f"cannot parse SMILES: {text!r}"))
    return Chem.AddHs(mol), text


def generate(mol, n: int, rmsd: float, *, seed: int = DEFAULT_SEED, force_field: str = "MMFF94", heavy_only: bool = True,
             max_iters: int = 200) -> list[Conformer]:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign

    if not 1 <= n <= MAX_CONFORMERS:
        raise ConformerError(L(f"作る配座の数は 1〜{MAX_CONFORMERS} にしてください", f"the number of conformers must be 1..{MAX_CONFORMERS}"))
    if rmsd < 0:
        raise ConformerError(L("重複とみなす RMSD は 0 以上にしてください", "the duplicate RMSD must be 0 or more"))
    if force_field not in FORCE_FIELDS:
        raise ConformerError(L(f"力場は {' / '.join(FORCE_FIELDS)} のどれかです", f"force field must be one of {' / '.join(FORCE_FIELDS)}"))
    if max_iters < 1:
        raise ConformerError(L("最適化の反復の上限は 1 以上にしてください", "max iterations must be at least 1"))
    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=int(n), params=params))
    if not ids:
        raise ConformerError(L("RDKit の ETKDG で 3 次元座標を 1 つも作れませんでした", "RDKit ETKDG could not embed any conformer"))
    if len(ids) > 1:
        x0 = mol.GetConformer(ids[0]).GetPositions()
        if max(float(np.abs(mol.GetConformer(i).GetPositions() - x0).max()) for i in ids[1:]) < 1e-6:
            raise ConformerError(L(f"RDKit が作った {len(ids)} 個の配座がすべて同じ座標です (乱数の種 {seed})。別の種にしてください "
                                   "(RDKit 2026.03 では種 0 でこうなりました)",
                                   f"all {len(ids)} conformers embedded by RDKit have identical coordinates (seed {seed}); use another seed "
                                   "(RDKit 2026.03 does this with seed 0)"))
    if force_field.startswith("MMFF"):
        if not AllChem.MMFFHasAllMoleculeParams(mol):
            raise ConformerError(L(f"この分子には {force_field} のパラメータがそろっていません (UFF を選ぶ方法があります)", f"{force_field} parameters are missing for this molecule (UFF is an alternative)"))
        res = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=int(max_iters), mmffVariant=force_field)
    else:
        if not AllChem.UFFHasAllMoleculeParams(mol):
            raise ConformerError(L("この分子には UFF のパラメータがそろっていません", "UFF parameters are missing for this molecule"))
        res = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=int(max_iters))
    confs = [Conformer(c.GetId(), float(e), not bool(nc)) for c, (nc, e) in zip(mol.GetConformers(), res)]
    confs.sort(key=lambda c: c.energy_kcal)
    ref = Chem.RemoveHs(mol) if heavy_only else Chem.Mol(mol)
    kept: list[Conformer] = []
    for c in confs:
        probe = Chem.Mol(ref)
        for k in kept:
            r = float(rdMolAlign.GetBestRMS(probe, ref, c.conf_id, k.conf_id))
            if r <= rmsd:
                c.duplicate_of, c.rmsd = k.conf_id, r
                break
        if c.duplicate_of is None:
            c.kept, c.rank = True, len(kept) + 1
            kept.append(c)
    return confs


def _atoms_of(mol, conf_id: int):
    from ase import Atoms

    conf = mol.GetConformer(conf_id)
    return Atoms(symbols=[a.GetSymbol() for a in mol.GetAtoms()], positions=conf.GetPositions())


def _crest_files(spec: CalculationSpec, atoms) -> dict[str, str] | None:
    import io

    from ase.io import write

    m = spec.method
    if m.code != "xtb" or m.gfn != "2" or m.solvation == "gbsa":
        return None
    args = ["crest", "struct.xyz", "--gfn2", f"--chrg {spec.structure.charge}", f"--uhf {spec.structure.multiplicity - 1}"]
    if m.solvation == "alpb":
        args.append(f"--alpb {m.solvent.strip().lower()}")
    args.append(f"-T {spec.runtime.omp_threads}")
    buf = io.StringIO()
    write(buf, atoms, format="xyz")
    run = ("#!/bin/bash\n" + L("# ADIT が生成。CREST (xtb によるメタダイナミクスの配座探索) を実行する。ADIT は実行しない\n",
                               "# generated by ADIT. Runs CREST (metadynamics conformer search with xtb); ADIT does not run it\n")
           + 'cd "$(dirname "$0")"\n' + " ".join(args) + " > crest.log 2>&1\n")
    readme = "\n".join([
        L("CREST の入力 (ADIT は実行しません)", "CREST input (ADIT does not run it)"), "",
        L("  struct.xyz   出発の構造 (力場のエネルギーが最も低い配座)", "  struct.xyz   starting structure (the conformer with the lowest force-field energy)"),
        f"  run_crest.sh {' '.join(args)}",
        L("               電荷・不対電子の数・溶媒・並列数は spec.json の xtb の設定から写したものです", "               charge, unpaired electrons, solvent and threads are copied from the xtb settings in spec.json"),
        "", L("インストール方法: conda install -c conda-forge crest   (xtb と同じ conda-forge)", "Install: conda install -c conda-forge crest   (conda-forge, like xtb)"),
        L(f"結果: crest_conformers.xyz (重複を除いた配座)、crest_best.xyz (最もエネルギーの低いもの)、crest.energies (相対エネルギー)。文書: {CREST_DOC}",
          f"Results: crest_conformers.xyz (unique conformers), crest_best.xyz (lowest), crest.energies (relative energies). Docs: {CREST_DOC}"),
    ]) + "\n"
    return {"struct.xyz": buf.getvalue(), "run_crest.sh": run, "README.txt": readme}


def write_conformers(spec: CalculationSpec, cfg, out_dir: Path | str, *, n: int, rmsd: float, seed: int = DEFAULT_SEED, force_field: str = "MMFF94",
                     heavy_only: bool = True, max_iters: int = 200, keep: int | None = None, smiles: str | None = None,
                     overwrite: bool = False) -> list[Path]:
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    from adit import __version__

    out = Path(out_dir).expanduser()
    if spec.structure.periodic:
        raise ConformerError(L("配座の候補は分子 (非周期) の構造で作ります", "conformers are generated for molecules (non-periodic)"))
    if keep is not None and keep < 1:
        raise ConformerError(L("残す配座の数は 1 以上にしてください", "the number of conformers to keep must be at least 1"))
    mol, source = molecule_from(spec, smiles)
    fc, rad = Chem.GetFormalCharge(mol), Descriptors.NumRadicalElectrons(mol)
    if fc != spec.structure.charge:
        raise ConformerError(L(f"分子の形式電荷の和 ({fc}) が spec.json の全電荷 ({spec.structure.charge}) と違います", f"the formal charge of the molecule ({fc}) differs from the total charge in spec.json ({spec.structure.charge})"))
    if rad + 1 != spec.structure.multiplicity:
        raise ConformerError(L(f"分子の不対電子 ({rad} 個、多重度 {rad + 1}) が spec.json の多重度 ({spec.structure.multiplicity}) と違います",
                               f"the molecule has {rad} radical electrons (multiplicity {rad + 1}) but spec.json has multiplicity {spec.structure.multiplicity}"))
    if spec.structure.fixed_atoms:
        mine = [a.GetSymbol() for a in mol.GetAtoms()]
        if list(spec.structure.atoms.symbols) != mine:
            raise ConformerError(L(f"固定原子の番号が、配座の原子の並びと合いません。spec.json は {len(spec.structure.atoms.symbols)} 原子 ("
                                   f"{''.join(list(spec.structure.atoms.symbols)[:8])}…)、配座は {len(mine)} 原子 ({''.join(mine[:8])}…) です。"
                                   "RDKit は SMILES から水素を足した順に原子を並べるので、番号が違うことがあります",
                                   f"the fixed-atom indices do not match the conformer atom order (spec.json has {len(spec.structure.atoms.symbols)} atoms "
                                   f"{''.join(list(spec.structure.atoms.symbols)[:8])}..., the conformers have {len(mine)} atoms {''.join(mine[:8])}...; "
                                   "RDKit orders atoms as SMILES plus added hydrogens, so the indices can differ)"))
    confs = generate(mol, n, rmsd, seed=seed, force_field=force_field, heavy_only=heavy_only, max_iters=max_iters)
    kept = [c for c in confs if c.kept][:keep] if keep else [c for c in confs if c.kept]
    e0 = kept[0].energy_kcal
    items = []
    base = spec.model_dump(mode="json")
    base["handoff"] = None
    for c in kept:
        data = dict(base)
        data["structure"] = dict(base["structure"], atoms=AtomsData.from_ase(_atoms_of(mol, c.conf_id)).model_dump(mode="json"), velocities=None)
        data["meta"] = dict(base["meta"], comment=(base["meta"].get("comment", "") + f" conformer {c.rank} (RDKit id {c.conf_id}, {force_field} {c.energy_kcal:.4f} kcal/mol)").strip())
        s = CalculationSpec.model_validate(data)
        items.append(batch.Item(f"conf_{c.rank:03d}", s, [
            L("== 配座の候補 ==", "== Conformer candidate =="),
            L(f"  配座 {c.rank} / {len(kept)} (RDKit の番号 {c.conf_id})。{force_field} のエネルギー {c.energy_kcal:.4f} kcal/mol "
              f"(最も低い配座との差 {c.energy_kcal - e0:.4f} kcal/mol)、力場の最適化は{'収束' if c.converged else '反復の上限で打ち切り'}。一覧は ../conformers.json",
              f"  Conformer {c.rank} of {len(kept)} (RDKit id {c.conf_id}); {force_field} energy {c.energy_kcal:.4f} kcal/mol "
              f"({c.energy_kcal - e0:.4f} kcal/mol above the lowest); force-field optimization {'converged' if c.converged else 'stopped at the iteration limit'}. List: ../conformers.json")]))
    crest = _crest_files(spec, _atoms_of(mol, kept[0].conf_id))
    dirs = batch.write_items(items, cfg, out, overwrite=overwrite)
    table = [{"rdkit_id": c.conf_id, "energy_kcal_mol": c.energy_kcal, "relative_kcal_mol": c.energy_kcal - e0, "converged": c.converged,
              "kept": c.kept and c in kept, "rank": c.rank if c in kept else None, "dir": f"conf_{c.rank:03d}" if c in kept else None,
              "duplicate_of_rdkit_id": c.duplicate_of, "rmsd_to_duplicate_ang": c.rmsd} for c in confs]
    batch.write_json(out / CONF_FILE, {
        "generated_by": f"adit {__version__}", "source": source, "requested": n, "embedded": len(confs), "kept": len(kept),
        "seed": seed, "force_field": force_field, "max_iters": max_iters, "rmsd_threshold_ang": rmsd,
        "rmsd_atoms": "heavy" if heavy_only else "all", "energy_unit": "kcal/mol (force field)", "conformers": table})
    w = Chem.SDWriter(str(out / SDF_FILE))
    for c in kept:
        mol.SetProp("_Name", f"conf_{c.rank:03d}")
        mol.SetProp("adit_energy_kcal_mol", f"{c.energy_kcal:.6f}")
        mol.SetProp("adit_force_field", force_field)
        w.write(mol, confId=c.conf_id)
    w.close()
    batch.write_scan_json(out, "conformer", [f"{c.rank:03d}" for c in kept], dirs)
    if crest:
        (out / "crest").mkdir(exist_ok=True)
        for name, text in crest.items():
            (out / "crest" / name).write_text(text, encoding="utf-8", newline="\n")
        (out / "crest" / "run_crest.sh").chmod(0o755)
    dup = sum(1 for c in confs if not c.kept)
    lines = [L(f"ADIT {__version__} が生成した、配座の候補の計算です ({spec.method.code})", f"Conformer candidates generated by ADIT {__version__} ({spec.method.code})"), "",
             L(f"  元の分子: {source}", f"  molecule: {source}"),
             L(f"  RDKit の ETKDG で {n} 個を頼み {len(confs)} 個を作り (乱数の種 {seed})、{force_field} で最適化 (反復の上限 {max_iters})。",
               f"  Requested {n} conformers from RDKit ETKDG, embedded {len(confs)} (seed {seed}), optimized with {force_field} (max {max_iters} iterations)."),
             L(f"  RMSD {rmsd:g} Å 以下 ({'水素を除いた原子。水素の向きだけが違う配座は重複になります' if heavy_only else '全原子'}) を重複として {dup} 個を外し、{len(kept)} 個を計算にしました。",
               f"  {dup} were dropped as duplicates (RMSD <= {rmsd:g} Å over {'heavy atoms; conformers differing only in hydrogen orientation count as duplicates' if heavy_only else 'all atoms'}); {len(kept)} became runs."),
             L("  力場のエネルギーは並べるだけです (どの配座を使うかは ADIT は判断しません)。", "  Force-field energies are listed only; ADIT does not choose among them."), "",
             L("  ディレクトリ  力場のエネルギー [kcal/mol]  最低との差", "  directory    force-field energy [kcal/mol]  above lowest")]
    lines += [f"  conf_{c.rank:03d}      {c.energy_kcal:12.4f}            {c.energy_kcal - e0:8.4f}" for c in kept]
    lines += ["", L("  conformers.json  全部の配座 (重複として外したものと、その相手と RMSD を含む)", "  conformers.json  all conformers (with the dropped duplicates, their partner and RMSD)"),
              L("  conformers.sdf   残した配座 (分子ビューアで開けます)", "  conformers.sdf   the kept conformers (opens in molecular viewers)")]
    if crest:
        lines.append(L("  crest/           CREST の入力 (ADIT は実行しません。crest/README.txt)", "  crest/           CREST input (ADIT does not run it; see crest/README.txt)"))
    else:
        lines.append(L("  CREST の入力は、計算コードが xtb で GFN2 (溶媒は ALPB かなし) のときだけ書きます (CREST の文書で確かめた引数の組のため)。",
                       "  CREST input is written only for xtb with GFN2 (solvent ALPB or none), the argument set verified in the CREST docs."))
    lines += ["", L("== 実行したあと ==", "== After running =="),
              L("  adit-analyze <このディレクトリ> --scan   (配座ごとの最終エネルギーの表 scan_energies.csv)", "  adit-analyze <this directory> --scan   (final energy per conformer in scan_energies.csv)")]
    batch.write_top(out, cfg, spec, dirs, lines, "配座の計算", "conformer runs")
    return dirs


__all__ = ["ConformerError", "Conformer", "FORCE_FIELDS", "DEFAULT_SEED", "molecule_from", "generate", "write_conformers"]
