
from __future__ import annotations

from typing import Annotated, Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from adit.lang import L
from adit.mixture import Component
from adit.structure import StructureError

Vec3 = tuple[float, float, float]


class RecipeError(StructureError):
    pass


class _M(BaseModel):
    model_config = ConfigDict(extra="forbid")


BASE_SOURCES = ("preset", "smiles", "file", "bulk", "surface", "mixture", "2d", "cluster", "polymer")


class Base(_M):

    source: Literal["preset", "smiles", "file", "bulk", "surface", "mixture", "2d", "cluster", "polymer"]
    ref: str | dict[str, Any]
    sha256: str = ""


class TwoD(_M):

    kind: Literal["graphene", "mx2", "nanoribbon", "nanotube"]
    formula: str = ""
    a: float | None = None
    thickness: float | None = None
    mx2_kind: Literal["2H", "1T"] = "2H"
    size: tuple[int, int, int] = (1, 1, 1)
    n: int = 6
    m: int = 0
    length: int = 1
    bond: float = 1.42
    symbol: str = "C"
    ribbon_type: Literal["zigzag", "armchair"] = "zigzag"
    saturated: bool = True
    vacuum: float = 10.0


class ClusterRef(_M):

    kind: Literal["icosahedron", "decahedron", "octahedron", "wulff"]
    symbol: str
    shells: int = 3
    p: int = 2
    q: int = 2
    r: int = 0
    length: int = 5
    cutoff: int = 0
    size: int = 100
    surfaces: list[tuple[int, int, int]] = Field(default_factory=lambda: [(1, 0, 0), (1, 1, 0), (1, 1, 1)])
    energies: list[float] = Field(default_factory=lambda: [1.0, 1.1, 0.9])
    structure: Literal["fcc", "bcc", "sc"] = "fcc"
    lattice_constant: float | None = None


class PolymerRef(_M):

    unit: str
    n: int = 10
    seed: int = 0


class Selection(_M):
    elements: list[str] = Field(default_factory=list)
    indices: list[int] = Field(default_factory=list)
    z_min: float | None = None
    z_max: float | None = None


class MoleculeRef(_M):

    kind: Literal["preset", "smiles", "file"] = "preset"
    ref: str


class _Step(_M):
    # A disabled step stays in the recipe but is skipped when building (OVITO-style pipeline toggle).
    enabled: bool = True


class Supercell(_Step):

    op: Literal["supercell"] = "supercell"
    repeat: tuple[int, int, int] | None = None
    matrix: list[list[int]] | None = None
    fit_components: list[Component] = Field(default_factory=list)
    fit_min_distance: float = 2.0

    @model_validator(mode="after")
    def _one(self) -> "Supercell":
        if (self.repeat is None) == (self.matrix is None):
            raise ValueError(L("repeat と matrix のどちらか一方を書いてください", "give exactly one of repeat or matrix"))
        if self.matrix is not None and (len(self.matrix) != 3 or any(len(r) != 3 for r in self.matrix)):
            raise ValueError(L("matrix は 3×3 の整数です", "matrix must be 3x3 integers"))
        if self.matrix is not None and self.fit_components:
            raise ValueError(L("溶液に合わせた断面の自動調整には、「変換行列」ではなく「各軸の繰り返し」を使ってください",
                               "automatic cross-section sizing for a solution requires axis repeats, not a transformation matrix"))
        return self


class Slab(_Step):

    op: Literal["slab"] = "slab"
    miller: tuple[int, int, int]
    layers: int = 3
    vacuum: float = 10.0
    termination: int = 0


class Vacuum(_Step):

    op: Literal["vacuum"] = "vacuum"
    axis: Literal[0, 1, 2] = 2
    thickness: float = 10.0


class Box(_Step):

    op: Literal["box"] = "box"
    padding: float = 5.0
    lengths: Vec3 | None = None
    max_multiple: int = 8


class Adsorb(_Step):

    op: Literal["adsorb"] = "adsorb"
    molecule: MoleculeRef
    site: str | None = None
    above_atom: int | None = None
    xy: tuple[float, float] | None = None
    height: float = 2.0
    down_atom: int = 0

    @model_validator(mode="after")
    def _one(self) -> "Adsorb":
        if sum(x is not None for x in (self.site, self.above_atom, self.xy)) != 1:
            raise ValueError(L("置き場所は site / above_atom / xy のどれか 1 つを書いてください",
                               "give exactly one of site / above_atom / xy"))
        return self


class _Pick(_Step):
    where: Selection = Field(default_factory=Selection)
    count: int | None = None
    fraction: float | None = None
    seed: int = 0

    @model_validator(mode="after")
    def _one(self):
        if (self.count is None) == (self.fraction is None):
            raise ValueError(L("count (個数) と fraction (割合) のどちらか一方を書いてください", "give exactly one of count or fraction"))
        return self


class Remove(_Pick):

    op: Literal["remove"] = "remove"


class Substitute(_Pick):

    op: Literal["substitute"] = "substitute"
    to: str


class SolventLayer(_Step):

    op: Literal["solvent_layer"] = "solvent_layer"
    components: list[Component]
    density_g_cm3: float = 1.0
    thickness: float | None = None
    gap: float = 2.0
    vacuum: float = 0.0
    min_distance: float = 2.0
    seed: int = 0
    max_tries: int = 2000


class Solvate(_Step):

    op: Literal["solvate"] = "solvate"
    components: list[Component]
    padding: float = 8.0
    density_g_cm3: float = 1.0
    min_distance: float = 2.0
    seed: int = 0
    max_tries: int = 2000


class Fix(_Step):

    op: Literal["fix"] = "fix"
    where: Selection = Field(default_factory=Selection)
    bottom_layers: int = 0
    layer_tolerance: float = 0.5


Step = Annotated[Supercell | Slab | Vacuum | Box | Adsorb | Remove | Substitute | SolventLayer | Solvate | Fix,
                 Field(discriminator="op")]

OP_LABELS: dict[str, tuple[str, str]] = {
    "base": ("土台", "base"), "supercell": ("超格子", "supercell"), "slab": ("面で切る", "slab"),
    "vacuum": ("真空層", "vacuum"), "box": ("直方体のセル", "box"), "adsorb": ("吸着", "adsorb"),
    "remove": ("原子を抜く", "remove"), "substitute": ("元素の置換", "substitute"),
    "solvent_layer": ("溶液の層 (界面)", "solvent layer (interface)"), "solvate": ("溶媒和", "solvate"), "fix": ("固定", "fix"),
}


def op_label(op: str) -> str:
    ja, en = OP_LABELS.get(op, (op, op))
    return L(ja, en)


def interface_steps(base_is_slab: bool, components: list[Component], *, miller: tuple[int, int, int] = (1, 0, 0), layers: int = 3,
                    min_distance: float = 2.0) -> list[Any]:
    comps = [c.model_copy() for c in components]
    steps: list[Any] = [] if base_is_slab else [Slab(miller=tuple(int(x) for x in miller), layers=layers, vacuum=0.0)]
    steps.append(Supercell(repeat=(1, 1, 1), fit_components=comps, fit_min_distance=min_distance))
    steps.append(SolventLayer(components=[c.model_copy() for c in comps], vacuum=0.0, min_distance=min_distance))
    steps.append(Fix(bottom_layers=1))
    return steps


class Recipe(_M):
    base: Base
    steps: list[Step] = Field(default_factory=list)

    def to_ref(self) -> str:
        import json

        data = self.model_dump(mode="json", exclude_none=True)
        for step in data["steps"]:
            if step.get("enabled", True):
                step.pop("enabled", None)     # only "enabled": false is written, so older refs compare equal
        return json.dumps(data, ensure_ascii=False, separators=(",", ":"))

    def active_steps(self) -> list[Any]:
        return [s for s in self.steps if s.enabled]

    @classmethod
    def from_ref(cls, ref: str) -> "Recipe":
        try:
            return cls.model_validate_json(ref)
        except ValueError as ex:
            errs = getattr(ex, "errors", None)
            if callable(errs):
                lines = "; ".join(f"{'.'.join(str(x) for x in e.get('loc', ()))}: {e.get('msg', '')}" for e in errs())
            else:
                lines = str(ex)
            raise RecipeError(L(f"組み立て手順を読めません: {lines}", f"cannot read the build recipe: {lines}")) from ex

    def total_charge(self) -> int:
        from adit.mixture import MixtureSpec

        q = 0
        if self.base.source == "mixture" and isinstance(self.base.ref, str):
            q += MixtureSpec.from_ref(self.base.ref).total_charge()
        for s in self.active_steps():
            if isinstance(s, (SolventLayer, Solvate)):
                q += sum(c.count * c.charge for c in s.components)
        return q


def has_op(steps: list, op: str) -> bool:
    """True when an enabled step of this kind is in the list."""
    return any(s.op == op and getattr(s, "enabled", True) for s in steps)
