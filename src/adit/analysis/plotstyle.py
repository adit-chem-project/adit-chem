
from __future__ import annotations

from dataclasses import dataclass

from adit.errors import AditValueError
from adit.lang import L

TICK_DIRECTIONS = ("in", "out", "inout")
GRID_CHOICES = ("both", "x", "y", "none")
SPINE_CHOICES = ("all", "left-bottom", "none")


class PlotStyleError(AditValueError):
    pass


@dataclass
class PlotStyle:
    colors: tuple[str, ...] = ()
    tick_direction: str = ""          # in / out / inout
    grid: str = ""                    # both / x / y / none
    spines: str = ""                  # all / left-bottom / none
    line_width: float = 0.0
    font_size: float = 0.0
    dpi: int = 0

    def apply(self) -> None:
        from matplotlib import rcParams
        from matplotlib.colors import is_color_like
        from cycler import cycler

        if self.colors:
            bad = [c for c in self.colors if not is_color_like(c)]
            if bad:
                raise PlotStyleError(L(f"色として読めません: {', '.join(bad)} "
                                       "(例: black、#1f77b4、tab:blue)",
                                       f"cannot read as a color: {', '.join(bad)} "
                                       "(e.g. black, #1f77b4, tab:blue)"))
            rcParams["axes.prop_cycle"] = cycler(color=list(self.colors))
        if self.tick_direction:
            if self.tick_direction not in TICK_DIRECTIONS:
                raise PlotStyleError(L(f"目盛りの向きは {', '.join(TICK_DIRECTIONS)} のどれかです: {self.tick_direction}",
                                       f"the tick direction must be one of {', '.join(TICK_DIRECTIONS)}: {self.tick_direction}"))
            rcParams["xtick.direction"] = rcParams["ytick.direction"] = self.tick_direction
        if self.grid:
            if self.grid not in GRID_CHOICES:
                raise PlotStyleError(L(f"目盛り線は {', '.join(GRID_CHOICES)} のどれかです: {self.grid}",
                                       f"the grid must be one of {', '.join(GRID_CHOICES)}: {self.grid}"))
            rcParams["axes.grid"] = self.grid != "none"
            if self.grid in ("x", "y"):
                rcParams["axes.grid.axis"] = self.grid
            else:
                rcParams["axes.grid.axis"] = "both"
        if self.spines:
            if self.spines not in SPINE_CHOICES:
                raise PlotStyleError(L(f"枠は {', '.join(SPINE_CHOICES)} のどれかです: {self.spines}",
                                       f"the spines must be one of {', '.join(SPINE_CHOICES)}: {self.spines}"))
            show_all = self.spines == "all"
            rcParams["axes.spines.top"] = rcParams["axes.spines.right"] = show_all
            keep = self.spines != "none"
            rcParams["axes.spines.left"] = rcParams["axes.spines.bottom"] = keep
        if self.line_width > 0:
            rcParams["lines.linewidth"] = self.line_width
        if self.font_size > 0:
            rcParams["font.size"] = self.font_size
        if self.dpi > 0:
            rcParams["savefig.dpi"] = self.dpi

    def summary(self) -> str:
        parts = []
        if self.colors:
            parts.append(L("線の色 ", "line colors ") + ", ".join(self.colors))
        if self.tick_direction:
            parts.append(L(f"目盛りの向き {self.tick_direction}", f"tick direction {self.tick_direction}"))
        if self.grid:
            parts.append(L(f"目盛り線 {self.grid}", f"grid {self.grid}"))
        if self.spines:
            parts.append(L(f"枠 {self.spines}", f"spines {self.spines}"))
        if self.line_width > 0:
            parts.append(L(f"線の太さ {self.line_width:g} pt", f"line width {self.line_width:g} pt"))
        if self.font_size > 0:
            parts.append(L(f"文字 {self.font_size:g} pt", f"font size {self.font_size:g} pt"))
        if self.dpi > 0:
            parts.append(L(f"解像度 {self.dpi} dpi", f"resolution {self.dpi} dpi"))
        return L("、", "; ").join(parts)


def from_text(colors: str = "", tick_direction: str = "", grid: str = "", spines: str = "",
              line_width: str | float = "", font_size: str | float = "", dpi: str | int = "") -> PlotStyle:
    def number(value, name: str) -> float:
        if value in ("", None):
            return 0.0
        try:
            out = float(value)
        except (TypeError, ValueError):
            raise PlotStyleError(L(f"{name} は数で書いてください: {value!r}", f"{name} must be a number: {value!r}")) from None
        if out < 0:
            raise PlotStyleError(L(f"{name} は 0 以上にしてください: {out:g}", f"{name} must be 0 or more: {out:g}"))
        return out

    return PlotStyle(
        colors=tuple(c.strip() for c in str(colors or "").split(",") if c.strip()),
        tick_direction=str(tick_direction or "").strip(),
        grid=str(grid or "").strip(),
        spines=str(spines or "").strip(),
        line_width=number(line_width, L("線の太さ", "line width")),
        font_size=number(font_size, L("文字の大きさ", "font size")),
        dpi=int(number(dpi, L("解像度", "resolution"))),
    )

_GRID_CHOICE = ""


def set_grid_choice(choice: str) -> None:
    global _GRID_CHOICE
    _GRID_CHOICE = choice or ""


def grid(ax, alpha: float = 0.3) -> None:
    if _GRID_CHOICE == "none":
        ax.grid(False)
        return
    if _GRID_CHOICE in ("x", "y"):
        ax.grid(True, axis=_GRID_CHOICE, alpha=alpha)
        return
    ax.grid(alpha=alpha)
