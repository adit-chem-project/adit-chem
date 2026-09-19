
from __future__ import annotations

import math

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QFont, QIcon, QPainter, QPainterPath, QPen, QPixmap

from adit.gui.style import GROUP_COLORS, Tokens, resolve_theme, tokens_for

_UNIT = 100.0
_CACHE: dict[tuple[str, int, str], QIcon] = {}
_TOKENS: Tokens | None = None
_THEME = "light"


def set_theme(theme: str) -> None:
    global _TOKENS, _THEME
    _THEME = "dark" if str(theme).endswith("dark") else "light"
    _TOKENS = tokens_for(_THEME)
    _CACHE.clear()


def _tokens() -> Tokens:
    global _TOKENS
    if _TOKENS is None:
        set_theme(resolve_theme())
    return _TOKENS  # type: ignore[return-value]


def _pen(p: QPainter, color: str, width: float = 8.0) -> None:
    pen = QPen(QColor(color)); pen.setWidthF(width)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap); pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
    p.setPen(pen); p.setBrush(Qt.BrushStyle.NoBrush)


def _fill(p: QPainter, color: str) -> None:
    p.setPen(Qt.PenStyle.NoPen); p.setBrush(QColor(color))


def _dot(p: QPainter, x: float, y: float, r: float, color: str) -> None:
    _fill(p, color); p.drawEllipse(QPointF(x, y), r, r)


def _line(p: QPainter, x1: float, y1: float, x2: float, y2: float, color: str, w: float = 8.0) -> None:
    _pen(p, color, w); p.drawLine(QPointF(x1, y1), QPointF(x2, y2))


def _poly(p: QPainter, pts: list[tuple[float, float]], color: str, w: float = 8.0, close: bool = True) -> None:
    path = QPainterPath(QPointF(*pts[0]))
    for x, y in pts[1:]:
        path.lineTo(QPointF(x, y))
    if close:
        path.closeSubpath()
    _pen(p, color, w); p.drawPath(path)


def _badge(p: QPainter, letter: str, color: str) -> None:
    _fill(p, color); p.drawRoundedRect(QRectF(8, 8, 84, 84), 24, 24)
    f = QFont(); f.setPixelSize(62); f.setBold(True); p.setFont(f)
    p.setPen(QColor("#FFFFFF"))
    p.drawText(QRectF(8, 6, 84, 88), Qt.AlignmentFlag.AlignCenter, letter)


def _draw(name: str, p: QPainter, t: Tokens) -> None:
    fg, accent, muted = t.fg, t.accent, t.muted

    if name == "preset":
        _line(p, 30, 68, 50, 40, muted, 9); _line(p, 70, 68, 50, 40, muted, 9)
        _dot(p, 50, 38, 20, accent); _dot(p, 26, 72, 14, fg); _dot(p, 74, 72, 14, fg)
    elif name == "smiles":
        _poly(p, [(14, 66), (34, 44), (54, 66), (74, 44), (90, 56)], accent, 9, close=False)
        _dot(p, 14, 66, 9, fg); _dot(p, 54, 66, 9, fg); _dot(p, 90, 56, 9, fg)
    elif name == "file":
        _poly(p, [(24, 10), (62, 10), (80, 30), (80, 90), (24, 90)], fg, 8)
        _poly(p, [(62, 10), (62, 30), (80, 30)], muted, 7, close=False)
        for y in (52, 68):
            _line(p, 36, y, 68, y, muted, 7)
    elif name == "bulk":
        _poly(p, [(20, 34), (50, 18), (80, 34), (80, 68), (50, 84), (20, 68)], fg, 8)
        _poly(p, [(20, 34), (50, 50), (80, 34)], muted, 7, close=False)
        _line(p, 50, 50, 50, 84, muted, 7)
    elif name == "slab":
        _fill(p, accent); p.drawRoundedRect(QRectF(14, 58, 72, 12), 5, 5)
        _fill(p, fg); p.drawRoundedRect(QRectF(14, 74, 72, 12), 5, 5)
        _pen(p, muted, 6); pen = p.pen(); pen.setStyle(Qt.PenStyle.DotLine); p.setPen(pen)
        p.drawLine(QPointF(14, 30), QPointF(86, 30)); p.drawLine(QPointF(14, 44), QPointF(86, 44))
    elif name == "mixture":
        _poly(p, [(30, 12), (30, 44), (16, 84), (84, 84), (70, 44), (70, 12)], fg, 8)
        _dot(p, 38, 70, 8, accent); _dot(p, 58, 72, 7, accent); _dot(p, 48, 56, 6, muted)
        _line(p, 26, 12, 74, 12, fg, 8)

    elif name == "2d":
        import math as _m
        pts = [(50 + 34 * _m.cos(_m.radians(30 + 60 * k)), 50 + 34 * _m.sin(_m.radians(30 + 60 * k))) for k in range(6)]
        _poly(p, pts, fg, 8)
        for k in (0, 2, 4):
            _dot(p, *pts[k], 9, accent)
    elif name == "cluster":
        for x, y in ((50, 30), (32, 46), (68, 46), (50, 50), (38, 68), (62, 68)):
            _dot(p, x, y, 12, accent if (x, y) == (50, 50) else fg)
    elif name == "polymer":
        _poly(p, [(26, 62), (40, 40), (54, 62), (68, 40)], accent, 8, close=False)
        _poly(p, [(20, 24), (12, 24), (12, 78), (20, 78)], fg, 7, close=False)
        _poly(p, [(80, 24), (88, 24), (88, 78), (80, 78)], fg, 7, close=False)

    elif name == "task_single_point":
        _pen(p, muted, 7); p.drawEllipse(QPointF(50, 50), 32, 32)
        _dot(p, 50, 50, 13, accent)
    elif name == "task_optimize":
        path = QPainterPath(QPointF(14, 24)); path.quadTo(QPointF(50, 96), QPointF(86, 40))
        _pen(p, muted, 8); p.drawPath(path)
        _dot(p, 47, 68, 13, accent)
    elif name == "task_md":
        path = QPainterPath(QPointF(12, 62)); path.cubicTo(QPointF(34, 14), QPointF(52, 92), QPointF(72, 44))
        _pen(p, muted, 8); p.drawPath(path)
        _dot(p, 78, 38, 14, accent)
    elif name == "task_vibrations":
        _poly(p, [(30, 50), (40, 32), (50, 68), (60, 32), (70, 50)], muted, 7, close=False)
        _dot(p, 18, 50, 14, accent); _dot(p, 82, 50, 14, accent)
    elif name == "task_bands":
        _pen(p, muted, 6); p.drawRect(QRectF(16, 16, 68, 68))
        path = QPainterPath(QPointF(20, 34)); path.quadTo(QPointF(50, 8), QPointF(80, 34))
        _pen(p, accent, 8); p.drawPath(path)
        path = QPainterPath(QPointF(20, 66)); path.quadTo(QPointF(50, 92), QPointF(80, 66))
        _pen(p, fg, 8); p.drawPath(path)

    elif name.startswith("code_"):
        letter, color = {"code_dftbplus": ("D", GROUP_COLORS["method"]), "code_vasp": ("V", GROUP_COLORS["runtime"]),
                         "code_xtb": ("x", GROUP_COLORS["task"]), "code_espresso": ("Q", GROUP_COLORS["kpoints"]),
                         "code_orca": ("O", GROUP_COLORS["analysis"]), "code_cp2k": ("C", GROUP_COLORS["structure"]),
                         "code_lammps": ("L", GROUP_COLORS["task"]), "code_gromacs": ("G", GROUP_COLORS["kpoints"]),
                         "code_mlip": ("M", GROUP_COLORS["analysis"])}[name]
        _badge(p, letter, color)

    elif name == "kpoints_gamma":
        _pen(p, muted, 6); p.drawRect(QRectF(16, 16, 68, 68))
        _dot(p, 50, 50, 14, accent)
    elif name == "kpoints_mesh":
        _pen(p, muted, 6); p.drawRect(QRectF(16, 16, 68, 68))
        for x in (30, 50, 70):
            for y in (30, 50, 70):
                _dot(p, x, y, 7, accent if (x == 50 and y == 50) else fg)
    elif name == "kpoints_density":
        _pen(p, muted, 6); p.drawRect(QRectF(16, 16, 68, 68))
        for x, r in ((30, 5), (50, 8), (70, 11)):
            for y in (32, 52, 72):
                _dot(p, x, y, r, accent if x == 70 else fg)

    elif name == "profile_direct":
        _pen(p, fg, 8); p.drawRoundedRect(QRectF(14, 20, 72, 48), 8, 8)
        _line(p, 40, 82, 60, 82, fg, 8); _line(p, 50, 68, 50, 82, fg, 8)
        _fill(p, accent); p.drawRoundedRect(QRectF(24, 30, 52, 28), 4, 4)
    elif name in ("profile_pbs", "profile_slurm"):
        for y in (16, 42, 68):
            _pen(p, fg, 7); p.drawRoundedRect(QRectF(16, y, 68, 18), 5, 5)
            _dot(p, 28, y + 9, 5, accent)

    elif name == "count":
        f = QFont(); f.setPixelSize(58); f.setBold(True); p.setFont(f)
        p.setPen(QColor(fg)); p.drawText(QRectF(0, 0, 100, 100), Qt.AlignmentFlag.AlignCenter, "×n")
    elif name == "box_density":
        _pen(p, muted, 6); p.drawRect(QRectF(16, 16, 68, 68))
        for x in (32, 50, 68):
            for y in (32, 50, 68):
                _dot(p, x, y, 6, accent)
    elif name == "box_edge":
        _pen(p, fg, 7); p.drawRect(QRectF(20, 20, 60, 60))
        _pen(p, accent, 7); p.drawLine(QPointF(20, 90), QPointF(80, 90))
        _line(p, 20, 84, 20, 96, accent, 6); _line(p, 80, 84, 80, 96, accent, 6)

    elif name == "generate":
        _poly(p, [(20, 12), (58, 12), (76, 32), (76, 88), (20, 88)], fg, 8)
        _line(p, 50, 44, 50, 68, accent, 9); _line(p, 38, 56, 62, 56, accent, 9)
    elif name == "run":
        _fill(p, accent)
        path = QPainterPath(QPointF(30, 18)); path.lineTo(QPointF(84, 50)); path.lineTo(QPointF(30, 82)); path.closeSubpath()
        p.drawPath(path)
    elif name == "draw":
        _poly(p, [(20, 80), (26, 60), (66, 20), (80, 34), (40, 74)], accent, 8)
        _line(p, 20, 80, 40, 74, accent, 8)
    elif name == "open":
        _poly(p, [(14, 26), (44, 26), (52, 38), (86, 38), (86, 78), (14, 78)], fg, 8)
    elif name == "save":
        _pen(p, fg, 8); p.drawRoundedRect(QRectF(16, 16, 68, 68), 8, 8)
        _fill(p, accent); p.drawRect(QRectF(34, 16, 32, 26)); p.drawRect(QRectF(30, 58, 40, 26))
    elif name == "settings":
        _fill(p, fg)
        for k in range(6):
            p.save(); p.translate(50, 50); p.rotate(k * 60)
            p.drawRoundedRect(QRectF(-9, -46, 18, 24), 4, 4); p.restore()
        _fill(p, fg); p.drawEllipse(QPointF(50, 50), 28, 28)
        _fill(p, t.card); p.drawEllipse(QPointF(50, 50), 13, 13)
    elif name == "reload":
        _pen(p, fg, 9)
        p.drawArc(QRectF(20, 20, 60, 60), 40 * 16, 280 * 16)
        _fill(p, fg)
        path = QPainterPath(QPointF(80, 22)); path.lineTo(QPointF(84, 52)); path.lineTo(QPointF(56, 40)); path.closeSubpath()
        p.drawPath(path)
    elif name == "quit":
        _pen(p, fg, 9); p.drawArc(QRectF(20, 20, 60, 60), 60 * 16, 240 * 16)
        _line(p, 50, 12, 50, 46, accent, 9)
    elif name in ("undo", "redo"):
        p.save()
        if name == "redo":
            p.translate(100, 0); p.scale(-1, 1)
        path = QPainterPath(QPointF(30, 40)); path.lineTo(QPointF(62, 40))
        path.cubicTo(QPointF(90, 40), QPointF(90, 80), QPointF(62, 80)); path.lineTo(QPointF(44, 80))
        _pen(p, fg, 10); p.drawPath(path)
        _poly(p, [(40, 24), (22, 40), (40, 56)], accent, 10, close=False)
        p.restore()
    elif name == "scan":
        for k, c in enumerate((muted, muted, fg)):
            x, y = 14 + k * 12, 10 + k * 12
            _fill(p, t.card); p.drawRoundedRect(QRectF(x, y, 50, 64), 6, 6)
            _pen(p, c, 7); p.drawRoundedRect(QRectF(x, y, 50, 64), 6, 6)
        _line(p, 48, 62, 76, 62, accent, 8); _line(p, 48, 74, 68, 74, accent, 8)
    elif name == "continue":
        _fill(p, t.card); p.drawRoundedRect(QRectF(10, 20, 36, 48), 5, 5)
        _pen(p, muted, 6); p.drawRoundedRect(QRectF(10, 20, 36, 48), 5, 5)
        _fill(p, t.card); p.drawRoundedRect(QRectF(54, 32, 36, 48), 5, 5)
        _pen(p, fg, 7); p.drawRoundedRect(QRectF(54, 32, 36, 48), 5, 5)
        _line(p, 30, 84, 62, 84, accent, 7); _poly(p, [(56, 76), (66, 84), (56, 92)], accent, 6, close=False)
    elif name == "stages":
        _poly(p, [(12, 84), (12, 66), (36, 66), (36, 46), (60, 46), (60, 26), (88, 26)], fg, 8, close=False)
        for x, y in ((24, 66), (48, 46), (74, 26)):
            _dot(p, x, y - 10, 6, accent)
    elif name == "batch_compare":
        _line(p, 50, 18, 50, 74, fg, 7); _line(p, 20, 30, 80, 30, fg, 7)
        _line(p, 26, 78, 74, 78, fg, 7)
        _poly(p, [(10, 52), (20, 30), (30, 52)], accent, 6); _poly(p, [(70, 52), (80, 30), (90, 52)], muted, 6)
    elif name == "batch_conformers":
        _poly(p, [(16, 70), (30, 44), (46, 66)], muted, 7, close=False)
        _poly(p, [(54, 66), (70, 38), (86, 60)], accent, 7, close=False)
        for x, y in ((16, 70), (46, 66), (54, 66), (86, 60)):
            _dot(p, x, y, 6, fg)
    elif name == "batch_neb":
        path = QPainterPath(QPointF(12, 78)); path.cubicTo(QPointF(38, 76), QPointF(40, 22), QPointF(50, 22))
        path.cubicTo(QPointF(60, 22), QPointF(62, 76), QPointF(88, 78))
        _pen(p, muted, 7); p.drawPath(path)
        for x, y in ((22, 77), (36, 52), (50, 22), (64, 52), (78, 77)):
            _dot(p, x, y, 6, accent if x == 50 else fg)
    elif name == "batch_phonons":
        for x in (22, 50, 78):
            for y in (24, 52, 80):
                _dot(p, x, y, 6, accent if (x, y) == (50, 52) else fg)
        _poly(p, [(28, 52), (34, 44), (40, 60), (46, 44), (50, 52)], muted, 5, close=False)
    elif name == "batch_elastic":
        _pen(p, muted, 6); p.drawRect(QRectF(30, 26, 40, 48))
        _pen(p, accent, 7); p.drawRect(QRectF(18, 34, 64, 32))
        _poly(p, [(6, 50), (14, 42), (14, 58)], accent, 5); _poly(p, [(94, 50), (86, 42), (86, 58)], accent, 5)
    elif name == "batch_ts":
        path = QPainterPath(QPointF(10, 76)); path.quadTo(QPointF(30, 24), QPointF(50, 44))
        path.quadTo(QPointF(70, 64), QPointF(90, 20))
        _pen(p, muted, 7); p.drawPath(path)
        _dot(p, 50, 44, 10, accent)
    elif name == "template_load":
        _pen(p, fg, 7); p.drawRoundedRect(QRectF(18, 12, 54, 70), 6, 6)
        for y in (30, 44, 58):
            _line(p, 30, y, 60, y, muted, 5)
        _line(p, 70, 64, 70, 92, accent, 8); _poly(p, [(58, 80), (70, 92), (82, 80)], accent, 7, close=False)
    elif name == "template_save":
        _pen(p, fg, 7); p.drawRoundedRect(QRectF(18, 12, 54, 70), 6, 6)
        for y in (30, 44, 58):
            _line(p, 30, y, 60, y, muted, 5)
        _line(p, 70, 64, 70, 92, accent, 8); _poly(p, [(58, 76), (70, 64), (82, 76)], accent, 7, close=False)
    elif name == "help":
        _pen(p, fg, 8); p.drawEllipse(QPointF(50, 50), 38, 38)
        f = QFont(); f.setPixelSize(56); f.setBold(True); p.setFont(f)
        p.setPen(QColor(accent)); p.drawText(QRectF(0, 0, 100, 100), Qt.AlignmentFlag.AlignCenter, "?")
    elif name == "about":
        _pen(p, fg, 8); p.drawEllipse(QPointF(50, 50), 38, 38)
        _dot(p, 50, 30, 6, accent); _line(p, 50, 46, 50, 72, accent, 10)
    elif name == "language":
        _pen(p, fg, 7); p.drawEllipse(QPointF(50, 50), 38, 38)
        p.drawEllipse(QPointF(50, 50), 16, 38)
        _line(p, 12, 50, 88, 50, muted, 6); _line(p, 18, 30, 82, 30, muted, 5); _line(p, 18, 70, 82, 70, muted, 5)
    elif name == "theme_light":
        _dot(p, 50, 50, 18, accent)
        for k in range(8):
            a = math.radians(k * 45)
            _line(p, 50 + 28 * math.cos(a), 50 + 28 * math.sin(a), 50 + 40 * math.cos(a), 50 + 40 * math.sin(a), fg, 7)
    elif name == "theme_dark":
        path = QPainterPath(); path.addEllipse(QPointF(50, 50), 36, 36)
        cut = QPainterPath(); cut.addEllipse(QPointF(66, 38), 30, 30)
        _fill(p, fg); p.drawPath(path.subtracted(cut))
    elif name == "theme_auto":
        _pen(p, fg, 7); p.drawEllipse(QPointF(50, 50), 36, 36)
        path = QPainterPath(QPointF(50, 14)); path.arcTo(QRectF(14, 14, 72, 72), 90, 180); path.closeSubpath()
        _fill(p, fg); p.drawPath(path)
    elif name in ("frame_auto", "frame_custom", "frame_native"):
        _pen(p, fg, 7); p.drawRoundedRect(QRectF(12, 18, 76, 64), 10, 10)
        _line(p, 12, 38, 88, 38, muted, 6)
        if name == "frame_custom":
            for k, c in enumerate(("#FEBC2E", "#28C840", "#FF5F57")):
                _dot(p, 56 + k * 11, 28, 4.5, c)
        elif name == "frame_native":
            _line(p, 62, 28, 68, 28, fg, 4); _pen(p, fg, 4); p.drawRect(QRectF(72, 25, 6, 6))
        else:
            f = QFont(); f.setPixelSize(34); f.setBold(True); p.setFont(f)
            p.setPen(QColor(accent)); p.drawText(QRectF(12, 40, 76, 42), Qt.AlignmentFlag.AlignCenter, "A")
    elif name == "show_run":
        path = QPainterPath(QPointF(10, 50)); path.quadTo(QPointF(50, 8), QPointF(90, 50)); path.quadTo(QPointF(50, 92), QPointF(10, 50))
        _pen(p, fg, 7); p.drawPath(path)
        _dot(p, 50, 50, 14, accent)
    elif name == "search":
        _pen(p, fg, 8); p.drawEllipse(QPointF(42, 42), 26, 26)
        _line(p, 62, 62, 84, 84, accent, 10)
    elif name == "keyboard":
        _pen(p, fg, 7); p.drawRoundedRect(QRectF(10, 28, 80, 44), 8, 8)
        for x in (26, 42, 58, 74):
            _dot(p, x, 42, 4, fg)
        _line(p, 30, 58, 70, 58, accent, 7)
    elif name in ("chevron_up", "chevron_down"):
        y0, y1 = (62, 38) if name == "chevron_up" else (38, 62)
        _poly(p, [(26, y0), (50, y1), (74, y0)], muted, 10, close=False)

    elif name == "tab_structure":
        _draw("preset", p, t)
    elif name == "tab_files":
        _draw("file", p, t)
    elif name == "tab_run":
        _draw("run", p, t)
    elif name == "tab_workspace":
        _pen(p, muted, 7); p.setBrush(Qt.BrushStyle.NoBrush); p.drawRoundedRect(QRectF(16, 24, 68, 52), 6, 6)
        _poly(p, [(32, 44), (44, 54), (32, 64)], accent, 7, close=False)
        _line(p, 52, 64, 70, 64, accent, 7)
    elif name == "tab_analysis":
        _fill(p, muted); p.drawRoundedRect(QRectF(18, 54, 16, 32), 4, 4)
        _fill(p, accent); p.drawRoundedRect(QRectF(42, 30, 16, 56), 4, 4)
        _fill(p, fg); p.drawRoundedRect(QRectF(66, 44, 16, 42), 4, 4)
    else:
        _dot(p, 50, 50, 22, muted)


def icon(name: str, size: int = 18) -> QIcon:
    key = (name, size, _THEME)
    got = _CACHE.get(key)
    if got is not None:
        return got
    t = _tokens()
    pix = QPixmap(size, size)
    pix.fill(Qt.GlobalColor.transparent)
    p = QPainter(pix)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setRenderHint(QPainter.RenderHint.TextAntialiasing)
    p.scale(size / _UNIT, size / _UNIT)
    try:
        _draw(name, p, t)
    finally:
        p.end()
    got = QIcon(pix)
    _CACHE[key] = got
    return got


def code_icon(code: str, size: int = 18) -> QIcon:
    return icon(f"code_{code}", size)


def source_icon(source: str, size: int = 18) -> QIcon:
    return icon({"surface": "slab", "fetch": "file"}.get(source, source), size)


def task_icon(task: str, size: int = 18) -> QIcon:
    return icon(f"task_{task}", size)


def deg(x: float) -> int:
    return int(x * 16)


__all__ = ["icon", "code_icon", "source_icon", "task_icon", "set_theme", "math"]
