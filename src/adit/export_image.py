
from __future__ import annotations

import math
from xml.sax.saxutils import escape

SVG_HEAD = ('<svg xmlns="http://www.w3.org/2000/svg" width="{w:.0f}" height="{h:.0f}" '
            'viewBox="{x0:.2f} {y0:.2f} {w:.2f} {h:.2f}">\n')
BOND_WIDTH = 2.0
DOUBLE_GAP = 3.5
LABEL_SIZE = 13.0


def _bounds(points, margin: float) -> tuple[float, float, float, float]:
    xs = [p[0] for p in points] or [0.0]
    ys = [p[1] for p in points] or [0.0]
    x0, x1 = min(xs) - margin, max(xs) + margin
    y0, y1 = min(ys) - margin, max(ys) + margin
    return x0, y0, max(x1 - x0, 1.0), max(y1 - y0, 1.0)


def sketch_to_svg(sketch, colors: dict[str, str] | None = None, background: str = "white") -> str:
    from adit.sketch import ELEMENT_COLORS

    colors = {**ELEMENT_COLORS, **(colors or {})}
    atoms, bonds = sketch.atoms, sketch.bonds
    if not atoms:
        raise ValueError("何も描かれていません")
    pts = [(a.x, a.y) for a in atoms]
    x0, y0, w, h = _bounds(pts, margin=28.0)
    out = [SVG_HEAD.format(w=w, h=h, x0=x0, y0=y0)]
    if background:
        out.append(f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{w:.2f}" height="{h:.2f}" fill="{background}"/>\n')
    for b in bonds:
        a, c = atoms[b.a], atoms[b.b]
        dx, dy = c.x - a.x, c.y - a.y
        length = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / length * DOUBLE_GAP, dx / length * DOUBLE_GAP
        offsets = {1: (0.0,), 2: (-1.0, 1.0), 3: (-1.0, 0.0, 1.0)}[b.order]
        for k in offsets:
            out.append(f'<line x1="{a.x + nx * k:.2f}" y1="{a.y + ny * k:.2f}" '
                       f'x2="{c.x + nx * k:.2f}" y2="{c.y + ny * k:.2f}" '
                       f'stroke="#202020" stroke-width="{BOND_WIDTH}" stroke-linecap="round"/>\n')
    for atom in atoms:
        if atom.elem == "C" and atom.charge == 0:
            continue
        text = atom.elem + ("+" * atom.charge if atom.charge > 0 else "−" * -atom.charge)
        out.append(f'<circle cx="{atom.x:.2f}" cy="{atom.y:.2f}" r="{LABEL_SIZE * 0.72:.2f}" fill="{background or "white"}"/>\n')
        out.append(f'<text x="{atom.x:.2f}" y="{atom.y:.2f}" font-family="Helvetica, Arial, sans-serif" '
                   f'font-size="{LABEL_SIZE}" fill="{colors.get(atom.elem, "#202020")}" '
                   f'text-anchor="middle" dominant-baseline="central">{escape(text)}</text>\n')
    out.append("</svg>\n")
    return "".join(out)


def sketch_to_mol(sketch) -> str:
    from rdkit import Chem
    from rdkit.Chem import rdDepictor

    mol = sketch.to_mol()
    rdDepictor.Compute2DCoords(mol)
    return Chem.MolToMolBlock(mol)


def scene_to_svg(scene, width: int = 640, height: int = 520, rx: float = -1.05, ry: float = 0.52,
                 background: str = "white", rotation=None) -> str:
    if scene.n_atoms == 0:
        raise ValueError("構造がありません")
    if rotation is not None:
        # the 3x3 matrix of the desktop viewer (screen = R @ p), so the file shows the view on screen
        R = [[float(x) for x in row] for row in rotation]

        def rotate(p):
            return [R[i][0] * p[0] + R[i][1] * p[1] + R[i][2] * p[2] for i in range(3)]
    else:
        def rotate(p):
            cy, sy = math.cos(ry), math.sin(ry)
            x, z = p[0] * cy + p[2] * sy, -p[0] * sy + p[2] * cy
            cx, sx = math.cos(rx), math.sin(rx)
            return [x, p[1] * cx - z * sx, p[1] * sx + z * cx]

    scale = min(width, height) * 0.42 / max(scene.scale, 1e-6)

    def project(p):
        r = rotate(p)
        return (width / 2 + r[0] * scale, height / 2 - r[1] * scale, r[2])

    out = [SVG_HEAD.format(w=width, h=height, x0=0.0, y0=0.0)]
    if background:
        out.append(f'<rect width="{width}" height="{height}" fill="{background}"/>\n')
    for a, b in scene.cell_lines:
        pa, pb = project(a), project(b)
        out.append(f'<line x1="{pa[0]:.2f}" y1="{pa[1]:.2f}" x2="{pb[0]:.2f}" y2="{pb[1]:.2f}" '
                   'stroke="#9a9a9a" stroke-width="1"/>\n')
    for ia, ib in scene.bonds:
        pa, pb = project(scene.positions[ia]), project(scene.positions[ib])
        mx, my = (pa[0] + pb[0]) / 2, (pa[1] + pb[1]) / 2
        width_b = max(1.5, 0.16 * scale)
        out.append(f'<line x1="{pa[0]:.2f}" y1="{pa[1]:.2f}" x2="{pb[0]:.2f}" y2="{pb[1]:.2f}" '
                   f'stroke="rgba(0,0,0,0.35)" stroke-width="{width_b + 2:.2f}" stroke-linecap="round"/>\n')
        for (start, end, color) in ((pa, (mx, my), scene.colors[ia]), ((mx, my), pb, scene.colors[ib])):
            out.append(f'<line x1="{start[0]:.2f}" y1="{start[1]:.2f}" x2="{end[0]:.2f}" y2="{end[1]:.2f}" '
                       f'stroke="{color}" stroke-width="{width_b:.2f}" stroke-linecap="round"/>\n')
    order = sorted(range(scene.n_atoms), key=lambda k: project(scene.positions[k])[2])
    for k in order:
        x, y, _ = project(scene.positions[k])
        r = max(2.0, scene.radii[k] * 0.55 * scale)
        out.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="{r:.2f}" fill="{scene.colors[k]}" '
                   'stroke="rgba(0,0,0,0.35)" stroke-width="1"/>\n')
    out.append("</svg>\n")
    return "".join(out)
