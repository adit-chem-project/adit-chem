
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QPoint, QSize
from PySide6.QtGui import QGuiApplication, QImage, QPainter
from PySide6.QtWidgets import QFileDialog, QMessageBox, QPushButton, QWidget

from adit.lang import L

COPY_TEXT = ("画像をコピー", "Copy image")
SAVE_TEXT = ("画像を保存…", "Save image…")


def widget_image(widget: QWidget, scale: int = 2) -> QImage:
    size = QSize(max(widget.width(), 1) * scale, max(widget.height(), 1) * scale)
    image = QImage(size, QImage.Format.Format_ARGB32)
    image.setDevicePixelRatio(scale)
    image.fill(0)
    painter = QPainter(image)
    painter.scale(scale, scale)
    widget.render(painter, QPoint(0, 0))
    painter.end()
    return image


def copy_widget(widget: QWidget, scale: int = 2) -> None:
    QGuiApplication.clipboard().setImage(widget_image(widget, scale))


def copy_file(path: Path) -> None:
    QGuiApplication.clipboard().setImage(QImage(str(path)))


def save_dialog(parent: QWidget, default_name: str, svg_text: str | None = None,
                widget: QWidget | None = None, mol_text: str | None = None) -> Path | None:
    filters = []
    if svg_text is not None:
        filters.append(L("SVG (拡大しても粗くならない) (*.svg)", "SVG, stays sharp when enlarged (*.svg)"))
    filters.append(L("PNG 画像 (*.png)", "PNG image (*.png)"))
    if mol_text is not None:
        filters.append(L("MOL 形式 (ChemDraw などで開けます) (*.mol)", "MOL file, opens in ChemDraw (*.mol)"))
    path_text, _ = QFileDialog.getSaveFileName(parent, L(*SAVE_TEXT), default_name, ";;".join(filters))
    if not path_text:
        return None
    path = Path(path_text)
    try:
        if path.suffix.lower() == ".svg" and svg_text is not None:
            path.write_text(svg_text, encoding="utf-8")
        elif path.suffix.lower() == ".mol" and mol_text is not None:
            path.write_text(mol_text, encoding="utf-8")
        else:
            if widget is None:
                raise ValueError(L("この形式では保存できません", "cannot save in this format"))
            target = path if path.suffix else path.with_suffix(".png")
            if not widget_image(widget).save(str(target)):
                raise OSError(L(f"画像を書けません: {target}", f"cannot write the image: {target}"))
            path = target
    except (OSError, ValueError) as ex:
        QMessageBox.warning(parent, L("保存できません", "Cannot save"), str(ex))
        return None
    return path


def copy_button(parent: QWidget | None = None) -> QPushButton:
    b = QPushButton(L(*COPY_TEXT), parent)
    b.setObjectName("link")
    b.setToolTip(L("Word や PowerPoint に貼れる形でクリップボードに入れます",
                   "Puts the image on the clipboard, ready to paste into Word or PowerPoint"))
    return b


def save_button(parent: QWidget | None = None) -> QPushButton:
    b = QPushButton(L(*SAVE_TEXT), parent)
    b.setObjectName("link")
    b.setToolTip(L("SVG (拡大しても粗くならない)、PNG、MOL で保存します",
                   "Saves as SVG (stays sharp), PNG, or MOL"))
    return b
