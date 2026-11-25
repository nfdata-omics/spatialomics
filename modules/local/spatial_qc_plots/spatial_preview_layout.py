"""Shared fixed-pixel layout helpers for spatial MultiQC previews."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

PREVIEW_WIDTH_PX = 1200
PREVIEW_DPI = 100
PREVIEW_MARGIN_PX = 48
PREVIEW_GUTTER_PX = 40
PREVIEW_PANEL_WIDTH_PX = (
    PREVIEW_WIDTH_PX - 2 * PREVIEW_MARGIN_PX - PREVIEW_GUTTER_PX
) / 2
PREVIEW_PANEL_LEFTS_PX = (
    PREVIEW_MARGIN_PX,
    PREVIEW_MARGIN_PX + PREVIEW_PANEL_WIDTH_PX + PREVIEW_GUTTER_PX,
)
PREVIEW_TITLE_STYLE = {
    "fontfamily": "DejaVu Sans",
    "fontsize": 10,
    "fontweight": "semibold",
}


def new_preview_figure(height_px):
    """Create a pixel-exact MultiQC preview figure with a white background."""
    return plt.figure(
        figsize=(PREVIEW_WIDTH_PX / PREVIEW_DPI, height_px / PREVIEW_DPI),
        dpi=PREVIEW_DPI,
        facecolor="white",
    )


def preview_add_axes(fig, left_px, top_px, width_px, height_px):
    """Reserve an Axes rectangle in output pixels, using a top-left origin."""
    height = fig.get_figheight() * fig.dpi
    return fig.add_axes([
        left_px / PREVIEW_WIDTH_PX,
        1 - (top_px + height_px) / height,
        width_px / PREVIEW_WIDTH_PX,
        height_px / height,
    ])


def set_axis_viewport(axis, bounds):
    """Set exact x/y bounds while preserving the plotted axis directions."""
    minx, miny, maxx, maxy = bounds
    axis.set_xlim((maxx, minx) if axis.xaxis_inverted() else (minx, maxx))
    axis.set_ylim((maxy, miny) if axis.yaxis_inverted() else (miny, maxy))


def preview_finish_map_axis(axis, box, viewport=None):
    """Restore a fixed panel box after spatialdata_plot and hide its title."""
    if viewport is not None:
        set_axis_viewport(axis, viewport)
    left_px, top_px, width_px, height_px = box
    fig_height_px = axis.figure.get_figheight() * axis.figure.dpi
    axis.set_position([
        left_px / PREVIEW_WIDTH_PX,
        1 - (top_px + height_px) / fig_height_px,
        width_px / PREVIEW_WIDTH_PX,
        height_px / fig_height_px,
    ])
    axis.set_aspect("equal", adjustable="box", anchor="C")
    axis.set_title("")
    axis.set_axis_off()


def preview_add_title(fig, left_px, title, top_px=35, width_px=PREVIEW_PANEL_WIDTH_PX):
    """Draw a title inside its reserved band, outside the map Axes."""
    height = fig.get_figheight() * fig.dpi
    fig.text(
        (left_px + width_px / 2) / PREVIEW_WIDTH_PX,
        1 - top_px / height,
        title,
        ha="center",
        va="center",
        **PREVIEW_TITLE_STYLE,
    )


def preview_remove_extra_axes(fig, plot_axes):
    """Drop unintended colorbar or legend axes created by a plotting producer."""
    for axis in list(fig.axes):
        if not any(axis is plot_axis for plot_axis in plot_axes):
            axis.remove()


def save_preview_webp(fig, output_path, quality=90):
    """Encode an Agg canvas as RGB WebP without an intermediate PNG file."""
    output_path = Path(output_path)
    if output_path.suffix.lower() != ".webp":
        raise ValueError(f"Expected a .webp MultiQC preview filename: {output_path}")
    fig.canvas.draw()
    image = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())).convert("RGB")
    image.save(output_path, format="WEBP", quality=quality, method=6)
