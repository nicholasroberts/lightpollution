#!/usr/bin/env python3
"""Shared, Illustrator-friendly vector figure defaults.

Presentation settings only: no changes to datasets, calibrations, statistics,
geometries used for analysis, or existing plotting scripts.

Typical use in a script run as `python src/my_plot.py`::

    from figure_style import apply_figure_style, figure_size, save_figure
    apply_figure_style()
    fig, ax = plt.subplots(figsize=figure_size("double"))
    # ... plot ...
    save_figure(fig, "figures/publication/my_figure", map_figure=True)
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable
import warnings

import matplotlib as mpl


MM_PER_INCH = 25.4
FIGURE_WIDTH_MM = {
    "single": 90.0,
    "double": 180.0,
    "slide": 240.0,
}

# Accessible distinction for lines/points; map colours are data-dependent.
SERIES_COLOURS = {
    "radiance": "#D55E00",
    "population": "#009E73",
    "spider": "#0072B2",
    "moth": "#CC79A7",
}

# Existing project convention: blue = decrease, red = increase.
CHANGE_CMAP = "RdBu_r"

DEFAULT_FONT = "DejaVu Sans"
DEFAULT_VECTOR_FORMATS = ("svg", "pdf")
STANDARD_SIZE_TARGET_MIB = 5.0
MAP_SIZE_TARGET_MIB = 10.0


def apply_figure_style(font: str = DEFAULT_FONT, font_size: float = 9.0) -> None:
    """Set consistent Matplotlib defaults for editable vector output.

    Install the same font on the machine running Adobe Illustrator to avoid
    substitution. Matplotlib does not make fonts available to Illustrator.
    """
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font, "Liberation Sans", "Arial", "sans-serif"],
        "font.size": font_size,
        "axes.titlesize": font_size + 1.0,
        "axes.labelsize": font_size,
        "xtick.labelsize": font_size - 1.0,
        "ytick.labelsize": font_size - 1.0,
        "legend.fontsize": font_size - 1.0,
        "figure.titlesize": font_size + 1.0,
        "lines.linewidth": 1.1,
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.unicode_minus": True,
        "text.usetex": False,
        "mathtext.fontset": "dejavusans",
        "svg.fonttype": "none",  # Keep ordinary text as editable SVG text.
        "pdf.fonttype": 42,     # Embed TrueType fonts instead of Type 3.
        "ps.fonttype": 42,
        "pdf.compression": 9,
        "savefig.facecolor": "white",
        "savefig.transparent": False,
    })


def figure_size(preset: str = "double", height_ratio: float = 0.65) -> tuple[float, float]:
    """Return (width, height) in inches from a standard width preset.

    `height_ratio` is the height divided by the width. Pass a custom
    `figsize` directly for unusually tall maps or multi-panel layouts.
    """
    if preset not in FIGURE_WIDTH_MM:
        raise ValueError(f"Unknown preset {preset!r}; choose from {tuple(FIGURE_WIDTH_MM)}")
    if height_ratio <= 0:
        raise ValueError("height_ratio must be positive")
    width = FIGURE_WIDTH_MM[preset] / MM_PER_INCH
    return (width, width * height_ratio)


def save_figure(
    fig,
    output_stem: str | Path,
    *,
    formats: Iterable[str] = DEFAULT_VECTOR_FORMATS,
    map_figure: bool = False,
    size_target_mib: float | None = None,
    pad_inches: float = 0.05,
) -> dict[str, Path]:
    """Export editable vector graphics and warn about oversized files.

    Outputs are never silently rasterised, simplified, or altered. For large
    native-pixel plots, explicitly rasterise *only* the heavy artist in the
    calling script, leaving text/axes/legends as vectors.

    `output_stem` must be a suffix-free filename; the return value maps
    output format to the created file. File-size thresholds are warnings,
    not hard limits that discard usable artwork.
    """
    stem = Path(output_stem)
    if stem.suffix:
        raise ValueError("output_stem must not have a file extension")
    if size_target_mib is None:
        size_target_mib = MAP_SIZE_TARGET_MIB if map_figure else STANDARD_SIZE_TARGET_MIB
    if size_target_mib <= 0:
        raise ValueError("size_target_mib must be positive")

    formats = tuple(str(fmt).lower().lstrip(".") for fmt in formats)
    if not formats or len(set(formats)) != len(formats):
        raise ValueError("Provide one or more distinct output formats")
    unsupported = set(formats) - {"svg", "pdf", "eps"}
    if unsupported:
        raise ValueError(f"Unsupported vector format(s): {sorted(unsupported)}")

    stem.parent.mkdir(parents=True, exist_ok=True)
    exported: dict[str, Path] = {}
    for fmt in formats:
        if fmt == "eps":
            warnings.warn(
                "EPS does not support transparency reliably; check any alpha layers "
                "in Illustrator. SVG is the preferred master.",
                UserWarning,
                stacklevel=2,
            )
        out = stem.with_suffix(f".{fmt}")
        fig.savefig(
            out,
            format=fmt,
            bbox_inches="tight",
            pad_inches=pad_inches,
            facecolor="white",
            transparent=False,
        )
        exported[fmt] = out
        size_mib = out.stat().st_size / (1024 * 1024)
        if size_mib > size_target_mib:
            warnings.warn(
                f"{out} is {size_mib:.1f} MiB; target is {size_target_mib:.1f} MiB. "
                "Inspect repeated artists and map path complexity; do not simplify "
                "underlying scientific values or silently rasterise text.",
                UserWarning,
                stacklevel=2,
            )
    return exported
