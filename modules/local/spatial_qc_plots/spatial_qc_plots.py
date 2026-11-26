"""Render per-sample spatial QC reports from the annotated SpatialData Zarr."""

import argparse
import importlib
import importlib.metadata
import sys

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import spatialdata
import spatialdata_plot  # noqa: F401  # pylint: disable=unused-import
import yaml
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

from spatial_preview_layout import (
    PREVIEW_PANEL_LEFTS_PX,
    PREVIEW_PANEL_WIDTH_PX,
    PREVIEW_WIDTH_PX,
    new_preview_figure,
    preview_add_axes,
    preview_add_title,
    preview_finish_map_axis,
    preview_remove_extra_axes,
    save_preview_webp,
)


def save_spatial_qc_preview(sdata, adata, sample_id, resolution, output_path):
    """Render the six QC maps with the original in-tissue and QC-pass sets."""
    in_tissue_mask = adata.obs["in_tissue"] == 1
    qc_pass_mask = adata.obs["qc_pass"].fillna(False).astype(bool)
    spatial_metrics = [
        ("Total counts", "total_counts", "UMI counts"),
        ("Detected genes", "n_genes_by_counts", "Genes per bin"),
        ("Mitochondrial fraction", "pct_counts_mt", "Mitochondrial counts (%)"),
    ]
    preview_columns = {}
    for label, metric, _colorbar_label in spatial_metrics:
        before_column = f"_spatial_qc_preview_{metric}"
        after_column = f"{before_column}_pass"
        adata.obs[before_column] = adata.obs[metric].where(in_tissue_mask)
        adata.obs[after_column] = adata.obs[metric].where(qc_pass_mask)
        preview_columns[label] = (before_column, after_column)

    capture_key = f"{sample_id}_capture_area"
    if capture_key not in sdata.shapes:
        raise KeyError(f"Capture-area shape '{capture_key}' was not found")
    minx, miny, maxx, maxy = sdata.shapes[capture_key].total_bounds
    padding = max(maxx - minx, maxy - miny) * 0.02
    spatial_bounds = (minx - padding, miny - padding, maxx + padding, maxy + padding)
    image_key = f"{sample_id}_hires_image"
    shape_key = f"{sample_id}_{resolution}"
    hires_image = sdata.images.get(image_key)
    qc_elements = spatialdata.SpatialData(
        images={image_key: hires_image} if hires_image is not None else {},
        shapes={shape_key: sdata.shapes[shape_key]},
        tables={resolution: adata},
    )
    sdata_crop = qc_elements.query.bounding_box(
        axes=("x", "y"),
        min_coordinate=list(spatial_bounds[:2]),
        max_coordinate=list(spatial_bounds[2:]),
        target_coordinate_system=sample_id,
    )

    canvas_height_px = 1600
    plot_top_px = 62
    plot_height_px = 405
    row_pitch_px = 520
    bar_width_px = 420
    bar_height_px = 16
    bar_left_px = (PREVIEW_WIDTH_PX - bar_width_px) / 2
    bar_offset_px = plot_top_px + plot_height_px + 18

    fig = new_preview_figure(canvas_height_px)
    plot_axes = []
    try:
        for row, (label, _metric, colorbar_label) in enumerate(spatial_metrics):
            before_column, after_column = preview_columns[label]
            values = pd.to_numeric(
                adata.obs[before_column], errors="coerce"
            ).dropna().to_numpy(dtype=float)
            if values.size:
                vmin = float(np.percentile(values, 1))
                vmax = float(np.percentile(values, 99))
                if vmin == vmax:
                    vmax = vmin + 1
            else:
                vmin, vmax = 0.0, 1.0
            norm = Normalize(vmin=vmin, vmax=vmax, clip=True)
            row_plot_top = plot_top_px + row * row_pitch_px

            for col, (column_title, column) in enumerate((
                ("In tissue · before QC", before_column),
                ("QC pass", after_column),
            )):
                left = PREVIEW_PANEL_LEFTS_PX[col]
                box = (left, row_plot_top, PREVIEW_PANEL_WIDTH_PX, plot_height_px)
                axis = preview_add_axes(fig, *box)
                plot_axes.append(axis)
                layers = sdata_crop
                if image_key in sdata_crop.images:
                    layers = layers.pl.render_images(image_key, cmap="gray", colorbar=False)
                layers.pl.render_shapes(
                    shape_key,
                    color=column,
                    cmap="viridis",
                    norm=norm,
                    fill_alpha=0.5,
                    outline_width=0,
                    na_color="#00000000",
                    colorbar=False,
                ).pl.show(
                    coordinate_systems=sample_id,
                    ax=axis,
                    colorbar=False,
                )
                preview_finish_map_axis(axis, box, viewport=spatial_bounds)
                preview_add_title(
                    fig,
                    left,
                    f"{label} · {column_title}",
                    top_px=35 + row * row_pitch_px,
                )

            preview_remove_extra_axes(fig, plot_axes)
            bar_axis = preview_add_axes(
                fig,
                bar_left_px,
                bar_offset_px + row * row_pitch_px,
                bar_width_px,
                bar_height_px,
            )
            colorbar = fig.colorbar(
                ScalarMappable(norm=norm, cmap="viridis"),
                cax=bar_axis,
                orientation="horizontal",
            )
            colorbar.set_label(colorbar_label, fontsize=9, labelpad=3)
            colorbar.ax.tick_params(labelsize=8, pad=2)
            plot_axes.append(bar_axis)

        save_preview_webp(fig, output_path, quality=90)
    finally:
        plt.close(fig)


def distribution_plots(adata, sample_id, thresholds):
    """Render the three QC distributions using their original full populations."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    sns.histplot(adata.obs["total_counts"], kde=True, color="skyblue", ax=axes[0])
    axes[0].set_title(f"{sample_id} - Total Counts per Spot")
    axes[0].axvline(
        thresholds["min_counts"],
        color="red",
        linestyle="--",
        label=f'Min counts ({thresholds["min_counts"]})',
    )
    axes[0].set_xlabel("Total counts")
    axes[0].set_ylabel("Number of spots")

    sns.histplot(
        adata.obs["n_genes_by_counts"],
        kde=True,
        bins=60,
        color="lightgreen",
        ax=axes[1],
    )
    axes[1].set_title(f"{sample_id} - Genes Detected per Spot")
    axes[1].axvline(
        thresholds["min_genes"],
        color="red",
        linestyle="--",
        label=f'Min genes ({thresholds["min_genes"]})',
    )
    axes[1].set_xlabel("Number of genes")
    axes[1].set_ylabel("Number of spots")

    sns.histplot(
        adata.obs["pct_counts_mt"], kde=True, bins=60, color="salmon", ax=axes[2]
    )
    axes[2].set_title(f"{sample_id} - Percent Mitochondrial Counts")
    axes[2].axvline(
        thresholds["max_mt"],
        color="red",
        linestyle="--",
        label=f'Max mt ({thresholds["max_mt"]})',
    )
    axes[2].set_xlabel("% mt counts")
    axes[2].set_ylabel("Number of spots")

    plt.tight_layout()
    plt.savefig(f"{sample_id}_qc_distributions.png")
    plt.close(fig)


def versions_yaml(process_name):
    """Return package versions used by the report renderer."""
    versions = {process_name: {"python": sys.version.split()[0]}}
    libraries = [
        "matplotlib",
        "numpy",
        "pillow",
        "pandas",
        "seaborn",
        "spatialdata",
        "spatialdata-plot",
    ]
    for library in libraries:
        try:
            version = importlib.metadata.version(library)
        except importlib.metadata.PackageNotFoundError:
            module = importlib.import_module(library.replace("-", "_"))
            version = getattr(module, "__version__", "unknown")
        versions[process_name][library] = version
    return yaml.dump(versions)


def main():
    parser = argparse.ArgumentParser(description="Render per-sample spatial QC reports")
    parser.add_argument("--zarr", help="QC-annotated SpatialData Zarr")
    parser.add_argument("--sample", help="Sample ID")
    parser.add_argument("--resolution", help="SpatialData table key")
    parser.add_argument("--versions-dict", help="Emit software versions and exit")
    args = parser.parse_args()

    if args.versions_dict:
        print(versions_yaml(args.versions_dict))
        return

    if not all((args.zarr, args.sample, args.resolution)):
        parser.error("--zarr, --sample and --resolution are required")

    sns.set_theme(style="white")
    sdata = spatialdata.read_zarr(args.zarr)
    if args.resolution not in sdata.tables:
        raise KeyError(f"Table '{args.resolution}' was not found in the QC SpatialData")
    adata = sdata.tables[args.resolution]
    qc_config = adata.uns.get("spatialomics_qc")
    if not qc_config or qc_config.get("resolution") != args.resolution:
        raise ValueError(
            "The QC table does not contain matching spatialomics_qc configuration"
        )

    distribution_plots(adata, args.sample, qc_config)
    save_spatial_qc_preview(
        sdata,
        adata,
        args.sample,
        args.resolution,
        f"{args.sample}_qc_mqc.webp",
    )


if __name__ == "__main__":
    main()
