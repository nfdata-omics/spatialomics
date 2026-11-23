
"""Spatial quality control module for spatial omics data."""

import importlib
import importlib.metadata
import argparse
import sys
from pathlib import Path
from PIL import Image
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import yaml

import numpy as np
import scipy as sp
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import spatialdata
import spatialdata_plot # noqa: F401 # pyright: ignore[reportUnusedImport] # pylint: disable=unused-import
import spotsweeper.local_outliers as lo
from spatialdata.models import Image2DModel
from spatialdata.transformations import get_transformation

matplotlib.use("Agg")

# The preview layout contract is intentionally local to this producer.
# The same pixel grid is used by the other independent Nextflow producers:
# 1200 px canvas, 48 px outer margins, 40 px inter-panel gutter, title
# center at y=35 px, 10 pt DejaVu Sans semibold.  Only MultiQC preview
# figures use this layout; analytical plots remain unchanged.
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
    """Reserve an Axes rectangle in output pixels (origin: top left)."""
    height = fig.get_figheight() * fig.dpi
    return fig.add_axes([
        left_px / PREVIEW_WIDTH_PX,
        1 - (top_px + height_px) / height,
        width_px / PREVIEW_WIDTH_PX,
        height_px / height,
    ])


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
    # spatialdata_plot can add the coordinate system as the Axes title.
    # Titles are drawn separately in a fixed band instead.
    axis.set_title("")
    axis.set_axis_off()


def preview_add_title(fig, left_px, title, top_px=35, width_px=PREVIEW_PANEL_WIDTH_PX):
    """Draw a title inside its own reserved band, outside the map Axes."""
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
    """Drop unintended colorbar/legend axes created by a plotting producer."""
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


def set_axis_viewport(axis, bounds):
    """Set exact x/y bounds while preserving the plotted axis directions."""
    minx, miny, maxx, maxy = bounds
    axis.set_xlim((maxx, minx) if axis.xaxis_inverted() else (minx, maxx))
    axis.set_ylim((maxy, miny) if axis.yaxis_inverted() else (miny, maxy))


def get_element_viewport(sdata, image_key, coordinate_system):
    """Get the full stored image extent in the plotting coordinate system."""
    extent = spatialdata.get_extent(
        sdata.images[image_key], coordinate_system=coordinate_system
    )
    return (
        float(extent["x"][0]),
        float(extent["y"][0]),
        float(extent["x"][1]),
        float(extent["y"][1]),
    )


def add_grayscale_preview_image(sdata, source_key, target_key):
    """Add a temporary single-channel luminance view with source transformations."""
    image = sdata.images[source_key]
    channel_names = set(image.coords["c"].values.tolist())
    if {"r", "g", "b"}.issubset(channel_names):
        grayscale = (
            0.2126 * image.sel(c="r")
            + 0.7152 * image.sel(c="g")
            + 0.0722 * image.sel(c="b")
        )
    else:
        grayscale = image.mean(dim="c")
    grayscale = grayscale.expand_dims(c=["gray"]).transpose("c", "y", "x")
    grayscale.attrs = {}
    sdata.images[target_key] = Image2DModel.parse(
        grayscale,
        transformations=get_transformation(image, get_all=True),
    )


def save_registration_tissue_preview(sdata, sample_id, resolution, output_path):
    """Save registered H&E/CytAssist and Space Ranger tissue-selection previews."""
    hires_key = f"{sample_id}_hires_image"
    cytassist_key = f"{sample_id}_cytassist_image"
    square_key = f"{sample_id}_{resolution}"
    capture_key = f"{sample_id}_capture_area"

    required = [hires_key, square_key, capture_key]
    missing = [key for key in required if key not in sdata]
    if missing:
        raise ValueError(
            "Cannot create the registration and tissue-selection preview; missing "
            f"SpatialData element(s): {', '.join(missing)}"
        )

    minx, miny, maxx, maxy = sdata.shapes[capture_key].total_bounds
    padding = max(maxx - minx, maxy - miny) * 0.02
    tissue_bounds = (minx - padding, miny - padding, maxx + padding, maxy + padding)

    # Use stored metadata for the H&E field of view, not the first image
    # artist produced by spatialdata_plot. This also makes it possible to
    # diagnose whether the imported hires and full H&E cover different areas.
    hires_viewport = get_element_viewport(sdata, hires_key, sample_id)
    print(f"Registration H&E source: {hires_key}; extent: {hires_viewport}")
    full_key = f"{sample_id}_full_image"
    if full_key in sdata.images:
        full_viewport = get_element_viewport(sdata, full_key, sample_id)
        print(f"Registration full H&E available: {full_key}; extent: {full_viewport}")
        if not np.allclose(hires_viewport, full_viewport, rtol=1e-3, atol=1e-3):
            print(
                "WARNING: hires H&E and full H&E have different spatial extents. "
                "The registration preview uses the hires element."
            )
    else:
        print("Registration: no *_full_image element found in SpatialData")

    hires_gray_key = "_registration_hires_gray"
    add_grayscale_preview_image(sdata, hires_key, hires_gray_key)

    has_cytassist = cytassist_key in sdata.images
    if has_cytassist:
        cytassist_gray_key = "_registration_cytassist_gray"
        add_grayscale_preview_image(sdata, cytassist_key, cytassist_gray_key)

    capture_crop = sdata.query.bounding_box(
        axes=("x", "y"),
        min_coordinate=list(tissue_bounds[:2]),
        max_coordinate=list(tissue_bounds[2:]),
        target_coordinate_system=sample_id,
    )

    # One row, using the same two-column boxes as other MultiQC producers.
    # With no CytAssist, the single panel spans the full two-column grid.
    canvas_height_px = 570
    plot_top_px = 62
    plot_height_px = 455
    panel_boxes = (
        [(left, plot_top_px, PREVIEW_PANEL_WIDTH_PX, plot_height_px)
         for left in PREVIEW_PANEL_LEFTS_PX]
        if has_cytassist else
        [(PREVIEW_MARGIN_PX, plot_top_px,
          PREVIEW_WIDTH_PX - 2 * PREVIEW_MARGIN_PX, plot_height_px)]
    )
    fig = new_preview_figure(canvas_height_px)
    axes = [preview_add_axes(fig, *box) for box in panel_boxes]
    try:
        if has_cytassist:
            cytassist_cmap = LinearSegmentedColormap.from_list(
                "cytassist_magenta",
                ["#000000", "#ff00ff"],
            )
            sdata.pl.render_images(
                hires_gray_key,
                cmap="gray",
                colorbar=False,
            ).pl.render_images(
                cytassist_gray_key,
                cmap=cytassist_cmap,
                alpha=0.45,
                colorbar=False,
            ).pl.render_shapes(
                capture_key,
                fill_alpha=0,
                outline_color="#39ff14",
                outline_width=2,
                outline_alpha=1,
            ).pl.show(
                coordinate_systems=sample_id,
                ax=axes[0],
                colorbar=False,
            )
            if axes[0].images:
                print(
                    "Registration first Matplotlib image artist extent: "
                    f"{axes[0].images[0].get_extent()}"
                )
            preview_finish_map_axis(
                axes[0], panel_boxes[0], viewport=hires_viewport
            )
            preview_add_title(
                fig, PREVIEW_PANEL_LEFTS_PX[0], "H&E / CytAssist registration"
            )
            tissue_axis = axes[1]
        else:
            print(
                f"CytAssist layer '{cytassist_key}' is unavailable; "
                "creating only the Space Ranger tissue-selection panel"
            )
            tissue_axis = axes[0]

        capture_crop.pl.render_images(
            hires_gray_key,
            cmap="gray",
            colorbar=False,
        ).pl.render_shapes(
            square_key,
            color="#00d5ff",
            fill_alpha=0.35,
            outline_width=0,
        ).pl.render_shapes(
            capture_key,
            fill_alpha=0,
            outline_color="#ffffff",
            outline_width=1.5,
            outline_alpha=0.9,
        ).pl.show(
            coordinate_systems=sample_id,
            ax=tissue_axis,
            colorbar=False,
        )
        set_axis_viewport(tissue_axis, tissue_bounds)
        preview_finish_map_axis(
            tissue_axis, panel_boxes[1] if has_cytassist else panel_boxes[0]
        )
        if has_cytassist:
            preview_add_title(
                fig, PREVIEW_PANEL_LEFTS_PX[1], "Space Ranger in-tissue bins"
            )
        else:
            preview_add_title(
                fig,
                PREVIEW_MARGIN_PX,
                "Space Ranger in-tissue bins",
                width_px=PREVIEW_WIDTH_PX - 2 * PREVIEW_MARGIN_PX,
            )

        preview_remove_extra_axes(fig, axes)
        save_preview_webp(fig, output_path, quality=90)
    finally:
        plt.close(fig)


def spatial_quality_control(
    zarr_folder,
    output_zarr,
    sample_id,
    resolution="square_016um",
    min_counts=100,
    min_genes=50,
    max_mt=20,
    novelty_thresh=0.0
):
    """
    Perform spatial quality control analysis on spatial omics data.

    Reads a spatialdata zarr object, calculates QC metrics, identifies global and
    local outliers, generates distribution plots and spatial visualizations, and
    saves the QC annotations in a new SpatialData Zarr and writes summary outputs.

    Parameters
    ----------
    zarr_folder : str
        Path to the input Zarr object containing spatialdata.
    output_zarr : str
        Path to the output SpatialData Zarr containing the QC annotations.
    sample_id : str
        Unique identifier for the sample being processed.
    resolution : str, optional
        Label of the table in the Zarr object pointing to the resolution to process. Default is "square_016um".
    min_counts : int, optional
        Minimum number of UMI counts per spot. Default is 100.
    min_genes : int, optional
        Minimum number of genes detected per spot. Default is 50.
    max_mt : float, optional
        Maximum percentage of mitochondrial counts allowed. Default is 20.
    novelty_thresh : float, optional
        Minimum novelty score (complexity) threshold. Default is 0.0.

    Returns
    -------
    None
    """

    sns.set_theme(style="white")
    print(f"\n=== Processing sample: {sample_id} ===\n")

    # Load AnnData at the specified resolution
    zarr_folder = Path(zarr_folder)
    output_zarr = Path(output_zarr)
    if zarr_folder.resolve() == output_zarr.resolve():
        raise ValueError("The input and output SpatialData Zarr paths must be different")

    sdata = spatialdata.read_zarr(zarr_folder)
    if resolution not in sdata.tables:
        raise KeyError(f"Table '{resolution}' was not found in the input SpatialData")
    adata = sdata.tables[resolution]
    adata.obs['sample'] = sample_id

    print(sdata)

    # Ensure the expression matrix is in sparse format
    if not sp.sparse.issparse(adata.X):
        adata.X = sp.sparse.csr_matrix(adata.X)

    # Flag genes for QC
    adata.var["gene_symbols_upper"] = adata.var_names.str.upper()

    adata.var["mt"] = adata.var["gene_symbols_upper"].str.startswith("MT-")
    adata.var["ribo"] = adata.var["gene_symbols_upper"].str.startswith(('RPS','RPL'))
    adata.var["hb"] = adata.var["gene_symbols_upper"].str.startswith("HB")

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=['mt','ribo','hb'],
        inplace=True,
        percent_top=None,
        log1p=False
    )

    # Generate distribution plots for QC metrics
    distribution_plots(adata, sample_id, min_counts, min_genes, max_mt)

    # Compute novelty (complexity) score
    adata.obs["novelty_score"] = adata.obs["n_genes_by_counts"] / adata.obs["total_counts"]

    # QC flags for global outliers
    adata.obs["qc_low_complexity"] = adata.obs["novelty_score"] <= novelty_thresh
    adata.obs["qc_lib_size"] = adata.obs["total_counts"] < min_counts
    adata.obs["qc_detected"] = adata.obs["n_genes_by_counts"] < min_genes
    adata.obs["qc_mito"] = adata.obs["pct_counts_mt"] > max_mt

    # Combine global outliers
    adata.obs["global_outliers"] = (
        adata.obs["qc_lib_size"] |
        adata.obs["qc_detected"] |
        adata.obs["qc_mito"] |
        adata.obs["qc_low_complexity"]
    )

    # Local outlier detection
    for metric in ["total_counts", "n_genes_by_counts", "pct_counts_mt"]:
        lo.local_outliers(adata, metric=metric, sample_key="region", n_neighbors=36)

    # Combine local outliers
    adata.obs["local_outliers"] = (
        adata.obs["total_counts_outliers"] |
        adata.obs["n_genes_by_counts_outliers"] |
        adata .obs["pct_counts_mt_outliers"]
    )

    qc_flag_columns = [
        "qc_low_complexity",
        "qc_lib_size",
        "qc_detected",
        "qc_mito",
        "global_outliers",
        "total_counts_outliers",
        "n_genes_by_counts_outliers",
        "pct_counts_mt_outliers",
        "local_outliers",
    ]

    # Mask all the filters for in-tissue spots only
    adata.obs[qc_flag_columns] = adata.obs[qc_flag_columns].astype("boolean")
    out_of_tissue_mask = adata.obs["in_tissue"] == 0
    adata.obs.loc[out_of_tissue_mask, qc_flag_columns] = pd.NA
    in_tissue_mask = adata.obs["in_tissue"] == 1
    global_outliers = adata.obs["global_outliers"].fillna(False).astype(bool)
    valid_bins = in_tissue_mask & ~global_outliers
    adata.obs["qc_pass"] = valid_bins.to_numpy(dtype=bool)
    valid_obs = adata.obs.loc[valid_bins]
    n_in_tissue = int(in_tissue_mask.sum())
    n_valid_bins = int(valid_bins.sum())
    valid_bins_fraction = n_valid_bins / n_in_tissue if n_in_tissue else float("nan")

    # save annotation to csv
    adata.obs.to_csv(f"{sample_id}_qc_annotated_obs.csv")

    # Record the QC configuration alongside the annotated table.
    adata.uns["spatialomics_qc"] = {
        "resolution": resolution,
        "min_counts": int(min_counts),
        "min_genes": int(min_genes),
        "max_mt": float(max_mt),
        "novelty_thresh": float(novelty_thresh),
        "pass_column": "qc_pass",
    }

    # Persist the annotated table as part of a new canonical SpatialData object.
    sdata.tables[resolution] = adata
    sdata.write(output_zarr)

    save_registration_tissue_preview(
        sdata,
        sample_id,
        resolution,
        f"{sample_id}_registration_tissue_mqc.webp",
    )

    # Collect QC summary
    summary_dict = {
        "Sample": sample_id,
        "Bin size": resolution,
        "Total spots": adata.n_obs,
        "In-tissue spots": n_in_tissue,
        "Valid bins": n_valid_bins,
        "Valid bins fraction": valid_bins_fraction,
        "Mean UMI counts in valid bins": valid_obs["total_counts"].mean(),
        "Mean genes in valid bins": valid_obs["n_genes_by_counts"].mean(),
        "Mean mitochondrial % in valid bins": valid_obs["pct_counts_mt"].mean(),
        "Mean novelty score in valid bins": valid_obs["novelty_score"].mean(),
        "Low number of UMI per spot": int(adata.obs["qc_lib_size"].sum()),
        "Low number of genes per spot": int(adata.obs["qc_detected"].sum()),
        "High % mitochondrial counts": int(adata.obs["qc_mito"].sum()),
        "Low complexity": int(adata.obs["qc_low_complexity"].sum()),
        "Global outliers": int(adata.obs["global_outliers"].sum()),
        "Local outliers": int(adata.obs["local_outliers"].sum())
    }

    # Save summary as CSV
    summary_df = pd.DataFrame([summary_dict])
    summary_df.to_csv(f"{sample_id}_qc_metrics.csv", index=False)

    adata.obs[qc_flag_columns] = adata.obs[qc_flag_columns].astype(str)

    # Create new annotation column for filtered QC metrics
    adata.obs["total_counts_filtered"] = adata.obs["total_counts"].where(valid_bins, pd.NA)
    adata.obs["n_genes_by_counts_filtered"] = adata.obs["n_genes_by_counts"].where(valid_bins, pd.NA)
    adata.obs["pct_counts_mt_filtered"] = adata.obs["pct_counts_mt"].where(valid_bins, pd.NA)

    spatial_metrics = [
        ("Total counts", "total_counts", "total_counts_filtered", "UMI counts"),
        (
            "Detected genes",
            "n_genes_by_counts",
            "n_genes_by_counts_filtered",
            "Genes per bin",
        ),
        (
            "Mitochondrial fraction",
            "pct_counts_mt",
            "pct_counts_mt_filtered",
            "Mitochondrial counts (%)",
        ),
    ]
    preview_columns = {}
    for label, metric, filtered_metric, _ in spatial_metrics:
        preview_column = f"_spatial_qc_preview_{metric}"
        adata.obs[preview_column] = adata.obs[metric].where(in_tissue_mask)
        preview_columns[label] = (preview_column, filtered_metric)

    # Use one capture-area viewport for all six panels. The absolute padding is
    # derived once in the plotted coordinate system and is not data-dependent.
    shape_key = f"{sample_id}_{resolution}"
    capture_key = f"{sample_id}_capture_area"
    if capture_key not in sdata.shapes:
        raise KeyError(f"Capture-area shape '{capture_key}' was not found")
    minx, miny, maxx, maxy = sdata.shapes[capture_key].total_bounds
    spatial_padding = max(maxx - minx, maxy - miny) * 0.02
    spatial_bounds = (
        minx - spatial_padding,
        miny - spatial_padding,
        maxx + spatial_padding,
        maxy + spatial_padding,
    )
    sdata_crop = sdata.query.bounding_box(
        axes=("x", "y"),
        min_coordinate=list(spatial_bounds[:2]),
        max_coordinate=list(spatial_bounds[2:]),
        target_coordinate_system=sample_id,
    )

    # MultiQC spatial QC preview: three deterministic two-panel rows.
    # The shared colorbar is deliberately narrower than a panel (420 px)
    # and centered below the pair, rather than spanning the entire canvas.
    save_spatial_qc_preview(
        sdata_crop=sdata_crop,
        adata=adata,
        sample_id=sample_id,
        resolution=resolution,
        spatial_bounds=spatial_bounds,
        spatial_metrics=spatial_metrics,
        preview_columns=preview_columns,
        output_path=f"{sample_id}_qc_spatial_plots.webp",
    )

    adata.obs.drop(columns=[column for column, _ in preview_columns.values()], inplace=True)

def save_spatial_qc_preview(
    sdata_crop,
    adata,
    sample_id,
    resolution,
    spatial_bounds,
    spatial_metrics,
    preview_columns,
    output_path,
):
    """Render the six QC maps on a fixed grid, with one short bar per row.

    All six maps preserve a common viewport and the same physical scale.
    Per metric, the same Normalize object is used for before/after QC.
    Titles and colorbars have dedicated pixel regions and cannot shift maps.
    """
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
        for row, (label, metric, _, colorbar_label) in enumerate(spatial_metrics):
            values = pd.to_numeric(
                adata.obs[preview_columns[label][0]], errors="coerce"
            ).dropna().to_numpy(dtype=float)
            if values.size:
                vmin = float(np.percentile(values, 1))
                vmax = float(np.percentile(values, 99))
                if vmin == vmax:
                    vmax = vmin + 1
            else:
                vmin, vmax = 0.0, 1.0
            norm = Normalize(vmin=vmin, vmax=vmax)

            row_plot_top = plot_top_px + row * row_pitch_px
            for col, (column_title, column) in enumerate(zip(
                ("In tissue · before QC", "QC pass"), preview_columns[label]
            )):
                left = PREVIEW_PANEL_LEFTS_PX[col]
                box = (left, row_plot_top, PREVIEW_PANEL_WIDTH_PX, plot_height_px)
                axis = preview_add_axes(fig, *box)
                plot_axes.append(axis)
                sdata_crop.pl.render_images(  # pylint: disable=no-member
                    f"{sample_id}_hires_image",
                    cmap="gray",
                    colorbar=False,
                ).pl.render_shapes(
                    f"{sample_id}_{resolution}",
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

            # Remove any incidental spatialdata_plot axes before adding a
            # short, centrally placed scale for both maps in this row.
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
            # The explicit colorbar axes must survive cleanup in later rows.
            plot_axes.append(bar_axis)

        save_preview_webp(fig, output_path, quality=90)
    finally:
        plt.close(fig)


def distribution_plots(adata, sample_id, min_counts, min_genes, max_mt):
    """
    Generate distribution plots for QC metrics.

    Parameters
    ----------
    adata : AnnData
        Annotated data object with QC metrics.
    sample_id : str
        Unique identifier for the sample.
    min_counts : int
        Minimum number of UMI counts per spot.
    min_genes : int
        Minimum number of genes detected per spot.
    max_mt : float
        Maximum percentage of mitochondrial counts allowed.
    """

    fig, axs = plt.subplots(1, 3, figsize=(18, 5))

    sns.histplot(adata.obs["total_counts"], kde=True, color="skyblue", ax=axs[0])
    axs[0].set_title(f"{sample_id} - Total Counts per Spot")
    axs[0].axvline(min_counts, color='red', linestyle='--', label=f'Min counts ({min_counts})')
    axs[0].set_xlabel("Total counts")
    axs[0].set_ylabel("Number of spots")

    sns.histplot(adata.obs["n_genes_by_counts"], kde=True, bins=60, color="lightgreen", ax=axs[1])
    axs[1].set_title(f"{sample_id} - Genes Detected per Spot")
    axs[1].axvline(min_genes, color='red', linestyle='--', label=f'Min genes ({min_genes})')
    axs[1].set_xlabel("Number of genes")
    axs[1].set_ylabel("Number of spots")

    sns.histplot(adata.obs["pct_counts_mt"], kde=True, bins=60, color="salmon", ax=axs[2])
    axs[2].set_title(f"{sample_id} - Percent Mitochondrial Counts")
    axs[2].axvline(max_mt, color='red', linestyle='--', label=f'Max mt ({max_mt})')
    axs[2].set_xlabel("% mt counts")
    axs[2].set_ylabel("Number of spots")

    plt.tight_layout()
    plt.savefig(f"{sample_id}_qc_distributions.png")
    plt.close(fig)


def versions_yaml(process_name, list_of_libs=None):
    """
    Generate YAML formatted string with versions of relevant libraries.

    Parameters
    ----------
    process_name : str
        Process name to use as key in the versions dictionary.
    list_of_libs : list of str, optional
        List of specific library names to include in the versions dictionary.

    Returns
    -------
    str
        YAML formatted string containing library versions and Python version.
    """

    versions = {}
    versions[process_name] = {}

    versions[process_name]['python'] = f"{sys.version_info.major}" \
        f".{sys.version_info.minor}.{sys.version_info.micro}"

    for lib in list_of_libs:
        try:
            version = importlib.metadata.version(lib)
        except importlib.metadata.PackageNotFoundError:
            try:
                module = importlib.import_module(lib)
                version = getattr(module, '__version__', 'unknown')
            except (ImportError, AttributeError):
                version = None
        if version is not None:
            versions[process_name][lib] = version

    return yaml.dump(versions)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Spatial quality control")
    parser.add_argument("--zarr", type=str,
                        help="Path to input Zarr object (zarr file)")
    parser.add_argument("--output-zarr", type=str,
                        help="Path to output SpatialData Zarr with QC annotations")
    parser.add_argument("--sample", type=str,
                        help="Sample ID to process")
    parser.add_argument("--resolution", type=str,
                        help="Label of the table in the Zarr object pointing to the resolution to process")
    parser.add_argument("--min-counts", type=int,
                        help="Minimum counts threshold for spot filtering")
    parser.add_argument("--min-genes", type=int,
                        help="Minimum genes threshold for spot filtering")
    parser.add_argument("--max-mt", type=int,
                        help="Maximum mitochondrial counts percentage threshold for spot filtering")
    parser.add_argument("--versions-dict", type=str,
                        help="Return dictionary of versions used by the module and exit")

    args = parser.parse_args()

    if args.versions_dict:
        libs = [
            "scanpy",
            "spatialdata",
            "spatialdata-plot",
            "spotsweeper",
            "seaborn",
            "matplotlib",
            "pandas",
            "scipy",
        ]
        print(versions_yaml(args.versions_dict, libs))
    else:
        if not args.zarr or not args.output_zarr or not args.sample:
            parser.error("--zarr, --output-zarr and --sample are required")
        spatial_quality_control(
            args.zarr,
            args.output_zarr,
            args.sample,
            resolution = args.resolution,
            min_counts = args.min_counts,
            min_genes = args.min_genes,
            max_mt = args.max_mt
        )
