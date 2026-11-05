
"""Spatial quality control module for spatial omics data."""

import importlib
import importlib.metadata
import argparse
import sys
from pathlib import Path
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
    capture_crop = sdata.query.bounding_box(
        axes=("x", "y"),
        min_coordinate=[minx - padding, miny - padding],
        max_coordinate=[maxx + padding, maxy + padding],
        target_coordinate_system=sample_id,
    )

    has_cytassist = cytassist_key in capture_crop.images
    n_panels = 2 if has_cytassist else 1
    fig, axes = plt.subplots(1, n_panels, figsize=(12, 5.5), squeeze=False)
    axes = axes.ravel()

    hires_gray_key = "_registration_hires_gray"
    add_grayscale_preview_image(capture_crop, hires_key, hires_gray_key)

    if has_cytassist:
        cytassist_gray_key = "_registration_cytassist_gray"
        add_grayscale_preview_image(capture_crop, cytassist_key, cytassist_gray_key)
        cytassist_cmap = LinearSegmentedColormap.from_list(
            "cytassist_magenta",
            ["#000000", "#ff00ff"],
        )
        capture_crop.pl.render_images(
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
            title="H&E / CytAssist registration",
        )
        axes[0].set_axis_off()
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
        title="Space Ranger in-tissue bins",
    )
    tissue_axis.set_axis_off()

    fig.subplots_adjust(left=0.01, right=0.99, top=0.94, bottom=0.01, wspace=0.05)
    fig.savefig(output_path, dpi=100, bbox_inches="tight", pad_inches=0.03)
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
        f"{sample_id}_registration_tissue_mqc.png",
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

    # Crop the spatialdata object to the visium area
    shape_key = f"{sample_id}_{resolution}"
    minx, miny, maxx, maxy = sdata[shape_key].total_bounds
    sdata_crop = sdata.query.bounding_box(
        axes=("x", "y"),
        min_coordinate=[minx-1000, miny-1000],
        max_coordinate=[maxx+1000, maxy+1000],
        target_coordinate_system=sample_id,
    )

    # Spatial QC plot
    with plt.rc_context({
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
    }):
        fig, axs = plt.subplots(3, 2, figsize=(11.8, 15.5), constrained_layout=False)
        cmap = "viridis"

        for row, (label, metric, _, colorbar_label) in enumerate(spatial_metrics):
            unfiltered_values = pd.to_numeric(
                adata.obs[preview_columns[label][0]], errors="coerce"
            ).dropna().to_numpy(dtype=float)
            if unfiltered_values.size > 0:
                vmin = np.percentile(unfiltered_values, 1)
                vmax = np.percentile(unfiltered_values, 99)
                if vmin == vmax:
                    vmax = vmin + 1
            else:
                vmin, vmax = 0, 1
            norm = Normalize(vmin=vmin, vmax=vmax)

            for col, (column_title, color) in enumerate(zip(
                ("In tissue · before QC", "QC pass"), preview_columns[label]
            )):
                axis = axs[row, col]
                sdata_crop.pl.render_images( # pylint: disable=no-member
                    f"{sample_id}_hires_image",
                    cmap="gray",
                ).pl.render_shapes(
                    f"{sample_id}_{resolution}",
                    color=color,
                    cmap=cmap,
                    norm=norm,
                    fill_alpha=0.5,
                    outline_width=0,
                    na_color="#00000000",
                    colorbar=False,
                ).pl.show(
                    coordinate_systems=sample_id,
                    dpi=100,
                    ax=axis,
                    colorbar=False,
                )
                axis.set_title(f"{label} · {column_title}", fontsize=8)
                axis.set_axis_off()

            colorbar = fig.colorbar(
                ScalarMappable(norm=norm, cmap=cmap),
                ax=axs[row, :],
                fraction=0.035,
                pad=0.025,
            )
            colorbar.set_label(colorbar_label, fontsize=7)
            colorbar.ax.tick_params(labelsize=7)

        fig.subplots_adjust(left=0.015, right=0.92, top=0.98, bottom=0.015,
                            wspace=0.06, hspace=0.12)

        fig.savefig(
            f"{sample_id}_qc_spatial_plots.png",
            dpi=100,
            bbox_inches="tight",
            pad_inches=0.05,
        )
        plt.close(fig)

    adata.obs.drop(columns=[column for column, _ in preview_columns.values()], inplace=True)

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
