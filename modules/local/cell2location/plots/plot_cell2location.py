#!/usr/bin/env python3
"""Plot Cell2location results stored in a SpatialData table."""

import argparse
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import spatialdata as sd
import spatialdata_plot  # noqa: F401  # pylint: disable=unused-import
from matplotlib.colors import LogNorm, Normalize, to_hex, to_rgb
from matplotlib.patches import Patch


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot Cell2location results stored in a SpatialData table"
    )
    parser.add_argument("--input-zarr", required=True, type=Path)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--output-prefix", help="Output filename prefix; defaults to sample name")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--table-key", default="square_016um")
    parser.add_argument(
        "--shape-key",
        help="SpatialData shapes element; defaults to <sample-name>_<table-key>",
    )
    parser.add_argument(
        "--image-key",
        help="Optional background image; defaults to <sample-name>_hires_image when present",
    )
    parser.add_argument(
        "--abundance-statistic",
        choices=("mean", "q05"),
        default="mean",
        help="Posterior summary used for biological plots and tables",
    )
    parser.add_argument(
        "--uncertainty-min-abundance",
        type=float,
        default=0.1,
        help="Minimum posterior mean used for relative uncertainty plots",
    )
    parser.add_argument(
        "--multiqc-cell-types",
        nargs="*",
        default=None,
        help="Up to three cell types for the compact MultiQC panel; default: top three",
    )
    return parser.parse_args()


def read_training_history(metadata):
    """Validate and return the training history stored by deconvolution."""
    history = metadata.get("training_history")
    if isinstance(history, dict):
        epochs = np.asarray(history.get("epoch", []))
        elbo = np.asarray(history.get("elbo_train", []), dtype=float)
    elif history is not None:
        # Backward compatibility with SpatialData produced before epochs were
        # stored explicitly alongside the ELBO values.
        elbo = np.asarray(history, dtype=float)
        epochs = np.arange(len(elbo))
    else:
        raise ValueError("table.uns['cell2location']['training_history'] is missing")
    if epochs.ndim != 1 or elbo.ndim != 1 or len(epochs) != len(elbo):
        raise ValueError("Cell2location training-history arrays are not aligned")
    if len(epochs) == 0:
        raise ValueError("Cell2location training history is empty")
    if not np.isfinite(elbo).all():
        raise ValueError("Cell2location training history contains non-finite ELBO values")
    return pd.DataFrame({"epoch": epochs, "elbo_train": elbo})


def plot_training_history(training_history, output_path):
    """Plot complete and post-warmup Cell2location training histories."""
    history_start = max(1, int(len(training_history) * 0.2))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(training_history["epoch"], training_history["elbo_train"])
    axes[0].set(title="Complete training history", xlabel="Epoch", ylabel="-ELBO loss")
    axes[1].plot(
        training_history["epoch"].iloc[history_start:],
        training_history["elbo_train"].iloc[history_start:],
    )
    axes[1].set(title="Training history after first 20%", xlabel="Epoch", ylabel="-ELBO loss")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_reconstruction_accuracy(reconstruction_qc, output_path):
    """Render the expected-vs-observed histogram computed during deconvolution."""
    if reconstruction_qc.get("status") != "success":
        return False

    histogram = np.asarray(reconstruction_qc.get("histogram"), dtype=float)
    observed_edges = np.asarray(
        reconstruction_qc.get("observed_log10_edges"), dtype=float
    )
    expected_edges = np.asarray(
        reconstruction_qc.get("expected_log10_edges"), dtype=float
    )
    if histogram.ndim != 2:
        raise ValueError("Reconstruction histogram is not two-dimensional")
    if len(observed_edges) != histogram.shape[0] + 1:
        raise ValueError("Observed reconstruction edges do not match the histogram")
    if len(expected_edges) != histogram.shape[1] + 1:
        raise ValueError("Expected reconstruction edges do not match the histogram")
    if not np.isfinite(histogram).all() or (histogram < 0).any():
        raise ValueError("Reconstruction histogram contains invalid counts")
    if not np.isfinite(observed_edges).all() or not np.isfinite(expected_edges).all():
        raise ValueError("Reconstruction histogram edges are not finite")
    positive_counts = histogram[histogram > 0]
    if positive_counts.size == 0:
        raise ValueError("Reconstruction histogram contains no observations")

    fig, axis = plt.subplots(figsize=(6, 5.5))
    mesh = axis.pcolormesh(
        observed_edges,
        expected_edges,
        np.ma.masked_less_equal(histogram.T, 0),
        cmap="viridis",
        norm=LogNorm(vmin=1, vmax=max(2, float(positive_counts.max()))),
        shading="auto",
    )
    lower = min(float(observed_edges[0]), float(expected_edges[0]))
    upper = max(float(observed_edges[-1]), float(expected_edges[-1]))
    axis.plot(
        [lower, upper],
        [lower, upper],
        color="black",
        linestyle="--",
        linewidth=0.8,
    )
    axis.set(
        title="Reconstruction accuracy",
        xlabel="Observed counts, log10(count + 1)",
        ylabel="Posterior expected counts, log10(count + 1)",
        xlim=(lower, upper),
        ylim=(lower, upper),
    )
    axis.set_aspect("equal", adjustable="box")
    fig.colorbar(mesh, ax=axis, label="Gene-location pairs per bin")
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    return True


def plot_posterior_mean_vs_q05(mean, q05, sample_name, output_path):
    """Plot faceted hexbins comparing posterior mean and q05."""
    cell_types = mean.columns.tolist()
    ncols = min(4, len(cell_types))
    nrows = math.ceil(len(cell_types) / ncols)
    gridsize = 40

    maximum = np.log1p(
        max(
            float(np.nanmax(mean.to_numpy())),
            float(np.nanmax(q05.to_numpy())),
        )
    )
    maximum = max(maximum, 1.0)

    max_bin_count = 2
    for cell_type in cell_types:
        x = np.log1p(mean[cell_type].dropna().to_numpy())
        y = np.log1p(q05.loc[mean[cell_type].dropna().index, cell_type].to_numpy())
        counts, _, _ = np.histogram2d(
            x,
            y,
            bins=gridsize,
            range=((0, maximum), (0, maximum)),
        )
        max_bin_count = max(max_bin_count, int(counts.max()))

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(3.4 * ncols, 3.2 * nrows),
        squeeze=False,
    )
    last_hexbin = None
    for axis, cell_type in zip(axes.ravel(), cell_types):
        valid = mean[cell_type].notna() & q05[cell_type].notna()
        x = np.log1p(mean.loc[valid, cell_type].to_numpy())
        y = np.log1p(q05.loc[valid, cell_type].to_numpy())
        last_hexbin = axis.hexbin(
            x,
            y,
            gridsize=gridsize,
            extent=(0, maximum, 0, maximum),
            mincnt=1,
            cmap="viridis",
            norm=LogNorm(vmin=1, vmax=max_bin_count),
        )
        axis.plot([0, maximum], [0, maximum], color="black", linestyle="--", linewidth=0.8)
        axis.set(
            title=cell_type,
            xlim=(0, maximum),
            ylim=(0, maximum),
            xlabel="log1p(posterior mean)",
            ylabel="log1p(posterior q05)",
        )
        axis.set_aspect("equal")

    for axis in axes.ravel()[len(cell_types) :]:
        axis.set_visible(False)

    if last_hexbin is not None:
        fig.colorbar(
            last_hexbin,
            ax=axes.ravel().tolist(),
            label="Spatial locations per hexagon",
            shrink=0.75,
        )
    fig.suptitle(f"{sample_name}: posterior mean versus q05", fontsize=14)
    fig.subplots_adjust(top=0.93, wspace=0.35, hspace=0.4)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_relative_posterior_uncertainty(
    sdata,
    shape_key,
    image_key,
    coordinate_system,
    cell_types,
    uncertainty_columns,
    uncertainty_vmax,
    minimum_abundance,
    output_path,
):
    """Plot relative posterior interval width for every cell type."""
    ncols = min(4, len(cell_types))
    nrows = math.ceil(len(cell_types) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 4.0 * nrows),
        squeeze=False,
    )

    for axis, cell_type, column in zip(axes.ravel(), cell_types, uncertainty_columns):
        if image_key:
            plot = sdata.pl.render_images(image_key).pl.render_shapes(
                shape_key,
                color=column,
                cmap="magma",
                norm=Normalize(vmin=0, vmax=uncertainty_vmax),
                fill_alpha=0.75,
                outline_width=0,
                na_color="#00000000",
            )
        else:
            plot = sdata.pl.render_shapes(
                shape_key,
                color=column,
                cmap="magma",
                norm=Normalize(vmin=0, vmax=uncertainty_vmax),
                fill_alpha=0.9,
                outline_width=0,
                na_color="#00000000",
            )
        plot.pl.show(coordinate_systems=coordinate_system, ax=axis)
        axis.set_title(cell_type, fontsize=9)
        axis.set_axis_off()

    for axis in axes.ravel()[len(cell_types) :]:
        axis.set_visible(False)

    fig.suptitle(
        f"Relative posterior interval width: (q95 - q05) / mean\n"
        f"Locations with mean < {minimum_abundance:g} are masked",
        fontsize=13,
    )
    fig.subplots_adjust(top=0.94, wspace=0.15, hspace=0.2)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_celltype_abundance(
    sdata,
    shape_key,
    image_key,
    coordinate_system,
    cell_types,
    abundance_columns,
    abundance,
    statistic,
    output_path,
):
    """Plot selected cell-abundance estimates for every cell type."""
    ncols = min(4, len(cell_types))
    nrows = math.ceil(len(cell_types) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 4.0 * nrows),
        squeeze=False,
    )

    for axis, cell_type, column in zip(axes.ravel(), cell_types, abundance_columns):
        values = abundance[cell_type].dropna().to_numpy()
        vmax = float(np.percentile(values, 99)) if values.size else 1.0
        vmax = max(vmax, np.finfo(float).eps)
        if image_key:
            plot = sdata.pl.render_images(image_key).pl.render_shapes(
                shape_key,
                color=column,
                cmap="viridis",
                norm=Normalize(vmin=0, vmax=vmax),
                fill_alpha=0.75,
                outline_width=0,
                na_color="#00000000",
            )
        else:
            plot = sdata.pl.render_shapes(
                shape_key,
                color=column,
                cmap="viridis",
                norm=Normalize(vmin=0, vmax=vmax),
                fill_alpha=0.9,
                outline_width=0,
                na_color="#00000000",
            )
        plot.pl.show(coordinate_systems=coordinate_system, ax=axis)
        axis.set_title(f"{cell_type}\ncolour maximum: p99 = {vmax:.3g}", fontsize=9)
        axis.set_axis_off()

    for axis in axes.ravel()[len(cell_types) :]:
        axis.set_visible(False)

    fig.suptitle(f"Cell-type abundance ({statistic})", fontsize=14)
    fig.subplots_adjust(top=0.95, wspace=0.15, hspace=0.25)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_celltype_fraction(
    sdata,
    shape_key,
    image_key,
    coordinate_system,
    cell_types,
    fraction_columns,
    statistic,
    output_path,
):
    """Plot within-location cell-type fractions for every cell type."""
    ncols = min(4, len(cell_types))
    nrows = math.ceil(len(cell_types) / ncols)
    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 4.0 * nrows),
        squeeze=False,
    )

    for axis, cell_type, column in zip(axes.ravel(), cell_types, fraction_columns):
        if image_key:
            plot = sdata.pl.render_images(image_key).pl.render_shapes(
                shape_key,
                color=column,
                cmap="viridis",
                norm=Normalize(vmin=0, vmax=1),
                fill_alpha=0.75,
                outline_width=0,
                na_color="#00000000",
            )
        else:
            plot = sdata.pl.render_shapes(
                shape_key,
                color=column,
                cmap="viridis",
                norm=Normalize(vmin=0, vmax=1),
                fill_alpha=0.9,
                outline_width=0,
                na_color="#00000000",
            )
        plot.pl.show(coordinate_systems=coordinate_system, ax=axis)
        axis.set_title(cell_type, fontsize=9)
        axis.set_axis_off()

    for axis in axes.ravel()[len(cell_types) :]:
        axis.set_visible(False)

    fig.suptitle(f"Cell-type fraction ({statistic})", fontsize=14)
    fig.subplots_adjust(top=0.95, wspace=0.15, hspace=0.2)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_total_abundance(
    sdata,
    shape_key,
    image_key,
    coordinate_system,
    total_column,
    total_abundance,
    statistic,
    sample_name,
    output_path,
):
    """Plot total-abundance distribution and spatial map."""
    values = total_abundance.dropna().to_numpy()
    mean_value = float(np.mean(values))
    median_value = float(np.median(values))
    vmax = max(float(np.percentile(values, 99)), np.finfo(float).eps)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    axes[0].hist(values, bins=50, color="#4c78a8", edgecolor="white", linewidth=0.3)
    axes[0].axvline(
        mean_value,
        color="#e45756",
        linestyle="--",
        linewidth=1.5,
        label=f"Mean: {mean_value:.2f}",
    )
    axes[0].axvline(
        median_value,
        color="#54a24b",
        linestyle=":",
        linewidth=1.8,
        label=f"Median: {median_value:.2f}",
    )
    axes[0].set(
        title="Distribution",
        xlabel="Total estimated cells per spatial location",
        ylabel="Number of spatial locations",
    )
    axes[0].legend(frameon=False)

    if image_key:
        plot = sdata.pl.render_images(image_key).pl.render_shapes(
            shape_key,
            color=total_column,
            cmap="viridis",
            norm=Normalize(vmin=0, vmax=vmax),
            fill_alpha=0.75,
            outline_width=0,
            na_color="#00000000",
        )
    else:
        plot = sdata.pl.render_shapes(
            shape_key,
            color=total_column,
            cmap="viridis",
            norm=Normalize(vmin=0, vmax=vmax),
            fill_alpha=0.9,
            outline_width=0,
            na_color="#00000000",
        )
    plot.pl.show(coordinate_systems=coordinate_system, ax=axes[1])
    axes[1].set_title(f"Spatial distribution\ncolour maximum: p99 = {vmax:.3g}")
    axes[1].set_axis_off()

    fig.suptitle(f"{sample_name}: total cell abundance ({statistic})", fontsize=14)
    fig.subplots_adjust(top=0.86, wspace=0.25)
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def plot_multiqc_deconvolution(
    sdata,
    shape_key,
    image_key,
    coordinate_system,
    total_column,
    total_abundance,
    selected_cell_types,
    mixture_column,
    mixture_palette,
    selected_coverage,
    statistic,
    sample_name,
    output_path,
):
    """Plot the compact two-panel deconvolution summary for MultiQC."""
    values = total_abundance.dropna().to_numpy()
    vmax = max(float(np.percentile(values, 99)), np.finfo(float).eps)
    colours = ["#d81b60", "#00a6d6", "#f4c430"]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    if image_key:
        total_plot = sdata.pl.render_images(image_key).pl.render_shapes(
            shape_key,
            color=total_column,
            cmap="viridis",
            norm=Normalize(vmin=0, vmax=vmax),
            fill_alpha=0.8,
            outline_width=0,
            na_color="#00000000",
            colorbar=False,
        )
    else:
        total_plot = sdata.pl.render_shapes(
            shape_key,
            color=total_column,
            cmap="viridis",
            norm=Normalize(vmin=0, vmax=vmax),
            fill_alpha=0.95,
            outline_width=0,
            na_color="#00000000",
            colorbar=False,
        )
    total_plot.pl.show(coordinate_systems=coordinate_system, ax=axes[0])
    axes[0].set_title(f"Total abundance ({statistic})\ncolour maximum: p99 = {vmax:.3g}")
    axes[0].set_axis_off()

    if image_key:
        mixture_plot = sdata.pl.render_images(image_key).pl.render_shapes(
            shape_key,
            color=mixture_column,
            groups=mixture_palette,
            palette=mixture_palette,
            fill_alpha=0.85,
            outline_width=0,
            na_color="#00000000",
            colorbar=False,
        )
    else:
        mixture_plot = sdata.pl.render_shapes(
            shape_key,
            color=mixture_column,
            groups=mixture_palette,
            palette=mixture_palette,
            fill_alpha=0.95,
            outline_width=0,
            na_color="#00000000",
            colorbar=False,
        )

    mixture_plot.pl.show(
        coordinate_systems=coordinate_system,
        ax=axes[1],
        legend_loc=None,
        na_in_legend=False,
    )
    axes[1].set_title(
        "Cell-type composition\n"
        f"selected types cover {selected_coverage:.1%} of slide abundance"
    )
    axes[1].set_axis_off()

    legend_handles = [
        Patch(facecolor=colour, label=cell_type)
        for cell_type, colour in zip(selected_cell_types, colours)
    ]
    legend_handles.append(Patch(facecolor="#9e9e9e", label="Other cell types"))
    axes[1].legend(
        handles=legend_handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.12),
        ncol=2,
        frameon=False,
        fontsize=8,
    )

    for extra_axis in list(fig.axes):
        if all(extra_axis is not axis for axis in axes.ravel()):
            extra_axis.remove()

    fig.suptitle(f"{sample_name}: Cell2location deconvolution", fontsize=14)
    fig.subplots_adjust(top=0.88, bottom=0.15, wspace=0.08)
    fig.savefig(output_path, dpi=160)
    plt.close(fig)


def main():
    """Read Cell2location results and create plots and MultiQC input tables."""
    args = parse_args()

    if not args.input_zarr.exists():
        raise FileNotFoundError(f"Input SpatialData does not exist: {args.input_zarr}")
    if args.uncertainty_min_abundance < 0:
        raise ValueError("--uncertainty-min-abundance must be greater than or equal to zero")

    requested_cell_types = []
    for value in args.multiqc_cell_types or []:
        requested_cell_types.extend(item.strip() for item in value.split(",") if item.strip())
    if len(requested_cell_types) > 3:
        raise ValueError("--multiqc-cell-types accepts at most three cell types")
    if len(requested_cell_types) != len(set(requested_cell_types)):
        raise ValueError("--multiqc-cell-types contains duplicate cell types")

    print(f"Reading SpatialData: {args.input_zarr}")
    sdata = sd.read_zarr(args.input_zarr)
    if args.table_key not in sdata.tables:
        raise KeyError(f"Table '{args.table_key}' was not found in the input SpatialData")

    shape_key = args.shape_key or f"{args.sample_name}_{args.table_key}"
    if shape_key not in sdata.shapes:
        raise KeyError(f"Shapes element '{shape_key}' was not found in the input SpatialData")
    if args.sample_name not in sdata.coordinate_systems:
        raise KeyError(f"Coordinate system '{args.sample_name}' was not found")

    default_image_key = f"{args.sample_name}_hires_image"
    image_key = args.image_key
    if image_key is None and default_image_key in sdata.images:
        image_key = default_image_key
    if image_key is not None and image_key not in sdata.images:
        raise KeyError(f"Image element '{image_key}' was not found in the input SpatialData")

    table = sdata.tables[args.table_key]
    posterior_keys = {
        "mean": "means_cell_abundance_w_sf",
        "std": "stds_cell_abundance_w_sf",
        "q05": "q05_cell_abundance_w_sf",
        "q95": "q95_cell_abundance_w_sf",
    }
    missing_keys = [key for key in posterior_keys.values() if key not in table.obsm]
    if missing_keys:
        raise KeyError(f"Missing Cell2location posterior matrices: {', '.join(missing_keys)}")

    cell2location_metadata = table.uns.get("cell2location")
    if not isinstance(cell2location_metadata, dict):
        raise ValueError("table.uns['cell2location'] is missing")
    metadata_cell_types = cell2location_metadata.get("cell_types")
    posterior = {}
    for summary_name, key in posterior_keys.items():
        values = table.obsm[key]
        if isinstance(values, pd.DataFrame):
            frame = values.copy()
            frame.index = table.obs_names
            frame.columns = frame.columns.astype(str)
        else:
            array = np.asarray(values)
            if metadata_cell_types is None:
                raise ValueError(
                    f"table.obsm['{key}'] has no column labels and "
                    "table.uns['cell2location']['cell_types'] is missing"
                )
            frame = pd.DataFrame(
                array,
                index=table.obs_names,
                columns=pd.Index(metadata_cell_types, dtype=str),
            )
        posterior[summary_name] = frame.astype(float)

    cell_types = posterior["mean"].columns.tolist()
    if not cell_types:
        raise ValueError("No Cell2location cell types were found")
    if len(cell_types) != len(set(cell_types)):
        raise ValueError("Cell2location cell-type names are not unique")
    invalid_cell_types = [name for name in cell_types if "\t" in name or "\n" in name]
    if invalid_cell_types:
        raise ValueError("Cell-type names cannot contain tabs or newlines")
    if "\t" in args.sample_name or "\n" in args.sample_name:
        raise ValueError("--sample-name cannot contain tabs or newlines")
    output_prefix = args.output_prefix or args.sample_name
    if "\t" in output_prefix or "\n" in output_prefix:
        raise ValueError("--output-prefix cannot contain tabs or newlines")
    for summary_name, frame in posterior.items():
        if not frame.index.equals(posterior["mean"].index):
            raise ValueError(f"Posterior matrix '{summary_name}' has different location labels")
        if frame.columns.tolist() != cell_types:
            raise ValueError(f"Posterior matrix '{summary_name}' has different cell-type labels")

    mean_present = posterior["mean"].notna()
    modelled = mean_present.all(axis=1)
    if (mean_present.any(axis=1) & ~modelled).any():
        raise ValueError("Posterior mean contains partially missing rows")
    if not modelled.any():
        raise ValueError("No spatial locations contain Cell2location posterior estimates")

    for summary_name, frame in posterior.items():
        present = frame.notna()
        if (present.any(axis=1) & ~present.all(axis=1)).any():
            raise ValueError(f"Posterior matrix '{summary_name}' contains partially missing rows")
        if not frame.notna().all(axis=1).equals(modelled):
            raise ValueError(
                f"Posterior matrix '{summary_name}' has a different missing-value mask"
            )
        modelled_values = frame.loc[modelled].to_numpy()
        if not np.isfinite(modelled_values).all():
            raise ValueError(f"Posterior matrix '{summary_name}' contains non-finite values")
        if (modelled_values < 0).any():
            raise ValueError(f"Posterior matrix '{summary_name}' contains negative values")

    missing_requested = [name for name in requested_cell_types if name not in cell_types]
    if missing_requested:
        raise ValueError(
            "Requested MultiQC cell types were not found: " + ", ".join(missing_requested)
        )

    abundance = posterior[args.abundance_statistic]
    total_abundance = abundance.sum(axis=1, min_count=1)
    fractions = abundance.div(total_abundance.replace(0, np.nan), axis=0)
    uncertainty = (posterior["q95"] - posterior["q05"]).div(posterior["mean"])
    uncertainty = uncertainty.where(posterior["mean"] > args.uncertainty_min_abundance)

    abundance_columns = []
    fraction_columns = []
    uncertainty_columns = []
    for index, cell_type in enumerate(cell_types):
        abundance_column = f"_c2l_plot_abundance_{index}"
        fraction_column = f"_c2l_plot_fraction_{index}"
        uncertainty_column = f"_c2l_plot_uncertainty_{index}"
        table.obs[abundance_column] = abundance[cell_type].to_numpy()
        table.obs[fraction_column] = fractions[cell_type].to_numpy()
        table.obs[uncertainty_column] = uncertainty[cell_type].to_numpy()
        abundance_columns.append(abundance_column)
        fraction_columns.append(fraction_column)
        uncertainty_columns.append(uncertainty_column)

    total_column = "_c2l_plot_total_abundance"
    table.obs[total_column] = total_abundance.to_numpy()

    slide_abundance = abundance.loc[modelled].sum(axis=0)
    slide_total = float(slide_abundance.sum())
    if slide_total <= 0:
        raise ValueError("Total estimated abundance across the slide is zero")
    slide_composition = slide_abundance / slide_total

    if requested_cell_types:
        selected_cell_types = requested_cell_types
    else:
        selected_cell_types = slide_abundance.sort_values(ascending=False).head(3).index.tolist()
    selected_coverage = float(slide_composition.loc[selected_cell_types].sum())

    selected_fractions = fractions[selected_cell_types]
    other_fraction = (1 - selected_fractions.sum(axis=1, min_count=1)).clip(lower=0, upper=1)
    basis_colours = np.asarray(
        [to_rgb(colour) for colour in ["#d81b60", "#00a6d6", "#f4c430"]]
    )
    mixture_rgb = (
        selected_fractions.fillna(0).to_numpy() @ basis_colours[: len(selected_cell_types)]
        + other_fraction.fillna(0).to_numpy()[:, None] * np.asarray(to_rgb("#9e9e9e"))
    )
    mixture_rgb = np.round(np.clip(mixture_rgb, 0, 1) * 15) / 15
    mixture_values = pd.Series(pd.NA, index=table.obs_names, dtype="object")
    mixture_values.loc[modelled] = [
        to_hex(colour) for colour in mixture_rgb[modelled.to_numpy()]
    ]
    mixture_column = "_c2l_plot_mixture_colour"
    table.obs[mixture_column] = pd.Categorical(mixture_values)
    mixture_palette = table.obs[mixture_column].cat.categories.tolist()

    finite_uncertainty = uncertainty.to_numpy()[np.isfinite(uncertainty.to_numpy())]
    uncertainty_vmax = (
        float(np.percentile(finite_uncertainty, 99)) if finite_uncertainty.size else 1.0
    )
    uncertainty_vmax = max(uncertainty_vmax, np.finfo(float).eps)

    sdata.tables[args.table_key] = table
    minx, miny, maxx, maxy = sdata.shapes[shape_key].total_bounds
    sdata_crop = sdata.query.bounding_box(
        axes=("x", "y"),
        min_coordinate=[minx, miny],
        max_coordinate=[maxx, maxy],
        target_coordinate_system=args.sample_name,
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    prefix = args.output_dir / output_prefix

    training_history = read_training_history(cell2location_metadata)
    training_history.to_csv(f"{prefix}_training_history.csv", index=False)
    plot_training_history(training_history, Path(f"{prefix}_training_history.png"))

    reconstruction_qc = cell2location_metadata.get("reconstruction_qc", {})
    reconstruction_plot_status = "failed"
    reconstruction_error = str(reconstruction_qc.get("error", ""))
    if plot_reconstruction_accuracy(
        reconstruction_qc,
        Path(f"{prefix}_reconstruction_accuracy.png"),
    ):
        reconstruction_plot_status = "success"
    else:
        if not reconstruction_error:
            reconstruction_error = "Reconstruction QC data are unavailable"
        Path(f"{prefix}_reconstruction_accuracy_error.txt").write_text(
            reconstruction_error + "\n",
            encoding="utf-8",
        )

    diagnostics = pd.DataFrame(
        [
            {
                "total_locations": int(table.n_obs),
                "modelled_locations": int(modelled.sum()),
                "shared_genes": int(
                    cell2location_metadata.get("n_shared_genes", table.n_vars)
                ),
                "cell_types": int(len(cell_types)),
                "epochs_completed": int(len(training_history)),
                "initial_elbo": float(training_history["elbo_train"].iloc[0]),
                "final_elbo": float(training_history["elbo_train"].iloc[-1]),
                "reconstruction_plot": reconstruction_plot_status,
                "reconstruction_error": reconstruction_error,
            }
        ]
    )
    diagnostics.to_csv(f"{prefix}_diagnostics.tsv", sep="\t", index=False)

    plot_posterior_mean_vs_q05(
        posterior["mean"].loc[modelled],
        posterior["q05"].loc[modelled],
        args.sample_name,
        Path(f"{prefix}_posterior_mean_vs_q05.pdf"),
    )
    plot_relative_posterior_uncertainty(
        sdata_crop,
        shape_key,
        image_key,
        args.sample_name,
        cell_types,
        uncertainty_columns,
        uncertainty_vmax,
        args.uncertainty_min_abundance,
        Path(f"{prefix}_relative_posterior_uncertainty_spatial.pdf"),
    )
    plot_celltype_abundance(
        sdata_crop,
        shape_key,
        image_key,
        args.sample_name,
        cell_types,
        abundance_columns,
        abundance,
        args.abundance_statistic,
        Path(f"{prefix}_celltype_abundance_{args.abundance_statistic}.pdf"),
    )
    plot_celltype_fraction(
        sdata_crop,
        shape_key,
        image_key,
        args.sample_name,
        cell_types,
        fraction_columns,
        args.abundance_statistic,
        Path(f"{prefix}_celltype_fraction_{args.abundance_statistic}.pdf"),
    )
    plot_total_abundance(
        sdata_crop,
        shape_key,
        image_key,
        args.sample_name,
        total_column,
        total_abundance,
        args.abundance_statistic,
        args.sample_name,
        Path(f"{prefix}_total_abundance_{args.abundance_statistic}.pdf"),
    )
    plot_multiqc_deconvolution(
        sdata_crop,
        shape_key,
        image_key,
        args.sample_name,
        total_column,
        total_abundance,
        selected_cell_types,
        mixture_column,
        mixture_palette,
        selected_coverage,
        args.abundance_statistic,
        args.sample_name,
        Path(f"{prefix}_deconvolution_mqc.png"),
    )

    composition_path = Path(
        f"{prefix}_celltype_composition_{args.abundance_statistic}_mqc.tsv"
    )
    with composition_path.open("w", encoding="utf-8") as handle:
        handle.write("# id: 'cell2location_celltype_composition'\n")
        handle.write("# section_name: 'Cell2location cell-type composition'\n")
        handle.write(
            "# description: 'Whole-slide cell-type fractions calculated from "
            f"posterior {args.abundance_statistic} abundance.'\n"
        )
        handle.write("# plot_type: 'bargraph'\n")
        handle.write("# pconfig:\n")
        handle.write("#   id: 'cell2location_celltype_composition_plot'\n")
        handle.write("#   title: 'Cell-type composition by sample'\n")
        handle.write("#   ylab: 'Fraction of estimated cells'\n")
        handle.write("Sample\t" + "\t".join(cell_types) + "\n")
        handle.write(
            args.sample_name
            + "\t"
            + "\t".join(f"{slide_composition[cell_type]:.10g}" for cell_type in cell_types)
            + "\n"
        )

    summary_path = Path(f"{prefix}_deconvolution_summary_{args.abundance_statistic}_mqc.tsv")
    with summary_path.open("w", encoding="utf-8") as handle:
        handle.write("# id: 'cell2location_deconvolution_summary'\n")
        handle.write("# section_name: 'Cell2location deconvolution summary'\n")
        handle.write(
            "# description: 'Summary statistics calculated from posterior "
            f"{args.abundance_statistic} abundance.'\n"
        )
        handle.write("# plot_type: 'table'\n")
        handle.write("# pconfig:\n")
        handle.write("#   id: 'cell2location_deconvolution_summary_table'\n")
        handle.write(
            "Sample\tSpatial locations\tTotal abundance\tMean abundance per location\t"
            "Median abundance per location\n"
        )
        handle.write(
            f"{args.sample_name}\t{int(modelled.sum())}\t{slide_total:.10g}\t"
            f"{total_abundance.loc[modelled].mean():.10g}\t"
            f"{total_abundance.loc[modelled].median():.10g}\n"
        )

    print(f"Wrote Cell2location plots and tables to: {args.output_dir}")
    print(f"Abundance statistic: {args.abundance_statistic}")
    print(f"Compact MultiQC cell types: {', '.join(selected_cell_types)}")


if __name__ == "__main__":
    main()
