
"""Spatial quality control module for spatial omics data."""

import importlib
import importlib.metadata
import argparse
import sys
import yaml

import scipy as sp
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import spatialdata
import spatialdata_plot # noqa: F401 # pyright: ignore[reportUnusedImport] # pylint: disable=unused-import
import spotsweeper.local_outliers as lo

matplotlib.use("Agg")

def qc_from_h5ad(
    zarr_folder,
    sample_id,
    resolution="016um",
    min_counts=100,
    min_genes=50,
    max_mt=20,
    novelty_thresh=0.0
):
    """
    Perform spatial quality control analysis on spatial omics data.

    Reads a spatialdata zarr object, calculates QC metrics, identifies global and
    local outliers, generates distribution plots and spatial visualizations, and
    saves annotated observations and QC summary statistics.

    Parameters
    ----------
    zarr_folder : str
        Path to the input Zarr object containing spatialdata.
    sample_id : str
        Unique identifier for the sample being processed.
    resolution : str, optional
        Resolution of the spatial data to process. Default is "016um".
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
    sdata = spatialdata.read_zarr(zarr_folder)
    # adata = data.tables[f'square_{resolution}']
    adata = sdata.tables[f'square_{resolution}']
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

    # save annotation to csv
    adata.obs.to_csv(f"{sample_id}_qc_annotated_obs.csv")

    # Save AnnData with quality control annotation
    adata.write(f"{sample_id}_{resolution}_qc.h5ad")

    # Collect QC summary
    summary_dict = {
        "Sample": sample_id,
        "Total spots": adata.n_obs,
        "In-tissue spots": (adata.obs['in_tissue'] == 1).sum(),
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

    # Spatial QC plot
    spatial_colors = [
        "total_counts",
        "n_genes_by_counts",
        "pct_counts_mt",
        "global_outliers",
        "local_outliers"
    ]
    with plt.rc_context({
        "font.size": 8,
        "axes.titlesize": 9,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
    }):
        fig, axs = plt.subplots(2, 3, figsize=(14, 8), constrained_layout=False)
        axs = axs.ravel()

        for i, color in enumerate(spatial_colors):
            sdata.pl.render_shapes(  # pylint: disable=no-member
                f"{sample_id}_square_{resolution}",
                color=color,
                cmap="viridis",
            ).pl.show(
                coordinate_systems=sample_id,
                title=color,
                ax=axs[i],
            )
            axs[i].set_title(color, fontsize=9)

        axs[5].axis("off")

        # More space between panels
        fig.subplots_adjust(wspace=0.35, hspace=0.35)

        for a in fig.axes:
            if a not in axs.flat:
                a.set_ylabel("")  # clear colorbar label
                a.set_title("")   # clear colorbar title if present
                a.tick_params(labelsize=7)

        fig.savefig(f"{sample_id}_qc_spatial_plots.png")


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
    parser.add_argument("--sample", type=str,
                        help="Sample ID to process")
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
        if not args.zarr or not args.sample:
            parser.error("--zarr and --sample are required")
        qc_from_h5ad(args.zarr, args.sample)
