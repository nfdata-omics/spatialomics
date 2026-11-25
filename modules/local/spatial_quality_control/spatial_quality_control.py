
"""Spatial quality control module for spatial omics data."""

import importlib
import importlib.metadata
import argparse
import sys
from pathlib import Path
import yaml

import scipy as sp
import pandas as pd
import scanpy as sc
import spatialdata
import spotsweeper.local_outliers as lo


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
    local outliers, and saves the QC annotations and summary outputs.

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
            "spotsweeper",
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
