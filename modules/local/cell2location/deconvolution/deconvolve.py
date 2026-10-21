#!/usr/bin/env python3
"""Deconvolve one SpatialData table using precomputed cell2location signatures."""

import argparse
import random
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run cell2location deconvolution on one SpatialData sample"
    )
    parser.add_argument("--input-zarr", type=Path)
    parser.add_argument("--signatures", type=Path)
    parser.add_argument("--output-zarr", type=Path)
    parser.add_argument("--output-abundance-mean", type=Path)
    parser.add_argument("--output-abundance-q05", type=Path)
    parser.add_argument("--output-diagnostics", type=Path)
    parser.add_argument("--table-key", default="square_016um")
    parser.add_argument("--qc-pass-key", default="qc_pass")
    parser.add_argument(
        "--counts-layer",
        default="X",
        help="Spatial layer containing raw counts; use X for the main matrix",
    )
    parser.add_argument(
        "--gene-key",
        default="gene_ids",
        help="Spatial var column containing gene IDs, or var_names",
    )
    parser.add_argument("--max-epochs", type=int, default=200)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--n-cells-per-location", type=float, default=3.0)
    parser.add_argument("--detection-alpha", type=float, default=20.0)
    parser.add_argument("--posterior-num-samples", type=int, default=1000)
    parser.add_argument("--posterior-batch-size", type=int, default=512)
    parser.add_argument("--qc-num-locations", type=int, default=1000)
    parser.add_argument("--minimum-shared-genes", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--accelerator",
        choices=("auto", "cpu", "gpu"),
        default="auto",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    missing = [
        name
        for name in (
            "input_zarr",
            "signatures",
            "output_zarr",
            "output_abundance_mean",
            "output_abundance_q05",
            "output_diagnostics",
        )
        if getattr(args, name) is None
    ]
    if missing:
        raise ValueError(f"Missing required arguments: {', '.join(missing)}")
    if not args.input_zarr.exists():
        raise FileNotFoundError(f"Input SpatialData does not exist: {args.input_zarr}")
    if not args.signatures.is_file():
        raise FileNotFoundError(f"Reference signatures do not exist: {args.signatures}")
    if args.input_zarr.resolve() == args.output_zarr.resolve():
        raise ValueError("The input and output SpatialData Zarr paths must be different")

    positive_parameters = {
        "--max-epochs": args.max_epochs,
        "--batch-size": args.batch_size,
        "--n-cells-per-location": args.n_cells_per_location,
        "--detection-alpha": args.detection_alpha,
        "--posterior-num-samples": args.posterior_num_samples,
        "--posterior-batch-size": args.posterior_batch_size,
        "--qc-num-locations": args.qc_num_locations,
    }
    for name, value in positive_parameters.items():
        if value <= 0:
            raise ValueError(f"{name} must be greater than zero")

    import numpy as np
    import pandas as pd
    import spatialdata as sd
    import torch
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    import cell2location
    from cell2location.models import Cell2location

    accelerator = args.accelerator
    if accelerator == "auto":
        accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    if accelerator == "gpu" and not torch.cuda.is_available():
        raise RuntimeError("--accelerator gpu was requested, but CUDA is not available")

    np.random.seed(args.seed)
    random.seed(args.seed)
    torch.manual_seed(args.seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(args.seed)

    print(f"Reading SpatialData: {args.input_zarr}")
    sdata = sd.read_zarr(args.input_zarr)
    if args.table_key not in sdata.tables:
        raise KeyError(f"Table '{args.table_key}' was not found in the input SpatialData")
    table = sdata.tables[args.table_key]
    if args.qc_pass_key not in table.obs:
        raise KeyError(f"Table obs does not contain '{args.qc_pass_key}'")
    if args.counts_layer != "X" and args.counts_layer not in table.layers:
        raise KeyError(f"Spatial table layers do not contain '{args.counts_layer}'")

    qc_pass = table.obs[args.qc_pass_key].fillna(False).astype(bool)
    if not qc_pass.any():
        raise ValueError(f"No spatial locations have {args.qc_pass_key}=True")
    adata = table[qc_pass].copy()
    print(f"Locations passing QC: {adata.n_obs} / {table.n_obs}")

    signatures = pd.read_csv(args.signatures, index_col=0)
    signatures.index = signatures.index.astype(str)
    if signatures.index.has_duplicates:
        duplicates = signatures.index[signatures.index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Signature gene IDs are not unique; examples: {duplicates}")

    if args.gene_key == "var_names":
        spatial_gene_ids = pd.Index(adata.var_names.astype(str))
    else:
        if args.gene_key not in adata.var:
            raise KeyError(f"Spatial table var does not contain '{args.gene_key}'")
        spatial_gene_ids = pd.Index(adata.var[args.gene_key].astype(str))
    if spatial_gene_ids.has_duplicates:
        duplicates = spatial_gene_ids[spatial_gene_ids.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Spatial gene IDs are not unique; examples: {duplicates}")

    shared_genes = spatial_gene_ids.intersection(signatures.index, sort=False)
    print(f"Shared genes: {len(shared_genes)}")
    if len(shared_genes) < args.minimum_shared_genes:
        raise ValueError(
            f"Only {len(shared_genes)} genes are shared between reference and spatial data; "
            f"expected at least {args.minimum_shared_genes}. Check the gene ID options."
        )

    keep_genes = spatial_gene_ids.isin(shared_genes)
    adata = adata[:, keep_genes].copy()
    adata.var_names = spatial_gene_ids[keep_genes]
    signatures = signatures.loc[adata.var_names].copy()

    layer = None if args.counts_layer == "X" else args.counts_layer
    Cell2location.setup_anndata(adata=adata, layer=layer)
    model = Cell2location(
        adata,
        cell_state_df=signatures,
        N_cells_per_location=args.n_cells_per_location,
        detection_alpha=args.detection_alpha,
    )

    print(f"Training spatial model on {accelerator}")
    model.train(
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        train_size=1,
        accelerator=accelerator,
    )

    training_history = pd.DataFrame(
        {
            "epoch": np.asarray(model.history_["elbo_train"].index),
            "elbo_train": np.asarray(model.history_["elbo_train"]).reshape(-1),
        }
    )
    training_history["elbo_train"] = training_history["elbo_train"].astype(float)

    adata = model.export_posterior(
        adata,
        sample_kwargs={
            "num_samples": args.posterior_num_samples,
            "batch_size": args.posterior_batch_size,
            "accelerator": accelerator,
        },
        add_to_obsm=["means", "stds", "q05", "q95"],
    )

    cell_types = signatures.columns.astype(str).tolist()
    abundance_keys = [
        "means_cell_abundance_w_sf",
        "stds_cell_abundance_w_sf",
        "q05_cell_abundance_w_sf",
        "q95_cell_abundance_w_sf",
    ]
    full_abundances = {}
    for abundance_key in abundance_keys:
        if abundance_key not in adata.obsm:
            raise KeyError(f"Cell2location did not create adata.obsm['{abundance_key}']")
        abundance_values = np.asarray(adata.obsm[abundance_key], dtype=np.float32)
        if abundance_values.shape != (adata.n_obs, len(cell_types)):
            raise ValueError(
                f"Unexpected shape for adata.obsm['{abundance_key}']: "
                f"{abundance_values.shape}, expected {(adata.n_obs, len(cell_types))}"
            )
        full_abundance = pd.DataFrame(
            np.full((table.n_obs, len(cell_types)), np.nan, dtype=np.float32),
            index=table.obs_names,
            columns=cell_types,
        )
        full_abundance.loc[adata.obs_names, :] = abundance_values
        table.obsm[abundance_key] = full_abundance
        full_abundances[abundance_key] = full_abundance

    posterior_metadata = adata.uns.get("mod", {})
    detection_sensitivity_keys = {}
    for summary_name, output_suffix in (
        ("means", "mean"),
        ("stds", "std"),
        ("q05", "q05"),
        ("q95", "q95"),
    ):
        posterior = posterior_metadata.get(f"post_sample_{summary_name}", {})
        if "detection_y_s" not in posterior:
            continue
        values = np.asarray(posterior["detection_y_s"], dtype=np.float32).reshape(-1)
        if values.size != adata.n_obs:
            raise ValueError(
                f"Unexpected shape for post_sample_{summary_name}['detection_y_s']: "
                f"{np.asarray(posterior['detection_y_s']).shape}"
            )
        output_key = f"c2l_detection_sensitivity_{output_suffix}"
        full_values = np.full(table.n_obs, np.nan, dtype=np.float32)
        full_values[qc_pass.to_numpy()] = values
        table.obs[output_key] = full_values
        detection_sensitivity_keys[summary_name] = output_key

    args.output_diagnostics.mkdir(parents=True, exist_ok=True)
    training_history.to_csv(args.output_diagnostics / "training_history.csv", index=False)

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
    fig.savefig(args.output_diagnostics / "training_history.png", dpi=200)
    plt.close(fig)

    reconstruction_status = "success"
    reconstruction_error = ""
    try:
        np.random.seed(args.seed)
        model.plot_QC(
            summary_name="means",
            use_n_obs=min(args.qc_num_locations, adata.n_obs),
        )
        fig = plt.gcf()
        fig.savefig(args.output_diagnostics / "reconstruction_accuracy.png", dpi=200)
        plt.close(fig)
    except Exception as error:  # Keep the expensive deconvolution result if diagnostic plotting fails.
        reconstruction_status = "failed"
        reconstruction_error = str(error).replace("\t", " ").replace("\n", " ")
        (args.output_diagnostics / "reconstruction_accuracy_error.txt").write_text(
            f"{type(error).__name__}: {error}\n",
            encoding="utf-8",
        )
        plt.close("all")

    final_elbo = float(training_history["elbo_train"].iloc[-1])
    diagnostics = pd.DataFrame(
        [
            {
                "total_locations": int(table.n_obs),
                "modelled_locations": int(adata.n_obs),
                "shared_genes": int(len(shared_genes)),
                "cell_types": int(len(cell_types)),
                "epochs_completed": int(len(training_history)),
                "initial_elbo": float(training_history["elbo_train"].iloc[0]),
                "final_elbo": final_elbo,
                "reconstruction_plot": reconstruction_status,
                "reconstruction_error": reconstruction_error,
            }
        ]
    )
    diagnostics.to_csv(args.output_diagnostics / "diagnostics.tsv", sep="\t", index=False)

    cell2location_metadata = {
        "model_name": posterior_metadata.get("model_name", "Cell2location"),
        "cell2location_version": getattr(cell2location, "__version__", "unknown"),
        "table_key": args.table_key,
        "qc_pass_key": args.qc_pass_key,
        "counts_layer": args.counts_layer,
        "gene_key": args.gene_key,
        "abundance_keys": abundance_keys,
        "detection_sensitivity_keys": detection_sensitivity_keys,
        "cell_types": cell_types,
        "n_total_locations": int(table.n_obs),
        "n_modelled_locations": int(adata.n_obs),
        "n_input_genes": int(table.n_vars),
        "n_shared_genes": int(len(shared_genes)),
        "n_cell_types": int(len(cell_types)),
        "n_cells_per_location": float(args.n_cells_per_location),
        "detection_alpha": float(args.detection_alpha),
        "max_epochs": int(args.max_epochs),
        "epochs_completed": int(len(training_history)),
        "batch_size": int(args.batch_size),
        "posterior_num_samples": int(args.posterior_num_samples),
        "posterior_batch_size": int(args.posterior_batch_size),
        "seed": int(args.seed),
        "accelerator": accelerator,
        "initial_elbo": float(training_history["elbo_train"].iloc[0]),
        "final_elbo": final_elbo,
        "training_history": training_history["elbo_train"].tolist(),
        "reconstruction_plot": reconstruction_status,
    }
    table.uns["cell2location"] = cell2location_metadata
    sdata.tables[args.table_key] = table

    args.output_zarr.parent.mkdir(parents=True, exist_ok=True)
    args.output_abundance_mean.parent.mkdir(parents=True, exist_ok=True)
    args.output_abundance_q05.parent.mkdir(parents=True, exist_ok=True)
    sdata.write(args.output_zarr)
    full_abundances["means_cell_abundance_w_sf"].to_csv(
        args.output_abundance_mean,
        index_label="location",
    )
    full_abundances["q05_cell_abundance_w_sf"].to_csv(
        args.output_abundance_q05,
        index_label="location",
    )
    print(f"Wrote deconvolved SpatialData: {args.output_zarr}")
    print(f"Wrote mean cell-type abundance table: {args.output_abundance_mean}")
    print(f"Wrote q05 cell-type abundance table: {args.output_abundance_q05}")
    print(f"Wrote model diagnostics: {args.output_diagnostics}")


if __name__ == "__main__":
    main()
