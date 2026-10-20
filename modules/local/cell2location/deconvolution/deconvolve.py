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
    parser.add_argument("--output-abundance", type=Path)
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
        for name in ("input_zarr", "signatures", "output_zarr", "output_abundance")
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
    }
    for name, value in positive_parameters.items():
        if value <= 0:
            raise ValueError(f"{name} must be greater than zero")

    import numpy as np
    import pandas as pd
    import spatialdata as sd
    import torch
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
    adata = model.export_posterior(
        adata,
        sample_kwargs={
            "num_samples": args.posterior_num_samples,
            "batch_size": args.posterior_batch_size,
            "accelerator": accelerator,
        },
    )

    abundance_key = "q05_cell_abundance_w_sf"
    if abundance_key not in adata.obsm:
        raise KeyError(f"Cell2location did not create adata.obsm['{abundance_key}']")
    abundance = adata.obsm[abundance_key]
    if not isinstance(abundance, pd.DataFrame):
        abundance = pd.DataFrame(
            abundance,
            index=adata.obs_names,
            columns=signatures.columns,
        )

    full_abundance = pd.DataFrame(
        np.nan,
        index=table.obs_names,
        columns=abundance.columns,
    )
    full_abundance.loc[adata.obs_names, :] = abundance.to_numpy()
    table.obsm[abundance_key] = full_abundance

    output_columns = [f"c2l_{str(name).replace(' ', '_')}" for name in abundance.columns]
    if len(output_columns) != len(set(output_columns)):
        raise ValueError("Cell-type names are not unique after replacing spaces with underscores")
    for source_name, output_name in zip(abundance.columns, output_columns):
        table.obs[output_name] = full_abundance[source_name].to_numpy()
    table.obs["c2l_total_abundance"] = full_abundance.sum(axis=1, min_count=1).to_numpy()

    cell2location_metadata = dict(adata.uns["mod"])
    cell2location_metadata.update(
        {
            "table_key": args.table_key,
            "qc_pass_key": args.qc_pass_key,
            "abundance_key": abundance_key,
            "n_cells_per_location": float(args.n_cells_per_location),
            "detection_alpha": float(args.detection_alpha),
        }
    )
    table.uns["cell2location"] = cell2location_metadata
    sdata.tables[args.table_key] = table

    args.output_zarr.parent.mkdir(parents=True, exist_ok=True)
    args.output_abundance.parent.mkdir(parents=True, exist_ok=True)
    sdata.write(args.output_zarr)
    full_abundance.to_csv(args.output_abundance, index_label="location")
    print(f"Wrote deconvolved SpatialData: {args.output_zarr}")
    print(f"Wrote cell-type abundance table: {args.output_abundance}")


if __name__ == "__main__":
    main()
