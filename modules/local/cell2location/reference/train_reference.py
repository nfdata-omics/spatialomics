#!/usr/bin/env python3
"""Train a cell2location reference model and export cell-type signatures."""

import argparse
import random
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train a cell2location RegressionModel on a single-cell reference"
    )
    parser.add_argument("--reference-h5ad", type=Path)
    parser.add_argument("--output-model", type=Path)
    parser.add_argument("--output-signatures", type=Path)
    parser.add_argument("--cell-type-key", default="celltype")
    parser.add_argument(
        "--batch-key",
        default="Batch_ID",
        help="Reference obs batch column; use 'none' to disable batch correction",
    )
    parser.add_argument(
        "--counts-layer",
        default="raw_counts",
        help="Reference layer containing raw counts; use X for reference.X",
    )
    parser.add_argument(
        "--gene-key",
        default="var_names",
        help="Reference var column containing gene IDs, or var_names",
    )
    parser.add_argument("--cell-count-cutoff", type=int, default=5)
    parser.add_argument("--cell-percentage-cutoff", type=float, default=0.03)
    parser.add_argument("--nonzero-mean-cutoff", type=float, default=1.12)
    parser.add_argument("--max-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=2500)
    parser.add_argument("--posterior-num-samples", type=int, default=1000)
    parser.add_argument("--posterior-batch-size", type=int, default=2500)
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
        for name in ("reference_h5ad", "output_model", "output_signatures")
        if getattr(args, name) is None
    ]
    if missing:
        raise ValueError(f"Missing required arguments: {', '.join(missing)}")
    if not args.reference_h5ad.is_file():
        raise FileNotFoundError(f"Reference AnnData does not exist: {args.reference_h5ad}")

    positive_parameters = {
        "--max-epochs": args.max_epochs,
        "--batch-size": args.batch_size,
        "--posterior-num-samples": args.posterior_num_samples,
        "--posterior-batch-size": args.posterior_batch_size,
    }
    for name, value in positive_parameters.items():
        if value <= 0:
            raise ValueError(f"{name} must be greater than zero")

    import numpy as np
    import pandas as pd
    import scanpy as sc
    import torch
    from cell2location.models import RegressionModel
    from cell2location.utils.filtering import filter_genes

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

    print(f"Reading reference: {args.reference_h5ad}")
    adata = sc.read_h5ad(args.reference_h5ad)

    if args.cell_type_key not in adata.obs:
        raise KeyError(f"Reference obs does not contain '{args.cell_type_key}'")
    batch_key = None if args.batch_key.lower() == "none" else args.batch_key
    if batch_key and batch_key not in adata.obs:
        raise KeyError(f"Reference obs does not contain '{batch_key}'")
    if args.counts_layer != "X":
        if args.counts_layer not in adata.layers:
            raise KeyError(f"Reference layers do not contain '{args.counts_layer}'")
        adata.X = adata.layers[args.counts_layer].copy()

    adata.obs[args.cell_type_key] = adata.obs[args.cell_type_key].astype("category")
    selected_genes = filter_genes(
        adata,
        cell_count_cutoff=args.cell_count_cutoff,
        cell_percentage_cutoff2=args.cell_percentage_cutoff,
        nonz_mean_cutoff=args.nonzero_mean_cutoff,
    )
    adata = adata[:, selected_genes].copy()
    print(f"Reference after gene filtering: {adata.shape}")

    setup_kwargs = {"labels_key": args.cell_type_key}
    if batch_key:
        setup_kwargs["batch_key"] = batch_key
    RegressionModel.setup_anndata(adata=adata, **setup_kwargs)

    print(f"Training reference model on {accelerator}")
    model = RegressionModel(adata)
    model.train(
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        train_size=1,
        lr=0.002,
        accelerator=accelerator,
    )
    adata = model.export_posterior(
        adata,
        sample_kwargs={
            "num_samples": args.posterior_num_samples,
            "batch_size": args.posterior_batch_size,
        },
    )

    factor_names = [str(name) for name in adata.uns["mod"]["factor_names"]]
    signature_columns = [f"means_per_cluster_mu_fg_{name}" for name in factor_names]
    if "means_per_cluster_mu_fg" in adata.varm:
        signatures = adata.varm["means_per_cluster_mu_fg"]
        if not isinstance(signatures, pd.DataFrame):
            signatures = pd.DataFrame(
                signatures,
                index=adata.var_names,
                columns=signature_columns,
            )
        signatures = signatures[signature_columns].copy()
    else:
        signatures = adata.var[signature_columns].copy()
    signatures.columns = factor_names

    if args.gene_key == "var_names":
        signatures.index = adata.var_names.astype(str)
    else:
        if args.gene_key not in adata.var:
            raise KeyError(f"Reference var does not contain '{args.gene_key}'")
        signatures.index = adata.var[args.gene_key].astype(str)
    if signatures.index.has_duplicates:
        duplicates = signatures.index[signatures.index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Reference gene IDs are not unique; examples: {duplicates}")

    args.output_model.parent.mkdir(parents=True, exist_ok=True)
    args.output_signatures.parent.mkdir(parents=True, exist_ok=True)
    model.save(str(args.output_model), save_anndata=True)
    signatures.to_csv(args.output_signatures, index_label="gene_id")
    print(f"Wrote reference model: {args.output_model}")
    print(f"Wrote reference signatures: {args.output_signatures}")


if __name__ == "__main__":
    main()
