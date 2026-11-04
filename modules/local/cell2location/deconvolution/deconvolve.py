#!/usr/bin/env python3
"""Train one Cell2location spatial model jointly across multiple SpatialData samples."""

import argparse
import json
import random
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run joint Cell2location deconvolution across SpatialData samples"
    )
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--signatures", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--output-manifest", required=True, type=Path)
    parser.add_argument("--output-training-history", required=True, type=Path)
    parser.add_argument("--output-diagnostics", required=True, type=Path)
    parser.add_argument("--table-key", default="square_016um")
    parser.add_argument("--qc-pass-key", default="qc_pass")
    parser.add_argument(
        "--counts-layer",
        default="X",
        help="Spatial layer containing raw counts; use X for the main matrix",
    )
    parser.add_argument(
        "--gene-key",
        default="var_names",
        help="Spatial var column containing gene IDs, or var_names",
    )
    parser.add_argument("--max-epochs", type=int, default=400)
    parser.add_argument("--batch-size", type=int, default=25_000)
    parser.add_argument("--n-cells-per-location", type=float, default=3.0)
    parser.add_argument("--detection-alpha", type=float, default=200.0)
    parser.add_argument("--posterior-num-samples", type=int, default=1000)
    parser.add_argument("--posterior-batch-size", type=int, default=2048)
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

    for path, label in (
        (args.samples, "Sample manifest"),
        (args.signatures, "Reference signatures"),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{label} does not exist: {path}")

    positive_parameters = {
        "--max-epochs": args.max_epochs,
        "--batch-size": args.batch_size,
        "--n-cells-per-location": args.n_cells_per_location,
        "--detection-alpha": args.detection_alpha,
        "--posterior-num-samples": args.posterior_num_samples,
        "--posterior-batch-size": args.posterior_batch_size,
        "--qc-num-locations": args.qc_num_locations,
        "--minimum-shared-genes": args.minimum_shared_genes,
    }
    for name, value in positive_parameters.items():
        if value <= 0:
            raise ValueError(f"{name} must be greater than zero")

    import anndata as ad
    import cell2location
    import numpy as np
    import pandas as pd
    import spatialdata as sd
    import torch
    from cell2location.models import Cell2location
    from scipy.sparse import issparse
    from scvi import REGISTRY_KEYS

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

    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    if samples.columns.tolist() != ["sample_id", "spatialdata"]:
        raise ValueError("--samples must contain exactly: sample_id, spatialdata")
    if samples.empty:
        raise ValueError("--samples does not contain any samples")
    if samples["sample_id"].isna().any() or samples["spatialdata"].isna().any():
        raise ValueError("--samples contains missing sample IDs or SpatialData paths")
    if samples["sample_id"].duplicated().any():
        duplicates = samples.loc[samples["sample_id"].duplicated(), "sample_id"].tolist()
        raise ValueError(f"Duplicate sample IDs: {duplicates}")
    for sample_id in samples["sample_id"]:
        if not sample_id or Path(sample_id).name != sample_id or "\t" in sample_id:
            raise ValueError(f"Sample ID is not safe for use as a filename: {sample_id!r}")
    samples = samples.sort_values("sample_id", kind="stable").reset_index(drop=True)

    signatures = pd.read_csv(args.signatures, index_col=0)
    signatures.index = signatures.index.astype(str)
    signatures.columns = signatures.columns.astype(str)
    if signatures.index.has_duplicates:
        duplicates = signatures.index[signatures.index.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Signature gene IDs are not unique; examples: {duplicates}")
    if signatures.columns.has_duplicates:
        duplicates = signatures.columns[signatures.columns.duplicated()].unique()[:5].tolist()
        raise ValueError(f"Signature cell-type names are not unique; examples: {duplicates}")

    spatial_adatas = []
    sample_metrics = []
    shared_gene_set = set(signatures.index)
    sample_key = "_c2l_sample_id"
    original_obs_key = "_c2l_original_obs_name"

    for row in samples.itertuples(index=False):
        input_zarr = Path(row.spatialdata)
        if not input_zarr.exists():
            raise FileNotFoundError(
                f"SpatialData for sample '{row.sample_id}' does not exist: {input_zarr}"
            )
        print(f"Reading SpatialData for {row.sample_id}: {input_zarr}")
        sdata = sd.read_zarr(input_zarr)
        if args.table_key not in sdata.tables:
            raise KeyError(
                f"Sample '{row.sample_id}' does not contain table '{args.table_key}'"
            )
        table = sdata.tables[args.table_key]
        if not table.obs_names.is_unique:
            duplicates = table.obs_names[table.obs_names.duplicated()].unique()[:5].tolist()
            raise ValueError(
                f"Sample '{row.sample_id}' has duplicate observation IDs; examples: {duplicates}"
            )
        if args.qc_pass_key not in table.obs:
            raise KeyError(
                f"Sample '{row.sample_id}' obs does not contain '{args.qc_pass_key}'"
            )
        if args.counts_layer != "X" and args.counts_layer not in table.layers:
            raise KeyError(
                f"Sample '{row.sample_id}' layers do not contain '{args.counts_layer}'"
            )

        qc_pass = table.obs[args.qc_pass_key].fillna(False).astype(bool)
        if not qc_pass.any():
            raise ValueError(
                f"Sample '{row.sample_id}' has no locations with {args.qc_pass_key}=True"
            )

        if args.gene_key == "var_names":
            spatial_gene_ids = pd.Index(table.var_names.astype(str))
        else:
            if args.gene_key not in table.var:
                raise KeyError(
                    f"Sample '{row.sample_id}' var does not contain '{args.gene_key}'"
                )
            spatial_gene_ids = pd.Index(table.var[args.gene_key].astype(str))
        if spatial_gene_ids.has_duplicates:
            duplicates = spatial_gene_ids[spatial_gene_ids.duplicated()].unique()[:5].tolist()
            raise ValueError(
                f"Sample '{row.sample_id}' gene IDs are not unique; examples: {duplicates}"
            )
        selected_overlap = len(spatial_gene_ids.intersection(signatures.index, sort=False))
        alternative_overlaps = {}
        for alternative_key in ("var_names", "gene_ids"):
            if alternative_key == args.gene_key:
                continue
            if alternative_key == "var_names":
                alternative_ids = pd.Index(table.var_names.astype(str))
            elif alternative_key in table.var:
                alternative_ids = pd.Index(table.var[alternative_key].astype(str))
            else:
                continue
            alternative_overlaps[alternative_key] = len(
                alternative_ids.intersection(signatures.index, sort=False)
            )
        print(f"Spatial gene ID examples for {row.sample_id}: {spatial_gene_ids[:5].tolist()}")
        print(
            f"Reference overlap for {row.sample_id}: {selected_overlap} genes"
            + (
                "; alternatives: "
                + ", ".join(
                    f"{key}={value}" for key, value in alternative_overlaps.items()
                )
                if alternative_overlaps
                else ""
            )
        )

        counts = table.X if args.counts_layer == "X" else table.layers[args.counts_layer]
        counts = counts[qc_pass.to_numpy(), :].copy()
        obs_names = table.obs_names[qc_pass].astype(str)
        obs = pd.DataFrame(index=obs_names)
        obs[original_obs_key] = obs_names.to_numpy()
        obs[sample_key] = row.sample_id
        obs.index = pd.Index([f"{row.sample_id}::{name}" for name in obs_names], dtype=str)
        sample_adata = ad.AnnData(
            X=counts,
            obs=obs,
            var=pd.DataFrame(index=spatial_gene_ids),
        )
        spatial_adatas.append(sample_adata)
        shared_gene_set.intersection_update(spatial_gene_ids)
        sample_metrics.append(
            {
                "sample_id": row.sample_id,
                "n_total_locations": int(table.n_obs),
                "n_modelled_locations": int(qc_pass.sum()),
                "n_input_genes": int(table.n_vars),
                "selected_gene_overlap": int(selected_overlap),
                "alternative_gene_overlaps": alternative_overlaps,
            }
        )
        print(
            f"Locations passing QC for {row.sample_id}: "
            f"{int(qc_pass.sum())} / {table.n_obs}"
        )
        del sdata, table, counts

    shared_genes = pd.Index(
        [gene for gene in signatures.index if gene in shared_gene_set],
        dtype=str,
    )
    print(f"Spatial gene key: {args.gene_key}")
    print(f"Reference gene ID examples: {signatures.index[:5].tolist()}")
    print(f"Joint shared genes: {len(shared_genes)}")
    if len(shared_genes) < args.minimum_shared_genes:
        overlap_summary = ", ".join(
            f"{row['sample_id']}={row['selected_gene_overlap']}"
            for row in sample_metrics
        )
        alternative_summary = "; ".join(
            f"{row['sample_id']}: "
            + ", ".join(
                f"{key}={value}"
                for key, value in row["alternative_gene_overlaps"].items()
            )
            for row in sample_metrics
            if row["alternative_gene_overlaps"]
        )
        raise ValueError(
            f"Only {len(shared_genes)} genes are shared by the reference and every spatial "
            f"sample; expected at least {args.minimum_shared_genes}. Check that all inputs "
            f"use the same gene identifier namespace. Per-sample overlaps: {overlap_summary}."
            + (
                f" Alternative gene-key overlaps: {alternative_summary}."
                if alternative_summary
                else ""
            )
        )

    spatial_adatas = [adata[:, shared_genes].copy() for adata in spatial_adatas]
    adata_joint = ad.concat(
        spatial_adatas,
        axis=0,
        join="inner",
        merge="same",
        index_unique=None,
    )
    if not adata_joint.obs_names.is_unique:
        raise ValueError("Joint observation IDs are not unique after sample prefixing")
    signatures = signatures.loc[shared_genes].copy()
    sample_ids = samples["sample_id"].tolist()

    Cell2location.setup_anndata(adata=adata_joint, batch_key=sample_key)
    model = Cell2location(
        adata_joint,
        cell_state_df=signatures,
        N_cells_per_location=args.n_cells_per_location,
        detection_alpha=args.detection_alpha,
    )

    print(
        f"Training one spatial model on {len(sample_ids)} samples, "
        f"{adata_joint.n_obs} locations and {adata_joint.n_vars} genes on {accelerator}"
    )
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
    args.output_training_history.parent.mkdir(parents=True, exist_ok=True)
    training_history.to_csv(args.output_training_history, index=False)

    adata_joint = model.export_posterior(
        adata_joint,
        sample_kwargs={
            "num_samples": args.posterior_num_samples,
            "batch_size": args.posterior_batch_size,
            "summary_frequency": 1,
            "accelerator": accelerator,
        },
        add_to_obsm=["means", "stds", "q05", "q95"],
    )

    cell_types = signatures.columns.tolist()
    abundance_keys = [
        "means_cell_abundance_w_sf",
        "stds_cell_abundance_w_sf",
        "q05_cell_abundance_w_sf",
        "q95_cell_abundance_w_sf",
    ]
    for abundance_key in abundance_keys:
        if abundance_key not in adata_joint.obsm:
            raise KeyError(f"Cell2location did not create adata.obsm['{abundance_key}']")
        values = np.asarray(adata_joint.obsm[abundance_key])
        if values.shape != (adata_joint.n_obs, len(cell_types)):
            raise ValueError(
                f"Unexpected shape for adata.obsm['{abundance_key}']: {values.shape}, "
                f"expected {(adata_joint.n_obs, len(cell_types))}"
            )

    posterior_metadata = adata_joint.uns.get("mod", {})
    detection_values = {}
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
        if values.size != adata_joint.n_obs:
            raise ValueError(f"Unexpected detection_y_s size for {summary_name}: {values.size}")
        detection_values[f"c2l_detection_sensitivity_{output_suffix}"] = values

    reconstruction_qc = {"status": "success", "error": ""}
    try:
        random_generator = np.random.default_rng(args.seed)
        reconstruction_indices = random_generator.choice(
            adata_joint.n_obs,
            min(args.qc_num_locations, adata_joint.n_obs),
            replace=False,
        )
        expected_parameters = model.module.model.compute_expected(
            model.samples["post_sample_means"],
            model.adata_manager,
            ind_x=reconstruction_indices,
        )
        observed_counts = model.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)[
            reconstruction_indices, :
        ]
        if issparse(observed_counts):
            observed_counts = observed_counts.toarray()
        observed_counts = np.asarray(observed_counts)
        expected_counts = np.asarray(expected_parameters["mu"])
        if observed_counts.shape != expected_counts.shape:
            raise ValueError(
                "Observed and expected reconstruction matrices have different shapes: "
                f"{observed_counts.shape} and {expected_counts.shape}"
            )
        if not np.isfinite(observed_counts).all() or not np.isfinite(expected_counts).all():
            raise ValueError("Observed or expected reconstruction values are not finite")
        reconstruction_histogram, observed_edges, expected_edges = np.histogram2d(
            np.log10(observed_counts.ravel() + 1),
            np.log10(expected_counts.ravel() + 1),
            bins=50,
        )
        reconstruction_qc.update(
            {
                "histogram": reconstruction_histogram.astype(np.int64).tolist(),
                "observed_log10_edges": observed_edges.tolist(),
                "expected_log10_edges": expected_edges.tolist(),
                "n_locations": int(len(reconstruction_indices)),
            }
        )
    except Exception as error:  # Preserve successful inference if optional QC fails.
        reconstruction_qc.update(
            {
                "status": "failed",
                "error": f"{type(error).__name__}: {error}".replace("\t", " ").replace(
                    "\n", " "
                ),
            }
        )

    diagnostics = {
        "model_name": posterior_metadata.get("model_name", "Cell2location"),
        "cell2location_version": getattr(cell2location, "__version__", "unknown"),
        "joint_model": True,
        "sample_key": sample_key,
        "sample_ids": sample_ids,
        "table_key": args.table_key,
        "qc_pass_key": args.qc_pass_key,
        "counts_layer": args.counts_layer,
        "gene_key": args.gene_key,
        "abundance_keys": abundance_keys,
        "detection_sensitivity_keys": list(detection_values),
        "cell_types": cell_types,
        "n_samples": len(sample_ids),
        "n_modelled_locations": int(adata_joint.n_obs),
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
        "final_elbo": float(training_history["elbo_train"].iloc[-1]),
        "reconstruction_qc": reconstruction_qc,
    }
    args.output_diagnostics.parent.mkdir(parents=True, exist_ok=True)
    args.output_diagnostics.write_text(
        json.dumps(diagnostics, indent=2) + "\n",
        encoding="utf-8",
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_rows = []
    metrics_by_sample = {row["sample_id"]: row for row in sample_metrics}
    for sample_id in sample_ids:
        sample_mask = adata_joint.obs[sample_key].astype(str).to_numpy() == sample_id
        sample_obs = adata_joint.obs.loc[
            sample_mask,
            [original_obs_key, sample_key],
        ].copy()
        sample_obs.index = pd.Index(sample_obs[original_obs_key].astype(str), dtype=str)
        posterior = ad.AnnData(obs=sample_obs)
        for abundance_key in abundance_keys:
            posterior.obsm[abundance_key] = pd.DataFrame(
                np.asarray(adata_joint.obsm[abundance_key])[sample_mask].astype(np.float32),
                index=posterior.obs_names,
                columns=cell_types,
            )
        for output_key, values in detection_values.items():
            posterior.obs[output_key] = values[sample_mask]
        posterior.uns["cell2location"] = {
            "joint_model": True,
            "sample_id": sample_id,
            "sample_key": sample_key,
            "cell_types": np.asarray(cell_types, dtype=str),
            "abundance_keys": np.asarray(abundance_keys, dtype=str),
            "detection_sensitivity_keys": np.asarray(list(detection_values), dtype=str),
        }
        posterior_name = f"{sample_id}_cell2location_posterior.h5ad"
        posterior.write_h5ad(args.output_dir / posterior_name)
        metrics = metrics_by_sample[sample_id]
        manifest_rows.append(
            {
                "sample_id": sample_id,
                "posterior_file": posterior_name,
                "n_total_locations": metrics["n_total_locations"],
                "n_modelled_locations": metrics["n_modelled_locations"],
                "n_input_genes": metrics["n_input_genes"],
            }
        )
        print(f"Wrote posterior result for {sample_id}: {posterior_name}")

    pd.DataFrame(manifest_rows).to_csv(args.output_manifest, sep="\t", index=False)
    print(f"Wrote joint posterior manifest: {args.output_manifest}")
    print(f"Wrote joint training history: {args.output_training_history}")
    print(f"Wrote joint diagnostics: {args.output_diagnostics}")


if __name__ == "__main__":
    main()
