#!/usr/bin/env python3
"""Reinsert joint Cell2location posterior results into per-sample SpatialData."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm


def parse_args():
    parser = argparse.ArgumentParser(
        description="Reinsert joint Cell2location posterior results into SpatialData"
    )
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--posterior-dir", required=True, type=Path)
    parser.add_argument("--posterior-manifest", required=True, type=Path)
    parser.add_argument("--training-history", required=True, type=Path)
    parser.add_argument("--diagnostics", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--table-key", default="square_016um")
    parser.add_argument("--qc-pass-key", default="qc_pass")
    return parser.parse_args()


def plot_training_history(training_history, output_path):
    history_start = max(1, int(len(training_history) * 0.2))
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    axes[0].plot(training_history["epoch"], training_history["elbo_train"])
    axes[0].set(title="Complete training history", xlabel="Epoch", ylabel="-ELBO loss")
    axes[1].plot(
        training_history["epoch"].iloc[history_start:],
        training_history["elbo_train"].iloc[history_start:],
    )
    axes[1].set(
        title="Training history after first 20%",
        xlabel="Epoch",
        ylabel="-ELBO loss",
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def plot_reconstruction_accuracy(reconstruction_qc, output_path):
    if reconstruction_qc.get("status") != "success":
        return False

    histogram = np.asarray(reconstruction_qc.get("histogram"), dtype=float)
    observed_edges = np.asarray(reconstruction_qc.get("observed_log10_edges"), dtype=float)
    expected_edges = np.asarray(reconstruction_qc.get("expected_log10_edges"), dtype=float)
    if histogram.ndim != 2:
        raise ValueError("Reconstruction histogram is not two-dimensional")
    if len(observed_edges) != histogram.shape[0] + 1:
        raise ValueError("Observed reconstruction edges do not match the histogram")
    if len(expected_edges) != histogram.shape[1] + 1:
        raise ValueError("Expected reconstruction edges do not match the histogram")
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
    axis.plot([lower, upper], [lower, upper], color="black", linestyle="--", linewidth=0.8)
    axis.set(
        title="Joint-model reconstruction accuracy",
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


def main():
    args = parse_args()

    import anndata as ad
    import spatialdata as sd

    for path, label in (
        (args.samples, "Sample manifest"),
        (args.posterior_manifest, "Posterior manifest"),
        (args.training_history, "Training history"),
        (args.diagnostics, "Joint diagnostics"),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{label} does not exist: {path}")
    if not args.posterior_dir.is_dir():
        raise FileNotFoundError(f"Posterior directory does not exist: {args.posterior_dir}")

    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    if samples.columns.tolist() != ["sample_id", "spatialdata"]:
        raise ValueError("--samples must contain exactly: sample_id, spatialdata")
    if samples["sample_id"].duplicated().any():
        raise ValueError("--samples contains duplicate sample IDs")
    samples = samples.sort_values("sample_id", kind="stable").reset_index(drop=True)

    posterior_manifest = pd.read_csv(
        args.posterior_manifest,
        sep="\t",
        dtype={"sample_id": str, "posterior_file": str},
    )
    required_manifest_columns = {
        "sample_id",
        "posterior_file",
        "n_total_locations",
        "n_modelled_locations",
        "n_input_genes",
    }
    if set(posterior_manifest.columns) != required_manifest_columns:
        raise ValueError(
            "Posterior manifest columns differ from the expected contract: "
            + ", ".join(sorted(required_manifest_columns))
        )
    if posterior_manifest["sample_id"].duplicated().any():
        raise ValueError("Posterior manifest contains duplicate sample IDs")
    posterior_manifest["sample_id"] = posterior_manifest["sample_id"].astype(str)
    posterior_manifest = posterior_manifest.set_index("sample_id")

    sample_ids = samples["sample_id"].tolist()
    if set(sample_ids) != set(posterior_manifest.index):
        missing = sorted(set(sample_ids) - set(posterior_manifest.index))
        extra = sorted(set(posterior_manifest.index) - set(sample_ids))
        raise ValueError(f"Posterior/sample mismatch; missing={missing}, extra={extra}")

    diagnostics = json.loads(args.diagnostics.read_text(encoding="utf-8"))
    if diagnostics.get("sample_ids") != sample_ids:
        raise ValueError(
            "Joint diagnostics sample order does not match the canonical reinsertion order"
        )
    if diagnostics.get("table_key") != args.table_key:
        raise ValueError(
            f"Joint model used table '{diagnostics.get('table_key')}', "
            f"but reinsertion requested '{args.table_key}'"
        )

    training_history = pd.read_csv(args.training_history)
    if training_history.columns.tolist() != ["epoch", "elbo_train"]:
        raise ValueError("Training history must contain epoch and elbo_train columns")
    if training_history.empty or not np.isfinite(training_history["elbo_train"]).all():
        raise ValueError("Training history is empty or contains non-finite ELBO values")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_training_history(
        training_history,
        args.output_dir / "cell2location_training_history.png",
    )

    reconstruction_qc = diagnostics.get("reconstruction_qc", {})
    reconstruction_status = "failed"
    reconstruction_error = str(reconstruction_qc.get("error", ""))
    try:
        if plot_reconstruction_accuracy(
            reconstruction_qc,
            args.output_dir / "cell2location_reconstruction_accuracy.png",
        ):
            reconstruction_status = "success"
        else:
            if not reconstruction_error:
                reconstruction_error = "Reconstruction QC data are unavailable"
    except Exception as error:
        reconstruction_error = f"{type(error).__name__}: {error}"
    if reconstruction_status != "success":
        (args.output_dir / "cell2location_reconstruction_accuracy_error.txt").write_text(
            reconstruction_error + "\n",
            encoding="utf-8",
        )

    pd.DataFrame(
        [
            {
                "samples": int(diagnostics["n_samples"]),
                "modelled_locations": int(diagnostics["n_modelled_locations"]),
                "shared_genes": int(diagnostics["n_shared_genes"]),
                "cell_types": int(diagnostics["n_cell_types"]),
                "epochs_completed": int(diagnostics["epochs_completed"]),
                "initial_elbo": float(diagnostics["initial_elbo"]),
                "final_elbo": float(diagnostics["final_elbo"]),
                "reconstruction_plot": reconstruction_status,
                "reconstruction_error": reconstruction_error,
            }
        ]
    ).to_csv(
        args.output_dir / "cell2location_joint_diagnostics.tsv",
        sep="\t",
        index=False,
    )

    abundance_keys = [str(value) for value in diagnostics["abundance_keys"]]
    detection_keys = [str(value) for value in diagnostics["detection_sensitivity_keys"]]
    cell_types = [str(value) for value in diagnostics["cell_types"]]

    for row in samples.itertuples(index=False):
        sample_id = row.sample_id
        input_zarr = Path(row.spatialdata)
        if not input_zarr.exists():
            raise FileNotFoundError(
                f"SpatialData for sample '{sample_id}' does not exist: {input_zarr}"
            )
        posterior_path = args.posterior_dir / posterior_manifest.loc[
            sample_id, "posterior_file"
        ]
        if not posterior_path.is_file():
            raise FileNotFoundError(
                f"Posterior result for sample '{sample_id}' does not exist: {posterior_path}"
            )

        print(f"Reinserting Cell2location posterior into {sample_id}: {input_zarr}")
        sdata = sd.read_zarr(input_zarr)
        if args.table_key not in sdata.tables:
            raise KeyError(f"Sample '{sample_id}' does not contain table '{args.table_key}'")
        table = sdata.tables[args.table_key]
        if args.qc_pass_key not in table.obs:
            raise KeyError(f"Sample '{sample_id}' obs does not contain '{args.qc_pass_key}'")
        if not table.obs_names.is_unique:
            raise ValueError(f"Sample '{sample_id}' has duplicate observation IDs")

        qc_pass = table.obs[args.qc_pass_key].fillna(False).astype(bool)
        posterior = ad.read_h5ad(posterior_path)
        if not posterior.obs_names.is_unique:
            raise ValueError(f"Posterior for sample '{sample_id}' has duplicate observations")
        posterior_sample_ids = posterior.obs["_c2l_sample_id"].astype(str).unique().tolist()
        if posterior_sample_ids != [sample_id]:
            raise ValueError(
                f"Posterior '{posterior_path.name}' contains sample labels {posterior_sample_ids}"
            )
        expected_obs = pd.Index(table.obs_names[qc_pass].astype(str))
        if set(posterior.obs_names) != set(expected_obs):
            missing = expected_obs.difference(posterior.obs_names)[:5].tolist()
            extra = posterior.obs_names.difference(expected_obs)[:5].tolist()
            raise ValueError(
                f"Posterior locations do not match QC-pass locations for '{sample_id}'; "
                f"missing examples={missing}, extra examples={extra}"
            )
        posterior = posterior[expected_obs].copy()

        full_abundances = {}
        for abundance_key in abundance_keys:
            if abundance_key not in posterior.obsm:
                raise KeyError(
                    f"Posterior for '{sample_id}' is missing obsm['{abundance_key}']"
                )
            values = posterior.obsm[abundance_key]
            if isinstance(values, pd.DataFrame):
                values = values.loc[expected_obs, cell_types].to_numpy()
            else:
                values = np.asarray(values)
            if values.shape != (len(expected_obs), len(cell_types)):
                raise ValueError(
                    f"Unexpected shape for '{sample_id}' {abundance_key}: {values.shape}"
                )
            full = pd.DataFrame(
                np.full((table.n_obs, len(cell_types)), np.nan, dtype=np.float32),
                index=table.obs_names,
                columns=cell_types,
            )
            full.loc[expected_obs, :] = values
            table.obsm[abundance_key] = full
            full_abundances[abundance_key] = full

        for detection_key in detection_keys:
            if detection_key not in posterior.obs:
                raise KeyError(
                    f"Posterior for '{sample_id}' is missing obs['{detection_key}']"
                )
            full_values = np.full(table.n_obs, np.nan, dtype=np.float32)
            full_values[qc_pass.to_numpy()] = posterior.obs[detection_key].to_numpy(
                dtype=np.float32
            )
            table.obs[detection_key] = full_values

        manifest_row = posterior_manifest.loc[sample_id]
        if int(manifest_row["n_total_locations"]) != table.n_obs:
            raise ValueError(f"Total location count changed for sample '{sample_id}'")
        if int(manifest_row["n_modelled_locations"]) != int(qc_pass.sum()):
            raise ValueError(f"QC-pass location count changed for sample '{sample_id}'")

        sample_metadata = dict(diagnostics)
        sample_metadata.update(
            {
                "sample_id": sample_id,
                "n_joint_modelled_locations": int(diagnostics["n_modelled_locations"]),
                "n_total_locations": int(table.n_obs),
                "n_modelled_locations": int(qc_pass.sum()),
                "n_input_genes": int(manifest_row["n_input_genes"]),
                "training_history": {
                    "epoch": training_history["epoch"].to_numpy(),
                    "elbo_train": training_history["elbo_train"].to_numpy(),
                },
            }
        )
        table.uns["cell2location"] = sample_metadata
        sdata.tables[args.table_key] = table

        output_zarr = args.output_dir / f"{sample_id}_cell2location.zarr"
        output_mean = (
            args.output_dir / f"{sample_id}_cell2location_abundance_mean.csv"
        )
        output_q05 = args.output_dir / f"{sample_id}_cell2location_abundance_q05.csv"
        sdata.write(output_zarr)
        full_abundances["means_cell_abundance_w_sf"].to_csv(
            output_mean,
            index_label="location",
        )
        full_abundances["q05_cell_abundance_w_sf"].to_csv(
            output_q05,
            index_label="location",
        )
        print(f"Wrote deconvolved SpatialData: {output_zarr}")


if __name__ == "__main__":
    main()
