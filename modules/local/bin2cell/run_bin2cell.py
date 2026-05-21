#!/usr/bin/env python3
"""Aggregate Visium HD bins into cell-level AnnData using Bin2Cell labels."""

import argparse
import importlib
import sys
from pathlib import Path

import bin2cell as b2c
import numpy as np
import pandas as pd
import scipy.sparse
import tifffile
import yaml


def package_version(package_name):
    """Return an installed package version without requiring importlib.metadata everywhere."""
    try:
        from importlib.metadata import version

        return version(package_name)
    except Exception:
        try:
            module = importlib.import_module(package_name)
            return getattr(module, "__version__", None)
        except Exception:
            return None


def write_versions(process_name):
    versions = {process_name: {}}
    for package_name in ("bin2cell", "scanpy", "numpy", "pandas", "scipy", "tifffile"):
        version = package_version(package_name)
        if version:
            versions[process_name][package_name] = version
    versions[process_name]["python"] = (
        f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    )
    yaml.safe_dump(versions, sys.stdout, sort_keys=True)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Convert a segmentation mask to sparse labels, assign Visium HD bins "
            "to labels with Bin2Cell, and write a cell-level AnnData object."
        )
    )
    parser.add_argument("--sample-name", help="Sample name used in outputs")
    parser.add_argument(
        "--spaceranger-outs",
        type=Path,
        help="Space Ranger outs directory",
    )
    parser.add_argument(
        "--segmentation-mask-tif",
        type=Path,
        help="2D segmentation mask TIFF",
    )
    parser.add_argument(
        "--source-image",
        type=Path,
        help="Image path passed to bin2cell.read_visium as source_image_path",
    )
    parser.add_argument(
        "--bin-size",
        default="square_002um",
        help="Space Ranger binned_outputs subdirectory to use",
    )
    parser.add_argument(
        "--volume-ratio",
        type=float,
        default=4.0,
        help="Volume ratio used by bin2cell.expand_labels",
    )
    parser.add_argument(
        "--output-h5ad",
        type=Path,
        help="Output cell-level AnnData .h5ad",
    )
    parser.add_argument(
        "--output-labels-npz",
        type=Path,
        help="Output sparse segmentation labels .npz",
    )
    parser.add_argument(
        "--output-summary",
        type=Path,
        help="Output one-row Bin2Cell assignment summary CSV",
    )
    parser.add_argument(
        "--versions-dict",
        help="Only print package versions for the provided Nextflow process name",
    )
    return parser.parse_args()


def validate_inputs(args):
    required_paths = (
        "spaceranger_outs",
        "segmentation_mask_tif",
        "source_image",
        "output_h5ad",
        "output_labels_npz",
        "output_summary",
    )
    missing = [field for field in required_paths if getattr(args, field) is None]
    if args.sample_name is None:
        missing.append("sample_name")
    if missing:
        raise ValueError(f"Missing required arguments: {', '.join(missing)}")

    binned_path = args.spaceranger_outs / "binned_outputs" / args.bin_size
    spaceranger_image_path = binned_path / "spatial"

    if not args.spaceranger_outs.exists():
        raise FileNotFoundError(f"Space Ranger outs directory does not exist: {args.spaceranger_outs}")
    if not binned_path.exists():
        raise FileNotFoundError(f"Space Ranger binned output does not exist: {binned_path}")
    if not spaceranger_image_path.exists():
        raise FileNotFoundError(f"Space Ranger spatial image directory does not exist: {spaceranger_image_path}")
    if not args.segmentation_mask_tif.exists():
        raise FileNotFoundError(f"Segmentation mask TIFF does not exist: {args.segmentation_mask_tif}")
    if not args.source_image.exists():
        raise FileNotFoundError(f"Source image does not exist: {args.source_image}")

    return binned_path, spaceranger_image_path


def read_mask_as_sparse(mask_path, output_labels_npz):
    mask = tifffile.imread(mask_path)
    if mask.ndim > 2:
        mask = np.squeeze(mask)
    if mask.ndim != 2:
        raise ValueError(f"Expected a 2D segmentation mask, got shape {mask.shape}")

    labels_npz = scipy.sparse.csr_matrix(mask)
    scipy.sparse.save_npz(output_labels_npz, labels_npz)
    segmentation_labels = int(np.unique(mask).size - (1 if np.any(mask == 0) else 0))
    return segmentation_labels


def write_summary(args, adata, adata_cells, segmentation_labels):
    labels = adata.obs["labels_expanded"].to_numpy()
    input_bins = int(adata.n_obs)
    assigned_bins = int(np.sum(labels > 0))
    unassigned_bins = int(np.sum(labels == 0))

    row = {
        "sample": args.sample_name,
        "input_bins": input_bins,
        "assigned_bins": assigned_bins,
        "unassigned_bins": unassigned_bins,
        "assigned_fraction": assigned_bins / input_bins if input_bins else 0.0,
        "unassigned_fraction": unassigned_bins / input_bins if input_bins else 0.0,
        "segmentation_labels": segmentation_labels,
        "output_cells": int(adata_cells.n_obs),
        "volume_ratio": args.volume_ratio,
        "bin_size": args.bin_size,
    }
    pd.DataFrame([row]).to_csv(args.output_summary, index=False)


def run_bin2cell(args):
    binned_path, spaceranger_image_path = validate_inputs(args)
    segmentation_labels = read_mask_as_sparse(args.segmentation_mask_tif, args.output_labels_npz)

    adata = b2c.read_visium(
        str(binned_path),
        source_image_path=str(args.source_image),
        spaceranger_image_path=str(spaceranger_image_path),
    )

    b2c.insert_labels(
        adata=adata,
        labels_npz_path=str(args.output_labels_npz),
        basis="spatial",
        spatial_key="spatial",
        mpp=None,
        labels_key="labels",
    )

    b2c.expand_labels(
        adata,
        labels_key="labels",
        expanded_labels_key="labels_expanded",
        volume_ratio=args.volume_ratio,
    )

    adata_cells = b2c.bin_to_cell(
        adata=adata,
        labels_key="labels_expanded",
        spatial_keys=["spatial"],
        diameter_scale_factor=None,
    )

    adata_cells.uns["bin2cell"] = {
        "sample": args.sample_name,
        "bin_size": args.bin_size,
        "volume_ratio": args.volume_ratio,
        "labels_key": "labels",
        "expanded_labels_key": "labels_expanded",
        "segmentation_labels": segmentation_labels,
    }
    adata_cells.write_h5ad(args.output_h5ad)
    write_summary(args, adata, adata_cells, segmentation_labels)


def main():
    args = parse_args()
    if args.versions_dict:
        write_versions(args.versions_dict)
        return
    run_bin2cell(args)


if __name__ == "__main__":
    main()
