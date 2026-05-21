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


def read_memmapped_2d_tiff(path, label):
    """Read a TIFF as a squeezed 2D memmap."""
    try:
        image = tifffile.memmap(path)
    except Exception as exc:  # noqa: BLE001 - preserve the memmap failure context
        raise ValueError(f"{label} TIFF is not memmappable: {path}") from exc

    image = image.squeeze()
    if image.ndim != 2:
        raise ValueError(f"Expected {label} TIFF to be 2D after squeeze, got shape {image.shape}")

    return image


def downsampled_size(size, factor):
    """Return the size produced by slicing an axis with [::factor]."""
    return 0 if size == 0 else ((size - 1) // factor) + 1


def infer_downsample_factor(full_shape, downsampled_shape):
    """Infer the integer stride that maps a full-resolution image to a mask."""
    if full_shape == downsampled_shape:
        return 1

    height_ratio = full_shape[0] / downsampled_shape[0]
    width_ratio = full_shape[1] / downsampled_shape[1]
    first_candidate = max(1, int(min(height_ratio, width_ratio)) - 2)
    last_candidate = int(max(height_ratio, width_ratio)) + 3

    for factor in range(first_candidate, last_candidate + 1):
        if tuple(downsampled_size(size, factor) for size in full_shape) == downsampled_shape:
            return factor

    raise ValueError(
        "Segmentation mask shape is not compatible with source image shape; "
        f"got mask shape {downsampled_shape} and source image shape {full_shape}"
    )


def read_mask_as_sparse(mask_path, output_labels_npz):
    mask = tifffile.imread(mask_path)
    if mask.ndim > 2:
        mask = np.squeeze(mask)
    if mask.ndim != 2:
        raise ValueError(f"Expected a 2D segmentation mask, got shape {mask.shape}")

    labels_npz = scipy.sparse.csr_matrix(mask)
    scipy.sparse.save_npz(output_labels_npz, labels_npz)
    segmentation_labels = int(np.unique(mask).size - (1 if np.any(mask == 0) else 0))
    return segmentation_labels, mask.shape


def store_labels_npz_path(adata, labels_npz_path, labels_key):
    """Store labels path in adata.uns using the same convention as Bin2Cell."""
    if "bin2cell" not in adata.uns:
        adata.uns["bin2cell"] = {}
    if "labels_npz_paths" not in adata.uns["bin2cell"]:
        adata.uns["bin2cell"]["labels_npz_paths"] = {}

    labels_npz_path = str(labels_npz_path)
    if labels_npz_path.startswith("/"):
        adata.uns["bin2cell"]["labels_npz_paths"][labels_key] = labels_npz_path
    else:
        adata.uns["bin2cell"]["labels_npz_paths"][labels_key] = str(Path.cwd() / labels_npz_path)


def insert_labels_compatible(
    adata,
    labels_npz_path,
    source_image_path,
    mask_shape,
    basis="spatial",
    spatial_key="spatial",
    labels_key="labels",
):
    """Insert labels while handling downsampled masks and SciPy sparse indexing quirks."""
    labels_sparse = scipy.sparse.load_npz(labels_npz_path)
    source_image = read_memmapped_2d_tiff(source_image_path, "source image")
    mask_downsample_factor = infer_downsample_factor(source_image.shape, mask_shape)

    store_labels_npz_path(adata, labels_npz_path, labels_key)

    coords = b2c.get_mpp_coords(adata, basis=basis, spatial_key=spatial_key, mpp=None)
    coords = np.asarray(coords, dtype=np.int64)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError(f"Expected 2D coordinates with two columns, got shape {coords.shape}")

    if mask_downsample_factor != 1:
        coords = coords // mask_downsample_factor

    adata.obs[labels_key] = 0
    mask = (
        (coords[:, 0] >= 0)
        & (coords[:, 0] < labels_sparse.shape[0])
        & (coords[:, 1] >= 0)
        & (coords[:, 1] < labels_sparse.shape[1])
    )

    if np.any(mask):
        rows = np.asarray(coords[mask, 0], dtype=np.int64).reshape(-1)
        cols = np.asarray(coords[mask, 1], dtype=np.int64).reshape(-1)
        values = labels_sparse[rows, cols]
        if scipy.sparse.issparse(values):
            values = values.toarray()
        values = np.asarray(values).reshape(-1)
        adata.obs.loc[mask, labels_key] = values.astype(np.int64, copy=False)

    return mask_downsample_factor


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


def labels_have_assignments(adata, labels_key):
    """Return True when at least one bin has a non-zero segmentation label."""
    labels = adata.obs[labels_key].to_numpy()
    return bool(np.any(labels > 0))


def run_bin2cell(args):
    binned_path, spaceranger_image_path = validate_inputs(args)
    segmentation_labels, mask_shape = read_mask_as_sparse(
        args.segmentation_mask_tif,
        args.output_labels_npz,
    )

    adata = b2c.read_visium(
        str(binned_path),
        source_image_path=str(args.source_image),
        spaceranger_image_path=str(spaceranger_image_path),
    )

    mask_downsample_factor = insert_labels_compatible(
        adata=adata,
        labels_npz_path=str(args.output_labels_npz),
        source_image_path=args.source_image,
        mask_shape=mask_shape,
        basis="spatial",
        spatial_key="spatial",
        labels_key="labels",
    )

    if labels_have_assignments(adata, "labels"):
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
    else:
        adata.obs["labels_expanded"] = 0
        adata_cells = adata[:0].copy()

    adata_cells.uns["bin2cell"] = {
        "sample": args.sample_name,
        "bin_size": args.bin_size,
        "volume_ratio": args.volume_ratio,
        "labels_key": "labels",
        "expanded_labels_key": "labels_expanded",
        "segmentation_labels": segmentation_labels,
        "segmentation_downsample_factor": mask_downsample_factor,
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
