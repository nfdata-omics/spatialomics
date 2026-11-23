#!/usr/bin/env python3
"""Aggregate raw Visium HD counts by ROI and combine annotated samples."""

import argparse
import gzip
import json
from pathlib import Path
import re

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from shapely import STRtree, points


NULLABLE_ANNOTATIONS = {"slide", "area", "classification_name"}


def require_unique(index, label):
    if not index.is_unique:
        duplicates = index[index.duplicated()].unique()
        raise ValueError(f"Duplicate {label}: {', '.join(map(str, duplicates[:5]))}")


def read_bin_table(path, downsample_factor):
    table = ad.read_zarr(path)
    require_unique(table.obs_names, f"bin barcodes in {path}")
    # The Visium HD reader stores full-resolution x and y in table row order.
    xy = np.asarray(table.obsm["spatial"], dtype=np.float64)
    if xy.shape != (table.n_obs, 2) or not np.isfinite(xy).all():
        raise ValueError(f"Expected finite bin centres in obsm['spatial'] in {path}")
    is_gene = table.var["feature_types"].eq("Gene Expression").to_numpy()
    gene_ids = pd.Index(table.var.loc[is_gene, "gene_ids"], name="gene_id")
    require_unique(gene_ids, f"gene IDs in {path}")
    var = pd.DataFrame({"gene_symbol": table.var_names[is_gene]}, index=gene_ids)
    counts = sparse.csr_matrix(table.X)
    if not is_gene.all():
        counts = counts[:, is_gene]
    return counts, xy / downsample_factor, var


def aggregate_counts(geometries, xy, counts):
    tree = STRtree(geometries)
    n_rois = len(geometries)
    sums = sparse.csr_matrix((n_rois, counts.shape[1]), dtype=np.int64)
    included = np.zeros(n_rois, dtype=np.int64)
    shared = np.zeros(n_rois, dtype=np.int64)

    # Bound temporary geometry and membership arrays on fine bin grids.
    chunk_size = 250_000
    for start in range(0, len(xy), chunk_size):
        stop = min(start + chunk_size, len(xy))
        bins, rois = tree.query(points(xy[start:stop]), predicate="covered_by")
        if not len(bins):
            continue
        membership = sparse.csr_matrix(
            (np.ones(len(bins), dtype=np.int64), (rois, bins)),
            shape=(n_rois, stop - start),
        )
        chunk = counts[start:stop]
        # The Space Ranger reader can store integral UMIs in a floating-point dtype.
        values = chunk.data
        if (
            not np.isfinite(values).all()
            or values.min(initial=0) < 0
            or values.max(initial=0).item() >= 2**63
            or (values.dtype.kind == "f" and not np.equal(values, np.floor(values)).all())
        ):
            raise ValueError("Expected non-negative integer UMI counts representable as int64 in the input table")
        sums += membership @ chunk.astype(np.int64, copy=False)
        included += np.bincount(rois, minlength=n_rois)
        shared += np.bincount(rois[np.bincount(bins)[bins] > 1], minlength=n_rois)
    return sums, included, shared


def aggregate(args, sample):
    import spatialdata as sd

    resolution = re.fullmatch(r"square_((?!000)[0-9]{3})um", args.bin_label)
    if not resolution:
        raise ValueError(f"Invalid bin label {args.bin_label!r}. Use a Space Ranger bin label such as square_008um")
    table_path = Path("tables") / args.bin_label
    if not (args.zarr / table_path).is_dir():
        raise ValueError(f"Missing SpatialData table {args.bin_label!r} in {args.zarr}")

    sdata = sd.read_zarr(args.zarr, selection=("shapes",))
    downsample_factor = sdata.attrs["roi_zarr_downsample_factor"]
    roi_key = args.zarr.stem
    shapes = sdata.shapes[roi_key]
    roi_ids = pd.Index(shapes["roi_id"], dtype=str)

    counts, xy, var = read_bin_table(args.zarr / table_path, downsample_factor)
    sums, included, shared = aggregate_counts(shapes.geometry, xy, counts)

    sample_id = sample["sample_id"]
    obs_names = [json.dumps([sample_id, roi_id], separators=(",", ":")) for roi_id in roi_ids]
    obs = pd.DataFrame(
        {
            "sample_id": sample_id,
            "roi_id": roi_ids,
            "classification_name": pd.array(shapes.get("classification_name", [None] * len(shapes)), dtype="string"),
            "n_bins_included": included,
            "n_bins_shared": shared,
            "total_counts": np.asarray(sums.sum(axis=1)).ravel(),
            "roi_status": np.where(included > 0, "ok", "no_bins"),
        },
        index=obs_names,
    )
    for key in ("slide", "area"):
        if key in sample:
            obs[key] = pd.array([sample[key]] * len(obs), dtype="string")
    for key, value in sdata.attrs.get("roi_sample_metadata", {}).items():
        if value is not None and not isinstance(value, str):
            value = json.dumps(value, sort_keys=True, separators=(",", ":"))
        # Escape path separators for HDF5 while keeping ordinary metadata keys readable.
        column = "sample_metadata." + key.replace("%", "%25").replace("/", "%2F")
        obs[column] = pd.array([value] * len(obs), dtype="string")
    uns = {
        "source_matrix": table_path.as_posix(),
        "bin_size_um": int(resolution[1]),
        "assignment_rule": "Assign each bin retained by Space Ranger to every ROI covering its centre",
        "samples": {
            sample_id: {
                "annotated_zarr": args.zarr.name,
                "roi_shape_key": roi_key,
                "zarr_downsample_factor": downsample_factor,
            }
        },
    }
    ad.AnnData(X=sums, obs=obs, var=var, uns=uns).write_h5ad(args.output)


def combine(args):
    samples = {}
    for path in sorted(args.inputs):
        data = ad.read_h5ad(path)
        sample_id = data.obs["sample_id"].iloc[0]
        if sample_id in samples:
            raise ValueError(f"Duplicate pseudobulk result for sample {sample_id!r}: {path}")
        samples[sample_id] = data

    reference_id, reference = next(iter(samples.items()))
    for sample_id, data in samples.items():
        missing = reference.var_names.difference(data.var_names)
        extra = data.var_names.difference(reference.var_names)
        if len(missing) or len(extra):
            raise ValueError(
                f"Gene IDs differ for sample {sample_id!r} relative to {reference_id!r}. "
                f"Missing ({len(missing)}): {', '.join(missing[:5])}. "
                f"Extra ({len(extra)}): {', '.join(extra[:5])}."
            )

    # Outer join retains all observation columns. Gene sets were checked above.
    combined = ad.concat(list(samples.values()), join="outer", merge="first", uns_merge="unique")
    for key in combined.obs:
        if key in NULLABLE_ANNOTATIONS or key.startswith("sample_metadata."):
            combined.obs[key] = combined.obs[key].astype("string")
    combined.write_h5ad(args.output)
    if args.export_tsv:
        export_tsv(combined, args.output)


def export_tsv(data, output):
    prefix = output.with_suffix("")
    data.obs.to_csv(f"{prefix}_rois.tsv", sep="\t", index_label="obs_id")
    data.var.to_csv(f"{prefix}_genes.tsv", sep="\t", index_label="gene_id")

    # Write a wide matrix without densifying all ROIs at once.
    with gzip.open(f"{prefix}_counts.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        for start in range(0, data.n_obs, 128):
            stop = start + 128
            counts = data.X[start:stop]
            if sparse.issparse(counts):
                counts = counts.toarray()
            frame = pd.DataFrame(counts, index=data.obs_names[start:stop], columns=data.var_names)
            frame.to_csv(handle, sep="\t", index_label="obs_id", header=start == 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    per_sample = commands.add_parser("aggregate")
    for name in ("zarr", "sample", "output"):
        per_sample.add_argument(f"--{name}", type=Path, required=True)
    per_sample.add_argument("--bin-label", required=True)
    merge = commands.add_parser("combine")
    merge.add_argument("--output", type=Path, required=True)
    merge.add_argument("--export-tsv", action="store_true", help="Also write count, ROI, and gene tables alongside the H5AD")
    merge.add_argument("inputs", type=Path, nargs="+")
    args = parser.parse_args()

    ad.settings.allow_write_nullable_strings = True
    if args.command == "combine":
        combine(args)
    else:
        sample = json.loads(args.sample.read_text(encoding="utf-8"))
        try:
            aggregate(args, sample)
        except (ValueError, KeyError, OSError) as error:
            raise SystemExit(f"Sample {sample['sample_id']!r}: {error}") from error


if __name__ == "__main__":
    main()
