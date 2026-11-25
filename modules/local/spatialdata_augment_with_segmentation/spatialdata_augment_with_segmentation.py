#!/usr/bin/env python3

"""Assemble the final Visium HD SpatialData Zarr with samplesheet H&E."""

import argparse
import shutil
from pathlib import Path

import anndata
import dask.array as da
import numpy as np
import spatialdata as sd
import tifffile
import xarray as xr
from spatialdata.models import Image2DModel, Labels2DModel
from spatialdata.transformations import Scale, Sequence, set_transformation


class TiffPixels:
    """Expose TIFF pixels for chunked reads without copying the memmap."""

    def __init__(self, image):
        self.image = image
        self.shape = image.shape
        self.dtype = image.dtype
        self.ndim = image.ndim

    def __getitem__(self, index):
        return self.image[index]


def read_memmapped_2d_tiff(path, label):
    """Read a TIFF as a squeezed 2D memmap."""
    try:
        image = tifffile.memmap(path, mode="r")
    except Exception as exc:  # noqa: BLE001 - preserve the memmap failure context
        raise ValueError(f"{label} TIFF is not memmappable: {path}") from exc

    if image.ndim != 2:
        image = image.squeeze()
    if image.ndim != 2:
        raise ValueError(f"Expected {label} TIFF to be 2D after squeeze, got shape {image.shape}")

    return TiffPixels(image)


def downsampled_size(size, factor):
    """Return the size produced by slicing an axis with [::factor]."""
    return 0 if size == 0 else ((size - 1) // factor) + 1


def infer_downsample_factor(full_shape, downsampled_shape):
    """Infer the integer stride that maps a full-resolution image to a mask."""
    if full_shape == downsampled_shape:
        return 1

    if len(full_shape) != 2 or len(downsampled_shape) != 2:
        raise ValueError(
            "Expected 2D shapes when inferring segmentation downsample factor; "
            f"got {full_shape} and {downsampled_shape}"
        )

    height_ratio = full_shape[0] / downsampled_shape[0]
    width_ratio = full_shape[1] / downsampled_shape[1]
    first_candidate = max(1, int(min(height_ratio, width_ratio)) - 2)
    last_candidate = int(max(height_ratio, width_ratio)) + 3

    for factor in range(first_candidate, last_candidate + 1):
        if tuple(downsampled_size(size, factor) for size in full_shape) == downsampled_shape:
            return factor

    raise ValueError(
        "Segmentation mask shape is not compatible with microscopy image shape; "
        f"got {downsampled_shape} and {full_shape}"
    )


def validate_positive_number(value, name):
    """Validate positive numeric CLI inputs."""
    if value <= 0:
        raise ValueError(f"{name} must be > 0, got {value}")


def check_layer_names(sdata, layer_names, overwrite):
    """Fail on existing output layers unless overwrite was requested."""
    existing = []
    for collection_name, key in layer_names:
        collection = getattr(sdata, collection_name)
        if key in collection:
            existing.append(f"{collection_name}['{key}']")

    if existing and not overwrite:
        joined = ", ".join(existing)
        raise ValueError(f"Output layer(s) already exist: {joined}. Use --overwrite to replace them.")


def remove_existing_layers(sdata, layer_names):
    """Remove existing layers that are about to be regenerated."""
    for collection_name, key in layer_names:
        collection = getattr(sdata, collection_name)
        if key in collection:
            del collection[key]


def add_layers(
    sdata,
    sample_name,
    segmentation_mask_tif,
    microscopy_tif,
    zarr_downsample_factor,
    microscopy_preview_max_dimension,
    chunk_size,
    overwrite,
):
    """Add segmentation and microscopy layers to a SpatialData object."""
    coordinate_system = sample_name
    if coordinate_system not in sdata.coordinate_systems:
        raise ValueError(f"Coordinate system '{coordinate_system}' not found in input zarr")

    segmentation_key = "segmentation_mask"
    microscopy_key = f"{sample_name}_microscopy"

    microscopy = read_memmapped_2d_tiff(microscopy_tif, "microscopy image")
    microscopy_downsample_factor = max(
        1,
        (max(microscopy.shape) + microscopy_preview_max_dimension - 1)
        // microscopy_preview_max_dimension,
    )
    microscopy_downsampled_key = (
        f"{sample_name}_microscopy_{microscopy_downsample_factor}x_downsampled"
    )
    print(
        "Adaptive microscopy preview: "
        f"source shape={microscopy.shape}, "
        f"maximum dimension={microscopy_preview_max_dimension}, "
        f"downsample factor={microscopy_downsample_factor}"
    )

    layer_names = [
        ("images", microscopy_key),
        ("images", microscopy_downsampled_key),
    ]
    if segmentation_mask_tif is not None:
        layer_names.append(("labels", segmentation_key))
    check_layer_names(sdata, layer_names, overwrite)
    if overwrite:
        remove_existing_layers(sdata, layer_names)

    chunks_2d = (chunk_size, chunk_size)
    chunks_3d = (1, chunk_size, chunk_size)

    image_to_sample = Scale(
        [1.0 / zarr_downsample_factor, 1.0 / zarr_downsample_factor],
        axes=("x", "y"),
    )

    if segmentation_mask_tif is not None:
        mask = read_memmapped_2d_tiff(segmentation_mask_tif, "segmentation mask")
        segmentation_downsample_factor = infer_downsample_factor(microscopy.shape, mask.shape)
        mask_dask = da.from_array(mask, chunks=chunks_2d, name=False).astype("uint32")
        labels = Labels2DModel.parse(mask_dask, dims=("y", "x"), chunks=chunks_2d)
        if segmentation_downsample_factor == 1:
            mask_to_sample = image_to_sample
        else:
            mask_to_sample = Sequence(
                [
                    Scale(
                        [segmentation_downsample_factor, segmentation_downsample_factor],
                        axes=("y", "x"),
                    ),
                    image_to_sample,
                ]
            )
        set_transformation(labels, mask_to_sample, to_coordinate_system=coordinate_system)
        sdata.labels[segmentation_key] = labels

    microscopy_dask = da.from_array(microscopy, chunks=chunks_2d, name=False)
    microscopy_xr = xr.DataArray(
        microscopy_dask[None, :, :],
        dims=("c", "y", "x"),
        coords={"c": ["microscopy"]},
    )
    microscopy_image = Image2DModel.parse(
        microscopy_xr,
        dims=("c", "y", "x"),
        chunks=chunks_3d,
    )
    set_transformation(microscopy_image, image_to_sample, to_coordinate_system=coordinate_system)
    sdata.images[microscopy_key] = microscopy_image

    # Sample the centers of a reduced grid spanning the native image field.
    # Dimension-specific scales account for rounding of the reduced dimensions.
    height, width = microscopy.shape
    preview_height = downsampled_size(height, microscopy_downsample_factor)
    preview_width = downsampled_size(width, microscopy_downsample_factor)
    y_scale = height / preview_height
    x_scale = width / preview_width
    downsampled = microscopy_dask
    if microscopy_downsample_factor > 1:
        y_indices = ((np.arange(preview_height) + 0.5) * y_scale).astype(int)
        x_indices = ((np.arange(preview_width) + 0.5) * x_scale).astype(int)
        downsampled = da.take(da.take(microscopy_dask, y_indices, axis=0), x_indices, axis=1)
    downsampled = downsampled.rechunk(chunks_2d)
    downsampled_xr = xr.DataArray(
        downsampled[None, :, :],
        dims=("c", "y", "x"),
        coords={"c": [f"microscopy_{microscopy_downsample_factor}x_downsampled"]},
    )
    downsampled_image = Image2DModel.parse(
        downsampled_xr,
        dims=("c", "y", "x"),
        chunks=chunks_3d,
    )
    downsampled_to_sample = Sequence(
        [
            Scale(
                [y_scale, x_scale],
                axes=("y", "x"),
            ),
            image_to_sample,
        ]
    )
    set_transformation(
        downsampled_image,
        downsampled_to_sample,
        to_coordinate_system=coordinate_system,
    )
    sdata.images[microscopy_downsampled_key] = downsampled_image


def save_zarr(sdata, output_zarr, overwrite):
    """Write the modified SpatialData object to a zarr folder."""
    if output_zarr.exists():
        if not overwrite:
            raise FileExistsError(f"Output zarr already exists: {output_zarr}")
        if output_zarr.is_dir():
            shutil.rmtree(output_zarr)
        else:
            output_zarr.unlink()

    output_zarr.parent.mkdir(parents=True, exist_ok=True)
    anndata.settings.allow_write_nullable_strings = True
    sdata.write(output_zarr)


def main():
    parser = argparse.ArgumentParser(description="Add samplesheet H&E and optional segmentation to SpatialData")
    parser.add_argument("--input-zarr", type=Path, required=True)
    parser.add_argument("--output-zarr", type=Path, required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--microscopy-tif", type=Path, required=True)
    parser.add_argument("--segmentation-mask-tif", type=Path)
    parser.add_argument("--zarr-downsample-factor", type=float, required=True)
    parser.add_argument("--microscopy-preview-max-dimension", type=int, required=True)
    parser.add_argument("--chunk-size", type=int, default=2048)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    validate_positive_number(args.zarr_downsample_factor, "--zarr-downsample-factor")
    validate_positive_number(args.microscopy_preview_max_dimension, "--microscopy-preview-max-dimension")
    validate_positive_number(args.chunk_size, "--chunk-size")
    if args.input_zarr.resolve() == args.output_zarr.resolve():
        parser.error("--output-zarr must be different from --input-zarr")

    sdata = sd.read_zarr(args.input_zarr)
    add_layers(
        sdata,
        args.sample_name,
        args.segmentation_mask_tif,
        args.microscopy_tif,
        args.zarr_downsample_factor,
        args.microscopy_preview_max_dimension,
        args.chunk_size,
        args.overwrite,
    )
    save_zarr(sdata, args.output_zarr, args.overwrite)


if __name__ == "__main__":
    main()
