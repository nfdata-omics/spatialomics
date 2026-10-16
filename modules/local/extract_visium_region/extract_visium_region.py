#!/usr/bin/env python3
"""Compute Visium bounds and extract the full-resolution capture-area image."""

import argparse
import importlib
import importlib.metadata
import math
import sys
from pathlib import Path

import numpy as np
import spatialdata as sd
import tifffile
import yaml
import zarr


def validate_positive_number(value, name):
    """Validate positive numeric CLI inputs."""
    if value <= 0:
        raise ValueError(f"{name} must be > 0, got {value}")


def extract_visium_region(
    input_zarr,
    sample_name,
    microscopy_tif,
    zarr_downsample_factor,
    crop_padding,
    output_bounds,
    output_crop,
):
    """Compute Visium bounds, write their metadata, and extract the padded TIFF crop."""
    validate_positive_number(zarr_downsample_factor, "--zarr-downsample-factor")
    if crop_padding < 0:
        raise ValueError(f"--crop-padding must be >= 0, got {crop_padding}")

    sdata = sd.read_zarr(input_zarr)
    square_key = f"{sample_name}_square_016um"
    if square_key not in sdata:
        raise ValueError(f"Required SpatialData element '{square_key}' not found")

    with tifffile.TiffFile(microscopy_tif) as tif:
        series = tif.series[0]
        shape = tuple(series.shape)
        axes = series.axes

    if len(shape) != len(axes) or "Y" not in axes or "X" not in axes:
        raise ValueError(f"Unsupported microscopy TIFF layout: shape={shape}, axes='{axes}'")
    unsupported = [
        (axis, size)
        for axis, size in zip(axes, shape)
        if axis not in {"Y", "X", "S", "C"} and size != 1
    ]
    channel_axes = [axis for axis in ("S", "C") if axis in axes and shape[axes.index(axis)] > 1]
    if unsupported or len(channel_axes) > 1:
        raise ValueError(f"Unsupported microscopy TIFF layout: shape={shape}, axes='{axes}'")

    height = shape[axes.index("Y")]
    width = shape[axes.index("X")]

    minx, miny, maxx, maxy = sdata[square_key].total_bounds
    x0 = max(0, math.floor(minx * zarr_downsample_factor))
    y0 = max(0, math.floor(miny * zarr_downsample_factor))
    x1 = min(width, math.ceil(maxx * zarr_downsample_factor))
    y1 = min(height, math.ceil(maxy * zarr_downsample_factor))

    if x0 >= x1 or y0 >= y1:
        raise ValueError(
            "Visium bounds do not overlap the microscopy image after coordinate conversion: "
            f"{x0}:{y0}:{x1}:{y1}"
        )

    crop_x0 = max(0, x0 - crop_padding)
    crop_y0 = max(0, y0 - crop_padding)
    crop_x1 = min(width, x1 + crop_padding)
    crop_y1 = min(height, y1 + crop_padding)

    output_bounds.parent.mkdir(parents=True, exist_ok=True)
    output_crop.parent.mkdir(parents=True, exist_ok=True)
    with output_bounds.open("w", encoding="utf-8") as handle:
        handle.write(
            "x0\ty0\tx1\ty1\tcrop_x0\tcrop_y0\tcrop_x1\tcrop_y1\t"
            "source_width\tsource_height\tpadding\n"
        )
        handle.write(
            f"{x0}\t{y0}\t{x1}\t{y1}\t{crop_x0}\t{crop_y0}\t{crop_x1}\t{crop_y1}\t"
            f"{width}\t{height}\t{crop_padding}\n"
        )

    store = tifffile.imread(microscopy_tif, series=0, level=0, aszarr=True)
    try:
        image = zarr.open(store, mode="r")
        slices = []
        crop_axes = []
        for axis, size in zip(axes, shape):
            if axis == "Y":
                slices.append(slice(crop_y0, crop_y1))
                crop_axes.append(axis)
            elif axis == "X":
                slices.append(slice(crop_x0, crop_x1))
                crop_axes.append(axis)
            elif axis in {"S", "C"} and size > 1:
                slices.append(slice(None))
                crop_axes.append(axis)
            else:
                slices.append(0)
        crop = np.asarray(image[tuple(slices)])
    finally:
        store.close()

    if channel_axes:
        crop = np.transpose(
            crop,
            [crop_axes.index("Y"), crop_axes.index("X"), crop_axes.index(channel_axes[0])],
        )
        if crop.shape[-1] != 3:
            raise ValueError(f"Expected an RGB microscopy TIFF, got {crop.shape[-1]} channels")
        photometric = "rgb"
        output_axes = "YXS"
    else:
        crop = np.transpose(crop, [crop_axes.index("Y"), crop_axes.index("X")])
        photometric = "minisblack"
        output_axes = "YX"

    tifffile.imwrite(
        output_crop,
        crop,
        bigtiff=True,
        tile=(1024, 1024),
        photometric=photometric,
        metadata={"axes": output_axes},
    )


def versions_yaml(process_name, list_of_libs):
    """Generate YAML formatted string containing relevant library versions."""
    versions = {process_name: {}}
    versions[process_name]["python"] = (
        f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    )

    for lib in list_of_libs:
        try:
            version = importlib.metadata.version(lib)
        except importlib.metadata.PackageNotFoundError:
            try:
                module = importlib.import_module(lib)
                version = getattr(module, "__version__", "unknown")
            except (ImportError, AttributeError):
                version = None
        if version is not None:
            versions[process_name][lib] = version

    return yaml.dump(versions)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Compute Visium bounds and extract the full-resolution capture-area image"
    )
    parser.add_argument("--input-zarr", type=Path, help="Path to input SpatialData zarr folder")
    parser.add_argument("--sample-name", help="Sample name used to find SpatialData elements")
    parser.add_argument("--microscopy-tif", type=Path, help="Full-resolution microscopy TIFF")
    parser.add_argument(
        "--zarr-downsample-factor",
        type=float,
        default=1.0,
        help="Downsampling factor between TIFF image coordinates and zarr coordinates",
    )
    parser.add_argument(
        "--crop-padding",
        type=int,
        default=0,
        help="Full-resolution pixels added around each side of the capture-area crop",
    )
    parser.add_argument("--output-bounds", type=Path, help="Output TSV bounds file")
    parser.add_argument("--output-crop", type=Path, help="Output full-resolution TIFF crop")
    parser.add_argument(
        "--versions-dict",
        help="If set, print versions of relevant libraries in YAML format and exit",
    )
    return parser.parse_args()


def main():
    """Run Visium region extraction."""
    args = parse_args()
    if args.versions_dict:
        print(versions_yaml(args.versions_dict, ["spatialdata", "tifffile"]))
        return

    missing = [
        name
        for name in ("input_zarr", "sample_name", "microscopy_tif", "output_bounds", "output_crop")
        if getattr(args, name) is None
    ]
    if missing:
        raise ValueError(f"Missing required arguments: {', '.join(missing)}")

    extract_visium_region(
        input_zarr=args.input_zarr,
        sample_name=args.sample_name,
        microscopy_tif=args.microscopy_tif,
        zarr_downsample_factor=args.zarr_downsample_factor,
        crop_padding=args.crop_padding,
        output_bounds=args.output_bounds,
        output_crop=args.output_crop,
    )


if __name__ == "__main__":
    main()
