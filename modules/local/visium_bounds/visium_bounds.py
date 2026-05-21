#!/usr/bin/env python3
"""Compute Visium capture-area bounds in full-resolution microscopy pixels."""

import argparse
import importlib
import importlib.metadata
import math
import sys
from pathlib import Path

import spatialdata as sd
import tifffile
import yaml


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


def validate_positive_number(value, name):
    """Validate positive numeric CLI inputs."""
    if value <= 0:
        raise ValueError(f"{name} must be > 0, got {value}")


def compute_visium_bounds(input_zarr, sample_name, microscopy_tif, zarr_downsample_factor):
    """Return clamped full-resolution microscopy bounds for the Visium square grid."""
    validate_positive_number(zarr_downsample_factor, "--zarr-downsample-factor")

    sdata = sd.read_zarr(input_zarr)
    square_key = f"{sample_name}_square_016um"
    if square_key not in sdata:
        raise ValueError(f"Required SpatialData element '{square_key}' not found")

    microscopy = read_memmapped_2d_tiff(microscopy_tif, "microscopy image")
    height, width = microscopy.shape

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

    return x0, y0, x1, y1


def write_bounds(output_bounds, bounds):
    """Write bounds as a small TSV file."""
    output_bounds.parent.mkdir(parents=True, exist_ok=True)
    x0, y0, x1, y1 = bounds
    with output_bounds.open("w") as handle:
        handle.write("x0\ty0\tx1\ty1\n")
        handle.write(f"{x0}\t{y0}\t{x1}\t{y1}\n")


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
        description="Compute Visium capture-area bounds in full-resolution microscopy pixels"
    )
    parser.add_argument("--input-zarr", type=Path, help="Path to input SpatialData zarr folder")
    parser.add_argument("--sample-name", help="Sample name used to find SpatialData elements")
    parser.add_argument("--microscopy-tif", type=Path, help="Memmappable 2D microscopy TIFF")
    parser.add_argument(
        "--zarr-downsample-factor",
        type=float,
        default=1.0,
        help="Downsampling factor between TIFF image coordinates and zarr coordinates",
    )
    parser.add_argument("--output-bounds", type=Path, help="Output TSV bounds file")
    parser.add_argument(
        "--versions-dict",
        help="If set, print versions of relevant libraries in YAML format and exit",
    )
    return parser.parse_args()


def main():
    """Run the Visium bounds calculation."""
    args = parse_args()
    if args.versions_dict:
        print(versions_yaml(args.versions_dict, ["spatialdata", "tifffile"]))
        return

    missing = [
        name
        for name in ("input_zarr", "sample_name", "microscopy_tif", "output_bounds")
        if getattr(args, name) is None
    ]
    if missing:
        raise ValueError(f"Missing required arguments: {', '.join(missing)}")

    bounds = compute_visium_bounds(
        input_zarr=args.input_zarr,
        sample_name=args.sample_name,
        microscopy_tif=args.microscopy_tif,
        zarr_downsample_factor=args.zarr_downsample_factor,
    )
    write_bounds(args.output_bounds, bounds)


if __name__ == "__main__":
    main()
