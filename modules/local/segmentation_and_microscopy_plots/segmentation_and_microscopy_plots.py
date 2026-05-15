#!/usr/bin/env python3
"""Add segmentation and microscopy layers to a Visium HD SpatialData object."""

import argparse
import csv
import importlib
import importlib.metadata
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path

import dask.array as da
import geopandas as gpd
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import spatialdata as sd
import spatialdata_plot  # noqa: F401  # pylint: disable=unused-import
import tifffile
import xarray as xr
from shapely.geometry import box
from spatialdata.models import Image2DModel, Labels2DModel, ShapesModel
from spatialdata.transformations import Identity, Scale, Sequence, set_transformation


@dataclass(frozen=True)
class CropArea:
    """Cropping area in target coordinate-system coordinates."""

    x0: float
    y0: float
    x1: float
    y1: float

    @property
    def bounds_for_query(self):
        """Return min/max coordinates in the axis order expected by SpatialData."""
        return [self.y0, self.x0], [self.y1, self.x1]


def parse_crop_areas(crop_areas):
    """Parse x0:y0:x1:y1 crop areas separated by semicolons."""
    if not crop_areas:
        return None

    parsed = []
    for raw_area in crop_areas.split(";"):
        raw_area = raw_area.strip()
        if not raw_area:
            continue

        parts = raw_area.split(":")
        if len(parts) != 4:
            raise ValueError(
                "Invalid crop area "
                f"'{raw_area}'. Expected syntax: x0:y0:x1:y1[;x0:y0:x1:y1]"
            )

        try:
            x0, y0, x1, y1 = [float(part) for part in parts]
        except ValueError as exc:
            raise ValueError(f"Crop area '{raw_area}' contains non-numeric coordinates") from exc

        if x0 >= x1 or y0 >= y1:
            raise ValueError(f"Crop area '{raw_area}' must satisfy x0 < x1 and y0 < y1")

        parsed.append(CropArea(x0=x0, y0=y0, x1=x1, y1=y1))

    if not parsed:
        raise ValueError("--crop-areas was provided but no valid crop areas were found")

    return parsed


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


def validate_paths(input_zarr, output_zarr):
    """Validate input and optional output zarr paths."""
    if not input_zarr.exists():
        raise FileNotFoundError(f"Input zarr does not exist: {input_zarr}")

    if output_zarr and input_zarr.resolve() == output_zarr.resolve():
        raise ValueError("--output-zarr must be different from --input-zarr")


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


def visium_bounds(sdata, square_key):
    """Return the full Visium square-grid bounds."""
    minx, miny, maxx, maxy = sdata[square_key].total_bounds
    return CropArea(x0=minx, y0=miny, x1=maxx, y1=maxy)


def fullres_bounds_to_mask_slices(area, microscopy_shape, mask_shape, zarr_downsample_factor, mask_downsample_factor):
    """Convert sample-coordinate bounds to full-resolution pixels and mask slices."""
    x0 = max(0.0, area.x0 * zarr_downsample_factor)
    y0 = max(0.0, area.y0 * zarr_downsample_factor)
    x1 = min(float(microscopy_shape[1]), area.x1 * zarr_downsample_factor)
    y1 = min(float(microscopy_shape[0]), area.y1 * zarr_downsample_factor)

    if x0 >= x1 or y0 >= y1:
        raise ValueError(
            "Visium bounds do not overlap the microscopy image after coordinate conversion: "
            f"{x0:g}:{y0:g}:{x1:g}:{y1:g}"
        )

    mask_x0 = max(0, int(np.floor(x0 / mask_downsample_factor)))
    mask_y0 = max(0, int(np.floor(y0 / mask_downsample_factor)))
    mask_x1 = min(mask_shape[1], int(np.ceil(x1 / mask_downsample_factor)))
    mask_y1 = min(mask_shape[0], int(np.ceil(y1 / mask_downsample_factor)))

    if mask_x0 >= mask_x1 or mask_y0 >= mask_y1:
        raise ValueError(
            "Visium bounds do not overlap the segmentation mask after coordinate conversion: "
            f"{mask_x0}:{mask_y0}:{mask_x1}:{mask_y1}"
        )

    return (
        CropArea(x0=x0, y0=y0, x1=x1, y1=y1),
        slice(mask_y0, mask_y1),
        slice(mask_x0, mask_x1),
    )


def label_counts(mask, chunk_size, y_slice=None, x_slice=None):
    """Count non-zero label pixels in chunks."""
    if y_slice is None:
        y_slice = slice(0, mask.shape[0])
    if x_slice is None:
        x_slice = slice(0, mask.shape[1])

    counts = {}
    y_start = 0 if y_slice.start is None else y_slice.start
    y_stop = mask.shape[0] if y_slice.stop is None else y_slice.stop
    x_start = 0 if x_slice.start is None else x_slice.start
    x_stop = mask.shape[1] if x_slice.stop is None else x_slice.stop

    for y0 in range(y_start, y_stop, chunk_size):
        y1 = min(y0 + chunk_size, y_stop)
        for x0 in range(x_start, x_stop, chunk_size):
            x1 = min(x0 + chunk_size, x_stop)
            chunk = np.asarray(mask[y0:y1, x0:x1])
            labels, label_counts_array = np.unique(chunk, return_counts=True)
            for label, count in zip(labels, label_counts_array):
                label = int(label)
                if label == 0:
                    continue
                counts[label] = counts.get(label, 0) + int(count)

    return counts


def summarize_segmentation_region(counts, total_mask_pixels, fullres_area_px, area_scale, prefix):
    """Summarize label counts for one region."""
    segmented_mask_pixels = int(sum(counts.values()))
    n_segments = len(counts)
    row = {
        f"{prefix}_n_segments": n_segments,
        f"{prefix}_segmented_fraction": (
            segmented_mask_pixels / total_mask_pixels if total_mask_pixels else 0.0
        ),
        f"{prefix}_segment_density_per_megapixel": (
            n_segments / (fullres_area_px / 1_000_000) if fullres_area_px else 0.0
        ),
    }

    area_columns = {
        "mean": f"{prefix}_mean_segment_area_fullres_px",
        "median": f"{prefix}_median_segment_area_fullres_px",
        "min": f"{prefix}_min_segment_area_fullres_px",
        "max": f"{prefix}_max_segment_area_fullres_px",
        "p05": f"{prefix}_p05_segment_area_fullres_px",
        "p95": f"{prefix}_p95_segment_area_fullres_px",
    }

    if counts:
        areas = np.asarray(list(counts.values()), dtype=np.float64) * area_scale
        row.update(
            {
                area_columns["mean"]: float(np.mean(areas)),
                area_columns["median"]: float(np.median(areas)),
                area_columns["min"]: float(np.min(areas)),
                area_columns["max"]: float(np.max(areas)),
                area_columns["p05"]: float(np.percentile(areas, 5)),
                area_columns["p95"]: float(np.percentile(areas, 95)),
            }
        )
    else:
        row.update({column: 0.0 for column in area_columns.values()})

    return row


def segmentation_statistics_row(
    sample_name,
    sdata,
    square_key,
    segmentation_mask_tif,
    microscopy_tif,
    zarr_downsample_factor,
    chunk_size,
):
    """Compute whole-image and Visium-area segmentation statistics."""
    mask = read_memmapped_2d_tiff(segmentation_mask_tif, "segmentation mask")
    microscopy = read_memmapped_2d_tiff(microscopy_tif, "microscopy image")
    mask_downsample_factor = infer_downsample_factor(microscopy.shape, mask.shape)
    visium_area = visium_bounds(sdata, square_key)
    visium_fullres_area, visium_y_slice, visium_x_slice = fullres_bounds_to_mask_slices(
        visium_area,
        microscopy.shape,
        mask.shape,
        zarr_downsample_factor,
        mask_downsample_factor,
    )

    area_scale = mask_downsample_factor * mask_downsample_factor
    whole_fullres_area_px = microscopy.shape[0] * microscopy.shape[1]
    visium_fullres_area_px = (
        (visium_fullres_area.x1 - visium_fullres_area.x0)
        * (visium_fullres_area.y1 - visium_fullres_area.y0)
    )

    whole_counts = label_counts(mask, chunk_size)
    visium_counts = label_counts(mask, chunk_size, visium_y_slice, visium_x_slice)

    row = {
        "sample": sample_name,
        "segmentation_downsample_factor": mask_downsample_factor,
        "microscopy_height_px": microscopy.shape[0],
        "microscopy_width_px": microscopy.shape[1],
        "segmentation_height_px": mask.shape[0],
        "segmentation_width_px": mask.shape[1],
        "visium_x0_px": visium_fullres_area.x0,
        "visium_y0_px": visium_fullres_area.y0,
        "visium_x1_px": visium_fullres_area.x1,
        "visium_y1_px": visium_fullres_area.y1,
        "visium_width_px": visium_fullres_area.x1 - visium_fullres_area.x0,
        "visium_height_px": visium_fullres_area.y1 - visium_fullres_area.y0,
    }
    row.update(
        summarize_segmentation_region(
            whole_counts,
            mask.shape[0] * mask.shape[1],
            whole_fullres_area_px,
            area_scale,
            "whole",
        )
    )
    row.update(
        summarize_segmentation_region(
            visium_counts,
            (visium_y_slice.stop - visium_y_slice.start)
            * (visium_x_slice.stop - visium_x_slice.start),
            visium_fullres_area_px,
            area_scale,
            "visium",
        )
    )
    return row


def save_segmentation_statistics(row, output_path):
    """Write one-row segmentation statistics CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def crop_areas_to_shapes(crop_areas, coordinate_system):
    """Create a SpatialData shapes model containing crop boxes."""
    crop_gdf = gpd.GeoDataFrame(
        {"crop": [f"crop_{i + 1}" for i in range(len(crop_areas))]},
        geometry=[box(area.x0, area.y0, area.x1, area.y1) for area in crop_areas],
    )
    return ShapesModel.parse(crop_gdf, transformations={coordinate_system: Identity()})


def add_layers(
    sdata,
    sample_name,
    segmentation_mask_tif,
    microscopy_tif,
    zarr_downsample_factor,
    microscopy_downsample_factor,
    chunk_size,
    overwrite,
):
    """Add segmentation and microscopy layers to a SpatialData object."""
    coordinate_system = sample_name
    if coordinate_system not in sdata.coordinate_systems:
        raise ValueError(f"Coordinate system '{coordinate_system}' not found in input zarr")

    segmentation_key = "segmentation_mask"
    microscopy_key = f"{sample_name}_microscopy"
    microscopy_downsampled_key = (
        f"{sample_name}_microscopy_{microscopy_downsample_factor}x_downsampled"
    )

    layer_names = [
        ("labels", segmentation_key),
        ("images", microscopy_key),
        ("images", microscopy_downsampled_key),
    ]
    check_layer_names(sdata, layer_names, overwrite)
    if overwrite:
        remove_existing_layers(sdata, layer_names)

    mask = read_memmapped_2d_tiff(segmentation_mask_tif, "segmentation mask")
    microscopy = read_memmapped_2d_tiff(microscopy_tif, "microscopy image")
    segmentation_downsample_factor = infer_downsample_factor(microscopy.shape, mask.shape)

    chunks_2d = (chunk_size, chunk_size)
    chunks_3d = (1, chunk_size, chunk_size)

    image_to_sample = Scale(
        [1.0 / zarr_downsample_factor, 1.0 / zarr_downsample_factor],
        axes=("x", "y"),
    )

    mask_dask = da.from_array(mask, chunks=chunks_2d).astype("uint32")
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

    microscopy_dask = da.from_array(microscopy, chunks=chunks_2d)
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

    downsampled = microscopy_dask[::microscopy_downsample_factor, ::microscopy_downsample_factor]
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
                [microscopy_downsample_factor, microscopy_downsample_factor],
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

    return segmentation_key, microscopy_key, microscopy_downsampled_key


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
    sdata.write(output_zarr)


def save_registration_plot(
    sdata,
    sample_name,
    square_key,
    hires_key,
    microscopy_downsampled_key,
    output_path,
):
    """Save side-by-side full-slide registration plot."""
    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    sdata.pl.render_images(
        hires_key,
        cmap="gray",
    ).pl.render_shapes(
        square_key,
        color="#00000000",
        fill_alpha=0,
        outline_color="yellow",
        outline_width=0.01,
        outline_alpha=0.4,
    ).pl.show(
        coordinate_systems=sample_name,
        ax=axes[0],
        title=f"{sample_name} - CytAssist image",
    )

    sdata.pl.render_images(
        microscopy_downsampled_key,
        cmap="gray",
    ).pl.render_shapes(
        square_key,
        color="#00000000",
        fill_alpha=0,
        outline_color="yellow",
        outline_width=0.01,
        outline_alpha=0.4,
    ).pl.show(
        coordinate_systems=sample_name,
        ax=axes[1],
        title=f"{sample_name} - microscopy image",
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_crop_area_plot(
    sdata,
    sample_name,
    crop_shape_key,
    microscopy_downsampled_key,
    visium_area,
    output_path,
):
    """Save downsampled microscopy plot with crop-area overlays."""
    min_coordinate, max_coordinate = visium_area.bounds_for_query
    visium_crop = sd.bounding_box_query(
        sdata,
        axes=("y", "x"),
        min_coordinate=min_coordinate,
        max_coordinate=max_coordinate,
        target_coordinate_system=sample_name,
    )

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))

    sdata.pl.render_images(
        microscopy_downsampled_key,
        cmap="gray",
    ).pl.render_shapes(
        crop_shape_key,
        fill_alpha=0,
        outline_color="red",
        outline_width=2,
        outline_alpha=1,
    ).pl.show(
        coordinate_systems=sample_name,
        ax=axes[0],
        title="Full downsampled microscopy",
    )

    visium_crop.pl.render_images(
        microscopy_downsampled_key,
        cmap="gray",
    ).pl.render_shapes(
        crop_shape_key,
        fill_alpha=0,
        outline_color="red",
        outline_width=2,
        outline_alpha=1,
    ).pl.show(
        coordinate_systems=sample_name,
        ax=axes[1],
        title="Downsampled microscopy, Visium area",
    )

    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def save_segmentation_crop_panels(
    sdata,
    sample_name,
    crop_areas,
    microscopy_key,
    segmentation_key,
    output_path,
):
    """Save high-resolution crop panels with and without segmentation overlay."""
    n_crops = len(crop_areas)
    fig_height = max(5, 4 * n_crops)
    fig, axes = plt.subplots(n_crops, 2, figsize=(12, fig_height), squeeze=False)

    for index, crop_area in enumerate(crop_areas):
        min_coordinate, max_coordinate = crop_area.bounds_for_query
        crop_sdata = sd.bounding_box_query(
            sdata,
            axes=("y", "x"),
            min_coordinate=min_coordinate,
            max_coordinate=max_coordinate,
            target_coordinate_system=sample_name,
        )

        crop_label = (
            f"Crop {index + 1}: "
            f"{crop_area.x0:g}:{crop_area.y0:g}:{crop_area.x1:g}:{crop_area.y1:g}"
        )

        crop_sdata.pl.render_images(
            microscopy_key,
            cmap="gray",
        ).pl.show(
            coordinate_systems=sample_name,
            ax=axes[index, 0],
            title=f"{crop_label} - microscopy",
        )

        crop_sdata.pl.render_images(
            microscopy_key,
            cmap="gray",
        ).pl.render_labels(
            segmentation_key,
            fill_alpha=0,
            outline_alpha=1,
            outline_color="red",
            na_color="red",
        ).pl.show(
            coordinate_systems=sample_name,
            ax=axes[index, 1],
            title=f"{crop_label} - segmentation",
        )

    fig.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def save_result_plots(
    sdata,
    sample_name,
    crop_areas,
    segmentation_key,
    microscopy_key,
    microscopy_downsampled_key,
    results_dir,
):
    """Save all result PNG plots."""
    square_key = f"{sample_name}_square_016um"
    hires_key = f"{sample_name}_hires_image"

    for key in [square_key, hires_key]:
        if key not in sdata:
            raise ValueError(f"Required SpatialData element '{key}' not found")

    visium_area = visium_bounds(sdata, square_key)
    if crop_areas is None:
        crop_areas = [visium_area]

    results_dir.mkdir(parents=True, exist_ok=True)

    crop_shape_key = "crop_areas_for_plotting"
    if crop_shape_key in sdata.shapes:
        del sdata.shapes[crop_shape_key]
    sdata.shapes[crop_shape_key] = crop_areas_to_shapes(crop_areas, sample_name)

    try:
        save_registration_plot(
            sdata,
            sample_name,
            square_key,
            hires_key,
            microscopy_downsampled_key,
            results_dir / f"{sample_name}_registration_full_slide_mqc.png",
        )
        save_crop_area_plot(
            sdata,
            sample_name,
            crop_shape_key,
            microscopy_downsampled_key,
            visium_area,
            results_dir / f"{sample_name}_crop_areas_downsampled_microscopy_mqc.png",
        )
        save_segmentation_crop_panels(
            sdata,
            sample_name,
            crop_areas,
            microscopy_key,
            segmentation_key,
            results_dir / f"{sample_name}_segmentation_crop_panels_mqc.png",
        )
    finally:
        if crop_shape_key in sdata.shapes:
            del sdata.shapes[crop_shape_key]


def run(args):
    """Run the conversion and plotting workflow."""
    input_zarr = Path(args.input_zarr)
    output_zarr = Path(args.output_zarr) if args.output_zarr else None
    segmentation_mask_tif = Path(args.segmentation_mask_tif)
    microscopy_tif = Path(args.microscopy_tif)
    results_dir = Path(args.results_dir)

    validate_positive_number(args.zarr_downsample_factor, "--zarr-downsample-factor")
    validate_positive_number(args.microscopy_downsample_factor, "--microscopy-downsample-factor")
    validate_positive_number(args.chunk_size, "--chunk-size")
    validate_paths(input_zarr, output_zarr)

    if not segmentation_mask_tif.exists():
        raise FileNotFoundError(f"Segmentation mask TIFF does not exist: {segmentation_mask_tif}")
    if not microscopy_tif.exists():
        raise FileNotFoundError(f"Microscopy TIFF does not exist: {microscopy_tif}")

    crop_areas = parse_crop_areas(args.crop_areas)
    sdata = sd.read_zarr(input_zarr)
    square_key = f"{args.sample_name}_square_016um"

    if square_key not in sdata:
        raise ValueError(f"Required SpatialData element '{square_key}' not found")

    segmentation_key, microscopy_key, microscopy_downsampled_key = add_layers(
        sdata=sdata,
        sample_name=args.sample_name,
        segmentation_mask_tif=segmentation_mask_tif,
        microscopy_tif=microscopy_tif,
        zarr_downsample_factor=args.zarr_downsample_factor,
        microscopy_downsample_factor=args.microscopy_downsample_factor,
        chunk_size=args.chunk_size,
        overwrite=args.overwrite,
    )

    stats_row = segmentation_statistics_row(
        sample_name=args.sample_name,
        sdata=sdata,
        square_key=square_key,
        segmentation_mask_tif=segmentation_mask_tif,
        microscopy_tif=microscopy_tif,
        zarr_downsample_factor=args.zarr_downsample_factor,
        chunk_size=args.chunk_size,
    )
    save_segmentation_statistics(
        stats_row,
        results_dir / f"{args.sample_name}_segmentation_stats.csv",
    )

    if output_zarr is not None:
        save_zarr(sdata, output_zarr, args.overwrite)

    save_result_plots(
        sdata=sdata,
        sample_name=args.sample_name,
        crop_areas=crop_areas,
        segmentation_key=segmentation_key,
        microscopy_key=microscopy_key,
        microscopy_downsampled_key=microscopy_downsampled_key,
        results_dir=results_dir,
    )


def versions_yaml(process_name, list_of_libs):
    """Generate YAML formatted string with versions of relevant libraries."""
    import yaml

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


def build_parser():
    """Build command-line parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Add segmentation and microscopy TIFF layers to a Visium HD SpatialData zarr "
            "and generate registration/segmentation result plots."
        )
    )
    parser.add_argument("--input-zarr", help="Path to input SpatialData zarr folder")
    parser.add_argument("--output-zarr", help="Optional output SpatialData zarr folder")
    parser.add_argument("--sample-name", help="Sample name and target coordinate system")
    parser.add_argument("--segmentation-mask-tif", help="Memmappable 2D segmentation mask TIFF")
    parser.add_argument("--microscopy-tif", help="Memmappable 2D microscopy image TIFF")
    parser.add_argument(
        "--zarr-downsample-factor",
        type=float,
        help="Factor by which the zarr coordinate space is downsampled relative to TIFFs",
    )
    parser.add_argument(
        "--microscopy-downsample-factor",
        type=int,
        help="Factor used to create the downsampled microscopy image layer",
    )
    parser.add_argument(
        "--crop-areas",
        help="Optional crop areas as x0:y0:x1:y1[;x0:y0:x1:y1] in sample coordinates",
    )
    parser.add_argument(
        "--results-dir",
        default=".",
        help="Directory for result PNG files. Default: current directory",
    )
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=2048,
        help="Chunk size for TIFF-backed dask arrays. Default: 2048",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output layers and output zarr path if present",
    )
    parser.add_argument(
        "--versions-dict",
        help="Return dictionary of versions used by the script and exit",
    )
    return parser


def main():
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()

    if args.versions_dict:
        libs = [
            "dask",
            "geopandas",
            "matplotlib",
            "shapely",
            "spatialdata",
            "spatialdata-plot",
            "tifffile",
            "xarray",
            "yaml",
        ]
        print(versions_yaml(args.versions_dict, libs))
        return

    required = [
        "input_zarr",
        "sample_name",
        "segmentation_mask_tif",
        "microscopy_tif",
        "zarr_downsample_factor",
        "microscopy_downsample_factor",
    ]
    missing = [f"--{name.replace('_', '-')}" for name in required if getattr(args, name) is None]
    if missing:
        parser.error(f"Missing required arguments: {', '.join(missing)}")

    run(args)


if __name__ == "__main__":
    main()
