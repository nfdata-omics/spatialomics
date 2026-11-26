#!/usr/bin/env python3

"""Render reports from final SpatialData and the matching Space Ranger output."""

import argparse
import json
from pathlib import Path

import geopandas as gpd
import h5py
import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import spatialdata as sd
import spatialdata_plot  # noqa: F401
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon
from PIL import Image
from shapely.geometry import Polygon as ShapelyPolygon
from skimage.transform import ProjectiveTransform
from spatialdata.models import Image2DModel, ShapesModel, get_table_keys
from spatialdata.transformations import Identity, get_transformation

from segmentation_report import (
    parse_crop_areas,
    save_result_plots,
    save_segmentation_statistics,
    segmentation_statistics_row,
)
from spatial_preview_layout import (
    PREVIEW_MARGIN_PX,
    PREVIEW_PANEL_LEFTS_PX,
    PREVIEW_PANEL_WIDTH_PX,
    PREVIEW_WIDTH_PX,
    new_preview_figure,
    preview_add_axes,
    preview_add_title,
    preview_finish_map_axis,
    preview_remove_extra_axes,
    save_preview_webp,
)


def grayscale_cytassist(image):
    """Convert CytAssist to luminance while retaining its stored transforms."""
    channel_names = set(image.coords["c"].values.tolist())
    if {"r", "g", "b"}.issubset(channel_names):
        gray = 0.2126 * image.sel(c="r") + 0.7152 * image.sel(c="g") + 0.0722 * image.sel(c="b")
    else:
        gray = image.mean(dim="c")
    gray = gray.expand_dims(c=["gray"]).transpose("c", "y", "x")
    gray.attrs = {}
    return Image2DModel.parse(gray, transformations=get_transformation(image, get_all=True))


def validated_polygon(sample_id, name, vertices):
    polygon = ShapelyPolygon(vertices)
    if not polygon.is_valid or polygon.area <= 0:
        raise ValueError(f"Invalid {name} polygon for {sample_id}")
    return polygon


def load_capture_area(spaceranger_outs, sample_id):
    positions_path = (
        spaceranger_outs
        / "binned_outputs"
        / "square_016um"
        / "spatial"
        / "tissue_positions.parquet"
    )
    if not positions_path.exists():
        raise FileNotFoundError(f"Cannot construct the Visium HD capture area without {positions_path}")

    positions = pd.read_parquet(
        positions_path,
        columns=["array_col", "array_row", "pxl_col_in_fullres", "pxl_row_in_fullres"],
    ).set_index(["array_col", "array_row"])
    if positions.empty:
        raise ValueError(f"Space Ranger tissue positions are empty: {positions_path}")

    cols = positions.index.get_level_values("array_col")
    rows = positions.index.get_level_values("array_row")
    left, right, top, bottom = cols.min(), cols.max(), rows.min(), rows.max()
    # Project the outer corners at half a bin beyond the corner bin centers.
    center_indices = [(left, top), (right, top), (right, bottom), (left, bottom)]
    corner_positions = positions.reindex(center_indices)
    if corner_positions.isna().any().any() or left == right or top == bottom:
        raise ValueError(f"Expected four distinct corner bin centers in {positions_path}")
    bins_to_image = ProjectiveTransform.from_estimate(
        center_indices, corner_positions.to_numpy()
    )
    if not bins_to_image:
        raise ValueError(f"Cannot project the capture-area corners from {positions_path}")
    vertices = bins_to_image(
        [
            (left - 0.5, top - 0.5), (right + 0.5, top - 0.5),
            (right + 0.5, bottom + 0.5), (left - 0.5, bottom + 0.5),
        ]
    )
    polygon = validated_polygon(sample_id, "capture_area", vertices)
    return ShapesModel.parse(
        gpd.GeoDataFrame(geometry=[polygon]),
        transformations={sample_id: Identity()},
    )


def load_cytassist_footprint(spaceranger_outs, sample_id):
    cytassist_path = spaceranger_outs / "spatial" / "cytassist_image.tiff"
    has_cytassist_image = cytassist_path.is_file()
    feature_slice_paths = list(spaceranger_outs.glob("*feature_slice.h5"))
    if not has_cytassist_image or not feature_slice_paths:
        print(
            f"Rendering the CytAssist overlay without clipping for {sample_id}: "
            f"image present={has_cytassist_image}, "
            f"feature_slice.h5 present={bool(feature_slice_paths)}"
        )
        return None
    if len(feature_slice_paths) != 1:
        raise ValueError(
            "Expected exactly one feature_slice.h5 file in the Space Ranger output, "
            f"found {len(feature_slice_paths)}"
        )

    with h5py.File(feature_slice_paths[0], "r") as feature_slice:
        metadata_json = feature_slice.attrs.get("metadata_json")
        metadata = json.loads(metadata_json) if metadata_json is not None else {}

    transform_matrices = metadata.get("transform_matrices", {})
    cytassist_to_array = transform_matrices.get("cytassist_colrow_to_spot_colrow")
    array_to_microscope = transform_matrices.get("spot_colrow_to_microscope_colrow")
    if cytassist_to_array is None or array_to_microscope is None:
        print(
            f"Rendering the CytAssist overlay without clipping for {sample_id}: "
            "registration matrices are incomplete"
        )
        return None

    with Image.open(cytassist_path) as image:
        width, height = image.size
    cytassist_to_array = np.asarray(cytassist_to_array).reshape(3, 3)
    array_to_microscope = np.asarray(array_to_microscope).reshape(3, 3)
    cytassist_to_microscope = ProjectiveTransform(
        array_to_microscope @ cytassist_to_array
    )
    corners = [
        (-0.5, -0.5), (width - 0.5, -0.5),
        (width - 0.5, height - 0.5), (-0.5, height - 0.5),
    ]
    return validated_polygon(
        sample_id, "cytassist_footprint", cytassist_to_microscope(corners)
    )


def save_registration_tissue_preview(
    sdata, sample_id, resolution, preview_key, capture, cytassist_footprint, output_path
):
    """Show the original Space Ranger selection and, where available, registration."""
    capture_key = f"{sample_id}_capture_area"
    cytassist_key = f"{sample_id}_cytassist_image"
    square_key = f"{sample_id}_{resolution}"
    preview = sdata.images[preview_key]
    cytassist = sdata.images.get(cytassist_key)

    capture_extent = sd.get_extent(capture, coordinate_system=sample_id)
    minx, maxx = capture_extent["x"]
    miny, maxy = capture_extent["y"]
    padding = max(maxx - minx, maxy - miny) * 0.02
    tissue_viewport = (minx - padding, miny - padding, maxx + padding, maxy + padding)

    table = sdata.tables[resolution]
    _, region_key, instance_key = get_table_keys(table)
    tissue_ids = table.obs.loc[
        (table.obs[region_key] == square_key) & (table.obs["in_tissue"] == 1),
        instance_key,
    ]
    squares = sdata.shapes[square_key]
    tissue_shapes = squares.loc[squares.index.isin(tissue_ids)]
    tissue = sd.SpatialData(images={preview_key: preview}, shapes={capture_key: capture})

    panel_boxes = [
        (left, 62, PREVIEW_PANEL_WIDTH_PX, 455) for left in PREVIEW_PANEL_LEFTS_PX
    ] if cytassist is not None else [
        (PREVIEW_MARGIN_PX, 62, PREVIEW_WIDTH_PX - 2 * PREVIEW_MARGIN_PX, 455)
    ]
    fig = new_preview_figure(570)
    axes = [preview_add_axes(fig, *box) for box in panel_boxes]
    try:
        if cytassist is not None:
            overview = sd.SpatialData(
                images={preview_key: preview, cytassist_key: grayscale_cytassist(cytassist)},
                shapes={capture_key: capture},
            )
            cytassist_cmap = LinearSegmentedColormap.from_list(
                "cytassist_magenta", ["#7600c8", "#d000ff", "#ff5ad9", "#fff0fa"]
            )
            overview.pl.render_images(
                preview_key, cmap="gray", colorbar=False,
            ).pl.render_images(
                cytassist_key, cmap=cytassist_cmap, alpha=0.55, colorbar=False,
            ).pl.render_shapes(
                capture_key, fill_alpha=0, outline_color="#39ff14", outline_width=2,
                outline_alpha=1,
            ).pl.show(coordinate_systems=sample_id, ax=axes[0], colorbar=False)
            if cytassist_footprint is not None:
                # SpatialData Plot draws the H&E and CytAssist image artists in
                # layer order. The footprint is in sample coordinates.
                axes[0].images[-1].set_clip_path(
                    Polygon(cytassist_footprint.exterior.coords, transform=axes[0].transData)
                )
            extent = sd.get_extent(preview, coordinate_system=sample_id)
            viewport = (extent["x"][0], extent["y"][0], extent["x"][1], extent["y"][1])
            preview_finish_map_axis(axes[0], panel_boxes[0], viewport)
            preview_add_title(fig, panel_boxes[0][0], "H&E / CytAssist registration")
        else:
            print(f"Registered CytAssist layer '{cytassist_key}' is unavailable for {sample_id}; rendering tissue selection only")

        if not tissue_shapes.empty:
            tissue.shapes[square_key] = tissue_shapes
        layers = tissue.pl.render_images(preview_key, cmap="gray", colorbar=False)
        if not tissue_shapes.empty:
            layers = layers.pl.render_shapes(
                square_key, color="#00d5ff", fill_alpha=0.45, outline_width=0,
            )
        layers.pl.render_shapes(
            capture_key, fill_alpha=0, outline_color="#ffffff", outline_width=1.5,
            outline_alpha=0.9,
        ).pl.show(coordinate_systems=sample_id, ax=axes[-1], colorbar=False)
        preview_finish_map_axis(axes[-1], panel_boxes[-1], tissue_viewport)
        preview_add_title(
            fig, panel_boxes[-1][0], "Space Ranger in-tissue bins", width_px=panel_boxes[-1][2]
        )
        preview_remove_extra_axes(fig, axes)
        save_preview_webp(fig, output_path)
    finally:
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Generate reports from final SpatialData")
    parser.add_argument("--zarr", required=True)
    parser.add_argument("--spaceranger-outs", type=Path, required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--resolution", required=True)
    parser.add_argument("--crop-areas")
    parser.add_argument("--segmentation", action="store_true")
    parser.add_argument("--zarr-downsample-factor", type=float, required=True)
    parser.add_argument("--chunk-size", type=int, required=True)
    args = parser.parse_args()

    sdata = sd.read_zarr(args.zarr)
    sample = args.sample
    capture = load_capture_area(args.spaceranger_outs, sample)
    cytassist_key = f"{sample}_cytassist_image"
    cytassist_footprint = None
    if cytassist_key in sdata.images:
        cytassist_footprint = load_cytassist_footprint(args.spaceranger_outs, sample)
    microscopy_key = f"{sample}_microscopy"
    (preview_key,) = [
        key for key in sdata.images
        if key.startswith(f"{microscopy_key}_") and key.endswith("x_downsampled")
    ]
    save_registration_tissue_preview(
        sdata, sample, args.resolution, preview_key, capture, cytassist_footprint,
        f"{sample}_registration_tissue_mqc.webp",
    )

    segmentation_key = "segmentation_mask"
    if args.segmentation:
        mask = sdata.labels[segmentation_key]
        mask_transform = get_transformation(mask, to_coordinate_system=sample)
        mask_scale = mask_transform.to_affine_matrix(input_axes=("x", "y"), output_axes=("x", "y"))[0, 0]
        mask_downsample_factor = round(mask_scale * args.zarr_downsample_factor)
        image = sdata.images[microscopy_key]
        stats = segmentation_statistics_row(
            sample, sdata, f"{sample}_square_016um", mask.data,
            (image.sizes["y"], image.sizes["x"]), mask_downsample_factor,
            args.zarr_downsample_factor, args.chunk_size,
        )
        save_segmentation_statistics(stats, f"{sample}_segmentation_stats.csv")
        save_result_plots(
            sdata, sample, parse_crop_areas(args.crop_areas), segmentation_key,
            microscopy_key, preview_key, capture.total_bounds,
        )


if __name__ == "__main__":
    main()
