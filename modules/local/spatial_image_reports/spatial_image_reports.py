#!/usr/bin/env python3

"""Render reports from the final SpatialData Zarr."""

import argparse

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import spatialdata as sd
import spatialdata_plot  # noqa: F401
from matplotlib.colors import LinearSegmentedColormap
from spatialdata.models import Image2DModel, get_table_keys
from spatialdata.transformations import get_transformation

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


def save_registration_tissue_preview(sdata, sample_id, resolution, preview_key, output_path):
    """Show the original Space Ranger selection and, where available, registration."""
    capture_key = f"{sample_id}_capture_area"
    cytassist_key = f"{sample_id}_cytassist_image"
    square_key = f"{sample_id}_{resolution}"
    capture = sdata.shapes[capture_key]
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
                "cytassist_magenta", ["#000000", "#ff00ff"]
            )
            overview.pl.render_images(
                preview_key, cmap="gray", colorbar=False,
            ).pl.render_images(
                cytassist_key, cmap=cytassist_cmap, alpha=0.45, colorbar=False,
            ).pl.render_shapes(
                capture_key, fill_alpha=0, outline_color="#39ff14", outline_width=2,
                outline_alpha=1,
            ).pl.show(coordinate_systems=sample_id, ax=axes[0], colorbar=False)
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
                square_key, color="#00d5ff", fill_alpha=0.35, outline_width=0,
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
    parser.add_argument("--sample", required=True)
    parser.add_argument("--resolution", required=True)
    parser.add_argument("--crop-areas")
    parser.add_argument("--segmentation", action="store_true")
    parser.add_argument("--zarr-downsample-factor", type=float, required=True)
    parser.add_argument("--chunk-size", type=int, required=True)
    args = parser.parse_args()

    sdata = sd.read_zarr(args.zarr)
    sample = args.sample
    microscopy_key = f"{sample}_microscopy"
    (preview_key,) = [
        key for key in sdata.images
        if key.startswith(f"{microscopy_key}_") and key.endswith("x_downsampled")
    ]
    save_registration_tissue_preview(
        sdata, sample, args.resolution, preview_key, f"{sample}_registration_tissue_mqc.webp",
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
            microscopy_key, preview_key,
        )


if __name__ == "__main__":
    main()
