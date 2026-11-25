"""Segmentation review figures and statistics from the final SpatialData layers."""

import csv
from dataclasses import dataclass

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import spatialdata as sd
from shapely.geometry import box
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity

PREVIEW_TITLE_STYLE = {
    "fontfamily": "DejaVu Sans",
    "fontsize": 9,
    "fontweight": "semibold",
    "pad": 6,
}


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
    mask,
    microscopy_shape,
    mask_downsample_factor,
    zarr_downsample_factor,
    chunk_size,
):
    """Compute whole-image and Visium-area segmentation statistics."""
    visium_area = visium_bounds(sdata, square_key)
    visium_fullres_area, visium_y_slice, visium_x_slice = fullres_bounds_to_mask_slices(
        visium_area,
        microscopy_shape,
        mask.shape,
        zarr_downsample_factor,
        mask_downsample_factor,
    )

    area_scale = mask_downsample_factor * mask_downsample_factor
    whole_fullres_area_px = microscopy_shape[0] * microscopy_shape[1]
    visium_fullres_area_px = (
        (visium_fullres_area.x1 - visium_fullres_area.x0)
        * (visium_fullres_area.y1 - visium_fullres_area.y0)
    )

    whole_counts = label_counts(mask, chunk_size)
    visium_counts = label_counts(mask, chunk_size, visium_y_slice, visium_x_slice)

    row = {
        "sample": sample_name,
        "segmentation_downsample_factor": mask_downsample_factor,
        "microscopy_height_px": microscopy_shape[0],
        "microscopy_width_px": microscopy_shape[1],
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
    with open(output_path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def save_crop_segmentation_review(
    sdata,
    sample_name,
    crop_areas,
    capture_area,
    microscopy_key,
    microscopy_downsampled_key,
    segmentation_key,
    output_path,
):
    """Save one location/segmentation row for every requested crop area."""
    sdata = sd.SpatialData(
        images={key: sdata.images[key] for key in (microscopy_key, microscopy_downsampled_key)},
        labels={segmentation_key: sdata.labels[segmentation_key]},
    )
    n_crops = len(crop_areas)
    fig, axes = plt.subplots(
        n_crops,
        2,
        figsize=(12, max(4, 4 * n_crops)),
        squeeze=False,
        facecolor="white",
    )

    capture_min, capture_max = capture_area.bounds_for_query
    capture_sdata = sd.bounding_box_query(
        sdata,
        axes=("y", "x"),
        min_coordinate=capture_min,
        max_coordinate=capture_max,
        target_coordinate_system=sample_name,
    )

    for index, crop_area in enumerate(crop_areas):
        selected_crop_key = "selected_crop_for_plotting"
        capture_sdata.shapes[selected_crop_key] = ShapesModel.parse(
            gpd.GeoDataFrame(
                geometry=[box(crop_area.x0, crop_area.y0, crop_area.x1, crop_area.y1)]
            ),
            transformations={sample_name: Identity()},
        )
        capture_sdata.pl.render_images(
            microscopy_downsampled_key,
            cmap="gray",
            colorbar=False,
        ).pl.render_shapes(
            selected_crop_key,
            color="#ff1744",
            fill_alpha=0,
            outline_width=2,
            outline_alpha=1,
            colorbar=False,
        ).pl.show(
            coordinate_systems=sample_name,
            ax=axes[index, 0],
            colorbar=False,
            legend_loc=None,
        )
        axes[index, 0].set_title(
            f"Crop {index + 1} — location",
            **PREVIEW_TITLE_STYLE,
        )
        axes[index, 0].set_axis_off()

        min_coordinate, max_coordinate = crop_area.bounds_for_query
        crop_sdata = sd.bounding_box_query(
            sdata,
            axes=("y", "x"),
            min_coordinate=min_coordinate,
            max_coordinate=max_coordinate,
            target_coordinate_system=sample_name,
        )

        crop_sdata.pl.render_images(
            microscopy_key,
            cmap="gray",
            colorbar=False,
        ).pl.render_labels(
            segmentation_key,
            fill_alpha=0,
            outline_alpha=1,
            outline_color="#ff1744",
            na_color="#ff1744",
        ).pl.show(
            coordinate_systems=sample_name,
            ax=axes[index, 1],
            colorbar=False,
        )
        axes[index, 1].set_title(
            f"Crop {index + 1} — segmentation",
            **PREVIEW_TITLE_STYLE,
        )
        axes[index, 1].set_axis_off()

    fig.subplots_adjust(
        left=0.04,
        right=0.96,
        top=0.94,
        bottom=0.04,
        wspace=0.08,
        hspace=0.22,
    )
    fig.savefig(output_path, dpi=100, facecolor="white")
    plt.close(fig)


def save_result_plots(
    sdata,
    sample_name,
    crop_areas,
    segmentation_key,
    microscopy_key,
    microscopy_downsampled_key,
):
    """Save the combined crop-location and segmentation-review preview."""
    capture_key = f"{sample_name}_capture_area"

    if capture_key not in sdata:
        raise ValueError(f"Required SpatialData element '{capture_key}' not found")

    minx, miny, maxx, maxy = sdata.shapes[capture_key].total_bounds
    capture_area = CropArea(x0=minx, y0=miny, x1=maxx, y1=maxy)
    if crop_areas is None:
        crop_areas = [capture_area]

    save_crop_segmentation_review(
        sdata,
        sample_name,
        crop_areas,
        capture_area,
        microscopy_key,
        microscopy_downsampled_key,
        segmentation_key,
        f"{sample_name}_crop_segmentation_mqc.png",
    )
