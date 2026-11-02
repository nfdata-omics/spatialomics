#!/usr/bin/env python3
"""Import crop-local ROI annotations into a SpatialData zarr store."""

import argparse
import csv
import json
import shutil
from pathlib import Path

import geopandas as gpd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from shapely import affinity
from shapely.geometry import MultiPolygon, Polygon, shape
import spatialdata as sd
from spatialdata.models import ShapesModel
from spatialdata.transformations import Identity
import tifffile


REQUIRED_BOUNDS = {
    "crop_x0",
    "crop_y0",
    "output_width",
    "output_height",
    "scale_x",
    "scale_y",
    "zarr_downsample_factor",
}


def read_bounds(path):
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1:
        raise ValueError(f"Expected exactly one bounds row in {path}, found {len(rows)}")
    missing = REQUIRED_BOUNDS.difference(rows[0])
    if missing:
        raise ValueError(f"Missing bounds columns in {path}: {', '.join(sorted(missing))}")
    bounds = {key: float(rows[0][key]) for key in REQUIRED_BOUNDS}
    if bounds["zarr_downsample_factor"] <= 0:
        raise ValueError("zarr_downsample_factor must be greater than zero")
    return bounds


def json_property(value):
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return json.dumps(value, sort_keys=True, separators=(",", ":"))


def read_annotations(path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("type") != "FeatureCollection":
        raise ValueError("ROI annotations must be a GeoJSON FeatureCollection")

    records = []
    geometries = []
    for index, feature in enumerate(payload.get("features", []), start=1):
        geometry_data = feature.get("geometry")
        if geometry_data is None:
            continue
        geometry = shape(geometry_data)
        if not geometry.is_valid:
            geometry = geometry.buffer(0)
        if geometry.is_empty:
            continue
        if not isinstance(geometry, (Polygon, MultiPolygon)):
            raise ValueError(
                f"Unsupported ROI geometry '{geometry.geom_type}' in feature {index}; "
                "only Polygon and MultiPolygon are supported"
            )

        properties = dict(feature.get("properties") or {})
        classification = properties.get("classification")
        if isinstance(classification, dict):
            properties["classification_name"] = classification.get("name")
            color = classification.get("color")
            if isinstance(color, list) and len(color) >= 3:
                properties["classification_color"] = "#{:02x}{:02x}{:02x}".format(
                    *[max(0, min(255, int(channel))) for channel in color[:3]]
                )
        properties = {key: json_property(value) for key, value in properties.items()}
        properties["roi_id"] = str(feature.get("id") or f"roi_{index}")
        records.append(properties)
        geometries.append(geometry)

    if not geometries:
        raise ValueError(f"No polygon ROI annotations found in {path}")
    return gpd.GeoDataFrame(records, geometry=geometries, crs=None)


def transform_annotations(annotations, bounds):
    factor = bounds["zarr_downsample_factor"]
    x_scale = bounds["scale_x"] / factor
    y_scale = bounds["scale_y"] / factor
    x_offset = bounds["crop_x0"] / factor
    y_offset = bounds["crop_y0"] / factor
    transformed = annotations.copy()
    transformed.geometry = transformed.geometry.map(
        lambda geometry: affinity.translate(
            affinity.scale(geometry, xfact=x_scale, yfact=y_scale, origin=(0, 0)),
            xoff=x_offset,
            yoff=y_offset,
        )
    )
    return transformed


def polygon_parts(geometry):
    return geometry.geoms if isinstance(geometry, MultiPolygon) else [geometry]


def save_plot(annotations, image_path, output_path):
    image = tifffile.imread(image_path)
    fig, axis = plt.subplots(figsize=(10, 10))
    axis.imshow(image, cmap="gray" if image.ndim == 2 else None)

    legend = {}
    for _, row in annotations.iterrows():
        label = row.get("classification_name")
        if not isinstance(label, str) or not label:
            label = "Unclassified"
        color = row.get("classification_color")
        if not isinstance(color, str) or not is_color_like(color):
            color = "#ffcc00"
        legend[label] = color
        for polygon in polygon_parts(row.geometry):
            x, y = polygon.exterior.xy
            axis.fill(x, y, facecolor=color, edgecolor=color, alpha=0.25, linewidth=1.5)
            for interior in polygon.interiors:
                x, y = interior.xy
                axis.fill(x, y, facecolor="none", edgecolor=color, linewidth=1.0)

    handles = [
        Line2D([0], [0], color=color, linewidth=3, label=label)
        for label, color in sorted(legend.items())
    ]
    if handles:
        axis.legend(handles=handles, loc="upper right", framealpha=0.9)
    axis.set_xlim(-0.5, image.shape[1] - 0.5)
    axis.set_ylim(image.shape[0] - 0.5, -0.5)
    axis.set_title("ROI annotations")
    axis.set_axis_off()
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def run(args):
    bounds = read_bounds(args.bounds)
    local_annotations = read_annotations(args.roi_annotations)
    minx, miny, maxx, maxy = local_annotations.total_bounds
    if minx < 0 or miny < 0 or maxx > bounds["output_width"] or maxy > bounds["output_height"]:
        print(
            "Warning: some ROI coordinates fall outside the capture-area image bounds: "
            f"ROI={minx:g}:{miny:g}:{maxx:g}:{maxy:g}, "
            f"image=0:0:{bounds['output_width']:g}:{bounds['output_height']:g}"
        )

    save_plot(local_annotations, args.capture_area_image, args.output_plot)
    annotations = transform_annotations(local_annotations, bounds)

    sdata = sd.read_zarr(args.input_zarr)
    if args.sample_name not in sdata.coordinate_systems:
        raise ValueError(f"Coordinate system '{args.sample_name}' not found in input zarr")
    shape_key = f"{args.sample_name}_roi_annotations"
    if shape_key in sdata.shapes:
        raise ValueError(f"SpatialData shape element '{shape_key}' already exists")
    sdata.shapes[shape_key] = ShapesModel.parse(
        annotations,
        transformations={args.sample_name: Identity()},
    )

    if args.output_zarr.exists():
        if args.output_zarr.is_dir():
            shutil.rmtree(args.output_zarr)
        else:
            args.output_zarr.unlink()
    sdata.write(args.output_zarr)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-zarr", required=True, type=Path)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--roi-annotations", required=True, type=Path)
    parser.add_argument("--bounds", required=True, type=Path)
    parser.add_argument("--capture-area-image", required=True, type=Path)
    parser.add_argument("--output-zarr", required=True, type=Path)
    parser.add_argument("--output-plot", required=True, type=Path)
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
