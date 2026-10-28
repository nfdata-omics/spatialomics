#!/usr/bin/env python3
"""Segment nuclei with thresholding and marker-controlled watershed."""

import argparse

import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage import color, exposure, feature, filters, measure, morphology, segmentation, util


def parse_args():
    parser = argparse.ArgumentParser(
        description="Segment nuclei in a Visium microscopy region using scikit-image watershed"
    )
    parser.add_argument("--input-image", required=True, help="Input grayscale or RGB TIFF")
    parser.add_argument("--visium-bounds", required=True, help="TSV containing x0, y0, x1, y1")
    parser.add_argument("--output-mask", required=True, help="Output instance-label TIFF")
    parser.add_argument("--downsample-factor", type=int, default=1)
    parser.add_argument("--crop-padding", type=int, default=256)
    parser.add_argument(
        "--threshold-scale",
        type=float,
        default=1.0,
        help="Multiplier applied to the automatically estimated Otsu threshold",
    )
    parser.add_argument(
        "--gaussian-sigma",
        type=float,
        default=2.0,
        help="Gaussian sigma in full-resolution pixels",
    )
    parser.add_argument(
        "--min-nucleus-distance",
        type=float,
        default=12.0,
        help="Minimum distance between watershed seeds in full-resolution pixels",
    )
    parser.add_argument(
        "--min-nucleus-area",
        type=float,
        default=80.0,
        help="Minimum nucleus area in full-resolution pixels squared",
    )
    parser.add_argument(
        "--max-nucleus-area",
        type=float,
        default=5000.0,
        help="Maximum nucleus area in full-resolution pixels squared; use 0 to disable",
    )
    parser.add_argument(
        "--hole-area",
        type=float,
        default=80.0,
        help="Maximum foreground-hole area to fill in full-resolution pixels squared",
    )
    parser.add_argument(
        "--keep-border-objects",
        action="store_true",
        help="Keep objects touching the boundary of the segmented crop",
    )
    return parser.parse_args()


def read_visium_bounds(path):
    with open(path, encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]

    if len(lines) < 2:
        raise ValueError(f"Visium bounds file must contain a header and one data row: {path}")

    row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    missing = [field for field in ("x0", "y0", "x1", "y1") if field not in row]
    if missing:
        raise ValueError(f"Visium bounds file is missing: {', '.join(missing)}")

    bounds = tuple(int(float(row[field])) for field in ("x0", "y0", "x1", "y1"))
    if bounds[0] >= bounds[2] or bounds[1] >= bounds[3]:
        raise ValueError(f"Invalid Visium bounds: {bounds}")
    return bounds


def downsampled_size(size, factor):
    return 0 if size == 0 else ((size - 1) // factor) + 1


def crop_bounds(bounds, image_shape, padding, factor):
    if padding < 0:
        raise ValueError("--crop-padding must be >= 0")
    if factor <= 0:
        raise ValueError("--downsample-factor must be > 0")

    height, width = image_shape[:2]
    x0, y0, x1, y1 = bounds
    x0 = max(0, x0 - padding)
    y0 = max(0, y0 - padding)
    x1 = min(width, x1 + padding)
    y1 = min(height, y1 + padding)
    x0 = (x0 // factor) * factor
    y0 = (y0 // factor) * factor

    if x0 >= x1 or y0 >= y1:
        raise ValueError("The padded Visium bounds do not overlap the microscopy image")
    return x0, y0, x1, y1


def segment_nuclei(image, args):
    factor = args.downsample_factor
    image = np.squeeze(image)

    if image.ndim == 3 and image.shape[0] in (3, 4) and image.shape[-1] not in (3, 4):
        image = np.moveaxis(image, 0, -1)

    if image.ndim == 3:
        if image.shape[-1] < 3:
            image = image[..., 0]
        else:
            rgb = util.img_as_float(image[..., :3])
            nuclei_signal = color.rgb2hed(rgb)[..., 0]
            print("Nuclear signal: hematoxylin channel from RGB H&E deconvolution")
    if image.ndim == 2:
        gray = exposure.rescale_intensity(util.img_as_float(image), out_range=(0.0, 1.0))
        nuclei_signal = 1.0 - gray
        print("Nuclear signal: inverted grayscale image")
    elif image.ndim != 3:
        raise ValueError(f"Expected a grayscale or RGB image, got shape {image.shape}")

    sigma = max(0.0, args.gaussian_sigma / factor)
    smoothed = filters.gaussian(nuclei_signal, sigma=sigma, preserve_range=True)
    threshold = filters.threshold_otsu(smoothed) * args.threshold_scale
    foreground = smoothed > threshold

    min_area = max(1, round(args.min_nucleus_area / factor**2))
    max_area = max(1, round(args.max_nucleus_area / factor**2)) if args.max_nucleus_area > 0 else None
    hole_area = max(1, round(args.hole_area / factor**2))
    min_distance = max(1, round(args.min_nucleus_distance / factor))

    foreground = morphology.remove_small_holes(foreground, area_threshold=hole_area)
    foreground = morphology.remove_small_objects(foreground, min_size=min_area)
    distance = ndi.distance_transform_edt(foreground)
    coordinates = feature.peak_local_max(
        distance,
        labels=foreground,
        min_distance=min_distance,
        exclude_border=False,
    )

    markers = np.zeros(foreground.shape, dtype=np.int32)
    if len(coordinates):
        markers[tuple(coordinates.T)] = np.arange(1, len(coordinates) + 1)
    else:
        markers, _ = ndi.label(foreground)

    labels = segmentation.watershed(-distance, markers, mask=foreground)
    sizes = np.bincount(labels.ravel())
    rejected = sizes < min_area
    rejected[0] = False
    if max_area is not None:
        rejected |= sizes > max_area
        rejected[0] = False
    labels[rejected[labels]] = 0

    if not args.keep_border_objects:
        labels = segmentation.clear_border(labels)

    labels, _, _ = segmentation.relabel_sequential(labels)
    print(
        f"Otsu threshold: {threshold:.6g}; seeds: {len(coordinates)}; "
        f"retained nuclei: {int(labels.max())}"
    )
    return labels.astype(np.uint32, copy=False)


def main():
    args = parse_args()

    try:
        image = tifffile.memmap(args.input_image, mode="r")
    except ValueError:
        image = tifffile.imread(args.input_image)

    image = np.squeeze(image)
    if image.ndim == 3 and image.shape[0] in (3, 4) and image.shape[-1] not in (3, 4):
        image = np.moveaxis(image, 0, -1)
    if image.ndim not in (2, 3):
        raise ValueError(f"Expected a 2D or RGB TIFF, got shape {image.shape}")

    full_shape = image.shape[:2]
    bounds = crop_bounds(
        read_visium_bounds(args.visium_bounds),
        full_shape,
        args.crop_padding,
        args.downsample_factor,
    )
    x0, y0, x1, y1 = bounds
    crop = np.asarray(image[y0:y1:args.downsample_factor, x0:x1:args.downsample_factor])
    print(f"Segmenting crop {x0}:{y0}:{x1}:{y1}; downsampled shape: {crop.shape[:2]}")
    crop_mask = segment_nuclei(crop, args)

    output_shape = tuple(downsampled_size(size, args.downsample_factor) for size in full_shape)
    output = np.zeros(output_shape, dtype=np.uint32)
    mask_y0 = y0 // args.downsample_factor
    mask_x0 = x0 // args.downsample_factor
    output[
        mask_y0 : mask_y0 + crop_mask.shape[0],
        mask_x0 : mask_x0 + crop_mask.shape[1],
    ] = crop_mask

    tifffile.imwrite(
        args.output_mask,
        output,
        bigtiff=True,
        compression=None,
        photometric="minisblack",
        contiguous=True,
        metadata={
            "axes": "YX",
            "SegmentationMethod": "scikit-image thresholding and watershed",
            "DownsampleFactor": args.downsample_factor,
            "VisiumCropFullRes": list(bounds),
            "CropPaddingFullRes": args.crop_padding,
        },
    )


if __name__ == "__main__":
    main()
