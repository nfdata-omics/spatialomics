#!/usr/bin/env python3
"""Refine an existing CytAssist-to-microscopy transform with a bounded local affine step."""

import argparse
import csv
import importlib
import importlib.metadata
import json
import math
import sys
import warnings
from pathlib import Path

import cv2
import numpy as np
import tifffile
import yaml
from PIL import Image
from scipy import optimize
from skimage import filters, morphology
from skimage.metrics import normalized_mutual_information


warnings.filterwarnings("ignore", category=FutureWarning)


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


def squeeze_supported_image(array, label):
    """Squeeze singleton dimensions while keeping common 2D/RGB/channel-first images."""
    image = np.squeeze(array)
    if image.ndim == 2:
        return image
    if image.ndim == 3 and (image.shape[-1] <= 10 or image.shape[0] <= 10):
        return image
    raise ValueError(f"Expected {label} to be 2D or have a small channel axis, got shape {image.shape}")


def spatial_shape(image):
    """Return YX spatial shape for a supported image array."""
    if image.ndim == 2:
        return image.shape
    if image.shape[-1] <= 10:
        return image.shape[:2]
    if image.shape[0] <= 10:
        return image.shape[1:]
    raise ValueError(f"Cannot infer spatial axes for shape {image.shape}")


def read_image_array(path):
    """Read an image, preferring TIFF memmap when possible."""
    try:
        return tifffile.memmap(path)
    except Exception:
        return tifffile.imread(path)


def downsample_to_preview(image, max_dim):
    """Create a 2D preview image with approximately max_dim pixels on the longest side."""
    image = squeeze_supported_image(image, "image")
    height, width = spatial_shape(image)
    stride = max(1, int(math.ceil(max(height, width) / max_dim)))

    color_preview = None
    if image.ndim == 2:
        preview = image[::stride, ::stride]
    elif image.shape[-1] <= 10:
        color_preview = np.asarray(image[::stride, ::stride, :3])
        preview = color_preview.astype(np.float32, copy=False).mean(axis=-1)
    else:
        color_preview = np.asarray(np.moveaxis(image[:3, ::stride, ::stride], 0, -1))
        preview = color_preview.astype(np.float32, copy=False).mean(axis=-1)

    return np.asarray(preview), color_preview, (height, width), stride


def normalize_float(image):
    """Robustly normalize an image to float32 in [0, 1]."""
    data = np.asarray(image, dtype=np.float32)
    data = np.nan_to_num(data, copy=False)
    low, high = np.percentile(data, [1, 99])
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        low = float(np.min(data))
        high = float(np.max(data))
    if high <= low:
        return np.zeros(data.shape, dtype=np.float32)
    return np.clip((data - low) / (high - low), 0, 1).astype(np.float32, copy=False)


def to_uint8(image):
    """Convert a normalized image to uint8."""
    return np.clip(image * 255, 0, 255).astype(np.uint8)


def normalize_rgb_uint8(image):
    """Robustly normalize an RGB preview to uint8 for color masking."""
    if image is None:
        return None
    rgb = np.asarray(image)
    if rgb.ndim != 3 or rgb.shape[-1] < 3:
        return None
    rgb = rgb[..., :3].astype(np.float32, copy=False)
    out = np.zeros(rgb.shape, dtype=np.uint8)
    for channel in range(3):
        values = rgb[..., channel]
        low, high = np.percentile(values, [1, 99])
        if not np.isfinite(low) or not np.isfinite(high) or high <= low:
            low = float(np.min(values))
            high = float(np.max(values))
        if high > low:
            out[..., channel] = np.clip((values - low) / (high - low) * 255, 0, 255).astype(np.uint8)
    return out


def clean_mask(mask):
    """Remove small mask speckles and holes."""
    min_size = max(16, int(mask.size * 0.0001))
    cleaned = morphology.remove_small_objects(mask.astype(bool), min_size=min_size)
    cleaned = morphology.remove_small_holes(cleaned, area_threshold=min_size)
    return cleaned.astype(bool)


def foreground_mask(image):
    """Estimate a tissue/content mask from a normalized preview image."""
    if np.max(image) <= np.min(image):
        return np.zeros(image.shape, dtype=bool)

    threshold = filters.threshold_otsu(image)
    high = image > threshold
    low = image <= threshold

    height, width = image.shape
    border_y = max(1, int(round(height * 0.05)))
    border_x = max(1, int(round(width * 0.05)))
    border = np.zeros(image.shape, dtype=bool)
    border[:border_y, :] = True
    border[-border_y:, :] = True
    border[:, :border_x] = True
    border[:, -border_x:] = True

    border_high_fraction = float(high[border].mean())
    preferred = low if border_high_fraction >= 0.5 else high
    fallback = high if border_high_fraction >= 0.5 else low
    for mask in (preferred, fallback):
        cleaned = clean_mask(mask)
        fraction = float(cleaned.mean())
        if 0.005 <= fraction <= 0.95:
            return cleaned
    return clean_mask(preferred)


def color_tissue_mask(color_preview, exclude_green=False):
    """Estimate tissue mask from RGB saturation/value, optionally excluding green fiducials."""
    rgb = normalize_rgb_uint8(color_preview)
    if rgb is None:
        return None

    hsv = cv2.cvtColor(rgb, cv2.COLOR_RGB2HSV)
    hue = hsv[..., 0]
    saturation = hsv[..., 1]
    value = hsv[..., 2]

    mask = (saturation > 18) & (value > 20) & (value < 250)
    if exclude_green:
        green = (hue >= 35) & (hue <= 95) & (saturation > 35)
        mask &= ~green

    mask = clean_mask(mask)
    fraction = float(mask.mean())
    if fraction < 0.005 or fraction > 0.90:
        return None
    return mask


def tissue_mask(gray_preview, color_preview=None, exclude_green=False):
    """Estimate tissue/content mask, preferring color information when available."""
    color_mask = color_tissue_mask(color_preview, exclude_green=exclude_green)
    if color_mask is not None:
        return color_mask
    return foreground_mask(gray_preview)


def suppress_border(mask, border_fraction):
    """Suppress mask content close to image borders."""
    if border_fraction <= 0:
        return mask
    height, width = mask.shape
    border_y = int(round(height * border_fraction))
    border_x = int(round(width * border_fraction))
    trimmed = mask.copy()
    if border_y > 0:
        trimmed[:border_y, :] = False
        trimmed[-border_y:, :] = False
    if border_x > 0:
        trimmed[:, :border_x] = False
        trimmed[:, -border_x:] = False
    return trimmed


def warp_image(image, matrix, output_shape, interpolation=cv2.INTER_LINEAR):
    """Warp an image into output_shape using an affine matrix."""
    fixed_h, fixed_w = output_shape
    return cv2.warpAffine(
        image.astype(np.float32),
        matrix,
        (fixed_w, fixed_h),
        flags=interpolation,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=0,
    )


def matrix_3x3(matrix):
    """Return a 3x3 homogeneous matrix from a 2x3 affine matrix."""
    return np.vstack([matrix, [0.0, 0.0, 1.0]])


def full_to_preview_matrix(full_matrix, cytassist_stride, microscopy_stride):
    """Convert full-resolution CytAssist-to-microscopy matrix into preview coordinates."""
    fixed_inv = np.diag([1.0 / microscopy_stride, 1.0 / microscopy_stride, 1.0])
    moving = np.diag([cytassist_stride, cytassist_stride, 1.0])
    return fixed_inv @ full_matrix @ moving


def preview_to_full_matrix(preview_matrix, cytassist_stride, microscopy_stride):
    """Convert preview-coordinate matrix to full-resolution coordinate matrix."""
    fixed = np.diag([microscopy_stride, microscopy_stride, 1.0])
    moving_inv = np.diag([1.0 / cytassist_stride, 1.0 / cytassist_stride, 1.0])
    return fixed @ matrix_3x3(preview_matrix) @ moving_inv


def transformed_center(matrix, moving_shape):
    """Transform moving-image center into fixed-image coordinates."""
    moving_h, moving_w = moving_shape
    center = np.array([moving_w / 2.0, moving_h / 2.0, 1.0], dtype=np.float64)
    transformed = matrix_3x3(matrix) @ center
    return transformed[:2] / transformed[2]


def mask_boundary(mask):
    """Return a thin binary boundary image from a mask."""
    eroded = morphology.erosion(mask)
    return (mask & ~eroded).astype(np.float32)


def normalize_nonzero(image):
    """Normalize a non-negative image to [0, 1]."""
    max_value = float(np.max(image))
    if max_value <= 0:
        return np.zeros(image.shape, dtype=np.float32)
    return (image / max_value).astype(np.float32, copy=False)


def gradient_signal(image):
    """Return normalized gradient magnitude."""
    image = image.astype(np.float32, copy=False)
    grad_x = cv2.Sobel(image, cv2.CV_32F, 1, 0, ksize=3)
    grad_y = cv2.Sobel(image, cv2.CV_32F, 0, 1, ksize=3)
    return normalize_nonzero(cv2.magnitude(grad_x, grad_y))


def refinement_signal(image, mask):
    """Create a local refinement signal from tissue texture and mask boundaries."""
    gradient = gradient_signal(image)
    boundary = cv2.GaussianBlur(mask_boundary(mask), (5, 5), 0)
    tissue_intensity = image * mask.astype(np.float32)
    signal = 0.55 * gradient + 0.30 * boundary + 0.15 * normalize_nonzero(tissue_intensity)
    signal *= morphology.binary_dilation(mask, morphology.disk(5)).astype(np.float32)
    return normalize_nonzero(signal)


def bbox_from_mask(mask, padding, shape):
    """Return y0, y1, x0, x1 for a mask bounding box plus padding."""
    coords = np.argwhere(mask)
    if coords.size == 0:
        return 0, shape[0], 0, shape[1]
    y0, x0 = coords.min(axis=0)
    y1, x1 = coords.max(axis=0) + 1
    y0 = max(0, int(y0) - padding)
    x0 = max(0, int(x0) - padding)
    y1 = min(shape[0], int(y1) + padding)
    x1 = min(shape[1], int(x1) + padding)
    return y0, y1, x0, x1


def crop_matrix(matrix, y0, x0):
    """Convert full fixed-preview coordinates to fixed-crop coordinates."""
    cropped = matrix.copy()
    cropped[0, 2] -= x0
    cropped[1, 2] -= y0
    return cropped


def normalized_product_score(fixed_signal, moving_signal, matrix):
    """Return normalized image product after warping moving signal into fixed space."""
    warped = warp_image(moving_signal, matrix, fixed_signal.shape)
    denominator = math.sqrt(float(np.sum(fixed_signal * fixed_signal)) * float(np.sum(warped * warped)))
    if denominator <= 0:
        return 0.0
    return float(np.sum(fixed_signal * warped) / denominator)


def compute_metrics(fixed, moving, matrix, fixed_mask, moving_mask, sample, label):
    """Compute registration QC metrics on preview images."""
    warped = warp_image(moving, matrix, fixed.shape)
    warped_mask = warp_image(
        moving_mask.astype(np.float32),
        matrix,
        fixed.shape,
        interpolation=cv2.INTER_NEAREST,
    ) > 0.5

    union = fixed_mask | warped_mask
    intersection = fixed_mask & warped_mask
    mask_iou = float(intersection.sum() / union.sum()) if union.any() else 0.0
    moving_overlap = float(intersection.sum() / warped_mask.sum()) if warped_mask.any() else 0.0
    fixed_overlap = float(intersection.sum() / fixed_mask.sum()) if fixed_mask.any() else 0.0

    overlap = intersection & (warped > 0)
    if overlap.sum() > 10 and np.std(fixed[overlap]) > 0 and np.std(warped[overlap]) > 0:
        ncc = float(np.corrcoef(fixed[overlap], warped[overlap])[0, 1])
    else:
        ncc = 0.0

    try:
        nmi = float(normalized_mutual_information(to_uint8(fixed), to_uint8(warped)))
    except Exception:
        nmi = 0.0

    a, b, tx = matrix[0]
    c, d, ty = matrix[1]
    scale_x = math.sqrt(a * a + c * c)
    scale_y = math.sqrt(b * b + d * d)
    determinant = a * d - b * c
    rotation = math.degrees(math.atan2(c, a))

    return {
        "sample": sample,
        "label": label,
        "mask_iou": mask_iou,
        "moving_mask_overlap_fraction": moving_overlap,
        "fixed_mask_overlap_fraction": fixed_overlap,
        "normalized_cross_correlation": ncc,
        "normalized_mutual_information": nmi,
        "preview_scale_x": scale_x,
        "preview_scale_y": scale_y,
        "preview_rotation_degrees": rotation,
        "preview_translation_x": float(tx),
        "preview_translation_y": float(ty),
        "preview_determinant": float(determinant),
    }, warped


def score_metrics(metrics, signal_score):
    """Score refinement while keeping mask overlap as a guardrail."""
    return (
        1.5 * signal_score
        + 1.0 * metrics["mask_iou"]
        + 0.5 * metrics["moving_mask_overlap_fraction"]
        + 0.25 * metrics["fixed_mask_overlap_fraction"]
        + 0.1 * max(metrics["normalized_cross_correlation"], 0.0)
    )


def linear_scale(matrix):
    """Return the mean linear scale of a 2x3 affine matrix."""
    a, b = matrix[0, 0], matrix[0, 1]
    c, d = matrix[1, 0], matrix[1, 1]
    return (math.sqrt(a * a + c * c) + math.sqrt(b * b + d * d)) / 2.0


def refine_affine(initial_matrix, fixed_signal, moving_signal, moving_shape, args):
    """Refine a transform with a bounded affine delta around its current center."""
    center = transformed_center(initial_matrix, moving_shape)
    origin_to_center = np.array([[1.0, 0.0, center[0]], [0.0, 1.0, center[1]], [0.0, 0.0, 1.0]])
    center_to_origin = np.array([[1.0, 0.0, -center[0]], [0.0, 1.0, -center[1]], [0.0, 0.0, 1.0]])
    initial_3x3 = matrix_3x3(initial_matrix)

    max_rotation = math.radians(args.max_rotation_degrees)
    min_log_scale = math.log(max(0.05, 1.0 - args.max_scale_change))
    max_log_scale = math.log(1.0 + args.max_scale_change)
    bounds = [
        (-max_rotation, max_rotation),
        (min_log_scale, max_log_scale),
        (min_log_scale, max_log_scale),
        (-args.max_shear, args.max_shear),
        (-args.max_shear, args.max_shear),
        (-args.max_translation, args.max_translation),
        (-args.max_translation, args.max_translation),
    ]

    def delta_matrix(params):
        theta, log_scale_x, log_scale_y, shear_x, shear_y, tx, ty = params
        cos_theta = math.cos(theta)
        sin_theta = math.sin(theta)
        rotation = np.array(
            [[cos_theta, -sin_theta, 0.0], [sin_theta, cos_theta, 0.0], [0.0, 0.0, 1.0]],
            dtype=np.float64,
        )
        affine_delta = np.array(
            [
                [math.exp(log_scale_x), shear_x, tx],
                [shear_y, math.exp(log_scale_y), ty],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        return origin_to_center @ affine_delta @ rotation @ center_to_origin @ initial_3x3

    def objective(params):
        return -normalized_product_score(fixed_signal, moving_signal, delta_matrix(params)[:2])

    result = optimize.minimize(
        objective,
        np.zeros(7, dtype=np.float64),
        method="Powell",
        bounds=bounds,
        options={"maxiter": args.iterations, "xtol": args.epsilon, "ftol": args.epsilon, "disp": False},
    )

    refined = delta_matrix(result.x)[:2]
    return refined, {
        "optimizer_success": bool(result.success),
        "optimizer_iterations": int(getattr(result, "nit", 0)),
        "objective_score": float(-result.fun),
        "delta_rotation_degrees": float(math.degrees(result.x[0])),
        "delta_scale_x": float(math.exp(result.x[1])),
        "delta_scale_y": float(math.exp(result.x[2])),
        "delta_shear_x": float(result.x[3]),
        "delta_shear_y": float(result.x[4]),
        "delta_translation_x": float(result.x[5]),
        "delta_translation_y": float(result.x[6]),
    }


def acceptance_status(initial_matrix, refined_matrix, initial_metrics, refined_metrics, initial_score, refined_score, args):
    """Decide whether to accept a refined transform."""
    initial_det = float(np.linalg.det(initial_matrix[:, :2]))
    refined_det = float(np.linalg.det(refined_matrix[:, :2]))
    if initial_det == 0 or refined_det == 0 or np.sign(initial_det) != np.sign(refined_det):
        return "rejected_determinant_sign"

    initial_scale = linear_scale(initial_matrix)
    refined_scale = linear_scale(refined_matrix)
    if initial_scale <= 0:
        return "rejected_invalid_initial_scale"
    if abs(refined_scale / initial_scale - 1.0) > args.max_scale_change:
        return "rejected_scale_change"

    if refined_metrics["mask_iou"] + args.max_mask_iou_drop < initial_metrics["mask_iou"]:
        return "rejected_mask_iou_drop"
    if refined_metrics["moving_mask_overlap_fraction"] + args.max_moving_overlap_drop < initial_metrics["moving_mask_overlap_fraction"]:
        return "rejected_moving_overlap_drop"
    if refined_score + args.min_score_improvement < initial_score:
        return "rejected_score_drop"
    return "accepted"


def write_overlay(path, fixed, warped):
    """Write a green/magenta registration overlay PNG."""
    fixed_u8 = to_uint8(fixed)
    warped_u8 = to_uint8(warped)
    rgb = np.zeros((*fixed.shape, 3), dtype=np.uint8)
    rgb[..., 1] = fixed_u8
    rgb[..., 0] = warped_u8
    rgb[..., 2] = warped_u8
    Image.fromarray(rgb).save(path)


def write_metrics(path, row):
    """Write a single-row metrics CSV."""
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
        writer.writeheader()
        writer.writerow(row)


def run(args):
    """Run local alignment refinement."""
    cytassist_raw, cytassist_color, cytassist_shape, cytassist_stride = downsample_to_preview(
        read_image_array(args.cytassist_image),
        args.max_preview_dim,
    )
    microscopy_raw, microscopy_color, microscopy_shape, microscopy_stride = downsample_to_preview(
        read_image_array(args.microscopy_tif),
        args.max_preview_dim,
    )

    moving = normalize_float(cytassist_raw)
    fixed = normalize_float(microscopy_raw)
    moving_mask = suppress_border(
        tissue_mask(moving, cytassist_color, exclude_green=args.exclude_green_fiducials),
        args.cytassist_border_fraction,
    )
    if not moving_mask.any():
        moving_mask = tissue_mask(moving, cytassist_color, exclude_green=False)
    fixed_mask = tissue_mask(fixed, microscopy_color, exclude_green=False)

    initial_payload = json.loads(args.initial_transform.read_text(encoding="utf-8"))
    initial_full = np.array(initial_payload["full_resolution_transform_3x3"], dtype=np.float64)
    initial_preview = full_to_preview_matrix(initial_full, cytassist_stride, microscopy_stride)[:2]

    fixed_signal = refinement_signal(fixed, fixed_mask)
    moving_signal = refinement_signal(moving, moving_mask)

    warped_initial_mask = warp_image(
        moving_mask.astype(np.float32),
        initial_preview,
        fixed.shape,
        interpolation=cv2.INTER_NEAREST,
    ) > 0.5
    overlap = fixed_mask & warped_initial_mask
    roi_mask = overlap if int(overlap.sum()) >= args.min_overlap_pixels else (fixed_mask | warped_initial_mask)
    y0, y1, x0, x1 = bbox_from_mask(roi_mask, args.roi_padding, fixed.shape)

    fixed_signal_roi = fixed_signal[y0:y1, x0:x1]
    initial_preview_roi = crop_matrix(initial_preview, y0, x0)
    initial_signal_score = normalized_product_score(fixed_signal_roi, moving_signal, initial_preview_roi)

    refined_preview_roi, refine_stats = refine_affine(
        initial_preview_roi,
        fixed_signal_roi,
        moving_signal,
        moving.shape,
        args,
    )
    refined_preview = refined_preview_roi.copy()
    refined_preview[0, 2] += x0
    refined_preview[1, 2] += y0

    initial_metrics, initial_warped = compute_metrics(
        fixed=fixed,
        moving=moving,
        matrix=initial_preview,
        fixed_mask=fixed_mask,
        moving_mask=moving_mask,
        sample=args.sample_name,
        label="initial",
    )
    refined_metrics, refined_warped = compute_metrics(
        fixed=fixed,
        moving=moving,
        matrix=refined_preview,
        fixed_mask=fixed_mask,
        moving_mask=moving_mask,
        sample=args.sample_name,
        label="refined",
    )
    refined_signal_score = normalized_product_score(fixed_signal_roi, moving_signal, crop_matrix(refined_preview, y0, x0))
    initial_score = score_metrics(initial_metrics, initial_signal_score)
    refined_score = score_metrics(refined_metrics, refined_signal_score)

    status = acceptance_status(
        initial_preview,
        refined_preview,
        initial_metrics,
        refined_metrics,
        initial_score,
        refined_score,
        args,
    )
    accepted = status == "accepted"
    output_preview = refined_preview if accepted else initial_preview
    output_full = preview_to_full_matrix(output_preview, cytassist_stride, microscopy_stride)
    output_metrics = refined_metrics if accepted else initial_metrics
    output_warped = refined_warped if accepted else initial_warped

    row = {
        "sample": args.sample_name,
        "status": status,
        "accepted": accepted,
        "roi_y0": y0,
        "roi_y1": y1,
        "roi_x0": x0,
        "roi_x1": x1,
        "initial_signal_score": initial_signal_score,
        "refined_signal_score": refined_signal_score,
        "initial_candidate_score": initial_score,
        "refined_candidate_score": refined_score,
        "initial_mask_iou": initial_metrics["mask_iou"],
        "refined_mask_iou": refined_metrics["mask_iou"],
        "initial_moving_mask_overlap_fraction": initial_metrics["moving_mask_overlap_fraction"],
        "refined_moving_mask_overlap_fraction": refined_metrics["moving_mask_overlap_fraction"],
        "initial_fixed_mask_overlap_fraction": initial_metrics["fixed_mask_overlap_fraction"],
        "refined_fixed_mask_overlap_fraction": refined_metrics["fixed_mask_overlap_fraction"],
        "initial_ncc": initial_metrics["normalized_cross_correlation"],
        "refined_ncc": refined_metrics["normalized_cross_correlation"],
        "output_mask_iou": output_metrics["mask_iou"],
        "output_moving_mask_overlap_fraction": output_metrics["moving_mask_overlap_fraction"],
        "output_fixed_mask_overlap_fraction": output_metrics["fixed_mask_overlap_fraction"],
        **refine_stats,
    }

    payload = {
        "sample": args.sample_name,
        "transform_direction": "cytassist_to_microscopy",
        "coordinate_units": "pixels",
        "method": "local_affine_refinement",
        "status": status,
        "accepted": accepted,
        "initial_transform_path": str(args.initial_transform),
        "preview_transform_3x3": matrix_3x3(output_preview).tolist(),
        "full_resolution_transform_3x3": output_full.tolist(),
        "initial_preview_transform_3x3": matrix_3x3(initial_preview).tolist(),
        "refined_preview_transform_3x3": matrix_3x3(refined_preview).tolist(),
        "cytassist_shape_yx": list(cytassist_shape),
        "microscopy_shape_yx": list(microscopy_shape),
        "cytassist_preview_stride": cytassist_stride,
        "microscopy_preview_stride": microscopy_stride,
        "metrics": row,
    }

    args.output_transform.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    write_metrics(args.output_metrics, row)
    write_overlay(args.output_overlay_before, fixed, initial_warped)
    write_overlay(args.output_overlay_after, fixed, output_warped)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Refine an existing CytAssist-to-microscopy transform")
    parser.add_argument("--cytassist-image", type=Path, help="Input CytAssist image")
    parser.add_argument("--microscopy-tif", type=Path, help="Input microscopy TIFF")
    parser.add_argument("--initial-transform", type=Path, help="Initial transform JSON")
    parser.add_argument("--sample-name", help="Sample name")
    parser.add_argument("--output-transform", type=Path, help="Output refined transform JSON")
    parser.add_argument("--output-metrics", type=Path, help="Output refinement metrics CSV")
    parser.add_argument("--output-overlay-before", type=Path, help="Output before-refinement overlay PNG")
    parser.add_argument("--output-overlay-after", type=Path, help="Output after-refinement overlay PNG")
    parser.add_argument("--max-preview-dim", type=int, default=2048, help="Maximum preview dimension")
    parser.add_argument("--cytassist-border-fraction", type=float, default=0.03, help="CytAssist border fraction to suppress")
    parser.add_argument(
        "--exclude-green-fiducials",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exclude green CytAssist fiducial/frame signal from tissue masking",
    )
    parser.add_argument("--roi-padding", type=int, default=96, help="Preview-pixel padding around overlap ROI")
    parser.add_argument("--min-overlap-pixels", type=int, default=500, help="Minimum overlap pixels before falling back to union ROI")
    parser.add_argument("--max-rotation-degrees", type=float, default=2.0, help="Maximum local rotation")
    parser.add_argument("--max-scale-change", type=float, default=0.025, help="Maximum fractional local scale change")
    parser.add_argument("--max-shear", type=float, default=0.025, help="Maximum local shear")
    parser.add_argument("--max-translation", type=float, default=60.0, help="Maximum local translation in preview pixels")
    parser.add_argument("--iterations", type=int, default=250, help="Maximum optimizer iterations")
    parser.add_argument("--epsilon", type=float, default=1e-6, help="Optimizer convergence epsilon")
    parser.add_argument("--max-mask-iou-drop", type=float, default=0.02, help="Maximum allowed mask IoU drop")
    parser.add_argument(
        "--max-moving-overlap-drop",
        type=float,
        default=0.05,
        help="Maximum allowed moving tissue overlap drop",
    )
    parser.add_argument(
        "--min-score-improvement",
        type=float,
        default=-0.02,
        help="Minimum refined score improvement; negative values allow small score drops",
    )
    parser.add_argument(
        "--versions-dict",
        help="If set, print versions of relevant libraries in YAML format and exit",
    )
    return parser.parse_args()


def main():
    """Run the CLI."""
    args = parse_args()
    if args.versions_dict:
        print(
            versions_yaml(
                args.versions_dict,
                ["opencv-python", "numpy", "Pillow", "scikit-image", "scipy", "tifffile"],
            )
        )
        return

    missing = [
        name
        for name in (
            "cytassist_image",
            "microscopy_tif",
            "initial_transform",
            "sample_name",
            "output_transform",
            "output_metrics",
            "output_overlay_before",
            "output_overlay_after",
        )
        if getattr(args, name) is None
    ]
    if missing:
        raise ValueError(f"Missing required arguments: {', '.join(missing)}")
    if args.max_preview_dim <= 0:
        raise ValueError("--max-preview-dim must be > 0")
    if not 0 <= args.cytassist_border_fraction < 0.5:
        raise ValueError("--cytassist-border-fraction must be >= 0 and < 0.5")
    if args.roi_padding < 0:
        raise ValueError("--roi-padding must be >= 0")
    if args.min_overlap_pixels < 0:
        raise ValueError("--min-overlap-pixels must be >= 0")
    if args.max_rotation_degrees < 0 or args.max_scale_change < 0 or args.max_shear < 0 or args.max_translation < 0:
        raise ValueError("Refinement bounds must be >= 0")
    if args.iterations <= 0:
        raise ValueError("--iterations must be > 0")
    if args.epsilon <= 0:
        raise ValueError("--epsilon must be > 0")
    run(args)


if __name__ == "__main__":
    main()
