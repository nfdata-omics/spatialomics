#!/usr/bin/env python3
"""Estimate a partial-FOV CytAssist-to-microscopy transform and write QC artifacts."""

import argparse
import csv
import importlib
import importlib.metadata
import json
import math
import sys
import warnings
from dataclasses import dataclass
from pathlib import Path

import cv2
import numpy as np
import SimpleITK as sitk
import tifffile
import yaml
from PIL import Image
from scipy import optimize
from skimage import filters, morphology
from skimage.metrics import normalized_mutual_information


warnings.filterwarnings("ignore", category=FutureWarning)


@dataclass
class Candidate:
    """A candidate transform from original CytAssist preview coordinates to microscopy preview coordinates."""

    matrix: np.ndarray
    method: str
    orientation: str
    scale: float
    translation_x: float
    translation_y: float
    coarse_score: float
    tested_candidates: int
    metrics: dict | None = None
    warped: np.ndarray | None = None


def versions_yaml(process_name, list_of_libs):
    """Generate YAML formatted string containing relevant library versions."""
    versions = {process_name: {}}
    versions[process_name]["python"] = (
        f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
    )

    for lib in list_of_libs:
        if lib == "SimpleITK":
            versions[process_name][lib] = sitk.Version_VersionString()
            continue

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
    if data.size == 0:
        raise ValueError("Cannot normalize an empty image")

    low, high = np.percentile(data, [1, 99])
    if not np.isfinite(low) or not np.isfinite(high) or high <= low:
        low = float(np.min(data))
        high = float(np.max(data))
    if high <= low:
        return np.zeros(data.shape, dtype=np.float32)

    data = (data - low) / (high - low)
    return np.clip(data, 0, 1).astype(np.float32, copy=False)


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

    # The border is usually background. Use it to choose Otsu polarity, so
    # large tissue sections are not rejected just because they occupy most of
    # the microscopy image.
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


def clean_mask(mask):
    """Remove small mask speckles and holes."""
    min_size = max(16, int(mask.size * 0.0001))
    cleaned = morphology.remove_small_objects(mask.astype(bool), min_size=min_size)
    cleaned = morphology.remove_small_holes(cleaned, area_threshold=min_size)
    return cleaned.astype(bool)


def suppress_border(mask, border_fraction):
    """Suppress mask content close to image borders, where CytAssist fiducials/frame can dominate."""
    if border_fraction <= 0:
        return mask
    height, width = mask.shape
    border_y = int(round(height * border_fraction))
    border_x = int(round(width * border_fraction))
    if border_y == 0 and border_x == 0:
        return mask
    trimmed = mask.copy()
    if border_y > 0:
        trimmed[:border_y, :] = False
        trimmed[-border_y:, :] = False
    if border_x > 0:
        trimmed[:, :border_x] = False
        trimmed[:, -border_x:] = False
    return trimmed


def mask_boundary(mask):
    """Return a thin binary boundary image from a mask."""
    eroded = morphology.erosion(mask)
    boundary = mask & ~eroded
    return boundary.astype(np.float32)


def registration_signal(mask):
    """Create a template-matching signal from tissue area and tissue boundary."""
    mask_f = mask.astype(np.float32)
    boundary_f = mask_boundary(mask)
    signal = 0.70 * mask_f + 0.30 * boundary_f
    if np.max(signal) > 0:
        signal = signal / np.max(signal)
    return signal.astype(np.float32)


def orientation_matrix(shape, orientation):
    """Return matrix mapping original moving preview coordinates into oriented-preview coordinates."""
    height, width = shape
    if orientation == "none":
        return np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    if orientation == "flip_y":
        return np.array([[1.0, 0.0, 0.0], [0.0, -1.0, height - 1.0], [0.0, 0.0, 1.0]])
    if orientation == "flip_x":
        return np.array([[-1.0, 0.0, width - 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    if orientation == "flip_xy":
        return np.array([[-1.0, 0.0, width - 1.0], [0.0, -1.0, height - 1.0], [0.0, 0.0, 1.0]])
    raise ValueError(f"Unsupported orientation: {orientation}")


def orient_image(image, orientation):
    """Apply a candidate axis flip."""
    if orientation == "none":
        return image
    if orientation == "flip_y":
        return np.flipud(image)
    if orientation == "flip_x":
        return np.fliplr(image)
    if orientation == "flip_xy":
        return np.flipud(np.fliplr(image))
    raise ValueError(f"Unsupported orientation: {orientation}")


def center_fit_transform(moving_shape, fixed_shape):
    """Return a similarity transform that scales the moving image into the fixed image center."""
    moving_h, moving_w = moving_shape
    fixed_h, fixed_w = fixed_shape
    scale = min(fixed_w / moving_w, fixed_h / moving_h)
    tx = (fixed_w - moving_w * scale) / 2.0
    ty = (fixed_h - moving_h * scale) / 2.0
    return np.array([[scale, 0.0, tx], [0.0, scale, ty]], dtype=np.float64)


def scale_grid(args, moving_mask=None, fixed_mask=None):
    """Return preview-coordinate scales to test."""
    if args.scale_min <= 0 or args.scale_max <= 0:
        raise ValueError("--scale-min and --scale-max must be > 0")
    if args.scale_min > args.scale_max:
        raise ValueError("--scale-min must be <= --scale-max")
    if args.scale_steps < 1:
        raise ValueError("--scale-steps must be >= 1")
    if args.scale_steps == 1:
        scales = [args.scale_min]
    else:
        scales = list(np.geomspace(args.scale_min, args.scale_max, args.scale_steps))

    if args.auto_scale_from_masks and moving_mask is not None and fixed_mask is not None:
        moving_area = float(np.sum(moving_mask))
        fixed_area = float(np.sum(fixed_mask))
        if moving_area > 0 and fixed_area > 0:
            area_scale = math.sqrt(fixed_area / moving_area)
            if area_scale < args.scale_min or area_scale > args.scale_max:
                adaptive_min = max(0.05, area_scale / args.auto_scale_factor)
                adaptive_max = area_scale * args.auto_scale_factor
                if adaptive_min <= adaptive_max:
                    scales.extend(np.geomspace(adaptive_min, adaptive_max, args.scale_steps))

    return np.array(sorted(set(round(float(scale), 8) for scale in scales)), dtype=np.float64)


def resize_signal(signal, scale):
    """Resize a registration signal for a candidate preview scale."""
    height, width = signal.shape
    new_w = max(1, int(round(width * scale)))
    new_h = max(1, int(round(height * scale)))
    return cv2.resize(signal, (new_w, new_h), interpolation=cv2.INTER_AREA)


def template_match_candidates(moving_signal, fixed_signal, moving_shape, moving_mask, fixed_mask, args):
    """Search translation, scale, and axis-flip candidates with mask template matching."""
    candidates = []
    tested = 0
    fixed_h, fixed_w = fixed_signal.shape
    max_template_dim = max(fixed_h, fixed_w) * args.max_template_fixed_ratio

    for orientation in args.orientations:
        oriented_signal = orient_image(moving_signal, orientation)
        orientation_3x3 = orientation_matrix(moving_shape, orientation)

        for scale in scale_grid(args, moving_mask=moving_mask, fixed_mask=fixed_mask):
            template = resize_signal(oriented_signal, scale)
            templ_h, templ_w = template.shape
            if templ_h < args.min_template_dim or templ_w < args.min_template_dim:
                continue
            if max(templ_h, templ_w) > max_template_dim:
                continue
            if float(template.sum()) <= 0:
                continue

            pad_y = templ_h
            pad_x = templ_w
            padded_fixed = np.pad(
                fixed_signal,
                ((pad_y, pad_y), (pad_x, pad_x)),
                mode="constant",
                constant_values=0,
            ).astype(np.float32, copy=False)
            response = cv2.matchTemplate(padded_fixed, template.astype(np.float32), cv2.TM_CCORR_NORMED)
            _, max_value, _, max_location = cv2.minMaxLoc(response)
            tx = float(max_location[0] - pad_x)
            ty = float(max_location[1] - pad_y)
            transform = np.array(
                [[scale, 0.0, tx], [0.0, scale, ty], [0.0, 0.0, 1.0]],
                dtype=np.float64,
            )
            original_to_fixed = transform @ orientation_3x3
            candidates.append(
                Candidate(
                    matrix=original_to_fixed[:2],
                    method="partial_fov_mask_search",
                    orientation=orientation,
                    scale=float(scale),
                    translation_x=tx,
                    translation_y=ty,
                    coarse_score=float(max_value),
                    tested_candidates=0,
                )
            )
            tested += 1

    for candidate in candidates:
        candidate.tested_candidates = tested
    return candidates


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


def transformed_corners(matrix, moving_shape):
    """Transform moving-image corners into fixed-image coordinates."""
    moving_h, moving_w = moving_shape
    corners = np.array(
        [[0.0, 0.0], [moving_w, 0.0], [moving_w, moving_h], [0.0, moving_h]],
        dtype=np.float64,
    )
    homogeneous = np.c_[corners, np.ones(corners.shape[0])]
    transformed = (np.vstack([matrix, [0.0, 0.0, 1.0]]) @ homogeneous.T).T
    return transformed[:, :2] / transformed[:, 2:3]


def footprint_inside_fraction(matrix, moving_shape, fixed_shape):
    """Approximate how much of the moving footprint bbox lies inside fixed image bounds."""
    fixed_h, fixed_w = fixed_shape
    corners = transformed_corners(matrix, moving_shape)
    min_xy = corners.min(axis=0)
    max_xy = corners.max(axis=0)
    width, height = np.maximum(max_xy - min_xy, 0)
    bbox_area = float(width * height)
    if bbox_area <= 0:
        return 0.0
    inter_min = np.maximum(min_xy, [0.0, 0.0])
    inter_max = np.minimum(max_xy, [float(fixed_w), float(fixed_h)])
    inter_width, inter_height = np.maximum(inter_max - inter_min, 0)
    return float((inter_width * inter_height) / bbox_area)


def compute_metrics(fixed, moving, matrix, feature_stats, method, sample, fixed_mask=None, moving_mask=None):
    """Compute registration QC metrics on preview images."""
    warped = warp_image(moving, matrix, fixed.shape)
    if fixed_mask is None:
        fixed_mask = foreground_mask(fixed)
    if moving_mask is None:
        moving_mask = foreground_mask(moving)
    warped_mask = warp_image(
        moving_mask.astype(np.float32),
        matrix,
        fixed.shape,
        interpolation=cv2.INTER_NEAREST,
    ) > 0.5

    union = fixed_mask | warped_mask
    intersection = fixed_mask & warped_mask
    mask_iou = float(intersection.sum() / union.sum()) if union.any() else 0.0
    moving_mask_overlap_fraction = (
        float(intersection.sum() / warped_mask.sum()) if warped_mask.any() else 0.0
    )
    fixed_mask_overlap_fraction = float(intersection.sum() / fixed_mask.sum()) if fixed_mask.any() else 0.0

    overlap = intersection & (warped > 0)
    if overlap.sum() > 10:
        fixed_values = fixed[overlap]
        warped_values = warped[overlap]
        fixed_std = float(np.std(fixed_values))
        warped_std = float(np.std(warped_values))
        if fixed_std > 0 and warped_std > 0:
            ncc = float(np.corrcoef(fixed_values, warped_values)[0, 1])
        else:
            ncc = 0.0
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
    rotation_degrees = math.degrees(math.atan2(c, a))
    determinant = a * d - b * c

    metrics = {
        "sample": sample,
        "status": "ok",
        "method": method,
        "moving_keypoints": feature_stats.get("moving_keypoints", 0),
        "fixed_keypoints": feature_stats.get("fixed_keypoints", 0),
        "matches": feature_stats.get("matches", 0),
        "inliers": feature_stats.get("inliers", 0),
        "inlier_ratio": (
            feature_stats.get("inliers", 0) / feature_stats.get("matches", 1)
            if feature_stats.get("matches", 0)
            else 0.0
        ),
        "mask_iou": mask_iou,
        "moving_mask_overlap_fraction": moving_mask_overlap_fraction,
        "fixed_mask_overlap_fraction": fixed_mask_overlap_fraction,
        "normalized_mutual_information": nmi,
        "normalized_cross_correlation": ncc,
        "preview_scale_x": scale_x,
        "preview_scale_y": scale_y,
        "preview_rotation_degrees": rotation_degrees,
        "preview_translation_x": float(tx),
        "preview_translation_y": float(ty),
        "preview_determinant": float(determinant),
        "footprint_inside_fraction": footprint_inside_fraction(matrix, moving.shape, fixed.shape),
    }
    return metrics, warped


def candidate_score(metrics, coarse_score):
    """Score a candidate for partial-FOV registration.

    CytAssist usually covers only part of the microscopy image, so explaining all
    fixed-image tissue is less important than retaining CytAssist tissue and
    matching its tissue shape. Over-weighting fixed overlap biases toward
    oversized transforms that cover more microscopy tissue but have the wrong
    geometry.
    """
    ncc = max(float(metrics.get("normalized_cross_correlation", 0.0)), 0.0)
    mask_iou = float(metrics.get("mask_iou", 0.0))
    moving_overlap = float(metrics.get("moving_mask_overlap_fraction", 0.0))
    fixed_overlap = float(metrics.get("fixed_mask_overlap_fraction", 0.0))
    return (
        (2.0 * coarse_score)
        + (1.5 * moving_overlap)
        + (1.0 * mask_iou)
        + (0.25 * fixed_overlap)
        + (0.1 * ncc)
    )


def choose_best_candidate(candidates, fixed, moving, fixed_mask, moving_mask, sample):
    """Evaluate candidates with QC metrics and return the best one."""
    best = None
    for candidate in candidates:
        feature_stats = {
            "moving_keypoints": 0,
            "fixed_keypoints": 0,
            "matches": 0,
            "inliers": 0,
            "coarse_score": candidate.coarse_score,
        }
        metrics, warped = compute_metrics(
            fixed=fixed,
            moving=moving,
            matrix=candidate.matrix,
            feature_stats=feature_stats,
            method=candidate.method,
            sample=sample,
            fixed_mask=fixed_mask,
            moving_mask=moving_mask,
        )
        metrics["candidate_orientation"] = candidate.orientation
        metrics["candidate_scale"] = candidate.scale
        metrics["candidate_translation_x"] = candidate.translation_x
        metrics["candidate_translation_y"] = candidate.translation_y
        metrics["coarse_template_score"] = candidate.coarse_score
        metrics["tested_candidates"] = candidate.tested_candidates
        metrics["candidate_score"] = candidate_score(metrics, candidate.coarse_score)
        candidate.metrics = metrics
        candidate.warped = warped
        if best is None or metrics["candidate_score"] > best.metrics["candidate_score"]:
            best = candidate
    return best


def matrix_3x3(matrix):
    """Return a 3x3 homogeneous matrix from a 2x3 affine matrix."""
    return np.vstack([matrix, [0.0, 0.0, 1.0]])


def transformed_center(matrix, moving_shape):
    """Transform the moving-image center into fixed-image coordinates."""
    moving_h, moving_w = moving_shape
    center = np.array([moving_w / 2.0, moving_h / 2.0, 1.0], dtype=np.float64)
    transformed = matrix_3x3(matrix) @ center
    return transformed[:2] / transformed[2]


def linear_scale(matrix):
    """Return the mean linear scale of a 2x3 affine matrix."""
    a, b = matrix[0, 0], matrix[0, 1]
    c, d = matrix[1, 0], matrix[1, 1]
    return (math.sqrt(a * a + c * c) + math.sqrt(b * b + d * d)) / 2.0


def normalized_product_score(fixed_signal, moving_signal, matrix):
    """Return normalized image product after warping moving signal into fixed space."""
    warped = warp_image(moving_signal, matrix, fixed_signal.shape)
    denominator = math.sqrt(float(np.sum(fixed_signal * fixed_signal)) * float(np.sum(warped * warped)))
    if denominator <= 0:
        return 0.0
    return float(np.sum(fixed_signal * warped) / denominator)


def blur_signal(signal, filter_size):
    """Blur a registration signal for smoother local optimization."""
    if filter_size <= 1:
        return signal.astype(np.float32, copy=False)
    return cv2.GaussianBlur(signal.astype(np.float32, copy=False), (filter_size, filter_size), 0)


def refine_affine_local(coarse_matrix, fixed_signal, moving_signal, args):
    """Refine a coarse moving-to-fixed transform with bounded local affine optimization."""
    fixed_refine = blur_signal(fixed_signal, args.refine_gaussian_filter_size)
    moving_refine = blur_signal(moving_signal, args.refine_gaussian_filter_size)
    coarse_3x3 = matrix_3x3(coarse_matrix)
    center = transformed_center(coarse_matrix, moving_signal.shape)
    center_to_origin = np.array([[1.0, 0.0, -center[0]], [0.0, 1.0, -center[1]], [0.0, 0.0, 1.0]])
    origin_to_center = np.array([[1.0, 0.0, center[0]], [0.0, 1.0, center[1]], [0.0, 0.0, 1.0]])
    max_rotation = math.radians(args.refine_max_rotation_degrees)
    min_log_scale = math.log(max(0.05, 1.0 - args.refine_max_scale_change))
    max_log_scale = math.log(1.0 + args.refine_max_scale_change)
    bounds = [
        (-max_rotation, max_rotation),
        (min_log_scale, max_log_scale),
        (min_log_scale, max_log_scale),
        (-args.refine_max_shear, args.refine_max_shear),
        (-args.refine_max_shear, args.refine_max_shear),
        (-args.refine_max_translation, args.refine_max_translation),
        (-args.refine_max_translation, args.refine_max_translation),
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
        return origin_to_center @ affine_delta @ rotation @ center_to_origin @ coarse_3x3

    def objective(params):
        matrix = delta_matrix(params)[:2]
        return -normalized_product_score(fixed_refine, moving_refine, matrix)

    try:
        result = optimize.minimize(
            objective,
            np.zeros(7, dtype=np.float64),
            method="Powell",
            bounds=bounds,
            options={
                "maxiter": args.refine_iterations,
                "xtol": args.refine_epsilon,
                "ftol": args.refine_epsilon,
                "disp": False,
            },
        )
    except Exception as error:
        return None, {"refinement_status": "failed", "refinement_error": str(error), "refinement_objective_score": 0.0}

    refined = delta_matrix(result.x)[:2]
    coarse_det = float(np.linalg.det(coarse_matrix[:, :2]))
    refined_det = float(np.linalg.det(refined[:, :2]))
    if coarse_det == 0 or refined_det == 0 or np.sign(coarse_det) != np.sign(refined_det):
        return None, {
            "refinement_status": "rejected_determinant_sign",
            "refinement_objective_score": float(-result.fun),
            "refined_preview_determinant": refined_det,
        }

    coarse_scale = linear_scale(coarse_matrix)
    refined_scale = linear_scale(refined)
    if coarse_scale <= 0:
        return None, {"refinement_status": "rejected_invalid_coarse_scale", "refinement_objective_score": float(-result.fun)}
    scale_change = abs(refined_scale / coarse_scale - 1.0)
    if scale_change > args.refine_max_scale_change:
        return None, {
            "refinement_status": "rejected_scale_change",
            "refinement_objective_score": float(-result.fun),
            "refinement_scale_change": float(scale_change),
        }

    center_shift = float(np.linalg.norm(transformed_center(refined, moving_signal.shape) - transformed_center(coarse_matrix, moving_signal.shape)))
    max_center_shift = args.refine_max_center_shift_fraction * max(fixed_signal.shape)
    if center_shift > max_center_shift:
        return None, {
            "refinement_status": "rejected_center_shift",
            "refinement_objective_score": float(-result.fun),
            "refinement_center_shift": center_shift,
        }

    return refined, {
        "refinement_status": "accepted",
        "refinement_objective_score": float(-result.fun),
        "refinement_optimizer_success": bool(result.success),
        "refinement_optimizer_iterations": int(getattr(result, "nit", 0)),
        "refinement_scale_change": float(scale_change),
        "refinement_center_shift": center_shift,
    }


def maybe_refine_candidate(best, fixed, moving, fixed_mask, moving_mask, fixed_signal, moving_signal, sample, args):
    """Optionally refine the selected candidate and keep it only if QC remains acceptable."""
    if best is None or not args.refine_affine:
        return best

    refined_matrix, refine_stats = refine_affine_local(best.matrix, fixed_signal, moving_signal, args)
    if refined_matrix is None:
        best.metrics.update(refine_stats)
        return best

    feature_stats = {
        "moving_keypoints": 0,
        "fixed_keypoints": 0,
        "matches": 0,
        "inliers": 0,
        "coarse_score": best.coarse_score,
    }
    refined_metrics, refined_warped = compute_metrics(
        fixed=fixed,
        moving=moving,
        matrix=refined_matrix,
        feature_stats=feature_stats,
        method="partial_fov_mask_search_affine_refined",
        sample=sample,
        fixed_mask=fixed_mask,
        moving_mask=moving_mask,
    )
    refined_metrics["candidate_orientation"] = best.orientation
    refined_metrics["candidate_scale"] = best.scale
    refined_metrics["candidate_translation_x"] = best.translation_x
    refined_metrics["candidate_translation_y"] = best.translation_y
    refined_metrics["coarse_template_score"] = best.coarse_score
    refined_metrics["tested_candidates"] = best.tested_candidates
    refined_metrics.update(refine_stats)
    refined_metrics["coarse_candidate_score"] = best.metrics["candidate_score"]
    refined_metrics["candidate_score"] = candidate_score(refined_metrics, best.coarse_score)

    score_drop = best.metrics["candidate_score"] - refined_metrics["candidate_score"]
    if score_drop > args.refine_max_score_drop:
        best.metrics.update(
            {
                "refinement_status": "rejected_score_drop",
                "refinement_objective_score": refine_stats.get("refinement_objective_score", 0.0),
                "refinement_score_drop": float(score_drop),
            }
        )
        return best

    return Candidate(
        matrix=refined_matrix,
        method="partial_fov_mask_search_affine_refined",
        orientation=best.orientation,
        scale=best.scale,
        translation_x=best.translation_x,
        translation_y=best.translation_y,
        coarse_score=best.coarse_score,
        tested_candidates=best.tested_candidates,
        metrics=refined_metrics,
        warped=refined_warped,
    )


def mark_candidate_status(metrics, args):
    """Set an explicit QC status for the selected candidate."""
    if metrics["method"] == "center_fit_fallback":
        metrics["status"] = "fallback"
        return metrics
    if metrics["mask_iou"] < args.min_mask_iou:
        metrics["status"] = "rejected_low_mask_iou"
        return metrics
    if metrics["coarse_template_score"] < args.min_coarse_score:
        metrics["status"] = "rejected_low_template_score"
        return metrics
    if metrics["moving_mask_overlap_fraction"] < args.min_moving_mask_overlap:
        metrics["status"] = "rejected_low_moving_mask_overlap"
        return metrics
    metrics["status"] = "ok"
    return metrics


def preview_to_full_matrix(preview_matrix, cytassist_stride, microscopy_stride):
    """Convert preview-coordinate matrix to full-resolution coordinate matrix."""
    preview_3x3 = np.vstack([preview_matrix, [0.0, 0.0, 1.0]])
    fixed_scale = np.diag([microscopy_stride, microscopy_stride, 1.0])
    moving_scale_inv = np.diag([1.0 / cytassist_stride, 1.0 / cytassist_stride, 1.0])
    return fixed_scale @ preview_3x3 @ moving_scale_inv


def write_overlay(output_overlay, fixed, warped):
    """Write a green/magenta registration overlay PNG."""
    fixed_u8 = to_uint8(fixed)
    warped_u8 = to_uint8(warped)
    rgb = np.zeros((*fixed.shape, 3), dtype=np.uint8)
    rgb[..., 1] = fixed_u8
    rgb[..., 0] = warped_u8
    rgb[..., 2] = warped_u8
    Image.fromarray(rgb).save(output_overlay)


def write_metrics(output_metrics, metrics):
    """Write a single-row metrics CSV."""
    with output_metrics.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(metrics.keys()))
        writer.writeheader()
        writer.writerow(metrics)


def write_transform(
    output_transform,
    sample,
    method,
    preview_matrix,
    full_matrix,
    cytassist_shape,
    microscopy_shape,
    cytassist_stride,
    microscopy_stride,
    metrics,
):
    """Write the transform and registration context as JSON."""
    payload = {
        "sample": sample,
        "transform_direction": "cytassist_to_microscopy",
        "coordinate_units": "pixels",
        "method": method,
        "preview_transform_3x3": np.vstack([preview_matrix, [0.0, 0.0, 1.0]]).tolist(),
        "full_resolution_transform_3x3": full_matrix.tolist(),
        "cytassist_shape_yx": list(cytassist_shape),
        "microscopy_shape_yx": list(microscopy_shape),
        "cytassist_preview_stride": cytassist_stride,
        "microscopy_preview_stride": microscopy_stride,
        "metrics": metrics,
    }
    output_transform.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def run_registration(args):
    """Run registration and write outputs."""
    cytassist_preview_raw, cytassist_color, cytassist_shape, cytassist_stride = downsample_to_preview(
        read_image_array(args.cytassist_image),
        args.max_preview_dim,
    )
    microscopy_preview_raw, microscopy_color, microscopy_shape, microscopy_stride = downsample_to_preview(
        read_image_array(args.microscopy_tif),
        args.max_preview_dim,
    )

    moving = normalize_float(cytassist_preview_raw)
    fixed = normalize_float(microscopy_preview_raw)

    moving_mask = suppress_border(
        tissue_mask(moving, cytassist_color, exclude_green=args.exclude_green_fiducials),
        args.cytassist_border_fraction,
    )
    if not moving_mask.any():
        moving_mask = tissue_mask(moving, cytassist_color, exclude_green=False)
    fixed_mask = tissue_mask(fixed, microscopy_color, exclude_green=False)
    moving_signal = registration_signal(moving_mask)
    fixed_signal = registration_signal(fixed_mask)

    candidates = template_match_candidates(
        moving_signal=moving_signal,
        fixed_signal=fixed_signal,
        moving_shape=moving.shape,
        moving_mask=moving_mask,
        fixed_mask=fixed_mask,
        args=args,
    )
    best = choose_best_candidate(
        candidates=candidates,
        fixed=fixed,
        moving=moving,
        fixed_mask=fixed_mask,
        moving_mask=moving_mask,
        sample=args.sample_name,
    )
    best = maybe_refine_candidate(
        best=best,
        fixed=fixed,
        moving=moving,
        fixed_mask=fixed_mask,
        moving_mask=moving_mask,
        fixed_signal=fixed_signal,
        moving_signal=moving_signal,
        sample=args.sample_name,
        args=args,
    )

    if best is None:
        matrix = center_fit_transform(moving.shape, fixed.shape)
        method = "center_fit_fallback"
        metrics, warped = compute_metrics(
            fixed=fixed,
            moving=moving,
            matrix=matrix,
            feature_stats={},
            method=method,
            sample=args.sample_name,
            fixed_mask=fixed_mask,
            moving_mask=moving_mask,
        )
        metrics["candidate_orientation"] = "none"
        metrics["candidate_scale"] = float("nan")
        metrics["candidate_translation_x"] = float(matrix[0, 2])
        metrics["candidate_translation_y"] = float(matrix[1, 2])
        metrics["coarse_template_score"] = 0.0
        metrics["tested_candidates"] = 0
        metrics["candidate_score"] = 0.0
    else:
        matrix = best.matrix
        method = best.method
        metrics = best.metrics
        warped = best.warped

    metrics = mark_candidate_status(metrics, args)
    full_matrix = preview_to_full_matrix(
        matrix,
        cytassist_stride=cytassist_stride,
        microscopy_stride=microscopy_stride,
    )

    write_overlay(args.output_overlay, fixed, warped)
    write_metrics(args.output_metrics, metrics)
    write_transform(
        output_transform=args.output_transform,
        sample=args.sample_name,
        method=method,
        preview_matrix=matrix,
        full_matrix=full_matrix,
        cytassist_shape=cytassist_shape,
        microscopy_shape=microscopy_shape,
        cytassist_stride=cytassist_stride,
        microscopy_stride=microscopy_stride,
        metrics=metrics,
    )


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Estimate a partial-FOV CytAssist-to-microscopy transform and write QC artifacts"
    )
    parser.add_argument("--cytassist-image", type=Path, help="Input CytAssist image")
    parser.add_argument("--microscopy-tif", type=Path, help="Input memmappable microscopy TIFF")
    parser.add_argument("--sample-name", help="Sample name")
    parser.add_argument("--output-transform", type=Path, help="Output transform JSON")
    parser.add_argument("--output-metrics", type=Path, help="Output registration metrics CSV")
    parser.add_argument("--output-overlay", type=Path, help="Output QC overlay PNG")
    parser.add_argument(
        "--max-preview-dim",
        type=int,
        default=2048,
        help="Maximum preview dimension used for registration",
    )
    parser.add_argument(
        "--scale-min",
        type=float,
        default=0.35,
        help="Minimum preview-coordinate scale to test",
    )
    parser.add_argument(
        "--scale-max",
        type=float,
        default=1.6,
        help="Maximum preview-coordinate scale to test",
    )
    parser.add_argument(
        "--scale-steps",
        type=int,
        default=24,
        help="Number of log-spaced preview scales to test",
    )
    parser.add_argument(
        "--auto-scale-from-masks",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Add an adaptive scale band estimated from fixed/moving tissue-mask area ratio",
    )
    parser.add_argument(
        "--auto-scale-factor",
        type=float,
        default=1.6,
        help="Multiplicative half-width around the mask-area-derived adaptive scale",
    )
    parser.add_argument(
        "--orientations",
        nargs="+",
        default=["none", "flip_y", "flip_x", "flip_xy"],
        choices=["none", "flip_y", "flip_x", "flip_xy"],
        help="Axis-flip candidate orientations to test",
    )
    parser.add_argument(
        "--cytassist-border-fraction",
        type=float,
        default=0.03,
        help="Fraction of the CytAssist preview border to suppress before mask matching",
    )
    parser.add_argument(
        "--exclude-green-fiducials",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Exclude green CytAssist fiducial/frame signal from RGB tissue masking",
    )
    parser.add_argument(
        "--max-template-fixed-ratio",
        type=float,
        default=2.5,
        help="Skip candidate templates whose largest dimension exceeds this multiple of the fixed preview",
    )
    parser.add_argument(
        "--min-template-dim",
        type=int,
        default=32,
        help="Skip candidate templates smaller than this many pixels in either dimension",
    )
    parser.add_argument(
        "--min-mask-iou",
        type=float,
        default=0.03,
        help="Minimum selected mask IoU required for status=ok",
    )
    parser.add_argument(
        "--min-coarse-score",
        type=float,
        default=0.05,
        help="Minimum selected template score required for status=ok",
    )
    parser.add_argument(
        "--min-moving-mask-overlap",
        type=float,
        default=0.20,
        help="Minimum fraction of transformed CytAssist tissue mask overlapping microscopy tissue required for status=ok",
    )
    parser.add_argument(
        "--refine-affine",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Refine the selected coarse transform with bounded local affine optimization on tissue registration signals",
    )
    parser.add_argument(
        "--refine-iterations",
        type=int,
        default=300,
        help="Maximum optimizer iterations for affine refinement",
    )
    parser.add_argument(
        "--refine-epsilon",
        type=float,
        default=1e-6,
        help="Optimizer convergence epsilon for affine refinement",
    )
    parser.add_argument(
        "--refine-gaussian-filter-size",
        type=int,
        default=5,
        help="Odd Gaussian filter size used to smooth registration signals during affine refinement",
    )
    parser.add_argument(
        "--refine-max-rotation-degrees",
        type=float,
        default=3.0,
        help="Maximum absolute local rotation allowed during affine refinement",
    )
    parser.add_argument(
        "--refine-max-shear",
        type=float,
        default=0.05,
        help="Maximum absolute local shear term allowed during affine refinement",
    )
    parser.add_argument(
        "--refine-max-translation",
        type=float,
        default=80.0,
        help="Maximum local translation in preview pixels allowed during affine refinement",
    )
    parser.add_argument(
        "--refine-max-scale-change",
        type=float,
        default=0.04,
        help="Maximum fractional scale change allowed during local affine refinement",
    )
    parser.add_argument(
        "--refine-max-center-shift-fraction",
        type=float,
        default=0.15,
        help="Maximum refined center shift as a fraction of the largest fixed preview dimension",
    )
    parser.add_argument(
        "--refine-max-score-drop",
        type=float,
        default=0.10,
        help="Maximum allowed candidate score drop after affine refinement",
    )
    parser.add_argument(
        "--min-feature-matches",
        type=int,
        default=12,
        help="Deprecated compatibility option; feature matching is not used by this implementation",
    )
    parser.add_argument(
        "--versions-dict",
        help="If set, print versions of relevant libraries in YAML format and exit",
    )
    return parser.parse_args()


def main():
    """Run the command-line interface."""
    args = parse_args()
    if args.versions_dict:
        print(
            versions_yaml(
                args.versions_dict,
                ["opencv-python", "numpy", "pandas", "Pillow", "scikit-image", "scipy", "SimpleITK", "tifffile"],
            )
        )
        return

    missing = [
        name
        for name in (
            "cytassist_image",
            "microscopy_tif",
            "sample_name",
            "output_transform",
            "output_metrics",
            "output_overlay",
        )
        if getattr(args, name) is None
    ]
    if missing:
        raise ValueError(f"Missing required arguments: {', '.join(missing)}")
    if args.max_preview_dim <= 0:
        raise ValueError("--max-preview-dim must be > 0")
    if args.auto_scale_factor <= 1:
        raise ValueError("--auto-scale-factor must be > 1")
    if not 0 <= args.cytassist_border_fraction < 0.5:
        raise ValueError("--cytassist-border-fraction must be >= 0 and < 0.5")
    if args.max_template_fixed_ratio <= 0:
        raise ValueError("--max-template-fixed-ratio must be > 0")
    if args.min_template_dim <= 0:
        raise ValueError("--min-template-dim must be > 0")
    if args.refine_iterations <= 0:
        raise ValueError("--refine-iterations must be > 0")
    if args.refine_epsilon <= 0:
        raise ValueError("--refine-epsilon must be > 0")
    if args.refine_gaussian_filter_size <= 0 or args.refine_gaussian_filter_size % 2 == 0:
        raise ValueError("--refine-gaussian-filter-size must be a positive odd integer")
    if args.refine_max_rotation_degrees < 0:
        raise ValueError("--refine-max-rotation-degrees must be >= 0")
    if args.refine_max_shear < 0:
        raise ValueError("--refine-max-shear must be >= 0")
    if args.refine_max_translation < 0:
        raise ValueError("--refine-max-translation must be >= 0")
    if args.refine_max_scale_change < 0:
        raise ValueError("--refine-max-scale-change must be >= 0")
    if args.refine_max_center_shift_fraction < 0:
        raise ValueError("--refine-max-center-shift-fraction must be >= 0")
    if args.refine_max_score_drop < 0:
        raise ValueError("--refine-max-score-drop must be >= 0")
    if args.min_feature_matches <= 0:
        raise ValueError("--min-feature-matches must be > 0")

    run_registration(args)


if __name__ == "__main__":
    main()
