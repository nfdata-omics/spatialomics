"""Tile-based Cellpose segmentation and stitching utilities for spatial omics images."""

import os
import sys
import importlib
import importlib.metadata
from typing import Tuple, Iterator
import argparse
import yaml
import numpy as np
import tifffile
import torch
from cellpose import models
import czifile


class UnionFind:
    """Disjoint-set data structure used to merge overlapping instance IDs."""

    def __init__(self):
        """Initialize an empty parent mapping."""
        self.parent = {}

    def find(self, x):
        """Return the representative for ``x`` with path compression."""
        if x not in self.parent:
            self.parent[x] = x
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a, b):
        """Merge the sets that contain ``a`` and ``b``."""
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra


def generate_tiles(
    img_shape: Tuple[int, int],
    tile_size: int,
    overlap: int
) -> Iterator[Tuple[int, int, int, int]]:
    """
    Yields (y0, y1, x0, x1) tile coordinates.
    """
    # Compute sliding-window step from tile size and overlap.
    h, w = img_shape
    stride = tile_size - overlap

    for y in range(0, h, stride):
        for x in range(0, w, stride):
            # Clamp tiles at image borders so the final tiles can be smaller.
            y1 = min(y + tile_size, h)
            x1 = min(x + tile_size, w)
            yield y, y1, x, x1


def remove_edge_instances(tile_mask: np.ndarray) -> np.ndarray:
    """
    Remove all instances that touch the border of the tile.

    Parameters
    ----------
    tile_mask : np.ndarray
        2D instance-labeled mask (0 = background, >0 = instance ID)

    Returns
    -------
    np.ndarray
        Instance mask with edge-touching instances removed
    """

    if tile_mask.ndim != 2:
        raise ValueError("tile_mask must be a 2D array")

    cleaned = tile_mask.copy()
    if cleaned.shape[0] < 3 or cleaned.shape[1] < 3:
        return np.zeros_like(cleaned)

    # Collect instance IDs that touch the tile boundary
    edge_ids = set()

    # Top and bottom rows
    edge_ids.update(np.unique(cleaned[1, :]))
    edge_ids.update(np.unique(cleaned[-2, :]))

    # Left and right columns
    edge_ids.update(np.unique(cleaned[:, 1]))
    edge_ids.update(np.unique(cleaned[:, -2]))

    # Remove background if present
    edge_ids.discard(0)

    # Delete edge-touching instances
    for inst_id in edge_ids:
        cleaned[cleaned == inst_id] = 0

    return cleaned


def tile_segmenter(tile, model):
    """Run Cellpose on a tile and remove instances touching tile borders."""
    # flows/styles are unused but returned by Cellpose.
    masks, _flows, _styles = model.eval(tile)
    tile_mask = remove_edge_instances(masks)
    return tile_mask


def stitch_instance_tiles_with_union(
    img2d,
    image_shape2d,
    tile_segmenter_function,
    tile_size=1024,
    overlap=128,
    dtype=np.uint32
):
    """
    Segment a 2D image in overlapping tiles and stitch instance masks.

    Overlapping instance IDs are merged via union-find and relabeled to a
    contiguous global ID space.
    """
    # Global canvas storing a provisional stitched instance map.
    stitched = np.zeros(image_shape2d, dtype=dtype)
    uf = UnionFind()

    current_max_id = 0
    stride = tile_size - overlap

    for y0 in range(0, image_shape2d[0], stride):
        for x0 in range(0, image_shape2d[1], stride):
            y1 = min(y0 + tile_size, image_shape2d[0])
            x1 = min(x0 + tile_size, image_shape2d[1])

            # Run segmentation on the current tile.
            tile = img2d[y0:y1, x0:x1]
            tile_mask = tile_segmenter_function(tile).astype(dtype)

            # Assign provisional global IDs to tile instances
            local_ids = np.unique(tile_mask)
            local_ids = local_ids[local_ids != 0]

            # Give every local tile label a unique provisional global ID.
            local_to_global = {}
            for lid in local_ids:
                current_max_id += 1
                local_to_global[lid] = current_max_id
                uf.find(current_max_id)

            # Replace tile-local IDs with provisional global IDs.
            remapped = np.zeros_like(tile_mask)
            for lid, gid in local_to_global.items():
                remapped[tile_mask == lid] = gid

            stitched_view = stitched[y0:y1, x0:x1]

            # Detect overlaps and UNION instances
            overlap_mask = (stitched_view > 0) & (remapped > 0)
            if overlap_mask.any():
                overlapping_pairs = np.stack(
                    (stitched_view[overlap_mask], remapped[overlap_mask]),
                    axis=1
                )
                for a, b in overlapping_pairs:
                    uf.union(int(a), int(b))

            # Fill only empty pixels; keep existing labels until final relabel.
            write_mask = (stitched_view == 0) & (remapped > 0)
            stitched_view[write_mask] = remapped[write_mask]

    # Compress all equivalent provisional IDs into contiguous IDs [1..N].
    relabel_map = {}
    next_id = 1

    flat = stitched.ravel()
    for i in range(flat.size):
        if flat[i] == 0:
            continue
        root = uf.find(int(flat[i]))
        if root not in relabel_map:
            relabel_map[root] = next_id
            next_id += 1
        flat[i] = relabel_map[root]

    return stitched


def downsampled_size(size, factor):
    """Return the size produced by slicing an axis with [::factor]."""
    return 0 if size == 0 else ((size - 1) // factor) + 1


def read_visium_bounds(bounds_path):
    """Read x0, y0, x1, y1 full-resolution bounds from a TSV file."""
    with open(bounds_path, "r", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle if line.strip()]

    if len(lines) < 2:
        raise ValueError(f"Visium bounds file must contain a header and one data row: {bounds_path}")

    header = lines[0].split("\t")
    values = lines[1].split("\t")
    if len(header) != len(values):
        raise ValueError(f"Visium bounds file has inconsistent header/data columns: {bounds_path}")

    row = dict(zip(header, values))
    missing = [field for field in ("x0", "y0", "x1", "y1") if field not in row]
    if missing:
        raise ValueError(f"Visium bounds file is missing required columns: {', '.join(missing)}")

    try:
        x0, y0, x1, y1 = [int(float(row[field])) for field in ("x0", "y0", "x1", "y1")]
    except ValueError as exc:
        raise ValueError(f"Visium bounds file contains non-numeric coordinates: {bounds_path}") from exc

    if x0 >= x1 or y0 >= y1:
        raise ValueError(f"Invalid Visium bounds; expected x0 < x1 and y0 < y1, got {x0}:{y0}:{x1}:{y1}")

    return x0, y0, x1, y1


def crop_bounds_for_downsample(bounds, image_shape, padding, downsample_factor):
    """Pad, clamp, and align crop starts to the downsample grid."""
    if padding < 0:
        raise ValueError(f"--crop-padding must be >= 0, got {padding}")
    if downsample_factor <= 0:
        raise ValueError(f"--downsample-factor must be > 0, got {downsample_factor}")

    height, width = image_shape
    x0, y0, x1, y1 = bounds
    x0 = max(0, x0 - padding)
    y0 = max(0, y0 - padding)
    x1 = min(width, x1 + padding)
    y1 = min(height, y1 + padding)

    x0 = (x0 // downsample_factor) * downsample_factor
    y0 = (y0 // downsample_factor) * downsample_factor

    if x0 >= x1 or y0 >= y1:
        raise ValueError(
            "Padded Visium crop does not overlap the microscopy image: "
            f"{x0}:{y0}:{x1}:{y1}"
        )

    return x0, y0, x1, y1


def paste_crop_mask(crop_mask, full_image_shape, crop_bounds, downsample_factor, dtype=np.uint32):
    """Paste a downsampled crop mask into a full-image downsampled mask."""
    full_mask_shape = tuple(downsampled_size(size, downsample_factor) for size in full_image_shape)
    full_mask = np.zeros(full_mask_shape, dtype=dtype)

    x0, y0, _x1, _y1 = crop_bounds
    mask_y0 = y0 // downsample_factor
    mask_x0 = x0 // downsample_factor
    mask_y1 = mask_y0 + crop_mask.shape[0]
    mask_x1 = mask_x0 + crop_mask.shape[1]

    if mask_y1 > full_mask.shape[0] or mask_x1 > full_mask.shape[1]:
        raise ValueError(
            "Cropped segmentation mask does not fit in full-image mask: "
            f"crop mask shape={crop_mask.shape}, full mask shape={full_mask.shape}, "
            f"paste slices y={mask_y0}:{mask_y1}, x={mask_x0}:{mask_x1}"
        )

    full_mask[mask_y0:mask_y1, mask_x0:mask_x1] = crop_mask
    return full_mask


def image_segmenter(
    image_path,
    output_path,
    dowsample_factor=2,
    tile_size=1024,
    overlap=50,
    visium_bounds_path=None,
    crop_padding=256,
):
    """
    Segment an image with Cellpose using tiled inference and save the mask.

    Parameters
    ----------
    image_path : str
        Input path to the source CZI/BTF/TIFF image.
    output_path : str
        Output path where the stitched segmentation mask is written as TIFF.
    dowsample_factor : int, optional
        Step used for spatial downsampling before segmentation.
    tile_size : int, optional
        Tile size used for tiled segmentation.
    overlap : int, optional
        Pixel overlap between neighboring tiles.
    visium_bounds_path : str, optional
        TSV file containing x0, y0, x1, y1 full-resolution Visium bounds.
    crop_padding : int, optional
        Full-resolution pixels added on each side of the Visium crop before segmentation.
    """

    # Load image at full resolution.
    ext = os.path.splitext(image_path)[1].lower()
    if ext == ".czi":
        img = czifile.imread(image_path)
    elif ext in {".btf", ".tif", ".tiff"}:
        img = tifffile.imread(image_path)
    else:
        raise ValueError(
            f"Unsupported image format '{ext}'. Supported formats: .czi, .btf, .tif, .tiff"
        )

    # Reduce to 2D for segmentation.
    img = np.squeeze(img)
    if img.ndim == 3:
        img = np.mean(img, axis=2)
    elif img.ndim != 2:
        raise ValueError(
            f"Expected 2D/3D image after squeeze, got shape {img.shape} (ndim={img.ndim})"
        )

    full_image_shape = img.shape
    crop_bounds = None
    if visium_bounds_path:
        visium_bounds = read_visium_bounds(visium_bounds_path)
        crop_bounds = crop_bounds_for_downsample(
            visium_bounds,
            full_image_shape,
            crop_padding,
            dowsample_factor,
        )
        x0, y0, x1, y1 = crop_bounds
        print(f"Segmenting padded Visium crop in full-resolution pixels: {x0}:{y0}:{x1}:{y1}")
        img_to_segment = img[y0:y1, x0:x1]
    else:
        img_to_segment = img

    # Downscale only the selected region.
    img_to_segment = img_to_segment[::dowsample_factor, ::dowsample_factor]
    print(img_to_segment.shape)

    # Use GPU only when available, otherwise fall back to CPU.
    use_gpu = torch.cuda.is_available()
    print(f"Cellpose GPU enabled: {use_gpu}")
    model = models.CellposeModel(gpu=use_gpu)

    stitched = stitch_instance_tiles_with_union(
        img_to_segment,
        img_to_segment.shape,
        lambda tile: tile_segmenter(tile, model),
        tile_size=tile_size,
        overlap=overlap
    )

    if crop_bounds is not None:
        stitched = paste_crop_mask(
            stitched,
            full_image_shape,
            crop_bounds,
            dowsample_factor,
        )

    metadata = {
        "axes": "YX",
        "DownsampleFactor": dowsample_factor,
    }
    if crop_bounds is not None:
        metadata.update(
            {
                "VisiumCropFullRes": list(crop_bounds),
                "CropPaddingFullRes": crop_padding,
            }
        )

    tifffile.imwrite(
        output_path,
        stitched,
        metadata=metadata,
    )


def versions_yaml(process_name, list_of_libs=[]):
    """
    Generate YAML formatted string with versions of relevant libraries.

    Parameters
    ----------
    process_name : str
        Process name to use as key in the versions dictionary.
    list_of_libs : list of str, optional
        List of specific library names to include in the versions dictionary.

    Returns
    -------
    str
        YAML formatted string containing library versions and Python version.
    """

    # Build output as {process_name: {lib: version, ...}}.
    versions = {}
    versions[process_name] = {}

    versions[process_name]['python'] = f"{sys.version_info.major}" \
        f".{sys.version_info.minor}.{sys.version_info.micro}"

    for lib in list_of_libs:
        try:
            version = importlib.metadata.version(lib)
        except importlib.metadata.PackageNotFoundError:
            try:
                module = importlib.import_module(lib)
                version = getattr(module, '__version__', 'unknown')
            except (ImportError, AttributeError):
                version = None
        if version is not None:
            versions[process_name][lib] = version

    return yaml.dump(versions)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Spatial quality control")
    parser.add_argument("--input-image", type=str, help="Path to input CZI image")
    parser.add_argument("--output-mask", type=str, help="Path to output segmentation mask TIFF")
    parser.add_argument("--downsample-factor", type=int, default=1,
                        help="Downsampling factor for input image")
    parser.add_argument("--tile-size", type=int, default=1024, help="Tile size for segmentation")
    parser.add_argument("--overlap", type=int, default=50, help="Overlap in pixels between tiles")
    parser.add_argument("--visium-bounds", type=str,
                        help="TSV file containing x0, y0, x1, y1 full-resolution Visium bounds")
    parser.add_argument("--crop-padding", type=int, default=256,
                        help="Full-resolution pixels to add around the Visium bounds before segmentation")
    parser.add_argument("--versions-dict", type=str,
                        help="If set, print versions of relevant libraries in YAML format and exit")
    args = parser.parse_args()

    if args.versions_dict:

        libs = [
            "tifffile",
            "cellpose",
            "czifile",
            "torch",
            "torchvision"
        ]
        print(versions_yaml(args.versions_dict, libs))

    else:

        image_segmenter(
            args.input_image,
            args.output_mask,
            dowsample_factor=args.downsample_factor,
            tile_size=args.tile_size,
            overlap=args.overlap,
            visium_bounds_path=args.visium_bounds,
            crop_padding=args.crop_padding,
        )
