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


def image_segmenter(image_path, output_path, dowsample_factor=2, tile_size=1024, overlap=50):
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

    # downscale
    img = img[::dowsample_factor, ::dowsample_factor]
    print(img.shape)

    # Use GPU only when available, otherwise fall back to CPU.
    use_gpu = torch.cuda.is_available()
    print(f"Cellpose GPU enabled: {use_gpu}")
    model = models.CellposeModel(gpu=use_gpu)

    stitched = stitch_instance_tiles_with_union(
        img,
        img.shape,
        lambda tile: tile_segmenter(tile, model),
        tile_size=tile_size,
        overlap=overlap
    )

    tifffile.imwrite(
        output_path,
        stitched,
        metadata={
            "axes": "YX",
            "DownsampleFactor": dowsample_factor,
        },
    )


def versions_yaml(process_name, list_of_libs=None):
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
            overlap=args.overlap
        )
