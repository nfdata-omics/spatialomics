#!/usr/bin/env python3

import argparse
import os
import shutil
import numpy as np
import tifffile
import czifile


TIFF_EXTENSIONS = {".btf", ".tif", ".tiff"}


def metadata(image_path: str, factor: int) -> dict:
    return {
        "axes": "YX",
        "DownsampleFactor": factor,
        "OriginalFile": os.path.basename(image_path),
    }


def downsampled_size(size: int, factor: int) -> int:
    return 0 if size == 0 else ((size - 1) // factor) + 1


def input_is_memmappable_2d_tiff(image_path: str) -> bool:
    """Return True if the TIFF can be used directly without pixel rewriting."""
    ext = os.path.splitext(image_path)[1].lower()
    if ext not in TIFF_EXTENSIONS:
        return False

    try:
        image = tifffile.memmap(image_path)
    except Exception:
        return False

    return image.squeeze().ndim == 2


def copy_image(image_path: str, output_path: str) -> None:
    """Copy an image without loading its pixel data into Python memory."""
    if os.path.abspath(image_path) == os.path.abspath(output_path):
        print("Input and output paths are identical; leaving file unchanged.")
        return

    shutil.copyfile(image_path, output_path)


def reduce_chunk_to_2d(chunk: np.ndarray) -> np.ndarray:
    """Reduce a TIFF strip/tile chunk to 2D without large temporary arrays."""
    if chunk.ndim == 2:
        return chunk

    if chunk.ndim == 3 and chunk.shape[-1] <= 10:
        if np.issubdtype(chunk.dtype, np.integer):
            return (
                chunk.astype(np.uint32, copy=False).sum(axis=-1) // chunk.shape[-1]
            ).astype(np.uint8)

        return np.clip(chunk.mean(axis=-1), 0, 255).astype(np.uint8)

    raise ValueError(f"Expected 2D image or YXS chunk, got {chunk.shape}")


def normalize_segment(data: np.ndarray, axes: str) -> np.ndarray:
    """Normalize tifffile segment output to either YX or YXS."""
    if data.ndim == 4 and data.shape[0] == 1:
        data = data[0]

    if axes == "YX" and data.ndim == 3 and data.shape[-1] == 1:
        data = data[..., 0]

    return data


def stream_tiff_to_memmappable_2d(
    image_path: str,
    output_path: str,
    factor: int,
) -> bool:
    """Convert a TIFF page strip-by-strip to avoid loading full images into RAM."""
    ext = os.path.splitext(image_path)[1].lower()
    if ext not in TIFF_EXTENSIONS:
        return False

    with tifffile.TiffFile(image_path) as tif:
        if len(tif.series) != 1 or len(tif.pages) != 1:
            return False

        series = tif.series[0]
        axes = series.axes

        if axes == "YX":
            height, width = series.shape
            output_dtype = series.dtype
        elif axes == "YXS" and series.shape[-1] <= 10:
            height, width = series.shape[:2]
            output_dtype = np.uint8
        else:
            return False

        print(f"Streaming TIFF conversion: shape={series.shape}, dtype={series.dtype}, axes={axes}")
        output = tifffile.memmap(
            output_path,
            shape=(downsampled_size(height, factor), downsampled_size(width, factor)),
            dtype=output_dtype,
            bigtiff=True,
            compression=None,
            photometric="minisblack",
            metadata=metadata(image_path, factor),
        )

        page = tif.pages[0]
        for data, indices, _shape in page.segments(maxworkers=1, sort=True):
            if data is None:
                continue

            y_start = indices[-3]
            x_start = indices[-2]
            chunk = normalize_segment(data, axes)

            y_indices = np.arange(y_start, y_start + chunk.shape[0])
            x_indices = np.arange(x_start, x_start + chunk.shape[1])
            y_local = np.flatnonzero(y_indices % factor == 0)
            x_local = np.flatnonzero(x_indices % factor == 0)
            if len(y_local) == 0 or len(x_local) == 0:
                continue

            reduced = reduce_chunk_to_2d(chunk[y_local][:, x_local])
            out_y = y_indices[y_local] // factor
            out_x = x_indices[x_local] // factor
            output[
                out_y[0] : out_y[-1] + 1,
                out_x[0] : out_x[-1] + 1,
            ] = reduced

        output.flush()

    return True


def read_image(image_path: str) -> np.ndarray:
    """Read supported microscopy image formats into a NumPy array."""
    ext = os.path.splitext(image_path)[1].lower()

    if ext == ".czi":
        return czifile.imread(image_path)

    if ext in TIFF_EXTENSIONS:
        return tifffile.imread(image_path)

    raise ValueError(
        f"Unsupported image format '{ext}'. Supported formats: .czi, .btf, .tif, .tiff"
    )


def reduce_to_2d(img: np.ndarray) -> np.ndarray:
    img = np.squeeze(img)

    if img.ndim == 2:
        return img

    if img.ndim == 3:
        # assume last axis is channel if small
        if img.shape[-1] <= 10:
            img = img.mean(axis=-1)
        else:
            img = img.mean(axis=0)

        if np.issubdtype(img.dtype, np.floating):
            img = np.clip(img, 0, 255).astype(np.uint8)

        return img

    raise ValueError(f"Expected 2D/3D image after squeeze, got {img.shape}")


def downsample_stride(img: np.ndarray, factor: int) -> np.ndarray:
    return img[::factor, ::factor]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input microscopy image (.czi, .btf, .tif, .tiff)")
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output TIFF path",
    )
    parser.add_argument(
        "-f",
        "--factor",
        type=int,
        default=8,
        help="Downsampling factor, e.g. 8 means 1/8 x 1/8",
    )
    parser.add_argument(
        "--compression",
        default="none",
        choices=["none"],
        help="TIFF compression. Only 'none' is supported because the output must be memmappable.",
    )
    args = parser.parse_args()

    image_path = args.input
    factor = args.factor

    if args.output is None:
        base = os.path.splitext(os.path.basename(image_path))[0]
        output_path = f"{base}_downsampled_{factor}x.ome.tif"
    else:
        output_path = args.output

    if factor == 1 and input_is_memmappable_2d_tiff(image_path):
        print(f"Input is already a memmappable 2D TIFF: {image_path}")
        print(f"Copying without loading pixel data: {output_path}")
        copy_image(image_path, output_path)
        print("Done.")
        return

    if stream_tiff_to_memmappable_2d(image_path, output_path, factor):
        print(f"Writing: {output_path}")
        print("Done.")
        return

    print(f"Reading: {image_path}")
    img = read_image(image_path)
    print(f"Raw shape: {img.shape}, dtype: {img.dtype}")

    img = reduce_to_2d(img)
    print(f"2D shape: {img.shape}, dtype: {img.dtype}")

    print(f"Downsampling by factor {factor}")
    img_down = downsample_stride(img, factor)
    print(f"Downsampled shape: {img_down.shape}, dtype: {img_down.dtype}")

    print(f"Writing: {output_path}")
    tifffile.imwrite(
        output_path,
        img_down,
        bigtiff=True,
        compression=None,
        photometric="minisblack",
        contiguous=True,
        metadata=metadata(image_path, factor),
    )

    print("Done.")


if __name__ == "__main__":
    main()
