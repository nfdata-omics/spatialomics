#!/usr/bin/env python3

import argparse
import os
import shutil

import numpy as np
import tifffile


TIFF_EXTENSIONS = {".btf", ".tif", ".tiff"}


def metadata(image_path: str, factor: int) -> dict:
    return {
        "axes": "YX",
        "DownsampleFactor": factor,
        "OriginalFile": os.path.basename(image_path),
    }


def downsampled_size(size: int, factor: int) -> int:
    return 0 if size == 0 else ((size - 1) // factor) + 1


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


def squeeze_singleton_axes(
    image: np.ndarray,
    axes: str,
) -> tuple[np.ndarray, str]:
    """Remove non-spatial singleton dimensions while keeping axes in sync."""
    if len(axes) != image.ndim:
        raise ValueError(
            f"TIFF axes '{axes}' do not match image shape {image.shape}"
        )

    squeeze_axes = tuple(
        index
        for index, (axis, size) in enumerate(zip(axes, image.shape))
        if axis not in {"Y", "X"} and size == 1
    )
    if squeeze_axes:
        image = np.squeeze(image, axis=squeeze_axes)
        axes = "".join(
            axis for index, axis in enumerate(axes) if index not in squeeze_axes
        )

    return image, axes


def reduce_loaded_tiff_to_2d(image: np.ndarray, axes: str) -> np.ndarray:
    """Reduce an already loaded grayscale or color TIFF to a 2D image."""
    image, axes = squeeze_singleton_axes(image, axes)

    if axes == "YX" and image.ndim == 2:
        return image

    channel_axes = [axis for axis in ("S", "C") if axis in axes]
    if (
        len(channel_axes) == 1
        and set(axes) == {"Y", "X", channel_axes[0]}
        and image.ndim == 3
    ):
        channel_axis = axes.index(channel_axes[0])
        channel_count = image.shape[channel_axis]
        if channel_count > 10:
            raise ValueError(
                f"Unsupported TIFF channel count {channel_count} for axes '{axes}'"
            )

        # Move channels last as a view. Process rows in chunks so the channel
        # reduction does not create an image-sized float64 temporary array.
        image = np.moveaxis(image, channel_axis, -1)
        output = np.empty(image.shape[:2], dtype=np.uint8)
        rows_per_chunk = 1024
        for y_start in range(0, image.shape[0], rows_per_chunk):
            y_stop = min(y_start + rows_per_chunk, image.shape[0])
            chunk = image[y_start:y_stop]
            if np.issubdtype(chunk.dtype, np.integer):
                reduced = np.sum(chunk, axis=-1, dtype=np.uint64)
                reduced //= channel_count
            else:
                reduced = np.mean(chunk, axis=-1, dtype=np.float32)

            output[y_start:y_stop] = np.clip(reduced, 0, 255).astype(np.uint8)

        return output

    raise ValueError(
        "Unsupported TIFF layout: "
        f"shape={image.shape}, axes='{axes}'. "
        "Supported layouts are YX grayscale and YXS/SYX/YXC/CYX color images; "
        "singleton non-spatial axes are allowed."
    )


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
    axes: str,
    shape: tuple[int, ...],
) -> None:
    """Convert a TIFF page strip-by-strip to avoid loading full images into RAM."""
    with tifffile.TiffFile(image_path) as tif:
        page = tif.pages[0]

        if axes == "YX":
            height, width = shape
            output_dtype = page.dtype
        else:  # YXS, selected in main
            height, width = shape[:2]
            output_dtype = np.uint8

        print(
            "Streaming TIFF conversion: "
            f"shape={shape}, dtype={page.dtype}, axes={axes}"
        )
        output = tifffile.memmap(
            output_path,
            shape=(downsampled_size(height, factor), downsampled_size(width, factor)),
            dtype=output_dtype,
            bigtiff=True,
            compression=None,
            photometric="minisblack",
            metadata=metadata(image_path, factor),
        )

        for data, indices, _shape in page.segments(maxworkers=1, sort=True):
            if data is None:
                continue

            y_start = indices[-3]
            x_start = indices[-2]
            chunk = normalize_segment(data, axes)

            # Decoded edge tiles can include padding beyond the logical image
            # dimensions. Discard it before calculating source/output indices.
            valid_height = min(chunk.shape[0], height - y_start)
            valid_width = min(chunk.shape[1], width - x_start)
            if valid_height <= 0 or valid_width <= 0:
                continue
            chunk = chunk[:valid_height, :valid_width]

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


def downsample_stride(img: np.ndarray, factor: int) -> np.ndarray:
    return img[::factor, ::factor]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Input microscopy image (.btf, .tif, .tiff)")
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

    if factor < 1:
        parser.error("--factor must be greater than or equal to 1")

    if args.output is None:
        base = os.path.splitext(os.path.basename(image_path))[0]
        output_path = f"{base}_downsampled_{factor}x.ome.tif"
    else:
        output_path = args.output

    extension = os.path.splitext(image_path)[1].lower()
    if extension not in TIFF_EXTENSIONS:
        raise ValueError(
            f"Unsupported image format '{extension}'. "
            "Supported formats: .btf, .tif, .tiff"
        )

    # Inspect the TIFF and collect the few values needed to choose one of the
    # four conversion paths below.
    with tifffile.TiffFile(image_path) as tif:
        page_count = len(tif.pages)
        series_count = len(tif.series)
        print(
            f"TIFF structure: pages={page_count}, series={series_count}, "
            f"bigtiff={tif.is_bigtiff}"
        )

        if page_count == 1:
            source = "page"
            page = tif.pages[0]
            shape = tuple(page.shape)
            axes = page.axes
            page_is_memmappable = page.is_memmappable
            print(
                "TIFF page: "
                f"shape={shape}, dtype={page.dtype}, axes={axes}, "
                f"photometric={page.photometric.name}, "
                f"planarconfig={getattr(page.planarconfig, 'name', None)}, "
                f"tiled={page.is_tiled}, compression={page.compression.name}"
            )
        elif series_count == 1:
            source = "series"
            series = tif.series[0]
            shape = tuple(series.shape)
            axes = series.axes
            page_is_memmappable = False
            print(
                "TIFF series: "
                f"shape={shape}, dtype={series.dtype}, axes={axes}"
            )
        else:
            raise ValueError(
                "Unsupported multi-series TIFF: "
                f"found {series_count} series and {page_count} pages"
            )

    # Ignore singleton dimensions such as T=1 or Z=1 when checking the layout.
    normalized = [
        (size, axis)
        for size, axis in zip(shape, axes)
        if axis in {"Y", "X"} or size != 1
    ]
    normalized_shape = tuple(size for size, _axis in normalized)
    normalized_axes = "".join(axis for _size, axis in normalized)

    grayscale = normalized_axes == "YX" and len(normalized_shape) == 2
    channel_axes = [axis for axis in ("S", "C") if axis in normalized_axes]
    color = False
    if len(channel_axes) == 1 and len(normalized_shape) == 3:
        channel_axis = channel_axes[0]
        channel_count = normalized_shape[normalized_axes.index(channel_axis)]
        color = (
            normalized_axes.replace(channel_axis, "") == "YX"
            and 1 <= channel_count <= 10
        )

    if not grayscale and not color:
        raise ValueError(
            "Unsupported TIFF layout: "
            f"shape={shape}, axes='{axes}'. "
            "Supported layouts are YX grayscale and YXS/SYX/YXC/CYX color "
            "images; singleton non-spatial axes are allowed."
        )

    image = None
    if source == "series":
        print("Selected conversion strategy: load-series")
        print(f"Reading TIFF series: {image_path}")
        with tifffile.TiffFile(image_path) as tif:
            image = tif.series[0].asarray()
    elif axes == "YX" and factor == 1 and page_is_memmappable:
        print("Selected conversion strategy: copy")
        print(f"Input is already a memmappable 2D TIFF: {image_path}")
        print(f"Copying without loading pixel data: {output_path}")
        copy_image(image_path, output_path)
    elif axes == "YX" or (axes == "YXS" and shape[-1] <= 10):
        print("Selected conversion strategy: stream-page")
        stream_tiff_to_memmappable_2d(
            image_path,
            output_path,
            factor,
            axes,
            shape,
        )
        print(f"Writing: {output_path}")
    else:
        print("Selected conversion strategy: load-page")
        print(f"Reading physical TIFF page: {image_path}")
        with tifffile.TiffFile(image_path) as tif:
            image = tif.pages[0].asarray()

    if image is not None:
        print(f"Raw shape: {image.shape}, dtype: {image.dtype}, axes: {axes}")
        image = reduce_loaded_tiff_to_2d(image, axes)
        print(f"2D shape: {image.shape}, dtype: {image.dtype}")

        print(f"Downsampling by factor {factor}")
        image = downsample_stride(image, factor)
        print(f"Downsampled shape: {image.shape}, dtype: {image.dtype}")

        print(f"Writing: {output_path}")
        tifffile.imwrite(
            output_path,
            image,
            bigtiff=True,
            compression=None,
            photometric="minisblack",
            contiguous=True,
            metadata=metadata(image_path, factor),
        )

    print("Done.")


if __name__ == "__main__":
    main()
