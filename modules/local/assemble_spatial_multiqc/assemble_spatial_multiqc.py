#!/usr/bin/env python3

import argparse
import base64
import html
import json
import mimetypes
import re
from pathlib import Path


IMAGE_SECTIONS = [
    (
        "registration_full_slide_mqc",
        "Registration overview",
        "Full-slide microscopy with the registered Visium capture area.",
    ),
    (
        "qc_mqc",
        "Spatial QC overlays",
        "Spatial QC metrics overlaid on the microscopy image.",
    ),
    (
        "crop_areas_downsampled_microscopy_mqc",
        "Selected crop areas",
        "Crop areas used for segmentation inspection, shown on the downsampled microscopy image.",
    ),
    (
        "segmentation_crop_panels_mqc",
        "Segmentation panels",
        "Microscopy and segmentation views for the selected crop areas.",
    ),
    (
        "deconvolution_mqc",
        "Deconvolution",
        "Total estimated cell abundance and condensed cell-type composition.",
    ),
]


def safe_id(value):
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", value).strip("_")
    return safe or "sample"


def image_sort_key(path):
    name = path.name
    for index, (token, _title, _description) in enumerate(IMAGE_SECTIONS):
        if token in name:
            return index, name
    return len(IMAGE_SECTIONS), name


def image_section(path):
    name = path.name
    for token, title, description in IMAGE_SECTIONS:
        if token in name:
            return title, description
    return path.stem.replace("_", " "), ""


def image_data_uri(path):
    mime_type = mimetypes.guess_type(path.name)[0] or "image/png"
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime_type};base64,{encoded}"


def render_html(sample_name, images):
    sample_id = safe_id(sample_name)
    sample_name_yaml = json.dumps(sample_name)
    header = f"""<!--
parent_id: spatial_overview
parent_name: "Spatial overview"
parent_description: |
  Per-sample spatial QC, imaging, segmentation, and deconvolution review outputs.
id: "spatial_overview_{sample_id}"
section_name: {sample_name_yaml}
plot_type: "html"
-->
"""

    blocks = [
        '<div class="spatialomics-sample-report" '
        'style="display:flex;flex-direction:column;gap:1.25rem;">'
    ]

    for path in sorted(images, key=image_sort_key):
        title, description = image_section(path)
        escaped_title = html.escape(title)
        escaped_description = html.escape(description)
        escaped_alt = html.escape(f"{sample_name} {title}")
        data_uri = image_data_uri(path)
        description_html = (
            f'<p style="margin:0 0 .5rem 0;">{escaped_description}</p>' if description else ""
        )
        blocks.append(
            '<section style="display:flex;flex-direction:column;gap:.35rem;">'
            f'<h4 style="margin:0;">{escaped_title}</h4>'
            f"{description_html}"
            '<figure style="margin:0;">'
            f'<img src="{data_uri}" alt="{escaped_alt}" '
            'style="max-width:100%;height:auto;border:1px solid #ddd;" />'
            "</figure>"
            "</section>"
        )

    blocks.append("</div>")
    return header + "\n".join(blocks) + "\n"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--images", required=True, nargs="+", type=Path)
    args = parser.parse_args()

    images = [path for path in args.images if path.exists() and path.stat().st_size > 0]
    if not images:
        raise SystemExit("No non-empty images were provided")

    args.output.write_text(render_html(args.sample_name, images), encoding="utf-8")


if __name__ == "__main__":
    main()
