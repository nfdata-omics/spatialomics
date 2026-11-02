# ROI annotation exchange and import

## Objective and scope

Support a two-stage manual ROI annotation workflow around the Visium capture
area:

1. export a downsampled H&E crop and its coordinate contract;
2. annotate that image with an external application and export GeoJSON;
3. import the ROI geometries into the final SpatialData;
4. include a compact ROI overlay in the per-sample MultiQC report.

The interchange format is intentionally independent of the annotation tool.
QuPath is the initial application used for testing, but the public input is
named `roi_annotations` and accepts GeoJSON rather than a QuPath-specific
format.

Relevant implementation files are:

- [`assets/schema_input.json`](../../assets/schema_input.json);
- [`workflows/spatialomics.nf`](../../workflows/spatialomics.nf);
- [`modules/local/extract_visium_region/`](../../modules/local/extract_visium_region/);
- [`modules/local/import_roi_annotations/`](../../modules/local/import_roi_annotations/);
- [`modules/local/assemble_spatial_multiqc/`](../../modules/local/assemble_spatial_multiqc/).

## Agreed behaviour

### Export and manual annotation

`EXTRACT_VISIUM_REGION` runs for every sample with a microscopy image,
independently of segmentation. It produces:

```text
<sample>_visium_capture_area.ome.tif
<sample>_visium_bounds.tsv
```

The image contains the Visium capture area after the configured full-resolution
padding and crop downsampling. `save_capture_area_image` controls publication of
both files under `annotation_inputs/`; it does not control process execution.

The default crop padding is zero full-resolution pixels. The default crop
downsampling factor is four, meaning that both image dimensions are reduced to
approximately one quarter. ROI coordinates exported by the annotation tool are
expected to be pixel coordinates local to this downsampled image.

Manual annotation creates a natural two-run workflow. The first run publishes
the annotation inputs. A later run supplies the GeoJSON and uses `-resume` to
reuse completed upstream tasks.

### Optional GeoJSON input

The samplesheet accepts an optional `roi_annotations` column containing an
existing `.geojson` file. Presence of this file enables ROI import for the
corresponding sample; there is no separate enable or skip parameter.

For samples represented by multiple samplesheet rows, all non-empty
`roi_annotations` values must resolve to the same file. Conflicting values are
an input error.

The annotation path is routed through a dedicated channel keyed by sample ID.
It is deliberately excluded from the shared `meta` map. Adding the path on a
subsequent run therefore does not change the inputs of upstream analysis tasks;
only the import and downstream reporting tasks should need to run when their
other inputs and code remain unchanged.

## Coordinate TSV contract

The bounds TSV contains one header row and exactly one data row. Together with
the matching crop, it is the complete coordinate contract required by the ROI
importer.

| Column                          | Meaning                                                                 |
| ------------------------------- | ----------------------------------------------------------------------- |
| `x0`, `y0`, `x1`, `y1`          | Unpadded Visium bounds in full-resolution microscopy pixels             |
| `crop_x0`, `crop_y0`            | Crop origin in full-resolution microscopy pixels                        |
| `crop_x1`, `crop_y1`            | Exclusive crop end in full-resolution microscopy pixels                 |
| `source_width`, `source_height` | Full-resolution microscopy dimensions                                   |
| `padding`                       | Requested full-resolution padding                                       |
| `downsample_factor`             | Requested capture-area crop downsampling factor                         |
| `output_width`, `output_height` | Actual dimensions of the exported crop                                  |
| `scale_x`, `scale_y`            | Actual crop-local pixel to full-resolution pixel scale                  |
| `zarr_downsample_factor`        | Full-resolution microscopy pixel to SpatialData sample-coordinate scale |

For a GeoJSON coordinate `(x_roi, y_roi)`, the importer calculates:

```text
x_sample = (crop_x0 + x_roi * scale_x) / zarr_downsample_factor
y_sample = (crop_y0 + y_roi * scale_y) / zarr_downsample_factor
```

`scale_x` and `scale_y` are stored independently because rounding the output
dimensions can make the realised scale differ slightly from the requested crop
downsampling factor.

TSV files generated before `zarr_downsample_factor` was added are not complete
contracts and are rejected by the importer.

## GeoJSON and SpatialData contract

`IMPORT_ROI_ANNOTATIONS` consumes the keyed combination of:

- the most advanced SpatialData available for the sample;
- `roi_annotations` GeoJSON;
- the matching bounds TSV;
- the matching capture-area image.

The GeoJSON must be a `FeatureCollection` containing non-empty `Polygon` or
`MultiPolygon` geometries. Invalid polygon geometries are repaired with a
zero-width buffer when possible. Other geometry types are rejected rather than
being silently converted. Coordinates outside the crop dimensions currently
produce a warning because annotations may legitimately touch or slightly cross
the image boundary.

GeoJSON coordinates are always interpreted as crop-local pixels. Any implicit
geographic CRS inference is ignored. After transformation, the geometries are
stored with an identity transformation in the sample coordinate system as:

```text
sdata.shapes["<sample>_roi_annotations"]
```

The feature ID is stored as `roi_id`. Scalar feature properties are retained.
Nested properties are serialised as compact JSON strings. A QuPath-style
`classification` object is additionally flattened into:

- `classification_name`;
- `classification_color`, represented as a hexadecimal RGB colour when
  available.

The importer currently fails if the target shape element already exists. It
does not silently overwrite an existing ROI annotation set.

Samples without `roi_annotations` retain their previous final SpatialData
unchanged. For annotated samples, the importer writes a new intermediate
`<sample>_roi_annotations.zarr`, which becomes the SpatialData published by the
main workflow under its normal final output name.

## Reporting contract

The importer produces:

```text
<sample>_roi_annotations_mqc.png
```

This image overlays the crop-local polygons on a greyscale rendering of the
capture-area H&E image, using the annotation classification colours when
available and a fallback colour otherwise. The greyscale background keeps the
histology visually distinct from the annotation classes. The legend displays
`classification_name`, or `Unclassified` when absent. Each GeoJSON feature is
labelled once with a progressive number placed inside the polygon; for a
`MultiPolygon`, the label is placed in its largest component. The plot title
also reports the total number of ROI features. Polygon exteriors and interior
rings are rendered as one compound path, so holes remain transparent and show
the underlying histology or separately annotated ROI islands.

`ASSEMBLE_SPATIAL_MULTIQC` includes the image in a per-sample `ROI annotations`
section. The plot is generated by the importer; the assembler only orders and
embeds the available per-sample images.

## Non-obvious decisions and rationale

- Keep the annotation input outside shared sample metadata so a downstream-only
  change does not invalidate expensive upstream tasks under `-resume`.
- Use the TSV as the explicit coordinate contract instead of relying on TIFF or
  GeoJSON metadata, because neither format currently carries the complete crop
  origin and SpatialData scale.
- Store ROI annotations as vector shapes rather than raster labels, preserving
  feature identity, classes, holes and overlaps.
- Keep the interface tool-independent so annotations produced outside QuPath
  remain valid when they satisfy the same GeoJSON contract.
- Run capture-area extraction even when segmentation is disabled because the
  exported image is an independent analysis input.
- Generate the ROI overlay in the import module and leave
  `ASSEMBLE_SPATIAL_MULTIQC` as a lightweight report assembler.

## Implementation and verification status

### Implemented

- `zarr_downsample_factor` in the bounds TSV and stub output;
- optional `roi_annotations` samplesheet field and dedicated keyed channel;
- duplicate annotation-path validation for repeated sample rows;
- `IMPORT_ROI_ANNOTATIONS` module and module metadata;
- polygon validation, property preservation and coordinate transformation;
- ROI shapes in the final SpatialData;
- classification-aware capture-area overlay;
- progressive ROI labels and total ROI count in the overlay;
- hole-aware rendering for `Polygon` and `MultiPolygon` annotations;
- per-sample MultiQC assembly of the ROI overlay;
- disabled direct publication of the importer outputs, because its Zarr is
  routed through final SpatialData selection and its plot is embedded in the
  assembled report;
- pass-through of unannotated samples.

### Checked statically

- Python syntax for extraction, import and report-assembly scripts;
- JSON and YAML parsing;
- Nextflow 25.10.4 DSL parsing in offline `--help` mode;
- whitespace and patch consistency with `git diff --check`.

### Functionally validated

The user has reported a successful pipeline test on the mock dataset. ROI
import, transformed SpatialData shapes, final Zarr publication and MultiQC
integration produced satisfactory results. The greyscale plot background and
disabled direct module publication were added after that test and remain to be
confirmed in the next pipeline run. Hole-aware rendering was checked directly
with the mock GeoJSON: its ROI 12 retained and displayed all 43 interior rings.

## Known limitation and open validation

- A sample that supplies `roi_annotations` but has no matching microscopy crop
  or bounds TSV currently does not reach `IMPORT_ROI_ANNOTATIONS`; the workflow
  join does not yet raise an explicit missing-match error. Before considering
  the interface fully validated, this path should produce an actionable input
  error instead of leaving the sample unannotated.
- Functional testing must confirm that the installed SpatialData version writes
  all imported properties and transformations as expected and that downstream
  readers recover the ROI shapes in the sample coordinate system.
- The visual treatment of overlapping ROIs and large numbers of annotation
  classes should be reviewed using representative user annotations.
