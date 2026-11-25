# MultiQC spatial previews

## Objective and scope

The H&E source and reporting task boundaries are
specified in [Final SpatialData assembly and spatial reporting](spatialdata_assembly_and_reporting.md).

Define the compact per-sample spatial previews included in MultiQC. These
figures support rapid quality review; they are not archival or detailed image
outputs. The underlying analytical results remain in the published SpatialData
object and in the existing analysis-specific outputs.

This specification covers:

- preview size and presentation conventions;
- registration, spatial QC, crop/segmentation and deconvolution panels;
- assembly of the per-sample spatial section in MultiQC;
- prevention of unnecessary image-payload duplication in MultiQC artifacts.

It does not change QC filtering, segmentation, registration, Cell2location
inference or other analytical results.

Relevant implementation files are:

- [`workflows/spatialomics.nf`](../../workflows/spatialomics.nf);
- [`modules/local/spaceranger_to_zarr/`](../../modules/local/spaceranger_to_zarr/);
- [`modules/local/assemble_spatial_multiqc/`](../../modules/local/assemble_spatial_multiqc/);
- [`modules/local/spatial_qc_plots/`](../../modules/local/spatial_qc_plots/);
- [`modules/local/spatial_image_reports/`](../../modules/local/spatial_image_reports/);
- [`modules/local/spatialdata_augment_with_segmentation/`](../../modules/local/spatialdata_augment_with_segmentation/);
- [`modules/local/cell2location/plots/`](../../modules/local/cell2location/plots/);
- [`assets/multiqc_config.yml`](../../assets/multiqc_config.yml).

## Motivation

A 24-sample production report generated with MultiQC 1.34 was approximately
948 MiB. It contained 120 embedded PNG images whose Base64 payloads occupied
approximately 404 MiB and whose decoded RGBA buffers required approximately
1.4 GiB. MultiQC also copied the custom HTML image payloads into
`llms-full.txt` and the report's `aiReportMetadataBase64`, causing a second
large representation of the same previews. The report opened but the browser
renderer crashed during navigation.

The existing figures are also visually inconsistent: layouts, margins, panel
sizes, titles and colorbars vary by producer, and the crop locator and
segmentation review are split into separate figures.

## Shared preview contract

Files named `*_mqc.*` are MultiQC previews. They are intentionally bounded
representations of data available elsewhere and must not be treated as
high-resolution deliverables.

The previews must follow these conventions:

- a shared outer canvas width, with a maximum of approximately 1200 pixels at
  the saved resolution;
- height determined by the specified layout, with a bounded height per row;
- a common two-column layout grid, with fixed outer margins and a fixed
  horizontal gutter. A two-panel row uses two equal-width panel boxes; a
  one-panel row spans both columns and the gutter, so its edges remain aligned
  with the rows above and below;
- consistent white background and vertical spacing between rows;
- a shared panel-title style: every title uses the same explicitly configured
  font family, font size and weight, and occupies the same title area rather
  than inheriting producer-specific plotting defaults;
- legends and colorbars are placed below the corresponding plot area and must
  not change its horizontal position or width. Colorbars are horizontal; a
  row-level colorbar shared by two panels spans the aligned row beneath them;
- image aspect ratios are preserved within the aligned panel boxes rather than
  stretching image content to fill them;
- no sample-level figure title because MultiQC already identifies the sample;
- no coordinate ticks or grayscale colorbars unless they convey information
  required to interpret the panel;
- human-readable labels rather than internal column or element names;
- image encoding selected for compact browser delivery while keeping labels,
  outlines and diagnostic image structure legible;
- a target encoded preview budget of approximately 1--2 MiB per sample for the
  standard registration, QC, segmentation-review and deconvolution set;
- lazy-loaded HTML images;
- no silent omission of user-requested crop areas.

The spatial previews do not provide useful LLM input. The supported upstream
mechanism that can prevent their inclusion is to disable all MultiQC AI
features with `no_ai: true`; there is no section-level option for excluding
only the spatial HTML from AI metadata. MultiQC 1.34, currently pinned by the
pipeline, still builds `llms-full.txt` and embeds `aiReportMetadataBase64` even
when `no_ai` is enabled. This is an upstream bug, not a custom-content
configuration problem.

Upstream pull request MultiQC/MultiQC#3715 makes `no_ai: true` skip generation
of the AI report metadata, pseudonym map and `llms-full.txt`. If removing the AI
toolbox is accepted, the implementation should use the first stable MultiQC
release containing that fix and enable `no_ai: true`; upgrading to a released
fix is preferred to vendoring an unreleased MultiQC patch or rewriting the
generated report. Until such a release is available, bounding the source
previews still reduces both copies substantially, but does not eliminate the
duplication in MultiQC 1.34.

## Per-sample panel design

### Registration and Space Ranger tissue detection

This preview uses one row with two equally sized panel boxes. The left panel
shows registration and the right panel shows the tissue area selected by Space
Ranger. The panels have distinct purposes and do not reuse rendered Space
Ranger QC images.

The registration panel is a composite in the sample coordinate system:

1. the reduced grayscale H&E from the samplesheet, read from the final Zarr;
2. the registered CytAssist image shown with a strong pseudocolor and partial
   transparency so tissue agreement and disagreement are immediately visible;
3. the capture-area boundary shown as an unfilled, high-contrast outline.

Its viewport is the complete transformed H&E image extent, including H&E
regions outside the capture area. The transformed CytAssist layer and the
capture-area outline are overlaid within that viewport and clipped only at the
H&E image boundary; neither the capture-area bounds nor the CytAssist bounds
may be used to crop the panel.

The panel must use image layers and transformations declared in SpatialData. It
must not read a rendered Space Ranger QC image or infer a transform from image
dimensions, filenames or visual similarity. The plotting code must remain
independent of the tool that originally produced the registered SpatialData.

The required general input contract is therefore:

- a microscopy image layer transformed into the sample coordinate system;
- a distinct CytAssist/reference image layer transformed into the same sample
  coordinate system;
- capture-area geometry transformed into that coordinate system.

Space Ranger stores the final registration in the `metadata_json` attribute of
`feature_slice.h5`. It provides the projective transformations between
CytAssist pixel coordinates, spot-array coordinates and microscope-image pixel
coordinates:

- `cytassist_colrow_to_spot_colrow` and its inverse;
- `spot_colrow_to_microscope_colrow` and its inverse.

The pinned `spatialdata_io.visium_hd()` reader already consumes these matrices.
With `load_all_images=True`, it loads `spatial/cytassist_image.tiff`, composes
the CytAssist-to-spot and spot-to-microscope transformations, performs the
required projective resampling and adds `<sample>_cytassist_image` to
SpatialData with a transformation into the sample coordinate system. The Space
Ranger ingestion path must enable this option. The downstream plotting code
then consumes only the resulting SpatialData layers and transformations; it
does not depend on Space Ranger paths or transformation metadata.

Because the filtered expression matrix contains only included bins, their
bounds do not necessarily represent the full capture area. Ingestion therefore
derives one compact `<sample>_capture_area` polygon from the complete 16 µm
`tissue_positions.parquet` grid and its Space Ranger bin scale. This preserves
the capture-area extent without loading the raw expression matrix or retaining
all out-of-tissue bins as SpatialData shapes.

The tissue-detection panel shows the H&E image cropped to the capture area,
with the bins classified as `in_tissue=1` by Space Ranger rendered as a
semi-transparent, high-contrast fill. Its opacity must keep the underlying
tissue morphology clearly visible. This panel is for checking whether Space
Ranger included background or excluded apparent tissue; it is not a continuous
pixel-wise tissue segmentation.

The tissue classification must come from the unmodified Space Ranger
`in_tissue` field at the resolution used by the pipeline. It is shown before
the pipeline's own QC filtering and must not be derived from or visually
conflated with `qc_pass`, expression intensity, global-outlier filtering or
segmentation results. When the filtered Space Ranger matrix exposes only the
positive class, the footprint of those included bins defines the overlay; the
full capture-area H&E remains visible as the reference against which omitted
regions are assessed.

For inputs without a CytAssist image or complete non-null registration
matrices, the ingestion path does not fabricate a CytAssist layer. The tissue
detection panel remains available when its H&E, capture-area and `in_tissue`
inputs are present, while the registration panel is omitted with an explicit
reason recorded in the process output or log. The pipeline must never present
two microscopy layers as a CytAssist-versus-H&E registration check. When only
the tissue-detection panel is rendered, it spans the full two-column grid while
retaining its capture-area viewport.

### Spatial QC overlays

The spatial QC preview uses three rows and two columns:

| Row                    | Left column                   | Right column    |
| ---------------------- | ----------------------------- | --------------- |
| Total counts           | In-tissue before QC filtering | QC-passing bins |
| Detected genes         | In-tissue before QC filtering | QC-passing bins |
| Mitochondrial fraction | In-tissue before QC filtering | QC-passing bins |

Each row uses one shared colormap, normalization and colorbar. Robust limits
are calculated from the unfiltered in-tissue values for that metric and applied
unchanged to both panels, making before/after colors directly comparable.
The shared colorbar is horizontal and is placed below both plot areas.

All six plot viewports use the same bounds: the capture-area bounding box
expanded on every side by an absolute padding equal to 2% of the larger of its
width and height. This padding is computed in the plotted coordinate system;
data outside the padded bounds do not enlarge the viewport. Out-of-tissue
observations are transparent. Existing QC flags, thresholds and the `qc_pass`
definition are unchanged.

### Crop location and segmentation review

The current crop-locator and segmentation figures are replaced by one preview.
Each requested crop occupies one row with two equally sized panel boxes:

- left: the H&E capture-area view with that crop shown as an unfilled outline;
- right: a magnified view of the crop with segmented-cell boundaries overlaid
  in a high-contrast color.

Rows use concise labels such as `Crop 1 — location` and
`Crop 1 — segmentation`. Numeric crop coordinates may appear in a compact
caption but not as the main panel title. Grayscale intensity colorbars and
coordinate ticks are omitted. Multiple requested crops add rows; they are not
silently truncated.

This combined preview replaces both
`*_crop_areas_downsampled_microscopy_mqc.png` and
`*_segmentation_crop_panels_mqc.png` with
`*_crop_segmentation_mqc.png` in the MultiQC input contract.

### Cell2location deconvolution

The compact deconvolution preview retains two aligned panels:

- total estimated abundance with its existing robust normalization and an
  explicit colorbar;
- selected cell-type composition with a compact categorical legend.

The sample-level suptitle is removed. Short panel titles remain, and margins,
axis boxes and typography follow the shared preview contract. The detailed PDF
plots and quantitative MultiQC tables remain unchanged.

### Optional ROI preview

When ROI annotations are present, their preview follows the same pixel,
typography and margin constraints. This work does not change the ROI coordinate
or data contracts.

## Output and assembly contract

The workflow continues to associate preview inputs by sample ID and must not
rely on channel emission order. The assembler accepts only non-empty previews,
orders them consistently and emits one custom-content HTML fragment per sample.

The standard order is:

1. registration and Space Ranger tissue detection, including whichever panels
   have their required inputs;
2. spatial QC overlays;
3. crop and segmentation review, when segmentation is available;
4. ROI annotations, when supplied;
5. Cell2location deconvolution, when enabled.

Removing or combining preview files must be propagated consistently through
module outputs, workflow channels, module metadata, MultiQC search patterns,
tests and user-facing output documentation where applicable.

## Implementation status

Implemented, except for the explicitly open raster-encoding decision and the
MultiQC AI-metadata behavior that requires an upstream release or a separate
decision. The implementation:

- conditionally imports the registered CytAssist layer and preserves a compact
  capture-area geometry during Space Ranger ingestion;
- renders registration and tissue selection from the assembled Zarr for
  samples with a samplesheet H&E;
- renders the 3-by-2 spatial QC preview from the QC-annotated Zarr;
- replaces the two crop/segmentation previews with one row-oriented preview;
- regularizes the Cell2location preview and adds its abundance colorbar;
- updates workflow channels, module metadata, MultiQC assembly and user-facing
  output documentation for the new contracts;
- lazy-loads spatial images in the assembled MultiQC custom HTML.

## Verification status

The production report was inspected statically to quantify report size,
embedded image payloads and decoded pixel buffers. Existing representative
figures were inspected to establish the presentation problems.

Before registered CytAssist ingestion was enabled, the raw and augmented mock SpatialData Zarr metadata were inspected directly:

- the raw conversion contains high- and low-resolution tissue images but no
  CytAssist image layer with the reader arguments used at that time;
- the augmentation adds full and downsampled microscopy layers in the sample
  coordinate system;
- the registration figure at that stage compared two representations of
  the microscopy image rather than microscopy and CytAssist modalities.

The Space Ranger mock `feature_slice.h5` metadata contains all four
CytAssist/spot/microscope transformation matrices. A read-only invocation of
the pinned `spatialdata_io` 0.6.0 container with `load_all_images=True`
successfully produced `<sample>_cytassist_image` and attached its
transformation into the sample coordinate system. The mock 16 µm table also
confirmed that `in_tissue` is retained in `obs`; because this reduced filtered
mock contains only included bins, all observed values are `1`.

These metadata-only checks establish the ingestion contract but do not by
themselves constitute visual validation of registration or tissue detection.
No functional pipeline or Nextflow stub run has been performed for this work.

A targeted containerized check, not a Nextflow run, used the local reduced
Space Ranger mock to exercise ingestion and the plotting scripts. It confirmed
that the registered CytAssist and capture-area elements survive Zarr
serialization and produced:

- a 1120-by-539 registration/tissue-selection PNG of approximately 1.1 MiB;
- a 1091-by-1522 spatial-QC PNG of approximately 1.8 MiB;
- a 992-by-412 single-row crop/segmentation PNG of approximately 520 KiB.

The newly generated registration, QC and crop/segmentation previews were
inspected visually. The registration check used the real mock registration;
the segmentation check used a synthetic temporary mask and therefore validates
rendering and layout only, not segmentation quality. The Cell2location change
was checked statically but its revised preview was not generated functionally.
The assembler placed these three previews under their expected headings and
added lazy image-decoding attributes; the resulting Base64 HTML fragment was
approximately 4.3 MiB.
The provisional PNG sizes do not yet establish the final per-sample encoded
budget, which depends on the open encoding decision and representative
production images.

The later layout refinement (shared fixed-width grid, full-H&E registration
viewport, capture-area QC bounds and below-panel colorbars/legends) has been
implemented and checked for Python syntax and imports in the pinned project
container. The previews have not been regenerated or visually inspected after
this refinement, and no Nextflow or stub run was performed.

The official MultiQC AI documentation and the source code in the pinned 1.34
container were inspected. The documented `no_ai` option removes the AI user
interface, but version 1.34 calls AI metadata and prompt generation
unconditionally. Upstream pull request MultiQC/MultiQC#3715 addresses this exact
behavior and is still open, so no released version containing the fix has yet
been identified.

## Open questions

- Confirm whether removing the MultiQC AI toolbox is acceptable. Retaining AI
  while excluding only spatial custom content is not supported by MultiQC 1.34
  or by the currently proposed upstream fix.
- Select the compact raster encoding and quality settings after comparing
  legibility and output size on representative microscopy and categorical
  panels.
