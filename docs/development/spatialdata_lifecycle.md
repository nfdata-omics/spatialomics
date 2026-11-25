# SpatialData lifecycle

## Objective and scope

Define how a per-sample SpatialData object is created, transformed, selected and
published across the pipeline. This specification records shared invariants for
all modules that read or write SpatialData, including serialization settings
that may otherwise appear to be local implementation details.

Feature-specific element contracts remain in their dedicated specifications:

- [`cell2location_deconvolution.md`](cell2location_deconvolution.md);
- [`roi_annotation_exchange.md`](roi_annotation_exchange.md).

Relevant implementation files are:

- [`workflows/spatialomics.nf`](../../workflows/spatialomics.nf);
- [`main.nf`](../../main.nf);
- [`modules/local/spaceranger_to_zarr/`](../../modules/local/spaceranger_to_zarr/);
- [`modules/local/spatial_quality_control/`](../../modules/local/spatial_quality_control/);
- [`modules/local/cell2location/deconvolution/`](../../modules/local/cell2location/deconvolution/);
- [`modules/local/cell2location/reinsert/`](../../modules/local/cell2location/reinsert/);
- [`modules/local/spatialdata_augment_with_segmentation/`](../../modules/local/spatialdata_augment_with_segmentation/);
- [`modules/local/import_roi_annotations/`](../../modules/local/import_roi_annotations/).

## Lifecycle and canonical selection

The final assembly and reporting contracts are specified in
[Final SpatialData assembly and spatial reporting](spatialdata_assembly_and_reporting.md).

Each sample has one logical SpatialData lineage. Transformations write new Zarr
directories rather than modifying an upstream task output in place:

```text
Space Ranger output
        |
        v
raw SpatialData Zarr -----------------------> raw .zarr.tar.gz archive
        |
        v
QC-annotated SpatialData
        |
        +--> Cell2location SpatialData, when a reference is supplied
        |              |
        +--------------+
                       v
        ROI annotation augmentation, when supplied
                       |
                       v
        samplesheet microscopy and optional segmentation assembly
                       |
                       v
        final published SpatialData directory
```

The canonical object at each decision point is the most advanced successful
version available for that sample:

1. `SPACERANGER_TO_ZARR` creates the raw SpatialData from Space Ranger output,
   including the registered CytAssist image when Space Ranger provides both the
   image and complete transformation matrices. It also preserves a compact
   capture-area polygon derived from the complete 16 µm position grid.
2. `SPATIAL_QUALITY_CONTROL` creates the mandatory QC-annotated version used by
   downstream analyses.
3. `CELL2LOCATION_DECONVOLUTION` trains one joint spatial model when a single-cell
   reference is supplied, and `CELL2LOCATION_REINSERT` creates the deconvolved
   per-sample SpatialData versions. Without a reference, the QC versions continue
   unchanged.
4. `IMPORT_ROI_ANNOTATIONS` augments the QC or Cell2location version for
   samples with a matching GeoJSON, capture-area crop and coordinate TSV.
   Unannotated samples retain the preceding version.
5. `SPATIALDATA_AUGMENT_WITH_SEGMENTATION` writes the final version for
   samples with a samplesheet H&E, adding the native-resolution microscopy
   and its reduced preview. It adds segmentation labels when a mask is
   generated. Samples without a samplesheet H&E retain the preceding Zarr.
6. The workflow publishes exactly the selected final version.

Optional branches are associated by sample ID. Channel emission order must
never be used to associate SpatialData versions with samples.

`EXTRACT_VISIUM_REGION` reads the raw SpatialData because it only needs the
original Visium geometry and runs independently of downstream analytical
branches. `BIN2CELL` currently produces separate cell-level outputs and does not
modify the final SpatialData lineage.

ROI aggregation reads the annotated SpatialData before microscopy augmentation.
The segmentation backend produces masks independently of ROI import and
pseudobulk. Final assembly consumes the selected Zarr, the prepared microscopy
image and any generated mask. Reporting reads the resulting Zarr as an
immutable input.

## Storage and publication contract

### Raw SpatialData

The raw uncompressed Zarr is an internal workflow object and remains available
through the Nextflow work directory and cache. It is archived as
`<sample>.zarr.tar.gz` and published under `count/` for completeness and
long-term recovery. The archive is not intended as the normal downstream
analysis interface. Registered image layers use the transformations produced by
Space Ranger through `spatialdata_io`; the pipeline does not reconstruct image
registration from rendered Space Ranger plots.

### Intermediate SpatialData

QC, Cell2location, segmentation and ROI Zarr directories are process outputs.
They are passed to subsequent processes and retained according to the Nextflow
work/cache lifecycle, but are not independently published as competing final
datasets. Associated CSV files, plots and reports may be published separately.

### Final SpatialData

The selected final object is published as an unpacked directory:

```text
spatialdata/<sample>_spatialdata.zarr/
```

The directory form is the supported downstream interface because SpatialData
can read it directly without first creating an additional unarchived copy. No
separate final-Zarr archive is produced.

Publication does not establish ownership of intermediate work directories.
Deleting the Nextflow work directory or cache can prevent `-resume` from reusing
intermediate SpatialData even when a published final copy still exists.

## Mutation and serialization contract

### New output for every transformation

A module that augments SpatialData must:

1. read its declared input Zarr;
2. preserve unrelated elements, coordinate systems, tables and metadata;
3. add or replace only the elements defined by that module's contract;
4. write a distinct output Zarr in its task work directory;
5. emit `[meta, zarr]` with the original sample identity.

In-place modification of an upstream Zarr is not allowed because it would break
Nextflow task isolation, caching and reproducibility. Rewriting a SpatialData
object may copy substantial data, but this cost is currently accepted to retain
immutable task outputs.

Output replacement inside a task must be explicit. Modules must not silently
overwrite a pre-existing analytical element unless overwrite behaviour is part
of their public contract.

### AnnData nullable strings

Calling `SpatialData.write()` serializes the complete object, including AnnData
tables that a module did not modify. A writer that only adds an image, label or
shape can therefore fail while rewriting a pre-existing table.

With the AnnData version used by the pipeline containers, Pandas nullable string
indices or columns require explicit opt-in before writing:

```python
import anndata

anndata.settings.allow_write_nullable_strings = True
sdata.write(output_zarr)
```

The setting is a process-local serialization requirement. It does not change
the analytical values or justify casting indices and columns to another dtype.
An error such as the following is therefore a writer configuration problem, not
evidence that the reported table or ROI element is scientifically invalid:

```text
RuntimeError: anndata.settings.allow_write_nullable_strings is False
Error raised while writing key '_index' ... to /tables/<table>/obs
```

Every custom Python module that calls `SpatialData.write()` on an object that
may contain AnnData tables should enable this setting immediately before the
write. This includes late augmentation modules because they rewrite all tables
inherited from earlier lifecycle stages.

### Coordinate systems and element identity

Writers must preserve existing coordinate systems and transformations. New
spatial elements must declare their transformation to the intended sample
coordinate system explicitly. Element names must be stable and collisions must
be handled deliberately; silently renaming an element would change the public
data contract.

Tables must retain row identity and alignment when new columns or matrices are
added. Filtering for an analysis must not silently remove rows from the
canonical table unless that behaviour is explicitly specified. Sample IDs in
element names, coordinate systems, channel metadata and published filenames
must remain consistent.

## Caching and optional augmentation

Downstream-only inputs must not be inserted into the shared sample `meta` map.
They should travel in independently keyed channels so adding an optional input
does not invalidate upstream SpatialData tasks under `-resume`.

An optional augmentation follows this pattern:

1. join its required inputs by sample ID;
2. run only for samples with a complete input set;
3. join its output back to the preceding canonical SpatialData channel with
   `remainder: true`;
4. select the augmented output when present and the preceding object otherwise.

Missing required inputs for an explicitly requested augmentation should produce
an actionable error. They must not silently turn the requested analysis into a
pass-through. Empty optional-input channels must still allow all unaffected
samples to reach final publication.

## Current implementation status

### Implemented

- raw conversion from Space Ranger and archival publication, including
  conditional registered CytAssist import and compact capture-area geometry;
- mandatory QC annotation as the downstream baseline;
- optional Cell2location replacement of the QC version;
- optional segmentation and microscopy augmentation with keyed fallback;
- optional ROI augmentation with keyed fallback;
- unpacked final Zarr publication through the workflow output block;
- explicit nullable-string write opt-in in segmentation/microscopy augmentation;
- explicit nullable-string write opt-in in ROI annotation import.

### Known consistency gaps

- `SPATIAL_QUALITY_CONTROL` and `CELL2LOCATION_REINSERT` call
  `SpatialData.write()` without explicitly enabling
  `anndata.settings.allow_write_nullable_strings`. A local Docker end-to-end run
  completed both stages successfully on 2026-11-04, including reinsertion of the
  Cell2location results. This validates the observed mock-dataset dtypes, but not
  every possible input dtype; an explicit shared writer invariant may still be
  preferable if nullable strings are introduced by other datasets.
- `SPACERANGER_TO_ZARR` also writes SpatialData without the explicit setting.
  Its object originates directly from `spatialdata_io`, so its exact AnnData
  dtypes and compatibility should be checked separately.
- ROI import currently uses inner joins for its required crop and bounds inputs.
  A requested annotation without a matching export contract can therefore be
  omitted instead of producing the required actionable error.
- The pipeline does not currently perform a lightweight reopen check after each
  Zarr write. Successful completion of `write()` is the present persistence
  check.

## Verification status

The raw, QC, Cell2location and segmentation lifecycle has been exercised during
pipeline development. Final publication as an unpacked Zarr directory has been
implemented and syntactically checked.

The user initially reported a nullable-string serialization failure in
`IMPORT_ROI_ANNOTATIONS`. After enabling the AnnData opt-in used by the existing
segmentation writer, a complete mock-dataset pipeline run successfully wrote
and published the ROI-augmented final Zarr. The result and MultiQC integration
were reported as satisfactory.

## Future validation

- Confirm that every lifecycle stage can write and reopen representative
  SpatialData containing nullable string indices and columns.
- Decide whether the nullable-string opt-in should remain explicit in every
  writer or be enforced through a small shared repository convention without
  obscuring the write operation.
- Add focused validation that the selected final Zarr contains the expected
  elements from every enabled upstream augmentation.
- Measure the storage and I/O cost of full-object rewrites before considering
  incremental or in-place SpatialData APIs; any alternative must preserve task
  isolation and resumability.
