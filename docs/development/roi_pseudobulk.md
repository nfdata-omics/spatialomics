# ROI pseudobulk from Visium HD counts

## Objective and scope

Produce an AnnData count matrix for downstream statistical analysis, with one
observation per annotated GeoJSON feature and one variable per gene. Each row
contains the sum of the Space Ranger UMI counts from bins covered by that ROI
at the selected resolution. The default uses native 2 µm bins. Features with
the same classification retain separate rows. A `MultiPolygon` is a single
feature and produces one row.

The existing annotation and coordinate contract applies. The pseudobulk branch
starts from samples that complete `IMPORT_ROI_ANNOTATIONS`, using their imported
geometries and the bin table from the annotated SpatialData.
[ROI annotation exchange](roi_annotation_exchange.md)

## Source counts and bin coordinates

Use the table at `tables/<bin_label>` in the annotated Zarr. The run-wide
`--roi_pseudobulk_bin_size` parameter selects the table independently of
`bin_analysis_bin_size`. It defaults to `square_002um` and accepts Space Ranger
bin labels such as `square_008um` and `square_016um`. If the table is absent,
report the sample ID and the requested table.

Use the original UMI counts from the Space Ranger filtered feature-barcode
matrix loaded by `SPACERANGER_TO_ZARR`. The table contains the bins retained
by Space Ranger's tissue selection. Sum its `Gene Expression` features using
an integer accumulator. Include these bins regardless of the pipeline's
`qc_pass` flag, without applying further expression thresholds or normalization.
[Space Ranger feature-barcode matrices](https://www.10xgenomics.com/support/software/space-ranger/latest/tutorials/outputs/space-ranger-feature-barcode-matrices)

Use the full-resolution (x, y) centres in the table's `obsm["spatial"]`,
already aligned to the count rows. Divide by `roi_zarr_downsample_factor`,
the factor that ROI import copied from the extraction TSV into the annotated
Zarr attributes. The importer has already applied the crop offset and
crop-local scales to the ROI geometries.

## Bin assignment and aggregation

- Assign the entire bin to every ROI whose geometry covers its centre. The
  polygon boundary is included, including the boundary of a hole. The interior
  of a hole is excluded.
- Coarser bins can move counts across ROI boundaries because their full counts
  are assigned by the centre. This can affect small or narrow ROIs in particular.
- Each bin contributes at most once to a given ROI. A bin covered by multiple
  ROI features contributes independently to each row, including when the
  centre is on a shared boundary. Record the number of shared bins per ROI.
- Keep every imported ROI in the result. A ROI with no measured bins has an
  all-zero count row, zero included bins and `roi_status = no_bins`. A ROI with
  measured bins and zero gene-expression UMIs can be identified separately by
  its positive bin count and zero `total_counts`.

Use sparse matrices throughout. A spatial index can restrict geometry checks to
candidate bins, and a sparse ROI-by-bin membership matrix can aggregate counts
without materialising a dense bin-by-gene matrix or dense membership masks.

## AnnData contract

`X` is a sparse ROI-by-gene matrix of integer UMI sums. Use the Space Ranger gene
ID as `var_names` and retain the gene label from the SpatialData table in `var`. Gene IDs must be unique.

Use a collision-free encoding of `(sample_id, roi_id)` for `obs_names` and
retain both identifiers as columns. Duplicate imported `roi_id` values within a
sample are an input error. The imported ROI ID links the row back to the
published spatial geometry.

Keep the sample annotations available from the input, including slide and
capture area where supplied. The current samplesheet schema has no dedicated
fields for biological covariates such as donor, condition, or batch.
[Samplesheet schema](../../assets/schema_input.json)

The GeoJSON supplies two kinds of analytical annotations for the pseudobulk
output.

- The ROI classification comes from the feature's `properties.classification`.
  The importer reads `name` for a simple class or joins the ordered `names`
  components with `: ` for a derived class. AnnData reads `classification_name`
  from the imported shape. An unclassified ROI has a missing value.
  [QuPath class serializer](https://github.com/qupath/qupath/blob/v0.7.0/qupath-core/src/main/java/qupath/lib/io/GsonTools.java#L433-L450)
- The QuPath project image metadata comes from the optional
  `sample_metadata` object described below.

| `obs` field       | Meaning                                                                        |
| ----------------- | ------------------------------------------------------------------------------ |
| `n_bins_included` | Number of table barcodes whose centres are covered by the ROI.                 |
| `n_bins_shared`   | Number of included bins also assigned to at least one other ROI in the sample. |
| `total_counts`    | Sum of gene-expression UMIs across the ROI row.                                |
| `roi_status`      | `ok` when at least one bin is included, otherwise `no_bins`.                   |

Record the source matrix, bin resolution, and spatial assignment rule in `uns`.
In the combined file, retain sample-specific provenance keyed by `sample_id`.

### QuPath project image metadata

Use the existing `roi_annotations` input to carry additional sample-level
covariates. Add an optional top-level `sample_metadata` object alongside
`features` in the GeoJSON `FeatureCollection`. GeoJSON allows such additional
members.
[GeoJSON RFC 7946 section 6.1](https://datatracker.ietf.org/doc/html/rfc7946#section-6.1)

For example, `sample_metadata` could contain `{"donor_id": "D01", "condition":
"treated", "batch": "B02"}`. `IMPORT_ROI_ANNOTATIONS` reads the object when
it imports the polygons and stores it once in
`sdata.attrs["roi_sample_metadata"]`. The same Zarr write persists the
polygons and the sample metadata, making both available in the published
SpatialData.
[SpatialData global attributes](https://spatialdata.scverse.org/en/stable/api/SpatialData.html#spatialdata.SpatialData.attrs)

The pseudobulk step reads the sample metadata and ROI classification directly
from the annotated Zarr. It repeats sample-level covariates on every ROI
row of the sample in namespaced `obs` columns such as `sample_metadata.donor_id`.
The samplesheet remains the source of the sample identifier. If the metadata
object includes `sample_id`, it must match the samplesheet. GeoJSON files
without the object remain valid, and missing covariates stay missing in the
combined dataset.

QuPath 0.7.0 exports ROI classifications through **File → Export objects as
GeoJSON**. The command does not include metadata from the project image entry.
To include those values without QuPath scripting, export as a `FeatureCollection`
and copy the values from the corresponding row in **File → Project… → Edit
project metadata** into the top-level `sample_metadata` object in the exported
file. Supply that GeoJSON through `roi_annotations`.
[QuPath GeoJSON export](https://github.com/qupath/qupath/blob/v0.7.0/qupath-gui-fx/src/main/java/qupath/lib/gui/commands/ExportObjectsCommand.java#L113-L160)
· [QuPath metadata editor](https://qupath.readthedocs.io/en/stable/docs/reference/commands.html#edit-project-metadata)

## Workflow integration and publication

Run per-sample aggregation downstream of `IMPORT_ROI_ANNOTATIONS` and
independently of segmentation. Read counts, bin coordinates, ROI geometries
and annotation metadata from the annotated Zarr. Write the aggregated data as
AnnData. The final published SpatialData also includes microscopy and
segmentation layers when that augmentation runs.

The per-sample publication requirement is
`roi_pseudobulk/<sample>_roi_pseudobulk.h5ad` for each annotated sample. The
current configuration keeps those files as workflow intermediates.

When at least one per-sample result exists, publish the combined dataset at
`pseudobulk/roi_pseudobulk.h5ad`. Concatenate observations and align genes by
ID, including when inputs list the same genes in different orders. Require
identical sets of gene IDs across samples for the combined file and report
mismatches explicitly. This avoids treating an unavailable feature as observed
zero expression. Only the combined output requires collecting per-sample
results. Adding an annotation on a resumed run should reuse completed upstream
tasks under the existing annotation input contract.

Publish the combined counts and annotations in the same directory as
`roi_pseudobulk_counts.tsv.gz` (`X`), `roi_pseudobulk_rois.tsv` (`obs`) and
`roi_pseudobulk_genes.tsv` (`var`). The counts and ROI annotations use the
same `obs_id` from `obs_names`. The counts use `var_names` for gene column
headers, matching `gene_id` in the gene annotations. Preserve the order and
integer counts from the combined H5AD. Write missing annotations as empty fields.

The resulting rows describe individual ROI features. For downstream inference,
ROI observations from the same biological replicate should be modelled as such.
When pooling overlapping ROIs into a single region, aggregate over the union
of their bins to count each bin once.
