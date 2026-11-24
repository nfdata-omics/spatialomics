# Single-cell preprocessing QC

## Objective and scope

Add a `qc` subcommand to [`bin/sc_preprocess.py`](../../bin/sc_preprocess.py)
for AnnData (`.h5ad`) and MuData (`.h5mu`) with GEX and/or ADT. QC computes
metrics and records which observations and features pass. Keep the original
counts, cell and feature identities, and axis order in the output dataset.
Subsequent preprocessing commands use the annotations to select their analysis
inputs. This specification covers the standalone script and its optional
report files.

## QC metrics and criteria

Compute Scanpy QC metrics from the full raw count matrix of each modality,
including cells annotated as doublets. Calculate metrics and apply thresholds
across the input dataset, with the same thresholds for every sample.

| Axis    | Modality | Criterion                                       | Optional bound          |
| ------- | -------- | ----------------------------------------------- | ----------------------- |
| Cell    | GEX      | Total counts per cell                           | Minimum and/or maximum  |
| Cell    | GEX      | Genes detected per cell                         | Minimum and/or maximum  |
| Cell    | GEX      | Percentage of mitochondrial counts              | Maximum                 |
| Cell    | GEX      | Doublet label in an existing observation column | Label to exclude        |
| Cell    | ADT      | Total counts per cell                           | Minimum                 |
| Cell    | ADT      | Antibodies detected per cell                    | Minimum                 |
| Feature | GEX      | Cells detecting each gene (`n_cells_by_counts`) | Minimum (`--min-cells`) |

Bounds are inclusive. Omitted bounds do not constrain a criterion. Record a
boolean per active criterion and a combined `qc_pass` on each axis. For a
criterion with both bounds, the boolean is true when both are satisfied.
`qc_pass` is the logical AND of the active criteria on that axis. With no
active criteria, every element on that axis passes.

Calculate `n_cells_by_counts` once across the full GEX matrix. Use that value
directly for the feature criterion, including detections in doublets and
cells that fail other QC criteria. All ADT features pass. Scanpy may report
ADT feature coverage as a metric.

Keep the GEX mitochondrial, ribosomal and hemoglobin percentages from the
existing QC, based on the same gene-symbol rules. Ribosomal and hemoglobin
percentages are descriptive. Retain the modality-specific MAD diagnostics
as individual outlier flags and a combined soft flag. Use the existing
five-MAD rule on log-transformed total counts, log-transformed detected
feature counts and the percentage of counts in the top features, using the
whole modality as the reference population. Use the top 20 features where
available and retain the historical ADT fallback to the top two. The GEX
combined soft flag also includes failure of the mitochondrial bound when
active, as in the legacy diagnostic. The soft flags do not contribute to
`qc_pass` or exclusion counts.

Enable doublet QC through a CLI option naming the existing GEX observation
column (proposed `--doublet-column`). A separate option sets the label to
exclude (proposed `--doublet-label`, default `doublet`). Reject a missing
column or missing labels when this criterion is requested. Each cell with a
matching label fails this criterion and retains its original row and metadata.

## Multimodal cell membership

Keep local `obs["qc_pass"]` and `var["qc_pass"]` in each modality. For MuData,
store an additional aggregate `qc_pass` in the top-level `obs`. Match cells
by observation identity. To pass the aggregate QC, a cell must be present in
and pass every modality participating in the analysis. A cell missing from
a participating modality fails the aggregate QC. Keep that case identifiable
in the reports.

Use the shared set of admitted cells for the preprocessing steps that feed
an integrated analysis, including the modality-specific PCA fits and WNN.
Use each modality's local feature mask. On a single-modality input, use its
local masks. Inputs without QC annotations follow existing preprocessing
behavior.

## Downstream calculations and storage

- **Normalization.** Normalize after QC on the admitted cells and features,
  using the existing method selection. Keep the full raw count matrix in
  `layers["count"]`. Store the normalized values at their original positions
  in `.X` and the corresponding normalized layer, using zeros outside the
  normalized rows and columns. Keep sparse inputs sparse where the
  normalization method permits it. Persist separate row and column masks
  describing the normalization that produced these values. Downstream
  commands must use those masks to distinguish placeholders from measured
  zeros and detect when changed QC requires normalization to be rerun.
- **HVG, PCA and cell cycle.** Fit or score only on admitted cells and
  features. GEX PCA uses HVGs within the admitted feature set. Reattach
  outputs to the original axes. Use `NaN` for uncomputed cell coordinates
  and numeric cell-cycle scores, and a missing value for uncomputed phase.
- **Neighbors and UMAP.** Compute neighbors and WNN on the admitted cells,
  then store the sparse graph on the original observation axis with no
  edges incident to excluded cells. Record the cells used to build the
  graph. Run UMAP on the corresponding subgraph and store coordinates on
  the original axis, with `NaN` for excluded cells.
- **Existing plot and coordinate outputs.** Use computed coordinates for
  plots. Keep the original cell identifiers and mark missing coordinates
  in coordinate exports. Report an actionable error if the admitted
  subset is too small for the requested calculation.

## Optional QC reports

Give each file below a dedicated CLI flag that takes an output path. Write
only explicitly requested files. All reports use the same underlying QC
metrics and masks.

An optional QC CLI argument (proposed `--sample-column`) names the observation
column used to group reports. When the option is supplied, require the named
column and usable sample identifiers. When the option is omitted, treat all
observations as one sample, even when a column named `sample` exists.

| Output                      | Content                                                                                                                                          |
| --------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------ |
| GEX distribution PDF        | One page per sample with histograms of total counts per cell, detected genes per cell, mitochondrial percentage and ribosomal percentage.        |
| ADT distribution PDF        | One page per sample with histograms of total counts per cell and detected antibodies per cell.                                                   |
| Standalone interactive HTML | A self-contained per-sample view of QC outcomes and criterion failures.                                                                          |
| MultiQC interactive HTML    | The same interactive content as an embeddable `*_mqc.html` fragment.                                                                             |
| MultiQC static violin plot  | A `*_mqc.png` figure with one panel per available modality and one violin per sample, with samples side by side. Plot raw total counts per cell. |
| MultiQC summary table       | One row per sample, with applied thresholds, numbers of cells evaluated and failing per criterion, and aggregate QC pass and fail counts.        |

For the PDFs, show the distribution of all cells present in the corresponding
sample and modality, including doublets and cells that fail QC. Overlay
density curves and draw a vertical line for each active bound of the plotted
metric. The GEX feature coverage threshold has no line on these per-cell
plots. The violin plots use the same population as the PDFs.

The interactive view contains two parts for every sample:

- A stacked bar of admitted and excluded cells according to aggregate
  `qc_pass`. Together they account for all cells in the sample.
- Separate bars for failure counts for each active criterion. Show GEX total
  counts, detected genes and mitochondrial percentage separately. Show ADT
  total counts and detected antibodies separately. Include a bar for doublets
  when enabled and a bar for cells missing from each required modality when
  applicable. For a criterion with both a lower and upper bound, include
  the below-minimum and above-maximum counts in the hover details.

Criterion bars can overlap in membership. A cell can fail several criteria
and contributes once to the aggregate fail count. Count criterion failures
among cells present and evaluated in the corresponding modality. Mark inactive
or unavailable criteria in the summary table so that zero denotes an evaluated
criterion with no failures. Report the GEX feature coverage threshold in the
table as a feature criterion, with no corresponding cell-failure count.

Build the standalone document and MultiQC fragment from the same summaries
and interactive view. Embed the fragment directly using MultiQC Custom
Content and include the assets it needs to remain interactive within the
report. Use a supported Custom Content table format for the summary table.
[MultiQC Custom Content](https://docs.seqera.io/multiqc/custom_content/)

## Acceptance criteria

- QC preserves counts and axis identities and supplies local boolean masks.
  A cell or gene at an active bound passes that bound. With no thresholds
  or doublet criterion, all local masks pass and metrics are still recorded.
- Doublet QC changes the relevant cell flags. The GEX feature coverage
  metric is identical with and without doublet QC and cellular thresholds.
  MAD flags have no effect on any combined pass flag.
- For multimodal inputs with different cell orders or missing cells, the
  aggregate mask is the identity-based intersection. Excluded cells retain
  their counts, have no graph edges and have missing derived coordinates.
- Normalized values on the admitted subset match normalization of that
  subset alone. Padded zeros remain identifiable through the saved masks.
- Request each QC report independently. Compare per-sample totals and
  failure counts against the stored flags, including a cell that fails
  multiple criteria. Check the standalone HTML and the embedded MultiQC
  fragment in a browser. With no reporting flags, QC writes no reports.
- With an explicit sample column, reports group observations by its values.
  Without the option, reports show one sample for the entire input, including
  when the dataset contains a `sample` column.
