# Cell2location deconvolution

## Objective and scope

Integrate Cell2location deconvolution into the SpatialData-centred workflow while:

- training the single-cell reference once per pipeline run;
- deconvolving every spatial sample independently;
- retaining posterior results and provenance in SpatialData;
- separating expensive inference from lightweight plotting and reporting;
- keeping detailed diagnostics outside MultiQC and a concise summary inside it.

The feature currently supports one single-cell reference dataset per pipeline
run. Support for multiple reference datasets is not yet implemented.

Relevant implementation files are:

- [`main.nf`](../../main.nf);
- [`workflows/spatialomics.nf`](../../workflows/spatialomics.nf);
- [`modules/local/cell2location/reference/`](../../modules/local/cell2location/reference/);
- [`modules/local/cell2location/deconvolution/`](../../modules/local/cell2location/deconvolution/);
- [`modules/local/cell2location/plots/`](../../modules/local/cell2location/plots/);
- [`modules/local/assemble_spatial_multiqc/`](../../modules/local/assemble_spatial_multiqc/).

## Agreed behaviour

### Activation and data flow

The presence of `single_cell_reference_h5ad` enables the complete Cell2location
branch; there is no separate skip flag. Without a reference, the QC-annotated
SpatialData continues directly downstream.

With a reference:

1. `CELL2LOCATION_REFERENCE` trains one reference model and exports signatures;
2. `CELL2LOCATION_DECONVOLUTION` uses those signatures for every spatial sample;
3. `CELL2LOCATION_PLOTS` produces diagnostics, biological plots and MultiQC inputs;
4. the deconvolved SpatialData becomes the input to segmentation augmentation;
5. the most advanced available SpatialData is emitted and published by the workflow.

The spatial model is not saved separately. Nextflow caching is relied upon to
avoid repeating successful GPU inference. Plotting remains a separate process
so plot changes do not invalidate deconvolution tasks.

### Reference contract

The reference is an H5AD containing raw counts and cell-type annotations. The
pipeline interface is intentionally tool-independent:

| Parameter                             | Default     | Contract                                  |
| ------------------------------------- | ----------- | ----------------------------------------- |
| `single_cell_reference_h5ad`          | `null`      | Existing reference H5AD                   |
| `single_cell_reference_cell_type_key` | `cell_type` | Cell-type column in `obs`                 |
| `single_cell_reference_batch_key`     | `none`      | Batch column in `obs`; `none` disables it |
| `single_cell_reference_counts_layer`  | `X`         | Raw counts in `X` or a named layer        |
| `single_cell_reference_gene_key`      | `var_names` | Gene IDs in `var_names` or a `var` column |

The reference model is trained once. Its signatures CSV is published under
`${params.outdir}/deconvolution/` and broadcast to all spatial samples. The
trained reference model remains in the process output/cache but is not currently
published.

### Spatial input and inference

The deconvolution input is the QC-annotated SpatialData table selected by
`bin_analysis_bin_size` (`square_016um` by default). Only locations with
`obs['qc_pass'] == true` are modelled; other locations remain in the table with
missing posterior values.

Spatial counts are read from `X`. `spatial_gene_key` defaults to `var_names`.
Reference and spatial gene IDs must be unique and at least 1,000 genes must be
shared by default. Failure messages report examples and alternative overlaps
when `var_names` or `gene_ids` are available.

Spatial-model training currently uses a default batch size of 35,000 locations,
chosen for the facility GPUs. It remains overridable with module-level
`--batch-size`. Posterior sampling uses the independent
`--posterior-batch-size` setting, currently 512, to control host-memory use.

## SpatialData contract

The selected table stores labelled cell-type abundance matrices in `obsm`:

- `means_cell_abundance_w_sf`;
- `stds_cell_abundance_w_sf`;
- `q05_cell_abundance_w_sf`;
- `q95_cell_abundance_w_sf`.

Rows are aligned to the complete spatial table. Locations excluded by QC contain
`NaN`. No permanent `c2l_<celltype>` or total-abundance columns are added to
`obs`, because choosing one posterior estimator is analysis-dependent.

When available, detection-sensitivity summaries are stored in `obs` as:

- `c2l_detection_sensitivity_mean`;
- `c2l_detection_sensitivity_std`;
- `c2l_detection_sensitivity_q05`;
- `c2l_detection_sensitivity_q95`.

`uns['cell2location']` stores provenance, model settings, cell-type names,
location and gene counts, training history, ELBO summaries, posterior-sampling
settings and compact reconstruction-QC data or its error state.

## Output and reporting contract

### Deconvolution outputs

`CELL2LOCATION_DECONVOLUTION` produces:

```text
<sample>_cell2location.zarr/
<sample>_cell2location_abundance_mean.csv
<sample>_cell2location_abundance_q05.csv
```

The CSV files are published under `${params.outdir}/deconvolution/`. The Zarr is
passed downstream and included in the final published SpatialData rather than
published independently by the module.

### Selected abundance statistic

`cell2location_abundance_statistic` selects the point estimate used for
biological plots and tables:

- `mean` (default): `means_cell_abundance_w_sf`;
- `q05`: `q05_cell_abundance_w_sf`.

Posterior diagnostic plots are independent of this selection. Total abundance
is the sum across cell types. Fractions are abundance divided by total abundance;
zero-total locations remain masked. Whole-slide composition is calculated from
summed abundance, not from the unweighted mean of location-level fractions.

When `q05` is selected, derived totals and fractions are normalised marginal
`q05` estimates, not posterior fifth percentiles of those derived quantities.

### Per-sample plot outputs

`CELL2LOCATION_PLOTS` publishes under
`${params.outdir}/deconvolution/<sample>_cell2location/`:

- training history CSV and PNG;
- reconstruction-accuracy PNG, or an error text file when unavailable;
- a compact diagnostics TSV;
- posterior mean-versus-`q05` hexbin PDF;
- spatial relative posterior interval-width PDF;
- cell-type abundance PDF;
- cell-type fraction PDF;
- total-abundance distribution and spatial-map PDF;
- compact deconvolution PNG for MultiQC;
- composition and summary custom-content TSV files for MultiQC.

Available tissue imagery is rendered in greyscale. Abundance panels share one
scale based on the largest cell-type 99th percentile. Fraction panels use an
independent 99th-percentile maximum. Relative posterior interval-width panels
share one robust scale and mask posterior means at or below
`cell2location_uncertainty_min_abundance` (`0.1` by default).

### MultiQC

MultiQC receives:

- a cross-sample bar plot of whole-slide cell-type composition;
- a cross-sample table with modelled locations, total abundance, and mean and
  median abundance per location;
- one compact per-sample image showing total abundance and the composition of up
  to three selected cell types.

`cell2location_multiqc_cell_types` optionally selects up to three types. Its
default is `null`, which selects the three most abundant types independently per
sample. Selected fractions are mixed using pink, cyan and yellow; all remaining
fractions contribute neutral grey and are not renormalised away.

`ASSEMBLE_SPATIAL_MULTIQC` groups independently available QC, imaging,
segmentation and deconvolution images by sample. It publishes
`<sample>_spatial_mqc.html` under `${params.outdir}/spatial_plots/`. Missing
segmentation or deconvolution content does not suppress other spatial content.

## Non-obvious decisions and rationale

- Keep all posterior summaries rather than selecting `mean` or `q05` at
  inference time, so downstream analyses can choose an estimator explicitly.
- Keep Cell2location matrices in the AnnData table inside SpatialData rather
  than passing a parallel AnnData representation through the workflow.
- Train the reference once and broadcast signatures, avoiding repeated training
  for every slide.
- Separate inference and plotting because GPU tasks are expensive and plot
  development should remain resumable independently.
- Store a compact reconstruction histogram rather than full expected-count
  matrices to limit SpatialData size.
- Keep detailed posterior diagnostics outside MultiQC; MultiQC is limited to
  cross-sample summaries and one compact image per sample.
- Add high-resolution microscopy and segmentation only after deconvolution, so
  expensive inference does not stage or copy the largest image layers.

## Implementation and verification status

### Implemented

- input-driven optional reference branch;
- single reference training and signature reuse;
- QC-aware per-sample deconvolution;
- posterior and provenance storage in SpatialData;
- diagnostic and biological plots;
- MultiQC TSV and compact-panel integration;
- propagation into the final SpatialData and optional segmentation augmentation.

### Checked statically

- Python syntax for the Cell2location scripts;
- Nextflow lint for the modified workflow and modules;
- YAML parsing for module metadata and MultiQC configuration;
- isolated generation of plotting and assembled spatial-report outputs.

### Functionally validated

The user has reported successful reference training, spatial inference and plot
generation on representative and mock inputs. Full regression testing and the
final MultiQC presentation remain separate validation work.

## Open questions

- Multiple reference H5ADs are not implemented. The unresolved contract is
  whether each sample/reference pair should produce its own final SpatialData,
  remain a separate deconvolution output, or be merged with namespaced keys.
- The 35,000-location training batch is still being benchmarked to determine
  whether facility GPUs can support a higher default.
- Figure pagination and MultiQC readability should be reviewed for references
  containing many cell types.
