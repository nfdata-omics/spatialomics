# nfdata-omics/spatialomics: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Raw SpatialData archive](#raw-spatialdata-archive) - SpatialData generated directly from Space Ranger output
- [Cell2location deconvolution](#cell2location-deconvolution) - Joint-model diagnostics, posterior tables, and per-sample spatial plots
- [Final SpatialData](#final-spatialdata) - Analysis-ready SpatialData selected from the most advanced completed workflow branch
- [FastQC](#fastqc) - Raw read QC
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Raw SpatialData archive

<details markdown="1">
<summary>Output files</summary>

- `count/`
  - `*.zarr.tar.gz`: Archived SpatialData generated directly from the Space Ranger output, before quality control and downstream analyses.

</details>

This archive is provided for completeness and long-term storage. Extract it before opening the contained Zarr directory with SpatialData. When Space Ranger provides a CytAssist image and complete registration matrices, the archive contains the registered CytAssist layer in the sample coordinate system. It also contains a compact capture-area geometry derived from the complete Space Ranger bin-position grid.

### Cell2location deconvolution

This section is produced when `--single_cell_reference_h5ad` is supplied. One spatial model is trained jointly across all samples in the pipeline run, using sample identity as the Cell2location batch.

<details markdown="1">
<summary>Output files</summary>

- `deconvolution/`
  - `*_reference_signatures.csv`: Cell-type expression signatures inferred from the single-cell reference.
  - `*_cell2location_abundance_mean.csv`: Per-location posterior mean abundance for each cell type.
  - `*_cell2location_abundance_q05.csv`: Per-location posterior q05 abundance for each cell type.
  - `cell2location_training_history.csv`: Joint-model training ELBO values by epoch.
  - `cell2location_training_history.png`: Complete and post-warmup training history.
  - `cell2location_joint_diagnostics.json`: Model settings, provenance and compact reconstruction data.
  - `cell2location_joint_diagnostics.tsv`: Selected joint-model diagnostic metrics.
  - `cell2location_posterior_manifest.tsv`: Sample-to-posterior association and input dimensions.
  - `cell2location_reconstruction_accuracy.png`: Joint expected-versus-observed reconstruction plot, when available.
  - `cell2location_reconstruction_accuracy_error.txt`: Diagnostic error recorded instead when the reconstruction plot cannot be generated.
  - `<sample>_cell2location/`
    - `*_posterior_mean_vs_q05.pdf`: Posterior mean and q05 comparison by cell type.
    - `*_relative_posterior_uncertainty_spatial.pdf`: Spatial posterior uncertainty maps.
    - `*_celltype_abundance_*.pdf`: Cell-type abundance maps.
    - `*_celltype_fraction_*.pdf`: Within-location cell-type fraction maps.
    - `*_dominant_celltype_*.pdf`: Dominant cell type and assignment-strength map.
    - `*_total_abundance_*.pdf`: Total-abundance distribution and spatial map.
    - `*_deconvolution_mqc.png`: Compact per-sample image included in MultiQC.
    - `*_diagnostics.tsv` and `*_mqc.tsv`: Per-sample diagnostics and MultiQC custom-content inputs.

</details>

The posterior matrices are also stored in the corresponding table of the final SpatialData. Locations excluded by spatial QC contain missing posterior values.

### Final SpatialData

<details markdown="1">
<summary>Output files</summary>

- `spatialdata/`
  - `*_spatialdata.zarr/`: Final SpatialData directory containing the most advanced result available for each sample.

</details>

The final output is published as an unpacked Zarr directory and can be opened directly:

```python
import spatialdata as sd

sdata = sd.read_zarr("results/spatialdata/SAMPLE_spatialdata.zarr")
```

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html`: FastQC report containing quality metrics.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data file and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your sequenced reads. It provides information about the quality score distribution across your reads, per base sequence content (%A/T/G/C), adapter contamination and overrepresented sequences. For further reading and documentation see the [FastQC help pages](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

The per-sample spatial section uses bounded preview images rather than archival-resolution figures. It includes H&E/CytAssist registration and Space Ranger tissue selection, three paired spatial QC rows, combined crop-location and segmentation review, optional ROI annotations, and Cell2location deconvolution when their corresponding inputs are available. Detailed data remain in the SpatialData and analysis-specific outputs.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
