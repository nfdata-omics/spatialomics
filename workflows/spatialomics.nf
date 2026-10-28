/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { SPACERANGER_COUNT           } from '../modules/nf-core/spaceranger/count/main'
include { COLLECT_SPACERANGER_METRICS } from '../modules/local/collect_spaceranger_metrics/main'
include { SPACERANGER_TO_ZARR         } from '../modules/local/spaceranger_to_zarr/main'
include { TAR as TAR_RAW               } from '../modules/nf-core/tar/main'
include { SPATIAL_QUALITY_CONTROL     } from '../modules/local/spatial_quality_control/main'
include { COLLECT_QC                  } from '../modules/local/collect_qc/main'
include { CELL2LOCATION_REFERENCE     } from '../modules/local/cell2location/reference/main'
include { CELL2LOCATION_DECONVOLUTION } from '../modules/local/cell2location/deconvolution/main'
include { CELL2LOCATION_PLOTS         } from '../modules/local/cell2location/plots/main'
include { CZI_TO_SPACERANGER_TIFF     } from '../modules/local/czi_to_spaceranger_tiff/main'
include { PREPARE_SEGMENTATION_IMAGE  } from '../modules/local/prepare_segmentation_image/main'
include { EXTRACT_VISIUM_REGION       } from '../modules/local/extract_visium_region/main'
include { CELLPOSE_SEGMENTATION       } from '../modules/local/cellpose_segmentation/main'
include { SKIMAGE_WATERSHED_SEGMENTATION } from '../modules/local/skimage_watershed_segmentation/main'
include { SPATIALDATA_AUGMENT_WITH_SEGMENTATION } from '../modules/local/spatialdata_augment_with_segmentation/main'
include { ASSEMBLE_SPATIAL_MULTIQC    } from '../modules/local/assemble_spatial_multiqc/main'
include { BIN2CELL                    } from '../modules/local/bin2cell/main'

include { PREPARE_REF            } from '../subworkflows/local/prepare_ref'
include { PREPARE_FASTQ          } from '../subworkflows/local/prepare_fastq'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_spatialomics_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SPATIALOMICS {

    take:
    ch_samplesheet        // channel: samplesheet read in from --input
    ch_spaceranger_outs   // channel: spaceranger output paths read in from --input
    ch_fasta              // value channel: path(fasta)
    ch_gtf                // value channel: path(gtf)
    ch_gff                // value channel: path(gff)
    ch_spaceranger_index  // value channel: path(spaceranger_index)
    ch_probeset           // value channel: path(probeset) (optional)
    ch_single_cell_reference // value channel: optional single-cell reference AnnData
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def ch_deconvolution_multiqc_plots = channel.empty()

    //
    // SUBWORKFLOW: Prepare reference genome
    //
    ch_prepared_spaceranger_index = channel.empty()
    def reference_requested = params.spaceranger_index || params.fasta || params.gtf || params.gff

    if (reference_requested) {
        PREPARE_REF (
            ch_fasta,
            ch_gtf,
            ch_gff,
            ch_spaceranger_index,
            params.reference_name
        )
        ch_prepared_spaceranger_index = PREPARE_REF.out.spaceranger_index
        ch_versions = ch_versions.mix(PREPARE_REF.out.versions)
    }

    //
    // SUBWORKFLOW: Prepare FastQ files
    //
    PREPARE_FASTQ (
        ch_samplesheet
    )
    ch_reads = PREPARE_FASTQ.out.reads
    ch_versions = ch_versions.mix(PREPARE_FASTQ.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_FASTQ.out.multiqc_files)

    // Space Ranger does not accept CZI directly. Convert CZI images once and reuse the
    // resulting full-resolution RGB BigTIFF throughout the rest of the workflow.
    ch_images_to_prepare = ch_reads.map { meta, _fastq -> [["id": meta.id], meta.image] }
    // Also prepare images supplied with pre-computed Space Ranger outputs for
    // capture-area extraction and any other image-dependent downstream processing.
    ch_images_to_prepare = ch_images_to_prepare.mix(
        ch_spaceranger_outs.map { meta, _out -> [["id": meta.id], meta.image] }
    )

    // Select only CZI images for conversion to Space Ranger-compatible TIFF, and pass through any other image formats unchanged
    ch_images_to_prepare
        .branch { _meta, image ->
            czi: image && image.toString().toLowerCase().endsWith('.czi')
            passthrough: true
        }
        .set { ch_images_by_format }

    CZI_TO_SPACERANGER_TIFF(
        ch_images_by_format.czi
    )

    CZI_TO_SPACERANGER_TIFF.out.tiff
        .mix(ch_images_by_format.passthrough)
        .set { ch_prepared_microscopy_images }

    ch_reads
        .map { meta, fastq -> [["id": meta.id], meta, fastq] }
        .join(ch_prepared_microscopy_images)
        .multiMap { sample_meta, meta, fastq, image ->
            reads: [ sample_meta, fastq ]
            slide_and_img: [ sample_meta, image, meta.slide, meta.area, meta.cytaimage, meta.darkimage, meta.colorizedimage, meta.alignment, meta.slidefile ]
        }
        .set { ch_reads_with_meta }

    //
    // MODULE: Align and quantify with Space Ranger
    //
    SPACERANGER_COUNT(
        ch_reads_with_meta.reads,
        ch_reads_with_meta.slide_and_img,
        ch_prepared_spaceranger_index,
        ch_probeset
    )
    ch_versions = ch_versions.mix(SPACERANGER_COUNT.out.versions)

    // Collect Space Ranger output paths for downstream processing
    SPACERANGER_COUNT.out.outs
        .mix(ch_spaceranger_outs) // Add any additional Space Ranger output paths provided via --input
        .set { ch_all_spaceranger_outs }

    //
    // MODULE: Collect Space Ranger metrics across samples
    //
    COLLECT_SPACERANGER_METRICS (
        ch_all_spaceranger_outs
            .collect{ _meta, folder -> folder }
    )
    ch_versions = ch_versions.mix(COLLECT_SPACERANGER_METRICS.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(COLLECT_SPACERANGER_METRICS.out.metrics)

    //
    // MODULE: Convert Space Ranger output to Zarr and compress it
    //
    SPACERANGER_TO_ZARR (
        ch_all_spaceranger_outs,
        "True"
    )
    ch_versions = ch_versions.mix(SPACERANGER_TO_ZARR.out.versions.first())

    TAR_RAW (
        SPACERANGER_TO_ZARR.out.zarr,
        '.gz'
    )
    ch_versions = ch_versions.mix(TAR_RAW.out.versions.first())

    //
    // MODULE: Spatial quality control
    //
    SPATIAL_QUALITY_CONTROL (
        SPACERANGER_TO_ZARR.out.zarr,
        params.bin_analysis_bin_size,
        params.bin_analysis_min_counts,
        params.bin_analysis_min_genes,
        params.bin_analysis_max_mt
    )
    ch_versions = ch_versions.mix(SPATIAL_QUALITY_CONTROL.out.versions.first())

    COLLECT_QC (
        SPATIAL_QUALITY_CONTROL.out.annotated_obs.collect{ _meta, path -> path },
        SPATIAL_QUALITY_CONTROL.out.metrics.collect{ _meta, path -> path },
        params.bin_analysis_min_counts,
        params.bin_analysis_min_genes,
        params.bin_analysis_max_mt
    )
    ch_versions = ch_versions.mix(COLLECT_QC.out.versions.first())

    ch_multiqc_files = ch_multiqc_files.mix(COLLECT_QC.out.metrics)
    ch_multiqc_files = ch_multiqc_files.mix(COLLECT_QC.out.distributions)

    // Use the QC-annotated SpatialData as the canonical input for downstream
    // transformations. If a single-cell reference is provided, replace it with
    // the deconvolved SpatialData; otherwise keep the QC output unchanged.
    if ( params.single_cell_reference_h5ad ) {
        //
        // MODULE: Train the Cell2location reference model once per pipeline run
        //
        CELL2LOCATION_REFERENCE (
            ch_single_cell_reference.map { h5ad -> [["id": h5ad.simpleName], h5ad] },
            params.single_cell_reference_cell_type_key,
            params.single_cell_reference_batch_key,
            params.single_cell_reference_counts_layer,
            params.single_cell_reference_gene_key
        )

        // Convert the single signatures output to a value channel so that the same
        // reference can be reused by every spatial sample.
        CELL2LOCATION_REFERENCE.out.signatures
            .map { _meta, signatures -> signatures }
            .first()
            .set { ch_cell2location_signatures }

        //
        // MODULE: Deconvolve each QC-annotated SpatialData sample
        //
        CELL2LOCATION_DECONVOLUTION (
            SPATIAL_QUALITY_CONTROL.out.zarr,
            ch_cell2location_signatures
        )

        // Generate complete per-sample result plots and MultiQC inputs.
        CELL2LOCATION_PLOTS (
            CELL2LOCATION_DECONVOLUTION.out.zarr,
            params.bin_analysis_bin_size,
            params.cell2location_abundance_statistic,
            params.cell2location_uncertainty_min_abundance,
            params.cell2location_multiqc_cell_types ?: ''
        )

        ch_deconvolution_multiqc_plots = CELL2LOCATION_PLOTS.out.mqc_plot
        ch_multiqc_files = ch_multiqc_files.mix(
            CELL2LOCATION_PLOTS.out.composition_mqc.map { _meta, path -> path }
        )
        ch_multiqc_files = ch_multiqc_files.mix(
            CELL2LOCATION_PLOTS.out.summary_mqc.map { _meta, path -> path }
        )

        ch_zarr_before_segmentation = CELL2LOCATION_DECONVOLUTION.out.zarr
    } else {
        ch_zarr_before_segmentation = SPATIAL_QUALITY_CONTROL.out.zarr
    }

    //
    // MODULE: Extract the full-resolution Visium capture area from the RGB microscopy image
    //

    ch_microscopy_images = ch_prepared_microscopy_images
        .filter { _meta, image -> image }

    SPACERANGER_TO_ZARR.out.zarr
        .map { meta, zarr -> [["id": meta.id], zarr] }
        .join(ch_microscopy_images)
        .set { ch_extract_visium_region_inputs }

    EXTRACT_VISIUM_REGION (
        ch_extract_visium_region_inputs,
        params.zarr_downsample_factor,
        params.capture_area_image_padding,
        params.capture_area_image_downsample_factor
    )
    ch_versions = ch_versions.mix(EXTRACT_VISIUM_REGION.out.versions.first())

    // Prepare segmentation-specific inputs unless segmentation is being skipped.
    if ( params.skip_segmentation ) {
        ch_segmentation_images = channel.empty()
        ch_crop_areas = channel.empty()
    } else {
        // A missing image is represented by an empty path value inside an emitted tuple,
        // not by an empty channel. Remove those tuples so image-dependent processes are
        // skipped idiomatically for samples without a microscopy image.
        ch_segmentation_images = ch_microscopy_images
        ch_reads.map { meta, _fastq -> [["id": meta.id], meta.crop_areas ] }
            .mix ( ch_spaceranger_outs.map { meta, _out -> [["id": meta.id], meta.crop_areas ] } )
            .set { ch_crop_areas }
    }

    // Select one segmentation backend while exposing the same mask and source-image
    // channels to all downstream modules.
    if ( params.segmentation_method == 'cellpose' ) {

        // Prepare grayscale memmappable OME-TIFF images for Cellpose and downstream use.
        PREPARE_SEGMENTATION_IMAGE(
            ch_segmentation_images
        )

        PREPARE_SEGMENTATION_IMAGE.out.tiff
            .join(EXTRACT_VISIUM_REGION.out.bounds)
            .set { ch_cellpose_inputs }

        CELLPOSE_SEGMENTATION (
            ch_cellpose_inputs
        )
        ch_versions = ch_versions.mix(CELLPOSE_SEGMENTATION.out.versions.first())
        ch_segmentation_masks = CELLPOSE_SEGMENTATION.out.mask
        ch_segmentation_source_images = PREPARE_SEGMENTATION_IMAGE.out.tiff

    } else if ( params.segmentation_method == 'skimage_watershed' ) {

        // Segment the RGB source directly using H&E deconvolution and watershed.
        // The module also emits the grayscale memmappable TIFF required downstream.
        ch_segmentation_images
            .join(EXTRACT_VISIUM_REGION.out.bounds)
            .set { ch_skimage_watershed_inputs }

        SKIMAGE_WATERSHED_SEGMENTATION (
            ch_skimage_watershed_inputs
        )
        ch_segmentation_masks = SKIMAGE_WATERSHED_SEGMENTATION.out.mask
        ch_segmentation_source_images = SKIMAGE_WATERSHED_SEGMENTATION.out.tiff

    } else {
        error "Unsupported segmentation method '${params.segmentation_method}'. Use 'cellpose' or 'skimage_watershed'."
    }

    //
    // MODULE: Augment SpatialData with segmentation and microscopy layers
    //

    ch_zarr_before_segmentation
        .map { meta, zarr -> [["id": meta.id], zarr] }
        .join(ch_segmentation_masks)
        .join(ch_segmentation_source_images)
        .join(ch_crop_areas)
        .set { ch_spatialdata_augment_with_segmentation_inputs }

    SPATIALDATA_AUGMENT_WITH_SEGMENTATION (
        ch_spatialdata_augment_with_segmentation_inputs,
        params.zarr_downsample_factor,
        params.microscopy_preview_max_dimension,
        2048
    )

    // Select the most advanced SpatialData available for every sample. Samples
    // without segmentation retain the QC or Cell2location Zarr selected above.
    ch_zarr_before_segmentation
        .map { meta, zarr -> [meta.id, meta, zarr] }
        .join(
            SPATIALDATA_AUGMENT_WITH_SEGMENTATION.out.zarr
                .map { meta, zarr -> [meta.id, zarr] },
            remainder: true
        )
        .map { _id, meta, input_zarr, augmented_zarr = null ->
            [meta, augmented_zarr ?: input_zarr]
        }
        .set { ch_final_zarr }

    //
    // MODULE: Assemble the available per-sample spatial plots for MultiQC
    //

    SPATIAL_QUALITY_CONTROL.out.mqc_plot
        .mix(SPATIALDATA_AUGMENT_WITH_SEGMENTATION.out.registration_plot)
        .mix(SPATIALDATA_AUGMENT_WITH_SEGMENTATION.out.crop_areas_plot)
        .mix(SPATIALDATA_AUGMENT_WITH_SEGMENTATION.out.segmentation_crop_panels_plot)
        .mix(ch_deconvolution_multiqc_plots)
        .map { meta, plot -> [meta.id, meta, plot] }
        .groupTuple(by: 0)
        .map { _id, metadata, plots -> [metadata.first(), plots] }
        .set { ch_spatial_multiqc_inputs }

    ASSEMBLE_SPATIAL_MULTIQC (
        ch_spatial_multiqc_inputs
    )
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLE_SPATIAL_MULTIQC.out.html.collect{ _meta, path -> path })
    ch_multiqc_files = ch_multiqc_files.mix(SPATIALDATA_AUGMENT_WITH_SEGMENTATION.out.segmentation_stats.collect{ _meta, paths -> paths })

    //
    // MODULE: Aggregate Visium HD bins into cell-level AnnData with Bin2Cell
    //

    if ( params.skip_bin2cell ) {
        ch_bin2cell_inputs = channel.empty()
    } else {
        ch_all_spaceranger_outs
            .map { meta, outs -> [["id": meta.id], outs] }
            .join(ch_segmentation_masks)
            .join(ch_segmentation_source_images)
            .set { ch_bin2cell_inputs }
    }

    BIN2CELL (
        ch_bin2cell_inputs,
        params.bin2cell_bin_size,
        params.bin2cell_volume_ratio
    )
    ch_versions = ch_versions.mix(BIN2CELL.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(BIN2CELL.out.summary.collect{ _meta, path -> path })

    //
    // Collate and save software versions
    //

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'spatialomics_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'spatialomics'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:
    multiqc_report    = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    final_spatialdata = ch_final_zarr                                                // channel: [ val(meta), path(zarr) ]
    versions          = ch_versions                                                   // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
