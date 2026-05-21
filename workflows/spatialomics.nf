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
include { TAR                         } from '../modules/nf-core/tar/main'
include { SPATIAL_QUALITY_CONTROL     } from '../modules/local/spatial_quality_control/main'
include { COLLECT_QC                  } from '../modules/local/collect_qc/main'
include { IMAGE_TO_TIFF               } from '../modules/local/image_to_tiff/main'
include { VISIUM_BOUNDS               } from '../modules/local/visium_bounds/main'
include { CELLPOSE_SEGMENTATION       } from '../modules/local/cellpose_segmentation/main'
include { SEGMENTATION_AND_MICROSCOPY_PLOTS } from '../modules/local/segmentation_and_microscopy_plots/main'
include { ASSEMBLE_IMAGING_MULTIQC    } from '../modules/local/assemble_imaging_multiqc/main'
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
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // SUBWORKFLOW: Prepare reference genome
    //
    PREPARE_REF (
        ch_fasta,
        ch_gtf,
        ch_gff,
        ch_spaceranger_index,
        params.reference_name
    )
    ch_versions = ch_versions.mix(PREPARE_REF.out.versions)

    //
    // SUBWORKFLOW: Prepare FastQ files
    //
    PREPARE_FASTQ (
        ch_samplesheet
    )
    ch_reads = PREPARE_FASTQ.out.reads
    ch_versions = ch_versions.mix(PREPARE_FASTQ.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(PREPARE_FASTQ.out.multiqc_files)

    ch_reads
        .multiMap { meta, fastq ->
            reads: [ ["id": meta.id], fastq ]
            slide_and_img: [ ["id": meta.id], meta.image, meta.slide, meta.area, meta.cytaimage, meta.darkimage, meta.colorizedimage, meta.alignment, meta.slidefile ]
        }
        .set { ch_reads_with_meta }

    //
    // MODULE: Align and quantify with Space Ranger
    //
    SPACERANGER_COUNT(
        ch_reads_with_meta.reads,
        ch_reads_with_meta.slide_and_img,
        PREPARE_REF.out.spaceranger_index,
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

    TAR (
        SPACERANGER_TO_ZARR.out.zarr,
        '.gz'
    )
    ch_versions = ch_versions.mix(TAR.out.versions.first())

    //
    // MODULE: Spatial quality control
    //
    SPATIAL_QUALITY_CONTROL (
        SPACERANGER_TO_ZARR.out.zarr
    )
    ch_versions = ch_versions.mix(SPATIAL_QUALITY_CONTROL.out.versions.first())

    COLLECT_QC (
        SPATIAL_QUALITY_CONTROL.out.annotated_obs.collect{ _meta, path -> path },
        SPATIAL_QUALITY_CONTROL.out.metrics.collect{ _meta, path -> path }
    )
    ch_versions = ch_versions.mix(COLLECT_QC.out.versions.first())

    ch_multiqc_files = ch_multiqc_files.mix(COLLECT_QC.out.metrics)
    ch_multiqc_files = ch_multiqc_files.mix(COLLECT_QC.out.distributions)

    // Prepare microscopy images and crop areas for downstream processing, unless segmentation is being skipped
    if ( params.skip_segmentation ) {
        ch_microscopy_images = channel.empty()
        ch_crop_areas = channel.empty()
    } else {
        ch_reads.map { meta, _fastq -> [["id": meta.id], meta.image] }
            .mix ( ch_spaceranger_outs.map { meta, _out -> [["id": meta.id], meta.image] } )
            .set { ch_microscopy_images }
        ch_reads.map { meta, _fastq -> [["id": meta.id], meta.crop_areas ] }
            .mix ( ch_spaceranger_outs.map { meta, _out -> [["id": meta.id], meta.crop_areas ] } )
            .set { ch_crop_areas }
    }

    //
    // MODULE: Convert microscopy images to memmappable OME-TIFF format
    //
    IMAGE_TO_TIFF(
        ch_microscopy_images
    )

    SPACERANGER_TO_ZARR.out.zarr
        .map { meta, zarr -> [["id": meta.id], zarr] }
        .join(IMAGE_TO_TIFF.out.tiff)
        .set { ch_visium_bounds_inputs }

    //
    // MODULE: Compute full-resolution microscopy bounds for the Visium capture area
    //
    VISIUM_BOUNDS (
        ch_visium_bounds_inputs,
        params.zarr_downsample_factor
    )
    ch_versions = ch_versions.mix(VISIUM_BOUNDS.out.versions.first())

    IMAGE_TO_TIFF.out.tiff
        .join(VISIUM_BOUNDS.out.bounds)
        .set { ch_cellpose_inputs }

    //
    // MODULE: Cell segmentation with Cellpose
    //
    CELLPOSE_SEGMENTATION (
        ch_cellpose_inputs
    )
    ch_versions = ch_versions.mix(CELLPOSE_SEGMENTATION.out.versions.first())

    SPACERANGER_TO_ZARR.out.zarr
        .map { meta, zarr -> [["id": meta.id], zarr] }
        .join(CELLPOSE_SEGMENTATION.out.mask)
        .join(IMAGE_TO_TIFF.out.tiff)
        .join(ch_crop_areas)
        .set { ch_segmentation_and_microscopy_inputs }

    //
    // MODULE: Generate segmentation and microscopy plots
    //
    SEGMENTATION_AND_MICROSCOPY_PLOTS (
        ch_segmentation_and_microscopy_inputs,
        params.zarr_downsample_factor,
        16,
        2048
    )

    if ( params.skip_segmentation ) {

        SPATIAL_QUALITY_CONTROL.out.mqc_plot
            .map { meta, plot -> [meta, [plot]]}
            .set { ch_imaging_multiqc_inputs }

    } else {

        SPATIAL_QUALITY_CONTROL.out.mqc_plot
            .map { meta, png -> [["id": meta.id], png] }
            .join(SEGMENTATION_AND_MICROSCOPY_PLOTS.out.registration_plot)
            .join(SEGMENTATION_AND_MICROSCOPY_PLOTS.out.crop_areas_plot)
            .join(SEGMENTATION_AND_MICROSCOPY_PLOTS.out.segmentation_crop_panels_plot)
            .map { meta, mqc_plot, registration_plot, crop_areas_plot, segmentation_crop_panels_plot
                     -> [meta, [mqc_plot, registration_plot, crop_areas_plot, segmentation_crop_panels_plot]]}
            .set { ch_imaging_multiqc_inputs }

    }

    ASSEMBLE_IMAGING_MULTIQC (
        ch_imaging_multiqc_inputs
    )
    ch_multiqc_files = ch_multiqc_files.mix(ASSEMBLE_IMAGING_MULTIQC.out.html.collect{ _meta, path -> path })
    ch_multiqc_files = ch_multiqc_files.mix(SEGMENTATION_AND_MICROSCOPY_PLOTS.out.segmentation_stats.collect{ _meta, paths -> paths })

    if ( params.skip_bin2cell ) {
        ch_bin2cell_inputs = channel.empty()
    } else {
        ch_all_spaceranger_outs
            .map { meta, outs -> [["id": meta.id], outs] }
            .join(CELLPOSE_SEGMENTATION.out.mask)
            .join(IMAGE_TO_TIFF.out.tiff)
            .set { ch_bin2cell_inputs }
    }

    //
    // MODULE: Aggregate Visium HD bins into cell-level AnnData with Bin2Cell
    //
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
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
