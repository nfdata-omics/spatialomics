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
include { CELLPOSE_SEGMENTATION       } from '../modules/local/cellpose_segmentation/main'

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

    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

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
    ch_multiqc_files = ch_multiqc_files.mix(SPATIAL_QUALITY_CONTROL.out.mqc_plot.collect{ _meta, path -> path })

    //
    // MODULE: Cell segmentation with Cellpose
    //
    ch_reads.map { meta, _fastq -> [meta, meta.image] }
        .mix { ch_spaceranger_outs.map { meta, _out -> [meta, meta.image] } }
        .set { ch_cellpose_input }

    CELLPOSE_SEGMENTATION (
        ch_cellpose_input
    )
    ch_versions = ch_versions.mix(CELLPOSE_SEGMENTATION.out.versions.first())

    //
    // Collate and save software versions
    //
    def topic_versions = Channel.topic("versions")
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

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'spatialomics_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
