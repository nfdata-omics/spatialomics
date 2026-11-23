#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nfdata-omics/spatialomics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nfdata-omics/spatialomics
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SPATIALOMICS  } from './workflows/spatialomics'
include { PREPARE_REF   } from './subworkflows/local/prepare_ref'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_spatialomics_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_spatialomics_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_spatialomics_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFDATAOMICS_SPATIALOMICS {

    take:
    samplesheet // channel: samplesheet read in from --input
    spaceranger_outs // channel: spaceranger output paths read in from --input
    roi_annotations // channel: optional per-sample ROI GeoJSON files

    main:

    // Explicit parameters take precedence over entries in the selected genome dictionary.
    def reference_spec = [
        fasta: params.fasta ?: getGenomeAttribute('fasta'),
        gtf: params.gtf ?: getGenomeAttribute('gtf'),
        gff: params.gff ?: getGenomeAttribute('gff'),
        spaceranger_index: params.spaceranger_index ?: getGenomeAttribute('spaceranger_index'),
        reference_name: params.reference_name
    ]

    // Validate references only for samples that require Space Ranger execution.
    ch_validated_samplesheet = samplesheet.map { meta, fastqs ->
        assert reference_spec.spaceranger_index || (reference_spec.fasta && (reference_spec.gtf || reference_spec.gff)):
            "FASTQ input requires --spaceranger_index or --fasta with --gtf/--gff"

        assert reference_spec.spaceranger_index || reference_spec.reference_name:
            "Building a Space Ranger reference requires --reference_name"

        return [ meta, fastqs ]
    }

    // Gate reference preparation on the presence of a validated FASTQ sample.
    ch_reference_trigger = ch_validated_samplesheet.map { true }.unique()

    // Resolve and prepare the Space Ranger reference required by FASTQ samples.
    PREPARE_REF(
        ch_reference_trigger,
        reference_spec
    )

    ch_probeset           = params.probeset          ? channel.value(file(params.probeset, checkIfExists: true))          : channel.value([])

    // Define channel for the optional single-cell reference dataset
    ch_single_cell_reference = params.single_cell_reference_h5ad ? channel.value(file(params.single_cell_reference_h5ad, checkIfExists: true)) : channel.empty()

    //
    // WORKFLOW: Run pipeline
    //
    SPATIALOMICS (
        ch_validated_samplesheet,
        spaceranger_outs,
        roi_annotations,
        PREPARE_REF.out.spaceranger_index.collect(),
        PREPARE_REF.out.versions,
        ch_probeset,
        ch_single_cell_reference,
        params.multiqc_config,
        params.multiqc_logo,
        params.multiqc_methods_description,
        params.outdir,
    )
    emit:
    multiqc_report    = SPATIALOMICS.out.multiqc_report    // channel: /path/to/multiqc_report.html
    final_spatialdata = SPATIALOMICS.out.final_spatialdata // channel: [ val(meta), path(zarr) ]
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.help,
        params.help_full,
        params.show_hidden
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFDATAOMICS_SPATIALOMICS (
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.spaceranger_outs,
        PIPELINE_INITIALISATION.out.roi_annotations
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        NFDATAOMICS_SPATIALOMICS.out.multiqc_report
    )

    publish:
    final_spatialdata = NFDATAOMICS_SPATIALOMICS.out.final_spatialdata
}

output {
    final_spatialdata {
        path { meta, zarr ->
            zarr >> "spatialdata/${meta.id}_spatialdata.zarr"
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
