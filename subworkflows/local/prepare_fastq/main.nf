
include { CAT_FASTQ } from '../../../modules/nf-core/cat/fastq/main'
include { FQ_LINT   } from '../../../modules/nf-core/fq/lint/main'
include { FASTQC    } from '../../../modules/nf-core/fastqc/main'

workflow PREPARE_FASTQ {

    take:
    ch_samplesheet             // channel: [ val(meta), [ reads ] ]

    main:

    ch_versions = Channel.empty()
    ch_reads    = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_lint_log = Channel.empty()

    // divide samplesheet into samples with single and multiple fastqs
    ch_samplesheet
        .branch { meta, fastqs ->
            single: fastqs.size() == 1
            return [meta, fastqs.flatten()]
            multiple: fastqs.size() > 1
            return [meta, fastqs.flatten()]
        }
        .set { ch_fastq }

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ(
        ch_fastq.multiple
    ).reads.mix(ch_fastq.single).set { ch_reads }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // MODULE: Lint FastQ files
    //
    FQ_LINT(
        ch_reads
    )
    ch_versions = ch_versions.mix(FQ_LINT.out.versions.first())
    ch_lint_log = ch_lint_log.mix(FQ_LINT.out.lint)

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())


    emit:
    reads           = ch_reads                // channel: [ val(meta), path(reads) ]
    lint_log        = ch_lint_log
    multiqc_files   = ch_multiqc_files
    versions        = ch_versions             // channel: [ versions.yml ]
}
