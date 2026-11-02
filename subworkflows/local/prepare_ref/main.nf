//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GTF              } from '../../../modules/nf-core/gunzip'
include { GUNZIP as GUNZIP_GFF              } from '../../../modules/nf-core/gunzip'
include { UNTAR as UNTAR_SPACERANGER_REF    } from "../../../modules/nf-core/untar"
include { GFFREAD                           } from '../../../modules/nf-core/gffread'
include { SPACERANGER_MKGTF                 } from '../../../modules/nf-core/spaceranger/mkgtf'
include { SPACERANGER_MKREF                 } from '../../../modules/nf-core/spaceranger/mkref'

workflow PREPARE_REF {

    take:
    reference_trigger        // channel: emits once when at least one FASTQ sample requires a reference
    reference_spec           // map: resolved reference paths and reference name

    main:

    // Versions collector
    ch_versions = channel.empty()

    if (reference_spec.spaceranger_index) {
        ch_fasta = channel.empty()
        ch_gtf   = channel.empty()

        ch_spaceranger_index_input = reference_trigger.map {
            file(reference_spec.spaceranger_index, checkIfExists: true)
        }

        // Extract an archived Space Ranger index; otherwise use the directory directly.
        if (reference_spec.spaceranger_index.toLowerCase().endsWith('.tar.gz')) {
            UNTAR_SPACERANGER_REF(
                ch_spaceranger_index_input.map { index ->
                    [[ id: index.name.replaceAll(/\.(tar)(\.gz)?$/, '') ], index]
                }
            )
            ch_spaceranger_index = UNTAR_SPACERANGER_REF.out.untar.map { _meta, index -> index }
            ch_versions = ch_versions.mix(UNTAR_SPACERANGER_REF.out.versions)
        } else {
            ch_spaceranger_index = ch_spaceranger_index_input
        }
    } else {
        // Uncompress GTF or convert GFF to GTF.
        ch_gtf = channel.empty()
        if (reference_spec.gtf) {
            ch_gtf_input = reference_trigger.map {
                file(reference_spec.gtf, checkIfExists: true)
            }
            if (reference_spec.gtf.toLowerCase().endsWith('.gz')) {
                GUNZIP_GTF(ch_gtf_input.map { file -> [ [:], file ] })
                ch_gtf = GUNZIP_GTF.out.gunzip.map { _meta, file -> file }
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions_gunzip)
            } else {
                ch_gtf = ch_gtf_input
            }
        } else if (reference_spec.gff) {
            ch_gff_input = reference_trigger.map {
                file(reference_spec.gff, checkIfExists: true)
            }
            if (reference_spec.gff.toLowerCase().endsWith('.gz')) {
                GUNZIP_GFF(ch_gff_input.map { file -> [ [:], file ] })
                ch_gff = GUNZIP_GFF.out.gunzip
                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions_gunzip)
            } else {
                ch_gff = ch_gff_input.map { file -> [ [:], file ] }
            }
            GFFREAD(ch_gff, [])
            ch_gtf = GFFREAD.out.gtf.map { _meta, file -> file }
            ch_versions = ch_versions.mix(GFFREAD.out.versions_gffread)
        }

        // Uncompress FASTA when needed.
        ch_fasta = channel.empty()
        if (reference_spec.fasta) {
            ch_fasta_input = reference_trigger.map {
                file(reference_spec.fasta, checkIfExists: true)
            }
            if (reference_spec.fasta.toLowerCase().endsWith('.gz')) {
                GUNZIP_FASTA(ch_fasta_input.map { file -> [ [:], file ] })
                ch_fasta = GUNZIP_FASTA.out.gunzip.map { _meta, file -> file }
                ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions_gunzip)
            } else {
                ch_fasta = ch_fasta_input
            }
        }

        // Filter supported annotation biotypes and build the Space Ranger index.
        SPACERANGER_MKGTF(ch_gtf)
        ch_versions = ch_versions.mix(SPACERANGER_MKGTF.out.versions)

        ch_reference_name = reference_spec.reference_name
            ? reference_trigger.map { reference_spec.reference_name }
            : channel.empty()

        SPACERANGER_MKREF(
            ch_fasta,
            SPACERANGER_MKGTF.out.gtf,
            ch_reference_name
        )
        ch_versions = ch_versions.mix(SPACERANGER_MKREF.out.versions)
        ch_spaceranger_index = SPACERANGER_MKREF.out.reference
    }

    emit:
    fasta             = ch_fasta                  // channel: path(genome.fasta)
    gtf               = ch_gtf                    // channel: path(genome.gtf)
    spaceranger_index = ch_spaceranger_index      // channel: path(/path/to/spaceranger/index/)
    versions          = ch_versions               // channel: [ versions.yml ]

}
