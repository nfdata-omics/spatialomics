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
    fasta                    // file: /path/to/genome.fasta (optional!)
    gtf                      // file: /path/to/genome.gtf
    gff                      // file: /path/to/genome.gff
    spaceranger_index        // directory: /path/to/spaceranger/index/ (optional!)
    reference_name           // string: name for the new reference (if building new index)

    main:

    assert (params.spaceranger_index) || (params.fasta && (params.gtf || params.gff)):
        "Must provide a the spaceranger index (--spaceranger_index) \
        or a fasta file ('--fasta') and a gtf/gff file ('--gtf'/'--gff') if no index is given!"

    assert (params.spaceranger_index) || (reference_name):
        "Must provide a reference name (--reference_name) when building a new spaceranger index!"

    // Versions collector
    ch_versions = Channel.empty()

    if (params.spaceranger_index) {

        ch_fasta = Channel.empty()
        ch_gtf   = Channel.empty()

        // Define spaceranger index channel from the user-provided one
        if (params.spaceranger_index ==~ /.*\.tar\.gz$/) {
            UNTAR_SPACERANGER_REF ([
                ["id": file(params.spaceranger_index).name.replaceAll(/\.(tar)(\.gz)?$/, '')],
                spaceranger_index
            ])
            ch_spaceranger_index = UNTAR_SPACERANGER_REF.out.untar.map{ it[1] }
            ch_versions = ch_versions.mix(UNTAR_SPACERANGER_REF.out.versions)
        } else {
            ch_spaceranger_index = spaceranger_index
        }

    } else {

        // Uncompress GTF or GFF, and convert GFF to GTF
        ch_gtf = Channel.empty()
        if (params.gtf) {
            if (params.gtf.endsWith('.gz')) {
                GUNZIP_GTF( gtf.map { [ [:], it ] } )
                ch_gtf      = GUNZIP_GTF.out.gunzip.map { it[1] }
                ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
            } else {
                ch_gtf = gtf
            }
        } else if (params.gff) {
            if (params.gff.endsWith('.gz')) {
                GUNZIP_GFF( gff.map { [ [:], it ] } )
                ch_gff      = GUNZIP_GFF.out.gunzip
                ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
            } else {
                ch_gff = gff.map { [ [:], it ] }
            }
            ch_gtf      = GFFREAD(ch_gff, []).gtf.map { it[1] }
            ch_versions = ch_versions.mix(GFFREAD.out.versions)
        }

        // Uncompress FASTA if needed
        ch_fasta = Channel.of([])
        if (params.fasta.endsWith('.gz')) {
            GUNZIP_FASTA( fasta.map { [ [:], it ] } )
            ch_fasta    = GUNZIP_FASTA.out.gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
        } else {
            ch_fasta = fasta
        }

        //
        // Prepare gft file by keeping specific biotypes
        //
        SPACERANGER_MKGTF(
            ch_gtf,
        )
        ch_gtf_filtered = SPACERANGER_MKGTF.out.gtf

        //
        // Create Spacer Ranger reference
        //
        SPACERANGER_MKREF(
            ch_fasta,
            ch_gtf_filtered,
            reference_name
        )

        // Channel to handle SPACERANGER_MKREF output
        ch_spaceranger_index = SPACERANGER_MKREF.out.reference.ifEmpty {
            file("no_spaceranger", checkIfExists: false)
        }
    }

    emit:
    fasta             = ch_fasta                  // channel: path(genome.fasta)
    gtf               = ch_gtf                    // channel: path(genome.gtf)
    spaceranger_index = ch_spaceranger_index      // channel: path(/path/to/spaceranger/index/)
    versions          = ch_versions               // channel: [ versions.yml ]

}
