process CZI_TO_SPACERANGER_TIFF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/bftools:8.0.0--hdfd78af_0'

    input:
    tuple val(meta), path(input_image)

    output:
    tuple val(meta), path("*_spaceranger.tif"), emit: tiff
    tuple val(meta), path("*_spaceranger_info.txt"), emit: info
    tuple val("${task.process}"), val('container'), val('quay.io/biocontainers/bftools:8.0.0--hdfd78af_0'), emit: versions_container, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def maxMemoryMb = task.memory ? (task.memory.toMega() * 0.75) as int : 4096
    """
    export BF_MAX_MEM=${maxMemoryMb}m

    # Keep Bio-Formats' default flattened-resolution mode: series 0 is then
    # the full-resolution image only, without the lower-resolution pyramid.
    bfconvert \
        -no-upgrade \
        -series 0 \
        -bigtiff \
        -tilex 1024 \
        -tiley 1024 \
        "${input_image}" \
        "${prefix}_spaceranger.tif"

    showinf \
        -no-upgrade \
        -nopix \
        -nometa \
        "${prefix}_spaceranger.tif" \
        > "${prefix}_spaceranger_info.txt"

    grep -q 'Series count = 1' "${prefix}_spaceranger_info.txt"
    grep -q 'Image count = 1' "${prefix}_spaceranger_info.txt"
    grep -q 'RGB = true (3)' "${prefix}_spaceranger_info.txt"
    grep -q 'SizeZ = 1' "${prefix}_spaceranger_info.txt"
    grep -q 'SizeT = 1' "${prefix}_spaceranger_info.txt"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_spaceranger.tif"
    touch "${prefix}_spaceranger_info.txt"
    """
}
