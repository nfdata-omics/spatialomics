process SKIMAGE_WATERSHED_SEGMENTATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/scikit-image:0.24.0'

    input:
    tuple val(meta), path(input_image), path(visium_bounds)

    output:
    tuple val(meta), path("*_mask.tif"), emit: mask
    tuple val("${task.process}"), val('container'), val('quay.io/biocontainers/scikit-image:0.24.0'), emit: versions_container, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat << 'END_SCRIPT' > skimage_watershed_segmentation.py
${file("${moduleDir}/skimage_watershed_segmentation.py").text}
END_SCRIPT

    python3 skimage_watershed_segmentation.py \
        ${args} \
        --input-image "${input_image}" \
        --visium-bounds "${visium_bounds}" \
        --output-mask "${prefix}_mask.tif"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_mask.tif"
    """
}
