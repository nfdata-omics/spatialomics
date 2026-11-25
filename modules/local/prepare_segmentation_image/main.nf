process PREPARE_SEGMENTATION_IMAGE {
    tag "$meta.id"
    label 'process_single', 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(input_image)

    output:
    tuple val(meta), path("*_memmappable.ome.tif"), emit: tiff
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'), emit: versions_container, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat << 'END_SCRIPT' > prepare_segmentation_image.py
${file("${moduleDir}/prepare_segmentation_image.py").text}
END_SCRIPT

    python3 prepare_segmentation_image.py \
        "${input_image}" \
        --output "${prefix}_memmappable.ome.tif" \
        --factor 1 \
        --compression none
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_memmappable.ome.tif"
    """
}
