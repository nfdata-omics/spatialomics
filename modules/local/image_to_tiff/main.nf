process IMAGE_TO_TIFF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/cellpose-cuda:v4.0.8-torch2.10.0-cuda12.1.1'

    input:
    tuple val(meta), path(input_image)

    output:
    tuple val(meta), path("*_memmappable.ome.tif"), emit: tiff
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/cellpose-cuda:v4.0.8-torch2.10.0-cuda12.1.1'), emit: versions_container, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat << 'END_SCRIPT' > convert_to_memmappable_tiff.py
${file("${moduleDir}/convert_to_memmappable_tiff.py").text}
END_SCRIPT

    python3 convert_to_memmappable_tiff.py \
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
