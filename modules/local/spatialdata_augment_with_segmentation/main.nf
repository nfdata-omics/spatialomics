process SPATIALDATA_AUGMENT_WITH_SEGMENTATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(zarr_folder), path(microscopy_image), path(segmentation_mask, arity: '0..*')
    val zarr_downsample_factor
    val microscopy_preview_max_dimension
    val chunk_size

    output:
    tuple val(meta), path("*_segmentation_microscopy.zarr"), emit: zarr
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'), emit: versions_container, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def mask_arg = segmentation_mask ? "--segmentation-mask-tif \"${segmentation_mask.first()}\"" : ''
    def overwrite_arg = task.ext.overwrite ? '--overwrite' : ''
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}

    cat << 'END_SCRIPT' > spatialdata_augment_with_segmentation.py
${file("${moduleDir}/spatialdata_augment_with_segmentation.py").text}
END_SCRIPT

    python3 spatialdata_augment_with_segmentation.py \\
        ${args} \\
        --input-zarr "${zarr_folder}" \\
        --output-zarr "${prefix}_segmentation_microscopy.zarr" \\
        --sample-name "${meta.id}" \\
        --microscopy-tif "${microscopy_image}" \\
        --zarr-downsample-factor "${zarr_downsample_factor}" \\
        --microscopy-preview-max-dimension "${microscopy_preview_max_dimension}" \\
        --chunk-size "${chunk_size}" \\
        ${mask_arg} ${overwrite_arg}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_segmentation_microscopy.zarr"
    """
}
