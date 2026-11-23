process IMPORT_ROI_ANNOTATIONS {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(zarr_folder), path(roi_annotations), path(bounds), path(capture_area_image)

    output:
    tuple val(meta), path("*_roi_annotations.zarr"), emit: zarr
    tuple val(meta), path("*_roi_annotations_mqc.png"), emit: roi_plot
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'),
        emit: versions_container, topic: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}
    export XDG_CACHE_HOME=\${TMPDIR:-/tmp}/.cache
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"

    cat << 'END_SCRIPT' > import_roi_annotations.py
${file("${moduleDir}/import_roi_annotations.py").text}
END_SCRIPT

    python3 import_roi_annotations.py \
        --input-zarr "${zarr_folder}" \
        --sample-name "${prefix}" \
        --sample-id "${meta.id}" \
        --roi-annotations "${roi_annotations}" \
        --bounds "${bounds}" \
        --capture-area-image "${capture_area_image}" \
        --output-zarr "${prefix}_roi_annotations.zarr" \
        --output-plot "${prefix}_roi_annotations_mqc.png"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_roi_annotations.zarr"
    touch "${prefix}_roi_annotations_mqc.png"
    """
}
