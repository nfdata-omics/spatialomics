process CELLPOSE_SEGMENTATION {
    tag "$meta.id"
    label 'process_medium', 'process_gpu', 'process_high_memory'

    container 'docker.io/nfdata/cellpose-cuda:v4.0.8-torch2.10.0-cuda12.1.1'

    input:
    tuple val(meta), path(input_image)

    output:
    tuple val(meta), path("*_mask.tif")
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export TORCHINDUCTOR_CACHE_DIR=\${TMPDIR:-/tmp}/torchinductor-cache
    export CELLPOSE_LOCAL_MODELS_PATH=\${TMPDIR:-/tmp}/models
    mkdir -p "\$TORCHINDUCTOR_CACHE_DIR"
    mkdir -p "\$CELLPOSE_LOCAL_MODELS_PATH"

    cat << 'END_SCRIPT' > cellpose_segmentation.py
${file("${moduleDir}/cellpose_segmentation.py").text}
END_SCRIPT

    python3 cellpose_segmentation.py \
    ${args} \
    --input-image "${input_image}" \
    --output-mask "${meta.id}_mask.tif"

    python3 cellpose_segmentation.py \
        --versions-dict "${task.process}" > versions.yml
    """

    stub:
    """
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export TORCHINDUCTOR_CACHE_DIR=\${TMPDIR:-/tmp}/torchinductor-cache
    export USER=\${USER:-nextflow}
    export LOGNAME=\${LOGNAME:-\${USER}}
    export HOME=\${HOME:-\${TMPDIR:-/tmp}}
    mkdir -p "\$TORCHINDUCTOR_CACHE_DIR"
    export CELLPOSE_LOCAL_MODELS_PATH=\${TMPDIR:-/tmp}/models

    cat << 'END_SCRIPT' > cellpose_segmentation.py
${file("${moduleDir}/cellpose_segmentation.py").text}
END_SCRIPT

    touch "${meta.id}_mask.tif"

    python3 cellpose_segmentation.py \
        --versions-dict "${task.process}" > versions.yml
    """
}
