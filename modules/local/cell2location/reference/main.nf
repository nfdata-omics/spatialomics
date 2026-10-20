process CELL2LOCATION_REFERENCE {
    tag "$meta.id"
    label 'process_single', 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/cell2location-cuda:v0.1.5-torch2.8.0-cuda12.9'

    input:
    tuple val(meta), path(reference_h5ad)
    val cell_type_key
    val batch_key
    val counts_layer
    val gene_key

    output:
    tuple val(meta), path("*_reference_model"),          emit: model
    tuple val(meta), path("*_reference_signatures.csv"), emit: signatures
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/cell2location-cuda:v0.1.5-torch2.8.0-cuda12.9'), emit: versions_container, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}
    export XDG_CACHE_HOME=\${TMPDIR:-/tmp}/.cache
    mkdir -p "\$XDG_CACHE_HOME"

    cat << 'END_SCRIPT' > train_reference.py
${file("${moduleDir}/train_reference.py").text}
END_SCRIPT

    python3 train_reference.py \\
        ${args} \\
        --reference-h5ad "${reference_h5ad}" \\
        --cell-type-key "${cell_type_key}" \\
        --batch-key "${batch_key}" \\
        --counts-layer "${counts_layer}" \\
        --gene-key "${gene_key}" \\
        --output-model "${prefix}_reference_model" \\
        --output-signatures "${prefix}_reference_signatures.csv"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_reference_model"
    touch "${prefix}_reference_signatures.csv"
    """
}
