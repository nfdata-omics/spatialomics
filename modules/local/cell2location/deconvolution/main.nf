process CELL2LOCATION_DECONVOLUTION {
    tag "$meta.id"
    label 'process_single', 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/cell2location-cuda:v0.1.5-torch2.8.0-cuda12.9'

    input:
    tuple val(meta), path(zarr_folder)
    path signatures

    output:
    tuple val(meta), path("*_cell2location.zarr"),          emit: zarr
    tuple val(meta), path("*_cell2location_abundance.csv"), emit: abundance
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

    cat << 'END_SCRIPT' > deconvolve.py
${file("${moduleDir}/deconvolve.py").text}
END_SCRIPT

    python3 deconvolve.py \\
        ${args} \\
        --input-zarr "${zarr_folder}" \\
        --signatures "${signatures}" \\
        --output-zarr "${prefix}_cell2location.zarr" \\
        --output-abundance "${prefix}_cell2location_abundance.csv"
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p "${prefix}_cell2location.zarr"
    touch "${prefix}_cell2location_abundance.csv"
    """
}
