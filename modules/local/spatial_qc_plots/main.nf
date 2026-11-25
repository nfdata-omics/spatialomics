process SPATIAL_QC_PLOTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(zarr_folder)
    val resolution

    output:
    path "versions.yml", emit: versions
    tuple val(meta), path("*_qc_distributions.png"), emit: distributions
    tuple val(meta), path("*_qc_mqc.webp"), emit: mqc_plot

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}
    export XDG_CACHE_HOME=\${TMPDIR:-/tmp}/.cache
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"

    cat << 'END_LAYOUT' > spatial_preview_layout.py
${file("${moduleDir}/spatial_preview_layout.py").text}
END_LAYOUT

    cat << 'END_SCRIPT' > spatial_qc_plots.py
${file("${moduleDir}/spatial_qc_plots.py").text}
END_SCRIPT

    python3 spatial_qc_plots.py \
        --zarr "${zarr_folder}" \
        --sample "${prefix}" \
        --resolution "${resolution}"

    python3 spatial_qc_plots.py \
        --versions-dict "${task.process}" > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    touch "${prefix}_qc_distributions.png"
    touch "${prefix}_qc_mqc.webp"
    cat << 'END_LAYOUT' > spatial_preview_layout.py
${file("${moduleDir}/spatial_preview_layout.py").text}
END_LAYOUT
    cat << 'END_SCRIPT' > spatial_qc_plots.py
${file("${moduleDir}/spatial_qc_plots.py").text}
END_SCRIPT
    python3 spatial_qc_plots.py --versions-dict "${task.process}" > versions.yml
    """
}
