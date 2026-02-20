process SPATIAL_QUALITY_CONTROL {
    tag "$meta.id"
    label 'process_low'

    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(zarr_folder)
    path script_file

    output:
    path "versions.yml",                             emit: versions
    tuple val(meta), path("*_qc.h5ad"),              emit: anndata
    tuple val(meta), path("*_qc_annotated_obs.csv"), emit: annotated_obs
    tuple val(meta), path("*_qc_metrics.csv"),       emit: metrics
    tuple val(meta), path("*_qc_distributions.png"), emit: distributions
    tuple val(meta), path("*_qc_mqc.png"),           emit: mqc_plot

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}

    python3 ${script_file} \
        --zarr "${zarr_folder}" \
        --sample "${prefix}"

    mv ${prefix}_qc_spatial_plots.png ${prefix}_qc_mqc.png

    python3 ${script_file} \
        --versions-dict "${task.process}" > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_qc.h5ad
    touch ${prefix}_qc_annotated_obs.csv
    touch ${prefix}_qc_metrics.csv
    touch ${prefix}_qc_distributions.png
    touch ${prefix}_qc_mqc.png

    python3 ${script_file} \
        --versions-dict "${task.process}" > versions.yml
    """
}
