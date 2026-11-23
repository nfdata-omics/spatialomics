process SC_PREPROCESS_PCA {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container 'docker.io/nfdata/muon-sc_rnaseq:v1.0.3'

    input:
    tuple val(meta), path(h5ad, stageAs: 'input/*', arity: '1')
    val color

    output:
    tuple val(meta), path("${meta.id}.pca.h5ad"), emit: h5ad
    tuple val(meta), path("${meta.id}.pca.pdf"), emit: plot
    tuple val(meta), path("${meta.id}.pca_variance.pdf"), emit: variance_plot
    tuple val(meta), path("${meta.id}.pca_scores.txt"), emit: scores
    tuple val(meta), path("${meta.id}.pca_explained_variance.txt"), emit: explained_variance

    script:
    def args = task.ext.args ?: ''
    def quoted_color = color.toString().replace("'", "'\"'\"'")
    """
    export NUMBA_CACHE_DIR=/tmp
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    sc_preprocess.py pca \\
        ${args} \\
        --input "${h5ad}" \\
        --output "${meta.id}.pca.h5ad" \\
        --color='${quoted_color}' \\
        --plot "${meta.id}.pca.pdf" \\
        --variance-plot "${meta.id}.pca_variance.pdf" \\
        --scores "${meta.id}.pca_scores.txt" \\
        --explained-variance "${meta.id}.pca_explained_variance.txt"
    """

    stub:
    """
    touch "${meta.id}.pca.h5ad" "${meta.id}.pca.pdf" "${meta.id}.pca_variance.pdf" "${meta.id}.pca_scores.txt" "${meta.id}.pca_explained_variance.txt"
    """
}
