process SC_PREPROCESS_HVG {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container 'docker.io/nfdata/muon-sc_rnaseq:v1.0.3'

    input:
    tuple val(meta), path(h5ad, stageAs: 'input/*', arity: '1')
    val n_top_genes

    output:
    tuple val(meta), path("${meta.id}.hvg.h5ad"), emit: h5ad

    script:
    def args = task.ext.args ?: ''
    def top_genes_arg = n_top_genes == '' ? '' : "--n-top-genes ${n_top_genes.toString().toBigInteger()}"
    """
    NUMBA_CACHE_DIR=\$(mktemp -d "\$PWD/.numba-cache.XXXXXX")
    export NUMBA_CACHE_DIR
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    sc_preprocess.py hvg \\
        ${args} \\
        ${top_genes_arg} \\
        --input "${h5ad}" \\
        --output "${meta.id}.hvg.h5ad"
    """

    stub:
    """
    touch "${meta.id}.hvg.h5ad"
    """
}
