process SC_PREPROCESS_NORMALIZE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container 'docker.io/nfdata/muon-sc_rnaseq:v1.0.3'

    input:
    tuple val(meta), path(h5ad, stageAs: 'input/*', arity: '1')

    output:
    tuple val(meta), path("${meta.id}.normalized.h5ad"), emit: h5ad

    script:
    def args = task.ext.args ?: ''
    """
    NUMBA_CACHE_DIR=\$(mktemp -d "\$PWD/.numba-cache.XXXXXX")
    export NUMBA_CACHE_DIR
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    sc_preprocess.py normalize \\
        ${args} \\
        --input "${h5ad}" \\
        --output "${meta.id}.normalized.h5ad"
    """

    stub:
    """
    touch "${meta.id}.normalized.h5ad"
    """
}
