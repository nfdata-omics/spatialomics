process SC_PREPROCESS_QC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/../environment.yml"
    container 'docker.io/nfdata/muon-sc_rnaseq:v1.0.3'

    input:
    tuple val(meta), path(h5ad, stageAs: 'input/*', arity: '1')
    val min_counts

    output:
    tuple val(meta), path("${meta.id}.qc.h5ad"), emit: h5ad
    tuple val(meta), path("${meta.id}.qc_gex.pdf"), emit: gex_pdf
    tuple val(meta), path("${meta.id}.qc.html"), emit: html
    tuple val(meta), path("${meta.id}_qc_outcomes_mqc.html"), emit: multiqc_html
    tuple val(meta), path("${meta.id}_qc_summary_mqc.json"), emit: multiqc_table

    script:
    """
    NUMBA_CACHE_DIR=\$(mktemp -d "\$PWD/.numba-cache.XXXXXX")
    export NUMBA_CACHE_DIR
    export MPLCONFIGDIR=/tmp
    export XDG_CONFIG_HOME=/tmp

    sc_preprocess.py qc \\
        --input "${h5ad}" \\
        --output "${meta.id}.qc.h5ad" \\
        --data-type gex \\
        --min-counts ${min_counts.toString().toBigInteger()} \\
        --report-observation-unit rois \\
        --gex-pdf "${meta.id}.qc_gex.pdf" \\
        --html "${meta.id}.qc.html" \\
        --multiqc-html "${meta.id}_qc_outcomes_mqc.html" \\
        --multiqc-table "${meta.id}_qc_summary_mqc.json"
    """

    stub:
    """
    touch "${meta.id}.qc.h5ad" "${meta.id}.qc_gex.pdf" "${meta.id}.qc.html" "${meta.id}_qc_outcomes_mqc.html" "${meta.id}_qc_summary_mqc.json"
    """
}
