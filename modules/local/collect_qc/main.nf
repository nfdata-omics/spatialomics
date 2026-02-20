process COLLECT_QC {
    tag "all samples"
    label 'process_low'

    container 'docker.io/nfdata/plotly:v6.5.2'

    input:
    path annotated_obs_files
    path metrics_files

    output:
    path "qc_distributions.html", emit: distributions
    path "qc_metrics.csv",  emit: metrics
    path "versions.yml",   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}

    cat << END_SCRIPT > collect_qc.py
${file("${moduleDir}/collect_qc.py").text}
END_SCRIPT

    python3 collect_qc.py \
        --qc-metrics ${metrics_files.join(" ")} \
        --annotated-obs ${annotated_obs_files.join(" ")}

    python3 collect_qc.py \
        --versions-dict "${task.process}" > versions.yml
    """

    stub:
    """

    cat << END_SCRIPT > collect_qc.py
${file("${moduleDir}/collect_qc.py").text}
END_SCRIPT

     touch qc_distributions.html
     touch qc_metrics.csv

    python3 collect_qc.py \
        --versions-dict "${task.process}" > versions.yml
    """
}
