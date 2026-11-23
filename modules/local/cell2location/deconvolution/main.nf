process CELL2LOCATION_DECONVOLUTION {
    tag "${sample_ids.size()} sample(s)"
    label 'process_medium', 'process_gpu'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/cell2location-cuda:v0.1.5-torch2.8.0-cuda12.9'

    input:
    tuple val(sample_ids), path(zarr_folders, arity: '1..*')
    path signatures

    output:
    path "cell2location_posteriors",             emit: posteriors
    path "cell2location_posterior_manifest.tsv", emit: manifest
    path "cell2location_training_history.csv",   emit: training_history
    path "cell2location_joint_diagnostics.json", emit: diagnostics
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/cell2location-cuda:v0.1.5-torch2.8.0-cuda12.9'), emit: versions_container, topic: versions

    script:
    def args = task.ext.args ?: ''
    if (sample_ids.size() != zarr_folders.size()) {
        error "CELL2LOCATION_DECONVOLUTION received ${sample_ids.size()} sample IDs and ${zarr_folders.size()} SpatialData paths"
    }
    def sample_rows = (0..<sample_ids.size()).collect { index ->
        "${sample_ids[index]}\t${zarr_folders[index]}"
    }.join('\n')
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}
    export XDG_CACHE_HOME=\${TMPDIR:-/tmp}/.cache
    mkdir -p "\$XDG_CACHE_HOME"

    cat << 'END_SAMPLES' > samples.tsv
sample_id	spatialdata
${sample_rows}
END_SAMPLES

    cat << 'END_SCRIPT' > deconvolve.py
${file("${moduleDir}/deconvolve.py").text}
END_SCRIPT

    python3 deconvolve.py \
        ${args} \
        --samples samples.tsv \
        --signatures "${signatures}" \
        --output-dir cell2location_posteriors \
        --output-manifest cell2location_posterior_manifest.tsv \
        --output-training-history cell2location_training_history.csv \
        --output-diagnostics cell2location_joint_diagnostics.json
    """

    stub:
    def posterior_files = sample_ids.collect { sample_id ->
        "touch cell2location_posteriors/${sample_id}_cell2location_posterior.h5ad"
    }.join('\n')
    """
    mkdir -p cell2location_posteriors
    ${posterior_files}
    touch cell2location_posterior_manifest.tsv
    touch cell2location_training_history.csv
    touch cell2location_joint_diagnostics.json
    """
}
