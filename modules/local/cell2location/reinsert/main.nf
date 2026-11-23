process CELL2LOCATION_REINSERT {
    tag "${sample_ids.size()} sample(s)"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(sample_ids), path(zarr_folders, arity: '1..*')
    path posterior_dir
    path posterior_manifest
    path training_history
    path diagnostics

    output:
    path "*_cell2location.zarr",               emit: zarr
    path "*_cell2location_abundance_mean.csv", emit: abundance_mean
    path "*_cell2location_abundance_q05.csv",  emit: abundance_q05
    path "cell2location_training_history.png", emit: training_history_plot
    path "cell2location_reconstruction_accuracy.png",       emit: reconstruction_accuracy, optional: true
    path "cell2location_reconstruction_accuracy_error.txt", emit: reconstruction_error, optional: true
    path "cell2location_joint_diagnostics.tsv", emit: joint_diagnostics
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'), emit: versions_container, topic: versions

    script:
    def args = task.ext.args ?: ''
    if (sample_ids.size() != zarr_folders.size()) {
        error "CELL2LOCATION_REINSERT received ${sample_ids.size()} sample IDs and ${zarr_folders.size()} SpatialData paths"
    }
    def sample_rows = (0..<sample_ids.size()).collect { index ->
        "${sample_ids[index]}\t${zarr_folders[index]}"
    }.join('\n')
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}
    export XDG_CACHE_HOME=\${TMPDIR:-/tmp}/.cache
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"

    cat << 'END_SAMPLES' > samples.tsv
sample_id	spatialdata
${sample_rows}
END_SAMPLES

    cat << 'END_SCRIPT' > reinsert.py
${file("${moduleDir}/reinsert.py").text}
END_SCRIPT

    python3 reinsert.py \
        ${args} \
        --samples samples.tsv \
        --posterior-dir "${posterior_dir}" \
        --posterior-manifest "${posterior_manifest}" \
        --training-history "${training_history}" \
        --diagnostics "${diagnostics}" \
        --output-dir .
    """

    stub:
    def zarr_outputs = sample_ids.collect { sample_id ->
        "mkdir -p ${sample_id}_cell2location.zarr"
    }.join('\n')
    def csv_outputs = sample_ids.collect { sample_id ->
        "touch ${sample_id}_cell2location_abundance_mean.csv ${sample_id}_cell2location_abundance_q05.csv"
    }.join('\n')
    """
    ${zarr_outputs}
    ${csv_outputs}
    touch cell2location_training_history.png
    touch cell2location_reconstruction_accuracy.png
    touch cell2location_joint_diagnostics.tsv
    """
}
