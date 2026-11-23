process ROI_PSEUDOBULK {
    tag "$meta.id"
    label 'process_medium'

    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(zarr)
    val bin_label

    output:
    tuple val(meta), path("${meta.id}_roi_pseudobulk.h5ad"), emit: h5ad
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'),
        emit: versions_container, topic: versions

    script:
    def sample = groovy.json.JsonOutput.toJson(
        [sample_id: meta.id] + meta.subMap(['slide', 'area']).findAll { _key, value -> value }
    )
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}

    cat << 'END_SAMPLE' > sample.json
${sample}
END_SAMPLE

    roi_pseudobulk.py aggregate \\
        --zarr "${zarr}" \\
        --bin-label "${bin_label}" \\
        --sample sample.json \\
        --output "${meta.id}_roi_pseudobulk.h5ad"
    """

    stub:
    """
    touch "${meta.id}_roi_pseudobulk.h5ad"
    """
}
