process BIN2CELL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'quay.io/biocontainers/bin2cell:0.3.4--pyhdfd78af_0'

    input:
    tuple val(meta), path(spaceranger_outs), path(segmentation_mask), path(source_image)
    val bin_size
    val volume_ratio

    output:
    tuple val(meta), path("*_bin2cell_cells.h5ad"), emit: h5ad
    tuple val(meta), path("*_bin2cell_labels.npz"), emit: labels_npz
    tuple val(meta), path("*_bin2cell_summary.csv"), emit: summary
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}
    export XDG_CACHE_HOME=\${TMPDIR:-/tmp}/.cache
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"

    cat << 'END_SCRIPT' > run_bin2cell.py
${file("${moduleDir}/run_bin2cell.py").text}
END_SCRIPT

    python3 run_bin2cell.py \\
        ${args} \\
        --sample-name "${prefix}" \\
        --spaceranger-outs "${spaceranger_outs}" \\
        --segmentation-mask-tif "${segmentation_mask}" \\
        --source-image "${source_image}" \\
        --bin-size "${bin_size}" \\
        --volume-ratio "${volume_ratio}" \\
        --output-h5ad "${prefix}_bin2cell_cells.h5ad" \\
        --output-labels-npz "${prefix}_bin2cell_labels.npz" \\
        --output-summary "${prefix}_bin2cell_summary.csv"

    python3 run_bin2cell.py \\
        --versions-dict "${task.process}" > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_bin2cell_cells.h5ad"
    touch "${prefix}_bin2cell_labels.npz"
    cat << 'END_SUMMARY' > "${prefix}_bin2cell_summary.csv"
sample,input_bins,assigned_bins,unassigned_bins,assigned_fraction,unassigned_fraction,segmentation_labels,output_cells,volume_ratio,bin_size
${prefix},0,0,0,0.0,0.0,0,0,${volume_ratio},${bin_size}
END_SUMMARY

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bin2cell: unknown
    END_VERSIONS
    """
}
