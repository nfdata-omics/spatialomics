process VISIUM_BOUNDS {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(zarr_folder), path(microscopy_image)
    val zarr_downsample_factor

    output:
    tuple val(meta), path("*_visium_bounds.tsv"), emit: bounds
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}

    cat << 'END_SCRIPT' > visium_bounds.py
${file("${moduleDir}/visium_bounds.py").text}
END_SCRIPT

    python3 visium_bounds.py \\
        --input-zarr "${zarr_folder}" \\
        --sample-name "${prefix}" \\
        --microscopy-tif "${microscopy_image}" \\
        --zarr-downsample-factor "${zarr_downsample_factor}" \\
        --output-bounds "${prefix}_visium_bounds.tsv"

    python3 visium_bounds.py \\
        --versions-dict "${task.process}" > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat << 'END_BOUNDS' > "${prefix}_visium_bounds.tsv"
x0	y0	x1	y1
0	0	1	1
END_BOUNDS

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    spatialdata: unknown
END_VERSIONS
    """
}
