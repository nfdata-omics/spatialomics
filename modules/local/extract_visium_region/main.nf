process EXTRACT_VISIUM_REGION {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(zarr_folder), path(microscopy_image)
    val zarr_downsample_factor
    val crop_padding
    val crop_downsample_factor

    output:
    tuple val(meta), path("*_visium_bounds.tsv"), emit: bounds
    tuple val(meta), path("*_visium_capture_area.ome.tif"), emit: crop
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}

    cat << 'END_SCRIPT' > extract_visium_region.py
${file("${moduleDir}/extract_visium_region.py").text}
END_SCRIPT

    python3 extract_visium_region.py \\
        --input-zarr "${zarr_folder}" \\
        --sample-name "${prefix}" \\
        --microscopy-tif "${microscopy_image}" \\
        --zarr-downsample-factor "${zarr_downsample_factor}" \\
        --crop-padding "${crop_padding}" \\
        --crop-downsample-factor "${crop_downsample_factor}" \\
        --output-bounds "${prefix}_visium_bounds.tsv" \\
        --output-crop "${prefix}_visium_capture_area.ome.tif"

    python3 extract_visium_region.py \\
        --versions-dict "${task.process}" > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat << 'END_BOUNDS' > "${prefix}_visium_bounds.tsv"
x0	y0	x1	y1	crop_x0	crop_y0	crop_x1	crop_y1	source_width	source_height	padding	downsample_factor	output_width	output_height	scale_x	scale_y
0	0	1	1	0	0	1	1	1	1	${crop_padding}	${crop_downsample_factor}	1	1	1	1
END_BOUNDS

    touch "${prefix}_visium_capture_area.ome.tif"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    spatialdata: unknown
END_VERSIONS
    """
}
