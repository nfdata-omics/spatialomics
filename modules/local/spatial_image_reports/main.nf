process SPATIAL_IMAGE_REPORTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(zarr_folder), val(crop_areas)
    val resolution
    val zarr_downsample_factor
    val chunk_size
    val segmentation_enabled

    output:
    tuple val(meta), path("*_registration_tissue_mqc.webp"), emit: registration_plot
    tuple val(meta), path("*_crop_segmentation_mqc.png"), optional: true, emit: crop_segmentation_plot
    tuple val(meta), path("*_segmentation_stats.csv"), optional: true, emit: segmentation_stats
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'), emit: versions_container, topic: versions

    script:
    def segmentation_arg = segmentation_enabled ? '--segmentation' : ''
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}
    export XDG_CACHE_HOME=\${TMPDIR:-/tmp}/.cache
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"

    cat << 'END_LAYOUT' > spatial_preview_layout.py
${file("${moduleDir}/../spatial_qc_plots/spatial_preview_layout.py").text}
END_LAYOUT

    cat << 'END_SEGMENTATION' > segmentation_report.py
${file("${moduleDir}/segmentation_report.py").text}
END_SEGMENTATION

    cat << 'END_SCRIPT' > spatial_image_reports.py
${file("${moduleDir}/spatial_image_reports.py").text}
END_SCRIPT

    python3 spatial_image_reports.py \\
        --zarr "${zarr_folder}" \\
        --sample "${meta.id}" \\
        --resolution "${resolution}" \\
        --zarr-downsample-factor "${zarr_downsample_factor}" \\
        --chunk-size "${chunk_size}" \\
        --crop-areas "${crop_areas}" ${segmentation_arg}
    """

    stub:
    def segmentation_outputs = segmentation_enabled ?
        "touch \"${meta.id}_crop_segmentation_mqc.png\" \"${meta.id}_segmentation_stats.csv\"" : ''
    """
    touch "${meta.id}_registration_tissue_mqc.webp"
    ${segmentation_outputs}
    """
}
