process SEGMENTATION_AND_MICROSCOPY_PLOTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(zarr_folder), path(segmentation_mask), path(microscopy_image), val(crop_areas)
    val zarr_downsample_factor
    val microscopy_downsample_factor
    val chunk_size

    output:
    tuple val(meta), path("*_segmentation_microscopy.zarr"), optional: true, emit: zarr
    tuple val(meta), path("*_mqc.png"), emit: mqc_plots
    tuple val(meta), path("*_segmentation_stats.csv"), emit: segmentation_stats
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'), emit: versions_container, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def output_zarr_arg = task.ext.save_zarr ? "--output-zarr \"${prefix}_segmentation_microscopy.zarr\"" : ''
    def crop_areas_arg = crop_areas ? "--crop-areas \"${crop_areas}\"" : ''
    def overwrite_arg = task.ext.overwrite ? '--overwrite' : ''
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}
    export XDG_CACHE_HOME=\${TMPDIR:-/tmp}/.cache
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"

    cat << 'END_SCRIPT' > segmentation_and_microscopy_plots.py
${file("${moduleDir}/segmentation_and_microscopy_plots.py").text}
END_SCRIPT

    python3 segmentation_and_microscopy_plots.py \
        ${args} \
        --input-zarr "${zarr_folder}" \
        ${output_zarr_arg} \
        --sample-name "${prefix}" \
        --segmentation-mask-tif "${segmentation_mask}" \
        --microscopy-tif "${microscopy_image}" \
        --zarr-downsample-factor "${zarr_downsample_factor}" \
        --microscopy-downsample-factor "${microscopy_downsample_factor}" \
        --chunk-size "${chunk_size}" \
        ${crop_areas_arg} \
        ${overwrite_arg} \
        --results-dir "."
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def save_zarr = task.ext.save_zarr
    """
    mkdir -p ${save_zarr ? "${prefix}_segmentation_microscopy.zarr" : "stub_no_zarr"}
    touch ${prefix}_registration_full_slide_mqc.png
    touch ${prefix}_crop_areas_downsampled_microscopy_mqc.png
    touch ${prefix}_segmentation_crop_panels_mqc.png
    touch ${prefix}_segmentation_stats.csv

    """
}
