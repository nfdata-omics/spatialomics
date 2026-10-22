process CELL2LOCATION_PLOTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(zarr_folder)
    val table_key
    val abundance_statistic
    val uncertainty_min_abundance
    val multiqc_cell_types

    output:
    tuple val(meta), path("*_posterior_mean_vs_q05.pdf"),                    emit: posterior_mean_vs_q05
    tuple val(meta), path("*_relative_posterior_uncertainty_spatial.pdf"),  emit: posterior_uncertainty
    tuple val(meta), path("*_celltype_abundance_*.pdf"),                   emit: celltype_abundance
    tuple val(meta), path("*_celltype_fraction_*.pdf"),                    emit: celltype_fraction
    tuple val(meta), path("*_total_abundance_*.pdf"),                      emit: total_abundance
    tuple val(meta), path("*_deconvolution_mqc.png"),                      emit: mqc_plot
    tuple val(meta), path("*_celltype_composition_*_mqc.tsv"),             emit: composition_mqc
    tuple val(meta), path("*_deconvolution_summary_*_mqc.tsv"),            emit: summary_mqc
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'), emit: versions_container, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def multiqc_cell_types_arg = multiqc_cell_types ? "--multiqc-cell-types \"${multiqc_cell_types}\"" : ''
    """
    export NUMBA_CACHE_DIR=\${TMPDIR:-/tmp}
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}
    export XDG_CACHE_HOME=\${TMPDIR:-/tmp}/.cache
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"

    cat << 'END_SCRIPT' > plot_cell2location.py
${file("${moduleDir}/plot_cell2location.py").text}
END_SCRIPT

    python3 plot_cell2location.py \
        ${args} \
        --input-zarr "${zarr_folder}" \
        --sample-name "${meta.id}" \
        --output-prefix "${prefix}" \
        --output-dir . \
        --table-key "${table_key}" \
        --abundance-statistic "${abundance_statistic}" \
        --uncertainty-min-abundance "${uncertainty_min_abundance}" \
        ${multiqc_cell_types_arg}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}_posterior_mean_vs_q05.pdf"
    touch "${prefix}_relative_posterior_uncertainty_spatial.pdf"
    touch "${prefix}_celltype_abundance_${abundance_statistic}.pdf"
    touch "${prefix}_celltype_fraction_${abundance_statistic}.pdf"
    touch "${prefix}_total_abundance_${abundance_statistic}.pdf"
    touch "${prefix}_deconvolution_mqc.png"
    touch "${prefix}_celltype_composition_${abundance_statistic}_mqc.tsv"
    touch "${prefix}_deconvolution_summary_${abundance_statistic}_mqc.tsv"
    """
}
