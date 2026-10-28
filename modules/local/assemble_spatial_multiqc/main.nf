process ASSEMBLE_SPATIAL_MULTIQC {
    tag "$meta.id"
    label 'process_single'

    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(images)

    output:
    tuple val(meta), path("*_spatial_mqc.html"), emit: html
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'), emit: versions_container, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat << 'END_SCRIPT' > assemble_spatial_multiqc.py
${file("${moduleDir}/assemble_spatial_multiqc.py").text}
END_SCRIPT

    python3 assemble_spatial_multiqc.py \\
        --sample-name "${prefix}" \\
        --output "${prefix}_spatial_mqc.html" \\
        --images ${images}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_id = prefix.replaceAll('[^A-Za-z0-9_]+', '_').replaceAll('^_+|_+$', '') ?: 'sample'
    def sample_name_yaml = groovy.json.JsonOutput.toJson(prefix)
    """
    cat << 'END_HTML' > ${prefix}_spatial_mqc.html
<!--
parent_id: spatial_overview
parent_name: "Spatial overview"
parent_description: |
  Per-sample spatial QC, imaging, segmentation, and deconvolution review outputs.
id: "spatial_overview_${sample_id}"
section_name: ${sample_name_yaml}
plot_type: "html"
-->
<div class="spatialomics-sample-report">
  <p>Stub spatial report for ${prefix}.</p>
</div>
END_HTML
    """
}
