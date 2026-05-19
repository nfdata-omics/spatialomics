process ASSEMBLE_IMAGING_MULTIQC {
    tag "$meta.id"
    label 'process_low'

    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(images)

    output:
    tuple val(meta), path("*_imaging_and_segmentation_mqc.html"), emit: html
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'), emit: versions_container, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat << 'END_SCRIPT' > assemble_imaging_multiqc.py
${file("${moduleDir}/assemble_imaging_multiqc.py").text}
END_SCRIPT

    python3 assemble_imaging_multiqc.py \\
        --sample-name "${prefix}" \\
        --output "${prefix}_imaging_and_segmentation_mqc.html" \\
        --images ${images}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_id = prefix.replaceAll('[^A-Za-z0-9_]+', '_').replaceAll('^_+|_+$', '') ?: 'sample'
    def sample_name_yaml = groovy.json.JsonOutput.toJson(prefix)
    """
    cat << 'END_HTML' > ${prefix}_imaging_and_segmentation_mqc.html
<!--
parent_id: imaging_and_segmentation
parent_name: "Imaging and segmentation"
parent_description: |
  Per-sample microscopy, registration, spatial QC, crop-area, and segmentation review outputs.
id: "imaging_and_segmentation_${sample_id}"
section_name: ${sample_name_yaml}
plot_type: "html"
-->
<div class="spatialomics-sample-report">
  <p>Stub imaging and segmentation report for ${prefix}.</p>
</div>
END_HTML
    """
}
