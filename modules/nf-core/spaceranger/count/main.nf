process SPACERANGER_COUNT {
    tag "$meta.id"
    label 'process_high'

    container "nf-core/spaceranger:3.1.3"

    input:
    tuple val(meta), path("fastqs/${meta.id}_S1_L001_R?_001.fastq.gz")
    tuple val(meta2), path(image), val(slide), val(area), path(cytaimage), path(darkimage), path(colorizedimage), path(alignment), path(slidefile)
    path(reference)
    path(probeset)

    output:
    tuple val(meta), path("*_web_summary.html"), emit: web_summary
    tuple val(meta), path("outs"), emit: outs
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPACERANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Add flags for optional inputs on demand.
    def probeset_arg = probeset ? "--probe-set=\"${probeset}\"" : ""
    def alignment_arg = alignment ? "--loupe-alignment=\"${alignment}\"" : ""
    def slidefile_arg = slidefile ? "--slidefile=\"${slidefile}\"" : ""
    def image_arg = image ? "--image=\"${image}\"" : ""
    def cytaimage_arg = cytaimage ? "--cytaimage=\"${cytaimage}\"" : ""
    def darkimage_arg = darkimage ? "--darkimage=\"${darkimage}\"" : ""
    def colorizedimage_arg = colorizedimage ? "--colorizedimage=\"${colorizedimage}\"" : ""
    if (slide.matches("visium-(.*)") && area == "" && slidefile == "") {
        slide_and_area = "--unknown-slide=\"${slide}\""
    } else {
        slide_and_area = "--slide=\"${slide}\" --area=\"${area}\""
    }
    """
    spaceranger count \\
        --id="${prefix}" \\
        --sample="${meta.id}" \\
        --fastqs=fastqs \\
        --transcriptome="${reference}" \\
        --localcores=${task.cpus} \\
        --localmem=${task.memory.toGiga()} \\
        $image_arg $cytaimage_arg $darkimage_arg $colorizedimage_arg \\
        $slide_and_area \\
        $probeset_arg \\
        $alignment_arg \\
        $slidefile_arg \\
        $args
    mv ${prefix}/outs outs
    mv outs/web_summary.html ${prefix}_web_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spaceranger: \$(spaceranger -V | sed -e "s/spaceranger spaceranger-//g")
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "SPACERANGER_COUNT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p outs/
    touch outs/fake_file.txt
    touch ${prefix}_web_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        spaceranger: \$(spaceranger -V | sed -e "s/spaceranger spaceranger-//g")
    END_VERSIONS
    """
}
