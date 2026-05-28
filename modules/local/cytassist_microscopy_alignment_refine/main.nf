process CYTASSIST_MICROSCOPY_ALIGNMENT_REFINE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/cytassist-microscopy-registration:local'

    input:
    tuple val(meta), path(microscopy_image), path(cytassist_image), path(initial_transform)

    output:
    tuple val(meta), path("*_refined_cytassist_to_microscopy_transform.json"), emit: transform
    tuple val(meta), path("*_cytassist_microscopy_alignment_refinement_metrics.csv"), emit: metrics
    tuple val(meta), path("*_cytassist_microscopy_alignment_refinement_before_mqc.png"), emit: overlay_before
    tuple val(meta), path("*_cytassist_microscopy_alignment_refinement_after_mqc.png"), emit: overlay_after
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MPLBACKEND=Agg
    export MPLCONFIGDIR=\${TMPDIR:-/tmp}
    export XDG_CONFIG_HOME=\${TMPDIR:-/tmp}
    export XDG_CACHE_HOME=\${TMPDIR:-/tmp}/.cache
    mkdir -p "\$XDG_CACHE_HOME/fontconfig"

    cat << 'END_SCRIPT' > cytassist_microscopy_alignment_refine.py
${file("${moduleDir}/cytassist_microscopy_alignment_refine.py").text}
END_SCRIPT

    python3 cytassist_microscopy_alignment_refine.py \\
        ${args} \\
        --cytassist-image "${cytassist_image}" \\
        --microscopy-tif "${microscopy_image}" \\
        --initial-transform "${initial_transform}" \\
        --sample-name "${prefix}" \\
        --output-transform "${prefix}_refined_cytassist_to_microscopy_transform.json" \\
        --output-metrics "${prefix}_cytassist_microscopy_alignment_refinement_metrics.csv" \\
        --output-overlay-before "${prefix}_cytassist_microscopy_alignment_refinement_before_mqc.png" \\
        --output-overlay-after "${prefix}_cytassist_microscopy_alignment_refinement_after_mqc.png"

    python3 cytassist_microscopy_alignment_refine.py \\
        --versions-dict "${task.process}" > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat > "${prefix}_refined_cytassist_to_microscopy_transform.json" <<END_JSON
{
  "sample": "${prefix}",
  "transform_direction": "cytassist_to_microscopy",
  "method": "local_affine_refinement",
  "status": "stub"
}
END_JSON

    cat > "${prefix}_cytassist_microscopy_alignment_refinement_metrics.csv" <<END_CSV
sample,status,method
${prefix},stub,local_affine_refinement
END_CSV

    touch "${prefix}_cytassist_microscopy_alignment_refinement_before_mqc.png"
    touch "${prefix}_cytassist_microscopy_alignment_refinement_after_mqc.png"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    python: unknown
END_VERSIONS
    """
}
