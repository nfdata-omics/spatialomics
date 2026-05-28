process CYTASSIST_MICROSCOPY_REGISTRATION {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/nfdata/cytassist-microscopy-registration:local'

    input:
    tuple val(meta), path(microscopy_image), path(cytassist_image)

    output:
    tuple val(meta), path("*_cytassist_to_microscopy_transform.json"), emit: transform
    tuple val(meta), path("*_cytassist_microscopy_registration_metrics.csv"), emit: metrics
    tuple val(meta), path("*_cytassist_microscopy_registration_overlay_mqc.png"), emit: overlay
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

    cat << 'END_SCRIPT' > cytassist_microscopy_registration.py
${file("${moduleDir}/cytassist_microscopy_registration.py").text}
END_SCRIPT

    python3 cytassist_microscopy_registration.py \\
        ${args} \\
        --cytassist-image "${cytassist_image}" \\
        --microscopy-tif "${microscopy_image}" \\
        --sample-name "${prefix}" \\
        --output-transform "${prefix}_cytassist_to_microscopy_transform.json" \\
        --output-metrics "${prefix}_cytassist_microscopy_registration_metrics.csv" \\
        --output-overlay "${prefix}_cytassist_microscopy_registration_overlay_mqc.png"

    python3 cytassist_microscopy_registration.py \\
        --versions-dict "${task.process}" > versions.yml
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat > "${prefix}_cytassist_to_microscopy_transform.json" <<END_JSON
{
  "sample": "${prefix}",
  "transform_direction": "cytassist_to_microscopy",
  "status": "stub"
}
END_JSON

    cat > "${prefix}_cytassist_microscopy_registration_metrics.csv" <<END_CSV
sample,status,method
${prefix},stub,stub
END_CSV

    touch "${prefix}_cytassist_microscopy_registration_overlay_mqc.png"

    cat > versions.yml <<END_VERSIONS
"${task.process}":
    python: unknown
END_VERSIONS
    """
}
