process SPACERANGER_TO_ZARR {
    tag "$meta.id"
    label 'process_low'

    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    tuple val(meta), path(spaceranger_output_dir)
    val filtered_counts_file

    output:
    tuple val(meta), path("*.zarr"), emit: zarr
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env python3

    import os

    os.environ["NUMBA_CACHE_DIR"] = os.environ.get("TMPDIR", "/tmp")
    os.environ["MPLCONFIGDIR"] = os.environ.get("TMPDIR", "/tmp")
    os.environ["XDG_CONFIG_HOME"] = os.environ.get("TMPDIR", "/tmp")

    import sys
    import importlib
    import pkg_resources
    import yaml
    import spatialdata_io

    # ----------------------------------
    # Load dataset using spatialdata_io
    # ----------------------------------

    # Here, the bin size is set to 16, but we should allow the option to also read 002 and 008
    data = spatialdata_io.visium_hd(
        "$spaceranger_output_dir",
        filtered_counts_file="${filtered_counts_file}",
        dataset_id="${prefix}"
    )

    # ----------------------------------
    # Save full dataset as Zarr
    # ----------------------------------
    data.write(f"${prefix}.zarr", overwrite=True)

    # ----------------------------------
    # Print versions of relevant libraries
    # ----------------------------------

    versions = {}
    versions["${task.process}"] = {}
    for lib in ['spatialdata_io', 'spatialdata', 'numpy', 'pandas', 'scipy']:
        try:
            version = pkg_resources.get_distribution(lib).version
        except Exception:
            try:
                module = importlib.import_module(lib)
                version = getattr(module, '__version__', 'unknown')
            except Exception:
                version = None
        if version is not None:
            versions["${task.process}"][lib] = version
    versions["${task.process}"]['python'] = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"

    with open('versions.yml', 'w') as f:
        yaml.dump(versions, f)
    """

    stub:
    """
    #!/usr/bin/env python3

    import os
    os.makedirs("${meta.id}.zarr", exist_ok=True)

    # ----------------------------------
    # Print versions of relevant libraries
    # ----------------------------------

    versions = {}
    versions["${task.process}"] = {}
    for lib in ['spatialdata_io', 'spatialdata', 'numpy', 'pandas', 'scipy']:
        try:
            version = pkg_resources.get_distribution(lib).version
        except Exception:
            try:
                module = importlib.import_module(lib)
                version = getattr(module, '__version__', 'unknown')
            except Exception:
                version = None
        if version is not None:
            versions["${task.process}"][lib] = version
    versions["${task.process}"]['python'] = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"

    with open('versions.yml', 'w') as f:
        yaml.dump(versions, f)

    """


}
