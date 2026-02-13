process COLLECT_SPACERANGER_METRICS {
    tag "all samples"
    label 'process_low'

    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    path "outs_*"

    output:
    path "spaceranger_metrics.csv",  emit: metrics
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python3

    import sys
    import importlib
    import pkg_resources
    import glob
    import yaml
    import pandas as pd

    input_files = glob.glob("outs_*/metrics_summary.csv")

    dfs = []
    for file in input_files:
        df = pd.read_csv(file)
        dfs.append(df)

    # Concatenate all DataFrames
    combined = pd.concat(dfs, ignore_index=True)

    # Save to a single CSV file
    combined.to_csv("spaceranger_metrics.csv", index=False)

    # ----------------------------------
    # Print versions of relevant libraries
    # ----------------------------------

    versions = {}
    versions["${task.process}"] = {}
    for lib in ['pandas']:
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

    import sys
    import importlib
    import pkg_resources
    import glob
    import yaml

    open('spaceranger_metrics.csv', 'w').close()

    # ----------------------------------
    # Print versions of relevant libraries
    # ----------------------------------

    versions = {}
    versions["${task.process}"] = {}
    for lib in ['pandas']:
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
