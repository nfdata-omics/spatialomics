process SPACERANGER_TO_ZARR {
    tag "$meta.id"
    label 'process_single'

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

    import json
    import os
    from pathlib import Path

    os.environ["NUMBA_CACHE_DIR"] = os.environ.get("TMPDIR", "/tmp")
    os.environ["MPLCONFIGDIR"] = os.environ.get("TMPDIR", "/tmp")
    os.environ["XDG_CONFIG_HOME"] = os.environ.get("TMPDIR", "/tmp")

    import sys
    import importlib
    import pkg_resources
    import geopandas as gpd
    import h5py
    import pandas as pd
    import yaml
    import spatialdata_io
    from shapely.geometry import box
    from spatialdata.models import ShapesModel
    from spatialdata.transformations import Identity

    # ----------------------------------
    # Load dataset using spatialdata_io
    # ----------------------------------

    # Here, the bin size is set to 16, but we should allow the option to also read 002 and 008
    spaceranger_path = Path("$spaceranger_output_dir")
    feature_slice_paths = list(spaceranger_path.glob("*feature_slice.h5"))
    if len(feature_slice_paths) != 1:
        raise FileNotFoundError(
            "Expected exactly one feature_slice.h5 file in the Space Ranger output, "
            f"found {len(feature_slice_paths)}"
        )
    with h5py.File(feature_slice_paths[0], "r") as feature_slice:
        metadata_json = feature_slice.attrs.get("metadata_json")
        metadata = json.loads(metadata_json) if metadata_json is not None else {}

    transform_matrices = metadata.get("transform_matrices", {})
    required_transform_keys = [
        "cytassist_colrow_to_spot_colrow",
        "spot_colrow_to_cytassist_colrow",
        "spot_colrow_to_microscope_colrow",
        "microscope_colrow_to_spot_colrow",
    ]

    def has_nine_values(matrix):
        if not isinstance(matrix, list):
            return False
        if len(matrix) == 9:
            return True
        return (
            len(matrix) == 3
            and all(isinstance(row, list) and len(row) == 3 for row in matrix)
        )

    has_registration_matrices = all(
        transform_matrices.get(key) is not None
        and has_nine_values(transform_matrices[key])
        for key in required_transform_keys
    )
    has_cytassist_image = any(spaceranger_path.rglob("cytassist_image.tiff"))
    load_registered_cytassist = has_cytassist_image and has_registration_matrices
    if not load_registered_cytassist:
        print(
            "Registered CytAssist image will not be loaded: "
            f"image present={has_cytassist_image}, "
            f"complete registration matrices={has_registration_matrices}"
        )

    data = spatialdata_io.visium_hd(
        "$spaceranger_output_dir",
        filtered_counts_file="${filtered_counts_file}",
        dataset_id="${prefix}",
        load_segmentations_only=False,
        load_all_images=load_registered_cytassist
    )

    # The filtered count matrix contains only bins selected as in-tissue, so
    # their bounds cannot define the complete capture area. Preserve one compact
    # capture-area geometry derived from all 16 um bin positions instead of
    # loading the raw expression matrix or millions of additional shapes.
    bin_spatial_dir = (
        spaceranger_path
        / "binned_outputs"
        / "square_016um"
        / "spatial"
    )
    positions_path = bin_spatial_dir / "tissue_positions.parquet"
    scalefactors_path = bin_spatial_dir / "scalefactors_json.json"
    if not positions_path.exists() or not scalefactors_path.exists():
        raise FileNotFoundError(
            "Cannot construct the Visium HD capture-area geometry: expected "
            f"{positions_path} and {scalefactors_path}"
        )

    positions = pd.read_parquet(
        positions_path,
        columns=["pxl_col_in_fullres", "pxl_row_in_fullres"],
    )
    if positions.empty:
        raise ValueError(f"Space Ranger tissue positions are empty: {positions_path}")

    with scalefactors_path.open() as handle:
        scalefactors = json.load(handle)
    half_bin_width = float(scalefactors["spot_diameter_fullres"]) / 2.0
    capture_area = box(
        positions["pxl_col_in_fullres"].min() - half_bin_width,
        positions["pxl_row_in_fullres"].min() - half_bin_width,
        positions["pxl_col_in_fullres"].max() + half_bin_width,
        positions["pxl_row_in_fullres"].max() + half_bin_width,
    )
    capture_area_gdf = gpd.GeoDataFrame(geometry=[capture_area])
    data.shapes["${prefix}_capture_area"] = ShapesModel.parse(
        capture_area_gdf,
        transformations={"${prefix}": Identity()},
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
    for lib in [
        'spatialdata_io',
        'spatialdata',
        'numpy',
        'pandas',
        'scipy',
        'h5py',
        'geopandas',
        'shapely',
    ]:
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
    import sys
    import importlib
    import pkg_resources
    import yaml
    os.makedirs("${meta.id}.zarr", exist_ok=True)

    # ----------------------------------
    # Print versions of relevant libraries
    # ----------------------------------

    versions = {}
    versions["${task.process}"] = {}
    for lib in [
        'spatialdata_io',
        'spatialdata',
        'numpy',
        'pandas',
        'scipy',
        'h5py',
        'geopandas',
        'shapely',
    ]:
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
