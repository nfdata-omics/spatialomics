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
    import numpy as np
    import pandas as pd
    import yaml
    import spatialdata_io
    from PIL import Image
    from shapely.geometry import Polygon
    from skimage.transform import ProjectiveTransform
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
    cytassist_paths = list(spaceranger_path.rglob("cytassist_image.tiff"))
    load_registered_cytassist = bool(cytassist_paths) and has_registration_matrices
    if load_registered_cytassist and len(cytassist_paths) != 1:
        raise ValueError(f"Expected one CytAssist TIFF for ${prefix}, found {len(cytassist_paths)}")
    if not load_registered_cytassist:
        print(
            "Registered CytAssist image will not be loaded: "
            f"image present={bool(cytassist_paths)}, "
            f"complete registration matrices={has_registration_matrices}"
        )

    data = spatialdata_io.visium_hd(
        "$spaceranger_output_dir",
        filtered_counts_file="${filtered_counts_file}",
        dataset_id="${prefix}",
        load_segmentations_only=False,
        load_all_images=load_registered_cytassist
    )

    def store_footprint(name, vertices):
        polygon = Polygon(vertices)
        if not polygon.is_valid or polygon.area <= 0:
            raise ValueError(f"Invalid {name} polygon for ${prefix}")
        data.shapes[f"${prefix}_{name}"] = ShapesModel.parse(
            gpd.GeoDataFrame(geometry=[polygon]),
            transformations={"${prefix}": Identity()},
        )

    if load_registered_cytassist:
        # Project the original TIFF support. The registered raster can include
        # fill pixels outside this polygon after resampling.
        with Image.open(cytassist_paths[0]) as image:
            width, height = image.size
        cytassist_to_array = np.asarray(transform_matrices["cytassist_colrow_to_spot_colrow"]).reshape(3, 3)
        array_to_microscope = np.asarray(transform_matrices["spot_colrow_to_microscope_colrow"]).reshape(3, 3)
        cytassist_to_microscope = ProjectiveTransform(array_to_microscope @ cytassist_to_array)
        corners = [
            (-0.5, -0.5), (width - 0.5, -0.5),
            (width - 0.5, height - 0.5), (-0.5, height - 0.5),
        ]
        store_footprint("cytassist_footprint", cytassist_to_microscope(corners))

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
    if not positions_path.exists():
        raise FileNotFoundError(f"Cannot construct the Visium HD capture area without {positions_path}")

    positions = pd.read_parquet(
        positions_path,
        columns=["array_col", "array_row", "pxl_col_in_fullres", "pxl_row_in_fullres"],
    ).set_index(["array_col", "array_row"])
    if positions.empty:
        raise ValueError(f"Space Ranger tissue positions are empty: {positions_path}")

    # Use the four corner bin centers to map array indices into image
    # coordinates, then extend by half a bin in array space.
    cols = positions.index.get_level_values("array_col")
    rows = positions.index.get_level_values("array_row")
    left, right, top, bottom = cols.min(), cols.max(), rows.min(), rows.max()
    center_indices = [(left, top), (right, top), (right, bottom), (left, bottom)]
    corner_positions = positions.reindex(center_indices)
    if corner_positions.isna().any().any() or left == right or top == bottom:
        raise ValueError(f"Expected four distinct corner bin centers in {positions_path}")
    bins_to_image = ProjectiveTransform.from_estimate(center_indices, corner_positions.to_numpy())
    if not bins_to_image:
        raise ValueError(f"Cannot project the capture-area corners from {positions_path}")
    store_footprint(
        "capture_area",
        bins_to_image([
            (left - 0.5, top - 0.5), (right + 0.5, top - 0.5),
            (right + 0.5, bottom + 0.5), (left - 0.5, bottom + 0.5),
        ]),
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
