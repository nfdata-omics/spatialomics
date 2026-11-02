# Functional test fixture

## Objective and scope

The default `test` profile exercises the downstream spatialomics workflow from
pre-computed Space Ranger output using the compact fixture published at
<https://doi.org/10.5281/zenodo.23888419>. The fixture is intended for pipeline
development and continuous integration.

## Agreed behaviour

- `assets/samplesheet_test_mock.csv` keeps the samplesheet in the repository
  and references the version-specific Zenodo record for all data files.
- The Space Ranger input is supplied as `spaceranger_outs.tar.gz` and extracted
  by the pipeline before downstream processing.
- The test uses `skimage_watershed` segmentation to avoid a GPU dependency.
- Cell2location is not enabled because this fixture has no suitable single-cell
  reference.
- The previous FASTQ-based test inputs, Space Ranger reference, probe set, and
  Cellpose test override are retained under the `test_full` profile.
- Both profiles cap resources at 4 CPUs, 12 GB memory, and one hour per task for
  compatibility with CI runners.

## Verification status

- The Zenodo record and published filenames have been checked.
- Configuration and samplesheet validation are checked statically as part of
  implementation.
- A realistic functional run remains to be performed separately.
