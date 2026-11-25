# Final SpatialData assembly and spatial reporting

## Status and scope

Implemented. Mixed-run, mock 1 geometry and `-resume` validation in the full pipeline are pending.

The final Zarr will contain the H&E supplied through the samplesheet whenever
it is available. Segmentation is an optional addition. Spatial QC and final
assembly will expose their results independently of rendering, so changes to
spatial figures can reuse cached analytical outputs and the assembled Zarr.
The existing Cell2location, ROI import and pseudobulk process boundaries stay
in place.

## Analytical spatial QC

- Keep `SPATIAL_QUALITY_CONTROL` upstream of Cell2location, ROI import and
  final assembly. It continues to emit the QC-annotated Zarr, per-bin
  annotations and per-sample QC metrics.
- Preserve the selected bin resolution, QC calculations, thresholds,
  `in_tissue`, `qc_pass` and missing-value semantics.
- Move the per-sample QC distributions PNG, spatial QC WebP and registration
  WebP out of the analytical job. Use the existing measurements, annotations
  and thresholds for the distributions and spatial QC maps. Preserve the
  distribution populations and the maps' comparison between unfiltered and
  QC-passing in-tissue bins. The cross-sample `COLLECT_QC` report already reads
  analytical outputs independently.

## Final SpatialData assembly

Use the existing late augmentation step as the final assembler. Its base is
the per-sample Zarr selected after QC, optional Cell2location and optional ROI
import. Associate it with the prepared samplesheet image and any segmentation
mask by sample ID.

| Inputs for the sample                          | Final SpatialData                                                                  |
| ---------------------------------------------- | ---------------------------------------------------------------------------------- |
| H&E and segmentation mask                      | Base Zarr plus native-resolution H&E, its reduced preview and segmentation labels. |
| H&E only, including with `--skip_segmentation` | Base Zarr plus native-resolution H&E and its reduced preview.                      |
| No H&E                                         | Publish the selected upstream Zarr directly.                                       |

- Use the prepared version of the samplesheet `image`, including the existing
  CZI conversion where applicable. Preserve the single-channel
  `<sample>_microscopy` contract at native spatial resolution and the existing
  reduced microscopy layer. Both layers cover the full image field. Keep the
  established `segmentation_mask` element and its coordinate semantics.
- Give the new layers explicit transforms into the sample coordinate system.
  Reuse the existing `zarr_downsample_factor` convention and account for the
  preview's reduction in resolution. The user remains responsible for
  matching an external H&E to precomputed Space Ranger coordinates, as stated
  in the existing warning.
- Preserve inherited elements, tables and metadata. Write the H&E and
  optional segmentation labels at final assembly, into a distinct Zarr
  output. Publish one final dataset at
  `spatialdata/<sample>_spatialdata.zarr/`.

## Spatial reporting

- Generate the registration and tissue-selection preview from the assembled
  Zarr. The left panel shows the full field of the reduced H&E derived from
  the samplesheet, with the registered CytAssist in pseudocolor and the
  capture-area outline. Use the stored SpatialData transformations and set
  the viewport to the transformed H&E extent.
- Keep the right panel focused on the capture area and show the original
  Space Ranger `in_tissue` selection over the same H&E. When H&E is available
  without a registered CytAssist layer, show only the tissue-selection panel.
  Omit this WebP for samples without a samplesheet H&E. Make the MultiQC
  description accurate for the panels present.
- Generate the crop/segmentation preview and segmentation statistics CSV
  outside the final Zarr writer, only for samples with a mask. Preserve the
  existing crop selection and statistics for the whole image and the Visium
  region.
- Keep the existing report filenames, publication destinations and MultiQC
  integration. Reporting tasks read upstream outputs as immutable inputs.
  Put rendering code and settings in reporting tasks so edits to them leave
  analytical and final-assembly task inputs and scripts unchanged. The
  number of reporting processes is an implementation choice.

## Acceptance criteria

- A mixed run publishes one correctly associated final Zarr per sample across
  all three input cases. QC, Cell2location and ROI results carry through as
  applicable. The new H&E and segmentation layers enter the Zarr lineage at
  final assembly.
- With `--skip_segmentation`, a sample with H&E retains both H&E layers and
  its tissue-selection preview, plus the registration panel when registered
  CytAssist is available. A sample without H&E still produces its analytical
  QC outputs.
- For mock 1, the four boundaries of the registration panel match the field
  of the H&E from the samplesheet. CytAssist is placed according to its stored
  transform and is clipped at the H&E boundary.
- Per-sample QC reports and segmentation statistics retain their existing
  numerical results for equivalent inputs. Changing only the rendering code
  permits `-resume` to reuse the analytical and assembly tasks. Runtime checks
  are left to the user.

## Related specifications

- [SpatialData lifecycle and serialization](spatialdata_lifecycle.md)
- [MultiQC spatial preview layout and encoding](multiqc_spatial_previews.md)
