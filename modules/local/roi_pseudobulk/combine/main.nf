process COMBINE_ROI_PSEUDOBULK {
    label 'process_single'

    container 'docker.io/nfdata/spatialdata:v0.7.2'

    input:
    path(samples, stageAs: 'samples/*')

    output:
    path("roi_pseudobulk.h5ad"), emit: h5ad
    tuple path("roi_pseudobulk_counts.tsv.gz"), path("roi_pseudobulk_rois.tsv"), path("roi_pseudobulk_genes.tsv"), emit: tsv
    tuple val("${task.process}"), val('container'), val('docker.io/nfdata/spatialdata:v0.7.2'),
        emit: versions_container, topic: versions

    script:
    """
    roi_pseudobulk.py combine --output roi_pseudobulk.h5ad --export-tsv samples/*.h5ad
    """

    stub:
    """
    touch roi_pseudobulk.h5ad roi_pseudobulk_counts.tsv.gz roi_pseudobulk_rois.tsv roi_pseudobulk_genes.tsv
    """
}
