include { SC_PREPROCESS_QC        } from '../../../modules/local/sc_preprocess/qc'
include { SC_PREPROCESS_NORMALIZE } from '../../../modules/local/sc_preprocess/normalize'
include { SC_PREPROCESS_HVG       } from '../../../modules/local/sc_preprocess/hvg'
include { SC_PREPROCESS_PCA       } from '../../../modules/local/sc_preprocess/pca'

workflow PSEUDOBULK_PREPROCESS {
    take:
    h5ad // tuple val(meta), path(h5ad), with raw RNA counts in .X
    min_counts // minimum total counts per ROI
    n_top_genes // optional, number of top highly variable genes to select for PCA
    color // required, name of the .obs column to use for coloring

    main:
    SC_PREPROCESS_QC(h5ad, min_counts)
    SC_PREPROCESS_NORMALIZE(SC_PREPROCESS_QC.out.h5ad)
    SC_PREPROCESS_HVG(SC_PREPROCESS_NORMALIZE.out.h5ad, n_top_genes)
    SC_PREPROCESS_PCA(SC_PREPROCESS_HVG.out.h5ad, color)

    emit:
    gex_pdf = SC_PREPROCESS_QC.out.gex_pdf
    html = SC_PREPROCESS_QC.out.html
    multiqc_html = SC_PREPROCESS_QC.out.multiqc_html
    multiqc_table = SC_PREPROCESS_QC.out.multiqc_table
    h5ad = SC_PREPROCESS_PCA.out.h5ad
    plot = SC_PREPROCESS_PCA.out.plot
    variance_plot = SC_PREPROCESS_PCA.out.variance_plot
    scores = SC_PREPROCESS_PCA.out.scores
    multiqc_scores = SC_PREPROCESS_PCA.out.multiqc_scores
    explained_variance = SC_PREPROCESS_PCA.out.explained_variance
}
