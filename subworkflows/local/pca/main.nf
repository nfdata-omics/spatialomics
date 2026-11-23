include { SC_PREPROCESS_NORMALIZE } from '../../../modules/sc_preprocess/normalize'
include { SC_PREPROCESS_HVG       } from '../../../modules/sc_preprocess/hvg'
include { SC_PREPROCESS_PCA       } from '../../../modules/sc_preprocess/pca'

workflow PCA {
    take:
    h5ad // tuple val(meta), path(h5ad), with raw RNA counts in .X
    n_top_genes // optional, number of top highly variable genes to select for PCA
    color // required, name of the .obs column to use for coloring

    main:
    SC_PREPROCESS_NORMALIZE(h5ad)
    SC_PREPROCESS_HVG(SC_PREPROCESS_NORMALIZE.out.h5ad, n_top_genes)
    SC_PREPROCESS_PCA(SC_PREPROCESS_HVG.out.h5ad, color)

    emit:
    h5ad = SC_PREPROCESS_PCA.out.h5ad
    plot = SC_PREPROCESS_PCA.out.plot
    variance_plot = SC_PREPROCESS_PCA.out.variance_plot
    scores = SC_PREPROCESS_PCA.out.scores
    explained_variance = SC_PREPROCESS_PCA.out.explained_variance
}
