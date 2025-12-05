
process PreprocessAnnData {

    cache 'lenient'

    input:
    path adata_rna
    path gname_rna
    val filter_outliers
    val min_genes
    val min_cells
    val n_mads
    val reference

    output:
    path "filtered_anndata.h5ad" , emit: filtered_anndata_rna
    path "rna_concatenated_adata.h5ad", emit: adata_rna
    path "figures", emit: figures_dir

    script:
        """
        preprocess_adata.py ${adata_rna} ${gname_rna} \\
            --filter_outliers ${filter_outliers} \\
            --min_genes ${min_genes} \\
            --min_cells ${min_cells} \\
            --n_mads ${n_mads} \\
            --reference ${reference} \\
            --output_dir figures
        mv concatenated_adata.h5ad rna_concatenated_adata.h5ad
        """

}
