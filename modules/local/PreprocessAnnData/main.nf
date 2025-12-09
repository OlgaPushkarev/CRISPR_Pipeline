
process PreprocessAnnData {

    cache 'lenient'

    input:
    path adata_rna
    path gname_rna
    val filter_outliers
    val min_genes_per_cell
    val min_cells_per_gene
    val min_counts_per_cell
    val min_counts_per_gene
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
            --min_genes_per_cell ${min_genes_per_cell} \\
            --min_cells_per_gene ${min_cells_per_gene} \\
            --min_counts_per_cell ${min_counts_per_cell} \\
            --min_counts_per_gene ${min_counts_per_gene} \\
            --n_mads ${n_mads} \\
            --reference ${reference} \\
            --output_dir figures
        mv concatenated_adata.h5ad rna_concatenated_adata.h5ad
        """

}
