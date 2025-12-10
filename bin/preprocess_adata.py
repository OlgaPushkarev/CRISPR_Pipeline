#!/usr/bin/env python
#!/usr/bin/env python
import argparse
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import scipy

matplotlib.rcParams["axes.spines.top"] = False
matplotlib.rcParams["axes.spines.right"] = False

def is_outlier(adata_rna, metric: str, nmads: int, verbose=False):
    """
    Taken from single cell best practices
    https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
    """
    M = adata_rna.obs[metric]
    lower_bound = np.median(M) - nmads * scipy.stats.median_abs_deviation(M)
    upper_bound = np.median(M) + nmads * scipy.stats.median_abs_deviation(M)
    if ("mt" in metric) or ("ribo" in metric):
        lower_bound = 0
    outlier = (M < lower_bound) | (upper_bound < M)

    if verbose:
        print("-" * 10)
        print(f"Non-outlier boundaries for {metric}:")
        print(f"Lower boundary: {lower_bound}")
        print(f"Upper boundary: {upper_bound}")
        return outlier, lower_bound, upper_bound
    else:
        return outlier

def main(
        adata_rna,
        gname_rna,
        filter_outliers,
        min_genes_per_cell,
        min_cells_per_gene,
        min_counts_per_cell,
        min_counts_per_gene,
        n_mads,
        reference,
        output_dir
    ):

    gname_rna_path = os.path.join(gname_rna, 'counts_unfiltered/cells_x_genes.genes.names.txt')
    gene_df = pd.read_csv(gname_rna_path, header=None)
    gene_names = gene_df[0].tolist()
    adata_rna = sc.read(adata_rna)

    if len(gene_names) == adata_rna.shape[1]:
        # keep ensembl id as var_names
        # matched gene names as column: gene_id
        # adata_rna.var_names = gene_names
        adata_rna.var["symbol"] = gene_names
        # modify the ensembl id in adata_rna.var_names
        adata_rna.var_names = adata_rna.var_names.str.split('.').str[0]
        adata_rna.var_names_make_unique()
    else: 
        raise ValueError("The number of gene names does not match the number of variables in adata_rna")

    os.makedirs('figures', exist_ok=True)

    # knee plots
    knee_df = pd.DataFrame({
        'sum': np.array(adata_rna.X.sum(axis=1)).flatten(),
        'barcodes': adata_rna.obs_names.values})
    knee_df = knee_df.sort_values('sum', ascending=False).reset_index(drop=True)
    knee_df['sum_log'] = np.log1p(knee_df['sum'])

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    ax.plot(knee_df.index, knee_df['sum_log'], marker='o', linestyle='-', markersize=3)
    ax.set_xlabel('Barcodes', size=12)
    ax.set_ylabel('$\log_{1+p}$( UMI Counts )', size=12)
    ax.set_title('Knee plot', size=12)
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, 'knee_plot_scRNA.png'),
        dpi=300,
        bbox_inches="tight",
        transparent=True
    )
    plt.close()

    # Add batch number
    adata_rna.obs['batch_number'] = adata_rna.obs['batch'].factorize()[0] + 1

    if reference == "human":
        mt_prefix = "MT-"
        ribo_prefix = ("RPS", "RPL")
    else:
        mt_prefix = "Mt-"
        ribo_prefix = ("Rps", "Rpl")

    adata_rna.var["mt"] = adata_rna.var['symbol'].str.startswith(mt_prefix)
    adata_rna.var["ribo"] = adata_rna.var['symbol'].str.startswith(ribo_prefix)

    # Calculate QC metrics
    adata_rna.X = adata_rna.X.astype(np.float32)
    sc.pp.calculate_qc_metrics(
        adata_rna,
        qc_vars=["mt", "ribo"],
        inplace=True,
        log1p=True
    )

    # Plot violin

    # 'n_genes_by_counts':# genes per cell
    # 'total_counts' = # UMIs per cell = sequencing depth = library size

    labeling_dict = {
        'n_genes_by_counts': '# genes per cell',
        'total_counts': '# UMIs per cell',
        'pct_counts_mt': '% mitochondrial UMIs',
        'pct_counts_ribo': '% ribosomal UMIs',
        "S_score": "S score",
        "G2M_score": "G2M score"
    }

    qc_metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    fig, ax = plt.subplots(
        1,
        len(qc_metrics),
        figsize=(len(qc_metrics) * 1.75, 5),
    )
    for i, qc_metric in enumerate(qc_metrics):
        outlier_mask, low, high = is_outlier(adata_rna, qc_metric, n_mads, verbose=True)
        low = max(0, low)  # Ensure non-negative lower bound

        sc.pl.violin(
            adata_rna,
            qc_metric,
            jitter=0.2,
            ax=ax[i],
            color="gray",
            show=False,
            size=0.5,
        )
        for coll in ax[i].collections:
            coll.set_zorder(2)
    
        data = adata_rna.obs[qc_metric].dropna()
        median_val = np.median(data)
        ax[i].hlines(
            y=median_val,
            xmin=-0.075,
            xmax=0.075,
            color="black",
            lw=1.5,
            zorder=4,
            label="median" if i == 0 else None,
        )
        for artist in ax[i].findobj(matplotlib.collections.PathCollection):
            artist.set_rasterized(True)
    
        data = adata_rna.obs[qc_metric].dropna()
        ax[i].text(
            0.5,
            1.05,
            f"threshold={str(round(high, 2))}\nmedian={str(round(median_val, 2))}\n",
            ha="center",
            va="bottom",
            transform=ax[i].transAxes,
            fontsize=8,
        )
        ax[i].axhline(low, ls="--", lw=0.5, color="red", zorder=3)
        ax[i].axhline(high, ls="--", lw=0.5, color="red", zorder=3)
        
        ax[i].set_ylabel(labeling_dict.get(qc_metric, qc_metric) if i == 0 else "", size=12)
        ax[i].set_xlabel("")
    plt.tight_layout()
    plt.savefig(
        os.path.join(output_dir, "plot_qcs_violin.png"),
        dpi=300,
        bbox_inches="tight",
        transparent=True
    )
    plt.close()

    # Plot scatter
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    x_label = "total_counts"
    y_label = "n_genes_by_counts"
    color_by = "pct_counts_mt"

    x = adata_rna.obs[x_label]
    y = adata_rna.obs[y_label]
    c = adata_rna.obs[color_by]

    mask = ~(x.isna() | y.isna() | c.isna())
    x, y, c = x[mask], y[mask], c[mask]
    scatter = ax.scatter(
        x,
        y,
        c=c,
        cmap="Blues",
        s=10,
        alpha=0.5,
        linewidth=0,
        rasterized=True
    )
    ax.set_xlabel(labeling_dict.get(x_label), size=12)
    ax.set_ylabel(labeling_dict.get(y_label), size=12)
    cbar = plt.colorbar(scatter, ax=ax)
    cbar.set_label(labeling_dict.get(color_by))
    plt.savefig(
        os.path.join(output_dir, "plot_total_gene_umis_vs_ngenes_by_counts_mt_color_before_filter_scatter.png"),
        dpi=300,
        bbox_inches="tight",
        transparent=True
    )
    plt.close()

    print(f"Before filtering: {adata_rna.n_obs} cells x {adata_rna.n_vars} genes")

    adata_rna.obs["outlier"] = (
        is_outlier(adata_rna, "log1p_total_counts", n_mads, verbose=False)
        | is_outlier(adata_rna, "log1p_n_genes_by_counts", n_mads, verbose=False)
        | is_outlier(adata_rna, "pct_counts_mt", n_mads, verbose=False)
        | is_outlier(adata_rna, "pct_counts_ribo", n_mads, verbose=False)
    )

    print(f"# outlier cells: {adata_rna.obs['outlier'].sum()}")
    # adata_rna.write(os.path.join(args.output_dir, 'unfiltered_scRNA_anndata.h5ad'))

    if filter_outliers:
        adata_rna = adata_rna[~adata_rna.obs.outlier].copy()
        print(f"After filtering: {adata_rna.n_obs} cells")
    else:
        print("Skipping outlier removal (filter_outliers=False)")

    # if add_extra_filter:
        # sc.pp.filter_cells(adata_rna, min_genes=min_genes_per_cell)
        # sc.pp.filter_cells(adata_rna, min_counts=min_counts_per_cell)
    
    sc.pp.filter_genes(adata_rna, min_cells=min_cells_per_gene)
    sc.pp.filter_genes(adata_rna, min_counts=min_counts_per_gene)
    print(f"\nFinal: {adata_rna.n_obs} cells x {adata_rna.n_vars} genes")
    adata_rna.write('filtered_anndata.h5ad')
    adata_rna.write(os.path.join(args.output_dir, 'filtered_anndata.h5ad'))

    print("Quality control completed and saved to 'filtered_anndata.h5ad'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Perform QC on AnnData.')
    parser.add_argument('adata_rna', type=str, help='Path to the AnnData file.')
    parser.add_argument('gname_rna', type=str, help='Path to the cells x genes txt file.')
    parser.add_argument('--filter_outliers', type=bool, default=True, help='Remove the outliers')
    parser.add_argument('--min_genes_per_cell', type=int, default=100, help='Minimum number of genes per cell')
    parser.add_argument('--min_cells_per_gene', type=int, default=10, help='Minimum number of cells per gene')
    parser.add_argument('--min_counts_per_cell', type=int, default=500, help='Minimum number of UMIs per cell')
    parser.add_argument('--min_counts_per_gene', type=int, default=20, help='Minimum number of UMIs per gene')
    parser.add_argument('--n_mads', type=float, default=5, help='n median absolute deviations to filteer outliers')
    parser.add_argument('--reference', type=str, required=True, help='Reference species')
    parser.add_argument('--output_dir', required=True, help='Directory where plots will be saved')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)


    args = parser.parse_args()
    main(
        args.adata_rna,
        args.gname_rna,
        args.filter_outliers,
        args.min_genes_per_cell,
        args.min_cells_per_gene,
        args.min_counts_per_cell,
        args.min_counts_per_gene,
        args.n_mads,
        args.reference,
        args.output_dir
    )
