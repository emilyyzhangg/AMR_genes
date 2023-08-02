import scanpy as sc
import utils
import pandas as pd
adata = sc.read('./data/MrVIoutputs/bacdrop.h5ad')

ranked_genes_gauss = utils.mrvi_identify_cell_states(adata, 'X_mrvi_z', neighbor_method='gauss')
ranked_genes_knn = utils.mrvi_identify_cell_states(adata, 'X_mrvi_z', neighbor_method='knn')

n_top_genes = 30

top_gauss = utils.filter_duplicates(ranked_genes_gauss, n_top_genes)
top_knn = utils.filter_duplicates(ranked_genes_knn, n_top_genes)

out_dir = "./out/"

top_gauss.to_csv(out_dir + "gauss.csv")
top_knn.to_csv(out_dir + "knn.csv")