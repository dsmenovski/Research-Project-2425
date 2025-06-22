import os
import argparse

import numpy as np
import scanpy as sc

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--input", required=True, type=str, help="path to input folder")
parser.add_argument("-o", "--output", required=True, type=str, help="path to output folder")

args = parser.parse_args()

adata_st = sc.read_h5ad(
    os.path.join(args.input,
                 "xenium_spatial_ad_adata.h5ad"))

adata_st = adata_st[adata_st.obs['predicted.celltype'].isin(['Astrocyte', 'Microglia', 'Oligodendrocyte'])].copy()

sc.pp.normalize_total(adata_st)
sc.pp.log1p(adata_st)

sc.pp.neighbors(adata_st, random_state=56)

sc.tl.umap(adata_st, random_state=56)

np.random.seed(0)
sc.tl.leiden(adata_st, flavor='leidenalg', random_state=56)

sc.tl.rank_genes_groups(adata_st, groupby='leiden', method='wilcoxon')

top_genes_df = sc.get.rank_genes_groups_df(adata_st, group=None)

top_genes_df.to_csv(os.path.join(args.output, "top_genes.csv"))
