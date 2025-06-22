import argparse
import os

import pandas as pd
import scanpy as sc
import scipy.stats as st
import tangram as tg

from pre_processing import prepare_sc_data_for_imputation

# Follwing the tutorial at https://www.sc-best-practices.org/spatial/imputation.html#using-tangram-in-practise

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--input", required=True, type=str, help="path to input folder")
parser.add_argument("-o", "--output", required=True, type=str, help="path to output folder")
parser.add_argument("-d", "--dataset", required=True, type=str, help="name of reference dataset (ROSMAP / SEA-AD)")
parser.add_argument("-e", "--epochs", type=int, help="number of epochs")

args = parser.parse_args()

if args.dataset == 'ROSMAP':
    adata_sc = sc.read_h5ad(
        os.path.join(args.input,
                     "rosmap_data_for_imputation.h5ad"))
if args.dataset == 'SEA-AD':
    adata_sc = sc.read_h5ad(
        os.path.join(args.input,
                     "sea-ad_data_for_imputation.h5ad"))

adata_st = sc.read_h5ad(
    os.path.join(args.input,
                 "xenium_spatial_ad_adata.h5ad"))

adata_st = adata_st[adata_st.obs['predicted.celltype'].isin(['Astrocyte', 'Microglia', 'Oligodendrocyte'])].copy()

sc.pp.normalize_total(adata_st)
sc.pp.log1p(adata_st)

adata_st.var["mt"] = adata_st.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata_st, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
adata_st = adata_st[adata_st.obs.pct_counts_mt < 5, :].copy()

common = list(set(adata_st.var_names) & set(adata_sc.var_names))

adata_st = adata_st[:, common]

top_genes = pd.read_csv(os.path.join(args.input, "top_genes.csv"))

top_gene_per_cluster = top_genes.sort_values('scores', ascending=False).groupby('group', observed=True).first()

markers = list(set(top_gene_per_cluster['names'].to_numpy()))

tg.pp_adatas(adata_sc, adata_st, genes=markers)

assert "training_genes" in adata_sc.uns
assert "training_genes" in adata_st.uns

print(f"Number of training_genes: {len(adata_sc.uns['training_genes'])}")

lc_markers = [x.lower() for x in markers]
cp_markers = lc_markers.copy()
correlations = pd.Series(index=markers)

for gene in markers:
    mask = adata_st.var_names != gene
    adata_st_subset = adata_st[:, mask].copy()

    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata_st_subset,
        mode="cells",
        density_prior="rna_count_based",
        random_state=56,
        num_epochs=args.epochs,
        device="cpu",  # or: cuda
    )
    ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)

    correlations[gene] = st.spearmanr(adata_st[:, gene.lower()].X.toarray().ravel().copy(),
                                      ad_ge[:, gene.lower()].X.toarray().ravel().copy())[0]
    print(correlations[gene])

print('entire correlations dataframe:')
print(correlations)
correlations.to_csv(os.path.join(args.output, str(args.dataset)+"_correlations.csv"))
