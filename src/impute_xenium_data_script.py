import argparse
import scanpy as sc
from pre_processing import *
from impute_missing_genes_spatial_data import raw_impute

parser = argparse.ArgumentParser(description="")

parser.add_argument("--spatial_input", required=True, type=str, help="path to Xenium h5ad data")
parser.add_argument("--reference_input", required=True, type=str, help="path to prepared reference h5ad data")
parser.add_argument("--dataset", required=True, type=str, help="name of reference dataset (ROSMAP / SEA-AD)")
parser.add_argument("--epochs", required=True, type=int, help="number of epochs for imputation (default=5)")
parser.add_argument("-t", "--temp", required=True, type=str, help="path to temporary files folder")
parser.add_argument("-o", "--output", required=True, type=str, help="path to output model folder")

args = parser.parse_args()

idata_sc = sc.read_h5ad(args.reference_input)

adata_st = sc.read_h5ad(args.spatial_input)
print(adata_st.n_obs)

# Do this when imputing using ROSMAP
cell_type = ''
if 'astrocytes' in str(args.reference_input):
    adata_st = adata_st[adata_st.obs['predicted.celltype'].isin(['Astrocyte'])]
    cell_type = '_astrocytes'
elif 'oligodendroglia' in str(args.reference_input):
    adata_st = adata_st[adata_st.obs['predicted.celltype'].isin(['Oligodendrocyte'])]
    cell_type = '_oligodendroglia'
elif 'microglia' in str(args.reference_input):
    adata_st = adata_st[adata_st.obs['predicted.celltype'].isin(['Microglia'])]
    cell_type = '_microglia'
print(adata_st.n_obs)

raw_impute(idata_sc, adata_st, args.output, args.dataset.lower()+"_imputed_spatial_data"+cell_type)

print("Succesfully imputed spatial data.")

del adata_st
del idata_sc
