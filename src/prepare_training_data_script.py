import argparse
import scanpy as sc
from pre_processing import *

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--input", required=True, type=str, help="path to input h5ad data")
parser.add_argument("--dataset", required=True, type=str, help="name of dataset (ROSMAP / SEA-AD)")
parser.add_argument("--metadata", required=False, type=str, help="path to ROSMAP metadata folder")
parser.add_argument("-t", "--temp", required=True, type=str, help="path to temporary files folder")
# parser.add_argument("-o", "--output", required=True, type=str, help="path to output folder")

args = parser.parse_args()

# subset data to 1k most hvg
# try:
#     filtered_adata_1k = sc.read_h5ad(
#         os.path.join(args.temp, str(args.dataset.lower()) + "_filtered_1k_genes_data.h5ad"))
# except FileNotFoundError as e:
# print(e)
adata_sc = sc.read_h5ad(args.input)

# use raw counts
if args.dataset == "SEA-AD":
    print(adata_sc.obs['Donor ID'].unique().tolist)
    adata_sc = adata_sc[adata_sc.obs['Subclass'].isin(['Astrocyte', 'Microglia-PVM', 'Oligodendrocyte'])]
    print(adata_sc.n_obs)
    adata_sc.X = adata_sc.layers['UMIs']

adata_1k = sc_data_to_1k(adata_sc, args.temp, args.dataset.lower()+"_top_1k_genes_data")
del adata_sc

if args.dataset == "ROSMAP":
    # label
    labeled_adata_1k = extract_labels_from_ROSMAP_metadata(
        adata_1k, args.temp, args.dataset.lower()+"_labeled_1k_genes_data", args.metadata)
    del adata_1k

    # remove other
    filtered_adata_1k = remove_unknown_class_records(
        labeled_adata_1k, args.temp, args.dataset.lower()+"_filtered_1k_genes_data", args.metadata)
    del labeled_adata_1k
elif args.dataset == "SEA-AD":
    adata_1k.obs['y'] = [0 if x is False else 1 for x in adata_1k.obs['Overall AD neuropathological Change'] == 'High']
    filtered_adata_1k = adata_1k
    del adata_1k
else:
    raise IOError("Unsupported dataset")

# convert to pca
# pca_training_data = pca_labeled_data(filtered_adata_1k, args.output, args.dataset.lower()+"_pca_training_data")
del filtered_adata_1k

# del pca_training_data
