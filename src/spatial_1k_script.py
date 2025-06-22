import argparse
import os
import scanpy as sc
from pre_processing import st_data_to_1k

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--input", required=True, type=str, help="path to input data folder")
parser.add_argument("-o", "--output", required=True, type=str, help="path to output data folder")

args = parser.parse_args()

scdata = sc.read_h5ad(os.path.join(args.output, "rosmap_filtered_1k_genes_data.h5ad"))

odata_st = sc.read_h5ad(os.path.join(args.input, "rosmap_imputed_spatial_data_oligodendroglia.h5ad"))

st_data_to_1k(odata_st, scdata, args.output, "rosmap_oligodendroglia_imputed_xenium_1k")

del odata_st

rdata_st = sc.read_h5ad(os.path.join(args.input, "rosmap_imputed_spatial_data.h5ad"))

st_data_to_1k(rdata_st, scdata, args.output, "rosmap_imputed_xenium_1k")

del rdata_st

del scdata

scdata = sc.read_h5ad(os.path.join(args.output, "sea-ad_top_1k_genes_data.h5ad"))

sdata_st = sc.read_h5ad(os.path.join(args.input, "sea-ad_imputed_spatial_data.h5ad"))

st_data_to_1k(sdata_st, scdata, args.output, "sea-ad_imputed_xenium_1k")

del sdata_st
