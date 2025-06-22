import argparse
import os
from pathlib import Path
import anndata as ad
import scanpy as sc

from pre_processing import pca_spatial_data

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--input", required=True, type=str, help="path to input folder")
parser.add_argument("-o", "--output", required=True, type=str, help="path to output folder")

args = parser.parse_args()

# only needed for the rosmap imputed data
adata = sc.read_h5ad(os.path.join(args.input, 'rosmap_imputed_spatial_data_'+"astrocytes.h5ad"))
mdata = sc.read_h5ad(os.path.join(args.input, 'rosmap_imputed_spatial_data_'+"microglia.h5ad"))
odata = sc.read_h5ad(os.path.join(args.input, 'rosmap_imputed_spatial_data_'+"oligodendroglia.h5ad"))

rdata = ad.concat([adata, mdata, odata], axis=0)

pca_spatial_data(adata, args.output, "a-rosmap_imputed_xenium_pca")
pca_spatial_data(mdata, args.output, "m-rosmap_imputed_xenium_pca")
pca_spatial_data(odata, args.output, "o-rosmap_imputed_xenium_pca")

del adata, mdata, odata

rdata.write(Path(os.path.join(args.output, "rosmap_imputed_spatial_data.h5ad")))

# pca_spatial_data(rdata, args.output, "rosmap_imputed_xenium_pca")

del rdata

sdata = sc.read_h5ad(os.path.join(args.output, "sea-ad_imputed_spatial_data.h5ad"))

# pca_spatial_data(sdata, args.output, "sea-ad_imputed_xenium_pca")

del sdata
