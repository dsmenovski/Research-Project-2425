import argparse
import os
from pathlib import Path
import anndata as ad
import scanpy as sc

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--input", required=True, type=str, help="path to input folder")
parser.add_argument("-o", "--output", required=True, type=str, help="path to output folder")

args = parser.parse_args()

adata = sc.read_h5ad(os.path.join(args.input, "astrocytes.h5ad"))
mdata = sc.read_h5ad(os.path.join(args.input, "microglia.h5ad"))
odata = sc.read_h5ad(os.path.join(args.input, "oligodendroglia.h5ad"))

fdata = ad.concat([adata, mdata, odata], axis=0)

cell_types = ["Astrocyte"] * adata.n_obs + ["Microglia"] * mdata.n_obs + ["Oligodendrocyte"] * odata.n_obs

fdata.obs["cell_type"] = cell_types

del adata, mdata, odata

print(fdata.obs.head())
print(fdata.obs["cell_type"].value_counts())

fdata.write(Path(os.path.join(args.output, "full_rosmap_data.h5ad")))
