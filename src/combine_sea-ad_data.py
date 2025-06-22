import argparse
import os
from pathlib import Path
import anndata as ad
import scanpy as sc

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--input", required=True, type=str, help="path to input folder")
parser.add_argument("-o", "--output", required=True, type=str, help="path to output folder")

args = parser.parse_args()

try:
    adata = sc.read_h5ad(os.path.join(args.input, "combined_objects_sea-ad.h5ad"))
except FileNotFoundError as e:
    print(e)
    adata1 = sc.read_h5ad(os.path.join(args.input, "SEA-AD", "donor_objects",
                                       "H19.33.004_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
    adata2 = sc.read_h5ad(os.path.join(args.input, "SEA-AD", "donor_objects",
                                       "H20.33.001_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
    adata3 = sc.read_h5ad(os.path.join(args.input, "SEA-AD", "donor_objects",
                                       "H20.33.002_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
    adata4 = sc.read_h5ad(os.path.join(args.input, "SEA-AD", "donor_objects",
                                       "H20.33.032_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
    adata5 = sc.read_h5ad(os.path.join(args.input, "SEA-AD", "donor_objects",
                                       "H20.33.036_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
    adata6 = sc.read_h5ad(os.path.join(args.input, "SEA-AD", "donor_objects",
                                       "H21.33.001_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
    adata7 = sc.read_h5ad(os.path.join(args.input, "SEA-AD", "donor_objects",
                                       "H21.33.003_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
    adata8 = sc.read_h5ad(os.path.join(args.input, "SEA-AD", "donor_objects",
                                       "H21.33.032_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))
    adata9 = sc.read_h5ad(os.path.join(args.input, "SEA-AD", "donor_objects",
                                       "H21.33.045_SEAAD_A9_RNAseq_final-nuclei.2024-02-13.h5ad"))

    adata = ad.concat([adata1, adata2, adata3, adata4, adata5, adata6, adata7, adata8, adata9], axis=0)
    del adata1, adata2, adata3, adata4, adata5, adata6, adata7, adata8, adata9

    adata = adata[adata.obs['Subclass'].isin(['Astrocyte', 'Microglia-PVM', 'Oligodendrocyte'])]

    adata.write(Path(os.path.join(args.output, "combined_objects_sea-ad.h5ad")))

del adata
