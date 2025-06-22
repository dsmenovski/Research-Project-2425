import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

# We do all the scanpy processing, following this tutorial:
# https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering-2017.html#preprocessing


def preprocess_scaled_data():
    adata = sc.read_h5ad(os.path.join("/tudelft.net", "staff-umbrella", "bachelorAD", "data", "ROSMAP",
                                      "microglia.h5ad"))

    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=200)

    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    adata = adata[adata.obs.pct_counts_mt < 5, :].copy()

    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, n_top_genes=1000)

    adata = adata[:, adata.var.highly_variable].copy()

    # Fix from https://github.com/theislab/scvelo/issues/255#issuecomment-739995301
    adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

    adata.write(Path(os.path.join(
        "/tudelft.net", "staff-umbrella", "bachelorAD", "dsmenovski", "bsc-rp-2425-dimitar-smenovski", "data",
        "dataset_1k.h5ad")),
        compression="gzip")


def sc_data_to_1k(adata_sc, output_folder: str, file_name: str):
    sc.pp.filter_cells(adata_sc, min_genes=100)
    sc.pp.filter_genes(adata_sc, min_cells=3)

    adata_sc.var["mt"] = adata_sc.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata_sc, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata_sc = adata_sc[adata_sc.obs.pct_counts_mt < 5, :].copy()

    sc.pp.highly_variable_genes(adata_sc, flavor="seurat_v3", n_top_genes=1000)
    adata_sc = adata_sc[:, adata_sc.var.highly_variable].copy()

    sc.pp.normalize_total(adata_sc)
    sc.pp.log1p(adata_sc)

    sc.pp.scale(adata_sc)

    adata_sc.write(Path(os.path.join(output_folder, file_name+".h5ad")), compression="gzip")

    return adata_sc


def st_data_to_1k(adata_st, adata_sc, output_folder: str, file_name: str):
    sc.pp.filter_cells(adata_st, min_genes=100)
    sc.pp.filter_genes(adata_st, min_cells=3)

    adata_st.var["mt"] = adata_st.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata_st, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata_st = adata_st[adata_st.obs.pct_counts_mt < 5, :].copy()

    lc = [var.lower() for var in adata_sc.var_names]
    adata_st = adata_st[:, adata_st.var_names.isin(lc)].copy()

    sc.pp.normalize_total(adata_st)
    sc.pp.log1p(adata_st)

    sc.pp.scale(adata_st)

    adata_st.write(Path(os.path.join(output_folder, file_name+".h5ad")), compression="gzip")

    return adata_st


def extract_labels_from_ROSMAP_metadata(adata, output_folder: str, file_name: str, metadata_folder: str):
    df = pd.read_csv(os.path.join(metadata_folder, "ROSMAP_clinical.csv"))

    # Determine AD and healthy individuals according to Wang (2021)
    ad = df.query('cogdx == 4 and braaksc >= 4 and ceradsc <= 2')

    df_labeled = df.copy()
    df_labeled['y'] = np.nan

    df_labeled.loc[df_labeled['individualID'].isin(ad['individualID']), 'y'] = 1
    df_labeled.loc[~df_labeled['individualID'].isin(ad['individualID']), 'y'] = 0

    # Generate mapping from metadata labels to AnnData records
    y = np.zeros(shape=(adata.n_obs, 1))

    for i in range(adata.n_obs):
        try:
            y[i][0] = df_labeled.loc[df_labeled['individualID'] == adata.obs['individualID'].iloc[i], 'y'].values[0]
        except IndexError:
            # Skip rows where individualID is missing
            y[i][0] = 0

    adata.obs['y'] = y

    adata.write(Path(os.path.join(output_folder, file_name+".h5ad")), compression="gzip")

    return adata


def remove_unknown_class_records(adata, output_folder: str, file_name: str, metadata_folder: str):
    df = pd.read_csv(os.path.join(metadata_folder, "ROSMAP_clinical.csv"))

    # Determine AD and healthy individuals according to Wang (2021)
    ad = df.query('cogdx == 4 and braaksc >= 4 and ceradsc <= 2')
    ct = df.query('cogdx == 1 and braaksc <= 3 and ceradsc >= 3')

    rows = [x for x in adata.obs['individualID'].to_numpy()
            if x in ad['individualID'].to_numpy() or x in ct['individualID'].to_numpy()]

    adata = adata[adata.obs['individualID'].isin(rows)].copy()

    adata.write(Path(os.path.join(output_folder, file_name+".h5ad")), compression="gzip")

    return adata


def prepare_sc_data_for_imputation(adata_sc, output_folder: str, file_name: str):
    adata_sc.var["mt"] = adata_sc.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata_sc, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata_sc = adata_sc[adata_sc.obs.pct_counts_mt < 5, :].copy()

    sc.pp.normalize_total(adata_sc)
    sc.pp.log1p(adata_sc)

    sc.tl.pca(adata_sc, random_state=56)

    sc.pp.neighbors(adata_sc, random_state=56)
    sc.tl.leiden(adata_sc, flavor="leidenalg", random_state=56, resolution=0.25)
    sc.tl.umap(adata_sc, random_state=56)
    # sc.pl.umap(adata_sc, color="leiden")

    if output_folder == '' and file_name == '':
        print('Not saving data...')
        return adata_sc

    adata_sc.write(Path(os.path.join(output_folder, file_name+".h5ad")), compression="gzip")

    return adata_sc


# first this
def filter_imputed_spatial_data(adata_sc, adata_st, output_folder: str, file_name: str):
    adata_genes = [gene.lower() for gene in set(adata_sc.var_names)]
    adata_st = adata_st[:, adata_genes].copy()
    adata_st.write(Path(os.path.join(output_folder, file_name+".h5ad")))


# then this
def pca_labeled_data(adata, output_folder: str, file_name: str):
    sc.tl.pca(adata, n_comps=50, random_state=56)
    adata.write(Path(os.path.join(output_folder, file_name+".h5ad")))

    return adata


def pca_spatial_data(adata, output_folder: str, file_name: str):
    sc.pp.highly_variable_genes(adata, n_top_genes=1000)
    adata = adata[:, adata.var.highly_variable].copy()
    sc.pp.scale(adata)
    adata.write(Path(os.path.join(output_folder, file_name+".h5ad")))

    return adata
