import os
from pathlib import Path
import scanpy as sc
import tangram as tg

sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor="white")

# Follwing the tutorial at https://www.sc-best-practices.org/spatial/imputation.html#using-tangram-in-practise


def dup(adata_sc, adata_st, num_epochs=5):

    markers = list(set(adata_st.var_names) & set(adata_sc.var_names))

    adata_st = adata_st[:, markers]

    tg.pp_adatas(adata_sc, adata_st, genes=markers)

    assert "training_genes" in adata_sc.uns
    assert "training_genes" in adata_st.uns

    print(f"Number of training_genes: {len(adata_sc.uns['training_genes'])}")

    ad_map = tg.map_cells_to_space(
        adata_sc,
        adata_st,
        mode="cells",
        density_prior="rna_count_based",
        num_epochs=num_epochs,
        random_state=56,
        device="cpu",  # or: cuda
    )

    ad_ge = tg.project_genes(adata_map=ad_map, adata_sc=adata_sc)

    print(ad_ge)

    return ad_ge


def raw_impute(adata_sc, adata_st, output_folder: str, file_name: str):
    adata_st = adata_st[adata_st.obs['predicted.celltype'].isin(['Astrocyte', 'Microglia', 'Oligodendrocyte'])]

    sc.pp.normalize_total(adata_st)
    sc.pp.log1p(adata_st)

    adata_st.var["mt"] = adata_st.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata_st, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    adata_st = adata_st[adata_st.obs.pct_counts_mt < 5, :].copy()

    ad_ge = dup(adata_sc, adata_st)

    ad_ge.write(
        Path(os.path.join(output_folder, file_name+".h5ad")),
        compression="gzip")
