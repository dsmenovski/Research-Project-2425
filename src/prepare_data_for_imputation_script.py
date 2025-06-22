import argparse
import os
import scanpy as sc
import numpy as np

from pre_processing import prepare_sc_data_for_imputation

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--input", required=True, type=str, help="path to input folder")
parser.add_argument("-o", "--output", required=True, type=str, help="path to output folder")
parser.add_argument("-d", "--dataset", required=True, type=str, help="dataset to use")
parser.add_argument("--oligs", required=False, type=bool, help="is the input data ROSMAP oligodendrocytes")
parser.add_argument("--astro", required=False, type=bool, help="is the input data ROSMAP astrocytes")

args = parser.parse_args()

if args.dataset == 'ROSMAP':
    fdata = sc.read_h5ad(os.path.join(args.input, 'full_rosmap_data.h5ad'))

    # Downsample the data
    if args.oligs or args.astro:
        downsample_rate = 4 if args.oligs else 2
        n_cells = fdata.n_obs
        print(n_cells)
        np.random.seed(56)
        selected_pos = np.random.choice(fdata.n_obs, size=int(n_cells/downsample_rate), replace=False)
        selected_obs_names = fdata.obs_names[selected_pos]
        fdata = fdata[fdata.obs_names.isin(selected_obs_names)].copy()
        print(fdata.n_obs)

    _ = prepare_sc_data_for_imputation(fdata,
                                   os.path.join(str(args.output)),
                                   "rosmap_data_for_imputation_"+args.input)

if args.dataset == 'SEA-AD':
    sdata = sc.read_h5ad(os.path.join(args.input, 'combined_objects_sea-ad.h5ad'))

    _ = prepare_sc_data_for_imputation(sdata,
                                       os.path.join(str(args.output)),
                                       "sea-ad_data_for_imputation")
