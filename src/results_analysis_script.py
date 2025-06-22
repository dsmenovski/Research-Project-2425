import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from matplotlib import pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
from torch import layout

from model_classes import *

parser = argparse.ArgumentParser(description="")

parser.add_argument("-i", "--input", required=True, type=str, help="path to folder with classified Xenium data")
parser.add_argument("-d", "--dataset", required=True, type=str, help="name of model dataset (ROSMAP / SEA-AD)")
parser.add_argument("-o", "--output", required=True, type=str, help="path to output model folder")
parser.add_argument("-c", "--celltype", required=False, type=str, help="a-/m-/o-", default="")

args = parser.parse_args()

modelLC = LinearClassificationModel()
modelMLP = MultilayerPerceptron()
modelDROP = DropoutMultilayerPerceptron()

models = [modelLC, modelMLP, modelDROP]
model_names = ["Linear Classification", "Multilayer Perceptron", "MLP + Dropout"]
short_names = ["LC", "MLP", "MLP+Dropout"]

celltype = ''
if args.celltype != '':
    celltype = args.celltype + '_'

data = []
for model in models:
    adata = sc.read_h5ad(
        os.path.join(str(args.input),
                     celltype + args.dataset.lower() + '_' + model.shortstr + '_labeled_spatial_data.h5ad'))
    data.append(adata.copy())
    del adata

xs, ys = [], []
for i, adata in enumerate(data):
    x_dists = adata.obs[['plaque_distance_micron']].values.flatten()
    y_label = adata.obs[['y']].values.flatten()
    xs.append(x_dists)
    ys.append(y_label)

df = pd.DataFrame({
    'plaque_distance_micron': np.concatenate(xs),
    'y': np.concatenate(ys).astype(str),
    'model': np.concatenate([
        np.repeat(name, len(arr))
        for name, arr in zip(model_names, xs)
    ])
})

if args.celltype != '':
    celltype = args.celltype + ' '

plt.style.use('default')
sns.set_theme(style=None, palette=plt.rcParams['axes.prop_cycle'].by_key()['color'], font_scale=1.2)
fig, axes = plt.subplots(figsize=(10, 8))
sns.violinplot(data=df, y='model', x='plaque_distance_micron', ax=axes,  hue='y', hue_order=['1', '0'], split=True,
               orient='h', inner='quart')
handles, _ = axes.get_legend_handles_labels()
axes.legend(handles, ['AD', 'CT'], title='Label')
plt.xlim(right=1000)
plt.tight_layout(rect=(0., 0., 1., 0.975))
plt.title(f"Violin Plots of Labels to Distance to Plaque ({celltype + args.dataset})")
plt.xlabel("Distance to Plaque in Microns")
plt.ylabel("Normalized Distribution of AD and CT Labels Per Model")
plt.savefig(os.path.join(args.output, celltype.replace(' ', '_') + args.dataset.lower() + '_' + "multi_violin_plot.png"))
plt.close()

plt.style.use('default')
sns.set_theme(style=None, palette=plt.rcParams['axes.prop_cycle'].by_key()['color'], font_scale=1)

if args.celltype == '':
    celltype = 'All Type'
else:
    celltype = str(args.celltype).capitalize()

# Map of cells to distance
plt.figure()  # figsize=(10, 8))
plt.scatter(data[0].obs['cell_centroid_x'], data[0].obs['cell_centroid_y'], s=4, alpha=0.5,
            c=np.log1p(data[0].obs['plaque_distance_micron']), cmap='RdBu')
plt.colorbar(label='Distance to nearest amyloid deposit (log1p(micron))')
plt.tight_layout(rect=(0., 0., 1., 0.945))
plt.title(f'{celltype} Cell Locations Within Tissue \n Colored By Distance to Plaque')  # (alpha=0.5)')

if args.celltype != '':
    celltype = args.celltype + ' '
else:
    celltype = ''

plt.axis('equal')
plt.savefig(os.path.join(args.output, celltype.replace(' ', '_') + args.dataset.lower() + '_' + "map_cells_tissue.png"))
plt.close()

# Map of cells to probability
subplot_width = 4
# fig, axes = plt.subplots(1, 3, figsize=(subplot_width*3, subplot_width), layout='constrained',
#                          subplot_kw={'box_aspect': 1.0})
#
# global_ymin = min(adata.obs['y'].min() for adata in data)
# global_ymax = max(adata.obs['y'].max() for adata in data)

alltypes = 'All Cells' if args.celltype == '' else celltype.capitalize().replace(' ', '')

for i, adata in enumerate(data):
    plt.figure()  # figsize=(10, 8))
    scatter = plt.scatter(adata.obs['cell_centroid_x'], adata.obs['cell_centroid_y'], s=4,
                          alpha=0.5, c=np.log1p(adata.obs['y_soft']), cmap='RdBu_r', vmin=0, vmax=1)
    # axes[i].set_title(f'{model_names[i]}')
    # axes[i].set_aspect('equal')

    cbar = plt.colorbar(scatter)
    cbar.set_label('Probability cell is AD')
    plt.tight_layout(rect=(0., 0., 1., 0.95))
    plt.title(f'Tissue Map Colored by AD Probability Using \n {args.dataset} {short_names[i]} Model on {alltypes} Data')
    plt.axis('equal')
    plt.savefig(os.path.join(args.output, celltype.replace(' ', '_') + args.dataset.lower() + '_' +
                             short_names[i].lower() + "_map_cells_tissue_probability.png"))
    plt.close()

spearman = pd.DataFrame(index=model_names, columns=['corr', 'p-value'])
mannwhitney = pd.DataFrame(index=model_names, columns=['stat', 'p-value'])

for adata, name in zip(data, model_names):
    corr, p_value = spearmanr(adata.obs['plaque_distance_micron'], adata.obs['y_soft'], alternative='less')
    spearman.loc[name] = [corr, p_value]

    ad_full = adata.obs['plaque_distance_micron'][(adata.obs['y'] == 1)]
    ct_full = adata.obs['plaque_distance_micron'][(adata.obs['y'] == 0)]
    if not len(ad_full.to_numpy()) == 0 and not len(ct_full.to_numpy()) == 0:
        stat, p_value = mannwhitneyu(ad_full.to_numpy(), ct_full.to_numpy(), alternative="less")
        mannwhitney.loc[name] = [stat, p_value]

if args.celltype != '':
    celltype = args.celltype + '_'

spearman.to_csv(os.path.join(str(args.output), celltype + args.dataset.lower() + '_spearman_distance_probability' + '.csv'))
mannwhitney.to_csv(os.path.join(str(args.output), celltype + args.dataset.lower() + '_mannwhitney' + '.csv'))
