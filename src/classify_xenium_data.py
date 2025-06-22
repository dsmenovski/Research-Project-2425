import os
from pathlib import Path

import torch
import numpy as np
import scanpy as sc
import torch.nn.functional as f
from matplotlib import pyplot as plt

from model_classes import *


def classify(model: Model, model_dataset: str, celltype: str, input_folder: str, output_folder: str):
    model.load_state_dict(torch.load(
        os.path.join(input_folder, model_dataset + "_" + str(model) + ".pth"), weights_only=True))

    if celltype != '':
        celltype = celltype + '_'

    if model_dataset.lower() == "sea-ad":
        ad_imp = sc.read_h5ad(os.path.join(input_folder, "rosmap_"+celltype+"imputed_xenium_1k.h5ad"))
    elif model_dataset.lower() == "rosmap":
        ad_imp = sc.read_h5ad(os.path.join(input_folder, "sea-ad_imputed_xenium_1k.h5ad"))
    else:
        raise IOError("Unsupported dataset")

    with torch.no_grad():
        x_tensor = torch.tensor(ad_imp.X, dtype=torch.float32)
        ad_pred = model(x_tensor)

    y_pred = torch.argmax(ad_pred, dim=1)
    ad_imp.obs['y'] = np.array(y_pred)

    print("AD: %d, total: %d" % (sum(y_pred), len(y_pred)))

    y_soft = f.softmax(ad_pred, dim=1)[:, 0]
    ad_imp.obs['y_soft'] = np.array(y_soft)

    y_total = sum(y_soft)
    y_mean = y_total / len(y_soft)
    print("total: %.2f, mean: %.2f" % (y_total, y_mean))

    ad_imp.write(
        Path(os.path.join(output_folder,
                          celltype + model_dataset.lower() + "_" + model.shortstr + "_labeled_spatial_data.h5ad")),
        compression="gzip")

    return y_soft.tolist(), y_pred.tolist()


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct * total / 100.0))
        return f'{val} ({pct:.1f}%)'

    return my_autopct


def plotBarChart(y_preds, dataset, model_names, output_folder):
    plt.figure()
    plt.style.use('tableau-colorblind10')

    yAD, yCT = [], []
    for y_pred in y_preds:
        yCT.append(y_pred.count(0))
        yAD.append(y_pred.count(1))

    plt.bar(x=model_names, height=yAD)
    plt.bar(x=model_names, height=yCT, bottom=yAD)
    plt.xlabel('Model')
    plt.ylabel('Count')
    plt.legend(['AD', 'CT'])
    plt.title(f'Predicted Labels Distribution of {dataset} Models')
    plt.savefig(output_folder + '/' + dataset + '_bar_chart.png')
    plt.close()


def plotPieChart(y_preds, dataset, model_names, output_folder):
    plt.figure()
    colors = [
        ['#a6cee3', '#1f78b4'],
        ['#b2df8a', '#33a02c'],
        ['#fdbf6f', '#ff7f00']
    ]

    sizes = []
    labels = []
    wedge_colors = []

    for i, y_pred in enumerate(y_preds):
        count_0 = y_pred.count(0)
        count_1 = y_pred.count(1)
        total = count_0 + count_1
        perc_0 = count_0 / total
        perc_1 = count_1 / total

        total_section = 360 / 3
        sizes.append(total_section * perc_0)
        sizes.append(total_section * perc_1)

        labels.append(f'{model_names[i]} - CT\n{count_0} ({perc_0 * 100:.1f}%)')
        labels.append(f'{model_names[i]} - AD\n{count_1} ({perc_1 * 100:.1f}%)')

        wedge_colors.extend(colors[i])

    fig, ax = plt.subplots(figsize=(8, 8))
    _, _ = ax.pie(
        sizes,
        labels=labels,
        colors=wedge_colors,
        startangle=90,
        wedgeprops=dict(edgecolor='white')
    )

    ax.set(aspect='equal')
    plt.tight_layout()
    plt.title(dataset + ' Label Distribution of All Models')
    plt.savefig(output_folder + '/' + dataset + '_pie_chart.png')
    plt.close()


def plotHistogram(data, bins, dataset, output_folder):
    plt.figure()
    plt.style.use('default')
    plt.hist(data[0], label='LR', alpha=0.5, bins=bins)
    plt.hist(data[1], label='MLP', alpha=0.5, bins=bins)
    plt.hist(data[2], label='MLP+Dropout', alpha=0.5, bins=bins)
    plt.legend()
    plt.title(dataset + ' AD Prediction Distribution of All Models')
    plt.xlabel('Probability of AD')
    plt.ylabel('Amount of Predictions')
    plt.savefig(output_folder + '/' + dataset + '_histogram.png')
    plt.close()
