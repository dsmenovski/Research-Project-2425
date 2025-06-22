import argparse
import random

from matplotlib import pyplot as plt
from model_classes import *
from classify_xenium_data import *


def classification(model_dataset, models, celltype, model_names, input, output):

    random.seed(56)
    np.random.seed(56)
    torch.manual_seed(56)

    ads, ys = [], []
    for model in models:
        ad_pred, y_pred = classify(model, model_dataset, celltype, input, output)
        ads.append(ad_pred)
        ys.append(y_pred)

    if celltype != '':
        celltype = celltype + ' '

    plotBarChart(ys, celltype.capitalize()+model_dataset, model_names, output)

    plotHistogram(ads, 15, celltype.capitalize()+model_dataset, output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--input", required=True, type=str, help="path to input folder")
    parser.add_argument("-o", "--output", required=True, type=str, help="path to output folder")

    args = parser.parse_args()

    modelLR = LinearClassificationModel()
    modelMLP = MultilayerPerceptron()
    modelDROP = DropoutMultilayerPerceptron()

    model_results = dict()
    models = [modelLR, modelMLP, modelDROP]
    model_names = ["Linear Classification", "Multilayer Perceptron", "MLP + Dropout"]

    # classification("SEA-AD", models, 'astrocytes', model_names, args.input, args.output)
    # classification("SEA-AD", models, 'microglia', model_names, args.input, args.output)
    classification("SEA-AD", models, 'oligodendroglia', model_names, args.input, args.output)
    classification("SEA-AD", models, '', model_names, args.input, args.output)
    classification("ROSMAP", models, '', model_names, args.input, args.output)
