import argparse
import random

import scanpy as sc
from model_classes import *
from model_methods import run_train, run_eval
from pre_processing import *
from model_methods import plot_training_accuracy, plot_testing_accuracy


def train_models(data, dataset: str, output: str, epochs: int):
    training_data = data

    # train
    modelLC = LinearClassificationModel()
    modelMLP = MultilayerPerceptron()
    modelDROP = DropoutMultilayerPerceptron()
    
    model_training_accuracies = dict()
    model_evaluation_accuracies = dict()
    model_testing_accuracies = dict()
    models = [modelLC, modelMLP, modelDROP]
    
    model_training_accuracies['LC'], model_evaluation_accuracies['LC'], model_testing_accuracies['LC'] = (
        run_train(modelLC, dataset, training_data, output, epochs))
    model_training_accuracies['MLP'], model_evaluation_accuracies['MLP'], model_testing_accuracies['MLP'] = (
        run_train(modelMLP, dataset, training_data, output, epochs))
    model_training_accuracies['DROP'], model_evaluation_accuracies['DROP'], model_testing_accuracies['DROP'] = (
        run_train(modelDROP, dataset, training_data, output, epochs))

    plot_training_accuracy(model_training_accuracies, models, dataset, output)

    plot_testing_accuracy(model_evaluation_accuracies, models, dataset, evaluation=True, model_folder=output)

    # plot_testing_accuracy(model_testing_accuracies, models, dataset, evaluation=False, model_folder=output)

    del training_data


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-i", "--input", required=True, type=str, help="path to prepared SEA-AD h5ad data")
    parser.add_argument("-r", "--rosmap", required=True, type=str, help="path to prepared ROSMAP h5ad data")
    parser.add_argument("-o", "--output", required=True, type=str, help="path to output model folder")
    parser.add_argument("-e", "--epochs", type=int, help="number of epochs")

    args = parser.parse_args()

    random.seed(56)
    np.random.seed(56)
    torch.manual_seed(56)

    rosmap_training_data = sc.read_h5ad(args.rosmap)

    train_models(rosmap_training_data, "ROSMAP", args.output, args.epochs)

    del rosmap_training_data

    random.seed(56)
    np.random.seed(56)
    torch.manual_seed(56)

    sea_ad_training_data = sc.read_h5ad(args.input)

    sea_ad_training_data.obs['y'] = [0 if x is False else 1 for x in
                                     sea_ad_training_data.obs['Overall AD neuropathological Change'] == 'High']

    train_models(sea_ad_training_data, "SEA-AD", args.output, args.epochs)

    del sea_ad_training_data
