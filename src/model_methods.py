import os
import torch
import scanpy as sc
from anndata.experimental import AnnLoader
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split


# Code written with the help of the following tutorials:
# https://www.geeksforgeeks.org/linear-regression-using-pytorch/
# https://anndata.readthedocs.io/en/latest/tutorials/notebooks/annloader.html

def train(dataloader, model, criterion, optimizer):
    epoch_loss = 0.0
    total_correct = 0
    model.train()
    for batch in dataloader:
        x = batch.X.float()
        y = batch.obs['y'].long().unsqueeze(1).squeeze(1)

        y_pred = model(x)
        loss = criterion(y_pred, y)

        optimizer.zero_grad()
        loss.backward()

        optimizer.step()

        with torch.no_grad():
            y_binp = torch.argmax(y_pred, dim=1)
            correct = (y_binp == y).sum().item()
            total_correct += correct

        epoch_loss += loss.item()

    normalizer_train = len(dataloader.dataset)
    training_loss = epoch_loss / normalizer_train
    accuracy = total_correct / normalizer_train
    return training_loss, accuracy


def evaluate(dataloader, model, criterion):
    loss = 0.0
    total_correct = 0
    batch_accuracies = []
    model.eval()
    with torch.no_grad():
        for batch in dataloader:
            x = batch.X.float()
            y = batch.obs['y'].long().unsqueeze(1).squeeze(1)

            y_pred = model(x)
            loss += criterion(y_pred, y).item()

            y_binp = torch.argmax(y_pred, dim=1)
            correct = (y_binp == y).sum().item()
            batch_accuracies.append(correct / len(batch))
            total_correct += correct

    avg_loss = loss / len(dataloader.dataset)
    accuracy = total_correct / len(dataloader.dataset)
    return avg_loss, accuracy, batch_accuracies


def run_train(model, dataset: str, adata, output_folder, NUM_EPOCHS=21):
    criterion = model.criterion
    optimizer = torch.optim.Adam(model.parameters(), lr=0.0001)  # , weight_decay=1e-3)

    test_dataset = 'sea-ad' if dataset.lower() == 'rosmap' else 'rosmap'
    if test_dataset == 'sea-ad':
        tdata = sc.read_h5ad(os.path.join(output_folder, test_dataset + "_top_1k_genes_data.h5ad"))
        tdata.obs['y'] = [0 if x is False else 1 for x in
                          tdata.obs['Overall AD neuropathological Change'] == 'High']
    else:
        tdata = sc.read_h5ad(os.path.join(output_folder, test_dataset + "_filtered_1k_genes_data.h5ad"))

    # idColName, train_ids, test_ids = split_data(adata, dataset)
    # adata_train = adata[adata.obs[idColName].isin(train_ids)].copy()
    # adata_test = adata[adata.obs[idColName].isin(test_ids)].copy()
    adata_train = adata

    use_cuda = torch.cuda.is_available()

    train_loader = AnnLoader([adata_train], batch_size=128, shuffle=True, use_cuda=use_cuda)

    # losses = []
    accuracies = []
    eval_accuracies = []
    test_accuracies = []
    for epoch in range(NUM_EPOCHS):
        total_loss, epoch_acc = train(train_loader, model, criterion, optimizer)
        # losses.append(total_loss)
        accuracies.append(epoch_acc)
        if epoch % 5 == 0 or epoch == NUM_EPOCHS - 1:
            print("[epoch %03d]  average training loss: %.4f / accuracy: %.4f" % (epoch, total_loss, epoch_acc))
        avg_loss, avg_acc, _ = run_eval(model, test_dataset, output_folder, tdata)
        # test_loss, test_acc, _ = run_eval(model, dataset, output_folder, adata_test)
        eval_accuracies.append(avg_acc)
        # test_accuracies.append(test_acc)

    torch.save(model.state_dict(), os.path.join(output_folder, dataset + "_" + str(model) + ".pth"))

    return accuracies, eval_accuracies, test_accuracies


def run_eval(model, dataset: str, model_folder, adata):  # data_file="pca_labeled_1k_filtered_v2.h5ad"):
    criterion = model.criterion

    adata_test = adata

    use_cuda = torch.cuda.is_available()
    test_loader = AnnLoader([adata_test], batch_size=128, shuffle=True, use_cuda=use_cuda)

    test_loss, test_acc, accuracies = evaluate(test_loader, model, criterion)
    # print(accuracies)
    print(f"({dataset}) avg test loss: %.4f / avg test accuracy: %.4f" % (test_loss, test_acc))

    return test_loss, test_acc, accuracies


def split_data(adata, dataset: str):
    idColName = 'individualID' if dataset == 'ROSMAP' else 'Donor ID' if dataset == 'SEA-AD' else 'Unsupported dataset'
    if idColName == 'Unsupported dataset':
        raise ValueError('Unsupported dataset')

    train_ids, test_ids = train_test_split(
        adata.obs[idColName].unique(),
        test_size=0.2,
        random_state=56
    )
    return idColName, train_ids, test_ids


def plot_training_accuracy(model_accuracies, models, dataset: str, output_folder):
    plt.figure()
    for i, model in enumerate(models, start=1):
        model_name, linestyle = ("LC", '-') if i == 1 else ("MLP", '--') if i == 2 else ("DROP", ':')
        label = model_name if i != 3 else "MLP+Dropout"
        plt.plot(model_accuracies[model_name], label=label, linestyle=linestyle)
    plt.legend()
    plt.title(dataset + " " + "Model Training Accuracy")
    plt.xlabel("Epoch")
    plt.ylabel("Accuracy")
    plt.savefig(output_folder + "/" + dataset + "_" + "full_training_accuracy.png")
    plt.close()


def plot_testing_accuracy(model_accuracies, models, dataset: str, evaluation: bool, model_folder):
    plt.figure()
    for i, model in enumerate(models, start=1):
        model_name, linestyle = ("LC", '-') if i == 1 else ("MLP", '--') if i == 2 else ("DROP", ':')
        label = model_name if i != 3 else "MLP+Dropout"
        plt.plot(model_accuracies[model_name], label=label, linestyle=linestyle)
    plt.legend()
    typeAcc = 'Test'
    if evaluation:
        typeAcc = 'Evaluation'
    plt.title(dataset + " " + f"Model {typeAcc} Accuracy")
    plt.xlabel("Epoch")
    plt.ylabel("Accuracy")
    typeAcc = 'testing'
    if evaluation:
        typeAcc = 'evaluation'
    plt.savefig(model_folder + "/" + dataset + "_" + f"full_{typeAcc}_accuracy.png")
    plt.close()
