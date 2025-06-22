import pandas as pd
import numpy as np
import random
import string
from matplotlib import pyplot as plt


def generate_random_label(existing_keys, length=6):
    while True:
        label = ''.join(random.choices(string.ascii_lowercase, k=length))
        if label not in existing_keys:
            return label


if __name__ == '__main__':

    # Insert key-value pairs of genes to correlations from the imputation_validation.py code csv file result into the
    # dictionary below.
    # E.g.
    # data = {'gene1': 0.5, 'gene2': 0.7}
    data = {}

    data = {gene.upper(): corr for gene, corr in data.items()}

    series = pd.Series(data)
    series.sort_values(inplace=True)

    plt.bar(series.index, series.values)
    plt.xticks(rotation=90)
    plt.xlabel('Genes')
    plt.yticks(np.arange(0, 1.1, 0.1))
    plt.ylabel('Spearman Correlation Coefficient')
    plt.title('Correlation between measured and predicted gene expressions')
    plt.tight_layout()
    plt.show()
