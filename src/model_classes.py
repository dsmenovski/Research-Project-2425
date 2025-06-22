from abc import abstractmethod

import torch


class Model(torch.nn.Module):
    def __init__(self, shortstr):
        super().__init__()
        self._shortstr = shortstr

    @property
    def shortstr(self):
        return self._shortstr


class LinearClassificationModel(Model):
    def __init__(self):
        super().__init__("lc")
        self.linear = torch.nn.Linear(1000, 2)
        self.criterion = torch.nn.CrossEntropyLoss(size_average=False)

    def forward(self, x):
        y_pred = self.linear(x)
        return y_pred

    def __str__(self):
        return "linear_classification"


class DropoutMultilayerPerceptron(Model):
    def __init__(self):
        super().__init__("dmlp")
        self.linear_relu_stack = torch.nn.Sequential(
            torch.nn.Linear(1000, 512),
            torch.nn.ReLU(),
            torch.nn.Dropout(0.3),
            torch.nn.Linear(512, 2)
        )
        self.criterion = torch.nn.CrossEntropyLoss(size_average=False)

    def forward(self, x):
        y_pred = self.linear_relu_stack(x)
        return y_pred

    def __str__(self):
        return "dropout_multilayer_perceptron"


class MultilayerPerceptron(Model):
    def __init__(self):
        super().__init__("mlp")
        self.linear_relu_stack = torch.nn.Sequential(
            torch.nn.Linear(1000, 512),
            torch.nn.ReLU(),
            torch.nn.Linear(512, 2)
        )
        self.criterion = torch.nn.CrossEntropyLoss(size_average=False)

    def forward(self, x):
        y_pred = self.linear_relu_stack(x)
        return y_pred

    def __str__(self):
        return "multilayer_perceptron"
