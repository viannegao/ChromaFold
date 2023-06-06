import numpy as np

import torch
import torch.nn as nn


class resblock(nn.Module):
    def __init__(self, ni):
        super(resblock, self).__init__()
        self.blocks = nn.Sequential(
            nn.Conv1d(ni, ni, 3, 1, 1),
            nn.BatchNorm1d(ni),
            nn.ReLU(),
            nn.Conv1d(ni, ni, 3, 1, 1),
            nn.BatchNorm1d(ni),
            nn.ReLU(),
        )

    def forward(self, x):
        residual = x
        out = self.blocks(x)
        out = out + residual
        return out


class symmetrize_bulk(nn.Module):
    def __init__(self):
        super(symmetrize_bulk, self).__init__()

    def forward(self, x):
        if len(x.shape) == 2:
            print("not implemented")
            return None
        else:
            if len(x.shape) == 3:
                a, b, c = x.shape
                x = x.reshape(a, b, 1, c)
                x = x.repeat(1, 1, c, 1)
                x_t = x.permute(0, 1, 3, 2)
                x_sym = torch.concat((x, x_t), axis=1)  # (x+x_t)/2
                return x_sym
            else:
                return None


# Pseudo-bulk branch
class branch_pbulk(nn.Module):
    def __init__(self):
        super(branch_pbulk, self).__init__()

        scatac_res = 50

        self.bulk_summed_2d = nn.Sequential(
            nn.AvgPool1d(kernel_size=np.int64(1e04 / scatac_res)), symmetrize_bulk()
        )

        self.bulk_extractor_2d = nn.Sequential(
            nn.Conv1d(
                in_channels=2,
                out_channels=16,
                kernel_size=11,
                stride=1,
                dilation=1,
                padding="same",
            ),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2),
            nn.Conv1d(
                in_channels=16,
                out_channels=32,
                kernel_size=7,
                stride=1,
                dilation=1,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2),
            nn.Conv1d(
                in_channels=32,
                out_channels=32,
                kernel_size=5,
                stride=1,
                dilation=1,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=2),
            nn.Conv1d(
                in_channels=32,
                out_channels=32,
                kernel_size=5,
                stride=1,
                dilation=1,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            #             nn.MaxPool1d(kernel_size = 2),
            nn.Conv1d(
                in_channels=32,
                out_channels=32,
                kernel_size=5,
                stride=1,
                dilation=1,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            #             nn.MaxPool1d(kernel_size = 2),
            nn.Conv1d(
                in_channels=32,
                out_channels=32,
                kernel_size=5,
                stride=1,
                dilation=2,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            #             nn.MaxPool1d(kernel_size = 2),
            nn.Conv1d(
                in_channels=32,
                out_channels=32,
                kernel_size=5,
                stride=1,
                dilation=3,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            #             nn.MaxPool1d(kernel_size = 2),
            nn.Conv1d(
                in_channels=32,
                out_channels=32,
                kernel_size=5,
                stride=1,
                dilation=5,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(
                in_channels=32,
                out_channels=32,
                kernel_size=5,
                stride=1,
                dilation=5,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(
                in_channels=32,
                out_channels=32,
                kernel_size=5,
                stride=1,
                dilation=7,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(
                in_channels=32,
                out_channels=32,
                kernel_size=5,
                stride=1,
                dilation=11,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.Conv1d(
                in_channels=32,
                out_channels=32,
                kernel_size=5,
                stride=1,
                dilation=11,
                padding="same",
            ),
            nn.BatchNorm1d(32),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=5),
            nn.Conv1d(
                in_channels=32,
                out_channels=16,
                kernel_size=3,
                stride=1,
                dilation=1,
                padding="same",
            ),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.MaxPool1d(kernel_size=5),
            nn.Conv1d(
                in_channels=16,
                out_channels=16,
                kernel_size=3,
                stride=1,
                dilation=1,
                padding="same",
            ),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            nn.Conv1d(
                in_channels=16,
                out_channels=16,
                kernel_size=3,
                stride=1,
                dilation=1,
                padding="same",
            ),
            nn.BatchNorm1d(16),
            nn.ReLU(),
            symmetrize_bulk(),
        )

        self.total_extractor_2d = nn.Sequential(
            nn.Conv2d(in_channels=36, out_channels=64, kernel_size=3, stride=2),
            nn.BatchNorm2d(64),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2),
            nn.Conv2d(in_channels=64, out_channels=32, kernel_size=3, stride=2),
            nn.BatchNorm2d(32),
            nn.ReLU(),
            nn.MaxPool2d(kernel_size=2),
            nn.Conv2d(in_channels=32, out_channels=16, kernel_size=3, stride=2),
            nn.BatchNorm2d(16),
            nn.ReLU(),
        )

        self.classifier = nn.Sequential(
            nn.Linear(in_features=(1936), out_features=512),
        )
        self.classifier2 = nn.Sequential(nn.Linear(in_features=(512), out_features=200))

    def forward(self, x2):

        x3_2d = self.bulk_summed_2d(x2)
        x2_2d = self.bulk_extractor_2d(x2)

        x4 = torch.cat((x3_2d, x2_2d), 1)
        x4 = self.total_extractor_2d(x4)
        x4 = torch.flatten(x4, 1)
        x4 = self.classifier(x4)
        out = self.classifier2(x4)

        return out
