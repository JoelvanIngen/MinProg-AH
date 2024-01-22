import torch
import torch.nn as nn

from dataset import FoldDataset
from utils import compute_bond_score_coordinates, fold_validity_quantified_loss


class FoldLoss(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, predictions, targets):
		# TODO: 
		# - add invalid loss
		# - add loss if valid but less good bond than optimal
		# - add loss of difference between float and rounded int coordinates (ie 1.25 is a bad output, should be 1.)

        return loss


def main():
	# WORK IN PROGRESS
	dataset = FoldDataset()
	sequence = dataset[5][0]
	coordinates = dataset[5][1] * dataset[5][3]

	print(compute_bond_score_coordinates(sequence, coordinates), dataset[5][2])

	sample_output = torch.tensor([[[ 0.0589, -0.1306, -0.0682],
								[ 0.0587, -0.0938, -0.0634],
								[ 0.0707, -0.1321, -0.0722],
								[ 0.0680, -0.1537, -0.0676],
								[ 0.0666, -0.1641, -0.0630]]])

	print(fold_validity_quantified_loss(torch.round(sample_output)))
	print(fold_validity_quantified_loss(coordinates))
	
if __name__ == "__main__":
	main()