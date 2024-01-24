import torch
import torch.nn as nn

from dataset import FoldDataset
from utils import compute_bond_score_coordinates, fold_validity_quantified_loss, StraightThroughEstimator


class FoldLoss(nn.Module):
	def __init__(
		self,
		rounding_scale = 1.0,
		validity_scale = 1.0,
		bond_qual_scale = 1.0
	):
		super().__init__()
		self.rounding_scale = rounding_scale
		self.validity_scale = validity_scale
		self.bond_qual_scale = bond_qual_scale

	def forward(
		self, 
		predictions, 
		norm_factor, 
		sequence,
		targets, 
		target_score,
	):
		# un-normalise predictions and targets
		predictions *= norm_factor
		targets *= norm_factor

		# round output with straight-through estimator (differentiable)
		rounded_predictions = StraightThroughEstimator.apply(predictions)

		# determine rounding error, scale and add to loss
		rounding_error = torch.mean((predictions - rounded_predictions)**2)

		# determine validity, scale and add to loss
		validity = fold_validity_quantified_loss(rounded_predictions)

		# if presented ordering is valid compute score difference
		if validity == 0:
			prediction_score = compute_bond_score_coordinates(sequence, predictions)
			score_diff = torch.abs(prediction_score - target_score)
			# unless similar score is achieved, compare to target coordinates
			if score_diff != 0:
				score_diff += torch.abs(torch.sum(rounded_predictions - targets))
		# to avoid creating minima going from invalid orders to valid orders
		# with bad scores, add the target score as minimum score loss
		else:
			score_diff = -target_score

		# calculate total loss
		loss = (
			rounding_error * self.rounding_scale +
			validity * self.validity_scale +
			score_diff * self.bond_qual_scale
		)
			
		return loss


def main():
	# WORK IN PROGRESS
	dataset = FoldDataset()
	sequence = dataset[5][0]
	coordinates = dataset[5][1]

	sample_output = torch.tensor([[[ 0.0589, -0.1306, -0.0682],
								[ 0.0587, -0.0938, -0.0634],
								[ 0.0707, -0.1321, -0.0722],
								[ 0.0680, -0.1537, -0.0676],
								[ 0.0666, -0.1641, -0.0630]]])

	loss = FoldLoss()
	l = loss.forward(
		sample_output, 
		dataset.norm_val,
		dataset[5][0],
		dataset[5][1],
		dataset[5][2]
	)
	print(f"loss: {l}")
	
if __name__ == "__main__":
	main()