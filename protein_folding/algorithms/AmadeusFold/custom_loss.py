import torch
import torch.nn as nn

from dataset import FoldDataset
from utils import compute_bond_score_coordinates, fold_validity_quantified_loss, StraightThroughEstimator


class FoldLoss(nn.Module):
	"""
	A custom loss function to create a general metric for the quality of a fold
	prediction model. There are three main components:

	Rounding error: the model outputs floats, but the grid is discrete. The
		difference between the nearest int and a float is added to total
		loss.

	Validity: a function has been created to quantify the "invalidness" of a
		model output, to allow for gradient descent. The more factors that make
		an output invalid (improper distancing, duplicate coordinates), the
		higher this score is.

	Score difference: if a valid fold has a worse score than the optimal score
		in the dataset, a loss is added that scales with the difference between
		the prediction and the order of the datapoint
	"""
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
			score_diff = (prediction_score - target_score).float()
			# if lower score is found, the dataset does not contain optimal folds
			# and should not be used
			if score_diff < 0:
				print(f"NOTICE: a score {prediction_score} has been found for a sequence that is lower than the optimal score {target_score} in the dataset. Consider the validity of your data!")
			# unless similar score is achieved, compare to target coordinates
			elif score_diff > 0:
				score_diff += torch.mean(torch.abs(rounded_predictions - targets))
		# to avoid creating minima going from invalid orders to valid orders
		# with bad scores, add the target score as maximum score loss
		else:
			score_diff = -target_score * len(sequence)
	
		# calculate total loss
		loss = (
			rounding_error * self.rounding_scale +
			validity * self.validity_scale +
			score_diff * self.bond_qual_scale
		)
			
		return loss


def main():
	dataset = FoldDataset(shuffle=False)
	sequence = dataset[5][0]
	coordinates = dataset[5][1]

	sample_output = torch.tensor([[[ 0.0589, -0.1306, -0.0682],
								[ 0.0587, -0.0938, -0.0634],
								[ 0.0707, -0.1321, -0.0722],
								[ 0.0680, -0.1537, -0.0676],
								[ 0.0666, -0.1641, -0.0630]]])

	sample_output_validity = torch.tensor([[[ 0.0000,  0.0000,  0.0000],
								[-0.0667,  0.0000,  0.0000],
								[-0.1333,  0.0000,  0.0000],
								[-0.1333, -0.0667,  0.0000],
								[-0.0667, -0.0667,  0.0000],
								[-0.0000, -0.0667,  0.0000],
								#[-0.0000, -0.0667,  0.0667],
								#[-0.1333, -0.0667,  0.0000],
								#[-0.1333, -0.0000,  0.0000],
								]])	

	sample_output_score = torch.tensor([[[ 0.0000,  0.0000,  0.0000],
								[-0.0667,  0.0000,  0.0000],
								[-0.1333,  0.0000,  0.0000],
								[-0.1333, -0.0667,  0.0000],
								[-0.2, -0.0667,  0.0000],
								[-0.2667, -0.0667,  0.0000],
								#[-0.0000, -0.0667,  0.0667],
								#[-0.1333, -0.0667,  0.0000],
								#[-0.1333, -0.0000,  0.0000],
								]])	

	loss = FoldLoss()
	l = loss.forward(
		sample_output, 
		dataset.norm_val,
		dataset[5][0],
		dataset[5][1],
		dataset[5][2]
	)
	print(f"loss: {l}\nTarget: {coordinates}\nSample output: {sample_output}")
	
if __name__ == "__main__":
	main()