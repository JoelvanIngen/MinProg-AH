from torch.utils.data import Dataset
from torch.nn.functional import normalize
from torchvision.transforms import ToTensor
import torch
import numpy as np
import pandas as pd
import ast

from utils import index_sequence


class FoldDataset(Dataset):
	"""
	Dataset class to accommodate protein folding data. To normalise, all
	coordinates are divided by the maximum distance from 0 found in any
	direction in the provided .csv.
	"""
	def __init__(
			self,
			csv_location='./data/test_data.csv',
			shuffle=True,
			normalize = True,
			norm_val = 8
			) -> None:
		self.normalize = normalize
		self.norm_val = norm_val
		
		# read in .csv file
		self.dataframe = pd.read_csv(csv_location)
		if shuffle:
			self.dataframe = self.dataframe.sample(frac = 1)

	def __len__(self) -> int:
		return len(self.dataframe)
	
	def __getitem__(self, idx: int) -> tuple:
		datapoint = self.dataframe.iloc[idx]
		sequence = index_sequence(datapoint['sequence'])[0]
		coordinates = np.asarray(ast.literal_eval(datapoint['coordinates']))
		# divide by maximum found value if normalized
		if self.normalize:
			coordinates/= self.norm_val
		score = datapoint['score']

		sequence = torch.tensor(sequence)
		coordinates = torch.tensor(coordinates, dtype=torch.float32, requires_grad=True) #self.transform(coordinates, dtype=torch.float32)
		score = torch.tensor(score)
	
				
		return sequence, coordinates, score


if __name__ == '__main__':
	dataset = FoldDataset()
	print(dataset[5])