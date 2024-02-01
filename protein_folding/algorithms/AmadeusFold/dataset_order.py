from torch.utils.data import Dataset
from torch.nn.functional import normalize
from torchvision.transforms import ToTensor
import torch
import numpy as np
import pandas as pd
import ast

from utils import index_sequence


order_index = {
	0: 0,
	1: 1,
	-1: 2,
	2: 3,
	-2: 4,
	3: 5,
	-3: 6
}


class FoldDatasetOrder(Dataset):
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
			n_dimensions = 2
			) -> None:
		self.normalize = normalize
		self.n_dimensions = n_dimensions
		
		# read in .csv file
		self.dataframe = pd.read_csv(csv_location)
		if shuffle:
			self.dataframe = self.dataframe.sample(frac = 1)

	def __len__(self) -> int:
		return len(self.dataframe)
	
	def __getitem__(self, idx: int) -> tuple:
		datapoint = self.dataframe.iloc[idx]
		sequence = index_sequence(datapoint['sequence'])[0]
		order = [0] + list(ast.literal_eval(datapoint['order']))
		order_indexed = [order_index[direction] for direction in order]
		score = datapoint['score']

		sequence = torch.tensor(sequence)
		order = torch.tensor(order_indexed)
		score = torch.tensor(score)
	
		return sequence, order, score


if __name__ == '__main__':
	dataset = FoldDataset()
	print(dataset[5])