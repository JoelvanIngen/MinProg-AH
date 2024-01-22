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
			transform = ToTensor(),
			normalize = True
			) -> None:
		self.transform = transform
		self.normalize = normalize
		
		# read in .csv file
		self.dataframe = pd.read_csv(csv_location)
		if shuffle:
			self.dataframe = self.dataframe.sample(frac = 1)

		# find maximum coordinate value in any direction
		self.max_coordinate_value = 0
		for idx in range(len(self.dataframe)):
			datapoint = self.dataframe.iloc[idx]
			coordinates_list = ast.literal_eval(datapoint['coordinates'])
			for coordinates in coordinates_list:
				for coordinate in coordinates:
					if abs(coordinate) > self.max_coordinate_value:
						self.max_coordinate_value = abs(coordinate)

	def __len__(self) -> int:
		return len(self.dataframe)
	
	def __getitem__(self, idx: int) -> tuple:
		datapoint = self.dataframe.iloc[idx]
		sequence = datapoint['sequence']
		coordinates = np.asarray(ast.literal_eval(datapoint['coordinates']))
		score = datapoint['score']

		if self.transform:
			coordinates = self.transform(coordinates)
			score = torch.tensor(score)
		# divide by maximum found value if normalized
		if self.normalize:
			coordinates/= self.max_coordinate_value
				
		return sequence, coordinates, score, self.max_coordinate_value


if __name__ == '__main__':
	dataset = FoldDataset()
	print(dataset[5])