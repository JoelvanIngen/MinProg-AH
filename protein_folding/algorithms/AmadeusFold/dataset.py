from torch.utils.data import Dataset
from torch.nn.functional import normalize
from torchvision.transforms import ToTensor
import numpy as np
import pandas as pd
import ast

class FoldDataset(Dataset):
	def __init__(
			self,
			csv_location='./data/test_data.csv',
			shuffle=True,
			transform = ToTensor(),
			normalize = True
			) -> None:
		self.dataframe = pd.read_csv(csv_location)
		self.transform = transform
		self.normalize = normalize
		if shuffle:
			self.dataframe = self.dataframe.sample(frac = 1)
	
	def __len__(self) -> int:
		return len(self.dataframe)
	
	def __getitem__(self, idx: int) -> tuple:
		datapoint = self.dataframe.iloc[idx]
		sequence = datapoint['sequence']
		coordinates = np.asarray(ast.literal_eval(datapoint['coordinates']))

		if self.transform:
			coordinates = self.transform(coordinates)
			if self.normalize:
				coordinates = normalize(coordinates, dim = 1)
				
		return sequence, coordinates


if __name__ == '__main__':
	dataset = FoldDataset()
	print(dataset[5])