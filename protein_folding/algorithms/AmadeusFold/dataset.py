from torch.utils.data import Dataset
import pandas as pd
import ast

class FoldDataset(Dataset):
	def __init__(self, csv_location='./data/test_data.csv', shuffle=True):
		self.dataframe = pd.read_csv(csv_location)
		if shuffle:
			self.dataframe = self.dataframe.sample(frac = 1)
	
	def __len__(self):
		return len(self.dataframe)
	
	def __getitem__(self, idx):
		datapoint = self.dataframe.iloc[idx]
		sequence = datapoint['sequence']
		coordinates = ast.literal_eval(datapoint['coordinates'])
		return sequence, coordinates


if __name__ == '__main__':
	dataset = FoldDataset()
	print(dataset[5])