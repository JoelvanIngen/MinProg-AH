import torch

from LSTM_model import AmadeusFold
from utils import index_sequence

def main():
	# settings
	model_path = './models/<your_model>.pt'
	sequence = 'HHPHHCCHPHCPP' # 'HHCHHCPCHCH'
	multiplier = 8 # check dataset settings

	# load model and weight
	print("Loading model...")
	model = AmadeusFold(n_molecules = 4)
	model.load_state_dict(torch.load(model_path))
	print("Model loaded.")

	# process and feed sequence to model
	sequence_indexed = torch.tensor(index_sequence(sequence))
	with torch.no_grad():
		model_output = model(sequence_indexed)
	coordinates = list(torch.round(model_output * multiplier))
	print(coordinates)

if __name__ == "__main__":
	main()
