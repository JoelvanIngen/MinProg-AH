import torch

from LSTM_model_classifier import AmadeusFoldClassifier
from utils import index_sequence

order_index_reverse = {
	0: 0,
	1: 1,
	2: -1,
	3: 2,
	4: -2,
	5: 3,
	6: -3
}

def main():
	# settings
	model_path = './models/<your_model>.pt'
	sequence = 'HHPHHCCHPHCPP' # 'HHCHHCPCHCH'
	multiplier = 8 # check dataset settings

	# load model and weight
	print("Loading model...")
	model = AmadeusFoldClassifier(n_molecules = 4)
	model.load_state_dict(torch.load(model_path))
	print("Model loaded.")

	# process and feed sequence to model
	sequence_indexed = torch.tensor(index_sequence(sequence))
	with torch.no_grad():
		model_output = model(sequence_indexed)
	indexed_order = list(torch.argmax(model_output, dim=2))[0]
	order = [order_index_reverse[int(direction)] for direction in indexed_order][1:]
	print(order)

if __name__ == "__main__":
	main()
