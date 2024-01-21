import torch
from torch import nn
from torch.utils.data import DataLoader

_molecule_indices = {'H': 0, 'C': 1, 'P': 2}


class AmadeusFold(nn.Module):
	"""
	A pytorch-based model to predict an optimal fold for a molecule sequence.
	Made up of:
		- an embedding layer: map letter connections to mathematical
		representation
		- lstm layers: a form of a recurrent neural network, best for
		sequenced data
		- a linear(dense) layer: for mapping sequential lstm output to
		3-dimensional representations
	
	Pre:
		- n_molecules is an int representing the number of different types
		of molecules that can exist in a sequenc, so 3 for ['H', 'C', 'P']
	Post:
		- a functioning model object is instantiated
	"""

	def __init__(self, n_molecules = 3):
		super().__init__()

		# settings
		self.embedding_size = 32
		self.lstm_hidden_size = 128
		self.lstm_n_hidden_layers = 1
		self.prediction_dimension = 3

		# embedding for vector mappings
		self.embedding = nn.Embedding(
			num_embeddings=n_molecules,
			embedding_dim=self.embedding_size
		)

		# rnn(lstm) or transformer? for sequential data
		self.lstm = nn.LSTM(
			input_size=self.embedding_size,
			hidden_size=self.lstm_hidden_size,
			num_layers=self.lstm_n_hidden_layers,
			batch_first=True
		)

		# dense layer to map to spatial data
		self.output_layer = nn.Linear(
			in_features=self.lstm_hidden_size,
			out_features=self.prediction_dimension
		)

	def forward(self, x):
		seq_embedded = self.embedding(x)
		output_lstm, _ = self.lstm(seq_embedded)
		output = self.output_layer(output_lstm)

		return output


def main():

	# check for available hardware to speed up calculations
	device = (
		"cuda"
		if torch.cuda.is_available()
		else "mps"
		if torch.backends.mps.is_available()
		else "cpu"
	)

	# instantiate model (untrained, so ouput doesn't make sense)
	model = AmadeusFold().to(device)
	
	# choose a sample sequence and convert to tensor of ints
	sequence = 'HCHHPHH'
	seq_indexed = [[_molecule_indices[letter] for letter in sequence]]
	
	# feed indexed input to model
	test_input = torch.tensor(seq_indexed)

	# determine and display output
	output = model(test_input)
	print(f"Input: {sequence}\nInput indexed: {test_input}\nOutput: {output}")


if __name__ == "__main__":
	main()

