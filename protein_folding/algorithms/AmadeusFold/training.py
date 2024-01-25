import random
import torch
from torch.utils.data import DataLoader, random_split, Subset

from utils import index_sequence
from dataset import FoldDataset
from LSTM_model import AmadeusFold
from custom_loss import FoldLoss


def train_loop(batch_size, norm_val, dataset, model, loss_fn, optimizer):
	"""
	Training loop for the AmadeusFold model. Loops over the dataset, divides it
	into batches, determines loss and performs gradient descent.
	"""
	# not necessary but good practice
	model.train()

	# loop over dataset
	for _ in range(int(len(dataset) / batch_size)):
		# separate dataset into batches
		rand_idx = random.randint(0, len(dataset) - batch_size)
		batch = Subset(
			dataset,
			range(rand_idx, rand_idx + batch_size)
		)

		# loop over batch and sum total loss
		loss = 0
		for (sequence, coordinates, score) in batch:
			indexed_sequence = torch.tensor(index_sequence(sequence))
			predicted_coords = model(indexed_sequence)

			loss += loss_fn(
				predictions = predicted_coords, 
				norm_factor = norm_val, 
				sequence = sequence,
				targets = coordinates, 
				target_score = score,
			)

		# perform gradient descent
		loss.backward()
		optimizer.step()
		optimizer.zero_grad()
	
	return loss/batch_size


def test_loop(batch_size, norm_val, dataset, model, loss_fn):
	"""
	Evaluation loop for the AmadeusFold model. Loops over the dataset, divides
	it into batches, determines loss.
	"""
	# not necessary but good practice
	model.eval()

	# loop over dataset
	with torch.no_grad():
		for _ in range(int(len(dataset) / batch_size)):
			# separate dataset into batches
			rand_idx = random.randint(0, len(dataset) - batch_size)
			batch = Subset(
				dataset,
				range(rand_idx, rand_idx + batch_size)
			)

			# loop over batch and sum total loss
			loss = 0
			for (sequence, coordinates, score) in batch:
				indexed_sequence = torch.tensor(index_sequence(sequence))
				predicted_coords = model(indexed_sequence)

				loss += loss_fn(
					predictions = predicted_coords, 
					norm_factor = norm_val, 
					sequence = sequence,
					targets = coordinates, 
					target_score = score,
				)

	return loss/batch_size


def main():
	# settings
	model_path = './models/first_test/model.pt' # where to save weights
	norm_val = 15 # normalisation value (see FoldDataset)
	train_test_split = .8 # ratio training/testing data split
	learning_rate = 1e-3 # model learning rate
	batch_size = 256 # amount of datapooints in
	epochs = 2 # number of epochs to train
	patience = 10 # number of epochs validation loss can rise before stopping

	# instantiate datasets for training and validation
	dataset = FoldDataset('./data/test_data.csv', norm_val = norm_val)
	train_len = int(train_test_split * len(dataset))
	test_len = len(dataset) - train_len
	train_dataset, test_dataset = random_split(dataset, [train_len, test_len])

	# instantiate model, loss function and optimizer
	model = AmadeusFold()
	loss_fn = FoldLoss()
	optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

	# variables for overfitting
	loss_prev_test = 999
	overfit_counter = 0

	# start training procedure
	for ep in range(epochs):
		# determine training loss
		loss_train = train_loop(
			batch_size, 
			norm_val, 
			train_dataset, 
			model, 
			loss_fn, 
			optimizer
		)
		# determine validation loss
		loss_test = test_loop(
			batch_size, 
			norm_val, 
			train_dataset, 
			model, 
			loss_fn
		)
		
		# display status update
		print(f"Epoch {ep + 1}/{epochs}. training loss: {loss_train} eval loss: {loss_test}")
		
		# check for overfitting and update validation loss
		if loss_test > loss_prev_test:
			overfit_counter += 1
		loss_prev_test = loss_test
		
		# abort training if overfitting
		if overfit_counter > patience:
			print(f"Patience reached. Aborting training.")
			break

	print(f"Training finished. Saving model at {model_path}")

	# save final weights
	torch.save(model.state_dict(), model_path)


if __name__ == "__main__":
	main()



	