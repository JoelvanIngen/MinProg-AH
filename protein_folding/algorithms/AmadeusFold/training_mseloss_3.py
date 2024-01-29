import random
from datetime import datetime
import sys
import torch
from torch.utils.data import DataLoader, random_split, Subset
from torchvision.transforms import ToTensor


from utils import index_sequence, collate_fn
from dataset import FoldDataset
from LSTM_model import AmadeusFold
from custom_loss import FoldLoss


def train_loop(batch_size, norm_val, dataloader, model, loss_fn, optimizer):
	"""
	Training loop for the AmadeusFold model. Loops over the dataset, divides it
	into batches, determines loss and performs gradient descent.
	"""
	# not necessary but good practice
	model.train()

	# loop over dataset
	loss_total = 0
	for (sequences, coordinates) in dataloader:

		# loop over batch and sum total loss
		predicted_coords = model(sequences)
		loss = loss_fn(predicted_coords, coordinates).float()
		loss_total += loss
		
		# perform gradient descent
		loss.backward()
		optimizer.step()
		optimizer.zero_grad()
	
	return loss_total


def test_loop(batch_size, norm_val, dataloader, model, loss_fn):
	"""
	Evaluation loop for the AmadeusFold model. Loops over the dataset, divides
	it into batches, determines loss.
	"""
	# not necessary but good practice
	model.eval()

	# loop over dataset
	loss_total = 0
	with torch.no_grad():
		# loop over dataset
		for (sequences, coordinates) in dataloader:

			# loop over batch and sum total loss
			predicted_coords = model(sequences)
			loss = loss_fn(predicted_coords, coordinates).float()
			loss_total += loss

	return loss_total


def main():
	# settings
	model_path = './models/first_test/model_2d_5k_mse.pt' # where to save weights
	norm_val = 8 # normalisation value (see FoldDataset)
	train_test_split = .8 # ratio training/testing data split
	learning_rate = 1e-3 # model learning rate
	batch_size = 64 # amount of datapoints in batch
	epochs = 2000 # number of epochs to train
	patience = 10 # number of epochs validation loss can rise before stopping
	n_datapoints = 50000
	print(f"Beginning training of model {model_path}")

	# instantiate datasets for training and validation
	dataset_full = FoldDataset('./data/50k_datapoints.csv', norm_val = norm_val)
	dataset = Subset(
			dataset_full,
			range(n_datapoints)
		)
	train_len = int(train_test_split * n_datapoints)
	test_len = n_datapoints - train_len
	train_dataset, test_dataset = random_split(dataset, [train_len, test_len])
	train_dataloader = DataLoader(train_dataset, batch_size=batch_size, collate_fn=collate_fn)
	test_dataloader = DataLoader(test_dataset, batch_size=batch_size, collate_fn=collate_fn)
	print("Datasets loaded")

	# instantiate model, loss function and optimizer
	model = AmadeusFold(n_molecules = 4) # for padding
	loss_fn = torch.nn.MSELoss() # FoldLoss()
	optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)

	# variables for overfitting
	curr_min_val_loss = 999
	overfit_counter = 0

	# start training procedure
	print(f"Started training at {datetime.now()}")
	for ep in range(epochs):

		# determine training loss
		loss_train = train_loop(
			batch_size, 
			norm_val, 
			train_dataloader, 
			model, 
			loss_fn, 
			optimizer
		)
		# determine validation loss
		loss_test = test_loop(
			batch_size, 
			norm_val, 
			test_dataloader, 
			model, 
			loss_fn
		)
		
		# display status update
		print(f"Epoch {ep + 1}/{epochs}. training loss: {float(loss_train/train_test_split)} eval loss: {float(loss_test/(1-train_test_split))} overfits: {overfit_counter}/{patience} time: {datetime.now()}")
		
		# check for overfitting and update validation loss
		if loss_test > curr_min_val_loss:
			overfit_counter += 1
		else:
			curr_min_val_loss = loss_test
			overfit_counter = 0
		
		# abort training if overfitting
		if overfit_counter > patience:
			print(f"Patience reached. Aborting training.")
			break

	print(f"Training finished. Saving model at {model_path}")

	# save final weights
	torch.save(model.state_dict(), model_path)


if __name__ == "__main__":
	main()



	