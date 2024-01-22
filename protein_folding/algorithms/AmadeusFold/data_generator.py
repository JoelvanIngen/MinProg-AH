import sys
sys.path.append('../../..')

import pandas as pd
import random
from tqdm import tqdm

from protein_folding.algorithms.bruteforce import BruteForce
from protein_folding.protein import Protein
from utils import sequence_generator

def data_generator(
	n_datapoints: int = 5000,
	seq_len_limits: list = [5, 10],
	max_bruteforce_iterations: int = 50000,
	csv_loc: str = './data/unnamed_data.csv'
	) -> None:

	dataset = list()

	# loop over number of desired sequences
	i = 0
	pbar = tqdm(total=n_datapoints)
	while i < n_datapoints:

		# generate a random sequence of random length within bounds
		n_molecules = random.randint(seq_len_limits[0], seq_len_limits[1])
		sequence = sequence_generator(n_molecules)
		protein = Protein(sequence)

		# generate brute force solutions
		algorithm = BruteForce(protein, dimensions=2, max_iterations=max_bruteforce_iterations)
		results = algorithm.run()

		# obtain data from solutions: sequence, score, order, node coordinates
		for	order in results:
			sequence_str = "".join(sequence)
			score = results[order]
			protein.set_order(order)
			coordinates = protein.get_node_coordinates()
			datapoint = [sequence_str, score, order, coordinates]
			dataset.append(datapoint)
			i += 1
			pbar.update(1)
	
	pbar.close()
	
	# convert to dataframe and save
	df = pd.DataFrame(dataset, columns=['sequence', 'score', 'order', 'coordinates'])
	df.to_csv(csv_loc)


if __name__ == "__main__":
	data_generator()