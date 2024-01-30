import time

from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import SimulatedAnnealing, SimulatedAnnealingMinDim


def main():
	"""
	This file runs the simulated annealing or simulated annealing with the
	minimise dimensions heuristic for a set amount of time.
	"""
	start = time.time()
	n_runs = 0
	n_sec_total = 60
	sequence = 'HHPCHHPCCPCPPHHHHPPHCHPHPHCHPP'
	create_experiment_folders()
	scores = list()

	while time.time() - start < n_sec_total:
		print(f"run: {n_runs} runtime: {time.time() - start}")
	
		protein = Protein(sequence)
		algorithm = SimulatedAnnealingMinDim(protein, dimensions=2, debug=True)
		score = algorithm.run()
	
		scores.append(score)
		n_runs += 1
		print(f"Score: {score} order: {protein.order}")
	
	print(f"scores = {scores}")

if __name__ == '__main__':
	main()
