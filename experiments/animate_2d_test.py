"""
Animates the simulated annealing algorithm in 2D.
Authors: Wolf
"""

from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import SimulatedAnnealing


def main():
	create_experiment_folders()
	sequence = 'HHPCHHHHHPPHCHPHPHCHPP'
	protein = Protein(sequence)

	algorithm = SimulatedAnnealing(protein, dimensions=2, debug=True)
	score = algorithm.run()
	print(f"Score: {score}")
	protein.animate_2d(algorithm.orders)


if __name__ == '__main__':
	main()