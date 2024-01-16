from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import Spiral


def main():
	create_experiment_folders()
	sequence = 'H' * 30
	protein = Protein(sequence)

	algorithm = Spiral(protein, dimensions=2, debug=True)
	score = algorithm.run()
	print(f"Score: {score}")

	protein.plot('./output/evaluate_spiral_protein_len30.png')


if __name__ == '__main__':
	main()
