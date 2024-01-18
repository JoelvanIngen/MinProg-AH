from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import Regression


def main():
	create_experiment_folders()
	sequence = 'HHPHHHPH'
	protein = Protein(sequence)

	algorithm = Regression(protein, dimensions=2, debug=True)
	score = algorithm.run()
	print(f"Score: {score}")

	protein.plot(f'./output/evaluate_regression_protein_len{len(sequence)}.png')

if __name__ == '__main__':
	main()
