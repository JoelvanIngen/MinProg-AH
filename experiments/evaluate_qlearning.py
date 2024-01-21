from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms.reinforcement import run_protein_folding


def main():
	create_experiment_folders()
	sequence = 'HHPCHHPCCPCPPHHHHPPHCHPHPHCHPP'

	algorithm = run_protein_folding(sequence, 1000)

	protein.plot(f'./output/{algorithm.get_name()}_len{len(sequence)}.png')


if __name__ == '__main__':
	main()
