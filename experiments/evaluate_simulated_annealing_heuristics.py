from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import SimulatedAnnealingHeuristics


def main():
	create_experiment_folders()
	sequence = 'HHPCHHPCCPCPPHHHHPPHCHPHPHCHPPHHPCHHPCCPCPPHHHHPPHCHPHPHCHPPHHPCHHPCCPCPPHHHHPPHCHPHPHCHPP'
	# sequence = 'HHHHHHHHHH'
	protein = Protein(sequence)

	algorithm = SimulatedAnnealingHeuristics(protein, dimensions=2, debug=True)
	score = algorithm.run()
	print(f"Score: {score}")

	protein.plot(f'./output/{algorithm.get_name()}_len{len(sequence)}.png')


if __name__ == '__main__':
	main()
