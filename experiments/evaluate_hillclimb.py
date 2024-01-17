from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms.hillclimb_1 import IterativeRandomHillClimb


def main():
    create_experiment_folders()

    sequence = sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"
    protein = Protein(sequence)

    algorithm = IterativeRandomHillClimb(protein, dimensions=2, debug=True)
    score = algorithm.run()

    protein.plot('./output/evaluate_purerandom_protein_len30.png')

    print(f'Final score: {score}')


if __name__ == '__main__':
    main()
