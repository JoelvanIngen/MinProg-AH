from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import IterativeRandom


def main():
    create_experiment_folders()

    sequence = 'HHHHHHHH'
    protein = Protein(sequence)

    algorithm = IterativeRandom(protein, dimensions=2)
    score = algorithm.run()

    protein.plot('./output/evaluate_iterativerandom_protein_len8.png')

    print(f'Final score: {score}')


if __name__ == '__main__':
    main()
