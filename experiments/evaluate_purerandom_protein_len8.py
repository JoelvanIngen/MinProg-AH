from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import PureRandom


def main():
    create_experiment_folders()

    sequence = 'HHHHHHHH'
    protein = Protein(sequence)

    algorithm = PureRandom(protein, dimensions=2)
    score = algorithm.run()

    protein.plot('./output/evaluate_purerandom_protein_len8.png')

    print(f'Final score: {score}')


if __name__ == '__main__':
    main()