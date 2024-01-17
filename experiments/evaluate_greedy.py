from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import Greedy
from random import shuffle


def main():
    create_experiment_folders()
    sequence = 'H' * 50
    sequence += 'C' * 10
    sequence += 'P' * 10
    l = list(sequence)
    shuffle(l)
    sequence = ''.join(l)

    protein = Protein(sequence)

    algorithm = Greedy(protein, dimensions=2, debug=True)
    score = algorithm.run()
    print(f"Score: {score}")

    protein.plot('./output/evaluate_greedy_protein_len300.png')


if __name__ == '__main__':
    main()