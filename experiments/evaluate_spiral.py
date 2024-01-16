from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import Spiral
from random import shuffle


def main():
    create_experiment_folders()
    sequence = 'H' * 100
    sequence += 'C' * 100
    sequence += 'P' * 100
    l = list(sequence)
    shuffle(l)
    sequence = ''.join(l)

    protein = Protein(sequence)

    algorithm = Spiral(protein, dimensions=2, debug=True)
    score = algorithm.run()
    print(f"Score: {score}")

    protein.plot('./output/evaluate_spiral_protein_len300.png')


if __name__ == '__main__':
    main()
