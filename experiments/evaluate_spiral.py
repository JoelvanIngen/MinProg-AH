from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import Spiral
from random import shuffle


def main():
    create_experiment_folders()
    # sequence = 'H' * 10
    # sequence += 'C' * 50
    # sequence += 'P' * 10
    # l = list(sequence)
    # shuffle(l)
    # sequence = ''.join(l)
    sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"

    protein = Protein(sequence)

    algorithm = Spiral(protein, dimensions=2, debug=True)
    score = algorithm.run()
    print(f"Score: {score}")

    protein.plot(f'./output/evaluate_spiral_protein_len{len(sequence)}.png')


if __name__ == '__main__':
    main()
