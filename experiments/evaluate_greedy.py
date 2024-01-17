from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import Greedy
from random import shuffle


def main():
    create_experiment_folders()
    sequence = 'H' * 40
    sequence += 'C' * 40
    sequence += 'P' * 0
    l = list(sequence)
    shuffle(l)
    sequence = ''.join(l)
    # sequence = "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"

    protein = Protein(sequence)
    dim = 2
    algorithm = Greedy(protein, dimensions=dim, debug=True)
    score = algorithm.run()
    print(f"Score: {score}")

    protein.plot(f'./output/evaluate_greedy_protein_len{len(sequence)}_dim{dim}.png')


if __name__ == '__main__':
    main()