from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.definitions import *


if __name__ == '__main__':
    # Create necessary folders
    create_experiment_folders()

    seq = "PHPPHPHH"
    order = [UP, RIGHT, UP, LEFT, LEFT, DOWN, DOWN]

    protein = Protein(seq)
    protein.set_order(order)

    protein.plot(filename='./output/plot_predetermined_protein.png')
