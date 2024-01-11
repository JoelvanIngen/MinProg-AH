import os
from protein_folding.protein import Protein
from protein_folding.definitions import *


if __name__ == '__main__':
    seq = "PHPPHPHH"
    order = [UP, RIGHT, UP, LEFT, LEFT, DOWN, DOWN]

    protein = Protein(seq)
    protein.set_order(order)

    # Create output directory if it does not exist yet
    if not os.path.exists('./output/'):
        os.makedirs('./output/')

    protein.plot(filename='./output/plot_predetermined_protein.png')
