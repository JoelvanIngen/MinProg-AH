from protein_folding.protein import Protein
from protein_folding.definitions import *


if __name__ == '__main__':
    seq = "PHPPHPHH"
    order = [UP, RIGHT, UP, LEFT, LEFT, DOWN, DOWN]

    protein = Protein(seq)
    protein.set_order(order)

    protein.plot()
