from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.definitions import *


if __name__ == '__main__':
    # Create necessary folders
    create_experiment_folders()

    seq = "HHPHHHPHPHHHPH"
    # order = [UP, RIGHT, UP, LEFT, LEFT, DOWN, DOWN]
    orders = {(-1, -1, 2, 1, 2, 1, -2, 1, -2, -2, -2, -1, 2): -6, (2, 2, 1, -2, 1, -2, -1, -2, -1, -1, -1, 2, 1): -6, (2, -1, -2, -2, 1, 1, -2, -2, -1, 2, -1, -1, 2): -6, (-2, -1, 2, 2, 1, 2, -1, 2, -1, -2, -1, -2, 1): -6, (-1, 2, 1, 1, -2, 1, 2, 2, -1, -1, 2, -1, -2): -6, (-1, -1, -2, 1, 1, 1, 2, 1, 2, -1, -1, 2, 1): -6, (-1, -2, 1, 1, -2, 1, 2, 1, 2, -1, -1, 2, -1): -6, (-2, -2, -1, 2, 2, 2, 1, 1, -2, -2, 1, -2, -1): -6, (2, -1, -2, -1, 2, -1, -2, -2, 1, 1, 1, -2, -1): -6, (-2, 1, 2, 2, 1, 2, -1, 2, -1, -2, -1, -2, 1): -6}
    print(min(orders.values()))
    for i, order in enumerate(orders.keys()):
        protein = Protein(seq)
        protein.set_order(list(order))

        protein.plot(filename=f'./output/plot_predetermined_protein{i}.png')
