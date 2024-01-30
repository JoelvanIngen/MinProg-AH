from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import Greedy
from random import shuffle


def main():
    create_experiment_folders()
    sequences = {"HHPHHHPHPHHHPH": 0,
                 "HPHPPHHPHPPHPHHPPHPH": 0,
                 "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP": 0,
                 "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH": 0,
                 "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP": 0,
                 "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC": 0,
                 "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH": 0,
                 "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH": 0}
    DIM = 3
    CYCLES = 1000

    for _ in range(CYCLES):
        print(f"Cycling {_}")
        for sequence in sequences.keys():
            protein = Protein(sequence)

            algorithm = Greedy(protein, dimensions=DIM, debug=False)
            algorithm.run()
            if sequences[sequence] > protein.get_bond_score():
                sequences[sequence] = protein.get_bond_score()
                print(f"Sequence: {sequence} Score: {protein.get_bond_score()}")

            # if dim == 2:
            #     protein.plot(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
            # elif dim == 3:
            #     protein.plot_3d(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')

    print(sequences)

if __name__ == '__main__':
    main()
