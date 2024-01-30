from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import Spiral
from random import shuffle


def main():
    create_experiment_folders()
    sequences = ["HHPHHHPHPHHHPH",
                 "HPHPPHHPHPPHPHHPPHPH",
                 "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP",
                 "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH",
                 "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP",
                 "CPPCHPPCHPPCPPHHHHHHCCPCHPPCPCHPPHPC",
                 "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH",
                 "HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH"]

    for sequence in sequences:
        protein = Protein(sequence)

        algorithm = Spiral(protein, dimensions=3, debug=True)
        score = algorithm.run()
        print(f"Sequence: {sequence} Score: {score}")

        # protein.plot(f'./output/{algorithm.get_name()}_len{len(sequence)}.png')


if __name__ == '__main__':
    main()
