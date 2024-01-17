from experiments_helper import create_experiment_folders, sequence1
from protein_folding.protein import Protein
from protein_folding.algorithms.hillclimb import HillClimb


def main():
    create_experiment_folders()

    sequence = sequence1
    protein = Protein(sequence)
    print(f"{protein.get_order}")
    

    algorithm = HillClimb(protein, dimensions=2, debug=True)
    score = algorithm.run()

    protein.plot('./output/evaluate_Hillclimb.png')

    print(f'Final score: {score}')


if __name__ == '__main__':
    main()
