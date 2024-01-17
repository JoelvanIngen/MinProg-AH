from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms.cystinhillclimbers import CystinHillClimber 

def main():
    create_experiment_folders()

    sequence = 'HHHHHCCCHHHHH'  
    protein = Protein(sequence)

    start_point = "straight_folded"  # Choose the starting point for folding ("straight_folded", "random_folded", "dept_chain")
    iterations = 1000  # Adjust the number of iterations as needed
    dimension = 2  

    algorithm = CystinHillClimber(protein, dimensions=2, debug=True)
    algorithm.execute(protein, start_point=start_point, iterations=iterations, dimension=2)


    score = protein.get_bond_score()

    protein.plot('./output/evaluate_cystinhillclimber_protein.png')

    print(f'Final score: {score}')

if __name__ == '__main__':
    main()