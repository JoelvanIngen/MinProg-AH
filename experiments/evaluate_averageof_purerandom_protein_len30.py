import matplotlib.pyplot as plt
import numpy as np
from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import PureRandom

N = 1000


def main():
    create_experiment_folders()

    sequence = 'H' * 30

    scores = []
    for _ in range(N):
        protein = Protein(sequence)

        algorithm = PureRandom(protein, dimensions=2, debug=False)
        score = algorithm.run()
        scores.append(score)

    score_avg = sum(scores) / N

    print(f'Average final score: {score_avg}')
    print(f"All scores: {scores}")

    plt.hist(scores, bins=10, color='blue', edgecolor='black')
    plt.axvline(np.mean(scores), color='red', linestyle='dashed', linewidth=2, label=f'Average: {np.mean(scores):.2f}')
    plt.axvline(min(scores), color='green', linestyle='dashed', linewidth=2, label=f'Min: {min(scores)}')
    plt.axvline(max(scores), color='orange', linestyle='dashed', linewidth=2, label=f'Max: {max(scores)}')

    # Adding labels and title
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Histogram of Values with Average, Min, and Max')
    plt.legend()
    plt.show()
if __name__ == '__main__':
    main()
