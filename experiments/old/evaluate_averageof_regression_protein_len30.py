import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm
from experiments_helper import create_experiment_folders, generate_random_sequence
from protein_folding.protein import Protein
from protein_folding.algorithms import Regression

N = 10000


def main():
    create_experiment_folders()

    sequence = generate_random_sequence(20)

    scores = []
    for _ in tqdm(range(N)):
        protein = Protein(sequence)

        algorithm = Regression(protein, dimensions=2, debug=False)
        score = algorithm.run()
        scores.append(score)

    std_deviation = np.std(scores)
    average_value = np.mean(scores)
    min_value = min(scores)
    max_value = max(scores)

    print(f'Average final score: {average_value:.2f}')
    print(f"All scores: {scores}")

    plt.hist(scores, bins=int(abs(min_value)), color='blue', edgecolor='black')
    plt.axvline(average_value, color='red', linestyle='dashed', linewidth=2, label=f'Average: {average_value:.2f}')
    plt.axvline(min_value, color='green', linestyle='dashed', linewidth=2, label=f'Min: {min_value}')
    # plt.axvline(max_value, color='orange', linestyle='dashed', linewidth=2, label=f'Max: {max_value}')
    # plt.axvline(average_value - std_deviation, color='purple', linestyle='dashed', linewidth=2,
    #             label=f'-1 Std Dev: {average_value - std_deviation:.2f}')
    # plt.axvline(average_value + std_deviation, color='purple', linestyle='dashed', linewidth=2,
    #             label=f'+1 Std Dev: {average_value + std_deviation:.2f}')

    # Adding labels and title
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Regression with Average, Min, and Max')
    plt.legend()
    plt.savefig('./output/RegressionStats')


if __name__ == '__main__':
    main()
