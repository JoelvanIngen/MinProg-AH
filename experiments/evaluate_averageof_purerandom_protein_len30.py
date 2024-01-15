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


if __name__ == '__main__':
    main()
