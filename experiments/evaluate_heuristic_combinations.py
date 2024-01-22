from experiments_helper import get_available_heuristics, generate_random_sequence
import itertools
from tqdm import tqdm

from protein_folding.protein import Protein
from protein_folding.algorithms import *
from protein_folding.algorithms.heuristics import Heuristic

PROTEIN_LENGTH = 10
N_ITERATIONS = 10


def create_heuristic_combinations(heuristics_list):
    """
    Creates a list containing all possible heuristics combinations for
    evaluation.

    Code adapted from Stackoverflow:
    https://stackoverflow.com/questions/8371887/making-all-possible-combinations-of-a-list
    """
    all_combinations = []

    for i in range(1, len(heuristics_list) + 1):
        els = [list(x) for x in itertools.combinations(heuristics_list, i)]
        all_combinations.extend(els)

    return all_combinations


def run_all_combinations():
    algorithms = (
        IterativeGreedy,
        SimulatedAnnealingHeuristics,
    )

    number_of_runs = len(algorithms) * len(create_heuristic_combinations(get_available_heuristics())) * N_ITERATIONS
    pbar = tqdm(range(number_of_runs))

    avg_scores = []
    for algorithm in algorithms:
        avg_scores.append(run_each_heuristic_for_algorithm(algorithm, pbar))

        print(list(zip(get_available_heuristics(), avg_scores)))


def run_each_heuristic_for_algorithm(algorithm: Algorithm, pbar):
    scores_per_heuristic = []
    for heuristic_combination in tqdm(
            create_heuristic_combinations(get_available_heuristics()),
            desc=f"Algorithm: {algorithm}"
    ):
        scores_per_heuristic.append(repeat_run_algorithm(algorithm, heuristic_combination, pbar))

    return [sum(scores) / len(scores) for scores in scores_per_heuristic]


def repeat_run_algorithm(algorithm: Algorithm, heuristic_combination: list[Heuristic], pbar):
    return [run_algorithm(algorithm, heuristic_combination, pbar) for _ in range(N_ITERATIONS)]


def run_algorithm(algorithm_to_run, heuristics_combination: list[Heuristic], pbar):
    pbar.update(1)
    pbar.desc = f"{algorithm_to_run, heuristics_combination}"
    protein = Protein(generate_random_sequence(PROTEIN_LENGTH))
    algorithm = algorithm_to_run(protein, dimensions=2, heuristics=heuristics_combination, show_progressbar=False)
    return algorithm.run()


def main():
    run_all_combinations()


if __name__ == '__main__':
    main()
