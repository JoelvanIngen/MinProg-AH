from experiments_helper import get_available_heuristics, generate_random_sequence
import itertools
from tqdm import tqdm

from protein_folding.protein import Protein
from protein_folding.algorithms import *
from protein_folding.algorithms.heuristics import Heuristic

PROTEIN_LENGTH = 10
N_ITERATIONS = 10


# TODO: Change name? What would we actually call something that makes combinations?
class Combinator:
    def __init__(self, algorithms=None, heuristics=None, show_progressbar=True):
        self.algorithms = algorithms if algorithms else get_algorithms()
        self.heuristics = heuristics if heuristics else get_heuristics()

        self.heuristic_combinations = create_heuristic_combinations(self.heuristics)

        self.avg_scores = []

        self.show_progressbar = show_progressbar
        self.pbar = None

    def run_all(self):
        if self.show_progressbar:
            number_of_runs = len(self.algorithms) * len(self.heuristic_combinations) * N_ITERATIONS
            self.pbar = tqdm(range(number_of_runs))

        for a in self.algorithms:
            self.run_all_heuristics(a)

    def run_all_heuristics(self, a):
        for combination in self.heuristic_combinations:
            runs_avg = self.repeat_run_algorithm(a, combination)
            self.avg_scores.append((f"{a.name} - {[h.name for h in combination]}", runs_avg))

    def repeat_run_algorithm(self, a, heuristics):
        return avg([self.run_algorithm_once(a, heuristics=heuristics) for _ in range(N_ITERATIONS)])

    def run_algorithm_once(self, a, heuristics):
        if self.pbar:
            self.pbar.update(1)
            self.pbar.desc = f"Algorithm: {a.name}, Heuristics: {[h.name for h in heuristics]}"

        protein = Protein(generate_random_sequence(PROTEIN_LENGTH))
        algorithm = a(protein, dimensions=2, heuristics=self.heuristic_combinations, show_progressbar=False)
        score = algorithm.run()
        return score

    def print_scores(self):
        print(self.avg_scores)


def get_algorithms():
    return (
        IterativeGreedy,
        SimulatedAnnealingHeuristics,
    )


def get_heuristics():
    return get_available_heuristics()


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


def avg(values: list[float]) -> float:
    return sum(values) / len(values)


def main():
    c = Combinator()
    c.run_all()

    c.print_scores()


if __name__ == '__main__':
    main()
