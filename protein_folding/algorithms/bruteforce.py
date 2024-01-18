from . import Algorithm
import itertools
from random import shuffle, sample
from protein_folding.protein import Protein
from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score


class ScoreTracker:
    def __init__(self):
        self.scores = {}
        self.max_scores = 10
        self.lowest_top_score = 0

    def add_score(self, order, score):
        # If the number of scores is already 10, remove the lowest score before adding a new one
        if len(self.scores) == self.max_scores:
            self.remove_worst_score()
        self.scores[order] = score
        self.lowest_top_score = self.get_worst_score()

    def remove_worst_score(self):
        if self.scores:
            lowest_order = max(self.scores, key=self.scores.get)
            del self.scores[lowest_order]

    def get_worst_score(self):
        if self.scores:
            return self.scores[max(self.scores, key=self.scores.get)]
        return None

    def get_best_score(self):
        if self.scores:
            return self.scores[min(self.scores, key=self.scores.get)]
        return None

    def get_scores(self):
        return self.scores


def generate_combinations(directions, n):
    return list(list(itertools.product(directions, repeat=n)))


def generate_new_combination(directions, prev_combination):
    new_combination = prev_combination
    for i, direction in enumerate(new_combination):
        if direction == directions[-1]:
            new_combination[i] = directions[0]
            continue
        else:
            new_combination[i] = directions[directions.index(direction) + 1]
            if fast_validate_protein(new_combination):
                return new_combination
    return None


class BruteForce(Algorithm):

    # TODO Implement generate_new_combination instead of generating all combinations
    def __init__(self, protein: 'Protein', dimensions: int, *args, max_iterations=0, **kwargs):
        super().__init__(protein, dimensions, *args, **kwargs)
        self.sequence = self.protein.sequence
        self.n = len(self.protein.sequence) - 1
        self.order_list = generate_combinations(self.directions, self.n)

        if 0 < max_iterations < len(self.order_list):
            self.order_list = sample(self.order_list, max_iterations)

        self.configs = len(self.order_list)

        self.score_tracker = ScoreTracker()

    def run(self) -> dict:
        for i, order in zip(range(self.configs), self.order_list):
            # print(f"Checking config {i + 1}/{self.configs} ({(i + 1) // self.configs * 100:.0f}%)")
            if fast_validate_protein(order):
                score = fast_compute_bond_score(seq=self.sequence, order=order)
                if score < self.score_tracker.lowest_top_score:
                    print(
                        f"Found score {score} on config {i + 1}/{self.configs} ({(i + 1) // self.configs * 100:.0f}%)")
                    self.score_tracker.add_score(order, score)
        return self.score_tracker.get_scores()
