from . import Algorithm
import itertools
from random import shuffle
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


class BruteForce(Algorithm):
    def __init__(self, protein: 'Protein', dimensions: int, *args, max_iterations=200000, **kwargs):
        super().__init__(protein, dimensions, *args, **kwargs)
        self.sequence = self.protein.sequence
        self.n = len(self.protein.sequence) - 1
        self.order_list = generate_combinations(self.directions, self.n)
        shuffle(self.order_list)

        self.configs = min(len(self.order_list), max_iterations)

        self.score_tracker = ScoreTracker()

    def run(self) -> dict:
        for i, order in zip(range(self.configs), self.order_list):
            print(f"Checking config {i}/{self.configs} ({i / self.configs * 100:.0f}%)")
            if fast_validate_protein(order):
                score = fast_compute_bond_score(seq=self.sequence, order=order)
                print(order, score)
                if score <= self.score_tracker.lowest_top_score:
                    self.score_tracker.add_score(order, score)
        return self.score_tracker.get_scores()
