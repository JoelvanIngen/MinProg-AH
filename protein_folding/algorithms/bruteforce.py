from . import Algorithm
from protein_folding.protein import Protein
from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score
from typing import Iterable


class ScoreTracker:
    def __init__(self, max_scores: int = 10):
        self.scores = {}
        self.max_scores = max_scores
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


def generate_new_combination(directions: list[int], prev_combination: Iterable[int]) -> Iterable[int] | None:
    new_combination = list(prev_combination)
    for i, direction in enumerate(new_combination):
        if direction == directions[-1]:
            new_combination[i] = directions[0]
            continue
        else:
            new_combination[i] = directions[directions.index(direction) + 1]
            # if fast_validate_protein(new_combination):
            #     return tuple(new_combination)
            return tuple(new_combination)
    return None


class BruteForce(Algorithm):

    def __init__(self, protein: 'Protein', dimensions: int, *args, max_iterations=0, verbose=False, **kwargs):
        super().__init__(protein, dimensions, *args, **kwargs)
        self.sequence = self.protein.sequence
        self.n = len(self.protein.sequence) - 1
        # self.order_list = generate_combinations(self.directions, self.n)

        # if 0 < max_iterations < len(self.order_list):
        #     self.order_list = sample(self.order_list, max_iterations)
        
        self.verbose = verbose

        self.configs = 1
        if max_iterations > 0:
            self.max_configs = max_iterations
        else:
            self.max_configs = 4 ** self.n
        self.score_tracker = ScoreTracker(max_scores=10)

    def run(self) -> dict:
        # for i, order in zip(range(self.configs), self.order_list):
        order = tuple([-2] * (len(self.protein.sequence) - 1))
        for _ in range(self.max_configs):
            if fast_validate_protein(order):
                score = fast_compute_bond_score(seq=self.sequence, order=order)
                if score < self.score_tracker.lowest_top_score:
                    if self.verbose:
                        print(f"Found score {score} on config {self.configs}/{self.max_configs}"
                              f"({100 * self.configs // self.max_configs}%)")
                    self.score_tracker.add_score(order, score)
            order = generate_new_combination(self.directions, order)
            self.configs += 1
            if self.verbose:
                print(f"Config {self.configs}/{self.max_configs}"
                      f"({100 * self.configs // self.max_configs}%)")
        print(f"Stopped at config {self.configs} with order {order}")
        return self.score_tracker.get_scores()
