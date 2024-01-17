import cProfile
import pstats

from experiments_helper import create_bruteforce_folders
from protein_folding.protein import Protein
from random import shuffle
import itertools

MAX_ITERATIONS = 1000000

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


def generate_combinations(n):
    directions = [1, -1, 2, -2]
    return list(list(itertools.product(directions, repeat=n)))


def run():
    # Create necessary folders
    create_bruteforce_folders()

    seq = "HHPHHHPHPHHHPH"

    n = len(seq) - 1
    order_list = generate_combinations(n)
    configs = min(len(order_list), MAX_ITERATIONS)
    shuffle(order_list)
    print(f"Amount of configs: {configs}")

    score_tracker = ScoreTracker()

    for i, order in enumerate(order_list):
        if i > MAX_ITERATIONS:
            break
        print(f"Checking config {i}/{configs} ({i / configs * 100:.0f}%)")
        protein = Protein(seq)
        protein.set_order(order)
        if protein.has_valid_order():
            score = protein.get_bond_score()
            print(order, score)
            if score <= score_tracker.lowest_top_score:
                score_tracker.add_score(order, score)
                score_tracker.lowest_top_score = score_tracker.get_worst_score()

    results = score_tracker.get_scores()
    print(results)
    print(len(results))
    print(f"Best score: {score_tracker.get_best_score()}")

    for i, order in enumerate(results.keys()):
        protein = Protein(seq)
        protein.set_order(order)
        score = results[order]
        protein.plot(filename=f'./bf_output/bf_#{i + 1}_score{score}.png')


if __name__ == '__main__':
    profiler = cProfile.Profile()
    profiler.runcall(run)
    stats = pstats.Stats(profiler)
    stats.strip_dirs().sort_stats('cumulative').print_stats(10)
