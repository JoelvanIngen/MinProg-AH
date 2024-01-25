import random
from tqdm import tqdm

from . import Algorithm
from .heuristics import *
from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein_folding.protein import Protein


class SimulatedAnnealingHeuristics(Algorithm):
    """
    A simulated annealing algorithm that compares the current protein ordering
    against a randomly permutated but legal state and sets the current state to
    said permutation if the bond score of that state is lower, and if a
    randomly chosen float is above a perpetually lowering threshold, inspired
    by the way metals anneal.
    """

    def __init__(self, protein: 'Protein', dimensions,
                 reset_threshold: int = 1000,
                 n_permutations: int = 5000, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

        decrease = 0.9997
        self.thresholds = [decrease ** i for i in range(n_permutations)]

        # Current iteration counter
        self.run_iteration = 0
        self._iteration = 0

        # Amount of valid permutations that the algorithm will perform in total
        self.n_permutations = n_permutations

        # Amount of iterations without improvement that are permitted before
        # the algorithm resets itself
        self.reset_threshold = reset_threshold
        self.iterations_until_reset = reset_threshold

        # Run stats
        self.overall_best_score = 1
        self.overall_best_order = []

    def _increment_iteration(self):
        self.run_iteration += 1
        self._iteration += 1
        if self.pbar:
            self.pbar.update(1)

    def _compare_score(self, order, run_best_score, run_best_order, threshold):
        score = fast_compute_bond_score(self.protein.sequence, order[1:])

        if score >= run_best_score:
            self.iterations_until_reset -= 1

        if score <= run_best_score or random.random() < threshold:
            self.protein.set_order(order[1:])
            return score, order
        else:
            return run_best_score, run_best_order

    def _find_valid_move(self, idx, order, directions):
        for direction in directions:
            order[idx] = direction
            if fast_validate_protein(order[1:]):
                return direction

        return None

    def _select_random_node(self):
        idx = random.randint(1, len(self.protein.sequence) - 1)
        node = self.protein.nodes[idx]

        return idx, node

    def _run_attempt(self):
        self.run_iteration = 0
        run_best_score = 1
        run_best_order = []

        self.iterations_until_reset = self.reset_threshold

        while True:
            self._increment_iteration()

            threshold = self.thresholds[self.run_iteration]

            node_idx, node = self._select_random_node()

            self.protein.unghost_all()
            free_directions = node.get_free_directions(self.directions)
            if not free_directions:
                continue

            free_directions_legal = []
            for direction in free_directions:
                test_order = self.protein.order[:]
                test_order[node_idx] = direction
                if fast_validate_protein(test_order[1:]):
                    free_directions_legal.append(direction)

            if not self.heuristics or random.random() < threshold:
                free_directions_sorted = free_directions
                random.shuffle(free_directions_sorted)
            else:
                _, free_directions_sorted = self._process_heuristics(
                    node_idx, free_directions_legal, unghost=True)

            # Copy protein order to experiment on
            test_order = self.protein.order[:]
            chosen_direction = self._find_valid_move(node_idx, test_order, free_directions_sorted)
            if not chosen_direction:
                continue

            run_best_score, run_best_order = self._compare_score(test_order, run_best_score, run_best_order, threshold)

            if self.iterations_until_reset < 0 or self._iteration > self.n_permutations:
                break

        if run_best_score < self.overall_best_score:
            self.overall_best_score = run_best_score
            self.overall_best_order = run_best_order
            return run_best_score, run_best_order

    def run(self) -> float:
        """
        Runs the algorithm up to the maximum amount of iterations, and resets
        if there has been no improvement for `self.reset_threshold` iterations.
        """
        if self._show_progress:
            self.pbar = tqdm(range(self.n_permutations))

        while True:
            self._run_attempt()

            if self._iteration > self.n_permutations:
                self.protein.set_order(self.overall_best_order[1:])
                return self.overall_best_score
