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
                 reset_threshold: int = 5000, **kwargs):
        super().__init__(protein, dimensions, **kwargs)
        # Decrease of threshold value per iteration
        self.decrease = .9997

        # Current iteration counter
        self._iteration = 0

        # Amount of valid permutations that the algorithm will perform in total
        self.n_permutations = 60000

        # Amount of iterations without improvement that are permitted before
        # the algorithm resets itself
        self.reset_threshold = reset_threshold

    def get_permutated_directions(self, node_idx: int):
        """
        TODO: should be implemented within protein
        A function to obtain a set of directions that differs from the current
        protein ordering starting from a direction change at node [node_idx].
        The returned directions are an ordering of the same sequence as
        self.protein, but not necessarily legal.

        Pre:
            - node_idx is an index between 1 and len(sequence) - 1
        Post:
            - returns a list of directions that can be fed to protein.set_order
        """
        dirs_total = [None] + self.protein.get_order()
        node = self.protein.nodes[node_idx]
        free_directions = node.get_free_directions(self.directions)

        if free_directions:
            direction = random.choice(free_directions)
            dirs_total[node_idx] = direction

        dirs_total = dirs_total[1:]

        return dirs_total

    def _compare_score(self, order, run_best_score, threshold):
        score = fast_compute_bond_score(self.protein.sequence, test_order)
        if score <= run_best_score or random.random() < threshold:
            return score,

    def _find_valid_move(self, idx, order, directions):
        for direction in directions:
            order[idx] = direction
            if fast_validate_protein(order):
                return direction

        return None

    def _select_node(self):
        idx = random.randint(1, len(self.protein.sequence) - 1)
        node = self.protein.nodes[idx]

        return idx, node

    def _run_attempt(self, overall_best_score: int):
        threshold = 1
        run_best_score = 1
        run_best_order = []

        while True:
            node_idx, node = self._select_node()

            free_directions = node.get_free_directions(self.directions)
            if not free_directions:
                continue

            if random.random() < threshold:
                free_directions_sorted = random.shuffle(free_directions)
            else:
                _, free_directions_sorted = self._process_heuristics(
                    node_idx, free_directions, self.heuristics)

            # Copy protein order to experiment on
            test_order = self.protein.order[:]
            chosen_direction = self._find_valid_move(node_idx, test_order, free_directions_sorted)
            if not chosen_direction:
                continue

            self._compare_score(test_order, run_best_score, threshold)

        if run_best_score < overall_best_score:
            return run_best_score, run_best_order

    def run(self) -> float:
        """
        Runs the algorithm up to the maximum amount of iterations, and resets
        if there has been no improvement for `self.reset_threshold` iterations.
        """
        pbar = tqdm(range(self.n_permutations))

        while True:
            self._run_attempt()


    def run(self) -> float:
        """
        Runs a simulated annealing algorithm for n_permutations iterations.

        Post:
            - self.protein is ordered in the way that the algorithm found to
            maximise the bond score.
        """
        progressbar = tqdm(range(self.n_permutations))

        best_order = []
        best_order_overall = []
        best_score = 1
        threshold = 1

        while True:
            # print(self.protein.order)
            # Get random node index to permutate and determine its free directions
            node_idx = random.randint(1, len(self.protein.sequence) - 1)
            node = self.protein.nodes[node_idx]
            free_directions = node.get_free_directions(self.directions)

            # print(f"Node: {node_idx}")

            if not free_directions:
                # if self._debug:
                #     print("No free directions found")
                continue

            if random.random() < threshold:
                free_directions_sorted = free_directions
            else:
                _, free_directions_sorted = self._process_heuristics(
                    node_idx, free_directions, self.heuristics)

            order_test = self.protein.order[:]
            for direction in free_directions_sorted:
                order_test[node_idx] = direction
                if fast_validate_protein(order_test[1:]):
                    # node.change_direction(direction)
                    # print(self.protein.pos_to_node)
                    break

            else:
                # If no direction results in a valid protein, continue to next permutation
                # if self._debug:
                #   print("No direction resulted in valid configuration")
                continue

            score = fast_compute_bond_score(self.protein.sequence, order_test[1:])
            if score <= best_score or threshold > random.random():
                if self._debug and score < best_score:
                    print(f"New best score: {score}, old best score: {best_score}")

                if score < best_score:
                    best_order_overall = order_test

                best_order = order_test
                best_score = score
                self.protein.unghost_all()
                node.change_direction(direction)
                self.protein.unghost_all()

            threshold *= self.decrease
            self._iteration += 1
            progressbar.update(1)
            if self._iteration > self.n_permutations:
                break

        self.protein.set_order(best_order_overall[1:])
        return best_score

        # best_order = []
        # score = 0
        # threshold = 1
        # for _ in tqdm(range(self.n_permutations)):
        #     # get random node idx to permutate from and new directions
        #     node_idx = random.randint(1, len(self.protein.sequence) - 1)
        #     dirs_total = self.get_permutated_directions(node_idx)
        #
        #     # Prevent computing score if order is not valid
        #     if fast_validate_protein(dirs_total):
        #
        #         comparison_score = fast_compute_bond_score(self.protein.sequence, dirs_total)
        #         decision_float = random.random()
        #         if comparison_score <= score or threshold > decision_float:
        #             self.protein.set_order(dirs_total)
        #             score = comparison_score
        #
        #         if comparison_score <= score:
        #             best_order = dirs_total
        #
        #     threshold *= self.decrease
        #
        # self.protein.set_order(best_order)
        # score = self.protein.get_bond_score()
        #
        # return score
