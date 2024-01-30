from tqdm import tqdm

from protein_folding.fast_protein import fast_compute_bond_score, fast_validate_protein
from . import Algorithm
from .pruning import *

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein


class QueueEmptyError(Exception):
    pass


class StateQueue:
    def __init__(self):
        self.scores: list[int] = [0]
        self.states: list[list[int]] = [[1]]

    def get_highest_score_state(self) -> list[int]:
        """
        Finds and returns state with the highest score so far
        """
        if len(self.scores) == 0:
            assert len(self.states) == 0
            raise QueueEmptyError

        max_score = max(self.scores)
        idx = self.scores.index(max_score)

        self.scores.pop(idx)
        return self.states.pop(idx)

    def push(self, order: list[int], score: int) -> None:
        self.scores.append(score)
        self.states.append(order)


class BeamSearch(Algorithm):
    """
    A greedy algorithm that applies a direction for each individual
        node, starting from the first node after the root node. It determines
        for each node what the available directions are and which directions
        are already occupied. Then it chooses the direction which yields the
        highest immediate score change. If there are no directions left i.e., it worked
        itself into a corner, it will restart from scratch.
    """

    def __init__(self, protein: 'Protein', dimensions: int, max_iterations: int = 5000, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

        self._iteration = 0
        self.max_iterations = max_iterations

        self.best_score: int = 1
        self.best_order: list[int] = []

        self.end_nodes_evaluated = 0
        self.amount_of_best_found = 0

        self.queue = StateQueue()

    def _increment_iteration(self):
        if self.pbar:
            self.pbar.update(1)

        # UNUSED
        if self.keep_score_history:
            self.best_score_history.append(self.best_score)

        # UNUSED
        if self.keep_order_history:
            if fast_validate_protein(self.protein.order):
                self.order_history.append(self.protein.order)

    def find_next_orders(self, order) -> list[list[int]]:
        valid_next_orders = []
        for direction in self.directions:
            if direction == order[-1]:
                continue

            test_order = order + [direction]
            if not fast_validate_protein(test_order):
                continue

            valid_next_orders.append(test_order)

        return valid_next_orders

    def get_best_order(self):
        return self.queue.get_highest_score_state()

    def save_orders(self, orders):
        for order in orders:
            self.save_order(order)

    def save_order(self, order):
        score = fast_compute_bond_score(self.protein.sequence[:len(order)], order)
        self.queue.push(order, score)

    def keep_best_scores(self, orders):
        for order in orders:
            score = fast_compute_bond_score(self.protein.sequence[:len(order)], order)
            if score < self.best_score:
                if self._debug:
                    print(f"Found new best score of {score} with order: {order}")

                self.best_score = score
                self.best_order = order[:]

    def process_order(self, order):
        next_orders = self.find_next_orders(order)
        if not next_orders:
            return

        depth = len(next_orders[0])
        if depth == len(self.protein) - 1:
            return self.keep_best_scores(next_orders)

        self.save_orders(next_orders)

        # if self._debug and self._iteration % 1000 == 0:
        #     print(depth)

    def run(self) -> int:
        if self._show_progress:
            self.pbar = tqdm(range(self.max_iterations))

        # Start at first node after root node
        for iteration in range(self.max_iterations):
            try:
                order = self.get_best_order()
            except QueueEmptyError:
                break

            self.process_order(order)

            self._increment_iteration()

        print(self.best_order)
        self.protein.set_order(self.best_order)
        return self.best_score

    # def _next_fold(self, depth: int):
    #     def _reached_protein_end() -> bool:
    #         # Check if we've reached the end of the protein and keep results
    #         if depth >= len(self.protein):
    #             self.end_nodes_evaluated += 1
    #
    #             bond_score = fast_compute_bond_score(self.protein.sequence, self.protein.order[1:])
    #             if bond_score < self.best_score:
    #                 self.best_score = bond_score
    #                 self.best_order = self.protein.order
    #                 self.amount_of_best_found = 1
    #             elif bond_score == self.best_score:
    #                 self.amount_of_best_found += 1
    #
    #             return True
    #
    #         else:
    #             return False
    #
    #     def _prune_score() -> bool:
    #         pruner = Score(self.protein, alpha=self.prune_alpha, beta=self.prune_beta)
    #         return pruner.run(best_score=self.best_score, depth=depth)
    #
    #     def _reached_max_iterations() -> bool:
    #         return self._iteration > self.max_iterations
    #
    #     if _prune_score():
    #         self.amount_pruned += 1
    #         return
    #
    #     if _reached_protein_end() or _reached_max_iterations():
    #         return
    #
    #     self._increment_iteration()
    #
    #     free_directions = self.protein.nodes[depth].get_free_directions(self.directions)
    #     if not free_directions:
    #         return
    #
    #     direction_scores, free_directions_sorted = self._process_heuristics(
    #         depth, free_directions)
    #
    #     for i, direction in enumerate(free_directions_sorted):
    #         self.protein.preserve()
    #         self.protein.nodes[depth].change_direction(direction)
    #         self._next_fold(depth + 1)
    #         self.protein.revert()
    #
    #         # Update lowest-reached node to see how much of the configuration has been tried
    #         if depth < self.lowest_callback:
    #             self.lowest_callback = depth
    #         if depth < self.batch_lowest_callback:
    #             self.batch_lowest_callback = depth
    #
    #     if self._debug and self._iteration % 1000 == 0:
    #         print(f'Iteration: {self._iteration}/{self.max_iterations}, Depth/batch lowest/lowest: '
    #               f'{depth}/{self.batch_lowest_callback}/{self.lowest_callback},'
    #               f' Best score: {self.best_score} ({self.amount_of_best_found} found), End nodes reached: '
    #               f'{self.end_nodes_evaluated}, Amount pruned: {self.amount_pruned}')
    #
    #         self.batch_lowest_callback = len(self.protein)
