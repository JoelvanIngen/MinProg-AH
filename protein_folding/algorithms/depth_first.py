from tqdm import tqdm

from protein_folding.fast_protein import fast_compute_bond_score, fast_validate_protein
from . import Algorithm
from .pruning import *

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein


class DepthFirst(Algorithm):
    """
    A greedy algorithm that applies a direction for each individual
        node, starting from the first node after the root node. It determines
        for each node what the available directions are and which directions
        are already occupied. Then it chooses the direction which yields the
        highest immediate score change. If there are no directions left i.e., it worked
        itself into a corner, it will restart from scratch.
    """

    def __init__(self, protein: 'Protein', dimensions: int, max_iterations: int = 5000,
                 prune_alpha=0.5, prune_beta=15, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

        self._iteration = 0
        self.max_iterations = max_iterations
        self.lowest_callback = len(protein)

        self.best_score: int = 1
        self.best_order: list[int] = []

        self.end_nodes_evaluated = 0
        self.amount_of_best_found = 0

        self.dimensions = dimensions

        self.prune_alpha = prune_alpha
        self.prune_beta = prune_beta

    def _increment_iteration(self):
        self._iteration += 1

        if self.pbar:
            self.pbar.update(1)

        if self.keep_score_history:
            self.best_score_history.append(self.best_score)

        if self.keep_order_history:
            if fast_validate_protein(self.protein.order[1:]):
                self.order_history.append(self.protein.order[1:])

    def _next_fold(self, depth: int):
        def _reached_protein_end() -> bool:
            # Check if we've reached the end of the protein and keep results
            if depth >= len(self.protein):
                self.end_nodes_evaluated += 1

                bond_score = fast_compute_bond_score(self.protein.sequence, self.protein.order[1:])
                if bond_score < self.best_score:
                    self.best_score = bond_score
                    self.best_order = self.protein.order
                    self.amount_of_best_found = 1
                elif bond_score == self.best_score:
                    self.amount_of_best_found += 1

                return True

            else:
                return False

        def _prune_score() -> bool:
            pruner = Score(self.protein, alpha=self.prune_alpha, beta=self.prune_beta)
            return pruner.run(best_score=self.best_score, depth=depth)

        def _reached_max_iterations() -> bool:
            return self._iteration > self.max_iterations

        if _reached_protein_end() or _prune_score() or _reached_max_iterations():
            return

        self._increment_iteration()

        free_directions = self.protein.nodes[depth].get_free_directions(self.directions)
        if not free_directions:
            return

        direction_scores, free_directions_sorted = self._process_heuristics(
            depth, free_directions)

        for i, direction in enumerate(free_directions_sorted):
            self.protein.preserve()
            self.protein.nodes[depth].change_direction(direction)
            self._next_fold(depth + 1)
            self.protein.revert()

            # Update lowest-reached node to see how much of the configuration has been tried
            if depth < self.lowest_callback:
                self.lowest_callback = depth

        if self._debug and self._iteration % 200 == 0:
            print(f'Iteration: {self._iteration}/{self.max_iterations}, Depth/lowest: {depth}/{self.lowest_callback},'
                  f' Best score: {self.best_score} ({self.amount_of_best_found} found), End nodes reached: '
                  f'{self.end_nodes_evaluated}')

    def run(self) -> float:
        if self._show_progress:
            self.pbar = tqdm(range(self.max_iterations))

        # Start at first node after root node
        self._next_fold(depth=1)

        if self._debug:
            print(self.best_order)

        self.protein.set_order(self.best_order[1:])
        return self.best_score
