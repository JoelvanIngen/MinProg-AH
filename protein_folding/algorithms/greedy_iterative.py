from . import Algorithm
from .heuristics import *
from protein_folding.fast_protein import fast_compute_bond_score

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein

direction_dict = {-1: "Left", 1: "Right", 2: "Up", -2: "Down", 3: "Forward", -3: "Backward"}


class IterativeGreedy(Algorithm):
    """
    A greedy algorithm that applies a direction for each individual
        node, starting from the first node after the root node. It determines
        for each node what the available directions are and which directions
        are already occupied. Then it chooses the direction which yields the
        highest immediate score change. If there are no directions left i.e., it worked
        itself into a corner, it will restart from scratch.
    """

    def __init__(self, protein: 'Protein', dimensions: int, max_iterations: int, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

        self._iteration = 0
        self.max_iterations = max_iterations
        self.lowest_callback = len(protein)

        self.best_score: int = 1
        self.best_order: list[int] = []

        self.amount_of_best_found = 0

        self.dimensions = dimensions

        self.budget = 0
        self.cost_per_iteration = max_iterations / dimensions ** (len(self.protein) - 1)

        self.heuristics = (
            MinimiseDimensions(self.protein),
            FoldAmount(self.protein),
            # Potential(self.protein),
        )

    def saving_by_cutting_branch(self, depth: int):
        return self.dimensions ** (len(self.protein) - depth - 1)

    def _next_fold(self, depth: int):
        # Check if we've reached the end
        if depth >= len(self.protein):

            bond_score = fast_compute_bond_score(self.protein.sequence, self.protein.order[1:])
            if bond_score < self.best_score:
                self.best_score = bond_score
                self.best_order = self.protein.order
                self.amount_of_best_found = 1
            elif bond_score == self.best_score:
                self.amount_of_best_found += 1

            return

        self._iteration += 1

        if self._iteration > self.max_iterations:
            return

        free_directions = self.protein.nodes[depth].get_free_directions(self.directions)
        if not free_directions:
            return

        direction_scores, free_directions_sorted = self._process_heuristics(
            depth, free_directions, self.heuristics)

        # Determine budget and decide on prining if necessary

        for i, direction in enumerate(free_directions_sorted):
            # if i == len(free_directions) - 1 and i != 0:
            #     # Return early to prevent searching for worst scored branch
            #     self._iteration -= 1
            #     return

            self.protein.preserve()
            self.protein.nodes[depth].change_direction(direction)
            self._next_fold(depth + 1)
            self.protein.revert()

            # Update lowest node to see how much of the configuration has been tried
            if depth < self.lowest_callback:
                self.lowest_callback = depth

        if self._debug and self._iteration % 100 == 0:
            print(f'Iteration: {self._iteration}/{self.max_iterations}, Depth/lowest: {depth}/{self.lowest_callback},'
                  f' Best score: {self.best_score} ({self.amount_of_best_found} found)')

    def run(self) -> float:
        # Start at first node after root node
        self._next_fold(depth=1)

        if self._debug:
            print(self.best_order)

        self.protein.set_order(self.best_order[1:])
        return self.best_score
