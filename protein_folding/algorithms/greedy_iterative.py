from copy import deepcopy

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

        self.best_score: int = 1
        self.best_order: list[int] = []

        # Test if deepcopy is actually necessary
        self.best_order_deepcopy: list[int] = []

        self.heuristics = (
            MinimiseDimensions(self.protein),
            FoldAmount(self.protein),
        )

    def _next_fold(self, depth: int):
        # Check if we've reached the end
        if depth >= len(self.protein):
            bond_score = fast_compute_bond_score(self.protein.sequence, self.protein.order[1:])
            if bond_score < self.best_score:
                self.best_score = bond_score
                self.best_order = self.protein.order
                self.best_order_deepcopy = deepcopy(self.protein.order)

            return

        self._iteration += 1

        if self._iteration > self.max_iterations:
            return

        if self._debug and self._iteration % 100 == 0:
            print(f'Iteration: {self._iteration}/{self.max_iterations}, Depth: {depth}, Best score: {self.best_score}')

        free_directions = self.protein.nodes[depth].get_free_directions(self.directions)

        if not free_directions:
            return

        direction_scores, free_directions_sorted = self._process_heurstics(
            depth, free_directions, self.heuristics)

        for direction in free_directions_sorted:
            self.protein.preserve()
            self.protein.nodes[depth].change_direction(direction)
            self._next_fold(depth + 1)
            self.protein.revert()

    def run(self) -> float:
        # Start at first node after root node
        self._next_fold(depth=1)

        if self._debug:
            print(self.best_order)
            print(self.best_order_deepcopy)

        self.protein.set_order(self.best_order[1:])
        return self.best_score
