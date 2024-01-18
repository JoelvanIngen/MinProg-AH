import random
from . import Algorithm
from typing import TYPE_CHECKING
from copy import deepcopy

from protein_folding.protein import Protein
from protein_folding.node import Node

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

    def _next_fold(self, node: Node | None, depth: int):
        self._iteration += 1

        if self._iteration > self.max_iterations:
            return

        if self._debug:
            print(f'Iteration: {self._iteration}/{self.max_iterations}, Depth: {depth}, Best score: {self.best_score}')

        if not node:
            # Cannot be false if code works correctly, since we only look for
            # free directions and revert if stuck, reaching the chain end means
            # having a valid order
            assert self.protein.has_valid_order()

            if self.protein.get_bond_score() < self.best_score:
                self.best_score = self.protein.get_bond_score()
                self.best_order = self.protein.order
                self.best_order_deepcopy = deepcopy(self.protein.order)

            return

        free_directions = node.get_free_directions(self.directions)

        for direction in free_directions:
            self.protein.preserve()
            node.change_direction(direction)
            self._next_fold(node.next, depth + 1)
            self.protein.revert()

    def run(self) -> float:
        start_node = self.protein.nodes[1]

        self._next_fold(start_node, depth=0)

        print(self.best_order)
        print(self.best_order_deepcopy)
        self.protein.set_order(self.best_order[1:])
        return self.best_score

