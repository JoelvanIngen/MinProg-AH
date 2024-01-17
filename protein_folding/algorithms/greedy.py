import random
from . import Algorithm
from typing import TYPE_CHECKING
from copy import deepcopy

if TYPE_CHECKING:
    from protein_folding.protein import Protein
    from protein_folding.node import Node

direction_dict = {1: "Left", -1: "Right", 2: "Up", -2: "Down"}


class Greedy(Algorithm):
    """
    A greedy algorithm that applies a direction for each individual
        node, starting from the first node after the root node. It determines
        for each node what the available directions are and which directions
        are already occupied. Then it chooses the direction which yields the
        highest immediate score change. If there are no directions left i.e., it worked
        itself into a corner, it will restart from scratch.
    """

    def __init__(self, protein: 'Protein', dimensions: int, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

    def _attempt_construct_order(self):
        # Order protein in straight line
        self.protein.straighten()

        # Iterate through nodes, check their available directions,
        # choose a random one and apply it
        for node in self.protein.nodes[1:]:
            current_direction = node.direction_from_previous
            direction_scores = {}
            free_directions = node.get_free_directions(self.directions)
            if not free_directions:
                print("No free directions found")
                return False

            for direction in free_directions:
                node.change_direction(direction)
                if self.protein.has_valid_order():
                    direction_scores[direction] = self.protein.get_bond_score()
                # print(f"Direction checked: {direction_dict[direction]}")
                self.protein.revert()

            print(f"Direction_scores: {direction_scores}")
            min_value = min(direction_scores.values())
            best_directions = [k for k, v in direction_scores.items() if v == min_value]
            print(f"Best directions: {best_directions}")
            direction = random.choice(best_directions)
            print(f"Direction selected: {direction_dict[direction]}")
            node.change_direction(direction)
        return True

    def run(self) -> float:
        failed_attempts = 0
        while True:
            if self._attempt_construct_order():
                break

            failed_attempts += 1

        if self._debug:
            print(f"Failed attempts: {failed_attempts}")

        # If the following fails, self._attempt_construct_order actually
        # failed and should not have returned
        assert self.protein.has_valid_order()

        # Compute and return score
        score = self.protein.get_bond_score()

        return score
