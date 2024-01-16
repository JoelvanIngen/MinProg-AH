import random

from . import Algorithm
from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein


class PureRandom(Algorithm):
    """
    A random algorithm that comes up with a random directions list, checks
        whether the configuration is valid, applies the configuration and
        evaluates the score of the configuration. If the order wasn't valid,
        it just comes up with a new random one and applies that.

    This algorithm should be quick when constructing a protein, but might
        need a lot of attempts before getting a valid protein, decreasing
        performance when working with large structures
    """

    def __init__(self, protein: 'Protein', dimensions, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

    def _get_new_direction(self, previous_direction) -> int:
        while True:
            direction = random.choice(self.directions)
            if direction != -previous_direction:
                return direction

    def _create_order_list(self) -> list[int]:
        order_length = len(self.protein.sequence) - 1

        # Initialise random order list
        order: list[int] = [random.choice(self.directions)]

        # Fill other spots with random directions
        for _ in range(1, order_length):
            order.append(self._get_new_direction(order[-1]))

        return order

    def run(self) -> float:
        failed_attempts = 0
        while True:
            random_order = self._create_order_list()

            if fast_validate_protein(random_order):
                break

            failed_attempts += 1

        if self._debug:
            print(f"Failed attempts: {failed_attempts}")

        score = fast_compute_bond_score(self.protein.sequence, random_order)
        self.protein.set_order(random_order)

        return score


class IterativeRandom(Algorithm):
    """
    A random algorithm that applies a random direction for each individual
        node, starting from the first node after the root node. It determines
        for each node what the available directions are and which directions
        are already occupied. If there are no directions left i.e., it worked
        itself into a corner, it will restart from scratch.

    I expect this algorithm to be slower to construct a protein, but to have a
        higher success rate (not overlapping) when working with large proteins.
    """

    def __init__(self, protein: 'Protein', dimensions: int, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

    def _attempt_construct_order(self):
        # Order protein in straight line
        self.protein.straighten()

        # Iterate through nodes, check their available directions,
        # choose a random one and apply it
        for node in self.protein.nodes[1:]:
            free_directions = node.get_free_directions(self.protein.node_positions, self.directions)
            if not free_directions:
                return False

            direction = random.choice(free_directions)
            node.change_direction(self.protein.node_positions, direction)

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
