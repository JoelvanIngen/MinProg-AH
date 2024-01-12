import random

from . import Algorithm

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein


class PureRandom(Algorithm):
    """
    A random algorithm that comes up with a random directions list, checks
        whether the configuration is valid, applies the configuration and
        evaluates the score of the configuration. If the order wasn't valid,
        it just comes up with a new random one and applies that.
    """

    def __init__(self, protein: 'Protein', dimensions):
        super().__init__(protein, dimensions)

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
        while True:
            random_order = self._create_order_list()
            self.protein.set_order(random_order)

            if self.protein.has_valid_order():
                break

        score = self.protein.get_bond_score()

        return score


class IterativeRandom(Algorithm):
    """
    A random algorithm that applies a random direction for each individual
        node, starting from the first node after the root node. It determines
        for each node what the available directions are and which directions
        are already occupied. If there are no directions left i.e., it worked
        itself into a corner, it will restart from scratch.
    """

    def __init__(self, protein: 'Protein', dimensions: int):
        super().__init__(protein, dimensions)

        raise NotImplementedError  # Remove when implemented

    def run(self) -> float:
        pass
