from math import floor
from collections import deque

from . import Algorithm

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein


class Spiral(Algorithm):
    def __init__(self, protein: 'Protein', dimensions, **kwargs):
        super().__init__(protein, dimensions, **kwargs)
        assert dimensions == 2, "Spiral only works in 2D"

    def create_repetitions(self, num_of_nodes: int) -> list:
        repetitions = []
        i = 2
        while sum(repetitions) < num_of_nodes:
            repetitions.append(floor(i / 2))
            i += 1
        if sum(repetitions) > num_of_nodes:
            repetitions[-1] = num_of_nodes - sum(repetitions[:-1])
        return repetitions

    def generate_order(self) -> list:
        rep = deque(self.create_repetitions(len(self.protein) - 1))
        order = []
        while rep:
            for index in [0, 1, 3, 2]:
                if rep:
                    repetition = rep.popleft()
                    order += (repetition * [self.directions[index]])
        return order

    def run(self):
        order = self.generate_order()
        self.protein.set_order(order)
        score = self.protein.get_bond_score()

        return score
