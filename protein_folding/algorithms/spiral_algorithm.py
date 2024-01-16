import random

from . import Algorithm

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein_folding.protein import Protein


class Circle(Algorithm):
    def __init__(self, protein: 'Protein', dimensions, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

    def run(self):
        iterations = len(self.protein) - 2
        order = [2]
        for i in range(iterations):
            i += 1
