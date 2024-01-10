from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .protein import Protein


class Algorithm:
    def __init__(self, protein: Protein):
        self._protein = protein

    def optimise(self):
        pass
