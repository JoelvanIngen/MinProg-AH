from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein


class Pruning:
    def __init__(self, protein: 'Protein'):
        self.protein = protein

    def run(self, **kwargs) -> bool:
        """
        Checks an order for specific requirements and returns a bool indicating
            whether the branch should be pruned or not.
        """

        pass
