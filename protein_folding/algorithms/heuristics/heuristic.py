from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein_folding.protein import Protein
    from protein_folding.node import Node


class Heuristic:
    def __init__(self, protein: Protein, weight: float = 1.):
        """
        Initialises the heuristic

        pre:
            - protein is a protein object to run on
            - weight is a float indicating how the heuristic should weigh its
            final output
        """
        self.protein = protein
        self.weight = weight

    def run(self, node: Node, directions: list[int]) -> list[float]:
        """
        Runs the heuristic and computes a score for each possible option

        pre:
            - node is a Node object to apply the heuristic on
            - directions is a list of directions that a node can take from
            its source (the previous node).
        post:
            - a float is returned representing a score for every direction.
        """
        pass
