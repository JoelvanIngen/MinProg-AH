from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein


class Heuristic:
    def __init__(self, protein: 'Protein'):
        """
        Initialises the heuristic

        pre:
            - protein is a protein object to run on
            - weight is a float indicating how the heuristic should weigh its
            final output
        """
        self.protein = protein

        self.score_per_direction = []

    @property
    def name(self) -> str:
        return self.__class__.__name__

    def run(self, **kwargs):
        """
        Runs the heuristic and computes a score for a protein configuration

        Post:
            - a float value representing the heuristic's evaluation of the
                current order is saved in self.score_per_direction
        """
        pass

    def interpret(self) -> list[float]:
        """
        Reads the scores the heuristic assigned to different positions, and
            interprets them for the algorithm. It can, for example, normalise
            values or invert them.

        In the return values, a higher score means a better position.
        """
        pass

    def reset(self) -> None:
        """
        Resets the heuristic by forgetting previous values per direction
        """

        self.score_per_direction = []
