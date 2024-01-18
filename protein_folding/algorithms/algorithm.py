from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein_folding.protein import Protein

# Maps the dimensionality to the list of available directions for the algorithm to take
_dimensions_to_directions_mapping = {
    2: [-2, -1, 1, 2],
    3: [-3, -2, -1, 1, 2, 3]
}


class Algorithm:
    def __init__(self, protein: 'Protein', dimensions: int, debug=False):
        self.protein = protein
        self.directions = _dimensions_to_directions_mapping[dimensions]

        self._debug = debug

    def get_name(self) -> str:
        return self.__class__.__name__

    def _process_heurstics(self, node_idx: int,
                           free_directions: list[int],
                           heuristics: callable) -> list[list[float]]:
        scores = []
        for direction in free_directions:
            self.protein.preserve()
            self.protein.nodes[node_idx].change_direction(direction)
            direction_scores = [heuristic.run() for heuristic in heuristics]
            scores.append(direction_scores)
            self.protein.revert()

        return scores

    def run(self) -> float:
        """
        Runs the algorithm on the protein and calculates and returns
            the final score of the optimised protein

        post:
            - A float is returned representing the score of the
            optimised protein configuration
        """
        pass
