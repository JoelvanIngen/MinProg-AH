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

    def _process_heuristics_single_direction(self,
                                             direction: int,
                                             node_idx: int,
                                             heuristics: callable) -> list[float]:
        """
        For a given node, computes and returns the heuristic scores for a
        single direction.
        """
        self.protein.preserve()
        self.protein.nodes[node_idx].change_direction(direction)
        direction_scores = [heuristic.run() for heuristic in heuristics]
        self.protein.revert()

        return direction_scores

    def _process_heurstics(self, free_directions: list[int],
                           *args) -> list[list[float]]:
        """
        For a given node, computes and returns the heuristic scores for each
        free direction.
        """
        scores = []
        for direction in free_directions:
            direction_scores = self._process_heuristics_single_direction(direction, *args)
            scores.append(direction_scores)

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
