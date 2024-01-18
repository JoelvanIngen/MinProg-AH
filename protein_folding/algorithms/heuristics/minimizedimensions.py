from . import Heuristic
from protein_folding.protein import Protein
from protein_folding.node import Node
from protein_folding.vector.tools import get_min_max  

class MinimizeDimensions(Heuristic):
    def __init__(self, *args):
        super().__init__(*args)

    def _calc_size_score(self, left_upper, right_lower):
        """
        Computes a score based on the physical dimensions of the protein. This
            value can be used as a penalty for a specific configuration.
        """
        box = right_lower - left_upper
        return box.len_sq()

    def run(self, node_idx: int, directions: list[int]) -> list[float]:
        # Evaluates each direction for a given node and calculates the dimension penalty.
        scores = []

        for direction in directions:
            # Create a test protein with a new direction for the node
            test_protein = Protein(self.protein.sequence)
            test_protein.set_order(self.protein.get_order())
            test_protein.nodes[node_idx].change_direction(direction)

            # Find bounds and calculate penalty
            left_upper, right_lower = get_min_max([node.pos for node in self.nodes])
            penalty = self._calc_size_score(left_upper, right_lower)
            scores.append(penalty)

        return scores