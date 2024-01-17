from . import Heuristic

from protein_folding.protein import Protein
from protein_folding.node import Node


class Potential(Heuristic):
    def __init__(self, *args):
        super().__init__(*args)

    def _find_new_positions(self, node_idx: int, direction: int):
        # Copy protein to experiment on copy
        test_protein = Protein(self.protein.sequence)
        test_protein.set_order(self.protein.get_order())

        # Set new direction
        test_protein.nodes[node_idx].change_direction(test_protein.node_positions, direction)

        return [node for node in test_protein.nodes[node_idx + 1:] if node.letter != 'P']

    def _find_sources(self, node_idx: int) -> list[Node]:
        return [node for node in self.protein.nodes[:node_idx] if node.letter != 'P']

    def _find_direction_score(self, sources: list[Node], receivers: list[Node]):
        # Higher score is better, because nodes closer to eachother
        score = 0.

        for source in sources:
            for receiver in receivers:
                delta_pos = source.pos - receiver.pos
                r_sq = delta_pos.len_sq()

                score += 1 / r_sq

        return score

    def run(self, node_idx: int, directions: list[int]) -> list[float]:
        sources = self._find_sources(node_idx)

        receivers_per_direction = []
        for direction in directions:
            receivers = self._find_new_positions(node_idx, direction)

            receivers_per_direction.append(receivers)

        scores = []
        for receivers in receivers_per_direction:
            scores.append(self._find_direction_score(sources, receivers))

        return scores
