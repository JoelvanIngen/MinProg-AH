from . import Heuristic

from protein_folding.protein import Protein
from protein_folding.node import Node


# Extra factor to multiply results with so they're not all in range {0.9, 1.]
_MULT_FACTOR = 7


class PotentialPlus(Heuristic):
    def __init__(self, *args):
        super().__init__(*args)

    def run(self):
        score = 0
        for i, node in enumerate(self.protein.nodes[3:], start=3):
            if node.ghost:
                break

            start_j = (i + 1) % 2
            for other_node in self.protein.nodes[start_j:i-1:2]:
                if other_node.ghost:
                    break

                match node.letter + other_node.letter:
                    case 'CC':
                        mult = 4
                    case 'HH' | 'CH' | 'HC':
                        mult = 1
                    case 'PP':
                        mult = 1
                    case 'PH' | 'HP':
                        mult = -1
                    case 'PC' | 'CP':
                        mult = -4
                    case _:
                        return Exception("We forgot to implement a letter combination")

                delta_vec = other_node.pos - node.pos
                len_sq = delta_vec.len_sq()
                score += (1 / len_sq) * mult

        self.score_per_direction.append(score)

    def interpret(self) -> list[float]:
        """
        Normalises all scores such that max(scores) == 1
        """
        _min = min(self.score_per_direction)

        # Ensure all scores are positive
        scores_positive = [val - _min for val in self.score_per_direction]  # [0, ->}

        _max = max(scores_positive)

        # Ensure no dividing by zero
        if _max == 0:
            return [0. for _ in self.score_per_direction]

        # Normalise all scores such that scores:  [0, 1]
        return [val / _max for val in scores_positive]

