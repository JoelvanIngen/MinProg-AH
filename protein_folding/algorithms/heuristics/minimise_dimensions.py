from . import Heuristic
from protein_folding.vector import Vec3D, get_min_max


class MinimiseDimensions(Heuristic):
    def __init__(self, *args):
        super().__init__(*args)

    def _calc_size_score(self, left_upper: Vec3D, right_lower: Vec3D):
        """
        Computes a score based on the physical dimensions of the protein. This
            value can be used as a penalty for a specific configuration.
        """
        box = right_lower - left_upper
        return box.len_sq()

    def run(self, **kwargs):
        # Evaluates the protein and calculates the dimension penalty.
        left_upper, right_lower = get_min_max([node.pos for node in self.protein.nodes if not node.ghost])

        self.score_per_direction.append(self._calc_size_score(left_upper, right_lower))

    def interpret(self) -> list[float]:
        """
        Normalises the values such that max(scores) = 1,
            and flips them such that a score of 1 is better than a score of 0.
        """
        _min = min(self.score_per_direction)
        _max = max(self.score_per_direction)
        delta = _max - _min

        if delta == 0:
            # Avoid division by zero if all scores are similar
            return [0. for _ in self.score_per_direction]

        scores_norm = [1 / (val / _min) for val in self.score_per_direction]
        return scores_norm
