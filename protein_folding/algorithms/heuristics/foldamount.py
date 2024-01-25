from . import Heuristic


class FoldAmount(Heuristic):
    """
    A heuristic that determines the amount of non-ordered links ('corners',
        'folds') in the protein, with the philosophy that less ordered links
        ('straight sections') will result in more compactness and thus a
        higher chance that a bond is formed.

    Returns the amount of non-ordered links for the protein.
    """
    def __init__(self, *args):
        super().__init__(*args)

    def run(self, **kwargs):
        # Evaluates the protein and calculates the number of non-ordered links.
        n_corners = 0

        for node in [node for node in self.protein.nodes[2:] if not node.ghost]:
            if node.direction_from_previous != node.prev.direction_from_previous:
                n_corners += 1

        self.score_per_direction.append(n_corners)

    def interpret(self) -> list[float]:
        """
        Normalises the values such that max(scores) = 1 and min(scores) > 0.
        """
        _min = min(self.score_per_direction)
        _max = max(self.score_per_direction)

        if _max == 0:
            # Avoid division by zero if all scores are similar
            return [0. for _ in self.score_per_direction]

        scores_norm = [value / _max for value in self.score_per_direction]
        return scores_norm
