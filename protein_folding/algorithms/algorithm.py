from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein_folding.protein import Protein

# Maps the dimensionality to the list of available directions for the algorithm to take
_dimensions_to_directions_mapping = {
    2: [-2, -1, 1, 2],
    3: [-3, -2, -1, 1, 2, 3]
}


class Algorithm:
    def __init__(self, protein: 'Protein', dimensions: int,
                 heuristics=(),
                 debug=False,
                 visualise=False,
                 show_progressbar=True,
                 keep_score_history=False):

        self.protein = protein
        self.directions = _dimensions_to_directions_mapping[dimensions]

        self.heuristics = [
            h(protein) for h in heuristics
        ]

        self.keep_score_history = keep_score_history
        self.best_score_history = []

        self._debug = debug

        self._visualise = visualise
        self._order_history = []

        self._show_progress = show_progressbar
        self.pbar = None

    @property
    def name(self):
        return self.get_name()

    def get_name(self) -> str:
        return self.__class__.__name__

    def _keep_order_history(self, order: list[int] | None):
        if self._visualise:
            self._order_history.append(order[:] if order else self.protein.order[:])

    def _process_heuristics_single_direction(self,
                                             node_idx: int,
                                             direction: int,
                                             heuristics: list[callable],
                                             unghost: bool) -> list[float]:
        """
        For a given node, computes and returns the heuristic scores for a
        single direction.
        """
        self.protein.preserve()
        self.protein.nodes[node_idx].change_direction(direction)
        if unghost:
            self.protein.unghost_all()
        direction_scores = [heuristic.run() for heuristic in heuristics]
        self.protein.revert()

        return direction_scores

    def _process_heuristics(self,
                            node_idx: int,
                            free_directions: list[int],
                            unghost: bool = False) -> tuple[list[list[float]], list[int]]:
        """
        For a given node, computes and returns the heuristic scores for each
            free direction. Sorts the scores and the corresponding directions
            and returns both in a tuple.
        """
        if not self.heuristics:
            return [0 for _ in free_directions], free_directions

        # Remove previously remembered heuristic scores
        for heuristic in self.heuristics:
            heuristic.reset()

        # Process heuristic for every free dimension
        for direction in free_directions:
            self._process_heuristics_single_direction(node_idx, direction, self.heuristics, unghost)

        # Process data for each heuristic
        heuristic_scores = [0 for _ in range(len(free_directions))]
        for heuristic in self.heuristics:
            res = heuristic.interpret()
            for i, val in enumerate(res):
                heuristic_scores[i] += val

        # Sort directions based on their score from best to worst
        sorted_pairs = sorted(zip(heuristic_scores, free_directions), reverse=True)
        scores_sorted = [x for x, _ in sorted_pairs]
        directions_sorted = [x for _, x in sorted_pairs]

        return scores_sorted, directions_sorted

    def run(self) -> float:
        """
        Runs the algorithm on the protein and calculates and returns
            the final score of the optimised protein

        post:
            - A float is returned representing the score of the
            optimised protein configuration
        """
        pass

    def plot_score_progress(self):
        import matplotlib.pyplot as plt
        plt.plot([*range(len(self.best_score_history))], self.best_score_history)
        plt.show()
