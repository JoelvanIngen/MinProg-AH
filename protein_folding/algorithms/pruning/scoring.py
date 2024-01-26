import math

from protein_folding.fast_protein import fast_compute_bond_score
from . import Pruning


class Score(Pruning):
    def __init__(self, protein, alpha=0.5, beta=15, use_linear=False):
        super().__init__(protein)

        self.protein_length = len(protein)

        # I called the parameters alpha and beta because I see other algorithms
        # use those letters. There is no further scientific reason.
        # Alpha: Affect scaling speed. Higher alpha: faster scaling
        #   -> more pruning
        # Beta: Affects fixed margin. Higher beta:
        #   -> less pruning
        self.alpha = alpha
        self.beta = beta
        self.use_linear = use_linear

    def run(self, *, best_score: int, depth: int) -> bool:
        """
        Counts the amount of nodes with no neighbours and decides whether the
            configuration should be pruned. Returns a bool indicating the
            decision. For this agorithm to return true there need to be a
            number of criteria:
                - The amount of non-ghosted nodes should be more than 6
                - The fraction of non-ghosted nodes with no neighbours other
                    than the nodes directly next to it in the chain must
                    be more than half
        """

        # Prevent pruning if we do not yet have a good score base, to prevent
        # pruning too aggressivel
        if depth < self.protein_length / 3 + 2 or best_score > -5:
            return False

        score = -fast_compute_bond_score(self.protein.sequence[:depth + 1], self.protein.order[1:depth + 1])

        if self.use_linear:
            threshold = self.protein_length * depth / (self.beta * -best_score)
        else:
            threshold = -best_score * math.exp(1 / self.protein_length * depth * self.alpha) - self.beta

        return score < threshold
