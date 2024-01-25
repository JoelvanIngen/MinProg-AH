from . import Pruning
from protein_folding.fast_protein import fast_compute_bond_score


class Score(Pruning):
    def __init__(self, protein, mult):
        super().__init__(protein)

        # Higher multiplier: less likely to prune
        self.mult = mult

    def run(self, *, best_score, depth) -> bool:
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

        non_ghosted = [node for node in self.protein.nodes if not node.ghost]
        amount_non_ghosted = len(non_ghosted)
        if depth < 8 or best_score > -2:
            return False

        score = -fast_compute_bond_score(self.protein.sequence[:depth + 1], self.protein.order[1:depth + 1])

        threshold = len(self.protein) * depth / (self.mult * -best_score)

        return score < threshold
