from . import Pruning


class Neighbours(Pruning):
    def __init__(self, *args):
        super().__init__(*args)

    def has_neighbours(self, node, others):
        for other in others:
            if abs(other.id - node.id) <= 1:
                continue

            if node.is_neighbour(other):
                return True

        return False

    def _count_nodes_with_no_neighbours(self, nodes):
        count = 0

        for node in nodes:
            if not self.has_neighbours(node, nodes):
                count += 1

        return count

    def run(self) -> bool:
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
        if amount_non_ghosted < 6:
            return False

        # print(self.protein.order, self._count_nodes_with_no_neighbours(non_ghosted), int(amount_non_ghosted / 2))

        if self._count_nodes_with_no_neighbours(non_ghosted) / amount_non_ghosted > 0.75:
            return True

        return False
