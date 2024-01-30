import random

from . import Algorithm

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein

direction_dict = {-1: "Left", 1: "Right", 2: "Up", -2: "Down", 3: "Forward", -3: "Backward"}


class Greedy(Algorithm):
    """
    A greedy algorithm that applies a direction for each individual
        node, starting from the first node after the root node. It determines
        for each node what the available directions are and which directions
        are already occupied. Then it chooses the direction which yields the
        highest immediate score change. If there are no directions left i.e., it worked
        itself into a corner, it will restart from scratch.
    """

    def __init__(self, protein: 'Protein', dimensions: int, **kwargs):
        super().__init__(protein, dimensions, **kwargs)

    def _attempt_construct_order(self):
        # Order protein in straight line
        self.protein.straighten()

        # Iterate through nodes, check their available directions,
        # choose a random one and apply it
        for node in self.protein.nodes[1:]:
            direction_scores = {}
            free_directions = node.get_free_directions(self.directions)

            if not free_directions:
                return False

            for direction in free_directions:
                self.protein.preserve()
                node.change_direction(direction)

                neighbours = []
                for pos in node.get_neighbouring_positions():
                    other_node = self.protein.pos_to_node.get(pos)
                    if other_node:
                        neighbours.append(node.bond_value(other_node))

                direction_scores[direction] = sum(neighbours)
                # direction_scores[direction] = self.protein.get_bond_score()
                self.protein.revert()

            # results = {}
            # for k, v in direction_scores.items():
            #     results[direction_dict[k]] = v
            # print(results)

            min_value = min(direction_scores.values())
            best_directions = [k for k, v in direction_scores.items() if v == min_value]
            chosen_direction = random.choice(best_directions)
            node.change_direction(chosen_direction)
        return True

    def run(self) -> float:
        failed_attempts = 0
        while True:
            if self._attempt_construct_order():
                break

            failed_attempts += 1

        if self._debug:
            print(f"Failed attempts: {failed_attempts}")

        # If the following fails, self._attempt_construct_order actually
        # failed and should not have returned
        assert self.protein.has_valid_order()

        # Compute and return score
        score = self.protein.get_bond_score()

        return score
