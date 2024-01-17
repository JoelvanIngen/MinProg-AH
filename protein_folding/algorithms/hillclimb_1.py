import random

from . import Algorithm
from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from protein_folding.protein import Protein

class PureRandomHillClimb(Algorithm):
    def __init__(self, protein: 'Protein', dimensions, max_iterations=1000, **kwargs):
        super().__init__(protein, dimensions, **kwargs)
        self.max_iterations = max_iterations
        
    def _create_order_list(self) -> list[int]:
        order_length = len(self.protein.sequence) - 1
        order = [random.choice(self.directions) for _ in range(order_length)]
        return order
    
    def _modify_order(self, order):
        # Randomly choose an index to modify
        index_to_modify = random.randint(0, len(order) - 1)
        # Choose a new direction, ensuring it's different from the current one
        new_direction = random.choice([d for d in self.directions if d != order[index_to_modify]])
        order[index_to_modify] = new_direction
        return order

    def run(self) -> float:
        order = self._create_order_list()
        best_score = fast_compute_bond_score(self.protein.sequence, order)
        self.protein.set_order(order)

        for _ in range(self.max_iterations):
            new_order = self._modify_order(order.copy())
            if fast_validate_protein(new_order):
                new_score = fast_compute_bond_score(self.protein.sequence, new_order)
                if new_score > best_score:
                    best_score = new_score
                    order = new_order
                    self.protein.set_order(order)

        return best_score


class IterativeRandomHillClimb(Algorithm):
    def __init__(self, protein: 'Protein', dimensions, max_iterations=1000, **kwargs):
        super().__init__(protein, dimensions, **kwargs)
        self.max_iterations = max_iterations

    def _try_improve_node(self, node):
        original_direction = node.direction
        best_direction = original_direction
        best_score = self.protein.get_bond_score()

        for direction in self.directions:
            if direction != original_direction:
                node.direction = direction
                self.protein.collect_node_positions()  # Update positions after direction change
                new_score = self.protein.get_bond_score()
                if new_score > best_score:
                    best_score = new_score
                    best_direction = direction

        # Revert to original direction if no improvement found
        if best_direction == original_direction:
            return False
        else:
            node.direction = best_direction
            self.protein.collect_node_positions()  # Update positions for new best direction
            return True

    def run(self) -> float:
        self.protein.straighten()
        best_score = self.protein.get_bond_score()

        for _ in range(self.max_iterations):
            improvement = False
            for node in self.protein.nodes[1:]:
                if self._try_improve_node(node):
                    improvement = True

            if not improvement:
                break  # Terminate if no improvements can be made

        return best_score
