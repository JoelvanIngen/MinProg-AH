import random
import math
import traceback
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein_folding.protein import Protein

class HillClimb:
    """
    Hill Climbing algorithm for protein folding, using set_order for generating neighboring configurations
    and existing evaluation methods for evaluating configurations.
    """

    def __init__(self, protein: 'Protein', dimensions: int, temperature: float = 1000.0, cooling_rate: float = 0.001):
        self.protein = protein
        self.temperature = temperature
        self.cooling_rate = cooling_rate
        self.dimensions = dimensions
                         
    def run(self):
        print("Starting hill climbing algorithm")
        # Initial random configuration
        self.protein.random_configuration()
        best_score = self.protein.calc_size_score() # or calc_area(), calc_volume() depending on the evaluation criteria

        while self.temperature > 1:
            print(f"Current temperature: {self.temperature}, Best score: {best_score}")
            # Generate a neighbor configuration
            previous_order = self.protein.get_order()
            neighbor_order = self._generate_neighbor_order()
            if self._is_valid_order(neighbor_order):
                try:
                    self.protein.set_order(neighbor_order)
                except KeyError as e:
                    print(f"Error encountered when setting order: {e}. Reverting to previous configuration.\nFull error details: {traceback.format_exc()}")
                    self.protein.set_order(previous_order) 

            # Evaluate the new configuration
            neighbor_score = self.protein.calc_size_score()  # or calc_area(), calc_volume()

            # Calculate probability of accepting the new configuration
            if neighbor_score < best_score:
                accept_probability = 1
            else:
                accept_probability = math.exp((best_score - neighbor_score) / self.temperature)

            # Determine if the neighbor configuration should be accepted
            if random.uniform(0, 1) < accept_probability:
                best_score = neighbor_score
            else:
                # Revert to previous configuration if not accepted
                self.protein.set_order(previous_order) 

            # Cooling down process
            self.temperature *= 1 - self.cooling_rate

        return self.protein

    def _is_valid_order(self, order):
        # Placeholder for order validation logic
        # For now, it just returns True. This can be expanded to check for specific criteria
        return self._validate_order(order)
    
    def _validate_order(self, order):
        # TODO: Implement order validation logic
        return True

    def _generate_neighbor_order(self):
        # Generate a slightly modified order based on the current order
        current_order = self.protein.get_order()
        neighbor_order = current_order.copy()

        # Example: Modify one random element in the order
        random_index = random.randint(0, len(neighbor_order) - 1)
        direction_options = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]  # 2D directions
        neighbor_order[random_index] = random.choice(direction_options)

        return neighbor_order
