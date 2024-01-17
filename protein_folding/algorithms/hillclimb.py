import random
import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein_folding.protein import Protein

class AdaptedHillClimb:
    """
    Hill Climbing algorithm for protein folding, using set_order for generating neighboring configurations
    and existing evaluation methods for evaluating configurations.
    """

    def __init__(self, protein: 'Protein', temperature: float = 1000.0, cooling_rate: float = 0.001):
        self.protein = protein
        self.temperature = temperature
        self.cooling_rate = cooling_rate
    
    def run(self):
        # Initial random configuration
        self.protein.random_configuration()
        best_score = self.protein.calc_size_score() # or calc_area(), calc_volume() depending on the evaluation criteria

        while self.temperature > 1:
            # Generate a neighbor configuration
            neighbor_order = self._generate_neighbor_order()
            self.protein.set_order(neighbor_order)

            # Evaluate the new configuration
            neighbor_score = self.protein.calc_size_score() # or calc_area(), calc_volume()

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
                self.protein.set_order_to_previous()

            # Cooling down process
            self.temperature *= 1 - self.cooling_rate

        return self.protein

    def _generate_neighbor_order(self):
        # Generate a slightly modified order based on the current order
        current_order = self.protein.get_current_order() # Assuming this method exists
        neighbor_order = current_order.copy()

        # Example: Modify one random element in the order
        random_index = random.randint(0, len(neighbor_order) - 1)
        direction_options = [1, -1, 2, -2, 3, -3] # Adjust as per your model
        neighbor_order[random_index] = random.choice(direction_options)

        return neighbor_order
