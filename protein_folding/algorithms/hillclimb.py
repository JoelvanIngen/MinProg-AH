import random
import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from protein_folding.protein import Protein

class HillClimb:
    """
    Hill Climbing algorithm for protein folding, utilizing compute_bond_score and self.directions
    for evaluating and generating protein configurations.
    """

    def __init__(self, protein: 'Protein', dimensions: int, temperature: float = 1000.0, cooling_rate: float = 0.001):
        self.protein = protein
        self.dimensions = dimensions
        self.temperature = temperature
        self.cooling_rate = cooling_rate
    
    def run(self):
        # Initial random configuration
        self.protein.random_configuration()
        best_score = self.protein.compute_bond_score()

        while self.temperature > 1:
            # Generate a neighbor configuration by modifying self.directions
            neighbor = self.protein.generate_neighbor_by_directions()

            # Evaluate the new configuration using compute_bond_score
            neighbor_score = neighbor.compute_bond_score()

            # Calculate probability of accepting the new configuration
            if neighbor_score < best_score:
                accept_probability = 1
            else:
                accept_probability = math.exp((best_score - neighbor_score) / self.temperature)

            # Determine if the neighbor configuration should be accepted
            if random.uniform(0, 1) < accept_probability:
                self.protein = neighbor
                best_score = neighbor_score

            # Cooling down process
            self.temperature *= 1 - self.cooling_rate

        return self.protein
